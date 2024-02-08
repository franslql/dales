!!> \file modsubgrid.f90
!!!  Calculates and applies the Sub Filter Scale diffusion
!
!>
!!  Calculates and applies the Sub Filter Scale diffusion
!>
!!  \author Pier Siebesma, K.N.M.I.
!!  \author Stephan de Roode,TU Delft
!!  \author Chiel van Heerwaarden, Wageningen U.R.
!!  \author Thijs Heus,MPI-M
!!  \par Revision list
!!  \todo Documentation
!  This file is part of DALES.
!
! DALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! DALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 1993-2009 Delft University of Technology, Wageningen University, Utrecht University, KNMI
!

module modsubgrid
use modsubgriddata
use modprecision, only: field_r
implicit none
save
  public :: subgrid, initsubgrid, exitsubgrid, subgridnamelist

contains
  subroutine initsubgrid
    use modglobal, only : ih,i1,jh,j1,k1,delta,deltai,dx,dy,zf,dzf,fkar,pi
    use modmpi, only : myid

    implicit none

    integer   :: k

    real :: ceps
    real :: mlen

    call subgridnamelist

    allocate(ekm(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(ekh(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(zlt(2-ih:i1+ih,2-jh:j1+jh,k1))
    allocate(sbdiss(2-ih:i1+ih,2-jh:j1+jh,1))
    allocate(sbshr(2-ih:i1+ih,2-jh:j1+jh,1))
    allocate(sbbuo(2-ih:i1+ih,2-jh:j1+jh,1))
    allocate(csz(k1))
    allocate(anis_fac(k1))

    ! Initialize variables to avoid problems when not using subgrid scheme JvdD
    ! Determination of subgrid constants is explained in De Roode et al. 2017
    ekm=0.; ekh=0.; zlt=0.; sbdiss=0.; sbshr=0.; sbbuo=0.; csz=0.
    anis_fac = 0.

    cm = cf / (2. * pi) * (1.5*alpha_kolm)**(-1.5)

!     ch   = 2. * alpha_kolm / beta_kolm
    ch   = 1.0/Prandtl
    ch2  = ch-ch1

    ceps = 2. * pi / cf * (1.5*alpha_kolm)**(-1.5)
    ce1  = (cn**2)* (cm/Rigc - ch1*cm )
    ce2  = ceps - ce1

    if(cs < 0.) then
      csz(:)  = (cm**3/ceps)**0.25   !< Smagorinsky constant
    else
      csz(:)  = cs
    end if

    if(lmason) then
      do k = 1,k1
        mlen   = (1. / (csz(k) * delta(k))**nmason + 1. / (fkar * zf(k))**nmason)**(-1./nmason)
        csz(k) = mlen / delta(k)
      end do
    end if

    if(lanisotrop) then
       ! Anisotropic diffusion scheme  https://doi.org/10.1029/2022MS003095
       ! length scale in TKE equation is delta z (private communication with Marat)
       if ((dx.ne.dy) .and. myid == 0) stop "The anisotropic diffusion assumes dx=dy."
       do k = 1,k1
          deltai   (k) = 1./dzf(k)       !overrules deltai (k) = 1/delta(k) as defined in initglobal
          anis_fac (k) = (dx/dzf(k))**2  !assumes dx=dy. is used to enhance horizontal diffusion
       end do
    else
       do k = 1,k1
          anis_fac (k) = 1.   !horizontal = vertical diffusion
       end do
    endif


    if (myid==0) then
      write (6,*) 'cf    = ',cf
      write (6,*) 'cm    = ',cm
      write (6,*) 'ch    = ',ch
      write (6,*) 'ch1   = ',ch1
      write (6,*) 'ch2   = ',ch2
      write (6,*) 'ceps  = ',ceps
      write (6,*) 'ceps1 = ',ce1
      write (6,*) 'ceps2 = ',ce2
      write (6,*) 'cs    = ',cs
      write (6,*) 'Rigc  = ',Rigc
    endif

  end subroutine initsubgrid

  subroutine subgridnamelist
    use modglobal, only : ifnamopt,fname_options,checknamelisterror
    use modmpi,    only : myid, comm3d, mpierr, mpi_logical, D_MPI_BCAST

    implicit none

    integer :: ierr

    namelist/NAMSUBGRID/ &
        ldelta,lmason,cf,cn,Rigc,Prandtl,lsmagorinsky,cs,nmason,sgs_surface_fix,ch1,lanisotrop

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMSUBGRID,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMSUBGRID')
      write(6 ,NAMSUBGRID)
      close(ifnamopt)

      if (lmason .and. .not. ldelta) stop "lmason = .true. requires ldelta = .true."
      if (lmason .and. lanisotrop) stop "lmason = .true. is not compatible with lanisotropic = .true."
    end if

    call D_MPI_BCAST(ldelta          ,1, 0,comm3d,mpierr)
    call D_MPI_BCAST(lmason          ,1, 0,comm3d,mpierr)
    call D_MPI_BCAST(nmason          ,1, 0,comm3d,mpierr)
    call D_MPI_BCAST(lsmagorinsky    ,1, 0,comm3d,mpierr)
    call D_MPI_BCAST(lanisotrop      ,1, 0,comm3d,mpierr)
    call D_MPI_BCAST(cs              ,1, 0,comm3d,mpierr)
    call D_MPI_BCAST(cf              ,1, 0,comm3d,mpierr)
    call D_MPI_BCAST(cn              ,1, 0,comm3d,mpierr)
    call D_MPI_BCAST(Rigc            ,1, 0,comm3d,mpierr)
    call D_MPI_BCAST(Prandtl         ,1, 0,comm3d,mpierr)
    call D_MPI_BCAST(sgs_surface_fix ,1, 0,comm3d,mpierr)
    call D_MPI_BCAST(ch1             ,1, 0,comm3d,mpierr)
  end subroutine subgridnamelist

  subroutine subgrid

 ! Diffusion subroutines
! Thijs Heus, Chiel van Heerwaarden, 15 June 2007

    use modglobal, only : nsv, lmoist
    use modfields, only : up,vp,wp,e12p,thl0,thlp,qt0,qtp,sv0,svp
    use modsurfdata,only : thlflux,qtflux,svflux
    implicit none
    integer n

    call closure
    call diffu(up)
    call diffv(vp)
    call diffw(wp)
    if (.not. lsmagorinsky) call diffe(e12p)
    call diffc(thl0,thlp,thlflux)
    if (lmoist) call diffc( qt0, qtp, qtflux)
    do n=1,nsv
      call diffc(sv0(:,:,:,n),svp(:,:,:,n),svflux(:,:,n))
    end do
    if (.not. lsmagorinsky) call sources
  end subroutine

  subroutine exitsubgrid
    implicit none
    deallocate(ekm,ekh,zlt,sbdiss,sbbuo,sbshr,csz)
  end subroutine exitsubgrid

   subroutine closure

!-----------------------------------------------------------------|
!                                                                 |
!*** *closure*  calculates K-coefficients                         |
!                                                                 |
!      Hans Cuijpers   I.M.A.U.   06/01/1995                      |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!     All the K-closure factors are calculated.                   |
!                                                                 |
!     ekm(i,j,k) = k sub m : for velocity-closure                 |
!     ekh(i,j,k) = k sub h : for temperture-closure               |
!     ekh(i,j,k) = k sub h = k sub c : for concentration-closure  |
!                                                                 |
!     We will use the next model for these factors:               |
!                                                                 |
!     k sub m = 0.12 * l * sqrt(E)                                |
!                                                                 |
!     k sub h = k sub c = ( 1 + (2*l)/D ) * k sub m               |
!                                                                 |
!           where : l = mixing length  ( in model = z2 )          |
!                   E = subgrid energy                            |
!                   D = grid-size distance                        |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!             *closure* is called from *program*.                 |
!                                                                 |
!-----------------------------------------------------------------|

  use modglobal,   only : i1,j1,kmax,k1,ih,jh,i2,j2,delta,ekmin,grav,zf,fkar,deltai, &
                          dxi,dyi,dzf,dzfi,dzhi
  use modfields,   only : dthvdz,e120,u0,v0,w0,thvf
  use modsurfdata, only : dudz,dvdz,z0m
  use modmpi,      only : excjs
  implicit none

  real    :: strain2,mlen
  integer :: i,j,k,kp,km,jp,jm

  if(lsmagorinsky) then
    do k = 1,kmax
      mlen        = csz(k) * delta(k)

      do i = 2,i1
        do j = 2,j1

          kp=k+1
          km=k-1
          jp=j+1
          jm=j-1

          if(k == 1) then
            strain2 =  ( &
              ((u0(i+1,j,k)-u0(i,j,k))   *dxi        )**2    + &
              ((v0(i,jp,k)-v0(i,j,k))    *dyi        )**2    + &
              ((w0(i,j,kp)-w0(i,j,k))    *dzfi(k)     )**2    )

            strain2 = strain2 + 0.5_field_r * ( &
              ( 0.25_field_r*(w0(i+1,j,kp)-w0(i-1,j,kp))*dxi + &
              dudz(i,j)   )**2 )

            strain2 = strain2 + 0.125_field_r * ( &
              ((u0(i,jp,k)-u0(i,j,k))     *dyi     + &
              (v0(i,jp,k)-v0(i-1,jp,k))  *dxi        )**2    + &
              ((u0(i,j,k)-u0(i,jm,k))     *dyi     + &
              (v0(i,j,k)-v0(i-1,j,k))    *dxi        )**2    + &
              ((u0(i+1,j,k)-u0(i+1,jm,k)) *dyi     + &
              (v0(i+1,j,k)-v0(i,j,k))    *dxi        )**2    + &
              ((u0(i+1,jp,k)-u0(i+1,j,k)) *dyi     + &
              (v0(i+1,jp,k)-v0(i,jp,k))  *dxi        )**2    )

            strain2 = strain2 + 0.5_field_r * ( &
              ( 0.25_field_r*(w0(i,jp,kp)-w0(i,jm,kp))*dyi + &
              dvdz(i,j)   )**2 )

          else

            strain2 =  ( &
              ((u0(i+1,j,k)-u0(i,j,k))   *dxi        )**2    + &
              ((v0(i,jp,k)-v0(i,j,k))    *dyi        )**2    + &
              ((w0(i,j,kp)-w0(i,j,k))    *dzfi(k)     )**2    )

            strain2 = strain2 + 0.125_field_r * ( &
              ((w0(i,j,kp)-w0(i-1,j,kp))  *dxi     + &
              (u0(i,j,kp)-u0(i,j,k))     * dzhi(kp)  )**2    + &
              ((w0(i,j,k)-w0(i-1,j,k))    *dxi     + &
              (u0(i,j,k)-u0(i,j,km))     * dzhi(k)   )**2    + &
              ((w0(i+1,j,k)-w0(i,j,k))    *dxi     + &
              (u0(i+1,j,k)-u0(i+1,j,km)) * dzhi(k)   )**2    + &
              ((w0(i+1,j,kp)-w0(i,j,kp))  *dxi     + &
              (u0(i+1,j,kp)-u0(i+1,j,k)) * dzhi(kp)  )**2    )

            strain2 = strain2 + 0.125_field_r * ( &
              ((u0(i,jp,k)-u0(i,j,k))     *dyi     + &
              (v0(i,jp,k)-v0(i-1,jp,k))  *dxi        )**2    + &
              ((u0(i,j,k)-u0(i,jm,k))     *dyi     + &
              (v0(i,j,k)-v0(i-1,j,k))    *dxi        )**2    + &
              ((u0(i+1,j,k)-u0(i+1,jm,k)) *dyi     + &
              (v0(i+1,j,k)-v0(i,j,k))    *dxi        )**2    + &
              ((u0(i+1,jp,k)-u0(i+1,j,k)) *dyi     + &
              (v0(i+1,jp,k)-v0(i,jp,k))  *dxi        )**2    )

            strain2 = strain2 + 0.125_field_r * ( &
              ((v0(i,j,kp)-v0(i,j,k))     * dzhi(kp) + &
              (w0(i,j,kp)-w0(i,jm,kp))   *dyi        )**2    + &
              ((v0(i,j,k)-v0(i,j,km))     * dzhi(k)+ &
              (w0(i,j,k)-w0(i,jm,k))     *dyi        )**2    + &
              ((v0(i,jp,k)-v0(i,jp,km))  * dzhi(k)+ &
              (w0(i,jp,k)-w0(i,j,k))     *dyi        )**2    + &
              ((v0(i,jp,kp)-v0(i,jp,k))   * dzhi(kp) + &
              (w0(i,jp,kp)-w0(i,j,kp))   *dyi        )**2    )
          end if

          ekm(i,j,k)  = mlen ** 2 * sqrt(2 * strain2)
          ekh(i,j,k)  = ekm(i,j,k) / Prandtl

          ekm(i,j,k) = max(ekm(i,j,k),ekmin)
          ekh(i,j,k) = max(ekh(i,j,k),ekmin)
        end do
      end do
    end do

  ! do TKE scheme
 else
    ! choose one of ldelta, ldelta+lmason, lanisotropic, or none of them for Deardorff length scale adjustment
    if (ldelta .and. .not. lmason) then
       do k=1,kmax
          do j=2,j1
             do i=2,i1
                zlt(i,j,k) = delta(k)

                ekm(i,j,k) = cm * zlt(i,j,k) * e120(i,j,k)
                ekh(i,j,k) = ch * ekm(i,j,k)

                ekm(i,j,k) = max(ekm(i,j,k),ekmin)
                ekh(i,j,k) = max(ekh(i,j,k),ekmin)
             end do
          end do
       end do
    else if (ldelta .and. lmason) then ! delta scheme with Mason length scale correction
       do k=1,kmax
          do j=2,j1
             do i=2,i1
                zlt(i,j,k) = delta(k)
                zlt(i,j,k) = (1 / zlt(i,j,k) ** nmason + 1 / ( fkar * (zf(k) + z0m(i,j)))**nmason) ** (-1./nmason)

                ekm(i,j,k) = cm * zlt(i,j,k) * e120(i,j,k)
                ekh(i,j,k) = ch * ekm(i,j,k)

                ekm(i,j,k) = max(ekm(i,j,k),ekmin)
                ekh(i,j,k) = max(ekh(i,j,k),ekmin)
             end do
          end do
       end do
    else if (lanisotrop) then ! Anisotropic diffusion,  https://doi.org/10.1029/2022MS003095
       do k=1,kmax
          do j=2,j1
             do i=2,i1
                zlt(i,j,k) = dzf(k)

                ekm(i,j,k) = cm * zlt(i,j,k) * e120(i,j,k)
                ekh(i,j,k) = ch * ekm(i,j,k)

                ekm(i,j,k) = max(ekm(i,j,k),ekmin)
                ekh(i,j,k) = max(ekh(i,j,k),ekmin)
             end do
          end do
       end do
    else ! Deardorff lengthscale correction
       do k=1,kmax
          do j=2,j1
             do i=2,i1
                zlt(i,j,k) = delta(k)

                !original
                !if (dthvdz(i,j,k) > 0) then
                !zlt(i,j,k) = min(delta(k),cn*e120(i,j,k)/sqrt(grav/thvf(k)*abs(dthvdz(i,j,k))))
                !end if

                ! alternative without if
                zlt(i,j,k) = min(delta(k), &
                     cn*e120(i,j,k) / sqrt( grav/thvf(k) * abs(dthvdz(i,j,k))) + &
                     delta(k) * (1-sign(1.0_field_r,dthvdz(i,j,k))))
                ! the final line is 0 if dthvdz(i,j,k) > 0, else 2*delta(k)
                ! ensuring that zlt(i,j,k) = delta(k) when dthvdz < 0, as
                ! in the original scheme.

                ekm(i,j,k) = cm * zlt(i,j,k) * e120(i,j,k)
                ekh(i,j,k) = (ch1 + ch2 * zlt(i,j,k)*deltai(k)) * ekm(i,j,k)

                ekm(i,j,k) = max(ekm(i,j,k),ekmin)
                ekh(i,j,k) = max(ekh(i,j,k),ekmin)
             end do
          end do
       end do
    end if
  end if

!*************************************************************
!     Set cyclic boundary condition for K-closure factors.
!*************************************************************

  call excjs( ekm           , 2,i1,2,j1,1,k1,ih,jh)
  call excjs( ekh           , 2,i1,2,j1,1,k1,ih,jh)

  do j=1,j2
    do i=1,i2
      ekm(i,j,k1)  = ekm(i,j,kmax)
      ekh(i,j,k1)  = ekh(i,j,kmax)
    end do
  end do

  return
  end subroutine closure
  subroutine sources


!-----------------------------------------------------------------|
!                                                                 |
!*** *sources*                                                    |
!      calculates various terms from the subgrid TKE equation     |
!                                                                 |
!     Hans Cuijpers   I.M.A.U.     06/01/1995                     |
!                                                                 |
!     purpose.                                                    |
!     --------                                                    |
!                                                                 |
!      Subroutine sources calculates all other terms in the       |
!      subgrid energy equation, except for the diffusion terms.   |
!      These terms are calculated in subroutine diff.             |
!                                                                 |
!**   interface.                                                  |
!     ----------                                                  |
!                                                                 |
!     *sources* is called from *program*.                         |
!                                                                 |
!-----------------------------------------------------------------|

  use modglobal,   only : i1,j1,kmax,dxi,dyi,dzfi,dzhi,grav,cu,cv,deltai
  use modfields,   only : u0,v0,w0,e120,e12p,dthvdz,thvf
  use modsurfdata,  only : dudz,dvdz,ustar,thlflux
  use modsubgriddata, only: sgs_surface_fix

  implicit none

  real    tdef2, uwflux, vwflux, local_dudz, local_dvdz, local_dthvdz, horv
  integer i,j,k,jm,jp,km,kp


  do k=2,kmax
  do j=2,j1
  do i=2,i1
    kp=k+1
    km=k-1
    jp=j+1
    jm=j-1

    tdef2 = 2 * ( &
             ((u0(i+1,j,k)-u0(i,j,k))   *dxi         )**2    + &
             ((v0(i,jp,k)-v0(i,j,k))    *dyi         )**2    + &
             ((w0(i,j,kp)-w0(i,j,k))    *dzfi(k)     )**2    )

    tdef2 = tdef2 + 0.25_field_r * ( &
              ((w0(i,j,kp)-w0(i-1,j,kp))  * dxi     + &
               (u0(i,j,kp)-u0(i,j,k))     * dzhi(kp)  )**2    + &
              ((w0(i,j,k)-w0(i-1,j,k))    * dxi     + &
               (u0(i,j,k)-u0(i,j,km))     * dzhi(k)   )**2    + &
              ((w0(i+1,j,k)-w0(i,j,k))    * dxi     + &
               (u0(i+1,j,k)-u0(i+1,j,km)) * dzhi(k)   )**2    + &
              ((w0(i+1,j,kp)-w0(i,j,kp))  * dxi     + &
               (u0(i+1,j,kp)-u0(i+1,j,k)) * dzhi(kp)  )**2    )

    tdef2 = tdef2 + 0.25_field_r * ( &
              ((u0(i,jp,k)-u0(i,j,k))     * dyi     + &
               (v0(i,jp,k)-v0(i-1,jp,k))  * dxi        )**2    + &
              ((u0(i,j,k)-u0(i,jm,k))     * dyi     + &
               (v0(i,j,k)-v0(i-1,j,k))    * dxi        )**2    + &
              ((u0(i+1,j,k)-u0(i+1,jm,k)) * dyi     + &
               (v0(i+1,j,k)-v0(i,j,k))    * dxi        )**2    + &
              ((u0(i+1,jp,k)-u0(i+1,j,k)) * dyi     + &
               (v0(i+1,jp,k)-v0(i,jp,k))  * dxi        )**2    )

    tdef2 = tdef2 + 0.25_field_r * ( &
              ((v0(i,j,kp)-v0(i,j,k))     * dzhi(kp) + &
               (w0(i,j,kp)-w0(i,jm,kp))   * dyi        )**2    + &
              ((v0(i,j,k)-v0(i,j,km))     * dzhi(k)+ &
               (w0(i,j,k)-w0(i,jm,k))     * dyi        )**2    + &
              ((v0(i,jp,k)-v0(i,jp,km))   * dzhi(k)+ &
               (w0(i,jp,k)-w0(i,j,k))     * dyi        )**2    + &
              ((v0(i,jp,kp)-v0(i,jp,k))   * dzhi(kp) + &
               (w0(i,jp,kp)-w0(i,j,kp))   * dyi        )**2    )

    e12p(i,j,k) = e12p(i,j,k) &
                + (ekm(i,j,k)*tdef2 - ekh(i,j,k)*grav/thvf(k)*dthvdz(i,j,k) ) / (2*e120(i,j,k)) &  !  sbshr and sbbuo
                - (ce1 + ce2*zlt(i,j,k)*deltai(k)) * e120(i,j,k)**2 /(2*zlt(i,j,k))               !  sbdiss

  end do
  end do
  end do
!     ----------------------------------------------end i,j,k-loop

!     --------------------------------------------
!     special treatment for lowest full level: k=1
!     --------------------------------------------


  do j=2,j1
    jp=j+1
    jm=j-1
  do i=2,i1



! **  Calculate "shear" production term: tdef2  ****************

    tdef2 =  2 * ( &
            ((u0(i+1,j,1)-u0(i,j,1))*dxi)**2 &
          + ((v0(i,jp,1)-v0(i,j,1))*dyi)**2 &
          + ((w0(i,j,2)-w0(i,j,1))*dzfi(1))**2   )

    if (sgs_surface_fix) then
          ! Use known surface flux and exchange coefficient to derive
          ! consistent gradient (such that correct flux will occur in
          ! shear production term)
          ! Make sure that no division by zero occurs in determination of the
          ! directional component; ekm should already be >= ekmin
          ! Replace the dudz by surface flux -uw / ekm
          horv = max(sqrt((u0(i,j,1)+cu)**2+(v0(i,j,1)+cv)**2),  0.01)
          uwflux = -ustar(i,j)*ustar(i,j)* ((u0(i,j,1)+cu)/horv)
          local_dudz = -uwflux / ekm(i,j,1)
          tdef2 = tdef2 + ( 0.25_field_r*(w0(i+1,j,2)-w0(i-1,j,2))*dxi + &
               local_dudz )**2
    else
          tdef2 = tdef2 + ( 0.25_field_r*(w0(i+1,j,2)-w0(i-1,j,2))*dxi + &
                                  dudz(i,j)   )**2
    endif

    tdef2 = tdef2 +   0.25_field_r *( &
          ((u0(i,jp,1)-u0(i,j,1))*dyi+(v0(i,jp,1)-v0(i-1,jp,1))*dxi)**2 &
         +((u0(i,j,1)-u0(i,jm,1))*dyi+(v0(i,j,1)-v0(i-1,j,1))*dxi)**2 &
         +((u0(i+1,j,1)-u0(i+1,jm,1))*dyi+(v0(i+1,j,1)-v0(i,j,1))*dxi)**2 &
         +((u0(i+1,jp,1)-u0(i+1,j,1))*dyi+ &
                                 (v0(i+1,jp,1)-v0(i,jp,1))*dxi)**2   )

    if (sgs_surface_fix) then
          ! Use known surface flux and exchange coefficient to derive
          ! consistent gradient (such that correct flux will occur in
          ! shear production term)
          ! Make sure that no division by zero occurs in determination of the
          ! directional component; ekm should already be >= ekmin
          ! Replace the dvdz by surface flux -vw / ekm
          horv = max(sqrt((u0(i,j,1)+cu)**2+(v0(i,j,1)+cv)**2),  0.01)
          vwflux = -ustar(i,j)*ustar(i,j)* ((v0(i,j,1)+cv)/horv)
          local_dvdz = -vwflux / ekm(i,j,1)
          tdef2 = tdef2 + ( 0.25_field_r*(w0(i,jp,2)-w0(i,jm,2))*dyi + &
                        local_dvdz  )**2
    else
         tdef2 = tdef2 + ( 0.25_field_r*(w0(i,jp,2)-w0(i,jm,2))*dyi + &
                                dvdz(i,j)   )**2
    endif

! **  Include shear and buoyancy production terms and dissipation **

    sbshr(i,j,1)  = ekm(i,j,1)*tdef2/ ( 2*e120(i,j,1))
    if (sgs_surface_fix) then
          ! Replace the -ekh *  dthvdz by the surface flux of thv
          ! (but we only have the thlflux , which seems at the surface to be
          ! equivalent
          local_dthvdz = -thlflux(i,j)/ekh(i,j,1)
          sbbuo(i,j,1)  = -ekh(i,j,1)*grav/thvf(1)*local_dthvdz/ ( 2*e120(i,j,1))
    else
          sbbuo(i,j,1)  = -ekh(i,j,1)*grav/thvf(1)*dthvdz(i,j,1)/ ( 2*e120(i,j,1))
    endif
    sbdiss(i,j,1) = - (ce1 + ce2*zlt(i,j,1)*deltai(1)) * e120(i,j,1)**2 /(2*zlt(i,j,1))
  end do
  end do

  e12p(2:i1,2:j1,1) = e12p(2:i1,2:j1,1) + &
            sbshr(2:i1,2:j1,1)+sbbuo(2:i1,2:j1,1)+sbdiss(2:i1,2:j1,1)

  return
  end subroutine sources

  subroutine diffc (a_in,a_out,flux)

    use modglobal, only : i1,ih,i2,j1,jh,j2,k1,kmax,dx2i,dzf,dzfi,dy2i,dzhi
    use modfields, only : rhobf,rhobh
    implicit none

    real(field_r), intent(in)    :: a_in(2-ih:i1+ih,2-jh:j1+jh,k1)
    real(field_r), intent(inout) :: a_out(2-ih:i1+ih,2-jh:j1+jh,k1)
    real, intent(in)    :: flux (i2,j2)

    integer i,j,k,jm,jp,km,kp

    do k=2,kmax
      kp=k+1
      km=k-1

      do j=2,j1
        jp=j+1
        jm=j-1

        do i=2,i1
          a_out(i,j,k) = a_out(i,j,k) &
                    +  0.5_field_r * ( &
                  ( (ekh(i+1,j,k)+ekh(i,j,k))*(a_in(i+1,j,k)-a_in(i,j,k)) &
                    -(ekh(i,j,k)+ekh(i-1,j,k))*(a_in(i,j,k)-a_in(i-1,j,k)))*dx2i * anis_fac(k) &
                    + &
                  ( (ekh(i,jp,k)+ekh(i,j,k)) *(a_in(i,jp,k)-a_in(i,j,k)) &
                    -(ekh(i,j,k)+ekh(i,jm,k)) *(a_in(i,j,k)-a_in(i,jm,k)) )*dy2i * anis_fac(k) &
                  + &
                  ( rhobh(kp)/rhobf(k) * (dzf(kp)*ekh(i,j,k) + dzf(k)*ekh(i,j,kp)) &
                    *  (a_in(i,j,kp)-a_in(i,j,k)) * dzhi(kp)**2 &
                    - &
                    rhobh(k)/rhobf(k) * (dzf(km)*ekh(i,j,k) + dzf(k)*ekh(i,j,km)) &
                    *  (a_in(i,j,k)-a_in(i,j,km)) * dzhi(k)**2           )*dzfi(k) &
                            )

        end do
      end do
    end do

    do j=2,j1
      do i=2,i1

        a_out(i,j,1) = a_out(i,j,1) &
                  + 0.5_field_r * ( &
                ( (ekh(i+1,j,1)+ekh(i,j,1))*(a_in(i+1,j,1)-a_in(i,j,1)) &
                  -(ekh(i,j,1)+ekh(i-1,j,1))*(a_in(i,j,1)-a_in(i-1,j,1)) )*dx2i * anis_fac(1) &
                  + &
                ( (ekh(i,j+1,1)+ekh(i,j,1))*(a_in(i,j+1,1)-a_in(i,j,1)) &
                  -(ekh(i,j,1)+ekh(i,j-1,1))*(a_in(i,j,1)-a_in(i,j-1,1)) )*dy2i * anis_fac(1) &
                  + &
                ( rhobh(2)/rhobf(1) * (dzf(2)*ekh(i,j,1) + dzf(1)*ekh(i,j,2)) &
                  *  (a_in(i,j,2)-a_in(i,j,1)) * dzhi(2)**2 &
                  + rhobh(1)/rhobf(1)*flux(i,j) *2                        )*dzfi(1) &
                          )

      end do
    end do

  end subroutine diffc



  subroutine diffe(a_out)

    use modglobal, only : i1,ih,j1,jh,k1,kmax,dx2i,dzf,dzfi,dy2i,dzhi
    use modfields, only : e120,rhobf,rhobh
    implicit none

    real(field_r), intent(inout) :: a_out(2-ih:i1+ih,2-jh:j1+jh,k1)
    integer             :: i,j,k,jm,jp,km,kp

    do k=2,kmax
      kp=k+1
      km=k-1

      do j=2,j1
        jp=j+1
        jm=j-1

        do i=2,i1

          a_out(i,j,k) = a_out(i,j,k) &
                  +  ( &
              ((ekm(i+1,j,k)+ekm(i,j,k))*(e120(i+1,j,k)-e120(i,j,k)) &
              -(ekm(i,j,k)+ekm(i-1,j,k))*(e120(i,j,k)-e120(i-1,j,k)))*dx2i * anis_fac(k) &
                  + &
              ((ekm(i,jp,k)+ekm(i,j,k)) *(e120(i,jp,k)-e120(i,j,k)) &
              -(ekm(i,j,k)+ekm(i,jm,k)) *(e120(i,j,k)-e120(i,jm,k)) )*dy2i * anis_fac(k) &
                  + &
              (rhobh(kp)/rhobf(k) * (dzf(kp)*ekm(i,j,k) + dzf(k)*ekm(i,j,kp)) &
              *(e120(i,j,kp)-e120(i,j,k)) * dzhi(kp)**2 &
              - rhobh(k)/rhobf(k) * (dzf(km)*ekm(i,j,k) + dzf(k)*ekm(i,j,km)) &
              *(e120(i,j,k)-e120(i,j,km)) * dzhi(k)**2        )*dzfi(k) &
                            )

        end do
      end do
    end do

  !     --------------------------------------------
  !     special treatment for lowest full level: k=1
  !     --------------------------------------------

    do j=2,j1
      do i=2,i1

        a_out(i,j,1) = a_out(i,j,1) + &
            ( (ekm(i+1,j,1)+ekm(i,j,1))*(e120(i+1,j,1)-e120(i,j,1)) &
              -(ekm(i,j,1)+ekm(i-1,j,1))*(e120(i,j,1)-e120(i-1,j,1)) )*dx2i * anis_fac(1) &
            + &
            ( (ekm(i,j+1,1)+ekm(i,j,1))*(e120(i,j+1,1)-e120(i,j,1)) &
              -(ekm(i,j,1)+ekm(i,j-1,1))*(e120(i,j,1)-e120(i,j-1,1)) )*dy2i * anis_fac(1) &
            + &
              ( rhobh(2)/rhobf(1) * (dzf(2)*ekm(i,j,1) + dzf(1)*ekm(i,j,2)) &
              *  (e120(i,j,2)-e120(i,j,1)) * dzhi(2)**2              )*dzfi(1)

      end do
    end do

  end subroutine diffe


  subroutine diffu (a_out)

    use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi,dx2i,dzf,dzfi,dyi,dzhi,cu,cv
    use modfields, only : u0,v0,w0,rhobf,rhobh
    use modsurfdata,only : ustar
    implicit none

    real(field_r), intent(inout) :: a_out(2-ih:i1+ih,2-jh:j1+jh,k1)
    real(field_r)                :: emmo,emom,emop,empo
    real(field_r)                :: fu
    real(field_r)                :: ucu, upcu
    integer             :: i,j,k,jm,jp,km,kp

    do k=2,kmax
      kp=k+1
      km=k-1

      do j=2,j1
        jp=j+1
        jm=j-1

        do i=2,i1

          emom = ( dzf(km) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i-1,j,km) ) ) * &
                    ( .25_field_r   * dzhi(k) )

          emop = ( dzf(kp) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
                      dzf(k)  * ( ekm(i,j,kp) + ekm(i-1,j,kp) ) ) * &
                    ( .25_field_r   * dzhi(kp) )

          empo = 0.25_field_r * ( &
                  ekm(i,j,k)+ekm(i,jp,k)+ekm(i-1,jp,k)+ekm(i-1,j,k)  )

          emmo = 0.25_field_r * ( &
                  ekm(i,j,k)+ekm(i,jm,k)+ekm(i-1,jm,k)+ekm(i-1,j,k)  )


          a_out(i,j,k) = a_out(i,j,k) &
                  + &
                  ( ekm(i,j,k)  * (u0(i+1,j,k)-u0(i,j,k)) &
                    -ekm(i-1,j,k)* (u0(i,j,k)-u0(i-1,j,k)) ) * 2 * dx2i * anis_fac(k) &
                  + &
                  ( empo * ( (u0(i,jp,k)-u0(i,j,k))   *dyi &
                            +(v0(i,jp,k)-v0(i-1,jp,k))*dxi) &
                    -emmo * ( (u0(i,j,k)-u0(i,jm,k))   *dyi &
                            +(v0(i,j,k)-v0(i-1,j,k))  *dxi)   ) * dyi * anis_fac(k) &
                  + &
                  ( rhobh(kp)/rhobf(k) * emop * ( (u0(i,j,kp)-u0(i,j,k))   *dzhi(kp) &
                            +(w0(i,j,kp)-w0(i-1,j,kp))*dxi) &
                    - rhobh(k)/rhobf(k) * emom * ( (u0(i,j,k)-u0(i,j,km))   *dzhi(k) &
                            +(w0(i,j,k)-w0(i-1,j,k))  *dxi)   ) *dzfi(k)

        end do
      end do
    end do

  !     --------------------------------------------
  !     special treatment for lowest full level: k=1
  !     --------------------------------------------

    do j=2,j1
      jp = j+1
      jm = j-1

      do i=2,i1

        empo = 0.25_field_r * ( &
              ekm(i,j,1)+ekm(i,jp,1)+ekm(i-1,jp,1)+ekm(i-1,j,1)  )

        emmo = 0.25_field_r * ( &
              ekm(i,j,1)+ekm(i,jm,1)+ekm(i-1,jm,1)+ekm(i-1,j,1)  )

        emop = ( dzf(2) * ( ekm(i,j,1) + ekm(i-1,j,1) )  + &
                    dzf(1) * ( ekm(i,j,2) + ekm(i-1,j,2) ) ) * &
                  ( .25_field_r   * dzhi(2) )


        ucu   = 0.5_field_r*(u0(i,j,1)+u0(i+1,j,1))+cu

        if(ucu >= 0.) then
          upcu  = max(ucu,1.e-10)
        else
          upcu  = min(ucu,-1.e-10)
        end if


        fu = ( 0.5_field_r*( ustar(i,j)+ustar(i-1,j) ) )**2  * &
                upcu/sqrt(upcu**2  + &
                ((v0(i,j,1)+v0(i-1,j,1)+v0(i,jp,1)+v0(i-1,jp,1))/4+cv)**2)

        a_out(i,j,1) = a_out(i,j,1) &
                + &
              ( ekm(i,j,1)  * (u0(i+1,j,1)-u0(i,j,1)) &
              -ekm(i-1,j,1)* (u0(i,j,1)-u0(i-1,j,1)) ) * 2 * dx2i * anis_fac(1) &
                + &
              ( empo * ( (u0(i,jp,1)-u0(i,j,1))   *dyi &
                        +(v0(i,jp,1)-v0(i-1,jp,1))*dxi) &
              -emmo * ( (u0(i,j,1)-u0(i,jm,1))   *dyi &
                        +(v0(i,j,1)-v0(i-1,j,1))  *dxi)   ) * dyi * anis_fac(1) &
               + &
              ( rhobh(2)/rhobf(1) * emop * ( (u0(i,j,2)-u0(i,j,1))    *dzhi(2) &
                        +(w0(i,j,2)-w0(i-1,j,2))  *dxi) &
                -rhobh(1)/rhobf(1)*fu   ) * dzfi(1)

      end do
    end do

  end subroutine diffu


  subroutine diffv (a_out)

    use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi,dzf,dzfi,dyi,dy2i,dzhi,cu,cv
    use modfields, only : u0,v0,w0,rhobf,rhobh
    use modsurfdata,only : ustar

    implicit none

    real(field_r), intent(inout) :: a_out(2-ih:i1+ih,2-jh:j1+jh,k1)
    real(field_r)                :: emmo, eomm,eomp,epmo
    real(field_r)                :: fv, vcv,vpcv
    integer             :: i,j,k,jm,jp,km,kp

    do k=2,kmax
      kp=k+1
      km=k-1

      do j=2,j1
        jp=j+1
        jm=j-1

        do i=2,i1

          eomm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i,jm,km) ) ) * &
                    ( .25_field_r   * dzhi(k) )

          eomp = ( dzf(kp) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
                      dzf(k)  * ( ekm(i,j,kp) + ekm(i,jm,kp) ) ) * &
                    ( .25_field_r   * dzhi(kp) )

          emmo = 0.25_field_r  * ( &
                ekm(i,j,k)+ekm(i,jm,k)+ekm(i-1,jm,k)+ekm(i-1,j,k)  )

          epmo = 0.25_field_r  * ( &
                ekm(i,j,k)+ekm(i,jm,k)+ekm(i+1,jm,k)+ekm(i+1,j,k)  )


        a_out(i,j,k) = a_out(i,j,k) &
                + &
              ( epmo * ( (v0(i+1,j,k)-v0(i,j,k))   *dxi &
                        +(u0(i+1,j,k)-u0(i+1,jm,k))*dyi) &
                -emmo * ( (v0(i,j,k)-v0(i-1,j,k))   *dxi &
                        +(u0(i,j,k)-u0(i,jm,k))    *dyi)   ) * dxi * anis_fac(k) &
                + &
              (ekm(i,j,k) * (v0(i,jp,k)-v0(i,j,k)) &
              -ekm(i,jm,k)* (v0(i,j,k)-v0(i,jm,k))  ) * 2 * dy2i * anis_fac(k) &
                + &
              ( rhobh(kp)/rhobf(k) * eomp * ( (v0(i,j,kp)-v0(i,j,k))    *dzhi(kp) &
                        +(w0(i,j,kp)-w0(i,jm,kp))  *dyi) &
                - rhobh(k)/rhobf(k) * eomm * ( (v0(i,j,k)-v0(i,j,km))    *dzhi(k) &
                        +(w0(i,j,k)-w0(i,jm,k))    *dyi)   ) * dzfi(k)

        end do
      end do
    end do

  !     --------------------------------------------
  !     special treatment for lowest full level: k=1
  !     --------------------------------------------

    do j=2,j1
      jp = j+1
      jm = j-1
      do i=2,i1

        emmo = 0.25_field_r * ( &
              ekm(i,j,1)+ekm(i,jm,1)+ekm(i-1,jm,1)+ekm(i-1,j,1)  )

        epmo = 0.25_field_r  * ( &
              ekm(i,j,1)+ekm(i,jm,1)+ekm(i+1,jm,1)+ekm(i+1,j,1)  )

        eomp = ( dzf(2) * ( ekm(i,j,1) + ekm(i,jm,1)  )  + &
                    dzf(1) * ( ekm(i,j,2) + ekm(i,jm,2) ) ) * &
                  ( .25_field_r   * dzhi(2) )

        vcv   = 0.5_field_r*(v0(i,j,1)+v0(i,j+1,1))+cv
        if(vcv >= 0.) then
          vpcv  = max(vcv,1.e-10)
        else
          vpcv  = min(vcv,-1.e-10)
        end if


        fv    = ( 0.5_field_r*( ustar(i,j)+ustar(i,j-1) ) )**2  * &
                    vpcv/sqrt(vpcv**2  + &
                ((u0(i,j,1)+u0(i+1,j,1)+u0(i,jm,1)+u0(i+1,jm,1))/4+cu)**2)

        a_out(i,j,1) = a_out(i,j,1) &
                  + &
                  ( epmo * ( (v0(i+1,j,1)-v0(i,j,1))   *dxi &
                            +(u0(i+1,j,1)-u0(i+1,jm,1))*dyi) &
                    -emmo * ( (v0(i,j,1)-v0(i-1,j,1))   *dxi &
                            +(u0(i,j,1)-u0(i,jm,1))    *dyi)   ) * dxi * anis_fac(1) &
                  + &
                ( ekm(i,j,1) * (v0(i,jp,1)-v0(i,j,1)) &
                  -ekm(i,jm,1)* (v0(i,j,1)-v0(i,jm,1))  ) * 2 * dy2i * anis_fac(1) &
                  + &
                ( rhobh(2)/rhobf(1) * eomp * ( (v0(i,j,2)-v0(i,j,1))     *dzhi(2) &
                          +(w0(i,j,2)-w0(i,jm,2))    *dyi) &
                  -rhobh(1)/rhobf(1)*fv   ) * dzfi(1)

      end do
    end do

  end subroutine diffv



  subroutine diffw(a_out)

    use modglobal, only : i1,ih,j1,jh,k1,kmax,dxi,dyi,dzf,dzfi,dzhi
    use modfields, only : u0,v0,w0,rhobh,rhobf
    implicit none

  !*****************************************************************

    real(field_r), intent(inout) :: a_out(2-ih:i1+ih,2-jh:j1+jh,k1)
    real(field_r)                :: emom, eomm, eopm, epom
    integer             :: i,j,k,jm,jp,km,kp

    do k=2,kmax
      kp=k+1
      km=k-1
      do j=2,j1
        jp=j+1
        jm=j-1
        do i=2,i1

          emom = ( dzf(km) * ( ekm(i,j,k)  + ekm(i-1,j,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i-1,j,km) ) ) * &
                    ( .25_field_r   * dzhi(k) )

          eomm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jm,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i,jm,km) ) ) * &
                    ( .25_field_r   * dzhi(k) )

          eopm = ( dzf(km) * ( ekm(i,j,k)  + ekm(i,jp,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i,jp,km) ) ) * &
                    ( .25_field_r   * dzhi(k) )

          epom = ( dzf(km) * ( ekm(i,j,k)  + ekm(i+1,j,k)  )  + &
                      dzf(k)  * ( ekm(i,j,km) + ekm(i+1,j,km) ) ) * &
                    ( .25_field_r   * dzhi(k) )


          a_out(i,j,k) = a_out(i,j,k) &
                + &
                  ( epom * ( (w0(i+1,j,k)-w0(i,j,k))    *dxi &
                            +(u0(i+1,j,k)-u0(i+1,j,km)) *dzhi(k) ) &
                    -emom * ( (w0(i,j,k)-w0(i-1,j,k))    *dxi &
                            +(u0(i,j,k)-u0(i,j,km))     *dzhi(k) ))*dxi * anis_fac(k) &
                + &
                  ( eopm * ( (w0(i,jp,k)-w0(i,j,k))     *dyi &
                            +(v0(i,jp,k)-v0(i,jp,km))   *dzhi(k) ) &
                    -eomm * ( (w0(i,j,k)-w0(i,jm,k))     *dyi &
                            +(v0(i,j,k)-v0(i,j,km))     *dzhi(k) ))*dyi * anis_fac(k) &
                + (1/rhobh(k))*&
                  ( rhobf(k) * ekm(i,j,k) * (w0(i,j,kp)-w0(i,j,k)) *dzfi(k) &
                  - rhobf(km) * ekm(i,j,km)* (w0(i,j,k)-w0(i,j,km)) *dzfi(km) ) * 2 &
                                                              * dzhi(k)

        end do
      end do
    end do

  end subroutine diffw

end module
