!> \file advec_kappa.f90
!!  Does advection with a kappa limiter scheme.
!! \par Revision list
!! \par Authors
!! \see Hundsdorfer et al 1995
!!
!! For advection of scalars that need to be strictly monotone (for example chemically reacting species)
!! the kappa scheme has been implemented:
!! \latexonly
!! \begin{eqnarray}
!!  F_{i-\frac{1}{2}}^{\kappa} &=& \fav{u}_{i-\frac{1}{2}}
!!  \left[\phi_{i-1}+\frac{1}{2}\kappa_{i-\frac{1}{2}}\left(\phi_{i-1}-\phi_{i-2}\right)\right],
!! \end{eqnarray}
!! in case $\fav{u}>0$. $\kappa_{i-\smfrac{1}{2}}$ serves as a switch between higher order advection and
!! first order upwind in case of strong upwind gradients of $\phi$.
!! \endlatexonly
!! This makes the scheme monotone, but also rather dissipative.
!!
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

  subroutine advecc_kappa(putin,putout)

  use modglobal, only : i1,i2,ih,j1,j2,jh,k1,kmax,dxi,dyi,dzf
  use modfields, only : u0, v0, w0, rhobf
  implicit none
  real,external :: rlim
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in) :: putin
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout
  real      d1,d2,cf

  integer   i,j,k

  do k=1,kmax
    do j=2,j1
      do i=2,i2
        if (u0(i,j,k)>0) then
          d1 = putin(i-1,j,k)-putin(i-2,j,k)
          d2 = putin(i  ,j,k)-putin(i-1,j,k)
          cf = putin(i-1,j,k)
        else
          d1 = putin(i  ,j,k)-putin(i+1,j,k)
          d2 = putin(i-1,j,k)-putin(i  ,j,k)
          cf = putin(i  ,j,k)
        end if
        cf = cf + rlim(d1,d2)
        putout(i-1,j,k) = putout(i-1,j,k) - cf * u0(i,j,k) * dxi
        putout(i,j,k)   = putout(i,j,k)   + cf * u0(i,j,k) * dxi
      end do
    end do
  end do

  do k=1,kmax
    do j=2,j2
      do i=2,i1
        if (v0(i,j,k)>0) then
          d1 = putin(i,j-1,k)-putin(i,j-2,k)
          d2 = putin(i,j  ,k)-putin(i,j-1,k)
          cf = putin(i,j-1,k)
        else
          d1 = putin(i,j  ,k)-putin(i,j+1,k)
          d2 = putin(i,j-1,k)-putin(i,j  ,k)
          cf = putin(i,j  ,k)
        end if
        cf = cf + rlim(d1,d2)
        putout(i,j-1,k) = putout(i,j-1,k) - cf * v0(i,j,k) * dyi
        putout(i,j,k)   = putout(i,j,k)   + cf * v0(i,j,k) * dyi
      end do
    end do
  end do

  do k=3,kmax
    do j=2,j1
      do i=2,i1
        if (w0(i,j,k)>0) then
          d1 = rhobf(k-1) * putin(i,j,k-1) - rhobf(k-2) * putin(i,j,k-2)
          d2 = rhobf(k)   * putin(i,j,k  ) - rhobf(k-1) * putin(i,j,k-1)
          cf = rhobf(k-1) * putin(i,j,k-1)
        else
          d1 = rhobf(k)   * putin(i,j,k  ) - rhobf(k+1) * putin(i,j,k+1)
          d2 = rhobf(k-1) * putin(i,j,k-1) - rhobf(k)   * putin(i,j,k  )
          cf = rhobf(k)   * putin(i,j,k  )
        end if
        cf = cf + rlim(d1,d2)
        putout(i,j,k-1) = putout(i,j,k-1) - (1./rhobf(k-1))*cf * w0(i,j,k) * (1.0 / dzf(k-1))
        putout(i,j,k)   = putout(i,j,k)   + (1./rhobf(k))*cf   * w0(i,j,k) * (1.0 / dzf(k))
      end do
    end do
  end do

  !special case for vertical advection at k=2
  do j=2,j1
    do i=2,i1
      if (w0(i,j,2)>0) then
        d1 = 0
        d2 = rhobf(2) * putin(i,j,2)   - rhobf(1) * putin(i,j,1)
        cf = rhobf(1) * putin(i,j,1)
      else
        d1 = rhobf(2) * putin(i,j,2)   - rhobf(3) * putin(i,j,3)
        d2 = rhobf(1) * putin(i,j,1)   - rhobf(2) * putin(i,j,2)
        cf = rhobf(2) * putin(i,j,2)
      end if
      cf = cf + rlim(d1,d2)
      putout(i,j,1) = putout(i,j,1) - (1./rhobf(1))*cf * w0(i,j,2) * (1.0 / dzf(1))
      putout(i,j,2) = putout(i,j,2) + (1./rhobf(2))*cf * w0(i,j,2) * (1.0 / dzf(2))
    end do
  end do

end subroutine advecc_kappa

  subroutine advecc_kappa_fast(putin,putout)

  use modglobal, only : i1,i2,ih,j1,j2,jh,k1,kmax,dxi,dyi,dzf
  use modfields, only : u0, v0, w0, rhobf
  implicit none
  real,external :: rlim
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in) :: putin
  real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout
   real      d1,d2,cf

  integer   i,j,k

  real :: d1m, d2m, d1p, cfm, cfp, work

  do k=1,kmax
    do j=2,j1
      do i=2,i2 ! YES
        d2m =  putin(i  ,j,k) -putin(i-1,j,k)
        ! d2p = -putin(i  ,j,k) +putin(i-1,j,k) ! d2p = -d2m

        d1m = putin(i-1,j,k)-putin(i-2,j,k)
        d1p = putin(i  ,j,k)-putin(i+1,j,k)

        cfm = putin(i-1,j,k)
        cfp = putin(i  ,j,k)

        !d1 = (0.5 + sign(0.5, u0(i,j,k))) * d1m + (0.5 - sign(0.5, u0(i,j,k))) * d1p
        !d2 = (0.5 + sign(0.5, u0(i,j,k))) * d2m + (0.5 - sign(0.5, u0(i,j,k))) * d2p
        !d2 = d2m * sign(1.0, u0(i,j,k))
        !cf = (0.5 + sign(0.5, u0(i,j,k))) * cfm + (0.5 - sign(0.5, u0(i,j,k))) * cfp

        if (u0(i,j,k) > 0) then
           d1 = d1m
           d2 = d2m
           cf = cfm
        else
           d1 = d1p
           d2 = -d2m
           cf = cfp
        end if
        
        
        work = cf + &
         min(abs(d1), abs(d2), abs((d1/6.0) + (d2/3.0))) * &
         (sign(0.5, d1) + sign(0.5, d2))

        work = work * u0(i,j,k) * dxi
        putout(i-1,j,k) = putout(i-1,j,k) - work
        putout(i,j,k)   = putout(i,j,k)   + work
      end do
    end do
!  end do

!  do k=1,kmax
    do j=2,j2
      do i=2,i1 ! YES
        !d1m = putin(i,j-1,k)-putin(i,j-2,k)
        !d1p = putin(i,j  ,k)-putin(i,j+1,k)

        !d2m = putin(i,j  ,k)-putin(i,j-1,k)
        ! d2p = putin(i,j-1,k)-putin(i,j  ,k) ! d2p = -d2m

        !cfm = putin(i,j-1,k)
        !cfp = putin(i,j  ,k)

        !d1 = (0.5 + sign(0.5, v0(i,j,k))) * d1m + (0.5 - sign(0.5, v0(i,j,k))) * d1p
        ! d2 = (0.5 + sign(0.5, v0(i,j,k))) * d2m + (0.5 - sign(0.5, v0(i,j,k))) * d2p
        !d2 = d2m * sign(1.0, v0(i,j,k))
        !cf = (0.5 + sign(0.5, v0(i,j,k))) * cfm + (0.5 - sign(0.5, v0(i,j,k))) * cfp
        if (v0(i,j,k) > 0) then
           d1 = putin(i,j-1,k)-putin(i,j-2,k)
           d2 = putin(i,j  ,k)-putin(i,j-1,k)
           cf = putin(i,j-1,k)
        else
           d1 = putin(i,j  ,k)-putin(i,j+1,k)
           d2 = -(putin(i,j  ,k)-putin(i,j-1,k))
           cf = putin(i,j  ,k)
        end if
        
        
        !cf = 0.5 *     v0(i,j,k) *(putin(i,j-1,k)+putin(i,j,k)) &
        !   + 0.5 * abs(v0(i,j,k))*(putin(i,j-1,k)-putin(i,j,k))


        
        work = cf + &
         min(abs(d1), abs(d2), abs((d1/6.0) + (d2/3.0))) * &
         (sign(0.5, d1) + sign(0.5, d2))

        work = work * v0(i,j,k) * dyi
        putout(i,j-1,k) = putout(i,j-1,k) - work
        putout(i,j,k)   = putout(i,j,k)   + work
      end do
    end do
!  end do

!  do k=3,kmax
    if (k >= 3) then
    do j=2,j1
      do i=2,i1 ! YES
        d1m = rhobf(k-1) * putin(i,j,k-1) - rhobf(k-2) * putin(i,j,k-2)
        d2m = rhobf(k)   * putin(i,j,k  ) - rhobf(k-1) * putin(i,j,k-1)

        d1p = rhobf(k)   * putin(i,j,k  ) - rhobf(k+1) * putin(i,j,k+1)
        ! d2p = rhobf(k-1) * putin(i,j,k-1) - rhobf(k)   * putin(i,j,k  ) ! d2p = -d2m

        cfm = rhobf(k-1) * putin(i,j,k-1)
        cfp = rhobf(k)   * putin(i,j,k  )

        !d1 = (0.5 + sign(0.5, w0(i,j,k))) * d1m + (0.5 - sign(0.5, w0(i,j,k))) * d1p
        ! d2 = (0.5 + sign(0.5, w0(i,j,k))) * d2m + (0.5 - sign(0.5, w0(i,j,k))) * d2p
        !d2 = d2m * sign(1.0, w0(i,j,k))
        !cf = (0.5 + sign(0.5, w0(i,j,k))) * cfm + (0.5 - sign(0.5, w0(i,j,k))) * cfp

        if (w0(i,j,k) > 0) then
           d1 = d1m
           d2 = d2m
           cf = cfm
        else
           d1 = d1p
           d2 = -d2m
           cf = cfp
        end if
        
        work = cf + &
         min(abs(d1), abs(d2), abs((d1/6.0) + (d2/3.0))) * &
         (sign(0.5, d1) + sign(0.5, d2))

        work = work * w0(i,j,k)
        putout(i,j,k-1) = putout(i,j,k-1) - (1./(rhobf(k-1)*dzf(k-1)))*work
        putout(i,j,k)   = putout(i,j,k)   + (1./(rhobf(k)  *dzf(k)  ))*work
      end do
   end do
   end if
  end do

  ! from layer 1 to 2, special case. k=2
  do j=2,j1
    do i=2,i1 ! YES
      d1m = 0
      d2m = rhobf(2) * putin(i,j,2) - rhobf(1) * putin(i,j,1)
      cfm = rhobf(1) * putin(i,j,1)
      
      d1p = rhobf(2) * putin(i,j,2) - rhobf(3) * putin(i,j,3)
     !d2p = rhobf(1) * putin(i,j,1) - rhobf(2) * putin(i,j,2  ) ! d2p = -d2m
      cfp = rhobf(2) * putin(i,j,2)

      d1 = (0.5 - sign(0.5, w0(i,j,2))) * d1p
      ! d2 = (0.5 + sign(0.5, w0(i,j,2))) * d2m + (0.5 - sign(0.5, w0(i,j,2))) * d2p
      d2 = d2m * sign(1.0, w0(i,j,2))
      cf = (0.5 + sign(0.5, w0(i,j,2))) * cfm + (0.5 - sign(0.5, w0(i,j,2))) * cfp

      work = cf + &
         min(abs(d1), abs(d2), abs((d1/6.0) + (d2/3.0))) * &
         (sign(0.5, d1) + sign(0.5, d2))

      work = work * w0(i,j,2)
      putout(i,j,1) = putout(i,j,1) - (1./(rhobf(1)*dzf(1)))*work
      putout(i,j,2) = putout(i,j,2) + (1./(rhobf(2)*dzf(2)))*work
    end do
  end do

  end subroutine advecc_kappa_fast

subroutine  halflev_kappa(putin,putout)

  use modglobal, only : i1,ih,j1,jh,k1
    use modfields, only : w0, rhobf, rhobh
    implicit none
    real,external :: rlim
    real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(in) :: putin
    real, dimension(2-ih:i1+ih,2-jh:j1+jh,k1), intent(inout) :: putout
    real      d1,d2,cf
    integer   i,j,k


    do k=3,k1
      do j=2,j1
        do i=2,i1
          if (w0(i,j,k)>=0) then
            d1 = rhobf(k-1) * putin(i,j,k-1) - rhobf(k-2) * putin(i,j,k-2)
            d2 = rhobf(k) * putin(i,j,k  )   - rhobf(k-1) * putin(i,j,k-1)
            cf = rhobf(k-1) * putin(i,j,k-1)
          else
            d1 = rhobf(k)   * putin(i,j,k  ) - rhobf(k+1) * putin(i,j,k+1)
            d2 = rhobf(k-1) * putin(i,j,k-1) - rhobf(k)   * putin(i,j,k  )
            cf = rhobf(k)   * putin(i,j,k  )
          end if
          putout(i,j,k) = (1./rhobh(k))*(cf + rlim(d1,d2))
        end do
      end do
    end do

    do j=2,j1
      do i=2,i1
        if (w0(i,j,2)>=0) then
          d1 = 0
          d2 = rhobf(2) * putin(i,j,2) - rhobf(1) * putin(i,j,1)
          cf = rhobf(1) * putin(i,j,1)
        else
          d1 = rhobf(2) * putin(i,j,2) - rhobf(3) * putin(i,j,3)
          d2 = rhobf(1) * putin(i,j,1) - rhobf(2) * putin(i,j,2)
          cf = rhobf(2) * putin(i,j,2)
        end if
        putout(i,j,2) = (1./rhobh(2))*(cf + rlim(d1,d2))
      end do
    end do

  end subroutine halflev_kappa

!> Determination of the limiter function
  real function rlim(d1,d2)
    use modglobal, only : eps1
    implicit none
    real, intent(in) :: d1 !< Scalar flux at 1.5 cells upwind
    real, intent(in) :: d2 !< Scalar flux at 0.5 cells upwind
    real :: eps2 = 1e-40
    real ri,phir

    ri    = (d2+eps2)/(d1+eps2)
    phir  = max(0.,min(2.*ri,min(1./3.+2./3.*ri,2.)))
    rlim  = 0.5*phir*d1
    end function rlim



