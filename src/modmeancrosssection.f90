!> \file modlsmcrosssection.f90
!! Calculates averages and first order statistics over the reduced dimension and time.
!>
!! Calculates averages and first order statistics over the reduced dimension and time.
!>
!!  \par Revision list
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
module modmeancrosssection
use netcdf
use modglobal, only : longint
implicit none
private
PUBLIC :: initmeancrosssection, exitmeancrosssection, meancrosssection
character(80) :: fnamexz = 'meancrossxz.xxx.xxx.nc', fnameyz = 'meancrossyz.xxx.xxx.nc'
logical :: lmeancross = .false., lreducex = .false., lreducey = .true.
real :: dtav,timeav
integer(kind=longint) :: idtav,itimeav,tnext,tnextwrite
integer :: ntime,nsamples
integer, dimension(5) :: dimid
integer :: ncidxz, ncidyz, varid
real, allocatable, dimension(:,:) :: usum,vsum,wsum,thlsum,qtsum,u2sum,v2sum,w2sum,&
  thl2sum,qt2sum,uwsum,vwsum,thlwsum,qtwsum,uxz,vxz,wxz,thlxz,qtxz,u2xz,v2xz,w2xz,&
  thl2xz,qt2xz,uwxz,vwxz,thlwxz,qtwxz,uxz_int,vxz_int,wxz_int,thlxz_int,qtxz_int,&
  u2xz_int,v2xz_int,w2xz_int,thl2xz_int,qt2xz_int,uwxz_int,vwxz_int,thlwxz_int,qtwxz_int,&
  usum_glob,vsum_glob,wsum_glob,thlsum_glob,qtsum_glob,u2sum_glob,v2sum_glob,w2sum_glob,&
  thl2sum_glob,qt2sum_glob,uwsum_glob,vwsum_glob,thlwsum_glob,qtwsum_glob
contains
  subroutine initmeancrosssection
    use modstat_nc, only : nc_fillvalue
    use modglobal, only : imax,i1,i2,kmax,k1,lboundary,lopenbc,dtav_glob,timeav_glob,&
      btime,dt_lim,zf,zh,cexpnr,fname_options,dx,dy,ifnamopt,tres,checknamelisterror
    use modmpi, only : myid,cmyidx,cmyidy,comm3d,mpierr,MPI_LOGICAL,MY_REAL,myidx,myidy
    implicit none
    integer :: ierr,i
    namelist/NAMMEANCROSSSECTION/ &
    lmeancross, dtav, timeav

    dtav=dtav_glob;timeav=timeav_glob

    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMMEANCROSSSECTION,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMMEANCROSSSECTION')
      write(6 ,NAMMEANCROSSSECTION)
      close(ifnamopt)
    end if

    call MPI_BCAST(timeav     ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(dtav       ,1,MY_REAL    ,0,comm3d,mpierr)
    call MPI_BCAST(lmeancross ,1,MPI_LOGICAL,0,comm3d,mpierr)

    if(.not.lmeancross) return
    idtav = dtav/tres
    itimeav = timeav/tres

    tnext      = idtav   +btime
    tnextwrite = itimeav +btime
    nsamples = itimeav/idtav
    dt_lim = min(dt_lim,tnext)

    if (abs(timeav/dtav-nsamples)>1e-4) then
      stop 'timeav must be a integer multiple of dtav'
    end if

    allocate(usum(2:i2,k1), vsum(2:i1,k1), wsum(1:i2,k1), thlsum(2:i1,k1), qtsum(2:i1,k1), &
            u2sum(2:i2,k1),v2sum(2:i1,k1),w2sum(1:i2,k1),thl2sum(2:i1,k1),qt2sum(2:i1,k1), &
            uwsum(2:i2,k1),vwsum(2:i1,k1),               thlwsum(2:i1,k1),qtwsum(2:i1,k1) )
    allocate(uxz(imax+merge(1,0,lopenbc.and.lboundary(2)),kmax), vxz(imax,kmax), wxz(imax,k1), thlxz(imax,kmax), qtxz(imax,kmax), &
            u2xz(imax+merge(1,0,lopenbc.and.lboundary(2)),kmax),v2xz(imax,kmax),w2xz(imax,k1),thl2xz(imax,kmax),qt2xz(imax,kmax), &
            uwxz(imax+merge(1,0,lopenbc.and.lboundary(2)),k1)  ,vwxz(imax,k1)  ,              thlwxz(imax,k1)  ,qtwxz(imax,k1) )
    allocate(uxz_int(imax+merge(1,0,lopenbc.and.lboundary(2)),kmax), vxz_int(imax,kmax), wxz_int(imax,kmax), thlxz_int(imax,kmax), qtxz_int(imax,kmax),&
            u2xz_int(imax+merge(1,0,lopenbc.and.lboundary(2)),kmax),v2xz_int(imax,kmax),w2xz_int(imax,kmax),thl2xz_int(imax,kmax),qt2xz_int(imax,kmax),&
            uwxz_int(imax+merge(1,0,lopenbc.and.lboundary(2)),k1),  vwxz_int(imax,k1) ,                    thlwxz_int(imax,k1)  ,qtwxz_int(imax,k1) )
    if(myidy==0) then
      allocate(usum_glob(2:i2,k1), vsum_glob(2:i1,k1), wsum_glob(1:i2,k1), thlsum_glob(2:i1,k1), qtsum_glob(2:i1,k1), &
              u2sum_glob(2:i2,k1),v2sum_glob(2:i1,k1),w2sum_glob(1:i2,k1),thl2sum_glob(2:i1,k1),qt2sum_glob(2:i1,k1), &
              uwsum_glob(2:i2,k1),vwsum_glob(2:i1,k1),                    thlwsum_glob(2:i1,k1),qtwsum_glob(2:i1,k1) )
    endif


    ntime = 1
    if(myidy==0.and.lreducey) then
      fnamexz(13:15) = cmyidx
      fnamexz(17:19) = cexpnr
      ! Create netcdf file
      call check( nf90_create(fnamexz, NF90_NETCDF4, ncidxz) )
      ! Specify dimensions
      call check( nf90_def_dim(ncidxz, 'time', NF90_UNLIMITED, dimid(1)) )
      call check( nf90_def_dim(ncidxz, 'zt', kmax, dimid(2)) )
      call check( nf90_def_dim(ncidxz, 'xt', imax, dimid(3)) )
      call check( nf90_def_dim(ncidxz, 'zm', k1, dimid(4)) )
      call check( nf90_def_dim(ncidxz, 'xm', imax+merge(1,0,lopenbc.and.lboundary(2)), dimid(5)) )
      ! Create variables
      call check( nf90_def_var(ncidxz, 'time', NF90_REAL, dimid(1), varid ) )
      call check( nf90_put_att(ncidxz,varid,'longname','Time') )
      call check( nf90_put_att(ncidxz,varid,'units','s') )
      call check( nf90_def_var(ncidxz, 'zt', NF90_REAL, dimid(2), varid ) )
      call check( nf90_put_att(ncidxz,varid,'longname','Vertical displacement of cell centers') )
      call check( nf90_put_att(ncidxz,varid,'units','m') )
      call check( nf90_def_var(ncidxz, 'xt', NF90_REAL, dimid(3), varid ) )
      call check( nf90_put_att(ncidxz,varid,'longname','West-East displacement of cell centers') )
      call check( nf90_put_att(ncidxz,varid,'units','m') )
      call check( nf90_def_var(ncidxz, 'zm', NF90_REAL, dimid(4), varid ) )
      call check( nf90_put_att(ncidxz,varid,'longname','Vertical displacement of cell edges') )
      call check( nf90_put_att(ncidxz,varid,'units','m') )
      call check( nf90_def_var(ncidxz, 'xm', NF90_REAL, dimid(5), varid ) )
      call check( nf90_put_att(ncidxz,varid,'longname','West-East displacement of cell edges') )
      call check( nf90_put_att(ncidxz,varid,'units','m') )
      call check( nf90_def_var(ncidxz, 'uxz', NF90_REAL, (/dimid(5),dimid(2),dimid(1)/), varid ) )
      call check( nf90_put_att(ncidxz,varid,'longname','y averaged west-east velocity') )
      call check( nf90_put_att(ncidxz,varid,'units','m/s') )
      call check( nf90_def_var_deflate(ncidxz,varid, 0, 1, deflate_level = 2) )
      call check( nf90_put_att(ncidxz, varid, '_FillValue',nc_fillvalue) )
      call check( nf90_def_var(ncidxz, 'vxz', NF90_REAL, (/dimid(3),dimid(2),dimid(1)/), varid ) )
      call check( nf90_put_att(ncidxz,varid,'longname','y averaged south-north velocity') )
      call check( nf90_put_att(ncidxz,varid,'units','m/s') )
      call check( nf90_def_var_deflate(ncidxz,varid, 0, 1, deflate_level = 2) )
      call check( nf90_put_att(ncidxz, varid, '_FillValue',nc_fillvalue) )
      call check( nf90_def_var(ncidxz, 'wxz', NF90_REAL, (/dimid(3),dimid(4),dimid(1)/), varid ) )
      call check( nf90_put_att(ncidxz,varid,'longname','y averaged vertical velocity') )
      call check( nf90_put_att(ncidxz,varid,'units','m/s') )
      call check( nf90_def_var_deflate(ncidxz,varid, 0, 1, deflate_level = 2) )
      call check( nf90_put_att(ncidxz, varid, '_FillValue',nc_fillvalue) )
      call check( nf90_def_var(ncidxz, 'thlxz', NF90_REAL, (/dimid(3),dimid(2),dimid(1)/), varid ) )
      call check( nf90_put_att(ncidxz,varid,'longname','y averaged liquid water potential temperature') )
      call check( nf90_put_att(ncidxz,varid,'units','K') )
      call check( nf90_def_var_deflate(ncidxz,varid, 0, 1, deflate_level = 2) )
      call check( nf90_put_att(ncidxz, varid, '_FillValue',nc_fillvalue) )
      call check( nf90_def_var(ncidxz, 'qtxz', NF90_REAL, (/dimid(3),dimid(2),dimid(1)/), varid ) )
      call check( nf90_put_att(ncidxz,varid,'longname','y averaged total water specific humidity') )
      call check( nf90_put_att(ncidxz,varid,'units','kg/kg') )
      call check( nf90_def_var_deflate(ncidxz,varid, 0, 1, deflate_level = 2) )
      call check( nf90_put_att(ncidxz, varid, '_FillValue',nc_fillvalue) )
      call check( nf90_def_var(ncidxz, 'u2xz', NF90_REAL, (/dimid(5),dimid(2),dimid(1)/), varid ) )
      call check( nf90_put_att(ncidxz,varid,'longname','west-east velocity variance in y direction') )
      call check( nf90_put_att(ncidxz,varid,'units','m^2/s^2') )
      call check( nf90_def_var_deflate(ncidxz,varid, 0, 1, deflate_level = 2) )
      call check( nf90_put_att(ncidxz, varid, '_FillValue',nc_fillvalue) )
      call check( nf90_def_var(ncidxz, 'v2xz', NF90_REAL, (/dimid(3),dimid(2),dimid(1)/), varid ) )
      call check( nf90_put_att(ncidxz,varid,'longname','south-north velocity variance in y direction') )
      call check( nf90_put_att(ncidxz,varid,'units','m^2/s^2') )
      call check( nf90_def_var_deflate(ncidxz,varid, 0, 1, deflate_level = 2) )
      call check( nf90_put_att(ncidxz, varid, '_FillValue',nc_fillvalue) )
      call check( nf90_def_var(ncidxz, 'w2xz', NF90_REAL, (/dimid(3),dimid(4),dimid(1)/), varid ) )
      call check( nf90_put_att(ncidxz,varid,'longname','vertical velocity variance in y direction') )
      call check( nf90_put_att(ncidxz,varid,'units','m^2/s^2') )
      call check( nf90_def_var_deflate(ncidxz,varid, 0, 1, deflate_level = 2) )
      call check( nf90_put_att(ncidxz, varid, '_FillValue',nc_fillvalue) )
      call check( nf90_def_var(ncidxz, 'thl2xz', NF90_REAL, (/dimid(3),dimid(2),dimid(1)/), varid ) )
      call check( nf90_put_att(ncidxz,varid,'longname','liquid water potential temperature variance in y direction') )
      call check( nf90_put_att(ncidxz,varid,'units','K^2') )
      call check( nf90_def_var_deflate(ncidxz,varid, 0, 1, deflate_level = 2) )
      call check( nf90_put_att(ncidxz, varid, '_FillValue',nc_fillvalue) )
      call check( nf90_def_var(ncidxz, 'qt2xz', NF90_REAL, (/dimid(3),dimid(2),dimid(1)/), varid ) )
      call check( nf90_put_att(ncidxz,varid,'longname','total water specific humidity variance in y direction') )
      call check( nf90_put_att(ncidxz,varid,'units','kg^2/kg^2') )
      call check( nf90_def_var_deflate(ncidxz,varid, 0, 1, deflate_level = 2) )
      call check( nf90_put_att(ncidxz, varid, '_FillValue',nc_fillvalue) )
      call check( nf90_def_var(ncidxz, 'uwxz', NF90_REAL, (/dimid(5),dimid(4),dimid(1)/), varid ) )
      call check( nf90_put_att(ncidxz,varid,'longname','momentum flux (uw) from y direction') )
      call check( nf90_put_att(ncidxz,varid,'units','m^2/s^2') )
      call check( nf90_def_var_deflate(ncidxz,varid, 0, 1, deflate_level = 2) )
      call check( nf90_put_att(ncidxz, varid, '_FillValue',nc_fillvalue) )
      call check( nf90_def_var(ncidxz, 'vwxz', NF90_REAL, (/dimid(3),dimid(4),dimid(1)/), varid ) )
      call check( nf90_put_att(ncidxz,varid,'longname','momentum flux (vw) from y direction') )
      call check( nf90_put_att(ncidxz,varid,'units','m^2/s^2') )
      call check( nf90_def_var_deflate(ncidxz,varid, 0, 1, deflate_level = 2) )
      call check( nf90_put_att(ncidxz, varid, '_FillValue',nc_fillvalue) )
      call check( nf90_def_var(ncidxz, 'thlwxz', NF90_REAL, (/dimid(3),dimid(4),dimid(1)/), varid ) )
      call check( nf90_put_att(ncidxz,varid,'longname','liquid water potential temperature flux (thlw) from y direction') )
      call check( nf90_put_att(ncidxz,varid,'units','K m/s') )
      call check( nf90_def_var_deflate(ncidxz,varid, 0, 1, deflate_level = 2) )
      call check( nf90_put_att(ncidxz, varid, '_FillValue',nc_fillvalue) )
      call check( nf90_def_var(ncidxz, 'qtwxz', NF90_REAL, (/dimid(3),dimid(4),dimid(1)/), varid ) )
      call check( nf90_put_att(ncidxz,varid,'longname','total water specific humidity flux (qtw) from y direction') )
      call check( nf90_put_att(ncidxz,varid,'units','kg/kg m/s') )
      call check( nf90_def_var_deflate(ncidxz,varid, 0, 1, deflate_level = 2) )
      call check( nf90_put_att(ncidxz, varid, '_FillValue',nc_fillvalue) )
      call check( nf90_enddef(ncidxz) )
      call check( nf90_inq_varid(ncidxz, 'zt', varid) )
      call check( nf90_put_var(ncidxz, varid, zf(1:kmax)) )
      call check( nf90_inq_varid(ncidxz, 'xt', varid) )
      call check( nf90_put_var(ncidxz, varid, (/(dx*i+myidx*imax*dx,i=0,imax-1)/)) )
      call check( nf90_inq_varid(ncidxz, 'zm', varid) )
      call check( nf90_put_var(ncidxz, varid, zh(1:k1)) )
      call check( nf90_inq_varid(ncidxz, 'xm', varid) )
      call check( nf90_put_var(ncidxz, varid, (/(dx*i+myidx*imax*dx,i=0,imax-1+merge(1,0,lopenbc.and.lboundary(2)))/)) )
    endif
    if(myidx==0.and.lreducex) then
      fnameyz(13:15) = cmyidy
      fnameyz(17:19) = cexpnr
    endif
  end subroutine initmeancrosssection

  subroutine exitmeancrosssection
    use modmpi, only: myidy
    use modglobal, only : lboundary,lopenbc
    implicit none
    if(.not.lmeancross) return
    deallocate(usum, vsum, wsum, thlsum, qtsum, &
              u2sum,v2sum,w2sum,thl2sum,qt2sum, &
              uwsum,vwsum,      thlwsum,qtwsum)
    deallocate(uxz, vxz, wxz, thlxz, qtxz, &
              u2xz,v2xz,w2xz,thl2xz,qt2xz, &
              uwxz,vwxz,     thlwxz,qtwxz )
    deallocate(uxz_int, vxz_int, wxz_int, thlxz_int, qtxz_int,&
              u2xz_int,v2xz_int,w2xz_int,thl2xz_int,qt2xz_int,&
              uwxz_int,vwxz_int,         thlwxz_int,qtwxz_int )
    if(myidy==0) then
      deallocate(usum_glob, vsum_glob, wsum_glob, thlsum_glob, qtsum_glob, &
                u2sum_glob,v2sum_glob,w2sum_glob,thl2sum_glob,qt2sum_glob, &
                uwsum_glob,vwsum_glob,           thlwsum_glob,qtwsum_glob)
    endif
    if(myidy==0.and.lreducey) call check( nf90_close(ncidxz) )
  end subroutine exitmeancrosssection

  subroutine meancrosssection
    use modglobal, only : timee, dt_lim, rk3step
    use modmpi, only : myid
    implicit none
    if(.not.lmeancross) return
    if (rk3step/=3) return
    if(timee<tnext .and. timee<tnextwrite) then
      dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
      return
    end if
    if (timee>=tnext) then ! Sample statistics
      if(lreducey) call reducey
      if(lreducex) call reducex
      ! Set next sample time
      tnext = tnext+idtav
    end if
    if (timee>=tnextwrite) then ! Write statistics
      call writemeancrosssection
      ! Set next write time
      tnextwrite = tnextwrite+itimeav
      ! Reset integration sums
      uxz_int    = 0.; vxz_int = 0.; wxz_int = 0.; thlxz_int = 0.; qtxz_int = 0.
      u2xz_int   = 0.;v2xz_int = 0.;w2xz_int = 0.;thl2xz_int = 0.;qt2xz_int = 0.
      uwxz_int   = 0.;vwxz_int = 0.;thlwxz_int = 0.;qtwxz_int = 0.
    end if
    dt_lim = minval((/dt_lim,tnext-timee,tnextwrite-timee/))
  end subroutine meancrosssection

  subroutine reducey
    use mpi
    use modfields, only : um,vm,wm,thlm,qtm
    use modglobal, only : imax,jmax,kmax,i1,k1,i2,j2,lopenbc,lboundary,j1,dzh,dzf
    use modmpi, only : commcol,myidy,MY_REAL,mpierr,nprocy
    implicit none
    integer :: i,j,k
    ! Set initial sums to 0
    usum  = 0.; vsum  = 0.; wsum  = 0.; thlsum  = 0.; qtsum  = 0.
    u2sum = 0.; v2sum = 0.; w2sum = 0.; thl2sum = 0.; qt2sum = 0.
    uwsum = 0.; vwsum = 0.; thlwsum = 0.; qtwsum = 0.
    ! Start spatial reduction in y direction
    do j = 2,j1
      do i = 2,i1
        ! k=1 only for full levels with respect to k
        usum(i,1) = usum(i,1)+um(i,j,1)
        vsum(i,1) = vsum(i,1)+vm(i,j,1)
        thlsum(i,1) = thlsum(i,1)+thlm(i,j,1)
        qtsum(i,1) = qtsum(i,1)+qtm(i,j,1)
        u2sum(i,1) = u2sum(i,1)+um(i,j,1)**2
        v2sum(i,1) = v2sum(i,1)+vm(i,j,1)**2
        thl2sum(i,1) = thl2sum(i,1)+thlm(i,j,1)**2
        qt2sum(i,1) = qt2sum(i,1)+qtm(i,j,1)**2
        do k = 2,k1
          usum(i,k) = usum(i,k)+um(i,j,k)
          vsum(i,k) = vsum(i,k)+vm(i,j,k)
          wsum(i,k) = wsum(i,k)+wm(i,j,k)
          thlsum(i,k) = thlsum(i,k)+thlm(i,j,k)
          qtsum(i,k) = qtsum(i,k)+qtm(i,j,k)
          u2sum(i,k) = u2sum(i,k)+um(i,j,k)**2
          v2sum(i,k) = v2sum(i,k)+vm(i,j,k)**2
          w2sum(i,k) = w2sum(i,k)+wm(i,j,k)**2
          thl2sum(i,k) = thl2sum(i,k)+thlm(i,j,k)**2
          qt2sum(i,k) = qt2sum(i,k)+qtm(i,j,k)**2
          uwsum(i,k) = uwsum(i,k)+(wm(i,j,k)+wm(i-1,j,k))*(um(i,j,k)*dzf(k-1)+um(i,j,k-1)*dzf(k))/(4*dzh(k))
          vwsum(i,k) = vwsum(i,k)+(wm(i,j,k)+wm(i,j-1,k))*(vm(i,j,k)*dzf(k-1)+vm(i,j,k-1)*dzf(k))/(4*dzh(k))
          thlwsum(i,k) = thlwsum(i,k)+wm(i,j,k)*(thlm(i,j,k)*dzf(k-1)+thlm(i,j,k-1)*dzf(k))/(2*dzh(k))
          qtwsum(i,k) = qtwsum(i,k)+wm(i,j,k)*(qtm(i,j,k)*dzf(k-1)+qtm(i,j,k-1)*dzf(k))/(2*dzh(k))
        end do
      end do
      ! i=1 and i=i2 for w, and i=i2 for u and u2
      do k = 2,k1
        wsum(1,k)   = wsum(1,k)+wm(1,j,k)
        wsum(i2,k)  = wsum(i2,k)+wm(i2,j,k)
        usum(i2,k)  = usum(i2,k)+um(i2,j,k)
        u2sum(i2,k) = u2sum(i2,k)+um(i2,j,k)**2
      end do
    end do
    if(lopenbc .and. lboundary(4)) then ! Add j = j2 for v fields
      j = j2
      do i = 2,i1
        ! Treat first height level
        vsum(i,1) = vsum(i,1)+vm(i,j,1)
        v2sum(i,1) = v2sum(i,1)+vm(i,j,1)**2
        do k =2,k1 ! Other height levels
          vsum(i,k) = vsum(i,k)+vm(i,j,k)
          v2sum(i,k) = v2sum(i,k)+vm(i,j,k)**2
          vwsum(i,k) = vwsum(i,k)+(wm(i,j,k)+wm(i,j-1,k))*(vm(i,j,k)*dzf(k-1)+vm(i,j,k-1)*dzf(k))/(4*dzh(k))
        end do
      end do
    endif
    ! Sum over processies in y direction
    call mpi_reduce(usum,usum_glob,i1*k1,MY_REAL,MPI_SUM,0,commcol,mpierr)
    call mpi_reduce(vsum,vsum_glob,imax*k1,MY_REAL,MPI_SUM,0,commcol,mpierr)
    call mpi_reduce(wsum,wsum_glob,i2*k1,MY_REAL,MPI_SUM,0,commcol,mpierr)
    call mpi_reduce(thlsum,thlsum_glob,imax*k1,MY_REAL,MPI_SUM,0,commcol,mpierr)
    call mpi_reduce(qtsum,qtsum_glob,imax*k1,MY_REAL,MPI_SUM,0,commcol,mpierr)
    call mpi_reduce(u2sum,u2sum_glob,i1*k1,MY_REAL,MPI_SUM,0,commcol,mpierr)
    call mpi_reduce(v2sum,v2sum_glob,imax*k1,MY_REAL,MPI_SUM,0,commcol,mpierr)
    call mpi_reduce(w2sum,w2sum_glob,imax*k1,MY_REAL,MPI_SUM,0,commcol,mpierr)
    call mpi_reduce(thl2sum,thl2sum_glob,imax*k1,MY_REAL,MPI_SUM,0,commcol,mpierr)
    call mpi_reduce(qt2sum,qt2sum_glob,imax*k1,MY_REAL,MPI_SUM,0,commcol,mpierr)
    call mpi_reduce(uwsum,uwsum_glob,i1*k1,MY_REAL,MPI_SUM,0,commcol,mpierr)
    call mpi_reduce(vwsum,vwsum_glob,imax*k1,MY_REAL,MPI_SUM,0,commcol,mpierr)
    call mpi_reduce(thlwsum,thlwsum_glob,imax*k1,MY_REAL,MPI_SUM,0,commcol,mpierr)
    call mpi_reduce(qtwsum,qtwsum_glob,imax*k1,MY_REAL,MPI_SUM,0,commcol,mpierr)
    if(myidy==0) then
      ! Calculate variables
      uxz    = usum_glob(2:i1+merge(1,0,lopenbc.and.lboundary(2)),1:kmax)/(jmax*nprocy)
      vxz    = vsum_glob(2:i1,1:kmax)/(jmax*nprocy)
      wxz    = wsum_glob(2:i1,1:k1)/(jmax*nprocy)
      thlxz  = thlsum_glob(2:i1,1:kmax)/(jmax*nprocy)
      qtxz   = qtsum_glob(2:i1,1:kmax)/(jmax*nprocy)
      u2xz   = max(u2sum_glob(2:i1+merge(1,0,lopenbc.and.lboundary(2)),1:kmax)/(jmax*nprocy)-uxz**2,0.)
      v2xz   = max(v2sum_glob(2:i1,1:kmax)/(jmax*nprocy)-vxz**2,0.)
      w2xz   = max(w2sum_glob(2:i1,1:k1)/(jmax*nprocy)-wxz**2,0.)
      thl2xz = max(thl2sum_glob(2:i1,1:kmax)/(jmax*nprocy)-thlxz**2,0.)
      qt2xz  = max(qt2sum_glob(2:i1,1:kmax)/(jmax*nprocy)-qtxz**2,0.)
      do k = 2,k1
        do i = 2,i1
          uwxz(i-1,k) = uwsum_glob(i,k)/(jmax*nprocy)- &
            (usum_glob(i,k)*dzf(k-1)+usum_glob(i,k-1)*dzf(k))*(wsum_glob(i,k)+wsum_glob(i-1,k))/ &
            (4*dzh(k)*(jmax*nprocy)**2)
          vwxz(i-1,k) = vwsum_glob(i,k)/(jmax*nprocy)- &
            (vsum_glob(i,k)*dzf(k-1)+vsum_glob(i,k-1)*dzf(k))*wsum_glob(i,k)/ &
            (2*dzh(k)*(jmax*nprocy)**2)
          thlwxz(i-1,k) = thlwsum_glob(i,k)/(jmax*nprocy)- &
            (thlsum_glob(i,k)*dzf(k-1)+thlsum_glob(i,k-1)*dzf(k))*wsum_glob(i,k)/ &
            (2*dzh(k)*(jmax*nprocy)**2)
          qtwxz(i-1,k) = qtwsum_glob(i,k)/(jmax*nprocy)-&
            (qtsum_glob(i,k)*dzf(k-1)+qtsum_glob(i,k-1)*dzf(k))*wsum_glob(i,k)/ &
            (2*dzh(k)*(jmax*nprocy)**2)
        end do
      end do
      if(lopenbc.and.lboundary(2)) then
        do k = 2,k1
          uwxz(i2-1,k) = uwsum_glob(i2,k)/(jmax*nprocy)- &
            (usum_glob(i2,k)*dzf(k-1)+usum_glob(i2,k-1)*dzf(k))*(wsum_glob(i2,k)+wsum_glob(i2-1,k))/ &
            (4*dzh(k)*(jmax*nprocy)**2)
        end do
      endif
      ! Add to time integration
      uxz_int    = uxz_int+uxz
      vxz_int    = vxz_int+vxz
      wxz_int    = wxz_int+wxz
      thlxz_int  = thlxz_int+thlxz
      qtxz_int   = qtxz_int+qtxz
      u2xz_int   = u2xz_int+u2xz
      v2xz_int   = v2xz_int+v2xz
      w2xz_int   = w2xz_int+w2xz
      thl2xz_int = thl2xz_int+thl2xz
      qt2xz_int  = qt2xz_int+qt2xz
      uwxz_int   = uwxz_int+uwxz
      vwxz_int   = vwxz_int+vwxz
      thlwxz_int = thlwxz_int+thlwxz
      qtwxz_int  = qtwxz_int+qtwxz
    endif
  end subroutine reducey

  subroutine reducex
    implicit none

  end subroutine reducex

  subroutine writemeancrosssection
    use modglobal, only : imax,kmax,lopenbc,lboundary,k1,rtimee,lopenbc
    use modmpi, only : myidy
    implicit none
    if(lreducey.and.myidy==0) then
      call check( nf90_inq_varid(ncidxz, 'time', varid) )
      call check( nf90_put_var(ncidxz, varid, (/rtimee/), start = (/ntime/), count=(/1/)) )
      call check( nf90_inq_varid(ncidxz, 'uxz', varid) )
      call check( nf90_put_var(ncidxz, varid, uxz_int/nsamples, start    = (/1,1,ntime/), count=(/imax+merge(1,0,lopenbc.and.lboundary(2)),kmax,1/)) )
      call check( nf90_inq_varid(ncidxz, 'vxz', varid) )
      call check( nf90_put_var(ncidxz, varid, vxz_int/nsamples, start    = (/1,1,ntime/), count=(/imax,kmax,1/)) )
      call check( nf90_inq_varid(ncidxz, 'wxz', varid) )
      call check( nf90_put_var(ncidxz, varid, wxz_int/nsamples, start    = (/1,1,ntime/), count=(/imax,k1,1/)) )
      call check( nf90_inq_varid(ncidxz, 'thlxz', varid) )
      call check( nf90_put_var(ncidxz, varid, thlxz_int/nsamples, start  = (/1,1,ntime/), count=(/imax,kmax,1/)) )
      call check( nf90_inq_varid(ncidxz, 'qtxz', varid) )
      call check( nf90_put_var(ncidxz, varid, qtxz_int/nsamples, start   = (/1,1,ntime/), count=(/imax,kmax,1/)) )
      call check( nf90_inq_varid(ncidxz, 'u2xz', varid) )
      call check( nf90_put_var(ncidxz, varid, u2xz_int/nsamples, start   = (/1,1,ntime/), count=(/imax+merge(1,0,lopenbc.and.lboundary(2)),kmax,1/)) )
      call check( nf90_inq_varid(ncidxz, 'v2xz', varid) )
      call check( nf90_put_var(ncidxz, varid, v2xz_int/nsamples, start   = (/1,1,ntime/), count=(/imax,kmax,1/)) )
      call check( nf90_inq_varid(ncidxz, 'w2xz', varid) )
      call check( nf90_put_var(ncidxz, varid, w2xz_int/nsamples, start   = (/1,1,ntime/), count=(/imax,k1,1/)) )
      call check( nf90_inq_varid(ncidxz, 'thl2xz', varid) )
      call check( nf90_put_var(ncidxz, varid, thl2xz_int/nsamples, start = (/1,1,ntime/), count=(/imax,kmax,1/)) )
      call check( nf90_inq_varid(ncidxz, 'qt2xz', varid) )
      call check( nf90_put_var(ncidxz, varid, qt2xz_int/nsamples, start  = (/1,1,ntime/), count=(/imax,kmax,1/)) )
      call check( nf90_inq_varid(ncidxz, 'uwxz', varid) )
      call check( nf90_put_var(ncidxz, varid, uwxz_int/nsamples, start   = (/1,1,ntime/), count=(/imax+merge(1,0,lopenbc.and.lboundary(2)),k1,1/)) )
      call check( nf90_inq_varid(ncidxz, 'vwxz', varid) )
      call check( nf90_put_var(ncidxz, varid, vwxz_int/nsamples, start   = (/1,1,ntime/), count=(/imax,k1,1/)) )
      call check( nf90_inq_varid(ncidxz, 'thlwxz', varid) )
      call check( nf90_put_var(ncidxz, varid, thlwxz_int/nsamples, start = (/1,1,ntime/), count=(/imax,k1,1/)) )
      call check( nf90_inq_varid(ncidxz, 'qtwxz', varid) )
      call check( nf90_put_var(ncidxz, varid, qtwxz_int/nsamples, start  = (/1,1,ntime/), count=(/imax,k1,1/)) )
    endif
    ntime = ntime+1
  end subroutine writemeancrosssection

  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
  end subroutine check
end module modmeancrosssection
