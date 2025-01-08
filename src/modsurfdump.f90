!> \file modsurfdump.f90
!!  Dumps instantaneous surface fields and fluxes
!>
!!  xy-cross sections leads the dcape.myid.expnr.nc output
!>
!!  \author Frans Liqui Lung
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
module modsurfdump
use modsurfdata
use modglobal, only : longint,kmax
implicit none
private
PUBLIC :: initsurfdump,dosurfdump,exitsurfdump
save
  !NetCDF variables
  integer,parameter :: nvar = 11
  integer :: ncid4 = 0
  integer :: nrec = 0
  character(80) :: fname = 'surf.xxxxyxxx.xxx.nc'
  character(80),dimension(nvar,4) :: ncname
  character(80),dimension(1,4) :: tncname
  real    :: dtav
  integer(kind=longint) :: idtav,tnext
  logical :: lsurfdump = .false. !< switch for doing the surfdump (on/off)
contains

subroutine initsurfdump
    use modmpi,   only :myid,mpierr,comm3d,mpi_logical,cmyid,D_MPI_BCAST
    use modglobal,only :imax,jmax,ifnamopt,fname_options,dtmax,dtav_glob,ladaptive,dt_lim,cexpnr,tres,btime,checknamelisterror,&
                        output_prefix
    use modstat_nc,only : lnetcdf,open_nc, define_nc, redefine_nc,ncinfo,nctiminfo,writestat_dims_nc
    implicit none
    integer :: ierr

    namelist/NAMSURFDUMP/ &
    lsurfdump, dtav

    dtav = dtav_glob
    if(myid==0)then
      open(ifnamopt,file=fname_options,status='old',iostat=ierr)
      read (ifnamopt,NAMSURFDUMP,iostat=ierr)
      call checknamelisterror(ierr, ifnamopt, 'NAMSURFDUMP')
      write(6 ,NAMSURFDUMP)
      close(ifnamopt)
    end if

    call D_MPI_BCAST(dtav    ,1,0,comm3d,mpierr)
    call D_MPI_BCAST(lsurfdump   ,1,0,comm3d,mpierr)

    idtav = dtav/tres
    tnext   = idtav+btime
    if(.not.(lsurfdump)) return
    dt_lim = min(dt_lim,tnext)

    if (.not. ladaptive .and. abs(dtav/dtmax-nint(dtav/dtmax))>1e-4) then
      stop 'surfdump: dtav should be a integer multiple of dtmax'
    end if
    if (lnetcdf) then
    fname(6:13) = cmyid
    fname(15:17) = cexpnr
    call nctiminfo(tncname(1,:))
    call ncinfo(ncname( 1,:),'z0m','roughness length for momentum','m','tt0t')
    call ncinfo(ncname( 2,:),'z0h','roughness length for heat','m','tt0t')
    call ncinfo(ncname( 3,:),'rs','composite resistance','s/m','tt0t')
    call ncinfo(ncname( 4,:),'ra','aerodynamic resistance','s/m','tt0t')
    call ncinfo(ncname( 5,:),'thlflux','kinematic temperature flux','K m/s','tt0t')
    call ncinfo(ncname( 6,:),'qtflux','kinematic specific humidity flux','kg/kg m/s','tt0t')
    call ncinfo(ncname( 7,:),'ustar','friction velocity','m/s','tt0t')
    call ncinfo(ncname( 8,:),'obl','obukhov length','m','tt0t')
    call ncinfo(ncname( 9,:),'Cm','drag coefficient for momentum','-','tt0t')
    call ncinfo(ncname( 10,:),'Cs','Drag coefficient for scalars','-','tt0t')
    call ncinfo(ncname( 11,:),'tskin','skin temperature','K','tt0t')
    call open_nc(trim(output_prefix)//fname,  ncid4,nrec,n1=imax,n2=jmax)
    if (nrec==0) then
      call define_nc( ncid4, 1, tncname)
      call writestat_dims_nc(ncid4)
    end if
    call define_nc( ncid4, NVar, ncname)
    end if
end subroutine initsurfdump

subroutine dosurfdump
!>Run crosssection.
    use modglobal, only : imax,jmax,i1,j1,rk3step,timee,rtimee,dt_lim,linit_out
    use modstat_nc, only : lnetcdf, writestat_nc
    use modmpi
    implicit none
    real, allocatable :: vars(:,:,:)

    if (.not. lsurfdump) return
    if (rk3step/=3 .and. .not. linit_out) return
    if(timee<tnext .and. .not. linit_out) then
      dt_lim = min(dt_lim,tnext-timee)
      return
    end if
    if(.not. linit_out) then 
      tnext = tnext+idtav
      dt_lim = minval((/dt_lim,tnext-timee/))
    endif

    if (lnetcdf) then
      allocate(vars(1:imax,1:jmax,nvar))
      vars(:,:,1) = z0m(2:i1,2:j1)
      vars(:,:,2) = z0h(2:i1,2:j1)
      vars(:,:,3) = rs(2:i1,2:j1)
      vars(:,:,4) = ra(2:i1,2:j1)
      vars(:,:,5) = thlflux(2:i1,2:j1)
      vars(:,:,6) = qtflux(2:i1,2:j1)
      vars(:,:,7) = ustar(2:i1,2:j1)
      vars(:,:,8) = obl(2:i1,2:j1)
      vars(:,:,9) = Cm(2:i1,2:j1)
      vars(:,:,10) = Cs(2:i1,2:j1)
      vars(:,:,11)= tskin(2:i1,2:j1)
      call writestat_nc(ncid4,1,tncname,(/rtimee/),nrec,.true.)
      call writestat_nc(ncid4,nvar,ncname(1:nvar,:),vars,nrec,imax,jmax)
      deallocate(vars)
    end if
end subroutine dosurfdump

subroutine exitsurfdump
    use modstat_nc, only : exitstat_nc,lnetcdf
    implicit none
    if(lsurfdump .and. lnetcdf) then
    call exitstat_nc(ncid4)
    end if
end subroutine exitsurfdump
end module modsurfdump