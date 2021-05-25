!> \file modopenboundary.f90
!!  Sets ghost cells at domain boundary and implements radiation
!!  conditions for the boundary normal velocity components
!>
!!  Sets ghost cells at domain boundary and implements radiation
!!  conditions for the boundary normal velocity components
!>
!!  \author Frans Liqui Lung
!  This file is part of DALES.
!  To do:
!  - Allow for different zint and adjust rhointi accordingly
!  - Correct non divergence free input
!  - Allow for non-homogeneous starting conditions
!  - Change definition uphase for division by 0
!  - When to use nextval and currentval for nudging and check rtimee
!  - How to handle vertical derivative in top boundary condition (full levels)
!  - Clean loop for full level boundary conditions
!  - Use correct velocity level to determine in or outflow in full levels
!  - Check rtimee and half level nudgin and correction term
!  - Use um or u0 in half levels
!  - Check starting points in modforces
!  - Add possibility for higher order integration schemes
!  - Adjust turbulent pertubation generation
!  - Tau0 to input variables
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
module modopenboundary
use modglobal, only : boundary_type,boundary,lopenbc,lboundary,lperiodic,lsynturb,dzint,dxint,dyint,ntboundary,tboundary,tau0
use modsynturb, only : synturb,initsynturb,exitsynturb
use netcdf
implicit none
integer :: nxpatch, nypatch, nzpatch
real, dimension(:,:), allocatable :: uturbtemp,vturbtemp,wturbtemp
real, dimension(:), allocatable :: rhointi
integer :: nx1max,nx2max

contains
  subroutine initopenboundary
    ! Initialisation routine for openboundaries
    use modmpi, only : myidx, myidy, nprocx, nprocy
    use modglobal, only : imax,jmax,kmax,i1,j1,k1,dx,dy,dzf,itot,jtot,zf,zh,solver_id
    use modfields, only : rhobf
    implicit none
    integer :: i,j
    real :: yt(jmax),ym(j1)

    if(.not.lopenbc) return
    ! Check if hypre solver is selected
    if(solver_id /= 1) stop 'Openboundaries only possible with HYPRE pressure solver, change solver_id to 1'
    ! Check if boundary is present on process
    if(myidx==0)        lboundary(1) = .true.
    if(myidx==nprocx-1) lboundary(2) = .true.
    if(myidy==0)        lboundary(3) = .true.
    if(myidy==nprocy-1) lboundary(4) = .true.
    lboundary(5) = .true.
    ! Set boundary names
    boundary(1)%name = "west"
    boundary(2)%name = "east"
    boundary(3)%name = "south"
    boundary(4)%name = "north"
    boundary(5)%name = "top"
    ! Set dimension for each boundary
    boundary(1)%nx1  = jmax; boundary(1)%nx2  = kmax
    boundary(2)%nx1  = jmax; boundary(2)%nx2  = kmax
    boundary(3)%nx1  = imax; boundary(3)%nx2  = kmax
    boundary(4)%nx1  = imax; boundary(4)%nx2  = kmax
    boundary(5)%nx1  = imax; boundary(5)%nx2  = jmax
    boundary(1)%nx1u = jmax; boundary(1)%nx2u = kmax
    boundary(2)%nx1u = jmax; boundary(2)%nx2u = kmax
    boundary(3)%nx1u = i1;   boundary(3)%nx2u = kmax
    boundary(4)%nx1u = i1;   boundary(4)%nx2u = kmax
    boundary(5)%nx1u = i1;   boundary(5)%nx2u = jmax
    boundary(1)%nx1v = j1;   boundary(1)%nx2v = kmax
    boundary(2)%nx1v = j1;   boundary(2)%nx2v = kmax
    boundary(3)%nx1v = imax; boundary(3)%nx2v = kmax
    boundary(4)%nx1v = imax; boundary(4)%nx2v = kmax
    boundary(5)%nx1v = imax; boundary(5)%nx2v = j1
    boundary(1)%nx1w = jmax; boundary(1)%nx2w = k1
    boundary(2)%nx1w = jmax; boundary(2)%nx2w = k1
    boundary(3)%nx1w = imax; boundary(3)%nx2w = k1
    boundary(4)%nx1w = imax; boundary(4)%nx2w = k1
    boundary(5)%nx1w = imax; boundary(5)%nx2w = jmax
    ! Set number of patches for correction factor for radiation boundary conditions
    if(dxint == -1.) dxint = real(itot)*dx ! Set dxint to entire width as default
    if(dyint == -1.) dxint = real(jtot)*dy ! Set dyint to entire width as default
    nxpatch = int(dx/dxint*real(itot));
    nypatch = int(dy/dyint*real(jtot));
    nzpatch = kmax ! For now vertical integration scale is set equal to dz
    if(mod(dxint,dx)/=0 .or. mod(dyint,dy)/=0) stop 'dxint and dyint should be multiples of dx and dy respectively.'
    boundary(1)%nx1patch = nypatch; boundary(1)%nx2patch = nzpatch
    boundary(2)%nx1patch = nypatch; boundary(2)%nx2patch = nzpatch
    boundary(3)%nx1patch = nxpatch; boundary(3)%nx2patch = nzpatch
    boundary(4)%nx1patch = nxpatch; boundary(4)%nx2patch = nzpatch
    boundary(5)%nx1patch = nxpatch; boundary(5)%nx2patch = nypatch
    ! Allocate phase velocity and correction term radiation boundaries
    do i = 1,5
      if(.not.lboundary(i) .or. lperiodic(i)) cycle ! Open boundary not present
      allocate(boundary(i)%radcorr(boundary(i)%nx1patch,boundary(i)%nx2patch), &
        boundary(i)%radcorrsingle(boundary(i)%nx1patch,boundary(i)%nx2patch), &
        boundary(i)%uphase(boundary(i)%nx1,boundary(i)%nx2), &
        boundary(i)%uphasesingle(boundary(i)%nx1,boundary(i)%nx2), &
        boundary(i)%uturb(boundary(i)%nx1u,boundary(i)%nx2u), &
        boundary(i)%vturb(boundary(i)%nx1v,boundary(i)%nx2v), &
        boundary(i)%wturb(boundary(i)%nx1w,boundary(i)%nx2w), &
        boundary(i)%thlturb(boundary(i)%nx1,boundary(i)%nx2), &
        boundary(i)%qtturb(boundary(i)%nx1,boundary(i)%nx2), &
        boundary(i)%e12turb(boundary(i)%nx1,boundary(i)%nx2))
        boundary(i)%uturb = 0.; boundary(i)%vturb = 0.; boundary(i)%wturb = 0.
        boundary(i)%thlturb = 0.; boundary(i)%qtturb = 0.; boundary(i)%e12turb = 0.
    end do
    call initsynturb
  end subroutine initopenboundary

  subroutine exitopenboundary
    ! Exit routine for openboundaries
    implicit none
    integer :: i
    if(.not.lopenbc) return
    deallocate(tboundary)
    do i = 1,5
      if(.not.lboundary(i) .or. lperiodic(i)) cycle
      deallocate(boundary(i)%thl,boundary(i)%qt,boundary(i)%e12, &
        boundary(i)%u,boundary(i)%v,boundary(i)%w,boundary(i)%uphasesingle,boundary(i)%uphase, &
        boundary(i)%radcorr,boundary(i)%radcorrsingle,boundary(i)%uturb,boundary(i)%vturb, &
          boundary(i)%wturb,boundary(i)%thlturb,boundary(i)%qtturb,boundary(i)%name,boundary(i)%e12turb)
    end do
    deallocate(rhointi)
    call exitsynturb
  end subroutine exitopenboundary

  subroutine openboundary_readboundary
    use mpi
    use modglobal, only : dzf,kmax,cexpnr,imax,jmax,itot,jtot,k1,ntboundary, &
      tboundary,dzh,dx,dy,i1,j1,i2,j2,kmax,xsize,ysize,zh
    use modfields, only : rhobf,rhobh,uprof,vprof,thlprof,qtprof,e12prof,u0,um,v0,vm,w0,wm
    use modmpi, only : myid,comm3d,myidy,myidx,MY_REAL,nprocs
    implicit none
    integer :: it,i,j,k,ib,sy,sx,ey,ex,ip
    character(len = nf90_max_name) :: RecordDimName
    integer :: VARID,STATUS,NCID,mpierr,timeID
    integer, dimension(3) :: istart
    real,dimension(:),allocatable :: sumdiv,sumdivtot

    if(.not.lopenbc) return
    do ip = 0,nprocs-1
    call MPI_Barrier(comm3d,mpierr)
    if(myid==ip) then ! Read netcdf file one by one
    !--- open nc file ---
    STATUS = NF90_OPEN('openboundaries.inp.'//cexpnr//'.nc', nf90_nowrite, NCID)
    if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
    !--- get time dimensions
    status = nf90_inq_dimid(ncid, "time", timeID)
    if (status /= nf90_noerr) call handle_err(status)
    status = nf90_inquire_dimension(NCID, timeID, len=ntboundary, name=RecordDimName)
    if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
    !--- read time
    allocate(tboundary(ntboundary))
    STATUS = NF90_INQ_VARID(NCID, 'time', VARID)
    if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
    STATUS = NF90_GET_VAR (NCID, VARID, tboundary, start=(/1/), count=(/ntboundary/) )
    if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
    do ib = 1,5 ! loop over boundaries
      ! Allocate fields
      if(.not.lboundary(ib) .or. lperiodic(ib)) cycle ! Open boundary not present
      allocate(boundary(ib)%thl(boundary(ib)%nx1,boundary(ib)%nx2,ntboundary), &
        boundary(ib)%qt(boundary(ib)%nx1,boundary(ib)%nx2,ntboundary),  &
        boundary(ib)%e12(boundary(ib)%nx1,boundary(ib)%nx2,ntboundary), &
        boundary(ib)%u(boundary(ib)%nx1u,boundary(ib)%nx2u,ntboundary), &
        boundary(ib)%v(boundary(ib)%nx1v,boundary(ib)%nx2v,ntboundary), &
        boundary(ib)%w(boundary(ib)%nx1w,boundary(ib)%nx2w,ntboundary), &
        )
      if(lsynturb) then ! Allocate turbulent fields
        allocate(boundary(ib)%u2(boundary(ib)%nx1patch,boundary(ib)%nx2patch,ntboundary), &
        & boundary(ib)%v2(boundary(ib)%nx1patch,boundary(ib)%nx2patch,ntboundary), &
        & boundary(ib)%w2(boundary(ib)%nx1patch,boundary(ib)%nx2patch,ntboundary), &
        & boundary(ib)%uv(boundary(ib)%nx1patch,boundary(ib)%nx2patch,ntboundary), &
        & boundary(ib)%uw(boundary(ib)%nx1patch,boundary(ib)%nx2patch,ntboundary), &
        & boundary(ib)%vw(boundary(ib)%nx1patch,boundary(ib)%nx2patch,ntboundary), &
        & boundary(ib)%thl2(boundary(ib)%nx1patch,boundary(ib)%nx2patch,ntboundary),&
        & boundary(ib)%qt2(boundary(ib)%nx1patch,boundary(ib)%nx2patch,ntboundary), &
        & boundary(ib)%wthl(boundary(ib)%nx1patch,boundary(ib)%nx2patch,ntboundary), &
        & boundary(ib)%wqt(boundary(ib)%nx1patch,boundary(ib)%nx2patch,ntboundary))
      endif
      ! Read fields
      select case(ib)
      case(1,2)
        istart = (/myidy*jmax+1,1,1/)
      case(3,4)
        istart = (/myidx*imax+1,1,1/)
      case(5)
        istart = (/myidx*imax+1,myidy*jmax+1,1/)
      end select
      ! Read u
      STATUS = NF90_INQ_VARID(NCID,'u'//boundary(ib)%name, VARID)
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%u, start=istart, &
        & count=(/boundary(ib)%nx1u,boundary(ib)%nx2u,ntboundary/))
      ! Read v
      STATUS = NF90_INQ_VARID(NCID, 'v'//boundary(ib)%name, VARID)
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%v, start=istart, &
        & count=(/boundary(ib)%nx1v,boundary(ib)%nx2v,ntboundary/))
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      ! Read w
      STATUS = NF90_INQ_VARID(NCID, 'w'//boundary(ib)%name, VARID)
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%w, start=istart, &
        & count=(/boundary(ib)%nx1w,boundary(ib)%nx2w,ntboundary/))
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      ! Read thl
      STATUS = NF90_INQ_VARID(NCID, 'thl'//boundary(ib)%name, VARID)
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%thl, start=istart, &
        & count=(/boundary(ib)%nx1,boundary(ib)%nx2,ntboundary/))
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      ! Read qt
      STATUS = NF90_INQ_VARID(NCID, 'qt'//boundary(ib)%name, VARID)
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%qt, start=istart, &
        & count=(/boundary(ib)%nx1,boundary(ib)%nx2,ntboundary/))
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      ! Read e12
      STATUS = NF90_INQ_VARID(NCID, 'e12'//boundary(ib)%name, VARID)
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%e12, start=istart, &
        & count=(/boundary(ib)%nx1,boundary(ib)%nx2,ntboundary/))
      if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      ! Read input for turbulent pertubations
      if(lsynturb) then
        ! Read u2
        STATUS = NF90_INQ_VARID(NCID, 'u2'//boundary(ib)%name, VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%u2, start=(/1,1,1/), &
          & count=(/boundary(ib)%nx1patch,boundary(ib)%nx2patch,ntboundary/))
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        ! Read v2
        STATUS = NF90_INQ_VARID(NCID, 'v2'//boundary(ib)%name, VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%v2, start=(/1,1,1/), &
          & count=(/boundary(ib)%nx1patch,boundary(ib)%nx2patch,ntboundary/))
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        ! Read w2
        STATUS = NF90_INQ_VARID(NCID, 'w2'//boundary(ib)%name, VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%w2, start=(/1,1,1/), &
          & count=(/boundary(ib)%nx1patch,boundary(ib)%nx2patch,ntboundary/))
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        ! Read uv
        STATUS = NF90_INQ_VARID(NCID, 'uv'//boundary(ib)%name, VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%uv, start=(/1,1,1/), &
          & count=(/boundary(ib)%nx1patch,boundary(ib)%nx2patch,ntboundary/))
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        ! Read uw
        STATUS = NF90_INQ_VARID(NCID, 'uw'//boundary(ib)%name, VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%uw, start=(/1,1,1/), &
          & count=(/boundary(ib)%nx1patch,boundary(ib)%nx2patch,ntboundary/))
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        ! Read vw
        STATUS = NF90_INQ_VARID(NCID, 'vw'//boundary(ib)%name, VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%vw, start=(/1,1,1/), &
          & count=(/boundary(ib)%nx1patch,boundary(ib)%nx2patch,ntboundary/))
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        ! Read thl2
        STATUS = NF90_INQ_VARID(NCID, 'thl2'//boundary(ib)%name, VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%thl2, start=(/1,1,1/), &
          & count=(/boundary(ib)%nx1patch,boundary(ib)%nx2patch,ntboundary/))
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        ! Read qt2
        STATUS = NF90_INQ_VARID(NCID, 'qt2'//boundary(ib)%name, VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%qt2, start=(/1,1,1/), &
          & count=(/boundary(ib)%nx1patch,boundary(ib)%nx2patch,ntboundary/))
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        ! Read wthl
        STATUS = NF90_INQ_VARID(NCID, 'wthl'//boundary(ib)%name, VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%wthl, start=(/1,1,1/), &
          & count=(/boundary(ib)%nx1patch,boundary(ib)%nx2patch,ntboundary/))
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        ! Read wqt
        STATUS = NF90_INQ_VARID(NCID, 'wqt'//boundary(ib)%name, VARID)
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%wqt, start=(/1,1,1/), &
          & count=(/boundary(ib)%nx1patch,boundary(ib)%nx2patch,ntboundary/))
        if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        ! ! Read randthl
        ! STATUS = NF90_INQ_VARID(NCID, 'randthl'//boundary(ib)%name, VARID)
        ! if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        ! STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%randthl, start=(/1,1/), &
        !   & count=(/boundary(ib)%nx1,boundary(ib)%nx2/))
        ! if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        ! ! Read randqt
        ! STATUS = NF90_INQ_VARID(NCID, 'randqt'//boundary(ib)%name, VARID)
        ! if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
        ! STATUS = NF90_GET_VAR (NCID, VARID, boundary(ib)%randqt, start=(/1,1/), &
        !   & count=(/boundary(ib)%nx1,boundary(ib)%nx2/))
        ! if (STATUS .ne. nf90_noerr) call handle_err(STATUS)
      endif
    end do
    status = nf90_close(ncid)
    if (status /= nf90_noerr) call handle_err(status)
    endif
    end do

    ! Check divergence of data
    allocate(sumdiv(ntboundary)); sumdiv = 0.
    allocate(sumdivtot(ntboundary)); sumdivtot = 0.
    do it = 1,ntboundary
      if(lboundary(1)) then
        do j = 1,jmax
          do k = 1,kmax
            sumdiv(it) = sumdiv(it) - rhobf(k)*boundary(1)%u(j,k,it)*dzf(k)*dy
          end do
        end do
      endif
      if(lboundary(2)) then
        do j = 1,jmax
          do k = 1,kmax
            sumdiv(it) = sumdiv(it) + rhobf(k)*boundary(2)%u(j,k,it)*dzf(k)*dy
          end do
        end do
      endif
      if(lboundary(3)) then
        do i = 1,imax
          do k = 1,kmax
            sumdiv(it) = sumdiv(it) - rhobf(k)*boundary(3)%v(i,k,it)*dzf(k)*dx
          end do
        end do
      endif
      if(lboundary(4)) then
        do i = 1,imax
          do k = 1,kmax
            sumdiv(it) = sumdiv(it) + rhobf(k)*boundary(4)%v(i,k,it)*dzf(k)*dx
          end do
        end do
      endif
      if(lboundary(5)) then
        do i = 1,imax
          do j = 1,jmax
            sumdiv(it) = sumdiv(it) + rhobh(k1)*boundary(5)%w(i,j,it)*dx*dy
          end do
        end do
      endif
    end do
    call MPI_ALLREDUCE(sumdiv,sumdivtot,ntboundary,MY_REAL,MPI_SUM,comm3d,mpierr)
    print *, 'mean and max abs divergence',sum(abs(sumdivtot))/ntboundary,maxval(abs(sumdivtot))
    if(lboundary(5)) then ! Any divergence is compensated at the top boundary
      do it = 1,ntboundary
        boundary(5)%w(:,:,it)=boundary(5)%w(:,:,it)-sumdivtot(it)/(rhobh(k1)*xsize*ysize)
      end do
    endif
    ! Copy data to boundary information
    if(lboundary(1).and..not.lperiodic(1)) then
      sy = myidy*jmax+1; ey = sy+jmax-1
      do j = 2,j1
        do k = 1,kmax
          u0(2,j,k) = boundary(1)%u(j-1,k,1)
          um(2,j,k) = boundary(1)%u(j-1,k,1)
        end do
      end do
    endif
    if(lboundary(2).and..not.lperiodic(2)) then
      sy = myidy*jmax+1; ey = sy+jmax-1
      do j = 2,j1
        do k = 1,kmax
          u0(i2,j,k) = boundary(2)%u(j-1,k,1)
          um(i2,j,k) = boundary(2)%u(j-1,k,1)
        end do
      end do
    endif
    if(lboundary(3).and..not.lperiodic(3)) then
      sx = myidx*imax+1; ex = sx+imax-1
      do i = 2,i1
        do k = 1,kmax
          v0(i,2,k) = boundary(3)%v(i-1,k,1)
          vm(i,2,k) = boundary(3)%v(i-1,k,1)
        end do
      end do
    endif
    if(lboundary(4).and..not.lperiodic(4)) then
      sx = myidx*imax+1; ex = sx+imax-1
      do i = 2,i1
        do k = 1,kmax
          v0(i,j2,1:k) = boundary(4)%v(i-1,k,1)
          vm(i,j2,1:k) = boundary(4)%v(i-1,k,1)
        end do
      end do
    endif
    if(lboundary(5).and..not.lperiodic(5)) then
      sx = myidx*imax+1; ex = sx+imax-1
      sy = myidy*jmax+1; ey = sy+jmax-1
      do i = 2,i1
        do j = 2,j1
          w0(i,j,k1) = boundary(5)%w(i-1,j-1,1)
          wm(i,j,k1) = boundary(5)%w(i-1,j-1,1)
        end do
      end do
    endif
    allocate(rhointi(k1))
    rhointi = 1./(rhobf*dzf)
    deallocate(sumdiv,sumdivtot)
  end subroutine openboundary_readboundary

  subroutine openboundary_ghost
    ! Subroutine that fills the ghost cells for the cell centred variables at the boundary
    use modglobal, only : i1,j1,k1,ih,jh,nsv
    use modfields, only : um,u0,vm,v0,wm,w0,e12m,e120,thlm,thl0,qtm,qt0,svm,sv0,thl0av,qt0av,u0av,v0av
    implicit none
    integer :: i,n
    if(.not.lopenbc) return
    ! Apply non domain boundaries and forced periodic boundaries
    call openboundary_excjs(um   , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    call openboundary_excjs(u0   , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    call openboundary_excjs(vm   , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    call openboundary_excjs(v0   , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    call openboundary_excjs(wm   , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    call openboundary_excjs(w0   , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    call openboundary_excjs(e12m , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    call openboundary_excjs(e120 , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    call openboundary_excjs(thlm , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    call openboundary_excjs(thl0 , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    call openboundary_excjs(qtm  , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    call openboundary_excjs(qt0  , 2,i1,2,j1,1,k1,ih,jh,.not.lboundary(1:4).or.lperiodic(1:4))
    do n = 1,nsv
      call openboundary_excjs(svm(:,:,:,n), 2,i1,2,j1,1,k1,ih,jh,(/.true.,.true.,.true.,.true./))
      call openboundary_excjs(sv0(:,:,:,n), 2,i1,2,j1,1,k1,ih,jh,(/.true.,.true.,.true.,.true./))
    end do
    ! Apply open boundaries for domain boundaries for full levels (ghost cells)
    do i = 1,5 ! Loop over boundaries
      if(.not.lboundary(i).or.lperiodic(i)) cycle
      call applyboundaryf(thlm ,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%thl,boundary(i)%nx1,boundary(i)%nx2,0,boundary(i)%thlturb,profile=thl0av)
      call applyboundaryf(thl0 ,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%thl,boundary(i)%nx1,boundary(i)%nx2,0,boundary(i)%thlturb,profile=thl0av)
      call applyboundaryf(qtm  ,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%qt,boundary(i)%nx1,boundary(i)%nx2,1,boundary(i)%qtturb,profile=qt0av)
      call applyboundaryf(qt0  ,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%qt,boundary(i)%nx1,boundary(i)%nx2,1,boundary(i)%qtturb,profile=qt0av)
      call applyboundaryf(e12m ,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%e12,boundary(i)%nx1,boundary(i)%nx2,1,boundary(i)%e12turb)
      call applyboundaryf(e120 ,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%e12,boundary(i)%nx1,boundary(i)%nx2,1,boundary(i)%e12turb)
      if(i/=1.and.i/=2) call applyboundaryf(um,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%u,boundary(i)%nx1u,boundary(i)%nx2u,0,boundary(i)%uturb,profile=u0av)
      if(i/=1.and.i/=2) call applyboundaryf(u0,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%u,boundary(i)%nx1u,boundary(i)%nx2u,0,boundary(i)%uturb,profile=u0av)
      if(i/=3.and.i/=4) call applyboundaryf(vm,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%v,boundary(i)%nx1v,boundary(i)%nx2v,0,boundary(i)%vturb,profile=v0av)
      if(i/=3.and.i/=4) call applyboundaryf(v0,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%v,boundary(i)%nx1v,boundary(i)%nx2v,0,boundary(i)%vturb,profile=v0av)
      if(i/=5) call applyboundaryf(wm,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%w,boundary(i)%nx1w,boundary(i)%nx2w,0,boundary(i)%wturb)
      if(i/=5) call applyboundaryf(w0,2,i1,2,j1,1,k1,ih,jh,i,boundary(i)%w,boundary(i)%nx1w,boundary(i)%nx2w,0,boundary(i)%wturb)
    end do
  end subroutine openboundary_ghost

  subroutine openboundary_tend
    ! Subroutine that handles the tendencies of the boundary normal velocity components
    ! Radiation boundary conditions are used for outflow boundaries and nudging
    ! boundary conditions for the inflow boundaries.
    ! Outflow:
    ! du/dt = -uphase*du/dx
    ! Inflow:
    ! du/dt = (u-uboundary)/tau
    ! There is a possibility to add synthetic turbulence at the inflow boundary.
    use modfields, only : up,vp,wp
    use modglobal, only : i1,j1,k1,kmax,i2,j2,jmax,imax
    implicit none
    integer :: i

    if(.not.lopenbc) return
    if(lboundary(1).and..not.lperiodic(1)) call applyboundaryh(1,boundary(1)%nx1,boundary(1)%nx2,boundary(1)%uturb)
    if(lboundary(2).and..not.lperiodic(2)) call applyboundaryh(2,boundary(2)%nx1,boundary(2)%nx2,boundary(2)%uturb)
    if(lboundary(3).and..not.lperiodic(3)) call applyboundaryh(3,boundary(3)%nx1,boundary(3)%nx2,boundary(3)%vturb)
    if(lboundary(4).and..not.lperiodic(4)) call applyboundaryh(4,boundary(4)%nx1,boundary(4)%nx2,boundary(4)%vturb)
    if(lboundary(5).and..not.lperiodic(5)) call applyboundaryh(5,boundary(5)%nx1,boundary(5)%nx2,boundary(5)%wturb)
    ! Calculate and add correction term to guarantee conservation of mass
    do i = 1,5
      if(.not. lboundary(i).or.lperiodic(i)) cycle
      call radcorrection(i)
    end do
  end subroutine openboundary_tend

  subroutine openboundary_phasevelocity
    ! Subroutine that calculates the phase velocity that is required for the
    ! radiation outflow boundary. The phase velocity is calculated from the phase
    ! velocity of one gridcell to the interior at the prior time and is averaged
    ! over the integration length scales.
    use mpi
    use modmpi, only : comm3d,commrow,commcol,myidx,myidy,mpierr,MY_REAL
    use modglobal, only : imax,jmax,kmax,i1,j1,dx,dy,dzf
    use modfields, only : u0,up,v0,vp,w0,wp,rhobh,rhobf
    implicit none
    integer :: ib,i,j,k,ipatch,jpatch,kpatch
    real :: ipos,jpos

    if(.not.lopenbc) return
    do ib = 1,5
      if(.not.lboundary(ib).or.lperiodic(ib)) cycle
      select case(ib) ! Select boundary
      case(1) ! West
        boundary(1)%uphasesingle=0.
        do j = 1,jmax
          jpos = j + (myidy * jmax) - 1
          jpatch = int((jpos-0.5)*dy/dyint)+1
          do k = 1,kmax
            kpatch = k
            boundary(1)%uphasesingle(jpatch,kpatch) = boundary(1)%uphasesingle(jpatch,kpatch) + &
              (-up(3,j+1,k)*dx/sign(max(abs(u0(4,j+1,k)-u0(3,j+1,k)),1e-10),u0(4,j+1,k)-u0(3,j+1,k))) &
              *dy*rhobf(k)*dzf(k)/dyint*rhointi(kpatch)
          end do
        end do
        ! Integrate over processes
        call MPI_ALLREDUCE(boundary(1)%uphasesingle,boundary(1)%uphase,nypatch*nzpatch,MY_REAL, &
                           MPI_SUM, commcol,mpierr)
      case(2) ! East
        boundary(2)%uphasesingle=0.
        do j = 1,jmax
          jpos = j + (myidy * jmax) - 1
          jpatch = int((jpos-0.5)*dy/dyint)+1
          do k = 1,kmax
            kpatch = k
            boundary(2)%uphasesingle(jpatch,kpatch) = boundary(2)%uphasesingle(jpatch,kpatch) + &
              (-up(i1,j+1,k)*dx/sign(max(abs(u0(i1,j+1,k)-u0(i1-1,j+1,k)),1e-10),u0(i1,j+1,k)-u0(i1-1,j+1,k))) &
              *dy*rhobf(k)*dzf(k)/dyint*rhointi(kpatch)
          end do
        end do
        ! Integrate over processes
        call MPI_ALLREDUCE(boundary(2)%uphasesingle,boundary(2)%uphase,nypatch*nzpatch,MY_REAL, &
                           MPI_SUM, commcol,mpierr)
      case(3) ! South
        boundary(3)%uphasesingle=0.
        do i = 1,imax
          ipos = i + (myidx * imax) - 1
          ipatch = int((ipos-0.5)*dx/dxint)+1
          do k = 1,kmax
            kpatch = k
            boundary(3)%uphasesingle(ipatch,kpatch) = boundary(3)%uphasesingle(ipatch,kpatch) + &
              (-vp(i+1,3,k)*dy/sign(max(abs(v0(i+1,4,k)-v0(i+1,3,k)),1e-10),v0(i+1,4,k)-v0(i+1,3,k))) &
              *dx*rhobf(k)*dzf(k)/dxint*rhointi(kpatch)
          end do
        end do
        ! Integrate over processes
        call MPI_ALLREDUCE(boundary(3)%uphasesingle,boundary(3)%uphase,nxpatch*nzpatch,MY_REAL, &
                           MPI_SUM, commrow,mpierr)
      case(4) ! North
        boundary(4)%uphasesingle=0.
        do i = 1,imax
          ipos = i + (myidx * imax) - 1
          ipatch = int((ipos-0.5)*dx/dxint)+1
          do k = 1,kmax
            kpatch = k
            boundary(4)%uphasesingle(ipatch,kpatch) = boundary(4)%uphasesingle(ipatch,kpatch) + &
              (-vp(i+1,j1,k)*dy/sign(max(abs(v0(i+1,j1,k)-v0(i+1,j1-1,k)),1e-10),v0(i+1,j1,k)-v0(i+1,j1-1,k))) &
              *dx*rhobf(k)/dxint*rhointi(kpatch)
          end do
        end do
        ! Integrate over processes
        call MPI_ALLREDUCE(boundary(4)%uphasesingle,boundary(4)%uphase,nxpatch*nzpatch,MY_REAL, &
                           MPI_SUM, commrow,mpierr)
      case(5) ! Top
        boundary(5)%uphasesingle=0.
        do i = 1,imax
          ipos = i + (myidx * imax) - 1
          ipatch = int((ipos-0.5)*dx/dxint)+1
          do j = 1,jmax
            jpos = j + (myidy * jmax) - 1
            jpatch = int((jpos-0.5)*dy/dyint)+1
            boundary(5)%uphasesingle(ipatch,jpatch) = boundary(5)%uphasesingle(ipatch,jpatch) + &
              (-wp(i+1,j+1,kmax)*dzf(kmax-1)/sign(max(abs(rhobh(kmax)*w0(i+1,j+1,kmax)-rhobh(kmax-1)*w0(i+1,j+1,kmax-1))/rhobh(kmax),1e-10), &
              & rhobh(kmax)*w0(i+1,j+1,kmax)-rhobh(kmax-1)*w0(i+1,j+1,kmax-1)))*dx*dy/(dxint*dyint)
          end do
        end do
        ! Integrate over processes
        call MPI_ALLREDUCE(boundary(5)%uphasesingle,boundary(5)%uphase,nxpatch*nypatch,MY_REAL, &
                           MPI_SUM, comm3d,mpierr)
      end select
    end do
  end subroutine openboundary_phasevelocity

  subroutine openboundary_turb
    ! Subroutine that calls the synthetic turbulence routine for the generation
    ! of synthetic turbulence at the dirichlet inflow boundaries.
    use modglobal, only : rk3step
    implicit none
    if(rk3step == 1) call synturb()
  end subroutine openboundary_turb

  subroutine openboundary_excjs(a,sx,ex,sy,ey,sz,ez,ih,jh,switch)
    ! Subroutine that handles periodic boundaries. Based on the excjs function
    ! present in modmpi
    use mpi
    use modmpi, only : comm3d, mpierr, my_real, nbreast, nbrnorth, nbrsouth, nbrwest, nprocx, nprocy, myidx, myidy,myid
    implicit none
    integer sx, ex, sy, ey, sz, ez, ih, jh
    real a(sx-ih:ex+ih, sy-jh:ey+jh, sz:ez)
    integer status(MPI_STATUS_SIZE)
    integer ii, i, j, k
    integer reqn, reqs, reqe, reqw
    integer nssize, ewsize
    logical, dimension(4), intent(in) :: switch
    real,allocatable, dimension(:) :: sendn,recvn
    real,allocatable, dimension(:) :: sends,recvs
    real,allocatable, dimension(:) :: sende,recve
    real,allocatable, dimension(:) :: sendw,recvw

    !   Calculate buffer size
    nssize = jh*(ex - sx + 1 + 2*ih)*(ez - sz + 1)
    ewsize = ih*(ey - sy + 1 + 2*jh)*(ez - sz + 1)
    !   Allocate send / receive buffers
    if(nprocx .gt. 1) then ! Processor communication required
      if(switch(1)) then ! Send west boundary
        allocate(sendw(ewsize),recvw(ewsize))
        ii = 0
        do i=1,ih
        do k=sz,ez
        do j=sy-jh,ey+jh
          ii = ii + 1
          sendw(ii) = a(sx+i-1,j,k)
        enddo
        enddo
        enddo
        call MPI_ISEND(sendw, ewsize, MY_REAL, nbrwest, 7, comm3d, reqw, mpierr)
      endif
      if(switch(2)) then ! Send east boundary
        allocate(sende(ewsize),recve(ewsize))
        ii = 0
        do i=1,ih
        do k=sz,ez
        do j=sy-jh,ey+jh
          ii = ii + 1
          sende(ii) = a(ex-i+1,j,k)
        enddo
        enddo
        enddo
        call MPI_ISEND(sende, ewsize, MY_REAL, nbreast, 6, comm3d, reqe, mpierr)
      endif
      if(switch(1)) then ! Receive west boundary
        call MPI_RECV(recvw, ewsize, MY_REAL, nbrwest, 6, comm3d, status, mpierr)
        ii = 0
        do i=1,ih
        do k=sz,ez
        do j=sy-jh,ey+jh
          ii = ii + 1
          a(sx-i,j,k) = recvw(ii)
        enddo
        enddo
        enddo
      endif
      if(switch(2)) then ! Receive east boundary
        call MPI_RECV(recve, ewsize, MY_REAL, nbreast, 7, comm3d, status, mpierr)
        ii = 0
        do i=1,ih
        do k=sz,ez
        do j=sy-jh,ey+jh
          ii = ii + 1
          a(ex+i,j,k) = recve(ii)
        enddo
        enddo
        enddo
      endif
      if(switch(1)) then
        call MPI_WAIT(reqw, status, mpierr)
        deallocate(sendw,recvw)
      endif
      if(switch(2)) then
        call MPI_WAIT(reqe, status, mpierr)
        deallocate(sende,recve)
      endif
    else ! Single processor
      if(switch(1)) then
        do i=1,ih
        do k=sz,ez
        do j=sy-jh,ey+jh
          a(sx-i,j,k) = a(ex-i+1,j,k)
          a(ex+i,j,k) = a(sx+i-1,j,k)
        enddo
        enddo
        enddo
      endif
    endif
    if(nprocy.gt.1) then ! Processor communication required
      if(switch(3)) then ! Send south boundary
        allocate(sends(nssize),recvs(nssize))
        ii = 0
        do j=1,jh
        do k=sz,ez
        do i=sx-ih,ex+ih
          ii = ii + 1
          sends(ii) = a(i,sy+j-1,k)
        enddo
        enddo
        enddo
        call MPI_ISEND(sends, nssize, MY_REAL, nbrsouth, 5, comm3d, reqs, mpierr)
      endif
      if(switch(4)) then ! Send north boundary
        allocate(sendn(nssize),recvn(nssize))
        ii = 0
        do j=1,jh
        do k=sz,ez
        do i=sx-ih,ex+ih
          ii = ii + 1
          sendn(ii) = a(i,ey-j+1,k)
        enddo
        enddo
        enddo
        call MPI_ISEND(sendn, nssize, MY_REAL, nbrnorth, 4, comm3d, reqn, mpierr)
      endif
      if(switch(3)) then ! Receive south boundary
        call MPI_RECV(recvs, nssize, MY_REAL, nbrsouth, 4, comm3d, status, mpierr)
        ii = 0
        do j=1,jh
        do k=sz,ez
        do i=sx-ih,ex+ih
          ii = ii + 1
          a(i,sy-j,k) = recvs(ii)
        enddo
        enddo
        enddo
      endif
      if(switch(4)) then ! Receive south boundary
        call MPI_RECV(recvn, nssize, MY_REAL, nbrnorth, 5, comm3d, status, mpierr)
        ii = 0
        do j=1,jh
        do k=sz,ez
        do i=sx-ih,ex+ih
          ii = ii + 1
          a(i,ey+j,k) = recvn(ii)
        enddo
        enddo
        enddo
      endif
      if(switch(3)) then
        call MPI_WAIT(reqs, status, mpierr)
        deallocate(sends,recvs)
      endif
      if(switch(4)) then
        call MPI_WAIT(reqn, status, mpierr)
        deallocate(sendn,recvn)
      endif
    else ! Single processor
      if(switch(3)) then
        do j=1,jh
        do k=sz,ez
        do i=sx-ih,ex+ih
          a(i,sy-j,k) = a(i,ey-j+1,k)
          a(i,ey+j,k) = a(i,sy+j-1,k)
        enddo
        enddo
        enddo
      endif
    endif
  end subroutine openboundary_excjs

  subroutine applyboundaryf(a,sx,ex,sy,ey,sz,ez,ih,jh,ib,val,nx1,nx2,lmax0,turb,profile)
    ! Routine fills ghost cells based on dirichlet (inflow) or
    ! homogeneous neumann (outflow) boundary conditions. Adds synthetic
    ! turbulent pertubations to dirichlet condition if lsynturb=true.
    use modglobal, only : dzh,dx,dy,imax,jmax,kmax,rtimee,rdt,i2,j2,k1,i1,j1
    use modfields, only : u0,v0,w0,e120
    use modmpi, only : myid
    implicit none
    integer, intent(in) :: sx,ex,sy,ey,sz,ez,ih,jh,ib,nx1,nx2,lmax0
    real, intent(in), dimension(nx1,nx2,ntboundary) :: val
    real, intent(in), dimension(nx1,nx2) :: turb
    real, intent(in), dimension(k1), optional :: profile ! optional for top boundary to take gradient into account
    real, intent(inout), dimension(sx-ih:ex+ih,sy-jh:ey+jh,sz:ez) :: a
    integer :: i,j,k,itp,itm,kav=5,itpn,itmn
    real :: coefdir,coefneu,tp,tm,fp,fm,fpn,fmn,ddz,valtarget,un,e

    ! Get interpolation coefficients for boundary input
    itm=1
    if(ntboundary>1) then
      do while(rtimee-rdt>tboundary(itm))
        itm=itm+1
      end do
      if (rtimee-rdt>tboundary(1)) then
        itm=itm-1
      end if
      itp = itm+1
      tm = tboundary(itm)
      tp = tboundary(itp)
      fm = (tp-rtimee+rdt)/(tp-tm)
      fp = (rtimee-rdt-tm)/(tp-tm)
    else
      itp = 1
      fp  = 0.
      fm  = 1.
    endif
    select case(ib) ! Select domain boundary
    case(1) ! West
      do k = 1,nx2
        do j = 1,nx1
          un = u0(sx,min(j+1,j1),min(k,kmax))
          if(un<=0) then ! Homogeneous Neumann outflow
            a(sx-1,j+1,k)=a(sx,j+1,k)
          else ! Robin inflow conditions
            e = e120(sx,min(j+1,j1),min(k,kmax))
            coefdir = 1.
            coefneu = -un*tau0-e/un*dx
            valtarget = fp*val(j,k,itp)+fm*val(j,k,itm)+turb(j,k)
            if(lmax0==1) valtarget = max(valtarget,0.)
            a(sx-1,j+1,k) = ( 2.*dx*valtarget - &
              a(sx,j+1,k)*(coefdir*dx+2.*coefneu) ) / (coefdir*dx-2.*coefneu)
          endif
        end do
      end do
    case(2) ! East
      do k = 1,nx2
        do j = 1,nx1
          un = u0(ex+1,min(j+1,j1),min(k,kmax))
          if(un>=0) then ! Homogeneous Neumann outflow
            a(ex+1,j+1,k)=a(ex,j+1,k)
          else ! Robin or Dirichlet inflow conditions
            e = e120(ex,min(j+1,j1),min(k,kmax))
            coefdir = 1.
            coefneu = -un*tau0-e/un*dx
            valtarget = fp*val(j,k,itp)+fm*val(j,k,itm)+turb(j,k)
            if(lmax0==1) valtarget = max(valtarget,0.)
            a(ex+1,j+1,k) = ( 2.*dx*valtarget - &
              a(ex,j+1,k)*(coefdir*dx-2.*coefneu) ) / (coefdir*dx+2.*coefneu)
          endif
        end do
      end do
    case(3) ! South
      do k = 1,nx2
        do i = 1,nx1
          un = v0(min(i+1,i1),sy,min(k,kmax))
          if(un<=0) then ! Homogeneous Neumann outflow
            a(i+1,sy-1,k)=a(i+1,sy,k)
          else ! Robin or Dirichlet inflow conditions
            e = e120(min(i+1,i1),sy,min(k,kmax))
            coefdir = 1.
            coefneu = -un*tau0-e/un*dy
            valtarget = fp*val(i,k,itp)+fm*val(i,k,itm)+turb(i,k)
            if(lmax0==1) valtarget = max(valtarget,0.)
            a(i+1,sy-1,k) = ( 2.*dy*valtarget - &
              a(i+1,sy,k)*(coefdir*dy+2.*coefneu) ) / (coefdir*dy-2.*coefneu)
          endif
        end do
      end do
    case(4) ! North
      do k = 1,nx2
        do i = 1,nx1
          un = v0(min(i+1,i1),ey+1,min(k,kmax))
          if(un>=0) then ! Homogeneous Neumann outflow
            a(i+1,ey+1,k)=a(i+1,ey,k)
          else ! Robin or Dirichlet inflow conditions
            e = e120(min(i+1,i1),ey,min(k,kmax))
            coefdir = 1.
            coefneu = -un*tau0-e/un*dy
            valtarget = fp*val(i,k,itp)+fm*val(i,k,itm)+turb(i,k)
            if(lmax0==1) valtarget = max(valtarget,0.)
            a(i+1,ey+1,k) = ( 2.*dy*valtarget - &
              a(i+1,ey,k)*(coefdir*dy-2.*coefneu) ) / (coefdir*dy+2.*coefneu)
          endif
        end do
      end do
    case(5) ! Top
      ! Obtain verticle gradient if slab averaged profile is given
      if(present(profile)) then
        ddz = sum((profile(kmax-kav+1:kmax)-profile(kmax-kav:kmax-1))/ &
                   dzh(kmax-kav+1:kmax))/kav
      else
        ddz = 0.
      endif
      do i = 1,nx1
        do j = 1,nx2
          un = w0(min(i+1,i1),min(j+1,j1),ez)
          if(un>=0) then ! Neumann outflow
            a(i+1,j+1,ez)=ddz*dzh(ez)+a(i+1,j+1,ez-1)
          else ! Robin inflow conditions
            e = e120(min(i+1,i1),min(j+1,j1),ez-1)
            coefdir = 1.
            coefneu = -un*tau0-e/un*dzh(ez)
            valtarget = fp*val(i,j,itp)+fm*val(i,j,itm)+turb(i,j)-(un*tau0+e/un*dzh(ez))*ddz
            if(lmax0==1) valtarget = max(valtarget,0.)
            a(i+1,j+1,ez) = ( 2.*dzh(ez)*valtarget - &
              a(i+1,j+1,ez-1)*(coefdir*dzh(ez)-2.*coefneu) ) / (coefdir*dzh(ez)+2.*coefneu)
          endif
        end do
      end do
    end select
  end subroutine applyboundaryf

  subroutine applyboundaryh(ib,nx1,nx2,turb)
    ! Subroutine that applies the radiation and dirichlet boundary conditions
    ! for the boundary normal velocity components. Adds synthetic turbulence to
    ! the inflow dirichlet boundaries if lsyntrub=.true.
    use mpi
    use modmpi, only : MY_REAL,myidx,myidy
    use modglobal, only : dx,dy,dzf,dxi,dyi,rdt,i2,j2,k1,i1,j1,kmax,rtimee,rdt,itot,jtot,imax,jmax,grav
    use modfields, only : um,u0,up,vm,v0,vp,wm,w0,wp,rhobf,rhobh,thvh,thv0h
    implicit none
    integer, intent(in) :: nx1,nx2,ib
    real, intent(in), dimension(nx1,nx2) :: turb
    integer :: i,j,k,itmc,itmn,itpc,itpn,ipatch,jpatch,kpatch
    real :: tm,tp,fpc,fmc,fpn,fmn,unext,uwallcurrent,ipos,jpos
    itmc=1
    itmn=1
    if(ntboundary>1) then
      do while(rtimee-rdt>tboundary(itmc))
        itmc=itmc+1
      end do
      if (rtimee-rdt>tboundary(1)) then
        itmc=itmc-1
      end if
      do while(tboundary(itmn)<rtimee)
        itmn=itmn+1
      end do
      if (rtimee>tboundary(1)) then
        itmn=itmn-1
      end if
      itpc = itmc+1
      itpn = itmn+1
      tm = tboundary(itmc)
      tp = tboundary(itpc)
      fmc = (tp-rtimee+rdt)/(tp-tm)
      fpc = (rtimee-rdt-tm)/(tp-tm)
      tm = tboundary(itmn)
      tp = tboundary(itpn)
      fmn = (tp-rtimee)/(tp-tm)
      fpn = (rtimee-tm)/(tp-tm)
    else
      itpc = 1
      itpn = 1
      fpc  = 0.
      fmc  = 1.
      fpn  = 0.
      fmn  = 1.
    endif
    ! Apply domain boundaries
    select case(ib) ! Select boundary
    case(1) ! West
      do j = 1,nx1
        jpos = j + (myidy * jmax) - 1
        jpatch = int((jpos-0.5)*dy/dyint)+1
        do k = 1,nx2
          kpatch = k
          uwallcurrent = fpc*boundary(1)%u(j,k,itpc)+fmc*boundary(1)%u(j,k,itmc)
          if(uwallcurrent<=0.) then ! Outflow (Radiation)
            up(2,j+1,k) = -max(min(boundary(1)%uphase(jpatch,kpatch),uwallcurrent),-dx/rdt) * &
              (u0(3,j+1,k)-u0(2,j+1,k))*dxi
          else ! Inflow (Dirichlet)
            unext = fpn*boundary(1)%u(j,k,itpn)+fmn*boundary(1)%u(j,k,itmn)
            up(2,j+1,k) = ((unext+turb(j,k)) - um(2,j+1,k))/rdt
          endif
        end do
      end do
    case(2) ! East
      do j = 1,nx1
        jpos = j + (myidy * jmax) - 1
        jpatch = int((jpos-0.5)*dy/dyint)+1
        do k = 1,nx2
          kpatch = k
          uwallcurrent = fpc*boundary(2)%u(j,k,itpc)+fmc*boundary(2)%u(j,k,itmc)
          if(uwallcurrent>=0.) then ! Outflow (Radiation)
            up(i2,j+1,k) = -min(max(boundary(2)%uphase(jpatch,kpatch),uwallcurrent),dx/rdt) * &
              (u0(i2,j+1,k)-u0(i1,j+1,k))*dxi
          else ! Inflow (Dirichlet)
            unext = fpn*boundary(2)%u(j,k,itpn)+fmn*boundary(2)%u(j,k,itmn)
            up(i2,j+1,k) = ((unext+turb(j,k)) - um(i2,j+1,k))/rdt
          endif
        end do
      end do
    case(3) ! South
      do i = 1,nx1
        ipos = i + (myidx * imax) - 1
        ipatch = int((ipos-0.5)*dx/dxint)+1
        do k = 1,nx2
          kpatch = k
          uwallcurrent = fpc*boundary(3)%v(i,k,itpc)+fmc*boundary(3)%v(i,k,itmc)
          if(uwallcurrent<=0.) then ! Outflow (Radiation)
            vp(i+1,2,k) = -max(min(boundary(3)%uphase(ipatch,kpatch),uwallcurrent),-dy/rdt) * &
              (v0(i+1,3,k)-v0(i+1,2,k))*dyi
          else ! Inflow (Dirichlet)
            unext = fpn*boundary(3)%v(i,k,itpn)+fmn*boundary(3)%v(i,k,itmn)
            vp(i+1,2,k) = ((unext+turb(i,k)) - vm(i+1,2,k))/rdt
          endif
        end do
      end do
    case(4) ! North
      do i = 1,nx1
        ipos = i + (myidx * imax) - 1
        ipatch = int((ipos-0.5)*dx/dxint)+1
        do k = 1,nx2
          kpatch = k
          uwallcurrent = fpc*boundary(4)%v(i,k,itpc)+fmc*boundary(4)%v(i,k,itmc)
          if(uwallcurrent>=0.) then ! Outflow (Radiation)
            vp(i+1,j2,k) = -min(max(boundary(4)%uphase(ipatch,kpatch),uwallcurrent),dy/rdt) * &
              (v0(i+1,j2,k)-v0(i+1,j1,k))*dyi
          else ! Inflow (Dirichlet)
            unext = fpn*boundary(4)%v(i,k,itpn)+fmn*boundary(4)%v(i,k,itmn)
            vp(i+1,j2,k) = ((unext+turb(i,k)) - vm(i+1,j2,k))/rdt
          endif
        end do
      end do
    case(5) ! Top
      do i = 1,nx1
        ipos = i + (myidx * imax) - 1
        ipatch = int((ipos-0.5)*dx/dxint)+1
        do j = 1,nx2
          jpos = j + (myidy * jmax) - 1
          jpatch = int((jpos-0.5)*dy/dyint)+1
          uwallcurrent = fpc*boundary(5)%w(i,j,itpc)+fmc*boundary(5)%w(i,j,itmc)
          if(uwallcurrent>=0.) then ! Outflow (Radiation)
            wp(i+1,j+1,k1) = -min(max(boundary(5)%uphase(ipatch,jpatch),uwallcurrent),dzf(kmax)/rdt) * &
              (rhobh(k1)*w0(i+1,j+1,k1)-rhobh(kmax)*w0(i+1,j+1,kmax))/(dzf(kmax)*rhobh(k1)) + &
              grav*(thv0h(i+1,j+1,k1)-thvh(k1))/thvh(k1)
          else ! Inflow (Dirichlet)
            unext = fpn*boundary(5)%w(i,j,itpn)+fmn*boundary(5)%w(i,j,itmn)
            wp(i+1,j+1,k1) = ((unext+turb(i,j)) - wm(i+1,j+1,k1))/rdt
          endif
        end do
      end do
    end select
  end subroutine applyboundaryh

  subroutine radcorrection(ib)
    ! Calculates the integrated mass correction term for the boundary normal
    ! velocity components
    use mpi
    use modmpi, only : comm3d,commrow,commcol,myidx,myidy,mpierr,MY_REAL
    use modglobal, only : jmax,imax,kmax,i1,j1,dx,dy,dzf,i2,j2,k1,dxi,dyi,rtimee,rdt
    use modfields, only : rhobf, up, vp, wp
    implicit none
    integer, intent(in) :: ib
    integer :: ipos,jpos,kpos,ipatch,jpatch,kpatch,i,j,k,itp,itm
    real :: sum,tp,tm,idtb,dubdt

    itm = 1
    if(ntboundary>1) then
      do while(rtimee-rdt>tboundary(itm))
        itm=itm+1
      end do
      if (rtimee-rdt>tboundary(1)) then
        itm=itm-1
      end if
      itp = itm+1
    else
      itp = 1
    endif
    tm = tboundary(itm)
    tp = tboundary(itp)
    idtb = 1./max(1e-6,tp-tm)
    select case(ib) ! Select boundary
    case(1) ! West
      ! Calculate correction term for each patch
      boundary(1)%radcorrsingle = 0.
      do j = 2,j1
        jpos = j + (myidy * jmax) - 1
        jpatch = int((jpos-0.5)*dy/dyint)+1
        do k = 1,kmax
          kpatch = k
          dubdt = (boundary(1)%u(j-1,k,itp)-boundary(1)%u(j-1,k,itm))*idtb
          boundary(1)%radcorrsingle(jpatch,kpatch) = boundary(1)%radcorrsingle(jpatch,kpatch) + &
            rhobf(k)*(-up(2,j,k)+dubdt)*dzf(k)*dy*rhointi(kpatch)/dyint
        end do
      end do
      ! Communicate integration between processes
      call MPI_ALLREDUCE(boundary(1)%radcorrsingle,boundary(1)%radcorr,nypatch*nzpatch,MY_REAL, &
                         MPI_SUM, commcol,mpierr)
      ! Apply correction term
      do j = 2,j1
        jpos = j + (myidy * jmax) - 1
        jpatch = int((jpos-0.5)*dy/dyint)+1
        do k = 1,kmax
          kpatch = k
          up(2,j,k) = up(2,j,k) + boundary(1)%radcorr(jpatch,kpatch)
        end do
      end do
    case(2) ! East
      ! Calculate correction term for each patch
      boundary(2)%radcorrsingle = 0.
      do j = 2,j1
        jpos = j + (myidy * jmax) - 1
        jpatch = int((jpos-0.5)*dy/dyint)+1
        do k = 1,kmax
          kpatch = k
          dubdt = (boundary(2)%u(j-1,k,itp)-boundary(2)%u(j-1,k,itm))*idtb
          boundary(2)%radcorrsingle(jpatch,kpatch) = boundary(2)%radcorrsingle(jpatch,kpatch) + &
            rhobf(k)*(-up(i2,j,k)+dubdt)*dzf(k)*dy*rhointi(kpatch)/dyint
        end do
      end do
      ! Communicate integration between processes
      call MPI_ALLREDUCE(boundary(2)%radcorrsingle,boundary(2)%radcorr,nypatch*nzpatch,MY_REAL, &
                         MPI_SUM, commcol,mpierr)
      ! Apply correction term
      do j = 2,j1
        jpos = j + (myidy * jmax) - 1
        jpatch = int((jpos-0.5)*dy/dyint)+1
        do k = 1,kmax
          kpatch = k
          up(i2,j,k) = up(i2,j,k) + boundary(2)%radcorr(jpatch,kpatch)
        end do
      end do
    case(3) ! South
      ! Calculate correction term for each patch
      boundary(3)%radcorrsingle= 0.
      do i = 2,i1
        ipos = i + (myidx * imax) - 1
        ipatch = int((ipos-0.5)*dx/dxint)+1
        do k = 1,kmax
          kpatch = k
          dubdt = (boundary(3)%v(i-1,k,itp)-boundary(3)%v(i-1,k,itm))*idtb
          boundary(3)%radcorrsingle(ipatch,kpatch) = boundary(3)%radcorrsingle(ipatch,kpatch) + &
            rhobf(k)*(-vp(i,2,k)+dubdt)*dzf(k)*dx*rhointi(kpatch)/dxint
        end do
      end do
      ! Communicate integration between processes
      call MPI_ALLREDUCE(boundary(3)%radcorrsingle,boundary(3)%radcorr,nxpatch*nzpatch,MY_REAL, &
                         MPI_SUM, commrow,mpierr)
      ! Apply correction term
      do i = 2,i1
        ipos = i + (myidx * imax) - 1
        ipatch = int((ipos-0.5)*dx/dxint)+1
        do k = 1,kmax
          kpatch = k
          vp(i,2,k) = vp(i,2,k) + boundary(3)%radcorr(ipatch,kpatch)
        end do
      end do
    case(4) ! North
      ! Calculate correction term for each patch
      boundary(4)%radcorrsingle = 0.
      do i = 2,i1
        ipos = i + (myidx * imax) - 1
        ipatch = int((ipos-0.5)*dx/dxint)+1
        do k = 1,kmax
          kpatch = k
          dubdt = (boundary(4)%v(i-1,k,itp)-boundary(4)%v(i-1,k,itm))*idtb
          boundary(4)%radcorrsingle(ipatch,kpatch) = boundary(4)%radcorrsingle(ipatch,kpatch) + &
            rhobf(k)*(-vp(i,j2,k)+dubdt)*dzf(k)*dx*rhointi(kpatch)/dxint
        end do
      end do
      ! Communicate integration between processes
      call MPI_ALLREDUCE(boundary(4)%radcorrsingle,boundary(4)%radcorr,nxpatch*nzpatch,MY_REAL, &
                         MPI_SUM, commrow,mpierr)
      ! Apply correction term
      do i = 2,i1
        ipos = i + (myidx * imax) - 1
        ipatch = int((ipos-0.5)*dx/dxint)+1
        do k = 1,kmax
          kpatch = k
          vp(i,j2,k) = vp(i,j2,k) + boundary(4)%radcorr(ipatch,kpatch)
        end do
      end do
    case(5) ! Top
      ! Calculate correction term for each patch
      boundary(5)%radcorrsingle = 0.
      do i = 2,i1
        ipos = i + (myidx * imax) - 1
        ipatch = int((ipos-0.5)*dx/dxint)+1
        do j = 2,j1
          jpos = j + (myidy * jmax) - 1
          jpatch = int((jpos-0.5)*dy/dyint)+1
          dubdt = (boundary(5)%w(i-1,j-1,itp)-boundary(5)%w(i-1,j-1,itm))*idtb
          boundary(5)%radcorrsingle(ipatch,jpatch) = boundary(5)%radcorrsingle(ipatch,jpatch) + &
            (-wp(i,j,k1)+dubdt)*dx*dy/(dxint*dyint)
        end do
      end do
      ! Communicate integration between processes
      call MPI_ALLREDUCE(boundary(5)%radcorrsingle,boundary(5)%radcorr,nxpatch*nypatch,MY_REAL, &
                         MPI_SUM, comm3d,mpierr)
      ! Apply correction term
      do i = 2,i1
        ipos = i + (myidx * imax) - 1
        ipatch = int((ipos-0.5)*dx/dxint)+1
        do j = 2,j1
          jpos = j + (myidy * jmax) - 1
          jpatch = int((jpos-0.5)*dy/dyint)+1
          wp(i,j,k1) = wp(i,j,k1) + boundary(5)%radcorr(ipatch,jpatch)
        end do
      end do
    end select
  end subroutine radcorrection
  subroutine handle_err(errcode)

  implicit none

  integer errcode

  write(6,*) 'Error: ', nf90_strerror(errcode)
  stop 2

  end subroutine handle_err
end module modopenboundary
