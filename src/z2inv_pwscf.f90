PROGRAM z2inv
  !
  !  Author: Chao Cao @ Hangzhou Normal University (2015)
  !
  USE iotk_module
  USE kinds,            ONLY : DP
  USE becmod,           ONLY : bec_type, becp, calbec, &
                               allocate_bec_type, deallocate_bec_type
  USE io_files,         ONLY : prefix, iunwfc, nwordwfc, pseudo_dir, psfile
  !
  implicit none
  !
  integer, parameter :: iun=10
  character(len=32) :: arg              ! Command line argument
  integer nb                            ! Electron fill to band
  integer  npw, ipw                     ! npw: number of gvecs for the kpt
  character(iotk_attlenx) :: attr
  real(dp) :: bvec(1:3, 1:3)            ! b1, b2, b3
  real(dp) :: xk(1:3)                   ! Kpt coordinate
  real(dp) :: xx(1:3)                   ! temporary
  real(dp), allocatable :: gvec(:,:)    ! gvec coordinate
  integer, allocatable :: grid(:)       ! gvec grid
  integer, allocatable :: rgidx(:)      ! The index of -k
  real(dp) int0, int1, ksi(8)           ! int0 : normal, int1: parity
  complex(dp), allocatable :: wfc(:)    ! Wave function (pseudo)
  complex(dp), allocatable :: swfc(:)   ! S|psi> useless here?
  !
  integer ik, ib, isp, ii, jj           ! local variables
  !
  if (iargc()<1) then
    write(*,*) "!!! Did you read the help message in the code?"
    stop
  endif
  !
  call getarg(1, arg)
  !
  read(arg,*) nb
  !
  if (mod(nb, 2)==1) then
    write(*,*) "Number of bands is incorrect", nb
    stop
  endif
  !
  call iotk_open_read(iun, "data-file.xml")
  !
  !  READ b1, b2, b3 (for constructing gvec later)
  !
  call iotk_scan_begin(iun, "CELL")
  call iotk_scan_begin(iun, "RECIPROCAL_LATTICE_VECTORS")
  call iotk_scan_dat(iun, "b1", bvec(1,:))
  call iotk_scan_dat(iun, "b2", bvec(2,:))
  call iotk_scan_dat(iun, "b3", bvec(3,:))
  call iotk_scan_end(iun, "RECIPROCAL_LATTICE_VECTORS")
  call iotk_scan_end(iun, "CELL")
  !
  do ik=1, 8
    !
    ! Initialize ksi
    !
    ksi(ik)=1.0d0
    !
    ! READ K-POINT coordinate (in unit of 2piba)
    !
    call iotk_scan_begin(iun, "EIGENVALUES")
    call iotk_scan_begin(iun, "K-POINT"//iotk_index(ik))
    call iotk_scan_dat(iun, "K-POINT_COORDS", xk)
    call iotk_scan_end(iun, "K-POINT"//iotk_index(ik))
    call iotk_scan_end(iun, "EIGENVALUES")
    !
    write (*, '(A,I2,A,3F12.9)' ) "Kpoint #", ik, ": ", xk(1:3)
    !
    ! READ Specification of wave functions (pseudo)
    !
    call iotk_scan_begin(iun, "EIGENVECTORS")
    call iotk_scan_begin(iun, "K-POINT"//iotk_index(ik))
    call iotk_scan_dat(iun, "NUMBER_OF_GK-VECTORS", npw)
    allocate(grid(1:npw*3), gvec(1:npw, 1:3), wfc(1:npw), rgidx(1:npw))
    call iotk_scan_begin(iun, "GK-VECTORS")
    call iotk_scan_dat(iun, "GRID", grid(1:npw*3))
    call iotk_scan_end(iun, "GK-VECTORS")
    !
    ! Find out the real coordinate of gvecs
    !
    do ii=1, npw
        gvec(ii,:)=xk(:)+bvec(1,:)*grid(ii*3-2)+bvec(2,:)*grid(ii*3-1)+bvec(3,:)*grid(ii*3)
    enddo
    !
    ! Determine the index of -k
    !
    do ii=0, npw-1
      rgidx(ii+1)=-1
      do jj=0, npw-1
        xx(:)=gvec(ii+1, :)+gvec(jj+1, :)
        xx(:)=xx(:)-nint(xx(:))
        if (sum(xx(:)*xx(:))< 1E-5) then
          rgidx(ii+1)=jj+1
        endif
      enddo
      if (rgidx(ii+1)<0) then
        write(*,*) "!!!! Cannot find -k for grid", ii+1, ": ", grid(ii*3+1:ii*3+3)
      endif
    enddo
    !
    ! Do the parity <psi_i|P|psi_i> for one of the Kramer pairs
    !
    do ib=1, nb, 2
      !
      int0=0.0
      int1=0.0
      !
      ! For spinors, there are two components
      !
      do isp=1, 2
        !
        call iotk_scan_begin(iun, "WFC"//iotk_index(isp))
        call iotk_scan_dat(iun, "evc"//iotk_index(ib), wfc(1:npw))
        call iotk_scan_end(iun, "WFC"//iotk_index(isp))
        !
        ! int0 is the integral of the pseudo wfc
        ! int1 is <psi_i|P|psi_i>
        ! Strictly speaking, S|psi> should be calculated
        !
        int0=int0+sum(conjg(wfc(:))*wfc(:))
        int1=int1+sum(conjg(wfc(rgidx(:)))*wfc(:))
        !
      enddo
      !
      write(*,'(A, 1I4, A, 1F14.9)') "BAND ", (ib+1)/2, ": ", int1/int0
      !
      ! The actual index should only be calculated for 
      !   one of the Kramer pairs
      !
      ksi(ik)=ksi(ik)*(int1/int0)
      !
    enddo
    !
    deallocate(grid, wfc, rgidx, gvec)
    !
    call iotk_scan_end(iun, "K-POINT"//iotk_index(ik))
    call iotk_scan_end(iun, "EIGENVECTORS")
    !
  enddo
  !
  write(*,*) "!!!!!!!! FINAL RESULTS !!!!!!!!"
  do ik=1, 8
    write(*,'(A, 1I2, A, 1F14.9)') "kpt: ", ik, " ksi: ", ksi(ik)
  enddo
  !
  call iotk_close_read(iun)
  !
END PROGRAM
