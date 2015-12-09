MODULE wvfunc
  !
  USE constants,  ONLY : dp
  !
  implicit none
  !
  integer nbnd          ! Total number of bands
  integer nkpt          ! Total number of K-points
  integer nsp           ! Number of spinors (should be 2)
  integer npwmax        ! Maximum number of plane-waves
  integer,allocatable :: npw(:)         ! # of pl wvs per ik
  !
  real(dp),dimension(3,3) :: bvec       ! Reciprocal lattice
  real(dp),allocatable :: xk(:, :)      ! Kpt coordinates
  real(dp),allocatable :: eig(:, :)     ! Eigenvalues
  real(dp),allocatable :: occ(:, :)     ! Occupations
  real(dp),allocatable :: gvec(:, :, :) ! Actual Gvectors...
  !
  integer,allocatable :: grid(:, :, :)  ! Gvectors... grid(1:3, ipw, ik)
  complex(dp),allocatable :: wfc(:, :, :, :)    ! Wavefunction |psi> (Pseudo)
  complex(dp),allocatable :: swfc(:, :, :, :)   ! S|psi>
  !
 CONTAINS
  !
  SUBROUTINE read_wfc()
    !
    USE constants,  ONLY : dp, fin, stdout, cmplx_0
    USE iotk_module
    !
    implicit none
    !
    logical lspinorb
    !
    integer ik, ib, isp, ii
    !
    CALL iotk_open_read(fin, "data-file.xml")
    !
    !  Read in basic infos (dimension, etc)
    !
    write(stdout, *) " # Reading basic info and dimensions..."
    !
    CALL iotk_scan_begin(fin, "CELL")
    CALL iotk_scan_begin(fin, "RECIPROCAL_LATTICE_VECTORS")
    CALL iotk_scan_dat(fin, "b1", bvec(:,1))
    CALL iotk_scan_dat(fin, "b2", bvec(:,2))
    CALL iotk_scan_dat(fin, "b3", bvec(:,3))
    CALL iotk_scan_end(fin, "RECIPROCAL_LATTICE_VECTORS")
    CALL iotk_scan_end(fin, "CELL")
    CALL iotk_scan_begin(fin, "SPIN")
    CALL iotk_scan_dat(fin, "SPIN-ORBIT_CALCULATION", lspinorb)
    CALL iotk_scan_dat(fin, "SPINOR_DIM", nsp)
    CALL iotk_scan_end(fin, "SPIN")
    CALL iotk_scan_begin(fin, "BRILLOUIN_ZONE")
    CALL iotk_scan_dat(fin, "NUMBER_OF_K-POINTS", nkpt)
    CALL iotk_scan_end(fin, "BRILLOUIN_ZONE")
    CALL iotk_scan_begin(fin, "BAND_STRUCTURE_INFO")
    CALL iotk_scan_dat(fin, "NUMBER_OF_BANDS", nbnd)
    CALL iotk_scan_end(fin, "BAND_STRUCTURE_INFO")
    CALL iotk_scan_begin(fin, "EIGENVECTORS")
    CALL iotk_scan_dat(fin, "MAX_NUMBER_OF_GK-VECTORS", npwmax)
    CALL iotk_scan_end(fin, "EIGENVECTORS")
    !
    ! Check
    !
    if (.not.lspinorb) then
      write(stdout, *) "  Fatal ERROR: Z2 calculation requires SOC!"
      stop
    else if (nsp.ne.2) then
      write(stdout, *) "  Fatal ERROR: Weird, spinor #:", nsp
      stop
    endif
    !
    ! allocate arrays...
    !
    allocate(npw(nkpt), xk(3, nkpt))
    allocate(eig(nbnd, nkpt), occ(nbnd, nkpt))
    allocate(grid(3, npwmax, nkpt), gvec(3, npwmax, nkpt))
    allocate(wfc(npwmax, nbnd, nkpt, nsp))
    !
    grid(:, :, :)=0
    gvec(:, :, :)=0.d0
    wfc(:, :, :, :)=cmplx_0
    !
    ! Now read in the actual data...
    !
    do ik=1, nkpt
      !
      write(stdout, *) "   # Reading KPT: ", ik
      !
      CALL iotk_scan_begin(fin, "EIGENVALUES")
      CALL iotk_scan_begin(fin, "K-POINT"//iotk_index(ik))
      !
      CALL iotk_scan_dat(fin, "K-POINT_COORDS", xk(:, ik))
      !
      CALL iotk_scan_begin(fin, "DATAFILE")
      CALL iotk_scan_dat(fin, "EIGENVALUES", eig(:, ik))
      CALL iotk_scan_dat(fin, "OCCUPATIONS", occ(:, ik))
      CALL iotk_scan_end(fin, "DATAFILE")
      !
      CALL iotk_scan_end(fin, "K-POINT"//iotk_index(ik))
      CALL iotk_scan_end(fin, "EIGENVALUES")
      !
      CALL iotk_scan_begin(fin, "EIGENVECTORS")
      CALL iotk_scan_begin(fin, "K-POINT"//iotk_index(ik))
      CALL iotk_scan_dat(fin, "NUMBER_OF_GK-VECTORS", npw(ik))
      CALL iotk_scan_begin(fin, "GK-VECTORS")
      CALL iotk_scan_dat(fin, "GRID", grid(:, 1:npw(ik), ik))
      CALL iotk_scan_end(fin, "GK-VECTORS")
      !
      ! Get the actual (cartesian) coordinates of Gvecs
      !
      do ii=1, npw(ik)
        gvec(:, ii, ik)=xk(:, ik)+bvec(:, 1)*grid(1, ii, ik)+ &
                                  bvec(:, 2)*grid(2, ii, ik)+ &
                                  bvec(:, 3)*grid(3, ii, ik)
      enddo
      !
      do ib=1, nbnd
        !
        ! Two spinors
        !
        do isp=1, nsp
          !
          call iotk_scan_begin(fin, "WFC"//iotk_index(isp))
          call iotk_scan_dat(fin, "evc"//iotk_index(ib), wfc(1:npw(ik), ib, ik, isp))
          call iotk_scan_end(fin, "WFC"//iotk_index(isp))
          !
        enddo ! isp
        !
      enddo ! ib
      !
      CALL iotk_scan_end(fin, "K-POINT"//iotk_index(ik))
      CALL iotk_scan_end(fin, "EIGENVECTORS")
      !
    enddo ! ik
    !
  END SUBROUTINE
  !
  SUBROUTINE finalize_wvfunc()
    !
    implicit none
    !
    if (allocated(npw)) deallocate(npw)
    if (allocated(xk)) deallocate(xk)
    if (allocated(eig)) deallocate(eig)
    if (allocated(occ)) deallocate(occ)
    if (allocated(gvec)) deallocate(gvec)
    if (allocated(grid)) deallocate(grid)
    if (allocated(wfc)) deallocate(wfc)
    if (allocated(swfc)) deallocate(swfc)
    !
  END SUBROUTINE
  !
END MODULE
