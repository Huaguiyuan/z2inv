PROGRAM z2inv
  !
  !  Author: Chao Cao @ Hangzhou Normal University (2015)
  !
  USE constants
  USE wvfunc
  !
  implicit none
  !
  character(len=32) :: arg              ! Command line argument
  integer nb                            ! Electron fill to band
  real(dp) :: xx(1:3)                   ! temporary
  integer, allocatable :: igrev(:)      ! Index of -k
  real(dp) int0, ksi(8)                 ! int0 : normal, int1: parity
  complex(dp) int1
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
  CALL read_wfc()
  !
  do ik=1, nkpt
    !
    write(stdout, *) 'Working on KPT', ik
    !
    ! Initialize
    !
    ksi(ik)=1.0d0
    allocate(igrev(npw(ik)))
    !
    ! Determine the index of -(k+G)
    !
    do ii=1, npw(ik)
      igrev(ii)=-1
      do jj=1, npw(ik)
        xx(:)=gvec(:, ii, ik)+gvec(:, jj, ik)
        if (sum(xx(:)*xx(:))< 1E-5) then
          igrev(ii)=jj
        endif
      enddo
      if (igrev(ii)<0) then
        write(*,*) "!!!! Cannot find -k for grid", ii, ": ", grid(:, ii, ik)
      endif
    enddo
    !
    ! Do the parity <psi_i|P|psi_i> for one of the Kramer pairs
    !
    do ib=1, nbnd, 2
      !
      int0=0.0
      int1=cmplx(0.d0,0.d0)
      !
      ! For spinors, there are two components
      !
      do isp=1, 2
        !
        ! int0 is the integral of the pseudo wfc
        ! int1 is <psi_i|P|psi_i>
        ! Strictly speaking, S|psi> should be calculated
        !
        int0=int0+sum(conjg(wfc(:, ib, ik, isp))*wfc(:, ib, ik, isp))
        int1=int1+sum(conjg(wfc(igrev(:), ib, ik, isp))*wfc(:, ib, ik, isp))
        !
      enddo
      !
      write(*,'(A, 1I4, A, 2F14.9)') "BAND ", (ib+1)/2, ": ", int1/int0
      !
      if ((nb>0.and.ib<nb).or.(nb==0.and.nint(occ(ib, ik)).eq.1)) &
        ksi(ik)=ksi(ik)*(int1/int0)
      !
    enddo !ib
    !
    deallocate(igrev)
    !
  enddo
  !
  write(*,*) "!!!!!!!! FINAL RESULTS !!!!!!!!"
  do ik=1, nkpt
    write(*,'(A, 1I2, A, 1F14.9)') "kpt: ", ik, " ksi: ", ksi(ik)
  enddo
  !
END PROGRAM
