!include 'lapack.f90'

PROGRAM z2inv
  !
  use constants
  use wvfunc
  !use lapack95,  only : heev
  !
  implicit none
  !
  real(dp) ksi(8), xx(3)
  real(dp), allocatable :: gvec(:, :)
  real(dp) int0
  complex(dp) int1
  !complex(dp),allocatable :: int1(2, 2)
  !real(dp),allocatable :: ksi_b(2)
  integer ii, jj, ik, ib, jb, jj1, jj2
  integer, allocatable :: igrev(:)      ! -k list
  character(len=32) :: arg              ! Command line argument
  integer nb                            ! Electron fill to band
  !integer info
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
  CALL read_wavecar
  !
  !allocate(int1(nbnd/2, nbnd/2), ksi_b(nbnd/2))
  !
  do ik=1, nkpt ! should be <=8
    !
    write(stdout, *) 'Working on KPT: ', ik
    write(stdout, '(3F14.9)') kvec(:, ik)
    !
    ksi(ik)=1.0d0
    !
    allocate(gvec(3, nplane(ik)))
    allocate(igrev(nplane(ik)))
    !
    do ii=1, nplane(ik)
      !
      gvec(:, ii)=kvec(:, ik)+gkvec(:, ii, ik)
      !
    enddo
    !
    igrev(:)=-1
    !
    do ii=0, nplane(ik)-1
      !
      ! For spin unpolarized
      !
      if (ii.lt.nplane(ik)/2) then
        jj1=0
        jj2=nplane(ik)/2-1
      else
        jj1=nplane(ik)/2
        jj2=nplane(ik)-1
      endif
      !
      ! For AFM
      !
      !if (ii.lt.nplane(ik)/2) then
      !  jj1=nplane(ik)/2
      !  jj2=nplane(ik)-1
      !else
      !  jj1=0
      !  jj2=nplane(ik)/2-1
      !endif
      !
      do jj=jj1, jj2
        xx(:)=gvec(:, ii+1)+gvec(:, jj+1)
        if (sum(xx(:)**2)<eps5) then
          igrev(ii+1)=jj+1
        endif
      enddo
      if (igrev(ii)<0) then
        write(*, *) '!!!! Cannot find -k for grid', ii+1, ":", gkvec(:, ii+1, ik)
        write(*, *) '!!!!  non TR k?'
        stop
      endif
    enddo
    !
    do ib=1, nbnd, 2
      int1=cmplx(0.d0,0.d0)
      !
      int0=sum(conjg(coeff(:, ib, ik))*coeff(:, ib, ik))
      !
      do ii=1, nplane(ik)
        int1=int1+conjg(coeff(igrev(ii), ib, ik))*coeff(ii, ib, ik)
      enddo
      !
      write(*, '(A, 1I4, A, 2F14.9)') "BAND ", (ib+1)/2, ": ", int1/int0
      !
      if ((nb>0.and.ib<nb).or.(nb==0.and.nint(occ(ib, ik)).eq.1)) &
         ksi(ik)=ksi(ik)*(int1/int0)
      !
    enddo
    !
    !CALL heev(int1, ksi_b, 'V', 'U', info)
    !
    !do ib=1, nbnd/2
    !  write(*, '(A, 1I4, A, 2F14.9)') "BAND ", ib, ": ", ksi_b(ib)
    !  if ((nb>0.and.ib<nb/2).or.(nb==0.and.nint(occ(ib, ik)).eq.1)) &
    !     ksi(ik)=ksi(ik)*ksi_b(ib)
    !enddo
    !
    deallocate(igrev, gvec)
    !
  enddo ! ik
  !
  !deallocate(int1, ksi_b)
  !
  write(*,*) "!!!!!!!! FINAL RESULTS !!!!!!!!"
  do ik=1, nkpt
    write(*,'(A, 1I2, A, 1F14.9)') "kpt: ", ik, " ksi: ", ksi(ik)
  enddo
  !
END PROGRAM


