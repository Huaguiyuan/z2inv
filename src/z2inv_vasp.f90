PROGRAM z2inv
  !
  use constants
  use wvfunc
  !
  implicit none
  !
  real(dp) ksi(8), xx(3)
  real(dp), allocatable :: gvec(:, :)
  real(dp) int0, int1
  integer ii, jj, ik, ib, jj1, jj2
  integer, allocatable :: igrev(:)      ! -k list
  character(len=32) :: arg              ! Command line argument
  integer nb                            ! Electron fill to band
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
      gvec(:, ii)=(kvec(1, ik)+gkvec(1, ii, ik))*bvec(:, 1) + &
                  (kvec(2, ik)+gkvec(2, ii, ik))*bvec(:, 2) + &
                  (kvec(3, ik)+gkvec(3, ii, ik))*bvec(:, 3)
      !
    enddo
    !
    igrev(:)=-1
    !
    do ii=0, nplane(ik)-1
      if (ii.lt.nplane(ik)/2) then
        jj1=0
        jj2=nplane(ik)/2-1
      else
        jj1=nplane(ik)/2
        jj2=nplane(ik)-1
      endif
      do jj=jj1, jj2
        xx(:)=gvec(:, ii+1)+gvec(:, jj+1)
        xx(:)=xx(:)-nint(xx(:))
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
      !
      int0=sum(conjg(coeff(:, ib, ik))*coeff(:, ib, ik))
      !
      int1=0.d0
      !
      do ii=1, nplane(ik)
        int1=int1+conjg(coeff(igrev(ii), ib, ik))*coeff(ii, ib, ik)
      enddo
      !
      write(*, '(A, 1I4, A, 1F14.9)') "BAND ", (ib+1)/2, ": ", int1/int0
      !
      if ((nb>0.and.ib<nb).or.(nb==0.and.nint(occ(ib, ik)).eq.1)) &
         ksi(ik)=ksi(ik)*(int1/int0)
      !
    enddo
    !
    deallocate(igrev, gvec)
    !
  enddo ! ik
  !
  write(*,*) "!!!!!!!! FINAL RESULTS !!!!!!!!"
  do ik=1, nkpt
    write(*,'(A, 1I2, A, 1F14.9)') "kpt: ", ik, " ksi: ", ksi(ik)
  enddo
  !
END PROGRAM


