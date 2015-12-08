MODULE wvfunc
  !
  use constants
  !
  implicit none
  !
  integer nsp           ! # of spin-polarization
  integer nkpt          ! # of kpoints
  integer nbnd          ! # of bands
  integer npmax         ! Max # of planewaves
  real(dp) omega        ! Volume of cell
  real(dp) ecut         ! Energy cutoff
  real(dp), dimension(3, 3) :: alat, bvec       ! Real/Reciprocal lattice
  integer, allocatable :: nplane(:)             ! # of planewaves nplane(ik)
  integer, allocatable :: gkvec(:, :, :)        ! Gvec of planewaves gkvec(3, np, ik)
  real(dp), allocatable :: occ(:, :)            ! Occupations occ(ib, ik)
  real(dp), allocatable :: kvec(:, :)           ! K points kvec(3, ik)
  complex*8, allocatable :: coeff(:, :, :)    ! Coefficients of WFC coeff(ig, ib, ik)
  !
 CONTAINS
  !
  SUBROUTINE read_wavecar()
    !
    use constants
    !
    implicit none
    !
    real(dp) pi
    integer nrecl
    real(dp) tmp1, tmp2
    complex(dp) ctmp
    real(dp), dimension(3) :: btmp
    real(dp), dimension(3) :: bvl
    integer nb1max, nb2max, nb3max
    integer nb1A, nb2A, nb3A
    integer nb1B, nb2B, nb3B
    integer nb1C, nb2C, nb3C
    integer ii, ik, ib, ig3, ig2, ig1, ig3p, ig2p, ig1p
    !
    pi=4.d0*atan(1.d0)
    !
    ! Open WAVECAR and read in the recl parameter
    nrecl=16
    open(unit=fin, file='WAVECAR', access='direct', recl=nrecl)
    read(unit=fin, rec=1) tmp1, tmp2
    close(unit=fin)
    nrecl=nint(tmp1)
    nsp=nint(tmp2)
    !
    ! Use the correct recl parameter and reopen WAVECAR
    open(unit=fin, file='WAVECAR', access='direct', recl=nrecl)
    read(unit=fin, rec=2) tmp1, tmp2, ecut, alat, bvec
    nkpt=nint(tmp1)
    nbnd=nint(tmp2)
    !
    write(stdout, *) '##### WAVECAR INFO #####'
    write(stdout, *) ' nkpt: ', nkpt
    write(stdout, *) ' nbnd: ', nbnd
    write(stdout, '(A, 1F14.9, A)') ' Energy cutoff: ', ecut, 'eV'
    !
    if ((nsp.ne.1).or.(nkpt.ne.8)) then
      write(stdout, *) ' WARNNING!!!'
      write(stdout, *) '  Make sure your WAVECAR file is for Z2 calculation'
    endif
    !
    allocate(nplane(nkpt), occ(nbnd, nkpt), kvec(3, nkpt))
    !
    call cross_prod(bvec(:, 1), alat(:, 2), alat(:, 3))
    call cross_prod(bvec(:, 2), alat(:, 3), alat(:, 1))
    call cross_prod(bvec(:, 3), alat(:, 1), alat(:, 2))
    !
    omega=dot_product(alat(:, 1), bvec(:, 1))
    !
    bvec(:, :)=2.d0*pi*bvec(:, :)/omega
    !
    write(stdout, '(A, 1F14.9)') ' System Volume: ', omega
    write(stdout, *) ' Reciprocal lattice: '
    do ii=1, 3
      write(stdout, '(3F14.9)') bvec(:, ii)
    enddo
    !
    do ii=1, 3
      bvl(ii)=dsqrt(sum(bvec(:, ii)**2))
    enddo
    !
    tmp1=acos(dot_product(bvec(:, 1), bvec(:, 2))/(bvl(1)*bvl(2)))
    call cross_prod(btmp, bvec(:, 1), bvec(:, 2))
    tmp2=dot_product(btmp, bvec(:, 3))/(bvl(3)*dsqrt(sum(btmp(:)**2)))
    nb1A=(dsqrt(ecut*c)/(bvl(1)*abs(sin(tmp1))))+1
    nb2A=(dsqrt(ecut*c)/(bvl(2)*abs(sin(tmp1))))+1
    nb3A=(dsqrt(ecut*c)/(bvl(3)*abs(tmp2)))+1
    !
    write(*,*) nb1A, nb2A, nb3A
    !
    tmp1=acos(dot_product(bvec(:, 1), bvec(:, 3))/(bvl(1)*bvl(3)))
    call cross_prod(btmp, bvec(:, 1), bvec(:, 3))
    tmp2=dot_product(btmp, bvec(:, 2))/(bvl(2)*dsqrt(sum(btmp(:)**2)))
    nb1B=(dsqrt(ecut*c)/(bvl(1)*abs(sin(tmp1))))+1
    nb2B=(dsqrt(ecut*c)/(bvl(2)*abs(tmp2)))+1
    nb3B=(dsqrt(ecut*c)/(bvl(3)*abs(sin(tmp1))))+1
    !
    write(*,*) nb1B, nb2B, nb3B
    !
    tmp1=acos(dot_product(bvec(:, 2), bvec(:, 3))/(bvl(2)*bvl(3)))
    call cross_prod(btmp, bvec(:, 2), bvec(:, 3))
    tmp2=dot_product(btmp, bvec(:, 1))/(bvl(1)*dsqrt(sum(btmp(:)**2)))
    nb1C=(dsqrt(ecut*c)/(bvl(1)*abs(tmp2)))+1
    nb2C=(dsqrt(ecut*c)/(bvl(2)*abs(sin(tmp1))))+1
    nb3C=(dsqrt(ecut*c)/(bvl(3)*abs(sin(tmp1))))+1
    !
    write(*,*) nb1C, nb2C, nb3C
    !
    nb1max=max0(nb1A, nb1B, nb1C)
    nb2max=max0(nb2A, nb2B, nb2C)
    nb3max=max0(nb3A, nb3B, nb3C)
    !
    nb1A=nint(4.d0*pi*nb1A*nb2A*nb3A/3.d0)
    nb1B=nint(4.d0*pi*nb1B*nb2B*nb3B/3.d0)
    nb1C=nint(4.d0*pi*nb1C*nb2C*nb3C/3.d0)
    !
    npmax=2*min0(nb1A, nb1B, nb1C)
    !
    write(stdout, *) 'Maximum number of planewaves: ', npmax
    !
    allocate(gkvec(3, npmax, nkpt), coeff(npmax, nbnd, nkpt))
    !
    coeff(:, :, :)=cmplx(0.0, 0.0)
    !
    do ik=1, nkpt
      !
      ii=3+(ik-1)*(nbnd+1)
      read(unit=fin, rec=ii) tmp1, kvec(:, ik), &
        (ctmp, occ(ib, ik), ib=1, nbnd)
      nplane(ik)=nint(tmp1)
      write(stdout, *) '   Reading KPT: ', ik
      !
      ! Find Gvec and check
      !
      ii=0
      do ig3=0, 2*nb3max
        if (ig3.gt.nb3max) then
          ig3p=ig3-2*nb3max-1
        else
          ig3p=ig3
        endif
        do ig2=0, 2*nb2max
          if (ig2.gt.nb2max) then
            ig2p=ig2-2*nb2max-1
          else
            ig2p=ig2
          endif
          do ig1=0, 2*nb1max
            if (ig1.gt.nb1max) then
              ig1p=ig1-2*nb1max-1
            else
              ig1p=ig1
            endif
            btmp(:)=(kvec(1, ik)+ig1p)*bvec(:, 1) + &
                    (kvec(2, ik)+ig2p)*bvec(:, 2) + &
                    (kvec(3, ik)+ig3p)*bvec(:, 3)
            if (sum(btmp(:)**2).lt.(ecut*c)) then
              ii=ii+1
              gkvec(1, ii, ik)=ig1p
              gkvec(2, ii, ik)=ig2p
              gkvec(3, ii, ik)=ig3p
            endif
          enddo ! ig1
        enddo ! ig2
      enddo ! ig3
      !
      if (2*ii.ne.nplane(ik)) then
        write(stdout, *) '!!! ERROR Wrong # of planewaves!'
        write(stdout, *) 'nplane: ', nplane(ik), 'computed: ', 2*ii
        stop
      endif
      !
      gkvec(:, ii+1:2*ii, ik)=gkvec(:, 1:ii, ik)
      !
      do ib=1, nbnd
        read(unit=fin, rec=(ik-1)*(nbnd+1)+ib+3) (coeff(ii, ib, ik), ii=1, nplane(ik))
      enddo
      !
    enddo
    !
    close(unit=fin)
    !
  END SUBROUTINE
  !
  SUBROUTINE finalize_wvfunc()
    !
    implicit none
    !
    if (allocated(gkvec)) deallocate(gkvec)
    if (allocated(coeff)) deallocate(coeff)
    if (allocated(nplane)) deallocate(nplane)
    if (allocated(kvec)) deallocate(kvec)
    if (allocated(occ)) deallocate(occ)
    !
  END SUBROUTINE
  !
END MODULE
