MODULE wvfunc
  !
  use constants,        ONLY: dp
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
    use constants,      ONLY : dp, twopi, cmplx_0, cc, stdout, fin
    !
    integer nrecl
    real(dp) tmp1, tmp2, ecut, omega
    real(dp) b1, b2, b3
    real(dp),dimension(3) :: btmp
    real(dp),dimension(3,3) :: alat
    complex(dp),allocatable :: ctmp(:)
    complex*8,allocatable :: wftmp(:)
    integer nb1A, nb1B, nb1C, nb2A, nb2B, nb2C, nb3A, nb3B, nb3C
    integer nb1max, nb2max, nb3max
    integer ii, ik, ib, isp
    integer ig1, ig2, ig3, ig1p, ig2p, ig3p
    !
    ! Open WAVECAR and read in the recl parameter
    !
    nrecl=16
    open(unit=fin, file='WAVECAR', access='direct', recl=nrecl)
    read(unit=fin, rec=1) tmp1, tmp2
    close(unit=fin)
    !
    nrecl=nint(tmp1)
    nsp=nint(tmp2)
    !
    nsp=2               ! VASP uses nsp=1 for noncollinear
    !
    ! Use the correct recl parameter and reopen WAVECAR
    !
    open(unit=fin, file='WAVECAR', access='direct', recl=nrecl)
    read(unit=fin, rec=2) tmp1, tmp2, ecut, alat(:, :)
    nkpt=nint(tmp1)
    nbnd=nint(tmp2)
    !
    allocate(npw(nkpt), occ(nbnd, nkpt), xk(3, nkpt))
    !
    call cross_prod(bvec(:, 1), alat(:, 2), alat(:, 3))
    call cross_prod(bvec(:, 2), alat(:, 3), alat(:, 1))
    call cross_prod(bvec(:, 3), alat(:, 1), alat(:, 2))
    !
    omega=dot_product(alat(:, 1), bvec(:, 1))
    !
    bvec(:, :)=twopi*bvec(:, :)/omega
    !
    b1=dsqrt(sum(bvec(:, 1)**2))
    b2=dsqrt(sum(bvec(:, 2)**2))
    b3=dsqrt(sum(bvec(:, 3)**2))
    !
    tmp1=acos(dot_product(bvec(:, 1), bvec(:, 2))/(b1*b2))
    call cross_prod(btmp, bvec(:, 1), bvec(:, 2))
    tmp2=dot_product(btmp, bvec(:, 3))/(b3*dsqrt(sum(btmp(:)**2)))
    nb1A=(dsqrt(ecut*cc)/(b1*abs(sin(tmp1))))+1
    nb2A=(dsqrt(ecut*cc)/(b2*abs(sin(tmp1))))+1
    nb3A=(dsqrt(ecut*cc)/(b3*abs(tmp2)))+1
    !
    tmp1=acos(dot_product(bvec(:, 1), bvec(:, 3))/(b1*b3))
    call cross_prod(btmp, bvec(:, 1), bvec(:, 3))
    tmp2=dot_product(btmp, bvec(:, 2))/(b2*dsqrt(sum(btmp(:)**2)))
    nb1B=(dsqrt(ecut*cc)/(b1*abs(sin(tmp1))))+1
    nb2B=(dsqrt(ecut*cc)/(b2*abs(tmp2)))+1
    nb3B=(dsqrt(ecut*cc)/(b3*abs(sin(tmp1))))+1
    !
    tmp1=acos(dot_product(bvec(:, 2), bvec(:, 3))/(b2*b3))
    call cross_prod(btmp, bvec(:, 2), bvec(:, 3))
    tmp2=dot_product(btmp, bvec(:, 1))/(b1*dsqrt(sum(btmp(:)**2)))
    nb1C=(dsqrt(ecut*cc)/(b1*abs(tmp2)))+1
    nb2C=(dsqrt(ecut*cc)/(b2*abs(sin(tmp1))))+1
    nb3C=(dsqrt(ecut*cc)/(b3*abs(sin(tmp1))))+1
    !
    nb1max=max0(nb1A, nb1B, nb1C)
    nb2max=max0(nb2A, nb2B, nb2C)
    nb3max=max0(nb3A, nb3B, nb3C)
    !
    nb1A=nint(2.d0*twopi*nb1A*nb2A*nb3A/3.d0)
    nb1B=nint(2.d0*twopi*nb1B*nb2B*nb3B/3.d0)
    nb1C=nint(2.d0*twopi*nb1C*nb2C*nb3C/3.d0)
    !
    npwmax=min0(nb1A, nb1B, nb1C)
    !
    allocate(grid(3, npwmax, nkpt), wfc(npwmax, nbnd, nkpt, nsp))
    allocate(gvec(3, npwmax, nkpt))
    !
    wfc(:, :, :, :)=cmplx_0
    !
    do ik=1, nkpt
      !
      allocate(ctmp(nbnd))
      ii=3+(ik-1)*(nbnd+1)
      read(unit=fin, rec=ii) tmp1, xk(:, ik), &
        (ctmp(ib), occ(ib, ik), ib=1, nbnd)
      eig(:, ik)=real(ctmp(:))
      deallocate(ctmp)
      !
      npw(ik)=nint(tmp1)/2              ! Noncollinear stores
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
          !
          do ig1=0, 2*nb1max
            !
            if (ig1.gt.nb1max) then
              ig1p=ig1-2*nb1max-1
            else
              ig1p=ig1
            endif
            !
            btmp(:)=(xk(1, ik)+ig1p)*bvec(:, 1) + &
                    (xk(2, ik)+ig2p)*bvec(:, 2) + &
                    (xk(3, ik)+ig3p)*bvec(:, 3)
            !
            if (sum(btmp(:)**2).lt.(ecut*cc)) then
              ii=ii+1
              grid(1, ii, ik)=ig1p
              grid(2, ii, ik)=ig2p
              grid(3, ii, ik)=ig3p
              gvec(:, ii, ik)=xk(:, ik)+grid(:, ii, ik)
            endif
          enddo ! ig1
        enddo ! ig2
      enddo ! ig3
      !
      if (ii.ne.npw(ik)) then
        write(stdout, *) '!!! ERROR Wrong # of planewaves!'
        write(stdout, *) 'npw: ', npw(ik), 'computed: ', ii
        stop
      endif
      !
      allocate(wftmp(npw(ik)*2))
      do ib=1, nbnd
        read(unit=fin, rec=(ik-1)*(nbnd+1)+ib+3) wftmp(:)
        wfc(1:npw(ik), ib, ik, 1)=wftmp(1:npw(ik))
        wfc(1:npw(ik), ib, ik, 2)=wftmp(npw(ik)+1:2*npw(ik))
      enddo
      deallocate(wftmp)
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
  SUBROUTINE cross_prod(v1, v2, v3)
    !
    use constants, only :dp
    !
    implicit none
    !
    real(dp), dimension(3) :: v1, v2, v3
    !
    v1(1)=v2(2)*v3(3)-v2(3)*v3(2)
    v1(2)=v2(3)*v3(1)-v2(1)*v3(3)
    v1(3)=v2(1)*v3(2)-v2(2)*v3(1)
    !
  END SUBROUTINE
  !
END MODULE
