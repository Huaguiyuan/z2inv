include 'lapack.f90'

PROGRAM z2inv
  !
  !  Author: Chao Cao @ Hangzhou Normal University (2015)
  !
  USE constants
  USE wvfunc
  USE parity
  USE lapack95,  only : heev
  !
  implicit none
  !
  character(len=32) :: arg              ! Command line argument
  integer nb                            ! Electron fill to band
  real(dp) ksi(8)                       ! ksi: z2 number @ KPTs
  real(dp),allocatable :: ksi_b(:)      ! ksi_b: parity for each bands
  integer ngrp                          ! # of Groups of degenerate eigenstates
  integer,allocatable :: grp(:)         ! Index of the group states
  complex(dp),allocatable :: p(:, :)    ! Parity matrix
  real(dp),allocatable :: tmp(:)        ! Temporary
  !
  integer ik, ib, jb, ii, nn, info      ! local variables
  !
  if (iargc()<1) then
    !  Use electron filling to determine
    nb=0
    !
  else
    !
    CALL getarg(1, arg)
    read(arg,*) nb
    !
  endif
  !
  if (mod(nb, 2)==1) then
    write(*,*) "Number of bands is incorrect", nb
    write(*,*) "It must be even due to Kramer degeneracy."
    stop
  endif
  !
  CALL read_wfc()
  !
  allocate(ksi_b(nbnd/2), grp(nbnd/2+1))
  !
  do ik=1, nkpt
    !
    write(*,*) " # Working on KPT", ik
    !
    ! Initialize
    !
    ksi(ik)=1.0d0
    !
    ! Determine the index of -(k+G)
    !
    CALL invglist(ik)
    !
    CALL groupstates(grp, ngrp, ik)
    !
    ! For all the groups...
    !
    do ii=1, ngrp-1
      !
      ! Dimension of extra symmetry
      !
      nn=(grp(ii+1)-grp(ii))/2
      !
      if (nn.eq.1) then
        !
        ! No extra symmetry
        !
        allocate(p(nn,nn))
        !
        CALL calc_pmat(p(1,1), grp(ii), grp(ii), ik)
        !
        ib=(grp(ii)+1)/2
        ksi_b(ib)=real(p(1,1))
        !
        deallocate(p)
        !
      else
        !
        ! With extra symmetry, P may not be good quantum number
        !   So we have to do the diagonalization
        !
        allocate(p(nn*2,nn*2))
        !
        do ib=0, 2*nn-1
          !
          ! Actual band index: grp(ii)+ib
          do jb=0, 2*nn-1
            !
            ! grp(ii)+jb
            CALL calc_pmat(p(ib+1, jb+1), grp(ii)+ib, grp(ii)+jb, ik)
            !
          enddo
          !
        enddo
        !
        ib=(grp(ii)-1)/2
        !
        allocate(tmp(2*nn))
        !
        CALL heev(p, tmp, 'N', 'U', info)
        !
        do jb=1, nn
          ksi_b(ib+jb)=tmp(2*jb-1)
        enddo
        !
        deallocate(tmp)
        !
        deallocate(p)
        !
      endif
      !
    enddo ! ii (group)
    !
    do ib=1, nbnd/2
      !
      write(*,'(A, 1I4, A, 2F14.9)') "BAND ", ib, ": ", ksi_b(ib)
      !
      if ((nb>0.and.ib*2-1<nb).or.(nb==0.and.nint(occ(ib*2-1, ik)).eq.1)) &
        ksi(ik)=ksi(ik)*ksi_b(ib)
      !
    enddo !ib
    !
  enddo
  !
  deallocate(ksi_b, grp)
  !
  write(*,*) "!!!!!!!! FINAL RESULTS !!!!!!!!"
  do ik=1, nkpt
    write(*,'(A, 1I2, A, 1F14.9)') "kpt: ", ik, " ksi: ", ksi(ik)
  enddo
  !
END PROGRAM

SUBROUTINE groupstates(grp, ngrp, ik)
  !
  use constants,        ONLY : dp, eps6
  use wvfunc,           ONLY : nbnd, eig
  !
  implicit none
  !
  integer,dimension(nbnd/2+1) ::  grp
  integer ngrp, ik
  !
  integer ib, ii
  real(dp) etmp
  !
  etmp=eig(1,ik)-1.0
  !
  ii=0
  !
  do ib=1, nbnd, 2
    !
    if (abs(etmp-eig(ib,ik))>eps6) then
      ii=ii+1
      grp(ii)=ib
      etmp=eig(ib, ik)
    endif
    !
  enddo
  !
  ngrp=ii+1
  grp(ngrp)=nbnd+1
  !
END SUBROUTINE

