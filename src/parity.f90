MODULE parity
  !
  implicit none
  !
  integer, allocatable :: igrev(:)
  !
 CONTAINS
  !
  SUBROUTINE invglist(ik)
    !
    use constants,        ONLY : dp, eps5
    use wvfunc,           ONLY : npw, gvec, grid
    !
    implicit none
    !
    integer       ik
    !
    ! Temporary ...
    !
    real(dp) :: xx(1:3)
    integer ii, jj           ! local variables
    !
    if (allocated(igrev)) deallocate(igrev)
    allocate(igrev(npw(ik)))
    !
    ! Determine the index of -(k+G)
    !
    do ii=1, npw(ik)
      igrev(ii)=-1
      do jj=1, npw(ik)
        xx(:)=gvec(:, ii, ik)+gvec(:, jj, ik)
        if (sum(xx(:)*xx(:))< eps5) then
          igrev(ii)=jj
        endif
      enddo
      if (igrev(ii)<0) then
        write(*,*) "!!!! Cannot find -k for grid", ii, ": ", grid(:, ii, ik)
        write(*,*) "  Please make sure last calculation is for Z2!"
        stop
      endif
    enddo
    !
  END SUBROUTINE
  !
  SUBROUTINE calc_pmat(pmat, ib, jb, ik)
    !
    use constants,        ONLY : dp
    use wvfunc,           ONLY : wfc
    !
    implicit none
    !
    complex(dp) pmat
    !
    real(dp) int0
    complex(dp) int1
    !
    integer ib, jb, ik
    !
    ! int0 is <psi_i | psi_j>
    !
    int0=sum(conjg(wfc(:, ib, ik, 1))*wfc(:, ib, ik, 1))+ &
         sum(conjg(wfc(:, ib, ik, 2))*wfc(:, ib, ik, 2))
    !
    if (ib.ne.jb) then
      int0=int0*(sum(conjg(wfc(:, jb, ik, 1))*wfc(:, jb, ik, 1))+ &
                 sum(conjg(wfc(:, jb, ik, 2))*wfc(:, jb, ik, 2)))
      int0=sqrt(int0)
    endif
    !
    ! int1 is <psi_i|P|psi_i>
    ! Strictly speaking, S|psi> should be calculated
    !
    int1=sum(conjg(wfc(igrev(:), ib, ik, 1))*wfc(:, jb, ik, 1))+ &
         sum(conjg(wfc(igrev(:), ib, ik, 2))*wfc(:, jb, ik, 2))
    !
    pmat=int1/int0
    !
  END SUBROUTINE
  !
END MODULE
