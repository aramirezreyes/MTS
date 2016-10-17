module DFLib

  USE SysLinCSR
  USE Precond
  USE BLASS

  IMPLICIT NONE

  CONTAINS

!cccccccccccccccccccccccccccccccccccc
  subroutine FillMatrixDF(N,AA,JA,IA)
!cccccccccccccccccccccccccccccccccccc
    IMPLICIT NONE
    integer, intent(in) :: N
    real (kind=8), intent(out), dimension(:) :: AA
    integer, intent(out), dimension(:) :: JA,IA
    integer :: valind,ind,line
! TO BE COMPLETED
    valind = 1
    IA(1) = 1
    AA(valind:valind+3-1) = (/4.,-1.,-1./)
    JA(valind:valind+3-1) = (/1,2,N+1/)
    valind = valind + 3
    DO ind=2,N-1
       IA(ind) = valind
       AA(valind:valind+4-1) = (/-1.,4.,-1.,-1./)
       JA(valind:valind+4-1) = (/ind - 1,ind,ind+1,ind+N/)
       valind = valind + 4
    END DO
    ind = N
    IA(ind) = N-1
    AA(valind:valind+3-1) = (/-1.,-4.,-1./)
    JA(valind:valind+3-1) = (/N-1,N,2*N/)
    valind = valind + 3



    DO sub=1,N-2
       
          IA(N+1) = 1
          AA(valind:valind+4-1) = (/-1,-4.,-1.,-1./)
          JA(valind:valind+4-1) = (/1,N+1,N+2,2*N+1/)
          valind = valind + 4
          
          DO ind=2,N-1
             IA(sub*N+ind) = valind
             AA(valind:valind+5-1) = (/-1.,-1.,4.,-1.,-1./)
             JA(valind:valind+5-1) = (/sub*N+ind-1,ind - 1,ind,ind+1,ind+N/)
             valind = valind + 5
          END DO
          ind = N
          IA(ind) = N-1
          AA(valind:valind+3-1) = (/-1.,-4.,-1./)
          JA(valind:valind+3-1) = (/N-1,N,2*N/)
          valind = valind + 3
       
    END DO

   
  end subroutine FillMatrixDF

!cccccccccccccccccccccccccccccccc
subroutine FillRHSDF(N,X,Y,RHS)
!cccccccccccccccccccccccccccccccc

  IMPLICIT NONE
  integer,intent(in) :: N
  real (kind=8), intent(in), dimension(:)  :: X,Y
  real (kind=8), intent(out), dimension(:) :: RHS
  
! TO BE COMPLETED

  end subroutine FillRHSDF
 
!ccccccccccccccccccccccccccccc
function source(x,y) result(z)
!ccccccccccccccccccccccccccccc

      implicit none
         real(kind=8), intent(in) :: x,y
         real(kind=8) :: r,z
      r=sqrt(x**2.+y**2.)
      z=0.
      if ((r>1.5).and.(r<2)) then
         z=z+10000*((r-1.5)**2.)*(r-2)**2.
      end if
      if (r<1.) then
         z=z-250*(r**2.)*((r-1)**2.)
      end if
      return

end function source


!ccccccccccccccccccccccccccccccccccccccccccccc
  subroutine ManageBC(N,AA,JA,IA,RHS,BDIND,CL)
!ccccccccccccccccccccccccccccccccccccccccccccc

  IMPLICIT NONE

    integer, intent(in) :: N
    real (kind=8), intent(out), dimension(:) :: AA, RHS
    integer, intent(in), dimension(:) :: JA,IA
    integer, intent(in), dimension(:) :: BDIND
    real (kind=8),intent(in) :: CL

! TO BE COMPLETED 

  end subroutine ManageBC

  
end module DFLib
