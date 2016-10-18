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
    integer :: valind,ind,line,sub
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


    DO sub=2,N-1
          write(*,*) "Holas"       
          IA((sub-1)*N+1) = 1
          write(*,*) "Hola2"
          AA(valind:valind+4-1) = (/-1.,-4.,-1.,-1./)
          write(*,*) "Hola3"
          JA(valind:valind+4-1) = (/1,N+1,N+2,2*N+1/)
          valind = valind + 4
          
          DO ind=2,N-1
             write(*,*) (sub-1)*N+ind,"ind",ind,"sub",sub
             IA((sub-1)*N+ind) = ind+(sub-1)*N
             write(*,*) "Hola"
              AA(valind:valind+5-1) = (/-1.,-1.,4.,-1.,-1./)
             JA(valind:valind+5-1) = (/N*(sub-1)+ind,N*(sub)+ind+1,N*(sub)+ind+2,N*(sub)+ind+3,N*(sub+1)+ind/) 
             valind = valind + 5
          END DO
          ind = (sub-1)*N+N
          IA(N+(sub*N)) = N-1
          AA(valind:valind+3-1) = (/-1.,-4.,-1./)
          JA(valind:valind+3-1) = (/N-1,N,2*N/)
          valind = valind + 3
       
    END DO

    
    IA((N-1)*N+1) = (N-1)*N-1
    AA(valind:valind+3-1) = (/-1.,-1.,4./)
    JA(valind:valind+3-1) = (/(N*(N-2))+1,N*(N-1)+1,N*(N-1)+2/)
    valind = valind + 3
    DO ind=2,N-1
       IA((N-1)*(N)+ind) = valind
       AA(valind:valind+4-1) = (/-1.,-1.,4.,-1./)
       JA(valind:valind+4-1) = (/(N*(N-2))+ind,(N*(N-1))+ind-1,(N*(N-2))+ind,(N*(N-2))+ind+1/)
       valind = valind + 4
    END DO
    ind = N*2
    IA((N-1)*N+N) =N*(N-1)+1
    AA(valind:valind+3-1) = (/-1.,-1.,-4./)
    JA(valind:valind+3-1) = (/N*(N-1),N*(N-1)+1,N*(N-1)+2/)
    valind = valind + 3
    write (*,*) "Adios"

    write(*,*) AA
    write(*,*) JA
    write(*,*) IA
   
  end subroutine FillMatrixDF

!cccccccccccccccccccccccccccccccc
subroutine FillRHSDF(N,X,Y,RHS)
!cccccccccccccccccccccccccccccccc

  IMPLICIT NONE
  integer,intent(in) :: N
  real (kind=8), intent(in), dimension(:)  :: X,Y
  real (kind=8), intent(out), dimension(:) :: RHS
  integer :: i
  write(*,*) "I will fill RHS"
  DO i=1,N*N
     write(*,*) source(X(i),Y(i))
write(*,*) "Written element ", i, "of RHS"
  END DO

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
    real (kind=8), dimension(:),ALLOCATABLE :: V1,V2
    real :: g1
    integer :: IADD, LIG,ind
! TO BE COMPLETED 

    ALLOCATE (V1(N*N),V2(N*N))

    V1 = 0.0; V2=0.0
    DO ind=1,2*N+2*(N-2)
       LIG = BDIND(ind)
       V1(LIG) = CL
       CALL AMUX(N*N,V1,V2,AA,JA,IA)
       RHS = RHS - V2
       V1 = 0.0
       WHERE(JA.EQ.LIG) AA=0.0
       AA(IA(LIG):IA(LIG+1)-1) = 0.0
       g1 = GETELM(LIG,LIG,AA,JA,IA,IADD,.true.)
       AA(IADD) = 1.0
       RHS(LIG) = CL
    END DO


  end subroutine ManageBC

  
end module DFLib
