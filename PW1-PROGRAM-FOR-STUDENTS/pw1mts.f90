program pw1mts

!cccccccccccccccccccccccccccc
! DECLARATION OF USED MODULUS
!cccccccccccccccccccccccccccc
  
  USE SysLinCSR
  USE Precond
  USE BLASS
  USE DFLib

!ccccccccccccccccccccccccc
! DECLARATION OF VARIABLES
!ccccccccccccccccccccccccc

  IMPLICIT NONE

! Domain dimension 
  REAL (KIND=8) :: XYMIN, XYMAX
! Discretization step
   REAL (KIND=8) :: H
! Value of the BC to be applied on the boundary
   REAL (KIND=8) :: CL
! Number of nodes in each direction of the mesh
  INTEGER :: N
! Number of nodes in the mesh
  INTEGER :: NNODES
! Number of non nul coefficients in the matrix
  INTEGER :: NNE
! COORDINATES OF THE MESH NODES
  REAL (KIND=8), DIMENSION(:), ALLOCATABLE :: X,Y
! VECTORS FOR THE CSR STORAGE OF THE MATRIX
  REAL (KIND=8), DIMENSION(:), ALLOCATABLE:: AA
  INTEGER, DIMENSION(:), ALLOCATABLE :: IA, JA
! SOME VECTORS
  REAL (KIND=8), DIMENSION(:), ALLOCATABLE:: SOL, RHS
  INTEGER, DIMENSION(:), ALLOCATABLE :: BDIND
  INTEGER, DIMENSION(:), ALLOCATABLE :: iwork
! MATRIX OF THE LINEAR SYSTEM
  TYPE(MatriceCSR) :: Mat2D
! SOME LOCAL VARIABLES
  INTEGER :: ind

!cccccccccccccccccccccc
! VARIABLES AFFECTATION
!cccccccccccccccccccccc
  XYMIN = -3.0
  XYMAX = 3.0
  N = 50
  H=(XYMAX-XYMIN)/(N-1)
  CL=20.

!cccccccccccccccccc
! NODES COORDINATES
!cccccccccccccccccc
  ALLOCATE(X(N),Y(N))
  X=(/ XYMIN, (XYMIN+(ind-1)*H, ind=2,N-1), XYMAX /)
  Y=X

!cccccccccccccccccccccccccccc
! NUMBER OF NODES IN THE MESH
!cccccccccccccccccccccccccccc
  NNODES =  N*N 

!ccccccccccccccccccccccccccccccccccccccccccccccc
! NUMBER OF NON ZEROS COEFFICIENTS IN THE MATRIX
!ccccccccccccccccccccccccccccccccccccccccccccccc
  NNE    =  N*(5*N-4)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
! ALLOCATION AND FILLING OF THE ARRAYS FOR CSR STORAGE
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ALLOCATE(AA(NNE),JA(NNE),IA(NNODES+1))
  AA=0.; JA=0 ; IA=0
  CALL FillMatrixDF(N,AA,JA,IA)

  allocate(iwork(max(NNODES+1,2*NNE)))
  call csort(NNODES,AA,JA,IA,iwork,.true.) 
  deallocate(iwork)

!ccccccccccccccccccccccccccccccccccccccc
! RIGHT-HAND-SIDE ALLOCATION AND FILLING 
!ccccccccccccccccccccccccccccccccccccccc
  allocate(RHS(NNODES))
  CALL FillRHSDF(N,X,Y,RHS)

!cccccccccccccccccccccccccccccc
! BOUNDARY CONDITIONS TREATMENT
!cccccccccccccccccccccccccccccc
! Numbers of the boundary nodes
  ALLOCATE(BDIND(2*N+2*(N-2)))
  BDIND=(/ (ind,ind=1,N) , (ind*N+1,ind=1,N-2), ((ind+1)*N,ind=1,N-2), (N*(N-1)+ind, ind=1,N)  /)
  CALL ManageBC(N,AA,JA,IA,RHS,BDIND,CL)

!cccccccccccccccccccccccccccccccccccccccccccccc
! CREATION AND AFFECTATION OF THE DATA IN Mat2D 
!cccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccc
  CALL CREATE(Mat2D,NNODES,2*NNE)
!cccccccccccccccccccccccccccccccc
  ALLOCATE (Mat2D%Guess(Mat2d%NbLine)) ! Needed for iterative method
! Copy data in the matrix structure Mat2D
  Mat2D%Resolution = 1
  Mat2D%Coeff(:) = AA(:)
  Mat2D%IA(:) = IA(:)
  Mat2D%JA(:) = JA(:)
  Mat2D%Guess(:) = 0.
! Set-up the preconditioner ILUT(15, 1E-4)  ! new definition of lfil
!
  Mat2D%Precond=3        ! Precodintioner type: ILUT
  Mat2D%LFill=3          ! Filling level
  Mat2D%DropTol=1.0D-4   ! Tolerance
  ! Computation of the preconditioner
  CALL CalPrecond(Mat2D)

!cccccccccccccccccccccccccccc
! COMPUTATION OF THE SOLUTION 
!cccccccccccccccccccccccccccc
  allocate(SOL(NNODES))
  CALL RESOL(Mat2D,RHS(1:NNODES),SOL(1:NNODES))

!ccccccccccccccccccccccccccccccccccccccccccccccccc
! WRITE IN OUTPUT FILES FOR MATLAB POST-PROCESSING
!ccccccccccccccccccccccccccccccccccccccccccccccccc

  open(unit=100,file='solution.dat',status='unknown')
  do ind=1,NNODES
        write(100,*) SOL(ind)
  end do
  close(100)
  open(unit=200,file='valX.dat',status='unknown')
  do ind=1,N
        write(200,*) X(ind)
  end do
  close(200)
  open(unit=300,file='valY.dat',status='unknown')
  do ind=1,N
        write(300,*) Y(ind)
  end do
  close(300)

!ccccccccccccc
! DEALLOCATION
!ccccccccccccc
  CALL DESTROY(Mat2D)
  DEALLOCATE(X,Y,AA,JA,IA,BDIND,RHS,SOL)

end program pw1mts
