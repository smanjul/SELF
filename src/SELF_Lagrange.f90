! SELF_Lagrange.f90
!
! Copyright 2020 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : self@higherordermethods.org
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !

module SELF_Lagrange

  use iso_fortran_env
  use iso_c_binding
  use hipfort
  use hipfort_check
  use hipfort_hipmalloc
  use hipfort_hipblas

  use SELF_Constants
  use SELF_SupportRoutines
  use SELF_Quadrature
  !use SELF_HDF5
  !use HDF5

  use hipfort_hipblas

  use iso_c_binding

  implicit none

! #ifdef DOUBLE_PRECISION

! #define hipblasXgemm(handle,opA,opB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) hipblasDgemm(handle,opA,opB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc)
! #define hipblasXgemvStridedBatched(handle,opA,m,n,alpha,A,lda,stridA,x,incx,stridex,beta,y,incy,stridey,batchCount) hipblasDgemvStridedBatched(handle,opA,m,n,alpha,A,lda,stridA,x,incx,stridex,beta,y,incy,stridey,batchCount)

! #else

! #define hipblasXgemm(handle,opA,opB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc) hipblasSgemm(handle,opA,opB,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc)
! #define hipblasXgemvStridedBatched(handle,opA,m,n,alpha,A,lda,stridA,x,incx,stridex,beta,y,incy,stridey,batchCount) hipblasSgemvStridedBatched(handle,opA,m,n,alpha,A,lda,stridA,x,incx,stridex,beta,y,incy,stridey,batchCount)

! #endif

  type,public :: Lagrange
    !! A data structure for working with Lagrange Interpolating Polynomials in one, two, and three dimensions.
    !! The Lagrange data-structure stores the information necessary to interpolate between two
    !! sets of grid-points and to estimate the derivative of data at native grid points. Routines for
    !! multidimensional interpolation are based on the tensor product of 1-D interpolants. It is
    !! assumed that the polynomial degree (and the interpolation nodes) are the same in each direction.
    !! This assumption permits the storage of only one array of interpolation nodes and barycentric
    !! weights and is what allows this data structure to be flexible.

    integer :: N
      !! The number of control points.

    integer :: controlNodeType

    integer :: M
      !! The number of target points.

    integer :: targetNodeType

    real(prec),pointer,dimension(:) :: controlPoints
      !! The set of nodes in one dimension where data is known.
      !! To create higher dimension interpolation and differentiation operators, structured grids in two and three
      !! dimensions are created by tensor products of the controlPoints. This design decision implies that all
      !! spectral element methods supported by the Lagrange class have the same polynomial degree in each
      !! computational/spatial dimension. In practice, the controlPoints are the Legendre-Gauss, Legendre-Gauss-Lobatto,
      !! Legendre-Gauss-Radau, Chebyshev-Gauss, Chebyshev-Gauss-Lobatto, or Chebyshev-Gauss-Radau quadrature points over
      !! the domain [-1,1] (computational space). The Init routine for this class restricts controlPoints to one of
      !! these quadrature types or uniform points on [-1,1].

    real(prec),pointer,dimension(:) :: targetPoints
      !! The set of nodes in one dimension where data is to be interpolated to. To create higher dimension interpolation
      !! and differentiation operators, structured grids in two and three dimensions are created by tensor products of
      !! the targetPoints. In practice, the targetPoints are set to a uniformly distributed set of points between [-1,1]
      !! (computational space) to allow for interpolation from unevenly spaced quadrature points to a plotting grid.

    real(prec),pointer,dimension(:) :: bWeights
      !! The barycentric weights that are calculated from the controlPoints and used for interpolation.

    real(prec),pointer,dimension(:) :: qWeights
      !! The quadrature weights for discrete integration. The quadradture weights depend on the type of controlPoints
      !! provided; one of Legendre-Gauss, Legendre-Gauss-Lobatto, Legendre-Gauss-Radau, Chebyshev-Gauss,
      !! Chebyshev-Gauss-Lobatto, Chebyshev-Gauss Radau, or Uniform. If Uniform, the quadrature weights are constant
      !! $$dx = \frac{2.0}{N+1}$$.

    real(prec),pointer,dimension(:,:) :: iMatrix
      !! The interpolation matrix (transpose) for mapping data from the control grid to the target grid.

    real(prec),pointer,dimension(:,:) :: dMatrix
      !! The derivative matrix for mapping function nodal values to a nodal values of the derivative estimate. The
      !! dMatrix is based on a strong form of the derivative.

    real(prec),pointer,dimension(:,:) :: dgMatrix
      !! The derivative matrix for mapping function nodal values to a nodal values of the derivative estimate. The dgMatrix is based
      !! on a weak form of the derivative. It must be used with bMatrix to account for boundary contributions in the weak form.

    real(prec),pointer,dimension(:,:) :: bMatrix
      !! The boundary interpolation matrix that is used to map a grid of nodal values at the control points to the element boundaries.

  contains

    procedure,public :: Init => Init_Lagrange
    procedure,public :: Free => Free_Lagrange

    procedure,public :: UpdateDevice => UpdateDevice_Lagrange

    generic,public :: ScalarGridInterp_1D => ScalarGridInterp_1D_cpu,ScalarGridInterp_1D_gpu
    procedure,private :: ScalarGridInterp_1D_cpu,ScalarGridInterp_1D_gpu

    generic,public :: ScalarGridInterp_2D => ScalarGridInterp_2D_cpu,ScalarGridInterp_2D_gpu
    procedure,private :: ScalarGridInterp_2D_cpu,ScalarGridInterp_2D_gpu

    ! GENERIC,PUBLIC :: VectorGridInterp_2D => VectorGridInterp_2D_cpu,VectorGridInterp_2D_gpu
    ! PROCEDURE,PRIVATE :: VectorGridInterp_2D_cpu,VectorGridInterp_2D_gpu

    ! GENERIC,PUBLIC :: ScalarGridInterp_3D => ScalarGridInterp_3D_cpu,ScalarGridInterp_3D_gpu
    ! PROCEDURE,PRIVATE :: ScalarGridInterp_3D_cpu,ScalarGridInterp_3D_gpu

    ! GENERIC,PUBLIC :: VectorGridInterp_3D => VectorGridInterp_3D_cpu,VectorGridInterp_3D_gpu
    ! PROCEDURE,PRIVATE :: VectorGridInterp_3D_cpu,VectorGridInterp_3D_gpu

    generic,public :: ScalarBoundaryInterp_1D => ScalarBoundaryInterp_1D_cpu,ScalarBoundaryInterp_1D_gpu
    procedure,private :: ScalarBoundaryInterp_1D_cpu,ScalarBoundaryInterp_1D_gpu

    ! GENERIC,PUBLIC :: ScalarBoundaryInterp_2D => ScalarBoundaryInterp_2D_cpu,ScalarBoundaryInterp_2D_gpu
    ! PROCEDURE,PRIVATE :: ScalarBoundaryInterp_2D_cpu,ScalarBoundaryInterp_2D_gpu

    ! GENERIC,PUBLIC :: VectorBoundaryInterp_2D => VectorBoundaryInterp_2D_cpu,VectorBoundaryInterp_2D_gpu
    ! PROCEDURE,PRIVATE :: VectorBoundaryInterp_2D_cpu,VectorBoundaryInterp_2D_gpu

    ! GENERIC,PUBLIC :: TensorBoundaryInterp_2D => TensorBoundaryInterp_2D_cpu,TensorBoundaryInterp_2D_gpu
    ! PROCEDURE,PRIVATE :: TensorBoundaryInterp_2D_cpu,TensorBoundaryInterp_2D_gpu

    ! GENERIC,PUBLIC :: ScalarBoundaryInterp_3D => ScalarBoundaryInterp_3D_cpu,ScalarBoundaryInterp_3D_gpu
    ! PROCEDURE,PRIVATE :: ScalarBoundaryInterp_3D_cpu,ScalarBoundaryInterp_3D_gpu

    ! GENERIC,PUBLIC :: VectorBoundaryInterp_3D => VectorBoundaryInterp_3D_cpu,VectorBoundaryInterp_3D_gpu
    ! PROCEDURE,PRIVATE :: VectorBoundaryInterp_3D_cpu,VectorBoundaryInterp_3D_gpu

    ! GENERIC,PUBLIC :: TensorBoundaryInterp_3D => TensorBoundaryInterp_3D_cpu,TensorBoundaryInterp_3D_gpu
    ! PROCEDURE,PRIVATE :: TensorBoundaryInterp_3D_cpu,TensorBoundaryInterp_3D_gpu

    generic,public :: Derivative_1D => Derivative_1D_cpu,Derivative_1D_gpu
    procedure,private :: Derivative_1D_cpu,Derivative_1D_gpu

    ! generic,public :: DGDerivative_1D => DGDerivative_1D_cpu,DGDerivative_1D_gpu
    ! procedure,private :: DGDerivative_1D_cpu,DGDerivative_1D_gpu

    ! GENERIC,PUBLIC :: ScalarGradient_2D => ScalarGradient_2D_cpu,ScalarGradient_2D_gpu
    ! PROCEDURE,PRIVATE :: ScalarGradient_2D_cpu,ScalarGradient_2D_gpu

    ! GENERIC,PUBLIC :: VectorGradient_2D => VectorGradient_2D_cpu,VectorGradient_2D_gpu
    ! PROCEDURE,PRIVATE :: VectorGradient_2D_cpu,VectorGradient_2D_gpu

    ! GENERIC,PUBLIC :: VectorDivergence_2D => VectorDivergence_2D_cpu,VectorDivergence_2D_gpu
    ! PROCEDURE,PRIVATE :: VectorDivergence_2D_cpu,VectorDivergence_2D_gpu

    ! GENERIC,PUBLIC :: VectorDGDivergence_2D => VectorDGDivergence_2D_cpu,VectorDGDivergence_2D_gpu
    ! PROCEDURE,PRIVATE :: VectorDGDivergence_2D_cpu,VectorDGDivergence_2D_gpu

    ! GENERIC,PUBLIC :: ScalarGradient_3D => ScalarGradient_3D_cpu,ScalarGradient_3D_gpu
    ! PROCEDURE,PRIVATE :: ScalarGradient_3D_cpu,ScalarGradient_3D_gpu

    ! GENERIC,PUBLIC :: VectorGradient_3D => VectorGradient_3D_cpu,VectorGradient_3D_gpu
    ! PROCEDURE,PRIVATE :: VectorGradient_3D_cpu,VectorGradient_3D_gpu

    ! GENERIC,PUBLIC :: VectorDivergence_3D => VectorDivergence_3D_cpu,VectorDivergence_3D_gpu
    ! PROCEDURE,PRIVATE :: VectorDivergence_3D_cpu,VectorDivergence_3D_gpu

    ! GENERIC,PUBLIC :: VectorDGDivergence_3D => VectorDGDivergence_3D_cpu,VectorDGDivergence_3D_gpu
    ! PROCEDURE,PRIVATE :: VectorDGDivergence_3D_cpu,VectorDGDivergence_3D_gpu

    !procedure,public :: WriteHDF5 => WriteHDF5_Lagrange
    procedure,private :: CalculateBarycentricWeights
    procedure,private :: CalculateInterpolationMatrix
    procedure,private :: CalculateDerivativeMatrix
    procedure,private :: CalculateLagrangePolynomials

  end type Lagrange

  ! /////////////// !
  ! Boundary Interpolation Routines

  ! interface
  !   subroutine ScalarBoundaryInterp_1D_gpu_wrapper(bMatrix_dev,f_dev,fTarget_dev,N,nVar,nEl) &
  !     bind(c,name="ScalarBoundaryInterp_1D_gpu_wrapper")
  !     use iso_c_binding
  !     implicit none
  !     type(c_ptr) :: bMatrix_dev,f_dev,fTarget_dev
  !     integer(c_int),value :: N,nVar,nEl
  !   end subroutine ScalarBoundaryInterp_1D_gpu_wrapper
  ! end interface

  ! INTERFACE
  !   SUBROUTINE ScalarBoundaryInterp_2D_gpu_wrapper(bMatrix_dev,f_dev,fTarget_dev,N,nVar,nEl) &
  !     bind(c,name="ScalarBoundaryInterp_2D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: bMatrix_dev,f_dev,fTarget_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE ScalarBoundaryInterp_2D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE VectorBoundaryInterp_2D_gpu_wrapper(bMatrix_dev,f_dev,fTarget_dev,N,nVar,nEl) &
  !     bind(c,name="VectorBoundaryInterp_2D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: bMatrix_dev,f_dev,fTarget_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE VectorBoundaryInterp_2D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE TensorBoundaryInterp_2D_gpu_wrapper(bMatrix_dev,f_dev,fTarget_dev,N,nVar,nEl) &
  !     bind(c,name="TensorBoundaryInterp_2D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: bMatrix_dev,f_dev,fTarget_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE TensorBoundaryInterp_2D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE ScalarBoundaryInterp_3D_gpu_wrapper(bMatrix_dev,f_dev,fTarget_dev,N,nVar,nEl) &
  !     bind(c,name="ScalarBoundaryInterp_3D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: bMatrix_dev,f_dev,fTarget_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE ScalarBoundaryInterp_3D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE VectorBoundaryInterp_3D_gpu_wrapper(bMatrix_dev,f_dev,fTarget_dev,N,nVar,nEl) &
  !     bind(c,name="VectorBoundaryInterp_3D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: bMatrix_dev,f_dev,fTarget_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE VectorBoundaryInterp_3D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE TensorBoundaryInterp_3D_gpu_wrapper(bMatrix_dev,f_dev,fTarget_dev,N,nVar,nEl) &
  !     bind(c,name="TensorBoundaryInterp_3D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: bMatrix_dev,f_dev,fTarget_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE TensorBoundaryInterp_3D_gpu_wrapper
  ! END INTERFACE

  ! /////////////// !

  ! interface
  !   subroutine Derivative_1D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
  !     bind(c,name="Derivative_1D_gpu_wrapper")
  !     use iso_c_binding
  !     implicit none
  !     type(c_ptr) :: dMatrixT_dev,f_dev,df_dev
  !     integer(c_int),value :: N,nVar,nEl
  !   end subroutine Derivative_1D_gpu_wrapper
  ! end interface

  ! interface
  !   subroutine DGDerivative_1D_gpu_wrapper(dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev,N,nVar,nEl) &
  !     bind(c,name="DGDerivative_1D_gpu_wrapper")
  !     use iso_c_binding
  !     implicit none
  !     type(c_ptr) :: dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev
  !     integer(c_int),value :: N,nVar,nEl
  !   end subroutine DGDerivative_1D_gpu_wrapper
  ! end interface

  ! INTERFACE
  !   SUBROUTINE ScalarGradient_2D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
  !     bind(c,name="ScalarGradient_2D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE ScalarGradient_2D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE VectorGradient_2D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
  !     bind(c,name="VectorGradient_2D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE VectorGradient_2D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE VectorDGGradient_2D_gpu_wrapper(dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev,N,nVar,nEl) &
  !     bind(c,name="VectorDGGradient_2D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE VectorDGGradient_2D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE VectorDivergence_2D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
  !     bind(c,name="VectorDivergence_2D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE VectorDivergence_2D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !  SUBROUTINE VectorDGDivergence_2D_gpu_wrapper(dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev,N,nVar,nEl) &
  !     bind(c,name="VectorDGDivergence_2D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE VectorDGDivergence_2D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE ScalarGradient_3D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
  !     bind(c,name="ScalarGradient_3D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE ScalarGradient_3D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE VectorGradient_3D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
  !     bind(c,name="VectorGradient_3D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE VectorGradient_3D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !   SUBROUTINE VectorDivergence_3D_gpu_wrapper(dMatrixT_dev,f_dev,df_dev,N,nVar,nEl) &
  !     bind(c,name="VectorDivergence_3D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: dMatrixT_dev,f_dev,df_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE VectorDivergence_3D_gpu_wrapper
  ! END INTERFACE

  ! INTERFACE
  !  SUBROUTINE VectorDGDivergence_3D_gpu_wrapper(dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev,N,nVar,nEl) &
  !     bind(c,name="VectorDGDivergence_3D_gpu_wrapper")
  !     USE iso_c_binding
  !     IMPLICIT NONE
  !     TYPE(c_ptr) :: dgMatrixT_dev,bMatrix_dev,qWeights_dev,f_dev,bf_dev,df_dev
  !     INTEGER(C_INT),VALUE :: N,nVar,nEl
  !   END SUBROUTINE VectorDGDivergence_3D_gpu_wrapper
  ! END INTERFACE

contains

  subroutine Init_Lagrange(this,N,controlNodeType,M,targetNodeType)
    !! Initialize an instance of the Lagrange class
    !! On output, all of the attributes for the Lagrange class are allocated and values are initialized according to the number of
    !! control points, number of target points, and the types for the control and target nodes.
    !! If a GPU is available, device pointers for the Lagrange attributes are allocated and initialized.
    implicit none
    class(Lagrange),intent(out) :: this
    !! Lagrange class instance
    integer,intent(in)          :: N
    !! The number of control points for interpolant
    integer,intent(in)          :: M
    !! The number of target points for the interpolant
    integer,intent(in)          :: controlNodeType
    !! The integer code specifying the type of control points. Parameters are defined in SELF_Constants.f90. One of GAUSS(=1),
    !! GAUSS_LOBATTO(=2), or UNIFORM(=3)
    integer,intent(in)          :: targetNodeType
    !! The integer code specifying the type of target points. Parameters are defined in SELF_Constants.f90. One of GAUSS(=1),
    !! GAUSS_LOBATTO(=2), or UNIFORM(=3)
    ! -------!
    ! Local
    real(prec) :: q(0:M)

    this % N = N
    this % M = M
    this % controlNodeType = controlNodeType
    this % targetNodeType = targetNodeType

    call hipcheck(hipMallocManaged(this % controlPoints,N + 1,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % targetPoints,M + 1,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % bWeights,N + 1,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % qWeights,N + 1,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % iMatrix,N + 1,M + 1,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % dMatrix,N + 1,N + 1,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % dgMatrix,N + 1,N + 1,hipMemAttachGlobal))
    call hipcheck(hipMallocManaged(this % bMatrix,N + 1,2,hipMemAttachGlobal))

    if (controlNodeType == GAUSS .or. controlNodeType == GAUSS_LOBATTO) then

      call LegendreQuadrature(N, &
                              this % controlPoints, &
                              this % qWeights, &
                              controlNodeType)

    elseif (controlNodeType == CHEBYSHEV_GAUSS .or. controlNodeType == CHEBYSHEV_GAUSS_LOBATTO) then

      call ChebyshevQuadrature(N, &
                               this % controlPoints, &
                               this % qWeights, &
                               controlNodeType)

    elseif (controlNodeType == UNIFORM) then

      this % controlPoints = UniformPoints(-1.0_prec,1.0_prec,0,N)
      this % qWeights = 2.0_prec/real(N,prec)

    end if

    ! Target Points
    if (targetNodeType == GAUSS .or. targetNodeType == GAUSS_LOBATTO) then

      call LegendreQuadrature(M, &
                              this % targetPoints, &
                              q, &
                              targetNodeType)

    elseif (targetNodeType == UNIFORM) then

      this % targetPoints = UniformPoints(-1.0_prec,1.0_prec,0,M)

    end if

    call this % CalculateBarycentricWeights()
    call this % CalculateInterpolationMatrix()
    call this % CalculateDerivativeMatrix()
    this % bMatrix(1:N + 1,1) = this % CalculateLagrangePolynomials(-1.0_prec)
    this % bMatrix(1:N + 1,2) = this % CalculateLagrangePolynomials(1.0_prec)

    call this % UpdateDevice()

  end subroutine Init_Lagrange

  subroutine Free_Lagrange(this)
    !! Frees all memory (host and device) associated with an instance of the Lagrange class
    implicit none
    class(Lagrange),intent(inout) :: this
    !! Lagrange class instance

    call hipcheck(hipFree(this % controlPoints))
    call hipcheck(hipFree(this % targetPoints))
    call hipcheck(hipFree(this % bWeights))
    call hipcheck(hipFree(this % qWeights))
    call hipcheck(hipFree(this % iMatrix))
    call hipcheck(hipFree(this % dMatrix))
    call hipcheck(hipFree(this % dgMatrix))
    call hipcheck(hipFree(this % bMatrix))

  end subroutine Free_Lagrange

  subroutine UpdateDevice_Lagrange(this)
    !! Copy the Lagrange attributes from the host (CPU) to the device (GPU)
    implicit none
    class(Lagrange),intent(inout) :: this
    !! Lagrange class instance

    call hipcheck(hipMemPrefetchAsync(c_loc(this % controlPoints),sizeof(this % controlPoints),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % targetPoints),sizeof(this % targetPoints),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % bWeights),sizeof(this % bWeights),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % qWeights),sizeof(this % qWeights),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % iMatrix),sizeof(this % iMatrix),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % dMatrix),sizeof(this % dMatrix),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % dgMatrix),sizeof(this % dgMatrix),0,c_null_ptr))
    call hipcheck(hipMemPrefetchAsync(c_loc(this % bMatrix),sizeof(this % bMatrix),0,c_null_ptr))

  end subroutine UpdateDevice_Lagrange

  subroutine ScalarGridInterp_1D_cpu(this,f,fTarget,nvars,nelems)
    !! Host (CPU) implementation of the ScalarGridInterp_1D interface.
    !! In most cases, you should use the `ScalarGridInterp_1D` generic interface,
    !! rather than calling this routine directly.
    !! Interpolate a scalar-1D (real) array from the control grid to the target grid.
    !! The control and target grids are the ones associated with an initialized
    !! Lagrange instance.
    !!
    !! Interpolation is applied using a series of matrix-vector multiplications, using
    !! the Lagrange class's interpolation matrix
    !!
    !! $$ \tilde{f}_{m,iel,ivar} = \sum_{i=0}^N f_{i,iel,ivar} I_{i,m} $$
    !!
    implicit none
    class(Lagrange),intent(in) :: this
    !! Lagrange class instance
    integer,intent(in)     :: nvars
    !! The number of variables/functions that are interpolated
    integer,intent(in)     :: nelems
    !! The number of spectral elements in the SEM grid
    real(prec),intent(in)  :: f(1:this % N + 1,1:nelems,1:nvars)
    !! (Input) Array of function values, defined on the control grid
    real(prec),intent(out) :: fTarget(1:this % M + 1,1:nelems,1:nvars)
    !! (Output) Array of function values, defined on the target grid
    ! Local
    integer :: iel,ivar,i,ii
    real(prec) :: floc

    do ivar = 1,nvars
      do iel = 1,nelems
        do i = 1,this % M + 1
          floc = 0.0_prec
          do ii = 1,this % N + 1
            floc = floc + this % iMatrix(ii,i)*f(ii,iel,ivar)
          end do
          fTarget(i,iel,ivar) = floc
        end do
      end do
    end do

  end subroutine ScalarGridInterp_1D_cpu

  subroutine ScalarGridInterp_1D_gpu(this,f,fTarget,nvars,nelems,hipblas_handle)
    !! Device (GPU) implementation of the ScalarGridInterp_1D interface.
    !! In most cases, you should use the `ScalarGridInterp_1D` generic interface,
    !! rather than calling this routine directly.
    !! This routine calls hip/SELF_Lagrange.cpp:ScalarGridInterp_1D_gpu_wrapper
    !! Interpolate a scalar-1D (real) array from the control grid to the target grid.
    !! The control and target grids are the ones associated with an initialized
    !! Lagrange instance.
    !!
    !! Interpolation is applied using a series of matrix-vector multiplications, using
    !! the Lagrange class's interpolation matrix
    !!
    !! $$ \tilde{f}_{m,iel,ivar} = \sum_{i=0}^N f_{i,iel,ivar} I_{i,m} $$
    !!
    implicit none
    class(Lagrange),intent(in) :: this
    !! Lagrange class instance
    integer,intent(in) :: nvars
    !! The number of variables/functions that are interpolated
    integer,intent(in) :: nelems
    !! The number of spectral elements in the SEM grid
    real(prec),pointer,intent(in)  :: f(:,:,:)
    !! (Input) Array of function values, defined on the control grid
    real(prec),pointer,intent(inout) :: fTarget(:,:,:)
    !! (Output) Array of function values, defined on the target grid
    type(c_ptr),intent(inout) :: hipblas_handle
    ! Local
    integer(c_int) :: m
    integer(c_int) :: n
    integer(c_int) :: k
    real(c_prec) :: alpha
    integer(c_int) :: lda
    integer(c_int) :: ldb
    integer(c_int) :: ldc
    real(c_prec) :: beta

    m = this % M + 1 ! number of rows of A^T
    n = nvars*nelems ! number of columns of B
    k = this % N + 1! number of columns of A^T
    alpha = 1.0_c_prec
    lda = k ! leading dimension of A (interoplation matrix)
    ldb = k ! leading dimension of B (f)
    ldc = m ! leading dimension of C (fTarget)
    beta = 0.0_c_prec

#ifdef DOUBLE_PRECISION
    call hipblasCheck(hipblasDgemm(hipblas_handle, &
                                   HIPBLAS_OP_T,HIPBLAS_OP_N, &
                                   m,n,k,alpha, &
                                   c_loc(this % iMatrix),lda, &
                                   c_loc(f),ldb, &
                                   beta, &
                                   c_loc(fTarget),ldc))
#else
    call hipblasCheck(hipblasSgemm(hipblas_handle,&
        HIPBLAS_OP_T, HIPBLAS_OP_N, &
        m, n, k, alpha, &
        c_loc(this % iMatrix), lda, &
        c_loc(f), ldb, &
        beta, &
        c_loc(fTarget), ldc))
#endif

  end subroutine ScalarGridInterp_1D_gpu

  subroutine ScalarGridInterp_2D_cpu(this,f,fTarget,nvars,nelems)
    !! Host (CPU) implementation of the ScalarGridInterp_2D interface.
    !! In most cases, you should use the `ScalarGridInterp_2D` generic interface,
    !! rather than calling this routine directly.
    !! Interpolate a scalar-2D (real) array from the control grid to the target grid.
    !! The control and target grids are the ones associated with an initialized
    !! Lagrange instance.
    !!
    !! Interpolation is applied using a series of matrix-vector multiplications, using
    !! the Lagrange class's interpolation matrix
    !!
    !! $$ \tilde{f}_{m,n,iel,ivar} = \sum_{j=0}^N \sum_{i=0}^N f_{i,j,iel,ivar} I_{i,m} I_{j,n} $$
    !!
    implicit none
    class(Lagrange),intent(in) :: this
    !! Lagrange class instance
    integer,intent(in)     :: nvars
    !! The number of variables/functions that are interpolated
    integer,intent(in)     :: nelems
    !! The number of spectral elements in the SEM grid
    real(prec),intent(in)  :: f(1:this % N+1,1:this % N+1,1:nelems,1:nvars)
    !! (Input) Array of function values, defined on the control grid
    real(prec),intent(inout) :: fTarget(1:this % M+1,1:this % M+1,1:nelems,1:nvars)
    !! (Output) Array of function values, defined on the target grid
    ! Local
    integer :: i,j,ii,jj,iel,ivar
    real(prec) :: fi,fij

    do ivar = 1,nvars
      do iel = 1,nelems
        do j = 1,this % M+1
          do i = 1,this % M+1

            fij = 0.0_prec
            do jj = 1,this % N+1
              fi = 0.0_prec
              do ii = 1,this % N+1
                fi = fi + f(ii,jj,iel,ivar)*this % iMatrix(ii,i)
              end do
              fij = fij + fi*this % iMatrix(jj,j)
            end do
            fTarget(i,j,iel,ivar) = fij

          end do
        end do
      end do
    end do

  end subroutine ScalarGridInterp_2D_cpu

  subroutine ScalarGridInterp_2D_gpu(this,f,fInt,fTarget,nvars,nelems,hipblas_handle)
    !! Device (GPU) implementation of the ScalarGridInterp_2D interface.
    !! In most cases, you should use the `ScalarGridInterp_2D` generic interface,
    !! rather than calling this routine directly.
    !! This routine calls hip/SELF_Lagrange.cpp:ScalarGridInterp_2D_gpu_wrapper
    !! Interpolate a scalar-2D (real) array from the control grid to the target grid.
    !! The control and target grids are the ones associated with an initialized
    !! Lagrange instance.
    !!
    !! Interpolation is applied using a series of matrix-vector multiplications, using
    !! the Lagrange class's interpolation matrix
    !!
    !! $$ \tilde{f}_{m,n,iel,ivar} = \sum_{j=0}^N \sum_{i=0}^N f_{i,j,iel,ivar} I_{i,m} I_{j,n} $$
    !!
    implicit none
    class(Lagrange),intent(in) :: this
    !! Lagrange class instance
    integer,intent(in) :: nvars
    !! The number of variables/functions that are interpolated
    integer,intent(in) :: nelems
    !! The number of spectral elements in the SEM grid
    real(prec),pointer,intent(in)  :: f(:,:,:,:)
    !! (Input) Array of function values, defined on the control grid
    real(prec),pointer,intent(inout) :: fInt(:,:,:,:)
    !! (Inout) workspace array for handling intermediate values interpolated in one direction
    real(prec),pointer,intent(inout) :: fTarget(:,:,:,:)
    !! (Output) Array of function values, defined on the target grid
    type(c_ptr),intent(inout) :: hipblas_handle
    ! Local
    integer(c_int) :: m
    integer(c_int) :: n
    integer(c_int) :: k
    real(c_prec) :: alpha
    integer(c_int) :: lda
    integer(c_int) :: ldb
    integer(c_int) :: ldc
    real(c_prec) :: beta
    ! for gemvstridedbatch
    integer :: i
    integer(c_int64_t) :: strideA
    integer(c_int) :: incx
    integer(c_int64_t) :: stridex
    integer(c_int) :: incy
    integer(c_int64_t) :: stridey
    integer(c_int) :: batchCount
    integer(kind(HIPBLAS_STATUS_SUCCESS)) :: status

    m = this % M + 1 ! number of rows of A^T
    n = nvars*nelems*(this % N + 1) ! number of columns of B
    k = this % N + 1! number of columns of A^T
    alpha = 1.0_c_prec
    lda = k ! leading dimension of A (interoplation matrix)
    ldb = k ! leading dimension of B (f)
    ldc = m ! leading dimension of C (fTarget)
    beta = 0.0_c_prec

#ifdef DOUBLE_PRECISION
    ! First pass interpolates in the first quadrature dimension
    call hipblasCheck(hipblasDgemm(hipblas_handle, &
                                   HIPBLAS_OP_T,HIPBLAS_OP_N, &
                                   m,n,k,alpha, &
                                   c_loc(this % iMatrix),lda, &
                                   c_loc(f),ldb,beta, &
                                   c_loc(fInt),ldc))
#else
    ! First pass interpolates in the first quadrature dimension
    call hipblasCheck(hipblasSgemm(hipblas_handle, &
                                   HIPBLAS_OP_T,HIPBLAS_OP_N, &
                                   m,n,k,alpha, &
                                   c_loc(this % iMatrix),lda, &
                                   c_loc(f),ldb,beta, &
                                   c_loc(fInt),ldc))
#endif

    m = this % N + 1 ! number of rows of A
    n = this % M + 1 ! number of columns of A
    alpha = 1.0_c_prec
    lda = m ! leading dimension of A
    strideA = 0 ! stride for the batches of A (no stride)
    incx = this % M + 1 !
    stridex = (this % N + 1)*(this % M + 1)
    beta = 0.0_c_prec
    incy = this % M + 1
    stridey = (this % M + 1)*(this % M + 1)
    batchCount = nvars*nelems
    do i = 0,this % M
#ifdef DOUBLE_PRECISION
      call hipblasCheck(hipblasDgemvStridedBatched(hipblas_handle, &
                                                   HIPBLAS_OP_T, &
                                                   m,n,alpha, &
                                                   c_loc(this % iMatrix),lda,strideA, &
                                                   c_loc(fInt(1 + i,1,1,1)),incx,stridex,beta, &
                                                   c_loc(fTarget(1 + i,1,1,1)),incy,stridey,batchCount))
#else
      call hipblasCheck(hipblasSgemvStridedBatched(hipblas_handle, &
                                                   HIPBLAS_OP_T, &
                                                   m,n,alpha, &
                                                   c_loc(this % iMatrix),lda,strideA, &
                                                   c_loc(fInt(1 + i,1,1,1)),incx,stridex,beta, &
                                                   c_loc(fTarget(1 + i,1,1,1)),incy,stridey,batchCount))
#endif
    end do

  end subroutine ScalarGridInterp_2D_gpu

!   subroutine VectorGridInterp_2D_cpu(this,f,fTarget,nvars,nelems)
!     !! Host (CPU) implementation of the VectorGridInterp_2D interface.
!     !! In most cases, you should use the `VectorGridInterp_2D` generic interface,
!     !! rather than calling this routine directly.
!     !! Interpolate a vector-2D (real) array from the control grid to the target grid.
!     !! The control and target grids are the ones associated with an initialized
!     !! Lagrange instance.
!     !!
!     !! Interpolation is applied using a series of matrix-vector multiplications, using
!     !! the Lagrange class's interpolation matrix
!     !!
!     !! $$ \tilde{f}_{dir,m,n,iel,ivar} = \sum_{j=0}^N \sum_{i=0}^N f_{dir,i,j,iel,ivar} I_{i,m} I_{j,n} $$
!     !!
!     implicit none
!     class(Lagrange),intent(in) :: this
!     !! Lagrange class instance
!     integer,intent(in)     :: nvars
!     !! The number of variables/functions that are interpolated
!     integer,intent(in)     :: nelems
!     !! The number of spectral elements in the SEM grid
!     real(prec),intent(in)  :: f(1:2,0:this % N,0:this % N,1:nelems,1:nvars)
!     !! (Input) Array of function values, defined on the control grid
!     real(prec),intent(out) :: fTarget(1:2,0:this % M,0:this % M,1:nelems,1:nvars)
!     !! (Output) Array of function values, defined on the target grid
!     ! Local
!     integer :: i,j,ii,jj,iel,ivar
!     real(prec) :: fi(1:2)

!     do iel = 1,nelems
!       do ivar = 1,nvars
!         do j = 0,this % M
!           do i = 0,this % M

!             fTarget(1,i,j,iel,ivar) = 0.0_prec
!             fTarget(2,i,j,iel,ivar) = 0.0_prec

!             do jj = 0,this % N

!               fi(1:2) = 0.0_prec
!               do ii = 0,this % N
!                 fi(1:2) = fi(1:2) + f(1:2,ii,jj,iel,ivar)*this % iMatrix (ii,i)
!               end do

!               fTarget(1:2,i,j,iel,ivar) = fTarget(1:2,i,j,iel,ivar) + fi(1:2)*this % iMatrix (jj,j)

!             end do

!           end do
!         end do
!       end do
!     end do

!   end subroutine VectorGridInterp_2D_cpu
! !
!   subroutine VectorGridInterp_2D_gpu(this,f_dev,fTarget_dev,nvars,nelems)
!     !! Device (GPU) implementation of the VectorGridInterp_2D interface.
!     !! In most cases, you should use the `VectorGridInterp_2D` generic interface,
!     !! rather than calling this routine directly.
!     !! This routine calls hip/SELF_Lagrange.cpp:VectorGridInterp_2D_gpu_wrapper
!     !! Interpolate a vector-2D (real) array from the control grid to the target grid.
!     !! The control and target grids are the ones associated with an initialized
!     !! Lagrange instance.
!     !!
!     !! Interpolation is applied using a series of matrix-vector multiplications, using
!     !! the Lagrange class's interpolation matrix
!     !!
!     !! $$ \tilde{f}_{dir,m,n,iel,ivar} = \sum_{j=0}^N \sum_{i=0}^N f_{dir,i,j,iel,ivar} I_{i,m} I_{j,n} $$
!     !!
!     implicit none
!     class(Lagrange),intent(in) :: this
!     !! Lagrange class instance
!     integer,intent(in) :: nvars
!     !! The number of variables/functions that are interpolated
!     integer,intent(in) :: nelems
!     !! The number of spectral elements in the SEM grid
!     type(c_ptr),intent(in)  :: f_dev
!     !! (Input) Array of function values, defined on the control grid
!     type(c_ptr),intent(out) :: fTarget_dev
!     !! (Output) Array of function values, defined on the target grid

!     call VectorGridInterp_2D_gpu_wrapper(this % iMatrix % deviceData, &
!                                          f_dev,fTarget_dev, &
!                                          this % N,this % M, &
!                                          nvars,nelems)

!   end subroutine VectorGridInterp_2D_gpu

!   subroutine ScalarGridInterp_3D_cpu(this,f,fTarget,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)     :: nvars,nelems
!     real(prec),intent(in)  :: f(0:this % N,0:this % N,0:this % N,1:nelems,1:nvars)
!     real(prec),intent(out) :: fTarget(0:this % M,0:this % M,0:this % M,1:nelems,1:nvars)
!     ! Local
!     integer :: iel,ivar,i,j,k,ii,jj,kk
!     real(prec) :: fi,fij,fijk

!     do iel = 1,nelems
!       do ivar = 1,nvars
!         do k = 0,this % M
!           do j = 0,this % M
!             do i = 0,this % M

!               fijk = 0.0_prec
!               do kk = 0,this % N
!                 fij = 0.0_prec
!                 do jj = 0,this % N
!                   fi = 0.0_prec
!                   do ii = 0,this % N
!                     fi = fi + f(ii,jj,kk,iel,ivar)*this % iMatrix (ii,i)
!                   end do
!                   fij = fij + fi*this % iMatrix (jj,j)
!                 end do
!                 fijk = fijk + fij*this % iMatrix (kk,k)
!               end do
!               fTarget(i,j,k,iel,ivar) = fijk

!             end do
!           end do
!         end do
!       end do
!     end do

!   end subroutine ScalarGridInterp_3D_cpu
! !
!   subroutine ScalarGridInterp_3D_gpu(this,f_dev,fTarget_dev,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in) :: nvars,nelems
!     type(c_ptr),intent(in)  :: f_dev
!     type(c_ptr),intent(out) :: fTarget_dev

!     call ScalarGridInterp_3D_gpu_wrapper(this % iMatrix % deviceData, &
!                                          f_dev,fTarget_dev, &
!                                          this % N,this % M, &
!                                          nvars,nelems)

!   end subroutine ScalarGridInterp_3D_gpu

!   subroutine VectorGridInterp_3D_cpu(this,f,fTarget,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)     :: nvars,nelems
!     real(prec),intent(in)  :: f(1:3,0:this % N,0:this % N,0:this % N,1:nelems,1:nvars)
!     real(prec),intent(out) :: fTarget(1:3,0:this % M,0:this % M,0:this % M,1:nelems,1:nvars)
!     ! Local
!     integer :: iel,ivar,i,j,k,ii,jj,kk
!     real(prec) :: fi(1:3),fij(1:3)

!     do iel = 1,nelems
!       do ivar = 1,nvars
!         do k = 0,this % M
!           do j = 0,this % M
!             do i = 0,this % M

!               fTarget(1:3,i,j,k,iel,ivar) = 0.0_prec
!               do kk = 0,this % N
!                 fij(1:3) = 0.0_prec
!                 do jj = 0,this % N
!                   fi(1:3) = 0.0_prec
!                   do ii = 0,this % N
!                     fi(1:3) = fi(1:3) + f(1:3,ii,jj,kk,iel,ivar)*this % iMatrix (ii,i)
!                   end do
!                   fij(1:3) = fij(1:3) + fi(1:3)*this % iMatrix (jj,j)
!                 end do
!                 fTarget(1:3,i,j,k,iel,ivar) = fTarget(1:3,i,j,k,iel,ivar) + fij(1:3)*this % iMatrix (kk,k)
!               end do

!             end do
!           end do
!         end do
!       end do
!     end do

!   end subroutine VectorGridInterp_3D_cpu
! !
!   subroutine VectorGridInterp_3D_gpu(this,f_dev,fTarget_dev,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in) :: nvars,nelems
!     type(c_ptr),intent(in)  :: f_dev
!     type(c_ptr),intent(out) :: fTarget_dev

!     call VectorGridInterp_3D_gpu_wrapper(this % iMatrix % deviceData, &
!                                          f_dev,fTarget_dev, &
!                                          this % N,this % M, &
!                                          nvars,nelems)

!   end subroutine VectorGridInterp_3D_gpu

! Derivative_1D
!
!   Calculates the derivative of the Lagrange interpolant given a set of nodal function values at
!   the native interpolation nodes
!
!   Given a set of nodal values at the interpolation nodes, the derivative of a function through
!   the interpolation nodes can be estimated by
!
!                       f'_a = \sum_{i=0}^N f_{i} l'_i(\xi_a),   a=0,1,2,...,N
!
!   where l_i(\xi) are the Lagrange interpolating polynomials through the interpolation points.
!   The derivative matrix is D_{a,i} = l'_i(\xi_a) maps an array of nodal values at the interpolation
!   nodes to its estimated derivative. This routine serves as a wrapper to call either the CUDA
!   kernel (if CUDA is enabled) or the CPU version.
!
!   Usage :
!
!     TYPE(Lagrange) :: interp
!     INTEGER        :: nvars, nelems
!     REAL(prec)     :: f(0:interp % N,1:nelems,1:nvars)
!     REAL(prec)     :: derF(0:interp % N,1:nelems,1:nvars)
!
!       CALL interp % Derivative_1D( f, derF, nvars, nelems )
!
!     * If CUDA is enabled, the fnative and ftarget arrays must be CUDA device variables.
!
!   Parameters :
!
!     interp (in)
!       A previously constructed Lagrange data-structure.
!
!     f (in)
!       Array of function nodal values at the native interpolation nodes.
!
!     nvars (in)
!
!     nelems (in)
!
!     derF (out)
!      Array of derivative values at the target interpolation nodes.
!
! ================================================================================================ !

  subroutine Derivative_1D_cpu(this,f,df,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)     :: nvars,nelems
    real(prec),intent(in)  :: f(1:this % N + 1,1:nelems,1:nvars)
    real(prec),intent(out) :: df(1:this % N + 1,1:nelems,1:nvars)
    ! Local
    integer :: i,ii,iel,ivar
    real(prec) :: dfloc

    do iel = 1,nelems
      do ivar = 1,nvars
        do i = 1,this % N + 1

          dfloc = 0.0_prec
          do ii = 1,this % N + 1
            dfloc = dfloc + this % dMatrix(ii,i)*f(ii,iel,ivar)
          end do
          df(i,iel,ivar) = dfloc

        end do
      end do
    end do

  end subroutine Derivative_1D_cpu

  subroutine Derivative_1D_gpu(this,f,df,nvars,nelems,hipblas_handle)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in) :: nvars,nelems
    real(prec),pointer,intent(in)  :: f(:,:,:)
    real(prec),pointer,intent(out) :: df(:,:,:)
    type(c_ptr),intent(inout) :: hipblas_handle
    ! Local
    integer(c_int) :: m
    integer(c_int) :: n
    integer(c_int) :: k
    real(c_prec) :: alpha
    integer(c_int) :: lda
    integer(c_int) :: ldb
    integer(c_int) :: ldc
    real(c_prec) :: beta

    m = this % N + 1 ! number of rows of A^T
    n = nvars*nelems ! number of columns of B
    k = this % N + 1! number of columns of A^T
    alpha = 1.0_c_prec
    lda = k ! leading dimension of A (interoplation matrix)
    ldb = k ! leading dimension of B (f)
    ldc = m ! leading dimension of C (fTarget)
    beta = 0.0_c_prec

#ifdef DOUBLE_PRECISION
    call hipblasCheck(hipblasDgemm(hipblas_handle, &
                                   HIPBLAS_OP_T,HIPBLAS_OP_N, &
                                   m,n,k,alpha, &
                                   c_loc(this % dMatrix),lda, &
                                   c_loc(f),ldb, &
                                   beta, &
                                   c_loc(df),ldc))
#else
    call hipblasCheck(hipblasSgemm(hipblas_handle, &
                                   HIPBLAS_OP_T,HIPBLAS_OP_N, &
                                   m,n,k,alpha, &
                                   c_loc(this % dMatrix),lda, &
                                   c_loc(f),ldb, &
                                   beta, &
                                   c_loc(df),ldc))
#endif

  end subroutine Derivative_1D_gpu

!   subroutine DGDerivative_1D_cpu(this,f,bf,df,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)     :: nvars,nelems
!     real(prec),intent(in)  :: f(0:this % N,1:nelems,1:nvars)
!     real(prec),intent(in)  :: bf(1:nvars,1:2,1:nelems)
!     real(prec),intent(out) :: df(0:this % N,1:nelems,1:nvars)
!     ! Local
!     integer :: i,ii,iel,ivar

!     do iel = 1,nelems
!       do ivar = 1,nvars
!         do i = 0,this % N

!           ! Interior Derivative Matrix Application
!           df(i,iel,ivar) = 0.0_prec
!           do ii = 0,this % N
!             df(i,iel,ivar) = df(i,iel,ivar) + this % dgMatrix (ii,i)*f(ii,iel,ivar)
!           end do

!           ! Boundary Contribution
!           df(i,iel,ivar) = df(i,iel,ivar) + (bf(ivar,2,iel)*this % bMatrix (i,1) + &
!                                              bf(ivar,1,iel)*this % bMatrix (i,0))/ &
!                            this % qWeights (i)

!         end do

!       end do
!     end do

!   end subroutine DGDerivative_1D_cpu

!   subroutine DGDerivative_1D_gpu(this,f_dev,bf_dev,df_dev,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in) :: nvars,nelems
!     type(c_ptr),intent(in)  :: f_dev
!     type(c_ptr),intent(in)  :: bf_dev
!     type(c_ptr),intent(out) :: df_dev

!     call DGDerivative_1D_gpu_wrapper(this % dgMatrix % deviceData, &
!                                      this % bMatrix % deviceData, &
!                                      this % qWeights % deviceData, &
!                                      f_dev,bf_dev,df_dev, &
!                                      this % N, &
!                                      nvars,nelems)

!   end subroutine DGDerivative_1D_gpu

! ! ================================================================================================ !
! !
! ! CalculateGradient_2D
! !
! !   Calculates the gradient of a 2-D function, represented by a 2-D array of nodal values.
! !
! !   Given a set of nodal values at the interpolation nodes, the gradient of a function through
! !   the interpolation nodes can be estimated by
! !
! !                       (df/dx)_{a,b} = \sum_{i=0}^N f_{i,b} l'_i(\xi_a),   a,b=0,1,2,...,N
! !                       (df/dy)_{a,b} = \sum_{j=0}^N f_{a,j} l'_j(\xi_b),   a,b=0,1,2,...,N
! !
! !   where l_i(\xi) are the Lagrange interpolating polynomials through the interpolation points.
! !   The derivative matrix is D_{a,i} = l'_i(\xi_a) maps an array of nodal values at the interpolation
! !   nodes to its estimated derivative. This routine serves as a wrapper to call either the CUDA
! !   kernel (if CUDA is enabled) or the CPU version.
! !
! !   Usage :
! !
! !     TYPE(Lagrange) :: interp
! !     INTEGER        :: nvars, nelems
! !     REAL(prec)     :: f(0:interp % N,0:interp % N,1:nelems,1:nvars)
! !     REAL(prec)     :: gradF(1:2,0:interp % N,0:interp % N,1:nelems,1:nvars)
! !
! !       CALL interp % CalculateGradient_2D( f, gradF, nvars, nelems )
! !
! !     * If CUDA is enabled, the fnative and ftarget arrays must be CUDA device variables.
! !
! !   Parameters :
! !
! !     interp (in)
! !       A previously constructed Lagrange data-structure.
! !
! !     f (in)
! !       Array of function nodal values at the native interpolation nodes.
! !
! !     nvars (in)
! !
! !     nelems (in)
! !
! !     gradF (out)
! !      Array of derivative values at the target interpolation nodes.
! !
! ! ================================================================================================ !
! !
!   subroutine ScalarGradient_2D_cpu(this,f,gradF,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)     :: nvars,nelems
!     real(prec),intent(in)  :: f(0:this % N,0:this % N,1:nelems,1:nvars)
!     real(prec),intent(out) :: gradF(1:2,0:this % N,0:this % N,1:nelems,1:nvars)
!     ! Local
!     integer    :: i,j,ii,iel,ivar

!     do iel = 1,nelems
!       do ivar = 1,nvars
!         do j = 0,this % N
!           do i = 0,this % N

!             gradF(1,i,j,iel,ivar) = 0.0_prec
!             gradF(2,i,j,iel,ivar) = 0.0_prec
!             do ii = 0,this % N
!               gradF(1,i,j,iel,ivar) = gradF(1,i,j,iel,ivar) + this % dMatrix (ii,i)*f(ii,j,iel,ivar)
!               gradF(2,i,j,iel,ivar) = gradF(2,i,j,iel,ivar) + this % dMatrix (ii,j)*f(i,ii,iel,ivar)
!             end do

!           end do
!         end do
!       end do
!     end do

!   end subroutine ScalarGradient_2D_cpu

!   subroutine ScalarGradient_2D_gpu(this,f_dev,gradF_dev,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)         :: nvars,nelems
!     type(c_ptr),intent(in)     :: f_dev
!     type(c_ptr),intent(out)    :: gradF_dev

!     call ScalarGradient_2D_gpu_wrapper(this % dMatrix % deviceData, &
!                                        f_dev,gradF_dev,this % N, &
!                                        nvars,nelems)

!   end subroutine ScalarGradient_2D_gpu
! !
! !
!   ! SUBROUTINE ScalarDGGradient_2D_cpu(this,f,bf,gradF,nvars,nelems)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)     :: nvars,nelems
!   !   REAL(prec),INTENT(in)  :: f(0:this % N,0:this % N,1:nelems,1:nvars)
!   !   REAL(prec),INTENT(in)  :: bf(0:this % N,1:nvars,1:4,1:nelems)
!   !   REAL(prec),INTENT(out) :: gradF(1:2,0:this % N,0:this % N,1:nelems,1:nvars)
!   !   ! Local
!   !   INTEGER    :: i,j,ii,iel,ivar

!   !   DO iel = 1,nelems
!   !     DO ivar = 1,nvars
!   !       DO j = 0,this % N
!   !         DO i = 0,this % N

!   !           gradF(1,i,j,iel,ivar) = 0.0_prec
!   !           gradF(2,i,j,iel,ivar) = 0.0_prec
!   !           DO ii = 0,this % N
!   !             gradF(1,i,j,iel,ivar) = gradF(1,i,j,iel,ivar) + this % dgMatrix (ii,i)*f(ii,j,iel,ivar)
!   !             gradF(2,i,j,iel,ivar) = gradF(2,i,j,iel,ivar) + this % dgMatrix (ii,j)*f(i,ii,iel,ivar)
!   !           END DO

!   !           ! Boundary Contribution
!   !           gradF(1,i,j,iel,ivar) = gradF(1,i,j,iel,ivar) + (bf(j,ivar,2,iel)*this % bMatrix (i,1) + &
!   !                                                            bf(j,ivar,4,iel)*this % bMatrix (i,0))/ &
!   !                                   this % qWeights (i)

!   !           gradF(2,i,j,iel,ivar) = gradF(2,i,j,iel,ivar) + (bf(i,ivar,3,iel)*this % bMatrix (j,1) + &
!   !                                                            bf(i,ivar,1,iel)*this % bMatrix (j,0))/ &
!   !                                   this % qWeights (j)

!   !         END DO
!   !       END DO
!   !     END DO
!   !   END DO

!   ! END SUBROUTINE ScalarDGGradient_2D_cpu

!   ! SUBROUTINE ScalarDGGradient_2D_gpu(this,f_dev,bf_dev,gradF_dev,nvars,nelems)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)         :: nvars,nelems
!   !   TYPE(c_ptr),INTENT(in)     :: f_dev
!   !   TYPE(c_ptr),INTENT(in)     :: bf_dev
!   !   TYPE(c_ptr),INTENT(out)    :: gradF_dev

!   !   CALL ScalarDGGradient_2D_gpu_wrapper(this % dgMatrix % deviceData, &
!   !                                        this % bMatrix % deviceData, &
!   !                                        this % qWeights % deviceData, &
!   !                                        f_dev,bf_dev,gradF_dev,this % N, &
!   !                                        nvars,nelems)

!   ! END SUBROUTINE ScalarDGGradient_2D_gpu

!   subroutine VectorGradient_2D_cpu(this,f,gradF,nvars,nelems)
!     !
!     ! Input : Vector(1:2,...)
!     ! Output : Tensor(1:2,1:2,....)
!     !          > Tensor(1,1) = d/ds1( Vector(1,...) )
!     !          > Tensor(2,1) = d/ds1( Vector(2,...) )
!     !          > Tensor(1,2) = d/ds2( Vector(1,...) )
!     !          > Tensor(2,2) = d/ds2( Vector(2,...) )
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)     :: nvars,nelems
!     real(prec),intent(in)  :: f(1:2,0:this % N,0:this % N,1:nelems,1:nvars)
!     real(prec),intent(out) :: gradF(1:2,1:2,0:this % N,0:this % N,1:nelems,1:nvars)
!     ! Local
!     integer    :: i,j,ii,iel,ivar
!     real(prec) :: gf(1:2,1:2)

!     do iel = 1,nelems
!       do ivar = 1,nvars
!         do j = 0,this % N
!           do i = 0,this % N

!             gf(1,1) = 0.0_prec
!             gf(2,1) = 0.0_prec
!             gf(1,2) = 0.0_prec
!             gf(2,2) = 0.0_prec
!             do ii = 0,this % N
!               gf(1,1) = gf(1,1) + this % dMatrix (ii,i)*f(1,ii,j,iel,ivar)
!               gf(2,1) = gf(2,1) + this % dMatrix (ii,i)*f(2,ii,j,iel,ivar)
!               gf(1,2) = gf(1,2) + this % dMatrix (ii,j)*f(1,i,ii,iel,ivar)
!               gf(2,2) = gf(2,2) + this % dMatrix (ii,j)*f(2,i,ii,iel,ivar)
!             end do
!             gradF(1:2,1:2,i,j,iel,ivar) = gf(1:2,1:2)

!           end do
!         end do
!       end do
!     end do

!   end subroutine VectorGradient_2D_cpu

!   subroutine VectorGradient_2D_gpu(this,f_dev,gradF_dev,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)         :: nvars,nelems
!     type(c_ptr),intent(in)     :: f_dev
!     type(c_ptr),intent(out)    :: gradF_dev

!     call VectorGradient_2D_gpu_wrapper(this % dMatrix % deviceData, &
!                                        f_dev,gradF_dev,this % N, &
!                                        nvars,nelems)

!   end subroutine VectorGradient_2D_gpu

!   ! SUBROUTINE VectorDGGradient_2D_cpu(this,f,bf,gradF,nvars,nelems)
!   !   !
!   !   ! Input : Vector(1:2,...)
!   !   ! Output : Tensor(1:2,1:2,....)
!   !   !          > Tensor(1,1) = d/ds1( Vector(1,...) )
!   !   !          > Tensor(2,1) = d/ds1( Vector(2,...) )
!   !   !          > Tensor(1,2) = d/ds2( Vector(1,...) )
!   !   !          > Tensor(2,2) = d/ds2( Vector(2,...) )
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)     :: nvars,nelems
!   !   REAL(prec),INTENT(in)  :: f(1:2,0:this % N,0:this % N,1:nelems,1:nvars)
!   !   REAL(prec),INTENT(in)  :: bf(1:2,0:this % N,1:nvars,1:4,1:nelems)
!   !   REAL(prec),INTENT(out) :: gradF(1:2,1:2,0:this % N,0:this % N,1:nelems,1:nvars)
!   !   ! Local
!   !   INTEGER    :: i,j,ii,iel,ivar

!   !   DO iel = 1,nelems
!   !     DO ivar = 1,nvars
!   !       DO j = 0,this % N
!   !         DO i = 0,this % N

!   !           gradF(1,1,i,j,iel,ivar) = 0.0_prec
!   !           gradF(2,1,i,j,iel,ivar) = 0.0_prec
!   !           gradF(1,2,i,j,iel,ivar) = 0.0_prec
!   !           gradF(2,2,i,j,iel,ivar) = 0.0_prec
!   !           DO ii = 0,this % N
!   !             gradF(1,1,i,j,iel,ivar) = gradF(1,1,i,j,iel,ivar) + this % dgMatrix (ii,i)*f(1,ii,j,iel,ivar)
!   !             gradF(2,1,i,j,iel,ivar) = gradF(2,1,i,j,iel,ivar) + this % dgMatrix (ii,i)*f(2,ii,j,iel,ivar)
!   !             gradF(1,2,i,j,iel,ivar) = gradF(1,2,i,j,iel,ivar) + this % dgMatrix (ii,j)*f(1,i,ii,iel,ivar)
!   !             gradF(2,2,i,j,iel,ivar) = gradF(2,2,i,j,iel,ivar) + this % dgMatrix (ii,j)*f(2,i,ii,iel,ivar)
!   !           END DO
!   !           gradF(1,1,i,j,iel,ivar) = gradF(1,1,i,j,iel,ivar) + (this % bMatrix (i,1)*bf(1,j,ivar,2,iel) + &
!   !                                                                this % bMatrix (i,0)*bf(1,j,ivar,4,iel))/ &
!   !                                     this % qWeights (i)

!   !           gradF(2,1,i,j,iel,ivar) = gradF(2,1,i,j,iel,ivar) + (this % bMatrix (i,1)*bf(2,j,ivar,2,iel) + &
!   !                                                                this % bMatrix (i,0)*bf(2,j,ivar,4,iel))/ &
!   !                                     this % qWeights (i)

!   !           gradF(1,2,i,j,iel,ivar) = gradF(1,2,i,j,iel,ivar) + (this % bMatrix (j,1)*bf(1,i,ivar,3,iel) + &
!   !                                                                this % bMatrix (j,0)*bf(1,i,ivar,1,iel))/ &
!   !                                     this % qWeights (j)

!   !           gradF(2,2,i,j,iel,ivar) = gradF(2,2,i,j,iel,ivar) + (this % bMatrix (j,1)*bf(2,i,ivar,3,iel) + &
!   !                                                                this % bMatrix (j,0)*bf(2,i,ivar,1,iel))/ &
!   !                                     this % qWeights (j)

!   !         END DO
!   !       END DO
!   !     END DO
!   !   END DO

!   ! END SUBROUTINE VectorDGGradient_2D_cpu

!   ! SUBROUTINE VectorDGGradient_2D_gpu(this,f_dev,bf_dev,gradF_dev,nvars,nelems)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)         :: nvars,nelems
!   !   TYPE(c_ptr),INTENT(in)     :: f_dev
!   !   TYPE(c_ptr),INTENT(in)     :: bf_dev
!   !   TYPE(c_ptr),INTENT(out)    :: gradF_dev

!   !   CALL VectorDGGradient_2D_gpu_wrapper(this % dMatrix % deviceData, &
!   !                                        this % bMatrix % deviceData, &
!   !                                        this % qWeights % deviceData, &
!   !                                        f_dev,bf_dev,gradF_dev,this % N, &
!   !                                        nvars,nelems)

!   ! END SUBROUTINE VectorDGGradient_2D_gpu

!   subroutine VectorDivergence_2D_cpu(this,f,dF,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)     :: nvars,nelems
!     real(prec),intent(in)  :: f(1:2,0:this % N,0:this % N,1:nelems,1:nvars)
!     real(prec),intent(out) :: dF(0:this % N,0:this % N,1:nelems,1:nvars)
!     ! Local
!     integer    :: i,j,ii,iel,ivar

!     do iel = 1,nelems
!       do ivar = 1,nvars
!         do j = 0,this % N
!           do i = 0,this % N

!             dF(i,j,iel,ivar) = 0.0_prec
!             do ii = 0,this % N
!               dF(i,j,iel,ivar) = dF(i,j,iel,ivar) + this % dMatrix (ii,i)*f(1,ii,j,iel,ivar) + &
!                                  this % dMatrix (ii,j)*f(2,i,ii,iel,ivar)
!             end do

!           end do
!         end do
!       end do
!     end do

!   end subroutine VectorDivergence_2D_cpu

!   subroutine VectorDivergence_2D_gpu(this,f_dev,dF_dev,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)         :: nvars,nelems
!     type(c_ptr),intent(in)     :: f_dev
!     type(c_ptr),intent(out)    :: dF_dev

!     call VectorDivergence_2D_gpu_wrapper(this % dMatrix % deviceData, &
!                                          f_dev,dF_dev,this % N, &
!                                          nvars,nelems)

!   end subroutine VectorDivergence_2D_gpu

!   subroutine VectorDGDivergence_2D_cpu(this,f,bF,dF,nvars,nelems)
!     ! Assumes bF is the vector component in the direction normal to the boundary
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)     :: nvars,nelems
!     real(prec),intent(in)  :: f(1:2,0:this % N,0:this % N,1:nelems,1:nvars)
!     real(prec),intent(in)  :: bF(0:this % N,1:nvars,1:4,1:nelems)
!     real(prec),intent(out) :: dF(0:this % N,0:this % N,1:nelems,1:nvars)
!     ! Local
!     real(prec) :: dfLoc
!     integer    :: i,j,ii,iel,ivar

!     do iel = 1,nelems
!       do ivar = 1,nvars
!         do j = 0,this % N
!           do i = 0,this % N

!             dfLoc = 0.0_prec
!             do ii = 0,this % N
!               dfLoc = dfLoc + this % dgMatrix (ii,i)*f(1,ii,j,iel,ivar) + &
!                       this % dgMatrix (ii,j)*f(2,i,ii,iel,ivar)
!             end do

!             dfLoc = dfLoc + (this % bMatrix (i,1)*bF(j,ivar,2,iel) + &
!                              this % bMatrix (i,0)*bF(j,ivar,4,iel))/ &
!                     this % qWeights (i) + &
!                     (this % bMatrix (j,1)*bF(i,ivar,3,iel) + &
!                      this % bMatrix (j,0)*bF(i,ivar,1,iel))/ &
!                     this % qWeights (j)
!             dF(i,j,iel,ivar) = dFLoc

!           end do
!         end do
!       end do
!     end do

!   end subroutine VectorDGDivergence_2D_cpu

!   subroutine VectorDGDivergence_2D_gpu(this,f_dev,bF_dev,dF_dev,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)         :: nvars,nelems
!     type(c_ptr),intent(in)     :: f_dev
!     type(c_ptr),intent(in)     :: bF_dev
!     type(c_ptr),intent(out)    :: dF_dev

!     call VectorDGDivergence_2D_gpu_wrapper(this % dgMatrix % deviceData, &
!                                            this % bMatrix % deviceData, &
!                                            this % qWeights % deviceData, &
!                                            f_dev,bF_dev,dF_dev,this % N, &
!                                            nvars,nelems)

!   end subroutine VectorDGDivergence_2D_gpu

!   ! SUBROUTINE VectorCurl_2D_cpu(this,f,dF,nvars,nelems)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)     :: nvars,nelems
!   !   REAL(prec),INTENT(in)  :: f(1:2,0:this % N,0:this % N,1:nelems,1:nvars)
!   !   REAL(prec),INTENT(out) :: dF(0:this % N,0:this % N,1:nelems,1:nvars)
!   !   ! Local
!   !   INTEGER    :: i,j,ii,iel,ivar

!   !   DO iel = 1,nelems
!   !     DO ivar = 1,nvars
!   !       DO j = 0,this % N
!   !         DO i = 0,this % N

!   !           dF(i,j,iel,ivar) = 0.0_prec
!   !           DO ii = 0,this % N
!   !             dF(i,j,iel,ivar) = dF(i,j,iel,ivar) + this % dMatrix (ii,j)*f(1,i,ii,iel,ivar) - &
!   !                                this % dMatrix (ii,i)*f(2,ii,j,iel,ivar)
!   !           END DO

!   !         END DO
!   !       END DO
!   !     END DO
!   !   END DO

!   ! END SUBROUTINE VectorCurl_2D_cpu

!   ! SUBROUTINE VectorCurl_2D_gpu(this,f_dev,dF_dev,nvars,nelems)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)         :: nvars,nelems
!   !   TYPE(c_ptr),INTENT(in)     :: f_dev
!   !   TYPE(c_ptr),INTENT(out)    :: dF_dev

!   !   CALL VectorCurl_2D_gpu_wrapper(this % dMatrix % deviceData, &
!   !                                  f_dev,dF_dev,this % N, &
!   !                                  nvars,nelems)

!   ! END SUBROUTINE VectorCurl_2D_gpu

!   ! SUBROUTINE P2VectorDivergence_2D_cpu(this,f,dF,nvars,nelems)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)     :: nvars,nelems
!   !   REAL(prec),INTENT(in)  :: f(1:2,0:this % N,0:this % N,0:this % N,1:nelems,1:nvars)
!   !   REAL(prec),INTENT(out) :: dF(0:this % N,0:this % N,1:nelems,1:nvars)
!   !   ! Local
!   !   INTEGER    :: i,j,n,iel,ivar
!   !   REAL(prec) :: dfloc

!   !   DO iel = 1,nelems
!   !     DO ivar = 1,nvars
!   !       DO j = 0,this % N
!   !         DO i = 0,this % N

!   !           dfloc = 0.0_prec
!   !           DO n = 0,this % N
!   !             dfloc = dfloc + this % dMatrix (n,i)*f(1,n,i,j,iel,ivar) + &
!   !                             this % dMatrix (n,j)*f(2,n,i,j,iel,ivar)
!   !           END DO

!   !           dF(i,j,iel,ivar) = 2.0_prec*dfloc

!   !         END DO
!   !       END DO
!   !     END DO
!   !   END DO

!   ! END SUBROUTINE P2VectorDivergence_2D_cpu

!   ! SUBROUTINE P2VectorDivergence_2D_gpu(this,f_dev,dF_dev,nvars,nelems)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)         :: nvars,nelems
!   !   TYPE(c_ptr),INTENT(in)     :: f_dev
!   !   TYPE(c_ptr),INTENT(out)    :: dF_dev

!   !   CALL P2VectorDivergence_2D_gpu_wrapper(this % dMatrix % deviceData, &
!   !                                        f_dev,dF_dev,this % N, &
!   !                                        nvars,nelems)

!   ! END SUBROUTINE P2VectorDivergence_2D_gpu

!   ! SUBROUTINE P2VectorDGDivergence_2D_cpu(this,f,bF,dF,nvars,nelems)
!   !   ! Assumes bF is the vector component in the direction normal to the boundary
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)     :: nvars,nelems
!   !   REAL(prec),INTENT(in)  :: f(1:2,0:this % N,0:this % N,0:this % N,1:nelems,1:nvars)
!   !   REAL(prec),INTENT(in)  :: bF(0:this % N,1:nvars,1:4,1:nelems)
!   !   REAL(prec),INTENT(out) :: dF(0:this % N,0:this % N,1:nelems,1:nvars)
!   !   ! Local
!   !   REAL(prec) :: dfLoc
!   !   INTEGER    :: i,j,n,iel,ivar

!   !   DO iel = 1,nelems
!   !     DO ivar = 1,nvars
!   !       DO j = 0,this % N
!   !         DO i = 0,this % N

!   !           dfLoc = 0.0_prec
!   !           DO n = 0,this % N
!   !             dfLoc = dfLoc + this % dgMatrix (n,i)*f(1,n,i,j,iel,ivar) + &
!   !                             this % dgMatrix (n,j)*f(2,n,i,j,iel,ivar)
!   !           END DO

!   !           dfLoc = dfLoc + (this % bMatrix (i,1)*bF(j,ivar,2,iel) + &
!   !                            this % bMatrix (i,0)*bF(j,ivar,4,iel))/ &
!   !                              this % qWeights (i) + &
!   !                           (this % bMatrix (j,1)*bF(i,ivar,3,iel) + &
!   !                            this % bMatrix (j,0)*bF(i,ivar,1,iel))/ &
!   !                              this % qWeights (j)
!   !           dF(i,j,iel,ivar) = dFLoc

!   !         END DO
!   !       END DO
!   !     END DO
!   !   END DO

!   ! END SUBROUTINE P2VectorDGDivergence_2D_cpu

!   ! SUBROUTINE P2VectorDGDivergence_2D_gpu(this,f_dev,bF_dev,dF_dev,nvars,nelems)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)         :: nvars,nelems
!   !   TYPE(c_ptr),INTENT(in)     :: f_dev
!   !   TYPE(c_ptr),INTENT(in)     :: bF_dev
!   !   TYPE(c_ptr),INTENT(out)    :: dF_dev

!   !   CALL P2VectorDGDivergence_2D_gpu_wrapper(this % dgMatrix % deviceData, &
!   !                                          this % bMatrix % deviceData, &
!   !                                          this % qWeights % deviceData, &
!   !                                          f_dev,bF_dev,dF_dev,this % N, &
!   !                                          nvars,nelems)

!   ! END SUBROUTINE P2VectorDGDivergence_2D_gpu

!   ! SUBROUTINE TensorDivergence_2D_cpu(this,f,dF,nvars,nelems)
!   !   ! Note that the divergence is taken over the first dimension (row dimension) of the tensor matrix
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)     :: nvars,nelems
!   !   REAL(prec),INTENT(in)  :: f(1:2,1:2,0:this % N,0:this % N,1:nelems,1:nvars)
!   !   REAL(prec),INTENT(out) :: dF(1:2,0:this % N,0:this % N,1:nelems,1:nvars)
!   !   ! Local
!   !   INTEGER    :: i,j,ii,iel,ivar

!   !   DO iel = 1,nelems
!   !     DO ivar = 1,nvars
!   !       DO j = 0,this % N
!   !         DO i = 0,this % N

!   !           dF(1,i,j,iel,ivar) = 0.0_prec
!   !           dF(2,i,j,iel,ivar) = 0.0_prec
!   !           DO ii = 0,this % N
!   !             dF(1,i,j,iel,ivar) = dF(1,i,j,iel,ivar) + this % dMatrix (ii,i)*f(1,1,ii,j,iel,ivar) + &
!   !                                  this % dMatrix (ii,j)*f(2,1,i,ii,iel,ivar)
!   !             dF(2,i,j,iel,ivar) = dF(2,i,j,iel,ivar) + this % dMatrix (ii,i)*f(1,2,ii,j,iel,ivar) + &
!   !                                  this % dMatrix (ii,j)*f(2,2,i,ii,iel,ivar)
!   !           END DO

!   !         END DO
!   !       END DO
!   !     END DO
!   !   END DO

!   ! END SUBROUTINE TensorDivergence_2D_cpu

!   ! SUBROUTINE TensorDivergence_2D_gpu(this,f_dev,dF_dev,nvars,nelems)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)         :: nvars,nelems
!   !   TYPE(c_ptr),INTENT(in)     :: f_dev
!   !   TYPE(c_ptr),INTENT(out)    :: dF_dev

!   !   CALL TensorDivergence_2D_gpu_wrapper(this % dMatrix % deviceData, &
!   !                                        f_dev,dF_dev,this % N, &
!   !                                        nvars,nelems)

!   ! END SUBROUTINE TensorDivergence_2D_gpu

!   ! SUBROUTINE TensorDGDivergence_2D_cpu(this,f,bF,dF,nvars,nelems)
!   !   ! Note that the divergence is taken over the first dimension (row dimension) of the tensor matrix
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)     :: nvars,nelems
!   !   REAL(prec),INTENT(in)  :: f(1:2,1:2,0:this % N,0:this % N,1:nelems,1:nvars)
!   !   REAL(prec),INTENT(in)  :: bf(1:2,1:2,0:this % N,1:nvars,1:4,1:nelems)
!   !   REAL(prec),INTENT(out) :: dF(1:2,0:this % N,0:this % N,1:nelems,1:nvars)
!   !   ! Local
!   !   INTEGER    :: i,j,ii,iel,ivar

!   !   DO iel = 1,nelems
!   !     DO ivar = 1,nvars
!   !       DO j = 0,this % N
!   !         DO i = 0,this % N

!   !           dF(1,i,j,iel,ivar) = 0.0_prec
!   !           dF(2,i,j,iel,ivar) = 0.0_prec
!   !           DO ii = 0,this % N
!   !             dF(1,i,j,iel,ivar) = dF(1,i,j,iel,ivar) + this % dgMatrix (ii,i)*f(1,1,ii,j,iel,ivar) + &
!   !                                  this % dgMatrix (ii,j)*f(2,1,i,ii,iel,ivar)
!   !             dF(2,i,j,iel,ivar) = dF(2,i,j,iel,ivar) + this % dgMatrix (ii,i)*f(1,2,ii,j,iel,ivar) + &
!   !                                  this % dgMatrix (ii,j)*f(2,2,i,ii,iel,ivar)
!   !           END DO

!   !           dF(1,i,j,iel,ivar) = dF(1,i,j,iel,ivar) + (this % bMatrix (i,1)*bf(1,1,j,ivar,2,iel) + &
!   !                                                      this % bMatrix (i,0)*bf(1,1,j,ivar,4,iel))/ &
!   !                                this % qWeights (i) + &
!   !                                (this % bMatrix (j,1)*bf(2,1,i,ivar,3,iel) + &
!   !                                 this % bMatrix (j,0)*bf(2,1,i,ivar,1,iel))/ &
!   !                                this % qWeights (j)

!   !           dF(2,i,j,iel,ivar) = dF(2,i,j,iel,ivar) + (this % bMatrix (i,1)*bf(1,2,j,ivar,2,iel) + &
!   !                                                      this % bMatrix (i,0)*bf(1,2,j,ivar,4,iel))/ &
!   !                                this % qWeights (i) + &
!   !                                (this % bMatrix (j,1)*bf(2,2,i,ivar,3,iel) + &
!   !                                 this % bMatrix (j,0)*bf(2,2,i,ivar,1,iel))/ &
!   !                                this % qWeights (j)
!   !         END DO
!   !       END DO
!   !     END DO
!   !   END DO

!   ! END SUBROUTINE TensorDGDivergence_2D_cpu

!   ! SUBROUTINE TensorDGDivergence_2D_gpu(this,f_dev,bF_dev,dF_dev,nvars,nelems)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)         :: nvars,nelems
!   !   TYPE(c_ptr),INTENT(in)     :: f_dev
!   !   TYPE(c_ptr),INTENT(in)     :: bf_dev
!   !   TYPE(c_ptr),INTENT(out)    :: dF_dev

!   !   CALL TensorDGDivergence_2D_gpu_wrapper(this % dgMatrix % deviceData, &
!   !                                          this % bMatrix % deviceData, &
!   !                                          this % qWeights % deviceData, &
!   !                                          f_dev,bF_dev,dF_dev,this % N, &
!   !                                          nvars,nelems)

!   ! END SUBROUTINE TensorDGDivergence_2D_gpu

!   subroutine ScalarGradient_3D_cpu(this,f,gradF,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)     :: nvars,nelems
!     real(prec),intent(in)  :: f(0:this % N,0:this % N,0:this % N,1:nelems,1:nvars)
!     real(prec),intent(out) :: gradF(1:3,0:this % N,0:this % N,0:this % N,1:nelems,1:nvars)
!     ! Local
!     integer    :: i,j,k,ii,iel,ivar
!     real(prec) :: gf(1:3)

!     do iel = 1,nelems
!       do ivar = 1,nvars
!         do k = 0,this % N
!           do j = 0,this % N
!             do i = 0,this % N

!               gF(1) = 0.0_prec
!               gF(2) = 0.0_prec
!               gF(3) = 0.0_prec
!               do ii = 0,this % N
!                 gF(1) = gF(1) + this % dMatrix (ii,i)*f(ii,j,k,iel,ivar)
!                 gF(2) = gF(2) + this % dMatrix (ii,j)*f(i,ii,k,iel,ivar)
!                 gF(3) = gF(3) + this % dMatrix (ii,k)*f(i,j,ii,iel,ivar)
!               end do

!               gradF(1,i,j,k,iel,ivar) = gF(1)
!               gradF(2,i,j,k,iel,ivar) = gF(2)
!               gradF(3,i,j,k,iel,ivar) = gF(3)

!             end do
!           end do
!         end do
!       end do
!     end do

!   end subroutine ScalarGradient_3D_cpu

!   subroutine ScalarGradient_3D_gpu(this,f_dev,gradF_dev,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)         :: nvars,nelems
!     type(c_ptr),intent(in)     :: f_dev
!     type(c_ptr),intent(out)    :: gradF_dev

!     call ScalarGradient_3D_gpu_wrapper(this % dMatrix % deviceData, &
!                                        f_dev,gradF_dev,this % N, &
!                                        nvars,nelems)

!   end subroutine ScalarGradient_3D_gpu
! !
!   subroutine VectorGradient_3D_cpu(this,f,gradF,nvars,nelems)
!     !
!     ! Input : Vector(1:3,...)
!     ! Output : Tensor(1:3,1:3,....)
!     !          > Tensor(1,1) = d/ds1( Vector(1,...) )
!     !          > Tensor(2,1) = d/ds1( Vector(2,...) )
!     !          > Tensor(3,1) = d/ds1( Vector(3,...) )
!     !          > Tensor(1,2) = d/ds2( Vector(1,...) )
!     !          > Tensor(2,2) = d/ds2( Vector(2,...) )
!     !          > Tensor(3,2) = d/ds2( Vector(3,...) )
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)     :: nvars,nelems
!     real(prec),intent(in)  :: f(1:3,0:this % N,0:this % N,0:this % N,1:nelems,1:nvars)
!     real(prec),intent(out) :: gradF(1:3,1:3,0:this % N,0:this % N,0:this % N,1:nelems,1:nvars)
!     ! Local
!     integer    :: i,j,k,ii,iel,ivar
!     real(prec) :: gF(1:3,1:3)

!     do iel = 1,nelems
!       do ivar = 1,nvars
!         do k = 0,this % N
!           do j = 0,this % N
!             do i = 0,this % N

!               gF = 0.0_prec
!               do ii = 0,this % N
!                 gF(1,1) = gF(1,1) + this % dMatrix (ii,i)*f(1,ii,j,k,iel,ivar)
!                 gF(2,1) = gF(2,1) + this % dMatrix (ii,i)*f(2,ii,j,k,iel,ivar)
!                 gF(3,1) = gF(3,1) + this % dMatrix (ii,i)*f(3,ii,j,k,iel,ivar)
!                 gF(1,2) = gF(1,2) + this % dMatrix (ii,j)*f(1,i,ii,k,iel,ivar)
!                 gF(2,2) = gF(2,2) + this % dMatrix (ii,j)*f(2,i,ii,k,iel,ivar)
!                 gF(3,2) = gF(3,2) + this % dMatrix (ii,j)*f(3,i,ii,k,iel,ivar)
!                 gF(1,3) = gF(1,3) + this % dMatrix (ii,k)*f(1,i,j,ii,iel,ivar)
!                 gF(2,3) = gF(2,3) + this % dMatrix (ii,k)*f(2,i,j,ii,iel,ivar)
!                 gF(3,3) = gF(3,3) + this % dMatrix (ii,k)*f(3,i,j,ii,iel,ivar)
!               end do

!               gradF(1,1,i,j,k,iel,ivar) = gF(1,1)
!               gradF(2,1,i,j,k,iel,ivar) = gF(2,1)
!               gradF(3,1,i,j,k,iel,ivar) = gF(3,1)
!               gradF(1,2,i,j,k,iel,ivar) = gF(1,2)
!               gradF(2,2,i,j,k,iel,ivar) = gF(2,2)
!               gradF(3,2,i,j,k,iel,ivar) = gF(3,2)
!               gradF(1,3,i,j,k,iel,ivar) = gF(1,3)
!               gradF(2,3,i,j,k,iel,ivar) = gF(2,3)
!               gradF(3,3,i,j,k,iel,ivar) = gF(3,3)

!             end do
!           end do
!         end do
!       end do
!     end do

!   end subroutine VectorGradient_3D_cpu

!   subroutine VectorGradient_3D_gpu(this,f_dev,gradF_dev,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)         :: nvars,nelems
!     type(c_ptr),intent(in)     :: f_dev
!     type(c_ptr),intent(out)    :: gradF_dev

!     call VectorGradient_3D_gpu_wrapper(this % dMatrix % deviceData, &
!                                        f_dev,gradF_dev,this % N, &
!                                        nvars,nelems)

!   end subroutine VectorGradient_3D_gpu

!   subroutine VectorDivergence_3D_cpu(this,f,dF,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)     :: nvars,nelems
!     real(prec),intent(in)  :: f(1:3,0:this % N,0:this % N,0:this % N,1:nelems,1:nvars)
!     real(prec),intent(out) :: dF(0:this % N,0:this % N,0:this % N,1:nelems,1:nvars)
!     ! Local
!     integer    :: i,j,k,ii,iel,ivar

!     do iel = 1,nelems
!       do ivar = 1,nvars
!         do k = 0,this % N
!           do j = 0,this % N
!             do i = 0,this % N

!               dF(i,j,k,iel,ivar) = 0.0_prec
!               do ii = 0,this % N
!                 dF(i,j,k,iel,ivar) = dF(i,j,k,iel,ivar) + this % dMatrix (ii,i)*f(1,ii,j,k,iel,ivar) + &
!                                      this % dMatrix (ii,j)*f(2,i,ii,k,iel,ivar) + &
!                                      this % dMatrix (ii,k)*f(3,i,j,ii,iel,ivar)
!               end do

!             end do
!           end do
!         end do
!       end do
!     end do

!   end subroutine VectorDivergence_3D_cpu

!   subroutine VectorDivergence_3D_gpu(this,f_dev,dF_dev,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)         :: nvars,nelems
!     type(c_ptr),intent(in)     :: f_dev
!     type(c_ptr),intent(out)    :: dF_dev

!     call VectorDivergence_3D_gpu_wrapper(this % dMatrix % deviceData, &
!                                          f_dev,dF_dev,this % N, &
!                                          nvars,nelems)

!   end subroutine VectorDivergence_3D_gpu

!   subroutine VectorDGDivergence_3D_cpu(this,f,bF,dF,nvars,nelems)
!     ! Assumes bF is the vector component in the direction normal to the element boundaries
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)     :: nvars,nelems
!     real(prec),intent(in)  :: f(1:3,0:this % N,0:this % N,0:this % N,1:nelems,1:nvars)
!     real(prec),intent(in)  :: bf(0:this % N,0:this % N,1:nvars,1:6,1:nelems)
!     real(prec),intent(out) :: dF(0:this % N,0:this % N,0:this % N,1:nelems,1:nvars)
!     ! Local
!     integer    :: i,j,k,ii,iel,ivar

!     do iel = 1,nelems
!       do ivar = 1,nvars
!         do k = 0,this % N
!           do j = 0,this % N
!             do i = 0,this % N

!               dF(i,j,k,iel,ivar) = 0.0_prec
!               do ii = 0,this % N
!                 dF(i,j,k,iel,ivar) = dF(i,j,k,iel,ivar) + this % dgMatrix (ii,i)*f(1,ii,j,k,iel,ivar) + &
!                                      this % dgMatrix (ii,j)*f(2,i,ii,k,iel,ivar) + &
!                                      this % dgMatrix (ii,k)*f(3,i,j,ii,iel,ivar)
!               end do

!               dF(i,j,k,iel,ivar) = dF(i,j,k,iel,ivar) + (this % bMatrix (i,1)*bF(j,k,ivar,3,iel) + & ! east
!                                                          this % bMatrix (i,0)*bF(j,k,ivar,5,iel))/ &  ! west
!                                    this % qWeights (i) + &
!                                    (this % bMatrix (j,1)*bF(i,k,ivar,4,iel) + & ! north
!                                     this % bMatrix (j,0)*bF(i,k,ivar,2,iel))/ &  ! south
!                                    this % qWeights (j) + &
!                                    (this % bMatrix (k,1)*bF(i,j,ivar,6,iel) + & ! top
!                                     this % bMatrix (k,0)*bF(i,j,ivar,1,iel))/ &  ! bottom
!                                    this % qWeights (k)

!             end do
!           end do
!         end do
!       end do
!     end do

!   end subroutine VectorDGDivergence_3D_cpu

!   subroutine VectorDGDivergence_3D_gpu(this,f_dev,bF_dev,dF_dev,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)         :: nvars,nelems
!     type(c_ptr),intent(in)     :: f_dev
!     type(c_ptr),intent(in)     :: bF_dev
!     type(c_ptr),intent(out)    :: dF_dev

!     call VectorDGDivergence_3D_gpu_wrapper(this % dgMatrix % deviceData, &
!                                            this % bMatrix % deviceData, &
!                                            this % qWeights % deviceData, &
!                                            f_dev,bF_dev,dF_dev,this % N, &
!                                            nvars,nelems)

!   end subroutine VectorDGDivergence_3D_gpu

!   ! SUBROUTINE VectorCurl_3D_cpu(this,f,dF,nvars,nelems)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)     :: nvars,nelems
!   !   REAL(prec),INTENT(in)  :: f(1:3,0:this % N,0:this % N,0:this % N,1:nelems,1:nvars)
!   !   REAL(prec),INTENT(out) :: dF(1:3,0:this % N,0:this % N,0:this % N,1:nelems,1:nvars)
!   !   ! Local
!   !   INTEGER    :: i,j,k,ii,iel,ivar

!   !   DO iel = 1,nelems
!   !     DO ivar = 1,nvars
!   !       DO k = 0,this % N
!   !         DO j = 0,this % N
!   !           DO i = 0,this % N

!   !             dF(1,i,j,k,iel,ivar) = 0.0_prec
!   !             dF(2,i,j,k,iel,ivar) = 0.0_prec
!   !             dF(3,i,j,k,iel,ivar) = 0.0_prec
!   !             DO ii = 0,this % N
!   !               dF(1,i,j,k,iel,ivar) = dF(1,i,j,k,iel,ivar) + this % dMatrix (ii,j)*f(3,i,ii,k,iel,ivar) - &
!   !                                      this % dMatrix (ii,k)*f(2,i,j,ii,iel,ivar)
!   !               dF(2,i,j,k,iel,ivar) = dF(2,i,j,k,iel,ivar) + this % dMatrix (ii,k)*f(1,i,j,ii,iel,ivar) - &
!   !                                      this % dMatrix (ii,i)*f(3,ii,j,k,iel,ivar)
!   !               dF(3,i,j,k,iel,ivar) = dF(3,i,j,k,iel,ivar) + this % dMatrix (ii,i)*f(2,ii,j,k,iel,ivar) - &
!   !                                      this % dMatrix (ii,j)*f(1,i,ii,k,iel,ivar)
!   !             END DO

!   !           END DO
!   !         END DO
!   !       END DO
!   !     END DO
!   !   END DO

!   ! END SUBROUTINE VectorCurl_3D_cpu

!   ! SUBROUTINE VectorCurl_3D_gpu(this,f_dev,dF_dev,nvars,nelems)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)         :: nvars,nelems
!   !   TYPE(c_ptr),INTENT(in)     :: f_dev
!   !   TYPE(c_ptr),INTENT(out)    :: dF_dev

!   !   CALL VectorCurl_3D_gpu_wrapper(this % dMatrix % deviceData, &
!   !                                  f_dev,dF_dev,this % N, &
!   !                                  nvars,nelems)

!   ! END SUBROUTINE VectorCurl_3D_gpu

!   ! SUBROUTINE TensorDivergence_3D_cpu(this,f,dF,nvars,nelems)
!   !   ! Note that the divergence is taken over the first dimension (row dimension) of the tensor matrix
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)     :: nvars,nelems
!   !   REAL(prec),INTENT(in)  :: f(1:3,1:3,0:this % N,0:this % N,0:this % N,1:nelems,1:nvars)
!   !   REAL(prec),INTENT(out) :: dF(1:3,0:this % N,0:this % N,0:this % N,1:nelems,1:nvars)
!   !   ! Local
!   !   INTEGER    :: i,j,k,ii,iel,ivar

!   !   DO iel = 1,nelems
!   !     DO ivar = 1,nvars
!   !       DO k = 0,this % N
!   !         DO j = 0,this % N
!   !           DO i = 0,this % N

!   !             dF(1,i,j,k,iel,ivar) = 0.0_prec
!   !             dF(2,i,j,k,iel,ivar) = 0.0_prec
!   !             dF(3,i,j,k,iel,ivar) = 0.0_prec
!   !             DO ii = 0,this % N
!   !               dF(1,i,j,k,iel,ivar) = dF(1,i,j,k,iel,ivar) + &
!   !                                      this % dMatrix (ii,i)*f(1,1,ii,j,k,iel,ivar) + &
!   !                                      this % dMatrix (ii,j)*f(2,1,i,ii,k,iel,ivar) + &
!   !                                      this % dMatrix (ii,k)*f(3,1,i,j,ii,iel,ivar)

!   !               dF(2,i,j,k,iel,ivar) = dF(2,i,j,k,iel,ivar) + &
!   !                                      this % dMatrix (ii,i)*f(1,2,ii,j,k,iel,ivar) + &
!   !                                      this % dMatrix (ii,j)*f(2,2,i,ii,k,iel,ivar) + &
!   !                                      this % dMatrix (ii,k)*f(3,2,i,j,ii,iel,ivar)

!   !               dF(3,i,j,k,iel,ivar) = dF(3,i,j,k,iel,ivar) + &
!   !                                      this % dMatrix (ii,i)*f(1,3,ii,j,k,iel,ivar) + &
!   !                                      this % dMatrix (ii,j)*f(2,3,i,ii,k,iel,ivar) + &
!   !                                      this % dMatrix (ii,k)*f(3,3,i,j,ii,iel,ivar)
!   !             END DO

!   !           END DO
!   !         END DO
!   !       END DO
!   !     END DO
!   !   END DO

!   ! END SUBROUTINE TensorDivergence_3D_cpu

!   ! SUBROUTINE TensorDivergence_3D_gpu(this,f_dev,dF_dev,nvars,nelems)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)         :: nvars,nelems
!   !   TYPE(c_ptr),INTENT(in)     :: f_dev
!   !   TYPE(c_ptr),INTENT(out)    :: dF_dev

!   !   CALL TensorDivergence_3D_gpu_wrapper(this % dMatrix % deviceData, &
!   !                                        f_dev,dF_dev,this % N, &
!   !                                        nvars,nelems)

!   ! END SUBROUTINE TensorDivergence_3D_gpu

!   ! SUBROUTINE TensorDGDivergence_3D_cpu(this,f,bF,dF,nvars,nelems)
!   !   ! Note that the divergence is taken over the first dimension (row dimension) of the tensor matrix
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)     :: nvars,nelems
!   !   REAL(prec),INTENT(in)  :: f(1:3,1:3,0:this % N,0:this % N,0:this % N,1:nelems,1:nvars)
!   !   REAL(prec),INTENT(in)  :: bF(1:3,1:3,0:this % N,0:this % N,1:nvars,1:6,1:nelems)
!   !   REAL(prec),INTENT(out) :: dF(1:3,0:this % N,0:this % N,0:this % N,1:nelems,1:nvars)
!   !   ! Local
!   !   INTEGER    :: i,j,k,ii,iel,ivar

!   !   DO iel = 1,nelems
!   !     DO ivar = 1,nvars
!   !       DO k = 0,this % N
!   !         DO j = 0,this % N
!   !           DO i = 0,this % N

!   !             dF(1,i,j,k,iel,ivar) = 0.0_prec
!   !             dF(2,i,j,k,iel,ivar) = 0.0_prec
!   !             dF(3,i,j,k,iel,ivar) = 0.0_prec
!   !             DO ii = 0,this % N
!   !               dF(1,i,j,k,iel,ivar) = dF(1,i,j,k,iel,ivar) + &
!   !                                      this % dgMatrix (ii,i)*f(1,1,ii,j,k,iel,ivar) + &
!   !                                      this % dgMatrix (ii,j)*f(2,1,i,ii,k,iel,ivar) + &
!   !                                      this % dgMatrix (ii,k)*f(3,1,i,j,ii,iel,ivar)

!   !               dF(2,i,j,k,iel,ivar) = dF(2,i,j,k,iel,ivar) + &
!   !                                      this % dgMatrix (ii,i)*f(1,2,ii,j,k,iel,ivar) + &
!   !                                      this % dgMatrix (ii,j)*f(2,2,i,ii,k,iel,ivar) + &
!   !                                      this % dgMatrix (ii,k)*f(3,2,i,j,ii,iel,ivar)

!   !               dF(3,i,j,k,iel,ivar) = dF(3,i,j,k,iel,ivar) + &
!   !                                      this % dgMatrix (ii,i)*f(1,3,ii,j,k,iel,ivar) + &
!   !                                      this % dgMatrix (ii,j)*f(2,3,i,ii,k,iel,ivar) + &
!   !                                      this % dgMatrix (ii,k)*f(3,3,i,j,ii,iel,ivar)
!   !             END DO

!   !             dF(1,i,j,k,iel,ivar) = dF(1,i,j,k,iel,ivar) + (this % bMatrix (i,1)*bF(1,1,j,k,ivar,3,iel) + & ! east
!   !                                                    this % bMatrix (i,0)*bF(1,1,j,k,ivar,5,iel))/ &  ! west
!   !                                    this % qWeights (i) + &
!   !                                    (this % bMatrix (j,1)*bF(2,1,i,k,ivar,4,iel) + & ! north
!   !                                     this % bMatrix (j,0)*bF(2,1,i,k,ivar,2,iel))/ &  ! south
!   !                                    this % qWeights (j) + &
!   !                                    (this % bMatrix (k,1)*bF(3,1,i,j,ivar,6,iel) + & ! top
!   !                                     this % bMatrix (k,0)*bF(3,1,i,j,ivar,1,iel))/ &  ! bottom
!   !                                    this % qWeights (k)

!   !             dF(2,i,j,k,iel,ivar) = dF(2,i,j,k,iel,ivar) + (this % bMatrix (i,1)*bF(1,2,j,k,ivar,3,iel) + & ! east
!   !                                                    this % bMatrix (i,0)*bF(1,2,j,k,ivar,5,iel))/ &  ! west
!   !                                    this % qWeights (i) + &
!   !                                    (this % bMatrix (j,1)*bF(2,2,i,k,ivar,4,iel) + & ! north
!   !                                     this % bMatrix (j,0)*bF(2,2,i,k,ivar,2,iel))/ &  ! south
!   !                                    this % qWeights (j) + &
!   !                                    (this % bMatrix (k,1)*bF(3,2,i,j,ivar,6,iel) + & ! top
!   !                                     this % bMatrix (k,0)*bF(3,2,i,j,ivar,1,iel))/ &  ! bottom
!   !                                    this % qWeights (k)

!   !             dF(3,i,j,k,iel,ivar) = dF(3,i,j,k,iel,ivar) + (this % bMatrix (i,1)*bF(1,3,j,k,ivar,3,iel) + & ! east
!   !                                                    this % bMatrix (i,0)*bF(1,3,j,k,ivar,5,iel))/ &  ! west
!   !                                    this % qWeights (i) + &
!   !                                    (this % bMatrix (j,1)*bF(2,3,i,k,ivar,4,iel) + & ! north
!   !                                     this % bMatrix (j,0)*bF(2,3,i,k,ivar,2,iel))/ &  ! south
!   !                                    this % qWeights (j) + &
!   !                                    (this % bMatrix (k,1)*bF(3,3,i,j,ivar,6,iel) + & ! top
!   !                                     this % bMatrix (k,0)*bF(3,3,i,j,ivar,1,iel))/ &  ! bottom
!   !                                    this % qWeights (k)

!   !           END DO
!   !         END DO
!   !       END DO
!   !     END DO
!   !   END DO

!   ! END SUBROUTINE TensorDGDivergence_3D_cpu

!   ! SUBROUTINE TensorDGDivergence_3D_gpu(this,f_dev,bF_dev,dF_dev,nvars,nelems)
!   !   IMPLICIT NONE
!   !   CLASS(Lagrange),INTENT(in) :: this
!   !   INTEGER,INTENT(in)         :: nvars,nelems
!   !   TYPE(c_ptr),INTENT(in)     :: f_dev
!   !   TYPE(c_ptr),INTENT(in)     :: bF_dev
!   !   TYPE(c_ptr),INTENT(out)    :: dF_dev

!   !   CALL TensorDGDivergence_3D_gpu_wrapper(this % dgMatrix % deviceData, &
!   !                                          this % bMatrix % deviceData, &
!   !                                          this % qWeights % deviceData, &
!   !                                          f_dev,bF_dev,dF_dev,this % N, &
!   !                                          nvars,nelems)

!   ! END SUBROUTINE TensorDGDivergence_3D_gpu
!   ! /////////////////////////////// !
!   ! Boundary Interpolation Routines !

  subroutine ScalarBoundaryInterp_1D_cpu(this,f,fTarget,nvars,nelems)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)         :: nvars,nelems
    real(prec),intent(in)      :: f(1:this % N+1,1:nelems,1:nvars)
    real(prec),intent(inout)   :: fTarget(1:2,1:nelems,1:nvars)
    ! Local
    integer :: ii,iel,ivar
    real(prec) :: fb(1:2)

    do iel = 1,nelems
      do ivar = 1,nvars
        fb(1:2) = 0.0_prec
        do ii = 1,this % N+1
          fb(1) = fb(1) + this % bMatrix (ii,1)*f(ii,iel,ivar) ! West
          fb(2) = fb(2) + this % bMatrix (ii,2)*f(ii,iel,ivar) ! East
        end do
        fTarget(1:2,iel,ivar) = fb(1:2)
      end do
    end do

  end subroutine ScalarBoundaryInterp_1D_cpu

  subroutine ScalarBoundaryInterp_1D_gpu(this,f,fTarget,nvars,nelems,hipblas_handle)
    implicit none
    class(Lagrange),intent(in) :: this
    integer,intent(in)  :: nvars,nelems
    real(prec),pointer,intent(in)      :: f(:,:,:)
    real(prec),pointer,intent(inout)     :: fTarget(:,:,:)
    type(c_ptr),intent(inout) :: hipblas_handle
    ! Local
    integer(c_int) :: m
    integer(c_int) :: n
    integer(c_int) :: k
    real(c_prec) :: alpha
    integer(c_int) :: lda
    integer(c_int) :: ldb
    integer(c_int) :: ldc
    real(c_prec) :: beta

    m = 2 ! number of rows of A^T
    n = nvars*nelems ! number of columns of B
    k = this % N + 1! number of columns of A^T
    alpha = 1.0_c_prec
    lda = k ! leading dimension of A (interoplation matrix)
    ldb = k ! leading dimension of B (f)
    ldc = m ! leading dimension of C (fTarget)
    beta = 0.0_c_prec

#ifdef DOUBLE_PRECISION
    call hipblasCheck(hipblasDgemm(hipblas_handle, &
                                   HIPBLAS_OP_T,HIPBLAS_OP_N, &
                                   m,n,k,alpha, &
                                   c_loc(this % bMatrix),lda, &
                                   c_loc(f),ldb, &
                                   beta, &
                                   c_loc(fTarget),ldc))
#else
    call hipblasCheck(hipblasSgemm(hipblas_handle,&
        HIPBLAS_OP_T, HIPBLAS_OP_N, &
        m, n, k, alpha, &
        c_loc(this % bMatrix), lda, &
        c_loc(f), ldb, &
        beta, &
        c_loc(fTarget), ldc))
#endif

  end subroutine ScalarBoundaryInterp_1D_gpu

!   subroutine ScalarBoundaryInterp_2D_cpu(this,f,fTarget,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)         :: nvars,nelems
!     real(prec),intent(in)      :: f(0:this % N,0:this % N,1:nelems,1:nvars)
!     real(prec),intent(out)     :: fTarget(0:this % N,1:nvars,1:4,1:nelems)
!     ! Local
!     integer :: i,ii,iel,ivar
!     real(prec) :: fb(1:4)

!     do iel = 1,nelems
!       do ivar = 1,nvars
!         do i = 0,this % N

!           fb(1:4) = 0.0_prec

!           do ii = 0,this % N
!             fb(1) = fb(1) + this % bMatrix (ii,0)*f(i,ii,iel,ivar) ! South
!             fb(2) = fb(2) + this % bMatrix (ii,1)*f(ii,i,iel,ivar) ! East
!             fb(3) = fb(3) + this % bMatrix (ii,1)*f(i,ii,iel,ivar) ! North
!             fb(4) = fb(4) + this % bMatrix (ii,0)*f(ii,i,iel,ivar) ! West
!           end do

!           fTarget(i,ivar,1:4,iel) = fb(1:4)

!         end do
!       end do
!     end do

!   end subroutine ScalarBoundaryInterp_2D_cpu

!   subroutine ScalarBoundaryInterp_2D_gpu(this,f,fTarget,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)  :: nvars,nelems
!     type(c_ptr),intent(in)  :: f
!     type(c_ptr),intent(out)  :: fTarget

!     call ScalarBoundaryInterp_2D_gpu_wrapper(this % bMatrix % deviceData, &
!                                              f,fTarget,this % N,nvars,nelems)

!   end subroutine ScalarBoundaryInterp_2D_gpu

!   subroutine VectorBoundaryInterp_2D_cpu(this,f,fTarget,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)  :: nvars,nelems
!     real(prec),intent(in)  :: f(1:2,0:this % N,0:this % N,1:nelems,1:nvars)
!     real(prec),intent(out)  :: fTarget(1:2,0:this % N,1:nvars,1:4,1:nelems)
!     ! Local
!     integer :: i,ii,idir,iel,ivar
!     real(prec) :: fb(1:2,1:4)

!     do iel = 1,nelems
!       do ivar = 1,nvars
!         do i = 0,this % N

!           fb(1:2,1:4) = 0.0_prec
!           do ii = 0,this % N
!             do idir = 1,2
!               fb(idir,1) = fb(idir,1) + this % bMatrix (ii,0)*f(idir,i,ii,iel,ivar) ! South
!               fb(idir,2) = fb(idir,2) + this % bMatrix (ii,1)*f(idir,ii,i,iel,ivar) ! East
!               fb(idir,3) = fb(idir,3) + this % bMatrix (ii,1)*f(idir,i,ii,iel,ivar) ! North
!               fb(idir,4) = fb(idir,4) + this % bMatrix (ii,0)*f(idir,ii,i,iel,ivar) ! West
!             end do
!           end do

!           do idir = 1,2
!             fTarget(idir,i,ivar,1:4,iel) = fb(idir,1:4)
!           end do

!         end do
!       end do
!     end do

!   end subroutine VectorBoundaryInterp_2D_cpu

!   subroutine VectorBoundaryInterp_2D_gpu(this,f,fTarget,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)  :: nvars,nelems
!     type(c_ptr),intent(in)  :: f
!     type(c_ptr),intent(out)  :: fTarget

!     call VectorBoundaryInterp_2D_gpu_wrapper(this % bMatrix % deviceData, &
!                                              f,fTarget,this % N,nvars,nelems)

!   end subroutine VectorBoundaryInterp_2D_gpu

!   subroutine TensorBoundaryInterp_2D_cpu(this,f,fTarget,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)  :: nvars,nelems
!     real(prec),intent(in)  :: f(1:2,1:2,0:this % N,0:this % N,1:nelems,1:nvars)
!     real(prec),intent(out)  :: fTarget(1:2,1:2,0:this % N,1:nvars,1:4,1:nelems)
!     ! Local
!     integer :: i,ii,idir,jdir,iel,ivar
!     real(prec) :: fb(1:2,1:2,1:4)

!     do iel = 1,nelems
!       do ivar = 1,nvars
!         do i = 0,this % N

!           fb(1:2,1:2,1:4) = 0.0_prec
!           do ii = 0,this % N
!             do jdir = 1,2
!               do idir = 1,2
!                 fb(idir,jdir,1) = fb(idir,jdir,1) + this % bMatrix (ii,0)*f(idir,jdir,i,ii,iel,ivar) ! South
!                 fb(idir,jdir,2) = fb(idir,jdir,2) + this % bMatrix (ii,1)*f(idir,jdir,ii,i,iel,ivar) ! East
!                 fb(idir,jdir,3) = fb(idir,jdir,3) + this % bMatrix (ii,1)*f(idir,jdir,i,ii,iel,ivar) ! North
!                 fb(idir,jdir,4) = fb(idir,jdir,4) + this % bMatrix (ii,0)*f(idir,jdir,ii,i,iel,ivar) ! West
!               end do
!             end do
!           end do

!           do jdir = 1,2
!             do idir = 1,2
!               fTarget(idir,jdir,i,ivar,1:4,iel) = fb(idir,jdir,1:4)
!             end do
!           end do

!         end do
!       end do
!     end do

!   end subroutine TensorBoundaryInterp_2D_cpu

!   subroutine TensorBoundaryInterp_2D_gpu(this,f,fTarget,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)  :: nvars,nelems
!     type(c_ptr),intent(in)  :: f
!     type(c_ptr),intent(out)  :: fTarget

!     call TensorBoundaryInterp_2D_gpu_wrapper(this % bMatrix % deviceData, &
!                                              f,fTarget,this % N,nvars,nelems)

!   end subroutine TensorBoundaryInterp_2D_gpu

!   subroutine ScalarBoundaryInterp_3D_cpu(this,f,fTarget,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)         :: nvars,nelems
!     real(prec),intent(in)      :: f(0:this % N,0:this % N,0:this % N,1:nelems,1:nvars)
!     real(prec),intent(out)     :: fTarget(0:this % N,0:this % N,1:nvars,1:6,1:nelems)
!     ! Local
!     integer :: i,j,ii,iel,ivar
!     real(prec) :: fb(1:6)

!     do iel = 1,nelems
!       do ivar = 1,nvars
!         do j = 0,this % N
!           do i = 0,this % N

!             fb(1:6) = 0.0_prec

!             do ii = 0,this % N
!               fb(1) = fb(1) + this % bMatrix (ii,0)*f(i,j,ii,iel,ivar) ! Bottom
!               fb(2) = fb(2) + this % bMatrix (ii,0)*f(i,ii,j,iel,ivar) ! South
!               fb(3) = fb(3) + this % bMatrix (ii,1)*f(ii,i,j,iel,ivar) ! East
!               fb(4) = fb(4) + this % bMatrix (ii,1)*f(i,ii,j,iel,ivar) ! North
!               fb(5) = fb(5) + this % bMatrix (ii,0)*f(ii,i,j,iel,ivar) ! West
!               fb(6) = fb(6) + this % bMatrix (ii,1)*f(i,j,ii,iel,ivar) ! Top
!             end do

!             fTarget(i,j,ivar,1:6,iel) = fb(1:6)

!           end do
!         end do
!       end do
!     end do

!   end subroutine ScalarBoundaryInterp_3D_cpu

!   subroutine ScalarBoundaryInterp_3D_gpu(this,f,fTarget,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)  :: nvars,nelems
!     type(c_ptr),intent(in)  :: f
!     type(c_ptr),intent(out)  :: fTarget

!     call ScalarBoundaryInterp_3D_gpu_wrapper(this % bMatrix % deviceData, &
!                                              f,fTarget,this % N,nvars,nelems)

!   end subroutine ScalarBoundaryInterp_3D_gpu

!   subroutine VectorBoundaryInterp_3D_cpu(this,f,fTarget,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)  :: nvars,nelems
!     real(prec),intent(in)  :: f(1:3,0:this % N,0:this % N,0:this % N,1:nelems,1:nvars)
!     real(prec),intent(out)  :: fTarget(1:3,0:this % N,0:this % N,1:nvars,1:6,1:nelems)
!     ! Local
!     integer :: i,j,ii,idir,iel,ivar
!     real(prec) :: fb(1:3,1:6)

!     do iel = 1,nelems
!       do ivar = 1,nvars
!         do j = 0,this % N
!           do i = 0,this % N

!             fb(1:3,1:6) = 0.0_prec
!             do ii = 0,this % N
!               do idir = 1,3
!                 fb(idir,1) = fb(idir,1) + this % bMatrix (ii,0)*f(idir,i,j,ii,iel,ivar) ! Bottom
!                 fb(idir,2) = fb(idir,2) + this % bMatrix (ii,0)*f(idir,i,ii,j,iel,ivar) ! South
!                 fb(idir,3) = fb(idir,3) + this % bMatrix (ii,1)*f(idir,ii,i,j,iel,ivar) ! East
!                 fb(idir,4) = fb(idir,4) + this % bMatrix (ii,1)*f(idir,i,ii,j,iel,ivar) ! North
!                 fb(idir,5) = fb(idir,5) + this % bMatrix (ii,0)*f(idir,ii,i,j,iel,ivar) ! West
!                 fb(idir,6) = fb(idir,6) + this % bMatrix (ii,1)*f(idir,i,j,ii,iel,ivar) ! Top
!               end do
!             end do

!             do idir = 1,3
!               fTarget(idir,i,j,ivar,1:6,iel) = fb(idir,1:6)
!             end do

!           end do
!         end do
!       end do
!     end do

!   end subroutine VectorBoundaryInterp_3D_cpu

!   subroutine VectorBoundaryInterp_3D_gpu(this,f,fTarget,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)  :: nvars,nelems
!     type(c_ptr),intent(in)  :: f
!     type(c_ptr),intent(out)  :: fTarget

!     call VectorBoundaryInterp_3D_gpu_wrapper(this % bMatrix % deviceData, &
!                                              f,fTarget,this % N,nvars,nelems)

!   end subroutine VectorBoundaryInterp_3D_gpu

!   subroutine TensorBoundaryInterp_3D_cpu(this,f,fTarget,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)  :: nvars,nelems
!     real(prec),intent(in)  :: f(1:3,1:3,0:this % N,0:this % N,0:this % N,1:nelems,1:nvars)
!     real(prec),intent(out)  :: fTarget(1:3,1:3,0:this % N,0:this % N,1:nvars,1:6,1:nelems)
!     ! Local
!     integer :: i,j,ii,idir,jdir,iel,ivar
!     real(prec) :: fb(1:3,1:3,1:6)

!     do iel = 1,nelems
!       do ivar = 1,nvars
!         do j = 0,this % N
!           do i = 0,this % N

!             fb(1:3,1:3,1:6) = 0.0_prec
!             do ii = 0,this % N
!               do jdir = 1,3
!                 do idir = 1,3
!                   fb(idir,jdir,1) = fb(idir,jdir,1) + this % bMatrix (ii,0)*f(idir,jdir,i,j,ii,iel,ivar) ! Bottom
!                   fb(idir,jdir,2) = fb(idir,jdir,2) + this % bMatrix (ii,0)*f(idir,jdir,i,ii,j,iel,ivar) ! South
!                   fb(idir,jdir,3) = fb(idir,jdir,3) + this % bMatrix (ii,1)*f(idir,jdir,ii,i,j,iel,ivar) ! East
!                   fb(idir,jdir,4) = fb(idir,jdir,4) + this % bMatrix (ii,1)*f(idir,jdir,i,ii,j,iel,ivar) ! North
!                   fb(idir,jdir,5) = fb(idir,jdir,5) + this % bMatrix (ii,0)*f(idir,jdir,ii,i,j,iel,ivar) ! West
!                   fb(idir,jdir,6) = fb(idir,jdir,6) + this % bMatrix (ii,1)*f(idir,jdir,i,j,ii,iel,ivar) ! Top
!                 end do
!               end do
!             end do

!             do jdir = 1,3
!               do idir = 1,3
!                 fTarget(idir,jdir,i,j,ivar,1:6,iel) = fb(idir,jdir,1:6)
!               end do
!             end do

!           end do
!         end do
!       end do
!     end do

!   end subroutine TensorBoundaryInterp_3D_cpu

!   subroutine TensorBoundaryInterp_3D_gpu(this,f,fTarget,nvars,nelems)
!     implicit none
!     class(Lagrange),intent(in) :: this
!     integer,intent(in)  :: nvars,nelems
!     type(c_ptr),intent(in)  :: f
!     type(c_ptr),intent(out)  :: fTarget

!     call TensorBoundaryInterp_3D_gpu_wrapper(this % bMatrix % deviceData, &
!                                              f,fTarget,this % N,nvars,nelems)

!   end subroutine TensorBoundaryInterp_3D_gpu

! ================================================================================================ !
!
! CalculateBarycentricWeights (PRIVATE)
!
!   A PRIVATE routine that calculates and stores the barycentric weights for the Lagrange
!   data-structure.
!
!   This routine is from Alg. 30 on pg. 75 of D.A. Kopriva, 2009.
!
! ================================================================================================ !

  subroutine CalculateBarycentricWeights(this)
    implicit none
    class(Lagrange),intent(inout) :: this
    ! Local
    integer :: i,j
    real(real64) :: bWeights(0:this % N)
    real(real64) :: controlPoints(0:this % N)

    do i = 0,this % N
      bWeights(i) = 1.0_real64
      controlPoints(i) = real(this % controlPoints(i + 1),real64)
    end do

    ! Computes the product w_k = w_k*(s_k - s_j), k /= j
    do j = 1,this % N
      do i = 0,j - 1

        bWeights(i) = bWeights(i)*(controlPoints(i) - controlPoints(j))
        bWeights(j) = bWeights(j)*(controlPoints(j) - controlPoints(i))

      end do
    end do

    do j = 0,this % N
      bWeights(j) = 1.0_prec/bWeights(j)
      this % bWeights(j + 1) = real(bWeights(j),prec)
    end do

  end subroutine CalculateBarycentricWeights

! ================================================================================================ !
!
! CalculateInterpolationMatrix (PRIVATE)
!
!   A PRIVATE routine that fills in the interpolation matrix for the Lagrange data structure.
!
!   This function is from Alg. 32 on pg. 76 of D.A. Kopriva, 2009.
!
! ================================================================================================ !

  subroutine CalculateInterpolationMatrix(this)
    implicit none
    class(Lagrange),intent(inout) :: this
    ! Local
    integer    :: row,col
    logical    :: rowHasMatch
    real(real64) :: temp1,temp2
    real(real64) :: iMatrix(0:this % M,0:this % N)
    real(real64) :: bWeights(0:this % N)
    real(real64) :: controlPoints(0:this % N)
    real(real64) :: targetPoints(0:this % M)

    do col = 0,this % N
      controlPoints(col) = real(this % controlPoints(col + 1),real64)
      bWeights(col) = real(this % bWeights(col + 1),real64)
    end do
    do row = 0,this % M
      targetPoints(row) = real(this % targetPoints(row + 1),real64)
    end do

    do row = 0,this % M

      rowHasMatch = .false.

      do col = 0,this % N

        iMatrix(row,col) = 0.0_real64

        if (AlmostEqual(targetPoints(row),controlPoints(col))) then
          rowHasMatch = .true.
          iMatrix(row,col) = 1.0_real64
        end if

      end do

      if (.not. (rowHasMatch)) then

        temp1 = 0.0_real64

        do col = 0,this % N
          temp2 = bWeights(col)/ &
                  (targetPoints(row) - &
                   controlPoints(col))
          iMatrix(row,col) = temp2
          temp1 = temp1 + temp2
        end do

        do col = 0,this % N
          iMatrix(row,col) = iMatrix(row,col)/temp1
        end do

      end if

    end do

    do row = 0,this % M
      do col = 0,this % N
        this % iMatrix(col + 1,row + 1) = real(iMatrix(row,col),prec)
      end do
    end do

  end subroutine CalculateInterpolationMatrix

! ================================================================================================ !
!
! CalculateDerivativeMatrix (PRIVATE)
!
!   Calculates and stores the derivative matrix and its transpose.
!   Generates a matrix that can be used to approximate derivatives at the interpolation nodes.
!
!   This function is from Alg. 37 on pg. 82 of D.A. Kopriva, 2009.
!
! ================================================================================================ !

  subroutine CalculateDerivativeMatrix(this)
    implicit none
    class(Lagrange),intent(inout) :: this
    ! Local
    integer      :: row,col
    real(real64) :: dmat(0:this % N,0:this % N)
    real(real64) :: dgmat(0:this % N,0:this % N)
    real(real64) :: bWeights(0:this % N)
    real(real64) :: qWeights(0:this % N)
    real(real64) :: controlPoints(0:this % N)

    do row = 0,this % N
      bWeights(row) = real(this % bWeights(row + 1),real64)
      qWeights(row) = real(this % qWeights(row + 1),real64)
      controlPoints(row) = real(this % controlPoints(row + 1),real64)
    end do

    do row = 0,this % N

      dmat(row,row) = 0.0_prec

      do col = 0,this % N

        if (.not. (col == row)) then

          dmat(row,col) = bWeights(col)/ &
                          (bWeights(row)* &
                           (controlPoints(row) - &
                            controlPoints(col)))

          dmat(row,row) = dmat(row,row) - dmat(row,col)

        end if

      end do

    end do

    do row = 0,this % N
      do col = 0,this % N
        dgmat(row,col) = -dmat(col,row)* &
                         qWeights(col)/ &
                         qWeights(row)
      end do
    end do

    do row = 0,this % N
      do col = 0,this % N
        this % dMatrix(row + 1,col + 1) = real(dmat(col,row),prec)
        this % dgMatrix(row + 1,col + 1) = real(dgmat(col,row),prec)
      end do
    end do

  end subroutine CalculateDerivativeMatrix

! ================================================================================================ !
!
! CalculateLagrangePolynomials
!
!   Evaluates each of the 1-D Lagrange interpolating polynomials at a specified point.
!
!   This function is from Alg. 34 on pg. 77 of D.A. Kopriva, 2009.
!
! ================================================================================================ !

  function CalculateLagrangePolynomials(this,sE) result(lAtS)
    implicit none
    class(Lagrange) :: this
    real(prec)      :: sE
    real(prec)      :: lAtS(0:this % N)
    ! Local
    integer    :: j
    logical    :: xMatchesNode
    real(real64) :: temp1,temp2
    real(real64) :: sELocal
    real(real64) :: controlPoints(0:this % N)
    real(real64) :: bWeights(0:this % N)
    real(real64) :: lS(0:this % N)

    sELocal = real(sE,real64)
    do j = 0,this % N
      controlPoints(j) = real(this % controlPoints(j + 1),real64)
      bWeights(j) = real(this % bWeights(j + 1),real64)
    end do

    xMatchesNode = .false.

    do j = 0,this % N

      lS(j) = 0.0_real64
      if (AlmostEqual(sELocal,controlPoints(j))) then
        lS(j) = 1.0_real64
        xMatchesNode = .true.
      end if

    end do

    if (xMatchesNode) then
      do j = 0,this % N
        lAtS(j) = real(lS(j),prec)
      end do
      return
    end if

    temp1 = 0.0_real64

    do j = 0,this % N
      temp2 = bWeights(j)/(sE - controlPoints(j))
      lS(j) = temp2
      temp1 = temp1 + temp2
    end do

    lS = lS/temp1

    do j = 0,this % N
      lAtS(j) = real(lS(j),prec)
    end do

  end function CalculateLagrangePolynomials

  ! subroutine WriteHDF5_Lagrange(this,fileId)
  !   implicit none
  !   class(Lagrange),intent(in) :: this
  !   integer(HID_T),intent(in) :: fileId

  !   call CreateGroup_HDF5(fileId,'/interp')

  !   call WriteArray_HDF5(fileId,'/interp/controlpoints', &
  !                        this % controlPoints)

  !   call WriteArray_HDF5(fileId,'/interp/qweights', &
  !                        this % qWeights)

  !   call WriteArray_HDF5(fileId,'/interp/dgmatrix', &
  !                        this % dgMatrix)

  !   call WriteArray_HDF5(fileId,'/interp/dmatrix', &
  !                        this % dMatrix)

  !   call WriteArray_HDF5(fileId,'/interp/bmatrix', &
  !                        this % bMatrix)

  !   call WriteArray_HDF5(fileId,'/interp/imatrix', &
  !                        this % iMatrix)

  ! end subroutine WriteHDF5_Lagrange

end module SELF_Lagrange
