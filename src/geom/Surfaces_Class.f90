! Surfaces_CLASS.f90
!
! Copyright 2017 Joseph Schoonover <joe@fluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !



MODULE Surfaces_CLASS

! src/COMMON/
  USE ModelPrecision
  USE ConstantsDictionary
  USE Lagrange_CLASS

  IMPLICIT NONE

!  The Surfaces CLASS provides attributes and TYPE-bound procedures for defining and manipulating
!  surfaces in multiple dimensions.
!
!  A surface is a geometric primitive that can be described two free parameters.
!
!  As an example, a surface in three-dimensions is represented as
!
!       \vec{x}(\xi^1,\xi^2) = x(\xi^1,\xi^2) \hat{x} + y(\xi^1,\xi^2) \hat{y} + z(\xi^1,\xi^2) \hat{z}
!
!  where (\xi^1,\xi^2)  are parameters defined on [-1,1]x[-1,1].
!  In the SELF, surfaces in 3-D are primarily USEd in the generation of mappings between physical
!  and computational space for hexahedral elements.
!
!  The Surfaces CLASS permits surfaces that reside in higher dimensions.
!
!

  TYPE Surfaces
    INTEGER                 :: N, nSurfaces
    TYPE(Lagrange)          :: interp
    REAL(prec), ALLOCATABLE :: x(:,:,:,:)
    REAL(prec), ALLOCATABLE :: dxds(:,:,:,:,:)

#ifdef HAVE_CUDA
    INTEGER, DEVICE, ALLOCATABLE    :: N_dev, nSurfaces_dev
    REAL(prec), DEVICE, ALLOCATABLE :: x_dev(:,:,:,:)
    REAL(prec), DEVICE, ALLOCATABLE :: dxds_dev(:,:,:,:,:)
#endif

  CONTAINS

    PROCEDURE :: Build => Build_Surfaces
    PROCEDURE :: Trash => Trash_Surfaces

#ifdef HAVE_CUDA
    PROCEDURE :: UpdateHost => UpdateHost_Surfaces
    PROCEDURE :: UpdateDevice => UpdateDevice_Surfaces
#endif

    PROCEDURE :: Set_Surfaces
    PROCEDURE :: CalculateSlope => CalculateSlope_Surfaces
    PROCEDURE :: Evaluate       => Evaluate_Surfaces
    PROCEDURE :: EvaluateSlope  => EvaluateSlope_Surfaces

  END TYPE Surfaces

#ifdef HAVE_CUDA
  INTEGER, DEVICE, ALLOCATABLE, PRIVATE :: nDim_dev
#endif

CONTAINS
!
! ================================================================================================ !
!  Build_Surfaces
! Initializes and assigns the attributes of the Surfaces CLASS.
!
!  This SUBROUTINE depends on
!   Module \ref Lagrange : S/R \ref Build_Lagrange
!   Module \ref Surfaces_CLASS : S/R \ref CalculateSlope_Surfaces
!
!  Usage : </H2>
! TYPE</B>(Surfaces) :: this
! INTEGER</B>       :: N, nDim
! REAL</B>(prec)    :: x(0:N,0:N,1:nDim), nodes(0:N)
!         ....
!     CALL</B> this % Build( x, nodes, N, nDim )
!
!   Parameters : </H2>
!
!     out <th> mySurfaces  Surfaces  On output, the surface interpolant, position vectors,
!                                              and derivatives are filled in.
!     in <th> x(0:N,0:N,1:nDim)  REAL(prec)  Surfaces position vectors
!     in <th> nodes(0:N)*  REAL(prec)  Discrete locations of the free parameter.
!     in <th> N  INTEGER  Polynomial degree of the interpolant that describes
!                                         the surface.
!     in <th> nDim  INTEGER  Number of spatial dimensions that the surface resides in.
!  </table>
!  * IF the surface is being USEd for Mapped-Geometry element construction, the nodes must be between
!    [-1,1]. A tensor product of "nodes" with itself is USEd to construct the 2-D free parameter
!    space
!
! ================================================================================================ !

  SUBROUTINE Build_Surfaces( mySurfaces, nodes, N, nSurfaces )
    IMPLICIT NONE
    CLASS( Surfaces ), INTENT(out) :: mySurfaces
    INTEGER, INTENT(in)            :: N, nSurfaces
    REAL(prec), INTENT(in)         :: nodes(0:N)


    mySurfaces % N         = N
    mySurfaces % nSurfaces = nSurfaces

    ALLOCATE( mySurfaces % x(0:N,0:N,1:3,1:nSurfaces), &
      mySurfaces % dxds(1:2,0:N,0:N,1:3,1:nSurfaces) )

    CALL mySurfaces % interp % Build( N, N, nodes, nodes )
    mySurfaces % x = 0.0_prec
    mySurfaces % dxds = 0.0_prec

#ifdef HAVE_CUDA

    ALLOCATE( mySurfaces % N_dev, mySurfaces % nSurfaces_dev )
    ALLOCATE( mySurfaces % x_dev(0:N,0:N,1:3,1:nSurfaces), &
      mySurfaces % dxds_dev(1:2,0:N,0:N,1:3,1:nSurfaces) )

    mySurfaces % N_dev         = N
    mySurfaces % nSurfaces_dev = nSurfaces

    ALLOCATE( nDim_dev )
    nDim_dev = 3

#endif


  END SUBROUTINE Build_Surfaces

! ================================================================================================ !
!  Trash_Surfaces
! Frees memory held by the attributes of the Surfaces CLASS.
!
!  This SUBROUTINE depends on
!   Module \ref Lagrange : S/R \ref Trash_Lagrange
!
!  Usage : </H2>
! TYPE</B>(Surfaces) :: this
!         ....
!     CALL</B> this % Trash( )
!
!   Parameters : </H2>
!
!     in/out <th> mySurfaces  Surfaces
!                         On input</B>, a previously constructed Surfaces DATA-structure,
!                         On output</B>, the memory held by its attributes is freed.
!  </table>
!
! ================================================================================================ !

  SUBROUTINE Trash_Surfaces( mySurfaces )
    IMPLICIT NONE
    CLASS( Surfaces ), INTENT(inout)     :: mySurfaces


    CALL mySurfaces % interp % Trash( )
    DEALLOCATE( mySurfaces % x, mySurfaces % dxds )

#ifdef HAVE_CUDA
    DEALLOCATE( mySurfaces % N_dev, mySurfaces % nSurfaces_dev )
    DEALLOCATE( mySurfaces % x_dev, mySurfaces % dxds_dev )

    DEALLOCATE( nDim_dev )
#endif


  END SUBROUTINE Trash_Surfaces

#ifdef HAVE_CUDA
  SUBROUTINE UpdateDevice_Surfaces( mySurfaces )
    IMPLICIT NONE
    CLASS( Surfaces ), INTENT(inout) :: mySurfaces


    mySurfaces % x_dev    = mySurfaces % x
    mySurfaces % dxds_dev = mySurfaces % dxds


  END SUBROUTINE UpdateDevice_Surfaces


  SUBROUTINE UpdateHost_Surfaces( mySurfaces )
    IMPLICIT NONE
    CLASS( Surfaces ), INTENT(inout) :: mySurfaces


    mySurfaces % x    = mySurfaces % x_dev
    mySurfaces % dxds = mySurfaces % dxds_dev


  END SUBROUTINE UpdateHost_Surfaces
#endif

! ================================================================================================ !
!  Set_Surfaces
!  Sets the Surfaces position DATA and re-calculates the surface slope DATA.
!
!  IF you want to re-USE an previously constructed surface, this routine overwrites the surface position
!  and slope information. Note that the number of position DATA points must not change.
!
!  This SUBROUTINE depends on
!   Module \ref Surfaces_CLASS : S/R \ref CalculateSlope_Surfaces
!
!  Usage : </H2>
! TYPE</B>(Surfaces) :: this
! REAL</B>(prec)    :: x(0:this % N,0:this % N,1:this % nDim)
!         ....
!     CALL</B> this % Reset( x )
!
!   Parameters : </H2>
!
!     inout <th> mySurfaces  Surfaces
!                         On input</B>, a previously constructed Surfaces,
!                         On output</B>, new surface positions and surface slopes are filled in
!     in <th> x(0:mySurfaces % N,0:mySurfaces % N, 1: mySurfaces % nDim)  REAL(prec)
!                     Surfaces position vectors
!  </table>
!
! ================================================================================================ !

  SUBROUTINE Set_Surfaces( mySurfaces, x )
    IMPLICIT NONE
    CLASS( Surfaces ), INTENT(inout) :: mySurfaces
    REAL(prec), INTENT(in)           :: x(0:mySurfaces % N, 0:mySurfaces % N, 1:3, 1:mySurfaces % nSurfaces)

    mySurfaces % x = x

#ifdef HAVE_CUDA

    CALL mySurfaces % UpdateDevice( )

#endif


    CALL mySurfaces % CalculateSlope( )


#ifdef HAVE_CUDA

    CALL mySurfaces % UpdateHost( )

#endif

  END SUBROUTINE Set_Surfaces

! ================================================================================================ !
!  CalculateSlope_Surfaces
! Calculates the derivative of the surface position wrt to the free parameter and stores the result
! within the DATA structure (attribute "dxds").
!
! The interpolant's derivative matrix is USEd to quickly compute the surface slope at each of the
! points where the surface position is known.
!
!  This SUBROUTINE depends on
!   Module \ref Lagrange : S/R \ref ApplyDerivativeMatrix_Lagrange
!
!  Usage : </H2>
! TYPE</B>(Surfaces) :: this
!         ....
!     CALL</B> this % CalculateSlope(  )
!
!   Parameters : </H2>
!
!     in/out <th> mySurfaces  Surfaces
!                         On input</B>, a Surfaces structure with the position ("x") attribute
!                         filled in,
!                         On output</B>, the surface slope ("dxds") attribute is filled in
!  </table>
!
! ================================================================================================ !

  SUBROUTINE CalculateSlope_Surfaces( mySurfaces )
    IMPLICIT NONE
    CLASS( Surfaces ), INTENT(inout) :: mySurfaces

#ifdef HAVE_CUDA

    CALL mySurfaces % interp % CalculateGradient_2D( mySurfaces % x_dev, mySurfaces % dxds_dev, &
      nDim_dev, mySurfaces % nSurfaces_dev )

#else

    CALL mySurfaces % interp % CalculateGradient_2D( mySurfaces % x, mySurfaces % dxds, &
      3, mySurfaces % nSurfaces )

#endif

  END SUBROUTINE CalculateSlope_Surfaces

! ================================================================================================ !
!  Evaluate_Surfaces
! Estimates the surface position at a given value of the surface parameter.
!
!  This SUBROUTINE depends on
!   Module \ref Lagrange : S/R \ref Interpolate_Lagrange
!
!  Usage : </H2>
! TYPE</B>(Surfaces) :: this
! REAL</B>(prec)    :: s(1:2)
! REAL</B>(prec)    :: x(1:this % nDim)
!         ....
!     x = this % Evaluate( s )
!
!   Parameters : </H2>
!
!     in <th> mySurfaces  Surfaces  A previously constructed Surfaces DATA structure
!     in <th> s(1:2)  REAL(prec)  Value of the surface parameters where the surface position
!                                            is desired.
!     in <th> x(1:mySurfaces % nDim)  REAL(prec)  Position of the surface at "s"
!  </table>
!
! ================================================================================================ !

  FUNCTION Evaluate_Surfaces( mySurfaces, s, j ) RESULT( x )
    IMPLICIT NONE
    CLASS( Surfaces ) :: mySurfaces
    REAL(prec)        :: s(1:2)
    REAL(prec)        :: x(1:3)
    INTEGER           :: j
    ! Local
    INTEGER :: i


    DO i = 1, 3

      x(i) = mySurfaces % interp % Interpolate_2D( mySurfaces % x(:,:,i,j), s )

    ENDDO


  END FUNCTION Evaluate_Surfaces

! ================================================================================================ !
!  EvaluateSlope_Surfaces
!
!   A surfaceEstimates the surface slope at a given value of the surface parameters.
!
!   This SUBROUTINE depends on
!    Module \ref Lagrange : S/R \ref DIFferentiate_Lagrange
!
!   Usage : </H2>
!  TYPE</B>(Surfaces) :: this
!  REAL</B>(prec)    :: s(1:2)
! REAL</B>(prec)    :: dxds(1:this % nDim,1:2)
!         ....
!     dxds = this % EvaluateSlope( s )
!
!   Parameters : </H2>
!
!     in <th> mySurfaces  Surfaces  A previously constructed Surfaces DATA structure
!     in <th> s(1:2)  REAL(prec)  Value of the surface parameters where the surface slope
!                                            is desired.
!     in <th> dxds(1:mySurfaces % nDim,1:2)  REAL(prec)  Slope of the surface at "s"
!  </table>
!
! ================================================================================================ !

  FUNCTION EvaluateSlope_Surfaces( mySurfaces, s, j ) RESULT( dxds )
    IMPLICIT NONE
    CLASS( Surfaces ) :: mySurfaces
    REAL(prec)        :: s(1:2)
    REAL(prec)        :: dxds(1:3,1:2)
    INTEGER           :: j
    ! Local
    INTEGER ::  i


    DO i = 1, 3

      dxds(i,1) = mySurfaces % interp % Interpolate_2D( mySurfaces % dxds(1,:,:,i,j), s )
      dxds(i,2) = mySurfaces % interp % Interpolate_2D( mySurfaces % dxds(2,:,:,i,j), s )

    ENDDO


  END FUNCTION EvaluateSlope_Surfaces

END MODULE Surfaces_CLASS