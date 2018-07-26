! SELF_Fluids_Driver.f90
!
! Copyright 2018 Joseph Schoonover <joe@myFluidnumerics.consulting>, Fluid Numerics LLC
! All rights reserved.
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !


PROGRAM SELF_Fluids_Driver

  USE ModelPrecision
  USE ModelParameters_Class
  USE HexMesh_Class
  USE Fluid_EquationParser_Class
  USE Fluid_Class

  IMPLICIT NONE

! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !

  TYPE( Fluid )                :: myFluid
  TYPE( Fluid_EquationParser ) :: myFluidConditions
  LOGICAL                      :: setupSuccess
  LOGICAL                      :: initializeFromScratch
  LOGICAL                      :: pickupFileExists
  LOGICAL                      :: run_MeshGenOnly, run_UpToInitOnly


! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !


  CALL Setup( )

  IF( setupSuccess )THEN

    IF( run_MeshGenOnly )THEN

      CALL MeshGen( )
    
    ELSE

      CALL Initialize( )

      IF( .NOT. run_UpToInitOnly )THEN
        CALL MainLoop( )
      ENDIF

      CALL Cleanup( )

    ENDIF


  ENDIF

CONTAINS

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !

  SUBROUTINE Setup( )
    IMPLICIT NONE
    ! Local
    INTEGER :: nArg, argID
    CHARACTER(500) :: argName
    LOGICAL :: helpNeeded

    helpNeeded       = .FALSE.
    run_MeshGenOnly  = .FALSE.
    run_UpToInitOnly = .FALSE.

    nArg = command_argument_count( )

    IF( nArg > 0 )THEN

      CALL get_command_argument( 1, argName )

      SELECT CASE( TRIM( argName ) )

        CASE( "meshgen" )

          run_MeshGenOnly  = .TRUE.
          run_UpToInitOnly = .FALSE.
          setupSuccess     = .TRUE.

          IF( nArg > 1 ) THEN
            PRINT*, '  Too many command line arguments for sfluid meshgen.'
            helpNeeded   = .TRUE.
            setupSuccess = .FALSE.
          ENDIF

        CASE( "init" )

          run_MeshGenOnly  = .FALSE.
          run_UpToInitOnly = .TRUE.
          setupSuccess     = .TRUE.

          IF( nArg > 1 ) THEN
            PRINT*, '  Too many command line arguments for sfluid init.'
            helpNeeded   = .TRUE.
            setupSuccess = .FALSE.
          ENDIF

         CASE( "help" )
            helpNeeded   = .TRUE.
            setupSuccess = .FALSE.

         CASE DEFAULT
            run_MeshGenOnly  = .FALSE.
            run_UpToInitOnly = .FALSE.
            helpNeeded       = .FALSE.

      END SELECT

    ENDIF

    IF( helpNeeded ) THEN

      PRINT*, 'SELF-Fluids (sfluid) Command Line Tool'      
      PRINT*, ' '
      PRINT*, ' A program for solving Compressible Navier-Stokes using the'
      PRINT*, ' Nodal Discontinuous Galerkin Spectral Element Method.'
      PRINT*, ' '
      PRINT*, '  sfluid [tool]'      
      PRINT*, ' '
      PRINT*, ' [tool] can be :'
      PRINT*, ' '
      PRINT*, '   help'
      PRINT*, '     Display this help message'
      PRINT*, ' '
      PRINT*, '   meshgen'
      PRINT*, '     Run only the mesh generator to generate a structured mesh.'
      PRINT*, '     The structured mesh is built using nXElem, nYElem, and nZElem'
      PRINT*, '     specified in runtime.params. Further, domain decomposition '
      PRINT*, '     for the structured mesh is done by setting the number of'
      PRINT*, '     of processes in each direction (nProcX, nProcY, nProcZ)'
      PRINT*, '     Topography can be set using an equation like'
      PRINT*, ' '
      PRINT*, '           h = exp( -(x-500.0)^2/200.0 )'
      PRINT*, ' '
      PRINT*, '     in the self.equations file. This will result in a terrain-following'
      PRINT*, '     structured mesh.'
      PRINT*, ' '
      PRINT*, '     Future releases of SELF-Fluids will offer more complete support'
      PRINT*, '     for working with typical unstructured mesh formats. '
      PRINT*, ' '
      PRINT*, '   init'
      PRINT*, '     Run up to the initial condition generation and do not forward'
      PRINT*, '     step the model. The initial conditions are read in from the '
      PRINT*, '     self.equations file. '
      PRINT*, ' '
      PRINT*, ' '
      PRINT*, ' '

      setupSuccess = .FALSE.
      RETURN
    ENDIF

    IF( .NOT. run_MeshGenOnly )THEN
      CALL myFluid % Build( setupSuccess )
    ENDIF


  END SUBROUTINE Setup

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !

  SUBROUTINE MeshGen( )


      CALL myFluid % ExtComm % SetRanks( )

      IF( myFluid % ExtComm % myRank == 0 )THEN
        PRINT*, '  Generating structured mesh...'
        CALL StructuredMeshGenerator_3D( )
      ENDIF

      CALL myFluid % ExtComm % Finalize( )

      PRINT*, '  Done'

  END SUBROUTINE MeshGen

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !

  SUBROUTINE Initialize( )

    ! Attempt to read the fluid pickup file. If it doesn't exist, this routine
    ! returns FALSE.
    CALL myFluid % ReadPickup( pickupFileExists )

    ! If the pickup file doesn't exist, then the initial conditions are generated
    ! from the equation parser.
    IF( .NOT. pickupFileExists )THEN

      PRINT(MsgFMT), 'Pickup file not found.'

    ENDIF

    IF( .NOT. pickupFileExists .OR. run_UpToInitOnly )THEN

      PRINT(MsgFMT), 'Attempting initial condition generation from self.equations'
      CALL myFluidConditions % Build( 'self.equations' )
      CALL InitialConditions( )

      CALL myFluid % WritePickup( )
      CALL myFluid % WriteTecplot( )

    ENDIF

#ifdef HAVE_DIAGNOSTICS
    CALL myFluid % WriteDiagnostics( )
#endif

  END SUBROUTINE Initialize

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !

  SUBROUTINE InitialConditions( )
    ! Local
    INTEGER    :: i, j, k, iEl
    INTEGER    :: iFace, bID, e1, s1, e2
    REAL(prec) :: x(1:3)
    REAL(prec) :: T, Tbar, u, v, w, rho, rhobar, s, s0


    myFluid % state % solution = 0.0_prec
    !$OMP PARALLEL
    CALL myFluid % CalculateStaticState( ) !! CPU Kernel
    !$OMP END PARALLEL

    DO iEl = 1, myFluid % mesh % elements % nElements
      DO k = 0, myFluid % params % polyDeg
        DO j = 0, myFluid % params % polyDeg
          DO i = 0, myFluid % params % polyDeg

            x(1:3) = myFluid % mesh % elements % x(i,j,k,1:3,iEl)

            IF( myFluidConditions % calculate_density_from_T )THEN
               
              u  = myFluidConditions % u % evaluate( x ) 
              v  = myFluidConditions % v % evaluate( x ) 
              w  = myFluidConditions % w % evaluate( x ) 
              T  = myFluidConditions % t % evaluate( x ) ! Potential temperature anomaly
              s  = myFluidConditions % tracer % evaluate( x )
              s0 = myFluidConditions % staticTracer % evaluate( x )

              Tbar = myFluid % static % solution(i,j,k,5,iEl)/myFluid % static % solution(i,j,k,4,iEl)

              myFluid % state % solution(i,j,k,4,iEl) = -myFluid % static % solution(i,j,k,4,iEl)*T/(Tbar + T)
              myFluid % state % solution(i,j,k,1,iEl) = ( myFluid % state % solution(i,j,k,4,iEl) + myFluid % static % solution(i,j,k,4,iEl) )*u
              myFluid % state % solution(i,j,k,2,iEl) = ( myFluid % state % solution(i,j,k,4,iEl) + myFluid % static % solution(i,j,k,4,iEl) )*v
              myFluid % state % solution(i,j,k,3,iEl) = ( myFluid % state % solution(i,j,k,4,iEl) + myFluid % static % solution(i,j,k,4,iEl) )*w
#ifdef PASSIVE_TRACERS
              myFluid % state % solution(i,j,k,6,iEl) = ( myFluid % state % solution(i,j,k,4,iEl) + myFluid % static % solution(i,j,k,4,iEl) )*s
              myFluid % static % solution(i,j,k,6,iEl) = ( myFluid % state % solution(i,j,k,4,iEl) + myFluid % static % solution(i,j,k,4,iEl) )*s0
#endif

            ELSE

              u   = myFluidConditions % u % evaluate( x ) 
              v   = myFluidConditions % v % evaluate( x ) 
              w   = myFluidConditions % w % evaluate( x ) 
              rho = myFluidConditions % rho % evaluate( x ) 
              T   = myFluidConditions % t % evaluate( x ) ! Potential temperature anomaly
              s = myFluidConditions % tracer % evaluate( x )


              myFluid % state % solution(i,j,k,4,iEl) = rho
              myFluid % state % solution(i,j,k,1,iEl) = ( myFluid % state % solution(i,j,k,4,iEl) + myFluid % static % solution(i,j,k,4,iEl) )*u
              myFluid % state % solution(i,j,k,2,iEl) = ( myFluid % state % solution(i,j,k,4,iEl) + myFluid % static % solution(i,j,k,4,iEl) )*v
              myFluid % state % solution(i,j,k,3,iEl) = ( myFluid % state % solution(i,j,k,4,iEl) + myFluid % static % solution(i,j,k,4,iEl) )*w
              myFluid % state % solution(i,j,k,5,iEl) = ( myFluid % state % solution(i,j,k,4,iEl) + myFluid % static % solution(i,j,k,4,iEl) )*T
#ifdef PASSIVE_TRACERS
              myFluid % state % solution(i,j,k,6,iEl) = ( myFluid % state % solution(i,j,k,4,iEl) + myFluid % static % solution(i,j,k,4,iEl) )*s
              myFluid % static % solution(i,j,k,6,iEl) = ( myFluid % state % solution(i,j,k,4,iEl) + myFluid % static % solution(i,j,k,4,iEl) )*s0
#endif

            ENDIF

            myFluid % sourceTerms % drag(i,j,k,iEl) = myFluidConditions % drag % evaluate( x )

          ENDDO
        ENDDO
      ENDDO
    ENDDO

#ifdef HAVE_CUDA
    myFluid % state % solution_dev = myFluid % state % solution
    CALL myFluid % sourceTerms % UpdateDevice( )
#endif

    !$OMP PARALLEL
    CALL myFluid % EquationOfState( )
    !$OMP END PARALLEL

#ifdef HAVE_CUDA
    myFluid % state % boundarySolution  = myFluid % state % boundarySolution_dev
#endif

    CALL myFluid % UpdateExternalStaticState( )

    DO bID = 1, myFluid % extComm % nBoundaries

       iFace = myFluid % extComm % boundaryIDs( bID )
       e1    = myFluid % mesh % faces % elementIDs(1,iFace)
       s1    = myFluid % mesh % faces % elementSides(1,iFace)
       e2    = myFluid % mesh % faces % elementIDs(2,iFace)

       IF( e2 == PRESCRIBED )THEN


         IF( myFluidConditions % calculate_density_from_T )THEN

           DO j = 0, myFluid % params % polyDeg
             DO i = 0, myFluid % params % polyDeg

               x(1:3) = myFluid % mesh % elements % xBound(i,j,1:3,s1,e1)

               u = myFluidConditions % u % evaluate( x ) 
               v = myFluidConditions % v % evaluate( x ) 
               w = myFluidConditions % w % evaluate( x ) 
               T = myFluidConditions % t % evaluate( x ) ! Potential temperature anomaly
               s = myFluidConditions % tracer % evaluate( x )
      
               Tbar = myFluid % static % boundarySolution(i,j,5,s1,e1)/myFluid % static % boundarySolution(i,j,4,s1,e1)
      
               myFluid % state % prescribedState(i,j,4,bID) = -myFluid % static % boundarySolution(i,j,4,s1,e1)*T/(Tbar + T)
               myFluid % state % prescribedState(i,j,1,bID) = ( myFluid % state % prescribedState(i,j,4,bID) + myFluid % static % boundarySolution(i,j,4,s1,e1) )*u
               myFluid % state % prescribedState(i,j,2,bID) = ( myFluid % state % prescribedState(i,j,4,bID) + myFluid % static % boundarySolution(i,j,4,s1,e1) )*v
               myFluid % state % prescribedState(i,j,3,bID) = ( myFluid % state % prescribedState(i,j,4,bID) + myFluid % static % boundarySolution(i,j,4,s1,e1) )*w
#ifdef PASSIVE_TRACERS
               myFluid % state % prescribedState(i,j,6,bID) = ( myFluid % state % prescribedState(i,j,4,bID) + myFluid % static % boundarySolution(i,j,4,s1,e1) )*s
#endif

             ENDDO
           ENDDO
          
         ELSE

           DO j = 0, myFluid % params % polyDeg
             DO i = 0, myFluid % params % polyDeg

               x(1:3) = myFluid % mesh % elements % xBound(i,j,1:3,s1,e1)

               u = myFluidConditions % u % evaluate( x ) 
               v = myFluidConditions % v % evaluate( x ) 
               w = myFluidConditions % w % evaluate( x ) 
               T = myFluidConditions % t % evaluate( x ) ! Potential temperature anomaly
               rho = myFluidConditions % rho % evaluate( x ) 
               s = myFluidConditions % tracer % evaluate( x )
  
               Tbar = myFluid % static % boundarySolution(i,j,5,s1,e1)/myFluid % static % boundarySolution(i,j,4,e1,s1)
  
               myFluid % state % prescribedState(i,j,4,bID) = rho
               myFluid % state % prescribedState(i,j,1,bID) = ( rho + myFluid % static % boundarySolution(i,j,4,s1,e1) )*u
               myFluid % state % prescribedState(i,j,2,bID) = ( rho + myFluid % static % boundarySolution(i,j,4,s1,e1) )*v
               myFluid % state % prescribedState(i,j,3,bID) = ( rho + myFluid % static % boundarySolution(i,j,4,s1,e1) )*w
               myFluid % state % prescribedState(i,j,5,bID) = ( rho + myFluid % static % boundarySolution(i,j,4,s1,e1) )*T
#ifdef PASSIVE_TRACERS
               myFluid % state % prescribedState(i,j,6,bID) = ( rho + myFluid % static % boundarySolution(i,j,4,s1,e1) )*s
#endif

             ENDDO
           ENDDO

         ENDIF

       ENDIF

    ENDDO


  END SUBROUTINE InitialConditions

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !

  SUBROUTINE Cleanup( )
    IMPLICIT NONE

    CALL myFluid % Trash( )

  END SUBROUTINE Cleanup

! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !

  SUBROUTINE MainLoop( )
    IMPLICIT NONE
    INTEGER    :: iT
! ------------------------------------------------------------------------------ !
! <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><> !
! ------------------------------------------------------------------------------ !
    DO iT = 1, myFluid % params % nDumps ! Loop over time-steps

      CALL myFluid % ForwardStepRK3( myFluid % params % nStepsPerDump ) ! Forward Step

#ifdef HAVE_CUDA
      myFluid % state % solution = myFluid % state % solution_dev ! Update the host from the GPU
#endif

      CALL myFluid % WritePickup( )
      CALL myFluid % WriteTecplot( )

#ifdef HAVE_DIAGNOSTICS
      CALL myFluid % Diagnostics( )
      CALL myFluid % WriteDiagnostics( )
#endif

    ENDDO


  END SUBROUTINE MainLoop


END PROGRAM SELF_Fluids_Driver

