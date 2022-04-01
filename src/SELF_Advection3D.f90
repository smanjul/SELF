!
! Copyright 2020-2022 Fluid Numerics LLC
! Author : Joseph Schoonover (joe@fluidnumerics.com)
! Support : support@fluidnumerics.com
!
! //////////////////////////////////////////////////////////////////////////////////////////////// !
MODULE SELF_Advection3D

  USE SELF_Metadata
  USE SELF_Mesh
  USE SELF_MappedData
  USE SELF_Model


  TYPE,EXTENDS(Model3D) :: Advection3D

    CONTAINS

    ! Concretized Methods
    PROCEDURE :: SourceMethod => Source_Advection3D
    PROCEDURE :: FluxMethod => Flux_Advection3D
    PROCEDURE :: RiemannSolver => RiemannSolver_Advection3D
    PROCEDURE :: SetBoundaryCondition => SetBoundaryCondition_Advection3D
    
    ! Overridden Methods
    PROCEDURE :: CalculateEntropy => CalculateEntropy_Advection3D

  END TYPE Advection3D

CONTAINS

  SUBROUTINE CalculateEntropy_Advection3D(this)
  !! Base method for calculating entropy of a model
  !! Calculates the entropy as the integration of the 
  !! squared tracer over the domain
    IMPLICIT NONE
    CLASS(Advection3D), INTENT(inout) :: this
    ! Local
    INTEGER :: i, j, k, iVar, iEl
    REAL(prec) :: Jacobian, s
    REAL(prec) :: wi,wj,wk

    ! TO DO : GPU reduction

    this % entropy = 0.0_prec

    DO iEl = 1, this % geometry % x % nElem
      DO iVar = 1, this % geometry % x % nVar
        DO k = 0, this % geometry % x % interp % N
          DO j = 0, this % geometry % x % interp % N
            DO i = 0, this % geometry % x % interp % N

              ! Coordinate mapping Jacobian
              Jacobian = this % geometry % J % interior % hostData(i,j,k,1,iEl)

              ! Quadrature weights
              wi = this % geometry % x % interp % qWeights % hostData(i) 
              wj = this % geometry % x % interp % qWeights % hostData(j) 
              wk = this % geometry % x % interp % qWeights % hostData(k) 

              ! Solution
              s = this % solution % interior % hostData(i,j,k,iVar,iEl)

              this % entropy = this % entropy + s*s*Jacobian*wi*wj*wk

            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE CalculateEntropy_Advection3D

  SUBROUTINE SetBoundaryCondition_Advection3D(this)
    IMPLICIT NONE
    CLASS(Advection3D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iSide,iEl,iVar,e2      

    DO iEl = 1, this % solution % nElem
      DO iSide = 1, 6
        DO iVar = 1, this % solution % nVar
          DO j = 0, this % solution % interp % N
            DO i = 0, this % solution % interp % N

              e2 = this % mesh % sideInfo % hostData(3,iSide,iEl)

              IF( e2 == 0 )THEN
                this % solution % extBoundary % hostData(i,j,iVar,iSide,iEl) = 0.0_prec
              ENDIF

            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO


  END SUBROUTINE SetBoundaryCondition_Advection3D 

  SUBROUTINE Source_Advection3D(this)
    IMPLICIT NONE
    CLASS(Advection3D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,k,iEl,iVar

    DO iEl = 1, this % source % nElem
      DO iVar = 1, this % source % nVar
        DO k = 0, this % source % interp % N
          DO j = 0, this % source % interp % N
            DO i = 0, this % source % interp % N

              this % source % interior % hostData(i,j,k,iVar,iEl) = 0.0_prec

            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE Source_Advection3D

  SUBROUTINE Flux_Advection3D(this)
    IMPLICIT NONE
    CLASS(Advection3D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,k,iEl,iVar

    DO iEl = 1, this % source % nElem
      DO iVar = 1, this % source % nVar
        DO k = 0, this % source % interp % N
          DO j = 0, this % source % interp % N
            DO i = 0, this % source % interp % N

              ! f1 = u*s (x-component)
              this % flux % interior % hostData(1,i,j,k,iVar,iEl) = &
                    this % velocity % interior % hostData(1,i,j,k,iVar,iEl)*&
                    this % solution % interior % hostData(i,j,k,iVar,iEl)
        
              ! f2 = v*s (y-component)
              this % flux % interior % hostData(2,i,j,k,iVar,iEl) = &
                    this % velocity % interior % hostData(2,i,j,k,iVar,iEl)*&
                    this % solution % interior % hostData(i,j,k,iVar,iEl)

              ! f3 = w*s (z-component)
              this % flux % interior % hostData(3,i,j,k,iVar,iEl) = &
                    this % velocity % interior % hostData(3,i,j,k,iVar,iEl)*&
                    this % solution % interior % hostData(i,j,k,iVar,iEl)

            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE Flux_Advection3D

  SUBROUTINE RiemannSolver_Advection3D(this)
    IMPLICIT NONE
    CLASS(Advection3D),INTENT(inout) :: this
    ! Local
    INTEGER :: i,j,iSide,iEl,iVar
    REAL(prec) :: extState, intState
    REAL(prec) :: nhat(1:3), nmag, un


    DO iEl = 1, this % solution % nElem
      DO iSide = 1, 6
        DO iVar = 1, this % solution % nVar
          DO j = 0, this % solution % interp % N
            DO i = 0, this % solution % interp % N

             ! Get the boundary normals on cell edges from the mesh geometry
             nhat(1:3) = this % geometry % nHat % boundary % hostData(1:3,i,j,1,iSide,iEl)

             ! Calculate the normal velocity at the cell edges
             un = this % velocity % boundary % hostData(1,i,j,1,iSide,iEl)*nHat(1)+&
                  this % velocity % boundary % hostData(2,i,j,1,iSide,iEl)*nHat(2)+&
                  this % velocity % boundary % hostData(3,i,j,1,iSide,iEl)*nHat(3)

             ! Pull external and internal state for the Riemann Solver (Lax-Friedrichs)
             extState = this % solution % extBoundary % hostData(i,j,iVar,iSide,iEl)
             intState = this % solution % boundary % hostData(i,j,iVar,iSide,iEl)
             nmag = this % geometry % nScale % boundary % hostData(i,j,1,iSide,iEl)

             ! Calculate the flux
             this % flux % boundaryNormal % hostData(i,j,iVar,iSide,iEl) = 0.5_prec*&
                 ( un*(intState + extState) - abs(un)*(extState - intState) )*nmag

            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE RiemannSolver_Advection3D

END MODULE SELF_Advection3D
