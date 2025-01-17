! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !
!
! Maintainers : support@fluidnumerics.com
! Official Repository : https://github.com/FluidNumerics/self/
!
! Copyright © 2024 Fluid Numerics LLC
!
! Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in
!    the documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS” AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
! THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// !

program test

  implicit none
  integer :: exit_code

  exit_code = mappedscalarbrderivative_1d_constant()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains
  integer function mappedscalarbrderivative_1d_constant() result(r)
    use SELF_Constants
    use SELF_Lagrange
    use SELF_MappedScalar_1D
    use SELF_Mesh_1D
    use SELF_Geometry_1D

    implicit none

    integer,parameter :: controlDegree = 7
    integer,parameter :: targetDegree = 16
    integer,parameter :: nvar = 1
    integer,parameter :: nelem = 100
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 10.0_prec**(-7)
#else
    real(prec),parameter :: tolerance = 10.0_prec**(-3)
#endif
    type(MappedScalar1D) :: f
    type(MappedScalar1D) :: df
    type(Lagrange),target :: interp
    type(Mesh1D),target :: mesh
    type(Geometry1D),target :: geometry

    call mesh%StructuredMesh(nElem=nelem, &
                             x=(/0.0_prec,10.0_prec/))
    ! Create an interpolant
    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    ! Generate geometry (metric terms) from the mesh elements
    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    ! Initialize scalars
    call f%Init(interp,nvar,nelem)
    call df%Init(interp,nvar,nelem)
    call f%AssociateGeometry(geometry)
    call df%AssociateGeometry(geometry)

    call f%SetEquation(1,'f = 1.0')
    call f%SetInteriorFromEquation(0.0_prec)
    print*,"min, max (interior)",minval(f%interior),maxval(f%interior)

    call f%BoundaryInterp()
    call f%UpdateHost()
    print*,"min, max (boundary)",minval(f%boundary),maxval(f%boundary)

    call f%SideExchange(mesh)
    ! Set boundary conditions
    f%extBoundary(1,1,1) = 1.0_prec ! Left most
    f%extBoundary(2,nelem,1) = 1.0_prec ! Right most
    print*,"min, max (extboundary)",minval(f%extBoundary),maxval(f%extBoundary)

    call f%UpdateDevice()
    call f%AverageSides()
    call f%UpdateHost()
    print*,"min, max (avgboundary)",minval(f%avgboundary),maxval(f%avgboundary)

    ! Compute "fluxes"
    f%boundarynormal(1,:,:) = -f%avgBoundary(1,:,:) ! Account for left facing normal
    f%boundarynormal(2,:,:) = f%avgBoundary(2,:,:) ! Account for right facing normal

    call f%UpdateDevice()
#ifdef ENABLE_GPU
    call f%MappedDGDerivative(df%interior_gpu)
#else
    call f%MappedDGDerivative(df%interior)
#endif
    call df%UpdateHost()

    ! Calculate diff from exact
    df%interior = abs(df%interior-0.0_prec)

    if(maxval(df%interior) <= tolerance) then
      r = 0
    else
      print*,"Max error : ",maxval(df%interior)
      r = 1
    endif

    ! Clean up
    call f%DissociateGeometry()
    call df%DissociateGeometry()

    call mesh%Free()
    call geometry%Free()
    call interp%free()
    call f%free()
    call df%free()

  endfunction mappedscalarbrderivative_1d_constant
endprogram test
