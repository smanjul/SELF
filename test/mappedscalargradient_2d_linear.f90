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

  exit_code = mappedscalargradient_2d_linear()
  if(exit_code /= 0) then
    stop exit_code
  endif

contains
  integer function mappedscalargradient_2d_linear() result(r)

    use SELF_constants
    use SELF_Lagrange
    use SELF_Mesh_2D
    use SELF_Geometry_2D
    use SELF_MappedScalar_2D
    use SELF_MappedVector_2D

    implicit none

    integer,parameter :: controlDegree = 7
    integer,parameter :: targetDegree = 16
    integer,parameter :: nvar = 1
#ifdef DOUBLE_PRECISION
    real(prec),parameter :: tolerance = 10.0_prec**(-7)
#else
    real(prec),parameter :: tolerance = 10.0_prec**(-2)
#endif
    type(Lagrange),target :: interp
    type(Mesh2D),target :: mesh
    type(SEMQuad),target :: geometry
    type(MappedScalar2D) :: f
    type(MappedVector2D) :: df
    character(LEN=255) :: WORKSPACE
    integer :: iel,j,i
    integer(HID_T) :: fileId

    ! Create an interpolant
    call interp%Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

    ! Create a uniform block mesh
    call get_environment_variable("WORKSPACE",WORKSPACE)
    call mesh%Read_HOPr(trim(WORKSPACE)//"/share/mesh/Block2D/Block2D_mesh.h5")

    ! Generate geometry (metric terms) from the mesh elements
    call geometry%Init(interp,mesh%nElem)
    call geometry%GenerateFromMesh(mesh)

    call f%Init(interp,nvar,mesh%nelem)
    call df%Init(interp,nvar,mesh%nelem)
    call f%AssociateGeometry(geometry)

    call f%SetEquation(1,'f = x*y')

    call f%SetInteriorFromEquation(geometry,0.0_prec)
    print*,"min, max (interior)",minval(f%interior),maxval(f%interior)

#ifdef ENABLE_GPU
    call f%MappedGradient(df%interior_gpu)
#else
    call f%MappedGradient(df%interior)
#endif
    call df%UpdateHost()

    ! Calculate diff from exact
    do iel = 1,mesh%nelem
      do j = 1,controlDegree+1
        do i = 1,controlDegree+1
          df%interior(i,j,iel,1,1) = abs(df%interior(i,j,iel,1,1)-geometry%x%interior(i,j,iel,1,2)) ! df/dx = y
          df%interior(i,j,iel,1,2) = abs(df%interior(i,j,iel,1,2)-geometry%x%interior(i,j,iel,1,1)) ! df/dy = x
        enddo
      enddo
    enddo

    print*,maxval(df%interior),tolerance
    if(maxval(df%interior) <= tolerance) then
      r = 0
    else
      r = 1
    endif

    ! Clean up
    call f%DissociateGeometry()
    call geometry%Free()
    call mesh%Free()
    call interp%Free()
    call f%free()
    call df%free()

    r = 0

  endfunction mappedscalargradient_2d_linear
endprogram test
