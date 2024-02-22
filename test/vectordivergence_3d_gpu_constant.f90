program test

  implicit none
  integer :: exit_code
  
  exit_code = vectordivergence_3d_gpu_constant()
  stop exit_code

contains
integer function vectordivergence_3d_gpu_constant() result(r)
  use SELF_Constants
  use SELF_Lagrange
  use SELF_Data
  use iso_c_binding
  use hipfort_hipblas

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
  type(Vector3D) :: f
  type(Scalar3D) :: df
  type(Lagrange),target :: interp
  type(c_ptr) :: handle

  call hipblasCheck(hipblasCreate(handle))

  ! Create an interpolant
  call interp % Init(N=controlDegree, &
                     controlNodeType=GAUSS, &
                     M=targetDegree, &
                     targetNodeType=UNIFORM)

  ! Initialize vectors
  call f % Init(interp,nvar,nelem)

  call df % Init(interp,nvar,nelem)

  ! Set the source vector (on the control grid) to a non-zero constant
  f % interior  = 1.0_prec

  call f % updatedevice()

  call f % Divergence(df, handle)

  call hipcheck(hipdevicesynchronize())

  ! Calculate diff from exact
  df % interior  = abs(df % interior  - 0.0_prec)

  if (maxval(df % interior ) <= tolerance) then
    r = 0
  else
    r = 1
  end if

  call f % free()
  call df % free()
  call interp % free()
  call hipblasCheck(hipblasDestroy(handle))

end function vectordivergence_3d_gpu_constant
end program test
