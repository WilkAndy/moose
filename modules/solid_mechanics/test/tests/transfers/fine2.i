[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 2
  xmin = 0
  xmax = 1
[]

[Variables]
  [v]
    order = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  [v]
    type = NullKernel
    variable = v
  []
[]

[AuxVariables]
  [coarse_disp_x]
    order = FIRST
    family = LAGRANGE
  []
  [fine_residual]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  [fine_residual]
    type = FunctionAux
    variable = fine_residual
    function = '-3 - x'
  []
[]

[Outputs]
  exodus = true
  csv = true
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  dt = 1
  end_time = 2
[]
