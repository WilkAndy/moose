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
  [aux]
    order = FIRST
    family = LAGRANGE
  []
[]

[Outputs]
  exodus = true
  csv = true
  file_base = fine_out
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  dt = 1
  end_time = 1
[]
