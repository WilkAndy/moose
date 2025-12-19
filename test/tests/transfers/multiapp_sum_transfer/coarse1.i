[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
  xmin = 0
  xmax = 1
[]

[Variables]
  [u]
    order = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  [u]
    type = NullKernel
    variable = u
  []
[]

[ICs]
  [linear]
    type = FunctionIC
    variable = u
    function = '1 + x'
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
  end_time = 1
[]

[MultiApps]
  [fine]
    type = TransientMultiApp
    input_files = fine1.i
    execute_on = TIMESTEP_END
  []
[]

[Transfers]
  [sum_T]
    type = MultiAppSumTransfer
    to_multi_app = fine
    fe_family = LAGRANGE
    fe_order = FIRST
    skip_outside_points = true
    verbose = true
    from_variable = u
    to_variable = aux
    execute_on = TIMESTEP_END
  []
[]
