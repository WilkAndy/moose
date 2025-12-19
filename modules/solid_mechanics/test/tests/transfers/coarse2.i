[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
  xmin = 0
  xmax = 1
[]

[Variables]
  [disp_x]
    order = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  [what_we_do_not_want]
    type = NullKernel
    variable = disp_x
  []
[]
[NodalKernels]
  inactive = 'what_we_want'
  [what_we_want]
    type = CoupledForceNodalKernel
    variable = disp_x
    v = residual_from_fine
  []
[]

[AuxVariables]
  [residual_from_fine]
    order = FIRST
    family = LAGRANGE
  []
[]

[ICs]
  [linear]
    type = FunctionIC
    variable = disp_x
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
  end_time = 2
[]

[MultiApps]
  [fine]
    type = TransientMultiApp
    input_files = fine2.i
    execute_on = TIMESTEP_BEGIN
  []
[]

[Transfers]
  [get_residual_from_fine] # TODO
    type = MultiAppSumTransfer
    from_multi_app = fine
    fe_family = LAGRANGE
    fe_order = FIRST
    skip_outside_points = true
    verbose = true
    from_variable = fine_residual
    to_variable = residual_from_fine
  []
  [disp_to_fine]
    type = MultiAppSumTransfer
    to_multi_app = fine
    fe_family = LAGRANGE
    fe_order = FIRST
    skip_outside_points = true
    verbose = true
    from_variable = disp_x
    to_variable = coarse_disp_x
  []
[]
