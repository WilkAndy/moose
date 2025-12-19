[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 1
  xmin = 0
  xmax = 1
[]

[Problem]
  extra_tag_matrices = mass
[]

[Variables]
  [disp_x]
    order = FIRST
    family = LAGRANGE
  []
[]

[Kernels]
  [massmatrix_x]
    type = MassMatrix
    density = 1
    matrix_tags = mass
    variable = disp_x
  []
[]

[NodalKernels]
  [residual_from_fine]
    type = CoupledForceNodalKernel
    variable = disp_x
    v = residual_from_fine
    coef = -1 # CoupledForceNodalKernel has a negative sign in it
  []
[]

[BCs]
  [left_fixed]
    type = ExplicitDirichletBC
    variable = disp_x
    boundary = left
    value = 0
  []
[]

[AuxVariables]
  [residual_from_fine]
    order = FIRST
    family = LAGRANGE
  []
[]

[Outputs]
  exodus = true
  csv = true
[]


[Executioner]
  type = Transient

  [TimeIntegrator]
    type = ExplicitMixedOrder
    mass_matrix_tag = mass
    second_order_vars = disp_x
    use_constant_mass = true
  []

  dt = 0.5
  end_time = 10
[]


[MultiApps]
  [fine]
    type = TransientMultiApp
    input_files = fine3.i
    execute_on = TIMESTEP_BEGIN
  []
[]

[Transfers]
  [get_residual_from_fine]
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
