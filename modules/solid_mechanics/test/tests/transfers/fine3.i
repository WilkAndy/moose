[Mesh]
  type = GeneratedMesh
  dim = 1
  nx = 10
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
  [divstress]
    type = DynamicStressDivergenceTensors
    component = 0
    use_displaced_mesh = false
    variable = disp_x
    implicit = false
    displacements = disp_x # this is a required parameter, but has no impact here
    save_in = fine_residual # note the BC residual does not get saved, hmmmm.  Fortunately, the coarse-mesh BC will zero its version
  []
  [massmatrix_x]
    type = MassMatrix
    density = 100 # 10 eles in fine3 mesh, 1 ele in coarse3 mesh, density = 1 in coarse3 mesh
    matrix_tags = mass
    variable = disp_x
  []
[]

[BCs]
  [left_fixed]
    type = ExplicitDirichletBC
    variable = disp_x
    boundary = left
    value = 0
    implicit = false
  []
[]

[ICs]
  [stretched]
    type = FunctionIC
    variable = disp_x
    function = '0.5 * x'
  []
[]

[Materials]
  [Elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1
    poissons_ratio = 0
    implicit = false
  []
  [strain]
    type = ComputeSmallStrain
#    implicit = false # Not sure about this - we want to use old disp_x but up-to-date coarse_disp_x.  true_disp_x gets updated on TIMESTEP_BEGIN
    displacements = true_disp_x # NOTE: different than normal solid mechanics
  []
  [stress]
    type = ComputeLinearElasticStress
    implicit = false
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
  [true_disp_x]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  [true_disp_x]
    type = ParsedAux
    variable = true_disp_x
    coupled_variables = 'coarse_disp_x disp_x'
    expression = 'coarse_disp_x + disp_x'
    execute_on = TIMESTEP_BEGIN
  []
[]

[Postprocessors]
  [end_disp_x]
    type = PointValue
    point = '1 0 0'
    variable = true_disp_x
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
