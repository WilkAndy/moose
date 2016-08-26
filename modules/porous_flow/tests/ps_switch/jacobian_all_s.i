# Newton cooling from a bar.  1-phase transient
[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1 # 10
  ny = 1
  xmin = 0
  xmax = 100
  ymin = 0
  ymax = 1
[]

[GlobalParams]
  PorousFlowDictator = dictator
[]

[UserObjects]
  [./dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'thevar'
    number_fluid_phases = 1
    number_fluid_components = 1
  [../]
  [./setter]
    type = PorousFlow1PhaseSet
    p_or_s_var = thevar
    encoder = encoder
    porepressure = porepressure
    saturation = saturation
    execute_on = 'timestep_begin'
  [../]
[]

[Variables]
  [./thevar]
  [../]
[]

[AuxVariables]
  [./porepressure]
  [../]
  [./saturation]
  [../]
  [./encoder]
  [../]
[]

[ICs]
  [./pp]
    type = FunctionIC
    variable = porepressure
    function = '1000*(-10-x)'
  [../]
  [./sat]
    type = FunctionIC
    variable = saturation
    function = 'if(x<-10,1,pow(1.0+pow(-0.1*(-10-x),1.0/(1.0-0.8)),-0.8))' # 'if(x<50,1,pow(1.0+pow(-0.1*(50-x),1.0/(1.0-0.8)),-0.8))'
  [../]
[]

[AuxKernels]
  [./updater_pp]
    type = PorousFlowVariableSwitchUpdater
    p_or_s_variable = thevar
    encoder_auxvar = encoder
    m = 0.8
    al = 1E-4
    updating_porepressure = true
    variable = porepressure
    execute_on = 'linear timestep_end'
  [../]
  [./updater_sat]
    type = PorousFlowVariableSwitchUpdater
    p_or_s_variable = thevar
    encoder_auxvar = encoder
    m = 0.8
    al = 1E-4
    updating_porepressure = false
    variable = saturation
    execute_on = 'linear timestep_end'
  [../]
[]

[Kernels]
  [./mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = thevar
  [../]
  [./flux]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    gravity = '0 0 0'
    variable = thevar
  [../]
[]

[Materials]
  [./temperature]
    type = PorousFlowTemperature
  [../]
  [./nnn]
    type = PorousFlowNodeNumber
    on_initial_only = true
  [../]
  [./ppss_new]
    type = PorousFlow1PhasePorS_VG
    primary_variable = thevar
    decoder = encoder
    porepressure = porepressure
    saturation = saturation
    m = 0.8
    al = 1E-4
  [../]
  [./massfrac]
    type = PorousFlowMassFraction
  [../]
  [./dens0]
    type = PorousFlowDensityConstBulk
    density_P0 = 1000
    bulk_modulus = 1.0E6
    phase = 0
  [../]
  [./dens_all]
    type = PorousFlowJoiner
    include_old = true
    material_property = PorousFlow_fluid_phase_density
  [../]
  [./dens_all_at_quadpoints]
    type = PorousFlowJoiner
    material_property = PorousFlow_fluid_phase_density_qp
    at_qps = true
  [../]
  [./porosity]
    type = PorousFlowPorosityConst
    porosity = 0.1
  [../]
  [./permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1E-15 0 0 0 1E-15 0 0 0 1E-15'
  [../]
  [./relperm]
    type = PorousFlowRelativePermeabilityCorey # irrelevant in this fully-saturated situation
    n_j = 2
    phase = 0
  [../]
  [./relperm_all]
    type = PorousFlowJoiner
    material_property = PorousFlow_relative_permeability
  [../]
  [./visc0]
    type = PorousFlowViscosityConst
    viscosity = 1E-3
    phase = 0
  [../]
  [./visc_all]
    type = PorousFlowJoiner
    material_property = PorousFlow_viscosity
  [../]
[]

[BCs]
  [./left]
    type = DirichletBC
    variable = thevar
    boundary = left
    value = 1E5
  [../]
[]

[VectorPostprocessors]
  [./porepressure]
    type = LineValueSampler
    variable = thevar
    start_point = '0 0.5 0'
    end_point = '100 0.5 0'
    sort_by = x
    num_points = 20
    execute_on = timestep_end
  [../]
[]

[Preconditioning]
  active = 'andy'
  [./andy]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason -snes_test_display'
    petsc_options_iname = '-ksp_type -pc_type -snes_atol -snes_rtol -snes_max_it -snes_type'
    petsc_options_value = 'bcgs bjacobi 1E-12 1E-5 20 test'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = 30E8
  dt = 1E8
[]

[Outputs]
  file_base = ps03_with_switching
  exodus = true
  print_linear_residuals = true
  print_perf_log = true
  [./along_line]
    type = CSV
    execute_vector_postprocessors_on = final
  [../]
[]
