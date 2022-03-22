[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 50
    xmin = 0
    xmax = 200 # km
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]

[Variables]
  [f0] # Ca
    initial_condition = 4.0084607901018e-05
  []
  [f1] # SO4
    initial_condition = 3.842745743203e-06
  []
  [f2] # CH3COO
    initial_condition = 5.9051288195564e-11
  []
  [f3] # HS
    initial_condition = 3.3071701740266e-11
  []
  [f4] # CH4
    initial_condition = 1.6044444376963e-11
  []
  [porepressure]
    initial_condition = 0
  []
[]

[Modules]
  [FluidProperties]
    [the_simple_fluid]
      type = SimpleFluidProperties
      thermal_expansion = 0
      bulk_modulus = 2E9
      viscosity = 1E-3
      density0 = 1000
    []
  []
[]

[PorousFlowFullySaturated]
  coupling_type = Hydro
  porepressure = porepressure
  temperature = 0
  mass_fraction_vars = 'f0 f1 f2 f3 f4'
  save_component_rate_in = 'rate_Ca rate_SO4 rate_CH3COO rate_HS rate_CH4 rate_H2O' # change in kg at every node / dt
  fp = the_simple_fluid
  temperature_unit = Celsius
[]

[AuxVariables]
  [rate_Ca]
  []
  [rate_SO4]
  []
  [rate_CH3COO]
  []
  [rate_HS]
  []
  [rate_CH4]
  []
  [rate_H2O]
  []
[]

[BCs]
  [inject_Ca]
    type = PorousFlowSink
    boundary = left
    variable = f0
    flux_function = -1
  []
  [remove_Ca]
    type = PorousFlowOutflowBC
    boundary = right
    include_relperm = false
    mass_fraction_component = 0
    variable = f0
  []
  [remove_SO4]
    type = PorousFlowOutflowBC
    boundary = right
    include_relperm = false
    mass_fraction_component = 1
    variable = f1
  []
  [remove_CH3COO]
    type = PorousFlowOutflowBC
    boundary = right
    include_relperm = false
    mass_fraction_component = 2
    variable = f2
  []
  [remove_HS]
    type = PorousFlowOutflowBC
    boundary = right
    include_relperm = false
    mass_fraction_component = 3
    variable = f3
  []
  [remove_CH4]
    type = PorousFlowOutflowBC
    boundary = right
    include_relperm = false
    mass_fraction_component = 4
    variable = f4
  []
  [remove_H2O]
    type = PorousFlowOutflowBC
    boundary = right
    include_relperm = false
    mass_fraction_component = 5
    variable = porepressure
  []
[]

[Materials]
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.3
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1E-12 0 0   0 1E-12 0   0 0 1E-12'
  []
[]

[Preconditioning]
  [typically_efficient]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_hypre_type'
    petsc_options_value = ' hypre    boomeramg'
  []
[]
  
[Executioner]
  type = Transient
  solve_type = Newton
  end_time = 1000
  dt = 10
  #nl_abs_tol = 1E-4
[]

[Outputs]
  csv = true
  exodus = true
[]
