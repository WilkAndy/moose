# groundwater velocity is 10m.yr^-1, porosity is 0.3, so injection rate is 0.3 * 10m^3.yr^-1.m^-2 = 0.3E4kg.yr^-1.m^-2
flow_rate = 0.3E4
# To support this flow rate, the porepressure gradient = flow_rate * viscosity / permeability / density
natural_gradient = ${fparse flow_rate * 1E-3 / (3600 * 24 * 365.0) / 1E-9 / 1E3}
# To the groundwater system, CH3COO- is added at rate 4E-6 mol.m^-3.yr^-1 = 4E-6 * 59.0445E-3 kg.m^-3.yr^-1
source_CH3COO = 236.178E-9
# To the groundwater system, Ca++ is added at rate 2E-6 mol.m^-3.yr^-1 = 2E-6 * 40.08E-3 kg.m^-3.yr^-1
source_Ca = 80.16E-9
# The following are the mass-fractions of the species in the groundwater
# The numerical values can be obtained by running the geochemistry simulation with a very small timestep so no kinetics are active
eqm_Ca = 4.0084607901018e-05
eqm_SO4 = 3.842745743203e-06
eqm_CH3COO = 5.9051288195564e-11
eqm_HS = 3.3071701740266e-11
eqm_CH4 = 1.6044444376963e-11
# The following are the injection rates of each species
inj_Ca = ${fparse eqm_Ca * flow_rate}
inj_SO4 = ${fparse eqm_SO4 * flow_rate}
inj_CH3COO = ${fparse eqm_CH3COO * flow_rate}
inj_HS = ${fparse eqm_HS * flow_rate}
inj_CH4 = ${fparse eqm_CH4 * flow_rate}
inj_H2O = ${fparse (1.0 - eqm_Ca - eqm_SO4 - eqm_CH3COO - eqm_HS - eqm_CH4) * flow_rate}
# The following are scalings used in calculating the residual.  Eg, because the concentration of CH3COO is so low, its residual is always tiny, so to get better accuracy it should be scaled
scale_Ca = ${fparse 1.0 / eqm_Ca}
scale_SO4 = ${fparse 1.0 / eqm_SO4}
scale_CH3COO = ${fparse 1.0 / eqm_CH3COO}
scale_HS = ${fparse 1.0 / eqm_HS}
scale_CH4 = ${fparse 1.0 / eqm_CH4}
[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 1
    nx = 50
    xmin = 0
    xmax = 200000
  []
[]

[GlobalParams]
  PorousFlowDictator = dictator
  gravity = '0 0 0'
[]

[Variables]
  [f0] # Ca
    initial_condition = ${eqm_Ca}
    scaling = ${scale_Ca}
  []
  [f1] # SO4
    initial_condition = ${eqm_SO4}
    scaling = ${scale_SO4}
  []
  [f2] # CH3COO
    initial_condition = ${eqm_CH3COO}
    scaling = ${scale_CH3COO}
  []
  [f3] # HS
    initial_condition = ${eqm_HS}
    scaling = ${scale_HS}
  []
  [f4] # CH4
    initial_condition = ${eqm_CH4}
    scaling = ${scale_CH4}
  []
  [porepressure]
  []
[]

[ICs]
  [porepressure]
    type = FunctionIC
    variable = porepressure
    function = '(200000 - x) * ${natural_gradient}'
  []
[]

[Modules]
  [FluidProperties]
    [the_simple_fluid]
      type = SimpleFluidProperties
      viscosity = 1E-3
      thermal_expansion = 0
    []
  []
[]

[PorousFlowFullySaturated]
  coupling_type = Hydro
  porepressure = porepressure
  mass_fraction_vars = 'f0 f1 f2 f3 f4'
  save_component_rate_in = 'rate_Ca rate_SO4 rate_CH3COO rate_HS rate_CH4 rate_H2O' # change in kg at every node / dt
  fp = the_simple_fluid
  time_unit = years
[]

[Kernels]
  [add_Ca]
    type = BodyForce
    variable = f0
    value = ${source_Ca}
  []
  [add_CH3COO]
    type = BodyForce
    variable = f2
    value = ${source_CH3COO}
  []
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
    flux_function = -${inj_Ca}
  []
  [inject_SO4]
    type = PorousFlowSink
    boundary = left
    variable = f1
    flux_function = -${inj_SO4}
  []
  [inject_CH3COO]
    type = PorousFlowSink
    boundary = left
    variable = f2
    flux_function = -${inj_CH3COO}
  []
  [inject_HS]
    type = PorousFlowSink
    boundary = left
    variable = f3
    flux_function = -${inj_HS}
  []
  [inject_CH4]
    type = PorousFlowSink
    boundary = left
    variable = f4
    flux_function = -${inj_CH4}
  []
  [inject_H2O]
    type = PorousFlowSink
    boundary = left
    variable = porepressure
    flux_function = -${inj_H2O}
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
    type = DirichletBC # Could employ a PorousFlowOutflowBC here, but then porepressure can rise or fall arbitrarily
    boundary = right
    variable = porepressure
    value = 0
  []
[]

[Materials]
  [porosity]
    type = PorousFlowPorosityConst
    porosity = 0.3 # irrelevant here
  []
  [permeability]
    type = PorousFlowPermeabilityConst
    permeability = '1E-9 0 0   0 1E-9 0   0 0 1E-9'
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
  [TimeStepper]
    type = FunctionDT
  function = '0.1 * (t + 1)'
  []
  end_time = 100000
  nl_abs_tol = 1E-5
[]

[Outputs]
  csv = true
  exodus = true
[]

[MultiApps]
  [react]
    type = TransientMultiApp
    input_files = bio_zoning.i
    clone_master_mesh = true
    execute_on = 'timestep_end'
  []
[]

[Transfers]
  [changes_due_to_flow]
    type = MultiAppCopyTransfer
    direction = to_multiapp
    source_variable = 'rate_Ca rate_SO4 rate_CH3COO rate_HS rate_CH4 rate_H2O' # change in kg at every node / dt
    variable = 'pf_rate_Ca pf_rate_SO4 pf_rate_CH3COO pf_rate_HS pf_rate_CH4 pf_rate_H2O'
    multi_app = react
  []
  [massfrac_from_geochem]
    type = MultiAppCopyTransfer
    direction = from_multiapp
    source_variable = 'massfrac_Ca massfrac_SO4 massfrac_CH3COO massfrac_HS massfrac_CH4'
    variable = 'f0 f1 f2 f3 f4'
    multi_app = react
  []
[]
  
