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
  point = '0 0 0'
  reactor = reactor
[]
  
[SpatialReactionSolver]
  model_definition = definition
  geochemistry_reactor_name = reactor
  swap_into_basis = 'Siderite'
  swap_out_of_basis = 'Fe++'
  charge_balance_species = "HCO3-"
  constraint_species = "H2O              Ca++             HCO3-            SO4--            CH3COO-          HS-              CH4(aq)          Siderite         H+" 
# ASSUME that 1 litre of solution initially contains:
  constraint_value = "  1.0              1E-3             2E-3             0.04E-3          1E-9             1E-9             1E-9             1               -7.5"
  constraint_meaning = "kg_solvent_water bulk_composition bulk_composition bulk_composition bulk_composition bulk_composition bulk_composition free_mineral     log10activity"
  constraint_unit = "   kg               moles            moles            moles            moles            moles            moles            cm3            dimensionless"
  controlled_activity_name = 'H+'
  controlled_activity_value = 3.16227E-8 # this is pH=7.5
  kinetic_species_name = "sulfate_reducer methanogen"
  kinetic_species_initial_value = '1E-6 1E-6' # molecular weight of the microbes = 1, so with 1kg solvent water, this is 1E-6 mmolal
  kinetic_species_unit = 'mg mg'
  source_species_names = "H2O              Ca++             SO4--            CH3COO-          HS-              CH4(aq)"
  source_species_rates = "rate_H2O_per_1l  rate_Ca_per_1l   rate_SO4_per_1l  rate_CH3COO_per_1l rate_HS_per_1l rate_CH4_per_1l"
  ramp_max_ionic_strength_initial = 0
  execute_console_output_on = 'INITIAL TIMESTEP_END'
  solver_info = true
  evaluate_kinetic_rates_always = true
  precision = 16
[]

[UserObjects]
  [rate_sulfate_reducer]
    type = GeochemistryKineticRate
    kinetic_species_name = "sulfate_reducer"
    intrinsic_rate_constant = 31.536 # 1E-9 mol/mg/s = 31.536 mol/g/year
    multiply_by_mass = true
    promoting_species_names = 'CH3COO- SO4--'
    promoting_indices = '1 1'
    promoting_monod_indices = '1 1'
    promoting_half_saturation = '70E-6 200E-6'
    direction = dissolution
    kinetic_biological_efficiency = 4.3
    energy_captured = 45E3
    theta = 0.2
    eta = 1
  []
  [death_sulfate_reducer]
    type = GeochemistryKineticRate
    kinetic_species_name = "sulfate_reducer"
    intrinsic_rate_constant = 0.031536 # 1E-9/s = 0.031536/year
    multiply_by_mass = true
    direction = death
  []
  [rate_methanogen]
    type = GeochemistryKineticRate
    kinetic_species_name = "methanogen"
    intrinsic_rate_constant = 63.072 # 2E-9 mol/mg/s = 63.072 mol/g/year
    multiply_by_mass = true
    promoting_species_names = 'CH3COO-'
    promoting_indices = '1'
    promoting_monod_indices = '1'
    promoting_half_saturation = '20E-3'
    direction = dissolution
    kinetic_biological_efficiency = 2.0
    energy_captured = 24E3
    theta = 0.5
    eta = 1
  []
  [death_methanogen]
    type = GeochemistryKineticRate
    kinetic_species_name = "methanogen"
    intrinsic_rate_constant = 0.031536 # 1E-9/s = 0.031536/year
    multiply_by_mass = true
    direction = death
  []
  [definition]
    type = GeochemicalModelDefinition
    database_file = "db_sulf_meth.json"
    basis_species = "H2O H+ CH3COO- CH4(aq) HS- Ca++ HCO3- SO4-- Fe++"
    kinetic_redox = "sulfate_reducer methanogen"
    equilibrium_minerals = "Siderite Mackinawite"
    kinetic_rate_descriptions = "rate_sulfate_reducer death_sulfate_reducer death_methanogen"
  []
  [nodal_void_volume_uo]
    type = NodalVoidVolume
    porosity = porosity
    execute_on = 'initial timestep_end' # "initial" means this is evaluated properly for the first timestep
  []
[]


[Executioner]
  type = Transient
  dt = 1 # years
  end_time = 1E-10
[]

[AuxVariables]
  [porosity]
    initial_condition = 0.3
  []
  [nodal_void_volume]
  []
  [pf_rate_H2O] # change in H2O mass (kg/s) at each node provided by the porous-flow simulation
  []
  [pf_rate_Ca] # change in Ca++ mass (kg/s) at each node provided by the porous-flow simulation
  []
  [pf_rate_SO4] # change in SO4-- mass (kg/s) at each node provided by the porous-flow simulation
  []
  [pf_rate_CH3COO] # change in CH3COO- (kg/s) at each node provided by the porous-flow simulation
  []
  [pf_rate_HS] # change in HS- (kg/s) at each node provided by the porous-flow simulation
  []
  [pf_rate_CH4] # change in CH4(aq) (kg/s) at each node provided by the porous-flow simulation
  []

  [rate_H2O_per_1l] # change in H2O per 1 litre of aqueous solution that we consider at each node
  []
  [rate_Ca_per_1l] # change in Ca++ per 1 litre of aqueous solution that we consider at each node
  []
  [rate_SO4_per_1l] # change in SO4-- per 1 litre of aqueous solution that we consider at each node
  []
  [rate_CH3COO_per_1l] # change in CH3COO- (kg/s) per 1 litre of aqueous solution that we consider at each node
  []
  [rate_HS_per_1l] # change in HS- (kg/s) per 1 litre of aqueous solution that we consider at each node
  []
  [rate_CH4_per_1l] # change in CH4(aq) (kg/s) per 1 litre of aqueous solution that we consider at each node
  []

  [transported_H2O]
  []
  [transported_Ca]
  []
  [transported_SO4]
  []
  [transported_CH3COO]
  []
  [transported_HS]
  []
  [transported_CH4]
  []
  [transported_mass]
  []

  [massfrac_H2O]
  []
  [massfrac_Ca]
  []
  [massfrac_SO4]
  []
  [massfrac_CH3COO]
  []
  [massfrac_HS]
  []
  [massfrac_CH4]
  []
[]

[AuxKernels]
  [nodal_void_volume_auxk]
    type = NodalVoidVolumeAux
    variable = nodal_void_volume
    nodal_void_volume_uo = nodal_void_volume_uo
    execute_on = 'initial timestep_end' # "initial" to ensure it is properly evaluated for the first timestep
  []
  [rate_H2O_per_1l]
    type = ParsedAux
    args = 'pf_rate_H2O nodal_void_volume'
    variable = rate_H2O_per_1l
# pf_rate = change in kg at every node
# pf_rate * 1000 / molar_mass_in_g_per_mole = change in moles at every node
# pf_rate * 1000 / molar_mass / (nodal_void_volume_in_m^3 * 1000) = change in moles per litre of aqueous solution
    function = 'pf_rate_H2O / 18.0152 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_Ca_per_1l]
    type = ParsedAux
    args = 'pf_rate_Ca nodal_void_volume'
    variable = rate_Ca_per_1l
    function = 'pf_rate_Ca / 40.08 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_SO4_per_1l]
    type = ParsedAux
    args = 'pf_rate_SO4 nodal_void_volume'
    variable = rate_SO4_per_1l
    function = 'pf_rate_SO4 / 96.0576 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_CH3COO_per_1l]
    type = ParsedAux
    args = 'pf_rate_CH3COO nodal_void_volume'
    variable = rate_CH3COO_per_1l
    function = 'pf_rate_CH3COO / 59.0445 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_HS_per_1l]
    type = ParsedAux
    args = 'pf_rate_HS nodal_void_volume'
    variable = rate_HS_per_1l
    function = 'pf_rate_HS / 33.0679 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [rate_CH4_per_1l]
    type = ParsedAux
    args = 'pf_rate_CH4 nodal_void_volume'
    variable = rate_CH4_per_1l
    function = 'pf_rate_CH4 / 16.0426 / nodal_void_volume'
    execute_on = 'timestep_begin'
  []
  [transported_H2O]
    type = GeochemistryQuantityAux
    variable = transported_H2O
    species = H2O
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_end'
  []
  [transported_Ca]
    type = GeochemistryQuantityAux
    variable = transported_Ca
    species = "Ca++"
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_end'
  []
  [transported_SO4]
    type = GeochemistryQuantityAux
    variable = transported_SO4
    species = "SO4--"
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_end'
  []
  [transported_CH3COO]
    type = GeochemistryQuantityAux
    variable = transported_CH3COO
    species = "CH3COO-"
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_end'
  []
  [transported_HS]
    type = GeochemistryQuantityAux
    variable = transported_HS
    species = "HS-"
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_end'
  []
  [transported_CH4]
    type = GeochemistryQuantityAux
    variable = transported_CH4
    species = "CH4(aq)"
    quantity = transported_moles_in_original_basis
    execute_on = 'timestep_end'
  []
  [transported_mass]
    type = ParsedAux
    variable = transported_mass
    args = 'transported_H2O transported_Ca transported_SO4 transported_CH3COO transported_HS transported_CH4'
    function = 'transported_H2O * 18.0152 + transported_Ca * 40.08 + transported_SO4 * 96.0576 +  transported_CH3COO * 59.0445 +  transported_HS * 33.0679 + transported_CH4 * 16.0426'
    execute_on = 'timestep_end'
  []
  [massfrac_H2O]
    type = ParsedAux
    args = 'transported_H2O transported_mass'
    variable = massfrac_H2O
    function = 'transported_H2O * 18.0152 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_Ca]
    type = ParsedAux
    args = 'transported_Ca transported_mass'
    variable = massfrac_Ca
    function = 'transported_Ca * 40.08 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_SO4]
    type = ParsedAux
    args = 'transported_SO4 transported_mass'
    variable = massfrac_SO4
    function = 'transported_SO4 * 96.0576 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_CH3COO]
    type = ParsedAux
    args = 'transported_CH3COO transported_mass'
    variable = massfrac_CH3COO
    function = 'transported_CH3COO * 59.0445 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_HS]
    type = ParsedAux
    args = 'transported_HS transported_mass'
    variable = massfrac_HS
    function = 'transported_HS * 33.0679 / transported_mass'
    execute_on = 'timestep_end'
  []
  [massfrac_CH4]
    type = ParsedAux
    args = 'transported_CH4 transported_mass'
    variable = massfrac_CH4
    function = 'transported_CH4 * 16.0426 / transported_mass'
    execute_on = 'timestep_end'
  []
[]

[Postprocessors]
  [massfrac_H2O]
    type = PointValue
    point = '0 0 0'
    variable = massfrac_H2O
  []
  [massfrac_Ca]
    type = PointValue
    point = '0 0 0'
    variable = massfrac_Ca
  []
  [massfrac_SO4]
    type = PointValue
    point = '0 0 0'
    variable = massfrac_SO4
  []
  [massfrac_CH3COO]
    type = PointValue
    point = '0 0 0'
    variable = massfrac_CH3COO
  []
  [massfrac_HS]
    type = PointValue
    point = '0 0 0'
    variable = massfrac_HS
  []
  [massfrac_CH4]
    type = PointValue
    point = '0 0 0'
    variable = massfrac_CH4
  []
[]

[Outputs]
  csv = true
  exodus = true
[]
