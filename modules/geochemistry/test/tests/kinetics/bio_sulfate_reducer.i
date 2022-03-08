# Example of a microbe-catalysed reaction:
# CH3COO- + SO4-- -> 2HCO3- + HS-
# at pH = 7.2
# at temperature = 25degC
# with log10(K) = 8.404 for this reaction (defined in the database: see bio_sulfate_reducer0.i)
[TimeDependentReactionSolver]
  model_definition = definition
  geochemistry_reactor_name = reactor
  charge_balance_species = "Cl-"
  constraint_species = "H2O              Na+              Ca++             Fe++             Cl-              SO4--            HCO3-            CH3COO-          HS-              H+" 
  constraint_value = "  1.0              501E-3           20E-3            2E-3             500E-3           20E-3            2E-3             1E-3             0.3E-6           -7.2"
  constraint_meaning = "kg_solvent_water bulk_composition bulk_composition bulk_composition bulk_composition bulk_composition bulk_composition bulk_composition free_concentration log10activity"
  constraint_unit = "   kg               moles            moles            moles            moles            moles            moles            moles            molal              dimensionless"
  controlled_activity_name = 'H+'
  controlled_activity_value = 6.30957E-8 # this is pH=7.2
  kinetic_species_name = "sulfate_reducer"
  kinetic_species_initial_value = 0.1 # molecular weight of sulfate_reducer = 1, so this is the amount of mmoles too
  kinetic_species_unit = mg
  ramp_max_ionic_strength_initial = 0
  stoichiometric_ionic_str_using_Cl_only = true # for comparison with GWB
  execute_console_output_on = 'INITIAL TIMESTEP_END'
  mol_cutoff = 1E-20
  solver_info = true
  evaluate_kinetic_rates_always = true
  precision = 16
  prevent_precipitation = 'Pyrite Troilite'
[]

[UserObjects]
  [rate_sulfate_reducer]
    type = GeochemistryKineticRate
    kinetic_species_name = "sulfate_reducer"
    intrinsic_rate_constant = 0.0864 # 1E-9 mol/mg/s = 0.0864 mol/g/day
    multiply_by_mass = true
    promoting_species_names = 'CH3COO-'
    promoting_indices = '1'
    promoting_monod_indices = '1'
    promoting_half_saturation = 70E-6
    direction = dissolution
    biological_efficiency = 4.3
    energy_captured = 45E3
    theta = 0.2
    eta = 1
  []
  [definition]
    type = GeochemicalModelDefinition
    database_file = "db.json"
    basis_species = "H2O              Na+              Ca++             Fe++             Cl-              SO4--            HCO3-            CH3COO-          HS-              H+"
    equilibrium_minerals = "*" #Mackinawite Siderite
    kinetic_redox = "sulfate_reducer"
    kinetic_rate_descriptions = "rate_sulfate_reducer"
  []
[]

[Functions]
  [timestepper]
    type = PiecewiseLinear
    x = '0 10 18  21'
    y = '1E-2 1E-2  1   1'
  []
[]

[Executioner]
  type = Transient
  [TimeStepper]
    type = FunctionDT
    function = timestepper
  []
  end_time = 7
[]

[AuxVariables]
  [diss]
  []
  [moles_acetate]
  []
[]
[AuxKernels]
#  [moles_acetate]
#    type = GeochemistryQuantityAux
#    species = 'CH3COO-'
#    reactor = reactor
#    variable = moles_acetate
#    quantity = kinetic_moles
#  []
  [diss]
    type = ParsedAux
    args = "moles_acetate"
    function = '1.0E-10 - moles_acetate'
    variable = diss
  []
[]
[Postprocessors]
  [dissolved_moles]
    type = PointValue
    point = '0 0 0'
    variable = diss
  []
[]
[Outputs]
  csv = true
[]
