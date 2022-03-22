# Example of a microbe-catalysed reaction:
# Lactate + 2HAsO4-- + 2H2O -> CH3COO- + CO3-- + 2As(OH)4-
# at pH = 9.8
# at temperature = 20degC
# The equation in the database involving lactate is
# Lactate- + 3O2(aq) -> 2H+ + 3HCO3-
# with log10(K) = 231.4 at 20degC
[TimeDependentReactionSolver]
  model_definition = definition
  geochemistry_reactor_name = reactor
  swap_into_basis = 'CO3--'
  swap_out_of_basis = 'HCO3-'
  charge_balance_species = "Cl-"
  constraint_species = "H2O              Na+              CO3--            Lactate-         Cl-              AsO4---          CH3COO-          As(OH)4-         H+" 
  constraint_value = "  1.0              1448E-3          24E-3            10E-3            1500E-3          10E-3            1E-6             1E-6             -9.8"
  constraint_meaning = "kg_solvent_water bulk_composition bulk_composition bulk_composition bulk_composition bulk_composition bulk_composition bulk_composition log10activity"
  constraint_unit = "   kg               moles            moles            moles            moles            moles            moles            moles            dimensionless"
  controlled_activity_name = 'H+'
  controlled_activity_value = 1.58489E-10 # this is pH=9.8
  kinetic_species_name = "arsenate_reducer"
  kinetic_species_initial_value = 0.5 # molecular weight of arsenate_reducer = 1, so this is the amount of mmoles too
  kinetic_species_unit = mg
  ramp_max_ionic_strength_initial = 0
  stoichiometric_ionic_str_using_Cl_only = true # for comparison with GWB
  execute_console_output_on = 'INITIAL TIMESTEP_END'
  mol_cutoff = 1E-20
  solver_info = true
  evaluate_kinetic_rates_always = true
  precision = 16
[]

[UserObjects]
  [rate_arsenate_reducer]
    type = GeochemistryKineticRate
    kinetic_species_name = "arsenate_reducer"
    intrinsic_rate_constant = 0.6048 # 7E-9 mol/mg/s = 0.6048 mol/g/day
    promoting_species_names = 'HAsO4--'
    promoting_indices = '1'
    promoting_monod_indices = '1'
    promoting_half_saturation = 10E-6
    multiply_by_mass = true
    direction = dissolution
    kinetic_biological_efficiency = 5
    energy_captured = 125E3
    theta = 0.25
    eta = 1
  []
  [definition]
    type = GeochemicalModelDefinition
    database_file = "db_lactate.json"
    basis_species = "H2O Na+ Cl- HCO3- H+ As(OH)4- Lactate- CH3COO- AsO4---"
    kinetic_redox = "arsenate_reducer"
    kinetic_rate_descriptions = "rate_arsenate_reducer"
  []
[]

[Functions]
  [timestepper]
    type = PiecewiseLinear
    x = '0 10 18  21'
    y = '1E-1 1E-1  1   1'
  []
[]

[Executioner]
  type = Transient
  [TimeStepper]
    type = FunctionDT
    function = timestepper
  []
  end_time = 2
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
