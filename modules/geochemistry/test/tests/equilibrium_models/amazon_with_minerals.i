[Mesh]
type = GeneratedMesh
  dim = 1
[]

  [Variables]
  [./u]
  [../]
  []

  [Kernels]
    [./u]
  type = Diffusion
  variable = u
  [../]
  []

  [Executioner]
  type = Steady
  []


[UserObjects]
  [./definition]
    type = GeochemicalModelDefinition
    database_file = "../../../database/moose_geochemdb.json"
    basis_species = "H2O H+ SiO2(aq) Al+++ Fe++ Ca++ Mg++ Na+ HCO3- SO4-- Cl- O2(aq)"
    equilibrium_minerals = "Kaolinite Hematite"
  [../]
  [./solver]
    type = EquilibriumReactionSolver
    model_definition = definition
  swap_out_of_basis = "Al+++ Fe++"
  swap_into_basis = "Kaolinite Hematite"
  verbose = true
    charge_balance_species = "Cl-"
    mass_solvent_water = 1
  free_moles_mineral_species = "Kaolinite Hematite"
  free_moles_mineral_values = " 1E-2      0.03"
# Approximately TDS = 44mg/kg
# Note that TDS = (mass_non-water) / (mass_solvent_water + mass_non-water),
# so with mass_solvent_water = 1kg, mass_non-water = 4.4E-5kg kg, and total_mass = 1.000044kg
# concentration of SiO2(aq) = 7mg/kg -> moles = 7E-3 * 1.000044 / 60.0843 = 0.0001165
# easy:
  free_molality_species = "O2(aq)"
  free_molality_values = "1.813E-4" 
    bulk_moles_species = "H+     SiO2(aq)  Ca++      Mg++     Na+     HCO3-     SO4--     Cl-      "
    bulk_moles_values = "11.9E-5 0.0001165 0.0001073 4.526E-5 7.83E-5 0.0003114 3.1233E-5 5.3594E-5"

# hard:
#  activity_species = "H+"
#  activity_values = "1.0"
#  free_molality_species = "H+      O2(aq)"
#  free_molality_values = "3.162E-7 1.813E-4" 
#    bulk_moles_species = "SiO2(aq)  Ca++      Mg++     Na+     HCO3-     SO4--     Cl-      "
#    bulk_moles_values = " 0.0001165 0.0001073 4.526E-5 7.83E-5 0.0003114 3.1233E-5 5.3594E-5"
  max_iter = 1E3
  [../]  
[]

