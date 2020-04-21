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
    type = ERS1
    model_definition = definition
  swap_out_of_basis = "Al+++ Fe++"
  swap_into_basis = "Kaolinite Hematite"
  verbose = true
    charge_balance_species = "Cl-"
# Approximately TDS = 44mg/kg
# Note that TDS = (mass_non-water) / (mass_solvent_water + mass_non-water),
# so with mass_solvent_water = 1kg, mass_non-water = 4.4E-5kg kg, and total_mass = 1.000044kg
# concentration of SiO2(aq) = 7mg/kg -> moles = 7E-3 * 1.000044 / 60.0843 = 0.0001165
    constraint_species = "H2O              H+                   O2(aq)        SiO2(aq)           Kaolinite                  Hematite               Ca++               Mg++               Na+                HCO3-              SO4--              Cl-      "
    constraint_value = "  1.0              3.162E-7             1.813E-4      0.0001165          1E-2                       0.03                       0.0001073          4.526E-5           7.83E-5            0.0003114          3.1233E-5          5.3594E-5"
    constraint_meaning = "kg_solvent_water activity             free_molality moles_bulk_species free_moles_mineral_species free_moles_mineral_species moles_bulk_species moles_bulk_species moles_bulk_species moles_bulk_species moles_bulk_species moles_bulk_species"
  max_iter = 1E3
  [../]  
[]

