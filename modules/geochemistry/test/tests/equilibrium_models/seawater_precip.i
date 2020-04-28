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
    basis_species = "H2O H+ Cl- Na+ SO4-- Mg++ Ca++ K+ HCO3- SiO2(aq) O2(aq)"
    equilibrium_minerals = "Antigorite Tremolite Talc Chrysotile Sepiolite Anthophyllite Dolomite Dolomite-ord Huntite Dolomite-dis Magnesite Calcite Aragonite Quartz"
  
  [../]
  [./solver]
    type = EquilibriumReactionSolver
    model_definition = definition
    swap_out_of_basis = "H+   "
    swap_into_basis = "  MgCO3"
    charge_balance_species = "Cl-" # this means the bulk moles of Cl- will not be exactly as set below
    mass_solvent_water = 1
# Bethke specifies TDS = 35080mg/kg.
# Note that TDS = (mass_non-water) / (mass_solvent_water + mass_non-water),
# so with mass_solvent_water = 1kg, mass_non-water = 0.036355346 kg, and total_mass = 1.036355346 kg
# concentration of Cl- = 19350mg/kg -> moles = 19350E-3 * 1.036355346 / 35.4530 = 0.565635516
  free_molality_species = "MgCO3    O2(aq)"
  free_molality_values = "0.1068E-3 0.2151E-3"
    #bulk_moles_species = "Cl-        Na+         SO4--       Mg++        Ca++        K+          HCO3-        SiO2(aq)"
    #bulk_moles_values = "0.565635516 0.485049175 0.029237905 0.055005077 0.010627297 0.010576055 0.0024118232 0.00010349"
    bulk_moles_species = "Cl-   Na+    SO4--   Mg++    Ca++    K+          HCO3-    SiO2(aq)"
    bulk_moles_values = "0.5656 0.4850 0.02924 0.05501 0.01063 0.010576055 0.002412 0.00010349"
    #prevent_precipitation = "Tremolite Talc Chrysotile Sepiolite Anthophyllite Dolomite-ord Huntite Dolomite-dis Magnesite Calcite Aragonite"
    prevent_precipitation = "Dolomite-dis Dolomite-ord"
  verbose = true
  [../]  
[]

