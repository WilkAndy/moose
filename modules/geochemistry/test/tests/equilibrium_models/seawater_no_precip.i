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
    equilibrium_gases = "O2(g) CO2(g)"
  
  [../]
  [./solver]
    type = EquilibriumReactionSolver
    model_definition = definition
    swap_out_of_basis = "H+     O2(aq)"
    swap_into_basis = "  CO2(g) O2(g)"
    charge_balance_species = "Cl-" # this means the bulk moles of Cl- will not be exactly as set below
    activity_species = "CO2(g)       O2(g)" # fix fugacity
    activity_values = "0.0003162278  0.2"
    mass_solvent_water = 1
# Bethke specifies TDS = 35080mg/kg.
# Note that TDS = (mass_non-water) / (mass_solvent_water + mass_non-water),
# so with mass_solvent_water = 1kg, mass_non-water = 0.036355346 kg, and total_mass = 1.036355346 kg
# concentration of Cl- = 19350mg/kg -> moles = 19350E-3 * 1.036355346 / 35.4530 = 0.565635516
    #bulk_moles_species = "Cl-        Na+         SO4--       Mg++        Ca++        K+          HCO3-        SiO2(aq)"
    #bulk_moles_values = "0.565635516 0.485049175 0.029237905 0.055005077 0.010627297 0.010576055 0.0024118232 0.00010349"
    bulk_moles_species = "Cl-   Na+    SO4--   Mg++    Ca++    K+          HCO3-    SiO2(aq)"
    bulk_moles_values = "0.5656 0.4850 0.02924 0.05501 0.01063 0.010576055 0.002412 0.00010349"
# no minerals are initially in the basis, so cannot specify their free mole number or bulk mole number
# no minerals will be swapped into the basis because of the following prevent_precipitation, but if that were not the case and minerals were swapped into the basis, due to supersaturation, their free mole number will assumed to be initially zero
    prevent_precipitation = "Antigorite Tremolite Talc Chrysotile Sepiolite Anthophyllite Dolomite Dolomite-ord Huntite Dolomite-dis Magnesite Calcite Aragonite Quartz"
  verbose = true
  [../]  
[]

