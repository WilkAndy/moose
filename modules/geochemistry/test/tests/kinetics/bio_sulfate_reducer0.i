# The purpose of this is to get the equilibrium constant 8.404
# CH3COO- = -1*SO4-- + 1*HS- + 2*HCO3-  .  log10(K) = 8.404
[GeochemicalModelInterrogator]
  model_definition = definition
  swap_into_basis = 'HS-'
  swap_out_of_basis = 'H+'
[]

[UserObjects]
  [definition]
    type = GeochemicalModelDefinition
    database_file = "../../../database/moose_geochemdb.json"
    basis_species = "H2O SO4-- H+ O2(aq) HCO3-"
    piecewise_linear_interpolation = true
  []
[]
