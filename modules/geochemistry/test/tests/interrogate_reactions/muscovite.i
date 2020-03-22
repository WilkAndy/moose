[GeochemicalModelInterrogator]
  model_root = root
  swap_out_of_basis = "Al+++"
  swap_into_basis = "  Kaolinite"
  activity_species = "H2O"
  activity_values = "1.0"
  equilibrium_species = Muscovite
  interrogation = activity
[]

[UserObjects]
  [./root]
    type = GeochemicalModelRoot
    database_file = "../../data/entire_gwb.txt"
    basis_species = "H2O K+ Al+++ SiO2(aq) H+"
    equilibrium_minerals = "Muscovite Kaolinite"
  [../]
[]

