[GeochemicalModelInterrogator]
  model_root = root
  swap_out_of_basis = "Al+++     SiO2(aq) H+"
  swap_into_basis = "  Muscovite Quartz   OH-"
  equilibrium_species = "Clinoptil-Ca"
[]

[UserObjects]
  [./root]
    type = GeochemicalModelRoot
    database_file = "../../data/entire_gwb.txt"
    basis_species = "H2O Ca++ Al+++ SiO2(aq) H+ K+"
    equilibrium_minerals = "Clinoptil-Ca Muscovite Quartz"
  [../]
[]

