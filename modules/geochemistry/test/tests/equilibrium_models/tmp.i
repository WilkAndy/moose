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
    basis_species = "H2O H+ Cl- Al+++"
  [../]
  [./solver]
    type = EquilibriumReactionSolver
    model_definition = definition
  verbose = true
    charge_balance_species = "Cl-"
    mass_solvent_water = 1
  rel_tol = 1E-200
  abs_tol = 1E-10

  # this is difficult because one huge aqueous complex gets an enormous molality because of [H+]^-32 
   # bulk_moles_species = "H+ Al+++   Cl-"
   # bulk_moles_values = " 1E-15  2.5945E-6 5.4E-5"

  # this does not work because one huge aqueous complex gets an enormous molality because of [H+]^-32
# and then it means the bulk conc of H+ is huge, which means the bulk of Cl- must be super negative to ensure charge balance, which means the free molality of Cl- must be as close to zero as possible, which then prevents anything happening due to the damping
  free_molality_species = "H+     "
  free_molality_values = "3.162E-7"
    bulk_moles_species = "Al+++   Cl-"
    bulk_moles_values = " 2.5945E-6 5.4E-5"

  max_iter = 100
  [../]  
[]

