[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 0
  []
[]

[Variables]
  [./u]
  [../]
[]

[Kernels]
  [./zero]
    type = CoefTimeDerivative
    variable = u
    Coefficient = 0.0
  [../]
[]

[Executioner]
  type = Transient
  end_time = 1
  dt = 1
[]
  
[Outputs]
  csv = true
[]
