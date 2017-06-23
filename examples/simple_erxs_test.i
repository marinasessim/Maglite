[Mesh]
 type = GeneratedMesh
 dim = 1
 xmin = 0
 xmax = 1
 nx = 2
[]

[Problem]
  type = FEProblem
  kernel_coverage_check = false
[]

[Variables]
  [./test_var]
  [../]
[]

[UserObjects]
  [./XSGenerator]
     type = ElasticRecoilCrossSectionUserObject

     # Inputs
     atomic_mass = 1
     recoil_energy_limits = '4882812.5 976562.5 195312.5 39062.5 7812.5 1562.5 312.5 62.5 12.5 2.5 0.5  0.1' # 11 groups
     neutron_energy_limits = '1e6 1e5 1e4 1e3 1e2 1e1' # 5 groups

     # Functions
     neutron_spectrum = neutron_spectrum
     scattering_law = scattering_law
     elastic_xs = elastic_xs

     execute_on = timestep_end
  [../]
[]

# t is equal to Ei
[Functions]
  [./neutron_spectrum]
     type = ParsedFunction
     value = '1 / t'
  [../]

  # t is equal to Ei
  [./elastic_xs]
     type = ConstantFunction
     value = '1e6'
  [../]

  # t is equal to mu_c
  [./scattering_law]
      type = ParsedFunction
     value = '0.5'
  [../]
[]

[Postprocessors]
  [./test]
    type = AverageNodalVariableValue
    variable = test_var
  [../]
[]

[Executioner]
  type = Steady
[]
