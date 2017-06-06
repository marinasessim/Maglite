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
     isotope_type = 1
     atomic_mass = 100
     recoil_energy_limits = '1e5 1e4 1e3 1e2 1e1'
     neutron_energy_limits = '20e6 1e6 1e3 1e1'

     # Functions
     neutron_spectrum = neutron_spectrum
     scattering_law = scattering_law
     simple_erxs = elastic_recoil_xs

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
  [./elastic_recoil_xs]
     type = ParsedFunction
     value = '1'
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

[Debug]
[]
