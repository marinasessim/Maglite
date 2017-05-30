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
     recoil_energy_limits = 1
     neutron_energy_limits = 1

     # Functions
     neutron_spectrum = spectrum_simple
     epithermal_law = epithermal_scattering_law
     fast_law = fast_scattering_law
     simple_erxs = elastic_cross_section

     execute_on = timestep_end
  [../]
[]

# t is equal to Ei
[Functions]
  [./spectrum_simple]
     type = ParsedFunction
     value = '1 / t'
  [../]

  # t is equal to Ei
  [./elastic_cross_section]
     type = ParsedFunction
     value = '1 + 0.0000001 * t'
  [../]

  # t is equal to mu_c, but it is always 1
  [./epithermal_scattering_law]
      type = ParsedFunction
     value = '1'
  [../]

  # t is equal to mu_c
  [./fast_scattering_law]
      type = ParsedFunction
     value = '1 - 0.5 * t'
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
