# t is equal to Ei
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
  [./dummy]
  [../]
[]

[UserObjects]
  [./XSGenerator]
    # Inputs
    # Functions
    type = ElasticRecoilCrossSectionUserObject
    isotope_type = '1'
    atomic_mass = '100'
    recoil_energy_limits = '1'
    neutron_energy_limits = '1'
    neutron_spectrum = spectrum_simple
    epithermal_law = epithermal_scattering_law
    fast_law = fast_scattering_law
    simple_erxs = elastic_cross_section
    execute_on = 'timestep_end'
  [../]
[]

[Functions]
  # t is equal to Ei
  # t is equal to mu_c, but it is always 1
  # t is equal to mu_c
  [./spectrum_simple]
    type = ParsedFunction
    function = '1 / t'
  [../]
  [./elastic_cross_section]
    type = ParsedFunction
    function = '1 + 0.0000001 * t'
  [../]
  [./epithermal_scattering_law]
    type = ParsedFunction
    function = '1'
  [../]
  [./fast_scattering_law]
    type = ParsedFunction
    function = '1 - 0.5 * t'
  [../]
[]

[Postprocessors]
  [./test]
    type = AverageNodalVariableValue
    variable = 'dummy'
  [../]
[]

[Executioner]
  type = Steady
[]

