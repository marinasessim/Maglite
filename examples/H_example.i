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
     legendre_order = 1
     recoil_energy_limits = '2 1.95 1.9 1.85 1.8 1.75 1.7 1.65 1.6 1.55 1.5 1.45 1.4 1.35 1.3
                             1.25 1.2 1.15 1.1 1.05 1.0 0.95 0.9 0.85 0.8 0.75 0.7 0.65 0.6 0.55 0.5 0.45
                             0.4 0.35 0.3 0.25 0.2 0.15 0.1 0.05 0.0' # 11 groups
     neutron_energy_limits = '2 1' # 5 groups

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
     value = '1'
  [../]

  # t is equal to Ei
  [./elastic_xs]
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
