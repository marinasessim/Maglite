#                               Case F                                      #
# neutron spectrum = 1/E; cross section = data; A = 1; scattering law = 1/2 #
#                                                                           #

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
     neutron_energy_limits = '2 1'
     recoil_energy_limits = '10.0 9.9 9.8 9.7 9.6 9.5 9.4 9.3 9.2 9.1 9.0 8.9 8.8 8.7 8.6 8.5 8.4 8.3 8.2 8.1 8.0 7.9 7.8 7.7 7.6 7.5 7.4 7.3 7.2 7.1 7.0 6.9 6.8 6.7 6.6 6.5 6.4 6.3 6.2 6.1 6.0 5.9 5.8 5.7 5.6 5.5 5.4 5.3 5.2 5.1 5.0 4.9 4.8 4.7 4.6 4.5 4.4 4.3 4.2 4.1 4.0 3.9 3.8 3.7 3.6 3.5 3.4 3.3 3.2 3.1 3.0 2.9 2.8 2.7 2.6 2.5 2.4 2.3 2.2 2.1 2.0 1.9 1.8 1.7 1.6 1.5 1.4 1.3 1.2 1.1 1.0 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0.0'

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
     value = '1/t'
  [../]

  # t is equal to Ei
  [./elastic_xs]
    type = PiecewiseLinear
    data_file = H1_input.csv
    xy_in_file_only = false
    format = columns
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
