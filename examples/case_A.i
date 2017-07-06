#                             Case A                                   #
# neutron spectrum = 1; cross section = 1; A = 1; scattering law = 1/2 #
#                                                                      #

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
     output_file_name = case_A_out.csv
     atomic_mass = 1
     legendre_order = 4
     neutron_energy_limits = '1e7 1e6 1e5 1e4 1e3 1e2 1e1 1e0'
     recoil_energy_limits = '2097152 1048576 524288 262144 131072 65536 32768 16384 8192 4096 2048 1024 512 256 128 64 32 16 8 4 2 1 0'

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
