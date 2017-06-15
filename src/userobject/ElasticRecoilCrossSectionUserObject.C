/****************************************************************/
/*            Elastic Recoil Cross Section                      */
/****************************************************************/

// MOOSE includes
#include "ElasticRecoilCrossSectionUserObject.h"
#include "Function.h"

// libmesh includes
#include "libmesh/quadrature.h"

// Output includes
#include <cassert>
#include <fstream>
#include <stdexcept>

template <>
InputParameters
validParams<ElasticRecoilCrossSectionUserObject>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addParam<Real>(
      "atomic_mass", 1, "Atomic Mass of the isotope. Defaults to Hydrogen A = 1");
  params.addRequiredParam<std::vector<Real>>(
      "neutron_energy_limits", "Energy limits of the incident neutron in [eV] and descending order");
  params.addRequiredParam<std::vector<Real>>(
      "recoil_energy_limits", "Energy limits of the recoil atom in [eV] and descending order");
  params.addRequiredParam<FunctionName>(
      "neutron_spectrum","Function representing the reactor neutron spectrum");
  params.addRequiredParam<FunctionName>(
      "scattering_law", "Function representing the scattering law for neutrons");
  params.addRequiredParam<FunctionName>(
      "elastic_xs", "Function representing the elastic cross section");
  params.addParam<unsigned int>(
      "legendre_order", 6, "Order of Legendre polynomials; Default to P5, where n = 0, ..., 5");
  return params;
}

ElasticRecoilCrossSectionUserObject::ElasticRecoilCrossSectionUserObject(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _neutron_spectrum(getFunction("neutron_spectrum")),
    _scattering_law(getFunction("scattering_law")),
    _elastic_xs(getFunction("elastic_xs")),
    _L(getParam<unsigned int>("legendre_order")),
    _atomic_mass(getParam<Real>("atomic_mass")),
    _neutron_energy_limits(getParam<std::vector<Real>>("neutron_energy_limits")),
    _recoil_energy_limits(getParam<std::vector<Real>>("recoil_energy_limits"))
{
   /// See "enum_order.h and enum_quadrature_type.h"
  _quadrature = QBase::build(QGAUSS, 1, FORTYTHIRD).release();

  /// Extract 1 dim quadrature rule from libmesh and store in points and weights
  const std::vector<Point> & temp_points = _quadrature->get_points();
  for (auto & pt : temp_points)
  {
    _quad_points.push_back(pt(0));
  }
  _quad_weights = _quadrature->get_weights();

  _alpha = pow(((_atomic_mass - 1)/(_atomic_mass + 1)),2);
  _gamma = 4 * _atomic_mass / std::pow(( _atomic_mass + 1),2);
  _G = _neutron_energy_limits.size() - 1;
  _T = _recoil_energy_limits.size() - 1;
}

Real
ElasticRecoilCrossSectionUserObject::legendreP(unsigned int  n, Real x)
{
  switch (n)
  {
    case 0: return 1;
    case 1: return x;
    case 2: return 0.5 * (3 * pow(x, 2) - 1);
    case 3: return  0.5 * (5 * pow(x, 3) - 3 * x);
    case 4: return 0.125 * (35 * pow(x, 4) - 30 * pow(x, 2) + 3);
    case 5: return 0.125 * (63 * pow(x, 5) - 70 * pow(x, 3) + 15 * x);
    default:
      mooseError("Implementation of Legendre polynomials goes up to n = 5");
  }
}

ElasticRecoilCrossSectionUserObject::~ElasticRecoilCrossSectionUserObject()
{
  delete _quadrature;
}

void
ElasticRecoilCrossSectionUserObject::initialize()
{
  /// Calculate neutron spectrum over group g
  _xi_g.resize(_G);
  for (unsigned int g = 0; g < _G; ++g)
  {
    Real E_max = _neutron_energy_limits[g];
    Real E_min = _neutron_energy_limits[g + 1];

    for (unsigned int p = 0; p < _quad_points.size(); ++p)
    {
      Real E = 0.5 * (E_max - E_min) * _quad_points[p] + 0.5 * (E_max + E_min);
      Real w_E = 0.5 * _quad_weights[p] * (E_max - E_min);

      _xi_g[g] += w_E * _neutron_spectrum.value(E, Point());
    }
  }
}

/// Find the neutron energy group given a neutron energy
unsigned int
ElasticRecoilCrossSectionUserObject::findNeutronEnergyGroup(Real energy)
{
  for (unsigned int g = 0; g < _G; ++g)
  {
    if (energy < _neutron_energy_limits[g] && energy > _neutron_energy_limits[g+1])
      return g;
  }
  /// C++ is stupid so this is necessary
  mooseError("Should never get here");
  return 0;
}

void
ElasticRecoilCrossSectionUserObject::execute()
{
  /// Size the cross section array
    _erxs_coeff.resize(_L);
    for (unsigned int l = 0; l < _L; ++l)
    {
      _erxs_coeff[l].resize(_T);
      for (unsigned int t = 0; t < _T; ++t)
      {
        _erxs_coeff[l][t].resize(_G);
      }
    }

    /// Loop over all recoil energy bins
    for (unsigned int l = 0; l < _L; ++l)
    {
      for (unsigned int t = 0; t < _T; ++t)
      {
        Real T_max = _recoil_energy_limits[t];
        Real T_min = _recoil_energy_limits[t + 1];

        /// Loop over quadrature points within the current recoil bin [T_min, T_max]
        for (unsigned int i_T = 0; i_T < _quad_points.size(); ++i_T)
        {
          Real T_i = 0.5 * (T_max - T_min) * _quad_points[i_T] + 0.5 * (T_max + T_min);
          Real w_T = 0.5 * _quad_weights[i_T] * (T_max - T_min);

          /// Loop over quadrature points to integrate over theta_c [-pi, pi]
          for (unsigned int i_theta = 0; i_theta < _quad_points.size(); ++i_theta)
          {
            Real theta_c = libMesh::pi * _quad_points[i_theta];
            Real w_theta = _quad_weights[i_theta] * libMesh::pi;

            /// Calculate the neutron incident energy associated with the combination of T and theta_c
            Real E_i = 2 * T_i / (_gamma * (1 - std::cos(theta_c)));

            /// Find neutron energy group that E_i belongs to
            Real threshold = (T_i / (1 - _alpha));
            if ( E_i > threshold && E_i < _neutron_energy_limits[0] && E_i > _neutron_energy_limits[_G])
            {
              unsigned int g = findNeutronEnergyGroup(E_i);

              /// Save this contribution (g -> t) to the summation of the erxs
              _erxs_coeff[l][t][g] += libMesh::pi * _elastic_xs.value(E_i, Point()) * _neutron_spectrum.value(E_i, Point()) *
                                            std::cos(0.5 * theta_c) * _scattering_law.value(E_i, Point()) *
                                            legendreP(l, std::sin(0.5 * theta_c)) * w_T * w_theta / _xi_g[g];

            }
          }
        }
      }
    }
  }

void
ElasticRecoilCrossSectionUserObject::finalize()
{
  /**
   * Write the elastic recoil cross section (g -> t) output file
   * G (neutron energy groups) rows
   * T (recoil energy bins) columns
   */
  std::ofstream output_file;
  output_file.open ("erxs_output.csv");
  for (unsigned int l = 0; l < _L; ++l)
  {
    for (unsigned int g = 0; g < _G; ++g)
    {
      for (unsigned int t = 0; t < _T; ++t)
      {
        if (t < _T - 1)
        output_file << _erxs_coeff[l][t][g] << ',';
        else
        output_file << _erxs_coeff[l][t][g];
      }
      output_file << std::endl;
    }
    output_file << std::endl << std::endl;
  }
  output_file.close();
}
