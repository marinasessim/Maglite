/******************************************************************************/
/*                          Elastic Recoil Cross Section                      */
/******************************************************************************/

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
    Real E_u = _neutron_energy_limits[g];
    Real E_l = _neutron_energy_limits[g + 1];

    for (unsigned int p = 0; p < _quad_points.size(); ++p)
    {
      Real E = 0.5 * (E_u - E_l) * _quad_points[p] + 0.5 * (E_u + E_l);
      Real w_E = 0.5 * _quad_weights[p] * (E_u - E_l);

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

  for (unsigned int l = 0; l < _L; ++l)
  {
    for (unsigned int g = 0; g < _G; ++g)
    {
      Real E_u = _neutron_energy_limits[g];
      Real E_l = _neutron_energy_limits[g + 1];

      for (unsigned int i_E = 0; i_E < _quad_points.size(); ++i_E)
      {
        Real E = 0.5 * (E_u - E_l) * _quad_points[i_E] + 0.5 * (E_u + E_l);
        Real w_E = 0.5 * _quad_weights[i_E] * (E_u - E_l);

        Real T_max = _gamma * E;
        for (unsigned int t = 0; t < _T; ++t)
        {
          Real T_u = _recoil_energy_limits[t];
          Real T_l = _recoil_energy_limits[t + 1];

          /// Case III: T_max < interval, stop
          if (T_max < T_l)
            continue; // Breaks for loop

          /// Case II: T_max inside interval, refit
          if (T_max > T_l && T_max < T_u)
            T_u = T_max; // Refit to new interval between T_l and T_max

          /// Case I: T_max > interval, ok
          /// Case II also utilizes this piece of code with T_u = T_max
          for (unsigned int i_T = 0; i_T < _quad_points.size(); ++i_T)
          {
            Real T = 0.5 * (T_u - T_l) * _quad_points[i_T] + 0.5 * (T_u + T_l);
            Real w_T = 0.5 * _quad_weights[i_T] * (T_u - T_l);

            /// Calculate cosine of recoil angle in CM frame
            Real mu_c = 1 - 2 * T / (E * _gamma);

            /// Calculate cosine of recoil angle in Lab frame
            Real mu_L =  sqrt((1 - mu_c) / 2);

            /// Calculate contribution to cross section coefficients
            _erxs_coeff[l][t][g] += 2 * libMesh::pi / _xi_g[g] *
                                    _elastic_xs.value(E,Point()) *
                                    _neutron_spectrum.value(E, Point()) *
                                    _scattering_law.value(mu_c, Point()) *
                                    legendreP(l, mu_L) * w_T * w_E;
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

    if (l < _L - 1)
      output_file << std::endl << std::endl;
  }
  output_file.close();
}
