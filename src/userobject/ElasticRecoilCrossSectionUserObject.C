/****************************************************************/
/*            Elastic Recoil Cross Section                      */
/****************************************************************/

// MOOSE includes
#include "ElasticRecoilCrossSectionUserObject.h"
#include "Function.h"

// libmesh includes
#include "libmesh/quadrature.h"

// Outputs
#include <cassert>
#include <fstream>
#include <stdexcept>

template <>
InputParameters
validParams<ElasticRecoilCrossSectionUserObject>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addRequiredParam<std::vector<Real>>(
      "atomic_mass", "Atomic Mass of the isotopes");
  params.addRequiredParam<std::vector<unsigned int>>(
      "isotope_type", "Type of the target isotopes");
  params.addRequiredParam<std::vector<Real>>(
      "neutron_energy_limits", "Energy limits of the incident neutron. [eV] and descending order");
  params.addRequiredParam<std::vector<Real>>(
      "recoil_energy_limits", "Energy limits of the recoil atom. [eV] and descending order");
  params.addRequiredParam<FunctionName>(
      "neutron_spectrum","Function representing the reactor neutron spectrum");
  params.addRequiredParam<FunctionName>(
      "scattering_law", "Function representing the scattering law for neutrons");
  params.addRequiredParam<FunctionName>(
      "simple_erxs", "Function that handles the cross section calculation");
  params.addParam<unsigned int>(
      "legendre_order", 5, "Order of Legendre polynomials; Default = 5");

  return params;
}

ElasticRecoilCrossSectionUserObject::ElasticRecoilCrossSectionUserObject(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _neutron_spectrum(getFunction("neutron_spectrum")),
    _scattering_law(getFunction("scattering_law")),
    _simple_erxs(getFunction("simple_erxs")),

    _atomic_mass(getParam<std::vector<Real>>("atomic_mass")),
    _isotope_type(getParam<std::vector<unsigned int>>("isotope_type")),
    _number_isotope(_isotope_type.size()),
    _neutron_energy_limits(getParam<std::vector<Real>>("neutron_energy_limits")),
    _recoil_energy_limits(getParam<std::vector<Real>>("recoil_energy_limits")),
    _L(getParam<unsigned int>("legendre_order"))
{
   /// See "enum_order.h and enum_quadrature_type.h"
  _quadrature = QBase::build(QGAUSS, 1, FORTYTHIRD).release();

  const std::vector<Point> & temp_points = _quadrature->get_points();
  for (auto & pt : temp_points)
  {
    _quad_points.push_back(pt(0));
    _console << pt(0) << std::endl;
  }
  _quad_weights = _quadrature->get_weights();

}

Real
ElasticRecoilCrossSectionUserObject::legendreP(unsigned int  n, Real x)
{
  switch (n)
  {
    case 0:
      return 1;

    case 1:
      return x;

    case 2:
      return 0.5 * (3 * pow(x, 2) - 1);

    case 3:
      return  0.5 * (5 * pow(x, 3) - 3 * x);

    case 4:
      return 0.125 * (35 * pow(x, 4) - 30 * pow(x, 2) + 3);

    case 5:
      return 0.125 * (63 * pow(x, 5) - 70 * pow(x, 3) + 15 * x);

    default:
      mooseError("Implementation of Legendre polynomials goes up to n = 5");
  }
}

/// Class destructor
ElasticRecoilCrossSectionUserObject::~ElasticRecoilCrossSectionUserObject()
{
  delete _quadrature;
}

void
ElasticRecoilCrossSectionUserObject::initialize()
{
  /// Loop to calculate neutron spectrum over group g xi(E)
  unsigned int G = _neutron_energy_limits.size() - 1;
  _xi_g.resize(G);
  for (unsigned int g = 0; g < G; ++g)
  {
    Real E_max = _neutron_energy_limits[g];
    Real E_min = _neutron_energy_limits[g + 1];

    unsigned int nqp = _quad_points.size();

    for (unsigned int p = 0; p < nqp; ++p)
    {
      Real E = 0.5 * (E_max - E_min) * _quad_points[p] + 0.5 * (E_max + E_min);
      Real w_E = 0.5 * _quad_weights[p] * (E_max - E_min);

      _xi_g[g] += w_E * _neutron_spectrum.value(E, Point());
    }
  }
}

/// Find the neutron energy group
unsigned int
ElasticRecoilCrossSectionUserObject::findNeutronEnergyGroup(Real energy)
{
  unsigned int G = _neutron_energy_limits.size() - 1;
  Real E_max = _neutron_energy_limits[0];
  unsigned int group;

  for (unsigned int g = 0; g < G; ++g)
  {
    if (energy > E_max)
    {
      group = 1;
    }
    else if (energy > _neutron_energy_limits[g+1] && energy < _neutron_energy_limits[g])
    {
      group = g + 1;
    }
    else
    {
      group = G;
    }
  }
  return group;
}

void
ElasticRecoilCrossSectionUserObject::execute()
{
  _gamma = 4 * _atomic_mass[0] / std::pow(( _atomic_mass[0] + 1),2);
  unsigned int T = _recoil_energy_limits.size() - 1;
  unsigned int G = _neutron_energy_limits.size() - 1;

  for (unsigned int g = 0; g < G; ++g)
  {
    // size the cross section array
    _elastic_recoil_xs.resize(_L);
    for (unsigned int l = 0; l < _L; ++l)
    {
      _elastic_recoil_xs[l].resize(T);
      for (unsigned int t = 0; t < T; ++t)
        _elastic_recoil_xs[l][t].resize(G);
    }

    // For all T bins
    // there should be a loop over l here but for now I set l = 0
    for (unsigned int l = 0; l < _L; ++l)
    {
      for (unsigned int t = 0; t < T; ++t)
      {
        Real T_max = _recoil_energy_limits[t];
        Real T_min = _recoil_energy_limits[t + 1];

        // loop over quadrature points within the recoil bin [T_min, T_max]
        unsigned int nqp = _quad_points.size();
        for (unsigned int i_T = 0; i_T < nqp; ++i_T)
        {
          Real T_i = 0.5 * (T_max - T_min) * _quad_points[i_T] + 0.5 * (T_max + T_min);
          Real w_T = 0.5 * _quad_weights[i_T] * (T_max - T_min);

          // loop over quadrature points to integrate over theta_c [-pi, pi]
          for (unsigned int i_theta = 0; i_theta < nqp; ++i_theta)
          {
            Real theta_c = libMesh::pi * _quad_points[i_theta];
            Real w_theta = _quad_weights[i_theta] * libMesh::pi;

            Real Ei = 2 * T_i / (_gamma * (1 - std::cos(theta_c)));

            unsigned int g = findNeutronEnergyGroup(Ei);

            _elastic_recoil_xs[l][t][g] += 0.5 * std::cos(0.5 * theta_c) * _scattering_law.value(Ei, Point()) *
                                          legendreP(l, std::sin(0.5 * theta_c)) * w_T * w_theta / _xi_g[g];

            _console << l << ',' << t << ',' << g << ':' << _elastic_recoil_xs[l][t][g] << std::endl;
          }
        }
      }
    }
  }
}

void
ElasticRecoilCrossSectionUserObject::finalize()
{
  std::ofstream output_file;
  output_file.open ("erxs_output.csv");

  unsigned int T = _recoil_energy_limits.size() - 1;
  unsigned int G = _neutron_energy_limits.size() - 1;

for (unsigned int l = 0; l < _L; ++l)
{
  for (unsigned int g = 0; g < G; ++g) // g rows
  {
    for (unsigned int t = 0; t < T; ++t) // t columns
    {
      output_file << _elastic_recoil_xs[0][t][g] << ',';
    }
      output_file << std::endl; // l matrices
  }
  output_file << std::endl << std::endl;
}
  output_file.close();
}
