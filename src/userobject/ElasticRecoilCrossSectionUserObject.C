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
      "legendre_order", 8, "Order of Legendre polynomials; Default to P5, where n = 0, ..., 10 fo LP and up to n = 8 for shifted LP");
  params.addParam<std::string>(
      "output_file_name", "erxs_output.csv", "Name of the output file (.csv)");
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
    _recoil_energy_limits(getParam<std::vector<Real>>("recoil_energy_limits")),
    _file_name(getParam<std::string>("output_file_name"))
{

  // See "enum_order.h and enum_quadrature_type.h"
  _quadrature = QBase::build(QGAUSS, 1, FORTYTHIRD).release();
  /*
  // Extract 1 dim quadrature rule from libmesh and store in points and weights
  const std::vector<Point> & temp_points = _quadrature->get_points();
  for (auto & pt : temp_points)
    _quad_points.push_back(pt(0));
  _quad_weights = _quadrature->get_weights();
  */

  _quad_points =   {-0.999714, -0.998492, -0.996295, -0.993125, -0.988984, -0.983878,
                    -0.977809, -0.970786, -0.962814, -0.953901, -0.944056, -0.933289,
                    -0.921609, -0.90903, -0.895562, -0.881219, -0.866015, -0.849965,
                    -0.833084, -0.815389, -0.796898, -0.777628, -0.757598, -0.736828,
                    -0.715338, -0.693149, -0.670283, -0.646762, -0.622609, -0.597847,
                    -0.572502, -0.546597, -0.520158, -0.493211, -0.465782, -0.437897,
                    -0.409585, -0.380873, -0.351789, -0.32236, -0.292617, -0.262588,
                    -0.232302, -0.20179, -0.17108, -0.140203, -0.109189, -0.0780686,
                    -0.0468717, -0.015629, 0.015629, 0.0468717, 0.0780686, 0.109189,
                    0.140203, 0.17108, 0.20179, 0.232302, 0.262588, 0.292617, 0.32236,
                    0.351789, 0.380873, 0.409585, 0.437897, 0.465782, 0.493211, 0.520158,
                    0.546597, 0.572502, 0.597847, 0.622609, 0.646762, 0.670283, 0.693149,
                    0.715338, 0.736828, 0.757598, 0.777628, 0.796898, 0.815389, 0.833084,
                    0.849965, 0.866015, 0.881219, 0.895562, 0.90903, 0.921609, 0.933289,
                    0.944056, 0.953901, 0.962814, 0.970786, 0.977809, 0.983878, 0.988984,
                    0.993125, 0.996295, 0.998492, 0.999714};

  _quad_weights =  {0.000734634, 0.00170939, 0.00268393, 0.00365596, 0.00462445,
                    0.00558843, 0.00654695, 0.00749907, 0.00844387, 0.00938042,
                    0.0103078, 0.0112251, 0.0121315, 0.0130259, 0.0139077, 0.0147759,
                    0.0156296, 0.0164681, 0.0172905, 0.0180959, 0.0188837, 0.0196531,
                    0.0204032, 0.0211334, 0.021843, 0.0225312, 0.0231974, 0.023841,
                    0.0244612, 0.0250575, 0.0256294, 0.0261762, 0.0266975, 0.0271926,
                    0.0276612, 0.0281028, 0.0285169, 0.0289031, 0.0292611, 0.0295905,
                    0.029891, 0.0301623, 0.0304041, 0.0306162, 0.0307984, 0.0309505,
                    0.0310723, 0.0311638, 0.0312249, 0.0312554, 0.0312554, 0.0312249,
                    0.0311638, 0.0310723, 0.0309505, 0.0307984, 0.0306162, 0.0304041,
                    0.0301623, 0.029891, 0.0295905, 0.0292611, 0.0289031, 0.0285169,
                    0.0281028, 0.0276612, 0.0271926, 0.0266975, 0.0261762, 0.0256294,
                    0.0250575, 0.0244612, 0.023841, 0.0231974, 0.0225312, 0.021843,
                    0.0211334, 0.0204032, 0.0196531, 0.0188837, 0.0180959, 0.0172905,
                    0.0164681, 0.0156296, 0.0147759, 0.0139077, 0.0130259, 0.0121315,
                    0.0112251, 0.0103078, 0.00938042, 0.00844387, 0.00749907, 0.00654695,
                    0.00558843, 0.00462445, 0.00365596, 0.00268393, 0.00170939,
                    0.000734634};

  _alpha = pow(((_atomic_mass - 1) / (_atomic_mass + 1)), 2);
  _gamma = 4 * _atomic_mass / std::pow(( _atomic_mass + 1), 2);
  _G = _neutron_energy_limits.size() - 1;
  _T = _recoil_energy_limits.size() - 1;
}

Real
ElasticRecoilCrossSectionUserObject::legendreP(unsigned int n, Real x)
{
  switch (n)
  {
    case 0: return 1;
    case 1: return x;
    case 2: return 0.5 * (3 * pow(x, 2) - 1);
    case 3: return 0.5 * (5 * pow(x, 3) - 3 * x);
    case 4: return 0.125 * (35 * pow(x, 4) - 30 * pow(x, 2) + 3);
    case 5: return 0.125 * (63 * pow(x, 5) - 70 * pow(x, 3) + 15 * x);
    case 6: return 0.0625 * (231 * pow(x, 6) - 315 * pow(x, 4) + 105 * pow(x, 2) - 5);
    case 7: return 0.0625 * (429 * pow(x, 7) - 693 * pow(x, 5) + 315 * pow(x, 3) - 35 * x);
    case 8: return 0.0625 * (6435 * pow(x,8) - 12012 * pow(x,6) + 6930 * pow(x, 4) - 1260 * pow(x, 2) + 35);
    case 9: return 0.0078125 * (12155 * pow(x, 9) - 25740 * pow(x, 7) + 18018 * pow(x,5) - 4620 * pow(x, 3) + 315 * x);
    case 10: return 0.00390625 * (46189 * pow(x, 10) - 109395 * pow(x, 8) + 90090 * pow(x, 6) - 30030 * pow(x, 4) + 3465 * pow(x, 2) - 63);
    default:
      mooseError("Implementation of Legendre polynomials goes up to n = 10");
  }
}

Real
ElasticRecoilCrossSectionUserObject::shiftedLP(unsigned int n, Real x)
{
  switch (n)
  {
    case 0: return 1;
    case 1: return 2 * x - 1;
    case 2: return 6 * pow(x, 2) - 6 * x + 1;
    case 3: return 20 * pow(x, 3) - 30 * pow(x, 2) + 12 * x - 1;
    case 4: return 70 * pow(x, 4) - 140 * pow(x, 3) + 90 * pow(x, 2) - 20 * x + 1;
    case 5: return -1 + 30 * x - 210 * pow(x, 2) + 560 * pow(x, 3) - 630 * pow(x, 4) + 252 * pow(x, 5);
    case 6: return 1 - 42 * x + 420 * pow(x, 2) - 1680 * pow(x, 3) + 3150 * pow(x, 4) - 2772 * pow(x, 5) + 924 * pow(x, 6);
    case 7: return -1 + 56 * x - 756 * pow(x, 2) + 4200 * pow(x, 3) - 11550 * pow(x, 4) + 16632 * pow(x, 5) - 12012 * pow(x, 6) + 3432 * pow(x, 7);
    case 8: return 1 - 72 * x + 1260 * pow(x, 2) - 9240 * pow(x, 3) + 34650 * pow(x, 4) - 72072 * pow(x, 5) + 84084 * pow(x, 6) - 51480 * pow(x, 7) + 12870 * pow(x, 8);
    default:
      mooseError("Implementation of Shifted Legendre polynomials goes up to n = 8");
  }
}

ElasticRecoilCrossSectionUserObject::~ElasticRecoilCrossSectionUserObject()
{
  delete _quadrature;
}

void
ElasticRecoilCrossSectionUserObject::initialize()
{
  // Calculate neutron spectrum over group g
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

// Find the neutron energy group given a neutron energy
unsigned int
ElasticRecoilCrossSectionUserObject::findNeutronEnergyGroup(Real energy)
{
  for (unsigned int g = 0; g < _G; ++g)
  {
    if (energy < _neutron_energy_limits[g] && energy > _neutron_energy_limits[g + 1])
      return g;
  }
  /// C++ is stupid so this is necessary
  mooseError("Should never get here");
  return 0;
}

void
ElasticRecoilCrossSectionUserObject::execute()
{
  // Size the cross section array
  _erxs_coeff.resize(_L + 1);
  for (unsigned int l = 0; l < _L + 1; ++l)
  {
    _erxs_coeff[l].resize(_T);
    for (unsigned int t = 0; t < _T; ++t)
      _erxs_coeff[l][t].resize(_G);
  }

  _save_mu_L.resize(_T);
  for (unsigned int t = 0; t < _T; ++t)
    _save_mu_L[t].resize(_G);


  for (unsigned int l = 0; l < _L + 1; ++l)
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

Real mu_C_max = 1 - 2 * T_u / E_l;
if (mu_C_max < -1)
  mu_C_max = -1;
Real mu_C_min = 1 - 2 * T_l / E_u;
Real mu_L_max = sqrt((1 - mu_C_max) / 2);
Real mu_L_min = sqrt((1 - mu_C_min) / 2);

          // Case III: T_max < interval, stop
          if (T_max < T_l)
            continue; // Breaks for loop

          // Case II: T_max inside interval, refit
          if (T_max > T_l && T_max < T_u)
            T_u = T_max; // Refit to new interval between T_l and T_max

          // Case I: T_max > interval, ok
          // Case II also utilizes this piece of code with T_u = T_max
          for (unsigned int i_T = 0; i_T < _quad_points.size(); ++i_T)
          {
            Real T = 0.5 * (T_u - T_l) * _quad_points[i_T] + 0.5 * (T_u + T_l);
            Real w_T = 0.5 * _quad_weights[i_T] * (T_u - T_l);

            // Calculate cosine of recoil angle in CM frame
            Real mu_c = 1 - 2 * T / (E * _gamma);

            // Calculate cosine of recoil angle in Lab frame
            Real mu_L = sqrt((1 - mu_c) / 2);
if (mu_L * 1.1 > mu_L_max  || mu_L / 1.1 < mu_L_min)
  _console << "Violated min/max mu_L; E " << E << " T " << E << " mu_L " << mu_L << std::endl;
else
  _console << "good" << std::endl;
            _save_mu_L[t][g] =  mu_L;

            //_console << E << ", "  << T << ", " << mu_c << ", " << mu_L << std::endl;

            /// Calculate contribution to cross section coefficients
            _erxs_coeff[l][t][g] += 1 / _xi_g[g] *
                                    _elastic_xs.value(E,Point()) *
                                    _neutron_spectrum.value(E, Point()) *
                                    _scattering_law.value(mu_c, Point()) *
                                    shiftedLP(l, mu_L) * w_T * w_E;
          }
        }
      }
    }
  }
}

void
ElasticRecoilCrossSectionUserObject::finalize()
{
  /*
   * Write the elastic recoil cross section (g -> t) output file
   *
   *            l = 0              l = 1              l = ...
   *       g = 0, 1, 2, ...   g = 0, 1, 2, ...   g = 0, 1, 2, ...
   *  t 0
   *    1
   *    2
   *    .
   *    .
   *    .
   *
   */

  std::ofstream output_file;
  output_file.open(_file_name);
  for (unsigned int t = 0; t < _T; ++t)
  {
    for (unsigned int l = 0; l < _L + 1; ++l)
    {
      for (unsigned int g = 0; g < _G; ++g)
      {
        if (l < _L || g < _G - 1)
        output_file << _erxs_coeff[l][t][g] << ',';
        else
        output_file << _erxs_coeff[l][t][g];
      }
    }
    output_file << std::endl;
  }
  output_file.close();

  std::ofstream output_file2;
  output_file2.open("mu_L_out.csv");
  for (unsigned int t = 0; t < _T; ++t)
  {
      for (unsigned int g = 0; g < _G; ++g)
      {
        if (g < _G - 1)
        output_file2 << _save_mu_L[t][g] << ',';
        else
        output_file2 << _save_mu_L[t][g];
      }
    output_file2 << std::endl;
  }
  output_file2.close();
}
