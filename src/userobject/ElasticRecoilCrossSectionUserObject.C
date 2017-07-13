
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
  params.addClassDescription("Calculate the recoil atom elastic cross section given the isotope and energy groups."
      "It outputs the coefficients for the Legendre expansion of the cross section according to the order of Legendre"
      "polynomials l, neutron energy groups g and recoil energy bins t. It also outputs the maximum and mininum"
      "cosines of the recoil atom scatter angle in the Laboratory frame.");
  params.addParam<Real>(
      "atomic_mass", 1, "Atomic Mass of the isotope. Default to Hydrogen A = 1");
  params.addRequiredParam<std::vector<Real>>(
      "neutron_energy_limits", "Energy limits of the incident neutron in [eV] and descending order");
  params.addRequiredParam<std::vector<Real>>(
      "recoil_energy_limits", "Energy limits of the recoil atom in [eV] and descending order");
  params.addRequiredParam<FunctionName>(
      "neutron_spectrum","Function representing the reactor neutron spectrum");
  params.addRequiredParam<FunctionName>(
      "scattering_law", "Function representing the scattering law for neutrons");
  params.addRequiredParam<FunctionName>(
      "elastic_xs", "Function representing the neutron elastic cross section");
  params.addParam<unsigned int>(
      "legendre_order", 10, "Order of Legendre polynomials where n = 0, ..., 10. Default to P10");
  params.addParam<std::string>(
      "erxs_output_file_name", "erxs_output.csv", "Name of the output file with the cross section coefficients (.csv)");
  params.addParam<std::string>(
      "mu_L_output_file_name", "mu_L_out.csv", "Name of the output file with the mininum and maximum mu_L (.csv)");
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
    _erxs_file_name(getParam<std::string>("erxs_output_file_name")),
    _mu_L_file_name(getParam<std::string>("mu_L_output_file_name"))
{
  _quad_points =   {-0.999982, -0.999905, -0.999767, -0.999567, -0.999305, -0.998982,
                    -0.998598, -0.998152, -0.997645, -0.997076, -0.996446, -0.995755,
                    -0.995002, -0.994188, -0.993313, -0.992377, -0.99138, -0.990322,
                    -0.989203, -0.988023, -0.986782, -0.98548, -0.984118, -0.982696,
                    -0.981213, -0.979669, -0.978065, -0.976401, -0.974677, -0.972893,
                    -0.971049, -0.969146, -0.967182, -0.96516, -0.963077, -0.960936,
                    -0.958735, -0.956476, -0.954157, -0.95178, -0.949345, -0.946851,
                    -0.944298, -0.941688, -0.939019, -0.936293, -0.933509, -0.930668,
                    -0.92777, -0.924814, -0.921802, -0.918732, -0.915607, -0.912425,
                    -0.909186, -0.905892, -0.902542, -0.899137, -0.895676, -0.89216,
                    -0.88859, -0.884964, -0.881284, -0.87755, -0.873762, -0.86992,
                    -0.866025, -0.862076, -0.858074, -0.85402, -0.849913, -0.845753,
                    -0.841542, -0.837279, -0.832964, -0.828598, -0.824181, -0.819713,
                    -0.815195, -0.810627, -0.806009, -0.801341, -0.796624, -0.791858,
                    -0.787043, -0.78218, -0.777269, -0.77231, -0.767303, -0.762249,
                    -0.757148, -0.752001, -0.746807, -0.741568, -0.736282, -0.730952,
                    -0.725576, -0.720156, -0.714692, -0.709183, -0.703631, -0.698036,
                    -0.692397, -0.686716, -0.680993, -0.675228, -0.669422, -0.663574,
                    -0.657685, -0.651756, -0.645787, -0.639778, -0.633729, -0.627642,
                    -0.621516, -0.615352, -0.60915, -0.60291, -0.596634, -0.59032,
                    -0.583971, -0.577585, -0.571164, -0.564708, -0.558217, -0.551691,
                    -0.545132, -0.538539, -0.531913, -0.525254, -0.518563, -0.51184,
                    -0.505086, -0.4983, -0.491484, -0.484638, -0.477761, -0.470856,
                    -0.463921, -0.456958, -0.449967, -0.442948, -0.435902, -0.428829,
                    -0.421729, -0.414604, -0.407453, -0.400277, -0.393076, -0.385851,
                    -0.378603, -0.371331, -0.364036, -0.356719, -0.34938, -0.34202,
                    -0.334638, -0.327236, -0.319814, -0.312372, -0.304911, -0.297431,
                    -0.289933, -0.282416, -0.274883, -0.267333, -0.259766, -0.252183,
                    -0.244585, -0.236971, -0.229343, -0.221701, -0.214046, -0.206377,
                    -0.198695, -0.191001, -0.183296, -0.175579, -0.167851, -0.160113,
                    -0.152366, -0.144608, -0.136842, -0.129068, -0.121285, -0.113495,
                    -0.105698, -0.0978951, -0.0900857, -0.0822707, -0.0744507,
                    -0.0666261, -0.0587973, -0.050965, -0.0431296, -0.0352914,
                    -0.0274511, -0.0196092, -0.011766, -0.00392208, 0.00392208, 0.011766,
                    0.0196092, 0.0274511, 0.0352914, 0.0431296, 0.050965, 0.0587973,
                    0.0666261, 0.0744507, 0.0822707, 0.0900857, 0.0978951, 0.105698,
                    0.113495, 0.121285, 0.129068, 0.136842, 0.144608, 0.152366, 0.160113,
                    0.167851, 0.175579, 0.183296, 0.191001, 0.198695, 0.206377, 0.214046,
                    0.221701, 0.229343, 0.236971, 0.244585, 0.252183, 0.259766, 0.267333,
                    0.274883, 0.282416, 0.289933, 0.297431, 0.304911, 0.312372, 0.319814,
                    0.327236, 0.334638, 0.34202, 0.34938, 0.356719, 0.364036, 0.371331,
                    0.378603, 0.385851, 0.393076, 0.400277, 0.407453, 0.414604, 0.421729,
                    0.428829, 0.435902, 0.442948, 0.449967, 0.456958, 0.463921, 0.470856,
                    0.477761, 0.484638, 0.491484, 0.4983, 0.505086, 0.51184, 0.518563,
                    0.525254, 0.531913, 0.538539, 0.545132, 0.551691, 0.558217, 0.564708,
                    0.571164, 0.577585, 0.583971, 0.59032, 0.596634, 0.60291, 0.60915,
                    0.615352, 0.621516, 0.627642, 0.633729, 0.639778, 0.645787, 0.651756,
                    0.657685, 0.663574, 0.669422, 0.675228, 0.680993, 0.686716, 0.692397,
                    0.698036, 0.703631, 0.709183, 0.714692, 0.720156, 0.725576, 0.730952,
                    0.736282, 0.741568, 0.746807, 0.752001, 0.757148, 0.762249, 0.767303,
                    0.77231, 0.777269, 0.78218, 0.787043, 0.791858, 0.796624, 0.801341,
                    0.806009, 0.810627, 0.815195, 0.819713, 0.824181, 0.828598, 0.832964,
                    0.837279, 0.841542, 0.845753, 0.849913, 0.85402, 0.858074, 0.862076,
                    0.866025, 0.86992, 0.873762, 0.87755, 0.881284, 0.884964, 0.88859,
                    0.89216, 0.895676, 0.899137, 0.902542, 0.905892, 0.909186, 0.912425,
                    0.915607, 0.918732, 0.921802, 0.924814, 0.92777, 0.930668, 0.933509,
                    0.936293, 0.939019, 0.941688, 0.944298, 0.946851, 0.949345, 0.95178,
                    0.954157, 0.956476, 0.958735, 0.960936, 0.963077, 0.96516, 0.967182,
                    0.969146, 0.971049, 0.972893, 0.974677, 0.976401, 0.978065, 0.979669,
                    0.981213, 0.982696, 0.984118, 0.98548, 0.986782, 0.988023, 0.989203,
                    0.990322, 0.99138, 0.992377, 0.993313, 0.994188, 0.995002, 0.995755,
                    0.996446, 0.997076, 0.997645, 0.998152, 0.998598, 0.998982, 0.999305,
                    0.999567, 0.999767, 0.999905, 0.999982};

  _quad_weights =  {0.0000462637, 0.00010769, 0.000169201, 0.00023071, 0.000292206,
                    0.000353684, 0.000415141, 0.000476572, 0.000537974, 0.000599343,
                    0.000660675, 0.000721967, 0.000783214, 0.000844413, 0.00090556,
                    0.000966651, 0.00102768, 0.00108865, 0.00114955, 0.00121038,
                    0.00127114, 0.00133182, 0.00139241, 0.00145292, 0.00151334,
                    0.00157367, 0.0016339, 0.00169403, 0.00175406, 0.00181398,
                    0.00187379, 0.00193348, 0.00199305, 0.0020525, 0.00211182,
                    0.00217102, 0.00223008, 0.002289, 0.00234778, 0.00240642, 0.00246491,
                    0.00252325, 0.00258143, 0.00263945, 0.00269731, 0.00275501,
                    0.00281254, 0.00286989, 0.00292706, 0.00298406, 0.00304087,
                    0.0030975, 0.00315393, 0.00321017, 0.00326622, 0.00332206, 0.0033777,
                    0.00343313, 0.00348835, 0.00354335, 0.00359814, 0.0036527,
                    0.00370704, 0.00376115, 0.00381503, 0.00386868, 0.00392209,
                    0.00397525, 0.00402817, 0.00408085, 0.00413327, 0.00418544,
                    0.00423735, 0.004289, 0.00434039, 0.0043915, 0.00444235, 0.00449293,
                    0.00454323, 0.00459325, 0.00464298, 0.00469244, 0.0047416,
                    0.00479047, 0.00483904, 0.00488732, 0.0049353, 0.00498297,
                    0.00503034, 0.0050774, 0.00512414, 0.00517058, 0.00521669,
                    0.00526248, 0.00530794, 0.00535309, 0.0053979, 0.00544238,
                    0.00548652, 0.00553033, 0.00557379, 0.00561692, 0.00565969,
                    0.00570212, 0.0057442, 0.00578593, 0.0058273, 0.00586831, 0.00590896,
                    0.00594924, 0.00598916, 0.00602871, 0.00606789, 0.0061067,
                    0.00614513, 0.00618318, 0.00622086, 0.00625815, 0.00629505,
                    0.00633157, 0.0063677, 0.00640343, 0.00643878, 0.00647372,
                    0.00650827, 0.00654242, 0.00657616, 0.0066095, 0.00664243,
                    0.00667496, 0.00670707, 0.00673877, 0.00677006, 0.00680093,
                    0.00683138, 0.00686141, 0.00689102, 0.00692021, 0.00694896,
                    0.0069773, 0.0070052, 0.00703267, 0.00705971, 0.00708631, 0.00711248,
                    0.00713821, 0.0071635, 0.00718835, 0.00721276, 0.00723672,
                    0.00726024, 0.00728332, 0.00730594, 0.00732811, 0.00734984,
                    0.00737111, 0.00739193, 0.00741229, 0.0074322, 0.00745165,
                    0.00747064, 0.00748917, 0.00750724, 0.00752485, 0.00754199,
                    0.00755868, 0.00757489, 0.00759064, 0.00760593, 0.00762074,
                    0.00763509, 0.00764896, 0.00766237, 0.00767531, 0.00768777,
                    0.00769976, 0.00771127, 0.00772231, 0.00773288, 0.00774297,
                    0.00775258, 0.00776172, 0.00777038, 0.00777856, 0.00778626,
                    0.00779349, 0.00780023, 0.00780649, 0.00781228, 0.00781758,
                    0.0078224, 0.00782674, 0.0078306, 0.00783398, 0.00783687, 0.00783928,
                    0.00784121, 0.00784266, 0.00784363, 0.00784411, 0.00784411,
                    0.00784363, 0.00784266, 0.00784121, 0.00783928, 0.00783687,
                    0.00783398, 0.0078306, 0.00782674, 0.0078224, 0.00781758, 0.00781228,
                    0.00780649, 0.00780023, 0.00779349, 0.00778626, 0.00777856,
                    0.00777038, 0.00776172, 0.00775258, 0.00774297, 0.00773288,
                    0.00772231, 0.00771127, 0.00769976, 0.00768777, 0.00767531,
                    0.00766237, 0.00764896, 0.00763509, 0.00762074, 0.00760593,
                    0.00759064, 0.00757489, 0.00755868, 0.00754199, 0.00752485,
                    0.00750724, 0.00748917, 0.00747064, 0.00745165, 0.0074322,
                    0.00741229, 0.00739193, 0.00737111, 0.00734984, 0.00732811,
                    0.00730594, 0.00728332, 0.00726024, 0.00723672, 0.00721276,
                    0.00718835, 0.0071635, 0.00713821, 0.00711248, 0.00708631,
                    0.00705971, 0.00703267, 0.0070052, 0.0069773, 0.00694896, 0.00692021,
                    0.00689102, 0.00686141, 0.00683138, 0.00680093, 0.00677006,
                    0.00673877, 0.00670707, 0.00667496, 0.00664243, 0.0066095,
                    0.00657616, 0.00654242, 0.00650827, 0.00647372, 0.00643878,
                    0.00640343, 0.0063677, 0.00633157, 0.00629505, 0.00625815,
                    0.00622086, 0.00618318, 0.00614513, 0.0061067, 0.00606789,
                    0.00602871, 0.00598916, 0.00594924, 0.00590896, 0.00586831,
                    0.0058273, 0.00578593, 0.0057442, 0.00570212, 0.00565969, 0.00561692,
                    0.00557379, 0.00553033, 0.00548652, 0.00544238, 0.0053979,
                    0.00535309, 0.00530794, 0.00526248, 0.00521669, 0.00517058,
                    0.00512414, 0.0050774, 0.00503034, 0.00498297, 0.0049353, 0.00488732,
                    0.00483904, 0.00479047, 0.0047416, 0.00469244, 0.00464298,
                    0.00459325, 0.00454323, 0.00449293, 0.00444235, 0.0043915,
                    0.00434039, 0.004289, 0.00423735, 0.00418544, 0.00413327, 0.00408085,
                    0.00402817, 0.00397525, 0.00392209, 0.00386868, 0.00381503,
                    0.00376115, 0.00370704, 0.0036527, 0.00359814, 0.00354335,
                    0.00348835, 0.00343313, 0.0033777, 0.00332206, 0.00326622,
                    0.00321017, 0.00315393, 0.0030975, 0.00304087, 0.00298406,
                    0.00292706, 0.00286989, 0.00281254, 0.00275501, 0.00269731,
                    0.00263945, 0.00258143, 0.00252325, 0.00246491, 0.00240642,
                    0.00234778, 0.002289, 0.00223008, 0.00217102, 0.00211182, 0.0020525,
                    0.00199305, 0.00193348, 0.00187379, 0.00181398, 0.00175406,
                    0.00169403, 0.0016339, 0.00157367, 0.00151334, 0.00145292,
                    0.00139241, 0.00133182, 0.00127114, 0.00121038, 0.00114955,
                    0.00108865, 0.00102768, 0.000966651, 0.00090556, 0.000844413,
                    0.000783214, 0.000721967, 0.000660675, 0.000599343, 0.000537974,
                    0.000476572, 0.000415141, 0.000353684, 0.000292206, 0.00023071,
                    0.000169201, 0.00010769, 0.0000462637};

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

  // Size the lab frame cosine array
  _save_mu_L.resize(_T);
  for (unsigned int t = 0; t < _T; ++t)
  {
    _save_mu_L[t].resize(_G);
    for (unsigned int g = 0; g < _G; ++g)
      _save_mu_L[t][g].resize(1);
  }

  // Loop over all Legendre orders to calculate the coefficients
  for (unsigned int l = 0; l < _L + 1; ++l)
  {
    // Loop over all the neutron energy groups
    for (unsigned int g = 0; g < _G; ++g)
    {
      // Gets the upper and lower energy values of that group
      Real E_u = _neutron_energy_limits[g];
      Real E_l = _neutron_energy_limits[g + 1];

      // Loop over all the neutron energies within group g
      for (unsigned int i_E = 0; i_E < _quad_points.size(); ++i_E)
      {
        // Get the energies in Gaussian quadrature points and their respective weights
        Real E = 0.5 * (E_u - E_l) * _quad_points[i_E] + 0.5 * (E_u + E_l);
        Real w_E = 0.5 * _quad_weights[i_E] * (E_u - E_l);

        // Calculate maximum amount of energy transferred to the recoil atom
        Real T_max = _gamma * E;

        // Loop over all the possible recoil energy bins
        for (unsigned int t = 0; t < _T; ++t)
        {
          // Gets the upper and lower energy values of that bin
          Real T_u = _recoil_energy_limits[t];
          Real T_l = _recoil_energy_limits[t + 1];

          /*
           * Calculate possible range of angles according to neutron energy group
           * and recoil energy bin. This approach avoids unphysical behavior (negative
           * values due to convergence issue of expansion) of recoil cross section.
           * The subscript C means the Center of Mass (CM) frame and L, the Lab frame.
           */
          Real mu_C_max = 1 - 2 * T_u / E_l;
          if (mu_C_max < -1)
            mu_C_max = -1;
          Real mu_C_min = 1 - 2 * T_l / E_u;
          Real mu_L_max = sqrt((1 - mu_C_max) / 2);
          Real mu_L_min = sqrt((1 - mu_C_min) / 2);

          /*
           * Elastic scaterring case III: T_max < T_l, stop
           * When the maximum recoil energy is lower than the lowest energy of the recoil energy bin t
           */
          if (T_max < T_l)
            continue;

          /*
           * Elastic scaterring case II: T_max inside T bin, refit
           * When the maximum recoil energy is whitin the recoil energy bin t, we need to refit
           * the quadrature rule for the new interval between T_l and T_max
           */
          if (T_max > T_l && T_max < T_u)
            T_u = T_max;

          /*
           * Elastic scaterring case I: T_max > interval, ok
           * When the maximum recoil energy is greater than the highest energy of the recoil energy bin t
           * Case II also utilizes this piece of code with T_u = T_max
           */
          for (unsigned int i_T = 0; i_T < _quad_points.size(); ++i_T)
          {
            // Get the energies in Gaussian quadrature points and their respective weights
            Real T = 0.5 * (T_u - T_l) * _quad_points[i_T] + 0.5 * (T_u + T_l);
            Real w_T = 0.5 * _quad_weights[i_T] * (T_u - T_l);

            // Calculate cosine of recoil angle in the CM frame given the neutron and recoil atom energies
            Real mu_C = 1 - 2 * T / (E * _gamma);

            // Calculate cosine of recoil angle in the Lab frame according to geometry rules
            Real mu_L = sqrt((1 - mu_C) / 2);

            // Save maximum and mininum Lab frame cosine values
            _save_mu_L[t][g][0] = mu_L_max;
            _save_mu_L[t][g][1] = mu_L_min;

            /*
             * Calculate contribution to cross section coefficients
             * mu_L is scaled from its possible range of values [mu_L_min, mu_L_max] to fit the interval [-1,1]
             * of the Legendre polynomials
             */
             Real scaled_mu_L = 2 * (mu_L - mu_L_min) / (mu_L_max - mu_L_min) - 1;
            _erxs_coeff[l][t][g] += 1 / _xi_g[g] *
                                    _elastic_xs.value(E,Point()) *
                                    _neutron_spectrum.value(E, Point()) *
                                    _scattering_law.value(mu_C, Point()) *
                                    legendreP(l, scaled_mu_L) * w_T * w_E;
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
   */
  std::ofstream output_file;
  output_file.open(_erxs_file_name);
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

  /*
   * Writes output file with maximum and mininum cosines in the Lab frame (mu_L)
   * It follows same structure as ERXS outfile file, but saves the mu_L_max and
   * mu_L_min per g -> t combination.
   * mu_L_max, mu_L_min
   */
  std::ofstream output_file2;
  output_file2.open(_mu_L_file_name);
  for (unsigned int t = 0; t < _T; ++t)
  {
      for (unsigned int g = 0; g < _G; ++g)
      {
        if (g < _G - 1)
          output_file2 << _save_mu_L[t][g][0] << ',' << _save_mu_L[t][g][1] << ',';
        else
          output_file2 << _save_mu_L[t][g][0] << ',' << _save_mu_L[t][g][1];
      }
    output_file2 << std::endl;
  }
  output_file2.close();
}
