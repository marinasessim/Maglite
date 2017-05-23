/****************************************************************/
/*            Elastic Recoil Cross Section Base                 */
/*          Base class for the further calculations             */
/*              Purely virtual (abstract) class                 */
/****************************************************************/

// MOOSE includes
#include "ElasticRecoilCrossSectionBaseUserObject.h"

// libmesh includes
#include "libmesh/quadrature.h"

template <>
InputParameters
validParams<ElasticRecoilCrossSectionBaseUserObject>()
{
  InputParameters params = validParams<GeneralUserObject>();
  params.addRequiredParam<std::vector<Real>>(
      "atomic_mass", "Atomic Mass of the isotopes");
  params.addRequiredParam<std::vector<unsigned int>>(
      "isotope_type", "Type of the target isotopes");

  //Later we will read these from some interface:
  params.addRequiredParam<std::vector<Real>>(
      "incident_energy_boundaries", "Energy boundaries of the incident neutron");
  params.addRequiredParam<std::vector<Real>>(
      "recoil_energy_boundaries", "Energy that was transferred to the recoil atom");

  return params;
}

ElasticRecoilCrossSectionBaseUserObject::ElasticRecoilCrossSectionBaseUserObject(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _atomic_mass(getParam<std::vector<Real>>("atomic_mass")),
    _isotope_type(getParam<std::vector<unsigned int>>("isotope_type")),
    _number_isotope(_isotope_type.size()),
    _incident_energy_boundaries(getParam<std::vector<Real>>("incident_energy")),
    _recoil_energy_boundaries(getParam<std::vector<Real>>("recoil_energy"))
{
}

Real
ElasticRecoilCrossSectionBaseUserObject::legendreP(unsigned int  n, Real x)
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

void
ElasticRecoilCrossSectionBaseUserObject::execute()
{
}
