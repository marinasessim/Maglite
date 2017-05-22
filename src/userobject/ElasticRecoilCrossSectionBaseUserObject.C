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

void
ElasticRecoilCrossSectionBaseUserObject::execute()
{
}
