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
  params.addRequiredParam<std::vector<unsigned int>>(
      "number_isotope", "Number of possible isotopes");

  //Later we will read these from some interface:
  params.addRequiredParam<std::vector<Real>>(
      "incident_energy", "Energy of the incident neutron");
  params.addRequiredParam<std::vector<Real>>(
      "recoil_energy", "Energy that was transferred to the recoil atom");

  return params;
}

ElasticRecoilCrossSectionBaseUserObject::ElasticRecoilCrossSectionBaseUserObject(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _atomic_mass(getParam<std::vector<Real>>("atomic_mass")),
    _isotope_type(getParam<std::vector<unsigned int>>("isotope_type")),
    _number_isotope(getParam<std::vector<unsigned int>>("number_isotope")),
    _incident_energy(getParam<std::vector<Real>>("incident_energy")),
    _recoil_energy(getParam<std::vector<Real>>("recoil_energy")),
    _elastic_xs(),
    _neutron_energy(),
    _scattering_law()
{
}

// Can I do this in a base class?
Real
ElasticRecoilCrossSectionBaseUserObject::elasticCrossSection()
{
  return _elastic_xs;
}

Real
ElasticRecoilCrossSectionBaseUserObject::neutronSpectrum()
{
  return _neutron_energy;
}

Real
ElasticRecoilCrossSectionBaseUserObject::scatteringLaw()
{
  return _scattering_law;
}

void
ElasticRecoilCrossSectionBaseUserObject::execute()
{
}
