/****************************************************************/
/*            Elastic Recoil Cross Section                      */
/****************************************************************/

// MOOSE includes
#include "ElasticRecoilCrossSectionUserObject.h"

// libmesh includes
#include "libmesh/quadrature.h"

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
      "neutron_energy_limits", "Energy boundaries of the incident neutron");
  params.addRequiredParam<std::vector<Real>>(
      "recoil_energy_limits", "Energy that was transferred to the recoil atom");
  params.addRequiredParam<FunctionName>(
      "neutron_spectrum","Function representing the reactor neutron spectrum");
  params.addRequiredParam<FunctionName>(
      "epithermal_law", "Function representing the scattering law for epithermal neutrons");
  params.addRequiredParam<FunctionName>(
      "fast_law", "Function representing the scattering law for fast neutrons");
  params.addRequiredParam<FunctionName>(
      "simple_erxs", "Function that handles the cross section calculation");

  return params;
}

ElasticRecoilCrossSectionUserObject::ElasticRecoilCrossSectionUserObject(const InputParameters & parameters)
  : GeneralUserObject(parameters),
    _neutron_spectrum(getFunction("neutron_spectrum")),
    _epithermal_law(getFunction("epithermal_law")),
    _fast_law(getFunction("fast_law")),
    _simple_erxs(getFunction("simple_erxs")),

    _atomic_mass(getParam<std::vector<Real>>("atomic_mass")),
    _isotope_type(getParam<std::vector<unsigned int>>("isotope_type")),
    _number_isotope(_isotope_type.size()),
    _neutron_energy_limits(getParam<std::vector<Real>>("neutron_energy_limits")),
    _recoil_energy_limits(getParam<std::vector<Real>>("recoil_energy_limits"))
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

// Class destructor
ElasticRecoilCrossSectionUserObject::~ElasticRecoilCrossSectionUserObject()
{
  delete _quadrature;
}
