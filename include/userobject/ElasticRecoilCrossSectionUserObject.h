/****************************************************************/
/*            Elastic Recoil Cross Section                      */
/****************************************************************/

#ifndef ELASTICRECOILCROSSSECTIONUSEROBJECT_H
#define ELASTICRECOILCROSSSECTIONUSEROBJECT_H

// MOOSE includes
#include "GeneralUserObject.h"

// Forward Declarations
class ElasticRecoilCrossSectionUserObject;

template <>
InputParameters validParams<ElasticRecoilCrossSectionUserObject>();

class ElasticRecoilCrossSectionUserObject : public GeneralUserObject
{
public:
  ElasticRecoilCrossSectionUserObject(const InputParameters & parameters);

  /**
   * Class destructor.
   */
  virtual ~ElasticRecoilCrossSectionUserObject();

  virtual void initialize() override;

  virtual void execute() override;

  virtual void finalize() override;

protected:

  /// Function representing the reactor neutron spectrum.
  Function & _neutron_spectrum;

  /// Function representing the epithermal spectrum scattering law.
  Function & _scattering_law;

  /// Function representing the simple elastic recoil cross section calculation.
  Function & _simple_erxs;

  Real legendreP(unsigned int  n, Real x);

  QBase * _quadrature;
  std::vector<Real> _quad_points;
  std::vector<Real> _quad_weights;

  std::vector<Real> _atomic_mass;
  std::vector<unsigned int> _isotope_type;
  unsigned int _number_isotope;
  std::vector<Real> _neutron_energy_limits;
  std::vector<Real> _recoil_energy_limits;

  Real _L;
  std::vector<Real> _xi_g;
  Real _gamma;
  unsigned int findNeutronEnergyGroup(Real energy);

  std::vector<std::vector<std::vector<Real>>> _elastic_recoil_xs;

};

#endif
