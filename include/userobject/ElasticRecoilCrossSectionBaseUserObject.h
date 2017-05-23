/****************************************************************/
/*            Elastic Recoil Cross Section Base                 */
/*          Base class for the further calculations             */
/*              Purely virtual (abstract) class                 */
/****************************************************************/

#ifndef ELASTICRECOILCROSSSECTIONBASEUSEROBJECT_H
#define ELASTICRECOILCROSSSECTIONBASEUSEROBJECT_H

// MOOSE includes
#include "GeneralUserObject.h"

// Forward Declarations
class ElasticRecoilCrossSectionBaseUserObject;

template <>
InputParameters validParams<ElasticRecoilCrossSectionBaseUserObject>();

class ElasticRecoilCrossSectionBaseUserObject : public GeneralUserObject
{
public:
  ElasticRecoilCrossSectionBaseUserObject(const InputParameters & parameters);

  /**
   * Purely virtual methods.
   */
  virtual Real elasticCrossSection() = 0;

  virtual Real neutronSpectrum() = 0;

  virtual Real scatteringLaw() = 0;

  virtual void execute() override;

  Real legendreP(unsigned int  n, Real x);

protected:

  std::vector<Real> _atomic_mass;
  std::vector<unsigned int> _isotope_type;
  unsigned int _number_isotope;
  std::vector<Real> _incident_energy_boundaries;
  std::vector<Real> _recoil_energy_boundaries;

};

#endif
