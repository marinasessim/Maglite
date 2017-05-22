/****************************************************************/
/*            Scattering Cross Section User Object              */
/*           Calculation of elastic cross section               */
/*                       Derived class                          */
/****************************************************************/

#ifndef SCATTERINGCROSSSECTIONUO_H
#define SCATTERINGCROSSSECTIONUO_H

// MOOSE includes
#include "ElasticRecoilCrossSectionBaseUserObject.h"

// Forward Declarations
class ScatteringCrossSectionUO;

template <>
InputParameters validParams<ScatteringCrossSectionUO>();

class ScatteringCrossSectionUO : public ElasticRecoilCrossSectionBaseUserObject
{
public:
  ScatteringCrossSectionUO(const InputParameters & parameters);

  /**
   * Overrides elasticCrossSection() from the Base class.
   */
  virtual Real elasticCrossSection() override;


protected:

};

#endif
