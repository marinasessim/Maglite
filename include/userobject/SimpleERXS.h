/****************************************************************/
/*            Scattering Cross Section User Object              */
/*           Calculation of elastic cross section               */
/*                       Derived class                          */
/****************************************************************/

#ifndef SIMPLEERXS_H
#define SIMPLEERXS_H

// MOOSE includes
#include "ElasticRecoilCrossSectionBaseUserObject.h"

// Forward Declarations
class SimpleERXS;

template <>
InputParameters validParams<SimpleERXS>();

class SimpleERXS : public ElasticRecoilCrossSectionBaseUserObject
{
public:
  SimpleERXS(const InputParameters & parameters);

  /**
   * Overrides from the base class.
   */
  virtual Real elasticCrossSection() override;

  virtual Real neutronSpectrum() override;

  virtual Real scatteringLaw() override;

protected:

};

#endif
