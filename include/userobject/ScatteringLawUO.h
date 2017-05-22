/****************************************************************/
/*            Elastic Recoil Cross Section Base                 */
/*          Base class for the further calculations             */
/*              Purely virtual (abstract) class                 */
/****************************************************************/

#ifndef SCATTERINGLAWUO_H
#define SCATTERINGLAWUO_H

// MOOSE includes
#include "ElasticRecoilCrossSectionBaseUserObject.h"

// Forward Declarations
class ScatteringLawUO;

template <>
InputParameters validParams<ScatteringLawUO>();

class ScatteringLawUO : public ElasticRecoilCrossSectionBaseUserObject
{
public:
  ScatteringLawUO(const InputParameters & parameters);

  /**
   * Overrides scatteringLaw() from the Base class.
   */
  virtual Real scatteringLaw() override;


protected:

};

#endif
