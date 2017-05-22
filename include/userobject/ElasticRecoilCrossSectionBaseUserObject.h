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
  virtual Real neutronSpectrum() = 0; //How do I represent parameters here?

  virtual Real elasticCrossSection() = 0;

  virtual Real scatteringLaw() = 0;

protected:

};

#endif
