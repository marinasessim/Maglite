/****************************************************************/
/*            Elastic Recoil Cross Section Base                 */
/*          Base class for the further calculations             */
/*              Purely virtual (abstract) class                 */
/****************************************************************/

#ifndef NEUTRONSPECTRUMUO_H
#define NEUTRONSPECTRUMUO_H

// MOOSE includes
#include "ElasticRecoilCrossSectionBaseUserObject.h"

// Forward Declarations
class NeutronSpectrumUO;

template <>
InputParameters validParams<NeutronSpectrumUO>();

class NeutronSpectrumUO : public ElasticRecoilCrossSectionBaseUserObject
{
public:
  NeutronSpectrumUO(const InputParameters & parameters);

  /**
   * Overrides neutronSpectrum() from the Base class.
   */
  virtual Real neutronSpectrum() override;


protected:

};

#endif
