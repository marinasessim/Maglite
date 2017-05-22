/****************************************************************/
/*            Elastic Recoil Cross Section Base                 */
/*          Base class for the further calculations             */
/*              Purely virtual (abstract) class                 */
/****************************************************************/

#include "NeutronSpectrumUO.h"
#include "ElasticRecoilCrossSectionBaseUserObject.h"

template <>
InputParameters
validParams<NeutronSpectrumUO>()
{
  InputParameters params = validParams<ElasticRecoilCrossSectionBaseUserObject>();

  return params;
}

NeutronSpectrumUO::NeutronSpectrumUO(const InputParameters & parameters)
  : ElasticRecoilCrossSectionBaseUserObject(parameters)
{
}

Real
NeutronSpectrumUO::neutronSpectrum()
{
 return 1;
}
