/****************************************************************/
/*            Elastic Recoil Cross Section Base                 */
/*          Base class for the further calculations             */
/*              Purely virtual (abstract) class                 */
/****************************************************************/

#include "ScatteringLawUO.h"
#include "ElasticRecoilCrossSectionBaseUserObject.h"

template <>
InputParameters
validParams<ScatteringLawUO>()
{
  InputParameters params = validParams<ElasticRecoilCrossSectionBaseUserObject>();

  return params;
}

ScatteringLawUO::ScatteringLawUO(const InputParameters & parameters)
  : ElasticRecoilCrossSectionBaseUserObject(parameters)
{
}

Real
ScatteringLawUO::scatteringLaw()
{
 return 1;
}
