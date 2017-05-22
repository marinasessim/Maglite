/****************************************************************/
/*            Scattering Cross Section User Object              */
/*           Calculation of elastic cross section               */
/*                       Derived class                          */
/****************************************************************/

#include "ScatteringCrossSectionUO.h"
#include "ElasticRecoilCrossSectionBaseUserObject.h"

template <>
InputParameters
validParams<ScatteringCrossSectionUO>()
{
  InputParameters params = validParams<ElasticRecoilCrossSectionBaseUserObject>();

  return params;
}

ScatteringCrossSectionUO::ScatteringCrossSectionUO(const InputParameters & parameters)
  : ElasticRecoilCrossSectionBaseUserObject(parameters)
{
}

Real
ScatteringCrossSectionUO::elasticCrossSection()
{
 return 1;
}
