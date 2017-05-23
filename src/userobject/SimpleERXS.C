/****************************************************************/
/*            Scattering Cross Section User Object              */
/*           Calculation of elastic cross section               */
/*                       Derived class                          */
/****************************************************************/

#include "SimpleERXS.h"
#include "ElasticRecoilCrossSectionBaseUserObject.h"

template <>
InputParameters
validParams<SimpleERXS>()
{
  InputParameters params = validParams<ElasticRecoilCrossSectionBaseUserObject>();

  return params;
}

SimpleERXS::SimpleERXS(const InputParameters & parameters)
  : ElasticRecoilCrossSectionBaseUserObject(parameters)
{
}

Real
SimpleERXS::elasticCrossSection()
{
  return 1;
}

Real
SimpleERXS::neutronSpectrum()
{
  return 1;
}

Real
SimpleERXS::scatteringLaw()
{
  return 1;
}
