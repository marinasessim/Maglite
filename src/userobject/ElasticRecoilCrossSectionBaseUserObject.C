/****************************************************************/
/*            Elastic Recoil Cross Section Base                 */
/*          Base class for the further calculations             */
/*              Purely virtual (abstract) class                 */
/****************************************************************/

// MOOSE includes
#include "ElasticRecoilCrossSectionBaseUserObject.h"

// libmesh includes
#include "libmesh/quadrature.h"

template <>
InputParameters
validParams<ElasticRecoilCrossSectionBaseUserObject>()
{
  InputParameters params = validParams<GeneralUserObject>();

  return params;
}

ElasticRecoilCrossSectionBaseUserObject::ElasticRecoilCrossSectionBaseUserObject(const InputParameters & parameters)
  : GeneralUserObject(parameters)
{
}
