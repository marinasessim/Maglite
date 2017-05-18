/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "GeneralUserObject.h"
#include "NumeratorUserObject.h"

template <>
InputParameters
validParams<NumeratorUserObject>()
{
  InputParameters params = validParams<UserObject>();
  params += validParams<MaterialPropertyInterface>();
  return params;
}

NumeratorUserObject::NumeratorUserObject(const InputParameters & parameters)
  : UserObject(parameters),
    MaterialPropertyInterface(this),
    TransientInterface(this),
    DependencyResolverInterface(),
    UserObjectInterface(this),
    PostprocessorInterface(this),
    VectorPostprocessorInterface(this)
{
  _supplied_vars.insert(name());
}

const std::set<std::string> &
NumeratorUserObject::getRequestedItems()
{
  return _depend_vars;
}

const std::set<std::string> &
NumeratorUserObject::getSuppliedItems()
{
  return _supplied_vars;
}

const PostprocessorValue &
NumeratorUserObject::getPostprocessorValue(const std::string & name)
{
  _depend_vars.insert(_pars.get<PostprocessorName>(name));
  return PostprocessorInterface::getPostprocessorValue(name);
}

const PostprocessorValue &
NumeratorUserObject::getPostprocessorValueByName(const PostprocessorName & name)
{
  _depend_vars.insert(name);
  return PostprocessorInterface::getPostprocessorValueByName(name);
}

const VectorPostprocessorValue &
NumeratorUserObject::getVectorPostprocessorValue(const std::string & name,
                                               const std::string & vector_name)
{
  _depend_vars.insert(_pars.get<VectorPostprocessorName>(name));
  return VectorPostprocessorInterface::getVectorPostprocessorValue(name, vector_name);
}

const VectorPostprocessorValue &
NumeratorUserObject::getVectorPostprocessorValueByName(const VectorPostprocessorName & name,
                                                     const std::string & vector_name)
{
  _depend_vars.insert(name);
  return VectorPostprocessorInterface::getVectorPostprocessorValueByName(name, vector_name);
}

void
NumeratorUserObject::threadJoin(const UserObject &)
{
  mooseError("GeneralUserObjects do not execute using threads, this function does nothing and "
             "should not be used.");
}

void
NumeratorUserObject::subdomainSetup()
{
  mooseError("GeneralUserObjects do not execute subdomainSetup method, this function does nothing "
             "and should not be used.");
}

void
NumeratorUserObject::execute()
{
}
