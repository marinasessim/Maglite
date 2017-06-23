#include "MagliteApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

#include "ElasticRecoilCrossSectionUserObject.h"

template <>
InputParameters
validParams<MagliteApp>()
{
  InputParameters params = validParams<MooseApp>();
  return params;
}

MagliteApp::MagliteApp(InputParameters parameters) : MooseApp(parameters)
{
  Moose::registerObjects(_factory);
  ModulesApp::registerObjects(_factory);
  MagliteApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
  ModulesApp::associateSyntax(_syntax, _action_factory);
  MagliteApp::associateSyntax(_syntax, _action_factory);
}

MagliteApp::~MagliteApp() {}

// External entry point for dynamic application loading
extern "C" void
MagliteApp__registerApps()
{
  MagliteApp::registerApps();
}
void
MagliteApp::registerApps()
{
  registerApp(MagliteApp);
}

// External entry point for dynamic object registration
extern "C" void
MagliteApp__registerObjects(Factory & factory)
{
  MagliteApp::registerObjects(factory);
}
void
MagliteApp::registerObjects(Factory & factory)
{
  registerUserObject(ElasticRecoilCrossSectionUserObject);
}

// External entry point for dynamic syntax association
extern "C" void
MagliteApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory)
{
  MagliteApp::associateSyntax(syntax, action_factory);
}
void
MagliteApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
