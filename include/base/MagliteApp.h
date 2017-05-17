#ifndef MAGLITEAPP_H
#define MAGLITEAPP_H

#include "MooseApp.h"

class MagliteApp;

template <>
InputParameters validParams<MagliteApp>();

class MagliteApp : public MooseApp
{
public:
  MagliteApp(InputParameters parameters);
  virtual ~MagliteApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* MAGLITEAPP_H */
