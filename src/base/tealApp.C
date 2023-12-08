#include "tealApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
tealApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

tealApp::tealApp(InputParameters parameters) : MooseApp(parameters)
{
  tealApp::registerAll(_factory, _action_factory, _syntax);
}

tealApp::~tealApp() {}

void 
tealApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<tealApp>(f, af, s);
  Registry::registerObjectsTo(f, {"tealApp"});
  Registry::registerActionsTo(af, {"tealApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
tealApp::registerApps()
{
  registerApp(tealApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
tealApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  tealApp::registerAll(f, af, s);
}
extern "C" void
tealApp__registerApps()
{
  tealApp::registerApps();
}
