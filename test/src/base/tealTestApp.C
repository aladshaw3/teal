//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "tealTestApp.h"
#include "tealApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
tealTestApp::validParams()
{
  InputParameters params = tealApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

tealTestApp::tealTestApp(InputParameters parameters) : MooseApp(parameters)
{
  tealTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

tealTestApp::~tealTestApp() {}

void
tealTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  tealApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"tealTestApp"});
    Registry::registerActionsTo(af, {"tealTestApp"});
  }
}

void
tealTestApp::registerApps()
{
  registerApp(tealApp);
  registerApp(tealTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
tealTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  tealTestApp::registerAll(f, af, s);
}
extern "C" void
tealTestApp__registerApps()
{
  tealTestApp::registerApps();
}
