/*!
 *  \file HeatSource.h
 *  \brief Kernel for creating a heat source or sink by volume
 *  \details This file creates a kernel for the coupling a heat source or sink
 *            into an energy balance equation as shown below:
 *                  Res = test * v
 *                          where v = coupled heat source (in W/m^3)
 *
 *
 *  \author Austin Ladshaw
 *  \date 12/09/2023
 *  \copyright This kernel was designed and built by Austin Ladshaw (2023)
 */

#include "HeatSource.h"

registerMooseObject("tealApp", HeatSource);

InputParameters
HeatSource::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredCoupledVar("coupled_source",
                               "Name of the coupled heat source variable (W/m^3)");
  return params;
}

HeatSource::HeatSource(const InputParameters & parameters)
  : Kernel(parameters),
    _coupled_source(coupledValue("coupled_source")),
    _coupled_source_var(coupled("coupled_source"))
{
}

Real
HeatSource::computeQpResidual()
{
  return -_test[_i][_qp] * _coupled_source[_qp];
}

Real
HeatSource::computeQpJacobian()
{
  return 0.0;
}

Real
HeatSource::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _coupled_source_var)
  {
    return -_test[_i][_qp] * _phi[_j][_qp];
  }
  return 0.0;
}
