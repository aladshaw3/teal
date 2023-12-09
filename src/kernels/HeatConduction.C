/*!
 *  \file HeatConduction.h
 *  \brief Kernel for creating a heat conduction
 *  \details This file creates a kernel for the conduction of heat
 *            in an energy balance equation as shown below:
 *                  Res = grad_test * grad_u * K * fv
 *                          where K = thermal conductivity (in W/m/K)
 *							and   fv = volume fraction (-)
 *
 *
 *  \author Austin Ladshaw
 *  \date 12/09/2023
 *  \copyright This kernel was designed and built by Austin Ladshaw (2023)
 */

#include "HeatConduction.h"

registerMooseObject("tealApp", HeatConduction);

InputParameters
HeatConduction::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredCoupledVar("thermal_conductivity",
                               "Name of the thermal conductivity variable (W/m/K)");
  params.addCoupledVar(
      "volume_frac", 1, "Variable for volume fraction (solid volume / total volume) (-)");
  return params;
}

HeatConduction::HeatConduction(const InputParameters & parameters)
  : Kernel(parameters),
    _conductivity(coupledValue("thermal_conductivity")),
    _conductivity_var(coupled("thermal_conductivity")),
    _volfrac(coupledValue("volume_frac")),
    _volfrac_var(coupled("volume_frac"))
{
}

Real
HeatConduction::computeQpResidual()
{
  return _volfrac[_qp] * _conductivity[_qp] * _grad_test[_i][_qp] * _grad_u[_qp];
}

Real
HeatConduction::computeQpJacobian()
{
  return _volfrac[_qp] * _conductivity[_qp] * _grad_test[_i][_qp] * _grad_phi[_j][_qp];
}

Real
HeatConduction::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _conductivity_var)
  {
    return _volfrac[_qp] * _phi[_j][_qp] * _grad_test[_i][_qp] * _grad_u[_qp];
  }
  if (jvar == _volfrac_var)
  {
    return _phi[_j][_qp] * _conductivity[_qp] * _grad_test[_i][_qp] * _grad_u[_qp];
  }
  return 0.0;
}
