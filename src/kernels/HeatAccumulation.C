/*!
 *  \file HeatAccumulation.h
 *	\brief Kernel to create a heat accumulation kernel for thermal dynamics
 *	\details This file creates a heat accumulation kernel for thermal dynamics
 *				and introduces the following phyiscs:
 *						Res = test * fv * rho * cp * dTdt
 *								where fv = volume fraction (-)
 *									  rho = material density (kg/m^3)
 *									  cp = heat capacity of the material (J/kg/K)
 *									  dTdt = internal heat rate change (K/s)
 *
 *  \author Austin Ladshaw
 *  \date 12/09/2023
 *  \copyright This kernel was designed and built by Austin Ladshaw (2023)
 */

#include "HeatAccumulation.h"

registerMooseObject("tealApp", HeatAccumulation);

InputParameters
HeatAccumulation::validParams()
{
  InputParameters params = CoefTimeDerivative::validParams();
  params.addRequiredCoupledVar("density",
                               "The name of the density variable for the material (kg/m^3)");
  params.addRequiredCoupledVar("heat_capacity",
                               "The name of the heat capacity variable for the material (J/kg/K)");
  params.addCoupledVar(
      "volume_frac", 1, "Variable for volume fraction (solid volume / total volume) (-)");
  return params;
}

HeatAccumulation::HeatAccumulation(const InputParameters & parameters)
  : CoefTimeDerivative(parameters),
    _density(coupledValue("density")),
    _density_var(coupled("density")),
    _heat_cap(coupledValue("heat_capacity")),
    _heat_cap_var(coupled("heat_capacity")),
    _volfrac(coupledValue("volume_frac")),
    _volfrac_var(coupled("volume_frac"))
{
}

Real
HeatAccumulation::computeQpResidual()
{
  _coef = _density[_qp] * _heat_cap[_qp] * _volfrac[_qp];
  return CoefTimeDerivative::computeQpResidual();
}

Real
HeatAccumulation::computeQpJacobian()
{
  _coef = _density[_qp] * _heat_cap[_qp] * _volfrac[_qp];
  return CoefTimeDerivative::computeQpJacobian();
}

Real
HeatAccumulation::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (jvar == _density_var)
  {
    return _phi[_j][_qp] * _heat_cap[_qp] * _volfrac[_qp] * _test[_i][_qp] * _u_dot[_qp];
  }

  if (jvar == _heat_cap_var)
  {
    return _density[_qp] * _phi[_j][_qp] * _volfrac[_qp] * _test[_i][_qp] * _u_dot[_qp];
  }

  if (jvar == _volfrac_var)
  {
    return _density[_qp] * _heat_cap[_qp] * _phi[_j][_qp] * _test[_i][_qp] * _u_dot[_qp];
  }

  return 0.0;
}
