/*!
 *  \file HeatConvection.h
 *  \brief Kernel for creating an exchange of thermal energy between two phases
 *  \details This file creates a kernel for the coupling a pair of heat variables in
 *            the same domain as a form of convective transfer:
 *                  Res = test * h * A * fv * (T - T_other)
 *                          where T = temperature of this heat variable's phase (K)
 *                          and T_other = temperature of the other heat variable's phase (K)
 *                          h = heat transfer coefficient (W/m^2/K)
 *                          A = specific contact area per volume between the phases (m^-1)
 *                              = area of solids per volume of solids
 *                          fv = volume fraction of the phases (volume solids / total volume)
 *
 *
 *  \author Austin Ladshaw
 *  \date 12/09/2023
 *  \copyright This kernel was designed and built by Austin Ladshaw (2023)
 */

#include "HeatConvection.h"

registerMooseObject("tealApp", HeatConvection);

InputParameters
HeatConvection::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addRequiredCoupledVar("convection_coeff",
                               "Variable for heat transfer coefficient (W/m^2/K)");
  params.addRequiredCoupledVar("coupled_temperature",
                               "Variable for the other phase temperature (K)");
  params.addCoupledVar(
      "volume_frac", 1, "Variable for volume fraction (solid volume / total volume) (-)");
  params.addRequiredCoupledVar(
      "specific_area",
      "Specific area for transfer [surface area of solids / volume solids] (m^-1)");
  return params;
}

HeatConvection::HeatConvection(const InputParameters & parameters)
  : Kernel(parameters),
    _hs(coupledValue("convection_coeff")),
    _hs_var(coupled("convection_coeff")),
    _other_temp(coupledValue("coupled_temperature")),
    _other_temp_var(coupled("coupled_temperature")),
    _volfrac(coupledValue("volume_frac")),
    _volfrac_var(coupled("volume_frac")),
    _specarea(coupledValue("specific_area")),
    _specarea_var(coupled("specific_area"))
{
}

Real
HeatConvection::computeQpResidual()
{
  return _test[_i][_qp] * _hs[_qp] * _specarea[_qp] * _volfrac[_qp] * (_u[_qp] - _other_temp[_qp]);
}

Real
HeatConvection::computeQpJacobian()
{
  return _test[_i][_qp] * _hs[_qp] * _specarea[_qp] * _volfrac[_qp] * _phi[_j][_qp];
}

Real
HeatConvection::computeQpOffDiagJacobian(unsigned int jvar)
{

  if (jvar == _other_temp_var)
  {
    return -_test[_i][_qp] * _hs[_qp] * _specarea[_qp] * _volfrac[_qp] * _phi[_j][_qp];
  }

  if (jvar == _hs_var)
  {
    return _test[_i][_qp] * _phi[_j][_qp] * _specarea[_qp] * _volfrac[_qp] *
           (_u[_qp] - _other_temp[_qp]);
  }

  if (jvar == _volfrac_var)
  {
    return _test[_i][_qp] * _hs[_qp] * _specarea[_qp] * _phi[_j][_qp] *
           (_u[_qp] - _other_temp[_qp]);
  }

  if (jvar == _specarea_var)
  {
    return _test[_i][_qp] * _hs[_qp] * _phi[_j][_qp] * _volfrac[_qp] * (_u[_qp] - _other_temp[_qp]);
  }

  return 0.0;
}
