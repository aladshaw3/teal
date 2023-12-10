/*!
 *  \file ThermalFluidFluxBC.h
 *	\brief Boundary Condition kernel for the thermal fluid flux across a boundary of the domain
 *	\details This file creates a generic boundary condition kernel for the flux of thermal fluids
 *			at a boundary.
 *
 *
 *  \author Austin Ladshaw
 *  \date 12/09/2023
 *  \copyright This kernel was modified from the ConservativeAdvection
 *				kernel provided by the MOOSE Framework (2023).
 */

#include "ThermalFluidFluxBC.h"

registerMooseObject("tealApp", ThermalFluidFluxBC);

InputParameters
ThermalFluidFluxBC::validParams()
{
  InputParameters params = IntegratedBC::validParams();
  params.addRequiredCoupledVar("density",
                               "The name of the density variable for the material (kg/m^3)");
  params.addRequiredCoupledVar("heat_capacity",
                               "The name of the heat capacity variable for the material (J/kg/K)");
  params.addCoupledVar(
      "volume_frac", 1, "Variable for volume fraction (solid volume / total volume) (-)");

  params.addRequiredCoupledVar("vel_x", "Variable for velocity in x-direction (m/s)");
  params.addRequiredCoupledVar("vel_y", "Variable for velocity in y-direction (m/s)");
  params.addRequiredCoupledVar("vel_z", "Variable for velocity in z-direction (m/s)");

  params.addRequiredCoupledVar("outside_temperature",
                               "Variable for the other phase temperature (K)");
  return params;
}

ThermalFluidFluxBC::ThermalFluidFluxBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _density(coupledValue("density")),
    _density_var(coupled("density")),
    _heat_cap(coupledValue("heat_capacity")),
    _heat_cap_var(coupled("heat_capacity")),
    _volfrac(coupledValue("volume_frac")),
    _volfrac_var(coupled("volume_frac")),

    _ux(coupledValue("vel_x")),
    _ux_var(coupled("vel_x")),
    _uy(coupledValue("vel_y")),
    _uy_var(coupled("vel_y")),
    _uz(coupledValue("vel_z")),
    _uz_var(coupled("vel_z")),
    _outside_temp(coupledValue("outside_temperature")),
    _outside_temp_var(coupled("outside_temperature"))
{
}

Real
ThermalFluidFluxBC::computeQpResidual()
{
  Real r = 0;

  _vec(0) = _ux[_qp];
  _vec(1) = _uy[_qp];
  _vec(2) = _uz[_qp];

  // Output
  if ((_vec)*_normals[_qp] > 0.0)
  {
    r += _test[_i][_qp] * (_vec * _normals[_qp]) * _u[_qp] * _density[_qp] * _heat_cap[_qp] *
         _volfrac[_qp];
  }
  // Input
  else
  {
    r += _test[_i][_qp] * (_vec * _normals[_qp]) * _outside_temp[_qp] * _density[_qp] *
         _heat_cap[_qp] * _volfrac[_qp];
  }

  return r;
}

Real
ThermalFluidFluxBC::computeQpJacobian()
{
  Real r = 0;

  _vec(0) = _ux[_qp];
  _vec(1) = _uy[_qp];
  _vec(2) = _uz[_qp];

  // Output
  if ((_vec)*_normals[_qp] > 0.0)
  {
    r += _test[_i][_qp] * (_vec * _normals[_qp]) * _phi[_j][_qp] * _density[_qp] * _heat_cap[_qp] *
         _volfrac[_qp];
  }
  // Input
  else
  {
    r += 0.0;
  }

  return r;
}

Real
ThermalFluidFluxBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  _vec(0) = _ux[_qp];
  _vec(1) = _uy[_qp];
  _vec(2) = _uz[_qp];

  Real r = 0;

  if (jvar == _ux_var)
  {
    // Output
    if ((_vec)*_normals[_qp] > 0.0)
    {
      r += _test[_i][_qp] * _u[_qp] * (_phi[_j][_qp] * _normals[_qp](0)) * _density[_qp] *
           _heat_cap[_qp] * _volfrac[_qp];
    }
    // Input
    else
    {
      r += _test[_i][_qp] * _outside_temp[_qp] * (_phi[_j][_qp] * _normals[_qp](0)) *
           _density[_qp] * _heat_cap[_qp] * _volfrac[_qp];
    }
    return r;
  }

  if (jvar == _uy_var)
  {
    // Output
    if ((_vec)*_normals[_qp] > 0.0)
    {
      r += _test[_i][_qp] * _u[_qp] * (_phi[_j][_qp] * _normals[_qp](1)) * _density[_qp] *
           _heat_cap[_qp] * _volfrac[_qp];
    }
    // Input
    else
    {
      r += _test[_i][_qp] * _outside_temp[_qp] * (_phi[_j][_qp] * _normals[_qp](1)) *
           _density[_qp] * _heat_cap[_qp] * _volfrac[_qp];
    }
    return r;
  }

  if (jvar == _uz_var)
  {
    // Output
    if ((_vec)*_normals[_qp] > 0.0)
    {
      r += _test[_i][_qp] * _u[_qp] * (_phi[_j][_qp] * _normals[_qp](2)) * _density[_qp] *
           _heat_cap[_qp] * _volfrac[_qp];
    }
    // Input
    else
    {
      r += _test[_i][_qp] * _outside_temp[_qp] * (_phi[_j][_qp] * _normals[_qp](2)) *
           _density[_qp] * _heat_cap[_qp] * _volfrac[_qp];
    }
    return r;
  }

  return 0.0;
}