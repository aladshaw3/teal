/*!
 *  \file HeatAdvectionConservative.h
 *	\brief Kernel to create heat advection physics with upwinding schemes
 *	\details This file creates a heat advection kernel with optional upwinding
 *				and introduces the following phyiscs:
 *						Res = -grad_test * fv * vel * rho * cp * T
 *								where fv = volume fraction (-)
 *									  rho = material density (kg/m^3)
 *									  cp = heat capacity of the material (J/kg/K)
 *									  T = temperature of the fluid (K)
 *									  vel = velocity of the fluid (m/s)
 *
 * 	\note This REQUIRES use with ThermalFluidFluxBC due to Gauss Divergence
 *
 *  \author Austin Ladshaw
 *  \date 12/09/2023
 *  \copyright This kernel was modified from the ConservativeAdvection
 *				kernel provided by the MOOSE Framework (2023).
 */

#pragma once

#include "Kernel.h"

/**
 * Advection of the variable by the velocity provided by the user.
 * Options for numerical stabilization are: none; full upwinding
 */
class HeatAdvectionConservative : public Kernel
{
public:
  static InputParameters validParams();

  HeatAdvectionConservative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual void computeResidual() override;
  virtual void computeJacobian() override;

  /// Adding the off-diagonal components for better convergence
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  const VariableValue & _density;   ///< Coupled material density variable (kg/m^3)
  const unsigned int _density_var;  ///< Variable identification for density
  const VariableValue & _heat_cap;  ///< Coupled material heat capacity variable (J/kg/K)
  const unsigned int _heat_cap_var; ///< Variable identification for heat capacity
  const VariableValue & _volfrac;   ///< Variable for volume fraction (-)
  const unsigned int _volfrac_var;  ///< Variable identification for volume fraction

  const VariableValue & _ux;  ///< Velocity in the x-direction (m/s)
  const unsigned int _ux_var; ///< Variable identification for ux
  const VariableValue & _uy;  ///< Velocity in the y-direction (m/s)
  const unsigned int _uy_var; ///< Variable identification for uy
  const VariableValue & _uz;  ///< Velocity in the z-direction (m/s)
  const unsigned int _uz_var; ///< Variable identification for uz

  /// advection vector (constructed from piece-wise components)
  RealVectorValue _vec;

  /// enum to make the code clearer
  enum class JacRes
  {
    CALCULATE_RESIDUAL = 0,
    CALCULATE_JACOBIAN = 1
  };

  /// Type of upwinding
  const enum class UpwindingType { none, full } _upwinding;

  /// Nodal value of u, used for full upwinding
  const VariableValue & _u_nodal;

  /// In the full-upwind scheme, whether a node is an upwind node
  std::vector<bool> _upwind_node;

  /// In the full-upwind scheme d(total_mass_out)/d(variable_at_node_i)
  std::vector<Real> _dtotal_mass_out;

  /// Returns - _grad_test * velocity
  Real negSpeedQp();

  /// Calculates the fully-upwind Residual and Jacobian (depending on res_or_jac)
  void fullUpwind(JacRes res_or_jac);
};
