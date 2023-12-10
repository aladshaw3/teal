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

#pragma once

#include "IntegratedBC.h"
#include "libmesh/vector_value.h"

/// ThermalFluidFluxBC class object inherits from IntegratedBC object
/** This class object inherits from the IntegratedBC object.

  The flux BC uses the velocity in the system to apply a boundary
  condition based on whether or not material is leaving or entering the boundary. */
class ThermalFluidFluxBC : public IntegratedBC
{
public:
  /// Required new syntax for InputParameters
  static InputParameters validParams();

  /// Required constructor for BC objects in MOOSE
  ThermalFluidFluxBC(const InputParameters & parameters);

protected:
  /// Required function override for BC objects in MOOSE
  /** This function returns a residual contribution for this object.*/
  virtual Real computeQpResidual() override;
  /// Required function override for BC objects in MOOSE
  /** This function returns a Jacobian contribution for this object. The Jacobian being
    computed is the associated diagonal element in the overall Jacobian matrix for the
    system and is used in preconditioning of the linear sub-problem. */
  virtual Real computeQpJacobian() override;
  /// Not Required, but aids in the preconditioning step
  /** This function returns the off diagonal Jacobian contribution for this object. By
    returning a non-zero value we will hopefully improve the convergence rate for the
    cross coupling of the variables. */
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

  const VariableValue & _outside_temp;  ///< Variable for other phase temperature (K)
  const unsigned int _outside_temp_var; ///< Variable identification for other phase temperature

private:
};
