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

#pragma once

#include "CoefTimeDerivative.h"

/// HeatAccumulation class object inherits from CoefTimeDerivative object
/** This class object inherits from the CoefTimeDerivative object in the MOOSE framework.

    The kernel adds the following physics:
      Res = test * fv * rho * cp * dTdt
*/
class HeatAccumulation : public CoefTimeDerivative
{
public:
  /// Required new syntax for InputParameters
  static InputParameters validParams();

  /// Required constructor for objects in MOOSE
  HeatAccumulation(const InputParameters & parameters);

protected:
  /// Required residual function for standard kernels in MOOSE
  /** This function returns a residual contribution for this object.*/
  virtual Real computeQpResidual() override;
  /// Required Jacobian function for standard kernels in MOOSE
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
  const unsigned int _volfrac_var;  ///< Variable identification for volume fractio
};
