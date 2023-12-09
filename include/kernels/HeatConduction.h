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

#pragma once

#include "Kernel.h"

/// HeatConduction class object inherits from Kernel object
/** This class object inherits from the Kernel object in the MOOSE framework.

    The kernel adds the following physics:
      Res = grad_test * grad_u * K * fv
*/
class HeatConduction : public Kernel
{
public:
  /// Required new syntax for InputParameters
  static InputParameters validParams();

  /// Required constructor for objects in MOOSE
  HeatConduction(const InputParameters & parameters);

protected:
  /// Required residual function for standard kernels in MOOSE
  /** This function returns a residual contribution for this object.*/
  virtual Real computeQpResidual();

  /// Required Jacobian function for standard kernels in MOOSE
  /** This function returns a Jacobian contribution for this object. The Jacobian being
   computed is the associated diagonal element in the overall Jacobian matrix for the
   system and is used in preconditioning of the linear sub-problem. */
  virtual Real computeQpJacobian();

  /// Not Required, but aids in the preconditioning step
  /** This function returns the off diagonal Jacobian contribution for this object. By
   returning a non-zero value we will hopefully improve the convergence rate for the
   cross coupling of the variables. */
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const VariableValue & _conductivity;  ///< Thermal Conductivity variable (W/m/K)
  const unsigned int _conductivity_var; ///< Variable identification for the thermal conductivity
  const VariableValue & _volfrac;       ///< Variable for volume fraction (-)
  const unsigned int _volfrac_var;      ///< Variable identification for volume fractio

private:
};
