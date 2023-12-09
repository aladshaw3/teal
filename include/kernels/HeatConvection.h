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

#pragma once

#include "Kernel.h"

/// HeatConvection class object inherits from Kernel object
/** This class object inherits from the Kernel object in the MOOSE framework.

    The kernel adds the following physics:
      Res = test * h * A * fv * (T - T_other)
*/
class HeatConvection : public Kernel
{
public:
  /// Required new syntax for InputParameters
  static InputParameters validParams();

  /// Required constructor for objects in MOOSE
  HeatConvection(const InputParameters & parameters);

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

  const VariableValue & _hs;          ///< Variable for Heat transfer coefficient (W/m^2/K)
  const unsigned int _hs_var;         ///< Variable identification for hw
  const VariableValue & _other_temp;  ///< Variable for other phase temperature (K)
  const unsigned int _other_temp_var; ///< Variable identification for other phase temperature
  const VariableValue & _volfrac;     ///< Variable for volume fraction (-)
  const unsigned int _volfrac_var;    ///< Variable identification for volume fraction
  const VariableValue & _specarea;    ///< Variable for specific area (m^-1)
  const unsigned int _specarea_var;   ///< Variable identification for specific area

private:
};
