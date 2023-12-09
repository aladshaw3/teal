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

#pragma once

#include "Kernel.h"

/// HeatSource class object inherits from Kernel object
/** This class object inherits from the Kernel object in the MOOSE framework.

    The kernel adds the following physics:
      Res = test * v
*/
class HeatSource : public Kernel
{
public:
  /// Required new syntax for InputParameters
  static InputParameters validParams();

  /// Required constructor for objects in MOOSE
  HeatSource(const InputParameters & parameters);

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

  const VariableValue & _coupled_source;  ///< Coupled variable (W/m^3)
  const unsigned int _coupled_source_var; ///< Variable identification for the coupled variable

private:
};
