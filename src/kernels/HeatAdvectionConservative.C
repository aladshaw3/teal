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

#include "HeatAdvectionConservative.h"
#include "SystemBase.h"

registerMooseObject("tealApp", HeatAdvectionConservative);

InputParameters
HeatAdvectionConservative::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Conservative form of $\\nabla \\cdot \\vec{v} u$ which in its weak "
                             "form is given by: $(-\\nabla \\psi_i, \\vec{v} u)$.");

  params.addRequiredCoupledVar("density",
                               "The name of the density variable for the material (kg/m^3)");
  params.addRequiredCoupledVar("heat_capacity",
                               "The name of the heat capacity variable for the material (J/kg/K)");
  params.addCoupledVar(
      "volume_frac", 1, "Variable for volume fraction (solid volume / total volume) (-)");

  params.addRequiredCoupledVar("vel_x", "Variable for velocity in x-direction (m/s)");
  params.addRequiredCoupledVar("vel_y", "Variable for velocity in y-direction (m/s)");
  params.addRequiredCoupledVar("vel_z", "Variable for velocity in z-direction (m/s)");

  MooseEnum upwinding_type("none full", "none");
  params.addParam<MooseEnum>("upwinding_type",
                             upwinding_type,
                             "Type of upwinding used.  None: Typically results in overshoots and "
                             "undershoots, but numerical diffusion is minimized.  Full: Overshoots "
                             "and undershoots are avoided, but numerical diffusion is large");
  return params;
}

HeatAdvectionConservative::HeatAdvectionConservative(const InputParameters & parameters)
  : Kernel(parameters),

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

    _upwinding(getParam<MooseEnum>("upwinding_type").getEnum<UpwindingType>()),
    _u_nodal(_var.dofValues()),
    _upwind_node(0),
    _dtotal_mass_out(0)
{
}

Real
HeatAdvectionConservative::negSpeedQp()
{
  _vec(0) = _ux[_qp];
  _vec(1) = _uy[_qp];
  _vec(2) = _uz[_qp];
  return -(_grad_test[_i][_qp] * _vec) * _density[_qp] * _heat_cap[_qp] * _volfrac[_qp];
}

Real
HeatAdvectionConservative::computeQpResidual()
{
  // This is the no-upwinded version
  // It gets called via Kernel::computeResidual()
  return negSpeedQp() * _u[_qp];
}

Real
HeatAdvectionConservative::computeQpJacobian()
{
  // This is the no-upwinded version
  // It gets called via Kernel::computeJacobian()
  return negSpeedQp() * _phi[_j][_qp];
}

void
HeatAdvectionConservative::computeResidual()
{
  switch (_upwinding)
  {
    case UpwindingType::none:
      Kernel::computeResidual();
      break;
    case UpwindingType::full:
      fullUpwind(JacRes::CALCULATE_RESIDUAL);
      break;
  }
}

void
HeatAdvectionConservative::computeJacobian()
{
  switch (_upwinding)
  {
    case UpwindingType::none:
      Kernel::computeJacobian();
      break;
    case UpwindingType::full:
      fullUpwind(JacRes::CALCULATE_JACOBIAN);
      break;
  }
}

Real
HeatAdvectionConservative::computeQpOffDiagJacobian(unsigned int jvar)
{
  _vec(0) = _ux[_qp];
  _vec(1) = _uy[_qp];
  _vec(2) = _uz[_qp];

  if (jvar == _ux_var)
  {
    return -_u[_qp] * (_phi[_j][_qp] * _grad_test[_i][_qp](0)) * _density[_qp] * _heat_cap[_qp] *
           _volfrac[_qp];
  }

  if (jvar == _uy_var)
  {
    return -_u[_qp] * (_phi[_j][_qp] * _grad_test[_i][_qp](1)) * _density[_qp] * _heat_cap[_qp] *
           _volfrac[_qp];
  }

  if (jvar == _uz_var)
  {
    return -_u[_qp] * (_phi[_j][_qp] * _grad_test[_i][_qp](2)) * _density[_qp] * _heat_cap[_qp] *
           _volfrac[_qp];
  }

  if (jvar == _density_var)
  {
    return -_u[_qp] * (_grad_test[_i][_qp] * _vec) * _phi[_j][_qp] * _heat_cap[_qp] * _volfrac[_qp];
  }

  if (jvar == _heat_cap_var)
  {
    return -_u[_qp] * (_grad_test[_i][_qp] * _vec) * _density[_qp] * _phi[_j][_qp] * _volfrac[_qp];
  }

  if (jvar == _volfrac_var)
  {
    return -_u[_qp] * (_grad_test[_i][_qp] * _vec) * _density[_qp] * _heat_cap[_qp] * _phi[_j][_qp];
  }

  return 0.0;
}

void
HeatAdvectionConservative::fullUpwind(JacRes res_or_jac)
{
  // The number of nodes in the element
  const unsigned int num_nodes = _test.size();

  // Even if we are computing the Jacobian we still need to compute the outflow from each node to
  // see which nodes are upwind and which are downwind
  prepareVectorTag(_assembly, _var.number());

  if (res_or_jac == JacRes::CALCULATE_JACOBIAN)
    prepareMatrixTag(_assembly, _var.number(), _var.number());

  // Compute the outflux from each node and store in _local_re
  // If _local_re is positive at the node, mass (or whatever the Variable represents) is flowing out
  // of the node
  _upwind_node.resize(num_nodes);
  for (_i = 0; _i < num_nodes; ++_i)
  {
    for (_qp = 0; _qp < _qrule->n_points(); _qp++)
      _local_re(_i) += _JxW[_qp] * _coord[_qp] * negSpeedQp();
    _upwind_node[_i] = (_local_re(_i) >= 0.0);
  }

  // Variables used to ensure mass conservation
  Real total_mass_out = 0.0;
  Real total_in = 0.0;
  if (res_or_jac == JacRes::CALCULATE_JACOBIAN)
    _dtotal_mass_out.assign(num_nodes, 0.0);

  for (unsigned int n = 0; n < num_nodes; ++n)
  {
    if (_upwind_node[n])
    {
      if (res_or_jac == JacRes::CALCULATE_JACOBIAN)
      {
        if (_test.size() == _phi.size())
          /* u at node=n depends only on the u at node=n, by construction.  For
           * linear-lagrange variables, this means that Jacobian entries involving the derivative
           * will only be nonzero for derivatives wrt variable at node=n.  Hence the
           * (n, n) in the line below.  The above "if" statement catches other variable types
           * (eg constant monomials)
           */
          _local_ke(n, n) += _local_re(n);

        _dtotal_mass_out[n] += _local_ke(n, n);
      }
      _local_re(n) *= _u_nodal[n];
      total_mass_out += _local_re(n);
    }
    else                        // downwind node
      total_in -= _local_re(n); // note the -= means the result is positive
  }

  // Conserve mass over all phases by proportioning the total_mass_out mass to the inflow nodes,
  // weighted by their local_re values
  for (unsigned int n = 0; n < num_nodes; ++n)
  {
    if (!_upwind_node[n]) // downwind node
    {
      if (res_or_jac == JacRes::CALCULATE_JACOBIAN)
        for (_j = 0; _j < _phi.size(); _j++)
          _local_ke(n, _j) += _local_re(n) * _dtotal_mass_out[_j] / total_in;
      _local_re(n) *= total_mass_out / total_in;
    }
  }

  // Add the result to the residual and jacobian
  if (res_or_jac == JacRes::CALCULATE_RESIDUAL)
  {
    accumulateTaggedLocalResidual();

    if (_has_save_in)
    {
      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
      for (const auto & var : _save_in)
        var->sys().solution().add_vector(_local_re, var->dofIndices());
    }
  }

  if (res_or_jac == JacRes::CALCULATE_JACOBIAN)
  {
    accumulateTaggedLocalMatrix();

    if (_has_diag_save_in)
    {
      unsigned int rows = _local_ke.m();
      DenseVector<Number> diag(rows);
      for (unsigned int i = 0; i < rows; i++)
        diag(i) = _local_ke(i, i);

      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
      for (const auto & var : _diag_save_in)
        var->sys().solution().add_vector(diag, var->dofIndices());
    }
  }
}
