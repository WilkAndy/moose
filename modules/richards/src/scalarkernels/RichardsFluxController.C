/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/
#include "RichardsFluxController.h"

template<>
InputParameters validParams<RichardsFluxController>()
{
  InputParameters params = validParams<ScalarKernel>();
  params.addRequiredParam<PostprocessorName>("pp0", "The pp0 postprocessor.  This kernel enforces total_flux = pp0 - pp1*variable.");
  params.addRequiredParam<PostprocessorName>("pp1", "The pp1 postprocessor.  This kernel enforces total_flux = pp0 - pp1*variable.");
  params.addRequiredParam<Real>("total_flux", "Desired total flux.  This kernel enforces total_flux = pp0 - pp1*variable.");

  return params;
}

RichardsFluxController::RichardsFluxController(const std::string & name, InputParameters parameters) :
    ScalarKernel(name, parameters),
    _total_flux(getParam<Real>("total_flux")),
    _pp0_value(getPostprocessorValue("pp0")),
    _pp1_value(getPostprocessorValue("pp1"))
{
}

RichardsFluxController::~RichardsFluxController()
{
}

void
RichardsFluxController::reinit()
{
}

void
RichardsFluxController::computeResidual()
{
  DenseVector<Number> & re = _assembly.residualBlock(_var.number());
  for (_i = 0; _i < re.size(); _i++)
    re(_i) += computeQpResidual();
}

Real
RichardsFluxController::computeQpResidual()
{
  return _total_flux - _pp0_value + _pp1_value*_u[0];
}

void
RichardsFluxController::computeJacobian()
{
  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
  for (_i = 0; _i < ke.m(); _i++)
    ke(_i, _i) += computeQpJacobian();
}

Real
RichardsFluxController::computeQpJacobian()
{
  return _pp1_value;
}

void
RichardsFluxController::computeOffDiagJacobian(unsigned int jvar)
{
  Moose::out << "RichardsFluxController offdiag, jvar = " << jvar << "\n"; //_pp0_value = " << _pp0_value << " _pp1_value = " << _pp1_value << "_u[0] = " << _u[0] << "\n";
}

Real
RichardsFluxController::computeQpOffDiagJacobian(unsigned int /*jvar*/)
{
  return 0.;
}
