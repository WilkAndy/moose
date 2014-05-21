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

#include "Awesomeness.h"


template<>
InputParameters validParams<Awesomeness>()
{
  InputParameters p = validParams<Kernel>();
  return p;
}


Awesomeness::Awesomeness(const std::string & name, InputParameters parameters) :
    Kernel(name, parameters)
{
  Moose::out << "Moose is awesome\n";
}

Awesomeness::~Awesomeness()
{
}

Real
Awesomeness::computeQpResidual()
{
  return 0.0;
}

Real
Awesomeness::computeQpJacobian()
{
  return 0.0;
}
