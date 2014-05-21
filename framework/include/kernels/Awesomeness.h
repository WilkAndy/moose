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

#ifndef AWESOMENESS_H
#define AWESOMENESS_H

#include "Kernel.h"

class Awesomeness;

template<>
InputParameters validParams<Awesomeness>();

/**
 * Class to prove Moose's awesomeness
 */
class Awesomeness : public Kernel
{
public:
  Awesomeness(const std::string & name, InputParameters parameters);
  virtual ~Awesomeness();

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
};


#endif /* AWESOMENESS_H */
