/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/


#ifndef POROUSFLOW1PHASESET_H
#define POROUSFLOW1PHASESET_H

#include "NodalUserObject.h"

//Forward Declarations
class PorousFlow1PhaseSet;

template<>
InputParameters validParams<PorousFlow1PhaseSet>();

class PorousFlow1PhaseSet :
  public NodalUserObject
{
public:
  PorousFlow1PhaseSet(const InputParameters & parameters);

  /**
   * This function will get called on each node.
   */
  virtual void execute();

  virtual void initialize(){ Moose::out << "PF1PS initialise\n";};
  virtual void threadJoin(const UserObject&){};
  virtual void finalize();

protected:
  const unsigned _p_or_s_num;
  const AuxVariableName _encoder_name;
  const AuxVariableName _porepressure_name;
  const AuxVariableName _saturation_name;
};

#endif // POROUSFLOW1PHASESET_H
