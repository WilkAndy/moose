/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/


#ifndef POROUSFLOWVARIABLESWITCHUPDATER_H
#define POROUSFLOWVARIABLESWITCHUPDATER_H

#include "AuxKernel.h"
#include "PorousFlowDictator.h"

//Forward Declarations
class PorousFlowVariableSwitchUpdater;

template<>
InputParameters validParams<PorousFlowVariableSwitchUpdater>();

/**
 * Updates in var switch
 */
class PorousFlowVariableSwitchUpdater: public AuxKernel
{
public:
  PorousFlowVariableSwitchUpdater(const InputParameters & parameters);

protected:
  virtual Real computeValue();

  /// van-Genuchten parameter
  const Real _al;

  /// van-Genuchten parameter
  const Real _m;

  /// van-Genuchten parameter
  const Real _p0;

  const PorousFlowDictator & _dictator;

  const VariableValue & _p_or_s;

  /// Whether the AuxVariable for this Kernel is in fact porepressure
  const bool _updating_pp;

  /// Encodes whether the _p_or_s at the node is porepressure (_encoder>0) or saturation (_encoder<0)
  const VariableValue & _encoder;
};

#endif // POROUSFLOWVARIABLESWITCHUPDATER_H
