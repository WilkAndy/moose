/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef POROUSFLOW1PHASEPORS_VG_H
#define POROUSFLOW1PHASEPORS_VG_H

#include "PorousFlowVariableBase.h"

//Forward Declarations
class PorousFlow1PhasePorS_VG;

template<>
InputParameters validParams<PorousFlow1PhasePorS_VG>();

/**
 * Material designed to calculate fluid-phase porepressure and saturation
 * for the single-phase situation, with variable-switching between
 * porepressure and saturation.  A van-Genuchten capillary expression is used
 */
class PorousFlow1PhasePorS_VG : public PorousFlowVariableBase
{
public:
  PorousFlow1PhasePorS_VG(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties();

  virtual void computeProperties();

  virtual void computeQpProperties();

  /// van-Genuchten parameter
  const Real _al;

  /// van-Genuchten parameter
  const Real _m;

  /// van-Genuchten parameter
  const Real _p0;

  /// Nodal value of the decoder
  const VariableValue & _decoder_nodal_var;

  /** Flag determining whether the whole element has porepressure
   * nodes, or saturation nodes, or whether there is a mixture
   * of nodes.
   */
  int _entire_elem_decoder;

  /// MOOSE variable number of the primary variable
  const unsigned _prime_varnum;

  /// PorousFlow variable number of the primary variable
  const unsigned _dvar;

  /// nodal porepressure
  const VariableValue & _pp_nodal_var;

  /// qp porepressure
  const VariableValue & _pp_qp_var;

  /// Gradient(porepressure)
  const VariableGradient & _gradpp_qp_var;

  /// nodal saturation
  const VariableValue & _sat_nodal_var;

  /// qp saturation
  const VariableValue & _sat_qp_var;

  /// Gradient(saturation)
  const VariableGradient & _gradsat_qp_var;

  virtual void buildPS();
};

#endif //POROUSFLOW1PHASEPORS_VG_H
