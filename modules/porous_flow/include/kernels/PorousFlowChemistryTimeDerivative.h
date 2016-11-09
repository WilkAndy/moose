/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#ifndef PORFLOWCHEMICALTIMEDERIVATIVE_H
#define PORFLOWCHEMICALTIMEDERIVATIVE_H

#include "TimeDerivative.h"
#include "PorousFlowDictator.h"

// Forward Declarations
class PorousFlowChemistryTimeDerivative;

template<>
InputParameters validParams<PorousFlowChemistryTimeDerivative>();

/**
 * Kernel = (conc - conc_old)/dt
 * where conc = porosity*sum_phases(saturation_phase*generalisedconc_phase^component)
 * It is lumped to the nodes
 */
class PorousFlowChemistryTimeDerivative : public TimeKernel
{
public:
  PorousFlowChemistryTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  virtual Real computeQpJacobian();

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /// the primary-species component index
  const unsigned _primary_species;

  /// holds info on the PorousFlow variables
  const PorousFlowDictator & _dictator;

  /// whether the Variable for this Kernel is a porous-flow variable according to the Dictator
  const bool _var_is_porflow_var;

  /// number of fluid phases
  const unsigned int _num_phases;

  /// porosity at the nodes, but it can depend on grad(variables) which are actually evaluated at the qps
  const MaterialProperty<Real> & _porosity;

  /// old value of porosity
  const MaterialProperty<Real> & _porosity_old;

  /// d(porosity)/d(porous-flow variable) - these derivatives will be wrt variables at the nodes
  const MaterialProperty<std::vector<Real> > & _dporosity_dvar;

  /// d(porosity)/d(grad porous-flow variable) - remember these derivatives will be wrt grad(vars) at qps
  const MaterialProperty<std::vector<RealGradient> > & _dporosity_dgradvar;

  /// nodal fluid saturation
  const MaterialProperty<std::vector<Real> > & _fluid_saturation_nodal;

  /// old value of fluid saturation
  const MaterialProperty<std::vector<Real> > & _fluid_saturation_nodal_old;

  /// d(nodal fluid saturation)/d(porous-flow variable)
  const MaterialProperty<std::vector<std::vector<Real> > > & _dfluid_saturation_nodal_dvar;

  /// Generalised concentrations in each phase
  const MaterialProperty<std::vector<std::vector<Real> > > & _gen_concentration;

  /// Generalised concentrations in each phase
  const MaterialProperty<std::vector<std::vector<Real> > > & _gen_concentration_old;

  /// Derivative of generalised concentrations in each phase
  const MaterialProperty<std::vector<std::vector<std::vector<Real> > > > & _dgen_concentration_dvar;  
  
  /**
   * Derivative of residual with respect to PorousFlow variable number pvar
   * This is used by both computeQpJacobian and computeQpOffDiagJacobian
   * @param pvar take the derivative of the residual wrt this PorousFlow variable
   */
  Real computeQpJac(unsigned int pvar);
};

#endif //PORFLOWCHEMICALTIMEDERIVATIVE_H
