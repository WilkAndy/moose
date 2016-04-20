/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/


#ifndef PORFLOWMATERIALMASSFRACTIONBUILDER_H
#define PORFLOWMATERIALMASSFRACTIONBUILDER_H

#include "DerivativeMaterialInterface.h"
#include "Material.h"

#include "PorousFlowDictator.h"

//Forward Declarations
class PorousFlowMaterialMassFractionBuilder;

template<>
InputParameters validParams<PorousFlowMaterialMassFractionBuilder>();

/**
 * Material designed to form a std::vector<std::vector>
 * of mass fractions from the individual mass fraction variables
 */
class PorousFlowMaterialMassFractionBuilder : public DerivativeMaterialInterface<Material>
{
public:
  PorousFlowMaterialMassFractionBuilder(const InputParameters & parameters);

protected:

  /// The variable names UserObject for the Porous-Flow variables
  const PorousFlowDictator & _porflow_name_UO;

  unsigned int _num_phases;

  unsigned int _num_components;


  MaterialProperty<std::vector<std::vector<Real> > > & _mass_frac;
  MaterialProperty<std::vector<std::vector<Real> > > & _mass_frac_old;
  MaterialProperty<std::vector<std::vector<RealGradient> > > & _grad_mass_frac;
  MaterialProperty<std::vector<std::vector<std::vector<Real> > > > & _dmass_frac_dvar;

  unsigned int _num_passed_mf_vars;

  std::vector<unsigned int> _mf_vars_num;
  std::vector<const VariableValue *> _mf_vars;
  std::vector<const VariableGradient *> _grad_mf_vars;

  virtual void initQpStatefulProperties();
  virtual void computeQpProperties();

 private:
  void build_mass_frac(unsigned int qp);
};

#endif //PORFLOWMATERIALMASSFRACTIONBUILDER_H
