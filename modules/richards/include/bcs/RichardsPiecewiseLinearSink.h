/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/


#ifndef RICHARDSPIECEWISELINEARSINK
#define RICHARDSPIECEWISELINEARSINK

#include "IntegratedBC.h"
#include "LinearInterpolation.h"
#include "Function.h"
#include "RichardsVarNames.h"
#include "RichardsDensity.h"
#include "RichardsRelPerm.h"
#include "RichardsSeff.h"


// Forward Declarations
class RichardsPiecewiseLinearSink;

template<>
InputParameters validParams<RichardsPiecewiseLinearSink>();

/**
 * Applies a flux sink to a boundary
 * The sink is a piecewise linear function of
 * porepressure (the "variable") at the quad points.
 * This is specified by _sink_func.
 * In addition, this sink can be multiplied by:
 *  (1) the relative permeability of the fluid at the quad point.
 *  (2) perm_nn*density/viscosity, where perm_nn is the
 *      permeability tensor projected in the normal direction.
 *  (3) a Function (which can be time-dependent, for instance)
 */
class RichardsPiecewiseLinearSink : public IntegratedBC
{
public:

  RichardsPiecewiseLinearSink(const std::string & name, InputParameters parameters);

protected:
  virtual void computeResidual();

  virtual Real computeQpResidual();

  virtual void computeJacobian();

  virtual Real computeQpJacobian();

  virtual void computeJacobianBlock(unsigned int jvar);

  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  /// whether to multiply the sink flux by permeability*density/viscosity
  bool _use_mobility;

  /// whether to multiply the sink flux by relative permeability
  bool _use_relperm;

  /// whether to use full upwinding
  bool _fully_upwind;

  /// piecewise-linear function of porepressure (this defines the strength of the sink)
  LinearInterpolation _sink_func;

  /// sink flux gets multiplied by this function
  Function & _m_func;

  /// holds info about the names and values of richards variable in the simulation
  const RichardsVarNames & _richards_name_UO;

  /// number of richards variables
  unsigned int _num_p;

  /// the moose internal variable number corresponding to the porepressure of this sink flux
  unsigned int _pvar;

  /// user object defining the density.  Only used if _fully_upwind = true
  const RichardsDensity * _density_UO;

  /// user object defining the effective saturation.  Only used if _fully_upwind = true
  const RichardsSeff * _seff_UO;

  /// user object defining the relative permeability.  Only used if _fully_upwind = true
  const RichardsRelPerm * _relperm_UO;

  /// number of nodes in this element.  Only used if _fully_upwind = true
  unsigned int _num_nodes;

  /**
   * nodal values of fluid density
   * These are used if _fully_upwind = true
   */
  std::vector<Real> _nodal_density;

  /**
   * d(_nodal_density)/d(variable_ph)  (variable_ph is the variable for phase=ph)
   * These are used in the jacobian calculations if _fully_upwind = true
   */
  std::vector<std::vector<Real> > _dnodal_density_dv;

  /**
   * nodal values of relative permeability
   * These are used if _fully_upwind = true
   */
  std::vector<Real> _nodal_relperm;

  /**
   * d(_nodal_relperm)/d(variable_ph)  (variable_ph is the variable for phase=ph)
   * These are used in the jacobian calculations if _fully_upwind = true
   */
  std::vector<std::vector<Real> > _dnodal_relperm_dv;

  /// porepressure values (only the _pvar component is used)
  MaterialProperty<std::vector<Real> > & _pp;

  /// d(porepressure_i)/d(variable_j)
  MaterialProperty<std::vector<std::vector<Real> > > & _dpp_dv;

  /// viscosity (only the _pvar component is used)
  MaterialProperty<std::vector<Real> > & _viscosity;

  /// permeability
  MaterialProperty<RealTensorValue> & _permeability;

  /**
   * derivative of effective saturation wrt variables
   * only _dseff_dv[_pvar][i] is used for i being all variables
   */
  MaterialProperty<std::vector<std::vector<Real> > > & _dseff_dv;

  /// relative permeability (only the _pvar component is used)
  MaterialProperty<std::vector<Real> > & _rel_perm;

  /// d(relperm_i)/d(variable_j)
  MaterialProperty<std::vector<std::vector<Real> > > & _drel_perm_dv;

  /// fluid density (only the _pvar component is used)
  MaterialProperty<std::vector<Real> > & _density;

  /// d(density_i)/d(variable_j)
  MaterialProperty<std::vector<std::vector<Real> > > & _ddensity_dv;

  /// optional Scalar variable coupling: this is the scalar variable's number
  unsigned int _p0_var;

  /// optional Scalar variable coupling: _p0[0] is the scalar variable's value
  VariableValue & _p0;

  /**
   * Holds the values of pressures at all the nodes of the element
   * Only used if _fully_upwind = true
   * Eg:
   * _ps_at_nodes[_pvar] is a pointer to this variable's nodal porepressure values
   * So: (*_ps_at_nodes[_pvar])[i] = _var.nodalSln()[i] = porepressure of pressure-variable _pvar at node i
   */
  std::vector<VariableValue *> _ps_at_nodes;


  /// calculates the nodal values of pressure, mobility, and derivatives thereof
  void prepareNodalValues();

  /// derivative of residual wrt the wrt_num Richards variable
  Real jac(unsigned int wrt_num);

  /// derivative of residual wrt to the optional scalar variable
  Real jacp0();


};

#endif //RichardsPiecewiseLinearSink
