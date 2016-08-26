/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlow1PhasePorS_VG.h"
#include "PorousFlowEffectiveSaturationVG.h"
#include "PorousFlowCapillaryPressureVG.h"
#include "libmesh/quadrature.h"
#include <limits>

template<>
InputParameters validParams<PorousFlow1PhasePorS_VG>()
{
  InputParameters params = validParams<PorousFlowVariableBase>();
  params.addRequiredCoupledVar("decoder", "Variable that specifies whether the node should use saturation or porepressure as its primary variable");
  params.addRequiredCoupledVar("primary_variable", "Variable that holds either the value of porepressure or saturation, depending on decoder");
  params.addRequiredCoupledVar("porepressure", "Aux Variable that is the porepressure");
  params.addRequiredCoupledVar("saturation", "Aux Variable that is the saturation");
  params.addRequiredRangeCheckedParam<Real>("al", "al > 0", "van-Genuchten alpha parameter.  Must be positive.  effectiveSaturation = (1 + (-al*c)^(1/(1-m)))^(-m)");
  params.addRequiredRangeCheckedParam<Real>("m", "m > 0 & m < 1", "van-Genuchten m parameter.  Must be between 0 and 1, and optimally should be set to >0.5   EffectiveSaturation = (1 + (-al*p)^(1/(1-m)))^(-m)");
  params.addClassDescription("This Material is used for the single-phase situation where the primary variable may be either porepressure or saturation.  A van-Genuchten capillary function is assumed");
  return params;
}

PorousFlow1PhasePorS_VG::PorousFlow1PhasePorS_VG(const InputParameters & parameters) :
    PorousFlowVariableBase(parameters),
    _al(getParam<Real>("al")),
    _m(getParam<Real>("m")),
    _p0(1.0/_al),

    _decoder_nodal_var(coupledNodalValue("decoder")),
    _entire_elem_decoder(0),

    _prime_varnum(coupled("primary_variable")),
    _dvar(_dictator.isPorousFlowVariable(_prime_varnum) ? _dictator.porousFlowVariableNum(_prime_varnum) : 0),

    _pp_nodal_var(coupledNodalValue("porepressure")),
    _pp_qp_var(coupledValue("porepressure")),
    _gradpp_qp_var(coupledGradient("porepressure")),

    _sat_nodal_var(coupledNodalValue("saturation")),
    _sat_qp_var(coupledValue("saturation")),
    _gradsat_qp_var(coupledGradient("saturation"))
{
  if (_dictator.numPhases() != 1)
    mooseError("The Dictator proclaims that the number of phases is " << _dictator.numPhases() << " whereas PorousFlow1PhasePorS_VG can only be used for 1-phase simulations.  Be aware that the Dictator has noted your mistake.");
}

void
PorousFlow1PhasePorS_VG::initQpStatefulProperties()
{
  PorousFlowVariableBase::initQpStatefulProperties();

  buildPS();
}

void
PorousFlow1PhasePorS_VG::computeProperties()
{
  unsigned num_pp_qps = 0;
  unsigned num_sat_qps = 0;
  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    if (_decoder_nodal_var[_node_number[_qp]] >= 0)
      num_pp_qps += 1;
    else
      num_sat_qps += 1;
  }
  if (num_pp_qps == _qrule->n_points())
    _entire_elem_decoder = 1;
  else if (num_sat_qps == _qrule->n_points())
    _entire_elem_decoder = -1;
  else
    _entire_elem_decoder = 0;

  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
    computeQpProperties();
}
  

void
PorousFlow1PhasePorS_VG::computeQpProperties()
{
  PorousFlowVariableBase::computeQpProperties();

  buildPS();

  if (_dictator.notPorousFlowVariable(_prime_varnum))
    return;

  if (_decoder_nodal_var[_node_number[_qp]] >= 0)
  {
    // the node that stores its info at this _qp is a porepressure node
    _dporepressure_nodal_dvar[_qp][0][_dvar] = 1.0;
    _dsaturation_nodal_dvar[_qp][0][_dvar] = PorousFlowEffectiveSaturationVG::dseff(_porepressure_nodal[_qp][0], _al, _m);
  }
  else
  {
    // the node that stores its info at this _qp is a saturation node
    _dsaturation_nodal_dvar[_qp][0][_dvar] = 1.0;
    _dporepressure_nodal_dvar[_qp][0][_dvar] = - PorousFlowCapillaryPressureVG::dCapillaryPressure(_saturation_nodal[_qp][0], _m, 0.0, 1.0, _p0, std::numeric_limits<Real>::max()); // minus sign comes because porepressure = -capillaryPressure
  }

  // Now do the quad-point derivatives, which are only possible
  // if all the nodes in the element are of the same type
  if (_entire_elem_decoder == 1)
  {
    // all the relevant nodes in this element are porepressure nodes
    _porepressure_qp[_qp][0] = _pp_qp_var[_qp];
    _gradp_qp[_qp][0] = _gradpp_qp_var[_qp];
    _dporepressure_qp_dvar[_qp][0][_dvar] = 1.0;
    _dgradp_qp_dgradv[_qp][0][_dvar] = 1.0;
    _dgradp_qp_dv[_qp][0][_dvar] = RealGradient();
    _saturation_qp[_qp][0] = PorousFlowEffectiveSaturationVG::seff(_porepressure_qp[_qp][0], _al, _m);
    _grads_qp[_qp][0] = _gradp_qp[_qp][0] * PorousFlowEffectiveSaturationVG::dseff(_porepressure_qp[_qp][0], _al, _m);
    _dsaturation_qp_dvar[_qp][0][_dvar] = PorousFlowEffectiveSaturationVG::dseff(_porepressure_qp[_qp][0], _al, _m);
    _dgrads_qp_dgradv[_qp][0][_dvar] = PorousFlowEffectiveSaturationVG::dseff(_porepressure_qp[_qp][0], _al, _m);
    _dgrads_qp_dv[_qp][0][_dvar] = _gradp_qp[_qp][0] * PorousFlowEffectiveSaturationVG::d2seff(_porepressure_qp[_qp][0], _al, _m);
  }
  else if (_entire_elem_decoder == -1)
  {
    // all the relevant nodes in this element are saturation nodes
    _saturation_qp[_qp][0] = _sat_qp_var[_qp];
    _grads_qp[_qp][0] = _gradsat_qp_var[_qp];
    _dsaturation_qp_dvar[_qp][0][_dvar] = 1.0;
    _dgrads_qp_dgradv[_qp][0][_dvar] = 1.0;
    _dgrads_qp_dv[_qp][0][_dvar] = RealGradient();
    // minus signs in following come because porepressure = -capillaryPressure
    _porepressure_qp[_qp][0] = - PorousFlowCapillaryPressureVG::capillaryPressure(_saturation_qp[_qp][0], _m, 0.0, 1.0, _p0, std::numeric_limits<Real>::max());
    _gradp_qp[_qp][0] = - _grads_qp[_qp][0] * PorousFlowCapillaryPressureVG::dCapillaryPressure(_saturation_qp[_qp][0], _m, 0.0, 1.0, _p0, std::numeric_limits<Real>::max());
    _dporepressure_qp_dvar[_qp][0][_dvar] = - PorousFlowCapillaryPressureVG::dCapillaryPressure(_saturation_qp[_qp][0], _m, 0.0, 1.0, _p0, std::numeric_limits<Real>::max());
    _dgradp_qp_dgradv[_qp][0][_dvar] = - PorousFlowCapillaryPressureVG::dCapillaryPressure(_saturation_qp[_qp][0], _m, 0.0, 1.0, _p0, std::numeric_limits<Real>::max());
    _dgradp_qp_dv[_qp][0][_dvar] = - _grads_qp[_qp][0] * PorousFlowCapillaryPressureVG::d2CapillaryPressure(_saturation_qp[_qp][0], _m, 0.0, 1.0, _p0, std::numeric_limits<Real>::max());
  }
  else
  {
    // the _qp versions below are just the defaults, which are as good
    // as I can get for "mixed" elements with some nodes being porepressure
    // and some nodes being saturation.
    // I can do better for "uniform" elements.  Say all nodes are saturation
    // nodes.  Then saturation varies linearly across the element, but
    // obviously porepressure doesn't.  But if i just use the _qp versions
    // below then i'll get linear variation in porepressure.  So, in the
    // "if" and "else if" above i address this by using the capillary equations
    _porepressure_qp[_qp][0] = _pp_qp_var[_qp];
    _gradp_qp[_qp][0] = _gradpp_qp_var[_qp];
    _saturation_qp[_qp][0] = _sat_qp_var[_qp];
    _grads_qp[_qp][0] = _gradsat_qp_var[_qp];
  }
}

void
PorousFlow1PhasePorS_VG::buildPS()
{
  _porepressure_nodal[_qp][0] = _pp_nodal_var[_node_number[_qp]];
  _saturation_nodal[_qp][0] = _sat_nodal_var[_node_number[_qp]];
  _porepressure_qp[_qp][0] = _pp_qp_var[_qp];
  _gradp_qp[_qp][0] = _gradpp_qp_var[_qp];
  _saturation_qp[_qp][0] = _sat_qp_var[_qp];
  _grads_qp[_qp][0] = _gradsat_qp_var[_qp];
}
