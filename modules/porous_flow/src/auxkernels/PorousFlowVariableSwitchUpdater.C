/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlowVariableSwitchUpdater.h"
#include "PorousFlowEffectiveSaturationVG.h"
#include "PorousFlowCapillaryPressureVG.h"

template<>
InputParameters validParams<PorousFlowVariableSwitchUpdater>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addRequiredRangeCheckedParam<Real>("al", "al > 0", "van-Genuchten alpha parameter.  Must be positive.  effectiveSaturation = (1 + (-al*c)^(1/(1-m)))^(-m)");
  params.addRequiredRangeCheckedParam<Real>("m", "m > 0 & m < 1", "van-Genuchten m parameter.  Must be between 0 and 1, and optimally should be set to >0.5   EffectiveSaturation = (1 + (-al*p)^(1/(1-m)))^(-m)");
  params.addRequiredParam<UserObjectName>("PorousFlowDictator", "The UserObject that holds the list of PorousFlow variable names");
  params.addRequiredCoupledVar("p_or_s_variable", "The Variable that represents porepressure or saturation");
  params.addRequiredParam<bool>("updating_porepressure", "Whether the variable for the AuxKernel is porepressure");
  params.addRequiredCoupledVar("encoder_auxvar", "Updater Aux Variable");
  return params;
}

PorousFlowVariableSwitchUpdater::PorousFlowVariableSwitchUpdater(const InputParameters & parameters) :
    AuxKernel(parameters),
    _al(getParam<Real>("al")),
    _m(getParam<Real>("m")),
    _p0(1.0/_al),
    _dictator(getUserObject<PorousFlowDictator>("PorousFlowDictator")),
    _p_or_s(coupledValue("p_or_s_variable")),
    _updating_pp(getParam<bool>("updating_porepressure")),
    _encoder(coupledValue("encoder_auxvar"))
{
  if (_dictator.numPhases() != 1)
    mooseError("PorousFlowVariableSwitchUpdater: num phases must be 1\n");
}

Real
PorousFlowVariableSwitchUpdater::computeValue()
{
  if (_updating_pp)
  {
    if (_encoder[_qp] >= 0) // node is porepressure
    {
      return _p_or_s[_qp];
    }
    else // node is saturation
    {
      return - PorousFlowCapillaryPressureVG::capillaryPressure(_p_or_s[_qp], _m, 0.0, 1.0, _p0, std::numeric_limits<Real>::max()); // minus sign comes because porepressure = -capillaryPressure
    }
  }
  else // updating saturation variable
  {
    if (_encoder[_qp] >= 0) // node is porepressure
    {
      return PorousFlowEffectiveSaturationVG::seff(_p_or_s[_qp], _al, _m);
    }
    else // node is saturation
    {
      return _p_or_s[_qp];
    }
  }
}
