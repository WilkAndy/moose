/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "PorousFlow1PhaseSet.h"
#include "NonlinearSystem.h"

template<>
InputParameters validParams<PorousFlow1PhaseSet>()
{
  InputParameters params = validParams<NodalUserObject>();
  params.addRequiredCoupledVar("p_or_s_var", "Nonlinear Variable whose value is either porepressure or saturation, depending upon encoder");
  params.addRequiredParam<AuxVariableName>("encoder", "Aux Variable whose value at a node determines the meaning of p_or_s_var");
  params.addRequiredParam<AuxVariableName>("porepressure", "Aux Variable that is the porepressure");
  params.addRequiredParam<AuxVariableName>("saturation", "Aux Variable that is the saturation");
  return params;
}

PorousFlow1PhaseSet::PorousFlow1PhaseSet(const InputParameters & parameters) :
    NodalUserObject(parameters),
    _p_or_s_num(coupled("p_or_s_var")),
    _encoder_name(parameters.get<AuxVariableName>("encoder")),
    _porepressure_name(parameters.get<AuxVariableName>("porepressure")),
    _saturation_name(parameters.get<AuxVariableName>("saturation"))
{
}

void
PorousFlow1PhaseSet::execute()
{
  NonlinearSystem & nl = _fe_problem.getNonlinearSystem();
  unsigned int sys_num = nl.number();
  AuxiliarySystem & aux = _fe_problem.getAuxiliarySystem();
  unsigned int aux_num = aux.number();

  const Node & node = *_current_node;
  unsigned dof_num_for_p_or_s =  node.dof_number(sys_num, _p_or_s_num, 0);
  unsigned encoder_num = node.dof_number(aux_num, _fe_problem.getAuxiliarySystem().getVariable(_tid, _encoder_name).number(), 0);

  const VariableValue pp = _fe_problem.getAuxiliarySystem().getVariable(_tid, _porepressure_name).nodalSln();
  if (pp[_qp] > -1000)
  {
    // set the variable to "porepressure"
    const VariableValue ppold = _fe_problem.getAuxiliarySystem().getVariable(_tid, _porepressure_name).nodalSlnOld();
    const VariableValue ppolder = _fe_problem.getAuxiliarySystem().getVariable(_tid, _porepressure_name).nodalSlnOlder();
    Moose::out << "porepressure_node pp = " << pp[_qp] << " ppolder = " << ppolder[_qp] << "\n";
    nl.solution().set(dof_num_for_p_or_s, pp[_qp]);
    nl.solutionOld().set(dof_num_for_p_or_s, ppold[_qp]);
    nl.solutionOlder().set(dof_num_for_p_or_s, ppolder[_qp]);
    aux.solution().set(encoder_num, 1.0);
  }
  else
  {
    // set the variable to "saturation"
    const VariableValue sat = _fe_problem.getAuxiliarySystem().getVariable(_tid, _saturation_name).nodalSln();
    const VariableValue satold = _fe_problem.getAuxiliarySystem().getVariable(_tid, _saturation_name).nodalSlnOld();
    const VariableValue satolder = _fe_problem.getAuxiliarySystem().getVariable(_tid, _saturation_name).nodalSlnOlder();
    Moose::out << "saturation_node sat = " << sat[_qp] << " satolder = " << satolder[_qp] << "\n";
    nl.solution().set(dof_num_for_p_or_s, sat[_qp]);
    nl.solutionOld().set(dof_num_for_p_or_s, satold[_qp]);
    nl.solutionOlder().set(dof_num_for_p_or_s, satolder[_qp]);
    aux.solution().set(encoder_num, -1.0);
  }
}

void
PorousFlow1PhaseSet::finalize()
{
  Moose::out << "PF1PS finalise\n";
  NonlinearSystem& nl = _fe_problem.getNonlinearSystem();
  nl.solution().close();
  nl.solutionOld().close();
  nl.solutionOlder().close();
  AuxiliarySystem & aux = _fe_problem.getAuxiliarySystem();
  aux.solution().close();
}

