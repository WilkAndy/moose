//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "MultiAppSumTransfer.h"
#include "Moose.h"
#include "MultiApp.h"
#include "MooseMesh.h"

// libMesh
#include "libmesh/point.h"
#include "libmesh/elem.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/point_locator_tree.h"
#include "libmesh/equation_systems.h"
#include "libmesh/system.h"
#include "libmesh/parallel.h"

registerMooseObject("MooseApp", MultiAppSumTransfer);

InputParameters
MultiAppSumTransfer::validParams()
{
  InputParameters params = MultiAppTransfer::validParams();

  // FE family: we target nodal LAGRANGE shapes
  MooseEnum fe_family("LAGRANGE", "LAGRANGE");
  params.addParam<MooseEnum>("fe_family",
                             fe_family,
                             "FE family for coarse nodal basis (use LAGRANGE for nodal shapes).");

  // FE order: FIRST, SECOND, THIRD,...
  MooseEnum fe_order("CONSTANT FIRST SECOND THIRD FOURTH", "FIRST");
  params.addParam<MooseEnum>("fe_order",
                             fe_order,
                             "FE order used for coarse nodal LAGRANGE shapes.");

  params.addParam<bool>("skip_outside_points",
                        true,
                        "If true, rows corresponding to fine nodes outside coarse mesh remain zeros.");

  params.addParam<bool>("verbose",
                        false,
                        "If true, prints summary of T assembly.");
  params.addRequiredParam<VariableName>("from_variable", "The coarse-mesh variable to transfer.");
  params.addRequiredParam<VariableName>("to_variable", "The fine-mesh variable to populate.");
  return params;
}

MultiAppSumTransfer::MultiAppSumTransfer(const InputParameters & params)
  : MultiAppTransfer(params),
    _fe_type(libMesh::FEType(
               static_cast<int>(getParam<MooseEnum>("fe_order").getEnum<libMesh::Order>()),
               getParam<MooseEnum>("fe_family").getEnum<libMesh::FEFamily>())),
    _skip_outside_points(getParam<bool>("skip_outside_points")),
    _verbose(getParam<bool>("verbose")),
    _from_var_name(getParam<VariableName>("from_variable")),
    _to_var_name(getParam<VariableName>("to_variable"))
{
}

void
MultiAppSumTransfer::initialSetup()
{
  MultiAppTransfer::initialSetup();

  getAppInfo();

 if (_from_meshes.size() != 1 || _to_meshes.size() != 1)
    mooseError("MultiAppSumTransfer expects exactly one local 'from' mesh (the coarse mesh) and one local 'to' mesh (the fine mesh).  Found from_meshes.size()=", _from_meshes.size(), "to_meshes.size()=", _to_meshes.size());

 // Get libMesh meshes
 const libMesh::MeshBase & coarse_mesh = (_current_direction == TO_MULTIAPP ? _from_meshes[0]->getMesh() : _to_meshes[0]->getMesh());
 const libMesh::MeshBase & fine_mesh   = (_current_direction == TO_MULTIAPP ? _to_meshes[0]->getMesh() : _from_meshes[0]->getMesh());

  // Resize T appropriately (rows=fine nodes, cols=coarse nodes)
  // Oh how we love dense matrices :-)
  // I'M ASSUMING NODE IDS GO FROM ZERO TO N_NODES()-1 !!
  const std::size_t n_rows = fine_mesh.n_nodes();
  const std::size_t n_cols = coarse_mesh.n_nodes();
  if (n_rows < n_cols)
    mooseError("MultiAppSumTransfer: it's assumed you run this from the coarse mesh.  In your case, n_nodes on the fine mesh is less than on the coarse mesh");
  _T.resize(n_rows, n_cols);
  _T.zero();

  // Assemble
  buildTransferMatrix(coarse_mesh, fine_mesh, _fe_type, _T);

  if (_verbose)
    std::cerr << "MultiAppSumTransfer: built T with size " <<  n_rows << " * " << n_cols << " using FEType(family=" << static_cast<int>(_fe_type.family) << ", order=" << static_cast<int>(_fe_type.order) << ")\n";
}


std::pair<libMesh::System *, unsigned int>
MultiAppSumTransfer::find_system_and_var(libMesh::EquationSystems & es, const std::string & var_name)
{
  for (unsigned int i = 0; i < es.n_systems(); ++i)
  {
    libMesh::System & sys = es.get_system(i);
    try
    {
      const unsigned int vn = sys.variable_number(var_name);
      return { &sys, vn };
    }
    catch (const std::exception &) { /* continue */ }
  }

  // Not found: emit useful diagnostics
  std::ostringstream oss;
  oss << "Variable '" << var_name << "' not found in any system of the inspected app."
      << " Systems present:";
  for (unsigned int i = 0; i < es.n_systems(); ++i)
  {
    libMesh::System & sys = es.get_system(i);
    oss << "\n  - " << sys.name() << " (n_vars=" << sys.n_vars() << ")";
  }
  throw std::runtime_error(oss.str());
}

void
MultiAppSumTransfer::execute()
{
  // Refresh app info (ensures _from_es _to_es _from_meshes _to_meshes are valid)
  getAppInfo();

  const std::size_t T_rows = _T.m();
  const std::size_t T_cols = _T.n();

  const libMesh::MeshBase & from_mesh = _from_meshes[0]->getMesh();
  const libMesh::MeshBase & to_mesh   = _to_meshes[0]->getMesh();


  // --- Source system & data ---
  libMesh::System * from_sys          = nullptr;
  unsigned int      from_var_number   = libMesh::invalid_uint;
  try
  {
    std::tie(from_sys, from_var_number) = find_system_and_var(*_from_es[0], _from_var_name);
  }
  catch (const std::exception & e)
  {
    mooseError("From variable lookup failed for '", _from_var_name, "': ", e.what());
  }

  const libMesh::DofMap & from_dof_map = from_sys->get_dof_map();
  libMesh::NumericVector<libMesh::Number> * from_vec = from_sys->solution.get();
  //        from_sys->current_local_solution.get() ? from_sys->current_local_solution.get()
  //                                              : from_sys->solution.get();
  if (!from_vec)
  {    
    mooseError("From system for variable '", _from_var_name, "' has no solution vector.");
  }

  // Build u_from (indexed by from_mesh node ID). Size guarded by T_cols or T_rows as appropriate.
  const std::size_t from_size = (_current_direction == TO_MULTIAPP ? T_cols : T_rows);
  std::vector<Real> u_from(from_size, 0.0);

  for (auto node_it = from_mesh.nodes_begin(); node_it != from_mesh.nodes_end(); ++node_it)
  {
    const libMesh::Node * from_node = *node_it;
    const std::size_t index = from_node->id();

    if (index >= from_size)
      mooseError("MultiAppSumTransfer: from node id ", index,
		 " exceeds T size of ", from_size,
		 ". Ensure T was sized to max node id + 1.");

    std::vector<libMesh::dof_id_type> dofs;
    from_dof_map.dof_indices(from_node, dofs, from_var_number);

    // Expect exactly one nodal DOF for LAGRANGE scalar variables
    if (dofs.empty())
    {
      // Leave as zero (e.g., non-nodal variable, or variable not defined at this node)
      continue;
    }
    else if (dofs.size() > 1)
    {
      // For higher-order or vector variables ????
      // Here we take the first DOF, assuming scalar nodal LAGRANGE.
    }

    u_from[index] = from_vec->el(dofs[0]);
  }

  // --- Destination system & target vector ---
  libMesh::System * to_sys          = nullptr;
  unsigned int      to_var_number   = libMesh::invalid_uint;
  try
  {
    std::tie(to_sys, to_var_number) = find_system_and_var(*_to_es[0], _to_var_name);
  }
  catch (const std::exception & e)
  {
    mooseError("Destination variable lookup failed for '", _to_var_name, "': ", e.what());
  }

  // We will write into the "to" multiapp transfer vector for this variable

  libMesh::NumericVector<Real> & tgt_vec = ((_current_direction == TO_MULTIAPP) ? getTransferVector(/* local to app index */ 0, _to_var_name) : *to_sys->solution.get());

  if (tgt_vec.size() == 0)
  {
    const auto n_global = to_sys->n_dofs();
    const auto n_local  = to_sys->n_local_dofs();
    if (n_global == 0)
      mooseError("Destination system for variable '", _to_var_name, "' has zero DOFs. INITIAL stage?");
    tgt_vec.init(n_global, n_local, /*fast=*/false, libMesh::AUTOMATIC);
    tgt_vec.zero(); // ensure a clean slate before setting values
  }

  const libMesh::DofMap & to_dof_map = to_sys->get_dof_map();
  const std::size_t to_size = (_current_direction == TO_MULTIAPP ? T_rows : T_cols);

  // --- Apply transfer: u_fine(row) = dot(T[row, :], u_from[:]) ---
  for (auto node_it = to_mesh.nodes_begin(); node_it != to_mesh.nodes_end(); ++node_it)
  {
    const libMesh::Node * to_node = *node_it;
    std::size_t index = to_node->id();

    if (index >= to_size)
      mooseError("MultiAppSumTransfer: to node id ", index,
		 " exceeds T size of ", to_size,
		 ". Ensure T was sized to max node id + 1.");

    // Dot product over all from node IDs
    libMesh::Number sum = 0.0;
    if (_current_direction == TO_MULTIAPP)
    {
      for (std::size_t col = 0; col < T_cols; ++col)
        sum += _T(index, col) * u_from[col];
    }
    else
    {
      for (std::size_t row = 0; row < T_rows; ++row)
        sum += _T(row, index) * u_from[row];
    }

    // Place into destination DOF(s) for this to node & variable
    std::vector<libMesh::dof_id_type> dofs;
    to_dof_map.dof_indices(to_node, dofs, to_var_number);

    if (dofs.empty())
    {
      // Variable not defined at this node; skip
      continue;
    }
    else if (dofs.size() > 1)
    {
      // For higher-order variables, distribute appropriately.
      // For now, assign to the first nodal DOF.
    }

    tgt_vec.set(dofs[0], static_cast<Real>(sum));
  }

  // Finalize target vector assembly for this variable
  tgt_vec.close();
  to_sys->update();


  std::cerr << "Transferred " << _from_var_name << " to " << _to_var_name << "\n";
}


void
MultiAppSumTransfer::buildTransferMatrix(const libMesh::MeshBase & coarse_mesh,
                                         const libMesh::MeshBase & fine_mesh,
                                         const libMesh::FEType & fe_type,
                                         libMesh::DenseMatrix<Real> & T_out)
{
  // Dimension inferred from mesh
  const unsigned int dim = coarse_mesh.mesh_dimension();

  // Point locator on the coarse mesh
  libMesh::PointLocatorTree locator(coarse_mesh);
  if (!locator.initialized())
    locator.init();

  // FE for shape evaluation on coarse elements
  std::unique_ptr<libMesh::FEBase> fe(libMesh::FEBase::build(dim, fe_type));

  // We will reinit with a single reference point each time; no quadrature rule is strictly necessary,
  // but FEBase expects a qrule for internal storage.  So, attach a 1-point Gauss in each dim, then
  // overwrite the qpoints on each reinit.
  std::unique_ptr<libMesh::QBase> qrule = libMesh::QBase::build(libMesh::QGAUSS, dim, libMesh::Order(1));
  // needed? correct?  qrule->init(1);
  fe->attach_quadrature_rule(qrule.get());


  // libMesh requires requesting what to compute prior to the first reinit()
  (void)fe->get_phi();

  // Iterate fine nodes (rows)
  // I love loops and loops

  for (auto node_it = fine_mesh.nodes_begin(); node_it != fine_mesh.nodes_end(); ++node_it)
  {
    const libMesh::Node * fine_node = *node_it;
    std::size_t row = fine_node->id();
    const libMesh::Point  p         = *fine_node; // physical coordinate

    // Find containing coarse element
    const libMesh::Elem * coarse_elem = locator(p);

    if (!coarse_elem)
    {
      if (_skip_outside_points)
        continue; // leave row as zeros
      // Optional fallback: nearest element search / projection could be added here.
      // but i think we probably won't need it, and could even throw an error
      continue;
    }

    // Compute reference coordinates xi for p in coarse_elem
    libMesh::Point xi = libMesh::FEInterface::inverse_map(dim, fe_type, coarse_elem, p);

    // Prepare one q-point at xi
    std::vector<libMesh::Point> qpoints(1);
    qpoints[0] = xi;

    // Evaluate shapes at xi
    fe->reinit(coarse_elem, &qpoints);

    // fe->get_phi() returns vector< vector<Real> > with size n_shape_functions x n_qpoints
    const std::vector<std::vector<libMesh::Real>> & phi = fe->get_phi();

    // For LAGRANGE, local shape functions correspond to local element nodes
    const unsigned int n_loc_shapes = phi.size();
    const unsigned int n_qp         = (n_loc_shapes > 0) ? phi[0].size() : 0;
    if (n_qp == 0)
      continue;

    // Map local shapes onto global coarse node columns
    // NOTE: For LAGRANGE of order >= 2, elem->n_nodes() includes edge/face/internal nodes.
    // We assume libMesh's local LAGRANGE shape ordering aligns with elem local node ordering.
    // If your version differs, adapt this loop accordingly.
    const unsigned int n_loc_nodes = coarse_elem->n_nodes();
    for (unsigned int j = 0; j < n_loc_nodes; ++j)
    {
      const libMesh::Node * coarse_node = coarse_elem->node_ptr(j);
      const libMesh::dof_id_type col    = coarse_node->id();

      // Shape value at our single qpoint
      const libMesh::Real val = phi[j][0];

      // Assign into T(row, col)
      T_out(row, col) = val;
    }

    // Note that for things that are not LAGRANGE, could have n_loc_shapes > n_loc_nodes
    // At the moment, restrict to nodal LAGRANGE by design.
  }
}
