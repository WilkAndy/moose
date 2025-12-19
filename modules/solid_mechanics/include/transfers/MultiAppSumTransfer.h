//* This file is part of the MOOSE framework
//* https://mooseframework.inl.gov
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE
#include "MultiAppTransfer.h"
#include "InputParameters.h"

// libMesh
#include "libmesh/mesh_base.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/fe_type.h"

/// TODO: Doco
// NOTE: assume a single pair of apps (index 0)
class MultiAppSumTransfer : public MultiAppTransfer
{
public:
  static InputParameters validParams();

  /// Standard constructor
  MultiAppSumTransfer(const InputParameters & parameters);

  /// Execute the transfer (here: build the T matrix)
  void initialSetup() override;

  /// Execute the transfer (here: build the T matrix)
  void execute() override;

  /// Accessor for the computed transfer matrix
  const libMesh::DenseMatrix<Real> & getTransferMatrix() const { return _T; }

protected:
  /**
   * Build the transfer matrix T using the provided coarse and fine meshes and FE type.
   *
   * Rows = fine-mesh nodes
   * Cols = coarse-mesh nodes
   * T(i,j) = coarse LAGRANGE shape function for coarse node j, evaluated at fine node i.
   *
   * @param coarse_mesh   the coarse application mesh
   * @param fine_mesh     the fine application mesh
   * @param fe_type       FE type (family, order) used for coarse nodal LAGRANGE shapes
   * @param T_out         output matrix to fill
   */
  void buildTransferMatrix(const libMesh::MeshBase & coarse_mesh,
                           const libMesh::MeshBase & fine_mesh,
                           const libMesh::FEType & fe_type,
                           libMesh::DenseMatrix<libMesh::Number> & T_out);

protected:
  /// Coarse FE type (family must be LAGRANGE for nodal shape functions)
  libMesh::FEType _fe_type;

  /// The assembled transfer matrix
  libMesh::DenseMatrix<libMesh::Number> _T;

  /// If a fine node lies outside the coarse mesh, leave row as zeros (optional future: nearest element)
  bool _skip_outside_points;

  /// Verbose logging
  bool _verbose;

  /// Variables involved in transfer
  const VariableName _from_var_name; // coarse mesh
  const VariableName _to_var_name;   // fine mesh

private:
  // Helper: find the libMesh System containing `var_name` in an EquationSystems.
  std::pair<libMesh::System *, unsigned int>
  find_system_and_var(libMesh::EquationSystems & es, const std::string & var_name);
};
