// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



// <h1> Systems Example 4 - Linear Elastic Cantilever </h1>
// \author David Knezevic
// \date 2012
//
// In this example we model a homogeneous isotropic cantilever
// using the equations of linear elasticity. We set the Poisson ratio to
// \nu = 0.3 and clamp the left boundary and apply a vertical load at the
// right boundary.


// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/perf_log.h"
#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/tecplot_io.h"

#include "libmesh/analytic_function.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;
Real exact_solution2 (const Real x,
                     const Real y,
                     const Real z = 0.)
{
  static const Real pi = acos(-1.);

  return sin(.5*pi*x)*cos(.5*pi*y)*cos(.5*pi*z);
}

Real exact_solution (const Real x,
                     const Real y,
                     const Real z = 0.)
{
  static const Real pi = acos(-1.);

  return cos(.5*pi*x)*sin(.5*pi*y)*cos(.5*pi*z);
}

void exact_solution_wrapper (DenseVector<Number> & output,
                             const Point & p,
                             const Real)
{
      output(0) = exact_solution(p(0), p(1), p(2));
}


// Matrix and right-hand side assemble
void assemble_elasticity(EquationSystems & es,
                         const std::string & system_name);

// Define the elasticity tensor, which is a fourth-order tensor
// i.e. it has four indices i, j, k, l
// Begin the main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh and any dependent libaries
  LibMeshInit init (argc, argv);

  // Initialize the cantilever mesh
  const unsigned int dim = 2;

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(dim <= LIBMESH_DIM, "2D support");

  // Create a 2D mesh distributed across the default MPI communicator.
  Mesh mesh(init.comm(), dim);
  MeshTools::Generation::build_square (mesh,
                                       15, 15,
                                       0., 1.,
                                       -1., 1.,
                                       QUAD9);


  // Print information about the mesh to the screen.
  mesh.print_info();

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system and its variables.
  // Create a system named "Elasticity"
  LinearImplicitSystem & system =
    equation_systems.add_system<LinearImplicitSystem> ("Elasticity");

  // Add two displacement variables, u and v, to the system
  unsigned int u_var = system.add_variable("u", SECOND, LAGRANGE);
  unsigned int v_var = system.add_variable("v", SECOND, LAGRANGE);

  system.attach_assemble_function (assemble_elasticity);
  std::set<boundary_id_type> boundary_ids;
    // the dim==1 mesh has two boundaries with IDs 0 and 1
  boundary_ids.insert(0);
  boundary_ids.insert(1);
  boundary_ids.insert(2);
  boundary_ids.insert(3);
  
  std::vector<unsigned int> variables(2);
  variables[0] = u_var; variables[1] = v_var;
  AnalyticFunction<> exact_solution_object(exact_solution_wrapper);

  DirichletBoundary dirichlet_bc(boundary_ids,
                                 variables,
                                 &exact_solution_object);
 
//  system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);

  equation_systems.init();

  // Print information about the system to the screen.
  equation_systems.print_info();

  // Solve the system
  system.solve();

  // Plot the solution
#ifdef LIBMESH_HAVE_EXODUS_API
  //ExodusII_IO (mesh).write_equation_systems("displacement.e", equation_systems);
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
  std::string tecplot_filename = "JV_coupled_system.plt";
  TecplotIO (mesh, true).write_equation_systems (tecplot_filename, equation_systems);

  // All done.
  return 0;
}


void assemble_elasticity(EquationSystems & es,
                         const std::string & system_name)
{
  libmesh_assert_equal_to (system_name, "Elasticity");

  const MeshBase & mesh = es.get_mesh();

  const unsigned int dim = mesh.mesh_dimension();
  Real pival = libMesh::pi;
  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>("Elasticity");

  const unsigned int u_var = system.variable_number ("u");
  const unsigned int v_var = system.variable_number ("v");

  const DofMap & dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(0);
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);
  

  UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));
  QGauss qface(dim-1, fe_type.default_quadrature_order());
  fe_face->attach_quadrature_rule (&qface);

  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();
  const std::vector<Point> & q_point = fe->get_xyz();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke),
    Kvu(Ke), Kvv(Ke);

  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe);

  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_v;

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      const Elem * elem = *el;

      dof_map.dof_indices (elem, dof_indices);
      dof_map.dof_indices (elem, dof_indices_u, u_var);
      dof_map.dof_indices (elem, dof_indices_v, v_var);

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_u_dofs = dof_indices_u.size();
      const unsigned int n_v_dofs = dof_indices_v.size();

      fe->reinit (elem);

      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
      Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);

      Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
      Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);

      Fu.reposition (u_var*n_u_dofs, n_u_dofs);
      Fv.reposition (v_var*n_u_dofs, n_v_dofs);

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          const Real x = q_point[qp](0);
          const Real y = q_point[qp](1);
          const Real eps = 1.e-3;

          for (unsigned int i=0; i<n_u_dofs; i++)
            for (unsigned int j=0; j<n_u_dofs; j++)
              {
                  Kuu(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp] )*x ; 
              }

          for (unsigned int i=0; i<n_u_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
              {
               Kuv(i,j) = 0.0;          
              }

          for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_u_dofs; j++)
              {
               Kvu(i,j) = 0.0;
              }

          for (unsigned int i=0; i<n_v_dofs; i++)
            for (unsigned int j=0; j<n_v_dofs; j++)
              {
               Kvv(i,j) +=JxW[qp]*(dphi[i][qp]*dphi[j][qp])*x ;            
              }
             Real fxyz = exact_solution(x,y)*pival*pival/2.0 + pival*0.5*sin(0.5*pival*x)*sin(0.5*pival*y)/x ;
             const Real fxy = -(exact_solution(x, y-eps) +
                        exact_solution(x, y+eps) +
                        exact_solution(x-eps, y) +
                        exact_solution(x+eps, y) -
                        4.*exact_solution(x, y))/eps/eps;


                  for (unsigned int i=0; i<n_u_dofs; i++)
                  Fu(i) += JxW[qp]*fxyz*phi[i][qp]*x;
                  for (unsigned int i = 0;i<n_v_dofs;i++)
                  Fv(i) += JxW[qp]*fxyz*phi[i][qp]*x;
        }
 	
        for (unsigned int s=0; s<elem->n_sides(); s++)
          if (elem->neighbor(s) == libmesh_nullptr)
            {
	const std::vector<std::vector<Real> > & phi_face = fe_face->get_phi();
	const std::vector<Real> & JxW_face = fe_face->get_JxW();
	const std::vector<Point> & qface_point = fe_face->get_xyz();
	const std::vector<Point>& qface_normals = fe_face->get_normals();
	fe_face->reinit(elem, s);
    	UniquePtr<Elem> side (elem->build_side(s));
	/*	  for (unsigned int qp=0; qp<qface.n_points(); qp++)
                {
                  // The location on the boundary of the current
                  // face quadrature point.
                  const Real xf = qface_point[qp](0);
                  const Real yf = qface_point[qp](1);

                  // The penalty value.  \frac{1}{\epsilon}
                  // in the discussion above.
                  const Real penalty = 1.e10;

                  // The boundary value.
                  const Real value = exact_solution(xf, yf);

                  // Matrix contribution of the L2 projection.
                  for (unsigned int i=0; i<phi_face.size(); i++)
                    for (unsigned int j=0; j<phi_face.size(); j++)
                      Ke(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp]*xf;

                  // Right-hand-side contribution of the L2
                  // projection.
                  for (unsigned int i=0; i<phi_face.size(); i++)
                    Fe(i) += JxW_face[qp]*penalty*value*phi_face[i][qp]*xf;
                }*/
      	double check = 0;
    		    {
    		      for (unsigned int ns=0; ns<side->n_nodes(); ns++)
    			{ const Real penalty = 1.e10;
    			  for (unsigned int n=0; n<elem->n_nodes(); n++)
    			    if (elem->node(n) == side->node(ns))
    			      { 
				Node *node = elem->get_node(n);

				Point poi = *node;
				const Real xf = poi(0);
		                const Real yf = poi(1);

    				for(unsigned int j=0; j<n_u_dofs; ++j)
    				  Kuu(n,j) = 0.;
                                for(unsigned int j=0; j<n_v_dofs; ++j)
                                  Kvv(n,j) = 0.;
				
    				Kuu(n,n) = 1.;
                                Kvv(n,n) = 1.;

    				Fu(n)   = exact_solution(xf,yf);			
                                Fv(n)   = exact_solution(xf,yf);
			      }
    			}
      		    }
            }
      
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
    }
}

