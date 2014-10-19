#include "ic.h"
#include "equation.h"
#include "claw.h"

using namespace dealii;

//--------------------------------------------------------------------------------------------
// Initial condition for Rayleigh-Taylor problem
// This is setup for 2-d case only
//--------------------------------------------------------------------------------------------
template <int dim>
void RayleighTaylor<dim>::vector_value (const Point<dim> &p,
                                        Vector<double>   &values) const
{
   // Density
   if(p[1] < 0.0)
      values[EulerEquations<dim>::density_component] = 1.0;
   else
      values[EulerEquations<dim>::density_component] = 2.0;
   
   // Momentum
   for(unsigned int d=0; d<dim; ++d)
      values[d] = 0.0;
   
   double vel = A *
                (1.0 + std::cos(2.0*numbers::PI*p[0]/Lx))/2.0 *
                (1.0 + std::cos(2.0*numbers::PI*p[1]/Ly))/2.0;
   
   values[1] = values[EulerEquations<dim>::density_component] * vel;
   
   double pressure = P0 - gravity * values[EulerEquations<dim>::density_component] * p[1];

   // Energy
   values[EulerEquations<dim>::energy_component] =
      pressure/(EulerEquations<dim>::gas_gamma - 1.0)
      + 0.5 * values[EulerEquations<dim>::density_component] * vel * vel;
}

//--------------------------------------------------------------------------------------------
// Initial condition for isentropic vortex problem
// This is setup for 2-d case only
//--------------------------------------------------------------------------------------------
template <int dim>
void IsentropicVortex<dim>::vector_value (const Point<dim> &p,
                                          Vector<double>   &values) const	
{
   const double gamma = EulerEquations<dim>::gas_gamma;
   
   double r2   = (p[0]-x0)*(p[0]-x0) + (p[1]-y0)*(p[1]-y0);
   
   double rho = std::pow(1.0 - a2*exp(1.0-r2), 1.0/(gamma-1.0));
   double vex =  - a1 * (p[1]-y0) * std::exp(0.5*(1.0-r2));
   double vey =  + a1 * (p[0]-x0) * std::exp(0.5*(1.0-r2));
   double pre = std::pow(rho, gamma) / gamma;
   
   values[0] = rho * vex;
   values[1] = rho * vey;
   values[EulerEquations<dim>::density_component] = rho;
   values[EulerEquations<dim>::energy_component] = pre/(gamma-1.0)
                                                   + 0.5 * rho * (vex*vex + vey*vey);
}

//--------------------------------------------------------------------------------------------
// Initial condition for system of three isentropic vortices.
// This is setup for 2-d case only
//--------------------------------------------------------------------------------------------
template <int dim>
void VortexSystem<dim>::vector_value (const Point<dim> &p,
                                      Vector<double>   &values) const
{
   const double gamma = EulerEquations<dim>::gas_gamma;
   
   double rho = 0, vex = 0, vey = 0, pre = 0;
   for(unsigned int i=0; i<3; ++i)
   {
      double r2   = (p[0]-x[i])*(p[0]-x[i]) + (p[1]-y[i])*(p[1]-y[i]);
      
      rho += std::pow(1.0 - a2*exp(1.0-r2), 1.0/(gamma-1.0));
      vex +=  - a1 * (p[1]-y[i]) * std::exp(0.5*(1.0-r2));
      vey +=  + a1 * (p[0]-x[i]) * std::exp(0.5*(1.0-r2));
   }
   
   rho -= 2.0;
   vex /= 3.0;
   vey /= 3.0;
   pre = std::pow(rho, gamma) / gamma;
   
   // Put large pressure in the region { |x| < 0.1 and |y| < 0.1 }
   if(std::fabs(p[0]) < 0.1 && std::fabs(p[1]) < 0.1) pre = 50.0;
   
   values[0] = rho * vex;
   values[1] = rho * vey;
   values[EulerEquations<dim>::density_component] = rho;
   values[EulerEquations<dim>::energy_component] = pre/(gamma-1.0)
                                                 + 0.5 * rho * (vex*vex + vey*vey);
}


//------------------------------------------------------------------------------
// Sets initial condition based on input file.
// For Qk basis we can just do interpolation.
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::set_initial_condition_Qk ()
{
   if(parameters.ic_function == "rt")
      VectorTools::interpolate(mapping(), dof_handler,
                               RayleighTaylor<dim>(parameters.gravity), old_solution);
   else if(parameters.ic_function == "isenvort")
      VectorTools::interpolate(mapping(), dof_handler,
                               IsentropicVortex<dim>(5.0, 0.0, 0.0), old_solution);
   else if(parameters.ic_function == "vortsys")
      VectorTools::interpolate(mapping(), dof_handler,
                               VortexSystem<dim>(), old_solution);
   else
      VectorTools::interpolate(mapping(), dof_handler,
                               parameters.initial_conditions, old_solution);
   
   current_solution = old_solution;
   predictor = old_solution;
}

//------------------------------------------------------------------------------
// Sets initial condition based on input file.
// For Pk basis we have to do an L2 projection.
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::set_initial_condition_Pk ()
{
   //QGauss<dim> quadrature_formula (fe.degree+1);
   //const unsigned int n_q_points = quadrature_formula.size();
   
   //FEValues<dim> fe_values (fe, quadrature_formula,
                            //update_values | update_q_points | update_JxW_values);
   
   //// Multiply by inverse mass matrix
   //const unsigned int   dofs_per_cell = fe.dofs_per_cell;
   //std::vector<unsigned int> dof_indices(fe.dofs_per_cell);
   //std::vector<unsigned int> local_dof_indices (dofs_per_cell);
   //Vector<double> initial_values (n_q_points);
   //Vector<double> cell_vector(dofs_per_cell);
   
   //typename DoFHandler<dim>::active_cell_iterator
      //cell = dof_handler.begin_active(),
      //endc = dof_handler.end();
      
   //if(parameters.ic_function == "rt")
   //{
	   //RayleighTaylor<dim> initial_condition(parameters.gravity);
	   //for (; cell!=endc; ++cell)
	   //if(cell->is_locally_owned())
	   //{
			//cell->get_dof_indices(local_dof_indices);
			//fe_values.reinit (cell);
			//const std::vector<Point<dim> > q_points(quadrature_formula.get_points());

			//old_solution = 0;
			//unsigned int c = cell_number(cell);
      
			//for(unsigned int i=0; i<dofs_per_cell; ++i)
			//{
				//for(unsigned int q=0; q<n_q_points; ++q)
				//{
					//initial_condition.vector_value (q_points[q], initial_values);
					//old_solution(dof_indices[i]) += old_solution(dof_indices[i]) * initial_values *
													//fe_values.shape_value(i,q) *
													//fe_values.JxW(q);									
				//}
				//old_solution(dof_indices[i]) *= inv_mass_matrix[c][i];
			//}
      
		//}
   //}

   //else if(parameters.ic_function == "isenvort")
   //{
	   //IsentropicVortex<dim> initial_condition(5.0, 0.0, 0.0);
	   //for (; cell!=endc; ++cell)
	   //if(cell->is_locally_owned())
	   //{
			//cell->get_dof_indices(local_dof_indices);
			//fe_values.reinit (cell);
      		//const std::vector<Point<dim> > q_points(quadrature_formula.get_points());

			//old_solution = 0;
			//unsigned int c = cell_number(cell);
      
			//for(unsigned int i=0; i<dofs_per_cell; ++i)
			//{
				//for(unsigned int q=0; q<n_q_points; ++q)
				//{
					//initial_condition.vector_value (q_points[q], initial_values);
					//old_solution(dof_indices[i]) += initial_values[q] *
													//fe_values.shape_value(i,q) *
													//fe_values.JxW(q);									
				//}
				//old_solution(dof_indices[i]) *= inv_mass_matrix[c][i];
			//}
      
		//}
   //}

   //else if(parameters.ic_function == "vortsys")
   //{
	   //VortexSystem<dim> initial_condition;
	   //for (; cell!=endc; ++cell)
	   //if(cell->is_locally_owned())
	   //{
			//cell->get_dof_indices(local_dof_indices);
			//fe_values.reinit (cell);
      		//const std::vector<Point<dim> > q_points(quadrature_formula.get_points());

			//old_solution = 0;
			//unsigned int c = cell_number(cell);
      
			//for(unsigned int i=0; i<dofs_per_cell; ++i)
			//{
				//for(unsigned int q=0; q<n_q_points; ++q)
				//{
					//initial_condition.vector_value (q_points[q], initial_values);
					//old_solution(dof_indices[i]) += initial_values[q] *
													//fe_values.shape_value(i,q) *
													//fe_values.JxW(q);									
				//}
				//old_solution(dof_indices[i]) *= inv_mass_matrix[c][i];
			//}
      
		//}
   //}

   //else
   //{
	   //dealii::FunctionParser<dim> initial_condition(parameters.initial_conditions);
	   //for (; cell!=endc; ++cell)
	   //if(cell->is_locally_owned())
	   //{
			//cell->get_dof_indices(local_dof_indices);
			//fe_values.reinit (cell);
      		//const std::vector<Point<dim> > q_points(quadrature_formula.get_points());
            
			//old_solution = 0;
			//unsigned int c = cell_number(cell);
      
			//for(unsigned int i=0; i<dofs_per_cell; ++i)
			//{
				//for(unsigned int q=0; q<n_q_points; ++q)
				//{
					//initial_condition.vector_value (q_points[q], initial_values);
					//old_solution(dof_indices[i]) += initial_values[q] *
													//fe_values.shape_value(i,q) *
													//fe_values.JxW(q);									
				//}
				//old_solution(dof_indices[i]) *= inv_mass_matrix[c][i];
			//}
      
		//}
   //}
   Vector<double> p_current_solution;
   p_current_solution.reinit 	(dof_handler.n_locally_owned_dofs());
   p_current_solution=current_solution;
   
   if(parameters.ic_function == "rt")
      VectorTools::create_right_hand_side (mapping(), dof_handler,
                                           QGauss<dim>(fe.degree+1),
                                           RayleighTaylor<dim>(parameters.gravity),
                                           p_current_solution);
   else if(parameters.ic_function == "isenvort")
      VectorTools::create_right_hand_side (mapping(), dof_handler,
                                           QGauss<dim>(fe.degree+1),
                                           IsentropicVortex<dim>(5.0, 0.0, 0.0),
                                           p_current_solution);
   else if(parameters.ic_function == "vortsys")
      VectorTools::create_right_hand_side (mapping(), dof_handler,
                                           QGauss<dim>(fe.degree+1),
                                           VortexSystem<dim>(),
                                           p_current_solution);
   else
      VectorTools::create_right_hand_side (mapping(), dof_handler,
                                           QGauss<dim>(fe.degree+1),
                                           parameters.initial_conditions,
                                           p_current_solution);
   
   std::vector<unsigned int> dof_indices(fe.dofs_per_cell);
   typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
   
   for (; cell!=endc; ++cell)
   if(cell->is_locally_owned())
   {
      cell->get_dof_indices(dof_indices);
      unsigned int c = cell_number(cell);
      
      for(unsigned int i=0; i<fe.dofs_per_cell; ++i)
         old_solution(dof_indices[i]) = p_current_solution(dof_indices[i]) *
                                        inv_mass_matrix[c][i];
   }
   
   //current_solution = old_solution;
   //predictor = old_solution;
   
   old_solution.compress(VectorOperation::insert);
   current_solution = old_solution;
   predictor = old_solution;

}

//------------------------------------------------------------------------------
// Set intitial condition by interpolation or projection depending on
// type of basis function.
//------------------------------------------------------------------------------
template <int dim>
void ConservationLaw<dim>::set_initial_condition ()
{
   if(parameters.basis == Parameters::AllParameters<dim>::Qk)
      set_initial_condition_Qk();
   else
      set_initial_condition_Pk();
}

template class RayleighTaylor<2>;
template class IsentropicVortex<2>;
template class VortexSystem<2>;
template class ConservationLaw<2>;
