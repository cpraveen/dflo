#ifndef __INTEGRATOR_H__
#define __INTEGRATOR_H__

/*
 *  integrator.h
 *  
 *
 *  Created by Praveen Chandrashekar on 25/02/11.
 *  Copyright 2011 TIFR-CAM. All rights reserved.
 *
 */

#include <dofs/dof_handler.h>

#include <lac/vector.h>

#include <meshworker/dof_info.h>
#include <meshworker/integration_info.h>
#include <meshworker/simple.h>
#include <meshworker/loop.h>

//------------------------------------------------------------------------------
// Class for integrating rhs using MeshWorker
// Explicit scheme: we only need rhs and there is no system matrix
//------------------------------------------------------------------------------
template <int dim>
class IntegratorExplicit
{
public:
   IntegratorExplicit (const dealii::DoFHandler<dim>& dof_handler)
      : 
      dof_info (dof_handler) {};
   
   dealii::MeshWorker::IntegrationInfoBox<dim> info_box;
   dealii::MeshWorker::DoFInfo<dim> dof_info;
   dealii::MeshWorker::Assembler::ResidualSimple< dealii::parallel::distributed::Vector<double> > assembler;

};

#endif
