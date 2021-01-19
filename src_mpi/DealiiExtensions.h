/*
 * DealiiExtensions.h
 *
 *  Created on: Dec 11, 2014
 *      Author: kraemer
 *
 * Modified on: March 25, 2016
 * 	    by: juanpablogallego
 */

#ifndef DEALIIEXTENSIONS_H_
#define DEALIIEXTENSIONS_H_

#include <vector>
#include <set>
#include <map>

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/table.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
//#include <deal.II/dofs/function_map.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/component_mask.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/distributed/tria.h>

/**
 * @short	Pair that includes face number and cell
 * @note	One pair to rule them all!
 */
template<size_t dim, size_t spacedim = dim>
using FaceCellPair = std::pair<typename dealii::parallel::distributed::Triangulation<dim>::cell_iterator, int >;


template<size_t dim, size_t spacedim = dim>
using FacePair = dealii::GridTools::PeriodicFacePair<
				    typename dealii::parallel::distributed::Triangulation<dim>::cell_iterator >;

/**
 * @short	A std::map that maps cells to face pairs
 * @note	The order of cells in the Face Pair is OF meaning
 */
template<size_t dim, size_t spacedim = dim>
using PeriodicCellMap = std::map< FaceCellPair<dim, spacedim>, FacePair<dim, spacedim> >;


/** @short
 * some extensions to deal.ii
 */
namespace DealIIExtensions {

//DEAL_II_NAMESPACE_OPEN

using namespace dealii;


/**
 *  @short 	Gathers cell pairs at a periodic boundary. This function starts at the coarsest level and
 *  		recursively visits the subcells of boundary cells. At the active level, cell pairs are added
 *  		to the cell_map.
 *  @param[in] cell_1 		cell at the one side of a periodic boundary
 *  @param[in] cell_2 		cell at the other side of a periodic boundary
 *  @param[in] face_nr_1	boundary face number of the cell_1
 *  @param[in] face_nr_2	boundary face number of the cell_2
 *  @param[out] cell_map	A map that stores cells and faces at the periodic boundary
 *  @param[in] face_orientation		see deal.ii Glossary
 *  @param[in] face_flip			see deal.ii Glossary
 *  @param[in] face_rotation		see deal.ii Glossary
 *  @note The implementation of this function is similar to dealii::make_periodicity_constraints
 */
template<typename DH>
void make_periodicity_map_dg(const typename DH::cell_iterator &cell_1,
			     const typename identity<typename DH::cell_iterator>::type &cell_2,
			     size_t face_nr_1,
			     size_t face_nr_2,
			     PeriodicCellMap<DH::dimension>& cell_map,
			     const bool face_orientation,
			     const bool face_flip,
			     const bool face_rotation);

/**
 * @short	High-level version of the first function, starting from PeriodicFacePairs
 * @param[in]	periodic_faces 	a vector of PeriodicFacePairs. They can be obtained by
 * 								calling collect_periodic_faces.
 * @param[out]	cell_map		A map that stores cells and faces at the periodic boundary
 */
template<typename DH>
void make_periodicity_map_dg(const std::vector<GridTools::PeriodicFacePair<typename DH::cell_iterator> > &periodic_faces,
			     PeriodicCellMap<DH::dimension>& cell_map);

/**
 * @short	Another high-level version of the first function, starting from boundary indicators
 * @param[in]	dof_handler 	a DoFHandler object
 * @param[in]	b_id1			boundary id at the first side of the periodic boundary
 * @param[in]	b_id2			boundary id at the second side of the periodic boundary
 * @param[out]	cell_map		A map that stores cells and faces at the periodic boundary
 * @note 	Creating periodic boundaries requires great care in the order of operations (at least for
 * 			Meshes of type parallel::distributed::Triangulation<dim>)
 * 			-# create unrefined mesh
 * 			-# call collect_periodic_faces<Mesh> to find out the periodic couplings at the coarsest level
 * 			-# call add_periodicity to create a layer of ghost nodes across the periodic boundary
 * 			-# call this function (make_periodicity_map_dg) to recursively find out the periodic couplings
 * 			   on the active levels
 */
template<typename DH>
void make_periodicity_map_dg(const DH &dof_handler,
			     size_t b_id1,
			     size_t b_id2,
			     const int direction,
			     PeriodicCellMap<DH::dimension>& cell_map);

//DEAL_II_NAMESPACE_CLOSE

} /* namespace DealIIExtensions */

#endif /* DEALIIEXTENSIONS_H_ */

