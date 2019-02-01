/*
 * DealiiExtensions.cpp
 *
 *  Created on: Dec 11, 2014
 *      Author: kraemer
 *
 * Modified on: March 25, 2016
 * 	    by: juanpablogallego
 */

#include "DealiiExtensions.h"

#include <deal.II/base/thread_management.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <algorithm>
#include <numeric>

using namespace dealii;
using namespace DoFTools;

namespace DealIIExtensions {

///////////////////////////////////////////
////////////// PERIODICITY ////////////////
///////////////////////////////////////////
// derived from make_periodicity constraints

template<typename DH>
void make_periodicity_map_dg(const typename DH::cell_iterator &cell_1,
			     const typename identity<typename DH::cell_iterator>::type &cell_2,
			     size_t face_nr_1,
			     size_t face_nr_2,
			     PeriodicCellMap<DH::dimension>& cell_map,
			     const bool face_orientation,
			     const bool face_flip,
			     const bool face_rotation)
{
	static const int dim = DH::dimension;
	static const int spacedim = DH::space_dimension;

	typedef typename DH::face_iterator FaceIterator;
	FaceIterator face_1 = cell_1->face(face_nr_1);
	FaceIterator face_2 = cell_2->face(face_nr_2);

	Assert(
			(dim != 1) || (face_orientation == true && face_flip == false && face_rotation == false),
			ExcMessage ("The supplied orientation " "(face_orientation, face_flip, face_rotation) " "is invalid for 1D"));

	Assert((dim != 2) || (face_orientation == true && face_rotation == false),
			ExcMessage ("The supplied orientation " "(face_orientation, face_flip, face_rotation) " "is invalid for 2D"));

	Assert(face_1 != face_2,
			ExcMessage ("face_1 and face_2 are equal! Cannot create boundary conditions DoFs " "on the very same face"));

	Assert(face_1->at_boundary() && face_2->at_boundary(),
			ExcMessage ("Faces for periodicity constraints must be on the " "boundary"));

	// A lookup table on how to go through the child cells depending on the
	// orientation:
	// see Documentation of GeometryInfo for details

	static const int lookup_table_2d[2][2] = {
	//          flip:
			{ 0, 1 }, //  false
			{ 1, 0 }, //  true
			};

	static const int lookup_table_3d[2][2][2][4] = {
	//                    orientation flip  rotation
			{ { { 0, 2, 1, 3 }, //  false       false false
					{ 2, 3, 0, 1 }, //  false       false true
					}, { { 3, 1, 2, 0 }, //  false       true  false
							{ 1, 0, 3, 2 }, //  false       true  true
					}, }, { { { 0, 1, 2, 3 }, //  true        false false
					{ 1, 3, 0, 2 }, //  true        false true
					}, { { 3, 2, 1, 0 }, //  true        true  false
							{ 2, 0, 3, 1 }, //  true        true  true
					}, }, };

	if (cell_1->has_children() && cell_2->has_children()) {
		// In the case that both faces have children, we loop over all
		// children and apply make_periodicty_constrains recursively:

		Assert(
				face_1->n_children() == GeometryInfo<dim>::max_children_per_face && face_2->n_children() == GeometryInfo<dim>::max_children_per_face,
				ExcNotImplemented());

		for (unsigned int i = 0; i < GeometryInfo<dim>::max_children_per_face;
				++i) {
			// Lookup the index for the second face
			unsigned int j;
			switch (dim) {
			case 2:
				j = lookup_table_2d[face_flip][i];
				break;
			case 3:
				j =
						lookup_table_3d[face_orientation][face_flip][face_rotation][i];
				break;
			default:
				AssertThrow(false, ExcNotImplemented())
				;
			}

			// find subcell ids that belong to the subface indices
			size_t child_cell_1 = GeometryInfo<dim>::child_cell_on_face(
					cell_1->refinement_case(), face_nr_1, i, face_orientation,
					face_flip, face_rotation, face_1->refinement_case());
			size_t child_cell_2 = GeometryInfo<dim>::child_cell_on_face(
					cell_2->refinement_case(), face_nr_2, j, face_orientation,
					face_flip, face_rotation, face_2->refinement_case());

			// precondition: subcell has the same orientation as cell (so that the face numbers coincide)
			// recursive call
			make_periodicity_map_dg<DH>(cell_1->child(child_cell_1),
					cell_2->child(child_cell_2), face_nr_1, face_nr_2, cell_map,
					face_orientation, face_flip, face_rotation);
		}
	} else {
		// Otherwise at least one of the two faces is active and
		// we need to do some work and enter the periodic face pairs!

		// This case could only be allowed with anisotropic refinement,
		// because the coarser cell is not allowed to have subfaces at the boundary.
		// Anisotropic refinement is also left out here for simplicity.
		// Consequently, opposite cells at periodic boundaries have to have the same
		// refinement level.
	/*	if (cell_1->has_children())
			PeriodicBoundaryNotPossible(
					"Opposite cells at periodic boundaries have to have the same refinement level.");
		if (cell_2->has_children())
			PeriodicBoundaryNotPossible(
					"Opposite cells at periodic boundaries have to have the same refinement level.");//*/

		// taken from make_periodicity_constraints: make sure faces are not artificial
		const unsigned int face_1_index = face_1->nth_active_fe_index(0);
		const unsigned int face_2_index = face_2->nth_active_fe_index(0);
		const unsigned int dofs_per_face = face_1->get_fe(face_1_index).dofs_per_face;
		std::vector<types::global_dof_index> dofs_1(dofs_per_face);
		std::vector<types::global_dof_index> dofs_2(dofs_per_face);
		face_1->get_dof_indices(dofs_1, face_1_index);
		face_2->get_dof_indices(dofs_2, face_2_index);
		for (unsigned int i = 0; i < dofs_per_face; i++) {
			if (dofs_1[i] == numbers::invalid_dof_index
					|| dofs_2[i] == numbers::invalid_dof_index) {
				/* If either of these faces have no indices, stop.  This is so
				 * that there is no attempt to match artificial cells of
				 * parallel distributed triangulations.
				 *
				 * While it seems like we ought to be able to avoid even calling
				 * set_periodicity_constraints for artificial faces, this
				 * situation can arise when a face that is being made periodic
				 * is only partially touched by the local subdomain.
				 * make_periodicity_constraints will be called recursively even
				 * for the section of the face that is not touched by the local
				 * subdomain.
				 *
				 * Until there is a better way to determine if the cells that
				 * neighbor a face are artificial, we simply test to see if the
				 * face does not have a valid dof initialization.
				 */
				return;
			}
		}

		// make sure cells are at least ghost cells (or even locally owned)
		if ((cell_1->has_children()) or (cell_2->has_children()))
			return;
		if ((cell_1->is_artificial()) or (cell_2->is_artificial()))
			return;
		/*Assert(not cell_1->is_artificial(),
				ExcMessage ("Cell at periodic boundary must not be artificial."));
		Assert(not cell_2->is_artificial(),
				ExcMessage ("Cell at periodic boundary must not be artificial."));
		 */

		// insert periodic face pair for both cells
		dealii::GridTools::PeriodicFacePair<
				dealii::TriaIterator<dealii::CellAccessor<dim, spacedim> > > face_pair1;//, face_pair2;
		face_pair1.cell[0] = cell_1;
		face_pair1.cell[1] = cell_2;
		face_pair1.face_idx[0] = face_nr_1;
		face_pair1.face_idx[1] = face_nr_2;
		face_pair1.orientation[0] = face_orientation;
		face_pair1.orientation[1] = face_flip;
		face_pair1.orientation[2] = face_rotation;
		
		/*
		 * 	We are interested in a map of the form
		 * 	(cell_1, face_1) --> FacePair
		 */
		
		FaceCellPair<dim,spacedim> Pair1 (cell_1,face_nr_1);
		cell_map.insert(std::pair<FaceCellPair<dim>,FacePair<dim,spacedim>>(Pair1,face_pair1));
		
		FaceCellPair<dim,spacedim> Pair2 (cell_2,face_nr_2);
		cell_map.insert(std::pair<FaceCellPair<dim>,FacePair<dim,spacedim>>(Pair2,face_pair1));
	}
}

template<typename DH>
void make_periodicity_map_dg(const std::vector<GridTools::PeriodicFacePair<typename DH::cell_iterator> > &periodic_faces,
			     PeriodicCellMap<DH::dimension>& cell_map)
{
	typedef std::vector<GridTools::PeriodicFacePair<typename DH::cell_iterator> > FaceVector;
	typename FaceVector::const_iterator it, end_periodic;
	it = periodic_faces.begin();
	end_periodic = periodic_faces.end();

	// Loop over all periodic faces...
	for (; it != end_periodic; ++it) {
		const typename DH::cell_iterator cell_1 = it->cell[0];
		size_t face_id_1 = it->face_idx[0];
		const typename DH::cell_iterator cell_2 = it->cell[1];
		size_t face_id_2 = it->face_idx[1];

		Assert(
				cell_1->face(face_id_1)->at_boundary() && cell_2->face(face_id_2)->at_boundary(),
				ExcInternalError());

		Assert(cell_1->face(face_id_1) != cell_2->face(face_id_2),
				ExcInternalError());

		// ... and apply the low level function to
		// every matching pair:
		make_periodicity_map_dg<DH>(cell_1, cell_2, face_id_1, face_id_2,
				cell_map, it->orientation[0], it->orientation[1],
				it->orientation[2]);
	}
}

// High level interface variants:

template<typename DH>
void make_periodicity_map_dg(const DH &dof_handler,
			     size_t b_id1,
			     size_t b_id2,
			     const int direction,
			     PeriodicCellMap<DH::dimension>& cell_map)
{
	static const int space_dim = DH::space_dimension;
	(void) space_dim;
	Assert(0<=direction && direction<space_dim,
			ExcIndexRange (direction, 0, space_dim));

	Assert(b_id1 != b_id2,
			ExcMessage ("The boundary indicators b_id1 and b_id2 must be" "different to denote different boundaries."));

	std::vector<GridTools::PeriodicFacePair<typename DH::cell_iterator> > matched_pairs;

	// Collect matching periodic cells on the coarsest level:
	GridTools::collect_periodic_faces(dof_handler, b_id1, b_id2, direction,
			matched_pairs);

	// call function to make the periodicity map
	make_periodicity_map_dg<DH>(matched_pairs, cell_map);
}

template
void make_periodicity_map_dg<DoFHandler<2> >(const typename DoFHandler<2>::cell_iterator &cell_1,
					     const identity<typename DoFHandler<2>::cell_iterator>::type &cell_2,
					     size_t face_nr_1,
					     size_t face_nr_2,
					     PeriodicCellMap<2>& cell_map,
					     const bool face_orientation,
					     const bool face_flip,
					     const bool face_rotation);

template
void make_periodicity_map_dg<DoFHandler<2> >(const std::vector<
						   GridTools::PeriodicFacePair<
						   typename DoFHandler<2>::cell_iterator> > &periodic_faces,
					     PeriodicCellMap<2>& cell_map);

template
void make_periodicity_map_dg<DoFHandler<2> >(const DoFHandler<2> &dof_handler,
					     size_t b_id1,
					     size_t b_id2,
					     const int direction,
					     PeriodicCellMap<2>& cell_map);

} /* namepace DealIIExtensions */
