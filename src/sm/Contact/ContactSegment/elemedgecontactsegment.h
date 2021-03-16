/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */



#pragma once
#include "boundarycontactsegment.h"
#include "contactsegment.h"
#include "node.h"
#include "inputrecord.h"
#include "Elements/structuralelement.h"
#include "set.h"
#include "feinterpol.h"
#include "feinterpol2d.h"
#include "Elements/nlstructuralelement.h"
#include "Materials/structuralmaterial.h"
#include "intarray.h"

#define _IFT_Linear2dElementEdgeContactSegment_Name "linear2delementedgecontactsegment"
 //#define _IFT_ElementEdgeContactSegment_elemSet "elemset"

namespace oofem {

	class Linear2dElementEdgeContactSegment : public BoundaryContactSegment
	{
	public:
		Linear2dElementEdgeContactSegment(int n, Domain *aDomain) : BoundaryContactSegment(n, aDomain) { ; }
		~Linear2dElementEdgeContactSegment() {};

		IRResultType initializeFrom(InputRecord * ir) override;
		//returns normalized n, which is an normal vector of contact
		void computeNormal(FloatArray& answer, Node * node, TimeStep* tstep) override;
		//returns non-normalized tangent vector at the point of contact
		void computeTangent(FloatArray &answer, Node *node, TimeStep *tstep) override;
		//returns an extended N (aka A) matrix, integrated at point of contact of given node
		void computeSegmentNMatrix(FloatMatrix& answer, Node* node, TimeStep * tStep) override;
		//@todo: the following two functions should be merged
		void computeSegmentBMatrix(FloatMatrix &answer, Node *node, TimeStep *tStep) override;
		void computedNdksi(FloatMatrix &answer, Node *node, TimeStep *tStep);
		//computes the penetration of the slave node to the given(closest) segment
		double computePenetration(Node * node, TimeStep * tStep) override;
		// computes metric tensor, m_ij = t_i \cdot t_j, for 2d linear element, this is equivalent to l^2
		void computeMetricTensor(FloatMatrix &answer, Node *node, TimeStep *tStep) override;
		//return location array of the slave node and master segment
		void giveLocationArray(const IntArray& dofIdArray, IntArray& s_loc, const UnknownNumberingScheme& c_s) const override;
		//this function is needed to allocate appropriate components of stiffness matrix
		void giveLocationArrays(const IntArray& dofIdArray, IntArray& s_loc, const UnknownNumberingScheme& c_s) override;

		const char *giveClassName() const override { return "Linear2delementedgecontactsegment"; }
		const char *giveInputRecordName() const override { return _IFT_Linear2dElementEdgeContactSegment_Name; }

	protected:


		//computes distance of a point given by a set of coordinates to a line given by two sets of coordinates
		//returns whether point of intersection is inbetween the line points
		bool computeContactPoint(FloatArray &ksi, Node* node, StructuralElement* element, int elemedge, TimeStep* tStep);


	};
} //end namespace oofem
