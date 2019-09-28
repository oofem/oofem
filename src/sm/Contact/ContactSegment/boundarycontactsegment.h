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
#include "contactsegment.h"
#include "Elements/nlstructuralelement.h"
#include "set.h"
#include "feinterpol2d.h"

#define _IFT_BoundaryContactSegment_edgeSet "edgeset"
#define _IFT_BoundaryContactSegment_pairUpdateMode "pairupdatemode"

namespace oofem {
	class BoundaryContactSegment : public ContactSegment
	{
	public:
		BoundaryContactSegment(int n, Domain *aDomain) : ContactSegment(n, aDomain) { ; }
		~BoundaryContactSegment() {};

		virtual IRResultType initializeFrom(InputRecord * ir) override;
		//returns normalized n, which is an normal vector of contact 
		virtual void computeNormal(FloatArray& answer, Node * node, TimeStep* tstep) override = 0;
		//returns non-normalized tangent vector at the point of contact
		virtual void computeTangent(FloatArray &answer, Node *node, TimeStep *tstep) override = 0;
		//returns an extended N (aka A) matrix, integrated at point of contact of given node
		virtual void computeSegmentNMatrix(FloatMatrix& answer, Node* node, TimeStep * tStep) override = 0;
		//@todo: the following two functions should be merged
		virtual void computeSegmentBMatrix(FloatMatrix &answer, Node *node, TimeStep *tStep) override = 0;
		virtual void computedNdksi(FloatMatrix &answer, Node *node, TimeStep *tStep) = 0;
		//computes the penetration of the slave node to the given(closest) segment
		virtual double computePenetration(Node * node, TimeStep * tStep) override = 0;
		// check for large strain mode      
		virtual bool hasNonLinearGeometry(Node *node, TimeStep *tStep) override;
		// computes metric tensor, m_ij = t_i \cdot t_j, for 2d linear element, this is equivalent to l^2
		virtual void computeMetricTensor(FloatMatrix &answer, Node *node, TimeStep *tStep) override = 0;
		//return location array of the slave node and master segment
		virtual void giveLocationArray(const IntArray& dofIdArray, IntArray& s_loc, const UnknownNumberingScheme& c_s) const override = 0;
		//this function is needed to allocate appropriate components of stiffness matrix
		virtual void giveLocationArrays(const IntArray& dofIdArray, IntArray& s_loc, const UnknownNumberingScheme& c_s) override = 0;
		virtual void updateYourself(TimeStep * tStep) override;

		virtual void postInitialize() override;

	protected:
		IntArray edges;
		IntArray lastEdge; //last edge worked with
		std::vector<Node*> knownNodes;
		std::vector<IntArray> knownClosestEdges;
		int setnum;

		enum UpdateMode { UM_Never = 0, UM_EachStep = 1, UM_EveryQuery = 2 };
		UpdateMode updateMode;

        //searches the array of nodes to which a closest edge was already computed
        //returns -1 if unsuccessful
        inline int giveIndexOfKnownNode( const Node *node ){
            for ( int i = knownNodes.size() - 1; i >= 0; i-- ) {
				if ( node == knownNodes.at( i ) ) return i;
            }
			return -1;
		}

		//gives the closest edge to a given node in the form of an IntArray(elempos,edgepos)
		//only computes it again if it wasn't determined for this node in this solution step yet
		//(i.e. the array of known closest edges is reset after step convergence)
		void giveClosestEdge(IntArray& answer, Node * node, TimeStep * tStep);

		//computes the closest point on a given element to a given node. Answer is in parametric coordinates of element boundary
		//specific per type of elements used - needs to be overriden by child class
		virtual bool computeContactPoint(FloatArray &ksi, Node* node, StructuralElement* element, int elemedge, TimeStep* tStep) = 0;
	};
} //end namespace oofem

