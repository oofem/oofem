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
#include "femcmpnn.h"
#include "node.h"
#include "inputrecord.h"

#define _IFT_ContactSegment_normMode "normmode"

namespace oofem {
    class ContactSegment : public FEMComponent
    {
    public:
		ContactSegment(int n, Domain *aDomain) : FEMComponent(n, aDomain){;}
		~ContactSegment() {};

		virtual IRResultType initializeFrom(InputRecord* ir) override {
			//@todo is this actually still necessary since the change of approach to normal computation?
			IRResultType result;

			normmode = NM_Never;
			int normmodeint = 0;
			IR_GIVE_OPTIONAL_FIELD(ir, normmodeint, _IFT_ContactSegment_normMode);
			if ( result == IRRT_OK ) {
              normmode = (NormalizationMode)normmodeint;
              if ( normmodeint < 0 || normmodeint > 2 ) OOFEM_ERROR("Contact segment normalization mode can be only 0, 1 or 2");
			}

			 return FEMComponent::initializeFrom(ir);
		};

		/**
		 * Computes the normal vector to the segment at the contact point of the given node.
		 * Contact point in this context is the point on the segment closest to the node,
		 * regardless of whether contact is actually occuring.
		 *
		 * @param answer [out] - the normal vector
		 * @param node - pointer to the node concerned
		 * @param tStep - the relevant time step
		 */
		virtual void computeNormal(FloatArray& answer, Node * node, TimeStep* tstep) = 0;
      
		/**
		 * Computes the tangential vector to the segment at the contact point of the given node.
		 * Contact point in this context is the point on the segment closest to the node,
		 * regardless of whether contact is actually occuring.
		 *
		 * @param answer [out] - the tangent vector
		 * @param node - pointer to the node concerned
		 * @param tStep - the relevant time step
		 */
		virtual void computeTangent( FloatArray &answer, Node *node, TimeStep *tstep ) = 0;

		/**
		 * Computes the N matrix (matrix of interpolation functions) of the segment at the contact point 
		 * of the given node. Contact point in this context is the point on the segment closest to the node,
		 * regardless of whether contact is actually occuring. 
		 *
		 * @param answer [out] - the N matrix
		 * @param node - pointer to the node concerned
		 * @param tStep - the relevant time step
		 */
		virtual void computeSegmentNMatrix(FloatMatrix& answer, Node* node, TimeStep * tStep) = 0;

		/**
		 * Computes the B matrix of the segment at the contact point
		 * of the given node. Contact point in this context is the point on the segment closest to the node,
		 * regardless of whether contact is actually occuring.
		 *
		 * For contact segments representing element boundaries, the B matrix is the derivative of the N matrix
		 * by the relevant parametric coordinate(s) of the element boundary
		 *
		 * @param answer [out] - the B matrix
		 * @param node - pointer to the node concerned
		 * @param tStep - the relevant time step
		 */
		virtual void computeSegmentBMatrix(FloatMatrix &answer, Node *node, TimeStep *tStep ) = 0;

		/**
		 * Gives the magnitude of the penetration of the given node into the contact segment, i.e. the Euclidian distace between
		 * node and the contact point on the segment.
		 *
		 * @returns the value of penetration, negative if penetrated, positive otherwise
		 *
		 * @param node - pointer to the node concerned
		 * @param tStep - the relevant time step
		 */
		virtual double computePenetration(Node * node, TimeStep * tStep) = 0;

		/**
		 * Answers whether the contact segment is comprised of elements with non-linear geometry,
		 * i.e. whether large strains shall be considered.
		 *
		 * @returns true or false
		 *
		 * @param node - pointer to the node concerned
		 * @param tStep - the relevant time step
		 */
		virtual bool hasNonLinearGeometry( Node *node, TimeStep *tStep ) = 0;

		/**
		 * Computes the metric tensor in contact of the given node with the segment.
		 *
		 * @param answer [out] the metric tensor.
		 * @param node - pointer to the node concerned
		 * @param tStep - the relevant time step
		 */
		virtual void computeMetricTensor( FloatMatrix &answer, Node *node, TimeStep *tStep ) = 0;

		/**
		 * Gives the location array of the segment relevant to the last contact considered.
		 * (differs from giveLocationArrays() for those segments which represent multiple elements.
		 *  In such case, this function gives only the location array for the relevant element and is to be
		 *  used for local-to-global matrix localization)
		 *
		 * @param dofIdArray the relevant DOFs
		 * @param s_loc [out] the answer
		 * @param tStep - the relevant time step
		 */
		virtual void giveLocationArray(const IntArray& dofIdArray, IntArray& s_loc, const UnknownNumberingScheme& c_s) const  = 0;

		/**
		 * Gives the location array of the segment relevant to all possible contact arrays.
		 * (differs from giveLocationArray() for those segments which represent multiple elements.
		 *  In such case, this function gives all the element locarrays and is to be used to build sparse matrices)
		 *
		 * @param dofIdArray the relevant DOFs
		 * @param s_loc [out] the answer
		 * @param tStep - the relevant time step
		 */
		virtual void giveLocationArrays(const IntArray& dofIdArray, IntArray& s_loc, const UnknownNumberingScheme& c_s) = 0;

		virtual void updateYourself(TimeStep * tStep) { ; };

		virtual void postInitialize(){;}

    protected:

        enum NormalizationMode { NM_Never = 0, NM_IfNotSmall = 1, NM_Always = 2 };
        NormalizationMode normmode;
    };
}

