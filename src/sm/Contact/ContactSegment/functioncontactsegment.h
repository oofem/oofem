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
#include "node.h"
#include "inputrecord.h"
#include "Elements/structuralelement.h"
#include "set.h"


namespace oofem {

    //virtual class parenting all analytical function contact segments
    //all children must implement computeContactPoint() - and need not implement anything else
    class FunctionContactSegment : public ContactSegment
    {
    public:
        FunctionContactSegment(int n, Domain *aDomain) : ContactSegment(n, aDomain) { ; }
        ~FunctionContactSegment() {};

        //returns normalized n, which is an normal vector of contact
        void computeNormal(FloatArray& answer, Node * node, TimeStep* tstep) override;

		void computeTangent( FloatArray &answer, Node *node, TimeStep *tstep ) override;

        //returns an extended N (aka A) matrix, integrated at point of contact of given node
        void computeSegmentNMatrix(FloatMatrix& answer, Node* node, TimeStep * tStep) override;

		//returns an extended N (aka A) matrix, integrated at point of contact of given node
        void computeSegmentBMatrix( FloatMatrix &answer, Node *node, TimeStep *tStep ) override;

		bool hasNonLinearGeometry( Node *node, TimeStep *tStep ) override;

		void computeMetricTensor( FloatMatrix &answer, Node *node, TimeStep *tStep ) override;

        //computes the penetration of node given 
        double computePenetration(Node * node, TimeStep * tStep) override;

        void giveLocationArray(const IntArray& dofIdArray, IntArray& s_loc, const UnknownNumberingScheme& c_s) const override;
        void giveLocationArrays(const IntArray& dofIdArray, IntArray& s_loc, const UnknownNumberingScheme& c_s) override;

    protected:

        virtual void computeContactPoint(FloatArray& answer, FloatArray& normal, const FloatArray& nodeCoords) = 0;
    };

}
