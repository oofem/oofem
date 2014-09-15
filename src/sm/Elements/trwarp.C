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

#include "../sm/Elements/trwarp.h"
#include "fei2dtrlin.h"
#include "node.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "classfactory.h"
#include "load.h"


namespace oofem {
REGISTER_Element(Tr_Warp);

FEI2dTrLin Tr_Warp :: interp(1, 2);

Tr_Warp :: Tr_Warp(int n, Domain *aDomain) :
    StructuralElement(n, aDomain), SpatialLocalizerInterface(this), ZZNodalRecoveryModelInterface(this)
    // Constructor.
{
    numberOfDofMans  = 3;
    numberOfGaussPoints = 1;
}

Tr_Warp :: ~Tr_Warp()
// Destructor
{ }


FEInterpolation *
Tr_Warp :: giveInterpolation() const { return & this->interp; }


void
Tr_Warp :: computeGaussPoints()
// Sets up the array containing the Gauss point of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}


IRResultType
Tr_Warp :: initializeFrom(InputRecord *ir)
{
    numberOfGaussPoints = 1;
    this->StructuralElement :: initializeFrom(ir);

    return IRRT_OK;
}


void
Tr_Warp :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer,
                                    int li, int ui)
// Returns the [3x6] strain-displacement matrix {B} of the receiver, eva-
// luated at gp.
{
    FloatMatrix dN;
    this->interp.evaldNdx( dN, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(2, 3);

    answer.at(1, 1) = dN.at(1, 1);
    answer.at(1, 2) = dN.at(2, 1);
    answer.at(1, 3) = dN.at(3, 1);

    answer.at(2, 1) = dN.at(1, 2);
    answer.at(2, 2) = dN.at(2, 2);
    answer.at(2, 3) = dN.at(3, 2);
}


double
Tr_Warp :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double determinant, weight, volume;
    determinant = fabs( this->interp.giveTransformationJacobian( * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
    weight = gp->giveWeight();
    volume = determinant * weight;

    return volume;
}

void
Tr_Warp :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */

    answer.resize(2);
    if ( iEdge == 1 ) { // edge between nodes 1,2
        answer.at(1) = 1;
        answer.at(2) = 2;
    } else if ( iEdge == 2 ) { // edge between nodes 2 3
        answer.at(1) = 2;
        answer.at(2) = 3;
    } else if ( iEdge == 3 ) { // edge between nodes 3 1
        answer.at(1) = 3;
        answer.at(2) = 1;
    } else {
        OOFEM_ERROR("wrong edge number");
    }
}


void
Tr_Warp :: computeEdgeLoadVectorAt(FloatArray &answer, Load *load,
                                             int iEdge, TimeStep *tStep, ValueModeType mode)
{
    // computes the edge load vector of the receiver corresponding to the inhomogeneous Neumann condition
    // the given load is a dummy variable because the boundary condition for the warping equation 
    // is determined by the geometry and thus the load intensity is not needed
    // (what is needed is just the indication that the given element edge is a part of the domain boundary)

    IntArray mask;
    this->giveEdgeDofMapping(mask, iEdge);

    // coordinates of the initial and final node of the edge
    FloatArray* coord1 = giveNode(mask.at(1)) -> giveCoordinates();
    FloatArray* coord2 = giveNode(mask.at(2)) -> giveCoordinates();
    // components of the edge vector (from the initial to the final node)
    double dx = coord2 -> at(1) - coord1 -> at(1); 
    double dy = coord2 -> at(2) - coord1 -> at(2); 
    // coordinates of the edge center
    double x = ( coord2 -> at(1) + coord1 -> at(1) ) / 2.; 
    double y = ( coord2 -> at(2) + coord1 -> at(2) ) / 2.; 
    
    // scalar product of the edge vector and the position vector of the edge center
    double f = dx * x + dy * y;

    // the load value has the meaning of relative twist
    FloatArray theta, dummy;
    load -> computeValueAt(theta, tStep, dummy, VM_Total);

    // equivalent nodal "loads" 
    FloatArray reducedAnswer;
    reducedAnswer.resize(2);
    reducedAnswer.at(1) = reducedAnswer.at(2) = theta.at(1) * f / 2.;
    answer.resize( this->computeNumberOfDofs() );
    answer.zero();
    answer.assemble(reducedAnswer, mask);
}


double
Tr_Warp :: giveThicknessAt(const FloatArray &gcoords)
{
    return 1.;
}


double
Tr_Warp :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double determinant = fabs( this->interp.edgeGiveTransformationJacobian( iEdge, * gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
    FloatArray gc;
    return determinant * gp->giveWeight();
}


void
Tr_Warp :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_w};
}

Interface *
Tr_Warp :: giveInterface(InterfaceType interface)
{
    if ( interface == SpatialLocalizerInterfaceType ) {
        return static_cast< SpatialLocalizerInterface * >(this);
	//    } else if ( interface == EIPrimaryFieldInterfaceType ) {
	//        return static_cast< EIPrimaryFieldInterface * >(this);
    } else if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    }

    return NULL;
}

double
Tr_Warp :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
{
    FloatArray lcoords(3), gcoords;
    lcoords.at(1) = lcoords.at(2) = lcoords.at(3) = 1. / 3.;
    this->computeGlobalCoordinates(gcoords, lcoords);
    return gcoords.distance(coords);
}
} // end namespace oofem
