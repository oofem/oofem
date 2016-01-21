/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
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

#include "trwarp.h"
#include "node.h"
#include "CrossSections/warpingcrosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "CrossSections/structuralcrosssection.h"
#include "classfactory.h"
#include "load.h"
#include "EngineeringModels/freewarping.h"
#include "engngm.h"
#include "dof.h"



namespace oofem {
REGISTER_Element(Tr_Warp);

FEI2dTrLin Tr_Warp :: interp(1, 2);

Tr_Warp :: Tr_Warp(int n, Domain *aDomain) :
    StructuralElement(n, aDomain), SpatialLocalizerInterface(this),  ZZNodalRecoveryModelInterface(this)
    // Constructor.
{
    numberOfDofMans  = 3;
    numberOfGaussPoints = 1;
}

Tr_Warp :: ~Tr_Warp()
// Destructor
{ }


void
Tr_Warp :: computeGaussPoints()
// Sets up the array containing the Gauss point of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 4) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}


IRResultType
Tr_Warp :: initializeFrom(InputRecord *ir)
{
    numberOfGaussPoints = 1;
    return StructuralElement :: initializeFrom(ir);
}


void
Tr_Warp :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveRealStress_Warping(answer, gp, strain, tStep);
}

void
Tr_Warp :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
{
    int numNodes = this->giveNumberOfDofManagers();
    FloatArray N(numNodes);

    //    int dim = this->giveSpatialDimension();

    answer.resize(1, numNodes);
    answer.zero();
    giveInterpolation()->evalN( N, iLocCoord, FEIElementGeometryWrapper(this) );

    answer.beNMatrixOf(N, 1);
}

void
Tr_Warp :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer,
                            int li, int ui)
// Returns the [2x4] strain-displacement matrix {B} of the receiver, eva-
// luated at gp.
{
    FloatMatrix dN;
    FloatArray tc(2);
    this->interp.evaldNdx( dN, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    // gausspoint coordinates
    FloatArray gcoords;
    Element *elem = gp->giveElement();
    elem->computeGlobalCoordinates( gcoords, gp->giveNaturalCoordinates() );

    this->transformCoordinates( tc, gcoords, this->giveCrossSection()->giveNumber() );

    answer.resize(2, 4);

    answer.at(1, 1) = dN.at(1, 1);
    answer.at(1, 2) = dN.at(2, 1);
    answer.at(1, 3) = dN.at(3, 1);
    answer.at(1, 4) = -tc.at(2);

    answer.at(2, 1) = dN.at(1, 2);
    answer.at(2, 2) = dN.at(2, 2);
    answer.at(2, 3) = dN.at(3, 2);
    answer.at(2, 4) = tc.at(1);
}


void
Tr_Warp :: transformCoordinates(FloatArray &answer, FloatArray &c, const int CGnumber) {
    answer.resize(2);
    FreeWarping *em = dynamic_cast< FreeWarping * >( this->giveDomain()->giveEngngModel() );
    if ( em ) {
        FloatMatrix CG;
        em->getCenterOfGravity(CG);
        answer.at(1) = c.at(1) - CG.at(CGnumber, 1);
        answer.at(2) = c.at(2) - CG.at(CGnumber, 2);
    } else   {
        OOFEM_ERROR("Error during transformCoordinates");
    }
}

void
Tr_Warp :: computeFirstMomentOfArea(FloatArray &answer)
// Returns the portion of the receiver which is attached to gp.
{
    answer.resize(2);

    FloatArray gcoords;
    GaussPoint *gp = this->giveDefaultIntegrationRulePtr()->getIntegrationPoint(0);
    double A = this->computeVolumeAround(gp);
    Element *elem = gp->giveElement();
    elem->computeGlobalCoordinates( answer, gp->giveNaturalCoordinates() );
    answer.times(A);
}

double
Tr_Warp :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double determinant, weight, volume;
    determinant = fabs( this->interp.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
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
Tr_Warp :: computeLocalForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
// computes the part of load vector, which is imposed by force loads acting

{
    this->computeEdgeLoadVectorAt(answer, NULL, tStep, mode);
}


void
Tr_Warp :: computeEdgeLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode)
{
    // computes the edge load vector of the receiver corresponding to the inhomogeneous Neumann condition
    // the given load is a dummy variable because the boundary condition for the warping equation
    // is determined by the geometry and thus the load intensity is not needed
    // (what is needed is just the indication that the given element edge is a part of the domain boundary)
    answer.resize(4);
    answer.zero();
    for ( int iEdge = 1; iEdge < 4; iEdge++ ) {
        IntArray mask;
        this->giveEdgeDofMapping(mask, iEdge);

        // coordinates of the initial and final node of the edge
        FloatArray *coord1 = giveNode( mask.at(1) )->giveCoordinates();
        FloatArray *coord2 = giveNode( mask.at(2) )->giveCoordinates();
        // components of the edge vector (from the initial to the final node)
        double dx = coord2->at(1) - coord1->at(1);
        double dy = coord2->at(2) - coord1->at(2);
        // coordinates of the initial node
        double x1 =  coord1->at(1);
        double y1 =  coord1->at(2);
        // coordinates of the final node
        double x2 = coord2->at(1);
        double y2 = coord2->at(2);

        // transform to coordinates w.r. center of gravity
        FloatArray tc1(2), c1(2), tc2(2), c2(2);
        c1.at(1) = x1;
        c1.at(2) = y1;
        this->transformCoordinates( tc1, c1, this->giveCrossSection()->giveNumber() );
        x1 = tc1.at(1);
        y1 = tc1.at(2);
        c2.at(1) = x2;
        c2.at(2) = y2;
        this->transformCoordinates( tc2, c2, this->giveCrossSection()->giveNumber() );
        x2 = tc2.at(1);
        y2 = tc2.at(2);

        // equivalent nodal "loads" (obtained by exact integration)
        double f1 = ( dx * x2 + dy * y2 ) / 3.0 + ( dx * x1 + dy * y1 ) / 6.0;
        double f2 = ( dx * x1 + dy * y1 ) / 3.0 + ( dx * x2 + dy * y2 ) / 6.0;

        // the load value has the meaning of relative twist
        double theta = this->giveDofManager(4)->giveDofWithID(24)->giveUnknown(VM_Total, tStep);
        FloatArray reducedAnswer, b;
        b.resize(4);
        b.zero();
        reducedAnswer.resize(2);
        reducedAnswer.at(1) =  f1 * theta;
        reducedAnswer.at(2) = f2 * theta;
        b.assemble(reducedAnswer, mask);
        answer.add(b);
    }
}


double
Tr_Warp :: giveThicknessAt(const FloatArray &gcoords)
{
    return 1.;
}


double
Tr_Warp :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double determinant = fabs( this->interp.edgeGiveTransformationJacobian( iEdge, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
    return determinant * gp->giveWeight();
}


void
Tr_Warp :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    // define new dof names (types) in dofiditem.h
    if ( inode > 0 && inode < 4 ) {
        answer = {
            Warp_PsiTheta
        };
    } else   {
        OOFEM_ERROR("Wrong numer of node");
    }
}


void
Tr_Warp :: giveInternalDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        Warp_Theta
    };
}

DofManager *
Tr_Warp ::  giveInternalDofManager(int i) const
{
    return this->giveDofManager(4);
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

int
Tr_Warp :: SpatialLocalizerI_containsPoint(const FloatArray &coords)
{
    FloatArray lcoords;
    return this->computeLocalCoordinates(lcoords, coords);
}


double
Tr_Warp :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
{
    FloatArray lcoords(3), gcoords;
    lcoords.at(1) = lcoords.at(2) = lcoords.at(3) = 1. / 3.;
    this->computeGlobalCoordinates(gcoords, lcoords);
    return gcoords.distance(coords);
}


void
Tr_Warp :: ZZNodalRecoveryMI_computeNNMatrix(FloatArray &answer, InternalStateType type)
{
    //
    // Returns NTN matrix (lumped) for Zienkiewicz-Zhu
    // The size of N mtrx is (nstresses, nnodes*nstreses)
    // Definition : sigmaVector = N * nodalSigmaVector
    //
    double volume = 0.0;
    FloatMatrix fullAnswer;
    FloatArray n;
    Element *elem  = this->ZZNodalRecoveryMI_giveElement();
    FEInterpolation *interpol = elem->giveInterpolation();
    IntegrationRule *iRule = elem->giveDefaultIntegrationRulePtr();

    if ( !interpol ) {
        OOFEM_ERROR( "ZZNodalRecoveryMI_computeNNMatrix: Element %d not providing interpolation", elem->giveNumber() );
    }

    int size = 3; //elem->giveNumberOfDofManagers();
    fullAnswer.resize(size, size);
    fullAnswer.zero();
    double pok = 0.0;

    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(i);
        double dV = elem->computeVolumeAround(gp);
        interpol->evalN( n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(elem) );
        fullAnswer.plusDyadSymmUpper(n, dV);
        pok += ( n.at(1) * dV ); ///@todo What is this? Completely unused.
        volume += dV;
    }


    fullAnswer.symmetrized();
    answer.resize(4);
    for ( int i = 1; i <= 3; i++ ) {
        double sum = 0.0;
        for ( int j = 1; j <= size; j++ ) {
            sum += fullAnswer.at(i, j);
        }

        answer.at(i) = sum;
    }
    answer.at(4) = 1.0;
}

bool
Tr_Warp :: ZZNodalRecoveryMI_computeNValProduct(FloatMatrix &answer, InternalStateType type,
                                                TimeStep *tStep)
{  // evaluates N^T sigma over element volume
   // N(nsigma, nsigma*nnodes)
   // Definition : sigmaVector = N * nodalSigmaVector
    FloatArray stressVector, n;
    Element *elem  = this->ZZNodalRecoveryMI_giveElement();
    FEInterpolation *interpol = elem->giveInterpolation();
    IntegrationRule *iRule = elem->giveDefaultIntegrationRulePtr();

    answer.clear();
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(i);
        double dV = elem->computeVolumeAround(gp);
        //this-> computeStressVector(stressVector, gp, tStep);
        if ( !elem->giveIPValue(stressVector, gp, type, tStep) ) {
            continue;
        }

        interpol->evalN( n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(elem) );
        answer.plusDyadUnsym(n, stressVector, dV);

        //  help.beTProductOf(n,stressVector);
        //  answer.add(help.times(dV));
    }
    answer.resizeWithData( 4, answer.giveNumberOfColumns() );
    for ( int i = 1; i <= answer.giveNumberOfColumns(); i++ ) {
        answer.at(4, i) = 0.0;
    }
    return true;
}

void
Tr_Warp :: postInitialize()
{
    Element :: postInitialize();
    dofManArray.resizeWithValues(4);
    WarpingCrossSection *wcs = dynamic_cast< WarpingCrossSection * >( this->giveCrossSection() );
    dofManArray.at(4) = wcs->giveWarpingNodeNumber();
}
} // end namespace oofem
