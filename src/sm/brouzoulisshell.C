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

#include "brouzoulisshell.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "structuralcrosssection.h"
#include "classfactory.h"
#include "node.h"
#include "dofmanager.h"
#include "layeredcrosssection.h"

namespace oofem {
REGISTER_Element(BrouzoulisShell);

FEI2dQuadLin BrouzoulisShell :: interpolation(1,1);

BrouzoulisShell :: BrouzoulisShell(int n, Domain *aDomain) : NLStructuralElement(n, aDomain)
{
    numberOfDofMans = 4;
    corotTriad.resize(3,3);
    corotTriad.beUnitMatrix();
    tempCorotTriad.resize(3,3);
    tempCorotTriad.beUnitMatrix();
}



IRResultType
BrouzoulisShell :: initializeFrom(InputRecord *ir)
{
    numberOfGaussPoints = -1;

    return this->NLStructuralElement :: initializeFrom(ir);
}




void
BrouzoulisShell :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    answer.setValues(6, D_u, D_v, D_w, R_u, R_v, R_w);
}



void
BrouzoulisShell :: computeStressVector(FloatArray &answer, const FloatArray &e, GaussPoint *gp, TimeStep *tStep)
{
    //this->giveStructuralCrossSection()->giveRealStress_Plate(answer, gp, e, tStep);
    
    //giveRealStressVector_PlateLayer(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedStrain, TimeStep *tStep)
}


void
BrouzoulisShell :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->give2dPlateStiffMtrx(answer, rMode, gp, tStep);

}


double
BrouzoulisShell :: computeVolumeAround(GaussPoint *gp)
{
    double determinant = fabs( this->interpolation.giveTransformationJacobian( * gp->giveCoordinates(), FEIElementGeometryWrapper(this) ) );
    double weight      = gp->giveWeight();
    return ( determinant * weight );
}


void
BrouzoulisShell :: computeGaussPoints()
// Sets up the array containing the Gauss points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ numberOfIntegrationRules ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 5);
        numberOfGaussPoints = 1;
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}


void
BrouzoulisShell :: computedNdX(FloatMatrix &answer)
{
    FloatMatrix coords;
    this->computeCoRotatedNodeCoords(coords);
    // A = 0.5 * (x31*y42 +x24*y31) = 4J
    double A = 0.5 * ( ( coords.at(1,3)-coords.at(1,1) )*( coords.at(2,4)-coords.at(2,2) ) + 
                       ( coords.at(1,2)-coords.at(1,4) )*( coords.at(2,3)-coords.at(2,1) ) ); 

    answer.resize(2,4);
    int x = 1; int y = 2;
    answer.at(1,1) = ( coords.at(y,2)-coords.at(y,4) ); 
    answer.at(1,2) = ( coords.at(y,3)-coords.at(y,1) );
    answer.at(1,3) = ( coords.at(y,4)-coords.at(y,2) );
    answer.at(1,4) = ( coords.at(y,1)-coords.at(y,3) );

    answer.at(2,1) = ( coords.at(x,4)-coords.at(x,2) ); 
    answer.at(2,2) = ( coords.at(x,1)-coords.at(x,3) );
    answer.at(2,3) = ( coords.at(x,2)-coords.at(x,4) );
    answer.at(2,4) = ( coords.at(x,3)-coords.at(x,1) );
    
    double fac1 = 1.0 / (2.0*A);
    answer.times(fac1);
}

void
BrouzoulisShell :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [6x60] strain-displacement matrix {B} of the receiver, eva-
// luated at gp.
// B matrix  -  3 rows : epsilon-X, epsilon-Y, gamma-YZ, gamma-ZX, gamma-XY  :
{

    FloatMatrix coords;
    this->computeCoRotatedNodeCoords(coords);
    // A = 0.5 * (x31*y42 +x24*y31) = 4J
    double A = 0.5 * ( ( coords.at(1,3)-coords.at(1,1) )*( coords.at(2,4)-coords.at(2,2) ) + 
                       ( coords.at(1,2)-coords.at(1,4) )*( coords.at(2,3)-coords.at(2,1) ) ); 

    FloatMatrix dNdx(2,4);
    int x = 1; int y = 2;
    dNdx.at(1,1) = ( coords.at(y,2)-coords.at(y,4) ); 
    dNdx.at(1,2) = ( coords.at(y,3)-coords.at(y,1) );
    dNdx.at(1,3) = ( coords.at(y,4)-coords.at(y,2) );
    dNdx.at(1,4) = ( coords.at(y,1)-coords.at(y,3) );

    dNdx.at(2,1) = ( coords.at(x,4)-coords.at(x,2) ); 
    dNdx.at(2,2) = ( coords.at(x,1)-coords.at(x,3) );
    dNdx.at(2,3) = ( coords.at(x,2)-coords.at(x,4) );
    dNdx.at(2,4) = ( coords.at(x,3)-coords.at(x,1) );
    
    //             
    double zeta = 0.0;
    double fac1 = 1.0 / (2.0*A);
    dNdx.times(fac1);
    answer.resize(5,20);

    double N = 0.25; //N(0,0) = 1/4
    for ( int i = 1, j = 0; i <= dNdx.giveNumberOfColumns(); i++, j += 5 ) {

        answer.at(1, j + 1) = dNdx.at(x, i);         // du/dx
        answer.at(1, j + 5) = zeta * dNdx.at(x, i);  // du/dx

        answer.at(2, j + 2) =         dNdx.at(y, i); // dv/dy
        answer.at(2, j + 4) = -zeta * dNdx.at(y, i); // dv/dy

        answer.at(3, j + 1) = dNdx.at(y, i);         // du/dy
        answer.at(3, j + 2) = dNdx.at(x, i);         // dv/dx
        answer.at(3, j + 4) = -zeta * dNdx.at(x, i);         
        answer.at(3, j + 5) =  zeta * dNdx.at(x, i);         

        answer.at(4, j + 3) = dNdx.at(x, i);         
        answer.at(4, j + 5) = N;

        answer.at(5, j + 3) = dNdx.at(y, i);         
        answer.at(5, j + 4) = -N;

        
    }
    
}


bool
BrouzoulisShell :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    answer.resize(24, 20);
    answer.zero();
    FloatMatrix temp;
    temp.beSubMatrixOf(this->tempCorotTriad,1,2,1,2);

    for ( int i = 0; i < 4; i++ ) { // Loops over nodes
        // In each node, transform global c.s. {D_u, D_v, D_w, R_u, R_v, R_w} into local c.s.
        answer.setSubMatrix(this->tempCorotTriad, 1 + i * 6, 1 + i * 5);     // Displacements
        answer.setSubMatrix(temp, 1 + i * 6 + 3, 1 + i * 5 + 3); // Rotations
    }
    temp = answer;
    answer.beTranspositionOf(temp);
    return true;
}

void
BrouzoulisShell :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    //// BH = [BH1 BH2]
    //// BH1 is for the lower order approximation
    //FloatMatrix dNdx1x2, dNdx3;

    //this->interpolation.evaldNdx( dNdx1x2, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) ); // Interpolation for x1 and x2
    //this->interpolation.evaldNdx( dNdx3, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );   // Interpolation for x3

    //answer.resize(9, 27); // 6*2 + 15*1 = 27
    //answer.zero();
    //// everything associated with 'w' should use the higher interpolation 
    //for ( int i = 1; i <= dNdx1x2.giveNumberOfRows(); i++ ) {
    //    answer.at(1, 3 * i - 2) = dNdx1x2.at(i, 1);     // du/dx
    //    answer.at(2, 3 * i - 1) = dNdx1x2.at(i, 2);     // dv/dy
    //    answer.at(3, 3 * i - 0) = dNdx3.  at(i, 3);     // dw/dz
    //    answer.at(4, 3 * i - 1) = dNdx1x2.at(i, 3);     // dv/dz
    //    answer.at(7, 3 * i - 0) = dNdx3  .at(i, 2);     // dw/dy
    //    answer.at(5, 3 * i - 2) = dNdx1x2.at(i, 3);     // du/dz
    //    answer.at(8, 3 * i - 0) = dNdx3  .at(i, 1);     // dw/dx
    //    answer.at(6, 3 * i - 2) = dNdx1x2.at(i, 2);     // du/dy
    //    answer.at(9, 3 * i - 1) = dNdx1x2.at(i, 1);     // dv/dx
    //}

    //// Add contributions from higher order interpolation
    //int pos = 18;
    //int nodePos = 6;
    //for ( int i = 1; i <= 15-6; i++ ) {
    //    answer.at(3, pos + i) = dNdx3.at(nodePos + i, 3);     // dw/dz
    //    answer.at(7, pos + i) = dNdx3.at(nodePos + i, 2);     // dw/dy
    //    answer.at(8, pos + i) = dNdx3.at(nodePos + i, 1);     // dw/dx

    //}

}


void
BrouzoulisShell :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at gp.
{
    FloatArray N;

    this->interpolation.evalN( N, iLocCoord, FEIElementGeometryWrapper(this) ); 

    answer.resize(2, 27);
    answer.zero();

    for ( int i = 1; i <= 6; i++ ) {
        answer.at(1, 3 * i - 2) = N.at(i);
        answer.at(2, 3 * i - 1) = N.at(i);
        answer.at(3, 3 * i - 0) = N.at(i);
    }
}


void
BrouzoulisShell :: computeCoRotatedNodeCoords(FloatMatrix &answer)
{
    FloatArray coordsNode, coords;
    FloatMatrix nodeCoords(3,4);

    for ( int i = 1; i <= 4; i++ ) {
        coords = *this->giveNode(i)->giveCoordinates();
        nodeCoords.setColumn(coords,i);
    }
    //nodeCoords.rotatedWith(this->tempCorotTriad,'t');
    answer.beTProductOf(this->tempCorotTriad, nodeCoords);
}

// ******************************
// ***  Surface load support  ***
// ******************************
#if 1
IntegrationRule *
BrouzoulisShell :: GetSurfaceIntegrationRule(int approxOrder)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, this, 1, 1);
    int npoints = iRule->getRequiredNumberOfIntegrationPoints(_Triangle, approxOrder);
    iRule->SetUpPointsOnTriangle(npoints, _Unknown);
    return iRule;
}

void
BrouzoulisShell :: computeSurfaceNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *sgp)
{
    FloatArray n;
    //interpolationUV.surfaceEvalN( n, iSurf, * sgp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    //answer.beNMatrixOf(n, 3);
}

void
BrouzoulisShell :: giveSurfaceDofMapping(IntArray &answer, int iSurf) const
{
    IntArray nodes;
    const int ndofsn = 3;

    //interpolation.computeLocalSurfaceMapping(nodes, iSurf); // Use lower order app for loads, ok?

    answer.resize(9);

    for ( int i = 1; i <= 3; i++ ) {
        answer.at(i * ndofsn - 2) = nodes.at(i) * ndofsn - 2;
        answer.at(i * ndofsn - 1) = nodes.at(i) * ndofsn - 1;
        answer.at(i * ndofsn) = nodes.at(i) * ndofsn;
    }
}



double
BrouzoulisShell :: computeSurfaceVolumeAround(GaussPoint *gp, int iSurf)
{
    double determinant, weight, volume;
    //determinant = fabs( interpolation.surfaceGiveTransformationJacobian( iSurf, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) ) );

    weight = gp->giveWeight();
    volume = determinant * weight;

    return volume;
}




void
BrouzoulisShell :: computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iSurf)
{
    //interpolation.surfaceLocal2global( answer, iSurf, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
}

int
BrouzoulisShell :: computeLoadLSToLRotationMatrix(FloatMatrix &answer, int iSurf, GaussPoint *gp)
{
    // returns transformation matrix from
    // surface local coordinate system
    // to element local c.s
    // (same as global c.s in this case)
    //
    // i.e. f(element local) = T * f(edge local)

    // definition of local c.s on surface:
    // local z axis - perpendicular to surface, pointing outwards from element
    // local x axis - is in global xy plane (perpendicular to global z axis)
    // local y axis - completes the righ hand side cs.

    /*
     * _error ("computeLoadLSToLRotationMatrix: surface local coordinate system not supported");
     * return 1;
     */
    FloatArray gc(3);
    FloatArray h1(3), h2(3), nn(3), n(3);
    IntArray snodes(4);

    answer.resize(3, 3);
    answer.zero();

    //this->interpolation.computeSurfaceMapping(snodes, dofManArray, iSurf);
    for ( int i = 1; i <= 3; i++ ) {
        gc.add( * domain->giveNode( snodes.at(i) )->giveCoordinates() );
    }

    gc.times(1. / 0.5);
    // determine "average normal"
    for ( int i = 1; i <= 4; i++ ) {
        int j = ( i ) % 4 + 1;
        h1.beDifferenceOf(* domain->giveNode( snodes.at(i) )->giveCoordinates(), gc);
        h2.beDifferenceOf(* domain->giveNode( snodes.at(j) )->giveCoordinates(), gc);
        n.beVectorProductOf(h1, h2);
        if ( n.computeSquaredNorm() > 1.e-6 ) {
            n.normalize();
        }

        nn.add(n);
    }

    nn.times(1. / 4.);
    if ( nn.computeSquaredNorm() < 1.e-6 ) {
        answer.zero();
    }

    nn.normalize();
    for ( int i = 1; i <= 3; i++ ) {
        answer.at(i, 3) = nn.at(i);
    }

    // determine lcs of surface
    // local x axis in xy plane
    double test = fabs(fabs( nn.at(3) ) - 1.0);
    if ( test < 1.e-5 ) {
        h1.at(1) = answer.at(1, 1) = 1.0;
        h1.at(2) = answer.at(2, 1) = 0.0;
    } else {
        h1.at(1) = answer.at(1, 1) = answer.at(2, 3);
        h1.at(2) = answer.at(2, 1) = -answer.at(1, 3);
    }

    h1.at(3) = answer.at(3, 1) = 0.0;
    // local y axis perpendicular to local x,z axes
    h2.beVectorProductOf(nn, h1);
    for ( int i = 1; i <= 3; i++ ) {
        answer.at(i, 2) = h2.at(i);
    }

    return 1;
}

Interface *
BrouzoulisShell :: giveInterface(InterfaceType interface)
{
    if( interface == NodalAveragingRecoveryModelInterfaceType ) {
//        return static_cast< NodalAveragingRecoveryModelInterface * >( this );
    } else {
        OOFEM_LOG_INFO("Interface on BrouzoulisShell element not supported");
        return NULL;
    }
}



#endif



void
BrouzoulisShell :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    // Computes internal forces as a summation of: sectional forces

    NLStructuralElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord);
    
}

void 
BrouzoulisShell :: computeSectionalForcesAt(FloatArray &sectionalForces, FloatArray &stress, GaussPoint *gp, TimeStep *tStep){

    FloatArray vStress;

    
    //FloatArray PG1(3), PG2(3), PG3(3);
    //FloatArray cartStressVector, contravarStressVector;
    //FloatMatrix lambda[3];

    //FloatMatrix Gcov, Gcon; 
    //FloatArray lCoords = *gp->giveCoordinates();
    ////this->evalInitialContravarBaseVectorsAt(lCoords, Gcon);
    //FloatMatrix P, Sig;
    //P.beMatrixFormOfStress(vStress);
    //Sig.beProductTOf(Gcon,P);       // Stress resultants stored in each column
    //PG1.beColumnOf(Sig,1);
    //PG2.beColumnOf(Sig,2);
    //PG3.beColumnOf(Sig,3);
    //        
    //double z = gp->giveCoordinate(3) * this->giveStructuralCrossSection()->give(CS_Thickness, gp);
    //this->computeLambdaMatrices(lambda, z); // associated with the variation of the test functions   
    //
    //// f = lambda_1^T * PG1 + lambda_2^T * PG2 + lambda_3^T * PG3
    //sectionalForces.resize(18);
    //FloatArray temp;
    //sectionalForces.zero();
    //temp.beTProductOf(lambda [0], PG1);
    //sectionalForces.add(temp);
    //temp.beTProductOf(lambda [1], PG2);
    //sectionalForces.add(temp);
    //temp.beTProductOf(lambda [2], PG3);
    //sectionalForces.add(temp);
}





void
BrouzoulisShell :: setupInitialNodeDirectors()
{
    FloatArray lcoords(2), nodeLocalXiCoords, nodeLocalEtaCoords;

    // Give the local coordinates for the element nodes - all at once

    this->initialNodeDirectors.resize(6);
    this->initialNodeMidSurface.resize(6);
    FloatMatrix N;
    FloatArray Xbar, M, Xtop, Xbottom;
    for ( int node = 1; node <= 3; node++ ) {
        this->initialNodeDirectors [ node - 1 ].resize(3);
        this->initialNodeDirectors [ node - 1 ].zero();
        this->initialNodeMidSurface [ node - 1 ].resize(3);
        this->initialNodeMidSurface [ node - 1 ].zero();

        Xbar = 0.5 * ( *this->giveNode(node)->giveCoordinates() + *this->giveNode(3+node)->giveCoordinates() );
        M =          ( *this->giveNode(node)->giveCoordinates() - *this->giveNode(3+node)->giveCoordinates() );

        this->initialNodeMidSurface [ node - 1 ] = Xbar;
        this->initialNodeDirectors [ node - 1 ] = M;
    }

    for ( int node = 1; node <= 3; node++ ) {
        this->initialNodeDirectors [ node + 3 - 1 ].resize(3);
        this->initialNodeDirectors [ node + 3 - 1 ].zero();
        this->initialNodeMidSurface [ node + 3 - 1 ].resize(3);
        this->initialNodeMidSurface [ node + 3 - 1 ].zero();

        Xbar = 0.5 * ( *this->giveNode(6 + node)->giveCoordinates() + *this->giveNode(9 + node)->giveCoordinates() );
        M =          ( *this->giveNode(6 + node)->giveCoordinates() - *this->giveNode(9 + node)->giveCoordinates() );

        this->initialNodeMidSurface [ node + 3 - 1 ] = Xbar;
        this->initialNodeDirectors [ node + 3 - 1 ] = M;
    }

}




} // end namespace oofem





// OLD


#if 0
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

#include "brouzoulisshell.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "structuralcrosssection.h"
#include "classfactory.h"
#include "node.h"
#include "dofmanager.h"

namespace oofem {
REGISTER_Element(BrouzoulisShell);

FEI3dWedgeLin BrouzoulisShell :: interpolationUV;
FEI3dWedgeQuad BrouzoulisShell :: interpolationW;

BrouzoulisShell :: BrouzoulisShell(int n, Domain *aDomain) : NLStructuralElement(n, aDomain)
{
    numberOfDofMans = 15;
    //numberOfDofMans = 12;
    
    //this->createInternalNodes();
    

}

void
BrouzoulisShell :: createInternalNodes()
{

    // Create the 3 mid nodes for the quadratic interpolation.
    // They should be placed in between the top and bottom ones

    Domain *domain = this->giveDomain();
    
    FloatArray coords, XTop, XBottom;
    for ( int i = 1; i <= 3; i++ ) {
        int number = domain->giveNumberOfDofManagers();
        DofManager *dMan = classFactory.createDofManager("Node", number + 1, domain );
        dMan->setGlobalNumber(number + 1);

 /*       computeDofManDofIdArray(dofIdArray, dMan);
        int nDofs = dMan->giveNumberOfDofs();
        
        if ( !dMan->hasDofID( ( DofIDItem ) ( dofIdArray.at(m) ) ) ) {
                    dMan->appendDof( new MasterDof( nDofs + m, dMan, ( DofIDItem ) ( dofIdArray.at(m) ) ) );
        }
        }*/


        XTop    = *this->giveNode(i+3)->giveCoordinates();
        XBottom = *this->giveNode(i  )->giveCoordinates();
        coords = 0.5 * ( XTop + XBottom ); 
        static_cast< Node *> (dMan)->setCoordinates(coords);    

        // Add the newly created node
        
        this->giveDomain()->setDofManager( number+1, dMan);
        this->addDofManager(dMan);
        this->dofManArray.printYourself();
    }

    numberOfDofMans = 15;
}

IRResultType
BrouzoulisShell :: initializeFrom(InputRecord *ir)
{
    numberOfGaussPoints = -1;

    return this->NLStructuralElement :: initializeFrom(ir);
}




void
BrouzoulisShell :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    if ( inode <= 6 ) { 
        answer.setValues(3, D_u, D_v, D_w);
    } else {
        answer.setValues(1, D_w);   // The last nodes only have extra u_3 dofs 
    }
}



void
BrouzoulisShell :: computeStressVector(FloatArray &answer, const FloatArray &e, GaussPoint *gp, TimeStep *tStep)
{
        this->giveStructuralCrossSection()->giveRealStress_3d(answer, gp, e, tStep);
}


void
BrouzoulisShell :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveStiffnessMatrix_3d(answer, rMode, gp, tStep);
}


double
BrouzoulisShell :: computeVolumeAround(GaussPoint *gp)
{
    double determinant = fabs( this->interpolationW.giveTransformationJacobian( * gp->giveCoordinates(), FEIElementGeometryWrapper(this) ) );
    double weight      = gp->giveWeight();

    return ( determinant * weight );
}


void
BrouzoulisShell :: computeGaussPoints()
// Sets up the array containing the Gauss points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ numberOfIntegrationRules ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 6);
        
        integrationRulesArray [ 0 ]->SetUpPointsOnWedge(6, 2, _3dMat);
        //this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}


void
BrouzoulisShell :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [6x60] strain-displacement matrix {B} of the receiver, eva-
// luated at gp.
// B matrix  -  6 rows : epsilon-X, epsilon-Y, epsilon-Z, gamma-YZ, gamma-ZX, gamma-XY  :
{
    FloatMatrix dNdx1x2, dNdx3;

    this->interpolationUV.evaldNdx( dNdx1x2, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) ); // Interpolation for x1 and x2
    this->interpolationW.evaldNdx( dNdx3, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );   // Interpolation for x3

    answer.resize(6, 27); // 6*2 + 15*1 = 27
    answer.zero();
    for ( int i = 1; i <= dNdx1x2.giveNumberOfRows(); i++ ) {

        answer.at(1, 3 * i - 2) = dNdx1x2.at(i, 1); // du/dx
        answer.at(2, 3 * i - 1) = dNdx1x2.at(i, 2); // dv/dy
        answer.at(3, 3 * i - 0) = dNdx3.at(i, 3);   // dw/dz

        answer.at(4, 3 * i - 1) = dNdx1x2.at(i, 3); // dv/dz
        answer.at(4, 3 * i - 0) = dNdx3.at(i, 2);   // dw/dy

        answer.at(5, 3 * i - 2) = dNdx1x2.at(i, 3); // du/dz
        answer.at(5, 3 * i - 0) = dNdx3.at(i, 1);   // dw/dx

        answer.at(6, 3 * i - 2) = dNdx1x2.at(i, 2); // du/dy
        answer.at(6, 3 * i - 1) = dNdx1x2.at(i, 1); // dv/dx
    }
    // Add contributions from higher order interpolation
    int pos = 18;
    int nodePos = 6;
    for ( int i = 1; i <= 15-6; i++ ) {
        answer.at(3, pos + i) = dNdx3.at(nodePos + i, 3);     // dw/dz
        answer.at(4, pos + i) = dNdx3.at(nodePos + i, 2);     // dw/dy
        answer.at(5, pos + i) = dNdx3.at(nodePos + i, 1);     // dw/dx

    }

}


void
BrouzoulisShell :: computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    // BH = [BH1 BH2]
    // BH1 is for the lower order approximation
    FloatMatrix dNdx1x2, dNdx3;

    this->interpolationUV.evaldNdx( dNdx1x2, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) ); // Interpolation for x1 and x2
    this->interpolationW.evaldNdx( dNdx3, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );   // Interpolation for x3

    answer.resize(9, 27); // 6*2 + 15*1 = 27
    answer.zero();
    // everything associated with 'w' should use the higher interpolation 
    for ( int i = 1; i <= dNdx1x2.giveNumberOfRows(); i++ ) {
        answer.at(1, 3 * i - 2) = dNdx1x2.at(i, 1);     // du/dx
        answer.at(2, 3 * i - 1) = dNdx1x2.at(i, 2);     // dv/dy
        answer.at(3, 3 * i - 0) = dNdx3.  at(i, 3);     // dw/dz
        answer.at(4, 3 * i - 1) = dNdx1x2.at(i, 3);     // dv/dz
        answer.at(7, 3 * i - 0) = dNdx3  .at(i, 2);     // dw/dy
        answer.at(5, 3 * i - 2) = dNdx1x2.at(i, 3);     // du/dz
        answer.at(8, 3 * i - 0) = dNdx3  .at(i, 1);     // dw/dx
        answer.at(6, 3 * i - 2) = dNdx1x2.at(i, 2);     // du/dy
        answer.at(9, 3 * i - 1) = dNdx1x2.at(i, 1);     // dv/dx
    }

    // Add contributions from higher order interpolation
    int pos = 18;
    int nodePos = 6;
    for ( int i = 1; i <= 15-6; i++ ) {
        answer.at(3, pos + i) = dNdx3.at(nodePos + i, 3);     // dw/dz
        answer.at(7, pos + i) = dNdx3.at(nodePos + i, 2);     // dw/dy
        answer.at(8, pos + i) = dNdx3.at(nodePos + i, 1);     // dw/dx

    }

}


void
BrouzoulisShell :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at gp.
{
    FloatArray Nx1x2, Nx3;

    this->interpolationUV.evalN( Nx1x2, iLocCoord, FEIElementGeometryWrapper(this) ); // Interpolation for x1 and x2
    this->interpolationW.evalN( Nx3, iLocCoord, FEIElementGeometryWrapper(this) );    // Interpolation for x3

    answer.resize(3, 27);
    answer.zero();

    for ( int i = 1; i <= 6; i++ ) {
        answer.at(1, 3 * i - 2) = Nx1x2.at(i);
        answer.at(2, 3 * i - 1) = Nx1x2.at(i);
        answer.at(3, 3 * i - 0) = Nx3.at(i);
    }
    for ( int i = 1; i <= 15-6; i++ ) {
        answer.at(3, 18 + i) = Nx3.at(6 + i);
    }
}




// ******************************
// ***  Surface load support  ***
// ******************************
#if 1
IntegrationRule *
BrouzoulisShell :: GetSurfaceIntegrationRule(int approxOrder)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, this, 1, 1);
    int npoints = iRule->getRequiredNumberOfIntegrationPoints(_Triangle, approxOrder);
    iRule->SetUpPointsOnTriangle(npoints, _Unknown);
    return iRule;
}

void
BrouzoulisShell :: computeSurfaceNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *sgp)
{
    FloatArray n;
    //interpolationUV.surfaceEvalN( n, iSurf, * sgp->giveCoordinates(), FEIElementGeometryWrapper(this) );
    //answer.beNMatrixOf(n, 3);
}

void
BrouzoulisShell :: giveSurfaceDofMapping(IntArray &answer, int iSurf) const
{
    IntArray nodes;
    const int ndofsn = 3;

    interpolationUV.computeLocalSurfaceMapping(nodes, iSurf); // Use lower order app for loads, ok?

    answer.resize(9);

    for ( int i = 1; i <= 3; i++ ) {
        answer.at(i * ndofsn - 2) = nodes.at(i) * ndofsn - 2;
        answer.at(i * ndofsn - 1) = nodes.at(i) * ndofsn - 1;
        answer.at(i * ndofsn) = nodes.at(i) * ndofsn;
    }
}



double
BrouzoulisShell :: computeSurfaceVolumeAround(GaussPoint *gp, int iSurf)
{
    double determinant, weight, volume;
    determinant = fabs( interpolationUV.surfaceGiveTransformationJacobian( iSurf, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) ) );

    weight = gp->giveWeight();
    volume = determinant * weight;

    return volume;
}




void
BrouzoulisShell :: computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iSurf)
{
    interpolationUV.surfaceLocal2global( answer, iSurf, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
}

int
BrouzoulisShell :: computeLoadLSToLRotationMatrix(FloatMatrix &answer, int iSurf, GaussPoint *gp)
{
    // returns transformation matrix from
    // surface local coordinate system
    // to element local c.s
    // (same as global c.s in this case)
    //
    // i.e. f(element local) = T * f(edge local)

    // definition of local c.s on surface:
    // local z axis - perpendicular to surface, pointing outwards from element
    // local x axis - is in global xy plane (perpendicular to global z axis)
    // local y axis - completes the righ hand side cs.

    /*
     * _error ("computeLoadLSToLRotationMatrix: surface local coordinate system not supported");
     * return 1;
     */
    FloatArray gc(3);
    FloatArray h1(3), h2(3), nn(3), n(3);
    IntArray snodes(4);

    answer.resize(3, 3);
    answer.zero();

    this->interpolationUV.computeSurfaceMapping(snodes, dofManArray, iSurf);
    for ( int i = 1; i <= 3; i++ ) {
        gc.add( * domain->giveNode( snodes.at(i) )->giveCoordinates() );
    }

    gc.times(1. / 0.5);
    // determine "average normal"
    for ( int i = 1; i <= 4; i++ ) {
        int j = ( i ) % 4 + 1;
        h1.beDifferenceOf(* domain->giveNode( snodes.at(i) )->giveCoordinates(), gc);
        h2.beDifferenceOf(* domain->giveNode( snodes.at(j) )->giveCoordinates(), gc);
        n.beVectorProductOf(h1, h2);
        if ( n.computeSquaredNorm() > 1.e-6 ) {
            n.normalize();
        }

        nn.add(n);
    }

    nn.times(1. / 4.);
    if ( nn.computeSquaredNorm() < 1.e-6 ) {
        answer.zero();
    }

    nn.normalize();
    for ( int i = 1; i <= 3; i++ ) {
        answer.at(i, 3) = nn.at(i);
    }

    // determine lcs of surface
    // local x axis in xy plane
    double test = fabs(fabs( nn.at(3) ) - 1.0);
    if ( test < 1.e-5 ) {
        h1.at(1) = answer.at(1, 1) = 1.0;
        h1.at(2) = answer.at(2, 1) = 0.0;
    } else {
        h1.at(1) = answer.at(1, 1) = answer.at(2, 3);
        h1.at(2) = answer.at(2, 1) = -answer.at(1, 3);
    }

    h1.at(3) = answer.at(3, 1) = 0.0;
    // local y axis perpendicular to local x,z axes
    h2.beVectorProductOf(nn, h1);
    for ( int i = 1; i <= 3; i++ ) {
        answer.at(i, 2) = h2.at(i);
    }

    return 1;
}

Interface *
BrouzoulisShell :: giveInterface(InterfaceType interface)
{
    if( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >( this );
    } else {
        OOFEM_LOG_INFO("Interface on BrouzoulisShell element not supported");
        return NULL;
    }
}



#endif

void BrouzoulisShell :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side, InternalStateType type, TimeStep *tStep)
{
    answer.resize(0);
}

void BrouzoulisShell :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep)
{
    // N*a
    FloatMatrix N;
    FloatArray lCoords, a, xi1, xi2, xi3;

    // Local coords for a quadratic wedge element (VTK cell type 26)
    double z = 1.0;
    xi1.setValues(15, 1., 0., 0., 1., 0., 0., .5, 0., .5, .5, 0., .5, 1., 0., 0.);      
    xi2.setValues(15, 0., 1., 0., 0., 1., 0., .5, .5, 0., .5, .5, 0., 0., 1., 0.);
    xi3.setValues(15, -z, -z, -z,  z,  z,  z, -z, -z, -z,  z,  z,  z, 0., 0., 0.);
    lCoords.setValues(3,  xi1.at(node), xi2.at(node), xi3.at(node) );
    this->computeNmatrixAt(lCoords, N);
    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, a);
    answer.beProductOf(N,a);





}





void
BrouzoulisShell :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    // Computes internal forces as a summation of: sectional forces

    NLStructuralElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord);
    return;

    int ndofs = this->giveNumberOfDofs();

    int numberOfLayers = 1; 
    FloatArray a, f(ndofs), ftemp, sectionalForces;
    FloatArray genEps, genEpsD, totalSolVec, lCoords;
    FloatMatrix B;
    //this->giveUpdatedSolutionVector(totalSolVec, tStep);    // => x, m, gam
    FloatArray vStress, vStrain;

    f.zero();
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRuleL = integrationRulesArray [ layer - 1 ];

        for ( int j = 1; j <= iRuleL->giveNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRuleL->getIntegrationPoint(j - 1);
            lCoords = *gp->giveCoordinates();
            this->computeBmatrixAt(gp, B);
            this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, a);
            vStrain.beProductOf(B, a);
            this->giveStructuralCrossSection()->giveRealStresses(vStress, gp, vStrain, tStep);



            this->computeSectionalForcesAt(sectionalForces, vStress, gp, tStep); // these are per unit volume
            

            // Computation of nodal forces: f = B^t*[N M T Ms Ts]^t
            ftemp.beTProductOf(B,sectionalForces);
            double dV = this->computeVolumeAround(gp);
            f.add(dV, ftemp);
        }
    }

    answer.resize( ndofs );
    answer.zero();

    //answer.assemble(f, ordering_all);
}

void 
BrouzoulisShell :: computeSectionalForcesAt(FloatArray &sectionalForces, FloatArray &stress, GaussPoint *gp, TimeStep *tStep){

    FloatArray vStress;

    
    FloatArray PG1(3), PG2(3), PG3(3);
    FloatArray cartStressVector, contravarStressVector;
    FloatMatrix lambda[3];

    FloatMatrix Gcov, Gcon; 
    FloatArray lCoords = *gp->giveCoordinates();
    //this->evalInitialContravarBaseVectorsAt(lCoords, Gcon);
    FloatMatrix P, Sig;
    P.beMatrixFormOfStress(vStress);
    Sig.beProductTOf(Gcon,P);       // Stress resultants stored in each column
    PG1.beColumnOf(Sig,1);
    PG2.beColumnOf(Sig,2);
    PG3.beColumnOf(Sig,3);
            
    double z = gp->giveCoordinate(3) * this->giveStructuralCrossSection()->give(CS_Thickness, gp);
    this->computeLambdaMatrices(lambda, z); // associated with the variation of the test functions   
    
    // f = lambda_1^T * PG1 + lambda_2^T * PG2 + lambda_3^T * PG3
    sectionalForces.resize(18);
    FloatArray temp;
    sectionalForces.zero();
    temp.beTProductOf(lambda [0], PG1);
    sectionalForces.add(temp);
    temp.beTProductOf(lambda [1], PG2);
    sectionalForces.add(temp);
    temp.beTProductOf(lambda [2], PG3);
    sectionalForces.add(temp);
}

void
BrouzoulisShell :: computeLambdaMatrices(FloatMatrix lambda [ 3 ], double zeta)
{
    // computes the lambda^g matrices associated with the variation and linearization 
    // of the base vectors g_i.

    // thickness coefficients
    double a = zeta ;
    double b = 0.5 * zeta * zeta ;

    // lambda1 =  ( I,   0,  a*I,   0 ,  b*I,   0 ,  0,  0 )
    // lambda2 =  ( 0,   I,   0 ,  a*I,   0 ,  b*I,  0,  0 )
    FloatMatrix eye(3,3), aEye(3,3), bEye(3,3);
    eye.beUnitMatrix();
    aEye=eye;  aEye.times(a);
    bEye=eye;  bEye.times(b);
    lambda[ 0 ].resize(3,24);   lambda[ 0 ].zero();  
    lambda[ 1 ].resize(3,24);   lambda[ 1 ].zero();
    lambda[ 0 ].setSubMatrix(eye,1,1);  lambda[ 0 ].setSubMatrix(aEye,1,7);  lambda[ 0 ].setSubMatrix(aEye,1,13);
    lambda[ 1 ].setSubMatrix(eye,1,4);  lambda[ 1 ].setSubMatrix(aEye,1,10); lambda[ 1 ].setSubMatrix(aEye,1,16);
    

    // lambda3 =  ( 0,  0,  0 , 0 , 0 , 0 , I , aI )
    lambda[ 2 ].resize(3,24);   lambda[ 2 ].zero();
    lambda[ 2 ].setSubMatrix(eye,1,19);  lambda[ 2 ].setSubMatrix(aEye,1,22); 

}


void
BrouzoulisShell :: evalInitialCovarBaseVectorsAt(FloatArray &lCoords, FloatMatrix &Gcov)
{
    double zeta = 0.0; //lcoords.at(3) * this->giveStructuralCrossSection()->give(CS_Thickness, gp);
    FloatArray M, Xbar;
    FloatMatrix dNdxi12, dNdxi3;

    // In plane base vectors
    this->interpolationUV.evaldNdxi(dNdxi12, lCoords, FEIElementGeometryWrapper(this));
    this->interpolationW.evaldNdxi(dNdxi3, lCoords, FEIElementGeometryWrapper(this));

    FloatArray G1(3), G2(3); 
    G1.zero();
    G2.zero();
    for ( int i = 1; i <= 6; i++ ) {
        FloatArray &xbar = * this->giveNode(i)->giveCoordinates();
        M = this->giveInitialNodeDirector(i);
        Xbar = this->giveInitialNodeMidSurface(i);
        FloatArray nodeCoords = (Xbar + zeta*M);
        G1.at(1) += dNdxi12.at(i, 1) * nodeCoords.at(1);
        G1.at(2) += dNdxi12.at(i, 1) * nodeCoords.at(2);
        G1.at(3) += dNdxi3.at(i, 1)  * nodeCoords.at(3);
        G2.at(1) += dNdxi12.at(i, 2) * nodeCoords.at(1);
        G2.at(2) += dNdxi12.at(i, 2) * nodeCoords.at(2);
        G2.at(3) += dNdxi3.at(i, 2)  * nodeCoords.at(3);
    }

    // Out of plane base vector = director
    FloatArray G3;
    //this->evalInitialDirectorAt(lcoords, G3);     // G3=M

    Gcov.resize(3,3);
    Gcov.setColumn(G1,1); Gcov.setColumn(G2,2); Gcov.setColumn(G3,3);
}

void 
BrouzoulisShell :: giveBmat(FloatArray &B, GaussPoint *gp)
{
    // Xbar = 0.5( X+ + X-)
    FloatMatrix dNdxibar;
    FloatArray lCoords = *gp->giveCoordinates();
    lCoords.at(3) = 0.0;  // Evaluate at the mid surface

    this->interpolationUV.evaldNdxi(dNdxibar, lCoords, FEIElementGeometryWrapper(this));
    this->interpolationW.evaldNdxi(dNdxibar, lCoords, FEIElementGeometryWrapper(this));
    //FloatMatrix N;

    //this->computeNmatrixAt(lCoords, N);


}



void
BrouzoulisShell :: setupInitialNodeDirectors()
{
    FloatArray lcoords(2), nodeLocalXiCoords, nodeLocalEtaCoords;

    // Give the local coordinates for the element nodes - all at once

    this->initialNodeDirectors.resize(6);
    this->initialNodeMidSurface.resize(6);
    FloatMatrix N;
    FloatArray Xbar, M, Xtop, Xbottom;
    for ( int node = 1; node <= 3; node++ ) {
        this->initialNodeDirectors [ node - 1 ].resize(3);
        this->initialNodeDirectors [ node - 1 ].zero();
        this->initialNodeMidSurface [ node - 1 ].resize(3);
        this->initialNodeMidSurface [ node - 1 ].zero();

        Xbar = 0.5 * ( *this->giveNode(node)->giveCoordinates() + *this->giveNode(3+node)->giveCoordinates() );
        M =          ( *this->giveNode(node)->giveCoordinates() - *this->giveNode(3+node)->giveCoordinates() );

        this->initialNodeMidSurface [ node - 1 ] = Xbar;
        this->initialNodeDirectors [ node - 1 ] = M;
    }

    for ( int node = 1; node <= 3; node++ ) {
        this->initialNodeDirectors [ node + 3 - 1 ].resize(3);
        this->initialNodeDirectors [ node + 3 - 1 ].zero();
        this->initialNodeMidSurface [ node + 3 - 1 ].resize(3);
        this->initialNodeMidSurface [ node + 3 - 1 ].zero();

        Xbar = 0.5 * ( *this->giveNode(6 + node)->giveCoordinates() + *this->giveNode(9 + node)->giveCoordinates() );
        M =          ( *this->giveNode(6 + node)->giveCoordinates() - *this->giveNode(9 + node)->giveCoordinates() );

        this->initialNodeMidSurface [ node + 3 - 1 ] = Xbar;
        this->initialNodeDirectors [ node + 3 - 1 ] = M;
    }

}


void
BrouzoulisShell :: giveLocalNodeCoords(FloatArray &nodeLocalXiCoords, FloatArray &nodeLocalEtaCoords)
{
    nodeLocalXiCoords.setValues(6, 1., 0., 0., .5, 0., .5);      // corner nodes then midnodes, uncertain of node numbering
    nodeLocalEtaCoords.setValues(6, 0., 1., 0., .5, .5, 0.);
}



} // end namespace oofem

#endif
