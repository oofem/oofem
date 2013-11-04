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

#include "intelline2.h"
#include "node.h"
#include "structuralinterfacecrosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "fei2dlinequad.h"
#include "classfactory.h"


namespace oofem {

REGISTER_Element( IntElLine2 );

FEI2dLineQuad IntElLine2 :: interp(1, 2);


IntElLine2 :: IntElLine2(int n, Domain *aDomain) :
    StructuralInterfaceElement(n, aDomain)
{
    numberOfDofMans = 6;
}


void
IntElLine2 :: computeNmatrixAt(GaussPoint *ip, FloatMatrix &answer)
{
    // Returns the modified N-matrix which multiplied with u give the spatial jump.

    FloatArray N;
    interp.evalN( N, * ip->giveCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(2, 12);
    answer.zero();

    answer.at(1, 2) = answer.at(2, 1) = -N.at(1);
    answer.at(1, 4) = answer.at(2, 3) = -N.at(2);
    answer.at(1, 6) = answer.at(2, 5) = -N.at(3);

    answer.at(1, 8) = answer.at(2, 7) = N.at(1);
    answer.at(1, 10) = answer.at(2, 9) = N.at(2);
    answer.at(1, 12) = answer.at(2, 11) = N.at(3);
}


void
IntElLine2 :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        //integrationRulesArray[0] = new LobattoIntegrationRule (1,domain, 1, 2);
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 2);
        integrationRulesArray [ 0 ]->SetUpPointsOnLine(4, _2dInterface); //@todo - should be a parameter with num of ip
    }
}

void
IntElLine2 :: computeCovarBaseVectorAt(IntegrationPoint *ip, FloatArray &G)
{
    FloatMatrix dNdxi;
    interp.evaldNdxi(dNdxi, * ip->giveCoordinates(), FEIElementGeometryWrapper(this) );
    FloatArray Xi;
    G.resize(2);
    G.zero();
    for ( int i = 1; i <= dNdxi.giveNumberOfRows(); i++ ) {
        Xi = 0.5* ( *this->giveNode(i)->giveCoordinates() + *this->giveNode(i+3)->giveCoordinates() ); // (mean) point on the fictious mid surface
        G.add(dNdxi.at(i,1), Xi); // Curvilinear base vector
    }
}

double
IntElLine2 :: computeAreaAround(IntegrationPoint *ip)
// Returns the length of the receiver. This method is valid only if 1
// Gauss point is used.
{
    FloatArray G;
    this->computeCovarBaseVectorAt(ip, G);

    double weight  = ip->giveWeight();
    double ds = sqrt( G.dotProduct(G) ) * weight;

    double thickness  = this->giveCrossSection()->give(CS_Thickness);    
    return ds * thickness;

#if 0
    // old code
    double weight  = ip->giveWeight();
    double ksi = ip->giveCoordinate(1);
    double dn1 = ksi - 0.5;
    double dn2 = ksi + 0.5;
    double dn3 = -2.0 * ksi;

    double x1 = this->giveNode(1)->giveCoordinate(1);
    double x2 = this->giveNode(2)->giveCoordinate(1);
    double x3 = this->giveNode(3)->giveCoordinate(1);

    double y1 = this->giveNode(1)->giveCoordinate(2);
    double y2 = this->giveNode(2)->giveCoordinate(2);
    double y3 = this->giveNode(3)->giveCoordinate(2);

    double dx = ( dn1 * x1 ) + ( dn2 * x2 ) + ( dn3 * x3 );
    double dy = ( dn1 * y1 ) + ( dn2 * y2 ) + ( dn3 * y3 );
    double thickness  = this->giveCrossSection()->give(CS_Thickness);
    return sqrt(dx * dx + dy * dy) * weight * thickness;
#endif
}


IRResultType
IntElLine2 :: initializeFrom(InputRecord *ir)
{
    return StructuralInterfaceElement :: initializeFrom(ir);
}


void
IntElLine2 :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    answer.setValues(2, D_u, D_v);
}

void 
IntElLine2 :: computeTransformationMatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    FloatArray G;
    this->computeCovarBaseVectorAt(gp, G);
    G.normalize();

    answer.resize(2,2);
    answer.setColumn(G,1);
    answer.at(1,2) = -G.at(2);
    answer.at(2,2) =  G.at(1);
}

bool
IntElLine2 :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    FloatArray grad(2);

    //double ksi = gp -> giveCoordinate(1) ;
    double ksi = 0.0; // compute tangent in the middle
    double dn1 = ksi - 0.5;
    double dn2 = ksi + 0.5;
    double dn3 = -2.0 * ksi;

    // tangent
    grad.at(1) = dn1 * this->giveNode(1)->giveCoordinate(1) + dn2 *this->giveNode(2)->giveCoordinate(1) + dn3 *this->giveNode(3)->giveCoordinate(1);
    grad.at(2) = dn1 * this->giveNode(1)->giveCoordinate(2) + dn2 *this->giveNode(2)->giveCoordinate(2) + dn3 *this->giveNode(3)->giveCoordinate(2);
    grad.normalize();

    answer.resize(12, 12);
    for ( int i = 0; i < 6; i++ ) {
        answer.at(i * 2 + 1, i * 2 + 1) = grad.at(1);
        answer.at(i * 2 + 1, i * 2 + 2) = grad.at(2);
        answer.at(i * 2 + 2, i * 2 + 1) = -grad.at(2);
        answer.at(i * 2 + 2, i * 2 + 2) = grad.at(1);
    }

    return 1;
}


FEInterpolation *
IntElLine2 :: giveInterpolation() const
{
    return &interp;
}


} // end namespace oofem
