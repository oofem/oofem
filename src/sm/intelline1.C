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

#include "intelline1.h"
#include "node.h"
#include "structuralinterfacecrosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "fei2dlinelin.h"
#include "classfactory.h"


namespace oofem {

REGISTER_Element( IntElLine1 );

FEI2dLineLin IntElLine1 :: interp(1, 1);


IntElLine1 :: IntElLine1(int n, Domain *aDomain) :
    StructuralInterfaceElement(n, aDomain)
{
    numberOfDofMans = 4;
}


void
IntElLine1 :: computeNmatrixAt(GaussPoint *ip, FloatMatrix &answer)
{
    // Returns the modified N-matrix which multiplied with u give the spatial jump.

    FloatArray N;
    FEInterpolation *interp = this->giveInterpolation(); 
    interp->evalN( N, * ip->giveCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(2, 8);
    answer.zero();
    answer.at(1, 1) = answer.at(2, 2) = -N.at(1);
    answer.at(1, 3) = answer.at(2, 4) = -N.at(2);

    answer.at(1, 5) = answer.at(2, 6) = N.at(1);
    answer.at(1, 7) = answer.at(2, 8) = N.at(2);
}


void
IntElLine1 :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        //integrationRulesArray[0] = new LobattoIntegrationRule (1,domain, 1, 2);
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 2);
        integrationRulesArray [ 0 ]->SetUpPointsOnLine(2, _2dInterface); //@todo - should be a parameter with num of ip
    }
}

void
IntElLine1 :: computeCovarBaseVectorAt(IntegrationPoint *ip, FloatArray &G)
{
    FloatMatrix dNdxi;
    FEInterpolation *interp = this->giveInterpolation(); 
    interp->evaldNdxi(dNdxi, * ip->giveCoordinates(), FEIElementGeometryWrapper(this) );
    G.resize(2);
    G.zero();
    int numNodes = this->giveNumberOfNodes();
    for ( int i = 1; i <= dNdxi.giveNumberOfRows(); i++ ) {
        double X1_i = 0.5* ( this->giveNode(i)->giveCoordinate(1) + this->giveNode(i+numNodes/2)->giveCoordinate(1) ); // (mean) point on the fictious mid surface
        double X2_i = 0.5* ( this->giveNode(i)->giveCoordinate(2) + this->giveNode(i+numNodes/2)->giveCoordinate(2) );
        G.at(1) += dNdxi.at(i,1) * X1_i;
        G.at(2) += dNdxi.at(i,1) * X2_i;
    }

    
}

double
IntElLine1 :: computeAreaAround(IntegrationPoint *ip)
{
    FloatArray G;
    this->computeCovarBaseVectorAt(ip, G);

    double weight  = ip->giveWeight();
    double ds = sqrt( G.dotProduct(G) ) * weight;

    double thickness  = this->giveCrossSection()->give(CS_Thickness);    
    return ds * thickness;
}


IRResultType
IntElLine1 :: initializeFrom(InputRecord *ir)
{
    return StructuralInterfaceElement :: initializeFrom(ir);
}


void
IntElLine1 :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    answer.setValues(2, D_u, D_v);
}

void 
IntElLine1 :: computeTransformationMatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    FloatArray G;
    this->computeCovarBaseVectorAt(gp, G);
    G.normalize();

    answer.resize(2,2);
    answer.at(1,1) =  G.at(1);
    answer.at(2,1) =  G.at(2);
    answer.at(1,2) = -G.at(2);
    answer.at(2,2) =  G.at(1);

    /*
    answer.resize(3,3);
    answer.at(1,1) =  G.at(1);
    answer.at(2,1) =  G.at(2);
    answer.at(3,1) =  0.0;
    answer.at(1,2) = -G.at(2);
    answer.at(2,2) =  G.at(1);
    answer.at(3,2) =  0.0;
    answer.at(3,3) =  1.0;
    */
}

FEInterpolation *
IntElLine1 :: giveInterpolation() const
{
    return &interp;
}


} // end namespace oofem
