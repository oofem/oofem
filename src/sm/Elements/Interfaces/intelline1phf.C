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

#include "Elements/Interfaces/intelline1phf.h"
#include "CrossSections/structuralinterfacecrosssection.h"
#include "Materials/InterfaceMaterials/structuralinterfacematerialphf.h"
#include "node.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "mathfem.h"
#include "fei2dlinelin.h"
#include "classfactory.h"


namespace oofem {
REGISTER_Element(IntElLine1PhF);

FEI2dLineLin IntElLine1PhF :: interp(1, 1);


IntElLine1PhF :: IntElLine1PhF(int n, Domain *aDomain) : StructuralInterfaceElementPhF(n, aDomain)
{
    numberOfDofMans = 4;
    axisymmode = false;

    numberOfGaussPoints = 4;
}


void
IntElLine1PhF :: computeNmatrixAt(GaussPoint *ip, FloatMatrix &answer)
{
    // Returns the modified N-matrix which multiplied with u give the spatial jump.

    FloatArray N;
    FEInterpolation *interp = this->giveInterpolation();
    interp->evalN( N, ip->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(2, 8);
    answer.zero();
    answer.at(1, 1) = answer.at(2, 2) = -N.at(1);
    answer.at(1, 3) = answer.at(2, 4) = -N.at(2);

    answer.at(1, 5) = answer.at(2, 6) = N.at(1);
    answer.at(1, 7) = answer.at(2, 8) = N.at(2);
}


void
IntElLine1PhF :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        //integrationRulesArray[ 0 ].reset( new LobattoIntegrationRule (1,domain, 1, 2) );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 2) );
        integrationRulesArray [ 0 ]->SetUpPointsOnLine(this->numberOfGaussPoints, _2dInterface);
    }
}

void
IntElLine1PhF :: computeCovarBaseVectorAt(IntegrationPoint *ip, FloatArray &G)
{
    FloatMatrix dNdxi;
    FEInterpolation *interp = this->giveInterpolation();
    interp->evaldNdxi( dNdxi, ip->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    G.resize(2);
    G.zero();
    int numNodes = this->giveNumberOfNodes();
    for ( int i = 1; i <= dNdxi.giveNumberOfRows(); i++ ) {
        double X1_i = 0.5 * ( this->giveNode(i)->giveCoordinate(1) + this->giveNode(i + numNodes / 2)->giveCoordinate(1) ); // (mean) point on the fictious mid surface
        double X2_i = 0.5 * ( this->giveNode(i)->giveCoordinate(2) + this->giveNode(i + numNodes / 2)->giveCoordinate(2) );
        G.at(1) += dNdxi.at(i, 1) * X1_i;
        G.at(2) += dNdxi.at(i, 1) * X2_i;
    }
}

double
IntElLine1PhF :: computeAreaAround(IntegrationPoint *ip)
{
    FloatArray G;
    this->computeCovarBaseVectorAt(ip, G);

    double weight  = ip->giveWeight();
    double ds = sqrt( G.dotProduct(G) ) * weight;
    if ( this->axisymmode ) {
        int numNodes = this->giveNumberOfNodes();
        FloatArray N;
        this->interp.evalN( N, ip->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
        // interpolate radius
        double r = 0.0;
        for ( int i = 1; i <= N.giveSize(); i++ ) {
            double X_i = 0.5 * ( this->giveNode(i)->giveCoordinate(1) + this->giveNode(i + numNodes / 2)->giveCoordinate(1) ); // X-coord of the fictious mid surface
            r += N.at(i) * X_i;
        }
        return ds * r;
        
    } else { // regular 2d
        double thickness  = this->giveCrossSection()->give(CS_Thickness, ip);
        return ds * thickness;
    }
}


IRResultType
IntElLine1PhF :: initializeFrom(InputRecord *ir)
{
    this->axisymmode = false;
    this->axisymmode = ir->hasField(_IFT_IntElLine1PhF_axisymmode);
    return StructuralInterfaceElement :: initializeFrom(ir);
    
}


void
IntElLine1PhF :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_u, D_v, T_f};
}

void
IntElLine1PhF :: computeTransformationMatrixAt(GaussPoint *gp, FloatMatrix &answer)
{
    // Transformation matrix to the local coordinate system
    FloatArray G;
    this->computeCovarBaseVectorAt(gp, G);
    G.normalize();

    answer.resize(2, 2);
    answer.at(1, 1) =  G.at(1);
    answer.at(2, 1) = -G.at(2);
    answer.at(1, 2) =  G.at(2);
    answer.at(2, 2) =  G.at(1);

}

FEInterpolation *
IntElLine1PhF :: giveInterpolation() const
{
    return & interp;
}


void
IntElLine1PhF :: giveDofManDofIDMask_u(IntArray &answer)
{
    answer = {D_u, D_v};
}

void
IntElLine1PhF :: giveDofManDofIDMask_d(IntArray &answer)
{
    // use temperature for now
    answer = {T_f};   
}

void
IntElLine1PhF :: getLocationArray_u( IntArray &answer )
{
    answer = {1, 2,  4, 5,  7, 8,  10, 11};
}
void
IntElLine1PhF :: getLocationArray_d( IntArray &answer )
{
    answer = {3, 6, 9, 12};
}


void
IntElLine1PhF :: giveEngTraction(FloatArray &answer, GaussPoint *gp, const FloatArray &jump, const double damage, TimeStep *tStep)
{
    
    StructuralInterfaceMaterialPhF *mat = static_cast< StructuralInterfaceMaterialPhF * >( this->giveInterfaceCrossSection()->giveInterfaceMaterial() );
    //double damage = this->computeDamageAt(gp, VM_Total, tStep);
    
    //StructuralInterfaceMaterial *mat = this->giveInterfaceCrossSection()->giveInterfaceMaterial();
    mat->giveEngTraction_2d(answer, gp, jump, damage, tStep);
    //this->giveInterfaceCrossSection()->giveEngTraction_2d(answer, gp, jump, tStep);
}


void
IntElLine1PhF::computeCovarBaseVectorsAt(GaussPoint *gp, FloatMatrix &G)
{
    
    FloatArray G1;
    this->computeCovarBaseVectorAt(gp, G1);
    
    G.resize(2, 1);
    G.at(1, 1) =  G1.at(1);
    G.at(2, 1) =  G1.at(2);
    
    
}


} // end namespace oofem
