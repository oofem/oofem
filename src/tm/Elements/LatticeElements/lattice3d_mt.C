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
 *               Copyright (C) 1993 - 2019   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "tm/Elements/LatticeElements/lattice3d_mt.h"
#include "tm/Materials/transportmaterial.h"
#include "tm/Materials/LatticeMaterials/latticetransmat.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "engngm.h"
#include "load.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif

namespace oofem {
REGISTER_Element(Lattice3d_mt);

Lattice3d_mt :: Lattice3d_mt(int n, Domain *aDomain, ElementMode em) :
    LatticeTransportElement(n, aDomain, em)
{
    this->numberOfDofMans  = 2;
}

double Lattice3d_mt :: giveLength()
{
    if ( this->geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    return this->length;
}

void Lattice3d_mt :: giveCrackLengths(FloatArray &lengths)
{
    if ( this->geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    lengths = this->crackLengths;
}


double
Lattice3d_mt :: giveArea()
{
    if ( this->geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    return this->area;
}

void
Lattice3d_mt :: computeNmatrixAt(FloatMatrix &answer, const FloatArray &coords)
{
    this->computeNSubMatrixAt(answer, coords);
}


void
Lattice3d_mt :: computeNSubMatrixAt(FloatMatrix &answer, const FloatArray &coords)
{
    double ksi = coords.at(1);
    double n1 = ( 1. - ksi ) * 0.5;
    double n2 = ( 1. + ksi ) * 0.5;

    answer.resize(1, 2);
    answer.zero();

    answer.at(1, 1) = n1;
    answer.at(1, 2) = n2;
}


void
Lattice3d_mt :: computeGradientMatrixAt(FloatMatrix &answer, const FloatArray &lcoords)
{
    double l = giveLength();

    answer.resize(1, 2);
    answer.zero();
    answer.at(1, 1) = -1.;
    answer.at(1, 2) = 1.;

    answer.times(1. / l);
}


void
Lattice3d_mt :: updateInternalState(TimeStep *tStep)
{
    FloatArray f, r;
    FloatMatrix n;
    TransportMaterial *mat = static_cast< TransportMaterial * >( this->giveMaterial() );

    // force updating ip values
    for ( auto &iRule: integrationRulesArray ) {
        for ( auto &gp: * iRule ) {
            this->computeNmatrixAt(n, gp->giveNaturalCoordinates() );
            this->computeVectorOf({ P_f }, VM_Total, tStep, r);
            f.beProductOf(n, r);
            mat->updateInternalState(f, gp, tStep);
        }
    }
}


void
Lattice3d_mt :: computeGaussPoints()
{
    integrationRulesArray.resize(1);
    integrationRulesArray [ 0 ] = std :: make_unique< GaussIntegrationRule >(1, this, 1, 2);
    integrationRulesArray [ 0 ]->SetUpPointsOnLine(1, _3dMTLattice);
}

void
Lattice3d_mt ::   giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = { P_f };
}

void
Lattice3d_mt :: initializeFrom(InputRecord &ir)
{
    this->Element :: initializeFrom(ir);

    numberOfGaussPoints = 1;

    minLength = 1.e-20;
    IR_GIVE_OPTIONAL_FIELD(ir, minLength, _IFT_Lattice3DMT_mlength);

    dimension = 3.;
    IR_GIVE_OPTIONAL_FIELD(ir, dimension, _IFT_Lattice3DMT_dim);


    IR_GIVE_OPTIONAL_FIELD(ir, area, _IFT_Lattice3DMT_area);


    polygonCoords.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, polygonCoords, _IFT_Lattice3DMT_polycoords);
    numberOfPolygonVertices = polygonCoords.giveSize() / 3.;

    crackWidths.resize(numberOfPolygonVertices);
    crackWidths.zero();
    IR_GIVE_OPTIONAL_FIELD(ir, crackWidths, _IFT_Lattice3DMT_crackwidths);

    couplingFlag = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, couplingFlag, _IFT_Lattice3DMT_couplingflag);

    couplingNumbers.resize(numberOfPolygonVertices);
    if ( couplingFlag == 1 ) {
        IR_GIVE_OPTIONAL_FIELD(ir, couplingNumbers, _IFT_Lattice3DMT_couplingnumber);
    }

    this->computeGaussPoints();
}

void
Lattice3d_mt :: computeFlow(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray r;
    IntArray dofid;
    double dV;
    double length = giveLength();
    double k;
    answer.resize(1);

    k = static_cast< TransportMaterial * >( this->giveMaterial() )->giveCharacteristicValue(Conductivity, gp, tStep);
    dV      = this->computeVolumeAround(gp);

    this->giveElementDofIDMask(dofid);
    this->computeVectorOf(dofid, VM_Total, tStep, r);
    double flow;
    double dP;

    dP = r.at(2) - r.at(1);
    flow = k * ( dP ) * dV / pow(length, 2.);
    answer.at(1) = fabs(flow);

    if ( !isActivated(tStep) ) {
        answer.zero();
    }
}

double
Lattice3d_mt :: computeVolumeAround(GaussPoint *aGaussPoint)
{
    return this->giveArea() * this->giveLength();
}


void
Lattice3d_mt :: computeGeometryProperties()
{
    //coordinates of the two nodes
    Node *nodeA, *nodeB;
    FloatArray coordsA(3), coordsB(3);

    nodeA  = this->giveNode(1);
    nodeB  = this->giveNode(2);

    for ( int i = 0; i < 3; i++ ) {
        coordsA.at(i + 1) =  nodeA->giveCoordinate(i + 1);
        coordsB.at(i + 1) =  nodeB->giveCoordinate(i + 1);
    }

    //Construct an initial temporary local coordinate system
    FloatArray s(3), t(3);

    //Calculate normal vector
    this->normal.resize(3);
    for ( int i = 0; i < 3; i++ ) {
        this->normal.at(i + 1) = coordsB.at(i + 1) - coordsA.at(i + 1);
    }

    // Compute midpoint
    this->midPoint.resize(3);
    for ( int i = 0; i < 3; i++ ) {
        this->midPoint.at(i + 1) = 0.5 * ( coordsB.at(i + 1) + coordsA.at(i + 1) );
    }

    this->length  = sqrt(pow(normal.at(1), 2.) + pow(normal.at(2), 2.) + pow(normal.at(3), 2.) );

    if ( this->length < this->minLength ) {
        this->length = this->minLength;
        printf("Length could be zero. Might not be possible to create local coordinate system. Need to calculate area differently. Crack lengths need to approximated.\n");
        //This will only work for cross-sections with three points at the moment.
        computeSpecialCrossSectionProperties();
        this->geometryFlag = 1;
        return;
    }

    for ( int i = 0; i < 3; i++ ) {
        this->normal.at(i + 1) /= length;
    }

    computeCrossSectionProperties();

    this->geometryFlag = 1;

    return;
}

void
Lattice3d_mt :: computeSpecialCrossSectionProperties()
{
    if ( this->numberOfPolygonVertices != 3 ) {
        return;
    }

    //Here, the area of the cross-section is calculated from the global vertex points without the generation of a local coordinate system.

    FloatArray pointOne(3);
    FloatArray pointTwo(3);
    FloatArray crackPoint(3);

    FloatArray l(3);
    this->crackLengths.resize(3);
    this->crackLengths.zero();

    for ( int m = 0; m < numberOfPolygonVertices; m++ ) {
        if ( m < numberOfPolygonVertices - 1 ) {
            for ( int n = 0; n < 3; n++ ) {
                pointOne.at(n + 1) = polygonCoords.at(3 * m + n + 1);
                pointTwo.at(n + 1) = polygonCoords.at(3 * ( m + 1 ) + n + 1);
            }
        } else {
            for ( int n = 0; n < 3; n++ ) {
                pointOne.at(n + 1) = polygonCoords.at(3 * m + n + 1);
                pointTwo.at(n + 1) = polygonCoords.at(n + 1);
            }
        }
        l.at(m + 1) = sqrt(pow(pointOne.at(1) - pointTwo.at(1), 2.) + pow(pointOne.at(2) - pointTwo.at(2), 2.) + pow(pointOne.at(3) - pointTwo.at(3), 2.) );
    }

    double s = ( l.at(1) + l.at(2) + l.at(3) ) / 2.;
    this->area = sqrt(s * ( s - l.at(1) ) * ( s - l.at(2) ) * ( s - l.at(3) ) );

    printf("area = %e, length =  %e\n", area, length);

    if ( this->area < pow(this->minLength, 2.) ) {
        this->area = pow(this->minLength, 2.);
        printf("small area fixed in special calculation\n");
    }

    //Approximate crack lengths assuming
    //an equilateral triangle of the same area as real triangle has.
    for ( int m = 0; m < numberOfPolygonVertices; m++ ) {
        this->crackLengths.at(m + 1) = pow(3, 0.25) * sqrt(this->area);
    }

    return;
}



// void
// Lattice3d_mt :: computeHomogenisedInternalForcesVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode, FloatArray &unknowns)
// {
//     IntArray dofids;
//     FloatArray tmp;

//     FloatMatrix s;
//     this->giveElementDofIDMask(dofids);

//     ///@todo Integrate and compute as a nonlinear problem instead of doing this tangent.
//     this->computeConductivityMatrix(s, Conductivity, tStep);
//     answer.beProductOf(s, unknowns);

//     this->computeInternalSourceRhsVectorAt(tmp, tStep, mode);
//     answer.subtract(tmp);

//     FloatMatrix bc_tangent;
//     this->computeBCMtrxAt(bc_tangent, tStep, VM_Total);
//     if ( bc_tangent.isNotEmpty() ) {
//         tmp.beProductOf(bc_tangent, unknowns);
//         answer.add(tmp);
//     }
// }


void
Lattice3d_mt :: computeCrossSectionProperties() {
    if ( this->numberOfPolygonVertices < 3 ) {
        return;
    }

    //Construct two perpendicular axis so that n is normal to the plane which they create
    //Check, if one of the components of the normal-direction is zero
    FloatArray s(3), t(3);
    if ( this->normal.at(1) == 0 ) {
        s.at(1) = 0.;
        s.at(2) = this->normal.at(3);
        s.at(3) = -this->normal.at(2);
    } else if ( this->normal.at(2) == 0 ) {
        s.at(1) = this->normal.at(3);
        s.at(2) = 0.;
        s.at(3) = -this->normal.at(1);
    } else {
        s.at(1) = this->normal.at(2);
        s.at(2) = -this->normal.at(1);
        s.at(3) = 0.;
    }

    s.normalize();

    t.beVectorProductOf(this->normal, s);
    t.normalize();

    //Set up rotation matrix
    FloatMatrix lcs(3, 3);

    for ( int i = 1; i <= 3; i++ ) {
        lcs.at(1, i) = this->normal.at(i);
        lcs.at(2, i) = s.at(i);
        lcs.at(3, i) = t.at(i);
    }


    //Calculate the local coordinates of the polygon vertices
    FloatArray help(3), test(3);
    FloatArray lpc(3 * numberOfPolygonVertices);
    for ( int k = 0; k < numberOfPolygonVertices; k++ ) {
        for ( int n = 0; n < 3; n++ ) {
            help(n) = polygonCoords(3 * k + n);
        }

        test.beProductOf(lcs, help);
        for ( int n = 0; n < 3; n++ ) {
            lpc(3 * k + n) = test(n);
        }
    }

    this->area = 0.;

    for ( int k = 0; k < numberOfPolygonVertices; k++ ) {
        if ( k < numberOfPolygonVertices - 1 ) {
            this->area += lpc(3 * k + 1) * lpc(3 * ( k + 1 ) + 2) - lpc(3 * ( k + 1 ) + 1) * lpc(3 * k + 2);
        } else {   //Back to zero for n+1
            this->area += lpc(3 * k + 1) * lpc(2) - lpc(1) * lpc(3 * k + 2);
        }
    }

    this->area *= 0.5;

    FloatArray tempCoords(3 * numberOfPolygonVertices);
    if ( this->area < 0 ) { //Set area to a positive value and rearrange the coordinate entries
        this->area *= -1.;

        for ( int k = 0; k < numberOfPolygonVertices; k++ ) {
            for ( int m = 0; m < 3; m++ ) {
                tempCoords.at(3 * k + m + 1) = polygonCoords.at(3 * ( numberOfPolygonVertices - k - 1 ) + m + 1);
            }
        }

        polygonCoords = tempCoords;

        // Calculate again local co-ordinate system for different order
        for ( int k = 0; k < numberOfPolygonVertices; k++ ) {
            for ( int n = 0; n < 3; n++ ) {
                help(n) = polygonCoords(3 * k + n);
            }

            test.beProductOf(lcs, help);
            for ( int n = 0; n < 3; n++ ) {
                lpc(3 * k + n) = test(n);
            }
        }
    }

    if ( this->area < pow(minLength, 2.) ) {
        this->area = pow(minLength, 2.);
        printf("Small area fixed.\n");
    }

    //Calculate centroids
    centroid.resize(3);
    for ( int k = 0; k < numberOfPolygonVertices; k++ ) {
        if ( k < numberOfPolygonVertices - 1 ) {
            centroid.at(2) += ( lpc(3 * k + 1) + lpc(3 * ( k + 1 ) + 1) ) * ( lpc(3 * k + 1) * lpc(3 * ( k + 1 ) + 2) - lpc(3 * ( k + 1 ) + 1) * lpc(3 * k + 2) );
            centroid.at(3) += ( lpc(3 * k + 2) + lpc(3 * ( k + 1 ) + 2) ) * ( lpc(3 * k + 1) * lpc(3 * ( k + 1 ) + 2) - lpc(3 * ( k + 1 ) + 1) * lpc(3 * k + 2) );
        } else { //Back to zero for n+1
            centroid.at(2) += ( lpc(3 * k + 1) + lpc(1) ) * ( lpc(3 * k + 1) * lpc(2) - lpc(1) * lpc(3 * k + 2) );
            centroid.at(3) += ( lpc(3 * k + 2) + lpc(2) ) * ( lpc(3 * k + 1) * lpc(2) - lpc(1) * lpc(3 * k + 2) );
        }
    }

    centroid.times(1. / ( 6. * this->area ) );

    centroid.at(1) = lpc.at(1); //The first component of all lpcs should be the same

    //Shift coordinates to centroi
    for ( int k = 0; k < numberOfPolygonVertices; k++ ) {
        for ( int l = 0; l < 3; l++ ) {
            lpc(3 * k + l) -= centroid(l);
        }
    }

    //Compute second moments of area.

    //This is for the temporary coordinate system
    double Ixx = 0.;
    double Iyy = 0.;
    double Ixy = 0.;
    double a;

    for ( int k = 0; k < numberOfPolygonVertices; k++ ) {
        if ( k < numberOfPolygonVertices - 1 ) {
            a = lpc(3 * k + 1) * lpc(3 * ( k + 1 ) + 2) - lpc(3 * ( k + 1 ) + 1) * lpc(3 * k + 2);

            Ixx += ( ( pow(lpc(3 * k + 2), 2.) + lpc(3 * k + 2) * lpc(3 * ( k + 1 ) + 2) + pow(lpc(3 * ( k + 1 ) + 2), 2.) ) * a ) / 12.;

            Iyy += ( ( pow(lpc(3 * k + 1), 2.) + lpc(3 * k + 1) * lpc(3 * ( k + 1 ) + 1) + pow(lpc(3 * ( k + 1 ) + 1), 2.) ) * a ) / 12.;

            Ixy += ( ( lpc(3 * k + 1) * lpc(3 * ( k + 1 ) + 2) + 2. * lpc(3 * k + 1) * lpc(3 * k + 2) +
                       2 * lpc(3 * ( k + 1 ) + 1) * lpc(3 * ( k + 1 ) + 2) + lpc(3 * ( k + 1 ) + 1) * lpc(3 * k + 2) ) * a ) / 24.;
        } else {   //Back to zero for n+1
            a = lpc(3 * k + 1) * lpc(2) - lpc(1) * lpc(3 * k + 2);

            Ixx += ( ( pow(lpc(3 * k + 2), 2.) + lpc(3 * k + 2) * lpc(2) + pow(lpc(2), 2.) ) * a ) / 12.;

            Iyy += ( ( pow(lpc(3 * k + 1), 2.) + lpc(3 * k + 1) * lpc(1) + pow(lpc(1), 2.) ) * a ) / 12.;

            Ixy += ( ( lpc(3 * k + 1) * lpc(2) + 2. * lpc(3 * k + 1) * lpc(3 * k + 2) +
                       2 * lpc(1) * lpc(2) + lpc(1) * lpc(3 * k + 2) ) * a ) / 24.;
        }
    }

    //Compute main axis of the cross-section
    double angleChange = 0.;
    double sum = fabs(Ixx + Iyy);
    double pi = 3.14159265;
    if ( ( fabs(Ixx - Iyy) / sum > 1.e-6 ) && fabs(Ixy) / sum > 1.e-6 ) {
        angleChange = 0.5 * atan(-2 * Ixy / ( Ixx - Iyy ) );
    } else if ( ( fabs(Ixx - Iyy) / sum < 1.e-6 ) && fabs(Ixy) / sum > 1.e-6 ) {
        angleChange = pi / 4.;
    }

    if ( Iyy > Ixx ) {
        angleChange = angleChange + pi / 2.;
    }

    //Moment of inertias saved in the element

    this->I1 = ( Ixx + Iyy ) / 2. + sqrt(pow( ( Ixx - Iyy ) / 2., 2. ) + pow(Ixy, 2.) );
    this->I2 = ( Ixx + Iyy ) / 2. - sqrt(pow( ( Ixx - Iyy ) / 2., 2. ) + pow(Ixy, 2.) );

    this->Ip = I1 + I2;

    //Rotation around normal axis by angleChange

    FloatMatrix rotationChange(3, 3);
    rotationChange.zero();

    rotationChange.at(1, 1) = 1.;
    rotationChange.at(2, 2) = cos(angleChange);
    rotationChange.at(2, 3) = -sin(angleChange);

    rotationChange.at(3, 2) = sin(angleChange);
    rotationChange.at(3, 3) = cos(angleChange);

    this->localCoordinateSystem.beProductOf(rotationChange, lcs);

    //Calculate the polygon vertices in the new coordinate system
    for ( int k = 0; k < numberOfPolygonVertices; k++ ) {
        for ( int n = 0; n < 3; n++ ) {
            help(n) = polygonCoords(3 * k + n);
        }

        test.beProductOf(this->localCoordinateSystem, help);
        for ( int n = 0; n < 3; n++ ) {
            lpc(3 * k + n) = test(n);
        }
    }

    //Calculate centroid again in local coordinate system
    centroid.zero();
    for ( int k = 0; k < numberOfPolygonVertices; k++ ) {
        if ( k < numberOfPolygonVertices - 1 ) {
            centroid.at(2) += ( lpc(3 * k + 1) + lpc(3 * ( k + 1 ) + 1) ) * ( lpc(3 * k + 1) * lpc(3 * ( k + 1 ) + 2) - lpc(3 * ( k + 1 ) + 1) * lpc(3 * k + 2) );
            centroid.at(3) += ( lpc(3 * k + 2) + lpc(3 * ( k + 1 ) + 2) ) * ( lpc(3 * k + 1) * lpc(3 * ( k + 1 ) + 2) - lpc(3 * ( k + 1 ) + 1) * lpc(3 * k + 2) );
        } else {   //Back to zero for n+1
            centroid.at(2) += ( lpc(3 * k + 1) + lpc(1) ) * ( lpc(3 * k + 1) * lpc(2) - lpc(1) * lpc(3 * k + 2) );
            centroid.at(3) += ( lpc(3 * k + 2) + lpc(2) ) * ( lpc(3 * k + 1) * lpc(2) - lpc(1) * lpc(3 * k + 2) );
        }
    }

    centroid.times(1. / ( 6. * area ) );

    centroid.at(1) = lpc.at(1); //The first component of all lpcs should be the same

    /*Express midpoint in local coordinate system.
     * Note that this point does not have to lie in the same plane as the cross-section.
     * However, we do not use the first componnent of the coordinate */
    FloatArray midPointLocal(3);
    midPointLocal.beProductOf(this->localCoordinateSystem, midPoint);

    //eccentricities stored in the element
    this->eccS = centroid.at(2) - midPointLocal.at(2);
    this->eccT = centroid.at(3) - midPointLocal.at(3);

    FloatMatrix transposeLCS;
    transposeLCS.beTranspositionOf(this->localCoordinateSystem);

    globalCentroid.beProductOf(transposeLCS, centroid);

    crackLengths.resize(numberOfPolygonVertices);
    double crackPointOne, crackPointTwo;
    for ( int m = 0; m < numberOfPolygonVertices; m++ ) {
        if ( m < numberOfPolygonVertices - 1 ) {
            crackPointOne = ( lpc(3 * ( m + 1 ) + 1) + lpc(3 * ( m ) + 1) ) / 2.;
            crackPointTwo = ( lpc(3 * ( m + 1 ) + 2) + lpc(3 * ( m ) + 2) ) / 2.;
        } else {
            crackPointOne = ( lpc(1) + lpc(3 * ( m ) + 1) ) / 2.;
            crackPointTwo = ( lpc(2) + lpc(3 * ( m ) + 2) ) / 2.;
        }

        crackLengths.at(m + 1) = sqrt(pow(crackPointOne - midPointLocal.at(2), 2.) + pow(crackPointTwo - midPointLocal.at(3), 2.) );
    }

    return;
}

void
Lattice3d_mt :: computeConductivityMatrix(FloatMatrix &answer, MatResponseMode rmode, TimeStep *tStep)
{
    GaussPoint *gp =  integrationRulesArray [ 0 ]->getIntegrationPoint(0);

    answer.resize(2, 2);
    answer.zero();
    answer.at(1, 1) = 1.;
    answer.at(1, 2) = -1.;
    answer.at(2, 1) = -1;
    answer.at(2, 2) = 1.;

    double length = giveLength();
    double k = static_cast< TransportMaterial * >( this->giveMaterial() )->giveCharacteristicValue(Conductivity, gp, tStep);
    double dV = this->computeVolumeAround(gp);
    double temp = k * dV / pow(length, 2.);
    answer.times(temp);

    return;
}


void
Lattice3d_mt :: computeCapacityMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    GaussPoint *gp =  integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    answer.resize(2, 2);
    answer.zero();
    answer.at(1, 1) = 2.;
    answer.at(1, 2) = 1.;
    answer.at(2, 1) = 1.;
    answer.at(2, 2) = 2.;
    double c = static_cast< TransportMaterial * >( this->giveMaterial() )->giveCharacteristicValue(Capacity, gp, tStep);
    double dV = this->computeVolumeAround(gp) / ( 6.0 * this->dimension );
    answer.times(c * dV);
}



void
Lattice3d_mt :: computeInternalSourceRhsVectorAt(FloatArray &answer, TimeStep *atTime, ValueModeType mode)
{
    int i, j, n, nLoads;
    double dV;
    bcGeomType ltype;
    Load *load;
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

    FloatArray deltaX(3), normalVector(3);
    FloatArray val, helpLoadVector, globalIPcoords;
    FloatMatrix nm;
    answer.resize(0);

    FloatArray gravityHelp(2);

    nLoads    = this->giveBodyLoadArray()->giveSize();
    for ( i = 1; i <= nLoads; i++ ) {
        n     = bodyLoadArray.at(i);
        load  = ( Load * ) domain->giveLoad(n);
        ltype = load->giveBCGeoType();

        if ( ltype == GravityPressureBGT ) {
            //Compute change of coordinates

            if ( this->geometryFlag == 0 ) {
                computeGeometryProperties();
            }

            deltaX.at(1) = this->normal.at(1) * this->length;
            deltaX.at(2) = this->normal.at(2) * this->length;
            deltaX.at(3) = this->normal.at(3) * this->length;

            gravityHelp.at(1) = 1.;
            gravityHelp.at(2) = -1.;

            dV  = this->computeVolumeAround(gp);
            load->computeValueAt(val, atTime, deltaX, mode);

            double k = static_cast< TransportMaterial * >( this->giveMaterial() )->giveCharacteristicValue(Conductivity, gp, atTime);

            double helpFactor = val.at(1) * k * dV;

            helpFactor /= pow(this->giveLength(), 2.);
            gravityHelp.times(helpFactor);

            if ( helpLoadVector.isEmpty() ) {
                helpLoadVector.resize(gravityHelp.giveSize() );
            }

            for ( j = 1; j <= gravityHelp.giveSize(); j++ ) {
                helpLoadVector.at(j) += gravityHelp.at(j);
            }
        }



        answer.add(helpLoadVector);
    }

    return;
}

int
Lattice3d_mt :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    if ( geometryFlag == 0 ) {
        computeGeometryProperties();
    }

    answer.resize(3);
    answer = midPoint;

    return 1;
}

bool
Lattice3d_mt :: computeLocalCoordinates(FloatArray &answer, const FloatArray &coords)
{
    answer.resize(1);
    answer.at(1) = 0.;

    return 1;
}


#ifdef __OOFEG

void
Lattice3d_mt :: drawYourself(oofegGraphicContext &gc, TimeStep *tStep)
{
    OGC_PlotModeType mode = gc.giveIntVarPlotMode();

    if ( mode == OGC_rawGeometry ) {
        this->drawRawGeometry(gc, tStep);
        this->drawRawCrossSections(gc, tStep);
    } else {
        OOFEM_ERROR("drawYourself : unsupported mode");
    }
}


void Lattice3d_mt :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    GraphicObj *go;

    WCRec p [ 2 ]; /* points */
    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor(gc.getElementColor() );
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);

    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(3);
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveCoordinate(3);

    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void Lattice3d_mt :: drawRawCrossSections(oofegGraphicContext &gc, TimeStep *tStep)
{
    GraphicObj *go;

    //  if (!go) { // create new one
    //Create as many points as we have polygon vertices
    WCRec p [ numberOfPolygonVertices ]; /* poin */

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor(gc.getElementColor() );
    EASValsSetLayer(OOFEG_RAW_CROSSSECTION_LAYER);

    for ( int i = 0; i < numberOfPolygonVertices; i++ ) {
        p [ i ].x = ( FPNum ) polygonCoords(3 * i);
        p [ i ].y = ( FPNum ) polygonCoords(3 * i + 1);
        p [ i ].z = ( FPNum ) polygonCoords(3 * i + 2);
    }

    WCRec pTemp [ 2 ]; /* points */
    for ( int i = 0; i < numberOfPolygonVertices; i++ ) {
        if ( i < numberOfPolygonVertices - 1 ) {
            pTemp [ 0 ] = p [ i ];
            pTemp [ 1 ] = p [ i + 1 ];
        } else {
            pTemp [ 0 ] = p [ i ];
            pTemp [ 1 ] = p [ 0 ];
        }

        go = CreateLine3D(pTemp);
        EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
        EGAttachObject(go, ( EObjectP ) this);
        EMAddGraphicsToModel(ESIModel(), go);
    }
}
#endif
} // end namespace oofem
