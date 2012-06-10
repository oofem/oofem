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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#include "rershell.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "domain.h"
#include "engngm.h"
#include "load.h"
#include "structuralcrosssection.h"
#include "mathfem.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
#endif

namespace oofem {
RerShell :: RerShell(int n, Domain *aDomain) :
    CCTPlate(n, aDomain)
{
    numberOfDofMans  = 3;
    GtoLRotationMatrix = NULL;
    Rx = 1.e+40;
    Ry = 1.e+40;
    Rxy = 1.e+40;

    numberOfGaussPoints = 1;
}


Interface *
RerShell :: giveInterface(InterfaceType interface)
{
    if ( interface == LayeredCrossSectionInterfaceType ) {
        return ( LayeredCrossSectionInterface * ) this;
    } else if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return ( ZZNodalRecoveryModelInterface * ) this;
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return ( NodalAveragingRecoveryModelInterface * ) this;
    }

    return NULL;
}


void
RerShell :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the [8x18] strain-displacement matrix {B} of the receiver,
// evaluated at aGaussPoint.
{
    int i, j, k, ii, ki;
    FloatArray b(3), c(3), nodeCoords;
    double x1, x2, x3, y1, y2, y3, area;

    this->giveLocalCoordinates( nodeCoords, * ( this->giveNode(1)->giveCoordinates() ) );
    x1 = nodeCoords.at(1);
    y1 = nodeCoords.at(2);

    this->giveLocalCoordinates( nodeCoords, * ( this->giveNode(2)->giveCoordinates() ) );
    x2 = nodeCoords.at(1);
    y2 = nodeCoords.at(2);

    this->giveLocalCoordinates( nodeCoords, * ( this->giveNode(3)->giveCoordinates() ) );
    x3 = nodeCoords.at(1);
    y3 = nodeCoords.at(2);


    b.at(1) = y2 - y3;
    b.at(2) = y3 - y1;
    b.at(3) = y1 - y2;

    c.at(1) = x3 - x2;
    c.at(2) = x1 - x3;
    c.at(3) = x2 - x1;

    area = this->giveArea();

    answer.resize(8, 18);
    answer.zero();

    for ( i = 1; i <= 3; i++ ) {
        j = i + 1 - i / 3 * 3;
        k = j + 1 - j / 3 * 3;

        ii = 6 * ( i - 1 );
        ki = 6 * ( k - 1 );

        answer.at(1, ii + 1) = b.at(i); // Eps_X
        answer.at(1, ki + 4) = -( b.at(i) - b.at(j) ) / Rx * area / 12.;
        answer.at(1, ki + 5) = -( c.at(i) - c.at(j) ) / Ry * area / 12.;

        answer.at(2, ii + 2) = c.at(i); // Eps_y
        answer.at(2, ki + 4) = -( b.at(i) - b.at(j) ) / Ry * area / 12.;
        answer.at(2, ki + 5) = -( c.at(i) - c.at(j) ) / Ry * area / 12.;

        answer.at(3, ii + 1) = c.at(i); // Gamma_xy;
        answer.at(3, ii + 2) = b.at(i);
        answer.at(3, ki + 4) = -( b.at(i) - b.at(j) ) / Rxy * area / 6.;
        answer.at(3, ki + 5) = -( c.at(i) - c.at(j) ) / Rxy * area / 6.;

        answer.at(4, ii + 5) = b.at(i); // Kappa_x

        answer.at(5, ii + 4) = -c.at(i); // Kappa_y

        answer.at(6, ii + 4) = -b.at(i); // Kappa_xy
        answer.at(6, ii + 5) = c.at(i);

        answer.at(7, ii + 3) = b.at(i); // Gamma_zx
        answer.at(7, ii + 5) += area * ( 2. / 3. );
        answer.at(7, ki + 5) += ( -2.0 * area + b.at(k) * ( c.at(i) - c.at(j) ) ) / 6.0;
        answer.at(7, ki + 4) = b.at(k) * ( b.at(i) - b.at(j) )            / 6.0;

        answer.at(8, ii + 3) = c.at(i); // Gamma_zy
        answer.at(8, ii + 4) += area * ( -2.0 / 3.0 );
        answer.at(8, ki + 4) += ( +2.0 * area + c.at(k) * ( b.at(i) - b.at(j) ) ) / 6.0;
        answer.at(8, ki + 5) = c.at(k) * ( c.at(i) - c.at(j) )            / 6.0;
    }

    answer.times( 1. / ( 2. * area ) );
}


void
RerShell :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 8);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Triangle, numberOfGaussPoints, _3dShell);
    }
}


void
RerShell :: computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver,
// evaluated at aGaussPoint.
{
    double x1, x2, x3, y1, y2, y3, l1, l2, l3, b1, b2, b3, c1, c2, c3;
    FloatArray nodeCoords;

    l1 = aGaussPoint->giveCoordinate(1);
    l2 = aGaussPoint->giveCoordinate(2);
    l3 = 1.0 - l1 - l2;

    this->giveLocalCoordinates( nodeCoords, * ( this->giveNode(1)->giveCoordinates() ) );
    x1 = nodeCoords.at(1);
    y1 = nodeCoords.at(2);

    this->giveLocalCoordinates( nodeCoords, * ( this->giveNode(2)->giveCoordinates() ) );
    x2 = nodeCoords.at(1);
    y2 = nodeCoords.at(2);

    this->giveLocalCoordinates( nodeCoords, * ( this->giveNode(3)->giveCoordinates() ) );
    x3 = nodeCoords.at(1);
    y3 = nodeCoords.at(2);

    /*
     * x1 = this -> giveNode(1) -> giveCoordinate(1);
     * x2 = this -> giveNode(2) -> giveCoordinate(1);
     * x3 = this -> giveNode(3) -> giveCoordinate(1);
     *
     * y1 = this -> giveNode(1) -> giveCoordinate(2);
     * y2 = this -> giveNode(2) -> giveCoordinate(2);
     * y3 = this -> giveNode(3) -> giveCoordinate(2);
     */
    b1 = y2 - y3;
    b2 = y3 - y1;
    b3 = y1 - y2;

    c1 = x3 - x2;
    c2 = x1 - x3;
    c3 = x2 - x1;

    answer.resize(5, 18);
    answer.zero();


    answer.at(1, 1) = l1;
    answer.at(1, 7) = l2;
    answer.at(1, 13) = l3;

    answer.at(2, 2) = l1;
    answer.at(2, 8) = l2;
    answer.at(2, 14) = l3;

    answer.at(3, 3) = l1;
    answer.at(3, 4) = -( l1 * l2 * b3 - l1 * l3 * b2 ) * 0.5;
    answer.at(3, 5) = -( l1 * l2 * c3 - l1 * l3 * c2 ) * 0.5;
    answer.at(3, 9) = l2;
    answer.at(3, 10) = -( l2 * l3 * b1 - l1 * l2 * b3 ) * 0.5;
    answer.at(3, 11) = -( l2 * l3 * c1 - l1 * l2 * c3 ) * 0.5;
    answer.at(3, 15) = l3;
    answer.at(3, 16) = -( l3 * l1 * b2 - l2 * l3 * b1 ) * 0.5;
    answer.at(3, 17) = -( l3 * l1 * c2 - l2 * l3 * c1 ) * 0.5;

    answer.at(4, 4) = l1;
    answer.at(4, 10) = l2;
    answer.at(4, 16) = l3;

    answer.at(5, 5) = l1;
    answer.at(5, 11) = l2;
    answer.at(5, 17) = l3;
}


double
RerShell :: giveArea()
// returns the area occupied by the receiver
{
    FloatArray nodeCoords;
    if ( area > 0 ) {
        return area;         // check if previously computed
    }

    double x1, x2, x3, y1, y2, y3;

    this->giveLocalCoordinates( nodeCoords, * ( this->giveNode(1)->giveCoordinates() ) );
    x1 = nodeCoords.at(1);
    y1 = nodeCoords.at(2);

    this->giveLocalCoordinates( nodeCoords, * ( this->giveNode(2)->giveCoordinates() ) );
    x2 = nodeCoords.at(1);
    y2 = nodeCoords.at(2);

    this->giveLocalCoordinates( nodeCoords, * ( this->giveNode(3)->giveCoordinates() ) );
    x3 = nodeCoords.at(1);
    y3 = nodeCoords.at(2);

    return ( area = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 ) );
}


IRResultType
RerShell :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    this->StructuralElement :: initializeFrom(ir);
    numberOfGaussPoints = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, IFT_RerShell_nip, "nip"); // Macro
    if ( numberOfGaussPoints != 1 ) {
        numberOfGaussPoints = 1;
    }

    this->computeGaussPoints();
    return IRRT_OK;
}


void
RerShell :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver.
{
    GaussPoint *gp;
    double dV, mss1;

    answer.resize(18, 18);
    answer.zero();
    gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

    dV = this->computeVolumeAround(gp);
    mss1 = dV * this->giveCrossSection()->give(CS_Thickness) * this->giveMaterial()->give('d', gp) / 3.;

    answer.at(1, 1) = mss1;
    answer.at(2, 2) = mss1;
    answer.at(3, 3) = mss1;

    answer.at(7, 7) = mss1;
    answer.at(8, 8) = mss1;
    answer.at(9, 9) = mss1;

    answer.at(13, 13) = mss1;
    answer.at(14, 14) = mss1;
    answer.at(15, 15) = mss1;
}


void
RerShell :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *stepN, ValueModeType mode)
// Computes numerically the load vector of the receiver due to the body
// loads, at stepN. - no support for momentum bodyload
{
    double dens, dV, load;
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    FloatArray f;
    FloatMatrix T;


    forLoad->computeComponentArrayAt(f, stepN, mode);
    //f.times( this->giveMaterial()->give('d') );

    if ( f.giveSize() == 0 ) {
        answer.resize(0);
        return;                                             // nil resultant
    } else {
        dens = this->giveMaterial()->give('d', gp);
        dV = this->computeVolumeAround(gp) * this->giveCrossSection()->give(CS_Thickness);

        answer.resize(18);

        load = f.at(1) * dens * dV / 3.0;
        answer.at(1) = load;
        answer.at(7) = load;
        answer.at(13) = load;

        load = f.at(2) * dens * dV / 3.0;
        answer.at(2) = load;
        answer.at(8) = load;
        answer.at(14) = load;

        load = f.at(3) * dens * dV / 3.0;
        answer.at(3) = load;
        answer.at(9) = load;
        answer.at(15) = load;

        // transform result from global cs to local element  cs.
        if ( this->computeGtoLRotationMatrix(T) ) {
            answer.rotatedWith(T, 'n');
        }

        return;
    }
}




FloatMatrix *
RerShell :: computeGtoLRotationMatrix()
// Returns the rotation matrix of the receiverof the size [3,3]
// coords(local) = T * coords(global)
//
// local coordinate (described by vector triplet e1',e2',e3') is defined as follows:
//
// e1'    : [N2-N1]    Ni - means i - th node
// help   : [N3-N1]
// e3'    : e1' x help
// e2'    : e3' x e1'
{
    int i;

    if ( GtoLRotationMatrix == NULL ) {
        FloatArray e1(3), e2, e3, help(3);

        for ( i = 1; i <= 3; i++ ) {
            e1.at(i) = ( this->giveNode(2)->giveCoordinate(i) - this->giveNode(1)->giveCoordinate(i) );
            help.at(i) = ( this->giveNode(3)->giveCoordinate(i) - this->giveNode(1)->giveCoordinate(i) );
        }

        // compute the norm of e1,help in order to normalize them
        e1.normalize();

        // compute vector product of e1' x help

        e3.beVectorProductOf(e1,help);

        // let us normalize e3'
        e3.normalize();

        // now from e3' x e1' compute e2'

        e2.beVectorProductOf(e3,e1);

        GtoLRotationMatrix = new FloatMatrix(3, 3);

        for ( i = 1; i <= 3; i++ ) {
            GtoLRotationMatrix->at(1, i) = e1.at(i);
            GtoLRotationMatrix->at(2, i) = e2.at(i);
            GtoLRotationMatrix->at(3, i) = e3.at(i);
        }
    }

    return GtoLRotationMatrix;
}


int
RerShell :: giveLocalCoordinateSystem(FloatMatrix &answer)
//
// returns a unit vectors of local coordinate system at element
// stored rowwise (mainly used by some materials with ortho and anisotrophy)
//
{
    this->computeGtoLRotationMatrix();
    answer = * GtoLRotationMatrix;
    return 1;
}


//converts global coordinates to local planar area coordinates, does not return a coordinate in the thickness direction, but
//does check that the point is in the element thickness
#define POINT_TOL 1.e-3
int
RerShell :: computeLocalCoordinates(FloatArray &answer, const FloatArray &coords)
{
    //set size of return value to 3 area coordinates
    answer.resize(3);

    //rotate the input point Coordinate System into the element CS
    FloatArray inputCoords_ElCS;
    this->giveLocalCoordinates( inputCoords_ElCS, const_cast< FloatArray & >(coords) );

    //Nodes are defined in the global CS, so they also need to be rotated into the element CS, therefore get the node points and
    //rotate them into the element CS
    FloatArray nodeCoords;
    double x1, x2, x3, y1, y2, y3, z1, z2, z3;

    this->giveLocalCoordinates( nodeCoords, * ( this->giveNode(1)->giveCoordinates() ) );
    x1 = nodeCoords.at(1);
    y1 = nodeCoords.at(2);
    z1 = nodeCoords.at(3);

    this->giveLocalCoordinates( nodeCoords, * ( this->giveNode(2)->giveCoordinates() ) );
    x2 = nodeCoords.at(1);
    y2 = nodeCoords.at(2);
    z2 = nodeCoords.at(3);

    this->giveLocalCoordinates( nodeCoords, * ( this->giveNode(3)->giveCoordinates() ) );
    x3 = nodeCoords.at(1);
    y3 = nodeCoords.at(2);
    z3 = nodeCoords.at(3);

    //Compute the area coordinates corresponding to this point
    double area;
    area = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 );

    answer.at(1) = ( ( x2 * y3 - x3 * y2 ) + ( y2 - y3 ) * inputCoords_ElCS.at(1) + ( x3 - x2 ) * inputCoords_ElCS.at(2) ) / 2. / area;
    answer.at(2) = ( ( x3 * y1 - x1 * y3 ) + ( y3 - y1 ) * inputCoords_ElCS.at(1) + ( x1 - x3 ) * inputCoords_ElCS.at(2) ) / 2. / area;
    answer.at(3) = ( ( x1 * y2 - x2 * y1 ) + ( y1 - y2 ) * inputCoords_ElCS.at(1) + ( x2 - x1 ) * inputCoords_ElCS.at(2) ) / 2. / area;

    //get midplane location at this point
    double midplZ;
    midplZ = z1 * answer.at(1) + z2 *answer.at(2) + z3 *answer.at(3);

    //check that the z is within the element
    StructuralCrossSection *cs;
    double elthick;

    cs = ( StructuralCrossSection * ) this->giveCrossSection();
    elthick = cs->give(CS_Thickness);

    if ( elthick / 2.0 + midplZ - fabs( inputCoords_ElCS.at(3) ) < -POINT_TOL ) {
        answer.zero();
        return 0;
    }

    //check that the point is in the element and set flag
    for ( int i = 1; i <= 3; i++ ) {
        if ( answer.at(i) < ( 0. - POINT_TOL ) ) {
            return 0;
        }

        if ( answer.at(i) > ( 1. + POINT_TOL ) ) {
            return 0;
        }
    }

    return 1;
}


bool
RerShell :: computeGtoLRotationMatrix(FloatMatrix &answer)
// Returns the rotation matrix of the receiver of the size [18,18]
// r(local) = T * r(global)
//
// local coordinate (described by vector triplet e1',e2',e3') is defined as follows:
//
// e1'    : [N2-N1]    Ni - means i - th node
// help   : [N3-N1]
// e3'    : e1' x help
// e2'    : e3' x e1'
{
    double val;
    int i, j;

    // test if pereviously computed
    if ( GtoLRotationMatrix == NULL ) {
        this->computeGtoLRotationMatrix();
    }


    answer.resize(18, 18);
    answer.zero();

    for ( i = 1; i <= 3; i++ ) {
        for ( j = 1; j <= 3; j++ ) {
            val = GtoLRotationMatrix->at(i, j);
            answer.at(i, j)     = val;
            answer.at(i + 3, j + 3) = val;
            answer.at(i + 6, j + 6) = val;
            answer.at(i + 9, j + 9) = val;
            answer.at(i + 12, j + 12) = val;
            answer.at(i + 15, j + 15) = val;
        }
    }

    return 1;
}


void
RerShell :: giveLocalCoordinates(FloatArray &answer, const FloatArray &global)
//
// Returns global coordinates given in global vector
// transformed into local coordinate system of the
// receiver
//
{
    // test the parameter
    if ( global.giveSize() != 3 ) {
        _error("GiveLocalCoordinate : cannot transform coordinates- size mismatch");
        exit(1);
    }

    // first ensure that receiver's GtoLRotationMatrix[3,3]
    // is defined

    this->computeGtoLRotationMatrix();
    answer.beProductOf(* GtoLRotationMatrix, global);
}


void
RerShell :: giveCharacteristicTensor(FloatMatrix &answer, CharTensor type, GaussPoint *gp, TimeStep *tStep)
//
// returns characteristic tensor of the receiver at given gp and tStep
//
{
    answer.resize(3, 3);
    answer.zero();

    if ( ( type == LocalForceTensor ) || ( type == GlobalForceTensor ) ) {
        FloatArray stress;
        this->computeStressVector(stress, gp, tStep);
        answer.at(1, 1) = stress.at(1);
        answer.at(2, 2) = stress.at(2);
        answer.at(1, 2) = stress.at(3);
        answer.at(2, 1) = stress.at(3);
        answer.at(1, 3) = stress.at(7);
        answer.at(3, 1) = stress.at(7);
        answer.at(2, 3) = stress.at(8);
        answer.at(3, 2) = stress.at(8);
    } else if ( ( type == LocalMomentumTensor ) || ( type == GlobalMomentumTensor ) ) {
        FloatArray stress;
        this->computeStressVector(stress, gp, tStep);
        answer.at(1, 1) = stress.at(4);
        answer.at(2, 2) = stress.at(5);
        answer.at(1, 2) = stress.at(6);
        answer.at(2, 1) = stress.at(6);
    } else if ( ( type == LocalStrainTensor ) || ( type == GlobalStrainTensor ) ) {
        FloatArray strain;
        this->computeStrainVector(strain, gp, tStep);

        answer.at(1, 1) = strain.at(1);
        answer.at(2, 2) = strain.at(2);
        answer.at(1, 2) = strain.at(3) / 2.;
        answer.at(2, 1) = strain.at(3) / 2.;
        answer.at(1, 3) = strain.at(7) / 2.;
        answer.at(3, 1) = strain.at(7) / 2.;
        answer.at(2, 3) = strain.at(8) / 2.;
        answer.at(3, 2) = strain.at(8) / 2.;
    } else if ( ( type == LocalCurvatureTensor ) || ( type == GlobalCurvatureTensor ) ) {
        FloatArray curv;
        this->computeStrainVector(curv, gp, tStep);
        answer.at(1, 1) = curv.at(4);
        answer.at(2, 2) = curv.at(5);
        answer.at(1, 2) = curv.at(6) / 2.;
        answer.at(2, 1) = curv.at(6) / 2.;
    } else {
        _error("GiveCharacteristicTensor: unsupported tensor mode");
        exit(1);
    }

    if ( ( type == GlobalForceTensor ) || ( type == GlobalMomentumTensor ) ||
        ( type == GlobalStrainTensor ) || ( type == GlobalCurvatureTensor ) ) {
        this->computeGtoLRotationMatrix();
        answer.rotatedWith(* GtoLRotationMatrix);
    }
}


void
RerShell :: computeStrainVectorInLayer(FloatArray &answer, GaussPoint *masterGp,
                                       GaussPoint *slaveGp, TimeStep *tStep)
//
// returns full 3d strain vector of given layer (whose z-coordinate from center-line is
// stored in slaveGp) for given tStep
//
{
    FloatArray masterGpStrain;
    double layerZeta, layerZCoord, top, bottom;

    this->computeStrainVector(masterGpStrain, masterGp, tStep);
    top    = masterGp->giveElement()->giveCrossSection()->give(CS_TopZCoord);
    bottom = masterGp->giveElement()->giveCrossSection()->give(CS_BottomZCoord);
    layerZeta = slaveGp->giveCoordinate(3);
    layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );

    answer.resize(6); // {Exx,Eyy,Ezz,GMyz,GMzx,GMxy}
    answer.zero();

    answer.at(1) = masterGpStrain.at(1) + masterGpStrain.at(4) * layerZCoord;
    answer.at(2) = masterGpStrain.at(2) + masterGpStrain.at(5) * layerZCoord;
    answer.at(6) = masterGpStrain.at(3) + masterGpStrain.at(6) * layerZCoord;
    answer.at(4) = masterGpStrain.at(8);
    answer.at(5) = masterGpStrain.at(7);
}


void
RerShell :: printOutputAt(FILE *file, TimeStep *stepN)
// Performs end-of-step operations.
{
    int i;
    GaussPoint *gp;
    FloatMatrix globTensorMembrane, globTensorPlate;

    fprintf(file, "element %d :\n", number);

    for ( i = 1; i <= numberOfGaussPoints; i++ ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(i - 1);
        if ( !domain->giveEngngModel()->isIncremental() ) {
            // delete this -> ComputeStrainVector(gp,stepN) ;
            // delete this -> ComputeStressVector(gp,stepN) ;
        }

        //gp   -> printOutputAt(file,stepN) ;


        fprintf( file, "  GP %d :", gp->giveNumber() );
        this->giveCharacteristicTensor(globTensorMembrane, GlobalStrainTensor, gp, stepN);
        this->giveCharacteristicTensor(globTensorPlate, GlobalCurvatureTensor, gp, stepN);
        fprintf(file, "  strains ");
        fprintf( file, " % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e ",
                globTensorMembrane.at(1, 1), globTensorMembrane.at(2, 2), globTensorMembrane.at(3, 3),
                2. * globTensorMembrane.at(2, 3), 2. * globTensorMembrane.at(3, 1), 2. * globTensorMembrane.at(1, 2),
                globTensorPlate.at(1, 1), globTensorPlate.at(2, 2), globTensorPlate.at(3, 3),
                2. * globTensorPlate.at(2, 3), 2. * globTensorPlate.at(1, 3), 2. * globTensorPlate.at(1, 2) );

        this->giveCharacteristicTensor(globTensorMembrane, GlobalForceTensor, gp, stepN);
        this->giveCharacteristicTensor(globTensorPlate, GlobalMomentumTensor, gp, stepN);
        fprintf(file, "\n          stresses");
        fprintf( file, " % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e % .4e ",
                globTensorMembrane.at(1, 1), globTensorMembrane.at(2, 2), globTensorMembrane.at(3, 3),
                globTensorMembrane.at(2, 3), globTensorMembrane.at(3, 1), globTensorMembrane.at(1, 2),
                globTensorPlate.at(1, 1), globTensorPlate.at(2, 2), globTensorPlate.at(3, 3),
                globTensorPlate.at(2, 3), globTensorPlate.at(1, 3), globTensorPlate.at(1, 2) );

        fprintf(file, "\n");
    }
}


void
RerShell :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    answer.setValues(6, D_u, D_v, D_w, R_u, R_v, R_w);
}


int
RerShell :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( type == IST_ShellForceMomentumTensor ) {
        return 12;
    }

    return 0;
}


void
RerShell :: ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatMatrix &answer, GaussPoint *aGaussPoint, InternalStateType type)
{
    // evaluates N matrix (interpolation estimated stress matrix)
    // according to Zienkiewicz & Zhu paper
    // N(nsigma, nsigma*nnodes)
    // Definition : sigmaVector = N * nodalSigmaVector
    double l1, l2, l3;

    l1 = aGaussPoint->giveCoordinate(1);
    l2 = aGaussPoint->giveCoordinate(2);
    l3 = 1.0 - l1 - l2;

    if ( type == IST_ShellForceMomentumTensor ) {
        answer.resize(1, 3);
    } else {
        return;
    }

    answer.zero();
    answer.at(1, 1) = l1;
    answer.at(1, 2) = l2;
    answer.at(1, 3) = l3;
}


void
RerShell :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                       InternalStateType type, TimeStep *tStep)
{
    this->giveIPValue(answer, integrationRulesArray [ 0 ]->getIntegrationPoint(0), type, tStep);
}


void
RerShell :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                      InternalStateType type, TimeStep *tStep)
{
    answer.resize(0);
}


int
RerShell :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type)
{
    if ( type == IST_ShellForceMomentumTensor ) {
        answer.resize(12);
        for ( int i = 1; i <= 12; i++ ) {
            answer.at(i) = i;
        }

        return 1;
    } else {
        return CCTPlate :: giveIntVarCompFullIndx(answer, type);
    }
}


int
RerShell :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
    FloatMatrix globTensor;
    CharTensor cht;

    answer.resize(12);

    if ( ( type == IST_ShellForceMomentumTensor ) || ( type == IST_ShellStrainCurvatureTensor ) ) {
        if ( type == IST_ShellForceMomentumTensor ) {
            cht = GlobalForceTensor;
        } else {
            cht = GlobalStrainTensor;
        }

        this->giveCharacteristicTensor(globTensor, cht, aGaussPoint, tStep);

        answer.at(1) = globTensor.at(1, 1);  //sxForce
        answer.at(2) = globTensor.at(2, 2);  //syForce
        answer.at(3) = globTensor.at(3, 3);  //szForce
        answer.at(4) = globTensor.at(2, 3);  //syzForce
        answer.at(5) = globTensor.at(1, 3);  //qxzForce
        answer.at(6) = globTensor.at(1, 2);  //qxyForce

        if ( type == IST_ShellForceMomentumTensor ) {
            cht = GlobalMomentumTensor;
        } else {
            cht = GlobalCurvatureTensor;
        }


        this->giveCharacteristicTensor(globTensor, cht, aGaussPoint, tStep);

        answer.at(7)  = globTensor.at(1, 1);  //mxForce
        answer.at(8)  = globTensor.at(2, 2);  //myForce
        answer.at(9)  = globTensor.at(3, 3);  //mzForce
        answer.at(10) = globTensor.at(2, 3);  //myzForce
        answer.at(11) = globTensor.at(1, 3);  //mxzForce
        answer.at(12) = globTensor.at(1, 2);  //mxyForce

        return 1;
    } else {
        answer.resize(0);
        return 0;
    }
}


#ifdef __OOFEG
/*
 * void
 * CCTPlate  :: drawRawGeometry ()
 * {
 * WCRec p[3];
 * GraphicObj *go;
 *
 * EASValsSetLineWidth(RAW_GEOMETRY_WIDTH);
 * EASValsSetColor(elementColor);
 * EASValsSetLayer(RAW_GEOMETRY_LAYER);
 * p[0].x = (FPNum) this->giveNode(1)->giveCoordinate(1);
 * p[0].y = (FPNum) this->giveNode(1)->giveCoordinate(2);
 * p[0].z = (FPNum) this->giveNode(1)->giveCoordinate(3);
 * p[1].x = (FPNum) this->giveNode(2)->giveCoordinate(1);
 * p[1].y = (FPNum) this->giveNode(2)->giveCoordinate(2);
 * p[1].z = (FPNum) this->giveNode(2)->giveCoordinate(3);
 * p[2].x = (FPNum) this->giveNode(3)->giveCoordinate(1);
 * p[2].y = (FPNum) this->giveNode(3)->giveCoordinate(2);
 * p[2].z = (FPNum) this->giveNode(3)->giveCoordinate(3);
 *
 * go =  CreateTriangle3D(p);
 * EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
 * EGAttachObject(go, (EObjectP) this);
 * EMAddGraphicsToModel(ESIModel(), go);
 * }
 *
 * void
 * CCTPlate  :: drawDeformedGeometry (UnknownType defMode)
 * {
 * WCRec p[3];
 * GraphicObj *go;
 * TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
 *
 * EASValsSetLineWidth(DEFORMED_GEOMETRY_WIDTH);
 * EASValsSetColor(deformedElementColor);
 * EASValsSetLayer(DEFORMED_GEOMETRY_LAYER);
 * p[0].x = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(1,tStep,defMode,defScale);
 * p[0].y = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(2,tStep,defMode,defScale);
 * p[0].z = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(3,tStep,defMode,defScale);
 * p[1].x = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(1,tStep,defMode,defScale);
 * p[1].y = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(2,tStep,defMode,defScale);
 * p[1].z = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(3,tStep,defMode,defScale);
 * p[2].x = (FPNum) this->giveNode(3)->giveUpdatedCoordinate(1,tStep,defMode,defScale);
 * p[2].y = (FPNum) this->giveNode(3)->giveUpdatedCoordinate(2,tStep,defMode,defScale);
 * p[2].z = (FPNum) this->giveNode(3)->giveUpdatedCoordinate(3,tStep,defMode,defScale);
 *
 * go =  CreateTriangle3D(p);
 * EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
 * EMAddGraphicsToModel(ESIModel(), go);
 * }
 *
 */


//void
//RerShell  :: drawScalar   (oofegGraphicContext& context)
//{}


/*
 * void
 * RerShell  :: drawInternalState (oofegGraphicContext& gc)
 * //
 * // Draws internal state graphics representation
 * //
 * {
 * int i;
 * WCRec p[3];
 * GraphicObj *tr;
 * double v1,v2,v3;
 * DrawMode mode = gc.getDrawMode();
 * TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
 * double defScale = gc.getDefScale();
 *
 * if (!gc.testElementGraphicActivity(this)) return;
 *
 * // check for valid DrawMode
 * if (!((mode == mxForce) || (mode == myForce) || (mode == mxyForce) ||
 * (mode == szxForce) || (mode == syzForce) || (mode == sxForce) ||
 * (mode == syForce) || (mode == sxyForce))) return;
 *
 * EASValsSetLayer(OOFEG_STRESS_CONTOUR_LAYER);
 * for (i=0; i< 3; i++) {
 * if (gc.getInternalVarsDefGeoFlag()) {
 *  // use deformed geometry
 *  p[i].x = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(1,tStep,DisplacementVector,defScale);
 *  p[i].y = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(2,tStep,DisplacementVector,defScale);
 *  p[i].z = (FPNum) this->giveNode(i+1)->giveUpdatedCoordinate(3,tStep,DisplacementVector,defScale);
 *
 * } else {
 *  p[i].x = (FPNum) this->giveNode(i+1)->giveCoordinate(1);
 *  p[i].y = (FPNum) this->giveNode(i+1)->giveCoordinate(2);
 *  p[i].z = (FPNum) this->giveNode(i+1)->giveCoordinate(3);
 * }
 * }
 *
 * int result = 0;
 * result+= this->giveInternalStateAtNode (gc, 1, &v1);
 * result+= this->giveInternalStateAtNode (gc, 2, &v2);
 * result+= this->giveInternalStateAtNode (gc, 3, &v3);
 *
 * if (result == 3) {
 *
 * tr = CreateTriangleWD3D (p,v1,v2,v3);
 * EGWithMaskChangeAttributes(LAYER_MASK, tr);
 * EMAddGraphicsToModel(ESIModel(), tr);
 * }
 * }
 */

/*
 * int
 * RerShell :: giveInternalStateAtNode (FloatArray& answer, InternalValueRequest& r, int node, TimeStep* atTime)
 *
 * //
 * // see eleemnt.h for description
 * //
 * {
 * opravit;
 *
 * DrawMode mode = gc.getDrawMode();
 *
 * if (!gc.testElementGraphicActivity(this)) return 0;
 *
 * const FloatArray* nodval;
 * int result = this->giveDomain()->giveSmoother()->giveNodalVector(nodval, this->giveNode(inode)->giveNumber(),
 *                       this->giveDomain()->giveSmoother()->giveElementRegion(this));
 * if (result) {
 * if (mode == sxForce ) {
 **val =  nodval->at(1);
 * return 1;
 * } else if (mode == syForce) {
 **val =  nodval->at(2);
 * return 1;
 * } else if (mode == sxyForce) {
 **val =  nodval->at(3);
 * return 1;
 * } else if (mode == mxForce ) {
 **val =  nodval->at(4);
 * return 1;
 * } else if (mode == myForce) {
 **val =  nodval->at(5);
 * return 1;
 * } else if (mode == mxyForce) {
 **val =  nodval->at(6);
 * return 1;
 * } else if (mode == szxForce ) {
 **val =  nodval->at(7);
 * return 1;
 * } else if (mode == syzForce) {
 **val =  nodval->at(8);
 * return 1;
 * } else return 0;
 * }
 * return 0;
 * }
 */
#endif
} // end namespace oofem
