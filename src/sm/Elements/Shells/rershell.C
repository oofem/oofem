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

#include "../sm/Elements/Shells/rershell.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "load.h"
#include "mathfem.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif

#include <cstdlib>

namespace oofem {
REGISTER_Element(RerShell);

RerShell :: RerShell(int n, Domain *aDomain) :
    CCTPlate3d(n, aDomain)
{
    Rx = 1.e+40;
    Ry = 1.e+40;
    Rxy = 1.e+40;

    numberOfGaussPoints = 1;
}


Interface *
RerShell :: giveInterface(InterfaceType interface)
{
    if ( interface == LayeredCrossSectionInterfaceType ) {
        return static_cast< LayeredCrossSectionInterface * >(this);
    } else if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    }

    return NULL;
}


void
RerShell :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [8x18] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
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

    for ( int i = 1; i <= 3; i++ ) {
        int j = i + 1 - i / 3 * 3;
        int k = j + 1 - j / 3 * 3;

        int ii = 6 * ( i - 1 );
        int ki = 6 * ( k - 1 );

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
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 8) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}


void
RerShell :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver,
// evaluated at gp.
{
    double x1, x2, x3, y1, y2, y3, l1, l2, l3, b1, b2, b3, c1, c2, c3;
    FloatArray nodeCoords;

    l1 = iLocCoord.at(1);
    l2 = iLocCoord.at(2);
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
    numberOfGaussPoints = 1;
    IRResultType result = StructuralElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    if ( numberOfGaussPoints != 1 ) {
        numberOfGaussPoints = 1;
    }

    return IRRT_OK;
}


void
RerShell :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->give3dShellStiffMtrx(answer, rMode, gp, tStep);
}

void
RerShell :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveGeneralizedStress_Shell(answer, gp, strain, tStep);
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
    mss1 = dV * this->giveCrossSection()->give(CS_Thickness, gp) * this->giveStructuralCrossSection()->give('d', gp) / 3.;

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
RerShell :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode)
// Computes numerically the load vector of the receiver due to the body
// loads, at tStep. - no support for momentum bodyload
{
    double dens, dV, load;
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    FloatArray f;
    FloatMatrix T;


    forLoad->computeComponentArrayAt(f, tStep, mode);

    if ( f.giveSize() == 0 ) {
        answer.clear();
        return;                                             // nil resultant
    } else {
        dens = this->giveStructuralCrossSection()->give('d', gp);
        dV = this->computeVolumeAround(gp) * this->giveCrossSection()->give(CS_Thickness, gp);

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

int
RerShell :: giveLocalCoordinateSystem(FloatMatrix &answer)
//
// returns a unit vectors of local coordinate system at element
// stored rowwise (mainly used by some materials with ortho and anisotrophy)
//
{
    this->computeGtoLRotationMatrix();
    answer = GtoLRotationMatrix;
    return 1;
}


//converts global coordinates to local planar area coordinates, does not return a coordinate in the thickness direction, but
//does check that the point is in the element thickness
#define POINT_TOL 1.e-3
bool
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
    StructuralCrossSection *cs = this->giveStructuralCrossSection();
    GaussPoint _gp(NULL, 1, answer, 1.0, _2dPlate);

    double elthick;

    elthick = cs->give(CS_Thickness, & _gp);

    if ( elthick / 2.0 + midplZ - fabs( inputCoords_ElCS.at(3) ) < -POINT_TOL ) {
        answer.zero();
        return false;
    }

    //check that the point is in the element and set flag
    for ( int i = 1; i <= 3; i++ ) {
        if ( answer.at(i) < ( 0. - POINT_TOL ) ) {
            return false;
        }

        if ( answer.at(i) > ( 1. + POINT_TOL ) ) {
            return false;
        }
    }

    return true;
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

    this->computeGtoLRotationMatrix();


    answer.resize(18, 18);
    answer.zero();

    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            val = GtoLRotationMatrix.at(i, j);
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
        OOFEM_ERROR("cannot transform coordinates- size mismatch");
        exit(1);
    }

    // first ensure that receiver's GtoLRotationMatrix[3,3]
    // is defined

    this->computeGtoLRotationMatrix();
    answer.beProductOf(GtoLRotationMatrix, global);
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
        FloatArray stress, strain;
        this->computeStrainVector(strain, gp, tStep);
        this->computeStressVector(stress, strain, gp, tStep);
        answer.at(1, 1) = stress.at(1);
        answer.at(2, 2) = stress.at(2);
        answer.at(1, 2) = stress.at(3);
        answer.at(2, 1) = stress.at(3);
        answer.at(1, 3) = stress.at(7);
        answer.at(3, 1) = stress.at(7);
        answer.at(2, 3) = stress.at(8);
        answer.at(3, 2) = stress.at(8);
    } else if ( ( type == LocalMomentumTensor ) || ( type == GlobalMomentumTensor ) ) {
        FloatArray stress, strain;
        this->computeStrainVector(strain, gp, tStep);
        this->computeStressVector(stress, strain, gp, tStep);
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
        OOFEM_ERROR("unsupported tensor mode");
        exit(1);
    }

    if ( ( type == GlobalForceTensor ) || ( type == GlobalMomentumTensor ) ||
        ( type == GlobalStrainTensor ) || ( type == GlobalCurvatureTensor ) ) {
        this->computeGtoLRotationMatrix();
        answer.rotatedWith(GtoLRotationMatrix);
    }
}


void
RerShell :: computeStrainVectorInLayer(FloatArray &answer, const FloatArray &masterGpStrain,
                                       GaussPoint *masterGp, GaussPoint *slaveGp, TimeStep *tStep)
//
// returns full 3d strain vector of given layer (whose z-coordinate from center-line is
// stored in slaveGp) for given tStep
//
{
    double layerZeta, layerZCoord, top, bottom;

    top    = this->giveCrossSection()->give(CS_TopZCoord, masterGp);
    bottom = this->giveCrossSection()->give(CS_BottomZCoord, masterGp);
    layerZeta = slaveGp->giveNaturalCoordinate(3);
    layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );

    answer.resize(5); // {Exx,Eyy,GMyz,GMzx,GMxy}

    answer.at(1) = masterGpStrain.at(1) + masterGpStrain.at(4) * layerZCoord;
    answer.at(2) = masterGpStrain.at(2) + masterGpStrain.at(5) * layerZCoord;
    answer.at(5) = masterGpStrain.at(3) + masterGpStrain.at(6) * layerZCoord;
    answer.at(3) = masterGpStrain.at(8);
    answer.at(4) = masterGpStrain.at(7);
}


void
RerShell :: printOutputAt(FILE *file, TimeStep *tStep)
{
    FloatArray v;

    fprintf(file, "element %d (%8d):\n", this->giveLabel(), number);

    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {

        fprintf( file, "  GP 1.%d :", gp->giveNumber() );
        this->giveIPValue(v, gp, IST_ShellStrainTensor, tStep);
        fprintf(file, "  strains    ");
        // eps_x, eps_y, eps_z, eps_yz, eps_xz, eps_xy (global)
        fprintf( file,
                " %.4e %.4e %.4e %.4e %.4e %.4e ",
                v.at(1), v.at(2), v.at(3), v.at(4), v.at(5), v.at(6) );

        this->giveIPValue(v, gp, IST_ShellCurvatureTensor, tStep);
        fprintf(file, "\n              curvatures ");
        // k_x, k_y, k_z, k_yz, k_xz, k_xy (global)
        fprintf( file,
                " %.4e %.4e %.4e %.4e %.4e %.4e ",
                v.at(1), v.at(2), v.at(3), v.at(4), v.at(5), v.at(6) );

        // Forces - Moments
        this->giveIPValue(v, gp, IST_ShellForceTensor, tStep);
        fprintf(file, "\n              stresses   ");
        // n_x, n_y, n_z, v_yz, v_xz, v_xy (global)
        fprintf( file,
                " %.4e %.4e %.4e %.4e %.4e %.4e ",
                v.at(1), v.at(2), v.at(3), v.at(4), v.at(5), v.at(6) );

        this->giveIPValue(v, gp, IST_ShellMomentumTensor, tStep);
        fprintf(file, "\n              moments    ");
        // m_x, m_y, m_z, m_yz, m_xz, m_xy (global)
        fprintf( file,
                " %.4e %.4e %.4e %.4e %.4e %.4e ",
                v.at(1), v.at(2), v.at(3), v.at(4), v.at(5), v.at(6) );

        fprintf(file, "\n");
    }
}


void
RerShell :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_u, D_v, D_w, R_u, R_v, R_w};
}


void
RerShell :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                       InternalStateType type, TimeStep *tStep)
{
    this->giveIPValue(answer, integrationRulesArray [ 0 ]->getIntegrationPoint(0), type, tStep);
}


int
RerShell :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    FloatMatrix globTensor;
    CharTensor cht;

    answer.resize(6);

    if (  type == IST_ShellCurvatureTensor || type == IST_ShellStrainTensor ) {
        if ( type == IST_ShellCurvatureTensor ) {
            cht = GlobalCurvatureTensor;
        } else {
            cht = GlobalStrainTensor;
        }

        this->giveCharacteristicTensor(globTensor, cht, gp, tStep);

        answer.at(1) = globTensor.at(1, 1); //xx
        answer.at(2) = globTensor.at(2, 2); //yy
        answer.at(3) = globTensor.at(3, 3); //zz
        answer.at(4) = 2*globTensor.at(2, 3); //yz
        answer.at(5) = 2*globTensor.at(1, 3); //xz
        answer.at(6) = 2*globTensor.at(2, 3); //yz

        return 1;
    } else if ( type == IST_ShellMomentumTensor || type == IST_ShellForceTensor ) {
        if ( type == IST_ShellMomentumTensor ) {
            cht = GlobalMomentumTensor;
        } else {
            cht = GlobalForceTensor;
        }

        this->giveCharacteristicTensor(globTensor, cht, gp, tStep);

        answer.at(1) = globTensor.at(1, 1); //xx
        answer.at(2) = globTensor.at(2, 2); //yy
        answer.at(3) = globTensor.at(3, 3); //zz
        answer.at(4) = globTensor.at(2, 3); //yz
        answer.at(5) = globTensor.at(1, 3); //xz
        answer.at(6) = globTensor.at(2, 3); //yz

        return 1;
    } else {
        return StructuralElement :: giveIPValue(answer, gp, type, tStep);
    }
}


#ifdef __OOFEG
/*
 * void
 * CCTPlate :: drawRawGeometry ()
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
 * CCTPlate :: drawDeformedGeometry (UnknownType defMode)
 * {
 * WCRec p[3];
 * GraphicObj *go;
 *
 * EASValsSetLineWidth(DEFORMED_GEOMETRY_WIDTH);
 * EASValsSetColor(deformedElementColor);
 * EASValsSetLayer(DEFORMED_GEOMETRY_LAYER);
 * p[0].x = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(1,tStep,defScale);
 * p[0].y = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(2,tStep,defScale);
 * p[0].z = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(3,tStep,defScale);
 * p[1].x = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(1,tStep,defScale);
 * p[1].y = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(2,tStep,defScale);
 * p[1].z = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(3,tStep,defScale);
 * p[2].x = (FPNum) this->giveNode(3)->giveUpdatedCoordinate(1,tStep,defScale);
 * p[2].y = (FPNum) this->giveNode(3)->giveUpdatedCoordinate(2,tStep,defScale);
 * p[2].z = (FPNum) this->giveNode(3)->giveUpdatedCoordinate(3,tStep,defScale);
 *
 * go =  CreateTriangle3D(p);
 * EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
 * EMAddGraphicsToModel(ESIModel(), go);
 * }
 *
 */


//void
//RerShell :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
//{}



/*
 * int
 * RerShell :: giveInternalStateAtNode (FloatArray& answer, InternalValueRequest& r, int node, TimeStep* tStep)
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
 ***val =  nodval->at(1);
 * return 1;
 * } else if (mode == syForce) {
 ***val =  nodval->at(2);
 * return 1;
 * } else if (mode == sxyForce) {
 ***val =  nodval->at(3);
 * return 1;
 * } else if (mode == mxForce ) {
 ***val =  nodval->at(4);
 * return 1;
 * } else if (mode == myForce) {
 ***val =  nodval->at(5);
 * return 1;
 * } else if (mode == mxyForce) {
 ***val =  nodval->at(6);
 * return 1;
 * } else if (mode == szxForce ) {
 ***val =  nodval->at(7);
 * return 1;
 * } else if (mode == syzForce) {
 ***val =  nodval->at(8);
 * return 1;
 * } else return 0;
 * }
 * return 0;
 * }
 */
#endif
} // end namespace oofem
