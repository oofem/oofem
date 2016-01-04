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

#include "../sm/Elements/Plates/dkt.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "fei2dtrlin.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "load.h"
#include "mathfem.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(DKTPlate);

FEI2dTrLin DKTPlate :: interp_lin(1, 2);

DKTPlate :: DKTPlate(int n, Domain *aDomain) :
    NLStructuralElement(n, aDomain),
    LayeredCrossSectionInterface(), ZZNodalRecoveryModelInterface(this),
    NodalAveragingRecoveryModelInterface(), SPRNodalRecoveryModelInterface(), ZZErrorEstimatorInterface(this)
{
    numberOfDofMans = 3;
    numberOfGaussPoints = 3;
    area = 0;
}


FEInterpolation *
DKTPlate :: giveInterpolation() const { return & interp_lin; }


FEInterpolation *
DKTPlate :: giveInterpolation(DofIDItem id) const
{
    return & interp_lin;
}


void
DKTPlate :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 5) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}


void
DKTPlate :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode)
// Computes numerically the load vector of the receiver due to the body loads, at tStep.
// load is assumed to be in global cs.
// load vector is then transformed to coordinate system in each node.
// (should be global coordinate system, but there may be defined
//  different coordinate system in each node)
{
    FloatArray force;
    FloatMatrix T;

    if ( ( forLoad->giveBCGeoType() != BodyLoadBGT ) || ( forLoad->giveBCValType() != ForceLoadBVT ) ) {
        OOFEM_ERROR("unknown load type");
    }

    GaussIntegrationRule irule(1, this, 1, 5);
    irule.SetUpPointsOnTriangle(1, _2dPlate);

    // note: force is assumed to be in global coordinate system.
    forLoad->computeComponentArrayAt(force, tStep, mode);

    if ( force.giveSize() ) {
        GaussPoint *gp = irule.getIntegrationPoint(0);
        double dens = this->giveStructuralCrossSection()->give('d', gp);
        double dV   = this->computeVolumeAround(gp) * this->giveCrossSection()->give(CS_Thickness, gp);

        answer.resize(9);
        answer.zero();

        double load = force.at(1) * dens * dV / 3.0;
        answer.at(1) = load;
        answer.at(4) = load;
        answer.at(7) = load;

        // transform result from global cs to local element cs.
        if ( this->computeGtoLRotationMatrix(T) ) {
            answer.rotatedWith(T, 'n');
        }
    } else {
        answer.clear();          // nil resultant
    }
}


void
DKTPlate :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [5x9] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
    // get node coordinates
    double x1, x2, x3, y1, y2, y3, z1, z2, z3;
    this->giveNodeCoordinates(x1, x2, x3, y1, y2, y3, z1, z2, z3);


    double ksi = gp->giveNaturalCoordinate(1);
    double eta = gp->giveNaturalCoordinate(2);

    double N1dk = 4.0 * ksi - 1.0;
    double N2dk = 0.0;
    double N3dk = -3.0 + 4.0 * ksi + 4.0 * eta;
    double N4dk = 4.0 * eta;
    double N5dk = -4.0 * eta;
    double N6dk = 4.0 * ( 1.0 - 2.0 * ksi - eta );

    double N1de = 0.0;
    double N2de = 4.0 * eta - 1.0;
    double N3de = -3.0 + 4.0 * eta + 4.0 * ksi;
    double N4de = 4.0 * ksi;
    double N5de = 4.0 * ( 1.0 - 2.0 * eta - ksi );
    double N6de = -4.0 * ksi;

    double A = N1dk * x1 + N2dk * x2 + N3dk * x3 + N4dk * ( x1 + x2 ) / 2 + N5dk * ( x2 + x3 ) / 2 + N6dk * ( x1 + x3 ) / 2;
    double B = N1dk * y1 + N2dk * y2 + N3dk * y3 + N4dk * ( y1 + y2 ) / 2 + N5dk * ( y2 + y3 ) / 2 + N6dk * ( y1 + y3 ) / 2;
    double C = N1de * x1 + N2de * x2 + N3de * x3 + N4de * ( x1 + x2 ) / 2 + N5de * ( x2 + x3 ) / 2 + N6de * ( x1 + x3 ) / 2;
    double D = N1de * y1 + N2de * y2 + N3de * y3 + N4de * ( y1 + y2 ) / 2 + N5de * ( y2 + y3 ) / 2 + N6de * ( y1 + y3 ) / 2;

    double dxdk = D / ( A * D - B * C );
    double dydk = -C / ( A * D - B * C );
    double dxde = -B / ( A * D - B * C );
    double dyde = A / ( A * D - B * C );

    double dN102 = N1dk * dxdk + N1de * dxde;
    double dN104 = N2dk * dxdk + N2de * dxde;
    double dN106 = N3dk * dxdk + N3de * dxde;
    double dN108 = N4dk * dxdk + N4de * dxde;
    double dN110 = N5dk * dxdk + N5de * dxde;
    double dN112 = N6dk * dxdk + N6de * dxde;

    double dN201 = -N1dk * dydk - N1de * dyde;
    double dN203 = -N2dk * dydk - N2de * dyde;
    double dN205 = -N3dk * dydk - N3de * dyde;
    double dN207 = -N4dk * dydk - N4de * dyde;
    double dN209 = -N5dk * dydk - N5de * dyde;
    double dN211 = -N6dk * dydk - N6de * dyde;

    // normals on element sides
    double dx4 = x2 - x1;
    double dy4 = y2 - y1;
    double l4 = sqrt(dx4 * dx4 + dy4 * dy4);

    double dx5 = x3 - x2;
    double dy5 = y3 - y2;
    double l5 = sqrt(dx5 * dx5 + dy5 * dy5);

    double dx6 = x1 - x3;
    double dy6 = y1 - y3;
    double l6 = sqrt(dx6 * dx6 + dy6 * dy6);

    double c4 = dy4 / l4;
    double s4 = -dx4 / l4;

    double c5 = dy5 / l5;
    double s5 = -dx5 / l5;

    double c6 = dy6 / l6;
    double s6 = -dx6 / l6;

    // transformation matrix between element DOFs (w, fi_x, fi_y)  and DOFs of element with qudratic rotations
    double T11 = -1.5 / l4 * c4;
    double T12 = -0.25 * c4 * c4 + 0.5 * s4 * s4;
    double T13 = -0.25 * c4 * s4 - 0.5 * c4 * s4;
    double T14 = 1.5 / l4 * c4;
    double T15 = -0.25 * c4 * c4 + 0.5 * s4 * s4;
    double T16 = -0.25 * c4 * s4 - 0.5 * c4 * s4;

    double T21 = -1.5 / l4 * s4;
    double T22 = -0.25 * c4 * s4 - 0.5 * c4 * s4;
    double T23 = -0.25 * s4 * s4 + 0.5 * c4 * c4;
    double T24 = 1.5 / l4 * s4;
    double T25 = -0.25 * c4 * s4 - 0.5 * c4 * s4;
    double T26 = -0.25 * s4 * s4 + 0.5 * c4 * c4;

    double T34 = -1.5 / l5 * c5;
    double T35 = -0.25 * c5 * c5 + 0.5 * s5 * s5;
    double T36 = -0.25 * c5 * s5 - 0.5 * c5 * s5;
    double T37 = 1.5 / l5 * c5;
    double T38 = -0.25 * c5 * c5 + 0.5 * s5 * s5;
    double T39 = -0.25 * c5 * s5 - 0.5 * c5 * s5;

    double T44 = -1.5 / l5 * s5;
    double T45 = -0.25 * c5 * s5 - 0.5 * c5 * s5;
    double T46 = -0.25 * s5 * s5 + 0.5 * c5 * c5;
    double T47 = 1.5 / l5 * s5;
    double T48 = -0.25 * c5 * s5 - 0.5 * c5 * s5;
    double T49 = -0.25 * s5 * s5 + 0.5 * c5 * c5;

    double T51 = 1.5 / l6 * c6;
    double T52 = -0.25 * c6 * c6 + 0.5 * s6 * s6;
    double T53 = -0.25 * c6 * s6 - 0.5 * c6 * s6;
    double T57 = -1.5 / l6 * c6;
    double T58 = -0.25 * c6 * c6 + 0.5 * s6 * s6;
    double T59 = -0.25 * c6 * s6 - 0.5 * c6 * s6;

    double T61 = 1.5 / l6 * s6;
    double T62 = -0.25 * c6 * s6 - 0.5 * c6 * s6;
    double T63 = -0.25 * s6 * s6 + 0.5 * c6 * c6;
    double T67 = -1.5 / l6 * s6;
    double T68 = -0.25 * c6 * s6 - 0.5 * c6 * s6;
    double T69 = -0.25 * s6 * s6 + 0.5 * c6 * c6;

    answer.resize(5, 9);
    answer.zero();

    answer.at(1, 1) = T21 * dN108 + T61 * dN112;
    answer.at(1, 2) = T22 * dN108 + T62 * dN112;
    answer.at(1, 3) = dN102 + T23 * dN108 + T63 * dN112;
    answer.at(1, 4) = T24 * dN108 + T44 * dN110;
    answer.at(1, 5) = T25 * dN108 + T45 * dN110;
    answer.at(1, 6) = dN104 + T26 * dN108 + T46 * dN110;
    answer.at(1, 7) = T47 * dN110 + T67 * dN112;
    answer.at(1, 8) = T48 * dN110 + T68 * dN112;
    answer.at(1, 9) = dN106 + T49 * dN110 + T69 * dN112;

    answer.at(2, 1) = T11 * dN207 + T51 * dN211;
    answer.at(2, 2) = dN201 + T12 * dN207 + T52 * dN211;
    answer.at(2, 3) = T13 * dN207 + T53 * dN211;
    answer.at(2, 4) = T14 * dN207 + T34 * dN209;
    answer.at(2, 5) = dN203 + T15 * dN207 + T35 * dN209;
    answer.at(2, 6) = T16 * dN207 + T36 * dN209;
    answer.at(2, 7) = T37 * dN209 + T57 * dN211;
    answer.at(2, 8) = dN205 + T38 * dN209 + T58 * dN211;
    answer.at(2, 9) = T39 * dN209 + T59 * dN211;

    answer.at(3, 1) = -T11 * dN108 - T51 * dN112 - T21 * dN207 - T61 * dN211;
    answer.at(3, 2) = -dN102 - T12 * dN108 - T52 * dN112 - T22 * dN207 - T62 * dN211;
    answer.at(3, 3) = -dN201 - T13 * dN108 - T53 * dN112 - T23 * dN207 - T63 * dN211;
    answer.at(3, 4) = -T14 * dN108 - T34 * dN110 - T24 * dN207 - T44 * dN209;
    answer.at(3, 5) = -dN104 - T15 * dN108 - T35 * dN110 - T25 * dN207 - T45 * dN209;
    answer.at(3, 6) = -dN203 - T16 * dN108 - T36 * dN110 - T26 * dN207 - T46 * dN209;
    answer.at(3, 7) = -T37 * dN110 - T57 * dN112 - T47 * dN209 - T67 * dN211;
    answer.at(3, 8) = -dN106 - T38 * dN110 - T58 * dN112 - T48 * dN209 - T68 * dN211;
    answer.at(3, 9) = -dN205 - T39 * dN110 - T59 * dN112 - T49 * dN209 - T69 * dN211;

    // Note: no shear strains, no shear forces => the 4th and 5th rows are zero
}


void
DKTPlate :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the [3x9] displacement interpolation matrix {N} of the receiver,
// evaluated at gp.
// Note: this interpolation is not available, as the deflection is cubic along the edges,
//       but not define in the interior of the element
// Note: the interpolation of rotations is quadratic
// NOTE: linear interpolation returned instead
{
    FloatArray N;

    answer.resize(3, 9);
    answer.zero();
    giveInterpolation()->evalN( N, iLocCoord, FEIElementGeometryWrapper(this) );

    answer.beNMatrixOf(N, 3);
}


void
DKTPlate :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveGeneralizedStress_Plate(answer, gp, strain, tStep);
}


void
DKTPlate :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->give2dPlateStiffMtrx(answer, rMode, gp, tStep);
}


void
DKTPlate :: giveNodeCoordinates(double &x1, double &x2, double &x3,
                                double &y1, double &y2, double &y3,
                                double &z1, double &z2, double &z3)
{
    FloatArray *nc1, *nc2, *nc3;
    nc1 = this->giveNode(1)->giveCoordinates();
    nc2 = this->giveNode(2)->giveCoordinates();
    nc3 = this->giveNode(3)->giveCoordinates();

    x1 = nc1->at(1);
    x2 = nc2->at(1);
    x3 = nc3->at(1);

    y1 = nc1->at(2);
    y2 = nc2->at(2);
    y3 = nc3->at(2);

    z1 = nc1->at(3);
    z2 = nc2->at(3);
    z3 = nc3->at(3);

}


IRResultType
DKTPlate :: initializeFrom(InputRecord *ir)
{
    return NLStructuralElement :: initializeFrom(ir);
}


void
DKTPlate :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_w, R_u, R_v};
}


void
DKTPlate :: computeMidPlaneNormal(FloatArray &answer, const GaussPoint *gp)
// returns normal vector to midPlane in GaussPoinr gp of receiver
{
    FloatArray u, v;
    u.beDifferenceOf( * this->giveNode(2)->giveCoordinates(), * this->giveNode(1)->giveCoordinates() );
    v.beDifferenceOf( * this->giveNode(3)->giveCoordinates(), * this->giveNode(1)->giveCoordinates() );

    answer.beVectorProductOf(u, v);
    answer.normalize();
}


double
DKTPlate :: giveCharacteristicLength(const FloatArray &normalToCrackPlane)
//
// returns receiver's characteristic length for crack band models
// for a crack formed in the plane with normal normalToCrackPlane.
//
{
    return this->giveCharacteristicLengthForPlaneElements(normalToCrackPlane);
}


double
DKTPlate :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double detJ, weight;
    std :: vector< FloatArray > lc = {FloatArray(3), FloatArray(3), FloatArray(3)};
    this->giveNodeCoordinates(lc[0].at(1), lc[1].at(1), lc[2].at(1), 
                              lc[0].at(2), lc[1].at(2), lc[2].at(2), 
                              lc[0].at(3), lc[1].at(3), lc[2].at(3));

    weight = gp->giveWeight();
    detJ = fabs( this->interp_lin.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIVertexListGeometryWrapper(lc) ) );
    return detJ * weight; ///@todo What about thickness?
}


void
DKTPlate :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver.
{
    answer.resize(9, 9);
    answer.zero();

    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

    double dV = this->computeVolumeAround(gp);
    double density = this->giveStructuralCrossSection()->give('d', gp);
    double mss1 = dV * this->giveCrossSection()->give(CS_Thickness, gp) * density / 3.;

    answer.at(1, 1) = mss1;
    answer.at(4, 4) = mss1;
    answer.at(7, 7) = mss1;
}


Interface *
DKTPlate :: giveInterface(InterfaceType interface)
{
    if ( interface == LayeredCrossSectionInterfaceType ) {
        return static_cast< LayeredCrossSectionInterface * >(this);
    } else if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return static_cast< NodalAveragingRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    } else if ( interface == ZZErrorEstimatorInterfaceType ) {
        return static_cast< ZZErrorEstimatorInterface * >(this);
    }


    return NULL;
}


#define POINT_TOL 1.e-3

bool
DKTPlate :: computeLocalCoordinates(FloatArray &answer, const FloatArray &coords)
//converts global coordinates to local planar area coordinates,
//does not return a coordinate in the thickness direction, but
//does check that the point is in the element thickness
{
    // get node coordinates
    double x1, x2, x3, y1, y2, y3, z1, z2, z3;
    this->giveNodeCoordinates(x1, x2, x3, y1, y2, y3, z1, z2, z3);

    // Fetch local coordinates.
    bool ok = this->interp_lin.global2local( answer, coords, FEIElementGeometryWrapper(this) ) > 0;

    //check that the point is in the element and set flag
    for ( int i = 1; i <= 3; i++ ) {
        if ( answer.at(i) < ( 0. - POINT_TOL ) ) {
            return false;
        }

        if ( answer.at(i) > ( 1. + POINT_TOL ) ) {
            return false;
        }
    }

    //get midplane location at this point
    double midplZ;
    midplZ = z1 * answer.at(1) + z2 * answer.at(2) + z3 * answer.at(3);

    //check that the z is within the element
    StructuralCrossSection *cs = this->giveStructuralCrossSection();
    double elthick = cs->give(CS_Thickness, answer, this);

    if ( elthick / 2.0 + midplZ - fabs( coords.at(3) ) < -POINT_TOL ) {
        answer.zero();
        return false;
    }


    return ok;
}


int
DKTPlate :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    FloatArray help;
    answer.resize(6);
    if ( type == IST_ShellForceTensor || type == IST_ShellStrainTensor ) {
        if ( type == IST_ShellForceTensor ) {
            help = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStressVector();
        } else {
            help = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStrainVector();
        }
        answer.at(1) = 0.0; // nx
        answer.at(2) = 0.0; // ny
        answer.at(3) = 0.0; // nz
        answer.at(4) = help.at(5); // vyz
        answer.at(5) = help.at(4); // vxz
        answer.at(6) = 0.0; // vxy
        return 1;
    } else if ( type == IST_ShellMomentumTensor || type == IST_ShellCurvatureTensor ) {
        if ( type == IST_ShellMomentumTensor ) {
            help = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStressVector();
        } else {
            help = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStrainVector();
        }
        answer.at(1) = help.at(1); // mx
        answer.at(2) = help.at(2); // my
        answer.at(3) = 0.0;        // mz
        answer.at(4) = 0.0;        // myz
        answer.at(5) = 0.0;        // mxz
        answer.at(6) = help.at(3); // mxy
        return 1;
    } else {
        return NLStructuralElement :: giveIPValue(answer, gp, type, tStep);
    }
}


//
// The element interface required by NodalAveragingRecoveryModel
//
void
DKTPlate :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                       InternalStateType type, TimeStep *tStep)
{
    GaussPoint *gp;
    if ( ( type == IST_ShellForceTensor ) || ( type == IST_ShellMomentumTensor ) ||
        ( type == IST_ShellStrainTensor )  || ( type == IST_ShellCurvatureTensor ) ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
        this->giveIPValue(answer, gp, type, tStep);
    } else {
        answer.clear();
    }
}


//
// The element interface required by SPRNodalRecoveryModelInterface
//
void
DKTPlate :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(3);
    pap.at(1) = this->giveNode(1)->giveNumber();
    pap.at(2) = this->giveNode(2)->giveNumber();
    pap.at(3) = this->giveNode(3)->giveNumber();
}


void
DKTPlate :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    answer.resize(1);
    if ( ( pap == this->giveNode(1)->giveNumber() ) ||
        ( pap == this->giveNode(2)->giveNumber() ) ||
        ( pap == this->giveNode(3)->giveNumber() ) ) {
        answer.at(1) = pap;
    } else {
        OOFEM_ERROR("node unknown");
    }
}


SPRPatchType
DKTPlate :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_2dxy;
}


//
// layered cross section support functions
//
void
DKTPlate :: computeStrainVectorInLayer(FloatArray &answer, const FloatArray &masterGpStrain,
                                       GaussPoint *masterGp, GaussPoint *slaveGp, TimeStep *tStep)
// returns full 3d strain vector of given layer (whose z-coordinate from center-line is
// stored in slaveGp) for given tStep
{
    double layerZeta, layerZCoord, top, bottom;

    top    = this->giveCrossSection()->give(CS_TopZCoord, masterGp);
    bottom = this->giveCrossSection()->give(CS_BottomZCoord, masterGp);
    layerZeta = slaveGp->giveNaturalCoordinate(3);
    layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );

    answer.resize(5); // {Exx,Eyy,GMyz,GMzx,GMxy}

    answer.at(1) = masterGpStrain.at(1) * layerZCoord;
    answer.at(2) = masterGpStrain.at(2) * layerZCoord;
    answer.at(5) = masterGpStrain.at(3) * layerZCoord;
    answer.at(3) = masterGpStrain.at(5);
    answer.at(4) = masterGpStrain.at(4);
}

// Edge load support
void
DKTPlate :: computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
{
    FloatArray n;

    this->interp_lin.edgeEvalN( n, iedge, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(3, 6);
    answer.at(1, 1) = n.at(1);
    answer.at(1, 4) = n.at(2);
    answer.at(2, 2) = answer.at(3, 3) = n.at(1);
    answer.at(2, 5) = answer.at(3, 6) = n.at(2);
}

void
DKTPlate :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    /*
     * provides dof mapping of local edge dofs (only nonzero are taken into account)
     * to global element dofs
     */

    answer.resize(6);
    if ( iEdge == 1 ) { // edge between nodes 1,2
        answer.at(1) = 1;
        answer.at(2) = 2;
        answer.at(3) = 3;
        answer.at(4) = 4;
        answer.at(5) = 5;
        answer.at(6) = 6;
    } else if ( iEdge == 2 ) { // edge between nodes 2 3
        answer.at(1) = 4;
        answer.at(2) = 5;
        answer.at(3) = 6;
        answer.at(4) = 7;
        answer.at(5) = 8;
        answer.at(6) = 9;
    } else if ( iEdge == 3 ) { // edge between nodes 3 1
        answer.at(1) = 7;
        answer.at(2) = 8;
        answer.at(3) = 9;
        answer.at(4) = 1;
        answer.at(5) = 2;
        answer.at(6) = 3;
    } else {
        OOFEM_ERROR("wrong edge number");
    }
}

double
DKTPlate :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double detJ = this->interp_lin.edgeGiveTransformationJacobian( iEdge, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    return detJ *gp->giveWeight();
}


void
DKTPlate :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    this->interp_lin.edgeLocal2global( answer, iEdge, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
}


int
DKTPlate :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
{
    // returns transformation matrix from
    // edge local coordinate system
    // to element local c.s
    // (same as global c.s in this case)
    //
    // i.e. f(element local) = T * f(edge local)
    //
    double dx, dy, length;
    IntArray edgeNodes;
    Node *nodeA, *nodeB;

    answer.resize(3, 3);
    answer.zero();

    this->interp_lin.computeLocalEdgeMapping(edgeNodes, iEdge);

    nodeA = this->giveNode( edgeNodes.at(1) );
    nodeB = this->giveNode( edgeNodes.at(2) );

    dx = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
    dy = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);
    length = sqrt(dx * dx + dy * dy);

    answer.at(1, 1) = 1.0;
    answer.at(2, 2) = dx / length;
    answer.at(2, 3) = -dy / length;
    answer.at(3, 2) = dy / length;
    answer.at(3, 3) = dx / length;

    return 1;
}



void
DKTPlate :: computeSurfaceNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *sgp)
{
    this->computeNmatrixAt(sgp->giveNaturalCoordinates(), answer);
}

void
DKTPlate :: giveSurfaceDofMapping(IntArray &answer, int iSurf) const
{
    answer.resize(9);
    answer.zero();
    if ( iSurf == 1 ) {
        for ( int i = 1; i <= 9; i++ ) {
            answer.at(i) = i;
        }
    } else {
        OOFEM_ERROR("wrong surface number");
    }
}

IntegrationRule *
DKTPlate :: GetSurfaceIntegrationRule(int approxOrder)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, this, 1, 1);
    int npoints = iRule->getRequiredNumberOfIntegrationPoints(_Triangle, approxOrder);
    iRule->SetUpPointsOnTriangle(npoints, _Unknown);
    return iRule;
}

double
DKTPlate :: computeSurfaceVolumeAround(GaussPoint *gp, int iSurf)
{
    return this->computeVolumeAround(gp);
}


void
DKTPlate :: computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int isurf)
{
    this->computeGlobalCoordinates( answer, gp->giveNaturalCoordinates() );
}


int
DKTPlate :: computeLoadLSToLRotationMatrix(FloatMatrix &answer, int isurf, GaussPoint *gp)
{
    return 0;
}

void
DKTPlate :: computeVertexBendingMoments(FloatMatrix &answer, TimeStep *tStep)
{
#ifdef DKT_EnableVertexMomentsCache
    if ( stateCounter == tStep->giveSolutionStateCounter() ) {
        answer = vertexMoments;
        return;
    }
#endif

    // the results should be cached somehow, as computing on the fly is highly inefficient
    // due to multiple requests
    answer.resize(5, 3);

    FloatArray eps, m;
    FloatArray coords [ 3 ]; // vertex local coordinates
    coords [ 0 ] = {
        1.0, 0.0
    };
    coords [ 1 ] = {
        0.0, 1.0
    };
    coords [ 2 ] = {
        0.0, 0.0
    };

    GaussIntegrationRule iRule = GaussIntegrationRule(1, this, 1, 1); // dummy rule used to evaluate B at vertices
    iRule.SetUpPointsOnTriangle(1, _Unknown);
    GaussPoint *vgp = iRule.getIntegrationPoint(0);

    for ( int i = 1; i <= this->numberOfDofMans; i++ ) {
        vgp->setNaturalCoordinates(coords [ i - 1 ]);
        this->computeStrainVector(eps, vgp, tStep);
        this->giveStructuralCrossSection()->giveGeneralizedStress_Plate(m, vgp, eps, tStep);
        answer.setColumn(m, i);
    }

#ifdef DKT_EnableVertexMomentsCache
    this->vertexMoments = answer;
    this->stateCounter = tStep->giveSolutionStateCounter();
#endif
}

void
DKTPlate :: computeShearForces(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    // as shear strains are enforced to be zero (art least on element edges) the shear forces are computed from equlibrium
    FloatMatrix m, dndx;
    answer.resize(5);

    this->computeVertexBendingMoments(m, tStep);
    this->interp_lin.evaldNdx( dndx, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    for ( int i = 1; i <= this->numberOfDofMans; i++ ) {
        answer.at(4) += m.at(1, i) * dndx.at(i, 1) + m.at(3, i) * dndx.at(i, 2); //dMxdx + dMxydy
        answer.at(5) += m.at(2, i) * dndx.at(i, 2) + m.at(3, i) * dndx.at(i, 1); //dMydy + dMxydx
    }
}


//
// io routines
//
#ifdef __OOFEG
void
DKTPlate :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    WCRec p [ 3 ];
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( this->giveMaterial()->isActivated(tStep) ) {
        EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
        EASValsSetColor( gc.getElementColor() );
        EASValsSetEdgeColor( gc.getElementEdgeColor() );
        EASValsSetEdgeFlag(true);
        EASValsSetFillStyle(FILL_SOLID);
        EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
        p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
        p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
        p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(3);
        p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
        p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
        p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveCoordinate(3);
        p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
        p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
        p [ 2 ].z = ( FPNum ) this->giveNode(3)->giveCoordinate(3);

        go =  CreateTriangle3D(p);
        EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
        EGAttachObject(go, ( EObjectP ) this);
        EMAddGraphicsToModel(ESIModel(), go);
    }
}


void
DKTPlate :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    WCRec p [ 3 ];
    GraphicObj *go;
    double defScale = gc.getDefScale();

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( this->giveMaterial()->isActivated(tStep) ) {
        EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
        EASValsSetColor( gc.getDeformedElementColor() );
        EASValsSetEdgeColor( gc.getElementEdgeColor() );
        EASValsSetEdgeFlag(true);
        EASValsSetFillStyle(FILL_SOLID);
        EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
        p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
        p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, defScale);
        p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, defScale);
        p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
        p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, defScale);
        p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, defScale);
        p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(1, tStep, defScale);
        p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(2, tStep, defScale);
        p [ 2 ].z = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(3, tStep, defScale);

        go =  CreateTriangle3D(p);
        EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
        EMAddGraphicsToModel(ESIModel(), go);
    }
}


void
DKTPlate :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int i, indx, result = 0;
    WCRec p [ 3 ];
    GraphicObj *tr;
    FloatArray v1, v2, v3;
    double s [ 3 ], defScale;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    if ( !this->giveMaterial()->isActivated(tStep) ) {
        return;
    }

    if ( gc.giveIntVarMode() == ISM_recovered ) {
        result += this->giveInternalStateAtNode(v1, gc.giveIntVarType(), gc.giveIntVarMode(), 1, tStep);
        result += this->giveInternalStateAtNode(v2, gc.giveIntVarType(), gc.giveIntVarMode(), 2, tStep);
        result += this->giveInternalStateAtNode(v3, gc.giveIntVarType(), gc.giveIntVarMode(), 3, tStep);
    } else if ( gc.giveIntVarMode() == ISM_local ) {
        GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
        result += giveIPValue(v1, gp, gc.giveIntVarType(), tStep);
        v2 = v1;
        v3 = v1;
        result *= 3;
    }

    if ( result != 3 ) {
        return;
    }

    indx = gc.giveIntVarIndx();

    s [ 0 ] = v1.at(indx);
    s [ 1 ] = v2.at(indx);
    s [ 2 ] = v3.at(indx);

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);

    if ( gc.getScalarAlgo() == SA_ISO_SURF ) {
        for ( i = 0; i < 3; i++ ) {
            if ( gc.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                defScale = gc.getDefScale();
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                p [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(3, tStep, defScale);
            } else {
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                p [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(3);
            }
        }

        //EASValsSetColor(gc.getYieldPlotColor(ratio));
        gc.updateFringeTableMinMax(s, 3);
        tr =  CreateTriangleWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ]);
        EGWithMaskChangeAttributes(LAYER_MASK, tr);
        EMAddGraphicsToModel(ESIModel(), tr);
    }
}

#endif
} // end namespace oofem
