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

#include "cct.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "structuralms.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "load.h"
#include "structuralcrosssection.h"
#include "mathfem.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {

FEI2dTrLin CCTPlate :: interp_lin(1, 2);
//FEI2dTrRot CCTPlate :: interp_rot(1, 2);

CCTPlate :: CCTPlate(int n, Domain *aDomain) :
    NLStructuralElement(n, aDomain),
    LayeredCrossSectionInterface(), ZZNodalRecoveryModelInterface(),
    NodalAveragingRecoveryModelInterface(), SPRNodalRecoveryModelInterface()
{
    numberOfDofMans = 3;
    numberOfGaussPoints = 1;
    area = 0;
}


FEInterpolation *
CCTPlate :: giveInterpolation(DofIDItem id)
{
    if (id == D_w) {
        return &interp_lin;
    } else {
        return NULL; //&interp_rot;
    }
}


void
CCTPlate :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 5);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Triangle, numberOfGaussPoints, _2dPlate);
    }
}


void
CCTPlate :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *stepN, ValueModeType mode)
// Computes numerically the load vector of the receiver due to the body loads, at stepN.
// load is assumed to be in global cs.
// load vector is then transformed to coordinate system in each node.
// (should be global coordinate system, but there may be defined
//  different coordinate system in each node)
{
    double dens, dV, load;
    GaussPoint *gp = NULL;
    FloatArray force;
    FloatMatrix T;

    if ( ( forLoad->giveBCGeoType() != BodyLoadBGT ) || ( forLoad->giveBCValType() != ForceLoadBVT ) ) {
        _error("computeBodyLoadVectorAt: unknown load type");
    }

    // note: force is assumed to be in global coordinate system.
    forLoad->computeComponentArrayAt(force, stepN, mode);

    if ( force.giveSize() ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

        dens = this->giveMaterial()->give('d', gp);
        dV   = this->computeVolumeAround(gp) * this->giveCrossSection()->give(CS_Thickness);

        answer.resize(9);
        answer.zero();

        load = force.at(1) * dens * dV / 3.0;
        answer.at(1) = load;
        answer.at(4) = load;
        answer.at(7) = load;

        // transform result from global cs to local element cs.
        if ( this->computeGtoLRotationMatrix(T) ) {
            answer.rotatedWith(T, 'n');
        }
    } else {
        answer.resize(0);          // nil resultant
    }
}


void
CCTPlate :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [5x9] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
    // get node coordinates
    double x1, x2, x3, y1, y2, y3;
    this->giveNodeCoordinates(x1, x2, x3, y1, y2, y3);

    //
    double area, b1, b2, b3, c1, c2, c3, l1, l2, l3;

    b1 = y2 - y3;
    b2 = y3 - y1;
    b3 = y1 - y2;

    c1 = x3 - x2;
    c2 = x1 - x3;
    c3 = x2 - x1;

    l1 = 1. / 3.;
    l2 = 1. / 3.;
    l3 = 1. / 3.;

    area = this->computeArea();

    answer.resize(5, 9);
    answer.zero();

    answer.at(1, 3) = b1;
    answer.at(1, 6) = b2;
    answer.at(1, 9) = b3;

    answer.at(2, 2) = -c1;
    answer.at(2, 5) = -c2;
    answer.at(2, 8) = -c3;

    answer.at(3, 2) = -b1;
    answer.at(3, 3) =  c1;
    answer.at(3, 5) = -b2;
    answer.at(3, 6) =  c2;
    answer.at(3, 8) = -b3;
    answer.at(3, 9) =  c3;

    answer.at(4, 1) =  b1;
    answer.at(4, 2) = ( -b1 * b3 * l2 + b1 * b2 * l3 ) * 0.5;
    answer.at(4, 3) = ( -b1 * c3 * l2 - b2 * c3 * l1 + b3 * c2 * l1 + b1 * c2 * l3 ) * 0.5 + l1 * 2. * area;
    answer.at(4, 4) =  b2;
    answer.at(4, 5) = ( -b2 * b1 * l3 + b2 * b3 * l1 ) * 0.5;
    answer.at(4, 6) = ( -b2 * c1 * l3 - b3 * c1 * l2 + b1 * c3 * l2 + b2 * c3 * l1 ) * 0.5 + l2 * 2. * area;
    answer.at(4, 7) =  b3;
    answer.at(4, 8) = ( -b3 * b2 * l1 + b3 * b1 * l2 ) * 0.5;
    answer.at(4, 9) = ( -b3 * c2 * l1 - b1 * c2 * l3 + b2 * c1 * l3 + b3 * c1 * l2 ) * 0.5 + l3 * 2. * area;

    answer.at(5, 1) =  c1;
    answer.at(5, 2) = ( -b3 * c1 * l2 - b3 * c2 * l1 + b2 * c3 * l1 + b2 * c1 * l3 ) * 0.5 - l1 * 2. * area;
    answer.at(5, 3) = ( -c1 * c3 * l2 + c1 * c2 * l3 ) * 0.5;
    answer.at(5, 4) =  c2;
    answer.at(5, 5) = ( -b1 * c2 * l3 - b1 * c3 * l2 + b3 * c1 * l2 + b3 * c2 * l1 ) * 0.5 - l2 * 2. * area;
    answer.at(5, 6) = ( -c2 * c1 * l3 + c2 * c3 * l1 ) * 0.5;
    answer.at(5, 7) =  c3;
    answer.at(5, 8) = ( -b2 * c3 * l1 - b2 * c1 * l3 + b1 * c2 * l3 + b1 * c3 * l2 ) * 0.5 - l3 * 2. * area;
    answer.at(5, 9) = ( -c3 * c2 * l1 + c3 * c1 * l2 ) * 0.5;

    answer.times( 1. / ( 2. * area ) );
}


void
CCTPlate :: computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the [3x9] displacement interpolation matrix {N} of the receiver,
// evaluated at gp.
{
    // get node coordinates
    double x1, x2, x3, y1, y2, y3;
    this->giveNodeCoordinates(x1, x2, x3, y1, y2, y3);

    //
    double l1, l2, l3, b1, b2, b3, c1, c2, c3;

    b1 = y2 - y3;
    b2 = y3 - y1;
    b3 = y1 - y2;

    c1 = x3 - x2;
    c2 = x1 - x3;
    c3 = x2 - x1;

    l1 = gp->giveCoordinate(1);
    l2 = gp->giveCoordinate(2);
    l3 = 1.0 - l1 - l2;

    //
    answer.resize(3, 9);
    answer.zero();

    answer.at(1, 1) = l1;
    answer.at(1, 2) = ( l1 * l2 * b3 - l3 * l1 * b2 ) * 0.5;
    answer.at(1, 3) = ( l1 * l2 * c3 - l3 * l1 * c2 ) * 0.5;
    answer.at(1, 4) = l2;
    answer.at(1, 5) = ( l2 * l3 * b1 - l1 * l2 * b3 ) * 0.5;
    answer.at(1, 6) = ( l2 * l3 * c1 - l1 * l2 * c3 ) * 0.5;
    answer.at(1, 7) = l3;
    answer.at(1, 8) = ( l3 * l1 * b2 - l2 * l3 * b1 ) * 0.5;
    answer.at(1, 9) = ( l3 * l1 * c2 - l2 * l3 * c1 ) * 0.5;

    answer.at(2, 2) = l1;
    answer.at(2, 5) = l2;
    answer.at(2, 8) = l3;

    answer.at(3, 3) = l1;
    answer.at(3, 6) = l2;
    answer.at(3, 9) = l3;
}


void
CCTPlate :: giveNodeCoordinates(double &x1, double &x2, double &x3,
                                double &y1, double &y2, double &y3,
                                double *z)
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

    if ( z ) {
        z [ 0 ] = nc1->at(3);
        z [ 1 ] = nc2->at(3);
        z [ 2 ] = nc3->at(3);
    }
}


IRResultType
CCTPlate :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                 // Required by IR_GIVE_FIELD macro

    this->NLStructuralElement :: initializeFrom(ir);
    numberOfGaussPoints = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, IFT_CCTPlate_nip, "nip"); // Macro
    if ( numberOfGaussPoints != 1 ) {
        numberOfGaussPoints = 1;
    }

    this->computeGaussPoints();
    return IRRT_OK;
}


void
CCTPlate :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    answer.setValues(3, D_w, R_u, R_v);
}


void
CCTPlate :: computeMidPlaneNormal(FloatArray &answer, const GaussPoint *gp)
// returns normal vector to midPlane in GaussPoinr gp of receiver
{
    FloatArray u, v;
    u.beDifferenceOf(*this->giveNode(2)->giveCoordinates(), *this->giveNode(1)->giveCoordinates());
    v.beDifferenceOf(*this->giveNode(3)->giveCoordinates(), *this->giveNode(1)->giveCoordinates());

    answer.beVectorProductOf(u,v);
    answer.normalize();
}


double
CCTPlate :: giveCharacteristicLenght(GaussPoint *gp, const FloatArray &normalToCrackPlane)
// returns receivers characteristic length in gp (for some material models)
// for crack formed in plane with normal normalToCrackPlane.
{
    return this->giveLenghtInDir(normalToCrackPlane);
}


double
CCTPlate :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double detJ, weight;

    weight = gp->giveWeight();
    detJ = fabs( this->interp_lin.giveTransformationJacobian(*gp->giveLocalCoordinates(), FEIElementGeometryWrapper(this)) );
    return detJ * weight; ///@todo What about thickness?
}


void
CCTPlate :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver.
{
    GaussPoint *gp;
    double dV, mss1;

    answer.resize(9, 9);
    answer.zero();

    gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

    dV = this->computeVolumeAround(gp);
    mss1 = dV * this->giveCrossSection()->give(CS_Thickness) * this->giveMaterial()->give('d', gp) / 3.;

    answer.at(1, 1) = mss1;
    answer.at(4, 4) = mss1;
    answer.at(7, 7) = mss1;
}


Interface *
CCTPlate :: giveInterface(InterfaceType interface)
{
    if ( interface == LayeredCrossSectionInterfaceType ) {
        return ( LayeredCrossSectionInterface * ) this;
    } else if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return ( ZZNodalRecoveryModelInterface * ) this;
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return ( NodalAveragingRecoveryModelInterface * ) this;
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return ( SPRNodalRecoveryModelInterface * ) this;
    }

    return NULL;
}


#define POINT_TOL 1.e-3

int
CCTPlate :: computeLocalCoordinates(FloatArray &answer, const FloatArray &coords)
//converts global coordinates to local planar area coordinates,
//does not return a coordinate in the thickness direction, but
//does check that the point is in the element thickness
{
    // get node coordinates
    double x1, x2, x3, y1, y2, y3, z [ 3 ];
    this->giveNodeCoordinates(x1, x2, x3, y1, y2, y3, z);

    // Fetch local coordinates.
    int ok = this->interp_lin.global2local(answer, coords, FEIElementGeometryWrapper(this));

    //get midplane location at this point
    double midplZ;
    midplZ = z [ 0 ] * answer.at(1) + z [ 1 ] * answer.at(2) + z [ 2 ] * answer.at(3);

    //check that the z is within the element
    StructuralCrossSection *cs;
    double elthick;

    cs = ( StructuralCrossSection * ) this->giveCrossSection();
    elthick = cs->give(CS_Thickness);

    if ( elthick / 2.0 + midplZ - fabs( coords.at(3) ) < -POINT_TOL ) {
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

    return ok;
}


int
CCTPlate :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *atTime)
{
    if ( type == IST_ShellForceMomentumTensor ) {
        answer = ( ( StructuralMaterialStatus * ) this->giveMaterial()->giveStatus(gp) )->giveStressVector();
        return 1;
    } else if ( type == IST_ShellStrainCurvatureTensor ) {
        answer = ( ( StructuralMaterialStatus * ) this->giveMaterial()->giveStatus(gp) )->giveStrainVector();
        return 1;
    } else {
        answer.resize(0);
        return 0;
    }
}

//
// The element interface required by ZZNodalRecoveryModel
//
int
CCTPlate :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( ( type == IST_ShellForceMomentumTensor || type == IST_ShellStrainCurvatureTensor) ) {
        return 5;
    }

    return 0;
}


void
CCTPlate :: ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatArray &answer, GaussPoint *gp, InternalStateType type)
// evaluates N matrix (interpolation estimated stress matrix)
// according to Zienkiewicz & Zhu paper
// N(nsigma, nsigma*nnodes)
// Definition : sigmaVector = N * nodalSigmaVector
{
    FloatArray n;

    if ( type == IST_ShellForceMomentumTensor || type == IST_ShellStrainCurvatureTensor ) {
        this->interp_lin.evalN(n, *gp->giveLocalCoordinates(), FEIElementGeometryWrapper(this));
        answer.resize(1, 3);
        answer.zero();
        answer.at(1) = n.at(1);
        answer.at(2) = n.at(2);
        answer.at(3) = n.at(3);
    } else {
        return;
    }
}


//
// The element interface required by NodalAveragingRecoveryModel
//
void
CCTPlate :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                       InternalStateType type, TimeStep *tStep)
{
    GaussPoint *gp;
    if ( type == IST_ShellForceMomentumTensor ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
        answer = ( ( StructuralMaterialStatus * ) this->giveMaterial()->giveStatus(gp) )->giveStressVector();
    } else if ( type == IST_ShellStrainCurvatureTensor ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
        answer = ( ( StructuralMaterialStatus * ) this->giveMaterial()->giveStatus(gp) )->giveStrainVector();
    } else {
        answer.resize(0);
    }
}


void
CCTPlate :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                      InternalStateType type, TimeStep *tStep)
{
    answer.resize(0);
}


//
// The element interface required by SPRNodalRecoveryModelInterface
//
void
CCTPlate :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(3);
    pap.at(1) = this->giveNode(1)->giveNumber();
    pap.at(2) = this->giveNode(2)->giveNumber();
    pap.at(3) = this->giveNode(3)->giveNumber();
}


void
CCTPlate :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    answer.resize(1);
    if ( ( pap == this->giveNode(1)->giveNumber() ) ||
        ( pap == this->giveNode(2)->giveNumber() ) ||
        ( pap == this->giveNode(3)->giveNumber() ) ) {
        answer.at(1) = pap;
    } else {
        _error("SPRNodalRecoveryMI_giveDofMansDeterminedByPatch: node unknown");
    }
}


void
CCTPlate :: SPRNodalRecoveryMI_computeIPGlobalCoordinates(FloatArray &coords, GaussPoint *gp)
{
    if ( gp == integrationRulesArray [ 0 ]->getIntegrationPoint(0) ) {
        this->computeGlobalCoordinates( coords, * gp->giveCoordinates() );
    } else {
        _error("SPRNodalRecoveryMI_computeIPGlobalCoordinates: unsupported ip num");
    }
}


SPRPatchType
CCTPlate :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_2dxy;
}


//
// layered cross section support functions
//
void
CCTPlate :: computeStrainVectorInLayer(FloatArray &answer, GaussPoint *masterGp,
                                       GaussPoint *slaveGp, TimeStep *tStep)
// returns full 3d strain vector of given layer (whose z-coordinate from center-line is
// stored in slaveGp) for given tStep
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

    answer.at(1) = masterGpStrain.at(1) * layerZCoord;
    answer.at(2) = masterGpStrain.at(2) * layerZCoord;
    answer.at(6) = masterGpStrain.at(3) * layerZCoord;
    answer.at(4) = masterGpStrain.at(5);
    answer.at(5) = masterGpStrain.at(4);
}


//
// io routines
//
#ifdef __OOFEG
void
CCTPlate  :: drawRawGeometry(oofegGraphicContext &gc)
{
    WCRec p [ 3 ];
    GraphicObj *go;
    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();

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
CCTPlate  :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
{
    WCRec p [ 3 ];
    GraphicObj *go;
    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
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
        p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
        p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
        p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);
        p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
        p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
        p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);
        p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
        p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
        p [ 2 ].z = ( FPNum ) this->giveNode(3)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);

        go =  CreateTriangle3D(p);
        EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
        EMAddGraphicsToModel(ESIModel(), go);
    }
}


void
CCTPlate  :: drawScalar(oofegGraphicContext &context)
{
    int i, indx, result = 0;
    WCRec p [ 3 ];
    GraphicObj *tr;
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    FloatArray v1, v2, v3;
    double s [ 3 ], defScale;
    IntArray map;

    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    if ( !this->giveMaterial()->isActivated(tStep) ) {
        return;
    }

    if ( context.giveIntVarMode() == ISM_recovered ) {
        result += this->giveInternalStateAtNode(v1, context.giveIntVarType(), context.giveIntVarMode(), 1, tStep);
        result += this->giveInternalStateAtNode(v2, context.giveIntVarType(), context.giveIntVarMode(), 2, tStep);
        result += this->giveInternalStateAtNode(v3, context.giveIntVarType(), context.giveIntVarMode(), 3, tStep);
    } else if ( context.giveIntVarMode() == ISM_local ) {
        GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
        result += giveIPValue(v1, gp, context.giveIntVarType(), tStep);
        v2 = v1;
        v3 = v1;
        result *= 3;
    }

    if ( result != 3 ) {
        return;
    }

    this->giveIntVarCompFullIndx( map, context.giveIntVarType() );

    if ( ( indx = map.at( context.giveIntVarIndx() ) ) == 0 ) {
        return;
    }

    s [ 0 ] = v1.at(indx);
    s [ 1 ] = v2.at(indx);
    s [ 2 ] = v3.at(indx);

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);

    if ( context.getScalarAlgo() == SA_ISO_SURF ) {
        for ( i = 0; i < 3; i++ ) {
            if ( context.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                defScale = context.getDefScale();
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
                p [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);
            } else {
                p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                p [ i ].z = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(3);
            }
        }

        //EASValsSetColor(gc.getYieldPlotColor(ratio));
        context.updateFringeTableMinMax(s, 3);
        tr =  CreateTriangleWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ]);
        EGWithMaskChangeAttributes(LAYER_MASK, tr);
        EMAddGraphicsToModel(ESIModel(), tr);
    }
}

/*
 * void
 * CCTPlate  :: drawInternalState (oofegGraphicContext & gc)
 * //
 * // Draws internal state graphics representation
 * //
 * {
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
 * (mode == szxForce) || (mode == syzForce))) return;
 *
 *
 *
 * EASValsSetLayer(OOFEG_STRESS_CONTOUR_LAYER);
 * if (gc.getInternalVarsDefGeoFlag()) {
 * // use deformed geometry
 * p[0].x = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(1,tStep,DisplacementVector,defScale);
 * p[0].y = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(2,tStep,DisplacementVector,defScale);
 * p[0].z = (FPNum) this->giveNode(1)->giveUpdatedCoordinate(3,tStep,DisplacementVector,defScale);
 * p[1].x = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(1,tStep,DisplacementVector,defScale);
 * p[1].y = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(2,tStep,DisplacementVector,defScale);
 * p[1].z = (FPNum) this->giveNode(2)->giveUpdatedCoordinate(3,tStep,DisplacementVector,defScale);
 * p[2].x = (FPNum) this->giveNode(3)->giveUpdatedCoordinate(1,tStep,DisplacementVector,defScale);
 * p[2].y = (FPNum) this->giveNode(3)->giveUpdatedCoordinate(2,tStep,DisplacementVector,defScale);
 * p[2].z = (FPNum) this->giveNode(3)->giveUpdatedCoordinate(3,tStep,DisplacementVector,defScale);
 * } else {
 * p[0].x = (FPNum) this->giveNode(1)->giveCoordinate(1);
 * p[0].y = (FPNum) this->giveNode(1)->giveCoordinate(2);
 * p[0].z = (FPNum) this->giveNode(1)->giveCoordinate(3);
 * p[1].x = (FPNum) this->giveNode(2)->giveCoordinate(1);
 * p[1].y = (FPNum) this->giveNode(2)->giveCoordinate(2);
 * p[1].z = (FPNum) this->giveNode(2)->giveCoordinate(3);
 * p[2].x = (FPNum) this->giveNode(3)->giveCoordinate(1);
 * p[2].y = (FPNum) this->giveNode(3)->giveCoordinate(2);
 * p[2].z = (FPNum) this->giveNode(3)->giveCoordinate(3);
 * }
 *
 * int result = 0;
 * result+= this->giveInternalStateAtNode (gc, 1, &v1);
 * result+= this->giveInternalStateAtNode (gc, 2, &v2);
 * result+= this->giveInternalStateAtNode (gc, 3, &v3);
 *
 * if (result == 3) {
 * tr = CreateTriangleWD3D (p,v1,v2,v3);
 * EGWithMaskChangeAttributes(LAYER_MASK, tr);
 * EMAddGraphicsToModel(ESIModel(), tr);
 * }
 * }
 */

#endif
} // end namespace oofem
