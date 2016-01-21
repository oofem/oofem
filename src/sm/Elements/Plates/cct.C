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

#include "../sm/Elements/Plates/cct.h"
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
REGISTER_Element(CCTPlate);

FEI2dTrLin CCTPlate :: interp_lin(1, 2);
//FEI2dTrRot CCTPlate :: interp_rot(1, 2);

CCTPlate :: CCTPlate(int n, Domain *aDomain) :
    NLStructuralElement(n, aDomain),
    LayeredCrossSectionInterface(), ZZNodalRecoveryModelInterface(this),
    NodalAveragingRecoveryModelInterface(), SPRNodalRecoveryModelInterface(), ZZErrorEstimatorInterface(this)
{
    numberOfDofMans = 3;
    numberOfGaussPoints = 1;
    area = 0;
}


FEInterpolation *
CCTPlate :: giveInterpolation() const { return & interp_lin; }


FEInterpolation *
CCTPlate :: giveInterpolation(DofIDItem id) const
{
    if ( id == D_w ) {
        return & interp_lin;
    } else {
        return NULL; //&interp_rot;
    }
}


void
CCTPlate :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 5) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}


void
CCTPlate :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode)
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
        // constant density and thickness assumed
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
CCTPlate :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [5x9] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
    // get node coordinates
  double x1, x2, x3, y1, y2, y3, z1, z2, z3;
  this->giveNodeCoordinates(x1, x2, x3, y1, y2, y3, z1, z2, z3);

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
CCTPlate :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the [3x9] displacement interpolation matrix {N} of the receiver,
// evaluated at gp.
{
    // get node coordinates
  double x1, x2, x3, y1, y2, y3, z1, z2, z3;
  this->giveNodeCoordinates(x1, x2, x3, y1, y2, y3, z1, z2, z3);

    //
    double l1, l2, l3, b1, b2, b3, c1, c2, c3;

    b1 = y2 - y3;
    b2 = y3 - y1;
    b3 = y1 - y2;

    c1 = x3 - x2;
    c2 = x1 - x3;
    c3 = x2 - x1;

    l1 = iLocCoord.at(1);
    l2 = iLocCoord.at(2);
    l3 = 1.0 - l1 - l2;

    //
    answer.resize(3, 9);
    answer.zero();

    answer.at(1, 1) = l1;
    answer.at(1, 2) = l1 * ( l2 * b3 - l3 * b2 ) * 0.5;
    answer.at(1, 3) = l1 * ( l2 * c3 - l3 * c2 ) * 0.5;
    answer.at(1, 4) = l2;
    answer.at(1, 5) = l2 * ( l3 * b1 - l1 * b3 ) * 0.5;
    answer.at(1, 6) = l2 * ( l3 * c1 - l1 * c3 ) * 0.5;
    answer.at(1, 7) = l3;
    answer.at(1, 8) = l3 * ( l1 * b2 - l2 * b1 ) * 0.5;
    answer.at(1, 9) = l3 * ( l1 * c2 - l2 * c1 ) * 0.5;

    answer.at(2, 2) = l1;
    answer.at(2, 5) = l2;
    answer.at(2, 8) = l3;

    answer.at(3, 3) = l1;
    answer.at(3, 6) = l2;
    answer.at(3, 9) = l3;
}


void
CCTPlate :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveGeneralizedStress_Plate(answer, gp, strain, tStep);
}


void
CCTPlate :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->give2dPlateStiffMtrx(answer, rMode, gp, tStep);
}


void
CCTPlate :: giveNodeCoordinates(double &x1, double &x2, double &x3,
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

double
CCTPlate :: computeArea ()
// returns the area occupied by the receiver
{
  // get node coordinates
  double x1, x2, x3, y1, y2, y3, z1, z2, z3;
  this->giveNodeCoordinates(x1, x2, x3, y1, y2, y3, z1, z2, z3);
  
  if (area > 0) return area;  // check if previously computed

  return (area = 0.5*(x2*y3+x1*y2+y1*x3-x2*y1-x3*y2-x1*y3)) ;
 
}

IRResultType
CCTPlate :: initializeFrom(InputRecord *ir)
{
    numberOfGaussPoints = 1;
    return NLStructuralElement :: initializeFrom(ir);
}


void
CCTPlate :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_w, R_u, R_v};
}


void
CCTPlate :: computeMidPlaneNormal(FloatArray &answer, const GaussPoint *gp)
// returns normal vector to midPlane in GaussPoinr gp of receiver
{
    FloatArray u, v;
    u.beDifferenceOf( * this->giveNode(2)->giveCoordinates(), * this->giveNode(1)->giveCoordinates() );
    v.beDifferenceOf( * this->giveNode(3)->giveCoordinates(), * this->giveNode(1)->giveCoordinates() );

    answer.beVectorProductOf(u, v);
    answer.normalize();
}


double
CCTPlate :: giveCharacteristicLength(const FloatArray &normalToCrackPlane)
// returns receiver's characteristic length (for crack band models)
// for a crack formed in the plane with normal normalToCrackPlane.
{
    return this->giveCharacteristicLengthForPlaneElements(normalToCrackPlane);
}


double
CCTPlate :: computeVolumeAround(GaussPoint *gp)
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
CCTPlate :: computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep)
// Returns the lumped mass matrix of the receiver.
{
    answer.resize(9, 9);
    answer.zero();

    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

    double dV = this->computeVolumeAround(gp);
    // constant thickness and density assumed
    double density = this->giveStructuralCrossSection()->give('d', gp);
    double mss1 = dV * this->giveCrossSection()->give(CS_Thickness, gp) * density / 3.;

    answer.at(1, 1) = mss1;
    answer.at(4, 4) = mss1;
    answer.at(7, 7) = mss1;
}


Interface *
CCTPlate :: giveInterface(InterfaceType interface)
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
CCTPlate :: computeLocalCoordinates(FloatArray &answer, const FloatArray &coords)
//converts global coordinates to local planar area coordinates,
//does not return a coordinate in the thickness direction, but
//does check that the point is in the element thickness
{
    // get node coordinates
    double x1, x2, x3, y1, y2, y3, z1, z2, z3;
    this->giveNodeCoordinates(x1, x2, x3, y1, y2, y3, z1, z2, z3);

    // Fetch local coordinates.
    bool ok = this->interp_lin.global2local( answer, coords, FEIElementGeometryWrapper(this) ) > 0;

    //get midplane location at this point
    double midplZ;
    midplZ = z1 * answer.at(1) + z2 * answer.at(2) + z3 * answer.at(3);

    //check that the z is within the element
    GaussPoint _gp(NULL, 1, answer, 1.0, _2dPlate);
    StructuralCrossSection *cs = this->giveStructuralCrossSection();
    double elthick = cs->give(CS_Thickness, & _gp);

    if ( elthick / 2.0 + midplZ - fabs( coords.at(3) ) < -POINT_TOL ) {
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

    return ok;
}


int
CCTPlate :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
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
CCTPlate :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
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
        OOFEM_ERROR("node unknown");
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
CCTPlate :: computeStrainVectorInLayer(FloatArray &answer, const FloatArray &masterGpStrain,
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
CCTPlate :: computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
{
    IntArray edgeNodes;
    FloatArray n;
    double b, c, n12;

    this->interp_lin.edgeEvalN( n, iedge, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    this->interp_lin.computeLocalEdgeMapping(edgeNodes, iedge);

    n12 = 0.5 * n.at(1) * n.at(2);
    b = this->giveNode( edgeNodes.at(1) )->giveCoordinate(2) - this->giveNode( edgeNodes.at(2) )->giveCoordinate(2);
    c = this->giveNode( edgeNodes.at(2) )->giveCoordinate(1) - this->giveNode( edgeNodes.at(1) )->giveCoordinate(1);


    answer.resize(3, 6);
    answer.at(1, 1) = n.at(1);
    answer.at(1, 2) = n12 * b;
    answer.at(1, 3) = n12 * c;
    answer.at(1, 4)  = n.at(2);
    answer.at(1, 5) = -n12 * b;
    answer.at(1, 6) = -n12 * c;
    //
    answer.at(2, 2) = answer.at(3, 3) = n.at(1);
    answer.at(2, 5) = answer.at(3, 6) = n.at(2);
}

void
CCTPlate :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
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
CCTPlate :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double detJ = this->interp_lin.edgeGiveTransformationJacobian( iEdge, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    return detJ *gp->giveWeight();
}


void
CCTPlate :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    this->interp_lin.edgeLocal2global( answer, iEdge, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
}


int
CCTPlate :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
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


//
// io routines
//
#ifdef __OOFEG
void
CCTPlate :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
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
CCTPlate :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
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
CCTPlate :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
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
