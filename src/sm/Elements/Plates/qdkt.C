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

#include "../sm/Elements/Plates/qdkt.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/CrossSections/structuralcrosssection.h"
#include "fei2dquadlin.h"
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
REGISTER_Element(QDKTPlate);

FEI2dQuadLin QDKTPlate :: interp_lin(1, 2);

QDKTPlate :: QDKTPlate(int n, Domain *aDomain) :
    NLStructuralElement(n, aDomain),
    ZZNodalRecoveryModelInterface(this),
    SPRNodalRecoveryModelInterface(),
    ZZErrorEstimatorInterface(this)
{
    numberOfDofMans = 4;
    numberOfGaussPoints = 4;
}


FEInterpolation *
QDKTPlate :: giveInterpolation() const { return & interp_lin; }


FEInterpolation *
QDKTPlate :: giveInterpolation(DofIDItem id) const
{
    return & interp_lin;
}


void
QDKTPlate :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 5) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);
    }
}

void
QDKTPlate :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the [5x12] strain-displacement matrix {B} of the receiver,
// evaluated at gp.
{
    // get node coordinates
  double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
  this->giveNodeCoordinates(x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4);

    // get gp coordinates
    double ksi = gp->giveNaturalCoordinate(1);
    double eta = gp->giveNaturalCoordinate(2);

    // geometrical characteristics of element sides
    double dx4 = x2-x1;
    double dy4 = y2-y1;
    double l4 = sqrt(dx4*dx4+dy4*dy4);
    
    double dx5 = x3-x2;
    double dy5 = y3-y2;
    double l5 = sqrt(dx5*dx5+dy5*dy5);
    
    double dx6 = x4-x3;
    double dy6 = y4-y3;
    double l6 = sqrt(dx6*dx6+dy6*dy6);

    double dx7 = x1-x4;
    double dy7 = y1-y4;
    double l7 = sqrt(dx7*dx7+dy7*dy7);

    double c4 = dy4/l4;
    double s4 = -dx4/l4;
    
    double c5 = dy5/l5;
    double s5 = -dx5/l5;
    
    double c6 = dy6/l6;
    double s6 = -dx6/l6;
    
    double c7 = dy7/l7;
    double s7 = -dx7/l7;

    // transformation matrix from vertex dofs (w,fi_x, fi_y) to quadratic rotation DOFs
    double T101 = -3./2./l4*c4;
    double T102 = -1./4.*c4*c4+1./2.*s4*s4;
    double T103 = -1./4.*c4*s4-1./2.*c4*s4;
    double T104 =  3./2./l4*c4;
    double T105 = -1./4.*c4*c4+1./2.*s4*s4;
    double T106 = -1./4.*c4*s4-1./2.*c4*s4;

    double T201 = -3./2./l4*s4;
    double T202 = -1./4.*c4*s4-1./2.*c4*s4;
    double T203 = -1./4.*s4*s4+1./2.*c4*c4;
    double T204 =  3./2./l4*s4;
    double T205 = -1./4.*c4*s4-1./2.*c4*s4;
    double T206 = -1./4.*s4*s4+1./2.*c4*c4;

    double T304 = -3./2./l5*c5;
    double T305 = -1./4.*c5*c5+1./2.*s5*s5;
    double T306 = -1./4.*c5*s5-1./2.*c5*s5;
    double T307 =  3./2./l5*c5;
    double T308 = -1./4.*c5*c5+1./2.*s5*s5;
    double T309 = -1./4.*c5*s5-1./2.*c5*s5;

    double T404 = -3./2./l5*s5;
    double T405 = -1./4.*c5*s5-1./2.*c5*s5;
    double T406 = -1./4.*s5*s5+1./2.*c5*c5;
    double T407 =  3./2./l5*s5;
    double T408 = -1./4.*c5*s5-1./2.*c5*s5;
    double T409 = -1./4.*s5*s5+1./2.*c5*c5;

    double T507 = -3./2./l6*c6;
    double T508 = -1./4.*c6*c6+1./2.*s6*s6;
    double T509 = -1./4.*c6*s6-1./2.*c6*s6;
    double T510 =  3./2./l6*c6;
    double T511 = -1./4.*c6*c6+1./2.*s6*s6;
    double T512 = -1./4.*c6*s6-1./2.*c6*s6;

    double T607 = -3./2./l6*s6;
    double T608 = -1./4.*c6*s6-1./2.*c6*s6;
    double T609 = -1./4.*s6*s6+1./2.*c6*c6;
    double T610 =  3./2./l6*s6;
    double T611 = -1./4.*c6*s6-1./2.*c6*s6;
    double T612 = -1./4.*s6*s6+1./2.*c6*c6;

    double T701 =  3./2./l7*c7;
    double T702 = -1./4.*c7*c7+1./2.*s7*s7;
    double T703 = -1./4.*c7*s7-1./2.*c7*s7;
    double T710 = -3./2./l7*c7;
    double T711 = -1./4.*c7*c7+1./2.*s7*s7;
    double T712 = -1./4.*c7*s7-1./2.*c7*s7;

    double T801 =  3./2./l7*s7;
    double T802 = -1./4.*c7*s7-1./2.*c7*s7;
    double T803 = -1./4.*s7*s7+1./2.*c7*c7;
    double T810 = -3./2./l7*s7;
    double T811 = -1./4.*c7*s7-1./2.*c7*s7;
    double T812 = -1./4.*s7*s7+1./2.*c7*c7;

    // derivatives of quadratic interpolation functions
    // we do not have "midside" nodes -> explicit here
    double N1dk =  0.25*(2.0*ksi+eta)*(1.0+eta);
    double N2dk =  0.25*(2.0*ksi-eta)*(1.0+eta);
    double N3dk =  0.25*(2.0*ksi+eta)*(1.0-eta);
    double N4dk =  0.25*(2.0*ksi-eta)*(1.0-eta);
    double N7dk = -ksi*(1.0-eta);
    double N8dk =  0.5*(1.0-eta*eta);
    double N5dk = -ksi*(1.0+eta);
    double N6dk = -0.5*(1.0-eta*eta);

    double N3de =  0.25*(2.0*eta+ksi)*(1.0-ksi);
    double N4de =  0.25*(2.0*eta-ksi)*(1.0+ksi);
    double N1de =  0.25*(2.0*eta+ksi)*(1.0+ksi);
    double N2de =  0.25*(2.0*eta-ksi)*(1.0-ksi);
    double N7de = -0.5*(1.0-ksi*ksi);
    double N8de = -eta*(1.0+ksi);
    double N5de =  0.5*(1.0-ksi*ksi);
    double N6de = -eta*(1.0-ksi);

    double detJ = 1./8.*((y4-y2)*(x3-x1)-(y3-y1)*(x4-x2))+
      ksi/8*((y3-y4)*(x2-x1)-(y2-y1)*(x3-x4))+
      eta/8*((y4-y1)*(x3-x2)-(y3-y2)*(x4-x1));
    
    double dxdk = -1.0/detJ * ((y3-y2)+(y4-y1+ksi*(y1-y2+y3-y4)))/4.0;
    double dxde =  1.0/detJ * ((y2-y1)+(y3-y4+eta*(y1-y2+y3-y4)))/4.0;
    double dydk =  1.0/detJ * ((x3-x2)+(x4-x1+ksi*(x1-x2+x3-x4)))/4.0;
    double dyde = -1.0/detJ * ((x2-x1)+(x3-x4+eta*(x1-x2+x3-x4)))/4.0;

    double dN102 = N1dk*dxdk+N1de*dxde;
    double dN104 = N2dk*dxdk+N2de*dxde;
    double dN106 = N3dk*dxdk+N3de*dxde;
    double dN108 = N4dk*dxdk+N4de*dxde;
    double dN110 = N5dk*dxdk+N5de*dxde;
    double dN112 = N6dk*dxdk+N6de*dxde;
    double dN114 = N7dk*dxdk+N7de*dxde;
    double dN116 = N8dk*dxdk+N8de*dxde;

    double dN201 = -N1dk*dydk-N1de*dyde;
    double dN203 = -N2dk*dydk-N2de*dyde;
    double dN205 = -N3dk*dydk-N3de*dyde;
    double dN207 = -N4dk*dydk-N4de*dyde; 
    double dN209 = -N5dk*dydk-N5de*dyde;
    double dN211 = -N6dk*dydk-N6de*dyde;
    double dN213 = -N7dk*dydk-N7de*dyde;
    double dN215 = -N8dk*dydk-N8de*dyde;

    answer.resize(5, 12);
    answer.zero();
    
    answer.at(1,1) = T201*dN110 + T801*dN116;
    answer.at(1,2) = T202*dN110 + T802*dN116;
    answer.at(1,3) = dN102 + T203*dN110 + T803*dN116;
    answer.at(1,4) = T204*dN110 + T404*dN112;
    answer.at(1,5) = T205*dN110 + T405*dN112;
    answer.at(1,6) = dN104 + T206*dN110 + T406*dN112;
    answer.at(1,7) = T407*dN112 + T607*dN114;
    answer.at(1,8) = T408*dN112 + T608*dN114;
    answer.at(1,9) = dN106 + T409*dN112 + T609*dN114;
    answer.at(1,10)= T610*dN114 + T810*dN116;
    answer.at(1,11)= T611*dN114 + T811*dN116;
    answer.at(1,12)= dN108 + T612*dN114 + T812*dN116;

    answer.at(2,1) = T101*dN209 + T701*dN215;
    answer.at(2,2) = dN201 + T102*dN209 + T702*dN215;
    answer.at(2,3) = T103*dN209 + T703*dN215;
    answer.at(2,4) = T104*dN209 + T304*dN211;
    answer.at(2,5) = dN203 + T105*dN209 + T305*dN211;
    answer.at(2,6) = T106*dN209 + T306*dN211;
    answer.at(2,7) = T307*dN211 + T507*dN213;
    answer.at(2,8) = dN205 + T308*dN211 + T508*dN213;
    answer.at(2,9) = T309*dN211 + T509*dN213;
    answer.at(2,10)= T510*dN213 + T710*dN215;
    answer.at(2,11)= dN207 + T511*dN213 + T711*dN215;
    answer.at(2,12)= T512*dN213 + T712*dN215;

    answer.at(3,1) = - T101*dN110 - T201*dN209 - T701*dN116 - T801*dN215;
    answer.at(3,2) = - dN102 - T102*dN110 - T202*dN209 - T702*dN116 - T802*dN215;
    answer.at(3,3) = - dN201 - T103*dN110 - T203*dN209 - T703*dN116 - T803*dN215;
    answer.at(3,4) = - T104*dN110 - T204*dN209 - T304*dN112 - T404*dN211;
    answer.at(3,5) = - dN104 - T105*dN110 - T205*dN209 - T305*dN112 - T405*dN211;
    answer.at(3,6) = - dN203 - T106*dN110 - T206*dN209 - T306*dN112 - T406*dN211;
    answer.at(3,7) = - T307*dN112 - T407*dN211 - T507*dN114 - T607*dN213;
    answer.at(3,8) = - dN106 - T308*dN112 - T408*dN211 - T508*dN114 - T608*dN213;
    answer.at(3,9) = - dN205 - T309*dN112 - T409*dN211 - T509*dN114 - T609*dN213;
    answer.at(3,10)= - T510*dN114 - T610*dN213 - T710*dN116 - T810*dN215;
    answer.at(3,11)= - dN108 - T511*dN114 - T611*dN213 - T711*dN116 - T811*dN215;
    answer.at(3,12)= - dN207 - T512*dN114 - T612*dN213 - T712*dN116 - T812*dN215;

    // Note: no shear strains, no shear forces => the 4th and 5th rows are zero
}


void
QDKTPlate :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the [3x9] displacement interpolation matrix {N} of the receiver,
// evaluated at gp.
// Note: this interpolation is not available, as the deflection is cubic along the edges,
//       but not define in the interior of the element
// Note: the interpolation of rotations is quadratic
// NOTE: linear interpolation returned instead
{
    FloatArray N;

    giveInterpolation()->evalN( N, iLocCoord, FEIElementGeometryWrapper(this) );

    answer.beNMatrixOf(N, 3);

}


void
QDKTPlate :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveGeneralizedStress_Plate(answer, gp, strain, tStep);
}


void
QDKTPlate :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->give2dPlateStiffMtrx(answer, rMode, gp, tStep);
}


void
QDKTPlate :: giveNodeCoordinates(double &x1, double &x2, double &x3, double &x4,
                                 double &y1, double &y2, double &y3, double &y4,
                                 double &z1, double &z2, double &z3, double &z4)
{
  FloatArray *nc1, *nc2, *nc3, *nc4;
    nc1 = this->giveNode(1)->giveCoordinates();
    nc2 = this->giveNode(2)->giveCoordinates();
    nc3 = this->giveNode(3)->giveCoordinates();
    nc4 = this->giveNode(4)->giveCoordinates();

    x1 = nc1->at(1);
    x2 = nc2->at(1);
    x3 = nc3->at(1);
    x4 = nc4->at(1);

    y1 = nc1->at(2);
    y2 = nc2->at(2);
    y3 = nc3->at(2);
    y4 = nc4->at(2);

    z1 = nc1->at(3);
    z2 = nc2->at(3);
    z3 = nc3->at(3);
    z4 = nc4->at(3);
    
}


IRResultType
QDKTPlate :: initializeFrom(InputRecord *ir)
{
    return NLStructuralElement :: initializeFrom(ir);
}


void
QDKTPlate :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_w, R_u, R_v};
}


void
QDKTPlate :: computeMidPlaneNormal(FloatArray &answer, const GaussPoint *gp)
// returns normal vector to midPlane in GaussPoinr gp of receiver
{
    FloatArray u, v;
    u.beDifferenceOf( * this->giveNode(2)->giveCoordinates(), * this->giveNode(1)->giveCoordinates() );
    v.beDifferenceOf( * this->giveNode(3)->giveCoordinates(), * this->giveNode(1)->giveCoordinates() );

    answer.beVectorProductOf(u, v);
    answer.normalize();
}


double
QDKTPlate :: giveCharacteristicLength(const FloatArray &normalToCrackPlane)
//
// returns receiver's characteristic length for crack band models
// for a crack formed in the plane with normal normalToCrackPlane.
//
{
    return this->giveCharacteristicLengthForPlaneElements(normalToCrackPlane);
}


double
QDKTPlate :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double detJ, weight;

    weight = gp->giveWeight();
    detJ = fabs( this->interp_lin.giveTransformationJacobian( gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) ) );
    return detJ * weight; ///@todo What about thickness?
}


Interface *
QDKTPlate :: giveInterface(InterfaceType interface)
{
    if ( interface == LayeredCrossSectionInterfaceType ) {
        return static_cast< LayeredCrossSectionInterface * >(this);
    } else if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return static_cast< ZZNodalRecoveryModelInterface * >(this);
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return static_cast< SPRNodalRecoveryModelInterface * >(this);
    } else if ( interface == ZZErrorEstimatorInterfaceType ) {
        return static_cast< ZZErrorEstimatorInterface * >(this);
    }


    return NULL;
}

void
QDKTPlate :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode)
// Computes numerically the load vector of the receiver due to the body
// loads, at tStep.
// load is assumed to be in global cs.
// load vector is then transformed to coordinate system in each node.
// (should be global coordinate system, but there may be defined
//  different coordinate system in each node)
{
    double dens, dV;
    FloatArray force, ntf;
    FloatMatrix n, T;

    if ( ( forLoad->giveBCGeoType() != BodyLoadBGT ) || ( forLoad->giveBCValType() != ForceLoadBVT ) ) {
        OOFEM_ERROR("unknown load type");
    }

    // note: force is assumed to be in global coordinate system.
    forLoad->computeComponentArrayAt(force, tStep, mode);
    // transform from global to element local c.s
    if ( this->computeLoadGToLRotationMtrx(T) ) {
        force.rotatedWith(T, 'n');
    }

    answer.clear();

    if ( force.giveSize() ) {
        for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {
            this->computeNmatrixAt(gp->giveSubPatchCoordinates(), n);
            dV  = this->computeVolumeAround(gp) * this->giveCrossSection()->give(CS_Thickness, gp);
            dens = this->giveCrossSection()->give('d', gp);
            ntf.beTProductOf(n, force);
            answer.add(dV * dens, ntf);
        }
    } else {
        return;
    }
}


#define POINT_TOL 1.e-3

bool
QDKTPlate :: computeLocalCoordinates(FloatArray &answer, const FloatArray &coords)
//converts global coordinates to local planar area coordinates,
//does not return a coordinate in the thickness direction, but
//does check that the point is in the element thickness
{
    // get node coordinates
  double x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4;
    this->giveNodeCoordinates(x1, x2, x3, x4, y1, y2, y3, y4, z1, z2, z3, z4);

    // Fetch local coordinates.
    bool ok = this->interp_lin.global2local( answer, coords, FEIElementGeometryWrapper(this) ) > 0;

    //check that the point is in the element and set flag
    for ( int i = 1; i <= 4; i++ ) {
        if ( answer.at(i) < ( 0. - POINT_TOL ) ) {
            return false;
        }

        if ( answer.at(i) > ( 1. + POINT_TOL ) ) {
            return false;
        }
    }

    //get midplane location at this point
    double midplZ;
    midplZ = z1 * answer.at(1) + z2 * answer.at(2) + z3 * answer.at(3) + z4 * answer.at(4);

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
QDKTPlate :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
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
// The element interface required by SPRNodalRecoveryModelInterface
//
void
QDKTPlate :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(4);
    pap.at(1) = this->giveNode(1)->giveNumber();
    pap.at(2) = this->giveNode(2)->giveNumber();
    pap.at(3) = this->giveNode(3)->giveNumber();
    pap.at(4) = this->giveNode(4)->giveNumber();
}


void
QDKTPlate :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    answer.resize(1);
    if ( ( pap == this->giveNode(1)->giveNumber() ) ||
        ( pap == this->giveNode(2)->giveNumber() ) ||
        ( pap == this->giveNode(3)->giveNumber() ) ||
        ( pap == this->giveNode(4)->giveNumber() ) ) {
      answer.at(1) = pap;
    } else {
        OOFEM_ERROR("node unknown");
    }
}


SPRPatchType
QDKTPlate :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_2dxy;
}


//
// layered cross section support functions
//
void
QDKTPlate :: computeStrainVectorInLayer(FloatArray &answer, const FloatArray &masterGpStrain,
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
QDKTPlate :: computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
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
QDKTPlate :: giveEdgeDofMapping(IntArray &answer, int iEdge) const
{
    if ( iEdge == 1 ) { // edge between nodes 1 2
      answer = {1, 2, 3, 4, 5, 6};
    } else if ( iEdge == 2 ) { // edge between nodes 2 3
      answer = {4, 5, 6, 7, 8, 9};
    } else if ( iEdge == 3 ) { // edge between nodes 3 4
      answer = {7, 8, 9, 10, 11, 12};
    } else if ( iEdge == 4 ) { // edge between nodes 4 1
      answer = {10, 11, 12, 1, 2, 3};
    } else {
        OOFEM_ERROR("wrong edge number");
    }
}

double
QDKTPlate :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
    double detJ = this->interp_lin.edgeGiveTransformationJacobian( iEdge, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    return detJ *gp->giveWeight();
}


void
QDKTPlate :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
    this->interp_lin.edgeLocal2global( answer, iEdge, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
}


int
QDKTPlate :: computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp)
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
QDKTPlate :: computeSurfaceNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *sgp)
{
    this->computeNmatrixAt(sgp->giveNaturalCoordinates(), answer);
}

void
QDKTPlate :: giveSurfaceDofMapping(IntArray &answer, int iSurf) const
{
    answer.resize(12);
    answer.zero();
    if ( iSurf == 1 ) {
        for (int i = 1; i<=12; i++) {
            answer.at(i) = i;
        }
    } else {
        OOFEM_ERROR("wrong surface number");
    }
}

IntegrationRule *
QDKTPlate :: GetSurfaceIntegrationRule(int approxOrder)
{
    IntegrationRule *iRule = new GaussIntegrationRule(1, this, 1, 1);
    int npoints = iRule->getRequiredNumberOfIntegrationPoints(_Square, approxOrder);
    iRule->SetUpPointsOnSquare(npoints, _Unknown);
    return iRule;
}

double
QDKTPlate :: computeSurfaceVolumeAround(GaussPoint *gp, int iSurf)
{
    return this->computeVolumeAround(gp);
}


void
QDKTPlate :: computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int isurf)
{
    this->computeGlobalCoordinates( answer, gp->giveNaturalCoordinates() );
}


int
QDKTPlate :: computeLoadLSToLRotationMatrix(FloatMatrix &answer, int isurf, GaussPoint *gp)
{
    return 0;
}


//
// io routines
//
#ifdef __OOFEG
void
QDKTPlate :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    WCRec p [ 4 ];
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
        p [ 3 ].x = ( FPNum ) this->giveNode(4)->giveCoordinate(1);
        p [ 3 ].y = ( FPNum ) this->giveNode(4)->giveCoordinate(2);
        p [ 3 ].z = 0.;

        go =  CreateQuad3D(p);
        EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
        EGAttachObject(go, ( EObjectP ) this);
        EMAddGraphicsToModel(ESIModel(), go);
    }
}


void
QDKTPlate :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    WCRec p [ 4 ];
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
        p [ 3 ].x = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(1, tStep, defScale);
        p [ 3 ].y = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(2, tStep, defScale);
        p [ 3 ].z = ( FPNum ) this->giveNode(4)->giveUpdatedCoordinate(3, tStep, defScale);

        go =  CreateQuad3D(p);
        EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
        EMAddGraphicsToModel(ESIModel(), go);
    }
}


void
QDKTPlate :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int i, indx, result = 0;
    WCRec p [ 4 ];
    GraphicObj *tr;
    FloatArray v [ 4 ];
    double s [ 4 ], defScale;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);
    if ( gc.giveIntVarMode() == ISM_recovered ) {
        for ( i = 1; i <= 4; i++ ) {
            result += this->giveInternalStateAtNode(v [ i - 1 ], gc.giveIntVarType(), gc.giveIntVarMode(), i, tStep);
        }

        if ( result != 4 ) {
            return;
        }

        indx = gc.giveIntVarIndx();

        for ( i = 1; i <= 4; i++ ) {
            s [ i - 1 ] = v [ i - 1 ].at(indx);
        }

        if ( gc.getScalarAlgo() == SA_ISO_SURF ) {
            for ( i = 0; i < 4; i++ ) {
                if ( gc.getInternalVarsDefGeoFlag() ) {
                    // use deformed geometry
                    defScale = gc.getDefScale();
                    p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                    p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                    p [ i ].z = 0.;
                } else {
                    p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                    p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                    p [ i ].z = 0.;
                }
            }

            //EASValsSetColor(gc.getYieldPlotColor(ratio));
            gc.updateFringeTableMinMax(s, 4);
            tr =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
            EGWithMaskChangeAttributes(LAYER_MASK, tr);
            EMAddGraphicsToModel(ESIModel(), tr);

        } else if ( ( gc.getScalarAlgo() == SA_ZPROFILE ) || ( gc.getScalarAlgo() == SA_COLORZPROFILE ) ) {
            double landScale = gc.getLandScale();

            for ( i = 0; i < 4; i++ ) {
                if ( gc.getInternalVarsDefGeoFlag() ) {
                    // use deformed geometry
                    defScale = gc.getDefScale();
                    p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                    p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                    p [ i ].z = s [ i ] * landScale;
                } else {
                    p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                    p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                    p [ i ].z = s [ i ] * landScale;
                }

                // this fixes a bug in ELIXIR
                if ( fabs(s [ i ]) < 1.0e-6 ) {
                    s [ i ] = 1.0e-6;
                }
            }

            if ( gc.getScalarAlgo() == SA_ZPROFILE ) {
                EASValsSetColor( gc.getDeformedElementColor() );
                EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
                tr =  CreateQuad3D(p);
                EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, tr);
            } else {
                gc.updateFringeTableMinMax(s, 4);
                tr =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
                EGWithMaskChangeAttributes(LAYER_MASK, tr);
            }

            EMAddGraphicsToModel(ESIModel(), tr);
        }
    } else if ( gc.giveIntVarMode() == ISM_local ) {
        if ( integrationRulesArray [ 0 ]->giveNumberOfIntegrationPoints() != 4 ) {
            return;
        }

        IntArray ind(4);
        WCRec pp [ 9 ];

        for ( i = 0; i < 4; i++ ) {
            if ( gc.getInternalVarsDefGeoFlag() ) {
                // use deformed geometry
                defScale = gc.getDefScale();
                pp [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(1, tStep, defScale);
                pp [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveUpdatedCoordinate(2, tStep, defScale);
                pp [ i ].z = 0.;
            } else {
                pp [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
                pp [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
                pp [ i ].z = 0.;
            }
        }

        for ( i = 0; i < 3; i++ ) {
            pp [ i + 4 ].x = 0.5 * ( pp [ i ].x + pp [ i + 1 ].x );
            pp [ i + 4 ].y = 0.5 * ( pp [ i ].y + pp [ i + 1 ].y );
            pp [ i + 4 ].z = 0.5 * ( pp [ i ].z + pp [ i + 1 ].z );
        }

        pp [ 7 ].x = 0.5 * ( pp [ 3 ].x + pp [ 0 ].x );
        pp [ 7 ].y = 0.5 * ( pp [ 3 ].y + pp [ 0 ].y );
        pp [ 7 ].z = 0.5 * ( pp [ 3 ].z + pp [ 0 ].z );

        pp [ 8 ].x = 0.25 * ( pp [ 0 ].x + pp [ 1 ].x + pp [ 2 ].x + pp [ 3 ].x );
        pp [ 8 ].y = 0.25 * ( pp [ 0 ].y + pp [ 1 ].y + pp [ 2 ].y + pp [ 3 ].y );
        pp [ 8 ].z = 0.25 * ( pp [ 0 ].z + pp [ 1 ].z + pp [ 2 ].z + pp [ 3 ].z );

        for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {
            const FloatArray &gpCoords = gp->giveNaturalCoordinates();
            if ( ( gpCoords.at(1) > 0. ) && ( gpCoords.at(2) > 0. ) ) {
                ind.at(1) = 0;
                ind.at(2) = 4;
                ind.at(3) = 8;
                ind.at(4) = 7;
            } else if ( ( gpCoords.at(1) < 0. ) && ( gpCoords.at(2) > 0. ) ) {
                ind.at(1) = 4;
                ind.at(2) = 1;
                ind.at(3) = 5;
                ind.at(4) = 8;
            } else if ( ( gpCoords.at(1) < 0. ) && ( gpCoords.at(2) < 0. ) ) {
                ind.at(1) = 5;
                ind.at(2) = 2;
                ind.at(3) = 6;
                ind.at(4) = 8;
            } else {
                ind.at(1) = 6;
                ind.at(2) = 3;
                ind.at(3) = 7;
                ind.at(4) = 8;
            }

            if ( giveIPValue(v [ 0 ], gp, gc.giveIntVarType(), tStep) == 0 ) {
                return;
            }

            indx = gc.giveIntVarIndx();

            for ( i = 1; i <= 4; i++ ) {
                s [ i - 1 ] = v [ 0 ].at(indx);
            }

            for ( i = 0; i < 4; i++ ) {
                p [ i ].x = pp [ ind.at(i + 1) ].x;
                p [ i ].y = pp [ ind.at(i + 1) ].y;
                p [ i ].z = pp [ ind.at(i + 1) ].z;
            }

            gc.updateFringeTableMinMax(s, 4);
            tr =  CreateQuadWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ], s [ 3 ]);
            EGWithMaskChangeAttributes(LAYER_MASK, tr);
            EMAddGraphicsToModel(ESIModel(), tr);
        }
    }
}

#endif
} // end namespace oofem
