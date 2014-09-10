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

#include "trplanrot.h"

#include "node.h"
#include "material.h"
#include "structuralcrosssection.h"
#include "structuralms.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "verbose.h"
#include "load.h"
#include "mathfem.h"
#include "classfactory.h"

namespace oofem {
REGISTER_Element(TrPlaneStrRot);

TrPlaneStrRot :: TrPlaneStrRot(int n, Domain *aDomain) :
    TrPlaneStress2d(n, aDomain)
{
    numberOfDofMans        = 3;
    numberOfGaussPoints    = 4;
    numberOfRotGaussPoints = 1;
}


void
TrPlaneStrRot :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(2);
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], numberOfGaussPoints, this);

        integrationRulesArray [ 1 ] = new GaussIntegrationRule(2, this, 4, 4);
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 1 ], numberOfRotGaussPoints, this);
    }
}


void
TrPlaneStrRot :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns part of strain-displacement matrix {B} of the receiver,
// (epsilon_x,epsilon_y,gamma_xy) = B . r
// type of this part is [3,9]  r=(u1,w1,fi1,u2,w2,fi2,u3,w3,fi3)
// evaluated at gp.
{
    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    //
    FloatArray b(3), c(3);

    for ( int i = 1; i <= 3; i++ ) {
        int j = i + 1 - i / 3 * 3;
        int k = j + 1 - j / 3 * 3;
        b.at(i) = y.at(j) - y.at(k);
        c.at(i) = x.at(k) - x.at(j);
    }

    //
    double area;
    area = this->giveArea();

    //
    int size;
    if ( ui == ALL_STRAINS ) {
        size = 4;
        ui = 4;
    } else {
        size = ui - li + 1;
    }

    if ( ( size < 0 ) || ( size > 4 ) ) {
        OOFEM_ERROR("ComputeBmatrixAt size mismatch");
    }

    //
    int ind = 1;

    answer.resize(size, 9);

    if ( ( li <= 2 ) ) {
        FloatArray nx = this->GiveDerivativeUX(gp);
        FloatArray ny = this->GiveDerivativeVY(gp);

        if ( ( li <= 1 ) && ( ui >= 1 ) ) {
            for ( int i = 1; i <= 3; i++ ) {
                answer.at(ind, 3 * i - 2) = b.at(i) * 1. / ( 2. * area );
                answer.at(ind, 3 * i - 0) = nx.at(i) * 1. / ( 2. * area );
            }

            ind++;
        }

        if ( ( li <= 2 ) && ( ui >= 2 ) ) {
            for ( int i = 1; i <= 3; i++ ) {
                answer.at(ind, 3 * i - 1) = c.at(i) * 1. / ( 2. * area );
                answer.at(ind, 3 * i - 0) = ny.at(i) * 1. / ( 2. * area );
            }

            ind++;
        }
    }

    if ( ( li <= 3 ) && ( ui >= 3 ) ) {
        GaussIntegrationRule ir(1, this, 1, 3);
        ir.SetUpPointsOnTriangle(1, _PlaneStress);

        FloatArray nx = this->GiveDerivativeVX( ir.getIntegrationPoint(0) );
        FloatArray ny = this->GiveDerivativeUY( ir.getIntegrationPoint(0) );

        for ( int i = 1; i <= 3; i++ ) {
            answer.at(ind, 3 * i - 2) = c.at(i) * 1. / ( 2. * area );
            answer.at(ind, 3 * i - 1) = b.at(i) * 1. / ( 2. * area );
            answer.at(ind, 3 * i - 0) = ( nx.at(i) + ny.at(i) ) * 1. / ( 2. * area );
        }

        ind++;
    }

    if ( ( li <= 4 ) && ( ui >= 4 ) ) {
        FloatArray shapeFunct(3);

        FloatArray nx = this->GiveDerivativeVX(gp);
        FloatArray ny = this->GiveDerivativeUY(gp);

        shapeFunct.at(1) = gp->giveNaturalCoordinate(1);
        shapeFunct.at(2) = gp->giveNaturalCoordinate(2);
        shapeFunct.at(3) = 1.0 - shapeFunct.at(1) - shapeFunct.at(2);

        for ( int i = 1; i <= 3; i++ ) {
            answer.at(ind, 3 * i - 2) = -1. * c.at(i) * 1.0 / 4.0 / area;
            answer.at(ind, 3 * i - 1) = b.at(i) * 1.0 / 4.0 / area;
            answer.at(ind, 3 * i - 0) = ( -4. * area * shapeFunct.at(i) + nx.at(i) - ny.at(i) ) * 1.0 / 4.0 / area;
        }
    }
}


void
TrPlaneStrRot :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver
// evaluated at gp.
{
    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    //
    FloatArray b(3), c(3), d(3);

    for ( int i = 1; i <= 3; i++ ) {
        int j = i + 1 - i / 3 * 3;
        int k = j + 1 - j / 3 * 3;
        b.at(i) = y.at(j) - y.at(k);
        c.at(i) = x.at(k) - x.at(j);
        d.at(i) = sqrt( b.at(i) * b.at(i) + c.at(i) * c.at(i) );
    }

    //
    FloatArray angles = this->GivePitch();

    //
    double l1, l2, l3, f11, f12, f13, f21, f22, f23;

    l1 = iLocCoord.at(1);
    l2 = iLocCoord.at(2);
    l3 = 1. - l1 - l2;

    f11 = d.at(2) / 2. *l1 *l3 *sin( angles.at(2) ) - d.at(3) / 2. *l2 *l1 *sin( angles.at(3) );
    f12 = d.at(3) / 2. *l2 *l1 *sin( angles.at(3) ) - d.at(1) / 2. *l3 *l2 *sin( angles.at(1) );
    f13 = d.at(1) / 2. *l3 *l2 *sin( angles.at(1) ) - d.at(2) / 2. *l1 *l3 *sin( angles.at(2) );

    f21 = d.at(3) / 2. *l2 *l1 *cos( angles.at(3) ) - d.at(2) / 2. *l1 *l3 *cos( angles.at(2) );
    f22 = d.at(1) / 2. *l3 *l2 *cos( angles.at(1) ) - d.at(3) / 2. *l2 *l1 *cos( angles.at(3) );
    f23 = d.at(2) / 2. *l1 *l3 *cos( angles.at(2) ) - d.at(1) / 2. *l3 *l2 *cos( angles.at(1) );

    //
    answer.resize(3, 9);
    answer.zero();

    answer.at(1, 1) = l1;
    answer.at(1, 3) = f11;
    answer.at(1, 4) = l2;
    answer.at(1, 6) = f12;
    answer.at(1, 7) = l3;
    answer.at(1, 9) = f13;

    answer.at(2, 2) = l1;
    answer.at(2, 3) = f21;
    answer.at(2, 5) = l2;
    answer.at(2, 6) = f22;
    answer.at(2, 8) = l3;
    answer.at(2, 9) = f23;

    answer.at(3, 3) = l1;
    answer.at(3, 6) = l2;
    answer.at(3, 9) = l3;
}


double
TrPlaneStrRot :: giveArea()
// returns the area occupied by the receiver
{
    if ( area > 0 ) { // check if previously computed
        return area;
    }

    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    //
    double x1, x2, x3, y1, y2, y3;
    x1 = x.at(1);
    x2 = x.at(2);
    x3 = x.at(3);

    y1 = y.at(1);
    y2 = y.at(2);
    y3 = y.at(3);

    return ( area = fabs( 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 ) ) );
    //    return ( area = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 ) );
}


void
TrPlaneStrRot :: giveNodeCoordinates(FloatArray &x, FloatArray &y)
{
    FloatArray *nc1, *nc2, *nc3;
    nc1 = this->giveNode(1)->giveCoordinates();
    nc2 = this->giveNode(2)->giveCoordinates();
    nc3 = this->giveNode(3)->giveCoordinates();

    x.at(1) = nc1->at(1);
    x.at(2) = nc2->at(1);
    x.at(3) = nc3->at(1);

    y.at(1) = nc1->at(2);
    y.at(2) = nc2->at(2);
    y.at(3) = nc3->at(2);

    //if (z) {
    //  z[0] = nc1->at(3);
    //  z[1] = nc2->at(3);
    //  z[2] = nc3->at(3);
    //}
}


FloatArray
TrPlaneStrRot :: GivePitch()
// Returns angles between each side and global x-axis
{
    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    //
    FloatArray angles(3);

    for ( int i = 0; i < 3; i++ ) {
        int j = i + 1 - i / 2 * 3;
        int k = j + 1 - j / 2 * 3;
        if ( x(k) == x(j) ) {
            if ( y(k) > y(j) ) {
                angles.at(i + 1) = 3.14159265358979 / 2.;
            } else {
                angles.at(i + 1) = 3.14159265358979 * 3. / 2.;
            }
        }

        if ( x(k) > x(j) ) {
            if ( y(k) >= y(j) ) {
                angles.at(i + 1) = atan( ( y(k) - y(j) ) / ( x(k) - x(j) ) );
            } else {
                angles.at(i + 1) = 2. * 3.14159265358979 - atan( ( y(j) - y(k) ) / ( x(k) - x(j) ) );
            }
        }

        if ( x(k) < x(j) ) {
            if ( y(k) >= y(j) ) {
                angles.at(i + 1) = 3.14159265358979 - atan( ( y(k) - y(j) ) / ( x(j) - x(k) ) );
            } else {
                angles.at(i + 1) = 3.14159265358979 + atan( ( y(j) - y(k) ) / ( x(j) - x(k) ) );
            }
        }
    }

    return angles;
}


FloatArray
TrPlaneStrRot :: GiveDerivativeUX(GaussPoint *gp)
{
    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    //
    FloatArray b(3), c(3), d(3);

    for ( int i = 1; i <= 3; i++ ) {
        int j = i + 1 - i / 3 * 3;
        int k = j + 1 - j / 3 * 3;
        b.at(i) = y.at(j) - y.at(k);
        c.at(i) = x.at(k) - x.at(j);
        d.at(i) = sqrt( b.at(i) * b.at(i) + c.at(i) * c.at(i) );
    }

    //
    FloatArray angles = this->GivePitch();

    //
    FloatArray shapeFunct(3);
    shapeFunct.at(1) = gp->giveNaturalCoordinate(1);
    shapeFunct.at(2) = gp->giveNaturalCoordinate(2);
    shapeFunct.at(3) = 1.0 - shapeFunct.at(1) - shapeFunct.at(2);

    //
    FloatArray nx(3);

    for ( int i = 1; i <= 3; i++ ) {
        int j = i + 1 - i / 3 * 3;
        int k = j + 1 - j / 3 * 3;
        nx.at(i) = ( d.at(j) / 2. * ( b.at(k) * shapeFunct.at(i) + shapeFunct.at(k) * b.at(i) ) * sin( angles.at(j) ) -
                     d.at(k) / 2. * ( b.at(i) * shapeFunct.at(j) + shapeFunct.at(i) * b.at(j) ) * sin( angles.at(k) ) );
    }

    return nx;
}


FloatArray
TrPlaneStrRot :: GiveDerivativeVX(GaussPoint *gp)
{
    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    //
    FloatArray b(3), c(3), d(3);

    for ( int i = 1; i <= 3; i++ ) {
        int j = i + 1 - i / 3 * 3;
        int k = j + 1 - j / 3 * 3;
        b.at(i) = y.at(j) - y.at(k);
        c.at(i) = x.at(k) - x.at(j);
        d.at(i) = sqrt( b.at(i) * b.at(i) + c.at(i) * c.at(i) );
    }

    //
    FloatArray angles = this->GivePitch();

    //
    FloatArray shapeFunct(3);
    shapeFunct.at(1) = gp->giveNaturalCoordinate(1);
    shapeFunct.at(2) = gp->giveNaturalCoordinate(2);
    shapeFunct.at(3) = 1.0 - shapeFunct.at(1) - shapeFunct.at(2);

    //
    FloatArray nx(3);

    for ( int i = 1; i <= 3; i++ ) {
        int j = i + 1 - i / 3 * 3;
        int k = j + 1 - j / 3 * 3;
        nx.at(i) = ( d.at(j) / 2. * ( b.at(k) * shapeFunct.at(i) + shapeFunct.at(k) * b.at(i) ) * cos( angles.at(j) ) -
                     d.at(k) / 2. * ( b.at(i) * shapeFunct.at(j) + shapeFunct.at(i) * b.at(j) ) * cos( angles.at(k) ) ) * ( -1.0 );
    }

    return nx;
}


FloatArray
TrPlaneStrRot :: GiveDerivativeUY(GaussPoint *gp)
{
    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    //
    FloatArray b(3), c(3), d(3);

    for ( int i = 1; i <= 3; i++ ) {
        int j = i + 1 - i / 3 * 3;
        int k = j + 1 - j / 3 * 3;
        b.at(i) = y.at(j) - y.at(k);
        c.at(i) = x.at(k) - x.at(j);
        d.at(i) = sqrt( b.at(i) * b.at(i) + c.at(i) * c.at(i) );
    }

    //
    FloatArray angles = this->GivePitch();

    //
    FloatArray shapeFunct(3);
    shapeFunct.at(1) = gp->giveNaturalCoordinate(1);
    shapeFunct.at(2) = gp->giveNaturalCoordinate(2);
    shapeFunct.at(3) = 1.0 - shapeFunct.at(1) - shapeFunct.at(2);

    //
    FloatArray ny(3);

    for ( int i = 1; i <= 3; i++ ) {
        int j = i + 1 - i / 3 * 3;
        int k = j + 1 - j / 3 * 3;
        ny.at(i) = ( d.at(j) / 2. * ( c.at(k) * shapeFunct.at(i) + shapeFunct.at(k) * c.at(i) ) * sin( angles.at(j) ) -
                     d.at(k) / 2. * ( c.at(i) * shapeFunct.at(j) + shapeFunct.at(i) * c.at(j) ) * sin( angles.at(k) ) );
    }

    return ny;
}


FloatArray
TrPlaneStrRot :: GiveDerivativeVY(GaussPoint *gp)
{
    // get node coordinates
    FloatArray x(3), y(3);
    this->giveNodeCoordinates(x, y);

    //
    FloatArray b(3), c(3), d(3);

    for ( int i = 1; i <= 3; i++ ) {
        int j = i + 1 - i / 3 * 3;
        int k = j + 1 - j / 3 * 3;
        b.at(i) = y.at(j) - y.at(k);
        c.at(i) = x.at(k) - x.at(j);
        d.at(i) = sqrt( b.at(i) * b.at(i) + c.at(i) * c.at(i) );
    }

    //
    FloatArray angles = this->GivePitch();

    //
    FloatArray shapeFunct(3);
    shapeFunct.at(1) = gp->giveNaturalCoordinate(1);
    shapeFunct.at(2) = gp->giveNaturalCoordinate(2);
    shapeFunct.at(3) = 1.0 - shapeFunct.at(1) - shapeFunct.at(2);

    //
    FloatArray ny(3);

    for ( int i = 1; i <= 3; i++ ) {
        int j = i + 1 - i / 3 * 3;
        int k = j + 1 - j / 3 * 3;
        ny.at(i) = ( d.at(j) / 2. * ( c.at(k) * shapeFunct.at(i) + shapeFunct.at(k) * c.at(i) ) * cos( angles.at(j) ) -
                     d.at(k) / 2. * ( c.at(i) * shapeFunct.at(j) + shapeFunct.at(i) * c.at(j) ) * cos( angles.at(k) ) ) * ( -1.0 );
    }

    return ny;
}


IRResultType
TrPlaneStrRot :: initializeFrom(InputRecord *ir)
{
    IRResultType result;              // Required by IR_GIVE_FIELD macro

    numberOfGaussPoints = 4;
    result = this->StructuralElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }

    numberOfRotGaussPoints = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfRotGaussPoints, _IFT_TrPlaneStrRot_niprot);

    if ( !( ( numberOfGaussPoints == 1 ) ||
           ( numberOfGaussPoints == 4 ) ||
           ( numberOfGaussPoints == 7 ) ) ) {
        numberOfGaussPoints = 4;
    }

    if ( !( ( numberOfRotGaussPoints == 1 ) ||
           ( numberOfRotGaussPoints == 4 ) ||
           ( numberOfRotGaussPoints == 7 ) ) ) {
        numberOfRotGaussPoints = 1;
    }

    return IRRT_OK;
}

#if 0
void
TrPlaneStrRot :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    FloatMatrix B, d;
    FloatArray vStress, vStrain, u;

    // This function can be quite costly to do inside the loops when one has many slave dofs.
    this->computeVectorOf(VM_Total, tStep, u);
    // subtract initial displacements, if defined
    if ( initialDisplacements ) {
        u.subtract(* initialDisplacements);
    }

    // zero answer will resize accordingly when adding first contribution
    answer.clear();

    for ( GaussPoint *gp: *this->giveDefaultIntegrationRulePtr() ) {

        // Engineering (small strain) stress
        if ( nlGeometry == 0 ) {
            this->computeBmatrixAt(gp, B);

            if ( !this->isActivated(tStep) ) {
                vStrain.resize( StructuralMaterial :: giveSizeOfVoigtSymVector( gp->giveMaterialMode() ) );
                vStrain.zero();
            }
            vStrain.beProductOf(B, u);
 #if 1
            this->computeStressVector(vStress, vStrain, gp, tStep);
 #else
            if ( useUpdatedGpRecord == 1 && false ) {
                ///@todo This is problematic, it only saves the plane stress components. Should we just keep the curvature field completely separate,
                /// or introduce the membrane+rotation mode directly into the cross-section.
                vStress = static_cast< StructuralMaterialStatus * >( gp->giveStatus() )->giveStressVector();
                /*
                 *              // Curvature field:
                 *              vStress.resizeWithValues(4);
                 *              //cs->givePlaneStressStiffMtrx(d, TangentStiffness, gp, tStep);
                 *              this->giveStructuralCrossSection()->giveCharMaterialStiffnessMatrix(d, ElasticStiffness, gp, tStep);
                 *              vStress.resizeWithValues(4, 0);
                 *              vStress.at(4) = vStrain.at(4) * d.at(3,3);
                 */
            } else {
                this->computeStressVector(vStress, vStrain, gp, tStep);
            }
 #endif
        } else if ( nlGeometry == 1 ) {  // First Piola-Kirchhoff stress
            OOFEM_ERROR("Only small strain mode is supported");
        }

        // Compute nodal internal forces at nodes as f = B^T*Stress dV
        double dV  = this->computeVolumeAround(gp);
        answer.plusProduct(B, vStress, dV);
    }

    // If inactive: update fields but do not give any contribution to the internal forces
    if ( !this->isActivated(tStep) ) {
        answer.zero();
        return;
    }
}
#endif


void
TrPlaneStrRot :: computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
// Computes the vector containing the strains at the Gauss point gp of
// the receiver, at time step tStep. The nature of these strains depends
// on the element's type.
{
    FloatMatrix b;
    FloatArray u, Epsilon;

    this->computeVectorOf(VM_Total, tStep, u);

    answer.resize(4);
    answer.zero();

    this->computeBmatrixAt(gp, b, 1, 3);
    Epsilon.beProductOf(b, u);
    answer.at(1) = Epsilon.at(1);
    answer.at(2) = Epsilon.at(2);
    answer.at(3) = Epsilon.at(3);

    if ( numberOfRotGaussPoints == 1 ) {
        //
        // if reduced integration in one gp only
        // force the evaluation of eps_fi in this gauss point
        // instead of evaluating in given gp
        //
        GaussPoint *helpGaussPoint;
        helpGaussPoint = integrationRulesArray [ 1 ]->getIntegrationPoint(0);

        this->computeBmatrixAt(helpGaussPoint, b, 4, 4);
    } else {
        OOFEM_ERROR("numberOfRotGaussPoints size mismatch");
    }

    Epsilon.beProductOf(b, u);
    answer.at(4) = Epsilon.at(1);
}


void
TrPlaneStrRot :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    StructuralCrossSection *cs = this->giveStructuralCrossSection();
    cs->giveMembraneRotStiffMtrx(answer, rMode, gp, tStep);
}


void
TrPlaneStrRot :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    StructuralCrossSection *cs = this->giveStructuralCrossSection();
    cs->giveGeneralizedStress_MembraneRot(answer, gp, strain, tStep);
}


void
TrPlaneStrRot :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {D_u, D_v, R_w};
}


void
TrPlaneStrRot :: computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode)
// Computes numerically the load vector of the receiver due to the body loads, at tStep.
// load is assumed to be in global cs.
// load vector is then transformed to coordinate system in each node.
// (should be global coordinate system, but there may be defined
//  different coordinate system in each node)
{
    double dens, dV, load;
    FloatArray force;
    FloatMatrix T;

    if ( ( forLoad->giveBCGeoType() != BodyLoadBGT ) || ( forLoad->giveBCValType() != ForceLoadBVT ) ) {
        OOFEM_ERROR("unknown load type");
    }

    // note: force is assumed to be in global coordinate system.
    forLoad->computeComponentArrayAt(force, tStep, mode);

    if ( force.giveSize() ) {
        GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

        dens = this->giveStructuralCrossSection()->give('d', gp);
        dV   = this->computeVolumeAround(gp) * this->giveCrossSection()->give(CS_Thickness, gp);

        answer.resize(9);
        answer.zero();

        load = force.at(1) * dens * dV / 3.0;
        answer.at(1) = load;
        answer.at(4) = load;
        answer.at(7) = load;

        load = force.at(2) * dens * dV / 3.0;
        answer.at(2) = load;
        answer.at(5) = load;
        answer.at(8) = load;

        // transform result from global cs to local element cs.
        if ( this->computeGtoLRotationMatrix(T) ) {
            answer.rotatedWith(T, 'n');
        }
    } else {
        answer.clear();          // nil resultant
    }
}


double
TrPlaneStrRot :: giveCharacteristicLength(const FloatArray &normalToCrackPlane)
//
// returns receiver's characteristic length for crack band models
// for a crack formed in the plane with normal normalToCrackPlane.
//
{
    return this->giveCharacteristicLengthForPlaneElements(normalToCrackPlane);
}

} // end namespace oofem
