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

#include "beam3d.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "flotarry.h"
#include "engngm.h"
#include "boundaryload.h"
#include "mathfem.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
Beam3d :: Beam3d(int n, Domain *aDomain) : StructuralElement(n, aDomain)
{
    numberOfDofMans = 2;
    referenceNode = 0;

    length = 0.;
    kappay = kappaz = -1.0;
    dofsToCondense = NULL;
}

Beam3d :: ~Beam3d()
{
    delete dofsToCondense;
}


void
Beam3d :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the strain matrix of the receiver.
// eeps = {\eps_x, \gamma_xz, \gamma_xy, \der{phi_x}{x}, \kappa_y, \kappa_z}^T
{
    double l, ksi, kappay, kappaz;

    l     = this->giveLength();
    ksi   = 0.5 + 0.5 * aGaussPoint->giveCoordinate(1);
    kappay = this->giveKappayCoeff();
    kappaz = this->giveKappazCoeff();

    answer.resize(6, 12);
    answer.zero();

    answer.at(1, 1) =  -1. / l;
    answer.at(1, 7) =   1. / l;
    answer.at(2, 2) =   ( -2. * kappaz ) / ( l * ( 1. + 2. * kappaz ) );
    answer.at(2, 6) =   kappaz / ( l * ( 1. + 2. * kappaz ) );
    answer.at(2, 8) =   2. * kappaz / ( l * ( 1. + 2. * kappaz ) );
    answer.at(2, 12) =   kappaz / ( l * ( 1. + 2. * kappaz ) );
    answer.at(3, 3) =   ( -2. * kappay ) / ( l * ( 1. + 2. * kappay ) );
    answer.at(3, 5) =   kappay / ( l * ( 1. + 2. * kappay ) );
    answer.at(3, 9) =   2. * kappay / ( l * ( 1. + 2. * kappay ) );
    answer.at(3, 11) =   kappay / ( l * ( 1. + 2. * kappay ) );

    answer.at(4, 4) =  -1. / l;
    answer.at(4, 10) =   1. / l;
    answer.at(5, 3) =   ( 6. - 12. * ksi ) / ( l * l * ( 1. + 2. * kappay ) );
    answer.at(5, 5) =   ( -2. * ( 2. + kappay ) + 6. * ksi ) / ( l * ( 1. + 2. * kappay ) );
    answer.at(5, 9) =   ( -6. + 12. * ksi ) / ( l * l * ( 1. + 2. * kappay ) );
    answer.at(5, 11) =   ( -2. * ( 1. - kappay ) + 6. * ksi ) / ( l * ( 1. + 2. * kappay ) );
    answer.at(6, 2) =   ( 6. - 12. * ksi ) / ( l * l * ( 1. + 2. * kappaz ) );
    answer.at(6, 6) =   ( -2. * ( 2. + kappaz ) + 6. * ksi ) / ( l * ( 1. + 2. * kappaz ) );
    answer.at(6, 8) =   ( -6. + 12. * ksi ) / ( l * l * ( 1. + 2. * kappaz ) );
    answer.at(6, 12) =   ( -2. * ( 1. - kappaz ) + 6. * ksi ) / ( l * ( 1. + 2. * kappaz ) );
}

void Beam3d :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( !integrationRulesArray ) {
        // the gauss point is used only when methods from crosssection and/or material
        // classes are requested
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 2);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Line, 3, _3dBeam);
    }
}


void
Beam3d :: computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint. Used for numerical calculation of consistent mass
// matrix. Must contain only interpolation for displacement terms,
// not for any rotations. (Inertia forces do not work on rotations).
// r = {u1,v1,w1,fi_x1,fi_y1,fi_z1,u2,v2,w2,fi_x2,fi_y2,fi_21}^T
{
    double l, ksi, ksi2, ksi3, kappay, kappaz, c1y, c1z;

    l     = this->giveLength();
    ksi =   0.5 + 0.5 * aGaussPoint->giveCoordinate(1);
    kappay = this->giveKappayCoeff();
    kappaz = this->giveKappazCoeff();
    c1y = 1. + 2. * kappay;
    c1z = 1. + 2. * kappaz;
    ksi2 = ksi * ksi;
    ksi3 = ksi2 * ksi;

    answer.resize(6, 12);
    answer.zero();

    answer.at(1, 1) = 1. - ksi;
    answer.at(1, 7) = ksi;
    answer.at(2, 2) = ( ( 1. + 2. * kappaz ) - 2. * kappaz * ksi - 3. * ksi2 + 2. * ksi3 ) / c1z;
    answer.at(2, 6) = -l * ( -( 1. + kappaz ) * ksi + ( 2. + kappaz ) * ksi2 - ksi3 ) / c1z;
    answer.at(2, 8) = ( 2. * kappaz * ksi + 3. * ksi2 - 2. * ksi3 ) / c1z;
    answer.at(2, 12) = -l * ( kappaz * ksi + ( 1. - kappaz ) * ksi2 - ksi3 ) / c1z;
    answer.at(3, 3) = ( ( 1. + 2. * kappay ) - 2. * kappay * ksi - 3. * ksi2 + 2. * ksi3 ) / c1y;
    answer.at(3, 5) = l * ( -( 1. + kappay ) * ksi + ( 2. + kappay ) * ksi2 - ksi3 ) / c1y;
    answer.at(3, 9) = ( 2. * kappay * ksi + 3. * ksi2 - 2. * ksi3 ) / c1y;
    answer.at(3, 11) = l * ( kappay * ksi + ( 1. - kappay ) * ksi2 - ksi3 ) / c1y;

    // rotations excluded
    answer.at(4, 4) = 1. - ksi;
    answer.at(4, 10) = ksi;
    answer.at(5, 3) = ( 6. * ksi - 6. * ksi2 ) / ( l * c1y );
    answer.at(5, 5) = ( ( 1. + 2. * kappay ) - 2. * ( 2. + kappay ) * ksi + 3. * ksi2 ) / c1y;
    answer.at(5, 9) = -( 6. * ksi + 6. * ksi2 ) / ( l * c1y );
    answer.at(5, 11) = ( -2. * ( 1. - kappay ) * ksi + 3. * ksi2 ) / c1y;
    answer.at(6, 2) = -( 6. * ksi - 6. * ksi2 ) / ( l * c1z );
    answer.at(6, 6) = ( ( 1. + 2. * kappaz ) - 2. * ( 2. + kappaz ) * ksi + 3. * ksi2 ) / c1z;
    answer.at(6, 8) = ( 6. * ksi + 6. * ksi2 ) / ( l * c1z );
    answer.at(6, 12) = ( -2. * ( 1. - kappaz ) * ksi + 3. * ksi2 ) / c1z;
}

void
Beam3d :: computeLocalStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
// Returns the stiffness matrix of the receiver, expressed in the global
// axes. No integration over volume done, beam with constant material and crosssection
// parameters assumed.
{
    // compute clamped stifness
    this->computeClampedStiffnessMatrix(answer, rMode, tStep);

    // condense requested dofs
    if ( dofsToCondense ) {
        this->condense(& answer, NULL, NULL, dofsToCondense);
    }
}


void
Beam3d :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
// Returns the stiffness matrix of the receiver, expressed in the global
// axes. No integration over volume done, beam with constant material and crosssection
// parameters assumed.
{
    // compute clamped stiffness
    this->computeLocalStiffnessMatrix(answer, rMode, tStep);
}


void
Beam3d :: computeClampedStiffnessMatrix(FloatMatrix &answer,
                                        MatResponseMode rMode, TimeStep *tStep)
// Returns the stiffness matrix of the receiver, expressed in the local
// axes. No integration over volume done, beam with constant material and crosssection
// parameters assumed.
{
    double l, l2, l3, eiy, eiz, kappay, kappaz, c1y, c1z;
    FloatMatrix d;

    this->computeConstitutiveMatrixAt(d, rMode, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);
    l     = this->giveLength();
    l2    = l * l;
    l3    = l2 * l;
    kappay = this->giveKappayCoeff();
    kappaz = this->giveKappazCoeff();
    c1y = 1. + 2. * kappay;
    c1z = 1. + 2. * kappaz;
    eiy = d.at(5, 5);
    eiz = d.at(6, 6);

    answer.resize(12, 12);
    answer.zero();

    answer.at(1, 1) =  d.at(1, 1) / l;
    answer.at(1, 7) = -d.at(1, 1) / l;
    answer.at(2, 2) =  eiz * 12. / ( l3 * c1z );
    answer.at(2, 6) =  eiz * 6.  / ( l2 * c1z );
    answer.at(2, 8) = -eiz * 12. / ( l3 * c1z );
    answer.at(2, 12) =  eiz * 6.  / ( l2 * c1z );
    answer.at(3, 3) =  eiy * 12. / ( l3 * c1y );
    answer.at(3, 5) = -eiy * 6.  / ( l2 * c1y );
    answer.at(3, 9) = -eiy * 12. / ( l3 * c1y );
    answer.at(3, 11) = -eiy * 6.  / ( l2 * c1y );

    answer.at(4, 4) =  d.at(4, 4) / l;
    answer.at(4, 10) = -d.at(4, 4) / l;
    answer.at(5, 5) =  eiy * 2. * ( 2. + kappay ) / ( l * c1y );
    answer.at(5, 9) =  eiy * 6.  / ( l2 * c1y );
    answer.at(5, 11) =  eiy * 2. * ( 1. - kappay ) / ( l * c1y );
    answer.at(6, 6) =  eiz * 2. * ( 2. + kappaz ) / ( l * c1z );
    answer.at(6, 8) = -eiz * 6.  / ( l2 * c1z );
    answer.at(6, 12) =  eiz * 2. * ( 1. - kappaz ) / ( l * c1z );

    answer.at(7, 7) =  d.at(1, 1) / l;
    answer.at(8, 8) =  eiz * 12. / ( l3 * c1z );
    answer.at(8, 12) = -eiz * 6.  / ( l2 * c1z );
    answer.at(9, 9) =  eiy * 12. / ( l3 * c1y );
    answer.at(9, 11) =  eiy * 6.  / ( l2 * c1y );
    answer.at(10, 10) = d.at(4, 4) / l;
    answer.at(11, 11) = eiy * 2. * ( 2. + kappay ) / ( l * c1y );
    answer.at(12, 12) = eiz * 2. * ( 2. + kappaz ) / ( l * c1z );

    answer.symmetrized();  // symmetrize answer
}


int
Beam3d :: computeLoadGToLRotationMtrx(FloatMatrix &answer)
{
    FloatMatrix lcs;
    int i, j;

    answer.resize(6, 6);
    answer.zero();

    this->giveLocalCoordinateSystem(lcs);
    for ( i = 1; i <= 3; i++ ) {
        for ( j = 1; j <= 3; j++ ) {
            answer.at(i, j) = lcs.at(i, j);
            answer.at(3 + i, 3 + j) = lcs.at(i, j);
        }
    }
    return 1;
}


bool
Beam3d :: computeGtoLRotationMatrix(FloatMatrix &answer)
// Returns the rotation matrix of the receiver.
{
    FloatMatrix lcs;
    int i, j;

    answer.resize(12, 12);
    answer.zero();

    this->giveLocalCoordinateSystem(lcs);
    for ( i = 1; i <= 3; i++ ) {
        for ( j = 1; j <= 3; j++ ) {
            answer.at(i, j) = lcs.at(i, j);
            answer.at(i + 3, j + 3) = lcs.at(i, j);
            answer.at(i + 6, j + 6) = lcs.at(i, j);
            answer.at(i + 9, j + 9) = lcs.at(i, j);
        }
    }
    return 1;
}


double
Beam3d :: computeVolumeAround(GaussPoint *aGaussPoint)
{
    double weight  = aGaussPoint->giveWeight();
    return weight * 0.5 * this->giveLength();
}


void
Beam3d :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    answer.setValues(6, D_u, D_v, D_w, R_u, R_v, R_w);
}


double
Beam3d :: giveLength()
// Returns the length of the receiver.
{
    double dx, dy, dz;
    Node *nodeA, *nodeB;

    if ( length == 0. ) {
        nodeA   = this->giveNode(1);
        nodeB   = this->giveNode(2);
        dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
        dy      = nodeB->giveCoordinate(2) - nodeA->giveCoordinate(2);
        dz      = nodeB->giveCoordinate(3) - nodeA->giveCoordinate(3);
        length  = sqrt(dx * dx + dy * dy + dz * dz);
    }

    return length;
}


void
Beam3d :: computeKappaCoeffs()
{
    // computes kappa coeff
    // kappa_y = (6*E*Iy)/(k*G*A*l^2)

    FloatMatrix d;
    double l = this->giveLength();

    this->computeConstitutiveMatrixAt( d, TangentStiffness, integrationRulesArray [ 0 ]->getIntegrationPoint(0), domain->giveEngngModel()->giveCurrentStep() );

    //  kappay = 6. * d.at(5, 5) / ( d.at(3, 3) * l * l );
    //  kappaz = 6. * d.at(6, 6) / ( d.at(2, 2) * l * l );
    if (d.at(3, 3) != 0.) {
        kappay = 6. * d.at(5, 5) / (d.at(3, 3) * l * l);
    } else {
        kappay = 0.;
    }
    if (d.at(2, 2) != 0.) {
        kappaz = 6. * d.at(6, 6) / (d.at(2, 2) * l * l);
    } else {
        kappaz = 0.;
    }
}


double
Beam3d :: giveKappayCoeff()
{
    if ( kappay < 0.0 ) {
        this->computeKappaCoeffs();
    }

    return kappay;
}


double
Beam3d :: giveKappazCoeff()
{
    if ( kappaz < 0.0 ) {
        this->computeKappaCoeffs();
    }

    return kappaz;
}


int
Beam3d :: giveLocalCoordinateSystem(FloatMatrix &answer)
//
// returns a unit vectors of local coordinate system at element
// stored rowwise (mainly used by some materials with ortho and anisotrophy)
//
{
    FloatArray lx(3), ly(3), lz(3), help(3);
    double length = this->giveLength();
    Node *nodeA, *nodeB, *refNode;
    int i;

    answer.resize(3, 3);
    answer.zero();
    nodeA  = this->giveNode(1);
    nodeB  = this->giveNode(2);
    refNode = this->giveDomain()->giveNode(this->referenceNode);

    for ( i = 1; i <= 3; i++ ) {
        lx.at(i) = ( nodeB->giveCoordinate(i) - nodeA->giveCoordinate(i) ) / length;
        help.at(i) = ( refNode->giveCoordinate(i) - nodeA->giveCoordinate(i) );
    }

    lz.beVectorProductOf(lx, help);
    lz.normalize();
    ly.beVectorProductOf(lz, lx);
    ly.normalize();

    for ( i = 1; i <= 3; i++ ) {
        answer.at(1, i) = lx.at(i);
        answer.at(2, i) = ly.at(i);
        answer.at(3, i) = lz.at(i);
    }

    return 1;
}


IRResultType
Beam3d :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    // first call parent
    StructuralElement :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, referenceNode, IFT_Beam3d_refnode, "refnode"); // Macro
    if ( referenceNode == 0 ) {
        _error("instanciateFrom: wrong reference node specified");
    }

    if ( ir->hasField(IFT_Beam3d_dofstocondense, "dofstocondense") ) {
        IntArray val;
        IR_GIVE_FIELD(ir, val, IFT_Beam3d_dofstocondense, "dofstocondense"); // Macro
        if ( val.giveSize() >= 12 ) {
            _error("instanciateFrom: wrong input data for condensed dofs");
        }

        dofsToCondense = new IntArray(val);
    } else {
        dofsToCondense = NULL;
    }

    this->computeGaussPoints();
    return IRRT_OK;
}


void
Beam3d :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    // stress equivalent vector in nodes (vector of internal forces)
    FloatArray prescStrainEndForces;
    FloatMatrix stiffness;
    FloatArray u;

    this->computeStiffnessMatrix(stiffness, SecantStiffness, tStep);
    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);
    answer.beProductOf(stiffness, u);

    this->computePrescribedStrainLoadVectorAt(prescStrainEndForces, tStep, VM_Total);
    if ( prescStrainEndForces.giveSize() ) {
        answer.subtract(prescStrainEndForces);
    }
}


void
Beam3d :: giveEndForcesVector(FloatArray &answer, TimeStep *tStep)
{
    // computes exact global end-forces vector
    FloatArray loadEndForces;

    this->giveInternalForcesVector(answer, tStep);

    // add exact end forces due to nonnodal loading
    this->computeForceLoadVector(loadEndForces, tStep, VM_Total);
    if ( loadEndForces.giveSize() ) {
        answer.subtract(loadEndForces);
    }
}


void
Beam3d :: computeEdgeLoadVectorAt(FloatArray &answer, Load *load, int iedge, TimeStep *tStep, ValueModeType mode)
{
    FloatArray coords, components, endComponents;
    FloatMatrix T;
    FloatArray floc(12);
    double l = this->giveLength();
    double kappay = this->giveKappayCoeff();
    double kappaz = this->giveKappazCoeff();
    double fx, fy, fz, fmx, fmy, fmz, dfx, dfy, dfz, dfmx, dfmy, dfmz;

    // evaluates the receivers edge load vector
    // for clamped beam
    //
    BoundaryLoad *edgeLoad = dynamic_cast< BoundaryLoad * >(load);
    if ( edgeLoad ) {
        if ( edgeLoad->giveNumberOfDofs() != 6 ) {
            _error("computeEdgeLoadVectorAt: load number of dofs mismatch");
        }

        answer.resize(12);
        answer.zero();

        switch ( edgeLoad->giveClassID() ) {
        case ConstantEdgeLoadClass:

            //edgeLoad->computeComponentArrayAt(components, tStep, mode);
            if ( edgeLoad->giveFormulationType() == BoundaryLoad :: BL_EntityFormulation ) {
                coords.resize(1);
                coords.at(1) = 0.0;
            } else {
                coords = * ( this->giveNode(1)->giveCoordinates() );
            }

            edgeLoad->computeValueAt(components, tStep, coords, mode);


            // prepare transformation coeffs
            if ( edgeLoad->giveCoordSystMode() == BoundaryLoad :: BL_GlobalMode ) {
                if ( this->computeLoadGToLRotationMtrx(T) ) {
                    components.rotatedWith(T, 'n');
                }
            }

            fx = components.at(1);
            fy = components.at(2);
            fz = components.at(3);
            fmx = components.at(4);
            fmy = components.at(5);
            fmz = components.at(6);


            answer.at(1) = fx * l / 2.;
            answer.at(2) = fy * l / 2. + fmz / ( 1. + 2. * kappaz );
            answer.at(3) = fz * l / 2. + fmy / ( 1. + 2. * kappay );
            answer.at(4) = fmx * l / 2.;
            answer.at(5) = ( -1. ) * fz * l * l / 12. + fmy * l * kappay / ( 1. + 2. * kappay );
            answer.at(6) = ( 1. ) * fy * l * l / 12. + fmz * l * kappaz / ( 1. + 2. * kappaz );

            answer.at(7) = fx * l / 2.;
            answer.at(8) = fy * l / 2. - fmz / ( 1. + 2. * kappaz );
            answer.at(9) = fz * l / 2. - fmy / ( 1. + 2. * kappay );
            answer.at(10) = fmx * l / 2.;
            answer.at(11) = ( 1. ) * fz * l * l / 12. + fmy * l * kappay / ( 1. + 2. * kappay );
            answer.at(12) = ( -1. ) * fy * l * l / 12. + fmz * l * kappaz / ( 1. + 2. * kappaz );
            break;
        case LinearEdgeLoadClass:

            /*
             * coords.at(1) = -1.;
             * edgeLoad->computeValueAt(components, tStep, coords, mode);
             */

            if ( edgeLoad->giveFormulationType() == BoundaryLoad :: BL_EntityFormulation ) {
                coords.resize(1);
                coords.at(1) = -1.0;
            } else {
                coords = * ( this->giveNode(1)->giveCoordinates() );
            }

            edgeLoad->computeValueAt(components, tStep, coords, mode);


            // prepare transformation coeffs
            if ( edgeLoad->giveCoordSystMode() == BoundaryLoad :: BL_GlobalMode ) {
                if ( this->computeLoadGToLRotationMtrx(T) ) {
                    components.rotatedWith(T, 'n');
                }
            }

            fx = components.at(1);
            fy = components.at(2);
            fz = components.at(3);
            fmx = components.at(4);
            fmy = components.at(5);
            fmz = components.at(6);

            /*
             * coords.at(1) = 1.;
             * edgeLoad->computeValueAt(endComponents, tStep, coords, mode);
             */
            if ( edgeLoad->giveFormulationType() == BoundaryLoad :: BL_EntityFormulation ) {
                coords.resize(1);
                coords.at(1) = 1.0;
            } else {
                coords = * ( this->giveNode(2)->giveCoordinates() );
            }

            edgeLoad->computeValueAt(endComponents, tStep, coords, mode);

            // prepare transformation coeffs
            if ( edgeLoad->giveCoordSystMode() == BoundaryLoad :: BL_GlobalMode ) {
                if ( T.isNotEmpty() ) {
                    endComponents.rotatedWith(T, 'n');
                }
            }

            // compute differences
            endComponents.subtract(components);

            dfx = endComponents.at(1);
            dfy = endComponents.at(2);
            dfz = endComponents.at(3);
            dfmx = endComponents.at(4);
            dfmy = endComponents.at(5);
            dfmz = endComponents.at(6);


            answer.at(1) = fx * l / 2. + dfx * l / 6.;
            answer.at(2) = fy * l / 2. + dfy * l * ( 20. * kappaz + 9 ) / ( 60. * ( 1. + 2. * kappaz ) ) +
                           fmz / ( 1. + 2. * kappaz ) + dfmz * ( 1. / 2. ) / ( 1. + 2. * kappaz );
            answer.at(3) = fz * l / 2. + dfz * l * ( 20. * kappay + 9 ) / ( 60. * ( 1. + 2. * kappay ) ) +
                           fmy / ( 1. + 2. * kappay ) + dfmy * ( 1. / 2. ) / ( 1. + 2. * kappay );
            answer.at(4) = fmx * l / 2. + dfmx * l / 6.;
            answer.at(5) = ( -1. ) * fz * l * l / 12. - dfz * l * l * ( 5. * kappay + 2. ) / ( 60. * ( 1. + 2. * kappay ) ) +
                           fmy * l * kappay / ( 1. + 2. * kappay ) + dfmy * l * ( 4. * kappay - 1. ) / ( 12. * ( 1. + 2. * kappay ) );
            answer.at(6) = ( 1. ) * fy * l * l / 12. + dfy * l * l * ( 5. * kappaz + 2. ) / ( 60. * ( 1. + 2. * kappaz ) ) +
                           fmz * l * kappaz / ( 1. + 2. * kappaz ) + dfmz * l * ( 4. * kappaz - 1. ) / ( 12. * ( 1. + 2. * kappaz ) );

            answer.at(7) = fx * l / 2. + dfx * l / 3.;
            answer.at(8) = fy * l / 2. + dfy * l * ( 40. * kappaz + 21 ) / ( 60. * ( 1. + 2. * kappaz ) ) -
                           fmz / ( 1. + 2. * kappaz ) - dfmz * ( 1. / 2. ) / ( ( 1. + 2. * kappaz ) );
            answer.at(9) = fz * l / 2. + dfz * l * ( 40. * kappay + 21 ) / ( 60. * ( 1. + 2. * kappay ) ) -
                           fmy / ( 1. + 2. * kappay ) - dfmy * ( 1. / 2. ) / ( ( 1. + 2. * kappay ) );
            answer.at(10) = fmx * l / 2. + dfmx * l / 3.;
            answer.at(11) = ( 1. ) * fz * l * l / 12. + dfz * l * l * ( 5. * kappay + 3. ) / ( 60. * ( 1. + 2. * kappay ) ) +
                            fmy * l * kappay / ( 1. + 2. * kappay ) + dfmy * l * ( 8. * kappay + 1. ) / ( 12. * ( 1. + 2. * kappay ) );
            answer.at(12) = ( -1. ) * fy * l * l / 12. - dfy * l * l * ( 5. * kappaz + 3. ) / ( 60. * ( 1. + 2. * kappaz ) ) +
                            fmz * l * kappaz / ( 1. + 2. * kappaz ) + dfmz * l * ( 8. * kappaz + 1. ) / ( 12. * ( 1. + 2. * kappaz ) );
            break;
        default:
            _error("computeEdgeLoadVectorAt: unsupported load type");
        }
    }
}


void
Beam3d :: printOutputAt(FILE *File, TimeStep *stepN)
{
    // Performs end-of-step operations.

    int i, n;
    FloatArray rl, Fl;
    FloatMatrix T;

    fprintf(File, "beam element %d :\n", number);

    // ask for global element displacement vector
    this->computeVectorOf(EID_MomentumBalance, VM_Total, stepN, rl);
    // ask for global element end forces vector
    this->giveEndForcesVector(Fl, stepN);

    fprintf(File, "  local displacements ");
    n = rl.giveSize();
    for ( i = 1; i <= n; i++ ) {
        fprintf( File, " % .4e", rl.at(i) );
    }

    fprintf(File, "\n  local end forces    ");
    n = Fl.giveSize();
    for ( i = 1; i <= n; i++ ) {
        fprintf( File, " % .4e", Fl.at(i) );
    }

    fprintf(File, "\n");
}


void
Beam3d :: computeLocalForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
// Computes the load vector of the receiver, at stepN.
{
    FloatMatrix stiff;

    StructuralElement :: computeLocalForceLoadVector(answer, stepN, mode); // in global c.s

    if ( answer.giveSize() && dofsToCondense ) {
        // condense requested dofs
        if ( answer.giveSize() != 0 ) {
            this->computeClampedStiffnessMatrix(stiff, TangentStiffness, stepN);
            this->condense(& stiff, NULL, & answer, dofsToCondense);
        }
    }
}


void
Beam3d :: computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode)
{
    StructuralElement::computeBodyLoadVectorAt(answer, load, tStep, mode);
    answer.times(this->giveCrossSection()->give(CS_Area));
}


/*
 * void
 * Beam3d :: computeForceLoadVector (FloatArray& answer, TimeStep* stepN, ValueModeType mode)
 * // Computes the load vector of the receiver, at stepN.
 * {
 * FloatMatrix stiff, *T;
 *
 * StructuralElement::computeForceLoadVector(answer, stepN, mode); // in global c.s
 *
 * if (answer.giveSize() && dofsToCondense) {
 * // condense requested dofs
 * if (answer.giveSize() != 0) {
 * this->computeClampedStiffnessMatrix (stiff, TangentStiffness, stepN) ;
 * this->condense (&stiff, NULL, &answer, dofsToCondense);
 * }
 *
 * }
 * }
 */


void
Beam3d :: computePrescribedStrainLocalLoadVectorAt(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
// Computes the load vector of the receiver, at stepN.
{
    StructuralElement :: computePrescribedStrainLocalLoadVectorAt(answer, stepN, mode); // ig g.c.s
    FloatMatrix stiff;

    if ( answer.giveSize() && dofsToCondense ) {
        // condense requested dofs
        if ( answer.giveSize() != 0 ) {
            this->computeClampedStiffnessMatrix(stiff, TangentStiffness, stepN);
            this->condense(& stiff, NULL, & answer, dofsToCondense);
        }
    }
}


void
Beam3d :: computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass)
{
    // computes mass matrix of the receiver

    FloatMatrix stiff;
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

    /*
     * StructuralElement::computeMassMatrix(answer, tStep);
     * answer.times(this->giveCrossSection()->give('A'));
     */
    double l = this->giveLength();
    double kappay = this->giveKappayCoeff();
    double kappaz = this->giveKappazCoeff();
    double kappay2 = kappay * kappay;
    double kappaz2 = kappaz * kappaz;
    double density = this->giveMaterial()->give('d', gp);
    double area = this->giveCrossSection()->give(CS_Area);
    double c2y = ( area * density ) / ( ( 1. + 2. * kappay ) * ( 1. + 2. * kappay ) );
    double c2z = ( area * density ) / ( ( 1. + 2. * kappaz ) * ( 1. + 2. * kappaz ) );
    double c1 = ( area * density );

    answer.resize(12, 12);
    answer.zero();

    answer.at(1, 1) = c1 * l / 3.0;
    answer.at(1, 7) = c1 * l / 6.0;
    answer.at(2, 2) = c2z * l * ( 13. / 35. + 7. * kappaz / 5. + 4. * kappaz2 / 3. );
    answer.at(2, 6) = -c2z * l * l * ( 11. / 210. + kappaz * 11. / 60. + kappaz2 / 6. );
    answer.at(2, 8) = c2z * l * ( 9. / 70. + kappaz * 3. / 5. + kappaz2 * 2. / 3. );
    answer.at(2, 12) = c2z * l * l * ( 13. / 420. + kappaz * 3. / 20. + kappaz2 / 6. );
    answer.at(3, 3) = c2y * l * ( 13. / 35. + 7. * kappay / 5. + 4. * kappay2 / 3. );
    answer.at(3, 5) = -c2y * l * l * ( 11. / 210. + kappay * 11. / 60. + kappay2 / 6. );
    answer.at(3, 9) = c2y * l * ( 9. / 70. + kappay * 3. / 5. + kappay2 * 2. / 3. );
    answer.at(3, 11) = c2y * l * l * ( 13. / 420. + kappay * 3. / 20. + kappay2 / 6. );
    answer.at(5, 5) = c2y * l * l * l * ( 1. / 105. + kappay / 30. + kappay2 / 30. );
    answer.at(5, 9) = -c2y * l * l * ( 13. / 420. + kappay * 3. / 20. + kappay2 / 6. );
    answer.at(5, 11) = -c2y * l * l * l * ( 1. / 140. + kappay / 30. + kappay2 / 30. );
    answer.at(6, 6) = c2z * l * l * l * ( 1. / 105. + kappaz / 30. + kappaz2 / 30. );
    answer.at(6, 8) = -c2z * l * l * ( 13. / 420. + kappaz * 3. / 20. + kappaz2 / 6. );
    answer.at(6, 12) = -c2z * l * l * l * ( 1. / 140. + kappaz / 30. + kappaz2 / 30. );


    answer.at(7, 7) = c1 * l / 3.0;
    answer.at(8, 8) = c2z * l * ( 13. / 35. + kappaz * 7. / 5. + kappaz2 * 4. / 3. );
    answer.at(8, 12) = c2z * l * l * ( 11. / 210. + kappaz * 11. / 60. + kappaz2 / 6. );
    answer.at(9, 9) = c2y * l * ( 13. / 35. + kappay * 7. / 5. + kappay2 * 4. / 3. );
    answer.at(9, 11) = c2y * l * l * ( 11. / 210. + kappay * 11. / 60. + kappay2 / 6. );
    answer.at(11, 11) = c2y * l * l * l * ( 1. / 105. + kappay / 30. + kappay2 / 30. );
    answer.at(12, 12) = c2z * l * l * l * ( 1. / 105. + kappaz / 30. + kappaz2 / 30. );

    answer.symmetrized();

    // condense requested dofs
    if ( dofsToCondense ) {
        this->computeClampedStiffnessMatrix(stiff, TangentStiffness, tStep);
        this->condense(& stiff, & answer, NULL, dofsToCondense);
    }

    mass = area * l * density;
}


int
Beam3d :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    double ksi, n1, n2;

    ksi = lcoords.at(1);
    n1  = ( 1. - ksi ) * 0.5;
    n2  = ( 1. + ksi ) * 0.5;

    answer.resize(3);
    answer.at(1) = n1 * this->giveNode(1)->giveCoordinate(1) + n2 *this->giveNode(2)->giveCoordinate(1);
    answer.at(2) = n1 * this->giveNode(1)->giveCoordinate(2) + n2 *this->giveNode(2)->giveCoordinate(2);
    answer.at(3) = n1 * this->giveNode(1)->giveCoordinate(3) + n2 *this->giveNode(2)->giveCoordinate(3);

    return 1;
}


void
Beam3d :: computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    // computes initial stress matrix of receiver (or geometric stiffness matrix)

    FloatMatrix stiff, T;
    FloatArray endForces;

    double l = this->giveLength();
    double kappay = this->giveKappayCoeff();
    double kappaz = this->giveKappazCoeff();
    double kappay2 = kappay * kappay;
    double kappaz2 = kappaz * kappaz;
    double minVal;
    double denomy = ( 1. + 2. * kappay ) * ( 1. + 2. * kappay ), denomz = ( 1. + 2. * kappaz ) * ( 1. + 2. * kappaz );
    double N;

    answer.resize(12, 12);
    answer.zero();

    answer.at(2, 2) = ( 4. * kappaz2 + 4. * kappaz + 6. / 5. ) / denomz;
    answer.at(2, 6) = ( l / 10. ) / denomz;
    answer.at(2, 8) = ( -4. * kappaz2 - 4. * kappaz - 6. / 5. ) / denomz;
    answer.at(2, 12) = ( l / 10. ) / denomz;

    answer.at(3, 3) = ( 4. * kappay2 + 4. * kappay + 6. / 5. ) / denomy;
    answer.at(3, 5) = ( -l / 10. ) / denomy;
    answer.at(3, 9) = ( -4. * kappay2 - 4. * kappay - 6. / 5. ) / denomy;
    answer.at(3, 11) = ( -l / 10. ) / denomy;

    answer.at(5, 5) = l * l * ( kappay2 / 3. + kappay / 3. + 2. / 15. ) / denomy;
    answer.at(5, 9) = ( l / 10. ) / denomy;
    answer.at(5, 11) = -l * l * ( kappay2 / 3. + kappay / 3. + 1. / 30. ) / denomy;

    answer.at(6, 6) = l * l * ( kappaz2 / 3. + kappaz / 3. + 2. / 15. ) / denomz;
    answer.at(6, 8) = ( -l / 10. ) / denomz;
    answer.at(6, 12) = -l * l * ( kappaz2 / 3. + kappaz / 3. + 1. / 30. ) / denomz;

    answer.at(8, 8) = ( 4. * kappaz2 + 4. * kappaz + 6. / 5. ) / denomz;
    answer.at(8, 12) = ( -l / 10. ) / denomz;

    answer.at(9, 9) = ( 4. * kappay2 + 4. * kappay + 6. / 5. ) / denomy;
    answer.at(9, 11) = ( l / 10. ) / denomy;

    answer.at(11, 11) = l * l * ( kappay2 / 3. + kappay / 3. + 2. / 15. ) / denomy;
    answer.at(12, 12) = l * l * ( kappaz2 / 3. + kappaz / 3. + 2. / 15. ) / denomz;

    minVal = min( fabs( answer.at(2, 2) ), fabs( answer.at(3, 3) ) );
    minVal = min( minVal, fabs( answer.at(5, 5) ) );
    minVal = min( minVal, fabs( answer.at(6, 6) ) );

    answer.at(1, 1) = minVal / 1000.;
    answer.at(1, 7) = -answer.at(1, 1);
    answer.at(7, 7) = answer.at(1, 1);

    answer.at(4, 4)   = minVal / 1000.;
    answer.at(4, 10)  = -answer.at(4, 4);
    answer.at(10, 10) = answer.at(4, 4);



    answer.symmetrized();
    // ask end forces in g.c.s
    this->giveEndForcesVector(endForces, tStep);

    N = ( -endForces.at(1) + endForces.at(7) ) / 2.;
    answer.times(N / l);

    // condense requested dofs
    if ( dofsToCondense ) {
        this->computeClampedStiffnessMatrix(stiff, TangentStiffness, tStep);
        this->condense(& stiff, & answer, NULL, dofsToCondense);
    }

    //answer.beLumpedOf (mass);
}


void
Beam3d :: FiberedCrossSectionInterface_computeStrainVectorInFiber(FloatArray &answer, GaussPoint *masterGp,
                                                                  GaussPoint *slaveGp, TimeStep *tStep)
{
    FloatArray masterGpStrain;
    double layerYCoord, layerZCoord;

    this->computeStrainVector(masterGpStrain, masterGp, tStep);
    layerZCoord = slaveGp->giveCoordinate(2);
    layerYCoord = slaveGp->giveCoordinate(1);

    answer.resize(6);  // {Exx,Eyy,Ezz,GMyz,GMzx,GMxy}
    answer.zero();

    answer.at(1) = masterGpStrain.at(1) + masterGpStrain.at(5) * layerZCoord - masterGpStrain.at(6) * layerYCoord;
    answer.at(5) = masterGpStrain.at(2);
    answer.at(6) = masterGpStrain.at(3);
}


Interface *
Beam3d :: giveInterface(InterfaceType interface)
{
    if ( interface == FiberedCrossSectionInterfaceType ) {
        return ( FiberedCrossSectionInterface * ) this;
    }

    return NULL;
}


#ifdef __OOFEG
void Beam3d :: drawRawGeometry(oofegGraphicContext &gc)
{
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    //  if (!go) { // create new one
    WCRec p [ 2 ];   /* poin */
    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
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


void Beam3d :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
{
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();
    double defScale = gc.getDefScale();
    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* poin */
    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);

    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(2, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}
#endif
} // end namespace oofem
