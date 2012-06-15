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

#include "beam2d.h"
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
// Set up interpolation coordinates
FEI2dLineLin Beam2d :: interp_geom(1, 3);
FEI2dLineHermite Beam2d :: interp_beam(1, 3);

Beam2d :: Beam2d(int n, Domain *aDomain) : StructuralElement(n, aDomain), LayeredCrossSectionInterface()
{
    numberOfDofMans = 2;

    kappa = -1; // set kappa to undef value (should be always > 0.)
    length = 0.;
    pitch = 10.;  // a dummy value

    dofsToCondense = NULL;
}


Beam2d :: ~Beam2d()
{
    delete dofsToCondense;
}


Interface *
Beam2d :: giveInterface(InterfaceType interface)
{
    if ( interface == LayeredCrossSectionInterfaceType ) {
        return ( LayeredCrossSectionInterface * ) this;
    }

    return NULL;
}


void
Beam2d :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the strain matrix of the receiver.
{
    double l, ksi, kappa;

    l     = this->giveLength();
    ksi   = 0.5 + 0.5 * gp->giveCoordinate(1);
    kappa = this->giveKappaCoeff();

    answer.resize(3, 6);
    answer.zero();

    answer.at(1, 1) =  -1. / l;
    answer.at(1, 4) =   1. / l;
    answer.at(2, 2) =   ( 6. - 12. * ksi ) / ( l * l * ( 1. + 2. * kappa ) );
    answer.at(2, 3) =   ( -2. * ( 2. + kappa ) + 6. * ksi ) / ( l * ( 1. + 2. * kappa ) );
    answer.at(2, 5) =   ( -6. + 12. * ksi ) / ( l * l * ( 1. + 2. * kappa ) );
    answer.at(2, 6) =   ( -2. * ( 1. - kappa ) + 6. * ksi ) / ( l * ( 1. + 2. * kappa ) );
    answer.at(3, 2) =   ( -2. * kappa ) / ( l * ( 1. + 2. * kappa ) );
    answer.at(3, 3) =   kappa / ( l * ( 1. + 2. * kappa ) );
    answer.at(3, 5) =   2. * kappa / ( l * ( 1. + 2. * kappa ) );
    answer.at(3, 6) =   kappa / ( l * ( 1. + 2. * kappa ) );
}


void
Beam2d :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( !integrationRulesArray ) {
        // the gauss point is used only when methods from crosssection and/or material
        // classes are requested
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Line, 3, _2dBeam);
    }
}


void
Beam2d :: computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at gp. Used for numerical calculation of consistent mass
// matrix. Must contain only interpolation for displacement terms,
// not for any rotations. (Inertia forces do not work on rotations).
// r = {u1,w1,fi_y1,u2,w2,fi_y2}^T

{
    double l, ksi, ksi2, ksi3, kappa, c1;

    l     = this->giveLength();
    ksi =   0.5 + 0.5 * gp->giveCoordinate(1);
    kappa = this->giveKappaCoeff();
    c1 = 1. + 2. * kappa;
    ksi2 = ksi * ksi;
    ksi3 = ksi2 * ksi;

    answer.resize(3, 6);
    answer.zero();

    answer.at(1, 1) = 1. - ksi;
    answer.at(1, 4) = ksi;
    answer.at(2, 2) = ( ( 1. + 2. * kappa ) - 2. * kappa * ksi - 3. * ksi2 + 2. * ksi3 ) / c1;
    answer.at(2, 3) = l * ( -( 1. + kappa ) * ksi + ( 2. + kappa ) * ksi2 - ksi3 ) / c1;
    answer.at(2, 5) = ( 2. * kappa * ksi + 3. * ksi2 - 2. * ksi3 ) / c1;
    answer.at(2, 6) = l * ( kappa * ksi + ( 1. - kappa ) * ksi2 - ksi3 ) / c1;
    answer.at(3, 2) = ( 6. * ksi - 6. * ksi2 ) / ( l * c1 );
    answer.at(3, 3) = ( ( 1. + 2. * kappa ) - 2. * ( 2. + kappa ) * ksi + 3. * ksi2 ) / c1;
    answer.at(3, 5) = ( -6. * ksi + 6. * ksi2 ) / ( l * c1 );
    answer.at(3, 6) = ( -2. * ( 1. - kappa ) * ksi + 3. * ksi2 ) / c1;
}


void
Beam2d :: computeLocalStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
// Returns the stiffness matrix of the receiver, expressed in the local
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
Beam2d :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
// Returns the stiffness matrix of the receiver, expressed in the global
// axes. No integration over volume done, beam with constant material and crosssection
// parameters assumed.
{
    // compute clamped stifness
    this->computeLocalStiffnessMatrix(answer, rMode, tStep);
}


void
Beam2d :: computeClampedStiffnessMatrix(FloatMatrix &answer,
                                        MatResponseMode rMode, TimeStep *tStep)
// Returns the stiffness matrix of the receiver, expressed in the global
// axes. No integration over volume done, beam with constant material and crosssection
// parameters assumed.
{
    double l, l2, l3, ei, kappa, c1;
    FloatMatrix d;

    this->computeConstitutiveMatrixAt(d, rMode, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);
    l     = this->giveLength();
    l2    = l * l;
    l3    = l2 * l;
    kappa = this->giveKappaCoeff();
    c1 = 1. + 2. * kappa;
    ei = d.at(2, 2);

    answer.resize(6, 6);
    answer.zero();

    answer.at(1, 1) =  d.at(1, 1) / l;
    answer.at(1, 4) = -d.at(1, 1) / l;
    answer.at(2, 2) =  ei * 12. / ( l3 * c1 );
    answer.at(2, 3) = -ei * 6.  / ( l2 * c1 );
    answer.at(2, 5) = -ei * 12. / ( l3 * c1 );
    answer.at(2, 6) = -ei * 6.  / ( l2 * c1 );
    answer.at(3, 3) =  ei * 2. * ( 2. + kappa ) / ( l * c1 );
    answer.at(3, 5) =  ei * 6.  / ( l2 * c1 );
    answer.at(3, 6) =  ei * 2. * ( 1. - kappa ) / ( l * c1 );
    answer.at(4, 4) =  d.at(1, 1) / l;
    answer.at(5, 5) =  ei * 12. / ( l3 * c1 );
    answer.at(5, 6) =  ei * 6.  / ( l2 * c1 );
    answer.at(6, 6) =  ei * 2. * ( 2. + kappa ) / ( l * c1 );

    answer.symmetrized();  // symmetrize answer
}


bool
Beam2d :: computeGtoLRotationMatrix(FloatMatrix &answer)
// Returns the rotation matrix of the receiver.
{
    double sine, cosine;

    answer.resize(6, 6);
    answer.zero();

    sine = sin( this->givePitch() );
    cosine  = cos(pitch);
    answer.at(1, 1) =  cosine;
    answer.at(1, 2) =  sine;
    answer.at(2, 1) = -sine;
    answer.at(2, 2) =  cosine;
    answer.at(3, 3) =  1.;
    answer.at(4, 4) =  cosine;
    answer.at(4, 5) =  sine;
    answer.at(5, 4) = -sine;
    answer.at(5, 5) =  cosine;
    answer.at(6, 6) =  1.;

    return true;
}


double
Beam2d :: computeVolumeAround(GaussPoint *gp)
{
    double weight  = gp->giveWeight();
    return weight * 0.5 * this->giveLength();
}


void
Beam2d :: computeStrainVectorInLayer(FloatArray &answer, GaussPoint *masterGp,
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

    answer.at(1) = masterGpStrain.at(1) + masterGpStrain.at(2) * layerZCoord;
    answer.at(5) = masterGpStrain.at(3);
}


void
Beam2d :: giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const
{
    answer.setValues(3, D_u, D_w, R_v);
}


double
Beam2d :: giveLength()
// Returns the length of the receiver.
{
    double dx, dy;
    Node *nodeA, *nodeB;

    if ( length == 0. ) {
        nodeA   = this->giveNode(1);
        nodeB   = this->giveNode(2);
        dx      = nodeB->giveCoordinate(1) - nodeA->giveCoordinate(1);
        dy      = nodeB->giveCoordinate(3) - nodeA->giveCoordinate(3);
        length  = sqrt(dx * dx + dy * dy);
    }

    return length;
}


double
Beam2d :: givePitch()
// Returns the pitch of the receiver.
{
    double xA, xB, yA, yB;
    Node *nodeA, *nodeB;

    if ( pitch == 10. ) {             // 10. : dummy initialization value
        nodeA  = this->giveNode(1);
        nodeB  = this->giveNode(2);
        xA     = nodeA->giveCoordinate(1);
        xB     = nodeB->giveCoordinate(1);
        yA     = nodeA->giveCoordinate(3);
        yB     = nodeB->giveCoordinate(3);
        pitch  = atan2(yB - yA, xB - xA);
    }

    return pitch;
}


double
Beam2d :: giveKappaCoeff()
{
    // returns kappa coeff
    // kappa = (6*E*I)/(k*G*A*l^2)

    if ( kappa < 0. ) {
        FloatMatrix d;
        double l = this->giveLength();

        this->computeConstitutiveMatrixAt( d, TangentStiffness, integrationRulesArray [ 0 ]->getIntegrationPoint(0), domain->giveEngngModel()->giveCurrentStep() );
        kappa = 6. * d.at(2, 2) / ( d.at(3, 3) * l * l );
    }

    return kappa;
}


int
Beam2d :: giveLocalCoordinateSystem(FloatMatrix &answer)
//
// returns a unit vectors of local coordinate system at element
// stored rowwise (mainly used by some materials with ortho and anisotrophy)
//
{
    double sine, cosine;

    answer.resize(3, 3);
    answer.zero();

    sine = sin( this->givePitch() );
    cosine = cos(pitch);

    answer.at(1, 1) = cosine;
    answer.at(1, 2) = sine;
    answer.at(2, 1) = -sine;
    answer.at(2, 2) = cosine;
    answer.at(3, 3) = 1.0;

    return 1;
}


IRResultType
Beam2d :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    // first call parent
    StructuralElement :: initializeFrom(ir);

    if ( ir->hasField(IFT_Beam2d_dofstocondense, "dofstocondense") ) {
        IntArray val;
        IR_GIVE_FIELD(ir, val, IFT_Beam2d_dofstocondense, "dofstocondense"); // Macro
        if ( val.giveSize() >= 6 ) {
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
Beam2d :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
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
Beam2d :: giveEndForcesVector(FloatArray &answer, TimeStep *tStep)
{
    // stress equivalent vector in nodes (vector of internal forces)
    FloatArray u, load;
    FloatMatrix stiffness;

    // compute stifness matrix in global cs
    this->computeLocalStiffnessMatrix(stiffness, SecantStiffness, tStep);

    // compute vector of unknowns in global cs
    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);
    answer.beProductOf(stiffness, u);

    // subtract prescribed strain load
    this->computePrescribedStrainLocalLoadVectorAt(load, tStep, VM_Total);
    if ( load.isNotEmpty() ) {
        answer.subtract(load);
    }

    // subtract exact end forces due to nonnodal loading
    this->computeLocalForceLoadVector(load, tStep, VM_Total);
    if ( load.isNotEmpty() ) {
        answer.subtract(load);
    }
}


void
Beam2d :: computeEdgeLoadVectorAt(FloatArray &answer, Load *load, int iedge, TimeStep *tStep, ValueModeType mode)
{
    int i;
    FloatArray coords, components, help;
    FloatMatrix T;
    double l = this->giveLength();
    double kappa = this->giveKappaCoeff();
    double fx, fz, fm, dfx, dfz, dfm;
    double cosine, sine;

    // evaluates the receivers edge load vector
    // for clamped beam
    //
    BoundaryLoad *edgeLoad = dynamic_cast< BoundaryLoad * >(load);
    if ( edgeLoad ) {
        if ( edgeLoad->giveNumberOfDofs() != 3 ) {
            _error("computeEdgeLoadVectorAt: load number of dofs mismatch");
        }

        answer.resize(6);
        answer.zero();
        //  edgeLoad->computeComponentArrayAt(components, tStep, mode);

        // prepare transformation coeffs
        sine = sin( this->givePitch() );
        cosine = cos(pitch);

        switch ( edgeLoad->giveClassID() ) {
        case ConstantEdgeLoadClass:
            coords.resize(1);
            if ( edgeLoad->giveFormulationType() == BoundaryLoad :: BL_EntityFormulation ) {
                coords.at(1) = 0.0;
            } else {
                coords = * ( this->giveNode(1)->giveCoordinates() );
            }

            edgeLoad->computeValueAt(components, tStep, coords, mode);

            if ( edgeLoad->giveCoordSystMode() == BoundaryLoad :: BL_GlobalMode ) {
                fx = cosine * components.at(1) + sine *components.at(2);
                fz = -sine *components.at(1) + cosine *components.at(2);
                fm = components.at(3);
            } else {
                fx = components.at(1);
                fz = components.at(2);
                fm = components.at(3);
            }

            answer.at(1) = fx * l / 2.;
            answer.at(2) = fz * l / 2. + fm / ( 1. + 2. * kappa );
            answer.at(3) = ( -1. ) * fz * l * l / 12. + fm * l * kappa / ( 1. + 2. * kappa );
            answer.at(4) = fx * l / 2.;
            answer.at(5) = fz * l / 2. - fm / ( 1. + 2. * kappa );
            answer.at(6) = fz * l * l / 12. + fm * l * kappa / ( 1. + 2. * kappa );
            break;
        case LinearEdgeLoadClass:

            components.resize(6);

            if ( edgeLoad->giveFormulationType() == BoundaryLoad :: BL_EntityFormulation ) {
                coords.resize(1);
                coords.at(1) = -1.0;
                edgeLoad->computeValueAt(help, tStep, coords, mode);
                for ( i = 1; i <= 3; i++ ) {
                    components.at(i) = help.at(i);
                }

                coords.at(1) = 1.0;
                edgeLoad->computeValueAt(help, tStep, coords, mode);
                for ( i = 1; i <= 3; i++ ) {
                    components.at(i + 3) = help.at(i);
                }
            } else {
                coords = * ( this->giveNode(1)->giveCoordinates() );
                edgeLoad->computeValueAt(help, tStep, coords, mode);
                for ( i = 1; i <= 3; i++ ) {
                    components.at(i) = help.at(i);
                }

                coords = * ( this->giveNode(2)->giveCoordinates() );
                edgeLoad->computeValueAt(help, tStep, coords, mode);
                for ( i = 1; i <= 3; i++ ) {
                    components.at(i + 3) = help.at(i);
                }
            }


            if ( edgeLoad->giveCoordSystMode() == BoundaryLoad :: BL_GlobalMode ) {
                fx = cosine * components.at(1) + sine *components.at(2);
                fz = -sine *components.at(1) + cosine *components.at(2);
                fm = components.at(3);

                dfx = cosine * ( components.at(4) - components.at(1) ) + sine * ( components.at(5) - components.at(2) );
                dfz = -sine * ( components.at(4) - components.at(1) ) + cosine * ( components.at(5) - components.at(2) );
                dfm = components.at(6) - components.at(3);
            } else {
                fx = components.at(1);
                fz = components.at(2);
                fm = components.at(3);
                dfx = components.at(4) - components.at(1);
                dfz = components.at(5) - components.at(2);
                dfm = components.at(6) - components.at(3);
            }

            answer.at(1) = fx * l / 2. + dfx * l / 6.;
            answer.at(2) = fz * l / 2. + dfz * l * ( 20. * kappa + 9 ) / ( 60. * ( 1. + 2. * kappa ) ) +
                           fm / ( 1. + 2. * kappa ) + dfm * ( 1. / 2. ) / ( 1. + 2. * kappa );
            answer.at(3) = ( -1. ) * fz * l * l / 12. - dfz * l * l * ( 5. * kappa + 2. ) / ( 60. * ( 1. + 2. * kappa ) ) +
                           fm * l * kappa / ( 1. + 2. * kappa ) + dfm * l * ( 4. * kappa - 1. ) / ( 12. * ( 1. + 2. * kappa ) );
            answer.at(4) = fx * l / 2. + dfx * l / 3.;
            answer.at(5) = fz * l / 2. + dfz * l * ( 40. * kappa + 21 ) / ( 60. * ( 1. + 2. * kappa ) ) -
                           fm / ( 1. + 2. * kappa ) - dfm * ( 1. / 2. ) / ( ( 1. + 2. * kappa ) );
            answer.at(6) = fz * l * l / 12. + dfz * l * l * ( 5. * kappa + 3. ) / ( 60. * ( 1. + 2. * kappa ) ) +
                           fm * l * kappa / ( 1. + 2. * kappa ) + dfm * l * ( 8. * kappa + 1. ) / ( 12. * ( 1. + 2. * kappa ) );
            break;
        default:
            _error("computeEdgeLoadVectorAt: unsupported load type");
        }
    }
}


void
Beam2d :: computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode)
{
    StructuralElement::computeBodyLoadVectorAt(answer, load, tStep, mode);
    answer.times(this->giveCrossSection()->give(CS_Area));
}


void
Beam2d :: printOutputAt(FILE *File, TimeStep *stepN)
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
Beam2d :: computeLocalForceLoadVector(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
// Computes the load vector of the receiver, at stepN.
{
    FloatMatrix stiff;

    StructuralElement :: computeLocalForceLoadVector(answer, stepN, mode);

    // condense requested dofs in local c.s
    if ( answer.giveSize() && dofsToCondense ) {
        if ( answer.giveSize() != 0 ) {
            this->computeClampedStiffnessMatrix(stiff, TangentStiffness, stepN);
            this->condense(& stiff, NULL, & answer, dofsToCondense);
        }
    }
}


void
Beam2d :: computePrescribedStrainLocalLoadVectorAt(FloatArray &answer, TimeStep *stepN, ValueModeType mode)
// Computes the load vector of the receiver, at stepN.
{
    StructuralElement :: computePrescribedStrainLocalLoadVectorAt(answer, stepN, mode);
    FloatMatrix stiff;

    // condense requested dofs
    if ( answer.giveSize() && dofsToCondense ) {
        if ( answer.giveSize() != 0 ) {
            this->computeClampedStiffnessMatrix(stiff, TangentStiffness, stepN);
            this->condense(& stiff, NULL, & answer, dofsToCondense);
        }
    }
}


void
Beam2d :: computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass)
{
    // computes mass matrix of the receiver

    FloatMatrix stiff;
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

    /*
     * StructuralElement::computeMassMatrix(answer, tStep);
     * answer.times(this->giveCrossSection()->give('A'));
     */
    double l = this->giveLength();
    double kappa = this->giveKappaCoeff();
    double kappa2 = kappa * kappa;
    double density = this->giveMaterial()->give('d', gp);
    double area = this->giveCrossSection()->give(CS_Area);
    double c2 = ( area * density ) / ( ( 1. + 2. * kappa ) * ( 1. + 2. * kappa ) );
    double c1 = ( area * density );

    answer.resize(6, 6);
    answer.zero();

    answer.at(1, 1) = c1 * l / 3.0;
    answer.at(1, 4) = c1 * l / 6.0;
    answer.at(2, 2) = c2 * l * ( 13. / 35. + 7. * kappa / 5. + 4. * kappa2 / 3. );
    answer.at(2, 3) = -c2 * l * l * ( 11. / 210. + kappa * 11. / 60. + kappa2 / 6. );
    answer.at(2, 5) = c2 * l * ( 9. / 70. + kappa * 3. / 5. + kappa2 * 2. / 3. );
    answer.at(2, 6) = c2 * l * l * ( 13. / 420. + kappa * 3. / 20. + kappa2 / 6. );
    answer.at(3, 3) = c2 * l * l * l * ( 1. / 105. + kappa / 30. + kappa2 / 30. );
    answer.at(3, 5) = -c2 * l * l * ( 13. / 420. + kappa * 3. / 20. + kappa2 / 6. );
    answer.at(3, 6) = -c2 * l * l * l * ( 1. / 140. + kappa / 30. + kappa2 / 30. );

    answer.at(4, 4) = c1 * l / 3.0;
    answer.at(5, 5) = c2 * l * ( 13. / 35. + kappa * 7. / 5. + kappa2 * 4. / 3. );
    answer.at(5, 6) = c2 * l * l * ( 11. / 210. + kappa * 11. / 60. + kappa2 / 6. );
    answer.at(6, 6) = c2 * l * l * l * ( 1. / 105. + kappa / 30. + kappa2 / 30. );

    answer.symmetrized();

    // condense requested dofs
    if ( dofsToCondense ) {
        this->computeClampedStiffnessMatrix(stiff, TangentStiffness, tStep);
        this->condense(& stiff, & answer, NULL, dofsToCondense);
    }

    mass = l * area * density;
}

/*
 * void
 * Beam2d :: giveMassMtrxIntegrationgMask (IntArray& answer)
 * {
 * answer.resize (3);
 *
 * answer.at(1) = 1;
 * answer.at(2) = 1;
 * answer.at(3) = 0;
 * }
 */

void
Beam2d :: computeInitialStressMatrix(FloatMatrix &answer, TimeStep *tStep)
{
    // computes initial stress matrix of receiver (or geometric stiffness matrix)

    FloatMatrix stiff;
    FloatArray endForces;

    double l = this->giveLength();
    double kappa = this->giveKappaCoeff();
    double kappa2 = kappa * kappa;
    double N;

    answer.resize(6, 6);
    answer.zero();

    answer.at(2, 2) = 4. * kappa2 + 4. * kappa + 6. / 5.;
    answer.at(2, 3) = -l / 10.;
    answer.at(2, 5) = -4. * kappa2 - 4. * kappa - 6. / 5.;
    answer.at(2, 6) = -l / 10.;
    answer.at(3, 3) = l * l * ( kappa2 / 3. + kappa / 3. + 2. / 15. );
    answer.at(3, 5) = l / 10.;
    answer.at(3, 6) = -l * l * ( kappa2 / 3. + kappa / 3. + 1. / 30. );

    answer.at(5, 5) = 4. * kappa2 + 4. * kappa + 6. / 5.;
    answer.at(5, 6) = l / 10.;
    answer.at(6, 6) = l * l * ( kappa2 / 3. + kappa / 3. + 2. / 15. );

    answer.at(1, 1) = min( fabs( answer.at(2, 2) ), fabs( answer.at(3, 3) ) ) / 1000.;
    answer.at(1, 4) = -answer.at(1, 1);
    answer.at(4, 4) = answer.at(1, 1);

    answer.symmetrized();
    // ask end forces in g.c.s
    this->giveEndForcesVector(endForces, tStep);

    N = ( -endForces.at(1) + endForces.at(4) ) / 2.;
    answer.times( N / ( l * ( 1. + 2. * kappa ) * ( 1. + 2. * kappa ) ) );

    // condense requested dofs
    if ( dofsToCondense ) {
        this->computeClampedStiffnessMatrix(stiff, TangentStiffness, tStep);
        this->condense(& stiff, & answer, NULL, dofsToCondense);
    }

    //answer.beLumpedOf (mass);
}


#ifdef __OOFEG
void Beam2d :: drawRawGeometry(oofegGraphicContext &gc)
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
    p [ 0 ].y = 0.;
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(3);
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = 0.;
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveCoordinate(3);
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}


void Beam2d :: drawDeformedGeometry(oofegGraphicContext &gc, UnknownType type)
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
    p [ 0 ].y = 0.;
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);

    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, EID_MomentumBalance, defScale);
    p [ 1 ].y = 0.;
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, EID_MomentumBalance, defScale);
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}
#endif
} // end namespace oofem
