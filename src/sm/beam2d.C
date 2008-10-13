/* $Header: /home/cvs/bp/oofem/sm/src/beam2d.C,v 1.4 2003/04/06 14:08:30 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   file libeam2d.cc

#include "beam2d.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "intarray.h"
#include "flotarry.h"
#include "dof.h"
#include "engngm.h"
#include "boundaryload.h"
#include "structuralcrosssection.h"
#include "structuralmaterial.h"
#include "mathfem.h"
#ifndef __MAKEDEPEND
#include <math.h>
#endif

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#endif

Beam2d :: Beam2d(int n, Domain *aDomain) : StructuralElement(n, aDomain), LayeredCrossSectionInterface()
    // Constructor.
{
    numberOfDofMans     = 2;
    rotationMatrix      = NULL;

    kappa = -1; // set kappa to undef value (should be always > 0.)
    length              = 0.;
    pitch               = 10.;  // a dummy value

    dofsToCondense = NULL;
}

Beam2d :: ~Beam2d()
{
    // destructor

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
Beam2d :: computeBmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer, int li, int ui)
// Returns the strain matrix of the receiver.
{
    double l, ksi, kappa;
    // FloatMatrix* answer ;

    l     = this->giveLength();
    ksi   = 0.5 + 0.5 * aGaussPoint->giveCoordinate(1);
    kappa = this->giveKappaCoeff();

    //answer = new FloatMatrix(3,6) ;
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

    return;
}

void Beam2d :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
  if (!integrationRulesArray) {
    
    // the gauss point is used only when methods from crosssection and/or material
    // classes are requested
    numberOfIntegrationRules = 1;
    integrationRulesArray = new IntegrationRule * [ 1 ];
    integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
    integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Line, 3, _2dBeam);
  }
}


void
Beam2d :: computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at aGaussPoint. Used for numerical calculation of consistent mass
// matrix. Must contain only interpolation for displacement terms,
// not for any rotations. (Inertia forces do not work on rotations).
// r = {u1,w1,fi_y1,u2,w2,fi_y2}^T

{
    double l, ksi, ksi2, ksi3, kappa, c1;
    // FloatMatrix* answer ;

    l     = this->giveLength();
    ksi =   0.5 + 0.5 * aGaussPoint->giveCoordinate(1);
    kappa = this->giveKappaCoeff();
    c1 = 1. + 2. * kappa;
    ksi2 = ksi * ksi;
    ksi3 = ksi2 * ksi;

    //answer = new FloatMatrix(3,6) ;
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

    return;
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

    // return result
    return;
}


void
Beam2d :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
// Returns the stiffness matrix of the receiver, expressed in the global
// axes. No integration over volume done, beam with constant material and crosssection
// parameters assumed.
{
    // compute clamped stifness
    this->computeLocalStiffnessMatrix(answer, rMode, tStep);
    // rotate answer to global coordinate system
    if ( this->updateRotationMatrix() ) {
        answer.rotatedWith(* this->rotationMatrix);
    }

    // return result
    return;
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

    // rotate answer to global coordinate system
    //  this -> giveRotationMatrix () ;
    //  if (rotationMatrix) answer.rotatedWith(*rotationMatrix) ;

    // delete d;
    // return result
    return;
}



int
Beam2d :: computeGtoLRotationMatrix(FloatMatrix &answer) // giveRotationMatrix ()
// Returns the rotation matrix of the receiver.
{
    double sine, cosine;

    answer.resize(6, 6);
    answer.zero();

    sine           = sin( this->givePitch() );
    cosine         = cos(pitch);
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

    return 1;
}


double
Beam2d :: computeVolumeAround(GaussPoint *aGaussPoint)
{
    double weight  = aGaussPoint->giveWeight();
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
    top    = masterGp->giveElement()->giveCrossSection()->give(TOPZCOORD);
    bottom = masterGp->giveElement()->giveCrossSection()->give(BOTTOMZCOORD);
    layerZeta = slaveGp->giveCoordinate(1);
    layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );

    answer.resize(6); // {Exx,Eyy,Ezz,GMyz,GMzx,GMxy}
    answer.zero();

    answer.at(1) = masterGpStrain.at(1) + masterGpStrain.at(2) * layerZCoord;
    answer.at(5) = masterGpStrain.at(3);

    //delete masterGpStrain;
    return;
}


void
Beam2d ::   giveDofManDofIDMask(int inode, EquationID, IntArray &answer) const {
    // returns DofId mask array for inode element node.
    // DofId mask array determines the dof ordering requsted from node.
    // DofId mask array contains the DofID constants (defined in cltypes.h)
    // describing physical meaning of particular DOFs.
    //IntArray* answer = new IntArray (3);
    answer.resize(3);

    answer.at(1) = D_u;
    answer.at(2) = D_w;
    answer.at(3) = R_v;

    return;
}



double Beam2d :: giveLength()
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


double Beam2d :: givePitch()
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
        // delete d;
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

    sine           = sin( this->givePitch() );
    cosine         = cos(pitch);

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
    // delete u;

    /* Substracted in PrintReaction Forces
     * if (loadEndForces = this-> ComputeLoadDependentPartOfLoadVector (tStep)) {
     * loadEndForces -> times(-1.0);
     * answer->add(loadEndForces);
     * }
     */
    this->computePrescribedStrainLoadVectorAt(prescStrainEndForces, tStep, VM_Total);
    if ( prescStrainEndForces.giveSize() ) {
        prescStrainEndForces.times(-1.0);
        answer.add(prescStrainEndForces);
    }

    return;
}

void
Beam2d :: giveEndForcesVector(FloatArray &answer, TimeStep *tStep)
{
    // stress equivalent vector in nodes (vector of internal forces)
    FloatArray u, load;
    FloatMatrix stiffness, T_GtoL, T_NtoG;

    this->computeGtoLRotationMatrix(T_GtoL);
    this->computeGNDofRotationMatrix(T_NtoG, _toGlobalCS);

    // compute stifness matrix in global cs
    this->computeLocalStiffnessMatrix(stiffness, SecantStiffness, tStep);
    stiffness.rotatedWith(T_GtoL);

    // compute vector of unknowns in global cs
    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);
    if ( T_NtoG.isNotEmpty() ) {
        u.rotatedWith(T_NtoG, 'n');
    }

    answer.beProductOf(stiffness, u);

    // substract prescribed strain load
    this->computePrescribedStrainLocalLoadVectorAt(load, tStep, VM_Total);
    if ( load.isNotEmpty() ) {
        load.rotatedWith(T_GtoL, 't');
        answer.substract(load);
    }

    // substract exact end forces due to nonnodal loading
    this->computeLocalForceLoadVector(load, tStep, VM_Total);
    if ( load.isNotEmpty() ) {
        load.rotatedWith(T_GtoL, 't');
        answer.substract(load);
    }

    return;
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
    BoundaryLoad *edgeLoad = dynamic_cast< BoundaryLoad * >( load );
    if ( edgeLoad ) {
        if ( edgeLoad->giveNumberOfDofs() != 3 ) {
            _error("computeEdgeLoadVectorAt: load number of dofs mismatch");
        }

        answer.resize(6);
        answer.zero();
        //  edgeLoad->computeComponentArrayAt(components, tStep, mode);

        // prepare transformation coeffs
        sine           = sin( this->givePitch() );
        cosine         = cos(pitch);

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

        //delete components;
    }

    return;
}


void Beam2d :: printOutputAt(FILE *File, TimeStep *stepN)
{
    // Performs end-of-step operations.

    int i, n;
    FloatArray rg, rl, Fg, Fl;
    FloatMatrix T;

    fprintf(File, "element %d :\n", number);

    //   for (i=0 ; i < numberOfIntegrationRules ; i++)
    //   integrationRulesArray[i]->printOutputAt(file,stepN);

    // ask for global element displacement vector
    this->computeVectorOf(EID_MomentumBalance, VM_Total, stepN, rg);
    if ( this->updateRotationMatrix() ) {
        rl.beProductOf(* this->rotationMatrix, rg);
        // delete rg;
    } else   {
        rl = rg;
    }

    // ask for global element end forces vector
    this->giveEndForcesVector(Fg, stepN);
    // if (computeGNLoadRotationMatrix (T, _toGlobalCS)) Fg.rotatedWith (T, 'n');
    if ( computeGtoLRotationMatrix(T) ) {
        Fg.rotatedWith(T, 'n');
    }

    Fl = Fg;

    fprintf(File, "  local displacements ");
    n = rl.giveSize();
    for ( i = 1; i <= n; i++ ) {
        fprintf( File, " % .4e", rl.at(i) );
    }

    // delete rl;

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

    return;
}


int
Beam2d :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    double ksi, n1, n2;

    ksi = lcoords.at(1);
    n1  = ( 1. - ksi ) * 0.5;
    n2  = ( 1. + ksi ) * 0.5;

    answer.resize(3);
    answer.at(1) = n1 * this->giveNode(1)->giveCoordinate(1) + n2 * this->giveNode(2)->giveCoordinate(1);
    answer.at(3) = n1 * this->giveNode(1)->giveCoordinate(3) + n2 * this->giveNode(2)->giveCoordinate(3);

    return 1;
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

    return;
}



void
Beam2d :: computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass)
{
    // computes mass matrix of the receiver

    FloatMatrix stiff;

    /*
     * StructuralElement::computeMassMatrix(answer, tStep);
     * answer.times(this->giveCrossSection()->give('A'));
     */
    double l = this->giveLength();
    double kappa = this->giveKappaCoeff();
    double kappa2 = kappa * kappa;
    double density = this->giveMaterial()->give('d');
    double area = this->giveCrossSection()->give('A');
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

    return;
}


/*
 * void
 * Beam2d :: computeTemperatureStrainVectorAt (FloatArray& answer, GaussPoint* gp, TimeStep* stepN, ValueModeType mode)
 * {
 * // computes temperature strain vector of the receiver
 * StructuralMaterial * mat = (StructuralMaterial*) this->giveMaterial();
 * StructuralCrossSection* cs = (StructuralCrossSection*) this->giveCrossSection();
 * FloatArray  et, e0 ;
 * double thick;
 *
 * if (this -> giveBodyLoadArray() -> isEmpty()) {answer.resize(0); return;}
 *
 * this -> computeResultingIPTemperatureAt (et, stepN, gp, mode);
 * if (et.giveSize() == 0) {answer.resize(0); return;}
 * if (et.giveSize() < 1) {
 * _error ("computeTemperatureStrainVectorAt - Bad format of TemperatureLoad");
 * exit (1);
 * }
 * mat->giveThermalDilatationVector (e0, gp,stepN);
 *
 * if (e0.giveSize()) {
 * answer.resize (3);
 * answer.zero();
 *
 * answer.at(1) = e0.at(1) * et.at(1);
 * if (et.giveSize() > 1) {
 * thick = cs->give(THICKNESS);
 * answer.at(2) = e0.at(1) * et.at(2)/ thick;   // kappa_x
 * }
 * }
 * //delete et;
 * //delete e0;
 *
 * return ;
 * }
 */

/*
 * int
 * Beam2d :: computeGtoNRotationMatrix (FloatMatrix& answer)
 * // returns transformation matrix from global coordinate set to
 * // nodal coordinate set
 * // return NULL if no trasformation necessary
 * {
 * FloatMatrix *triplet;
 * int i,flag=0,ii;
 *
 * for (i=1; i<= numberOfNodes; i++)
 * flag += this->giveNode(i)->hasLocalCS ();
 * if (flag == 0) {answer.beEmptyMtrx(); return 0 ;}
 *
 * answer.resize(6,6); answer.zero();
 * // loop over nodes
 * for (i=1; i<= numberOfNodes; i++) {
 * ii = (i-1)*3+1 ;
 * if (this->giveNode(i)->hasLocalCS ()) {
 * triplet = this->giveNode(i)->giveLocalCoordinateTriplet();
 * answer.at(ii,ii)     = triplet->at(1,1);
 * answer.at(ii,ii+1)   = triplet->at(1,2);
 * answer.at(ii+1,ii)   = triplet->at(2,1);
 * answer.at(ii+1,ii+1) = triplet->at(2,2);
 * } else {
 * // no transformation - unit matrix as
 * // transformation submatrix for node i
 * answer.at(ii,ii)     = 1.0;
 * answer.at(ii+1,ii+1) = 1.0;
 * }
 * // rotation in plane
 * answer.at(ii+2,ii+2)     = 1.0;
 * }
 *
 * return 1;
 * }
 */

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

    FloatMatrix stiff, T;
    FloatArray endForces;

    double l = this->giveLength();
    double kappa = this->giveKappaCoeff();
    double kappa2 = kappa * kappa;
    double N;

    answer.resize(6, 6);
    answer.zero();

    //answer.at (1,1) = 0.;
    //answer.at (1,4) = 0.;
    answer.at(2, 2) = 4. * kappa2 + 4. * kappa + 6. / 5.;
    answer.at(2, 3) = -l / 10.;
    answer.at(2, 5) = -4. * kappa2 - 4. * kappa - 6. / 5.;
    answer.at(2, 6) = -l / 10.;
    answer.at(3, 3) = l * l * ( kappa2 / 3. + kappa / 3. + 2. / 15. );
    answer.at(3, 5) = l / 10.;
    answer.at(3, 6) = -l * l * ( kappa2 / 3. + kappa / 3. + 1. / 30. );

    //answer.at (4,4) = 0.;
    answer.at(5, 5) = 4. * kappa2 + 4. * kappa + 6. / 5.;
    answer.at(5, 6) = l / 10.;
    answer.at(6, 6) = l * l * ( kappa2 / 3. + kappa / 3. + 2. / 15. );

    answer.at(1, 1) = min( fabs( answer.at(2, 2) ), fabs( answer.at(3, 3) ) ) / 1000.;
    answer.at(1, 4) = -answer.at(1, 1);
    answer.at(4, 4) = answer.at(1, 1);

    answer.symmetrized();
    // ask end forces in g.c.s
    this->giveEndForcesVector(endForces, tStep);
    //if (computeGNLoadRotationMatrix (T, _toGlobalCS)) endForces.rotatedWith (T, 'n');
    if ( computeGtoLRotationMatrix(T) ) {
        endForces.rotatedWith(T, 'n');
    }

    /*
     * this -> giveRotationMatrix () ;
     * if (rotationMatrix) {
     * endForces.rotatedWith (rotationMatrix, 'n');
     * }
     */
    N = ( -endForces.at(1) + endForces.at(4) ) / 2.;
    answer.times( N / ( l * ( 1. + 2. * kappa ) * ( 1. + 2. * kappa ) ) );

    // condense requested dofs
    if ( dofsToCondense ) {
        this->computeClampedStiffnessMatrix(stiff, TangentStiffness, tStep);
        this->condense(& stiff, & answer, NULL, dofsToCondense);
    }

    //answer.beLumpedOf (mass);

    if ( this->updateRotationMatrix() ) {
        answer.rotatedWith(* this->rotationMatrix);
    }

    return;
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
