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

#include "../sm/Elements/Beams/beam2d.h"
#include "../sm/Materials/structuralms.h"
#include "fei2dlinelin.h"
#include "fei2dlinehermite.h"
#include "node.h"
#include "material.h"
#include "crosssection.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "intarray.h"
#include "floatarray.h"
#include "engngm.h"
#include "boundaryload.h"
#include "mathfem.h"
#include "classfactory.h"
#include "elementinternaldofman.h"
#include "masterdof.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(Beam2d);

// Set up interpolation coordinates
FEI2dLineLin Beam2d :: interp_geom(1, 3);
FEI2dLineHermite Beam2d :: interp_beam(1, 3);

Beam2d :: Beam2d(int n, Domain *aDomain) : StructuralElement(n, aDomain), LayeredCrossSectionInterface()
{
    numberOfDofMans = 2;

    kappa = -1; // set kappa to undef value (should be always > 0.)
    length = 0.;
    pitch = 10.;  // a dummy value
    numberOfGaussPoints = 4;

    ghostNodes [ 0 ] = ghostNodes [ 1 ] = NULL;
    numberOfCondensedDofs = 0;
}


Beam2d :: ~Beam2d()
{
    if ( ghostNodes [ 0 ] ) {
        delete ghostNodes [ 0 ];
    }
    if ( ghostNodes [ 1 ] ) {
        delete ghostNodes [ 1 ];
    }
}


FEInterpolation *Beam2d :: giveInterpolation() const { return & interp_geom; }


Interface *
Beam2d :: giveInterface(InterfaceType interface)
{
    if ( interface == LayeredCrossSectionInterfaceType ) {
        return static_cast< LayeredCrossSectionInterface * >(this);
    }

    return NULL;
}


void
Beam2d :: computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui)
// Returns the strain matrix of the receiver.
{
    double l, ksi, kappa, c1;
    TimeStep *tStep = domain->giveEngngModel()->giveCurrentStep();

    l     = this->computeLength();
    ksi   = 0.5 + 0.5 * gp->giveNaturalCoordinate(1);
    kappa = this->giveKappaCoeff(tStep);
    c1 = 1. + 2. * kappa;

    answer.resize(3, 6);
    answer.zero();

    answer.at(1, 1) =  -1. / l;
    answer.at(1, 4) =   1. / l;
    answer.at(2, 2) =   ( 6. - 12. * ksi ) / ( l * l * c1 );
    answer.at(2, 3) =   ( -2. * ( 2. + kappa ) + 6. * ksi ) / ( l * c1 );
    answer.at(2, 5) =   ( -6. + 12. * ksi ) / ( l * l * c1 );
    answer.at(2, 6) =   ( -2. * ( 1. - kappa ) + 6. * ksi ) / ( l * c1 );
    answer.at(3, 2) =   ( -2. * kappa ) / ( l * c1 );
    answer.at(3, 3) =   kappa / ( c1 );
    answer.at(3, 5) =   2. * kappa / ( l * c1 );
    answer.at(3, 6) =   kappa / ( c1 );
}


void
Beam2d :: computeGaussPoints()
// Sets up the array of Gauss Points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        // the gauss point is used only when methods from crosssection and/or material
        // classes are requested
        integrationRulesArray.resize(1);
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], this->numberOfGaussPoints, this);
    }
}


void
Beam2d :: computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer)
// Returns the displacement interpolation matrix {N} of the receiver, eva-
// luated at gp. Used for numerical calculation of consistent mass
// matrix. Must contain only interpolation for displacement terms,
// not for any rotations. (Inertia forces do not work on rotations).
// r = {u1,w1,fi_y1,u2,w2,fi_y2}^T

{
    double l, ksi, ksi2, ksi3, kappa, c1;
    TimeStep *tStep = this->domain->giveEngngModel()->giveCurrentStep();

    l     = this->computeLength();
    ksi =   0.5 + 0.5 * iLocCoord.at(1);
    kappa = this->giveKappaCoeff(tStep);
    c1 = 1. + 2. * kappa;
    ksi2 = ksi * ksi;
    ksi3 = ksi2 * ksi;

    answer.resize(3, 6);
    answer.zero();

    answer.at(1, 1) = 1. - ksi;
    answer.at(1, 4) = ksi;
    answer.at(2, 2) = ( c1 - 2. * kappa * ksi - 3. * ksi2 + 2. * ksi3 ) / c1;
    answer.at(2, 3) = l * ( -( 1. + kappa ) * ksi + ( 2. + kappa ) * ksi2 - ksi3 ) / c1;
    answer.at(2, 5) = ( 2. * kappa * ksi + 3. * ksi2 - 2. * ksi3 ) / c1;
    answer.at(2, 6) = l * ( kappa * ksi + ( 1. - kappa ) * ksi2 - ksi3 ) / c1;
    answer.at(3, 2) = ( 6. * ksi - 6. * ksi2 ) / ( l * c1 );
    answer.at(3, 3) = ( c1 - 2. * ( 2. + kappa ) * ksi + 3. * ksi2 ) / c1;
    answer.at(3, 5) = ( -6. * ksi + 6. * ksi2 ) / ( l * c1 );
    answer.at(3, 6) = ( -2. * ( 1. - kappa ) * ksi + 3. * ksi2 ) / c1;
}

void
Beam2d :: computeLocalStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
// Returns the stiffness matrix of the receiver, expressed in the local
// axes. No integration over volume done, beam with constant material and crosssection
// parameters assumed.
{
    this->computeClampedStiffnessMatrix(answer, rMode, tStep);
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
}


void
Beam2d :: computeClampedStiffnessMatrix(FloatMatrix &answer,
                                        MatResponseMode rMode, TimeStep *tStep)
// Returns the stiffness matrix of the receiver, expressed in the global
// axes. No integration over volume done, beam with constant material and crosssection
// parameters assumed.
{
    double l = this->computeLength();
    FloatMatrix B, d, DB;
    answer.clear();
    for ( GaussPoint *gp : *this->giveDefaultIntegrationRulePtr() ) {
        this->computeBmatrixAt(gp, B);
        this->computeConstitutiveMatrixAt(d, rMode, gp, tStep);
        double dV = gp->giveWeight() * 0.5 * l;
        DB.beProductOf(d, B);
        answer.plusProductSymmUpper(B, DB, dV);
    }
    answer.symmetrized();
}


void
Beam2d :: computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->give2dBeamStiffMtrx(answer, rMode, gp, tStep);
}


void
Beam2d :: computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    this->giveStructuralCrossSection()->giveGeneralizedStress_Beam2d(answer, gp, strain, tStep);
}


bool
Beam2d :: computeGtoLRotationMatrix(FloatMatrix &answer)
// Returns the rotation matrix of the receiver.
{
    double sine, cosine;

    int ndofs = computeNumberOfGlobalDofs();
    answer.resize(ndofs, ndofs);
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

    for ( int i = 7; i <= ndofs; i++ ) {
        answer.at(i, i) = 1.0;
    }

    if ( this->hasDofs2Condense() ) {
        int condensedDofCounter = 0;
        DofIDItem dofids[] = {
            D_u, D_w, R_v
        };
        FloatMatrix l2p(6, ndofs); // local -> primary
        l2p.zero();
        // loop over nodes
        for ( int inode = 0; inode < 2; inode++ ) {
            // loop over DOFs
            for ( int idof = 0; idof < 3; idof++ ) {
                int eq = inode * 3 + idof + 1;
                if ( ghostNodes [ inode ] ) {
                    if ( ghostNodes [ inode ]->hasDofID(dofids [ idof ]) ) {
                        condensedDofCounter++;
                        l2p.at(eq, 6 + condensedDofCounter) = 1.0;
                        continue;
                    }
                }
                l2p.at(eq, eq) = 1.0;
            }
        }

        FloatMatrix g2l(answer);
        answer.beProductOf(l2p, g2l);
    }

    return true;
}


double
Beam2d :: computeVolumeAround(GaussPoint *gp)
{
    return 0.5 * this->computeLength() * gp->giveWeight();
}


void
Beam2d :: computeStrainVectorInLayer(FloatArray &answer, const FloatArray &masterGpStrain, GaussPoint *masterGp, GaussPoint *slaveGp, TimeStep *tStep)
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

    answer.resize(2); // {Exx,GMzx}

    answer.at(1) = masterGpStrain.at(1) + masterGpStrain.at(2) * layerZCoord;
    answer.at(2) = masterGpStrain.at(3);
}


void
Beam2d :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {
        D_u, D_w, R_v
    };
}


double
Beam2d :: computeLength()
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
Beam2d :: giveKappaCoeff(TimeStep *tStep)
{
    // returns kappa coeff
    // kappa = (6*E*I)/(k*G*A*l^2)

    if ( kappa < 0. ) {
        FloatMatrix d;
        double l = this->computeLength();

        this->computeConstitutiveMatrixAt(d, ElasticStiffness, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);
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
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    // first call parent
    StructuralElement :: initializeFrom(ir);

    if ( ir->hasField(_IFT_Beam2d_dofstocondense) ) {
        IntArray val;
        IR_GIVE_FIELD(ir, val, _IFT_Beam2d_dofstocondense);
        if ( val.giveSize() >= 6 ) {
            OOFEM_WARNING("wrong input data for condensed dofs");
            return IRRT_BAD_FORMAT;
        }

        DofIDItem mask[] = {
            D_u, D_w, R_v
        };
        this->numberOfCondensedDofs = val.giveSize();
        for ( int i = 1; i <= val.giveSize(); i++ ) {
            if ( val.at(i) <= 3 ) {
                if ( ghostNodes [ 0 ] == NULL ) {
                    ghostNodes [ 0 ] = new ElementDofManager(1, giveDomain(), this);
                }
                ghostNodes [ 0 ]->appendDof( new MasterDof(ghostNodes [ 0 ], mask [ val.at(i) - 1 ]) );
            } else {
                if ( ghostNodes [ 1 ] == NULL ) {
                    ghostNodes [ 1 ] = new ElementDofManager(1, giveDomain(), this);
                }
                ghostNodes [ 1 ]->appendDof( new MasterDof(ghostNodes [ 1 ], mask [ val.at(i) - 4 ]) );
            }
        }

    }
    return IRRT_OK;
}



void
Beam2d :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    StructuralElement :: giveInternalForcesVector(answer, tStep, useUpdatedGpRecord);

    ///@todo Pretty sure this won't work for nonlinear problems. I think dofsToCondense should just be replaced by an extra slave node.
    /*
     * if ( this->dofsToCondense ) {
     *  FloatMatrix stiff;
     *  this->computeClampedStiffnessMatrix(stiff, TangentStiffness, tStep);
     *  this->condense(& stiff, NULL, & answer, this->dofsToCondense);
     * }
     */
}


void
Beam2d :: giveEndForcesVector(FloatArray &answer, TimeStep *tStep)
{
    // stress equivalent vector in nodes (vector of internal forces)
    FloatArray load;

    this->giveInternalForcesVector(answer, tStep, false);

    // subtract exact end forces due to nonnodal loading
    this->computeLocalForceLoadVector(load, tStep, VM_Total);
    if ( load.isNotEmpty() ) {
        answer.subtract(load);
    }
}


void
Beam2d :: computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int edge, CharType type, ValueModeType mode, TimeStep *tStep)
{
    answer.clear();

    if ( edge != 1 ) {
        OOFEM_ERROR("Beam2D only has 1 edge (the midline) that supports loads. Attempted to apply load to edge %d", edge);
    }

    if ( type != ExternalForcesVector ) {
        return;
    }

    double l = this->computeLength();
    FloatArray coords, t;
    FloatMatrix N, T;

    answer.clear();
    for ( GaussPoint *gp : *this->giveDefaultIntegrationRulePtr() ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        this->computeNmatrixAt(lcoords, N);
        if ( load ) {
            this->computeGlobalCoordinates(coords, lcoords);
            load->computeValues(t, tStep, coords, { D_u, D_w, R_v }, mode);
        } else {
            load->computeValues(t, tStep, lcoords, { D_u, D_w, R_v }, mode);
        }

        if ( load->giveCoordSystMode() == Load :: CST_Global ) {
            if ( this->computeLoadGToLRotationMtrx(T) ) {
                t.rotatedWith(T, 'n');
            }
        }

        double dl = gp->giveWeight() * 0.5 * l;
        answer.plusProduct(N, t, dl);
    }

    // Loads from sets expects global c.s.
    this->computeGtoLRotationMatrix(T);
    answer.rotatedWith(T, 't');
    ///@todo Decide if we want local or global c.s. for loads over sets.
}


void
Beam2d :: computeEdgeLoadVectorAt(FloatArray &answer, Load *load, int iedge, TimeStep *tStep, ValueModeType mode)
{
    FloatArray coords, components, components2;
    double l = this->computeLength();
    double kappa = this->giveKappaCoeff(tStep);
    double fx, fz, fm, dfx, dfz, dfm;
    double cosine, sine;

    // evaluates the receivers edge load vector
    // for clamped beam
    //
    BoundaryLoad *edgeLoad = dynamic_cast< BoundaryLoad * >(load);
    if ( edgeLoad ) {
        answer.resize(6);
        answer.zero();

        // prepare transformation coeffs
        sine = sin( this->givePitch() );
        cosine = cos(pitch);


        switch ( edgeLoad->giveApproxOrder() ) {
        case 0:
            coords.resize(1);
            if ( edgeLoad->giveFormulationType() == Load :: FT_Entity ) {
                coords.at(1) = 0.0;
            } else {
                coords = * ( this->giveNode(1)->giveCoordinates() );
            }

            edgeLoad->computeValues(components, tStep, coords, { D_u, D_w, R_v }, mode);

            if ( edgeLoad->giveCoordSystMode() == Load :: CST_Global ) {
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

        case 1:
            components.resize(6);

            if ( edgeLoad->giveFormulationType() == Load :: FT_Entity ) {
                edgeLoad->computeValues(components, tStep, { -1.0 }, { D_u, D_w, R_v }, mode);
                edgeLoad->computeValues(components2, tStep, { 1.0 }, { D_u, D_w, R_v }, mode);
            } else {
                edgeLoad->computeValues(components, tStep, * this->giveNode(1)->giveCoordinates(), { D_u, D_w, R_v }, mode);
                edgeLoad->computeValues(components2, tStep, * this->giveNode(2)->giveCoordinates(), { D_u, D_w, R_v }, mode);
            }


            if ( edgeLoad->giveCoordSystMode() == Load :: CST_Global ) {
                fx = cosine * components.at(1) + sine *components.at(2);
                fz = -sine *components.at(1) + cosine *components.at(2);
                fm = components.at(3);

                dfx = cosine * ( components2.at(1) - components.at(1) ) + sine * ( components2.at(2) - components.at(2) );
                dfz = -sine * ( components2.at(1) - components.at(1) ) + cosine * ( components2.at(2) - components.at(2) );
                dfm = components2.at(3) - components.at(3);
            } else {
                fx = components.at(1);
                fz = components.at(2);
                fm = components.at(3);
                dfx = components2.at(1) - components.at(1);
                dfz = components2.at(2) - components.at(2);
                dfm = components2.at(3) - components.at(3);
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
            OOFEM_ERROR("unsupported load type");
        }
    }
}


void
Beam2d :: computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode)
{
    FloatArray lc(1);
    StructuralElement :: computeBodyLoadVectorAt(answer, load, tStep, mode);
    answer.times( this->giveCrossSection()->give(CS_Area, lc, this) );
}


int
Beam2d :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_BeamForceMomentumTensor ) {
        answer = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStressVector();
        return 1;
    } else if ( type == IST_BeamStrainCurvatureTensor ) {
        answer = static_cast< StructuralMaterialStatus * >( gp->giveMaterialStatus() )->giveStrainVector();
        return 1;
    } else {
        return StructuralElement :: giveIPValue(answer, gp, type, tStep);
    }
}


void
Beam2d :: printOutputAt(FILE *File, TimeStep *tStep)
{
    FloatArray rl, Fl;

    fprintf(File, "beam element %d (%8d) :\n", this->giveLabel(), this->giveNumber() );

    // ask for global element displacement vector
    this->computeVectorOf(VM_Total, tStep, rl);

    // ask for global element end forces vector
    this->giveEndForcesVector(Fl, tStep);

    fprintf(File, "  local displacements ");
    for ( auto &val : rl ) {
        fprintf(File, " %.4e", val);
    }

    fprintf(File, "\n  local end forces    ");
    for ( auto &val : Fl ) {
        fprintf(File, " %.4e", val);
    }

    fprintf(File, "\n");

    for ( auto &iRule : integrationRulesArray ) {
        iRule->printOutputAt(File, tStep);
    }
}


void
Beam2d :: computeLocalForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode)
// Computes the load vector of the receiver, at tStep.
{
    StructuralElement :: computeLocalForceLoadVector(answer, tStep, mode);
}


void
Beam2d :: computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass, const double *ipDensity)
{
    // computes mass matrix of the receiver

    FloatMatrix stiff;
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);

    /*
     * StructuralElement::computeMassMatrix(answer, tStep);
     * answer.times(this->giveCrossSection()->give('A'));
     */
    double l = this->computeLength();
    double kappa = this->giveKappaCoeff(tStep);
    double kappa2 = kappa * kappa;

    double density = this->giveStructuralCrossSection()->give('d', gp); // constant density assumed
    if ( ipDensity != NULL ) {
        // Override density if desired
        density = * ipDensity;
    }

    double area = this->giveCrossSection()->give(CS_Area, gp); // constant area assumed
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

    double l = this->computeLength();
    double kappa = this->giveKappaCoeff(tStep);
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

    //answer.beLumpedOf (mass);
}


#ifdef __OOFEG
void Beam2d :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
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


void Beam2d :: drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType type)
{
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    double defScale = gc.getDefScale();
    //  if (!go) { // create new one
    WCRec p [ 2 ]; /* poin */
    EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getDeformedElementColor() );
    EASValsSetLayer(OOFEG_DEFORMED_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 0 ].y = 0.;
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveUpdatedCoordinate(3, tStep, defScale);

    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(1, tStep, defScale);
    p [ 1 ].y = 0.;
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveUpdatedCoordinate(3, tStep, defScale);
    go = CreateLine3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}
#endif
} // end namespace oofem
