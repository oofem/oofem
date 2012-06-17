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

#include "tr1_2d_supg.h"
#include "node.h"
#include "material.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "engngm.h"
#include "fluiddynamicmaterial.h"
#include "load.h"
#include "timestep.h"
#include "boundaryload.h"
#include "geotoolbox.h"
#include "materialinterface.h"
#include "contextioerr.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
#endif

namespace oofem {
#define TRSUPG_ZERO_VOF 1.e-8

FEI2dTrLin TR1_2D_SUPG :: interp(1,2);

TR1_2D_SUPG :: TR1_2D_SUPG(int n, Domain *aDomain) :
    SUPGElement(n, aDomain), LEPlicElementInterface()
    // Constructor.
{
    numberOfDofMans  = 3;
}

TR1_2D_SUPG :: ~TR1_2D_SUPG()
// Destructor
{ }



int
TR1_2D_SUPG :: computeNumberOfDofs(EquationID ut)
{
    if ( ut == EID_MomentumBalance ) {
        return 6;
    } else if ( ut == EID_ConservationEquation ) {
        return 3;
    } else if ( ut == EID_MomentumBalance_ConservationEquation ) {
        return 9;
    } else {
        _error("computeNumberOfDofs: Unknown equation id encountered");
    }

    return 0;
}

void
TR1_2D_SUPG :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    if ( ( ut == EID_MomentumBalance ) || ( ut == EID_AuxMomentumBalance ) ) {
        answer.setValues(2, V_u, V_v);
    } else if ( ut == EID_ConservationEquation ) {
        answer.setValues(1, P_f);
    } else if ( ut == EID_MomentumBalance_ConservationEquation ) {
        answer.setValues(3, V_u, V_v, P_f);
    } else {
        _error("giveDofManDofIDMask: Unknown equation id encountered");
    }
}

void
TR1_2D_SUPG :: giveElementDofIDMask(EquationID ut, IntArray &answer) const
{
    this->giveDofManDofIDMask(1, ut, answer);
}


IRResultType
TR1_2D_SUPG :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;               // Required by IR_GIVE_FIELD macro

    this->SUPGElement :: initializeFrom(ir);

    this->vof = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, vof, IFT_TR12DSUPG_pvof, "pvof");
    if ( vof > 0.0 ) {
        setPermanentVolumeFraction(vof);
        this->temp_vof = vof;
    } else {
        this->vof = 0.0;
        IR_GIVE_OPTIONAL_FIELD(ir, vof, IFT_TR12DSUPG_vof, "vof");
        this->temp_vof = this->vof;
    }

    this->computeGaussPoints();
    this->initGeometry();
    return IRRT_OK;
}

void
TR1_2D_SUPG :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Triangle, 1, _2dFlow);
    }
}


void
TR1_2D_SUPG :: computeAccelerationTerm_MB(FloatMatrix &answer, TimeStep *atTime)
{
    answer.resize(6, 6);
    answer.zero();
    FloatArray un;
    double rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), atTime);

    double ar6 = rho * area / 6.0;
    double ar12 = rho * area / 12.0;
    double usum, vsum, coeff;

    /* consistent mass */

    answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = ar6;
    answer.at(4, 4) = answer.at(5, 5) = answer.at(6, 6) = ar6;

    answer.at(1, 3) = answer.at(1, 5) = ar12;
    answer.at(3, 1) = answer.at(3, 5) = ar12;
    answer.at(5, 1) = answer.at(5, 3) = ar12;

    answer.at(2, 4) = answer.at(2, 6) = ar12;
    answer.at(4, 2) = answer.at(4, 6) = ar12;
    answer.at(6, 2) = answer.at(6, 4) = ar12;

    /* SUPG stabilization term */
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime->givePreviousStep(), un);

    usum = un.at(1) + un.at(3) + un.at(5);
    vsum = un.at(2) + un.at(4) + un.at(6);
    coeff = rho * t_supg * area / 12.;

    answer.at(1, 1) += coeff * ( b [ 0 ] * ( usum + un.at(1) ) + c [ 0 ] * ( vsum + un.at(2) ) );
    answer.at(1, 3) += coeff * ( b [ 0 ] * ( usum + un.at(3) ) + c [ 0 ] * ( vsum + un.at(4) ) );
    answer.at(1, 5) += coeff * ( b [ 0 ] * ( usum + un.at(5) ) + c [ 0 ] * ( vsum + un.at(6) ) );

    answer.at(2, 2) += coeff * ( b [ 0 ] * ( usum + un.at(1) ) + c [ 0 ] * ( vsum + un.at(2) ) );
    answer.at(2, 4) += coeff * ( b [ 0 ] * ( usum + un.at(3) ) + c [ 0 ] * ( vsum + un.at(4) ) );
    answer.at(2, 6) += coeff * ( b [ 0 ] * ( usum + un.at(5) ) + c [ 0 ] * ( vsum + un.at(6) ) );

    answer.at(3, 1) += coeff * ( b [ 1 ] * ( usum + un.at(1) ) + c [ 1 ] * ( vsum + un.at(2) ) );
    answer.at(3, 3) += coeff * ( b [ 1 ] * ( usum + un.at(3) ) + c [ 1 ] * ( vsum + un.at(4) ) );
    answer.at(3, 5) += coeff * ( b [ 1 ] * ( usum + un.at(5) ) + c [ 1 ] * ( vsum + un.at(6) ) );

    answer.at(4, 2) += coeff * ( b [ 1 ] * ( usum + un.at(1) ) + c [ 1 ] * ( vsum + un.at(2) ) );
    answer.at(4, 4) += coeff * ( b [ 1 ] * ( usum + un.at(3) ) + c [ 1 ] * ( vsum + un.at(4) ) );
    answer.at(4, 6) += coeff * ( b [ 1 ] * ( usum + un.at(5) ) + c [ 1 ] * ( vsum + un.at(6) ) );

    answer.at(5, 1) += coeff * ( b [ 2 ] * ( usum + un.at(1) ) + c [ 2 ] * ( vsum + un.at(2) ) );
    answer.at(5, 3) += coeff * ( b [ 2 ] * ( usum + un.at(3) ) + c [ 2 ] * ( vsum + un.at(4) ) );
    answer.at(5, 5) += coeff * ( b [ 2 ] * ( usum + un.at(5) ) + c [ 2 ] * ( vsum + un.at(6) ) );

    answer.at(6, 2) += coeff * ( b [ 2 ] * ( usum + un.at(1) ) + c [ 2 ] * ( vsum + un.at(2) ) );
    answer.at(6, 4) += coeff * ( b [ 2 ] * ( usum + un.at(3) ) + c [ 2 ] * ( vsum + un.at(4) ) );
    answer.at(6, 6) += coeff * ( b [ 2 ] * ( usum + un.at(5) ) + c [ 2 ] * ( vsum + un.at(6) ) );
}


void
TR1_2D_SUPG :: computeAdvectionTerm_MB(FloatArray &answer, TimeStep *atTime)
{
    answer.resize(6);
    answer.zero();

    FloatArray u, un;
    double rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), atTime);
    double dudx, dudy, dvdx, dvdy, usum, vsum, coeff;
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime->givePreviousStep(), un);
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);

    dudx = b [ 0 ] * u.at(1) + b [ 1 ] * u.at(3) + b [ 2 ] * u.at(5);
    dudy = c [ 0 ] * u.at(1) + c [ 1 ] * u.at(3) + c [ 2 ] * u.at(5);
    dvdx = b [ 0 ] * u.at(2) + b [ 1 ] * u.at(4) + b [ 2 ] * u.at(6);
    dvdy = c [ 0 ] * u.at(2) + c [ 1 ] * u.at(4) + c [ 2 ] * u.at(6);

    usum = un.at(1) + un.at(3) + un.at(5);
    vsum = un.at(2) + un.at(4) + un.at(6);
    coeff = rho * area / 12.0;

    // standard galerkin term
    answer.at(1) = coeff * ( dudx * ( usum + un.at(1) ) + dudy * ( vsum + un.at(2) ) );
    answer.at(3) = coeff * ( dudx * ( usum + un.at(3) ) + dudy * ( vsum + un.at(4) ) );
    answer.at(5) = coeff * ( dudx * ( usum + un.at(5) ) + dudy * ( vsum + un.at(6) ) );
    answer.at(2) = coeff * ( dvdx * ( usum + un.at(1) ) + dvdy * ( vsum + un.at(2) ) );
    answer.at(4) = coeff * ( dvdx * ( usum + un.at(3) ) + dvdy * ( vsum + un.at(4) ) );
    answer.at(6) = coeff * ( dvdx * ( usum + un.at(5) ) + dvdy * ( vsum + un.at(6) ) );

    // supg stabilization term
    coeff = t_supg * rho * area / 12.0;
    double u1u1 = usum * usum + un.at(1) * un.at(1) + un.at(3) * un.at(3) + un.at(5) * un.at(5);
    double u1u2 = usum * vsum + un.at(1) * un.at(2) + un.at(3) * un.at(4) + un.at(5) * un.at(6);
    double u2u2 = vsum * vsum + un.at(2) * un.at(2) + un.at(4) * un.at(4) + un.at(6) * un.at(6);
    answer.at(1) += coeff * ( b [ 0 ] * ( dudx * u1u1 + dudy * u1u2 ) + c [ 0 ] * ( dudx * u1u2 + dudy * u2u2 ) );
    answer.at(3) += coeff * ( b [ 1 ] * ( dudx * u1u1 + dudy * u1u2 ) + c [ 1 ] * ( dudx * u1u2 + dudy * u2u2 ) );
    answer.at(5) += coeff * ( b [ 2 ] * ( dudx * u1u1 + dudy * u1u2 ) + c [ 2 ] * ( dudx * u1u2 + dudy * u2u2 ) );

    answer.at(2) += coeff * ( b [ 0 ] * ( dvdx * u1u1 + dvdy * u1u2 ) + c [ 0 ] * ( dvdx * u1u2 + dvdy * u2u2 ) );
    answer.at(4) += coeff * ( b [ 1 ] * ( dvdx * u1u1 + dvdy * u1u2 ) + c [ 1 ] * ( dvdx * u1u2 + dvdy * u2u2 ) );
    answer.at(6) += coeff * ( b [ 2 ] * ( dvdx * u1u1 + dvdy * u1u2 ) + c [ 2 ] * ( dvdx * u1u2 + dvdy * u2u2 ) );
}


void
TR1_2D_SUPG :: computeAdvectionDerivativeTerm_MB(FloatMatrix &answer, TimeStep *atTime)
{
    answer.resize(6, 6);
    answer.zero();

    FloatArray u, un;
    double rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), atTime);
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime->givePreviousStep(), un);

    double dudx [ 2 ] [ 2 ], usum [ 2 ];
    double coeff, ar12 = area / 12.;
    int i, j, k, m, w_dof_addr, u_dof_addr, d1j, d2j, dij;

    dudx [ 0 ] [ 0 ] = b [ 0 ] * u.at(1) + b [ 1 ] * u.at(3) + b [ 2 ] * u.at(5);
    dudx [ 0 ] [ 1 ] = c [ 0 ] * u.at(1) + c [ 1 ] * u.at(3) + c [ 2 ] * u.at(5);
    dudx [ 1 ] [ 0 ] = b [ 0 ] * u.at(2) + b [ 1 ] * u.at(4) + b [ 2 ] * u.at(6);
    dudx [ 1 ] [ 1 ] = c [ 0 ] * u.at(2) + c [ 1 ] * u.at(4) + c [ 2 ] * u.at(6);
    usum [ 0 ] = un.at(1) + un.at(3) + un.at(5);
    usum [ 1 ] = un.at(2) + un.at(4) + un.at(6);

    // dN(v)/dv
    for ( i = 1; i <= 2; i++ ) { // test function index
        for ( k = 1; k <= 3; k++ ) { // nodal val of function w
            for ( j = 1; j <= 2; j++ ) { // velocity vector component
                for ( m = 1; m <= 3; m++ ) { //  nodal components
                    w_dof_addr = ( k - 1 ) * 2 + i;
                    u_dof_addr = ( m - 1 ) * 2 + j;
                    d1j = ( j == 1 );
                    d2j = ( j == 2 );
                    dij = ( i == j );
                    coeff = ( m == k ) ? area / 6. : area / 12.;
                    answer.at(w_dof_addr, u_dof_addr) = rho * ( 0.0 * d1j * dudx [ i - 1 ] [ 0 ] * coeff + dij * b [ m - 1 ] * ar12 * ( usum [ 0 ] + un.at( ( k - 1 ) * 2 + 1 ) ) +
                                                               0.0 * d2j * dudx [ i - 1 ] [ 1 ] * coeff + dij * c [ m - 1 ] * ar12 * ( usum [ 1 ] + un.at( ( k - 1 ) * 2 + 2 ) ) );
                }
            }
        }
    }

    // stabilization term dN_delta/du
    for ( i = 1; i <= 2; i++ ) { // test function index
        for ( k = 1; k <= 3; k++ ) { // nodal val of function w
            for ( j = 1; j <= 2; j++ ) { // velocity vector component
                for ( m = 1; m <= 3; m++ ) { //  nodal components
                    w_dof_addr = ( k - 1 ) * 2 + i;
                    u_dof_addr = ( m - 1 ) * 2 + j;
                    d1j = ( j == 1 );
                    d2j = ( j == 2 );
                    dij = ( i == j );
                    answer.at(w_dof_addr, u_dof_addr) += t_supg * rho *
                                                         (
                        0.0 * d1j * b [ k - 1 ] * dudx [ i - 1 ] [ 0 ] * ar12 * ( usum [ 0 ] + un.at( ( m - 1 ) * 2 + 1 ) ) +
                        0.0 * d1j * b [ k - 1 ] * dudx [ i - 1 ] [ 1 ] * ar12 * ( usum [ 1 ] + un.at( ( m - 1 ) * 2 + 2 ) ) +
                        0.0 * d2j * c [ k - 1 ] * dudx [ i - 1 ] [ 0 ] * ar12 * ( usum [ 0 ] + un.at( ( m - 1 ) * 2 + 1 ) ) +
                        0.0 * d2j * c [ k - 1 ] * dudx [ i - 1 ] [ 1 ] * ar12 * ( usum [ 1 ] + un.at( ( m - 1 ) * 2 + 2 ) ) +

                        0.0 * d1j * b [ k - 1 ] * dudx [ i - 1 ] [ 0 ] * ar12 * ( usum [ 0 ] + un.at( ( m - 1 ) * 2 + 1 ) ) +
                        dij * b [ k - 1 ] * b [ m - 1 ] * ar12 * ( usum [ 0 ] * usum [ 0 ] + un.at(1) * un.at(1) + un.at(3) * un.at(3) + un.at(5) * un.at(5) ) +
                        0.0 * d2j * b [ k - 1 ] * dudx [ i - 1 ] [ 1 ] * ar12 * ( usum [ 0 ] + un.at( ( m - 1 ) * 2 + 1 ) ) +
                        dij * b [ k - 1 ] * c [ m - 1 ] * ar12 * ( usum [ 0 ] * usum [ 1 ] + un.at(1) * un.at(2) + un.at(3) * un.at(4) + un.at(5) * un.at(6) ) +

                        0.0 * d1j * c [ k - 1 ] * dudx [ i - 1 ] [ 0 ] * ar12 * ( usum [ 1 ] + un.at( ( m - 1 ) * 2 + 2 ) ) +
                        dij * c [ k - 1 ] * b [ m - 1 ] * ar12 * ( usum [ 0 ] * usum [ 1 ] + un.at(1) * un.at(2) + un.at(3) * un.at(4) + un.at(5) * un.at(6) ) +
                        0.0 * d2j * c [ k - 1 ] * dudx [ i - 1 ] [ 1 ] * ar12 * ( usum [ 1 ] + un.at( ( m - 1 ) * 2 + 2 ) ) +
                        dij * c [ k - 1 ] * c [ m - 1 ] * ar12 * ( usum [ 1 ] * usum [ 1 ] + un.at(2) * un.at(2) + un.at(4) * un.at(4) + un.at(6) * un.at(6) )
                                                         );
                }
            }
        }
    }
}


void
TR1_2D_SUPG :: computeDiffusionTerm_MB(FloatArray &answer, TimeStep *atTime)
{
    int i;
    answer.resize(6);
    answer.zero();
    FloatArray u, eps(3), stress;
    double Re = domain->giveEngngModel()->giveUnknownComponent(ReynoldsNumber, VM_Unknown, atTime, domain, NULL);

    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);

    eps.at(1) = ( b [ 0 ] * u.at(1) + b [ 1 ] * u.at(3) + b [ 2 ] * u.at(5) );
    eps.at(2) = ( c [ 0 ] * u.at(2) + c [ 1 ] * u.at(4) + c [ 2 ] * u.at(6) );
    eps.at(3) = ( b [ 0 ] * u.at(2) + b [ 1 ] * u.at(4) + b [ 2 ] * u.at(6) + c [ 0 ] * u.at(1) + c [ 1 ] * u.at(3) + c [ 2 ] * u.at(5) );
    ( ( FluidDynamicMaterial * ) this->giveMaterial() )->computeDeviatoricStressVector(stress, integrationRulesArray [ 0 ]->getIntegrationPoint(0),
                                                                                       eps, atTime);
    stress.times(1. / Re);

    // \int dNu/dxj \Tau_ij
    answer.resize(6);
    for ( i = 0; i < 3; i++ ) {
        //rh1p(1,lok) = -geome(7,ia)*( sigxx(ia)*geome(lok,ia) + sigxy(ia)*geome(lok1,ia) )*0.5d+00;
        //rh1p(2,lok) = -geome(7,ia)*( sigxy(ia)*geome(lok,ia) + sigyy(ia)*geome(lok1,ia) )*0.5d+00;
        answer.at( ( i ) * 2 + 1 ) = area * ( stress.at(1) * b [ i ] + stress.at(3) * c [ i ] );
        answer.at( ( i + 1 ) * 2 ) = area * ( stress.at(3) * b [ i ] + stress.at(2) * c [ i ] );
    }
}


void
TR1_2D_SUPG :: computeDiffusionDerivativeTerm_MB(FloatMatrix &answer, MatResponseMode mode, TimeStep *atTime)
{
    answer.resize(6, 6);
    answer.zero();
    FloatMatrix _db, _d, _b(3, 6);
    double Re = domain->giveEngngModel()->giveUnknownComponent(ReynoldsNumber, VM_Unknown, atTime, domain, NULL);

    _b.at(1, 1) = b [ 0 ];
    _b.at(1, 2) = 0.;
    _b.at(1, 3) = b [ 1 ];
    _b.at(1, 4) = 0.;
    _b.at(1, 5) = b [ 2 ];
    _b.at(1, 6) = 0.;
    _b.at(2, 1) = 0.;
    _b.at(2, 2) = c [ 0 ];
    _b.at(2, 3) = 0.;
    _b.at(2, 4) = c [ 1 ];
    _b.at(2, 5) = 0.;
    _b.at(2, 6) = c [ 2 ];
    _b.at(3, 1) = c [ 0 ];
    _b.at(3, 2) = b [ 0 ];
    _b.at(3, 3) = c [ 1 ];
    _b.at(3, 4) = b [ 1 ];
    _b.at(3, 5) = c [ 2 ];
    _b.at(3, 6) = b [ 2 ];

    ( ( FluidDynamicMaterial * ) this->giveMaterial() )->giveDeviatoricStiffnessMatrix(_d, mode,
                                                                                       integrationRulesArray [ 0 ]->getIntegrationPoint(0),
                                                                                       atTime);
    _db.beProductOf(_d, _b);
    answer.plusProductUnsym(_b, _db, area); //answer.plusProduct (_b,_db,area);
    //answer.symmetrized() ;
    answer.times(1. / Re);
}

void
TR1_2D_SUPG :: computePressureTerm_MB(FloatMatrix &answer, TimeStep *atTime)
{
    answer.resize(6, 3);
    answer.zero();
    FloatArray p, un;
    double usum, vsum;
    double ar3 = area / 3.0, coeff;

    this->computeVectorOf(EID_ConservationEquation, VM_Total, atTime, p);

    // G matrix
    answer.at(1, 1) = answer.at(1, 2) = answer.at(1, 3) = -b [ 0 ] * ar3;
    answer.at(3, 1) = answer.at(3, 2) = answer.at(3, 3) = -b [ 1 ] * ar3;
    answer.at(5, 1) = answer.at(5, 2) = answer.at(5, 3) = -b [ 2 ] * ar3;

    answer.at(2, 1) = answer.at(2, 2) = answer.at(2, 3) = -c [ 0 ] * ar3;
    answer.at(4, 1) = answer.at(4, 2) = answer.at(4, 3) = -c [ 1 ] * ar3;
    answer.at(6, 1) = answer.at(6, 2) = answer.at(6, 3) = -c [ 2 ] * ar3;

    // stabilization term (G_\delta mtrx)
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime->givePreviousStep(), un);
    usum = un.at(1) + un.at(3) + un.at(5);
    vsum = un.at(2) + un.at(4) + un.at(6);
    coeff = ar3 * t_supg;

    answer.at(1, 1) += coeff * ( usum * b [ 0 ] * b [ 0 ] + vsum * c [ 0 ] * b [ 0 ] );
    answer.at(1, 2) += coeff * ( usum * b [ 0 ] * b [ 1 ] + vsum * c [ 0 ] * b [ 1 ] );
    answer.at(1, 3) += coeff * ( usum * b [ 0 ] * b [ 2 ] + vsum * c [ 0 ] * b [ 2 ] );

    answer.at(3, 1) += coeff * ( usum * b [ 1 ] * b [ 0 ] + vsum * c [ 1 ] * b [ 0 ] );
    answer.at(3, 2) += coeff * ( usum * b [ 1 ] * b [ 1 ] + vsum * c [ 1 ] * b [ 1 ] );
    answer.at(3, 3) += coeff * ( usum * b [ 1 ] * b [ 2 ] + vsum * c [ 1 ] * b [ 2 ] );

    answer.at(5, 1) += coeff * ( usum * b [ 2 ] * b [ 0 ] + vsum * c [ 2 ] * b [ 0 ] );
    answer.at(5, 2) += coeff * ( usum * b [ 2 ] * b [ 1 ] + vsum * c [ 2 ] * b [ 1 ] );
    answer.at(5, 3) += coeff * ( usum * b [ 2 ] * b [ 2 ] + vsum * c [ 2 ] * b [ 2 ] );

    answer.at(2, 1) += coeff * ( usum * b [ 0 ] * c [ 0 ] + vsum * c [ 0 ] * c [ 0 ] );
    answer.at(2, 2) += coeff * ( usum * b [ 0 ] * c [ 1 ] + vsum * c [ 0 ] * c [ 1 ] );
    answer.at(2, 3) += coeff * ( usum * b [ 0 ] * c [ 2 ] + vsum * c [ 0 ] * c [ 2 ] );

    answer.at(4, 1) += coeff * ( usum * b [ 1 ] * c [ 0 ] + vsum * c [ 1 ] * c [ 0 ] );
    answer.at(4, 2) += coeff * ( usum * b [ 1 ] * c [ 1 ] + vsum * c [ 1 ] * c [ 1 ] );
    answer.at(4, 3) += coeff * ( usum * b [ 1 ] * c [ 2 ] + vsum * c [ 1 ] * c [ 2 ] );

    answer.at(6, 1) += coeff * ( usum * b [ 2 ] * c [ 0 ] + vsum * c [ 2 ] * c [ 0 ] );
    answer.at(6, 2) += coeff * ( usum * b [ 2 ] * c [ 1 ] + vsum * c [ 2 ] * c [ 1 ] );
    answer.at(6, 3) += coeff * ( usum * b [ 2 ] * c [ 2 ] + vsum * c [ 2 ] * c [ 2 ] );
}


void
TR1_2D_SUPG :: computeLSICStabilizationTerm_MB(FloatMatrix &answer, TimeStep *atTime)
{
    answer.resize(6, 6);
    answer.zero();
    double rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), atTime);
    double coeff = area * t_lsic * rho;
    double n[] = {
        b [ 0 ], c [ 0 ], b [ 1 ], c [ 1 ], b [ 2 ], c [ 2 ]
    };
    int i, j;

    for ( i = 1; i <= 6; i++ ) {
        for ( j = 1; j <= 6; j++ ) {
            answer.at(i, j) = coeff * n [ i - 1 ] * n [ j - 1 ];
        }
    }
}


void
TR1_2D_SUPG :: computeLinearAdvectionTerm_MC(FloatMatrix &answer, TimeStep *atTime)
{
    answer.resize(3, 6);
    answer.zero();
    double ar3 = area / 3.0;

    // G^T matrix
    answer.at(1, 1) = answer.at(2, 1) = answer.at(3, 1) = b [ 0 ] * ar3;
    answer.at(1, 3) = answer.at(2, 3) = answer.at(3, 3) = b [ 1 ] * ar3;
    answer.at(1, 5) = answer.at(2, 5) = answer.at(3, 5) = b [ 2 ] * ar3;

    answer.at(1, 2) = answer.at(2, 2) = answer.at(3, 2) = c [ 0 ] * ar3;
    answer.at(1, 4) = answer.at(2, 4) = answer.at(3, 4) = c [ 1 ] * ar3;
    answer.at(1, 6) = answer.at(2, 6) = answer.at(3, 6) = c [ 2 ] * ar3;
}

void
TR1_2D_SUPG :: computeAdvectionTerm_MC(FloatArray &answer, TimeStep *atTime)
{
    // N_epsilon (due to PSPG stabilization)
    double coeff = t_pspg * area / 3.0;
    double dudx, dudy, dvdx, dvdy, usum, vsum;
    FloatArray u, un;

    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime->givePreviousStep(), un);
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);

    dudx = b [ 0 ] * u.at(1) + b [ 1 ] * u.at(3) + b [ 2 ] * u.at(5);
    dudy = c [ 0 ] * u.at(1) + c [ 1 ] * u.at(3) + c [ 2 ] * u.at(5);
    dvdx = b [ 0 ] * u.at(2) + b [ 1 ] * u.at(4) + b [ 2 ] * u.at(6);
    dvdy = c [ 0 ] * u.at(2) + c [ 1 ] * u.at(4) + c [ 2 ] * u.at(6);

    usum = un.at(1) + un.at(3) + un.at(5);
    vsum = un.at(2) + un.at(4) + un.at(6);

    answer.resize(3);

    answer.at(1) = coeff * ( b [ 0 ] * ( dudx * usum + dudy * vsum ) + c [ 0 ] * ( dvdx * usum + dvdy * vsum ) );
    answer.at(2) = coeff * ( b [ 1 ] * ( dudx * usum + dudy * vsum ) + c [ 1 ] * ( dvdx * usum + dvdy * vsum ) );
    answer.at(3) = coeff * ( b [ 2 ] * ( dudx * usum + dudy * vsum ) + c [ 2 ] * ( dvdx * usum + dvdy * vsum ) );
}


void
TR1_2D_SUPG :: computeAdvectionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *atTime)
{
    answer.resize(3, 6);
    answer.zero();
    int k, j, m, w_dof_addr, u_dof_addr, d1j, d2j, km1, mm1;
    //double rho = this->giveMaterial()->giveCharacteristicValue(Density, integrationRulesArray[0]->getIntegrationPoint(0), atTime);
    FloatArray u, un;

    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime->givePreviousStep(), un);

    double dudx [ 2 ] [ 2 ], usum [ 2 ];
    double coeff;

    dudx [ 0 ] [ 0 ] = b [ 0 ] * u.at(1) + b [ 1 ] * u.at(3) + b [ 2 ] * u.at(5);
    dudx [ 0 ] [ 1 ] = c [ 0 ] * u.at(1) + c [ 1 ] * u.at(3) + c [ 2 ] * u.at(5);
    dudx [ 1 ] [ 0 ] = b [ 0 ] * u.at(2) + b [ 1 ] * u.at(4) + b [ 2 ] * u.at(6);
    dudx [ 1 ] [ 1 ] = c [ 0 ] * u.at(2) + c [ 1 ] * u.at(4) + c [ 2 ] * u.at(6);
    usum [ 0 ] = un.at(1) + un.at(3) + un.at(5);
    usum [ 1 ] = un.at(2) + un.at(4) + un.at(6);

    // dN_epsilon(v)/dv
    coeff = t_pspg * area / 3.;
    for ( k = 1; k <= 3; k++ ) { // nodal val of function w
        km1 = k - 1;
        for ( j = 1; j <= 2; j++ ) { // velocity vector component
            for ( m = 1; m <= 3; m++ ) { //  nodal components
                w_dof_addr = k;
                u_dof_addr = ( m - 1 ) * 2 + j;
                mm1 = m - 1;
                d1j = ( j == 1 );
                d2j = ( j == 2 );
                answer.at(w_dof_addr, u_dof_addr) = coeff * ( 0.0 * d1j * b [ km1 ] * dudx [ 0 ] [ 0 ] + d1j * b [ km1 ] * b [ mm1 ] * usum [ 0 ] +
                                                              0.0 * d2j * b [ km1 ] * dudx [ 0 ] [ 1 ] + d1j * b [ km1 ] * c [ mm1 ] * usum [ 1 ] +
                                                              0.0 * d1j * c [ km1 ] * dudx [ 1 ] [ 0 ] + d2j * c [ km1 ] * b [ mm1 ] * usum [ 0 ] +
                                                              0.0 * d2j * c [ km1 ] * dudx [ 1 ] [ 1 ] + d2j * c [ km1 ] * c [ mm1 ] * usum [ 1 ] );
            }
        }
    }
}

void
TR1_2D_SUPG :: computeAccelerationTerm_MC(FloatMatrix &answer, TimeStep *atTime)
{
    answer.resize(3, 6);
    answer.zero();
    double coeff = t_pspg * area / 3.0;
    // M_\epsilon

    answer.at(1, 1) = answer.at(1, 3) = answer.at(1, 5) = coeff * b [ 0 ];
    answer.at(1, 2) = answer.at(1, 4) = answer.at(1, 6) = coeff * c [ 0 ];
    answer.at(2, 1) = answer.at(2, 3) = answer.at(2, 5) = coeff * b [ 1 ];
    answer.at(2, 2) = answer.at(2, 4) = answer.at(2, 6) = coeff * c [ 1 ];
    answer.at(3, 1) = answer.at(3, 3) = answer.at(3, 5) = coeff * b [ 2 ];
    answer.at(3, 2) = answer.at(3, 4) = answer.at(3, 6) = coeff * c [ 2 ];
}

void
TR1_2D_SUPG :: computePressureTerm_MC(FloatMatrix &answer, TimeStep *atTime)
{
    double rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), atTime);
    int i, j;
    double coeff = t_pspg * area / rho;
    answer.resize(3, 3);


    for ( i = 1; i <= 3; i++ ) {
        for ( j = 1; j <= 3; j++ ) {
            answer.at(i, j) = coeff * ( b [ i - 1 ] * b [ j - 1 ] + c [ i - 1 ] * c [ j - 1 ] );
        }
    }
}


void
TR1_2D_SUPG :: computeSlipWithFrictionBCTerm_MB(FloatMatrix &answer, Load *load, int side, TimeStep *atTime)
{
    //answer.resize(6, 6);
    //answer.zero();

    int node1, node2, node3;
    int d1 = 1;
    int d2 = 1;
    int d3 = 1;
    double l, t1, t2, _t1, _t2;
    double beta;

    answer.resize(6, 6);
    answer.zero();
    //beta
    //area
    BoundaryLoad *edgeLoad = static_cast< BoundaryLoad * >(load);
    beta = edgeLoad->giveProperty('a');
    node1 = side;
    node2 = ( node1 == 3 ? 1 : node1 + 1 );

    node3 = ( node2 == 3 ? 1 : node2 + 1 );

    switch ( node3 ) {
    case 1:
        d1 = 0;
    case 2:
        d2 = 0;
    case 3:
        d3 = 0;
    }

    _t1 = giveNode(node2)->giveCoordinate(1) - giveNode(node1)->giveCoordinate(1);
    _t2 = giveNode(node2)->giveCoordinate(2) - giveNode(node1)->giveCoordinate(2);
    l = sqrt(_t1 * _t1 + _t2 * _t2);

    t1 = _t1 / l;
    t2 = _t2 / l;

    answer.at(1, 1) += ( l / 3 ) * beta * d1 * t1 * t1;
    answer.at(1, 2) += ( l / 3 ) * beta * d1 * t1 * t2;
    answer.at(2, 1) += ( l / 3 ) * beta * d1 * t1 * t2;
    answer.at(2, 2) += ( l / 3 ) * beta * d1 * t2 * t2;

    answer.at(1, 3) += ( l / 6 ) * beta * d1 * d2 * t1 * t1;
    answer.at(1, 4) += ( l / 6 ) * beta * d1 * d2 * t1 * t2;
    answer.at(2, 3) += ( l / 6 ) * beta * d1 * d2 * t1 * t2;
    answer.at(2, 4) += ( l / 6 ) * beta * d1 * d2 * t2 * t2;

    answer.at(1, 5) += ( l / 6 ) * beta * d1 * d3 * t1 * t1;
    answer.at(1, 6) += ( l / 6 ) * beta * d1 * d3 * t1 * t2;
    answer.at(2, 5) += ( l / 6 ) * beta * d1 * d3 * t1 * t2;
    answer.at(2, 6) += ( l / 6 ) * beta * d1 * d3 * t2 * t2;

    answer.at(3, 1) += ( l / 6 ) * beta * d1 * d2 * t1 * t1;
    answer.at(3, 2) += ( l / 6 ) * beta * d1 * d2 * t1 * t2;
    answer.at(4, 1) += ( l / 6 ) * beta * d1 * d2 * t1 * t2;
    answer.at(4, 2) += ( l / 6 ) * beta * d1 * d2 * t2 * t2;

    answer.at(3, 3) += ( l / 3 ) * beta * d2 * t1 * t1;
    answer.at(3, 4) += ( l / 3 ) * beta * d2 * t1 * t2;
    answer.at(4, 3) += ( l / 3 ) * beta * d2 * t1 * t2;
    answer.at(4, 4) += ( l / 3 ) * beta * d2 * t2 * t2;

    answer.at(3, 5) += ( l / 6 ) * beta * d2 * d3 * t1 * t1;
    answer.at(3, 6) += ( l / 6 ) * beta * d2 * d3 * t1 * t2;
    answer.at(4, 5) += ( l / 6 ) * beta * d2 * d3 * t1 * t2;
    answer.at(4, 6) += ( l / 6 ) * beta * d2 * d3 * t2 * t2;

    answer.at(5, 1) += ( l / 6 ) * beta * d3 * d1 * t1 * t1;
    answer.at(5, 2) += ( l / 6 ) * beta * d3 * d1 * t1 * t2;
    answer.at(6, 1) += ( l / 6 ) * beta * d3 * d1 * t1 * t2;
    answer.at(6, 2) += ( l / 6 ) * beta * d3 * d1 * t2 * t2;

    answer.at(5, 3) += ( l / 6 ) * beta * d3 * d2 * t1 * t1;
    answer.at(5, 4) += ( l / 6 ) * beta * d3 * d2 * t1 * t2;
    answer.at(6, 3) += ( l / 6 ) * beta * d3 * d2 * t1 * t2;
    answer.at(6, 4) += ( l / 6 ) * beta * d3 * d2 * t2 * t2;

    answer.at(5, 5) += ( l / 3 ) * beta * d3 * t1 * t1;
    answer.at(5, 6) += ( l / 3 ) * beta * d3 * t1 * t2;
    answer.at(6, 5) += ( l / 3 ) * beta * d3 * t1 * t2;
    answer.at(6, 6) += ( l / 3 ) * beta * d3 * t2 * t2;

    //answer.negated();
}


void
TR1_2D_SUPG :: computePenetrationWithResistanceBCTerm_MB(FloatMatrix &answer, Load *load, int side, TimeStep *atTime)
{
    //answer.resize(6, 6);
    //answer.zero();

    int node1, node2, node3;
    int d1 = 1;
    int d2 = 1;
    int d3 = 1;
    double l, n1, n2, _t1, _t2;
    double alpha;

    answer.resize(6, 6);
    answer.zero();

    BoundaryLoad *edgeLoad = static_cast< BoundaryLoad * >(load);
    alpha = edgeLoad->giveProperty('a');
    node1 = side;
    node2 = ( node1 == 3 ? 1 : node1 + 1 );

    node3 = ( node2 == 3 ? 1 : node2 + 1 );

    switch ( node3 ) {
    case 1:
        d1 = 0;
    case 2:
        d2 = 0;
    case 3:
        d3 = 0;
    }

    _t1 = giveNode(node2)->giveCoordinate(1) - giveNode(node1)->giveCoordinate(1);
    _t2 = giveNode(node2)->giveCoordinate(2) - giveNode(node1)->giveCoordinate(2);
    l = sqrt(_t1 * _t1 + _t2 * _t2);

    //t1 = _t1 / l;
    //t2 = _t2 / l;

    n1 = _t2 / l;
    n2 = -_t1 / l;

    answer.at(1, 1) += ( l / 3 ) * ( 1 / alpha ) * d1 * n1 * n1;
    answer.at(1, 2) += ( l / 3 ) * ( 1 / alpha ) * d1 * n1 * n2;
    answer.at(2, 1) += ( l / 3 ) * ( 1 / alpha ) * d1 * n1 * n2;
    answer.at(2, 2) += ( l / 3 ) * ( 1 / alpha ) * d1 * n2 * n2;

    answer.at(1, 3) += ( l / 6 ) * ( 1 / alpha ) * d1 * d2 * n1 * n1;
    answer.at(1, 4) += ( l / 6 ) * ( 1 / alpha ) * d1 * d2 * n1 * n2;
    answer.at(2, 3) += ( l / 6 ) * ( 1 / alpha ) * d1 * d2 * n1 * n2;
    answer.at(2, 4) += ( l / 6 ) * ( 1 / alpha ) * d1 * d2 * n2 * n2;

    answer.at(1, 5) += ( l / 6 ) * ( 1 / alpha ) * d1 * d3 * n1 * n1;
    answer.at(1, 6) += ( l / 6 ) * ( 1 / alpha ) * d1 * d3 * n1 * n2;
    answer.at(2, 5) += ( l / 6 ) * ( 1 / alpha ) * d1 * d3 * n1 * n2;
    answer.at(2, 6) += ( l / 6 ) * ( 1 / alpha ) * d1 * d3 * n2 * n2;

    answer.at(3, 1) += ( l / 6 ) * ( 1 / alpha ) * d1 * d2 * n1 * n1;
    answer.at(3, 2) += ( l / 6 ) * ( 1 / alpha ) * d1 * d2 * n1 * n2;
    answer.at(4, 1) += ( l / 6 ) * ( 1 / alpha ) * d1 * d2 * n1 * n2;
    answer.at(4, 2) += ( l / 6 ) * ( 1 / alpha ) * d1 * d2 * n2 * n2;

    answer.at(3, 3) += ( l / 3 ) * ( 1 / alpha ) * d2 * n1 * n1;
    answer.at(3, 4) += ( l / 3 ) * ( 1 / alpha ) * d2 * n1 * n2;
    answer.at(4, 3) += ( l / 3 ) * ( 1 / alpha ) * d2 * n1 * n2;
    answer.at(4, 4) += ( l / 3 ) * ( 1 / alpha ) * d2 * n2 * n2;

    answer.at(3, 5) += ( l / 6 ) * ( 1 / alpha ) * d2 * d3 * n1 * n1;
    answer.at(3, 6) += ( l / 6 ) * ( 1 / alpha ) * d2 * d3 * n1 * n2;
    answer.at(4, 5) += ( l / 6 ) * ( 1 / alpha ) * d2 * d3 * n1 * n2;
    answer.at(4, 6) += ( l / 6 ) * ( 1 / alpha ) * d2 * d3 * n2 * n2;

    answer.at(5, 1) += ( l / 6 ) * ( 1 / alpha ) * d3 * d1 * n1 * n1;
    answer.at(5, 2) += ( l / 6 ) * ( 1 / alpha ) * d3 * d1 * n1 * n2;
    answer.at(6, 1) += ( l / 6 ) * ( 1 / alpha ) * d3 * d1 * n1 * n2;
    answer.at(6, 2) += ( l / 6 ) * ( 1 / alpha ) * d3 * d1 * n2 * n2;

    answer.at(5, 3) += ( l / 6 ) * ( 1 / alpha ) * d3 * d2 * n1 * n1;
    answer.at(5, 4) += ( l / 6 ) * ( 1 / alpha ) * d3 * d2 * n1 * n2;
    answer.at(6, 3) += ( l / 6 ) * ( 1 / alpha ) * d3 * d2 * n1 * n2;
    answer.at(6, 4) += ( l / 6 ) * ( 1 / alpha ) * d3 * d2 * n2 * n2;

    answer.at(5, 5) += ( l / 3 ) * ( 1 / alpha ) * d3 * n1 * n1;
    answer.at(5, 6) += ( l / 3 ) * ( 1 / alpha ) * d3 * n1 * n2;
    answer.at(6, 5) += ( l / 3 ) * ( 1 / alpha ) * d3 * n1 * n2;
    answer.at(6, 6) += ( l / 3 ) * ( 1 / alpha ) * d3 * n2 * n2;

    //answer.negated();
}

void
TR1_2D_SUPG :: computeOutFlowBCTerm_MB(FloatMatrix &answer, int side, TimeStep *atTime)
{
    int node1, node2, node3 ;
    int d1 = 1;
    int d2 = 1;
    int d3 = 1;
    double l, n1, n2, t1, t2;

    answer.resize(6, 3);
    answer.zero();
    //beta
    //area
    node1 = side;
    node2 = ( node1 == 3 ? 1 : node1 + 1 );

    node3 = ( node2 == 3 ? 1 : node2 + 1 );

    switch ( node3 ) {
    case 1:
        d1 = 0;
    case 2:
        d2 = 0;
    case 3:
        d3 = 0;
    }


    t1 = giveNode(node2)->giveCoordinate(1) - giveNode(node1)->giveCoordinate(1);
    t2 = giveNode(node2)->giveCoordinate(2) - giveNode(node1)->giveCoordinate(2);
    l = sqrt(t1 * t1 + t2 * t2);

    n1 = t2 / l;
    n2 = -t1 / l;

    answer.at(1, 1) += ( l / 3 ) * d1 * d1 * n1;
    answer.at(1, 2) += ( l / 6 ) * d1 * d2 * n1;
    answer.at(1, 3) += ( l / 6 ) * d1 * d3 * n1;
    answer.at(2, 1) += ( l / 3 ) * d1 * d1 * n2;
    answer.at(2, 2) += ( l / 6 ) * d1 * d2 * n2;
    answer.at(2, 3) += ( l / 6 ) * d1 * d3 * n2;

    answer.at(3, 1) += ( l / 6 ) * d2 * d1 * n1;
    answer.at(3, 2) += ( l / 3 ) * d2 * d2 * n1;
    answer.at(3, 3) += ( l / 6 ) * d2 * d3 * n1;
    answer.at(4, 1) += ( l / 6 ) * d2 * d1 * n2;
    answer.at(4, 2) += ( l / 3 ) * d2 * d2 * n2;
    answer.at(4, 3) += ( l / 6 ) * d2 * d3 * n2;

    answer.at(5, 1) += ( l / 6 ) * d3 * d1 * n1;
    answer.at(5, 2) += ( l / 6 ) * d3 * d2 * n1;
    answer.at(5, 3) += ( l / 3 ) * d3 * d3 * n1;
    answer.at(6, 1) += ( l / 6 ) * d3 * d1 * n2;
    answer.at(6, 2) += ( l / 6 ) * d3 * d2 * n2;
    answer.at(6, 3) += ( l / 3 ) * d3 * d3 * n2;

    answer.negated();
}


void
TR1_2D_SUPG :: computeBCRhsTerm_MB(FloatArray &answer, TimeStep *atTime)
{
    answer.resize(6);
    answer.zero();

    int i, nLoads;
    double usum [ 2 ];
    Load *load;
    bcGeomType ltype;
    double rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), atTime);
    FloatArray un, gVector;

    // add body load (gravity) termms
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime->givePreviousStep(), un);

    usum [ 0 ] = un.at(1) + un.at(3) + un.at(5);
    usum [ 1 ] = un.at(2) + un.at(4) + un.at(6);
    double coeff = rho * area / 3.0;
    nLoads    = this->giveBodyLoadArray()->giveSize();
    for ( i = 1; i <= nLoads; i++ ) {
        load  = domain->giveLoad( bodyLoadArray.at(i) );
        ltype = load->giveBCGeoType();
        if ( ( ltype == BodyLoadBGT ) && ( load->giveBCValType() == ForceLoadBVT ) ) {
            load->computeComponentArrayAt(gVector, atTime, VM_Total);
            if ( gVector.giveSize() ) {
                answer.at(1) += coeff * ( gVector.at(1) * ( 1.0 + t_supg * ( b [ 0 ] * usum [ 0 ] + c [ 0 ] * usum [ 1 ] ) ) );
                answer.at(2) += coeff * ( gVector.at(2) * ( 1.0 + t_supg * ( b [ 0 ] * usum [ 0 ] + c [ 0 ] * usum [ 1 ] ) ) );
                answer.at(3) += coeff * ( gVector.at(1) * ( 1.0 + t_supg * ( b [ 1 ] * usum [ 0 ] + c [ 1 ] * usum [ 1 ] ) ) );
                answer.at(4) += coeff * ( gVector.at(2) * ( 1.0 + t_supg * ( b [ 1 ] * usum [ 0 ] + c [ 1 ] * usum [ 1 ] ) ) );
                answer.at(5) += coeff * ( gVector.at(1) * ( 1.0 + t_supg * ( b [ 2 ] * usum [ 0 ] + c [ 2 ] * usum [ 1 ] ) ) );
                answer.at(6) += coeff * ( gVector.at(2) * ( 1.0 + t_supg * ( b [ 2 ] * usum [ 0 ] + c [ 2 ] * usum [ 1 ] ) ) );
            }
        }
    }

    // loop over sides
    int n1, n2;
    double tx, ty, l;
    //IntArray nodecounter (3);

    if (  1 ) {
        FloatArray t, coords(1);
        int nLoads, n, id;
        BoundaryLoad *load;
        bcType loadtype;
        // integrate tractions

        // if no traction bc applied but side marked as with traction load
        // then zero traction is assumed !!!

        // loop over boundary load array
        nLoads    = this->giveBoundaryLoadArray()->giveSize() / 2;
        for ( i = 1; i <= nLoads; i++ ) {
            n     = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
            id    = boundaryLoadArray.at(i * 2);

            n1 = id;
            n2 = ( n1 == 3 ? 1 : n1 + 1 );

            tx = giveNode(n2)->giveCoordinate(1) - giveNode(n1)->giveCoordinate(1);
            ty = giveNode(n2)->giveCoordinate(2) - giveNode(n1)->giveCoordinate(2);
            l = sqrt(tx * tx + ty * ty);

            load  = dynamic_cast< BoundaryLoad * >( domain->giveLoad(n) );
            loadtype = load->giveType();
            if ( loadtype == TransmissionBC ) {
                load->computeValueAt(t, atTime, coords, VM_Total);

                // here it is assumed constant traction, one point integration only
                // n1 (u,v)
                answer.at( ( n1 - 1 ) * 2 + 1 ) += t.at(1) * l / 2.;
                answer.at(n1 * 2)       += t.at(2) * l / 2.;
                // n2 (u,v)
                answer.at( ( n2 - 1 ) * 2 + 1 ) += t.at(1) * l / 2.;
                answer.at(n2 * 2)       += t.at(2) * l / 2.;

                //answer.at(n1)+= (t.at(1)*nx + t.at(2)*ny) * l/2.;
                //answer.at(n2)+= (t.at(1)*nx + t.at(2)*ny) * l/2.;
            }
        }
    }
}

void
TR1_2D_SUPG :: computeBCRhsTerm_MC(FloatArray &answer, TimeStep *atTime)
{
    int i, nLoads;
    double coeff;
    Load *load;
    bcGeomType ltype;
    FloatArray gVector;

    answer.resize(3);
    answer.zero();
    coeff = t_pspg * area;
    nLoads    = this->giveBodyLoadArray()->giveSize();
    for ( i = 1; i <= nLoads; i++ ) {
        load  = domain->giveLoad( bodyLoadArray.at(i) );
        ltype = load->giveBCGeoType();
        if ( ( ltype == BodyLoadBGT ) && ( load->giveBCValType() == ForceLoadBVT ) ) {
            load->computeComponentArrayAt(gVector, atTime, VM_Total);
            if ( gVector.giveSize() ) {
                answer.at(1) += coeff * ( b [ 0 ] * gVector.at(1) + c [ 0 ] * gVector.at(2) );
                answer.at(2) += coeff * ( b [ 1 ] * gVector.at(1) + c [ 1 ] * gVector.at(2) );
                answer.at(3) += coeff * ( b [ 2 ] * gVector.at(1) + c [ 2 ] * gVector.at(2) );
            }
        }
    }
}


void
TR1_2D_SUPG :: updateStabilizationCoeffs(TimeStep *atTime)
{
    //TR1_2D_SUPG :: updateStabilizationCoeffs (atTime);
#if 0
    int i, j, k, m, km1, mm1, d1j, d2j, dij, w_dof_addr, u_dof_addr;
    double __g_norm, __gamma_norm, __gammav_norm, __beta_norm, __betav_norm, __c_norm, __ctilda_norm, __e_norm, __k_norm, __Re;
    double __t_p1, __t_p2, __t_p3, __t_s1, __t_s2, __t_s3;
    double nu, rho;
    double dudx [ 2 ] [ 2 ], usum [ 2 ];
    FloatArray u, un, a;

    // compute averaged viscosity based on rule of mixture
    GaussPoint *gp;
    if ( integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints() ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    } else {
        gp = integrationRulesArray [ 1 ]->getIntegrationPoint(0);
    }

    nu = this->giveMaterial()->giveCharacteristicValue( MRM_Viscosity, gp, atTime->givePreviousStep() );
    rho = this->giveMaterial()->giveCharacteristicValue( MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), atTime->givePreviousStep() );

    //this -> computeVectorOf(EID_MomentumBalance,VM_Total,atTime->givePreviousStep(),un) ;
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime->givePreviousStep(), u);
    this->computeVectorOf(EID_MomentumBalance, VM_Acceleration, atTime->givePreviousStep(), a);

    un = u;
    usum [ 0 ] = un.at(1) + un.at(3) + un.at(5);
    usum [ 1 ] = un.at(2) + un.at(4) + un.at(6);

    FloatMatrix __tmp;
    FloatArray __tmpvec, n;
    // assemble g matrix
    __tmp.resize(6, 3);
    double ar3 = area / 3.0;

    __tmp.at(1, 1) = __tmp.at(1, 2) = __tmp.at(1, 3) = b [ 0 ] * ar3;
    __tmp.at(3, 1) = __tmp.at(3, 2) = __tmp.at(3, 3) = b [ 1 ] * ar3;
    __tmp.at(5, 1) = __tmp.at(5, 2) = __tmp.at(5, 3) = b [ 2 ] * ar3;

    __tmp.at(2, 1) = __tmp.at(2, 2) = __tmp.at(2, 3) = c [ 0 ] * ar3;
    __tmp.at(4, 1) = __tmp.at(4, 2) = __tmp.at(4, 3) = c [ 1 ] * ar3;
    __tmp.at(6, 1) = __tmp.at(6, 2) = __tmp.at(6, 3) = c [ 2 ] * ar3;

    __g_norm = __tmp.computeFrobeniusNorm();

    // assemble \gamma matrix (advectionTerm of mass conservation eq)
    __tmp.resize(3, 6);
    dudx [ 0 ] [ 0 ] = b [ 0 ] * u.at(1) + b [ 1 ] * u.at(3) + b [ 2 ] * u.at(5);
    dudx [ 0 ] [ 1 ] = c [ 0 ] * u.at(1) + c [ 1 ] * u.at(3) + c [ 2 ] * u.at(5);
    dudx [ 1 ] [ 0 ] = b [ 0 ] * u.at(2) + b [ 1 ] * u.at(4) + b [ 2 ] * u.at(6);
    dudx [ 1 ] [ 1 ] = c [ 0 ] * u.at(2) + c [ 1 ] * u.at(4) + c [ 2 ] * u.at(6);


    // dN_epsilon(v)/dv
    double coeff = area / 3.;
    for ( k = 1; k <= 3; k++ ) { // nodal val of function w
        km1 = k - 1;
        for ( j = 1; j <= 2; j++ ) { // velocity vector component
            for ( m = 1; m <= 3; m++ ) { //  nodal components
                w_dof_addr = k;
                u_dof_addr = ( m - 1 ) * 2 + j;
                mm1 = m - 1;
                d1j = ( j == 1 );
                d2j = ( j == 2 );
                __tmp.at(w_dof_addr, u_dof_addr) = coeff * ( 0.0 * d1j * b [ km1 ] * dudx [ 0 ] [ 0 ] + d1j * b [ km1 ] * b [ mm1 ] * usum [ 0 ] +
                                                             0.0 * d2j * b [ km1 ] * dudx [ 0 ] [ 1 ] + d1j * b [ km1 ] * c [ mm1 ] * usum [ 1 ] +
                                                             0.0 * d1j * c [ km1 ] * dudx [ 1 ] [ 0 ] + d2j * c [ km1 ] * b [ mm1 ] * usum [ 0 ] +
                                                             0.0 * d2j * c [ km1 ] * dudx [ 1 ] [ 1 ] + d2j * c [ km1 ] * c [ mm1 ] * usum [ 1 ] );
            }
        }
    }

    __gamma_norm = __tmp.computeFrobeniusNorm();
    __tmpvec.beProductOf(__tmp, u);
    __gammav_norm = __tmpvec.computeNorm();
    // compute beta mtrx (acceleration term of mass conservation eq)
    __tmp.resize(3, 6);
    __tmp.zero();
    __tmp.at(1, 1) = __tmp.at(1, 3) = __tmp.at(1, 5) = ar3 * b [ 0 ];
    __tmp.at(1, 2) = __tmp.at(1, 4) = __tmp.at(1, 6) = ar3 * c [ 0 ];
    __tmp.at(2, 1) = __tmp.at(2, 3) = __tmp.at(2, 5) = ar3 * b [ 1 ];
    __tmp.at(2, 2) = __tmp.at(2, 4) = __tmp.at(2, 6) = ar3 * c [ 1 ];
    __tmp.at(3, 1) = __tmp.at(3, 3) = __tmp.at(3, 5) = ar3 * b [ 2 ];
    __tmp.at(3, 2) = __tmp.at(3, 4) = __tmp.at(3, 6) = ar3 * c [ 2 ];
    __beta_norm = __tmp.computeFrobeniusNorm();
    __tmpvec.beProductOf(__tmp, a);
    __betav_norm = __tmpvec.computeNorm();
    // compute c mtrx (advection term of momentum balance)
    // standard galerkin term
    __tmp.resize(6, 6);
    __tmp.zero();

    for ( i = 1; i <= 2; i++ ) { // test function index
        for ( k = 1; k <= 3; k++ ) { // nodal val of function w
            for ( j = 1; j <= 2; j++ ) { // velocity vector component
                for ( m = 1; m <= 3; m++ ) { //  nodal components
                    w_dof_addr = ( k - 1 ) * 2 + i;
                    u_dof_addr = ( m - 1 ) * 2 + j;
                    d1j = ( j == 1 );
                    d2j = ( j == 2 );
                    dij = ( i == j );
                    coeff = ( m == k ) ? area / 6. : area / 12.;
                    __tmp.at(w_dof_addr, u_dof_addr) = rho * ( 0.0 * d1j * dudx [ i - 1 ] [ 0 ] * coeff + dij * b [ m - 1 ] * ( area / 12. ) * ( usum [ 0 ] + un.at( ( k - 1 ) * 2 + 1 ) ) +
                                                              0.0 * d2j * dudx [ i - 1 ] [ 1 ] * coeff + dij * c [ m - 1 ] * ( area / 12. ) * ( usum [ 1 ] + un.at( ( k - 1 ) * 2 + 2 ) ) );
                }
            }
        }
    }

    __c_norm = __tmp.computeFrobeniusNorm();
    // compute ctilda mtrx (stabilization acceleration term of momentum balance)
    __tmp.resize(6, 6);
    __tmp.zero();
    coeff = rho * area / 12.;

    __tmp.at(1, 1) = coeff * ( b [ 0 ] * ( usum [ 0 ] + un.at(1) ) + c [ 0 ] * ( usum [ 1 ] + un.at(2) ) );
    __tmp.at(1, 3) = coeff * ( b [ 0 ] * ( usum [ 0 ] + un.at(3) ) + c [ 0 ] * ( usum [ 1 ] + un.at(4) ) );
    __tmp.at(1, 5) = coeff * ( b [ 0 ] * ( usum [ 0 ] + un.at(5) ) + c [ 0 ] * ( usum [ 1 ] + un.at(6) ) );

    __tmp.at(2, 2) = coeff * ( b [ 0 ] * ( usum [ 0 ] + un.at(1) ) + c [ 0 ] * ( usum [ 1 ] + un.at(2) ) );
    __tmp.at(2, 4) = coeff * ( b [ 0 ] * ( usum [ 0 ] + un.at(3) ) + c [ 0 ] * ( usum [ 1 ] + un.at(4) ) );
    __tmp.at(2, 6) = coeff * ( b [ 0 ] * ( usum [ 0 ] + un.at(5) ) + c [ 0 ] * ( usum [ 1 ] + un.at(6) ) );

    __tmp.at(3, 1) = coeff * ( b [ 1 ] * ( usum [ 0 ] + un.at(1) ) + c [ 1 ] * ( usum [ 1 ] + un.at(2) ) );
    __tmp.at(3, 3) = coeff * ( b [ 1 ] * ( usum [ 0 ] + un.at(3) ) + c [ 1 ] * ( usum [ 1 ] + un.at(4) ) );
    __tmp.at(3, 5) = coeff * ( b [ 1 ] * ( usum [ 0 ] + un.at(5) ) + c [ 1 ] * ( usum [ 1 ] + un.at(6) ) );

    __tmp.at(4, 2) = coeff * ( b [ 1 ] * ( usum [ 0 ] + un.at(1) ) + c [ 1 ] * ( usum [ 1 ] + un.at(2) ) );
    __tmp.at(4, 4) = coeff * ( b [ 1 ] * ( usum [ 0 ] + un.at(3) ) + c [ 1 ] * ( usum [ 1 ] + un.at(4) ) );
    __tmp.at(4, 6) = coeff * ( b [ 1 ] * ( usum [ 0 ] + un.at(5) ) + c [ 1 ] * ( usum [ 1 ] + un.at(6) ) );

    __tmp.at(5, 1) = coeff * ( b [ 2 ] * ( usum [ 0 ] + un.at(1) ) + c [ 2 ] * ( usum [ 1 ] + un.at(2) ) );
    __tmp.at(5, 3) = coeff * ( b [ 2 ] * ( usum [ 0 ] + un.at(3) ) + c [ 2 ] * ( usum [ 1 ] + un.at(4) ) );
    __tmp.at(5, 5) = coeff * ( b [ 2 ] * ( usum [ 0 ] + un.at(5) ) + c [ 2 ] * ( usum [ 1 ] + un.at(6) ) );

    __tmp.at(6, 2) = coeff * ( b [ 2 ] * ( usum [ 0 ] + un.at(1) ) + c [ 2 ] * ( usum [ 1 ] + un.at(2) ) );
    __tmp.at(6, 4) = coeff * ( b [ 2 ] * ( usum [ 0 ] + un.at(3) ) + c [ 2 ] * ( usum [ 1 ] + un.at(4) ) );
    __tmp.at(6, 6) = coeff * ( b [ 2 ] * ( usum [ 0 ] + un.at(5) ) + c [ 2 ] * ( usum [ 1 ] + un.at(6) ) );

    __ctilda_norm = __tmp.computeFrobeniusNorm();

    // compute e mtrx (advection term of momentum balance)
    __tmp.resize(6, 6);
    __tmp.zero();
    double __n[] = {
        b [ 0 ], c [ 0 ], b [ 1 ], c [ 1 ], b [ 2 ], c [ 2 ]
    };

    for ( i = 1; i <= 6; i++ ) {
        for ( j = 1; j <= 6; j++ ) {
            __tmp.at(i, j) = coeff * __n [ i - 1 ] * __n [ j - 1 ];
        }
    }

    __e_norm = __tmp.computeFrobeniusNorm();
    // compute element level Reynolds number
    // compute k~ norm first
    __tmp.resize(6, 6);
    __tmp.zero();
    FloatMatrix _b(3, 6), _d, _db;
    double ar12 = area / 12.;
    for ( i = 1; i <= 2; i++ ) { // test function index
        for ( k = 1; k <= 3; k++ ) { // nodal val of function w
            for ( j = 1; j <= 2; j++ ) { // velocity vector component
                for ( m = 1; m <= 3; m++ ) { //  nodal components
                    w_dof_addr = ( k - 1 ) * 2 + i;
                    u_dof_addr = ( m - 1 ) * 2 + j;
                    d1j = ( j == 1 );
                    d2j = ( j == 2 );
                    dij = ( i == j );
                    __tmp.at(w_dof_addr, u_dof_addr) += rho *
                                                        (
                        0.0 * d1j * b [ k - 1 ] * dudx [ i - 1 ] [ 0 ] * ar12 * ( usum [ 0 ] + un.at( ( m - 1 ) * 2 + 1 ) ) +
                        0.0 * d1j * b [ k - 1 ] * dudx [ i - 1 ] [ 1 ] * ar12 * ( usum [ 1 ] + un.at( ( m - 1 ) * 2 + 2 ) ) +
                        0.0 * d2j * c [ k - 1 ] * dudx [ i - 1 ] [ 0 ] * ar12 * ( usum [ 0 ] + un.at( ( m - 1 ) * 2 + 1 ) ) +
                        0.0 * d2j * c [ k - 1 ] * dudx [ i - 1 ] [ 1 ] * ar12 * ( usum [ 1 ] + un.at( ( m - 1 ) * 2 + 2 ) ) +

                        0.0 * d1j * b [ k - 1 ] * dudx [ i - 1 ] [ 0 ] * ar12 * ( usum [ 0 ] + un.at( ( m - 1 ) * 2 + 1 ) ) +
                        dij * b [ k - 1 ] * b [ m - 1 ] * ar12 * ( usum [ 0 ] * usum [ 0 ] + un.at(1) * un.at(1) + un.at(3) * un.at(3) + un.at(5) * un.at(5) ) +
                        0.0 * d2j * b [ k - 1 ] * dudx [ i - 1 ] [ 1 ] * ar12 * ( usum [ 0 ] + un.at( ( m - 1 ) * 2 + 1 ) ) +
                        dij * b [ k - 1 ] * c [ m - 1 ] * ar12 * ( usum [ 0 ] * usum [ 1 ] + un.at(1) * un.at(2) + un.at(3) * un.at(4) + un.at(5) * un.at(6) ) +

                        0.0 * d1j * c [ k - 1 ] * dudx [ i - 1 ] [ 0 ] * ar12 * ( usum [ 1 ] + un.at( ( m - 1 ) * 2 + 2 ) ) +
                        dij * c [ k - 1 ] * b [ m - 1 ] * ar12 * ( usum [ 0 ] * usum [ 1 ] + un.at(1) * un.at(2) + un.at(3) * un.at(4) + un.at(5) * un.at(6) ) +
                        0.0 * d2j * c [ k - 1 ] * dudx [ i - 1 ] [ 1 ] * ar12 * ( usum [ 1 ] + un.at( ( m - 1 ) * 2 + 2 ) ) +
                        dij * c [ k - 1 ] * c [ m - 1 ] * ar12 * ( usum [ 1 ] * usum [ 1 ] + un.at(2) * un.at(2) + un.at(4) * un.at(4) + un.at(6) * un.at(6) )
                                                        );
                }
            }
        }
    }

    __k_norm = __tmp.computeFrobeniusNorm();
    double u_1, u_2, vnorm = 0.;
    int im1;
    for ( i = 1; i <= 3; i++ ) {
        im1 = i - 1;
        u_1 = u.at( ( im1 ) * 2 + 1 );
        u_2 = u.at( ( im1 ) * 2 + 2 );
        vnorm = max( vnorm, sqrt(u_1 * u_1 + u_2 * u_2) );
    }

    if ( vnorm == 0.0 ) {
        //t_sugn1 = inf;
        double t_sugn2 = atTime->giveTimeIncrement() / 2.0;
        //t_sugn3 = inf;
        this->t_supg = 1. / sqrt( 1. / ( t_sugn2 * t_sugn2 ) );
        this->t_pspg = this->t_supg;
        this->t_lsic = 0.0;
        //printf ("%e %e\n",this->t_supg,this->t_pspg);
    } else {
        __Re = vnorm * vnorm * __c_norm / __k_norm / nu;

        __t_p1 = __g_norm / __gamma_norm;
        __t_p2 = atTime->giveTimeIncrement() * __g_norm / 2.0 / __beta_norm;
        __t_p3 = __t_p1 * __Re;
        this->t_pspg = 1. / sqrt( 1. / ( __t_p1 * __t_p1 ) + 1. / ( __t_p2 * __t_p2 ) + 1. / ( __t_p3 * __t_p3 ) );


        __t_s1 = __c_norm / __k_norm;
        __t_s2 = atTime->giveTimeIncrement() * __c_norm / __ctilda_norm / 2.0;
        __t_s3 = __t_s1 * __Re;
        this->t_supg = 1. / sqrt( 1. / ( __t_s1 * __t_s1 ) + 1. / ( __t_s2 * __t_s2 ) + 1. / ( __t_s3 * __t_s3 ) );

        //printf ("%e %e\n",this->t_supg,this->t_pspg);
        this->t_lsic = __c_norm / __e_norm;
        this->t_lsic = 0.0;
    }

#else
    /* UGN-Based Stabilization */
    double h_ugn, sum = 0.0, vnorm, t_sugn1, t_sugn2, t_sugn3, u_1, u_2, z, Re_ugn;
    double dscale, uscale, lscale, tscale, dt;
    //bool zeroFlag = false;
    int i, im1;
    FloatArray u;

    uscale = domain->giveEngngModel()->giveVariableScale(VST_Velocity);
    lscale = domain->giveEngngModel()->giveVariableScale(VST_Length);
    tscale = domain->giveEngngModel()->giveVariableScale(VST_Time);
    dscale = domain->giveEngngModel()->giveVariableScale(VST_Density);

    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime->givePreviousStep(), u);
    u.times(uscale);
    double nu;

    // compute averaged viscosity based on rule of mixture
    GaussPoint *gp;
    if ( integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints() ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    } else {
        gp = integrationRulesArray [ 1 ]->getIntegrationPoint(0);
    }

    nu = this->giveMaterial()->giveCharacteristicValue( MRM_Viscosity, gp, atTime->givePreviousStep() );
    nu *= domain->giveEngngModel()->giveVariableScale(VST_Viscosity);

    dt = atTime->giveTimeIncrement() * tscale;

    for ( i = 1; i <= 3; i++ ) {
        im1 = i - 1;
        sum += fabs(u.at( ( im1 ) * 2 + 1 ) * b [ im1 ] / lscale + u.at(im1 * 2 + 2) * c [ im1 ] / lscale);
    }

    /*
     * u_1=(u.at(1)+u.at(3)+u.at(5))/3.0;
     * u_2=(u.at(2)+u.at(4)+u.at(6))/3.0;
     * vnorm=sqrt(u_1*u_1+u_2*u_2);
     */
    vnorm = 0.;
    for ( i = 1; i <= 3; i++ ) {
        im1 = i - 1;
        u_1 = u.at( ( im1 ) * 2 + 1 );
        u_2 = u.at( ( im1 ) * 2 + 2 );
        vnorm = max( vnorm, sqrt(u_1 * u_1 + u_2 * u_2) );
    }

    if ( ( vnorm == 0.0 ) || ( sum <  vnorm * 1e-10 ) ) {
        //t_sugn1 = inf;
        t_sugn2 = dt / 2.0;
        //t_sugn3 = inf;
        this->t_supg = 1. / sqrt( 1. / ( t_sugn2 * t_sugn2 ) );
        this->t_pspg = this->t_supg;
        this->t_lsic = 0.0;
    } else {
        h_ugn = 2.0 * vnorm / sum;
        t_sugn1 = 1. / sum;
        t_sugn2 = dt / 2.0;
        t_sugn3 = h_ugn * h_ugn / 4.0 / nu;

        this->t_supg = 1. / sqrt( 1. / ( t_sugn1 * t_sugn1 ) + 1. / ( t_sugn2 * t_sugn2 ) + 1. / ( t_sugn3 * t_sugn3 ) );
        this->t_pspg = this->t_supg;

        Re_ugn = vnorm * h_ugn / ( 2. * nu );
        z = ( Re_ugn <= 3. ) ? Re_ugn / 3. : 1.0;
        this->t_lsic = h_ugn * vnorm * z / 2.0;
    }

    // if (this->number == 1) {
    //  printf ("t_supg %e t_pspg %e t_lsic %e\n", t_supg, t_pspg, t_lsic);
    // }


    this->t_supg *= uscale / lscale;
    this->t_pspg *= 1. / ( lscale * dscale );
    this->t_lsic *= ( dscale * uscale ) / ( lscale * lscale );

    this->t_lsic = 0.0;

#endif

    //this->t_lsic=0.0;
    //this->t_pspg=0.0;
}

/*
 * void
 * TR1_2D_SUPG :: updateStabilizationCoeffs (TimeStep* atTime)
 * {
 * double h_ugn, sum=0.0, snorm, vnorm, t_sugn1, t_sugn2, t_sugn3, u_1, u_2, z, Re_ugn;
 * bool zeroFlag = false;
 * int i,im1;
 * FloatArray u;
 * this -> computeVectorOf(EID_MomentumBalance,VM_Total,atTime, u) ;
 * double nu = this->giveMaterial()->giveCharacteristicValue(MRM_Viscosity, integrationRulesArray[0]->getIntegrationPoint(0), atTime);
 *
 * // UGN-Based Stabilization
 *
 * for (i=1;i<=3;i++) {
 *  im1=i-1;
 *  snorm = sqrt(u.at(im1*2+1)*u.at(im1*2+1) + u.at(im1*2+2)*u.at(im1*2+2));
 *  if (snorm == 0.0) zeroFlag = true;
 *  sum+= fabs(u.at((im1)*2+1)*b[im1] + u.at(im1*2+2)*c[im1])/snorm;
 * }
 * u_1=(u.at(1)+u.at(3)+u.at(5))/3.0;
 * u_2=(u.at(2)+u.at(4)+u.at(6))/3.0;
 * vnorm=sqrt(u_1*u_1+u_2*u_2);
 *
 * if (zeroFlag) {
 *
 *  //t_sugn1 = inf;
 *  t_sugn2 = atTime->giveTimeIncrement()/2.0;
 *  //t_sugn3 = inf;
 *  this->t_supg=1./sqrt(1./(t_sugn2*t_sugn2));
 *  this->t_pspg=this->t_supg;
 *  this->t_lsic=0.0;
 *
 * } else {
 *  h_ugn = 2.0*vnorm/sum;
 *  //t_sugn1 = h_ugn/(2.*vnorm);
 *  t_sugn1 = 1./sum;
 *  t_sugn2 = atTime->giveTimeIncrement()/2.0;
 *  t_sugn3 = h_ugn*h_ugn/4.0/nu;
 *
 *  this->t_supg=1./sqrt(1./(t_sugn1*t_sugn1) + 1./(t_sugn2*t_sugn2) + 1./(t_sugn3*t_sugn3));
 *  this->t_pspg=this->t_supg;
 *
 *  Re_ugn = vnorm*h_ugn/(2.*nu);
 *  z = (Re_ugn <= 3.)?Re_ugn/3. : 1.0 ;
 *  this->t_lsic=h_ugn*vnorm*z/2.0;
 *  this->t_lsic = 0.0;
 * }
 * //printf ("Element %d ( t_supg %e, t_pspg %e, t_lsic %e)\n", this->giveNumber(), t_supg, t_pspg, t_lsic);
 * }
 *
 */

double
TR1_2D_SUPG :: computeCriticalTimeStep(TimeStep *tStep)
{
    return 1.e6;
}


Interface *
TR1_2D_SUPG :: giveInterface(InterfaceType interface)
{
    if ( interface == ZZNodalRecoveryModelInterfaceType ) {
        return ( ZZNodalRecoveryModelInterface * ) this;
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType ) {
        return ( NodalAveragingRecoveryModelInterface * ) this;
    } else if ( interface == SPRNodalRecoveryModelInterfaceType ) {
        return ( SPRNodalRecoveryModelInterface * ) this;
    } else if ( interface == SpatialLocalizerInterfaceType ) {
        return ( SpatialLocalizerInterface * ) this;
    } else if ( interface == EIPrimaryFieldInterfaceType ) {
        return ( EIPrimaryFieldInterface * ) this;
    } else if ( interface == LEPlicElementInterfaceType ) {
        return ( LEPlicElementInterface * ) this;
    } else if ( interface == LevelSetPCSElementInterfaceType ) {
        return ( LevelSetPCSElementInterface * ) this;
    }

    return NULL;
}


int
TR1_2D_SUPG :: SpatialLocalizerI_containsPoint(const FloatArray &coords) {
    FloatArray lcoords;
    return this->computeLocalCoordinates(lcoords, coords);
}

double
TR1_2D_SUPG :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
{
    FloatArray lcoords(3), gcoords;
    double dist;
    int size, gsize;

    lcoords.at(1) = lcoords.at(2) = lcoords.at(3) = 1. / 3.;
    this->computeGlobalCoordinates(gcoords, lcoords);

    if ( ( size = coords.giveSize() ) < ( gsize = gcoords.giveSize() ) ) {
        _error("SpatialLocalizerI_giveDistanceFromParametricCenter: coordinates size mismatch");
    }

    if ( size == gsize ) {
        dist = coords.distance(gcoords);
    } else {
        FloatArray helpCoords = coords;

        helpCoords.resize(gsize);
        dist = helpCoords.distance(gcoords);
    }

    return dist;
}


void
TR1_2D_SUPG :: computeDeviatoricStrain(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    /* one should call material driver instead */
    FloatArray u(6);
    answer.resize(3);


    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);

    answer.at(1) = ( b [ 0 ] * u.at(1) + b [ 1 ] * u.at(3) + b [ 2 ] * u.at(5) );
    answer.at(2) = ( c [ 0 ] * u.at(2) + c [ 1 ] * u.at(4) + c [ 2 ] * u.at(6) );
    answer.at(3) = ( b [ 0 ] * u.at(2) + b [ 1 ] * u.at(4) + b [ 2 ] * u.at(6) + c [ 0 ] * u.at(1) + c [ 1 ] * u.at(3) + c [ 2 ] * u.at(5) );
}

void
TR1_2D_SUPG :: initGeometry()
{
    Node *node1, *node2, *node3;
    double x1, x2, x3, y1, y2, y3;

    node1 = giveNode(1);
    node2 = giveNode(2);
    node3 = giveNode(3);

    // init geometry data
    x1 = node1->giveCoordinate(1);
    x2 = node2->giveCoordinate(1);
    x3 = node3->giveCoordinate(1);

    y1 = node1->giveCoordinate(2);
    y2 = node2->giveCoordinate(2);
    y3 = node3->giveCoordinate(2);

    this->area = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 );

    if ( area < 0.0 ) {
        _error("Area is negative, check element numbering orientation");
    }

    b [ 0 ] = ( y2 - y3 ) / ( 2. * area );
    c [ 0 ] = ( x3 - x2 ) / ( 2. * area );
    b [ 1 ] = ( y3 - y1 ) / ( 2. * area );
    c [ 1 ] = ( x1 - x3 ) / ( 2. * area );
    b [ 2 ] = ( y1 - y2 ) / ( 2. * area );
    c [ 2 ] = ( x2 - x1 ) / ( 2. * area );
}


int
TR1_2D_SUPG :: checkConsistency()
{
    return SUPGElement :: checkConsistency();
}

void
TR1_2D_SUPG :: computeNMtrx(FloatArray &answer, GaussPoint *gp)
{
    double l1, l2;
    answer.resize(3);

    answer.at(1) = l1 = gp->giveCoordinate(1);
    answer.at(2) = l2 = gp->giveCoordinate(2);
    answer.at(3) = 1.0 - l1 - l2;
}

/*
 * double
 * TR1_2D_SUPG :: computeCriticalTimeStep (TimeStep* tStep)
 * {
 * FloatArray u;
 * double dt1, dt2, dt;
 * double Re = domain->giveEngngModel()->giveUnknownComponent(ReynoldsNumber, VM_Unknown, tStep, domain, NULL);
 *
 * this -> computeVectorOf(EID_MomentumBalance,VM_Total,tStep, u) ;
 *
 * double vn1 = sqrt(u.at(1)*u.at(1)+u.at(2)*u.at(2));
 * double vn2 = sqrt(u.at(3)*u.at(3)+u.at(4)*u.at(4));
 * double vn3 = sqrt(u.at(5)*u.at(5)+u.at(6)*u.at(6));
 * double veln = max (vn1, max(vn2,vn3));
 *
 * double l1 = 1.0/(sqrt(b[0]*b[0]+c[0]*c[0]));
 * double l2 = 1.0/(sqrt(b[1]*b[1]+c[1]*c[1]));
 * double l3 = 1.0/(sqrt(b[2]*b[2]+c[2]*c[2]));
 *
 * double ln = min (l1, min (l2,l3));
 *
 * // viscous limit
 * dt2 = 0.5*ln*ln*Re;
 * if (veln != 0.0) {
 *  dt1 = ln/veln;
 *  dt = dt1*dt2/(dt1+dt2);
 * } else {
 *  dt = dt2;
 * }
 * return dt;
 * }
 */


/* Note this LEPLIC implementation will not work for axi elements as the true volume
   is not taken into account (the effect of radius) !
*/


double
TR1_2D_SUPG :: computeLEPLICVolumeFraction(const FloatArray &n, const double p, LEPlic *matInterface, bool updFlag)
{
    Polygon pg;
    double answer, volume = computeMyVolume(matInterface, updFlag);
    this->formVolumeInterfacePoly(pg, matInterface, n, p, updFlag);
    answer = fabs(pg.computeVolume() / volume);
    if ( answer > 1.0 ) {
        return 1.0;
    } else {
        return answer;
    }
}

void
TR1_2D_SUPG :: formMaterialVolumePoly(Polygon &matvolpoly, LEPlic *matInterface,
                                      const FloatArray &normal, const double p, bool updFlag)
{
    int i;
    double x, y;
    Vertex v;

    matvolpoly.clear();

    if ( this->vof <= TRSUPG_ZERO_VOF ) {
        return;
    } else if ( this->vof >= ( 1 - TRSUPG_ZERO_VOF ) ) {
        for ( i = 1; i <= 3; i++ ) {
            if ( updFlag ) {
                x = matInterface->giveUpdatedXCoordinate( this->giveNode(i)->giveNumber() );
                y = matInterface->giveUpdatedYCoordinate( this->giveNode(i)->giveNumber() );
            } else {
                x = this->giveNode(i)->giveCoordinate(1);
                y = this->giveNode(i)->giveCoordinate(2);
            }

            v.setCoords(x, y);
            matvolpoly.addVertex(v);
        }

        return;
    }

    this->formVolumeInterfacePoly(matvolpoly, matInterface, normal, p, updFlag);
}


void
TR1_2D_SUPG :: formVolumeInterfacePoly(Polygon &matvolpoly, LEPlic *matInterface,
                                       const FloatArray &normal, const double p, bool updFlag)
{
    int i, next;
    bool nodeIn [ 3 ];
    double nx = normal.at(1), ny = normal.at(2), x, y;
    double tx, ty;
    Vertex v;

    matvolpoly.clear();

    for ( i = 1; i <= 3; i++ ) {
        if ( updFlag ) {
            x = matInterface->giveUpdatedXCoordinate( this->giveNode(i)->giveNumber() );
            y = matInterface->giveUpdatedYCoordinate( this->giveNode(i)->giveNumber() );
        } else {
            x = this->giveNode(i)->giveCoordinate(1);
            y = this->giveNode(i)->giveCoordinate(2);
        }

        if ( ( nx * x + ny * y + p ) >= 0. ) {
            nodeIn [ i - 1 ] = true;
        } else {
            nodeIn [ i - 1 ] = false;
        }
    }

    if ( nodeIn [ 0 ] && nodeIn [ 1 ] && nodeIn [ 2 ] ) { // all nodes inside
        for ( i = 1; i <= 3; i++ ) {
            if ( updFlag ) {
                x = matInterface->giveUpdatedXCoordinate( this->giveNode(i)->giveNumber() );
                y = matInterface->giveUpdatedYCoordinate( this->giveNode(i)->giveNumber() );
            } else {
                x = this->giveNode(i)->giveCoordinate(1);
                y = this->giveNode(i)->giveCoordinate(2);
            }

            v.setCoords(x, y);
            matvolpoly.addVertex(v);
        }

        return;
    } else if ( !( nodeIn [ 0 ] || nodeIn [ 1 ] || nodeIn [ 2 ] ) ) { // all nodes outside
        return;
    } else {
        for ( i = 1; i <= 3; i++ ) {
            next = i < 3 ? i + 1 : 1;
            if ( nodeIn [ i - 1 ] ) {
                if ( updFlag ) {
                    v.setCoords( matInterface->giveUpdatedXCoordinate( this->giveNode(i)->giveNumber() ),
                                matInterface->giveUpdatedYCoordinate( this->giveNode(i)->giveNumber() ) );
                } else {
                    v.setCoords( this->giveNode(i)->giveCoordinate(1),
                                this->giveNode(i)->giveCoordinate(2) );
                }

                matvolpoly.addVertex(v);
            }

            if ( nodeIn [ next - 1 ] ^ nodeIn [ i - 1 ] ) {
                // compute intersection with (i,next) edge
                if ( updFlag ) {
                    x = matInterface->giveUpdatedXCoordinate( this->giveNode(i)->giveNumber() );
                    y = matInterface->giveUpdatedYCoordinate( this->giveNode(i)->giveNumber() );
                    tx = matInterface->giveUpdatedXCoordinate( this->giveNode(next)->giveNumber() ) - x;
                    ty = matInterface->giveUpdatedYCoordinate( this->giveNode(next)->giveNumber() ) - y;
                } else {
                    x = this->giveNode(i)->giveCoordinate(1);
                    y = this->giveNode(i)->giveCoordinate(2);
                    tx = this->giveNode(next)->giveCoordinate(1) - x;
                    ty = this->giveNode(next)->giveCoordinate(2) - y;
                }

                double s, sd = nx * tx + ny * ty;
                if ( fabs(sd) > 1.e-6 ) {
                    s = ( -p - ( nx * x + ny * y ) ) / sd;
                    v.setCoords(x + tx * s, y + ty * s);
                    matvolpoly.addVertex(v);
                } else {
                    // pathological case - lines are parallel
                    if ( nodeIn [ i - 1 ] ) {
                        if ( updFlag ) {
                            v.setCoords( matInterface->giveUpdatedXCoordinate( this->giveNode(next)->giveNumber() ),
                                        matInterface->giveUpdatedYCoordinate( this->giveNode(next)->giveNumber() ) );
                        } else {
                            v.setCoords( this->giveNode(next)->giveCoordinate(1), this->giveNode(next)->giveCoordinate(2) );
                        }

                        matvolpoly.addVertex(v);
                    } else {
                        v.setCoords(x, y);
                        matvolpoly.addVertex(v);
                        if ( updFlag ) {
                            v.setCoords( matInterface->giveUpdatedXCoordinate( this->giveNode(next)->giveNumber() ),
                                        matInterface->giveUpdatedYCoordinate( this->giveNode(next)->giveNumber() ) );
                        } else {
                            v.setCoords( this->giveNode(next)->giveCoordinate(1), this->giveNode(next)->giveCoordinate(2) );
                        }

                        matvolpoly.addVertex(v);
                    }
                }
            }
        } // end loop over elem nodes

    }
}


double
TR1_2D_SUPG :: truncateMatVolume(const Polygon &matvolpoly, double &volume)
{
    Polygon me, clip;
    Graph g;

    this->formMyVolumePoly(me, NULL, false);
    g.clip(clip, me, matvolpoly);
#ifdef __OOFEG
    EASValsSetColor( gc [ 0 ].getActiveCrackColor() );
    //GraphicObj *go = clip.draw(::gc[OOFEG_DEBUG_LAYER],true);
    clip.draw(gc [ OOFEG_DEBUG_LAYER ], true);
    //EVFastRedraw(myview);
#endif
    volume = clip.computeVolume();
    return volume / area;
}

void
TR1_2D_SUPG :: formMyVolumePoly(Polygon &me, LEPlic *matInterface, bool updFlag)
{
    double x, y;
    Vertex v;

    me.clear();

    for ( int i = 1; i <= 3; i++ ) {
        if ( updFlag ) {
            x = matInterface->giveUpdatedXCoordinate( this->giveNode(i)->giveNumber() );
            y = matInterface->giveUpdatedYCoordinate( this->giveNode(i)->giveNumber() );
        } else {
            x = this->giveNode(i)->giveCoordinate(1);
            y = this->giveNode(i)->giveCoordinate(2);
        }

        v.setCoords(x, y);
        me.addVertex(v);
    }
}


double
TR1_2D_SUPG :: computeMyVolume(LEPlic *matInterface, bool updFlag)
{
    double x1, x2, x3, y1, y2, y3;
    if ( updFlag ) {
        x1 = matInterface->giveUpdatedXCoordinate( this->giveNode(1)->giveNumber() );
        x2 = matInterface->giveUpdatedXCoordinate( this->giveNode(2)->giveNumber() );
        x3 = matInterface->giveUpdatedXCoordinate( this->giveNode(3)->giveNumber() );

        y1 = matInterface->giveUpdatedYCoordinate( this->giveNode(1)->giveNumber() );
        y2 = matInterface->giveUpdatedYCoordinate( this->giveNode(2)->giveNumber() );
        y3 = matInterface->giveUpdatedYCoordinate( this->giveNode(3)->giveNumber() );
        return 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 );
    } else {
        return area;
    }
}

double
TR1_2D_SUPG :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the portion of the receiver which is attached to aGaussPoint.
{
    double determinant, weight, volume;

    determinant = fabs( 4 * area * area * ( c [ 1 ] * b [ 0 ] - c [ 0 ] * b [ 1 ] ) );


    weight      = aGaussPoint->giveWeight();
    volume      = determinant * weight;

    return volume;
}


double
TR1_2D_SUPG :: computeCriticalLEPlicTimeStep(TimeStep *tStep)
{
    FloatArray u;
    double Re = domain->giveEngngModel()->giveUnknownComponent(ReynoldsNumber, VM_Unknown, tStep, domain, NULL);

    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);

    double vn1 = sqrt( u.at(1) * u.at(1) + u.at(2) * u.at(2) );
    double vn2 = sqrt( u.at(3) * u.at(3) + u.at(4) * u.at(4) );
    double vn3 = sqrt( u.at(5) * u.at(5) + u.at(6) * u.at(6) );
    double veln = max( vn1, max(vn2, vn3) );

    double l1 = 1.0 / ( sqrt(b [ 0 ] * b [ 0 ] + c [ 0 ] * c [ 0 ]) );
    double l2 = 1.0 / ( sqrt(b [ 1 ] * b [ 1 ] + c [ 1 ] * c [ 1 ]) );
    double l3 = 1.0 / ( sqrt(b [ 2 ] * b [ 2 ] + c [ 2 ] * c [ 2 ]) );

    double ln = min( l1, min(l2, l3) );

    if ( veln != 0.0 ) {
        return ln / veln;
    } else {
        return 0.5 * ln * ln * Re;
    }
}


void
TR1_2D_SUPG :: giveElementCenter(LEPlic *mat_interface, FloatArray &center, bool upd)
{
    FloatArray coords;

    center.resize(2);
    center.zero();
    if ( upd ) {
        for ( int i = 1; i <= 3; i++ ) {
            mat_interface->giveUpdatedCoordinate( coords, giveNode(i)->giveNumber() );
            center.add(coords);
        }
    } else {
        for ( int i = 1; i <= 3; i++ ) {
            center.at(1) += this->giveNode(i)->giveCoordinate(1);
            center.at(2) += this->giveNode(i)->giveCoordinate(2);
        }
    }

    center.times(1. / 3.);
}

int
TR1_2D_SUPG :: EIPrimaryFieldI_evaluateFieldVectorAt(FloatArray &answer, PrimaryField &pf,
                                                     FloatArray &coords, IntArray &dofId, ValueModeType mode,
                                                     TimeStep *atTime)
{
    int i, j, indx, es;
    double sum;
    FloatArray elemvector, f, lc;
    //FloatMatrix n;
    IntArray elemdofs;
    // determine element dof ids
    this->giveElementDofIDMask(pf.giveEquationID(), elemdofs);
    es = elemdofs.giveSize();
    // first evaluate element unknown vector
    this->computeVectorOf(pf, mode, atTime, elemvector);
    // determine corresponding local coordinates
    if ( this->computeLocalCoordinates(lc, coords) ) {
        // compute interpolation matrix
        // this->computeNmatrixAt(n, &lc);
        // compute answer
        answer.resize( dofId.giveSize() );
        for ( i = 1; i <= dofId.giveSize(); i++ ) {
            if ( ( indx = elemdofs.findFirstIndexOf( dofId.at(i) ) ) ) {
                for ( j = 1, sum = 0.0; j <= 3; j++ ) {
                    sum += lc.at(j) * elemvector.at(es * ( j - 1 ) + indx);
                }

                answer.at(i) = sum;
            } else {
                // _error("EIPrimaryFieldI_evaluateFieldVectorAt: unknown dof id encountered");
                answer.at(i) = 0.0;
            }
        }

        return 0; // ok
    } else {
        _error("EIPrimaryFieldI_evaluateFieldVectorAt: target point not in receiver volume");
        return 1; // fail
    }
}


void
TR1_2D_SUPG :: updateYourself(TimeStep *tStep)
{
    SUPGElement :: updateYourself(tStep);
    LEPlicElementInterface :: updateYourself(tStep);
}

int
TR1_2D_SUPG :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    if ( type == IST_VOFFraction ) {
        MaterialInterface *mi = domain->giveEngngModel()->giveMaterialInterface( domain->giveNumber() );
        if ( mi ) {
            FloatArray val;
            mi->giveElementMaterialMixture( val, aGaussPoint->giveElement()->giveNumber() );
            answer.resize(1);
            answer.at(1) = val.at(1);
            return 1;
        } else {
            answer.resize(1);
            answer.at(1) = 1.0;
            return 1;
        }
    } else if ( type == IST_Density ) {
        answer.resize(1);
        answer.at(1) = this->giveMaterial()->giveCharacteristicValue(MRM_Density, aGaussPoint, atTime);
        return 1;
    } else {
        return SUPGElement :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}

int
TR1_2D_SUPG :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type)
{
    if ( ( type == IST_VOFFraction ) || ( type == IST_Density ) ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else {
        return SUPGElement :: giveIntVarCompFullIndx(answer, type);
    }
}


InternalStateValueType
TR1_2D_SUPG :: giveIPValueType(InternalStateType type)
{
    if ( ( type == IST_VOFFraction ) || ( type == IST_Density ) ) {
        return ISVT_SCALAR;
    } else {
        return SUPGElement :: giveIPValueType(type);
    }
}


int
TR1_2D_SUPG :: giveIPValueSize(InternalStateType type, GaussPoint *gp)
{
    if ( ( type == IST_VOFFraction ) || ( type == IST_Density ) ) {
        return 1;
    } else {
        return SUPGElement :: giveIPValueSize(type, gp);
    }
}


int
TR1_2D_SUPG :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( ( type == IST_StressTensor ) || ( type == IST_StrainTensor ) ) {
        return 4;
    }

    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    return this->giveIPValueSize(type, gp);
}


void
TR1_2D_SUPG :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                          InternalStateType type, TimeStep *tStep)
{
    GaussPoint *gp;
    gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    this->giveIPValue(answer, gp, type, tStep);
}

void
TR1_2D_SUPG :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                         InternalStateType type, TimeStep *tStep)
{
    answer.resize(0);
}

void
TR1_2D_SUPG :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(3);
    pap.at(1) = this->giveNode(1)->giveNumber();
    pap.at(2) = this->giveNode(2)->giveNumber();
    pap.at(3) = this->giveNode(3)->giveNumber();
}

void
TR1_2D_SUPG :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
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

int
TR1_2D_SUPG :: SPRNodalRecoveryMI_giveNumberOfIP()
{ return 1; }


void
TR1_2D_SUPG :: SPRNodalRecoveryMI_computeIPGlobalCoordinates(FloatArray &coords, GaussPoint *gp)
{
    if ( gp == integrationRulesArray [ 0 ]->getIntegrationPoint(0) ) {
        this->computeGlobalCoordinates( coords, * gp->giveCoordinates() );
    } else {
        _error("SPRNodalRecoveryMI_computeIPGlobalCoordinates: unsupported ip num");
    }
}

SPRPatchType
TR1_2D_SUPG :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_2dxy;
}



void
TR1_2D_SUPG :: printOutputAt(FILE *file, TimeStep *stepN)
// Performs end-of-step operations.
{
    SUPGElement :: printOutputAt(file, stepN);
    double rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), stepN);
    fprintf(file, "VOF %e, density %e\n\n", this->giveVolumeFraction(), rho);
}



contextIOResultType TR1_2D_SUPG :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full element context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;

    if ( ( iores = SUPGElement :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = LEPlicElementInterface :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}



contextIOResultType TR1_2D_SUPG :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full element context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;

    if ( ( iores = SUPGElement :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( iores = LEPlicElementInterface :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }


    return CIO_OK;
}


double
TR1_2D_SUPG :: LS_PCS_computeF(LevelSetPCS *ls, TimeStep *atTime)
{
    int i;
    double answer;
    FloatArray fi(3), un(6);

    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, un);

    for ( i = 1; i <= 3; i++ ) {
        fi.at(i) = ls->giveLevelSetDofManValue( dofManArray.at(i) );
    }

    double fix = b [ 0 ] * fi.at(1) + b [ 1 ] * fi.at(2) + b [ 2 ] * fi.at(3);
    double fiy = c [ 0 ] * fi.at(1) + c [ 1 ] * fi.at(2) + c [ 2 ] * fi.at(3);
    double norm = sqrt(fix * fix + fiy * fiy);

    answer = ( 1. / 3. ) * ( fix * ( un.at(1) + un.at(3) + un.at(5) ) + fiy * ( un.at(2) + un.at(4) + un.at(6) ) ) / norm;
    return answer;
}


double
TR1_2D_SUPG :: LS_PCS_computeS(LevelSetPCS *ls, TimeStep *atTime)
{
#if 0
    int i;
    double fi, answer, eps = 0.0;
    FloatArray s(3);

    for ( i = 1; i <= 3; i++ ) {
        fi = ls->giveLevelSetDofManValue( dofManArray.at(i) );
        s.at(i) = fi / sqrt(fi * fi + eps * eps);
    }

    answer = ( 1. / 3. ) * ( s.at(1) + s.at(2) + s.at(3) );
    return answer;

#else

    int i, neg = 0, pos = 0, zero = 0, si = 0;
    double x1, x2, x3, y1, y2, y3;
    FloatArray fi(3);

    for ( i = 1; i <= 3; i++ ) {
        fi.at(i) = ls->giveLevelSetDofManValue( dofManArray.at(i) );
    }


    for ( i = 1; i <= 3; i++ ) {
        if ( fi.at(i) >= 0. ) {
            pos++;
        }

        if ( fi.at(i) < 0.0 ) {
            neg++;
        }

        if ( fi.at(i) == 0. ) {
            zero++;
        }
    }

    if ( zero == 3 ) {
        return 0.0;
    } else if ( neg == 0 ) { // all level set values positive
        return 1.0; //return area;
    } else if ( pos == 0 ) { // all level set values negative
        return -1.0; //return -area;
    } else {
        // zero level set inside
        // find the vertex vith level set sign different from other two
        for ( i = 1; i <= 3; i++ ) {
            if ( ( pos > neg ) && ( fi.at(i) < 0.0 ) ) {
                si = i;
                break;
            }

            if ( ( pos < neg ) && ( fi.at(i) >= 0.0 ) ) {
                si = i;
                break;
            }
        }

        if ( si ) {
            x1 = this->giveNode(si)->giveCoordinate(1);
            y1 = this->giveNode(si)->giveCoordinate(2);

            // compute intersections
            int prev_node = ( si > 1 ) ? si - 1 : 3;
            int next_node = ( si < 3 ) ? si + 1 : 1;

            //double l = this->giveNode(si)->giveCoordinates()->distance(this->giveNode(next_node)->giveCoordinates());
            double t = fi.at(si) / ( fi.at(si) - fi.at(next_node) );
            x2 = x1 + t * ( this->giveNode(next_node)->giveCoordinate(1) - x1 );
            y2 = y1 + t * ( this->giveNode(next_node)->giveCoordinate(2) - y1 );

            //l = this->giveNode(si)->giveCoordinates()->distance(this->giveNode(prev_node)->giveCoordinates());
            t = fi.at(si) / ( fi.at(si) - fi.at(prev_node) );
            x3 = x1 + t * ( this->giveNode(prev_node)->giveCoordinate(1) - x1 );
            y3 = y1 + t * ( this->giveNode(prev_node)->giveCoordinate(2) - y1 );

            // compute area
            double __area = fabs( 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 ) );

            if ( pos > neg ) {
                // negative area computed
                return ( ( area - __area ) - __area ) / area;
            } else {
                // postive area computed
                return ( __area - ( area - __area ) ) / area;
            }
        } else {
            OOFEM_ERROR("TR1_2D_SUPG::LS_PCS_computeVOFFractions: internal consistency error");
            return 0.0;
        }
    }

#endif
}


void
TR1_2D_SUPG :: LS_PCS_computedN(FloatMatrix &answer)
{
    int i;
    answer.resize(3, 2);

    for ( i = 1; i <= 3; i++ ) {
        answer.at(i, 1) = b [ i - 1 ];
        answer.at(i, 2) = c [ i - 1 ];
    }
}


void
TR1_2D_SUPG :: LS_PCS_computeVOFFractions(FloatArray &answer, FloatArray &fi)
{
    int i, neg = 0, pos = 0, zero = 0, si = 0;
    double x1, x2, x3, y1, y2, y3;

    answer.resize(2);
    for ( i = 1; i <= 3; i++ ) {
        if ( fi.at(i) >= 0. ) {
            pos++;
        }

        if ( fi.at(i) < 0.0 ) {
            neg++;
        }

        if ( fi.at(i) == 0. ) {
            zero++;
        }
    }

    if ( neg == 0 ) { // all level set values positive
        answer.at(1) = 1.0;
        answer.at(2) = 0.0;
    } else if ( pos == 0 ) { // all level set values negative
        answer.at(1) = 0.0;
        answer.at(2) = 1.0;
    } else if ( zero == 3 ) {
        // ???????
        answer.at(1) = 1.0;
        answer.at(2) = 0.0;
    } else {
        // zero level set inside
        // find the vertex with level set sign different from other two
        for ( i = 1; i <= 3; i++ ) {
            if ( ( pos > neg ) && ( fi.at(i) < 0.0 ) ) {
                si = i;
                break;
            }

            if ( ( pos < neg ) && ( fi.at(i) >= 0.0 ) ) {
                si = i;
                break;
            }
        }

        if ( si && ( ( pos + neg ) == 3 ) ) {
            x1 = this->giveNode(si)->giveCoordinate(1);
            y1 = this->giveNode(si)->giveCoordinate(2);

            // compute intersections
            int prev_node = ( si > 1 ) ? si - 1 : 3;
            int next_node = ( si < 3 ) ? si + 1 : 1;

            double l = this->giveNode(si)->giveCoordinates()->distance( this->giveNode(next_node)->giveCoordinates() );
            double t = fi.at(si) / ( fi.at(si) - fi.at(next_node) );
            x2 = x1 + t * ( this->giveNode(next_node)->giveCoordinate(1) - this->giveNode(si)->giveCoordinate(1) );
            y2 = y1 + t * ( this->giveNode(next_node)->giveCoordinate(2) - this->giveNode(si)->giveCoordinate(2) );

            l = this->giveNode(si)->giveCoordinates()->distance( this->giveNode(prev_node)->giveCoordinates() );
            t = fi.at(si) / ( fi.at(si) - fi.at(prev_node) );
            x3 = x1 + t * ( this->giveNode(prev_node)->giveCoordinate(1) - this->giveNode(si)->giveCoordinate(1) );
            y3 = y1 + t * ( this->giveNode(prev_node)->giveCoordinate(2) - this->giveNode(si)->giveCoordinate(2) );

            // compute area
            double __area = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 );
            if ( fabs(__area) / area > 1.00001 ) {
                OOFEM_ERROR("TR1_2D_SUPG::LS_PCS_computeVOFFractions: internal consistency error");
            }

            // prevent some roundoff errors
            if ( fabs(__area) > area ) {
                __area = sgn(__area) * area;
            }

            if ( pos > neg ) {
                // negative area computed
                answer.at(2) = fabs(__area) / area;
                answer.at(1) = 1.0 - answer.at(2);
            } else {
                // postive area computed
                answer.at(1) = fabs(__area) / area;
                answer.at(2) = 1.0 - answer.at(1);
            }
        } else {
            OOFEM_ERROR("TR1_2D_SUPG::LS_PCS_computeVOFFractions: internal consistency error");
        }
    }
}


void
TR1_2D_SUPG :: giveLocalVelocityDofMap(IntArray &map)
{
    map.resize(6);
    map.at(1) = 1;
    map.at(2) = 2;
    map.at(3) = 4;
    map.at(4) = 5;
    map.at(5) = 7;
    map.at(6) = 8;
}

void
TR1_2D_SUPG :: giveLocalPressureDofMap(IntArray &map)
{
    map.resize(3);
    map.at(1) = 3;
    map.at(2) = 6;
    map.at(3) = 9;
}


#ifdef __OOFEG
int
TR1_2D_SUPG :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                       int node, TimeStep *atTime)
{
    /*
     * if (type == IST_VOFFraction) {
     * answer.resize(1);
     * answer.at(1) = this->giveTempVolumeFraction();
     * return 1;
     * } else if (type == IST_Density) {
     * answer.resize(1);
     * answer.at(1) = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray[0]-> getIntegrationPoint(0), atTime);
     * return 1;
     *
     * } else
     */return SUPGElement :: giveInternalStateAtNode(answer, type, mode, node, atTime);
}

void
TR1_2D_SUPG :: drawRawGeometry(oofegGraphicContext &gc)
{
    WCRec p [ 3 ];
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = 0.;
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = 0.;
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
    p [ 2 ].z = 0.;

    go =  CreateTriangle3D(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}

void TR1_2D_SUPG :: drawScalar(oofegGraphicContext &context)
{
    int i, indx, result = 0;
    WCRec p [ 3 ];
    GraphicObj *tr;
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    FloatArray v1, v2, v3;
    double s [ 3 ];
    IntArray map;

    if ( !context.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);

    // if ((context.giveIntVarMode() == ISM_local) && (context.giveIntVarType() ==  IST_VOFFraction)) {
    if ( ( context.giveIntVarType() ==  IST_VOFFraction ) && ( context.giveIntVarMode() == ISM_local ) ) {
        Polygon matvolpoly;
        this->formMaterialVolumePoly(matvolpoly, NULL, temp_normal, temp_p, false);
        EASValsSetColor( context.getStandardSparseProfileColor() );
        //GraphicObj *go = matvolpoly.draw(context,true,OOFEG_VARPLOT_PATTERN_LAYER);
        matvolpoly.draw(context, true, OOFEG_VARPLOT_PATTERN_LAYER);
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
            p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
            p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
            p [ i ].z = 0.;
        }

        //EASValsSetColor(gc.getYieldPlotColor(ratio));
        context.updateFringeTableMinMax(s, 3);
        tr =  CreateTriangleWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ]);
        EGWithMaskChangeAttributes(LAYER_MASK, tr);
        EMAddGraphicsToModel(ESIModel(), tr);
    } else if ( ( context.getScalarAlgo() == SA_ZPROFILE ) || ( context.getScalarAlgo() == SA_COLORZPROFILE ) ) {
        double landScale = context.getLandScale();

        for ( i = 0; i < 3; i++ ) {
            p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
            p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
            p [ i ].z = s [ i ] * landScale;
        }

        if ( context.getScalarAlgo() == SA_ZPROFILE ) {
            EASValsSetColor( context.getDeformedElementColor() );
            EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
            EASValsSetFillStyle(FILL_SOLID);
            tr =  CreateTriangle3D(p);
            EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | FILL_MASK | LAYER_MASK, tr);
        } else {
            context.updateFringeTableMinMax(s, 3);
            EASValsSetFillStyle(FILL_SOLID);
            tr =  CreateTriangleWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ]);
            EGWithMaskChangeAttributes(FILL_MASK | LAYER_MASK, tr);
        }

        EMAddGraphicsToModel(ESIModel(), tr);
    }
}



#endif
} // end namespace oofem
