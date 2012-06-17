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

#include "tr1_2d_supg2.h"
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
#include "fei2dtrlin.h"
#include "fei2dquadlin.h"
#include "geotoolbox.h"
#include "contextioerr.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
#endif

namespace oofem {
#define TRSUPG_ZERO_VOF 1.e-8
#define POINT_TOL 1.e-8


//#define TR1_2D_SUPG2_DEBUG

TR1_2D_SUPG2 :: TR1_2D_SUPG2(int n, Domain *aDomain) :
    TR1_2D_SUPG(n, aDomain)
    // Constructor.
{
    numberOfDofMans  = 3;
}

TR1_2D_SUPG2 :: ~TR1_2D_SUPG2()
// Destructor
{ }


void
TR1_2D_SUPG2 :: computeNMtrx(FloatArray &answer, GaussPoint *gp)
{
    double l1, l2;
    answer.resize(3);

    answer.at(1) = l1 = gp->giveCoordinate(1);
    answer.at(2) = l2 = gp->giveCoordinate(2);
    answer.at(3) = 1.0 - l1 - l2;
}


int
TR1_2D_SUPG2 :: computeNumberOfDofs(EquationID ut)
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
TR1_2D_SUPG2 :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
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
TR1_2D_SUPG2 :: giveElementDofIDMask(EquationID ut, IntArray &answer) const
{
    this->giveDofManDofIDMask(1, ut, answer);
}


IRResultType
TR1_2D_SUPG2 :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;               // Required by IR_GIVE_FIELD macro

    this->SUPGElement :: initializeFrom(ir);

    this->vof = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, vof, IFT_TR12DSUPG_pvof, "pvof");
    if ( vof > 0.0 ) {
        setPermanentVolumeFraction(vof);
        this->temp_vof = this->vof;
    } else {
        this->vof = 0.0;
        IR_GIVE_OPTIONAL_FIELD(ir, vof, IFT_TR12DSUPG_vof, "vof");
        this->temp_vof = this->vof;
    }

    this->mat [ 0 ] = this->mat [ 1 ] = this->material;
    IR_GIVE_OPTIONAL_FIELD(ir, mat [ 0 ], IFT_TR12DSUPG2_mat0, "mat0");
    IR_GIVE_OPTIONAL_FIELD(ir, mat [ 1 ], IFT_TR12DSUPG2_mat1, "mat1");
    this->material = this->mat [ 0 ];

    this->computeGaussPoints();
    this->initGeometry();
    this->updateIntegrationRules();
    return IRRT_OK;
}

void
TR1_2D_SUPG2 :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 2;
        integrationRulesArray = new IntegrationRule * [ 2 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3, true);
        integrationRulesArray [ 1 ] = new GaussIntegrationRule(2, this, 1, 3, true);
    }
}


/*
 * Integration template
 * // loop over each fluid
 * for (ifluid = 0; ifluid< 2; ifluid++) {
 * for (ip=0 ; ip < sub_IPRule[ifluid]->getNumberOfIntegrationPoints() ; ip++) {
 *  gp = sub_IPRule[ifluid]->getIntegrationPoint(ip) ;
 *  mapped_gp = sub_mapped_IPRule[ifluid]->getIntegrationPoint(ip) ;
 *  this->computeNMtrx (n, mapped_gp);
 *  dV = this->computeVolumeAround(gp,id[ifluid], vcoords[ifluid]) ;
 *    // compute integral here
 * }
 * }
 */



void
TR1_2D_SUPG2 :: computeAccelerationTerm_MB(FloatMatrix &answer, TimeStep *atTime)
{
    answer.resize(6, 6);
    answer.zero();
    FloatArray n, un;
    int ip, i, j, ifluid;
    double dV, val, rho;

    double u1, u2;

    GaussPoint *gp;
    //IntegrationRule* iRule = integrationRulesArray[giveDefaultIntegrationRule()];

    // consistent mass
    // loop over each fluid
    for ( ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( ip = 0; ip < integrationRulesArray [ ifluid ]->getNumberOfIntegrationPoints(); ip++ ) {
            gp = integrationRulesArray [ ifluid ]->getIntegrationPoint(ip);
            rho = this->_giveMaterial(ifluid)->giveCharacteristicValue(MRM_Density, gp, atTime);
            this->computeNMtrx(n, gp);
            dV = this->computeVolumeAround(gp, id [ ifluid ], vcoords [ ifluid ]);
            // compute integral here
            for ( i = 1; i <= 3; i++ ) {
                for ( j = 1; j <= 3; j++ ) {
                    val = rho * n.at(i) * n.at(j) * dV;
                    answer.at(i * 2 - 1, j * 2 - 1) += val;
                    answer.at(i * 2, j * 2)     += val;
                }
            }
        }
    }

    // SUPG stabilization term
    // loop over each fluid
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime->givePreviousStep(), un);
    for ( ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( ip = 0; ip < integrationRulesArray [ ifluid ]->getNumberOfIntegrationPoints(); ip++ ) {
            gp = integrationRulesArray [ ifluid ]->getIntegrationPoint(ip);
            rho = this->_giveMaterial(ifluid)->giveCharacteristicValue(MRM_Density, gp, atTime);
            this->computeNMtrx(n, gp);
            dV = this->computeVolumeAround(gp, id [ ifluid ], vcoords [ ifluid ]);
            u1 = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
            u2 = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);
            for ( i = 1; i <= 3; i++ ) {
                for ( j = 1; j <= 3; j++ ) {
                    val = dV * rho * t_supg * ( u1 * b [ i - 1 ] + u2 * c [ i - 1 ] ) * n.at(j);
                    answer.at(i * 2, j * 2) += val;
                    answer.at(i * 2 - 1, j * 2 - 1) += val;
                }
            }
        }
    }

    /*
     * // orginal version
     * // consistent mass
     * for (ip=0 ; ip < iRule->getNumberOfIntegrationPoints() ; ip++) {
     * gp = iRule->getIntegrationPoint(ip) ;
     * this->computeNMtrx (n, gp);
     * dV = 2.0*area*gp->giveWeight();
     * for (i=1;i<=3;i++) {
     *  for (j=1; j<=3; j++) {
     *    val = rho*n.at(i)*n.at(j)*dV;
     *    answer.at(i*2-1,j*2-1) += val;
     *    answer.at(i*2,j*2)     += val;
     *  }
     * }
     * }
     *
     * // SUPG stabilization term
     * this -> computeVectorOf(EID_MomentumBalance,VM_Total,atTime->givePreviousStep(), un) ;
     * for (ip=0 ; ip < iRule->getNumberOfIntegrationPoints() ; ip++) {
     * gp = iRule->getIntegrationPoint(ip) ;
     * this->computeNMtrx (n, gp);
     * dV = 2.0*area*gp->giveWeight();
     * u1=n.at(1)*un.at(1)+n.at(2)*un.at(3)+n.at(3)*un.at(5);
     * u2=n.at(1)*un.at(2)+n.at(2)*un.at(4)+n.at(3)*un.at(6);
     * for (i=1;i<=3;i++) {
     *  for (j=1; j<=3; j++) {
     *    val = dV*rho*t_supg*(u1*b[i-1]+u2*c[i-1])*n.at(j);
     *    answer.at(i*2,j*2) += val;
     *    answer.at(i*2-1, j*2-1) += val;
     *  }
     * }
     * }
     */
#ifdef TR1_2D_SUPG2_DEBUG
    /* test */
    FloatMatrix test;
    IntegrationRule *__ir0 = integrationRulesArray [ 0 ];
    integrationRulesArray [ 0 ] = integrationRulesArray [ 1 ];
    TR1_2D_SUPG :: computeAccelerationTerm_MB(test, atTime);
    integrationRulesArray [ 0 ] = __ir0;
    for ( i = 1; i <= 6; i++ ) {
        for ( j = 1; j <= 6; j++ ) {
            if ( fabs( ( answer.at(i, j) - test.at(i, j) ) / test.at(i, j) ) >= 1.e-10 ) {
                _error("computeAccelerationTerm_MB: test failure");
            }
        }
    }

#endif
}


void
TR1_2D_SUPG2 :: computeAdvectionTerm_MB(FloatArray &answer, TimeStep *atTime)
{
    answer.resize(6);
    answer.zero();

    FloatArray n, u, un;
    int ip, i, ifluid;
    double rho;
    double dV, dudx, dudy, dvdx, dvdy, u1, u2;
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime->givePreviousStep(), un);
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);

    GaussPoint *gp;
    //IntegrationRule* iRule = integrationRulesArray[giveDefaultIntegrationRule()];


    dudx = b [ 0 ] * u.at(1) + b [ 1 ] * u.at(3) + b [ 2 ] * u.at(5);
    dudy = c [ 0 ] * u.at(1) + c [ 1 ] * u.at(3) + c [ 2 ] * u.at(5);
    dvdx = b [ 0 ] * u.at(2) + b [ 1 ] * u.at(4) + b [ 2 ] * u.at(6);
    dvdy = c [ 0 ] * u.at(2) + c [ 1 ] * u.at(4) + c [ 2 ] * u.at(6);

    // standard galerkin term
    for ( ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( ip = 0; ip < integrationRulesArray [ ifluid ]->getNumberOfIntegrationPoints(); ip++ ) {
            gp = integrationRulesArray [ ifluid ]->getIntegrationPoint(ip);
            rho = this->_giveMaterial(ifluid)->giveCharacteristicValue(MRM_Density, gp, atTime);
            this->computeNMtrx(n, gp);
            dV = this->computeVolumeAround(gp, id [ ifluid ], vcoords [ ifluid ]);

            u1 = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
            u2 = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);
            for ( i = 1; i <= 3; i++ ) {
                answer.at(i * 2 - 1) += rho * dV * n.at(i) * ( u1 * dudx + u2 * dudy );
                answer.at(i * 2)   += rho * dV * n.at(i) * ( u1 * dvdx + u2 * dvdy );
            }
        }
    }

    // supg stabilization term
    for ( ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( ip = 0; ip < integrationRulesArray [ ifluid ]->getNumberOfIntegrationPoints(); ip++ ) {
            gp = integrationRulesArray [ ifluid ]->getIntegrationPoint(ip);
            rho = this->_giveMaterial(ifluid)->giveCharacteristicValue(MRM_Density, gp, atTime);
            this->computeNMtrx(n, gp);
            dV = this->computeVolumeAround(gp, id [ ifluid ], vcoords [ ifluid ]);

            u1 = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
            u2 = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);
            for ( i = 1; i <= 3; i++ ) {
                answer.at(i * 2 - 1) += t_supg * rho * dV * ( u1 * b [ i - 1 ] + u2 * c [ i - 1 ] ) * ( u1 * dudx + u2 * dudy );
                answer.at(i * 2)   += t_supg * rho * dV * ( u1 * b [ i - 1 ] + u2 * c [ i - 1 ] ) * ( u1 * dvdx + u2 * dvdy );
            }
        }
    }

#if 0
    // test of linearization
    FloatMatrix _h(6, 6);
    FloatArray _t(6);
    this->computeAdvectionDerivativeTerm_MB(_h, atTime);
    _t.beProductOf(_h, u);
    for ( i = 1; i <= 6; i++ ) {
        if ( ( fabs( answer.at(i) - _t.at(i) ) >= 1.e-6 ) ) {
            _error3( "computeAdvectionTerm_MB: test failure (elem %d, error=%e)", this->number, fabs( answer.at(i) - _t.at(i) ) );
        }
    }

#endif


#ifdef TR1_2D_SUPG2_DEBUG
    /* test */
    FloatArray test;
    IntegrationRule *__ir0 = integrationRulesArray [ 0 ];
    integrationRulesArray [ 0 ] = integrationRulesArray [ 1 ];
    TR1_2D_SUPG :: computeAdvectionTerm_MB(test, atTime);
    integrationRulesArray [ 0 ] = __ir0;
    for ( i = 1; i <= 6; i++ ) {
        if ( fabs( ( answer.at(i) - test.at(i) ) / test.at(i) ) >= 1.e-10 ) {
            _error("computeAdvectionTerm_MB: test failure");
        }
    }

#endif
}


void
TR1_2D_SUPG2 :: computeAdvectionDerivativeTerm_MB(FloatMatrix &answer, TimeStep *atTime)
{
    answer.resize(6, 6);
    answer.zero();

    FloatArray u, un, n;
    double rho;
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime->givePreviousStep(), un);

    double u1, u2, dV;
    int i, j, k, m, w_dof_addr, u_dof_addr, dij, ip, ifluid;

    GaussPoint *gp;
    //IntegrationRule* iRule = integrationRulesArray[giveDefaultIntegrationRule()];

    // dN(v)/dv
    for ( ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( ip = 0; ip < integrationRulesArray [ ifluid ]->getNumberOfIntegrationPoints(); ip++ ) {
            gp = integrationRulesArray [ ifluid ]->getIntegrationPoint(ip);
            rho = this->_giveMaterial(ifluid)->giveCharacteristicValue(MRM_Density, gp, atTime);
            this->computeNMtrx(n, gp);
            dV = this->computeVolumeAround(gp, id [ ifluid ], vcoords [ ifluid ]);

            u1 = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
            u2 = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);
            for ( i = 1; i <= 2; i++ ) { // test function index
                for ( k = 1; k <= 3; k++ ) { // nodal val of test function w
                    for ( j = 1; j <= 2; j++ ) { // velocity vector component
                        for ( m = 1; m <= 3; m++ ) { // nodal component
                            w_dof_addr = ( k - 1 ) * 2 + i;
                            u_dof_addr = ( m - 1 ) * 2 + j;
                            dij = ( i == j );
                            answer.at(w_dof_addr, u_dof_addr) += rho * dV * n.at(k) * ( u1 * dij * b [ m - 1 ] + u2 * dij * c [ m - 1 ] );
                        }
                    }
                }
            }
        }
    }

    // stabilization term dN_delta/du
    for ( ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( ip = 0; ip < integrationRulesArray [ ifluid ]->getNumberOfIntegrationPoints(); ip++ ) {
            gp = integrationRulesArray [ ifluid ]->getIntegrationPoint(ip);
            rho = this->_giveMaterial(ifluid)->giveCharacteristicValue(MRM_Density, gp, atTime);
            this->computeNMtrx(n, gp);
            dV = this->computeVolumeAround(gp, id [ ifluid ], vcoords [ ifluid ]);

            u1 = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
            u2 = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);
            for ( i = 1; i <= 2; i++ ) { // test function index
                for ( k = 1; k <= 3; k++ ) { // nodal val of function w
                    for ( j = 1; j <= 2; j++ ) { // velocity vector component
                        for ( m = 1; m <= 3; m++ ) { //  nodal components
                            w_dof_addr = ( k - 1 ) * 2 + i;
                            u_dof_addr = ( m - 1 ) * 2 + j;
                            dij = ( i == j );
                            answer.at(w_dof_addr, u_dof_addr) += t_supg * rho * dV * ( u1 * b [ k - 1 ] + u2 * c [ k - 1 ] ) * ( u1 * dij * b [ m - 1 ] + u2 * dij * c [ m - 1 ] );
                        }
                    }
                }
            }
        }
    }

#ifdef TR1_2D_SUPG2_DEBUG
    /* test */
    FloatMatrix test;
    IntegrationRule *__ir0 = integrationRulesArray [ 0 ];
    integrationRulesArray [ 0 ] = integrationRulesArray [ 1 ];
    TR1_2D_SUPG :: computeAdvectionDerivativeTerm_MB(test, atTime);
    integrationRulesArray [ 0 ] = __ir0;
    for ( i = 1; i <= 6; i++ ) {
        for ( j = 1; j <= 6; j++ ) {
            if ( fabs( ( answer.at(i, j) - test.at(i, j) ) / test.at(i, j) ) >= 1.e-8 ) {
                _error2( "computeAdvectionDerivativeTerm_MB: test failure (err=%e)", ( answer.at(i, j) - test.at(i, j) ) / test.at(i, j) );
            }
        }
    }

#endif
}


void
TR1_2D_SUPG2 :: computeDiffusionTerm_MB(FloatArray &answer, TimeStep *atTime)
{
    int i, ip, ifluid;
    answer.resize(6);
    answer.zero();
    FloatArray u, un, eps(3), stress;
    double dV, Re = domain->giveEngngModel()->giveUnknownComponent(ReynoldsNumber, VM_Unknown, atTime, domain, NULL);
    //double dudx,dudy,dvdx,dvdy;
    GaussPoint *gp;

    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);

    eps.at(1) = ( b [ 0 ] * u.at(1) + b [ 1 ] * u.at(3) + b [ 2 ] * u.at(5) );
    eps.at(2) = ( c [ 0 ] * u.at(2) + c [ 1 ] * u.at(4) + c [ 2 ] * u.at(6) );
    eps.at(3) = ( b [ 0 ] * u.at(2) + b [ 1 ] * u.at(4) + b [ 2 ] * u.at(6) + c [ 0 ] * u.at(1) + c [ 1 ] * u.at(3) + c [ 2 ] * u.at(5) );

    /*
     * this -> computeVectorOf(EID_MomentumBalance,VM_Total,atTime->givePreviousStep(), un) ;
     * dudx = b[0]*un.at(1)+b[1]*un.at(3)+b[2]*un.at(5);
     * dudy = c[0]*un.at(1)+c[1]*un.at(3)+c[2]*un.at(5);
     * dvdx = b[0]*un.at(2)+b[1]*un.at(4)+b[2]*un.at(6);
     * dvdy = c[0]*un.at(2)+c[1]*un.at(4)+c[2]*un.at(6);
     */

    for ( ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( ip = 0; ip < integrationRulesArray [ ifluid ]->getNumberOfIntegrationPoints(); ip++ ) {
            gp = integrationRulesArray [ ifluid ]->getIntegrationPoint(ip);
            dV = this->computeVolumeAround(gp, id [ ifluid ], vcoords [ ifluid ]);

            ( ( FluidDynamicMaterial * ) this->_giveMaterial(ifluid) )->computeDeviatoricStressVector(stress, gp, eps, atTime);
            stress.times(1. / Re);

            // \int dNu/dxj \Tau_ij
            answer.resize(6);
            for ( i = 0; i < 3; i++ ) {
                //rh1p(1,lok) = -geome(7,ia)*( sigxx(ia)*geome(lok,ia) + sigxy(ia)*geome(lok1,ia) )*0.5d+00;
                //rh1p(2,lok) = -geome(7,ia)*( sigxy(ia)*geome(lok,ia) + sigyy(ia)*geome(lok1,ia) )*0.5d+00;
                answer.at( ( i ) * 2 + 1 ) += dV * ( stress.at(1) * b [ i ] + stress.at(3) * c [ i ] );
                answer.at( ( i + 1 ) * 2 ) += dV * ( stress.at(3) * b [ i ] + stress.at(2) * c [ i ] );

                /*
                 * // stabilization term k_delta
                 * answer.at((i)*2+1) += t_supg* dV * (stress.at(1)*(dudx*b[i]+dvdx*c[i]) + stress.at(3)*(dudy*b[i]+dvdy*c[i]));
                 * answer.at((i+1)*2) += t_supg* dV * (stress.at(3)*(dudx*b[i]+dvdx*c[i]) + stress.at(2)*(dudy*b[i]+dvdy*c[i]));
                 */
            }
        }
    }

#ifdef TR1_2D_SUPG2_DEBUG
    /* test */
    FloatArray test;
    IntegrationRule *__ir0 = integrationRulesArray [ 0 ];
    integrationRulesArray [ 0 ] = integrationRulesArray [ 1 ];
    TR1_2D_SUPG :: computeDiffusionTerm_MB(test, atTime);
    integrationRulesArray [ 0 ] = __ir0;
    for ( i = 1; i <= 6; i++ ) {
        if ( fabs( ( answer.at(i) - test.at(i) ) / test.at(i) ) >= 1.e-10 ) {
            _error("computeDiffusionTerm_MB: test failure");
        }
    }

#endif
}


void
TR1_2D_SUPG2 :: computeDiffusionDerivativeTerm_MB(FloatMatrix &answer, MatResponseMode mode, TimeStep *atTime)
{
    int ifluid, ip;
    double dV;
    //double dudx, dudy, dvdx, dvdy;
    answer.resize(6, 6);
    answer.zero();
    FloatMatrix _db, _d, _b(3, 6), _bs(3, 6);
    //FloatArray un;
    GaussPoint *gp;
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

    /*
     * // stabilization terms
     * this -> computeVectorOf(EID_MomentumBalance,VM_Total,atTime->givePreviousStep(), un) ;
     * dudx = b[0]*un.at(1)+b[1]*un.at(3)+b[2]*un.at(5);
     * dudy = c[0]*un.at(1)+c[1]*un.at(3)+c[2]*un.at(5);
     * dvdx = b[0]*un.at(2)+b[1]*un.at(4)+b[2]*un.at(6);
     * dvdy = c[0]*un.at(2)+c[1]*un.at(4)+c[2]*un.at(6);
     *
     * _bs.at(1,1)=dudx*b[0]+dvdx*c[0]; _bs.at(1,2)=0.; _bs.at(1,3)=dudx*b[1]+dvdx*c[1]; _bs.at(1,4)=0.; _bs.at(1,5)=dudx*b[2]+dvdx*c[2]; _bs.at(1,6)=0.;
     * _bs.at(2,1)=0.; _bs.at(2,2)=dudy*b[0]+dvdy*c[0]; _bs.at(2,3)=0.; _bs.at(2,4)=dudy*b[1]+dvdy*c[1]; _bs.at(2,5)=0.; _bs.at(2,6)=dudy*b[2]+dvdy*c[2];
     * _bs.at(3,1)=dudy*b[0]+dvdy*c[0]; _bs.at(3,2)=dudx*b[0]+dvdx*c[0];
     * _bs.at(3,3)=dudy*b[1]+dvdy*c[1]; _bs.at(3,4)=dudx*b[1]+dvdx*c[1];
     * _bs.at(3,5)=dudy*b[2]+dvdy*c[2]; _bs.at(3,6)=dudx*b[2]+dvdx*c[2];
     */

    for ( ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( ip = 0; ip < integrationRulesArray [ ifluid ]->getNumberOfIntegrationPoints(); ip++ ) {
            gp = integrationRulesArray [ ifluid ]->getIntegrationPoint(ip);
            dV = this->computeVolumeAround(gp, id [ ifluid ], vcoords [ ifluid ]);


            ( ( FluidDynamicMaterial * ) this->_giveMaterial(ifluid) )->giveDeviatoricStiffnessMatrix(_d, mode,
                                                                                                      gp, atTime);
            _db.beProductOf(_d, _b);
            answer.plusProductSymmUpper(_b, _db, dV);
            //answer.plusProductSymmUpper (_bs,_db,dV*t_supg);

            answer.symmetrized();
        }
    }

    answer.times(1. / Re);

#ifdef TR1_2D_SUPG2_DEBUG
    /* test */
    int i, j;
    FloatMatrix test;
    IntegrationRule *__ir0 = integrationRulesArray [ 0 ];
    integrationRulesArray [ 0 ] = integrationRulesArray [ 1 ];
    TR1_2D_SUPG :: computeDiffusionDerivativeTerm_MB(test, mode, atTime);
    integrationRulesArray [ 0 ] = __ir0;
    for ( i = 1; i <= 6; i++ ) {
        for ( j = 1; j <= 6; j++ ) {
            if ( fabs( ( answer.at(i, j) - test.at(i, j) ) / test.at(i, j) ) >= 1.e-8 ) {
                _error2( "computeDiffusionDerivativeTerm_MB: test failure (err=%e)", ( answer.at(i, j) - test.at(i, j) ) / test.at(i, j) );
            }
        }
    }

#endif
}


void
TR1_2D_SUPG2 :: computePressureTerm_MB(FloatMatrix &answer, TimeStep *atTime)
{
    // TR1_2D_SUPG :: computePressureTerm_MB(answer, atTime);
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
TR1_2D_SUPG2 :: computeLSICStabilizationTerm_MB(FloatMatrix &answer, TimeStep *atTime)
{
    answer.resize(6, 6);
    answer.zero();
    double dV, rho;
    //double coeff = area*t_lsic*rho;
    double n[] = {
        b [ 0 ], c [ 0 ], b [ 1 ], c [ 1 ], b [ 2 ], c [ 2 ]
    };
    int i, j, ip, ifluid;
    GaussPoint *gp;

    for ( ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( ip = 0; ip < integrationRulesArray [ ifluid ]->getNumberOfIntegrationPoints(); ip++ ) {
            gp = integrationRulesArray [ ifluid ]->getIntegrationPoint(ip);
            rho = this->_giveMaterial(ifluid)->giveCharacteristicValue(MRM_Density, gp, atTime);
            dV = this->computeVolumeAround(gp, id [ ifluid ], vcoords [ ifluid ]);


            for ( i = 1; i <= 6; i++ ) {
                for ( j = 1; j <= 6; j++ ) {
                    answer.at(i, j) += dV * t_lsic * rho * n [ i - 1 ] * n [ j - 1 ];
                }
            }
        }
    }

#ifdef TR1_2D_SUPG2_DEBUG
    /* test */
    FloatMatrix test;
    IntegrationRule *__ir0 = integrationRulesArray [ 0 ];
    integrationRulesArray [ 0 ] = integrationRulesArray [ 1 ];
    TR1_2D_SUPG :: computeLSICStabilizationTerm_MB(test, atTime);
    integrationRulesArray [ 0 ] = __ir0;
    for ( i = 1; i <= 6; i++ ) {
        for ( j = 1; j <= 6; j++ ) {
            if ( fabs( ( answer.at(i, j) - test.at(i, j) ) / test.at(i, j) ) >= 1.e-8 ) {
                _error2( "computeLSICStabilizationTerm_MB: test failure (err=%e)", ( answer.at(i, j) - test.at(i, j) ) / test.at(i, j) );
            }
        }
    }

#endif
}


void
TR1_2D_SUPG2 :: computeLinearAdvectionTerm_MC(FloatMatrix &answer, TimeStep *atTime)
{
    TR1_2D_SUPG :: computeLinearAdvectionTerm_MC(answer, atTime);
}

void
TR1_2D_SUPG2 :: computeAdvectionTerm_MC(FloatArray &answer, TimeStep *atTime)
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
TR1_2D_SUPG2 :: computeAdvectionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *atTime)
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
TR1_2D_SUPG2 :: computeAccelerationTerm_MC(FloatMatrix &answer, TimeStep *atTime)
{
    TR1_2D_SUPG :: computeAccelerationTerm_MC(answer, atTime);
}

void
TR1_2D_SUPG2 :: computePressureTerm_MC(FloatMatrix &answer, TimeStep *atTime)
{
    double dV, rho;
    int i, j, ip, ifluid;
    //double coeff = t_pspg*area/rho;
    answer.resize(3, 3);
    answer.zero();

    GaussPoint *gp;

    for ( ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( ip = 0; ip < integrationRulesArray [ ifluid ]->getNumberOfIntegrationPoints(); ip++ ) {
            gp = integrationRulesArray [ ifluid ]->getIntegrationPoint(ip);
            rho = this->_giveMaterial(ifluid)->giveCharacteristicValue(MRM_Density, gp, atTime);
            dV = this->computeVolumeAround(gp, id [ ifluid ], vcoords [ ifluid ]);


            for ( i = 1; i <= 3; i++ ) {
                for ( j = 1; j <= 3; j++ ) {
                    answer.at(i, j) += t_pspg * dV * ( b [ i - 1 ] * b [ j - 1 ] + c [ i - 1 ] * c [ j - 1 ] ) / rho;
                }
            }
        }
    }

#ifdef TR1_2D_SUPG2_DEBUG
    /* test */
    FloatMatrix test;
    IntegrationRule *__ir0 = integrationRulesArray [ 0 ];
    integrationRulesArray [ 0 ] = integrationRulesArray [ 1 ];
    TR1_2D_SUPG :: computePressureTerm_MC(test, atTime);
    integrationRulesArray [ 0 ] = __ir0;
    for ( i = 1; i <= 3; i++ ) {
        for ( j = 1; j <= 3; j++ ) {
            if ( fabs( ( answer.at(i, j) - test.at(i, j) ) / test.at(i, j) ) >= 1.e-8 ) {
                _error2( "computePressureTerm_MC: test failure (err=%e)", ( answer.at(i, j) - test.at(i, j) ) / test.at(i, j) );
            }
        }
    }

#endif
}

/////////////////////////////////////////////////////////////////////////////////////////////
void
TR1_2D_SUPG2 :: computeBCRhsTerm_MB(FloatArray &answer, TimeStep *atTime)
{
    answer.resize(6);
    answer.zero();

    int i, l, ip, ifluid, nLoads;
    Load *load;
    bcGeomType ltype;
    double usum [ 2 ], rho;
    FloatArray un, n, gVector;
    double u1, u2, dV;

    GaussPoint *gp;
    //IntegrationRule* iRule = integrationRulesArray[giveDefaultIntegrationRule()];

    // add body load (gravity) termms
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime->givePreviousStep(), un);

    usum [ 0 ] = un.at(1) + un.at(3) + un.at(5);
    usum [ 1 ] = un.at(2) + un.at(4) + un.at(6);

    nLoads    = this->giveBodyLoadArray()->giveSize();
    for ( l = 1; l <= nLoads; l++ ) {
        load  = domain->giveLoad( bodyLoadArray.at(l) );
        ltype = load->giveBCGeoType();
        if ( ( ltype == BodyLoadBGT ) && ( load->giveBCValType() == ForceLoadBVT ) ) {
            load->computeComponentArrayAt(gVector, atTime, VM_Total);
            if ( gVector.giveSize() ) {
                for ( ifluid = 0; ifluid < 2; ifluid++ ) {
                    for ( ip = 0; ip < integrationRulesArray [ ifluid ]->getNumberOfIntegrationPoints(); ip++ ) {
                        gp = integrationRulesArray [ ifluid ]->getIntegrationPoint(ip);

                        rho = this->_giveMaterial(ifluid)->giveCharacteristicValue(MRM_Density, gp, atTime);
                        dV = this->computeVolumeAround(gp, id [ ifluid ], vcoords [ ifluid ]);
                        this->computeNMtrx(n, gp);

                        u1 = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
                        u2 = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);

                        for ( i = 1; i <= 3; i++ ) {
                            answer.at(i * 2 - 1) += rho * dV * gVector.at(1) * ( n.at(i) + t_supg * ( u1 * b [ i - 1 ] + u2 * c [ i - 1 ] ) );
                            answer.at(i * 2)   += rho * dV * gVector.at(2) * ( n.at(i) + t_supg * ( u1 * b [ i - 1 ] + u2 * c [ i - 1 ] ) );
                        }
                    }
                }
            }
        }
    }

    // loop over sides
    int n1, n2;
    int lnum, id;
    double tx, ty, length, nx, ny;
    FloatArray t, coords(1);
    BoundaryLoad *bload;

    // loop over boundary load array
    nLoads    = this->giveBoundaryLoadArray()->giveSize() / 2;
    for ( i = 1; i <= nLoads; i++ ) {
        if ( load->giveBCValType() == ForceLoadBVT ) {
            lnum  = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
            id    = boundaryLoadArray.at(i * 2);

            // integrate tractions
            n1 = id;
            n2 = ( n1 == 3 ? 1 : n1 + 1 );

            tx = giveNode(n2)->giveCoordinate(1) - giveNode(n1)->giveCoordinate(1);
            ty = giveNode(n2)->giveCoordinate(2) - giveNode(n1)->giveCoordinate(2);
            length = sqrt(tx * tx + ty * ty);
            nx = ty / length;
            ny = -tx / length;

            bload  = dynamic_cast< BoundaryLoad * >( domain->giveLoad(lnum) );
            if ( bload ) {
                bload->computeValueAt(t, atTime, coords, VM_Total);

                // here it is assumed constant traction, one point integration only
                // n1 (u,v)
                answer.at( ( n1 - 1 ) * 2 + 1 ) += t.at(1)  * l / 2.;
                answer.at(n1 * 2)       += t.at(2)  * l / 2.;
                // n2 (u,v)
                answer.at( ( n2 - 1 ) * 2 + 1 ) += t.at(1)  * l / 2.;
                answer.at(n2 * 2)       += t.at(2)  * l / 2.;

                //answer.at(n1)+= (t.at(1)*nx + t.at(2)*ny) * length/2.;
                //answer.at(n2)+= (t.at(1)*nx + t.at(2)*ny) * length/2.;
            }
        }
    }

#ifdef TR1_2D_SUPG2_DEBUG
    /* test */
    FloatArray test;
    IntegrationRule *__ir0 = integrationRulesArray [ 0 ];
    integrationRulesArray [ 0 ] = integrationRulesArray [ 1 ];
    TR1_2D_SUPG :: computeBCRhsTerm_MB(test, atTime);
    integrationRulesArray [ 0 ] = __ir0;
    for ( i = 1; i <= 6; i++ ) {
        if ( fabs( ( answer.at(i) - test.at(i) ) / test.at(i) ) >= 1.e-10 ) {
            _error("computeBCRhsTerm_MB: test failure");
        }
    }

#endif
}

void
TR1_2D_SUPG2 :: computeBCRhsTerm_MC(FloatArray &answer, TimeStep *atTime)
{
    TR1_2D_SUPG :: computeBCRhsTerm_MC(answer, atTime);
}


void
TR1_2D_SUPG2 :: updateStabilizationCoeffs(TimeStep *atTime)
{
    //TR1_2D_SUPG :: updateStabilizationCoeffs (atTime);
#if 0
    int i, j, k, l, w_dof_addr, u_dof_addr, ip, ifluid;
    double __g_norm, __gamma_norm, __gammav_norm, __beta_norm, __betav_norm, __c_norm, __e_norm, __k_norm, __Re;
    double __t_p1, __t_p2, __t_p3, __t_pv1, __t_pv2, __t_pv3;
    double nu, nu0, nu1, usum, vsum, rho, dV, u1, u2;
    FloatArray u, un, a;

    // compute averaged viscosity based on rule of mixture
    GaussPoint *gp;
    if ( integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints() ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    } else {
        gp = integrationRulesArray [ 1 ]->getIntegrationPoint(0);
    }

    nu0 = this->_giveMaterial(0)->giveCharacteristicValue(MRM_Viscosity, gp, atTime);
    nu1 = this->_giveMaterial(1)->giveCharacteristicValue(MRM_Viscosity, gp, atTime);
    nu = vof * nu0 + ( 1. - vof ) * nu1;

    //this -> computeVectorOf(EID_MomentumBalance,VM_Total,atTime->givePreviousStep(),un) ;
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);
    this->computeVectorOf(EID_MomentumBalance, VM_Acceleration, atTime, a);

    un = u;
    usum = un.at(1) + un.at(3) + un.at(5);
    vsum = un.at(2) + un.at(4) + un.at(6);

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
    for ( k = 1; k <= 3; k++ ) {
        for ( l = 1; l <= 3; l++ ) {
            __tmp.at(k, l * 2 - 1) = ar3 * b [ k - 1 ] * ( usum * b [ l - 1 ] + vsum * c [ l - 1 ] );
            __tmp.at(k, l * 2)  = ar3 * c [ k - 1 ] * ( usum * b [ l - 1 ] + vsum * c [ l - 1 ] );
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
    for ( ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( ip = 0; ip < integrationRulesArray [ ifluid ]->getNumberOfIntegrationPoints(); ip++ ) {
            gp = integrationRulesArray [ ifluid ]->getIntegrationPoint(ip);
            rho = this->_giveMaterial(ifluid)->giveCharacteristicValue(MRM_Density, gp, atTime);
            this->computeNMtrx(n, gp);
            dV = this->computeVolumeAround(gp, id [ ifluid ], vcoords [ ifluid ]);

            u1 = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
            u2 = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);
            for ( i = 1; i <= 2; i++ ) {
                for ( k = 1; k <= 3; k++ ) {
                    for ( l = 1; l <= 3; l++ ) {
                        w_dof_addr = ( k - 1 ) * 2 + i;
                        u_dof_addr = ( l - 1 ) * 2 + i;
                        __tmp.at(w_dof_addr, u_dof_addr) += rho * dV * n.at(k) * ( u1 * b [ l - 1 ] + u2 * c [ l - 1 ] );
                    }
                }
            }
        }
    }

    __c_norm = __tmp.computeFrobeniusNorm();
    // compute e mtrx (advection term of momentum balance)
    __tmp.resize(6, 6);
    __tmp.zero();
    double __n[] = {
        b [ 0 ], c [ 0 ], b [ 1 ], c [ 1 ], b [ 2 ], c [ 2 ]
    };
    for ( ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( ip = 0; ip < integrationRulesArray [ ifluid ]->getNumberOfIntegrationPoints(); ip++ ) {
            gp = integrationRulesArray [ ifluid ]->getIntegrationPoint(ip);
            rho = this->_giveMaterial(ifluid)->giveCharacteristicValue(MRM_Density, gp, atTime);
            dV = this->computeVolumeAround(gp, id [ ifluid ], vcoords [ ifluid ]);

            for ( i = 1; i <= 6; i++ ) {
                for ( j = 1; j <= 6; j++ ) {
                    __tmp.at(i, j) += dV * rho * __n [ i - 1 ] * __n [ j - 1 ];
                }
            }
        }
    }

    __e_norm = __tmp.computeFrobeniusNorm();
    // compute element level Reynolds number
    // compute k norm first
    __tmp.resize(6, 6);
    __tmp.zero();
    for ( ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( ip = 0; ip < integrationRulesArray [ ifluid ]->getNumberOfIntegrationPoints(); ip++ ) {
            gp = integrationRulesArray [ ifluid ]->getIntegrationPoint(ip);
            rho = this->_giveMaterial(ifluid)->giveCharacteristicValue(MRM_Density, gp, atTime);
            this->computeNMtrx(n, gp);
            dV = this->computeVolumeAround(gp, id [ ifluid ], vcoords [ ifluid ]);

            u1 = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
            u2 = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);
            for ( k = 1; k <= 3; k++ ) {
                for ( l = 1; l <= 3; l++ ) {
                    __tmp.at(k * 2 - 1, l * 2 - 1) += rho * dV * ( u1 * b [ k - 1 ] + u2 * c [ k - 1 ] ) * ( u1 * b [ l - 1 ] + u2 * c [ l - 1 ] );
                    __tmp.at(k * 2, l * 2)     += rho * dV * ( u1 * b [ k - 1 ] + u2 * c [ k - 1 ] ) * ( u1 * b [ l - 1 ] + u2 * c [ l - 1 ] );
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
    } else {
        __Re = vnorm * vnorm * __c_norm / __k_norm / nu;

        __t_p1 = __g_norm / __gamma_norm;
        __t_p2 = atTime->giveTimeIncrement() * __g_norm / 2.0 / __beta_norm;
        __t_p3 = __t_p1 * __Re;
        this->t_pspg = 1. / sqrt( 1. / ( __t_p1 * __t_p1 ) + 1. / ( __t_p2 * __t_p2 ) + 1. / ( __t_p3 * __t_p3 ) );

        __t_pv1 = __t_p1;
        __t_pv2 = __t_pv1 * __gammav_norm / __betav_norm;
        __t_pv3 = __t_pv1 * __Re;
        this->t_supg = 1. / sqrt( 1. / ( __t_pv1 * __t_pv1 ) + 1. / ( __t_pv2 * __t_pv2 ) + 1. / ( __t_pv3 * __t_pv3 ) );

        this->t_lsic = __c_norm / __e_norm;
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
    double nu, nu0, nu1;

    // compute averaged viscosity based on rule of mixture
    GaussPoint *gp;
    if ( integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints() ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    } else {
        gp = integrationRulesArray [ 1 ]->getIntegrationPoint(0);
    }

    nu0 = this->_giveMaterial(0)->giveCharacteristicValue(MRM_Viscosity, gp, atTime->givePreviousStep());
    nu1 = this->_giveMaterial(1)->giveCharacteristicValue(MRM_Viscosity, gp, atTime->givePreviousStep());
    nu = vof * nu0 + ( 1. - vof ) * nu1;
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

    if ( ( vnorm == 0.0 ) || ( sum == 0.0 ) ) {
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



double
TR1_2D_SUPG2 :: computeCriticalTimeStep(TimeStep *tStep)
{
    return 1.e3;
}


Interface *
TR1_2D_SUPG2 :: giveInterface(InterfaceType interface)
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
    }

    return NULL;
}


int
TR1_2D_SUPG2 :: SpatialLocalizerI_containsPoint(const FloatArray &coords) {
    FloatArray lcoords;
    return this->computeLocalCoordinates(lcoords, coords);
}

double
TR1_2D_SUPG2 :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
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
TR1_2D_SUPG2 :: computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    /* one computes here average deviatoric stress, based on rule of mixture (this is used only for postprocessing) */
    int i;
    FloatArray eps(3), s0(3), s1(3);
    answer.resize(3);

    this->computeDeviatoricStrain(eps, gp, tStep);

    ( ( FluidDynamicMaterial * ) this->_giveMaterial(0) )->computeDeviatoricStressVector(s0, gp, eps, tStep);
    ( ( FluidDynamicMaterial * ) this->_giveMaterial(1) )->computeDeviatoricStressVector(s1, gp, eps, tStep);

    for ( i = 1; i <= 3; i++ ) {
        answer.at(i) = ( temp_vof ) * s0.at(i) + ( 1. - temp_vof ) * s1.at(i);
    }
}


/*
 * double
 * TR1_2D_SUPG2 :: computeCriticalTimeStep (TimeStep* tStep)
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


double
TR1_2D_SUPG2 :: computeLEPLICVolumeFraction(const FloatArray &n, const double p, LEPlic *matInterface, bool updFlag)
{
    Polygon pg;
    double answer, volume = computeMyVolume(matInterface, updFlag);
    this->formVolumeInterfacePoly(pg, matInterface, n, p, updFlag);
    answer = fabs(pg.computeVolume() / volume);
    if ( answer > 1.000000001 ) {
        _warning2("VOF fraction out of bounds, vof = %e\n", answer);
        return 1.0;
    } else {
        return answer;
    }
}

void
TR1_2D_SUPG2 :: formMaterialVolumePoly(Polygon &matvolpoly, LEPlic *matInterface,
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
TR1_2D_SUPG2 :: formVolumeInterfacePoly(Polygon &matvolpoly, LEPlic *matInterface,
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
                if ( fabs(sd) > 1.e-10 ) {
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


void
TR1_2D_SUPG2 :: updateVolumePolygons(Polygon &referenceFluidPoly, Polygon &secondFluidPoly, int &rfPoints, int &sfPoints,
                                     const FloatArray &normal, const double p, bool updFlag)
{
    /*
     * this method updates two polygons, one filled with reference fluid and second filled with
     * other fluid (air). These two polygons are used in integrating element contributions.
     */
    int i, next;
    bool nodeIn [ 3 ];
    double nx = normal.at(1), ny = normal.at(2), x, y;
    double tx, ty;
    Vertex v;

    rfPoints = sfPoints = 0;
    referenceFluidPoly.clear();
    secondFluidPoly.clear();

    for ( i = 1; i <= 3; i++ ) {
        x = this->giveNode(i)->giveCoordinate(1);
        y = this->giveNode(i)->giveCoordinate(2);

        if ( ( nx * x + ny * y + p ) >= 0. ) {
            nodeIn [ i - 1 ] = true;
        } else {
            nodeIn [ i - 1 ] = false;
        }
    }

    if ( nodeIn [ 0 ] && nodeIn [ 1 ] && nodeIn [ 2 ] ) { // all nodes inside
        for ( i = 1; i <= 3; i++ ) {
            x = this->giveNode(i)->giveCoordinate(1);
            y = this->giveNode(i)->giveCoordinate(2);

            v.setCoords(x, y);
            referenceFluidPoly.addVertex(v);
            rfPoints++;
        }

        return;
    } else if ( !( nodeIn [ 0 ] || nodeIn [ 1 ] || nodeIn [ 2 ] ) ) { // all nodes outside
        for ( i = 1; i <= 3; i++ ) {
            x = this->giveNode(i)->giveCoordinate(1);
            y = this->giveNode(i)->giveCoordinate(2);

            v.setCoords(x, y);
            secondFluidPoly.addVertex(v);
            sfPoints++;
        }

        return;
    } else {
        for ( i = 1; i <= 3; i++ ) {
            next = i < 3 ? i + 1 : 1;
            if ( nodeIn [ i - 1 ] ) {
                v.setCoords( this->giveNode(i)->giveCoordinate(1),
                            this->giveNode(i)->giveCoordinate(2) );

                referenceFluidPoly.addVertex(v);
                rfPoints++;
            } else {
                v.setCoords( this->giveNode(i)->giveCoordinate(1),
                            this->giveNode(i)->giveCoordinate(2) );

                secondFluidPoly.addVertex(v);
                sfPoints++;
            }

            if ( nodeIn [ next - 1 ] ^ nodeIn [ i - 1 ] ) {
                // compute intersection with (i,next) edge
                x = this->giveNode(i)->giveCoordinate(1);
                y = this->giveNode(i)->giveCoordinate(2);
                tx = this->giveNode(next)->giveCoordinate(1) - x;
                ty = this->giveNode(next)->giveCoordinate(2) - y;

                double s, sd = nx * tx + ny * ty;
                if ( fabs(sd) > 1.e-10 ) {
                    s = ( -p - ( nx * x + ny * y ) ) / sd;
                    v.setCoords(x + tx * s, y + ty * s);
                    referenceFluidPoly.addVertex(v);
                    secondFluidPoly.addVertex(v);
                    rfPoints++;
                    sfPoints++;
                } else {
                    // pathological case - lines are parallel
                    if ( nodeIn [ i - 1 ] ) {
                        v.setCoords( this->giveNode(next)->giveCoordinate(1), this->giveNode(next)->giveCoordinate(2) );

                        referenceFluidPoly.addVertex(v);
                        secondFluidPoly.addVertex(v);
                        rfPoints++;
                        sfPoints++;
                    } else {
                        //v.setCoords(x, y);
                        //referenceFluidPoly.addVertex(v);
                        v.setCoords( this->giveNode(next)->giveCoordinate(1), this->giveNode(next)->giveCoordinate(2) );

                        referenceFluidPoly.addVertex(v);
                        secondFluidPoly.addVertex(v);
                        rfPoints++;
                        sfPoints++;
                    }
                }
            }
        } // end loop over elem nodes

    }
}

double
TR1_2D_SUPG2 :: truncateMatVolume(const Polygon &matvolpoly, double &volume)
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
TR1_2D_SUPG2 :: formMyVolumePoly(Polygon &me, LEPlic *matInterface, bool updFlag)
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
TR1_2D_SUPG2 :: computeMyVolume(LEPlic *matInterface, bool updFlag)
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

void
TR1_2D_SUPG2 :: giveElementCenter(LEPlic *mat_interface, FloatArray &center, bool upd)
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
TR1_2D_SUPG2 :: EIPrimaryFieldI_evaluateFieldVectorAt(FloatArray &answer, PrimaryField &pf,
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
                //_error("EIPrimaryFieldI_evaluateFieldVectorAt: unknown dof id encountered");
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
TR1_2D_SUPG2 :: updateYourself(TimeStep *tStep)
{
    SUPGElement :: updateYourself(tStep);
    LEPlicElementInterface :: updateYourself(tStep);
    //this->updateIntegrationRules ();
}

void
TR1_2D_SUPG2 :: updateIntegrationRules()
{
    int c [ 2 ];
    int i, j, ip;
    double x, y;
    Vertex v;

    for ( i = 0; i < 2; i++ ) {
        myPoly [ i ].clear();
        if ( vcoords [ i ] ) {
            delete vcoords [ i ];
        }

        vcoords [ i ] = NULL;
    }

    if ( this->temp_vof <= TRSUPG_ZERO_VOF ) {
        for ( i = 1; i <= 3; i++ ) {
            x = this->giveNode(i)->giveCoordinate(1);
            y = this->giveNode(i)->giveCoordinate(2);
            v.setCoords(x, y);
            myPoly [ 1 ].addVertex(v);
            c [ 1 ] = 3;
            c [ 0 ] = 0;
        }
    } else if ( this->temp_vof >= ( 1. - TRSUPG_ZERO_VOF ) ) {
        for ( i = 1; i <= 3; i++ ) {
            x = this->giveNode(i)->giveCoordinate(1);
            y = this->giveNode(i)->giveCoordinate(2);
            v.setCoords(x, y);
            myPoly [ 0 ].addVertex(v);
            c [ 0 ] = 3;
            c [ 1 ] = 0;
        }
    } else {
        this->updateVolumePolygons(myPoly [ 0 ], myPoly [ 1 ], c [ 0 ], c [ 1 ], temp_normal,  temp_p, false);
    }

    integrationRulesArray [ 0 ]->clear();
    integrationRulesArray [ 1 ]->clear();

    FloatArray gc, lc;
    const Vertex *p;
    double a;
    FEI2dTrLin triaApprox(1, 2);
    FEI2dQuadLin quadApprox(1, 2);
    FEInterpolation *approx = NULL;
    GaussPoint *gp;
    // set up integration points
    for ( i = 0; i < 2; i++ ) {
        if ( c [ i ] == 3 ) {
            id [ i ] = _Triangle;
            approx = & triaApprox;
        } else if ( c [ i ] == 4 ) {
            id [ i ] = _Square;
            approx = & quadApprox;
        } else if ( c [ i ] == 0 ) {
            continue;
        } else {
            _error2("updateYourself: cannot set up integration domain for %d vertex polygon", c [ i ]);
        }

        if ( i == 0 ) {
            a = area * this->temp_vof;
        } else {
            a = area * ( 1. - this->temp_vof );
        }

        if ( c [ i ] ) {
            vcoords [ i ] = new const FloatArray * [ c [ i ] ];
        }

        // set up vertex coords
        Polygon :: PolygonVertexIterator it(myPoly + i);
        j = 0;
        while ( it.giveNext(& p) ) {
            vcoords [ i ] [ j ] = p->getCoords();
            j++;
        }

        integrationRulesArray [ i ]->setUpIntegrationPoints(id [ i ], 4, _2dFlow);

        // remap ip coords into area coords of receiver
        for ( ip = 0; ip < integrationRulesArray [ i ]->getNumberOfIntegrationPoints(); ip++ ) {
            gp = integrationRulesArray [ i ]->getIntegrationPoint(ip);
            approx->local2global(gc, * gp->giveCoordinates(), FEIVertexListGeometryWrapper(c [ i ], vcoords [ i ]));
            triaApprox.global2local(lc, gc, FEIElementGeometryWrapper(this));
            // modify original ip coords to target ones
            gp->setLocalCoordinates( * gp->giveCoordinates() );
            gp->setCoordinates(lc);
            //gp->setWeight (gp->giveWeight()*a/area);
        }
    }

    // internal test -> compute receiver area
    int ifluid;
    double dV, __area = 0.0;
    for ( ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( ip = 0; ip < integrationRulesArray [ ifluid ]->getNumberOfIntegrationPoints(); ip++ ) {
            gp = integrationRulesArray [ ifluid ]->getIntegrationPoint(ip);
            dV = this->computeVolumeAround(gp, id [ ifluid ], vcoords [ ifluid ]);
            // compute integral here
            __area += dV;
        }
    }

    double __err = fabs(__area - area) / area;
    if ( __err > 1.e-6 ) {
        _warning2("updateIntegrationRules: volume inconsistency (%5.2f)", __err * 100);

        __area = 0.0;
        for ( ifluid = 0; ifluid < 2; ifluid++ ) {
            for ( ip = 0; ip < integrationRulesArray [ ifluid ]->getNumberOfIntegrationPoints(); ip++ ) {
                gp = integrationRulesArray [ ifluid ]->getIntegrationPoint(ip);
                dV = this->computeVolumeAround(gp, id [ ifluid ], vcoords [ ifluid ]);
                // compute integral here
                __area += dV;
            }
        }
    }
}


double
TR1_2D_SUPG2 :: computeVolumeAround(GaussPoint *gp, integrationDomain id, const FloatArray **idpoly)
{
    double weight = gp->giveWeight();

    if ( id == _Triangle ) {
        FEI2dTrLin __interpolation(1, 2);
        return weight *fabs( __interpolation.giveTransformationJacobian(* gp->giveLocalCoordinates(), FEIVertexListGeometryWrapper(3, idpoly)) );
    } else {
        FEI2dQuadLin __interpolation(1, 2);
        double det = fabs( __interpolation.giveTransformationJacobian(* gp->giveLocalCoordinates(), FEIVertexListGeometryWrapper(4, idpoly)) );
        return det * weight;
    }
}


/*
 * double
 * TR1_2D_SUPG2::computeVolumeAround(GaussPoint* gp, integrationDomain id, const Polygon& matvolpoly) {
 * }
 */

int
TR1_2D_SUPG2 :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    if ( type == IST_VOFFraction ) {
        answer.resize(1);
        answer.at(1) = this->giveTempVolumeFraction();
        return 1;
    } else if ( type == IST_Density ) {
        answer.resize(1);
        answer.at(1) = this->giveMaterial()->giveCharacteristicValue(MRM_Density, aGaussPoint, atTime);
        return 1;
    } else {
        return TR1_2D_SUPG :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}

int
TR1_2D_SUPG2 :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type)
{
    if ( ( type == IST_VOFFraction ) || ( type == IST_Density ) ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else {
        return TR1_2D_SUPG :: giveIntVarCompFullIndx(answer, type);
    }
}


InternalStateValueType
TR1_2D_SUPG2 :: giveIPValueType(InternalStateType type)
{
    if ( ( type == IST_VOFFraction ) || ( type == IST_Density ) ) {
        return ISVT_SCALAR;
    } else {
        return TR1_2D_SUPG :: giveIPValueType(type);
    }
}


int
TR1_2D_SUPG2 :: giveIPValueSize(InternalStateType type, GaussPoint *gp)
{
    if ( ( type == IST_VOFFraction ) || ( type == IST_Density ) ) {
        return 1;
    } else {
      return TR1_2D_SUPG::giveIPValueSize(type, gp);
    }
}


int
TR1_2D_SUPG2 :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( ( type == IST_StressTensor ) || ( type == IST_StrainTensor ) ) {
        return 4;
    }

    GaussPoint *gp;
    if ( integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints() ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    } else {
        gp = integrationRulesArray [ 1 ]->getIntegrationPoint(0);
    }

    return this->giveIPValueSize(type, gp);
}


void
TR1_2D_SUPG2 :: ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type)
{
    // evaluates N matrix (interpolation estimated stress matrix)
    // according to Zienkiewicz & Zhu paper
    // N(nsigma, nsigma*nnodes)
    // Definition : sigmaVector = N * nodalSigmaVector
    double l1, l2, l3;

    l1 = aGaussPoint->giveCoordinate(1);
    l2 = aGaussPoint->giveCoordinate(2);
    l3 = 1.0 - l1 - l2;

    if ( this->giveIPValueSize(type, aGaussPoint) ) {
        answer.resize(1, 3);
    } else {
        return;
    }

    answer.at(1) = l1;
    answer.at(2) = l2;
    answer.at(3) = l3;
}

void
TR1_2D_SUPG2 :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                           InternalStateType type, TimeStep *tStep)
{
    GaussPoint *gp;
    if ( integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints() ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    } else {
        gp = integrationRulesArray [ 1 ]->getIntegrationPoint(0);
    }

    this->giveIPValue(answer, gp, type, tStep);
}

void
TR1_2D_SUPG2 :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                          InternalStateType type, TimeStep *tStep)
{
    answer.resize(0);
}

void
TR1_2D_SUPG2 :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(3);
    pap.at(1) = this->giveNode(1)->giveNumber();
    pap.at(2) = this->giveNode(2)->giveNumber();
    pap.at(3) = this->giveNode(3)->giveNumber();
}

void
TR1_2D_SUPG2 :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
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
TR1_2D_SUPG2 :: SPRNodalRecoveryMI_giveNumberOfIP()
{ return 1; }


void
TR1_2D_SUPG2 :: SPRNodalRecoveryMI_computeIPGlobalCoordinates(FloatArray &coords, GaussPoint *gp)
{
    this->computeGlobalCoordinates( coords, * gp->giveCoordinates() );
}

SPRPatchType
TR1_2D_SUPG2 :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_2dxy;
}


void
TR1_2D_SUPG2 :: printOutputAt(FILE *file, TimeStep *stepN)
// Performs end-of-step operations.
{
    SUPGElement :: printOutputAt(file, stepN);

    GaussPoint *gp;
    if ( integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints() ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    } else {
        gp = integrationRulesArray [ 1 ]->getIntegrationPoint(0);
    }

    double rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, gp, stepN);
    fprintf(file, "VOF %e, density %e\n\n", this->giveVolumeFraction(), rho);
}


contextIOResultType TR1_2D_SUPG2 :: saveContext(DataStream *stream, ContextMode mode, void *obj)
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



contextIOResultType TR1_2D_SUPG2 :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
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




#ifdef __OOFEG
int
TR1_2D_SUPG2 :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
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
TR1_2D_SUPG2 :: drawRawGeometry(oofegGraphicContext &gc)
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

void TR1_2D_SUPG2 :: drawScalar(oofegGraphicContext &context)
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

    /*
     * if ((context.giveIntVarMode() == ISM_local) && (context.giveIntVarType() ==  IST_VOFFraction)) {
     * Polygon matvolpoly;
     * FloatArray _n;
     * double _p;
     * this->doCellDLS (_n, this->giveNumber());
     * this->findCellLineConstant (_p, _n, this->giveNumber());
     * this->formMaterialVolumePoly(matvolpoly, NULL, _n, _p, false);
     * GraphicObj *go = matvolpoly.draw(context,true,OOFEG_VARPLOT_PATTERN_LAYER);
     * return;
     * }
     */

    if ( context.giveIntVarMode() == ISM_recovered ) {
        result += this->giveInternalStateAtNode(v1, context.giveIntVarType(), context.giveIntVarMode(), 1, tStep);
        result += this->giveInternalStateAtNode(v2, context.giveIntVarType(), context.giveIntVarMode(), 2, tStep);
        result += this->giveInternalStateAtNode(v3, context.giveIntVarType(), context.giveIntVarMode(), 3, tStep);
    } else if ( context.giveIntVarMode() == ISM_local ) {
        GaussPoint *gp;
        if ( integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints() ) {
            gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
        } else {
            gp = integrationRulesArray [ 1 ]->getIntegrationPoint(0);
        }

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
