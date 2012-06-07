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

#include "tr1_2d_supg_axi.h"
#include "node.h"
#include "material.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "domain.h"
#include "mathfem.h"
#include "engngm.h"
#include "timestep.h"
#include "load.h"
#include "boundaryload.h"

#ifndef __MAKEDEPEND
 #include <math.h>
 #include <stdio.h>
#endif

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
#endif

namespace oofem {
TR1_2D_SUPG_AXI :: TR1_2D_SUPG_AXI(int n, Domain *aDomain) : TR1_2D_SUPG(n, aDomain)
    // Constructor.
{ }

TR1_2D_SUPG_AXI :: ~TR1_2D_SUPG_AXI()
// Destructor
{ }

void
TR1_2D_SUPG_AXI :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 1;
        integrationRulesArray = new IntegrationRule * [ 1 ];
        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Triangle, 7, _2dAxiFlow);
    }
}


void
TR1_2D_SUPG_AXI :: computeAccelerationTerm_MB(FloatMatrix &answer, TimeStep *atTime)
{
    answer.resize(6, 6);
    answer.zero();
    FloatArray un;
    double rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), atTime);
    FloatArray n(3);

    int i, j, ip;
    double _val, u, v, dV;
    GaussPoint *gp;
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime->givePreviousStep(), un);

    for ( ip = 0; ip < integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints(); ip++ ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(ip);
        dV = this->computeVolumeAround(gp);
        computeNVector(n, gp);

        for ( i = 1; i <= 3; i++ ) {
            for ( j = 1; j <= 3; j++ ) {
                /* consistent mass */
                _val = n.at(i) * n.at(j) * rho * dV;
                answer.at(2 * i - 1, 2 * j - 1) += _val;
                answer.at(2 * i, 2 * j)     += _val;

                /* SUPG stabilization term */
                u = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
                v = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);
                answer.at(2 * i - 1, 2 * j - 1) += ( u * b [ i - 1 ] + v * c [ j - 1 ] ) * n.at(j) * rho * t_supg * dV;
                answer.at(2 * i, 2 * j)     += ( u * b [ i - 1 ] + v * c [ j - 1 ] ) * n.at(j) * rho * t_supg * dV;
            }
        }
    }
}

void
TR1_2D_SUPG_AXI :: computeAdvectionTerm_MB(FloatArray &answer, TimeStep *atTime)
{
    answer.resize(6);
    answer.zero();

    double dV;
    double dudx, dudy, dvdx, dvdy, _u, _v;
    GaussPoint *gp;
    FloatArray u, un, n(3);
    int ip, i;
    double rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), atTime);

    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime->givePreviousStep(), un);
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);

    dudx = b [ 0 ] * u.at(1) + b [ 1 ] * u.at(3) + b [ 2 ] * u.at(5);
    dudy = c [ 0 ] * u.at(1) + c [ 1 ] * u.at(3) + c [ 2 ] * u.at(5);
    dvdx = b [ 0 ] * u.at(2) + b [ 1 ] * u.at(4) + b [ 2 ] * u.at(6);
    dvdy = c [ 0 ] * u.at(2) + c [ 1 ] * u.at(4) + c [ 2 ] * u.at(6);

    for ( ip = 0; ip < integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints(); ip++ ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(ip);
        dV = this->computeVolumeAround(gp);
        computeNVector(n, gp);

        _u = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
        _v = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);

        // standard galerkin term
        for ( i = 1; i <= 3; i++ ) {
            answer.at(2 * i - 1) += rho * n.at(i) * ( _u * dudx + _v * dudy ) * dV;
            answer.at(2 * i)   += rho * n.at(i) * ( _u * dvdx + _v * dvdy ) * dV;
        }

        // supg stabilization term
        for ( i = 1; i <= 3; i++ ) {
            answer.at(2 * i - 1) += ( _u * b [ i - 1 ] + _v * c [ i - 1 ] ) * ( _u * dudx + _v * dudy ) * rho * t_supg * dV;
            answer.at(2 * i)   += ( _u * b [ i - 1 ] + _v * c [ i - 1 ] ) * ( _u * dvdx + _v * dvdy ) * rho * t_supg * dV;
        }
    }
}

void
TR1_2D_SUPG_AXI :: computeAdvectionDerivativeTerm_MB(FloatMatrix &answer, TimeStep *atTime) {
    answer.resize(6, 6);
    answer.zero();

    FloatArray u, un, n(3);
    double rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), atTime);
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime->givePreviousStep(), un);
    double dudx [ 2 ] [ 2 ];
    int i, j, k, m, w_dof_addr, u_dof_addr, d1j, d2j, dij;

    int ip;
    double _u, _v, dV;
    GaussPoint *gp;

    dudx [ 0 ] [ 0 ] = b [ 0 ] * u.at(1) + b [ 1 ] * u.at(3) + b [ 2 ] * u.at(5);
    dudx [ 0 ] [ 1 ] = c [ 0 ] * u.at(1) + c [ 1 ] * u.at(3) + c [ 2 ] * u.at(5);
    dudx [ 1 ] [ 0 ] = b [ 0 ] * u.at(2) + b [ 1 ] * u.at(4) + b [ 2 ] * u.at(6);
    dudx [ 1 ] [ 1 ] = c [ 0 ] * u.at(2) + c [ 1 ] * u.at(4) + c [ 2 ] * u.at(6);

    for ( ip = 0; ip < integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints(); ip++ ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(ip);
        dV = this->computeVolumeAround(gp);
        computeNVector(n, gp);

        _u = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
        _v = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);

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
                        answer.at(w_dof_addr, u_dof_addr) += dV * rho * n.at(k) * ( 0.0 * d1j * n.at(m) * dudx [ i - 1 ] [ 0 ] + dij * _u * b [ m - 1 ] +
                                                                                    0.0 * d2j * n.at(m) * dudx [ i - 1 ] [ 0 ] + dij * _v * c [ m - 1 ] );
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
                        answer.at(w_dof_addr, u_dof_addr) += dV * t_supg * rho *
                                                             ( _u * b [ k - 1 ] + _v * c [ k - 1 ] ) * ( dij * _u * b [ m - 1 ] + dij * _v * c [ m - 1 ] );
                    }
                }
            }
        }
    }
}



void
TR1_2D_SUPG_AXI :: computeDiffusionTerm_MB(FloatArray &answer, TimeStep *atTime)
{
    answer.resize(6);
    answer.zero();
    FloatArray u, eps, stress;
    double Re = domain->giveEngngModel()->giveUnknownComponent(ReynoldsNumber, VM_Unknown, atTime, domain, NULL);
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);

    FloatMatrix _b(4, 6);

    int ip;
    double dV;
    GaussPoint *gp;
    FloatArray bs;

    // stabilization term K_delta
    FloatArray un(6), n(3);
    double _u, _v, _r;
    int i;
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime->givePreviousStep(), un);
    // end k_delta declaration

    for ( ip = 0; ip < integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints(); ip++ ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(ip);
        dV = this->computeVolumeAround(gp);
        this->computeBMtrx(_b, gp);
        eps.beProductOf(_b, u);
        ( ( FluidDynamicMaterial * ) this->giveMaterial() )->computeDeviatoricStressVector(stress, gp,
                                                                                           eps, atTime);
        stress.times(dV / Re);
        bs.beTProductOf(_b, stress);
        answer.add(bs);

#if 1
        // stabilization term k_delta
        _r = this->computeRadiusAt(gp);
        computeNVector(n, gp);
        _u = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
        _v = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);

        for ( i = 1; i <= 3; i++ ) {
            answer.at(2 * i - 1) -= t_supg * ( _u * b [ i - 1 ] + _v * c [ i - 1 ] ) * ( stress.at(1) / _r );
            answer.at(2 * i)   -= t_supg * ( _u * b [ i - 1 ] + _v * c [ i - 1 ] ) * ( stress.at(4) / _r );
        }

#endif
    }
}


void
TR1_2D_SUPG_AXI :: computeDiffusionDerivativeTerm_MB(FloatMatrix &answer, MatResponseMode mode, TimeStep *atTime)
{
    answer.resize(6, 6);
    answer.zero();
    FloatMatrix _db, _d, _b(4, 6);
    double Re = domain->giveEngngModel()->giveUnknownComponent(ReynoldsNumber, VM_Unknown, atTime, domain, NULL);

    int ip;
    double dV;
    GaussPoint *gp;

    // stabilization term K_delta
    FloatArray un(6), u(6), n(3), eps, stress;
    double _u, _v, _r;
    int i, j;
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime->givePreviousStep(), un);
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);


    for ( ip = 0; ip < integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints(); ip++ ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(ip);
        dV = this->computeVolumeAround(gp);
        this->computeBMtrx(_b, gp);
        ( ( FluidDynamicMaterial * ) this->giveMaterial() )->giveDeviatoricStiffnessMatrix(_d, mode,
                                                                                           gp, atTime);
        _db.beProductOf(_d, _b);
        answer.plusProductUnsym(_b, _db, dV);
        //answer.plusProductSymmUpper (_bs,_db,dV*t_supg);
        // }

#if 1
        _r = this->computeRadiusAt(gp);
        computeNVector(n, gp);
        eps.beProductOf(_b, u);
        ( ( FluidDynamicMaterial * ) this->giveMaterial() )->computeDeviatoricStressVector(stress, gp,
                                                                                           eps, atTime);
        //_mu = this->giveMaterial()->giveCharacteristicValue(MRM_Viscosity, gp, atTime);

        _u = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
        _v = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);

        for ( i = 1; i <= 3; i++ ) {
            for ( j = 1; j <= 6; j++ ) {
                //answer.at(2*i-1,j) -= t_supg*(_u*b[i-1]+_v*c[i-1])*(_d.at(1,1)*_b.at(1,j)/_r - stress.at(1)/_r/_r)*dV;
                //answer.at(2*i,j)   -= t_supg*(_u*b[i-1]+_v*c[i-1])*(_d.at(4,4)*_b.at(4,j)/_r - stress.at(4)/_r/_r)*dV;

                answer.at(2 * i - 1, j) -= t_supg * ( _u * b [ i - 1 ] + _v * c [ i - 1 ] ) * ( _db.at(1, j) / _r ) * dV;
                answer.at(2 * i, j)   -= t_supg * ( _u * b [ i - 1 ] + _v * c [ i - 1 ] ) * ( _db.at(4, j) / _r ) * dV;
            }
        }

#endif
    }



    answer.times(1. / Re);
}


void
TR1_2D_SUPG_AXI :: computePressureTerm_MB(FloatMatrix &answer, TimeStep *atTime)
{
    answer.resize(6, 3);
    answer.zero();
    FloatArray un, n(3);
    double _u, _v;

    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime->givePreviousStep(), un);
    int i, j, ip;
    double dV, _r;
    GaussPoint *gp;
    for ( ip = 0; ip < integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints(); ip++ ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(ip);
        dV = this->computeVolumeAround(gp);
        _r = this->computeRadiusAt(gp);
        computeNVector(n, gp);
        // G matrix
        for ( i = 1; i <= 3; i++ ) {
            for ( j = 1; j <= 3; j++ ) {
                answer.at(2 * i - 1, j) -= b [ i - 1 ] * n.at(j) * dV;
                answer.at(2 * i, j)   -= c [ i - 1 ] * n.at(j) * dV;

                answer.at(2 * i - 1, j) -= n.at(i) * n.at(j) * dV / _r;
            }
        }

        // stabilization term (G_\delta mtrx)
        _u = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
        _v = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);
        for ( i = 1; i <= 3; i++ ) {
            for ( j = 1; j <= 3; j++ ) {
                answer.at(2 * i - 1, j) += ( _u * b [ i - 1 ] + _v * c [ i - 1 ] ) * b [ j - 1 ] * dV * t_supg;
                answer.at(2 * i, j)   += ( _u * b [ i - 1 ] + _v * c [ i - 1 ] ) * c [ j - 1 ] * dV * t_supg;
            }
        }
    }
}

void
TR1_2D_SUPG_AXI :: computeLSICStabilizationTerm_MB(FloatMatrix &answer, TimeStep *atTime)
{
    answer.resize(6, 6);
    answer.zero();
    double rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), atTime);
    double n[] = {
        b [ 0 ], c [ 0 ], b [ 1 ], c [ 1 ], b [ 2 ], c [ 2 ]
    };
    int i, j;

    int ip;
    double dV;
    GaussPoint *gp;
    for ( ip = 0; ip < integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints(); ip++ ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(ip);
        dV = this->computeVolumeAround(gp);

        for ( i = 1; i <= 6; i++ ) {
            for ( j = 1; j <= 6; j++ ) {
                answer.at(i, j) += dV * t_lsic * rho * n [ i - 1 ] * n [ j - 1 ];
            }
        }
    }
}



void
TR1_2D_SUPG_AXI :: computeLinearAdvectionTerm_MC(FloatMatrix &answer, TimeStep *atTime)
{
    answer.resize(3, 6);
    answer.zero();

    int ip, i, j;
    double dV, _r;
    GaussPoint *gp;
    FloatArray n(3);
    for ( ip = 0; ip < integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints(); ip++ ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(ip);
        dV = this->computeVolumeAround(gp);
        _r = this->computeRadiusAt(gp);
        n.at(1) = gp->giveCoordinate(1);
        n.at(2) = gp->giveCoordinate(2);
        n.at(3) = 1. - n.at(1) - n.at(2);
        for ( i = 1; i <= 3; i++ ) {
            for ( j = 1; j <= 3; j++ ) {
                answer.at(j, 2 * i - 1) += b [ i - 1 ] * n.at(j) * dV;
                answer.at(j, 2 * i)   += c [ i - 1 ] * n.at(j) * dV;

                answer.at(i, 1 + ( j - 1 ) * 2) += n.at(i) * n.at(j) * dV / _r;
            }
        }
    }
}

void
TR1_2D_SUPG_AXI :: computeAdvectionTerm_MC(FloatArray &answer, TimeStep *atTime)
{
    // N_epsilon (due to PSPG stabilization)
    double dudx, dudy, dvdx, dvdy, _u, _v;
    FloatArray u, un, n(3);
    int ip, i;
    double dV;
    GaussPoint *gp;

    answer.resize(3);
    answer.zero();

    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime->givePreviousStep(), un);
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);
    dudx = b [ 0 ] * u.at(1) + b [ 1 ] * u.at(3) + b [ 2 ] * u.at(5);
    dudy = c [ 0 ] * u.at(1) + c [ 1 ] * u.at(3) + c [ 2 ] * u.at(5);
    dvdx = b [ 0 ] * u.at(2) + b [ 1 ] * u.at(4) + b [ 2 ] * u.at(6);
    dvdy = c [ 0 ] * u.at(2) + c [ 1 ] * u.at(4) + c [ 2 ] * u.at(6);

    for ( ip = 0; ip < integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints(); ip++ ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(ip);
        dV = this->computeVolumeAround(gp);
        computeNVector(n, gp);

        _u = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
        _v = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);

        for ( i = 1; i <= 3; i++ ) {
            answer.at(i) += t_pspg * dV * ( b [ i - 1 ] * ( _u * dudx + _v * dudy ) + c [ i - 1 ] * ( _u * dvdx + _v * dvdy ) );
        }
    }
}


void
TR1_2D_SUPG_AXI :: computeAdvectionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *atTime)
{
    answer.resize(3, 6);
    answer.zero();
    int k, j, m, w_dof_addr, u_dof_addr, d1j, d2j, km1, mm1;
    //double rho = this->giveMaterial()->giveCharacteristicValue(Density, integrationRulesArray[0]->getIntegrationPoint(0), atTime);
    FloatArray u, un, n(3);
    int ip;
    double dV, _u, _v;
    GaussPoint *gp;

    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime->givePreviousStep(), un);
    double dudx [ 2 ] [ 2 ];

    dudx [ 0 ] [ 0 ] = b [ 0 ] * u.at(1) + b [ 1 ] * u.at(3) + b [ 2 ] * u.at(5);
    dudx [ 0 ] [ 1 ] = c [ 0 ] * u.at(1) + c [ 1 ] * u.at(3) + c [ 2 ] * u.at(5);
    dudx [ 1 ] [ 0 ] = b [ 0 ] * u.at(2) + b [ 1 ] * u.at(4) + b [ 2 ] * u.at(6);
    dudx [ 1 ] [ 1 ] = c [ 0 ] * u.at(2) + c [ 1 ] * u.at(4) + c [ 2 ] * u.at(6);

    for ( ip = 0; ip < integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints(); ip++ ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(ip);
        dV = this->computeVolumeAround(gp);
        computeNVector(n, gp);

        _u = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
        _v = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);

        for ( k = 1; k <= 3; k++ ) { // nodal val of function w
            km1 = k - 1;
            for ( j = 1; j <= 2; j++ ) { // velocity vector component
                for ( m = 1; m <= 3; m++ ) { //  nodal components
                    w_dof_addr = k;
                    u_dof_addr = ( m - 1 ) * 2 + j;
                    mm1 = m - 1;
                    d1j = ( j == 1 );
                    d2j = ( j == 2 );
                    answer.at(w_dof_addr, u_dof_addr) += t_pspg * dV * ( b [ km1 ] * ( _u * d1j * b [ mm1 ] + _v * d1j * c [ mm1 ] ) + c [ km1 ] * ( _u * d2j * b [ mm1 ] + _v * d2j * c [ mm1 ] ) );
                }
            }
        }
    }
}

void TR1_2D_SUPG_AXI :: computeDiffusionTerm_MC(FloatArray &answer, TimeStep *atTime)
{
    answer.resize(3);
    answer.zero();

#if 1
    double Re = domain->giveEngngModel()->giveUnknownComponent(ReynoldsNumber, VM_Unknown, atTime, domain, NULL);
    double rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), atTime);
    int ip, i;
    double dV, _r;
    GaussPoint *gp;
    FloatArray eps, stress, u(6);
    FloatMatrix _b;

    // stabilization term K_eps
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);

    for ( ip = 0; ip < integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints(); ip++ ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(ip);
        dV = this->computeVolumeAround(gp);
        _r = this->computeRadiusAt(gp);
        this->computeBMtrx(_b, gp);
        eps.beProductOf(_b, u);
        ( ( FluidDynamicMaterial * ) this->giveMaterial() )->computeDeviatoricStressVector(stress, gp,
                                                                                           eps, atTime);
        stress.times(1. / Re);
        for ( i = 1; i <= 3; i++ ) {
            answer.at(i) -= t_pspg * ( b [ i - 1 ] * stress.at(1) + c [ i - 1 ] * stress.at(4) ) * dV / rho / _r;
        }
    }

#endif
}

void TR1_2D_SUPG_AXI :: computeDiffusionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *atTime)
{
    answer.resize(3, 6);
    answer.zero();

#if 1
    double Re = domain->giveEngngModel()->giveUnknownComponent(ReynoldsNumber, VM_Unknown, atTime, domain, NULL);
    double rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), atTime);
    int ip, i, j;
    double dV, _r, _mu;
    GaussPoint *gp;
    FloatMatrix _d, _b, _db;
    //FloatArray eps, stress,u;

    // stabilization term K_eps
    //this -> computeVectorOf(EID_MomentumBalance,VM_Total,atTime, u) ;
    for ( ip = 0; ip < integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints(); ip++ ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(ip);
        dV = this->computeVolumeAround(gp);
        _r = this->computeRadiusAt(gp);
        _mu = this->giveMaterial()->giveCharacteristicValue(MRM_Viscosity, gp, atTime);
        this->computeBMtrx(_b, gp);
        ( ( FluidDynamicMaterial * ) this->giveMaterial() )->giveDeviatoricStiffnessMatrix(_d, TangentStiffness, gp, atTime);
        _db.beProductOf(_d, _b);
        //eps.beProductOf (_b, u);
        //((FluidDynamicMaterial*) this->giveMaterial())->computeDeviatoricStressVector (stress,gp,eps,atTime);

        for ( i = 1; i <= 3; i++ ) {
            for ( j = 1; j <= 6; j++ ) {
                answer.at(i, j) -= t_pspg * ( b [ i - 1 ] * ( _db.at(1, j) / _r ) + c [ i - 1 ] * ( _db.at(4, j) / _r ) ) * dV / rho;
            }
        }
    }

    answer.times(1. / Re);
#endif
}


void
TR1_2D_SUPG_AXI :: computeAccelerationTerm_MC(FloatMatrix &answer, TimeStep *atTime)
{
    answer.resize(3, 6);
    answer.zero();
    int ip, i, j;
    double dV;
    GaussPoint *gp;
    FloatArray n(3);
    // M_\epsilon

    for ( ip = 0; ip < integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints(); ip++ ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(ip);
        dV = this->computeVolumeAround(gp);
        computeNVector(n, gp);

        for ( i = 1; i <= 3; i++ ) {
            for ( j = 1; j <= 3; j++ ) {
                answer.at(i, 2 * j - 1) += t_pspg * dV * b [ i - 1 ] * n.at(j);
                answer.at(i, 2 * j)   += t_pspg * dV * c [ i - 1 ] * n.at(j);
            }
        }
    }
}


void
TR1_2D_SUPG_AXI :: computePressureTerm_MC(FloatMatrix &answer, TimeStep *atTime)
{
    double rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), atTime);
    int ip, i, j;
    double coeff = t_pspg / rho;
    double dV;
    GaussPoint *gp;

    answer.resize(3, 3);
    answer.zero();

    for ( ip = 0; ip < integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints(); ip++ ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(ip);
        dV = this->computeVolumeAround(gp);

        for ( i = 1; i <= 3; i++ ) {
            for ( j = 1; j <= 3; j++ ) {
                answer.at(i, j) += dV * coeff * ( b [ i - 1 ] * b [ j - 1 ] + c [ i - 1 ] * c [ j - 1 ] );
            }
        }
    }
}

void
TR1_2D_SUPG_AXI :: computeDeviatoricStrain(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    /* one should call material driver instead */
    FloatArray u(6);
    FloatMatrix _b(4, 6);
    answer.resize(4);


    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);
    this->computeBMtrx(_b, gp);
    answer.beProductOf(_b, u);
}

void
TR1_2D_SUPG_AXI :: computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    /* one should call material driver instead */
    FloatArray eps(4);
    answer.resize(3);

    this->computeDeviatoricStrain(eps, gp, tStep);
    ( ( FluidDynamicMaterial * ) this->giveMaterial() )->computeDeviatoricStressVector(answer, gp, eps, tStep);
}



void
TR1_2D_SUPG_AXI :: computeBCRhsTerm_MB(FloatArray &answer, TimeStep *atTime)
{
    answer.resize(6);
    answer.zero();

    int i, ip, nLoads;
    Load *load;
    bcGeomType ltype;
    double rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), atTime);
    FloatArray un, gVector;
    GaussPoint *gp;
    FloatArray n(3);
    double dV, _r, coeff, u, v;

    // add body load (gravity) termms
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime->givePreviousStep(), un);
    nLoads    = this->giveBodyLoadArray()->giveSize();
    for ( i = 1; i <= nLoads; i++ ) {
        load  = domain->giveLoad( bodyLoadArray.at(i) );
        ltype = load->giveBCGeoType();
        if ( ( ltype == BodyLoadBGT ) && ( load->giveBCValType() == ForceLoadBVT ) ) {
            load->computeComponentArrayAt(gVector, atTime, VM_Total);
            if ( gVector.giveSize() ) {
                for ( ip = 0; ip < integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints(); ip++ ) {
                    gp = integrationRulesArray [ 0 ]->getIntegrationPoint(ip);
                    dV = this->computeVolumeAround(gp);
                    _r = this->computeRadiusAt(gp);
                    computeNVector(n, gp);
                    coeff = rho * dV;
                    u = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
                    v = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);

                    answer.at(1) += coeff * ( gVector.at(1) * ( n.at(1) + t_supg * ( b [ 0 ] * u + c [ 0 ] * v ) ) );
                    answer.at(2) += coeff * ( gVector.at(2) * ( n.at(1) + t_supg * ( b [ 0 ] * u + c [ 0 ] * v ) ) );
                    answer.at(3) += coeff * ( gVector.at(1) * ( n.at(2) + t_supg * ( b [ 1 ] * u + c [ 1 ] * v ) ) );
                    answer.at(4) += coeff * ( gVector.at(2) * ( n.at(2) + t_supg * ( b [ 1 ] * u + c [ 1 ] * v ) ) );
                    answer.at(5) += coeff * ( gVector.at(1) * ( n.at(3) + t_supg * ( b [ 2 ] * u + c [ 2 ] * v ) ) );
                    answer.at(6) += coeff * ( gVector.at(2) * ( n.at(3) + t_supg * ( b [ 2 ] * u + c [ 2 ] * v ) ) );
                }
            }
        }
    }

    //answer.times(rc);

    // loop over sides
    int j, n1, n2, code, sid;
    double tx, ty, l, nx, ny, side_r;
    //IntArray nodecounter (3);
    for ( j = 1; j <= boundarySides.giveSize(); j++ ) {
        code = boundaryCodes.at(j);
        sid = boundarySides.at(j);
        if ( ( code & FMElement_PrescribedTractionBC ) ) {
            FloatArray t, coords(1);
            int nLoads, n, id;
            BoundaryLoad *load;
            // integrate tractions
            n1 = sid;
            n2 = ( n1 == 3 ? 1 : n1 + 1 );

            tx = giveNode(n2)->giveCoordinate(1) - giveNode(n1)->giveCoordinate(1);
            ty = giveNode(n2)->giveCoordinate(2) - giveNode(n1)->giveCoordinate(2);
            l = sqrt(tx * tx + ty * ty);
            nx = ty / l;
            ny = -tx / l;
            // radius at side center
            side_r = 0.5 * ( giveNode(n2)->giveCoordinate(1) + giveNode(n1)->giveCoordinate(1) );

            // if no traction bc applied but side marked as with traction load
            // then zero traction is assumed !!!

            // loop over boundary load array
            nLoads    = this->giveBoundaryLoadArray()->giveSize() / 2;
            for ( i = 1; i <= nLoads; i++ ) {
                n     = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
                id    = boundaryLoadArray.at(i * 2);
                if ( id != sid ) {
                    continue;
                }

                load  = dynamic_cast< BoundaryLoad * >( domain->giveLoad(n) );
                if ( load ) {
                    load->computeValueAt(t, atTime, coords, VM_Total);

                    // here it is assumed constant traction, one point integration only
                    // n1 (u,v)
                    answer.at( ( n1 - 1 ) * 2 + 1 ) += t.at(1) * side_r * l / 2.;
                    answer.at(n1 * 2)       += t.at(2) * side_r * l / 2.;
                    // n2 (u,v)
                    answer.at( ( n2 - 1 ) * 2 + 1 ) += t.at(1) * side_r * l / 2.;
                    answer.at(n2 * 2)       += t.at(2) * side_r * l / 2.;

                    //answer.at(n1)+= (t.at(1)*nx + t.at(2)*ny) * side_r * l/2.;
                    //answer.at(n2)+= (t.at(1)*nx + t.at(2)*ny) * side_r * l/2.;
                }
            }
        }
    }
}

void
TR1_2D_SUPG_AXI :: computeBCRhsTerm_MC(FloatArray &answer, TimeStep *atTime)
{
    int i, ip, nLoads;
    Load *load;
    bcGeomType ltype;
    FloatArray gVector;
    GaussPoint *gp;
    double coeff, _r, dV;

    answer.resize(3);
    answer.zero();
    nLoads    = this->giveBodyLoadArray()->giveSize();
    for ( i = 1; i <= nLoads; i++ ) {
        load  = domain->giveLoad( bodyLoadArray.at(i) );
        ltype = load->giveBCGeoType();
        if ( ( ltype == BodyLoadBGT ) && ( load->giveBCValType() == ForceLoadBVT ) ) {
            load->computeComponentArrayAt(gVector, atTime, VM_Total);
            if ( gVector.giveSize() ) {
                for ( ip = 0; ip < integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints(); ip++ ) {
                    gp = integrationRulesArray [ 0 ]->getIntegrationPoint(ip);
                    dV = this->computeVolumeAround(gp);
                    _r = this->computeRadiusAt(gp);
                    coeff = t_pspg * dV;

                    answer.at(1) += coeff * ( b [ 0 ] * gVector.at(1) + c [ 0 ] * gVector.at(2) );
                    answer.at(2) += coeff * ( b [ 1 ] * gVector.at(1) + c [ 1 ] * gVector.at(2) );
                    answer.at(3) += coeff * ( b [ 2 ] * gVector.at(1) + c [ 2 ] * gVector.at(2) );
                }
            }
        }
    }
}


void
TR1_2D_SUPG_AXI :: computeBCLhsTerm_MB(FloatMatrix &answer, TimeStep *atTime)
{
    bcType boundarytype;
    int i, n, side;
    int nLoads = 0;
    Load *load;
    //bcType loadtype;
    FloatMatrix helpMatrix;
    // loop over boundary load array
    helpMatrix.resize(6, 6);
    helpMatrix.zero();

    answer.resize(6, 6);
    answer.zero();

    nLoads    = this->giveBoundaryLoadArray()->giveSize() / 2;

    if ( nLoads ) {
        for ( i = 1; i <= nLoads; i++ ) {
            n     = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
            side    = boundaryLoadArray.at(i * 2);
            load  = dynamic_cast< Load * >( domain->giveLoad(n) );
            boundarytype = load->giveType();
            if ( boundarytype == SlipWithFriction ) {
                this->computeSlipWithFrictionBCTerm_MB(helpMatrix, ( Load * ) load, side, atTime);
            } else if ( boundarytype == PenetrationWithResistance ) {
                this->computePenetrationWithResistanceBCTerm_MB(helpMatrix, ( Load * ) load, side, atTime);
            } else {
                answer.resize(6, 6);
                answer.zero();
                // _error("computeForceLoadVector : unsupported load type class");
            }

            answer.add(helpMatrix);
        }
    }
}


void
TR1_2D_SUPG_AXI :: computeSlipWithFrictionBCTerm_MB(FloatMatrix &answer, Load *load, int side, TimeStep *atTime)
{
    //answer.resize(6, 6);
    //answer.zero();

    int node1, node2, node3, i, j, ip;
    double l, t1, t2, _t1, _t2, dV;
    double beta;
    GaussPoint *gp;
    FloatArray nt(6), n(3);
    answer.resize(6, 6);
    answer.zero();

    BoundaryLoad *edgeLoad = static_cast< BoundaryLoad * >(load);
    beta = edgeLoad->giveProperty('a');
    node1 = side;
    node2 = ( node1 == 3 ? 1 : node1 + 1 );

    node3 = ( node2 == 3 ? 1 : node2 + 1 );

    _t1 = giveNode(node2)->giveCoordinate(1) - giveNode(node1)->giveCoordinate(1);
    _t2 = giveNode(node2)->giveCoordinate(2) - giveNode(node1)->giveCoordinate(2);
    l = sqrt(_t1 * _t1 + _t2 * _t2);

    t1 = _t1 / l;
    t2 = _t2 / l;

    for ( ip = 0; ip < integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints(); ip++ ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(ip);
        dV = this->computeVolumeAround(gp);
        computeNVector(n, gp);

        nt.at(1) = n.at(1) * t1;
        nt.at(2) = n.at(1) * t2;
        nt.at(3) = n.at(2) * t1;
        nt.at(4) = n.at(2) * t2;
        nt.at(5) = n.at(3) * t1;
        nt.at(6) = n.at(3) * t2;

        for ( i = 1; i <= 6; i++ ) {
            for ( j = 1; j <= 6; j++ ) {
                answer.at(i, j) += nt.at(i) * nt.at(j);
            }
        }

        answer.times(beta * dV);
    }
}

void
TR1_2D_SUPG_AXI :: computePenetrationWithResistanceBCTerm_MB(FloatMatrix &answer, Load *load, int side, TimeStep *atTime)
{
    //answer.resize(6, 6);
    //answer.zero();

    int node1, node2, node3, i, j, ip;
    double l, n1, n2, _t1, _t2, dV;
    double alpha;
    FloatArray nt(6), n(3);
    GaussPoint *gp;
    answer.resize(6, 6);
    answer.zero();

    BoundaryLoad *edgeLoad = static_cast< BoundaryLoad * >(load);
    alpha = edgeLoad->giveProperty('a');
    node1 = side;
    node2 = ( node1 == 3 ? 1 : node1 + 1 );

    node3 = ( node2 == 3 ? 1 : node2 + 1 );


    _t1 = giveNode(node2)->giveCoordinate(1) - giveNode(node1)->giveCoordinate(1);
    _t2 = giveNode(node2)->giveCoordinate(2) - giveNode(node1)->giveCoordinate(2);
    l = sqrt(_t1 * _t1 + _t2 * _t2);

    //t1 = _t1 / l;
    //t2 = _t2 / l;

    n1 = _t2 / l;
    n2 = -_t1 / l;


    for ( ip = 0; ip < integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints(); ip++ ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(ip);
        dV = this->computeVolumeAround(gp);
        this->computeNVector(n, gp);

        nt.at(1) = n.at(1) * n1;
        nt.at(2) = n.at(1) * n2;
        nt.at(3) = n.at(2) * n1;
        nt.at(4) = n.at(2) * n2;
        nt.at(5) = n.at(3) * n1;
        nt.at(6) = n.at(3) * n2;

        for ( i = 1; i <= 6; i++ ) {
            for ( j = 1; j <= 6; j++ ) {
                answer.at(i, j) += nt.at(i) * nt.at(j);
            }
        }

        answer.times( ( 1 / alpha ) * dV );
    }
}

void
TR1_2D_SUPG_AXI :: computeOutFlowBCTerm_MB(FloatMatrix &answer, int side, TimeStep *atTime)
{
    int node1, node2, node3, j, i, ip;
    double l, n1, n2, t1, t2, dV;
    FloatArray nt(6), n(3);
    GaussPoint *gp;
    answer.resize(6, 3);
    answer.zero();
    //beta
    //area
    node1 = side;
    node2 = ( node1 == 3 ? 1 : node1 + 1 );

    node3 = ( node2 == 3 ? 1 : node2 + 1 );


    t1 = giveNode(node2)->giveCoordinate(1) - giveNode(node1)->giveCoordinate(1);
    t2 = giveNode(node2)->giveCoordinate(2) - giveNode(node1)->giveCoordinate(2);
    l = sqrt(t1 * t1 + t2 * t2);

    n1 = t2 / l;
    n2 = -t1 / l;


    for ( ip = 0; ip < integrationRulesArray [ 0 ]->getNumberOfIntegrationPoints(); ip++ ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(ip);
        dV = this->computeVolumeAround(gp);
        computeNVector(n, gp);

        nt.at(1) = n.at(1) * n1;
        nt.at(2) = n.at(1) * n2;
        nt.at(3) = n.at(2) * n1;
        nt.at(4) = n.at(2) * n2;
        nt.at(5) = n.at(3) * n1;
        nt.at(6) = n.at(3) * n2;

        for ( i = 1; i <= 6; i++ ) {
            for ( j = 1; j <= 3; j++ ) {
                answer.at(i, j) += nt.at(i) * n.at(j);
            }
        }

        answer.times(dV);
    }

    answer.times(-1.0);
}



void
TR1_2D_SUPG_AXI :: computeBCLhsPressureTerm_MB(FloatMatrix &answer, TimeStep *atTime)
{
    bcType boundarytype;
    int i, n, side;
    int nLoads = 0;
    Load *load;
    //bcType loadtype;
    FloatMatrix helpMatrix;
    // loop over boundary load array
    helpMatrix.resize(6, 3);
    helpMatrix.zero();

    answer.resize(6, 3);
    answer.zero();

    nLoads    = this->giveBoundaryLoadArray()->giveSize() / 2;

    if ( nLoads ) {
        for ( i = 1; i <= nLoads; i++ ) {
            n     = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
            side    = boundaryLoadArray.at(i * 2);
            load  = dynamic_cast< Load * >( domain->giveLoad(n) );
            boundarytype = load->giveType();
            if ( boundarytype == OutFlowBC ) {
                this->computeOutFlowBCTerm_MB(helpMatrix, side, atTime);
            } else {
                answer.resize(6, 3);
                answer.zero();
                // _error("computeForceLoadVector : unsupported load type class");
            }

            answer.add(helpMatrix);
        }
    }
}


double
TR1_2D_SUPG_AXI :: computeRadiusAt(GaussPoint *gp)
{
    double r1, r2, r3;
    double n1, n2, n3;

    r1 = this->giveNode(1)->giveCoordinate(1);
    r2 = this->giveNode(2)->giveCoordinate(1);
    r3 = this->giveNode(3)->giveCoordinate(1);

    n1 = gp->giveCoordinate(1);
    n2 = gp->giveCoordinate(2);
    n3 = 1. - n1 - n2;

    return n1 * r1 + n2 * r2 + n3 * r3;
}

void TR1_2D_SUPG_AXI :: computeBMtrx(FloatMatrix &_b, GaussPoint *gp)
{
    _b.resize(4, 6);
    double _r = this->computeRadiusAt(gp);


    _b.at(1, 1) = b [ 0 ];
    _b.at(1, 2) = 0.;
    _b.at(1, 3) = b [ 1 ];
    _b.at(1, 4) = 0.;
    _b.at(1, 5) = b [ 2 ];
    _b.at(1, 6) = 0.;
    _b.at(2, 1) = 1. / _r;
    _b.at(2, 2) = 0.;
    _b.at(2, 3) = 1. / _r;
    _b.at(2, 4) = 0.;
    _b.at(2, 5) = 1. / _r;
    _b.at(2, 6) = 0.;
    _b.at(3, 1) = 0.0;
    _b.at(3, 2) = c [ 0 ];
    _b.at(3, 3) = 0.0;
    _b.at(3, 4) = c [ 1 ];
    _b.at(3, 5) = 0.0;
    _b.at(3, 6) = c [ 2 ];
    _b.at(4, 1) = c [ 0 ];
    _b.at(4, 2) = b [ 0 ];
    _b.at(4, 3) = c [ 1 ];
    _b.at(4, 4) = b [ 1 ];
    _b.at(4, 5) = c [ 2 ];
    _b.at(4, 6) = b [ 2 ];

    /*
     *
     * double _c = (1./3.);
     * double u1 = _c*(b[0]+1./_r), u2 = _c*(b[1]+1./_r), u3 = _c*(b[2]+1./_r);
     * double w1 = _c*c[0], w2 = _c*c[1], w3 = _c*c[2];
     *
     * //u1=u2=u3=w1=w2=w3=0.0;
     *
     * _b.at(1,1) = b[0]-u1; _b.at(1,2)=0.-w1;   _b.at(1,3)= b[1]-u2; _b.at(1,4)=0.-w2;   _b.at(1,5)=b[2]-u3; _b.at(1,6)=0.-w3;
     * _b.at(2,1) = 1./_r-u1;_b.at(2,2)=0.-w1;   _b.at(2,3)= 1./_r-u2;_b.at(2,4)=0.-w2;   _b.at(2,5)=1./_r-u3;_b.at(2,6)=0.-w3;
     * _b.at(3,1) = 0.0-u1;  _b.at(3,2)=c[0]-w1; _b.at(3,3)= 0.0-u2;  _b.at(3,4)=c[1]-w2; _b.at(3,5)=0.0-u3;  _b.at(3,6)=c[2]-w3;
     * _b.at(4,1) = c[0];    _b.at(4,2)=b[0];    _b.at(4,3)= c[1];    _b.at(4,4)=b[1];    _b.at(4,5)=c[2];    _b.at(4,6)=b[2];
     */
}

void
TR1_2D_SUPG_AXI :: computeNVector(FloatArray &n, GaussPoint *gp)
{
    n.resize(3);

    n.at(1) = gp->giveCoordinate(1);
    n.at(2) = gp->giveCoordinate(2);
    n.at(3) = 1. - n.at(1) - n.at(2);
}

double
TR1_2D_SUPG_AXI :: computeVolumeAround(GaussPoint *aGaussPoint)
{
    double _r, weight;

    weight  = aGaussPoint->giveWeight();
    _r = computeRadiusAt(aGaussPoint);
    //_r = 1.0;

    return 2.0 * area * weight * _r;
}

void
TR1_2D_SUPG_AXI :: initGeometry()
{
    TR1_2D_SUPG :: initGeometry();

    this->rc = ( this->giveNode(1)->giveCoordinate(1) +
                this->giveNode(2)->giveCoordinate(1) +
                this->giveNode(3)->giveCoordinate(1) ) / 3.0;
    //this->rc = 1.0;
}



void
TR1_2D_SUPG_AXI :: updateStabilizationCoeffs(TimeStep *atTime)
{
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
        //this->t_lsic = h_ugn * vnorm * z / 2.0;
        this->t_lsic = 0.0;
    }

    this->t_supg *= uscale / lscale;
    this->t_pspg *= 1. / ( lscale * dscale );
    this->t_lsic *= ( dscale * uscale ) / ( lscale * lscale );

    //this->t_lsic = 0.0;
}

void
TR1_2D_SUPG_AXI :: LS_PCS_computeVOFFractions(FloatArray &answer, FloatArray &fi)
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


            // compute volume associated to traingle (x1,y1; x2,y2; x3,y3)
            double __volume = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 ) * ((x1+x2+x3)/3.);
            double volume = this->area*((this->giveNode(1)->giveCoordinate(1)+this->giveNode(2)->giveCoordinate(1)+this->giveNode(3)->giveCoordinate(1))/3.);
            if ( fabs(__volume) / volume > 1.00001 ) {
                OOFEM_ERROR("TR1_2D_SUPG::LS_PCS_computeVOFFractions: internal consistency error");
            }

            // prevent some roundoff errors
            if ( fabs(__volume) > volume ) {
                __volume = sgn(__volume) * volume;
            }

            if ( pos > neg ) {
                // negative area computed
                answer.at(2) = fabs(__volume) / volume;
                answer.at(1) = 1.0 - answer.at(2);
            } else {
                // postive area computed
                answer.at(1) = fabs(__volume) / volume;
                answer.at(2) = 1.0 - answer.at(1);
            }
        } else {
            OOFEM_ERROR("TR1_2D_SUPG::LS_PCS_computeVOFFractions: internal consistency error");
        }
    }
}


} // end namespace oofem
