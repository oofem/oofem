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

#include "tr1_2d_supg_axi.h"
#include "fei2dtrlin.h"
#include "fluidmodel.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "domain.h"
#include "mathfem.h"
#include "engngm.h"
#include "timestep.h"
#include "load.h"
#include "boundaryload.h"
#include "fluiddynamicmaterial.h"
#include "fluidcrosssection.h"
#include "crosssection.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif

namespace oofem {
REGISTER_Element(TR1_2D_SUPG_AXI);

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
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 1 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], 7, this);
    }
}


void
TR1_2D_SUPG_AXI :: computeAccelerationTerm_MB(FloatMatrix &answer, TimeStep *tStep)
{
    answer.resize(6, 6);
    answer.zero();
    FloatArray un, n;
    double _val, u, v, dV, rho;
    this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), un);

    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        dV = this->computeVolumeAround(gp);
        rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity(gp);
        this->computeNVector(n, gp);

        for ( int i = 1; i <= 3; i++ ) {
            for ( int j = 1; j <= 3; j++ ) {
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
TR1_2D_SUPG_AXI :: computeAdvectionTerm_MB(FloatArray &answer, TimeStep *tStep)
{
    answer.resize(6);
    answer.zero();

    double dudx, dudy, dvdx, dvdy, _u, _v;
    FloatArray u, un, n;

    this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), un);
    this->computeVectorOfVelocities(VM_Total, tStep, u);

    dudx = b [ 0 ] * u.at(1) + b [ 1 ] * u.at(3) + b [ 2 ] * u.at(5);
    dudy = c [ 0 ] * u.at(1) + c [ 1 ] * u.at(3) + c [ 2 ] * u.at(5);
    dvdx = b [ 0 ] * u.at(2) + b [ 1 ] * u.at(4) + b [ 2 ] * u.at(6);
    dvdy = c [ 0 ] * u.at(2) + c [ 1 ] * u.at(4) + c [ 2 ] * u.at(6);

    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        double dV = this->computeVolumeAround(gp);
        double rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity(gp);
        this->computeNVector(n, gp);

        _u = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
        _v = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);

        // standard galerkin term
        for ( int i = 1; i <= 3; i++ ) {
            answer.at(2 * i - 1) += rho * n.at(i) * ( _u * dudx + _v * dudy ) * dV;
            answer.at(2 * i)   += rho * n.at(i) * ( _u * dvdx + _v * dvdy ) * dV;
        }

        // supg stabilization term
        for ( int i = 1; i <= 3; i++ ) {
            answer.at(2 * i - 1) += ( _u * b [ i - 1 ] + _v * c [ i - 1 ] ) * ( _u * dudx + _v * dudy ) * rho * t_supg * dV;
            answer.at(2 * i)   += ( _u * b [ i - 1 ] + _v * c [ i - 1 ] ) * ( _u * dvdx + _v * dvdy ) * rho * t_supg * dV;
        }
    }
}

void
TR1_2D_SUPG_AXI :: computeAdvectionDerivativeTerm_MB(FloatMatrix &answer, TimeStep *tStep)
{
    answer.resize(6, 6);
    answer.zero();

    FloatArray u, un, n;
    this->computeVectorOfVelocities(VM_Total, tStep, u);
    this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), un);
    double dudx [ 2 ] [ 2 ];
    int w_dof_addr, u_dof_addr, d1j, d2j, dij;
    double _u, _v, dV, rho;

    dudx [ 0 ] [ 0 ] = b [ 0 ] * u.at(1) + b [ 1 ] * u.at(3) + b [ 2 ] * u.at(5);
    dudx [ 0 ] [ 1 ] = c [ 0 ] * u.at(1) + c [ 1 ] * u.at(3) + c [ 2 ] * u.at(5);
    dudx [ 1 ] [ 0 ] = b [ 0 ] * u.at(2) + b [ 1 ] * u.at(4) + b [ 2 ] * u.at(6);
    dudx [ 1 ] [ 1 ] = c [ 0 ] * u.at(2) + c [ 1 ] * u.at(4) + c [ 2 ] * u.at(6);

    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        dV = this->computeVolumeAround(gp);
        rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity(gp);
        this->computeNVector(n, gp);

        _u = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
        _v = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);

        // dN(v)/dv
        for ( int i = 1; i <= 2; i++ ) { // test function index
            for ( int k = 1; k <= 3; k++ ) { // nodal val of function w
                for ( int j = 1; j <= 2; j++ ) { // velocity vector component
                    for ( int m = 1; m <= 3; m++ ) { //  nodal components
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
        for ( int i = 1; i <= 2; i++ ) { // test function index
            for ( int k = 1; k <= 3; k++ ) { // nodal val of function w
                for ( int j = 1; j <= 2; j++ ) { // velocity vector component
                    for ( int m = 1; m <= 3; m++ ) { //  nodal components
                        w_dof_addr = ( k - 1 ) * 2 + i;
                        u_dof_addr = ( m - 1 ) * 2 + j;
                        //d1j = ( j == 1 );
                        //d2j = ( j == 2 );
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
TR1_2D_SUPG_AXI :: computeDiffusionTerm_MB(FloatArray &answer, TimeStep *tStep)
{
    answer.resize(6);
    answer.zero();
    FloatArray u, eps, stress;
    double Re = static_cast< FluidModel * >( domain->giveEngngModel() )->giveReynoldsNumber();
    this->computeVectorOfVelocities(VM_Total, tStep, u);
    FluidDynamicMaterial *mat = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();
    FloatMatrix _b;

    // stabilization term K_delta
    FloatArray un, n;
    double _u, _v, _r;
    this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), un);
    // end k_delta declaration

    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        double dV = this->computeVolumeAround(gp);
        this->computeBMtrx(_b, gp);
        eps.beProductOf(_b, u);
        mat->computeDeviatoricStressVector(stress, gp, eps, tStep);
        answer.plusProduct(_b, stress, dV / Re);

#if 1
        // stabilization term k_delta
        _r = this->computeRadiusAt(gp);
        this->computeNVector(n, gp);
        _u = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
        _v = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);

        for ( int i = 1; i <= 3; i++ ) {
            answer.at(2 * i - 1) -= t_supg * ( _u * b [ i - 1 ] + _v * c [ i - 1 ] ) * ( stress.at(1) / _r ) * dV / Re;
            answer.at(2 * i)     -= t_supg * ( _u * b [ i - 1 ] + _v * c [ i - 1 ] ) * ( stress.at(4) / _r ) * dV / Re;
        }

#endif
    }
}


void
TR1_2D_SUPG_AXI :: computeDiffusionDerivativeTerm_MB(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep)
{
    answer.resize(6, 6);
    answer.zero();
    FloatMatrix _db, _d, _b;
    double Re = static_cast< FluidModel * >( domain->giveEngngModel() )->giveReynoldsNumber();
    FluidDynamicMaterial *mat = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();

    // stabilization term K_delta
    FloatArray un, u, n, eps, stress;
    double _u, _v, _r;
    this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), un);
    this->computeVectorOfVelocities(VM_Total, tStep, u);


    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        double dV = this->computeVolumeAround(gp);
        this->computeBMtrx(_b, gp);
        mat->giveDeviatoricStiffnessMatrix(_d, mode, gp, tStep);
        _db.beProductOf(_d, _b);
        answer.plusProductUnsym(_b, _db, dV);
        //answer.plusProductSymmUpper (_bs,_db,dV*t_supg);
        // }

#if 1
        _r = this->computeRadiusAt(gp);
        this->computeNVector(n, gp);
        eps.beProductOf(_b, u);
        mat->computeDeviatoricStressVector(stress, gp, eps, tStep);
        //_mu = mat->giveCharacteristicValue(MRM_Viscosity, gp, tStep);

        _u = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
        _v = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);

        for ( int i = 1; i <= 3; i++ ) {
            for ( int j = 1; j <= 6; j++ ) {
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
TR1_2D_SUPG_AXI :: computePressureTerm_MB(FloatMatrix &answer, TimeStep *tStep)
{
    answer.resize(6, 3);
    answer.zero();
    FloatArray un, n;
    double _u, _v;

    this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), un);

    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        double dV = this->computeVolumeAround(gp);
        double _r = this->computeRadiusAt(gp);
        this->computeNVector(n, gp);
        // G matrix
        for ( int i = 1; i <= 3; i++ ) {
            for ( int j = 1; j <= 3; j++ ) {
                answer.at(2 * i - 1, j) -= b [ i - 1 ] * n.at(j) * dV;
                answer.at(2 * i, j)   -= c [ i - 1 ] * n.at(j) * dV;

                answer.at(2 * i - 1, j) -= n.at(i) * n.at(j) * dV / _r;
            }
        }

        // stabilization term (G_\delta mtrx)
        _u = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
        _v = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);
        for ( int i = 1; i <= 3; i++ ) {
            for ( int j = 1; j <= 3; j++ ) {
                answer.at(2 * i - 1, j) += ( _u * b [ i - 1 ] + _v * c [ i - 1 ] ) * b [ j - 1 ] * dV * t_supg;
                answer.at(2 * i, j)   += ( _u * b [ i - 1 ] + _v * c [ i - 1 ] ) * c [ j - 1 ] * dV * t_supg;
            }
        }
    }
}

void
TR1_2D_SUPG_AXI :: computeLSICStabilizationTerm_MB(FloatMatrix &answer, TimeStep *tStep)
{
    answer.resize(6, 6);
    answer.zero();
    double n[] = {
        b [ 0 ], c [ 0 ], b [ 1 ], c [ 1 ], b [ 2 ], c [ 2 ]
    };

    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        double dV = this->computeVolumeAround(gp);
        double rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity(gp);

        for ( int i = 1; i <= 6; i++ ) {
            for ( int j = 1; j <= 6; j++ ) {
                answer.at(i, j) += dV * t_lsic * rho * n [ i - 1 ] * n [ j - 1 ];
            }
        }
    }
}



void
TR1_2D_SUPG_AXI :: computeLinearAdvectionTerm_MC(FloatMatrix &answer, TimeStep *tStep)
{
    answer.resize(3, 6);
    answer.zero();

    FloatArray n;
    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        double dV = this->computeVolumeAround(gp);
        double _r = this->computeRadiusAt(gp);
        this->computeNVector(n, gp);
        for ( int i = 1; i <= 3; i++ ) {
            for ( int j = 1; j <= 3; j++ ) {
                answer.at(j, 2 * i - 1) += b [ i - 1 ] * n.at(j) * dV;
                answer.at(j, 2 * i)   += c [ i - 1 ] * n.at(j) * dV;

                answer.at(i, 1 + ( j - 1 ) * 2) += n.at(i) * n.at(j) * dV / _r;
            }
        }
    }
}

void
TR1_2D_SUPG_AXI :: computeAdvectionTerm_MC(FloatArray &answer, TimeStep *tStep)
{
    // N_epsilon (due to PSPG stabilization)
    double dudx, dudy, dvdx, dvdy, _u, _v;
    FloatArray u, un, n;

    answer.resize(3);
    answer.zero();

    this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), un);
    this->computeVectorOfVelocities(VM_Total, tStep, u);
    dudx = b [ 0 ] * u.at(1) + b [ 1 ] * u.at(3) + b [ 2 ] * u.at(5);
    dudy = c [ 0 ] * u.at(1) + c [ 1 ] * u.at(3) + c [ 2 ] * u.at(5);
    dvdx = b [ 0 ] * u.at(2) + b [ 1 ] * u.at(4) + b [ 2 ] * u.at(6);
    dvdy = c [ 0 ] * u.at(2) + c [ 1 ] * u.at(4) + c [ 2 ] * u.at(6);

    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        double dV = this->computeVolumeAround(gp);
        this->computeNVector(n, gp);

        _u = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
        _v = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);

        for ( int i = 1; i <= 3; i++ ) {
            answer.at(i) += t_pspg * dV * ( b [ i - 1 ] * ( _u * dudx + _v * dudy ) + c [ i - 1 ] * ( _u * dvdx + _v * dvdy ) );
        }
    }
}


void
TR1_2D_SUPG_AXI :: computeAdvectionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *tStep)
{
    answer.resize(3, 6);
    answer.zero();
    int w_dof_addr, u_dof_addr, d1j, d2j, km1, mm1;
    FloatArray u, un, n;
    double _u, _v;

    this->computeVectorOfVelocities(VM_Total, tStep, u);
    this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), un);

    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        double dV = this->computeVolumeAround(gp);
        this->computeNVector(n, gp);

        _u = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
        _v = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);

        for ( int k = 1; k <= 3; k++ ) { // nodal val of function w
            km1 = k - 1;
            for ( int j = 1; j <= 2; j++ ) { // velocity vector component
                for ( int m = 1; m <= 3; m++ ) { //  nodal components
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

void TR1_2D_SUPG_AXI :: computeDiffusionTerm_MC(FloatArray &answer, TimeStep *tStep)
{
    answer.resize(3);
    answer.zero();

#if 1
    double Re = static_cast< FluidModel * >( domain->giveEngngModel() )->giveReynoldsNumber();
    FloatArray eps, stress, u;
    FloatMatrix _b;
    FluidDynamicMaterial *mat = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();

    // stabilization term K_eps
    this->computeVectorOfVelocities(VM_Total, tStep, u);

    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        double dV = this->computeVolumeAround(gp);
        double rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity(gp);
        double _r = this->computeRadiusAt(gp);
        this->computeBMtrx(_b, gp);
        eps.beProductOf(_b, u);
        mat->computeDeviatoricStressVector(stress, gp, eps, tStep);
        stress.times(1. / Re);
        for ( int i = 1; i <= 3; i++ ) {
            answer.at(i) -= t_pspg * ( b [ i - 1 ] * stress.at(1) + c [ i - 1 ] * stress.at(4) ) * dV / rho / _r;
        }
    }

#endif
}

void TR1_2D_SUPG_AXI :: computeDiffusionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *tStep)
{
    answer.resize(3, 6);
    answer.zero();

#if 1
    double Re = static_cast< FluidModel * >( domain->giveEngngModel() )->giveReynoldsNumber();
    FluidDynamicMaterial *mat = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial();
    FloatMatrix _d, _b, _db;
    //FloatArray eps, stress,u;

    // stabilization term K_eps
    //this -> computeVectorOfVelocities(VM_Total,tStep, u) ;
    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        double dV = this->computeVolumeAround(gp);
        double rho = mat->give('d', gp);
        double _r = this->computeRadiusAt(gp);
        this->computeBMtrx(_b, gp);
        mat->giveDeviatoricStiffnessMatrix(_d, TangentStiffness, gp, tStep);
        _db.beProductOf(_d, _b);
        //eps.beProductOf (_b, u);
        //mat->computeDeviatoricStressVector (stress,gp,eps,tStep);

        for ( int i = 1; i <= 3; i++ ) {
            for ( int j = 1; j <= 6; j++ ) {
                answer.at(i, j) -= t_pspg * ( b [ i - 1 ] * ( _db.at(1, j) / _r ) + c [ i - 1 ] * ( _db.at(4, j) / _r ) ) * dV / rho;
            }
        }
    }

    answer.times(1. / Re);
#endif
}


void
TR1_2D_SUPG_AXI :: computeAccelerationTerm_MC(FloatMatrix &answer, TimeStep *tStep)
{
    answer.resize(3, 6);
    answer.zero();
    FloatArray n;
    // M_\epsilon

    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        double dV = this->computeVolumeAround(gp);
        this->computeNVector(n, gp);

        for ( int i = 1; i <= 3; i++ ) {
            for ( int j = 1; j <= 3; j++ ) {
                answer.at(i, 2 * j - 1) += t_pspg * dV * b [ i - 1 ] * n.at(j);
                answer.at(i, 2 * j)     += t_pspg * dV * c [ i - 1 ] * n.at(j);
            }
        }
    }
}


void
TR1_2D_SUPG_AXI :: computePressureTerm_MC(FloatMatrix &answer, TimeStep *tStep)
{
    answer.resize(3, 3);
    answer.zero();

    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        double dV = this->computeVolumeAround(gp);
        double rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity(gp);
        double coeff = t_pspg / rho;

        for ( int i = 1; i <= 3; i++ ) {
            for ( int j = 1; j <= 3; j++ ) {
                answer.at(i, j) += dV * coeff * ( b [ i - 1 ] * b [ j - 1 ] + c [ i - 1 ] * c [ j - 1 ] );
            }
        }
    }
}

void
TR1_2D_SUPG_AXI :: computeDeviatoricStrain(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    /* one should call material driver instead */
    FloatArray u;
    FloatMatrix _b;

    this->computeVectorOfVelocities(VM_Total, tStep, u);
    this->computeBMtrx(_b, gp);
    answer.beProductOf(_b, u);
}

void
TR1_2D_SUPG_AXI :: computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    /* one should call material driver instead */
    FloatArray eps;

    this->computeDeviatoricStrain(eps, gp, tStep);
    static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial()->computeDeviatoricStressVector(answer, gp, eps, tStep);
}



void
TR1_2D_SUPG_AXI :: computeBCRhsTerm_MB(FloatArray &answer, TimeStep *tStep)
{
    answer.resize(6);
    answer.zero();

    int nLoads;
    FloatArray un, gVector;
    FloatArray nV;

    // add body load (gravity) termms
    this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), un);
    nLoads = this->giveBodyLoadArray()->giveSize();
    for ( int i = 1; i <= nLoads; i++ ) {
        Load *load = domain->giveLoad( bodyLoadArray.at(i) );
        bcGeomType ltype = load->giveBCGeoType();
        if ( ltype == BodyLoadBGT && load->giveBCValType() == ForceLoadBVT ) {
            load->computeComponentArrayAt(gVector, tStep, VM_Total);
            if ( gVector.giveSize() ) {
                for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
                    double dV = this->computeVolumeAround(gp);
                    double rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity(gp);
                    double coeff = rho * dV;
                    this->computeNVector(nV, gp);
                    double u = nV.at(1) * un.at(1) + nV.at(2) * un.at(3) + nV.at(3) * un.at(5);
                    double v = nV.at(1) * un.at(2) + nV.at(2) * un.at(4) + nV.at(3) * un.at(6);

                    answer.at(1) += coeff * ( gVector.at(1) * ( nV.at(1) + t_supg * ( b [ 0 ] * u + c [ 0 ] * v ) ) );
                    answer.at(2) += coeff * ( gVector.at(2) * ( nV.at(1) + t_supg * ( b [ 0 ] * u + c [ 0 ] * v ) ) );
                    answer.at(3) += coeff * ( gVector.at(1) * ( nV.at(2) + t_supg * ( b [ 1 ] * u + c [ 1 ] * v ) ) );
                    answer.at(4) += coeff * ( gVector.at(2) * ( nV.at(2) + t_supg * ( b [ 1 ] * u + c [ 1 ] * v ) ) );
                    answer.at(5) += coeff * ( gVector.at(1) * ( nV.at(3) + t_supg * ( b [ 2 ] * u + c [ 2 ] * v ) ) );
                    answer.at(6) += coeff * ( gVector.at(2) * ( nV.at(3) + t_supg * ( b [ 2 ] * u + c [ 2 ] * v ) ) );
                }
            }
        }
    }

    //answer.times(rc);

    // loop over sides
    int n1, n2, code, sid;
    double tx, ty, l, side_r;
    //IntArray nodecounter (3);
    for ( int j = 1; j <= boundarySides.giveSize(); j++ ) {
        code = boundaryCodes.at(j);
        sid = boundarySides.at(j);
        if ( ( code & FMElement_PrescribedTractionBC ) ) {
            FloatArray t, coords(1);
            int n, id;
            // integrate tractions
            n1 = sid;
            n2 = ( n1 == 3 ? 1 : n1 + 1 );

            tx = giveNode(n2)->giveCoordinate(1) - giveNode(n1)->giveCoordinate(1);
            ty = giveNode(n2)->giveCoordinate(2) - giveNode(n1)->giveCoordinate(2);
            l = sqrt(tx * tx + ty * ty);
            // radius at side center
            side_r = 0.5 * ( giveNode(n2)->giveCoordinate(1) + giveNode(n1)->giveCoordinate(1) );

            // if no traction bc applied but side marked as with traction load
            // then zero traction is assumed !!!

            // loop over boundary load array
            int numLoads = this->giveBoundaryLoadArray()->giveSize() / 2;
            for ( int i = 1; i <= numLoads; i++ ) {
                n = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
                id = boundaryLoadArray.at(i * 2);
                if ( id != sid ) {
                    continue;
                }

                auto load = dynamic_cast< BoundaryLoad * >( domain->giveLoad(n) );
                if ( load ) {
                    load->computeValueAt(t, tStep, coords, VM_Total);

                    // here it is assumed constant traction, one point integration only
                    // n1 (u,v)
                    answer.at(n1 * 2 - 1) += t.at(1) * side_r * l / 2.;
                    answer.at(n1 * 2)     += t.at(2) * side_r * l / 2.;
                    // n2 (u,v)
                    answer.at(n2 * 2 - 1) += t.at(1) * side_r * l / 2.;
                    answer.at(n2 * 2)     += t.at(2) * side_r * l / 2.;

                    //answer.at(n1)+= (t.at(1)*nx + t.at(2)*ny) * side_r * l/2.;
                    //answer.at(n2)+= (t.at(1)*nx + t.at(2)*ny) * side_r * l/2.;
                }
            }
        }
    }
}


void
TR1_2D_SUPG_AXI :: computeBCRhsTerm_MC(FloatArray &answer, TimeStep *tStep)
{
    int nLoads;
    FloatArray gVector;

    answer.resize(3);
    answer.zero();
    nLoads = this->giveBodyLoadArray()->giveSize();
    for ( int i = 1; i <= nLoads; i++ ) {
        Load *load  = domain->giveLoad( bodyLoadArray.at(i) );
        bcGeomType ltype = load->giveBCGeoType();
        if ( ( ltype == BodyLoadBGT ) && ( load->giveBCValType() == ForceLoadBVT ) ) {
            load->computeComponentArrayAt(gVector, tStep, VM_Total);
            if ( gVector.giveSize() ) {
                for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
                    double dV = this->computeVolumeAround(gp);
                    double coeff = t_pspg * dV;

                    answer.at(1) += coeff * ( b [ 0 ] * gVector.at(1) + c [ 0 ] * gVector.at(2) );
                    answer.at(2) += coeff * ( b [ 1 ] * gVector.at(1) + c [ 1 ] * gVector.at(2) );
                    answer.at(3) += coeff * ( b [ 2 ] * gVector.at(1) + c [ 2 ] * gVector.at(2) );
                }
            }
        }
    }
}


void
TR1_2D_SUPG_AXI :: computeLoadVector(FloatArray &answer, Load *load, CharType type, ValueModeType mode, TimeStep *tStep)
{
    if ( type != ExternalForcesVector ) {
        answer.clear();
        return;
    }

    FloatArray un, nV;

    this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), un);

    answer.resize(9);

    if ( load->giveBCValType() == ForceLoadBVT ) {
        FloatArray gVector;
        load->computeComponentArrayAt(gVector, tStep, VM_Total);

        for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
            double dV = this->computeVolumeAround(gp);
            double rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity(gp);
            double coeff;
            this->computeNVector(nV, gp);
            double u = nV.at(1) * un.at(1) + nV.at(2) * un.at(3) + nV.at(3) * un.at(5);
            double v = nV.at(1) * un.at(2) + nV.at(2) * un.at(4) + nV.at(3) * un.at(6);

            coeff = rho * dV;
            answer.at(1) += coeff * gVector.at(1) * ( nV.at(1) + t_supg * ( b [ 0 ] * u + c [ 0 ] * v ) );
            answer.at(2) += coeff * gVector.at(2) * ( nV.at(1) + t_supg * ( b [ 0 ] * u + c [ 0 ] * v ) );
            answer.at(4) += coeff * gVector.at(1) * ( nV.at(2) + t_supg * ( b [ 1 ] * u + c [ 1 ] * v ) );
            answer.at(5) += coeff * gVector.at(2) * ( nV.at(2) + t_supg * ( b [ 1 ] * u + c [ 1 ] * v ) );
            answer.at(7) += coeff * gVector.at(1) * ( nV.at(3) + t_supg * ( b [ 2 ] * u + c [ 2 ] * v ) );
            answer.at(8) += coeff * gVector.at(2) * ( nV.at(3) + t_supg * ( b [ 2 ] * u + c [ 2 ] * v ) );

            coeff = t_pspg * dV;
            answer.at(3) += coeff * ( b [ 0 ] * gVector.at(1) + c [ 0 ] * gVector.at(2) );
            answer.at(6) += coeff * ( b [ 1 ] * gVector.at(1) + c [ 1 ] * gVector.at(2) );
            answer.at(9) += coeff * ( b [ 2 ] * gVector.at(1) + c [ 2 ] * gVector.at(2) );
        }
    }
}


void
TR1_2D_SUPG_AXI :: computeBCLhsTerm_MB(FloatMatrix &answer, TimeStep *tStep)
{
    int nLoads;
    FloatMatrix helpMatrix;

    answer.clear();

    nLoads = this->giveBoundaryLoadArray()->giveSize() / 2;

    for ( int i = 1; i <= nLoads; i++ ) {
        int n = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
        int side = boundaryLoadArray.at(i * 2);
        Load *load = domain->giveLoad(n);
        bcType boundarytype = load->giveType();
        if ( boundarytype == SlipWithFriction ) {
            this->computeSlipWithFrictionBCTerm_MB(helpMatrix, load, side, tStep);
        } else if ( boundarytype == PenetrationWithResistance ) {
            this->computePenetrationWithResistanceBCTerm_MB(helpMatrix, load, side, tStep);
        } else {
            helpMatrix.clear();
            // OOFEM_ERROR("unsupported load type class");
        }

        answer.add(helpMatrix);
    }
}


void
TR1_2D_SUPG_AXI :: computeSlipWithFrictionBCTerm_MB(FloatMatrix &answer, Load *load, int side, TimeStep *tStep)
{
    int node1, node2;
    double l, t1, t2, _t1, _t2;
    double beta;
    FloatArray nt(6), n;
    answer.clear();

    BoundaryLoad *edgeLoad = static_cast< BoundaryLoad * >(load);
    beta = edgeLoad->giveProperty('a', tStep);
    node1 = side;
    node2 = ( node1 == 3 ? 1 : node1 + 1 );

    _t1 = giveNode(node2)->giveCoordinate(1) - giveNode(node1)->giveCoordinate(1);
    _t2 = giveNode(node2)->giveCoordinate(2) - giveNode(node1)->giveCoordinate(2);
    l = sqrt(_t1 * _t1 + _t2 * _t2);

    t1 = _t1 / l;
    t2 = _t2 / l;

    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        double dV = this->computeVolumeAround(gp);
        this->computeNVector(n, gp);

        nt.at(1) = n.at(1) * t1;
        nt.at(2) = n.at(1) * t2;
        nt.at(3) = n.at(2) * t1;
        nt.at(4) = n.at(2) * t2;
        nt.at(5) = n.at(3) * t1;
        nt.at(6) = n.at(3) * t2;

        answer.plusDyadSymmUpper(nt, beta * dV);
    }
    answer.symmetrized();
}

void
TR1_2D_SUPG_AXI :: computePenetrationWithResistanceBCTerm_MB(FloatMatrix &answer, Load *load, int side, TimeStep *tStep)
{
    int node1, node2;
    double l, n1, n2, _t1, _t2;
    double alpha;
    FloatArray nt(6), n;
    answer.clear();

    BoundaryLoad *edgeLoad = static_cast< BoundaryLoad * >(load);
    alpha = edgeLoad->giveProperty('a', tStep);
    node1 = side;
    node2 = ( node1 == 3 ? 1 : node1 + 1 );


    _t1 = giveNode(node2)->giveCoordinate(1) - giveNode(node1)->giveCoordinate(1);
    _t2 = giveNode(node2)->giveCoordinate(2) - giveNode(node1)->giveCoordinate(2);
    l = sqrt(_t1 * _t1 + _t2 * _t2);

    n1 = _t2 / l;
    n2 = -_t1 / l;


    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        double dV = this->computeVolumeAround(gp);
        this->computeNVector(n, gp);

        nt.at(1) = n.at(1) * n1;
        nt.at(2) = n.at(1) * n2;
        nt.at(3) = n.at(2) * n1;
        nt.at(4) = n.at(2) * n2;
        nt.at(5) = n.at(3) * n1;
        nt.at(6) = n.at(3) * n2;

        answer.plusDyadSymmUpper(nt, dV / alpha);
    }
    answer.symmetrized();
}

void
TR1_2D_SUPG_AXI :: computeOutFlowBCTerm_MB(FloatMatrix &answer, int side, TimeStep *tStep)
{
    int node1, node2;
    double l, n1, n2, t1, t2, dV;
    FloatArray nt(6), n;
    answer.clear();
    //beta
    //area
    node1 = side;
    node2 = ( node1 == 3 ? 1 : node1 + 1 );

    t1 = giveNode(node2)->giveCoordinate(1) - giveNode(node1)->giveCoordinate(1);
    t2 = giveNode(node2)->giveCoordinate(2) - giveNode(node1)->giveCoordinate(2);
    l = sqrt(t1 * t1 + t2 * t2);

    n1 = t2 / l;
    n2 = -t1 / l;


    for ( GaussPoint *gp: *integrationRulesArray [ 0 ] ) {
        dV = this->computeVolumeAround(gp);
        this->computeNVector(n, gp);

        nt.at(1) = n.at(1) * n1;
        nt.at(2) = n.at(1) * n2;
        nt.at(3) = n.at(2) * n1;
        nt.at(4) = n.at(2) * n2;
        nt.at(5) = n.at(3) * n1;
        nt.at(6) = n.at(3) * n2;

        answer.plusDyadUnsym(nt, n, dV);
    }

    answer.negated();
}



void
TR1_2D_SUPG_AXI :: computeBCLhsPressureTerm_MB(FloatMatrix &answer, TimeStep *tStep)
{
    int nLoads = 0;
    FloatMatrix helpMatrix;
    // loop over boundary load array

    answer.clear();

    nLoads = this->giveBoundaryLoadArray()->giveSize() / 2;

    for ( int i = 1; i <= nLoads; i++ ) {
        int n = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
        int side = boundaryLoadArray.at(i * 2);
        Load *load = domain->giveLoad(n);
        bcType boundarytype = load->giveType();
        if ( boundarytype == OutFlowBC ) {
            this->computeOutFlowBCTerm_MB(helpMatrix, side, tStep);
        } else {
            helpMatrix.clear();
            // OOFEM_ERROR("unsupported load type class");
        }

        answer.add(helpMatrix);
    }
}


double
TR1_2D_SUPG_AXI :: computeRadiusAt(GaussPoint *gp)
{
    FloatArray n;
    this->computeNVector(n, gp);

    return n.at(1) * this->giveNode(1)->giveCoordinate(1)
           + n.at(2) * this->giveNode(2)->giveCoordinate(1)
           + n.at(3) * this->giveNode(3)->giveCoordinate(1);
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
}

void
TR1_2D_SUPG_AXI :: computeNVector(FloatArray &n, GaussPoint *gp)
{
    this->interp.evalN( n, gp->giveSubPatchCoordinates(), FEIElementGeometryWrapper(this) );
}

double
TR1_2D_SUPG_AXI :: computeVolumeAround(GaussPoint *gp)
{
    double _r, weight, detJ;

    detJ = fabs( this->interp.giveTransformationJacobian( gp->giveSubPatchCoordinates(), FEIElementGeometryWrapper(this) ) );
    weight = gp->giveWeight();
    _r = computeRadiusAt(gp);

    return detJ * weight * _r;
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
TR1_2D_SUPG_AXI :: updateStabilizationCoeffs(TimeStep *tStep)
{
    /* UGN-Based Stabilization */
    double h_ugn, sum = 0.0, vnorm, t_sugn1, t_sugn2, t_sugn3, u_1, u_2;
    double dscale, uscale, lscale, tscale, dt;
    //bool zeroFlag = false;
    int im1;
    FloatArray u;

    uscale = domain->giveEngngModel()->giveVariableScale(VST_Velocity);
    lscale = domain->giveEngngModel()->giveVariableScale(VST_Length);
    tscale = domain->giveEngngModel()->giveVariableScale(VST_Time);
    dscale = domain->giveEngngModel()->giveVariableScale(VST_Density);

    this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), u);
    u.times(uscale);
    double nu;

    // compute averaged viscosity based on rule of mixture
    GaussPoint *gp;
    if ( integrationRulesArray [ 0 ]->giveNumberOfIntegrationPoints() ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    } else {
        gp = integrationRulesArray [ 1 ]->getIntegrationPoint(0);
    }

    nu = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial()->giveEffectiveViscosity( gp, tStep->givePreviousStep() );
    nu *= domain->giveEngngModel()->giveVariableScale(VST_Viscosity);

    dt = tStep->giveTimeIncrement() * tscale;

    for ( int i = 1; i <= 3; i++ ) {
        im1 = i - 1;
        sum += fabs(u.at( ( im1 ) * 2 + 1 ) * b [ im1 ] / lscale + u.at(im1 * 2 + 2) * c [ im1 ] / lscale);
    }

    /*
     * u_1=(u.at(1)+u.at(3)+u.at(5))/3.0;
     * u_2=(u.at(2)+u.at(4)+u.at(6))/3.0;
     * vnorm=sqrt(u_1*u_1+u_2*u_2);
     */
    vnorm = 0.;
    for ( int i = 1; i <= 3; i++ ) {
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
    int neg = 0, pos = 0, zero = 0, si = 0;
    double x1, x2, x3, y1, y2, y3;

    answer.resize(2);
    for ( double ifi: fi ) {
        if ( ifi >= 0. ) {
            pos++;
        } else if ( ifi < 0.0 ) {
            neg++;
        } else {
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
        for ( int i = 1; i <= 3; i++ ) {
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

            double t = fi.at(si) / ( fi.at(si) - fi.at(next_node) );
            x2 = x1 + t * ( this->giveNode(next_node)->giveCoordinate(1) - this->giveNode(si)->giveCoordinate(1) );
            y2 = y1 + t * ( this->giveNode(next_node)->giveCoordinate(2) - this->giveNode(si)->giveCoordinate(2) );

            t = fi.at(si) / ( fi.at(si) - fi.at(prev_node) );
            x3 = x1 + t * ( this->giveNode(prev_node)->giveCoordinate(1) - this->giveNode(si)->giveCoordinate(1) );
            y3 = y1 + t * ( this->giveNode(prev_node)->giveCoordinate(2) - this->giveNode(si)->giveCoordinate(2) );


            // compute volume associated to triangle (x1,y1; x2,y2; x3,y3)
            double __volume = 0.5 * ( x2 * y3 + x1 * y2 + y1 * x3 - x2 * y1 - x3 * y2 - x1 * y3 ) * ( ( x1 + x2 + x3 ) / 3. );
            double volume = this->area * ( ( this->giveNode(1)->giveCoordinate(1) + this->giveNode(2)->giveCoordinate(1) + this->giveNode(3)->giveCoordinate(1) ) / 3. );
            if ( fabs(__volume) / volume > 1.00001 ) {
                OOFEM_ERROR("internal consistency error");
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
            OOFEM_ERROR("internal consistency error");
        }
    }
}
} // end namespace oofem
