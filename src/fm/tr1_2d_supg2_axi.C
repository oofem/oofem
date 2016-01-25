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

#include "tr1_2d_supg2_axi.h"
#include "fluidmodel.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "engngm.h"
#include "fluiddynamicmaterial.h"
#include "fluidcrosssection.h"
#include "load.h"
#include "timestep.h"
#include "boundaryload.h"
#include "fei2dtrlin.h"
#include "fei2dquadlin.h"
#include "geotoolbox.h"
#include "crosssection.h"
#include "dynamicinputrecord.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif

namespace oofem {
#define TRSUPG_ZERO_VOF 1.e-8

REGISTER_Element(TR1_2D_SUPG2_AXI);

//#define TR1_2D_SUPG2_AXI_DEBUG

TR1_2D_SUPG2_AXI :: TR1_2D_SUPG2_AXI(int n, Domain *aDomain) :
    TR1_2D_SUPG(n, aDomain)
    // Constructor.
{
    numberOfDofMans = 3;
}

TR1_2D_SUPG2_AXI :: ~TR1_2D_SUPG2_AXI()
// Destructor
{ }


IRResultType
TR1_2D_SUPG2_AXI :: initializeFrom(InputRecord *ir)
{
    IRResultType result;               // Required by IR_GIVE_FIELD macro

    this->vof = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, vof, _IFT_Tr1SUPG_pvof);
    if ( vof > 0.0 ) {
        setPermanentVolumeFraction(vof);
        this->temp_vof = this->vof;
    } else {
        this->vof = 0.0;
        IR_GIVE_OPTIONAL_FIELD(ir, vof, _IFT_Tr1SUPG_vof);
        this->temp_vof = this->vof;
    }

    this->mat [ 0 ] = this->mat [ 1 ] = this->material;
    IR_GIVE_OPTIONAL_FIELD(ir, mat [ 0 ], _IFT_Tr1SUPG2_mat0);
    IR_GIVE_OPTIONAL_FIELD(ir, mat [ 1 ], _IFT_Tr1SUPG2_mat1);
    this->material = this->mat [ 0 ];

    result = SUPGElement :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }
    this->initGeometry();
    return IRRT_OK;
}


void
TR1_2D_SUPG2_AXI :: giveInputRecord(DynamicInputRecord &input)
{
    SUPGElement :: giveInputRecord(input);
    if ( this->permanentVofFlag ) {
        input.setField(this->vof, _IFT_Tr1SUPG_pvof);
    } else {
        input.setField(this->vof, _IFT_Tr1SUPG_vof);
    }

    input.setField(this->mat [ 0 ], _IFT_Tr1SUPG2_mat0);
    input.setField(this->mat [ 1 ], _IFT_Tr1SUPG2_mat1);
}

void
TR1_2D_SUPG2_AXI :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize( 2 );
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3, true) );
        integrationRulesArray [ 1 ].reset( new GaussIntegrationRule(2, this, 1, 3, true) );
    }
}


/*
 * Integration template
 * // loop over each fluid
 * for (ifluid = 0; ifluid< 2; ifluid++) {
 * for (GaussPoint *gp: *integrationRulesArray [ ifluid ] ) {
 *  mapped_gp = sub_mapped_IPRule[ifluid]->getIntegrationPoint(ip) ;
 *  this->computeNMtrx (n, mapped_gp);
 *  dV = this->computeVolumeAroundID(gp,id[ifluid], vcoords[ifluid]) ;
 *    // compute integral here
 * }
 * }
 */



void
TR1_2D_SUPG2_AXI :: computeAccelerationTerm_MB(FloatMatrix &answer, TimeStep *tStep)
{
    answer.resize(6, 6);
    answer.zero();
    FloatArray n, un;

    double _val, u, v;

    this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), un);

    // consistent mass
    // loop over each fluid
    for ( int ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( GaussPoint *gp: *integrationRulesArray [ ifluid ] ) {
            double rho = this->_giveMaterial(ifluid)->give('d', gp);
            double dV = this->computeVolumeAroundID(gp, id [ ifluid ], vcoords [ ifluid ]);
            this->computeNMtrx(n, gp);

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
}


void
TR1_2D_SUPG2_AXI :: computeAdvectionTerm_MB(FloatArray &answer, TimeStep *tStep)
{
    answer.resize(6);
    answer.zero();

    FloatArray n, u, un;
    double dudx, dudy, dvdx, dvdy, _u, _v;
    this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), un);
    this->computeVectorOfVelocities(VM_Total, tStep, u);


    dudx = b [ 0 ] * u.at(1) + b [ 1 ] * u.at(3) + b [ 2 ] * u.at(5);
    dudy = c [ 0 ] * u.at(1) + c [ 1 ] * u.at(3) + c [ 2 ] * u.at(5);
    dvdx = b [ 0 ] * u.at(2) + b [ 1 ] * u.at(4) + b [ 2 ] * u.at(6);
    dvdy = c [ 0 ] * u.at(2) + c [ 1 ] * u.at(4) + c [ 2 ] * u.at(6);

    // standard galerkin term
    for ( int ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( GaussPoint *gp: *integrationRulesArray [ ifluid ] ) {
            double rho = this->_giveMaterial(ifluid)->give('d', gp);
            double dV = this->computeVolumeAroundID(gp, id [ ifluid ], vcoords [ ifluid ]);
            this->computeNMtrx(n, gp);

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
}


void
TR1_2D_SUPG2_AXI :: computeAdvectionDerivativeTerm_MB(FloatMatrix &answer, TimeStep *tStep)
{
    answer.resize(6, 6);
    answer.zero();

    FloatArray u, un, n;
    this->computeVectorOfVelocities(VM_Total, tStep, u);
    this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), un);
    double _u, _v;
    int w_dof_addr, u_dof_addr, dij;

    // dN(v)/dv
    for ( int ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( GaussPoint *gp: *integrationRulesArray [ ifluid ] ) {
            double rho = this->_giveMaterial(ifluid)->give('d', gp);
            double dV = this->computeVolumeAroundID(gp, id [ ifluid ], vcoords [ ifluid ]);
            this->computeNMtrx(n, gp);

            _u = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
            _v = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);

            // dN(v)/dv
            for ( int i = 1; i <= 2; i++ ) { // test function index
                for ( int k = 1; k <= 3; k++ ) { // nodal val of function w
                    for ( int j = 1; j <= 2; j++ ) { // velocity vector component
                        for ( int m = 1; m <= 3; m++ ) { //  nodal components
                            w_dof_addr = ( k - 1 ) * 2 + i;
                            u_dof_addr = ( m - 1 ) * 2 + j;
                            //d1j = (j==1); d2j=(j==2);
                            dij = ( i == j );
                            answer.at(w_dof_addr, u_dof_addr) += dV * rho * n.at(k) * ( dij * _u * b [ m - 1 ] + dij * _v * c [ m - 1 ] );
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
                            //d1j = (j==1); d2j=(j==2);
                            dij = ( i == j );
                            answer.at(w_dof_addr, u_dof_addr) += dV * t_supg * rho *
                                                                 ( _u * b [ k - 1 ] + _v * c [ k - 1 ] ) * ( dij * _u * b [ m - 1 ] + dij * _v * c [ m - 1 ] );
                        }
                    }
                }
            }
        }
    }
}


void
TR1_2D_SUPG2_AXI :: computeDiffusionTerm_MB(FloatArray &answer, TimeStep *tStep)
{
    answer.resize(6);
    answer.zero();
    FloatArray u, un, eps, stress;
    double Re = static_cast< FluidModel * >( domain->giveEngngModel() )->giveReynoldsNumber();
    //double dudx,dudy,dvdx,dvdy;

    this->computeVectorOfVelocities(VM_Total, tStep, u);
    FloatArray n;
    double _u, _v, _r;
    this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), un);
    FloatMatrix _b(4, 6);

    for ( int ifluid = 0; ifluid < 2; ifluid++ ) {
        FluidDynamicMaterial *mat = static_cast< FluidDynamicMaterial * >( this->_giveMaterial(ifluid) );
        for ( GaussPoint *gp: *integrationRulesArray [ ifluid ] ) {
            double dV = this->computeVolumeAroundID(gp, id [ ifluid ], vcoords [ ifluid ]);

            this->computeBMtrx(_b, gp);
            eps.beProductOf(_b, u);
            mat->computeDeviatoricStressVector(stress, gp, eps, tStep);
            answer.plusProduct(_b, stress, dV / Re);

#if 1
            // stabilization term k_delta
            _r = this->computeRadiusAt(gp);
            computeNVector(n, gp);
            _u = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
            _v = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);

            for ( int i = 1; i <= 3; i++ ) {
                answer.at(2 * i - 1) -= t_supg * ( _u * b [ i - 1 ] + _v * c [ i - 1 ] ) * ( stress.at(1) / _r ) * dV / Re;
                answer.at(2 * i)     -= t_supg * ( _u * b [ i - 1 ] + _v * c [ i - 1 ] ) * ( stress.at(4) / _r ) * dV / Re;
            }

#endif
        }
    }
}


void
TR1_2D_SUPG2_AXI :: computeDiffusionDerivativeTerm_MB(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep)
{
    //double dudx, dudy, dvdx, dvdy;
    answer.resize(6, 6);
    answer.zero();
    FloatMatrix _db, _d, _b;
    //FloatArray un;
    double Re = static_cast< FluidModel * >( domain->giveEngngModel() )->giveReynoldsNumber();
    FloatArray un, u, n, eps, stress;
    double _u, _v, _r;

    this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), un);
    this->computeVectorOfVelocities(VM_Total, tStep, u);

    for ( int ifluid = 0; ifluid < 2; ifluid++ ) {
        FluidDynamicMaterial *mat = static_cast< FluidDynamicMaterial * >( this->_giveMaterial(ifluid) );
        for ( GaussPoint *gp: *integrationRulesArray [ ifluid ] ) {
            double dV = this->computeVolumeAroundID(gp, id [ ifluid ], vcoords [ ifluid ]);

            this->computeBMtrx(_b, gp);
            mat->giveDeviatoricStiffnessMatrix(_d, mode, gp, tStep);
            _db.beProductOf(_d, _b);
            answer.plusProductUnsym(_b, _db, dV);
            //answer.plusProductSymmUpper (_bs,_db,dV*t_supg);
            // }

#if 1
            _r = this->computeRadiusAt(gp);
            computeNVector(n, gp);
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
    }

    answer.times(1. / Re);
}


void
TR1_2D_SUPG2_AXI :: computePressureTerm_MB(FloatMatrix &answer, TimeStep *tStep)
{
    answer.resize(6, 3);
    answer.zero();
    FloatArray un, n;
    double _u, _v;

    this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), un);
    double dV, _r;
    for ( int ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( GaussPoint *gp: *integrationRulesArray [ ifluid ] ) {
            dV = this->computeVolumeAroundID(gp, id [ ifluid ], vcoords [ ifluid ]);
            _r = this->computeRadiusAt(gp);
            computeNVector(n, gp);
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
}


void
TR1_2D_SUPG2_AXI :: computeLSICStabilizationTerm_MB(FloatMatrix &answer, TimeStep *tStep)
{
    answer.resize(6, 6);
    answer.zero();
    //double coeff = area*t_lsic*rho;
    double n[] = {
        b [ 0 ], c [ 0 ], b [ 1 ], c [ 1 ], b [ 2 ], c [ 2 ]
    };

    for ( int ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( GaussPoint *gp: *integrationRulesArray [ ifluid ] ) {
            double rho = this->_giveMaterial(ifluid)->give('d', gp);
            double dV = this->computeVolumeAroundID(gp, id [ ifluid ], vcoords [ ifluid ]);

            for ( int i = 1; i <= 6; i++ ) {
                for ( int j = 1; j <= 6; j++ ) {
                    answer.at(i, j) += dV * t_lsic * rho * n [ i - 1 ] * n [ j - 1 ];
                }
            }
        }
    }
}


void
TR1_2D_SUPG2_AXI :: computeLinearAdvectionTerm_MC(FloatMatrix &answer, TimeStep *tStep)
{
    answer.resize(3, 6);
    answer.zero();

    FloatArray n;

    for ( int ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( auto &gp : *integrationRulesArray [ ifluid ] ) {
            double dV = this->computeVolumeAroundID(gp, id [ ifluid ], vcoords [ ifluid ]);
            double _r = this->computeRadiusAt(gp);
            computeNVector(n, gp);
            for ( int i = 1; i <= 3; i++ ) {
                for ( int j = 1; j <= 3; j++ ) {
                    answer.at(j, 2 * i - 1) += b [ i - 1 ] * n.at(j) * dV;
                    answer.at(j, 2 * i)   += c [ i - 1 ] * n.at(j) * dV;

                    answer.at(i, 1 + ( j - 1 ) * 2) += n.at(i) * n.at(j) * dV / _r;
                }
            }
        }
    }
}

void
TR1_2D_SUPG2_AXI :: computeAdvectionTerm_MC(FloatArray &answer, TimeStep *tStep)
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

    for ( int ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( auto &gp : *integrationRulesArray [ ifluid ] ) {
            double dV = this->computeVolumeAroundID(gp, id [ ifluid ], vcoords [ ifluid ]);
            computeNVector(n, gp);

            _u = n.at(1) * un.at(1) + n.at(2) * un.at(3) + n.at(3) * un.at(5);
            _v = n.at(1) * un.at(2) + n.at(2) * un.at(4) + n.at(3) * un.at(6);

            for ( int i = 1; i <= 3; i++ ) {
                answer.at(i) += t_pspg * dV * ( b [ i - 1 ] * ( _u * dudx + _v * dudy ) + c [ i - 1 ] * ( _u * dvdx + _v * dvdy ) );
            }
        }
    }
}


void
TR1_2D_SUPG2_AXI :: computeAdvectionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *tStep)
{
    answer.resize(3, 6);
    answer.zero();
    int w_dof_addr, u_dof_addr, d1j, d2j, km1, mm1;
    FloatArray u, un, n;
    double _u, _v;

    this->computeVectorOfVelocities(VM_Total, tStep, u);
    this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), un);

    for ( int ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( auto &gp : *integrationRulesArray [ ifluid ] ) {
            double dV = this->computeVolumeAroundID(gp, id [ ifluid ], vcoords [ ifluid ]);
            computeNVector(n, gp);

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
}

void TR1_2D_SUPG2_AXI :: computeDiffusionTerm_MC(FloatArray &answer, TimeStep *tStep)
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

    for ( int ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( auto &gp : *integrationRulesArray [ ifluid ] ) {
            double rho = this->_giveMaterial(ifluid)->give('d', gp);
            double dV = this->computeVolumeAroundID(gp, id [ ifluid ], vcoords [ ifluid ]);
            double _r = this->computeRadiusAt(gp);
            this->computeBMtrx(_b, gp);
            eps.beProductOf(_b, u);
            mat->computeDeviatoricStressVector(stress, gp, eps, tStep);
            stress.times(1. / Re);
            for ( int i = 1; i <= 3; i++ ) {
                answer.at(i) -= t_pspg * ( b [ i - 1 ] * stress.at(1) + c [ i - 1 ] * stress.at(4) ) * dV / rho / _r;
            }
        }
    }

#endif
}

void TR1_2D_SUPG2_AXI :: computeDiffusionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *tStep)
{
    answer.resize(3, 6);
    answer.zero();

#if 1
    double Re = static_cast< FluidModel * >( domain->giveEngngModel() )->giveReynoldsNumber();
    FloatMatrix _d, _b, _db;
    //FloatArray eps, stress,u;

    // stabilization term K_eps
    //this -> computeVectorOfVelocities(VM_Total,tStep, u) ;
    for ( int ifluid = 0; ifluid < 2; ifluid++ ) {
        FluidDynamicMaterial *mat = static_cast< FluidDynamicMaterial * >( this->_giveMaterial(ifluid) );
        for ( GaussPoint *gp: *integrationRulesArray [ ifluid ] ) {
            double dV = this->computeVolumeAroundID(gp, id [ ifluid ], vcoords [ ifluid ]);
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
    }

    answer.times(1. / Re);

#endif
}




void
TR1_2D_SUPG2_AXI :: computeAccelerationTerm_MC(FloatMatrix &answer, TimeStep *tStep)
{
    answer.resize(3, 6);
    answer.zero();
    FloatArray n;
    // M_\epsilon

    for ( int ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( GaussPoint *gp: *integrationRulesArray [ ifluid ] ) {
            double dV = this->computeVolumeAroundID(gp, id [ ifluid ], vcoords [ ifluid ]);
            computeNVector(n, gp);

            for ( int i = 1; i <= 3; i++ ) {
                for ( int j = 1; j <= 3; j++ ) {
                    answer.at(i, 2 * j - 1) += t_pspg * dV * b [ i - 1 ] * n.at(j);
                    answer.at(i, 2 * j)   += t_pspg * dV * c [ i - 1 ] * n.at(j);
                }
            }
        }
    }
}

void
TR1_2D_SUPG2_AXI :: computePressureTerm_MC(FloatMatrix &answer, TimeStep *tStep)
{
    answer.resize(3, 3);
    answer.zero();

    for ( int ifluid = 0; ifluid < 2; ifluid++ ) {
        FluidDynamicMaterial *mat = static_cast< FluidDynamicMaterial * >( this->_giveMaterial(ifluid) );
        for ( GaussPoint *gp: *integrationRulesArray [ ifluid ] ) {
            double rho = mat->give('d', gp);
            double dV = this->computeVolumeAroundID(gp, id [ ifluid ], vcoords [ ifluid ]);


            for ( int i = 1; i <= 3; i++ ) {
                for ( int j = 1; j <= 3; j++ ) {
                    answer.at(i, j) += t_pspg * dV * ( b [ i - 1 ] * b [ j - 1 ] + c [ i - 1 ] * c [ j - 1 ] ) / rho;
                }
            }
        }
    }

#ifdef TR1_2D_SUPG2_DEBUG
    /* test */
    FloatMatrix test;
    std :: swap(integrationRulesArray [ 0 ], integrationRulesArray [ 1 ]);
    TR1_2D_SUPG :: computePressureTerm_MC(test, tStep);
    std :: swap(integrationRulesArray [ 0 ], integrationRulesArray [ 1 ]);
    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            if ( fabs( ( answer.at(i, j) - test.at(i, j) ) / test.at(i, j) ) >= 1.e-8 ) {
                OOFEM_ERROR("test failure (err=%e)", ( answer.at(i, j) - test.at(i, j) ) / test.at(i, j));
            }
        }
    }

#endif
}

/////////////////////////////////////////////////////////////////////////////////////////////
void
TR1_2D_SUPG2_AXI :: computeBCRhsTerm_MB(FloatArray &answer, TimeStep *tStep)
{
    answer.resize(6);
    answer.zero();

    int nLoads;
    FloatArray un, nV, gVector;
    double u, v;

    // add body load (gravity) termms
    this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), un);

    nLoads = this->giveBodyLoadArray()->giveSize();
    for ( int i = 1; i <= nLoads; i++ ) {
        auto load  = domain->giveLoad( bodyLoadArray.at(i) );
        bcGeomType ltype = load->giveBCGeoType();
        if ( ltype == BodyLoadBGT && load->giveBCValType() == ForceLoadBVT ) {
            load->computeComponentArrayAt(gVector, tStep, VM_Total);
            if ( gVector.giveSize() ) {
                for ( int ifluid = 0; ifluid < 2; ifluid++ ) {
                    for ( GaussPoint *gp: *integrationRulesArray [ ifluid ] ) {
                        double rho = this->_giveMaterial(ifluid)->give('d', gp);
                        double dV = this->computeVolumeAroundID(gp, id [ ifluid ], vcoords [ ifluid ]);
                        double coeff = rho * dV;
                        computeNVector(nV, gp);
                        u = nV.at(1) * un.at(1) + nV.at(2) * un.at(3) + nV.at(3) * un.at(5);
                        v = nV.at(1) * un.at(2) + nV.at(2) * un.at(4) + nV.at(3) * un.at(6);

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
    }

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

                auto bLoad = dynamic_cast< BoundaryLoad * >( domain->giveLoad(n) );
                if ( bLoad ) {
                    bLoad->computeValueAt(t, tStep, coords, VM_Total);

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
TR1_2D_SUPG2_AXI :: computeBCRhsTerm_MC(FloatArray &answer, TimeStep *tStep)
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
                for ( int ifluid = 0; ifluid < 2; ifluid++ ) {
                    for ( GaussPoint *gp: *integrationRulesArray [ ifluid ] ) {
                        double dV = this->computeVolumeAroundID(gp, id [ ifluid ], vcoords [ ifluid ]);
                        double coeff = t_pspg * dV;

                        answer.at(1) += coeff * ( b [ 0 ] * gVector.at(1) + c [ 0 ] * gVector.at(2) );
                        answer.at(2) += coeff * ( b [ 1 ] * gVector.at(1) + c [ 1 ] * gVector.at(2) );
                        answer.at(3) += coeff * ( b [ 2 ] * gVector.at(1) + c [ 2 ] * gVector.at(2) );
                    }
                }
            }
        }
    }
}



void
TR1_2D_SUPG2_AXI :: updateStabilizationCoeffs(TimeStep *tStep)
{
    //TR1_2D_SUPG :: updateStabilizationCoeffs (tStep);
#if 0
    int i, j, k, l, w_dof_addr, u_dof_addr, ip, ifluid;
    double __g_norm, __gamma_norm, __gammav_norm, __beta_norm, __betav_norm, __c_norm, __e_norm, __k_norm, __Re;
    double __t_p1, __t_p2, __t_p3, __t_pv1, __t_pv2, __t_pv3;
    double nu, nu0, nu1, usum, vsum, rho, dV, u1, u2;
    FloatArray u, un, a;

    // compute averaged viscosity based on rule of mixture
    GaussPoint *gp;
    if ( integrationRulesArray [ 0 ]->giveNumberOfIntegrationPoints() ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    } else {
        gp = integrationRulesArray [ 1 ]->getIntegrationPoint(0);
    }

    nu0 = this->_giveMaterial(0)->giveCharacteristicValue( MRM_Viscosity, gp, tStep->givePreviousStep() );
    nu1 = this->_giveMaterial(1)->giveCharacteristicValue( MRM_Viscosity, gp, tStep->givePreviousStep() );
    nu = vof * nu0 + ( 1. - vof ) * nu1;

    //this -> computeVectorOfVelocities(VM_Total,tStep->givePreviousStep(),un) ;
    this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), u);
    this->computeVectorOfVelocities(VM_Acceleration, tStep, a);
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
        for ( GaussPoint *gp: integrationRulesArray [ ifluid ] ) {
            rho = this->_giveMaterial(ifluid)->give('d', gp);
            this->computeNMtrx(n, gp);
            dV = this->computeVolumeAroundID(gp, id [ ifluid ], vcoords [ ifluid ]);

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
        for ( GaussPoint *gp: *integrationRulesArray [ ifluid ] ) {
            rho = this->_giveMaterial(ifluid)->give('d', gp);
            dV = this->computeVolumeAroundID(gp, id [ ifluid ], vcoords [ ifluid ]);

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
        for ( GaussPoint *gp: *integrationRulesArray [ ifluid ] ) {
            rho = this->_giveMaterial(ifluid)->give('d', gp);
            this->computeNMtrx(n, gp);
            dV = this->computeVolumeAroundID(gp, id [ ifluid ], vcoords [ ifluid ]);

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
        double t_sugn2 = tStep->giveTimeIncrement() / 2.0;
        //t_sugn3 = inf;
        this->t_supg = 1. / sqrt( 1. / ( t_sugn2 * t_sugn2 ) );
        this->t_pspg = this->t_supg;
        this->t_lsic = 0.0;
    } else {
        __Re = vnorm * vnorm * __c_norm / __k_norm / nu;

        __t_p1 = __g_norm / __gamma_norm;
        __t_p2 = tStep->giveTimeIncrement() * __g_norm / 2.0 / __beta_norm;
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
    int im1;
    FloatArray u;

    uscale = domain->giveEngngModel()->giveVariableScale(VST_Velocity);
    lscale = domain->giveEngngModel()->giveVariableScale(VST_Length);
    tscale = domain->giveEngngModel()->giveVariableScale(VST_Time);
    dscale = domain->giveEngngModel()->giveVariableScale(VST_Density);

    this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), u);
    u.times(uscale);
    double nu, nu0, nu1;

    // compute averaged viscosity based on rule of mixture
    GaussPoint *gp;
    if ( integrationRulesArray [ 0 ]->giveNumberOfIntegrationPoints() ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    } else {
        gp = integrationRulesArray [ 1 ]->getIntegrationPoint(0);
    }

    nu0 = static_cast< FluidDynamicMaterial * >( this->_giveMaterial(0) )->giveEffectiveViscosity( gp, tStep->givePreviousStep() );
    nu1 = static_cast< FluidDynamicMaterial * >( this->_giveMaterial(1) )->giveEffectiveViscosity( gp, tStep->givePreviousStep() );
    nu = vof * nu0 + ( 1. - vof ) * nu1;
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
TR1_2D_SUPG2_AXI :: computeCriticalTimeStep(TimeStep *tStep)
{
    return 1.e3;
}


void
TR1_2D_SUPG2_AXI :: computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    /* one computes here average deviatoric stress, based on rule of mixture (this is used only for postprocessing) */
    FloatArray u, eps, s0, s1;
    FloatMatrix _b;
    answer.resize(3);


    this->computeVectorOfVelocities(VM_Total, tStep, u);
    this->computeBMtrx(_b, gp);
    eps.beProductOf(_b, u);

    static_cast< FluidDynamicMaterial * >( this->_giveMaterial(0) )->computeDeviatoricStressVector(s0, gp, eps, tStep);
    static_cast< FluidDynamicMaterial * >( this->_giveMaterial(1) )->computeDeviatoricStressVector(s1, gp, eps, tStep);

    for ( int i = 1; i <= 3; i++ ) {
        answer.at(i) = ( temp_vof ) * s0.at(i) + ( 1. - temp_vof ) * s1.at(i);
    }
}


/*
 * double
 * TR1_2D_SUPG2 :: computeCriticalTimeStep (TimeStep* tStep)
 * {
 * FloatArray u;
 * double dt1, dt2, dt;
 * double Re = static_cast<FluidModel*>(domain->giveEngngModel())->giveReynoldsNumber();
 *
 * this -> computeVectorOfVelocities(VM_Total,tStep, u) ;
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
TR1_2D_SUPG2_AXI :: computeLEPLICVolumeFraction(const FloatArray &n, const double p, LEPlic *matInterface, bool updFlag)
{
    Polygon pg;
    double answer, volume = computeMyVolume(matInterface, updFlag);
    this->formVolumeInterfacePoly(pg, matInterface, n, p, updFlag);
    answer = fabs(pg.computeVolume() / volume);
    if ( answer > 1.000000001 ) {
        OOFEM_WARNING("VOF fraction out of bounds, vof = %e\n", answer);
        return 1.0;
    } else {
        return answer;
    }
}

void
TR1_2D_SUPG2_AXI :: formMaterialVolumePoly(Polygon &matvolpoly, LEPlic *matInterface,
                                           const FloatArray &normal, const double p, bool updFlag)
{
    double x, y;
    Vertex v;

    matvolpoly.clear();

    if ( this->vof <= TRSUPG_ZERO_VOF ) {
        return;
    } else if ( this->vof >= ( 1 - TRSUPG_ZERO_VOF ) ) {
        for ( int i = 1; i <= 3; i++ ) {
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
TR1_2D_SUPG2_AXI :: formVolumeInterfacePoly(Polygon &matvolpoly, LEPlic *matInterface,
                                            const FloatArray &normal, const double p, bool updFlag)
{
    int next;
    bool nodeIn [ 3 ];
    double nx = normal.at(1), ny = normal.at(2), x, y;
    double tx, ty;
    Vertex v;

    matvolpoly.clear();

    for ( int i = 1; i <= 3; i++ ) {
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
        for ( int i = 1; i <= 3; i++ ) {
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
        for ( int i = 1; i <= 3; i++ ) {
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
TR1_2D_SUPG2_AXI :: updateVolumePolygons(Polygon &referenceFluidPoly, Polygon &secondFluidPoly, int &rfPoints, int &sfPoints,
                                         const FloatArray &normal, const double p, bool updFlag)
{
    /*
     * this method updates two polygons, one filled with reference fluid and second filled with
     * other fluid (air). These two polygons are used in integrating element contributions.
     */
    int next;
    bool nodeIn [ 3 ];
    double nx = normal.at(1), ny = normal.at(2), x, y;
    double tx, ty;
    Vertex v;

    rfPoints = sfPoints = 0;
    referenceFluidPoly.clear();
    secondFluidPoly.clear();

    for ( int i = 1; i <= 3; i++ ) {
        x = this->giveNode(i)->giveCoordinate(1);
        y = this->giveNode(i)->giveCoordinate(2);

        if ( ( nx * x + ny * y + p ) >= 0. ) {
            nodeIn [ i - 1 ] = true;
        } else {
            nodeIn [ i - 1 ] = false;
        }
    }

    if ( nodeIn [ 0 ] && nodeIn [ 1 ] && nodeIn [ 2 ] ) { // all nodes inside
        for ( int i = 1; i <= 3; i++ ) {
            x = this->giveNode(i)->giveCoordinate(1);
            y = this->giveNode(i)->giveCoordinate(2);

            v.setCoords(x, y);
            referenceFluidPoly.addVertex(v);
            rfPoints++;
        }

        return;
    } else if ( !( nodeIn [ 0 ] || nodeIn [ 1 ] || nodeIn [ 2 ] ) ) { // all nodes outside
        for ( int i = 1; i <= 3; i++ ) {
            x = this->giveNode(i)->giveCoordinate(1);
            y = this->giveNode(i)->giveCoordinate(2);

            v.setCoords(x, y);
            secondFluidPoly.addVertex(v);
            sfPoints++;
        }

        return;
    } else {
        for ( int i = 1; i <= 3; i++ ) {
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
TR1_2D_SUPG2_AXI :: truncateMatVolume(const Polygon &matvolpoly, double &volume)
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
TR1_2D_SUPG2_AXI :: formMyVolumePoly(Polygon &me, LEPlic *matInterface, bool updFlag)
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
TR1_2D_SUPG2_AXI :: computeMyVolume(LEPlic *matInterface, bool updFlag)
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
TR1_2D_SUPG2_AXI :: giveElementCenter(LEPlic *mat_interface, FloatArray &center, bool upd)
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

void
TR1_2D_SUPG2_AXI :: updateIntegrationRules()
{
    int c [ 2 ];
    int nip;
    double x, y;
    Vertex v;

    for ( int i = 0; i < 2; i++ ) {
        myPoly [ i ].clear();
    }

    if ( this->temp_vof <= TRSUPG_ZERO_VOF ) {
        for ( int i = 1; i <= 3; i++ ) {
            x = this->giveNode(i)->giveCoordinate(1);
            y = this->giveNode(i)->giveCoordinate(2);
            v.setCoords(x, y);
            myPoly [ 1 ].addVertex(v);
            c [ 1 ] = 3;
            c [ 0 ] = 0;
        }
    } else if ( this->temp_vof >= ( 1. - TRSUPG_ZERO_VOF ) ) {
        for ( int i = 1; i <= 3; i++ ) {
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
    FEI2dTrLin triaApprox(1, 2);
    FEI2dQuadLin quadApprox(1, 2);
    FEInterpolation *approx = NULL;
    // set up integration points
    for ( int i = 0; i < 2; i++ ) {
        if ( c [ i ] == 3 ) {
            id [ i ] = _Triangle;
            approx = & triaApprox;
        } else if ( c [ i ] == 4 ) {
            id [ i ] = _Square;
            approx = & quadApprox;
        } else if ( c [ i ] == 0 ) {
            continue;
        } else {
            OOFEM_ERROR("cannot set up integration domain for %d vertex polygon", c [ i ]);
        }

        vcoords [ i ].clear();
        vcoords [ i ].reserve( c [ i ] );

        // set up vertex coords
        Polygon :: PolygonVertexIterator it(myPoly + i);
        while ( it.giveNext(& p) ) {
            vcoords [ i ].push_back( *p->getCoords() );
        }

        if ( id [ i ] == _Triangle ) {
            nip = 4;
        } else {
            nip = 4;
        }

        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ i ], nip, this);

        // remap ip coords into area coords of receiver
        for ( GaussPoint *gp: *integrationRulesArray [ i ] ) {
            approx->local2global( gc, gp->giveNaturalCoordinates(), FEIVertexListGeometryWrapper(vcoords [ i ]) );
            triaApprox.global2local( lc, gc, FEIElementGeometryWrapper(this) );
            // modify original ip coords to target ones
            gp->setSubPatchCoordinates( gp->giveNaturalCoordinates() );
            gp->setNaturalCoordinates(lc);
            //gp->setWeight (gp->giveWeight()*a/area);
        }
    }

    // internal test -> compute receiver area
    double dV, __area = 0.0;
    for ( int ifluid = 0; ifluid < 2; ifluid++ ) {
        for ( GaussPoint *gp: *integrationRulesArray [ ifluid ] ) {
            dV = this->computeVolumeAroundID(gp, id [ ifluid ], vcoords [ ifluid ]);
            // compute integral here
            __area += dV;
        }
    }

    /*
     * double __err = fabs(__area-area)/area;
     * if (__err > 1.e-6) {
     * OOFEM_WARNING("volume inconsistency (%5.2f\%)", __err*100);
     *
     * __area=0.0;
     * for (ifluid = 0; ifluid< 2; ifluid++) {
     *  for (GaussPoint *gp: *integrationRulesArray[ifluid]) {
     *    dV = this->computeVolumeAroundID(gp,id[ifluid], vcoords[ifluid]) ;
     *    // compute integral here
     *    __area += dV;
     *  }
     * }
     * }
     */
}


double
TR1_2D_SUPG2_AXI :: computeRadiusAt(GaussPoint *gp)
{
    double r1, r2, r3;
    double n1, n2, n3;

    r1 = this->giveNode(1)->giveCoordinate(1);
    r2 = this->giveNode(2)->giveCoordinate(1);
    r3 = this->giveNode(3)->giveCoordinate(1);

    n1 = gp->giveNaturalCoordinate(1);
    n2 = gp->giveNaturalCoordinate(2);
    n3 = 1. - n1 - n2;

    return n1 * r1 + n2 * r2 + n3 * r3;
}



double
TR1_2D_SUPG2_AXI :: computeVolumeAroundID(GaussPoint *gp, integrationDomain id, const std::vector< FloatArray > &idpoly)
{
    double weight = gp->giveWeight();
    double _r = computeRadiusAt(gp);

    if ( id == _Triangle ) {
        FEI2dTrLin __interpolation(1, 2);
        return _r *weight *fabs( __interpolation.giveTransformationJacobian ( gp->giveSubPatchCoordinates(), FEIVertexListGeometryWrapper(idpoly) ) );
    } else {
        FEI2dQuadLin __interpolation(1, 2);
        double det = fabs( __interpolation.giveTransformationJacobian( gp->giveSubPatchCoordinates(), FEIVertexListGeometryWrapper(idpoly) ) );
        return _r * det * weight;
    }
}

void TR1_2D_SUPG2_AXI :: computeBMtrx(FloatMatrix &_b, GaussPoint *gp)
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

/*
 * double
 * TR1_2D_SUPG2::computeVolumeAroundID(GaussPoint* gp, integrationDomain id, const Polygon& matvolpoly) {
 * }
 */
void
TR1_2D_SUPG2_AXI :: computeNVector(FloatArray &n, GaussPoint *gp)
{
    this->interp.evalN( n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
}

void
TR1_2D_SUPG2_AXI :: printOutputAt(FILE *file, TimeStep *tStep)
// Performs end-of-step operations.
{
    SUPGElement :: printOutputAt(file, tStep);

    GaussPoint *gp;
    if ( integrationRulesArray [ 0 ]->giveNumberOfIntegrationPoints() ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    } else {
        gp = integrationRulesArray [ 1 ]->getIntegrationPoint(0);
    }

    double rho = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity(gp);
    fprintf(file, "VOF %e, density %e\n\n", this->giveVolumeFraction(), rho);
}


void
TR1_2D_SUPG2_AXI :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                               InternalStateType type, TimeStep *tStep)
{
    GaussPoint *gp;
    if ( integrationRulesArray [ 0 ]->giveNumberOfIntegrationPoints() ) {
        gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    } else {
        gp = integrationRulesArray [ 1 ]->getIntegrationPoint(0);
    }

    this->giveIPValue(answer, gp, type, tStep);
}

void
TR1_2D_SUPG2_AXI :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(3);
    pap.at(1) = this->giveNode(1)->giveNumber();
    pap.at(2) = this->giveNode(2)->giveNumber();
    pap.at(3) = this->giveNode(3)->giveNumber();
}

void
TR1_2D_SUPG2_AXI :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
{
    answer.resize(1);
    if ( ( pap == this->giveNode(1)->giveNumber() ) ||
        ( pap == this->giveNode(2)->giveNumber() ) ||
        ( pap == this->giveNode(3)->giveNumber() ) ) {
        answer.at(1) = pap;
    } else {
        OOFEM_ERROR("node unknown");
    }
}

int
TR1_2D_SUPG2_AXI :: SPRNodalRecoveryMI_giveNumberOfIP()
{ return 1; }


SPRPatchType
TR1_2D_SUPG2_AXI :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_2dxy;
}



#ifdef __OOFEG
int
TR1_2D_SUPG2_AXI :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                            int node, TimeStep *tStep)
{
    /*
     * if (type == IST_VOFFraction) {
     * answer.resize(1);
     * answer.at(1) = this->giveTempVolumeFraction();
     * return 1;
     * } else if (type == IST_Density) {
     * answer.resize(1);
     * answer.at(1) = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveDensity(gp);
     * return 1;
     *
     * } else
     */return SUPGElement :: giveInternalStateAtNode(answer, type, mode, node, tStep);
}

void
TR1_2D_SUPG2_AXI :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
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

void TR1_2D_SUPG2_AXI :: drawScalar(oofegGraphicContext &gc, TimeStep *tStep)
{
    int i, indx, result = 0;
    WCRec p [ 3 ];
    GraphicObj *tr;
    FloatArray v1, v2, v3;
    double s [ 3 ];

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLayer(OOFEG_VARPLOT_PATTERN_LAYER);

    // if ((gc.giveIntVarMode() == ISM_local) && (gc.giveIntVarType() ==  IST_VOFFraction)) {
    if ( ( gc.giveIntVarType() ==  IST_VOFFraction ) && ( gc.giveIntVarMode() == ISM_local ) ) {
        Polygon matvolpoly;
        this->formMaterialVolumePoly(matvolpoly, NULL, temp_normal, temp_p, false);
        EASValsSetColor( gc.getStandardSparseProfileColor() );
        //GraphicObj *go = matvolpoly.draw(gc,true,OOFEG_VARPLOT_PATTERN_LAYER);
        matvolpoly.draw(gc, true, OOFEG_VARPLOT_PATTERN_LAYER);
        return;
    }

    /*
     * if ((gc.giveIntVarMode() == ISM_local) && (gc.giveIntVarType() ==  IST_VOFFraction)) {
     * Polygon matvolpoly;
     * FloatArray _n;
     * double _p;
     * this->doCellDLS (_n, this->giveNumber());
     * this->findCellLineConstant (_p, _n, this->giveNumber());
     * this->formMaterialVolumePoly(matvolpoly, NULL, _n, _p, false);
     * GraphicObj *go = matvolpoly.draw(gc,true,OOFEG_VARPLOT_PATTERN_LAYER);
     * return;
     * }
     */

    if ( gc.giveIntVarMode() == ISM_recovered ) {
        result += this->giveInternalStateAtNode(v1, gc.giveIntVarType(), gc.giveIntVarMode(), 1, tStep);
        result += this->giveInternalStateAtNode(v2, gc.giveIntVarType(), gc.giveIntVarMode(), 2, tStep);
        result += this->giveInternalStateAtNode(v3, gc.giveIntVarType(), gc.giveIntVarMode(), 3, tStep);
    } else if ( gc.giveIntVarMode() == ISM_local ) {
        GaussPoint *gp;
        if ( integrationRulesArray [ 0 ]->giveNumberOfIntegrationPoints() ) {
            gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
        } else {
            gp = integrationRulesArray [ 1 ]->getIntegrationPoint(0);
        }

        result += giveIPValue(v1, gp, gc.giveIntVarType(), tStep);
        v2 = v1;
        v3 = v1;
        result *= 3;
    }

    if ( result != 3 ) {
        return;
    }

    indx = gc.giveIntVarIndx();

    s [ 0 ] = v1.at(indx);
    s [ 1 ] = v2.at(indx);
    s [ 2 ] = v3.at(indx);

    if ( gc.getScalarAlgo() == SA_ISO_SURF ) {
        for ( i = 0; i < 3; i++ ) {
            p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
            p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
            p [ i ].z = 0.;
        }

        //EASValsSetColor(gc.getYieldPlotColor(ratio));
        gc.updateFringeTableMinMax(s, 3);
        tr =  CreateTriangleWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ]);
        EGWithMaskChangeAttributes(LAYER_MASK, tr);
        EMAddGraphicsToModel(ESIModel(), tr);
    } else if ( ( gc.getScalarAlgo() == SA_ZPROFILE ) || ( gc.getScalarAlgo() == SA_COLORZPROFILE ) ) {
        double landScale = gc.getLandScale();

        for ( i = 0; i < 3; i++ ) {
            p [ i ].x = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(1);
            p [ i ].y = ( FPNum ) this->giveNode(i + 1)->giveCoordinate(2);
            p [ i ].z = s [ i ] * landScale;
        }

        if ( gc.getScalarAlgo() == SA_ZPROFILE ) {
            EASValsSetColor( gc.getDeformedElementColor() );
            EASValsSetLineWidth(OOFEG_DEFORMED_GEOMETRY_WIDTH);
            EASValsSetFillStyle(FILL_SOLID);
            tr =  CreateTriangle3D(p);
            EGWithMaskChangeAttributes(WIDTH_MASK | COLOR_MASK | FILL_MASK | LAYER_MASK, tr);
        } else {
            gc.updateFringeTableMinMax(s, 3);
            EASValsSetFillStyle(FILL_SOLID);
            tr =  CreateTriangleWD3D(p, s [ 0 ], s [ 1 ], s [ 2 ]);
            EGWithMaskChangeAttributes(FILL_MASK | LAYER_MASK, tr);
        }

        EMAddGraphicsToModel(ESIModel(), tr);
    }
}



#endif
} // end namespace oofem
