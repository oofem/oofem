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

#include "tet1_3d_supg.h"
#include "fei3dtetlin.h"
#include "fluidmodel.h"
#include "node.h"
#include "material.h"
#include "gausspoint.h"
#include "gaussintegrationrule.h"
#include "fluiddynamicmaterial.h"
#include "fluidcrosssection.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "engngm.h"
#include "timestep.h"
#include "crosssection.h"
#include "classfactory.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
REGISTER_Element(Tet1_3D_SUPG);

FEI3dTetLin Tet1_3D_SUPG :: interpolation;

Tet1_3D_SUPG :: Tet1_3D_SUPG(int n, Domain *aDomain) :
    SUPGElement2(n, aDomain)
    // Constructor.
{
    numberOfDofMans  = 4;
}

Tet1_3D_SUPG :: ~Tet1_3D_SUPG()
// Destructor
{ }

int
Tet1_3D_SUPG :: computeNumberOfDofs()
{
    return 16;
}

void
Tet1_3D_SUPG :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    answer = {V_u, V_v, V_w, P_f};
}


void
Tet1_3D_SUPG :: computeGaussPoints()
// Sets up the array containing the integration points of the receiver.
{
    if ( integrationRulesArray.size() == 0 ) {
        integrationRulesArray.resize(3);
        integrationRulesArray [ 0 ].reset( new GaussIntegrationRule(1, this, 1, 3) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 0 ], 1, this);

        integrationRulesArray [ 1 ].reset( new GaussIntegrationRule(2, this, 1, 3) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 1 ], 4, this);

        integrationRulesArray [ 2 ].reset( new GaussIntegrationRule(3, this, 1, 3) );
        this->giveCrossSection()->setupIntegrationPoints(* integrationRulesArray [ 2 ], 4, this);
    }
}


Interface *
Tet1_3D_SUPG :: giveInterface(InterfaceType interface)
{
    if ( interface == LevelSetPCSElementInterfaceType ) {
        return static_cast< LevelSetPCSElementInterface * >(this);
    }

    return NULL;
}


void
Tet1_3D_SUPG :: computeNuMatrix(FloatMatrix &answer, GaussPoint *gp)
{
    FloatArray n;
    this->interpolation.evalN( n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    answer.beNMatrixOf(n, 3);
}

void
Tet1_3D_SUPG :: computeUDotGradUMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatMatrix n, dn;
    FloatArray u, un;
    interpolation.evaldNdx( dn, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    this->computeNuMatrix(n, gp);
    this->computeVectorOfVelocities(VM_Total, tStep, un);

    u.beProductOf(n, un);

    answer.resize(3, 12);
    answer.zero();
    for ( int i = 1; i <= 4; i++ ) {
        answer.at(1, 3 * i - 2) = u.at(1) * dn.at(i, 1) + u.at(2) * dn.at(i, 2) + u.at(3) * dn.at(i, 3);
        answer.at(2, 3 * i - 1) = u.at(1) * dn.at(i, 1) + u.at(2) * dn.at(i, 2) + u.at(3) * dn.at(i, 3);
        answer.at(3, 3 * i - 0) = u.at(1) * dn.at(i, 1) + u.at(2) * dn.at(i, 2) + u.at(3) * dn.at(i, 3);
    }
}

void
Tet1_3D_SUPG :: computeDivTauMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep)
{
    answer.resize(3, 12);
    answer.zero();
}


void
Tet1_3D_SUPG :: computeGradUMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray u;
    FloatMatrix dn, um(3, 4);

    this->computeVectorOfVelocities(VM_Total, tStep, u);

    interpolation.evaldNdx( dn, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
    for ( int i = 1; i <= 4; i++ ) {
        um.at(1, i) = u.at(3 * i - 2);
        um.at(2, i) = u.at(3 * i - 1);
        um.at(3, i) = u.at(3 * i - 0);
    }

    answer.beProductOf(um, dn);
}


void
Tet1_3D_SUPG :: computeBMatrix(FloatMatrix &answer, GaussPoint *gp)
{
    FloatMatrix dn;
    interpolation.evaldNdx( dn, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(6, 12);
    answer.zero();

    for ( int i = 1; i <= 4; i++ ) {
        answer.at(1, 3 * i - 2) = dn.at(i, 1);
        answer.at(2, 3 * i - 1) = dn.at(i, 2);
        answer.at(3, 3 * i - 0) = dn.at(i, 3);

        answer.at(4, 3 * i - 1) = dn.at(i, 3);
        answer.at(4, 3 * i - 0) = dn.at(i, 2);

        answer.at(5, 3 * i - 2) = dn.at(i, 3);
        answer.at(5, 3 * i - 0) = dn.at(i, 1);

        answer.at(6, 3 * i - 2) = dn.at(i, 2);
        answer.at(6, 3 * i - 1) = dn.at(i, 1);
    }
}

void
Tet1_3D_SUPG :: computeDivUMatrix(FloatMatrix &answer, GaussPoint *gp)
{
    FloatMatrix dn;
    interpolation.evaldNdx( dn, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(1, 12);
    answer.zero();

    for ( int i = 1; i <= 4; i++ ) {
        answer.at(1, 3 * i - 2) = dn.at(i, 1);
        answer.at(1, 3 * i - 1) = dn.at(i, 2);
        answer.at(1, 3 * i - 0) = dn.at(i, 3);
    }
}

void
Tet1_3D_SUPG :: computeNpMatrix(FloatMatrix &answer, GaussPoint *gp)
{
    FloatArray n;
    this->interpolation.evalN( n, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.resize(1, 4);
    answer.zero();

    answer.at(1, 1)  = n.at(1);
    answer.at(1, 2)  = n.at(2);
    answer.at(1, 3)  = n.at(3);
    answer.at(1, 4)  = n.at(4);
}


void
Tet1_3D_SUPG :: computeGradPMatrix(FloatMatrix &answer, GaussPoint *gp)
{
    FloatMatrix dn;
    interpolation.evaldNdx( dn, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );

    answer.beTranspositionOf(dn);
}



void
Tet1_3D_SUPG :: updateStabilizationCoeffs(TimeStep *tStep)
{
    //TR1_2D_SUPG :: updateStabilizationCoeffs (tStep);
    /* UGN-Based Stabilization */
    double h_ugn, sum = 0.0, vnorm, t_sugn1, t_sugn2, t_sugn3, u_1, u_2, u_3, z, Re_ugn;
    double dscale, uscale, lscale, tscale, dt;
    //bool zeroFlag = false;
    int im1;
    FloatArray u, divu;
    FloatMatrix du;

    uscale = domain->giveEngngModel()->giveVariableScale(VST_Velocity);
    lscale = domain->giveEngngModel()->giveVariableScale(VST_Length);
    tscale = domain->giveEngngModel()->giveVariableScale(VST_Time);
    dscale = domain->giveEngngModel()->giveVariableScale(VST_Density);

    this->computeVectorOfVelocities(VM_Total, tStep->givePreviousStep(), u);

    u.times(uscale);
    double nu;

    // compute averaged viscosity based on rule of mixture

    dt = tStep->giveTimeIncrement() * tscale;

    std :: unique_ptr< IntegrationRule > &iRule = this->integrationRulesArray [ 1 ];
    GaussPoint *gp_first = iRule->getIntegrationPoint(0);
    nu = static_cast< FluidCrossSection * >( this->giveCrossSection() )->giveFluidMaterial()->giveEffectiveViscosity(gp_first, tStep->givePreviousStep() );
    nu *= domain->giveEngngModel()->giveVariableScale(VST_Viscosity);

    for ( auto *gp: *iRule ) {
        this->computeDivUMatrix(du, gp);
        divu.beProductOf(du, u);
        sum += divu.at(1);
    }

    sum *= ( 1. / lscale / iRule->giveNumberOfIntegrationPoints() );

    /*
     * for (i=1; i<=3;i++) {
     * im1=i-1;
     * sum+= fabs(u.at((im1)*2+1)*b[im1]/lscale + u.at(im1*2+2)*c[im1]/lscale);
     * }
     */
    vnorm = 0.;
    int nsd = this->giveNumberOfSpatialDimensions();
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        im1 = i - 1;
        u_1 = u.at( ( im1 ) * nsd + 1 );
        u_2 = u.at( ( im1 ) * nsd + 2 );
        if ( nsd > 2 ) {
            u_3 = u.at( ( im1 ) * nsd + 3 );
        } else {
            u_3 = 0.;
        }

        vnorm = max( vnorm, sqrt(u_1 * u_1 + u_2 * u_2 + u_3 * u_3) );
    }

    if ( ( vnorm == 0.0 ) || ( sum <  vnorm * 1e-10  ) ) {
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

    //this->t_lsic=0.0;
    //this->t_pspg=0.0;
}

int
Tet1_3D_SUPG :: giveNumberOfSpatialDimensions()
{
    return 3;
}

double
Tet1_3D_SUPG :: computeCriticalTimeStep(TimeStep *tStep)
{
    FloatArray u;
    double Re = static_cast< FluidModel * >( domain->giveEngngModel() )->giveReynoldsNumber();

    this->computeVectorOfVelocities(VM_Total, tStep, u);

    double vn1 = sqrt( u.at(1) * u.at(1) + u.at(2) * u.at(2) + u.at(3) * u.at(3) );
    double vn2 = sqrt( u.at(4) * u.at(4) + u.at(5) * u.at(5) + u.at(6) * u.at(6) );
    double vn3 = sqrt( u.at(7) * u.at(7) + u.at(8) * u.at(8) + u.at(9) * u.at(9) );
    double vn4 = sqrt( u.at(10) * u.at(10) + u.at(11) * u.at(11) + u.at(12) * u.at(12) );
    double veln = max( vn1, max( vn2, max(vn3, vn4) ) );

    double ln = 1.e6;
    Node *inode, *jnode, *knode, *lnode;
    FloatArray t1, t2, n3, n;
    for ( int l = 1; l <= 4; l++ ) {
        int i = ( l > 3 ) ? 1 : l + 1;
        int j = ( i > 3 ) ? 1 : i + 1;
        int k = ( j > 3 ) ? 1 : j + 1;

        inode = this->giveNode(i);
        jnode = this->giveNode(j);
        knode = this->giveNode(k);
        lnode = this->giveNode(l);
        t1.beDifferenceOf(*inode->giveCoordinates(), *jnode->giveCoordinates());

        t2.beDifferenceOf(*knode->giveCoordinates(), *jnode->giveCoordinates());

        n.beVectorProductOf(t1, t2);
        n.normalize();

        n3.beDifferenceOf(*lnode->giveCoordinates(), *jnode->giveCoordinates());

        ln = min( ln, sqrt( fabs( n.dotProduct(n3) ) ) );
    }

    if ( veln != 0.0 ) {
        return ln / veln;
    } else {
        return 0.5 * ln * ln * Re;
    }
}


double
Tet1_3D_SUPG :: computeVolumeAround(GaussPoint *gp)
// Returns the portion of the receiver which is attached to gp.
{
    double determinant, weight, volume;
    determinant = fabs( this->interpolation.giveTransformationJacobian( gp->giveNaturalCoordinates(),
                                                                       FEIElementGeometryWrapper(this) ) );

    weight = gp->giveWeight();
    volume = determinant * weight;

    return volume;
}


double
Tet1_3D_SUPG :: LS_PCS_computeF(LevelSetPCS *ls, TimeStep *tStep)
{
    double answer = 0.0, norm, dV, vol = 0.0;
    FloatMatrix n, dn;
    FloatArray fi(4), u, un, gfi;

    this->computeVectorOfVelocities(VM_Total, tStep, un);

    for ( int i = 1; i <= 4; i++ ) {
        fi.at(i) = ls->giveLevelSetDofManValue( dofManArray.at(i) );
    }

    for ( GaussPoint *gp: *this->integrationRulesArray [ 0 ] ) {
        dV  = this->computeVolumeAround(gp);
        interpolation.evaldNdx( dn, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
        this->computeNuMatrix(n, gp);
        u.beProductOf(n, un);
        gfi.beTProductOf(dn, fi);
        norm = gfi.computeNorm();

        vol += dV;
        answer += dV * u.dotProduct(gfi) / norm;
    }

    return answer / vol;
}

void
Tet1_3D_SUPG :: LS_PCS_computedN(FloatMatrix &answer)
{
    GaussPoint *gp = this->integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    interpolation.evaldNdx( answer, gp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this) );
}


double
Tet1_3D_SUPG :: LS_PCS_computeVolume()
{
    double answer = 0.0;
 
    for ( GaussPoint *gp: *this->integrationRulesArray [ 0 ] ) {
        answer += this->computeVolumeAround(gp);
    }

    return answer;
}

double
Tet1_3D_SUPG :: LS_PCS_computeS(LevelSetPCS *ls, TimeStep *tStep)
{
    FloatArray voff(2), fi(4);
    for ( int i = 1; i <= 4; i++ ) {
        fi.at(i) = ls->giveLevelSetDofManValue( dofManArray.at(i) );
    }

    this->LS_PCS_computeVOFFractions(voff, fi);
    return ( voff.at(1) - voff.at(2) );
}




void
Tet1_3D_SUPG :: LS_PCS_computeVOFFractions(FloatArray &answer, FloatArray &fi)
{
    int neg = 0, pos = 0, zero = 0, si = 0;
    double x1, y1, z1;
    answer.resize(2);

    for ( double ifi: fi ) {
        if ( ifi >= 0. ) {
            pos++;
        } else if ( ifi < 0. ) {
            neg++;
        } else {
            zero++;
        }
    }

    if ( neg == 0 ) { // all level set values positive
        answer.at(1) = 1.0;
        answer.at(2) = 0.0;
        return; //return area;
    } else if ( pos == 0 ) { // all level set values negative
        answer.at(1) = 0.0;
        answer.at(2) = 1.0;
        return; //return -area;
    } else if ( zero == 4 ) {
        // ???????
        answer.at(1) = 1.0;
        answer.at(2) = 0.0;
        return;
    } else {
        // zero level set inside
        // distinguish two poosible cases
        if ( max(pos, neg) == 3 ) {
            // find the vertex vith level set sign different from other three
            for ( int i = 1; i <= 4; i++ ) {
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
                z1 = this->giveNode(si)->giveCoordinate(3);


                int ii;
                double t, xi [ 3 ], yi [ 3 ], zi [ 3 ];
                // compute intersections with element sides originating from this vertex
                for ( int i = 0; i < 3; i++ ) {
                    ii = ( si + i ) % 4 + 1;
                    t = fi.at(si) / ( fi.at(si) - fi.at(ii) );
                    xi [ i ] = x1 + t * ( this->giveNode(ii)->giveCoordinate(1) - x1 );
                    yi [ i ] = y1 + t * ( this->giveNode(ii)->giveCoordinate(2) - y1 );
                    zi [ i ] = z1 + t * ( this->giveNode(ii)->giveCoordinate(3) - z1 );
                }

                // compute volume of this pyramid (si, xi[0], xi[1], xi[2])
                double __vol = fabs( ( 1. / 6. ) * ( ( x1 - xi [ 0 ] ) * ( ( yi [ 1 ] - yi [ 0 ] ) * ( zi [ 2 ] - zi [ 0 ] ) - ( zi [ 1 ] - zi [ 0 ] ) * ( yi [ 2 ] - yi [ 0 ] ) ) +
                                                    ( y1 - yi [ 0 ] ) * ( ( zi [ 1 ] - zi [ 0 ] ) * ( xi [ 2 ] - xi [ 0 ] ) - ( xi [ 1 ] - xi [ 0 ] ) * ( zi [ 2 ] - zi [ 0 ] ) ) +
                                                    ( z1 - zi [ 0 ] ) * ( ( xi [ 1 ] - xi [ 0 ] ) * ( yi [ 2 ] - yi [ 0 ] ) - ( yi [ 1 ] - yi [ 0 ] ) * ( xi [ 2 ] - xi [ 0 ] ) ) ) );


                double vol = LS_PCS_computeVolume();
                if ( ( fabs(__vol) - vol ) < 0.0000001 ) {
                    __vol = sgn(__vol) * vol;
                }

                if ( ( __vol < 0 ) || ( fabs(__vol) / vol > 1.0000001 ) ) {
                    OOFEM_ERROR("internal consistency error");
                }

                if ( pos > neg ) {
                    // negative vol computed
                    answer.at(2) = min(fabs(__vol) / vol, 1.0);
                    answer.at(1) = 1.0 - answer.at(2);
                } else {
                    // postive vol computed
                    answer.at(1) = min(fabs(__vol) / vol, 1.0);
                    answer.at(2) = 1.0 - answer.at(1);
                }
            } else {
                OOFEM_ERROR("internal consistency error");
            }
        } else if ( max(pos, neg) == 2 ) {
            // two vertices positive; two negative; compute positive volume
            int p1 = 0, p2 = 0;
            // find the vertex vith level set sign different from other three
            for ( int i = 1; i <= 4; i++ ) {
                if ( fi.at(i) >= 0.0 ) {
                    if ( p1 ) {
                        p2 = i;
                        break;
                    } else {
                        p1 = i;
                    }
                }
            }

            int _ind, ii;
            double t;
            double p1i_x [ 3 ], p1i_y [ 3 ], p1i_z [ 3 ];
            double p2i_x [ 3 ], p2i_y [ 3 ], p2i_z [ 3 ];

            if ( p1 && p2 ) {
                // find the two intersections sharing edge with p1 and p2
                p1i_x [ 0 ] = this->giveNode(p1)->giveCoordinate(1);
                p1i_y [ 0 ] = this->giveNode(p1)->giveCoordinate(2);
                p1i_z [ 0 ] = this->giveNode(p1)->giveCoordinate(3);

                p2i_x [ 0 ] = this->giveNode(p2)->giveCoordinate(1);
                p2i_y [ 0 ] = this->giveNode(p2)->giveCoordinate(2);
                p2i_z [ 0 ] = this->giveNode(p2)->giveCoordinate(3);

                _ind = 1;
                for ( int i = 0; i < 3; i++ ) {
                    ii = ( p1 + i ) % 4 + 1;
                    if ( ( ii == p2 ) || ( ii == p1 ) ) {
                        continue;
                    }

                    t = fi.at(p1) / ( fi.at(p1) - fi.at(ii) );
                    p1i_x [ _ind ] = p1i_x [ 0 ] + t * ( this->giveNode(ii)->giveCoordinate(1) - p1i_x [ 0 ] );
                    p1i_y [ _ind ] = p1i_y [ 0 ] + t * ( this->giveNode(ii)->giveCoordinate(2) - p1i_y [ 0 ] );
                    p1i_z [ _ind ] = p1i_z [ 0 ] + t * ( this->giveNode(ii)->giveCoordinate(3) - p1i_z [ 0 ] );

                    t = fi.at(p2) / ( fi.at(p2) - fi.at(ii) );
                    p2i_x [ _ind ] =   p2i_x [ 0 ] + t * ( this->giveNode(ii)->giveCoordinate(1) - p2i_x [ 0 ] );
                    p2i_y [ _ind ] =   p2i_y [ 0 ] + t * ( this->giveNode(ii)->giveCoordinate(2) - p2i_y [ 0 ] );
                    p2i_z [ _ind++ ] = p2i_z [ 0 ] + t * ( this->giveNode(ii)->giveCoordinate(3) - p2i_z [ 0 ] );
                }

                // compute volume of this wedge as a sum of volumes of three
                // pyramids
                double __v1 = ( ( p2i_x [ 0 ] - p1i_x [ 0 ] ) * ( p1i_y [ 1 ] - p1i_y [ 0 ] ) * ( p1i_z [ 2 ] - p1i_z [ 0 ] ) -
                               ( p2i_x [ 0 ] - p1i_x [ 0 ] ) * ( p1i_y [ 2 ] - p1i_y [ 0 ] ) * ( p1i_z [ 1 ] - p1i_z [ 0 ] ) +
                               ( p1i_x [ 2 ] - p1i_x [ 0 ] ) * ( p2i_y [ 0 ] - p1i_y [ 0 ] ) * ( p1i_z [ 1 ] - p1i_z [ 0 ] ) -
                               ( p1i_x [ 1 ] - p1i_x [ 0 ] ) * ( p2i_y [ 0 ] - p1i_y [ 0 ] ) * ( p1i_z [ 2 ] - p1i_z [ 0 ] ) +
                               ( p1i_x [ 1 ] - p1i_x [ 0 ] ) * ( p1i_y [ 2 ] - p1i_y [ 0 ] ) * ( p2i_z [ 0 ] - p1i_z [ 0 ] ) -
                               ( p1i_x [ 2 ] - p1i_x [ 0 ] ) * ( p1i_y [ 1 ] - p1i_y [ 0 ] ) * ( p2i_z [ 0 ] - p1i_z [ 0 ] ) ) / 6.0;

                double __v2 = ( ( p2i_x [ 0 ] - p1i_x [ 1 ] ) * ( p1i_y [ 2 ] - p1i_y [ 1 ] ) * ( p2i_z [ 1 ] - p1i_z [ 1 ] ) -
                               ( p2i_x [ 0 ] - p1i_x [ 1 ] ) * ( p2i_y [ 1 ] - p1i_y [ 1 ] ) * ( p1i_z [ 2 ] - p1i_z [ 1 ] ) +
                               ( p2i_x [ 1 ] - p1i_x [ 1 ] ) * ( p2i_y [ 0 ] - p1i_y [ 1 ] ) * ( p1i_z [ 2 ] - p1i_z [ 1 ] ) -
                               ( p1i_x [ 2 ] - p1i_x [ 1 ] ) * ( p2i_y [ 0 ] - p1i_y [ 1 ] ) * ( p2i_z [ 1 ] - p1i_z [ 1 ] ) +
                               ( p1i_x [ 2 ] - p1i_x [ 1 ] ) * ( p2i_y [ 1 ] - p1i_y [ 1 ] ) * ( p2i_z [ 0 ] - p1i_z [ 1 ] ) -
                               ( p2i_x [ 1 ] - p1i_x [ 1 ] ) * ( p1i_y [ 2 ] - p1i_y [ 1 ] ) * ( p2i_z [ 0 ] - p1i_z [ 1 ] ) ) / 6.0;

                double __v3 = ( ( p1i_x [ 2 ] - p2i_x [ 0 ] ) * ( p2i_y [ 1 ] - p2i_y [ 0 ] ) * ( p2i_z [ 2 ] - p2i_z [ 0 ] ) -
                               ( p1i_x [ 2 ] - p2i_x [ 0 ] ) * ( p2i_y [ 2 ] - p2i_y [ 0 ] ) * ( p2i_z [ 1 ] - p2i_z [ 0 ] ) +
                               ( p2i_x [ 2 ] - p2i_x [ 0 ] ) * ( p1i_y [ 2 ] - p2i_y [ 0 ] ) * ( p2i_z [ 1 ] - p2i_z [ 0 ] ) -
                               ( p2i_x [ 1 ] - p2i_x [ 0 ] ) * ( p1i_y [ 2 ] - p2i_y [ 0 ] ) * ( p2i_z [ 2 ] - p2i_z [ 0 ] ) +
                               ( p2i_x [ 1 ] - p2i_x [ 0 ] ) * ( p2i_y [ 2 ] - p2i_y [ 0 ] ) * ( p1i_z [ 2 ] - p2i_z [ 0 ] ) -
                               ( p2i_x [ 2 ] - p2i_x [ 0 ] ) * ( p2i_y [ 1 ] - p2i_y [ 0 ] ) * ( p1i_z [ 2 ] - p2i_z [ 0 ] ) ) / 6.0;

                double __vol = fabs(__v1) + fabs(__v2) + fabs(__v3);
                double vol = LS_PCS_computeVolume();

                if ( ( __vol < 0 ) || ( fabs(__vol) / vol > 1.0000001 ) ) {
                    OOFEM_ERROR("internal consistency error");
                }

                answer.at(1) = min(fabs(__vol) / vol, 1.0);
                answer.at(2) = 1.0 - answer.at(1);
            } else {
                OOFEM_ERROR("internal consistency error");
            }
        } else {
            OOFEM_ERROR("internal consistency error");
        }
    }
}


#ifdef __OOFEG

void Tet1_3D_SUPG :: drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep)
{
    WCRec p [ 4 ];
    GraphicObj *go;

    if ( !gc.testElementGraphicActivity(this) ) {
        return;
    }

    EASValsSetLineWidth(OOFEG_RAW_GEOMETRY_WIDTH);
    EASValsSetColor( gc.getElementColor() );
    EASValsSetEdgeColor( gc.getElementEdgeColor() );
    EASValsSetEdgeFlag(true);
    EASValsSetLayer(OOFEG_RAW_GEOMETRY_LAYER);
    EASValsSetFillStyle(FILL_SOLID);
    p [ 0 ].x = ( FPNum ) this->giveNode(1)->giveCoordinate(1);
    p [ 0 ].y = ( FPNum ) this->giveNode(1)->giveCoordinate(2);
    p [ 0 ].z = ( FPNum ) this->giveNode(1)->giveCoordinate(3);
    p [ 1 ].x = ( FPNum ) this->giveNode(2)->giveCoordinate(1);
    p [ 1 ].y = ( FPNum ) this->giveNode(2)->giveCoordinate(2);
    p [ 1 ].z = ( FPNum ) this->giveNode(2)->giveCoordinate(3);
    p [ 2 ].x = ( FPNum ) this->giveNode(3)->giveCoordinate(1);
    p [ 2 ].y = ( FPNum ) this->giveNode(3)->giveCoordinate(2);
    p [ 2 ].z = ( FPNum ) this->giveNode(3)->giveCoordinate(3);
    p [ 3 ].x = ( FPNum ) this->giveNode(4)->giveCoordinate(1);
    p [ 3 ].y = ( FPNum ) this->giveNode(4)->giveCoordinate(2);
    p [ 3 ].z = ( FPNum ) this->giveNode(4)->giveCoordinate(3);

    go =  CreateTetra(p);
    EGWithMaskChangeAttributes(WIDTH_MASK | FILL_MASK | COLOR_MASK | EDGE_COLOR_MASK | EDGE_FLAG_MASK | LAYER_MASK, go);
    EGAttachObject(go, ( EObjectP ) this);
    EMAddGraphicsToModel(ESIModel(), go);
}

#endif
} // end namespace oofem
