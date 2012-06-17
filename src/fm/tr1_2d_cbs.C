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

#include "tr1_2d_cbs.h"
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
#include "contextioerr.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "conTable.h"
#endif

namespace oofem {
#define TRSUPG_ZERO_VOF 1.e-8


FEI2dTrLin TR1_2D_CBS :: interp(1,2);

TR1_2D_CBS :: TR1_2D_CBS(int n, Domain *aDomain) :
    CBSElement(n, aDomain)
    //<RESTRICTED_SECTION>
    , LEPlicElementInterface()
    //</RESTRICTED_SECTION>
{
    // Constructor.
    numberOfDofMans  = 3;
}

TR1_2D_CBS :: ~TR1_2D_CBS()
// Destructor
{ }



int
TR1_2D_CBS :: computeNumberOfDofs(EquationID ut)
{
    if ( ut == EID_MomentumBalance || ut == EID_AuxMomentumBalance ) {
        return 6;
    } else if ( ut == EID_ConservationEquation ) {
        return 3;
    } else {
        _error("computeNumberOfDofs: Unknown equation id encountered");
    }

    return 0;
}

void
TR1_2D_CBS :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    if ( ( ut == EID_MomentumBalance ) || ( ut == EID_AuxMomentumBalance ) ) {
        answer.setValues(2, V_u, V_v);
    } else if ( ut == EID_ConservationEquation ) {
        answer.setValues(1, P_f);
    } else {
        _error("giveDofManDofIDMask: Unknown equation id encountered");
    }
}

void
TR1_2D_CBS :: giveElementDofIDMask(EquationID ut, IntArray &answer) const
{
    this->giveDofManDofIDMask(1, ut, answer);
}


IRResultType
TR1_2D_CBS :: initializeFrom(InputRecord *ir)
{
    //<RESTRICTED_SECTION>
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                // Required by IR_GIVE_FIELD macro
    //</RESTRICTED_SECTION>

    this->CBSElement :: initializeFrom(ir);

    //<RESTRICTED_SECTION>
    this->vof = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, vof, IFT_TR12DCBS_pvof, "pvof");
    if ( vof > 0.0 ) {
        setPermanentVolumeFraction(vof);
        this->temp_vof = vof;
    } else {
        this->vof = 0.0;
        IR_GIVE_OPTIONAL_FIELD(ir, vof, IFT_TR12DCBS_vof, "vof");
        this->temp_vof = this->vof;
    }

    //</RESTRICTED_SECTION>

    this->computeGaussPoints();
    return IRRT_OK;
}

void
TR1_2D_CBS :: computeGaussPoints()
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
TR1_2D_CBS :: computeConsistentMassMtrx(FloatMatrix &answer, TimeStep *atTime)
{
    answer.resize(6, 6);
    answer.zero();
    //double rho = this->giveMaterial()->give('d');
    double rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), atTime);

    double ar6 = rho * area / 6.0;
    double ar12 = rho * area / 12.0;

    answer.at(1, 1) = answer.at(2, 2) = answer.at(3, 3) = ar6;
    answer.at(4, 4) = answer.at(5, 5) = answer.at(6, 6) = ar6;

    answer.at(1, 3) = answer.at(1, 5) = ar12;
    answer.at(3, 1) = answer.at(3, 5) = ar12;
    answer.at(5, 1) = answer.at(5, 3) = ar12;

    answer.at(2, 4) = answer.at(2, 6) = ar12;
    answer.at(4, 2) = answer.at(4, 6) = ar12;
    answer.at(6, 2) = answer.at(6, 4) = ar12;
}


void
TR1_2D_CBS :: computeDiagonalMassMtrx(FloatArray &answer, TimeStep *atTime)
{
    int i;
    answer.resize(6);

    //double rho = this->giveMaterial()->give('d');
    double rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), atTime);
    double mm = rho * this->area / 3.0;
    for ( i = 1; i <= 6; i++ ) {
        answer.at(i) = mm;
    }
}

void
TR1_2D_CBS :: computeConvectionTermsI(FloatArray &answer, TimeStep *stepN)
{
    // calculates advection component for (*) velocities

    double ar12, ar3, dudx, dudy, dvdx, dvdy;
    double usum, vsum;
    double adu11, adu21, adu31, adv11, adv21, adv31;
    double adu12, adu22, adu32, adv12, adv22, adv32;
    double dt = stepN->giveTimeIncrement();
    //double rho = this->giveMaterial()->give('d');
    double rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), stepN);
    int nLoads, i;
    bcGeomType ltype;
    Load *load;
    FloatArray gVector;

    FloatArray u;

    ar12 = area / 12.;
    ar3 = area / 3.0;

    this->computeVectorOf(EID_MomentumBalance, VM_Total, stepN->givePreviousStep(), u);

    dudx = b [ 0 ] * u.at(1) + b [ 1 ] * u.at(3) + b [ 2 ] * u.at(5);
    dudy = c [ 0 ] * u.at(1) + c [ 1 ] * u.at(3) + c [ 2 ] * u.at(5);
    dvdx = b [ 0 ] * u.at(2) + b [ 1 ] * u.at(4) + b [ 2 ] * u.at(6);
    dvdy = c [ 0 ] * u.at(2) + c [ 1 ] * u.at(4) + c [ 2 ] * u.at(6);

    usum = u.at(1) + u.at(3) + u.at(5);
    vsum = u.at(2) + u.at(4) + u.at(6);

    // compute Cu*U term

    adu11 = ar12 * ( dudx * ( usum + u.at(1) ) + dudy * ( vsum + u.at(2) ) );
    adu21 = ar12 * ( dudx * ( usum + u.at(3) ) + dudy * ( vsum + u.at(4) ) );
    adu31 = ar12 * ( dudx * ( usum + u.at(5) ) + dudy * ( vsum + u.at(6) ) );
    adv11 = ar12 * ( dvdx * ( usum + u.at(1) ) + dvdy * ( vsum + u.at(2) ) );
    adv21 = ar12 * ( dvdx * ( usum + u.at(3) ) + dvdy * ( vsum + u.at(4) ) );
    adv31 = ar12 * ( dvdx * ( usum + u.at(5) ) + dvdy * ( vsum + u.at(6) ) );


    // compute Ku*U term
    double uu = ( usum * usum + u.at(1) * u.at(1) + u.at(3) * u.at(3) + u.at(5) * u.at(5) );
    double uv = ( usum * vsum + u.at(1) * u.at(2) + u.at(3) * u.at(4) + u.at(5) * u.at(6) );
    double vv = ( vsum * vsum + u.at(2) * u.at(2) + u.at(4) * u.at(4) + u.at(6) * u.at(6) );

    adu12 = dt * 0.5 * ar12 * ( b [ 0 ] * dudx * uu + b [ 0 ] * dudy * uv + c [ 0 ] * dudx * uv + c [ 0 ] * dudy * vv );
    adu22 = dt * 0.5 * ar12 * ( b [ 1 ] * dudx * uu + b [ 1 ] * dudy * uv + c [ 1 ] * dudx * uv + c [ 1 ] * dudy * vv );
    adu32 = dt * 0.5 * ar12 * ( b [ 2 ] * dudx * uu + b [ 2 ] * dudy * uv + c [ 2 ] * dudx * uv + c [ 2 ] * dudy * vv );
    adv12 = dt * 0.5 * ar12 * ( b [ 0 ] * dvdx * uu + b [ 0 ] * dvdy * uv + c [ 0 ] * dvdx * uv + c [ 0 ] * dvdy * vv );
    adv22 = dt * 0.5 * ar12 * ( b [ 1 ] * dvdx * uu + b [ 1 ] * dvdy * uv + c [ 1 ] * dvdx * uv + c [ 1 ] * dvdy * vv );
    adv32 = dt * 0.5 * ar12 * ( b [ 2 ] * dvdx * uu + b [ 2 ] * dvdy * uv + c [ 2 ] * dvdx * uv + c [ 2 ] * dvdy * vv );

    answer.resize(6);

    answer.at(1) = -adu11 - adu12;
    answer.at(3) = -adu21 - adu22;
    answer.at(5) = -adu31 - adu32;
    answer.at(2) = -adv11 - adv12;
    answer.at(4) = -adv21 - adv22;
    answer.at(6) = -adv31 - adv32;


    // body load (gravity effects)
    nLoads    = this->giveBodyLoadArray()->giveSize();
    for ( i = 1; i <= nLoads; i++ ) {
        load  = domain->giveLoad( bodyLoadArray.at(i) );
        ltype = load->giveBCGeoType();
        if ( ( ltype == BodyLoadBGT ) && ( load->giveBCValType() == ForceLoadBVT ) ) {
            load->computeComponentArrayAt(gVector, stepN, VM_Total);
            if ( gVector.giveSize() ) {
                answer.at(1) -= dt * 0.5 * ar3 * ( b [ 0 ] * usum + c [ 0 ] * vsum ) * gVector.at(1);
                answer.at(2) -= dt * 0.5 * ar3 * ( b [ 0 ] * usum + c [ 0 ] * vsum ) * gVector.at(2);
                answer.at(3) -= dt * 0.5 * ar3 * ( b [ 1 ] * usum + c [ 1 ] * vsum ) * gVector.at(1);
                answer.at(4) -= dt * 0.5 * ar3 * ( b [ 1 ] * usum + c [ 1 ] * vsum ) * gVector.at(2);
                answer.at(5) -= dt * 0.5 * ar3 * ( b [ 2 ] * usum + c [ 2 ] * vsum ) * gVector.at(1);
                answer.at(6) -= dt * 0.5 * ar3 * ( b [ 2 ] * usum + c [ 2 ] * vsum ) * gVector.at(2);
            }
        }
    }


    answer.times(rho);
}

void
TR1_2D_CBS :: computeDiffusionTermsI(FloatArray &answer, TimeStep *tStep)
{
    FloatArray stress;
    int i, j, nLoads;
    double Re = domain->giveEngngModel()->giveUnknownComponent(ReynoldsNumber, VM_Unknown, tStep, domain, NULL);
    double coeff, rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    bcGeomType ltype;
    Load *load;
    FloatArray gVector;
    double ar3 = area / 3.0;

    stress = ( ( FluidDynamicMaterialStatus * ) this->giveMaterial()->giveStatus(gp) )->giveDeviatoricStressVector();
    stress.times(1. / Re);

    // \int dNu/dxj \Tau_ij
    answer.resize(6);
    for ( i = 0; i < 3; i++ ) {
        answer.at( ( i ) * 2 + 1 ) = -area * ( stress.at(1) * b [ i ] + stress.at(3) * c [ i ] );
        answer.at( ( i + 1 ) * 2 ) = -area * ( stress.at(3) * b [ i ] + stress.at(2) * c [ i ] );
    }

    // add boundary termms
    coeff = ar3 * rho;
    // body load (gravity effects)
    nLoads    = this->giveBodyLoadArray()->giveSize();
    for ( i = 1; i <= nLoads; i++ ) {
        load  = domain->giveLoad( bodyLoadArray.at(i) );
        ltype = load->giveBCGeoType();
        if ( ( ltype == BodyLoadBGT ) && ( load->giveBCValType() == ForceLoadBVT ) ) {
            load->computeComponentArrayAt(gVector, tStep, VM_Total);
            if ( gVector.giveSize() ) {
                answer.at(1) += coeff * gVector.at(1);
                answer.at(2) += coeff * gVector.at(2);
                answer.at(3) += coeff * gVector.at(1);
                answer.at(4) += coeff * gVector.at(2);
                answer.at(5) += coeff * gVector.at(1);
                answer.at(6) += coeff * gVector.at(2);
            }
        }
    }

    // terms now computed even for boundary nodes that have prescribed velocities! (but they will not be localized)
    // loop over sides
    int n1, n2, code;
    double tx, ty, l, nx, ny;
    nLoads = boundarySides.giveSize();
    for ( j = 1; j <= nLoads; j++ ) {
        code = boundaryCodes.at(j);
        i = boundarySides.at(j);
        if ( ( code & FMElement_PrescribedTractionBC ) ) {
            //printf ("TR1_2D_CBS :: computeDiffusionTermsI tractions detected\n");

            //_error ("computeDiffusionTermsI: traction bc not supported");
            FloatArray t, coords(1);
            int nLoads, n, id;
            BoundaryLoad *load;
            // integrate tractions
            n1 = i;
            n2 = ( n1 == 3 ? 1 : n1 + 1 );

            tx = giveNode(n2)->giveCoordinate(1) - giveNode(n1)->giveCoordinate(1);
            ty = giveNode(n2)->giveCoordinate(2) - giveNode(n1)->giveCoordinate(2);
            l = sqrt(tx * tx + ty * ty);
            nx = ty / l;
            ny = -tx / l;

            // loop over boundary load array
            nLoads    = this->giveBoundaryLoadArray()->giveSize() / 2;
            for ( i = 1; i <= nLoads; i++ ) {
                n     = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
                id    = boundaryLoadArray.at(i * 2);
                load  = dynamic_cast< BoundaryLoad * >( domain->giveLoad(n) );
                if ( load ) {
                    load->computeValueAt(t, tStep, coords, VM_Total);

                    //printf ("TR1_2D_CBS :: computeDiffusionTermsI traction (%e,%e) detected\n", t.at(1), t.at(2));

                    answer.at( ( n1 - 1 ) * 2 + 1 ) += 0.5 * l * ( t.at(1) * nx );
                    answer.at( ( n1 ) * 2 )    += 0.5 * l * ( t.at(2) * ny );

                    answer.at( ( n2 - 1 ) * 2 + 1 ) += 0.5 * l * ( t.at(1) * nx );
                    answer.at( ( n2 ) * 2 )    += 0.5 * l * ( t.at(2) * ny );
                }
            }
        } else if ( !( ( code & FMElement_PrescribedUnBC ) && ( code & FMElement_PrescribedUsBC ) ) ) {
            /*
             * n1 = i;
             * n2 = (n1==3?n2=1:n2=n1+1);
             *
             * //if (giveNode(n1)->isBoundary() && giveNode(n2)->isBoundary()) {
             *
             * tx = giveNode(n2)->giveCoordinate(1)-giveNode(n1)->giveCoordinate(1);
             * ty = giveNode(n2)->giveCoordinate(2)-giveNode(n1)->giveCoordinate(2);
             * l = sqrt(tx*tx+ty*ty);
             * nx = ty/l; ny = -tx/l;
             *
             * // normal displacement precribed
             * answer.at((n1-1)*2+1)+=l*0.5*(stress.at(1)*nx + stress.at(3)*ny);
             * answer.at((n1)*2)    +=l*0.5*(stress.at(3)*nx + stress.at(2)*ny);
             *
             * answer.at((n2-1)*2+1)+=l*0.5*(stress.at(1)*nx + stress.at(3)*ny);
             * answer.at((n2)*2)    +=l*0.5*(stress.at(3)*nx + stress.at(2)*ny);
             * //}
             */
        }
    }
}

void
TR1_2D_CBS :: computeDensityRhsVelocityTerms(FloatArray &answer, TimeStep *tStep)
{
    // computes velocity terms on RHS for density equation

    int i, j, uadr, vadr;
    double velu = 0.0, velv = 0.0; // dudx=0.0, dvdy=0.0;
    double theta1 = domain->giveEngngModel()->giveUnknownComponent(Theta_1, VM_Unknown, tStep, domain, NULL);
    //double rho = this->giveMaterial()->give('d');
    double rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), tStep);
    FloatArray u(6), ustar(6);

    answer.resize(3);

    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep->givePreviousStep(), u);
    this->computeVectorOf(EID_AuxMomentumBalance, VM_Incremental, tStep, ustar);
    ustar.times(theta1);
    u.add(ustar);

    for ( i = 0; i < 3; i++ ) {
        uadr = i * 2 + 1;
        vadr = ( i + 1 ) * 2;
        velu += u.at(uadr);
        velv += u.at(vadr);
    }

    for ( i = 0; i < 3; i++ ) {
        answer.at(i + 1) = area * ( ( b [ i ] * velu + c [ i ] * velv ) ) / 3.0;
    }

    /* account for normal prescribed velocity on boundary*/
    /* on the rest of the boundary the tractions or pressure is prescribed -> to be implemented later */
    int n1, n2;
    double tx, ty, l, nx, ny, un1, un2;


    // loop over sides
    int code;
    for ( j = 1; j <= boundarySides.giveSize(); j++ ) {
        code = boundaryCodes.at(j);
        if ( ( code & FMElement_PrescribedPressureBC ) ) {
            continue;
        } else if ( ( code & FMElement_PrescribedUnBC ) ) {
            this->computeVectorOfPrescribed(EID_MomentumBalance, VM_Total, tStep, u);

            i = boundarySides.at(j);
            n1 = i;
            n2 = ( n1 == 3 ? 1 : n1 + 1 );

            //if (giveNode(n1)->isBoundary() && giveNode(n2)->isBoundary()) {

            tx = giveNode(n2)->giveCoordinate(1) - giveNode(n1)->giveCoordinate(1);
            ty = giveNode(n2)->giveCoordinate(2) - giveNode(n1)->giveCoordinate(2);
            l = sqrt(tx * tx + ty * ty);
            nx = ty / l;
            ny = -tx / l;
            un1 = nx * u.at( ( n1 - 1 ) * 2 + 1 ) + ny *u.at(n1 * 2);
            un2 = nx * u.at( ( n2 - 1 ) * 2 + 1 ) + ny *u.at(n2 * 2);
            //if ((un1 != 0.) && (un2 != 0.)) {
            // normal displacement precribed
            answer.at(n1) -= ( un1 * l / 3. + un2 * l / 6. );
            answer.at(n2) -= ( un2 * l / 3. + un1 * l / 6. );
            //}
            //}
        } else if ( ( code & FMElement_PrescribedTractionBC ) ) {
            continue;
        }
    }

    answer.times(rho);
}


void
TR1_2D_CBS :: computePrescribedTractionPressure(FloatArray &answer, TimeStep *tStep)
{
    /*
     * this method computes the prescribed pressure due to applied traction
     * p = tau(i,j)*n(i)*n(j) - traction(i)*n(i)
     * this pressure is enforced as dirichlet bc in density/pressure equation
     */
    FloatArray stress;
    int i, j;
    double Re = domain->giveEngngModel()->giveUnknownComponent(ReynoldsNumber, VM_Unknown, tStep, domain, NULL);
    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    stress = ( ( FluidDynamicMaterialStatus * ) this->giveMaterial()->giveStatus(gp) )->giveDeviatoricStressVector();
    stress.times(1. / Re);

    answer.resize(3);
    answer.zero();

    // loop over sides
    int n1, n2, code, sid;
    double tx, ty, l, nx, ny, pcoeff;
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


            //nodecounter.at(n1)++;
            //nodecounter.at(n2)++;
            pcoeff = stress.at(1) * nx * nx + stress.at(2) * ny * ny + 2.0 * stress.at(3) * nx * ny;
            answer.at(n1) += pcoeff;
            answer.at(n2) += pcoeff;

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
                    load->computeValueAt(t, tStep, coords, VM_Total);

                    answer.at(n1) -= t.at(1) * nx + t.at(2) * ny;
                    answer.at(n2) -= t.at(1) * nx + t.at(2) * ny;
                }
            }
        }
    }

    //for (i=1; i<=3; i++) {
    //  if (nodecounter.at(i)) answer.at(i) /= nodecounter.at(i);
    //}
}


void
TR1_2D_CBS :: computeNumberOfNodalPrescribedTractionPressureContributions(FloatArray &answer, TimeStep *tStep)
{
    /*
     * this method computes the prescribed pressure due to applied traction
     * p = tau(i,j)*n(i)*n(j) - traction(i)*n(i)
     * this pressure is enforced as dirichlet bc in density/pressure equation
     */
    int i, j;

    answer.resize(3);
    answer.zero();

    // loop over sides
    int n1, n2, code;
    IntArray nodecounter(3);
    for ( j = 1; j <= boundarySides.giveSize(); j++ ) {
        code = boundaryCodes.at(j);
        i = boundarySides.at(j);
        if ( ( code & FMElement_PrescribedTractionBC ) ) {
            n1 = i;
            n2 = ( n1 == 3 ? 1 : n1 + 1 );
            answer.at(n1)++;
            answer.at(n2)++;
        }
    }
}




void
TR1_2D_CBS :: computeDensityRhsPressureTerms(FloatArray &answer, TimeStep *tStep)
{
    // computes pressure terms on RHS for density equation
    FloatArray p(3);
    int i;
    double theta1 = domain->giveEngngModel()->giveUnknownComponent(Theta_1, VM_Unknown, tStep, domain, NULL);

    this->computeVectorOf(EID_ConservationEquation, VM_Total, tStep->givePreviousStep(), p);
    answer.resize(3);

    double dpdx = 0.0, dpdy = 0.0;
    for ( i = 0; i < 3; i++ ) {
        dpdx += b [ i ] * p.at(i + 1);
        dpdy += c [ i ] * p.at(i + 1);
    }

    for ( i = 0; i < 3; i++ ) {
        answer.at(i + 1) = -theta1 *tStep->giveTimeIncrement() * area * ( b [ i ] * dpdx + c [ i ] * dpdy );
    }
}


void
TR1_2D_CBS :: computePressureLhs(FloatMatrix &answer, TimeStep *tStep)
{
    // calculates the pressure LHS

    answer.resize(3, 3);

    answer.at(1, 1) = area * ( b [ 0 ] * b [ 0 ] + c [ 0 ] * c [ 0 ] );
    answer.at(2, 2) = area * ( b [ 1 ] * b [ 1 ] + c [ 1 ] * c [ 1 ] );
    answer.at(3, 3) = area * ( b [ 2 ] * b [ 2 ] + c [ 2 ] * c [ 2 ] );

    answer.at(1, 2) = answer.at(2, 1) = area * ( b [ 0 ] * b [ 1 ] + c [ 0 ] * c [ 1 ] );
    answer.at(1, 3) = answer.at(3, 1) = area * ( b [ 0 ] * b [ 2 ] + c [ 0 ] * c [ 2 ] );
    answer.at(2, 3) = answer.at(3, 2) = area * ( b [ 1 ] * b [ 2 ] + c [ 1 ] * c [ 2 ] );
}


void
TR1_2D_CBS :: computeCorrectionRhs(FloatArray &answer, TimeStep *tStep)
{
    //Evaluates the RHS of velocity correction step
    FloatArray p(3), u(6);
    int i;
    double pn1, ar3;
    double usum, vsum, coeff;


    this->computeVectorOf(EID_ConservationEquation, VM_Total, tStep, p);

    double dpdx = 0.0, dpdy = 0.0;
    for ( i = 0; i < 3; i++ ) {
        pn1 = p.at(i + 1);
        dpdx += b [ i ] * pn1;
        dpdy += c [ i ] * pn1;
    }

    answer.resize(6);

    ar3 = area / 3.0;
    answer.at(1) = answer.at(3) = answer.at(5) = -ar3 * dpdx;
    answer.at(2) = answer.at(4) = answer.at(6) = -ar3 * dpdy;


    this->computeVectorOf(EID_ConservationEquation, VM_Total, tStep->givePreviousStep(), p);
    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep->givePreviousStep(), u);
    dpdx = 0.0, dpdy = 0.0;
    for ( i = 0; i < 3; i++ ) {
        pn1 = p.at(i + 1);
        dpdx += b [ i ] * pn1;
        dpdy += c [ i ] * pn1;
    }

    usum = u.at(1) + u.at(3) + u.at(5);
    vsum = u.at(2) + u.at(4) + u.at(6);
    coeff = ar3 * tStep->giveTimeIncrement() / 2.0;
    answer.at(1) -= coeff * dpdx * ( b [ 0 ] * usum + c [ 0 ] * vsum );
    answer.at(3) -= coeff * dpdx * ( b [ 1 ] * usum + c [ 1 ] * vsum );
    answer.at(5) -= coeff * dpdx * ( b [ 2 ] * usum + c [ 2 ] * vsum );

    answer.at(2) -= coeff * dpdy * ( b [ 0 ] * usum + c [ 0 ] * vsum );
    answer.at(4) -= coeff * dpdy * ( b [ 1 ] * usum + c [ 1 ] * vsum );
    answer.at(6) -= coeff * dpdy * ( b [ 2 ] * usum + c [ 2 ] * vsum );
}


Interface *
TR1_2D_CBS :: giveInterface(InterfaceType interface)
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
    }
    //<RESTRICTED_SECTION>
    else if ( interface == LEPlicElementInterfaceType ) {
        return ( LEPlicElementInterface * ) this;
    }

    //</RESTRICTED_SECTION>
    return NULL;
}


int
TR1_2D_CBS :: SpatialLocalizerI_containsPoint(const FloatArray &coords) {
    FloatArray lcoords;
    return this->computeLocalCoordinates(lcoords, coords);
}

double
TR1_2D_CBS :: SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords)
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
TR1_2D_CBS :: computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    /* one should call material driver instead */
    FloatArray u(6), eps(3);
    answer.resize(3);


    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);

    eps.at(1) = ( b [ 0 ] * u.at(1) + b [ 1 ] * u.at(3) + b [ 2 ] * u.at(5) );
    eps.at(2) = ( c [ 0 ] * u.at(2) + c [ 1 ] * u.at(4) + c [ 2 ] * u.at(6) );
    eps.at(3) = ( b [ 0 ] * u.at(2) + b [ 1 ] * u.at(4) + b [ 2 ] * u.at(6) + c [ 0 ] * u.at(1) + c [ 1 ] * u.at(3) + c [ 2 ] * u.at(5) );
    ( ( FluidDynamicMaterial * ) this->giveMaterial() )->computeDeviatoricStressVector(answer, gp, eps, tStep);
}

int
TR1_2D_CBS :: checkConsistency()
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

    b [ 0 ] = ( y2 - y3 ) / ( 2. * area );
    c [ 0 ] = ( x3 - x2 ) / ( 2. * area );
    b [ 1 ] = ( y3 - y1 ) / ( 2. * area );
    c [ 1 ] = ( x1 - x3 ) / ( 2. * area );
    b [ 2 ] = ( y1 - y2 ) / ( 2. * area );
    c [ 2 ] = ( x2 - x1 ) / ( 2. * area );

    return CBSElement :: checkConsistency();
}

double
TR1_2D_CBS :: computeCriticalTimeStep(TimeStep *tStep)
{
    FloatArray u;
    double dt1, dt2, dt;
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

    // viscous limit
    dt2 = 0.5 * ln * ln * Re;
    if ( veln != 0.0 ) {
        dt1 = ln / veln;
        dt = dt1 * dt2 / ( dt1 + dt2 );
    } else {
        dt = dt2;
    }

    return dt;
}

//<RESTRICTED_SECTION>
double
TR1_2D_CBS :: computeLEPLICVolumeFraction(const FloatArray &n, const double p, LEPlic *matInterface, bool updFlag)
{
    Polygon pg;
    double answer, volume = computeMyVolume(matInterface, updFlag);
    this->formVolumeInterfacePoly(pg, matInterface, n, p, updFlag);
    answer = fabs(pg.computeVolume() / volume);
    if ( answer > 1.000000001 ) {
        return 1.0;
    } else {
        return answer;
    }
}

void
TR1_2D_CBS :: formMaterialVolumePoly(Polygon &matvolpoly, LEPlic *matInterface,
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
TR1_2D_CBS :: formVolumeInterfacePoly(Polygon &matvolpoly, LEPlic *matInterface,
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


double
TR1_2D_CBS :: truncateMatVolume(const Polygon &matvolpoly, double &volume)
{
    Polygon me, clip;
    Graph g;

    this->formMyVolumePoly(me, NULL, false);
    g.clip(clip, me, matvolpoly);
#ifdef __OOFEG
    EASValsSetColor( gc [ 0 ].getActiveCrackColor() );
    clip.draw(gc [ OOFEG_DEBUG_LAYER ], true);
    //EVFastRedraw(myview);
#endif
    volume = clip.computeVolume();
    return volume / area;
}

void
TR1_2D_CBS :: formMyVolumePoly(Polygon &me, LEPlic *matInterface, bool updFlag)
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
TR1_2D_CBS :: computeMyVolume(LEPlic *matInterface, bool updFlag)
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
TR1_2D_CBS :: giveElementCenter(LEPlic *mat_interface, FloatArray &center, bool upd)
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

//</RESTRICTED_SECTION>


int
TR1_2D_CBS :: EIPrimaryFieldI_evaluateFieldVectorAt(FloatArray &answer, PrimaryField &pf,
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
TR1_2D_CBS :: updateYourself(TimeStep *tStep)
{
    CBSElement :: updateYourself(tStep);
    //<RESTRICTED_SECTION>
    LEPlicElementInterface :: updateYourself(tStep);
    //</RESTRICTED_SECTION>
}

int
TR1_2D_CBS :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    //<RESTRICTED_SECTION>
    if ( type == IST_VOFFraction ) {
        answer.resize(1);
        answer.at(1) = this->giveTempVolumeFraction();
        return 1;
    } else
    //</RESTRICTED_SECTION>
    if ( type == IST_Density ) {
        answer.resize(1);
        answer.at(1) = this->giveMaterial()->giveCharacteristicValue(MRM_Density, aGaussPoint, atTime);
        return 1;
    } else {
        return CBSElement :: giveIPValue(answer, aGaussPoint, type, atTime);
    }
}

int
TR1_2D_CBS :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type)
{
    if ( ( type == IST_VOFFraction ) || ( type == IST_Density ) ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else {
        return CBSElement :: giveIntVarCompFullIndx(answer, type);
    }
}


InternalStateValueType
TR1_2D_CBS :: giveIPValueType(InternalStateType type)
{
    if ( ( type == IST_VOFFraction ) || ( type == IST_Density ) ) {
        return ISVT_SCALAR;
    } else {
        return CBSElement :: giveIPValueType(type);
    }
}


int
TR1_2D_CBS :: giveIPValueSize(InternalStateType type, GaussPoint *gp)
{
    return CBSElement::giveIPValueSize(type, gp);
}


int
TR1_2D_CBS :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( ( type == IST_StressTensor ) || ( type == IST_StrainTensor ) ) {
        return 4;
    }

    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    return this->giveIPValueSize(type, gp);
}


void
TR1_2D_CBS :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                         InternalStateType type, TimeStep *tStep)
{
    GaussPoint *gp;
    gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    this->giveIPValue(answer, gp, type, tStep);
}

void
TR1_2D_CBS :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                        InternalStateType type, TimeStep *tStep)
{
    answer.resize(0);
}

void
TR1_2D_CBS :: SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap)
{
    pap.resize(3);
    pap.at(1) = this->giveNode(1)->giveNumber();
    pap.at(2) = this->giveNode(2)->giveNumber();
    pap.at(3) = this->giveNode(3)->giveNumber();
}

void
TR1_2D_CBS :: SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap)
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
TR1_2D_CBS :: SPRNodalRecoveryMI_giveNumberOfIP()
{ return 1; }


void
TR1_2D_CBS :: SPRNodalRecoveryMI_computeIPGlobalCoordinates(FloatArray &coords, GaussPoint *gp)
{
    if ( gp == integrationRulesArray [ 0 ]->getIntegrationPoint(0) ) {
        this->computeGlobalCoordinates( coords, * gp->giveCoordinates() );
    } else {
        _error("SPRNodalRecoveryMI_computeIPGlobalCoordinates: unsupported ip num");
    }
}

SPRPatchType
TR1_2D_CBS :: SPRNodalRecoveryMI_givePatchType()
{
    return SPRPatchType_2dxy;
}



void
TR1_2D_CBS :: printOutputAt(FILE *file, TimeStep *stepN)
// Performs end-of-step operations.
{
    CBSElement :: printOutputAt(file, stepN);
    //<RESTRICTED_SECTION>
    double rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), stepN);
    fprintf(file, "VOF %e, density %e\n\n", this->giveVolumeFraction(), rho);
    //</RESTRICTED_SECTION>
}



contextIOResultType TR1_2D_CBS :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full element context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;

    if ( ( iores = CBSElement :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    //<RESTRICTED_SECTION>
    if ( ( iores = LEPlicElementInterface :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    //</RESTRICTED_SECTION>

    return CIO_OK;
}



contextIOResultType TR1_2D_CBS :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full element context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;

    if ( ( iores = CBSElement :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    //<RESTRICTED_SECTION>
    if ( ( iores = LEPlicElementInterface :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    //</RESTRICTED_SECTION>


    return CIO_OK;
}




#ifdef __OOFEG
int
TR1_2D_CBS :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                      int node, TimeStep *atTime)
{
    //<RESTRICTED_SECTION>
    if ( type == IST_VOFFraction ) {
        answer.resize(1);
        answer.at(1) = this->giveTempVolumeFraction();
        return 1;
    } else
    //</RESTRICTED_SECTION>
    if ( type == IST_Density ) {
        answer.resize(1);
        answer.at(1) = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), atTime);
        return 1;
    } else {
        return CBSElement :: giveInternalStateAtNode(answer, type, mode, node, atTime);
    }
}

void
TR1_2D_CBS :: drawRawGeometry(oofegGraphicContext &gc)
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

void TR1_2D_CBS :: drawScalar(oofegGraphicContext &context)
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
