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

#include "quad10_2d_supg.h"
#include "fei2dquadlin.h"
#include "fei2dquadconst.h"
#include "dofmanager.h"
#include "material.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "flotmtrx.h"
#include "flotarry.h"
#include "intarray.h"
#include "domain.h"
#include "mathfem.h"
#include "engngm.h"
#include "timestep.h"
#include "materialinterface.h"
#include "contextioerr.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

namespace oofem {
FEI2dQuadLin Quad10_2D_SUPG :: velocityInterpolation(1, 2);
FEI2dQuadConst Quad10_2D_SUPG :: pressureInterpolation(1, 2);


Quad10_2D_SUPG :: Quad10_2D_SUPG(int n, Domain *aDomain) :
    SUPGElement2(n, aDomain), pressureNode(1, aDomain, this)
{
    numberOfDofMans = 4;
}

Quad10_2D_SUPG :: ~Quad10_2D_SUPG()
{ }

FEInterpolation *
Quad10_2D_SUPG :: giveInterpolation()
{
    return & this->velocityInterpolation;
}

FEInterpolation *
Quad10_2D_SUPG :: giveInterpolation(DofIDItem id)
{
    if (id == P_f) {
        return & this->pressureInterpolation;
    } else {
        return & this->velocityInterpolation;
    }
}

DofManager *
Quad10_2D_SUPG :: giveInternalDofManager(int i) const
{
    //_error2("No such DOF available on Element %d", number);
    return ( DofManager * ) & pressureNode;
}


void Quad10_2D_SUPG :: giveLocationArray(IntArray &locationArray, EquationID ut, const UnknownNumberingScheme &s) const
// Returns the location array of the receiver. This array is obtained by
// simply appending the location array of every node of the receiver.
{
    IntArray nodeDofIDMask;
    IntArray nodalArray;
    int i;


    locationArray.resize(0);
    for ( i = 1; i <= numberOfDofMans; i++ ) {
        this->giveDofManDofIDMask(i, ut, nodeDofIDMask);
        this->giveDofManager(i)->giveLocationArray(nodeDofIDMask, nodalArray, s);
        locationArray.followedBy(nodalArray);
    }

    for ( i = 1; i <= 1; i++ ) {
        this->giveInternalDofManDofIDMask(i, ut, nodeDofIDMask);
        this->giveInternalDofManager(i)->giveLocationArray(nodeDofIDMask, nodalArray, s);
        locationArray.followedBy(nodalArray);
    }
}


int
Quad10_2D_SUPG :: giveTermIntergationRuleIndex(CharType termType)
{
    if ( ( termType == AccelerationTerm_MB ) || ( termType == AdvectionTerm_MB ) ||
        ( termType == AdvectionDerivativeTerm_MB ) ) {
        return 0;
    } else if ( ( termType == DiffusionTerm_MB ) || ( termType == DiffusionDerivativeTerm_MB ) ||
               ( termType == PressureTerm_MB ) || ( termType == AdvectionTerm_MC ) || ( termType == AdvectionDerivativeTerm_MC ) ||
               ( termType == DiffusionDerivativeTerm_MC ) || ( termType == BCRhsTerm_MC ) ) {
        return 0;
    } else if ( ( termType == LSICStabilizationTerm_MB ) || ( termType == LinearAdvectionTerm_MC ) ||
               ( termType == DiffusionTerm_MC ) || ( termType == AccelerationTerm_MC ) ||
               ( termType == PressureTerm_MC ) || ( termType ==  BCRhsTerm_MB ) ) {
        return 0;
    } else {
        _error2( "giveTermIntergationRuleIndex: Unknown CharType encountered [%s]", __CharTypeToString(termType) );
    }

    return 0;
}


int
Quad10_2D_SUPG :: computeNumberOfDofs(EquationID ut)
{
    if ( ut == EID_MomentumBalance ) {
        return 8;
    } else if ( ut == EID_ConservationEquation ) {
        return 1;
    } else if ( ut == EID_MomentumBalance_ConservationEquation ) {
        return 9;
    } else {
        _error("computeNumberOfDofs: Unknown equation id encountered");
    }

    return 0;
}


void
Quad10_2D_SUPG :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    if ( ( ut == EID_MomentumBalance ) || ( ut == EID_AuxMomentumBalance ) || ( ut == EID_MomentumBalance_ConservationEquation ) ) {
        answer.setValues(2, V_u, V_v);
    } else if ( ut == EID_ConservationEquation ) {
        answer.resize(0);
    } else  {
        _error("giveDofManDofIDMask: Unknown equation id encountered");
    }
}


void
Quad10_2D_SUPG ::   giveInternalDofManDofIDMask(int i, EquationID ut, IntArray &answer) const
{
    if ( ( ut == EID_MomentumBalance ) || ( ut == EID_AuxMomentumBalance ) ) {
        answer.resize(0);
    } else if ( ( ut == EID_ConservationEquation ) || ( ut == EID_MomentumBalance_ConservationEquation ) ) {
        answer.setValues(1, P_f);
    } else  {
        _error("giveDofManDofIDMask: Unknown equation id encountered");
    }
}


IRResultType
Quad10_2D_SUPG :: initializeFrom(InputRecord *ir)
{
    SUPGElement2 :: initializeFrom(ir);
    this->pressureNode.initializeFrom(ir);

    return IRRT_OK;
}


void
Quad10_2D_SUPG :: computeGaussPoints()
// Sets up the array containing the four Gauss points of the receiver.
{
    if ( !integrationRulesArray ) {
        numberOfIntegrationRules = 2;
        integrationRulesArray = new IntegrationRule * [ 2 ];


        integrationRulesArray [ 0 ] = new GaussIntegrationRule(1, this, 1, 3);
        integrationRulesArray [ 0 ]->setUpIntegrationPoints(_Square, 4, _2dFlow);

        //seven point Gauss integration
        integrationRulesArray [ 1 ] = new GaussIntegrationRule(2, this, 1, 3);
        integrationRulesArray [ 1 ]->setUpIntegrationPoints(_Square, 4, _2dFlow);
    }
}


void
Quad10_2D_SUPG :: computeNuMatrix(FloatMatrix &answer, GaussPoint *gp)
{
    int i;
    FloatArray n;

    this->velocityInterpolation.evalN(n, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));
    answer.resize(2, 8);
    answer.zero();

    for ( i = 1; i <= 4; i++ ) {
        answer.at(1, 2 * i - 1)  = n.at(i);
        answer.at(2, 2 * i - 0)  = n.at(i);
    }
}


void
Quad10_2D_SUPG :: computeUDotGradUMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *atTime)
{
    int i;
    FloatMatrix n, dn;
    FloatArray u, un;
    this->velocityInterpolation.evaldNdx( dn, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));
    this->computeNuMatrix(n, gp);
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, un);

    u.beProductOf(n, un);

    answer.resize(2, 8);
    answer.zero();
    for ( i = 1; i <= 4; i++ ) {
        answer.at(1, 2 * i - 1) = dn.at(i, 1) * u.at(1) + dn.at(i, 2) * u.at(2);
        answer.at(2, 2 * i)   = dn.at(i, 1) * u.at(1) + dn.at(i, 2) * u.at(2);
    }
}

void
Quad10_2D_SUPG :: computeBMatrix(FloatMatrix &answer, GaussPoint *gp)
{
    int i;
    FloatMatrix dn(4, 2);
    this->velocityInterpolation.evaldNdx(dn, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));

    answer.resize(3, 8);
    answer.zero();

    for ( i = 1; i <= 4; i++ ) {
        answer.at(1, 2 * i - 1) = dn.at(i, 1);
        answer.at(2, 2 * i)   = dn.at(i, 2);
        answer.at(3, 2 * i - 1) = dn.at(i, 2);
        answer.at(3, 2 * i)   = dn.at(i, 1);
    }
}

void
Quad10_2D_SUPG :: computeDivUMatrix(FloatMatrix &answer, GaussPoint *gp)
{
    int i;
    FloatMatrix dn(4, 2);
    velocityInterpolation.evaldNdx(dn, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));

    answer.resize(1, 8);
    answer.zero();

    for ( i = 1; i <= 4; i++ ) {
        answer.at(1, 2 * i - 1) = dn.at(i, 1);
        answer.at(1, 2 * i) = dn.at(i, 2);
    }
}

void
Quad10_2D_SUPG :: computeNpMatrix(FloatMatrix &answer, GaussPoint *gp)
{
    FloatArray n(1);
    pressureInterpolation.evalN(n, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));

    answer.resize(1, 1);
    answer.zero();

    answer.at(1, 1)  = n.at(1);
}


void
Quad10_2D_SUPG :: computeGradUMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *atTime)
{
    int i;
    FloatArray dnx(4), dny(4), u, u1(4), u2(4);
    FloatMatrix dn;

    answer.resize(2, 2);
    answer.zero();

    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);

    velocityInterpolation.evaldNdx(dn, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));
    for ( i = 1; i <= 4; i++ ) {
        dnx.at(i) = dn.at(i, 1);
        dny.at(i) = dn.at(i, 2);

        u1.at(i) = u.at(2 * i - 1);
        u2.at(i) = u.at(2 * i);
    }


    answer.at(1, 1) =  u1.dotProduct(dnx);
    answer.at(1, 2) =  u1.dotProduct(dny);
    answer.at(2, 1) =  u2.dotProduct(dnx);
    answer.at(2, 2) =  u2.dotProduct(dny);
}

void
Quad10_2D_SUPG :: computeGradPMatrix(FloatMatrix &answer, GaussPoint *gp)
{
    FloatMatrix dn(1, 2);
    pressureInterpolation.evaldNdx(dn, * gp->giveCoordinates(), FEIElementGeometryWrapper(this));

    answer.beTranspositionOf(dn);
}



void
Quad10_2D_SUPG :: computeDivTauMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *atTime)
{
    answer.resize(2, 8);
    answer.zero();
}


void
Quad10_2D_SUPG :: updateStabilizationCoeffs(TimeStep *atTime)
{
    double Re, norm_un, mu, mu_min, nu, norm_N, norm_N_d, norm_M_d, norm_LSIC, norm_G_c, norm_M_c, norm_N_c, t_p1, t_p2, t_p3, t_s1, t_s2, t_s3, rho;
    FloatMatrix dn, N, N_d, M_d, LSIC, G_c, M_c, N_c;
    FloatArray dN, s, lcoords_nodes, u, lcn, dn_a(2), n, u1(8), u2(8);
    GaussPoint *gp;
    int j;
    IntegrationRule *iRule;


    iRule = integrationRulesArray [ 1 ];
    mu_min = 1;
    rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, integrationRulesArray [ 0 ]->getIntegrationPoint(0), atTime);
    for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
        gp = iRule->getIntegrationPoint(j);
        mu = this->giveMaterial()->giveCharacteristicValue(MRM_Viscosity, gp, atTime);
        if ( mu_min > mu ) {
            mu_min = mu;
        }

        nu = mu_min / rho;
    }

    nu = mu_min / rho;

    //this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime->givePreviousStep(), un);
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);

    norm_un = u.computeNorm();

    this->computeAdvectionTerm(N, atTime);
    this->computeAdvectionDeltaTerm(N_d, atTime);
    this->computeMassDeltaTerm(M_d, atTime);
    this->computeLSICTerm(LSIC, atTime);
    this->computeLinearAdvectionTerm_MC(G_c, atTime);
    this->computeMassEpsilonTerm(M_c, atTime);
    this->computeAdvectionEpsilonTerm(N_c, atTime);


    norm_N = N.computeFrobeniusNorm();
    norm_N_d = N_d.computeFrobeniusNorm();
    norm_M_d = M_d.computeFrobeniusNorm();
    norm_LSIC = LSIC.computeFrobeniusNorm();
    norm_G_c = G_c.computeFrobeniusNorm();
    norm_M_c = M_c.computeFrobeniusNorm();
    norm_N_c = N_c.computeFrobeniusNorm();

    if ( ( norm_N == 0 ) || ( norm_N_d == 0 ) || ( norm_M_d == 0 ) ) {
        t_supg = 0;
    } else   {
        Re = ( norm_un / nu ) * ( norm_N / norm_N_d );

        t_s1 = norm_N / norm_N_d;

        t_s2 = atTime->giveTimeIncrement() * ( norm_N / norm_M_d ) * 0.5;

        t_s3 = t_s1 * Re;

        t_supg =  1. / sqrt( 1. / ( t_s1 * t_s1 ) + 1. / ( t_s2 * t_s2 ) + 1. / ( t_s3 * t_s3 ) );
        // t_supg = 0;
    }

    if ( norm_LSIC == 0 ) {
        t_lsic = 0;
    } else   {
        t_lsic = norm_N / norm_LSIC;

        // t_lsic = 0;
    }

    if ( ( norm_G_c == 0 ) || ( norm_N_c == 0 ) || ( norm_M_c == 0 ) ) {
        t_pspg = 0;
    } else  {
        Re = ( norm_un / nu ) * ( norm_N / norm_N_d );

        t_p1 = norm_G_c / norm_N_c;

        t_p2 = atTime->giveTimeIncrement() * ( norm_G_c / norm_M_c ) * 0.5;

        t_p3 = t_p1 * Re;

        t_pspg =  1. / sqrt( 1. / ( t_p1 * t_p1 ) + 1. / ( t_p2 * t_p2 ) + 1. / ( t_p3 * t_p3 ) );
        // t_pspg = 0;
    }

    //t_pspg = 0;
}




void
Quad10_2D_SUPG :: computeAdvectionTerm(FloatMatrix &answer, TimeStep *atTime)
{
    FloatMatrix n, b;
    double dV, rho;
    int k, undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    GaussPoint *gp;

    answer.resize(undofs, undofs);
    answer.zero();

    IntegrationRule *iRule = this->integrationRulesArray [ 1 ];
    /* consistent part + supg stabilization term */
    for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        this->computeNuMatrix(n, gp);
        this->computeUDotGradUMatrix(b, gp, atTime);
        dV  = this->computeVolumeAround(gp);
        rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, gp, atTime);
        answer.plusProductUnsym(n, b, rho * dV);
    }
}


void
Quad10_2D_SUPG :: computeAdvectionDeltaTerm(FloatMatrix &answer, TimeStep *atTime)
{
    FloatMatrix n, b;
    double dV, rho;
    int k, undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    GaussPoint *gp;

    answer.resize(undofs, undofs);
    answer.zero();

    IntegrationRule *iRule = this->integrationRulesArray [ 1 ];
    /* consistent part + supg stabilization term */
    for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        this->computeNuMatrix(n, gp);
        this->computeUDotGradUMatrix(b, gp, atTime);
        dV  = this->computeVolumeAround(gp);
        rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, gp, atTime);

        answer.plusProductUnsym(b, b, rho * dV);
    }
}



void
Quad10_2D_SUPG :: computeMassDeltaTerm(FloatMatrix &answer, TimeStep *atTime)
{
    FloatMatrix n, b;
    double dV, rho;
    int k, undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    GaussPoint *gp;

    answer.resize(undofs, undofs);
    answer.zero();
    IntegrationRule *iRule = this->integrationRulesArray [ 1 ];
    /* mtrx for computing t_supg, norm of this mtrx is computed */
    for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        this->computeNuMatrix(n, gp);
        this->computeUDotGradUMatrix(b, gp, atTime);
        dV  = this->computeVolumeAround(gp);
        rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, gp, atTime);

        answer.plusProductUnsym(b, n, rho * dV);
    }
}

void
Quad10_2D_SUPG :: computeLSICTerm(FloatMatrix &answer, TimeStep *atTime)
{
    int k, undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    double dV, rho;
    GaussPoint *gp;
    FloatMatrix b;

    answer.resize(undofs, undofs);
    answer.zero();

    IntegrationRule *iRule = this->integrationRulesArray [ 1 ];
    for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        dV  = this->computeVolumeAround(gp);
        rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, gp, atTime);
        this->computeDivUMatrix(b, gp);

        answer.plusProductSymmUpper(b, b, dV * rho);
    }

    answer.symmetrized();
}


void
Quad10_2D_SUPG :: computeAdvectionEpsilonTerm(FloatMatrix &answer, TimeStep *atTime)
{
    //  to compute t_pspg
    int pndofs = this->computeNumberOfDofs(EID_ConservationEquation);
    int undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    GaussPoint *gp;
    FloatMatrix g, b;
    double dV;
    int k;

    answer.resize(pndofs, undofs);
    answer.zero();

    IntegrationRule *iRule = this->integrationRulesArray [ 1 ];

    for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        this->computeGradPMatrix(g, gp);
        this->computeUDotGradUMatrix(b, gp, atTime);
        dV  = this->computeVolumeAround(gp);

        answer.plusProductUnsym(g, b, dV);
    }
}

void
Quad10_2D_SUPG :: computeMassEpsilonTerm(FloatMatrix &answer, TimeStep *atTime)
{
    // to compute t_pspg
    int pndofs = this->computeNumberOfDofs(EID_ConservationEquation);
    int undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    int k;
    double dV;
    FloatMatrix g, n;
    GaussPoint *gp;

    answer.resize(pndofs, undofs);
    answer.zero();

    IntegrationRule *iRule = this->integrationRulesArray [ 1 ];

    for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        this->computeGradPMatrix(g, gp);
        this->computeNuMatrix(n, gp);
        dV  = this->computeVolumeAround(gp);

        answer.plusProductUnsym(g, n, dV);
    }
}

int
Quad10_2D_SUPG :: giveNumberOfSpatialDimensions()
{
    return 2;
}



double
Quad10_2D_SUPG :: LS_PCS_computeF(LevelSetPCS *ls, TimeStep *atTime)
{
    /*
     * int i;
     * double answer;
     * FloatArray fi(3), un(6);
     *
     * this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, un);
     * for ( i = 1; i <= 3; i++ ) {
     *    fi.at(i) = ls->giveLevelSetDofManValue( dofManArray.at(i) );
     * }
     *
     * double fix = b [ 0 ] * fi.at(1) + b [ 1 ] * fi.at(2) + b [ 2 ] * fi.at(3);
     * double fiy = c [ 0 ] * fi.at(1) + c [ 1 ] * fi.at(2) + c [ 2 ] * fi.at(3);
     * double norm = sqrt(fix * fix + fiy * fiy);
     *
     * answer = ( 1. / 3. ) * ( fix * ( un.at(1) + un.at(3) + un.at(5) ) + fiy * ( un.at(2) + un.at(4) + un.at(6) ) ) / norm;
     * return answer;
     */
    return 0.0;
}


double
Quad10_2D_SUPG :: LS_PCS_computeS(LevelSetPCS *ls, TimeStep *atTime)
{
    return 0.0;
}


void
Quad10_2D_SUPG :: LS_PCS_computedN(FloatMatrix &answer)
{ }


void
Quad10_2D_SUPG :: LS_PCS_computeVOFFractions(FloatArray &answer, FloatArray &fi)
{ }



double
Quad10_2D_SUPG :: computeCriticalTimeStep(TimeStep *tStep)
{
    return 1.e6;
}



int
Quad10_2D_SUPG :: ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type)
{
    if ( ( type == IST_StressTensor ) || ( type == IST_StrainTensor ) ) {
        return 4;
    }

    GaussPoint *gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    return this->giveIPValueSize(type, gp);
}


void
Quad10_2D_SUPG :: NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node,
                                                          InternalStateType type, TimeStep *tStep)
{
    GaussPoint *gp;
    gp = integrationRulesArray [ 0 ]->getIntegrationPoint(0);
    this->giveIPValue(answer, gp, type, tStep);
}


void
Quad10_2D_SUPG :: NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side,
                                                         InternalStateType type, TimeStep *tStep)
{
    answer.resize(0);
}


int
Quad10_2D_SUPG :: checkConsistency()
{
    return SUPGElement :: checkConsistency();
}


void
Quad10_2D_SUPG :: updateYourself(TimeStep *tStep)
{
    SUPGElement :: updateYourself(tStep);
}

int
Quad10_2D_SUPG :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
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
Quad10_2D_SUPG :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type)
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
Quad10_2D_SUPG :: giveIPValueType(InternalStateType type)
{
    if ( ( type == IST_VOFFraction ) || ( type == IST_Density ) ) {
        return ISVT_SCALAR;
    } else {
        return SUPGElement :: giveIPValueType(type);
    }
}


int
Quad10_2D_SUPG :: giveIPValueSize(InternalStateType type, GaussPoint *gp)
{
    if ( ( type == IST_VOFFraction ) || ( type == IST_Density ) ) {
        return 1;
    } else {
      return SUPGElement::giveIPValueSize(type, gp);
    }
}


//void
//Quad10_2D_SUPG :: printOutputAt(FILE *file, TimeStep *stepN)
// Performs end-of-step operations.
//{
//    SUPGElement :: printOutputAt(file, stepN);
//}




contextIOResultType
Quad10_2D_SUPG :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full element context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;

    if ( ( iores = SUPGElement :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}



contextIOResultType Quad10_2D_SUPG :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full element context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;

    if ( ( iores = SUPGElement :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


double
Quad10_2D_SUPG :: computeVolumeAround(GaussPoint *aGaussPoint)
// Returns the portion of the receiver which is attached to aGaussPoint.
{
    double determinant, weight, volume;

    determinant = fabs( this->velocityInterpolation.giveTransformationJacobian(* aGaussPoint->giveCoordinates(), FEIElementGeometryWrapper(this)) );


    weight      = aGaussPoint->giveWeight();
    volume      = determinant * weight;

    return volume;
}

#if 0
double
Quad10_2D_SUPG :: computeVolumeAroundPressure(FEInterpolation2d& interpol, GaussPoint *aGaussPoint)
// Returns the portion of the receiver which is attached to aGaussPoint.
{
    double determinant, weight, volume;

    determinant = fabs( interpol.giveTransformationJacobian(domain, pressureDofManArray, * aGaussPoint->giveCoordinates(), 0.0) );

    weight = aGaussPoint->giveWeight();
    volume = determinant * weight;

    return volume;
}
#endif

Interface *
Quad10_2D_SUPG :: giveInterface(InterfaceType interface)
{
    if ( interface == LevelSetPCSElementInterfaceType ) {
        return ( LevelSetPCSElementInterface * ) this;
    } else if ( interface == ZZNodalRecoveryModelInterfaceType )  {
        return ( ZZNodalRecoveryModelInterface * ) this;
    } else if ( interface == NodalAveragingRecoveryModelInterfaceType )  {
        return ( NodalAveragingRecoveryModelInterface * ) this;
    }

    return NULL;
}


void
Quad10_2D_SUPG :: printOutputAt(FILE *file, TimeStep *stepN)
// Performs end-of-step operations.
{
    int i;

#ifdef __PARALLEL_MODE
    fprintf( file, "element %d [%8d] :\n", this->giveNumber(), this->giveGlobalNumber() );
#else
    fprintf(file, "element %d :\n", number);
#endif
    pressureNode.printOutputAt(file, stepN);

    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        integrationRulesArray [ i ]->printOutputAt(file, stepN);
    }
}



void
Quad10_2D_SUPG :: giveLocalVelocityDofMap(IntArray &map)
{
    map.resize(8);
    int i;

    for ( i = 1; i <= 8; i++ ) {
        map.at(i) = i;
    }
}
void
Quad10_2D_SUPG :: giveLocalPressureDofMap(IntArray &map)
{
    map.resize(1);
    map.at(1) = 9;
}



#ifdef __OOFEG
int
Quad10_2D_SUPG :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                       int node, TimeStep *atTime)
{
    return SUPGElement :: giveInternalStateAtNode(answer, type, mode, node, atTime);
}



void
Quad10_2D_SUPG :: drawRawGeometry(oofegGraphicContext &gc)
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

void Quad10_2D_SUPG :: drawScalar(oofegGraphicContext &context)
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
        //this->formMaterialVolumePoly(matvolpoly, NULL, temp_normal, temp_p, false);
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
