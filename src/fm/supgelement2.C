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

#include "supgelement2.h"
#include "domain.h"
#include "timestep.h"
#include "load.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "fluidmodel.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif

namespace oofem {
SUPGElement2 :: SUPGElement2(int n, Domain *aDomain) :
    SUPGElement(n, aDomain)
    // Constructor. Creates an element with number n, belonging to aDomain.
{}


SUPGElement2 :: ~SUPGElement2()
// Destructor.
{ }

IRResultType
SUPGElement2 :: initializeFrom(InputRecord *ir)
{
    //const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    //IRResultType result;                               // Required by IR_GIVE_FIELD macro

    return SUPGElement :: initializeFrom(ir);
}


void
SUPGElement2 :: giveInputRecord(DynamicInputRecord &input)
{
    SUPGElement :: giveInputRecord(input);
}


void
SUPGElement2 :: giveCharacteristicMatrix(FloatMatrix &answer,
                                         CharType mtrx, TimeStep *tStep)
//
// returns characteristics matrix of receiver according to mtrx
//
{
#if 0
    if ( mtrx == StiffnessMatrix ) {
        // support for stokes solver
        IntArray vloc, ploc;
        FloatMatrix h;
        int size = this->computeNumberOfDofs();
        this->giveLocalVelocityDofMap(vloc);
        this->giveLocalPressureDofMap(ploc);
        answer.resize(size, size);
        answer.zero();
        this->computeAdvectionDerivativeTerm_MB(h, tStep);
        answer.assemble(h, vloc);
        this->computeDiffusionDerivativeTerm_MB(h, TangentStiffness, tStep);
        answer.assemble(h, vloc);
        this->computePressureTerm_MB(h, tStep);
        answer.assemble(h, vloc, ploc);
        this->computeLinearAdvectionTerm_MC(h, tStep);
        answer.assemble(h, ploc, vloc);
        this->computeAdvectionDerivativeTerm_MC(h, tStep);
        answer.assemble(h, ploc, vloc);
        this->computeDiffusionDerivativeTerm_MC(h, tStep);
        answer.assemble(h, ploc, vloc);
        this->computePressureTerm_MC(h, tStep);
        answer.assemble(h, ploc);
        this->computeLSICStabilizationTerm_MB(h, tStep);
        answer.assemble(h, vloc);
    } else
#endif
    {
        _error("giveCharacteristicMatrix: Unknown Type of characteristic mtrx.");
    }
}


void
SUPGElement2 :: giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
                                         TimeStep *tStep)
//
// returns characteristics vector of receiver according to requested type
//
{
    if ( mtrx == ExternalForcesVector ) {
        // stokes flow
        IntArray vloc, ploc;
        FloatArray h;
        int size = this->computeNumberOfDofs();
        this->giveLocalVelocityDofMap(vloc);
        this->giveLocalPressureDofMap(ploc);
        answer.resize(size);
        answer.zero();
        this->computeBCRhsTerm_MB(h, tStep);
        answer.assemble(h, vloc);
        this->computeBCRhsTerm_MC(h, tStep);
        answer.assemble(h, ploc);
    } else
#if 0
    if ( mtrx == InternalForcesVector ) {
        // stokes flow
        IntArray vloc, ploc;
        FloatArray h;
        int size = this->computeNumberOfDofs();
        this->giveLocalVelocityDofMap(vloc);
        this->giveLocalPressureDofMap(ploc);
        answer.resize(size);
        answer.zero();
        this->computeAdvectionTerm_MB(h, tStep);
        answer.assemble(h, vloc);
        this->computeAdvectionTerm_MC(h, tStep);
        answer.assemble(h, ploc);
        this->computeDiffusionTerm_MB(h, tStep);
        answer.assemble(h, vloc);
        this->computeDiffusionTerm_MC(h, tStep);
        answer.assemble(h, ploc);

        FloatMatrix m1;
        FloatArray v, p;
        // add lsic stabilization term
        this->giveCharacteristicMatrix(m1, LSICStabilizationTerm_MB, tStep);
        //m1.times( lscale / ( dscale * uscale * uscale ) );
        this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, v);
        h.beProductOf(m1, v);
        answer.assemble(h, vloc);
        this->giveCharacteristicMatrix(m1, LinearAdvectionTerm_MC, tStep);
        //m1.times( 1. / ( dscale * uscale ) );
        h.beProductOf(m1, v);
        answer.assemble(h, ploc);


        // add pressure term
        this->giveCharacteristicMatrix(m1, PressureTerm_MB, tStep);
        this->computeVectorOf(EID_ConservationEquation, VM_Total, tStep, v);
        h.beProductOf(m1, v);
        answer.assemble(h, vloc);

        // pressure term
        this->giveCharacteristicMatrix(m1, PressureTerm_MC, tStep);
        this->computeVectorOf(EID_ConservationEquation, VM_Total, tStep, v);
        h.beProductOf(m1, v);
        answer.assemble(h, ploc);
    } else
#endif
    {
        _error("giveCharacteristicVector: Unknown Type of characteristic mtrx.");
    }
}


double
SUPGElement2 :: giveCharacteristicValue(CharType mtrx, TimeStep *tStep)
{
    if ( mtrx == CriticalTimeStep ) {
        return this->computeCriticalTimeStep(tStep);
    } else {
        _error("giveCharacteristicValue: Unknown Type of characteristic mtrx.");
    }

    return 0.0;
}


int
SUPGElement2 :: checkConsistency()
//
// check internal consistency
// mainly tests, whether material and crossSection data
// are safe for conversion to "Structural" versions
//
{
    int result = 1;
    /*
     * if (!this->giveMaterial()->testMaterialExtension(Material_TransportCapability)) {
     * _warning("checkConsistency : material without support for transport problems");
     * result =0;
     * }
     */
    return result;
}


void
SUPGElement2 :: updateInternalState(TimeStep *tStep)
{
    FloatArray stress;

    // force updating strains & stresses
    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        IntegrationRule *iRule = integrationRulesArray [ i ];
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            computeDeviatoricStress(stress, iRule->getIntegrationPoint(j), tStep);
        }
    }
}

void
SUPGElement2 :: printOutputAt(FILE *file, TimeStep *tStep)
// Performs end-of-step operations.
{
    int i;

#ifdef __PARALLEL_MODE
    fprintf( file, "element %d [%8d] :\n", this->giveNumber(), this->giveGlobalNumber() );
#else
    fprintf(file, "element %d :\n", number);
#endif

    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        integrationRulesArray [ i ]->printOutputAt(file, tStep);
    }
}


#ifdef __OOFEG
int
SUPGElement2 :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                        int node, TimeStep *tStep)
{
    int indx = 1;
    Node *n = this->giveNode(node);

    if ( type == IST_Velocity ) {
        answer.resize( this->giveSpatialDimension() );
        int dofindx;
        if ( ( dofindx = n->findDofWithDofId(V_u) ) ) {
            answer.at(indx++) = n->giveDof(dofindx)->giveUnknown(VM_Total, tStep);
        }

        if ( ( dofindx = n->findDofWithDofId(V_v) ) ) {
            answer.at(indx++) = n->giveDof(dofindx)->giveUnknown(VM_Total, tStep);
        }

        if ( ( dofindx = n->findDofWithDofId(V_w) ) ) {
            answer.at(indx++) = n->giveDof(dofindx)->giveUnknown(VM_Total, tStep);
        }

        return 1;
    } else if ( type == IST_Pressure ) {
        int dofindx;
        if ( ( dofindx = n->findDofWithDofId(P_f) ) ) {
            answer.resize(1);
            answer.at(1) = n->giveDof(dofindx)->giveUnknown(VM_Total, tStep);
            return 1;
        } else {
            return 0;
        }
    } else {
        return Element :: giveInternalStateAtNode(answer, type, mode, node, tStep);
    }
}

#endif

void
SUPGElement2 :: computeAccelerationTerm_MB(FloatMatrix &answer, TimeStep *tStep)
{
    FloatMatrix n, b;

    answer.resize(0, 0);

    int rule = 2;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];
    /* consistent part + supg stabilization term */
    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(k);
        this->computeNuMatrix(n, gp);
        this->computeUDotGradUMatrix( b, gp, tStep->givePreviousStep() );
        double dV = this->computeVolumeAround(gp);
        double rho = this->giveMaterial()->give('d', gp);
        /* consistent part */
        answer.plusProductUnsym(n, n, rho * dV);
        /* supg stabilization */
        answer.plusProductUnsym(b, n, rho * t_supg * dV);
    }
}

void
SUPGElement2 :: computeAdvectionTerm_MB(FloatArray &answer, TimeStep *tStep)
{
    FloatMatrix n, b, bn;
    FloatArray u, v(3);

    answer.resize(0);

    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);

    int rule = 2;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];
    /* consistent part + supg stabilization term */
    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(k);
        this->computeNuMatrix(n, gp);
        this->computeUDotGradUMatrix( bn, gp, tStep->givePreviousStep() );
        this->computeUDotGradUMatrix(b, gp, tStep);
        v.beProductOf(b, u);
        double dV = this->computeVolumeAround(gp);
        double rho = this->giveMaterial()->give('d', gp);
        /* consistent part */
        answer.plusProduct(n, v, rho * dV);

        /* supg stabilization */
        answer.plusProduct(bn, v, t_supg * rho * dV);
    }
}

void
SUPGElement2 :: computeAdvectionDerivativeTerm_MB(FloatMatrix &answer, TimeStep *tStep)
{
    FloatMatrix n, b, bn, grad_u, grad_uN, N;

    answer.resize(0, 0);

    int rule = 2;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];
    /* consistent part + supg stabilization term */
    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(k);
        this->computeNuMatrix(n, gp);
        this->computeUDotGradUMatrix( bn, gp, tStep->givePreviousStep() );
        this->computeUDotGradUMatrix(b, gp, tStep);
        double dV  = this->computeVolumeAround(gp);
        double rho = this->giveMaterial()->give('d', gp);

        this->computeGradUMatrix(grad_u, gp, tStep);

        /* consistent part */
        answer.plusProductUnsym(n, b, rho * dV);

        grad_uN.beProductOf(grad_u, n);
        answer.plusProductUnsym(n, grad_uN, rho * dV);
        /* supg stabilization */
        answer.plusProductUnsym(bn, b, t_supg * rho * dV);
        answer.plusProductUnsym(bn, grad_uN, t_supg * rho * dV);
    }
}

void
SUPGElement2 :: computeDiffusionTerm_MB(FloatArray &answer, TimeStep *tStep)
{
    FloatArray u, eps, stress, bs, dDB_u;
    FloatMatrix b, un_gu, dDB;
    double Re = static_cast< FluidModel * >( domain->giveEngngModel() )->giveReynoldsNumber();

    answer.resize(0);

    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);

    int rule = 1;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];
    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(k);
        double dV = this->computeVolumeAround(gp);
        this->computeBMatrix(b, gp);
        this->computeDivTauMatrix(dDB, gp, tStep);
        this->computeUDotGradUMatrix( un_gu, gp, tStep->givePreviousStep() );
        eps.beProductOf(b, u);
        static_cast< FluidDynamicMaterial * >( this->giveMaterial() )->computeDeviatoricStressVector(stress, gp, eps, tStep);
        dDB_u.beProductOf(dDB, u);
        /* consistent part */
        answer.plusProduct(b, stress, dV / Re);

        /* SUPG term */
        answer.plusProduct(un_gu, dDB_u, ( -1.0 ) * t_supg * dV);
    }
}

void
SUPGElement2 :: computeDiffusionDerivativeTerm_MB(FloatMatrix &answer, MatResponseMode mode, TimeStep *tStep)
{
    FloatMatrix _db, _d, _b, dDB, un_gu;
    double Re = static_cast< FluidModel * >( domain->giveEngngModel() )->giveReynoldsNumber();
    FloatArray dDB_u;

    answer.resize(0, 0);

    int rule = 1;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];
    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(k);
        double dV = this->computeVolumeAround(gp);
        this->computeBMatrix(_b, gp);
        static_cast< FluidDynamicMaterial * >( this->giveMaterial() )->giveDeviatoricStiffnessMatrix(_d, mode, gp, tStep);

        this->computeDivTauMatrix(dDB, gp, tStep);
        this->computeUDotGradUMatrix( un_gu, gp, tStep->givePreviousStep() );
        /* standard term */
        _db.beProductOf(_d, _b);
        answer.plusProductUnsym(_b, _db, dV / Re); //answer.plusProduct (_b,_db,area);

        /* SUPG term */
        answer.plusProductUnsym( un_gu, dDB, t_supg * dV * ( -1.0 ) * ( 1. / Re ) );
    }
}


void
SUPGElement2 :: computePressureTerm_MB(FloatMatrix &answer, TimeStep *tStep)
{
    FloatMatrix gu, np, b;

    answer.resize(0, 0);

    int rule = 1;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];
    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(k);
        double dV = this->computeVolumeAround(gp);
        this->computeDivUMatrix(gu, gp);
        this->computeNpMatrix(np, gp);

        /*alternative computing*/
        //this->computeNuMatrix(gu, gp);
        //this->computeGradPMatrix(np, gp);

        /* standard term */
        answer.plusProductUnsym(gu, np, ( -1.0 ) * dV);
    }

    iRule = this->integrationRulesArray [ 1 ];
    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(k);
        double dV = this->computeVolumeAround(gp);
        this->computeUDotGradUMatrix( b, gp, tStep->givePreviousStep() );
        this->computeGradPMatrix(np, gp);

        // supg term
        answer.plusProductUnsym(b, np, t_supg * dV);
    }
}
void
SUPGElement2 :: computeLSICStabilizationTerm_MB(FloatMatrix &answer, TimeStep *tStep)
{
    FloatMatrix b;

    answer.resize(0, 0);

    int rule = 0;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];
    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(k);
        double dV = this->computeVolumeAround(gp);
        double rho = this->giveMaterial()->give('d', gp);
        this->computeDivUMatrix(b, gp);

        answer.plusProductSymmUpper(b, b, dV * rho * t_lsic);
    }

    answer.symmetrized();
}

void
SUPGElement2 :: computeLinearAdvectionTerm_MC(FloatMatrix &answer, TimeStep *tStep)
{
    FloatMatrix gu, np;

    answer.resize(0, 0);

    int rule = 0;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];
    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(k);
        double dV = this->computeVolumeAround(gp);
        this->computeDivUMatrix(gu, gp);
        this->computeNpMatrix(np, gp);

        /* standard term */
        answer.plusProductUnsym(np, gu, dV);
    }
}

void
SUPGElement2 :: computeAdvectionTerm_MC(FloatArray &answer, TimeStep *tStep)
{
    // N_epsilon (due to PSPG stabilization)
    FloatMatrix g, b;
    FloatArray u, v;

    answer.resize(0);

    int rule = 1;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];
    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);

    /* pspg stabilization term */
    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(k);
        this->computeGradPMatrix(g, gp);
        this->computeUDotGradUMatrix(b, gp, tStep);
        v.beProductOf(b, u);
        double dV = this->computeVolumeAround(gp);
        answer.plusProduct(g, v, t_pspg * dV);
    }
}


void
SUPGElement2 :: computeAdvectionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *tStep)
{
    FloatMatrix g, b;

    answer.resize(0, 0);

    int rule = 1;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];
    /* pspg stabilization term */
    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(k);
        this->computeGradPMatrix(g, gp);
        this->computeUDotGradUMatrix(b, gp, tStep);
        double dV = this->computeVolumeAround(gp);
        answer.plusProductUnsym(g, b, dV * t_pspg);
    }
}

void
SUPGElement2 :: computeDiffusionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *tStep)
{
    FloatMatrix dDB, _d, g;

    answer.resize(0, 0);

    int rule = 1;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];

    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(k);
        double dV = this->computeVolumeAround(gp);
        double rho = this->giveMaterial()->give('d', gp);

        //( ( FluidDynamicMaterial * ) this->giveMaterial() )->giveDeviatoricStiffnessMatrix(_d, TangentStiffness,gp, tStep);

        this->computeDivTauMatrix(dDB, gp, tStep);
        this->computeGradPMatrix(g, gp);

        answer.plusProductUnsym(g, dDB, ( -1.0 ) * dV * t_pspg / rho);
    }

    answer.zero(); ///@todo Is this really correct? It makes the whole loop meaningless. If it is, please describe why with a comment.
}

void
SUPGElement2 :: computeDiffusionTerm_MC(FloatArray &answer, TimeStep *tStep)
{
    answer.resize(0);

    /*
     * FloatMatrix dDB, _d, g;
     * double sum, dV, coeff, rho, isd, nsd = this->giveNumberOfSpatialDimensions();
     * GaussPoint *gp;
     * IntegrationRule *iRule = this->integrationRulesArray [ 1 ];
     * int i,k;
     * FloatArray u, dDB_u;
     *
     * this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);
     *
     * for ( k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
     *  gp = iRule->getIntegrationPoint(k);
     *  dV  = this->computeVolumeAround(gp);
     *  rho = this->giveMaterial()->give('d', gp);
     *
     *  coeff = (-1.0) * dV * t_pspg / rho;
     *
     *  this->computeGradPMatrix(g, gp);
     *  this->computeDivTauMatrix(dDB, gp, tStep);
     *
     *  dDB_u.beProductOf(dDB, u);
     *
     *  for ( i = 1; i <= pndofs; i++ ) {
     *      for ( sum = 0, isd = 1; isd <= nsd; isd++ ) {
     *          sum += g.at(isd, i) * dDB_u.at(isd);
     *      }
     *
     *      answer.at(i) += coeff * sum;
     *  }
     *
     *
     *
     * }
     */
}


void
SUPGElement2 :: computeAccelerationTerm_MC(FloatMatrix &answer, TimeStep *tStep)
{
    FloatMatrix g, n;

    answer.resize(0, 0);
    // pspg stabilization term: M_\epsilon term


    int rule = 0;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];

    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(k);
        this->computeGradPMatrix(g, gp);
        this->computeNuMatrix(n, gp);
        double dV = this->computeVolumeAround(gp);
        answer.plusProductUnsym(g, n, dV * t_pspg);
    }
}

void
SUPGElement2 :: computePressureTerm_MC(FloatMatrix &answer, TimeStep *tStep)
{
    FloatMatrix g;

    answer.resize(0, 0);

    int rule = 0;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];

    for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
        GaussPoint *gp = iRule->getIntegrationPoint(k);
        this->computeGradPMatrix(g, gp);
        double dV  = this->computeVolumeAround(gp);
        double rho = this->giveMaterial()->give('d', gp);
        answer.plusProductSymmUpper(g, g, dV * t_pspg / rho);
    }

    answer.symmetrized();
}


void
SUPGElement2 :: computeBCRhsTerm_MB(FloatArray &answer, TimeStep *tStep)
{
    int nLoads;

    answer.resize(0);

    int rule = 0;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];
    FloatArray un, gVector, s, helpLoadVector;
    FloatMatrix b, nu;

    // add body load (gravity) termms
    nLoads = this->giveBodyLoadArray()->giveSize();
    for ( int i = 1; i <= nLoads; i++ ) {
        Load *load = domain->giveLoad( bodyLoadArray.at(i) );
        bcGeomType ltype = load->giveBCGeoType();
        if ( ( ltype == BodyLoadBGT ) && ( load->giveBCValType() == ForceLoadBVT ) ) {
            load->computeComponentArrayAt(gVector, tStep, VM_Total);
            if ( gVector.giveSize() ) {
                for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
                    GaussPoint *gp = iRule->getIntegrationPoint(k);
                    this->computeUDotGradUMatrix( b, gp, tStep->givePreviousStep() );
                    this->computeNuMatrix(nu, gp);
                    double dV  = this->computeVolumeAround(gp);
                    double rho = this->giveMaterial()->give('d', gp);
                    answer.plusProduct(b, gVector, t_supg * rho * dV);
                    answer.plusProduct(nu, gVector, rho * dV);
                }
            }
        }
    }

    // integrate tractions
    // if no traction bc applied but side marked as with traction load
    // then zero traction is assumed !!!

    // loop over boundary load array
    nLoads = this->giveBoundaryLoadArray()->giveSize() / 2;
    for ( int i = 1; i <= nLoads; i++ ) {
        int n = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
        int id = boundaryLoadArray.at(i * 2);
        Load *load  = domain->giveLoad(n);
        bcGeomType ltype = load->giveBCGeoType();
        if ( ltype == EdgeLoadBGT ) {
            this->computeEdgeLoadVector_MB(helpLoadVector, load, id, tStep);
            if ( helpLoadVector.giveSize() ) {
                answer.add(helpLoadVector);
            }
        } else if ( ltype == SurfaceLoadBGT ) {
            this->computeSurfaceLoadVector_MB(helpLoadVector, load, id, tStep);
            if ( helpLoadVector.giveSize() ) {
                answer.add(helpLoadVector);
            }
        } else {
            _error("computeForceLoadVector : unsupported load type class");
        }
    }
}



void
SUPGElement2 :: computeBCRhsTerm_MC(FloatArray &answer, TimeStep *tStep)
{
    int nLoads;
    FloatArray s, gVector, helpLoadVector;
    FloatMatrix g;

    int rule = 1;
    IntegrationRule *iRule = this->integrationRulesArray [ rule ];

    answer.resize(0);

    nLoads = this->giveBodyLoadArray()->giveSize();
    for ( int i = 1; i <= nLoads; i++ ) {
        Load *load  = domain->giveLoad( bodyLoadArray.at(i) );
        bcGeomType ltype = load->giveBCGeoType();
        if ( ( ltype == BodyLoadBGT ) && ( load->giveBCValType() == ForceLoadBVT ) ) {
            load->computeComponentArrayAt(gVector, tStep, VM_Total);
            if ( gVector.giveSize() ) {
                for ( int k = 0; k < iRule->giveNumberOfIntegrationPoints(); k++ ) {
                    GaussPoint *gp = iRule->getIntegrationPoint(k);
                    this->computeGradPMatrix(g, gp);
                    double dV = this->computeVolumeAround(gp);
                    answer.plusProduct(g, gVector, t_pspg * dV);
                }
            }
        }
    }

    // integrate tractions
    // if no traction bc applied but side marked as with traction load
    // then zero traction is assumed !!!

    // loop over boundary load array
    nLoads = this->giveBoundaryLoadArray()->giveSize() / 2;
    for ( int i = 1; i <= nLoads; i++ ) {
        int n = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
        int id = boundaryLoadArray.at(i * 2);
        Load *load = domain->giveLoad(n);
        bcGeomType ltype = load->giveBCGeoType();
        if ( ltype == EdgeLoadBGT ) {
            this->computeEdgeLoadVector_MC(helpLoadVector, load, id, tStep);
            if ( helpLoadVector.giveSize() ) {
                answer.add(helpLoadVector);
            }
        } else if ( ltype == SurfaceLoadBGT ) {
            this->computeSurfaceLoadVector_MC(helpLoadVector, load, id, tStep);
            if ( helpLoadVector.giveSize() ) {
                answer.add(helpLoadVector);
            }
        } else {
            _error("computeForceLoadVector : unsupported load type class");
        }
    }
}

void
SUPGElement2 :: computeEdgeLoadVector_MB(FloatArray &answer, Load *load, int id, TimeStep *tStep)
{
    _error("computeEdgeLoadVectorAt_MB: not implemented");
}

void
SUPGElement2 :: computeSurfaceLoadVector_MB(FloatArray &answer, Load *load, int id, TimeStep *tStep)
{
    _error("computeSurfaceLoadVectorAt_MB: not implemented");
}

void
SUPGElement2 :: computeEdgeLoadVector_MC(FloatArray &answer, Load *load, int id, TimeStep *tStep)
{
    _error("computeEdgeLoadVectorAt_MC: not implemented");
}

void
SUPGElement2 :: computeSurfaceLoadVector_MC(FloatArray &answer, Load *load, int id, TimeStep *tStep)
{
    _error("computeEdgeLoadVectorAt_MC: not implemented");
}


void
SUPGElement2 :: computeDeviatoricStrain(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray u;
    FloatMatrix b;
    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);

    this->computeBMatrix(b, gp);
    answer.beProductOf(b, u);
}

void
SUPGElement2 :: computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray eps;

    // compute deviatoric strain
    this->computeDeviatoricStrain(eps, gp, tStep);
    // call material to compute stress
    static_cast< FluidDynamicMaterial * >( this->giveMaterial() )->computeDeviatoricStressVector(answer, gp, eps, tStep);
}
} // end namespace oofem
