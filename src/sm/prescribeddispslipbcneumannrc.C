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
 *               Copyright (C) 1993 - 2021   Borek Patzak
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

#include "prescribeddispslipbcneumannrc.h"
#include "classfactory.h"
#include "node.h"
#include "masterdof.h"
#include "element.h"
#include "feinterpol.h"
#include "feinterpol2d.h"
#include "gausspoint.h"
#include "sparsemtrx.h"
#include "timestep.h"
#include "function.h"
#include "sparselinsystemnm.h"
#include "unknownnumberingscheme.h"
#include "engngm.h"
#include "mathfem.h"
#include "crosssection.h"
#include "Elements/structuralelement.h"


namespace oofem {
REGISTER_BoundaryCondition(PrescribedDispSlipBCNeumannRC);

PrescribedDispSlipBCNeumannRC :: PrescribedDispSlipBCNeumannRC(int n, Domain *d) :
    ActiveBoundaryCondition(n, d),
    PrescribedDispSlipHomogenization(),
    mpSigmaHom( new Node(0, d) ),
    lmTauHom( new Node(0, d) ),
    lmSigmaSHom( new Node(0, d) )
{
    int nsd = d->giveNumberOfSpatialDimensions();
    for ( int i = 0; i < nsd * nsd; i++ ) {
        // Just putting in X_i id-items since they don't matter.
        int dofId = d->giveNextFreeDofID();
        mSigmaIds.followedBy(dofId);
        mpSigmaHom->appendDof( new MasterDof( mpSigmaHom.get(), ( DofIDItem ) ( dofId ) ) );
    }

    for ( int i = 0; i < nsd; i++ ) {
        // Just putting in X_i id-items since they don't matter.
        int dofId = d->giveNextFreeDofID();
        lmTauIds.followedBy(dofId);
        lmTauHom->appendDof( new MasterDof( lmTauHom.get(), ( DofIDItem ) ( dofId ) ) );
    }

    for ( int i = 0; i < nsd ; i++ ) {
        // Just putting in X_i id-items since they don't matter.
        int dofId = d->giveNextFreeDofID();
        lmSigmaSIds.followedBy(dofId);
        lmSigmaSHom->appendDof( new MasterDof( lmSigmaSHom.get(), ( DofIDItem ) ( dofId ) ) );
    }
}

PrescribedDispSlipBCNeumannRC :: ~PrescribedDispSlipBCNeumannRC()
{
}


void PrescribedDispSlipBCNeumannRC :: initializeFrom(InputRecord &ir)
{
    ActiveBoundaryCondition :: initializeFrom(ir);
    PrescribedDispSlipHomogenization :: initializeFrom(ir);

    if ( this->dispGradient.isNotEmpty() ) {
        this->dispGradON = true;
    }

    if ( this->slipField.isNotEmpty() ) {
        this->slipON = true;
    }

    if ( this->slipGradient.isNotEmpty() ) {
        this->slipGradON = true;
    }

    if ( slipON || slipGradON) { //either slip or slip gradient (or both) will be prescribed weakly
        IR_GIVE_OPTIONAL_FIELD( ir, concreteVolSet, _IFT_PrescribedDispSlipBCNeumannRC_ConcreteVolSet );
        IR_GIVE_OPTIONAL_FIELD( ir, rebarSets, _IFT_PrescribedDispSlipBCNeumannRC_RebarSets );
    }
}


void PrescribedDispSlipBCNeumannRC :: giveInputRecord(DynamicInputRecord &input)
{
    ActiveBoundaryCondition :: giveInputRecord(input);
    PrescribedDispSlipHomogenization :: giveInputRecord(input);

    if ( slipON || slipGradON ) {
        input.setField(this->concreteVolSet, _IFT_PrescribedDispSlipBCNeumannRC_ConcreteVolSet);
        input.setField(this->rebarSets, _IFT_PrescribedDispSlipBCNeumannRC_RebarSets);
    }
}


int PrescribedDispSlipBCNeumannRC::giveNumberOfInternalDofManagers()
{
    int lmTot = (int)dispGradON + (int)slipON + (int)slipGradON;
    return lmTot;
}

DofManager *PrescribedDispSlipBCNeumannRC :: giveInternalDofManager(int i)
{
    if ( this->giveNumberOfInternalDofManagers() == 3 ) {
        if ( i == 1 ) { return mpSigmaHom.get(); }
        else if ( i == 2 ) { return lmTauHom.get(); }
        else if ( i == 3) { return lmSigmaSHom.get(); }
    } else if ( this->giveNumberOfInternalDofManagers() == 2) {
        if ( dispGradON && slipON ) {
            if ( i == 1 ) { return mpSigmaHom.get(); }
            else if ( i == 2 ) { return lmTauHom.get(); }
        } else if ( dispGradON && slipGradON ) {
            if ( i == 1 ) { return mpSigmaHom.get(); }
            else if ( i ==2 ) { return lmSigmaSHom.get(); }
        } else if ( slipON && slipGradON ) {
            if ( i == 1 ) { return lmTauHom.get(); }
            else if ( i == 2 ) { return lmSigmaSHom.get(); }
        }
    } else if ( this->giveNumberOfInternalDofManagers() == 1 ) {
        if ( dispGradON ) { return mpSigmaHom.get(); }
        else if ( slipON ) { return lmTauHom.get(); }
        else if ( slipGradON ) {return lmSigmaSHom.get(); }
    }

    return nullptr;
}

void PrescribedDispSlipBCNeumannRC :: scale(double s)
{
    this->dispGradient.times(s);
    this->slipField.times(s);
    this->slipGradient.times(s);
}


void PrescribedDispSlipBCNeumannRC :: assembleVector(FloatArray &answer, TimeStep *tStep,
                                                     CharType type, ValueModeType mode,
                                                     const UnknownNumberingScheme &s,
                                                     FloatArray *eNorm,
                                                     void* lock)
{
    if ( dispGradON ) {
        this->assembleVectorStress( answer, tStep, type, mode, s, eNorm );
    }

    if ( slipON ) {
        this->assembleVectorBStress(answer, tStep, type, mode, s);
    }

    if ( slipGradON ) {
        this->assembleVectorRStress( answer, tStep, type, mode, s);
    }
}


void PrescribedDispSlipBCNeumannRC::assembleVectorStress( FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode, const UnknownNumberingScheme &s, FloatArray *eNorm )
{
    IntArray sigma_loc;  // For the displacements and stress respectively
    mpSigmaHom->giveLocationArray(mSigmaIds, sigma_loc, s);

    if ( type == ExternalForcesVector ) {
        // The external forces have two contributions. On the additional equations for sigma, the load is simply the prescribed gradient.
        double rve_size = PrescribedDispSlipHomogenization::domainSize(this->giveDomain(), this->giveSetNumber());
        FloatArray stressLoad;
        FloatArray dispGrad;
        giveDispGradient(dispGrad);

        double loadLevel = this->giveTimeFunction()->evaluateAtTime(tStep->giveTargetTime());
        stressLoad.beScaled(-rve_size*loadLevel, dispGrad);

        answer.assemble(stressLoad, sigma_loc);
    } else if ( type == InternalForcesVector ) {
        FloatMatrix Ke;
        FloatArray fe_v, fe_s;
        FloatArray sigmaHom, e_u;
        IntArray loc, masterDofIDs, sigmaMasterDofIDs;

        // Fetch the current values of the stress;
        mpSigmaHom->giveUnknownVector(sigmaHom, mSigmaIds, mode, tStep);
        // and the master dof ids for sigmadev used for the internal norms
        mpSigmaHom->giveMasterDofIDArray(mSigmaIds, sigmaMasterDofIDs);

        // Assemble
        Set *setPointer = this->giveDomain()->giveSet(this->set);
        const IntArray &boundaries = setPointer->giveBoundaryList();
        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
            int boundary = boundaries.at(pos * 2);

            // Fetch the element information;
            e->giveLocationArray(loc, s, & masterDofIDs);
            e->computeVectorOf(mode, tStep, e_u);
            this->integrateTangentStress(Ke, e, boundary);

            // We just use the tangent, less duplicated code (the addition of sigmaDev is linear).
            fe_v.beProductOf(Ke, e_u);
            fe_s.beTProductOf(Ke, sigmaHom);

            // Note: The terms appear negative in the equations:
            fe_v.negated();
            fe_s.negated();

            answer.assemble(fe_s, loc); // Contributions to delta_v equations
            answer.assemble(fe_v, sigma_loc); // Contribution to delta_s_i equations
            if ( eNorm != NULL ) {
                eNorm->assembleSquared(fe_s, masterDofIDs);
                eNorm->assembleSquared(fe_v, sigmaMasterDofIDs);
            }
        }
    }
}


void PrescribedDispSlipBCNeumannRC::assembleVectorBStress( FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode, const UnknownNumberingScheme &s )
{
    IntArray tau_loc;
    lmTauHom->giveLocationArray(lmTauIds, tau_loc, s);

    if ( type == ExternalForcesVector ) {
        // The external forces have two contributions. On the additional equations for tau, the load is equal to the prescribed slip value.
        FloatArray bStressLoad;
        FloatArray slipVec;
        this->giveSlipField(slipVec);

        double loadLevel = this->giveTimeFunction()->evaluateAtTime(tStep->giveTargetTime());
        bStressLoad.beScaled(loadLevel, slipVec);

        answer.assemble(bStressLoad, tau_loc);
    } else if ( type == InternalForcesVector ) {
        FloatMatrix Ktau, Ktau_el;
        FloatArray fe_tau, fe_u;
        FloatArray tau, u_s, u_c;
        IntArray loc, loc_s, loc_c, masterDofIDs;

        // Fetch the current values of the Lagrange multiplier;
        lmTauHom->giveUnknownVector(tau, lmTauIds, mode, tStep);

        // CONTRIBUTION FROM REINFORCEMENT
        FloatMatrix C;
        double gammaBoxInt = computeInterfaceLength(rebarSets);
        this->computeWeightMatrix(C, rebarSets);
        Ktau.clear();
        Ktau_el.clear();
        fe_tau.clear();
        fe_u.clear();

        for (int i = 0; i < rebarSets.giveSize(); ++i) {
            Set *steelSet = this->giveDomain()->giveSet(rebarSets.at(i+1));

            for (int pos = 1; pos <= steelSet->giveElementList().giveSize(); ++pos) {
                Element *es = this->giveDomain()->giveElement( steelSet->giveElementList().at(pos) );

                // Fetch the element information;
                es->giveLocationArray(loc_s, s, &masterDofIDs);

                // Fetch the nodal displacements
                // since computeVectorOf gives the unknown vector in local coordinate system, it needs to be rotated to global coordinates
                FloatMatrix G2L;
                es->computeGtoLRotationMatrix(G2L);
                es->computeVectorOf(mode, tStep, u_s);
                u_s.rotatedWith(G2L, 't');

                //Compute the stiffness matrix expansion
                this->integrateTangentBStressSteel(Ktau_el, es, rebarSets.at(i+1));
                Ktau.beTProductOf(C, Ktau_el);
                Ktau.times(1/gammaBoxInt);

                //Compute the contribution to internal force vector
                fe_u.beTProductOf(Ktau, tau);
                fe_tau.beProductOf(Ktau, u_s);

                //Assemble
                answer.assemble(fe_u, loc_s);
                answer.assemble(fe_tau, tau_loc);
            }
        }
        Ktau.clear();
        fe_tau.clear();
        fe_u.clear();

        // CONTRIBUTION FROM CONCRETE
        Set *concreteSet = this->giveDomain()->giveSet(concreteVolSet);
        double omegaBox = PrescribedDispSlipHomogenization::domainSize(this->giveDomain(), this->giveSetNumber());

        for ( int pos = 1; pos <= concreteSet->giveElementList().giveSize(); ++pos) {
            Element* ec = this->giveDomain()->giveElement( concreteSet->giveElementList().at(pos) );

            // Fetch the element information;
            ec->giveLocationArray(loc_c, s, &masterDofIDs);

            // Fetch the nodal displacements
            ec->computeVectorOf(mode, tStep, u_c);

            //Compute the stiffness matrix expansion
            this->integrateTangentBStressConcrete(Ktau, ec);
            Ktau.times(1/omegaBox);

            //Compute the contribution to internal force vector
            fe_u.beTProductOf(Ktau, tau);
            fe_tau.beProductOf(Ktau, u_c);

            // Note: The terms appear negative in the equations:
            fe_u.negated();
            fe_tau.negated();

            //Assemble
            answer.assemble(fe_u, loc_c);
            answer.assemble(fe_tau, tau_loc);
        }
    }
}


void PrescribedDispSlipBCNeumannRC::assembleVectorRStress( FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode, const UnknownNumberingScheme &s )
{
    IntArray sigmaS_loc;
    lmSigmaSHom->giveLocationArray(lmSigmaSIds, sigmaS_loc, s);

    if ( type == ExternalForcesVector ) {
        // The external forces have two contributions. On the additional equations for sigmaS, the load is equal to the prescribed slip gradient value.
        FloatArray rStressLoad;
        FloatArray slipGrad;
        this->giveSlipGradient(slipGrad);
        slipGrad.resizeWithValues(2);

        double loadLevel = this->giveTimeFunction()->evaluateAtTime(tStep->giveTargetTime());
        rStressLoad.beScaled(loadLevel, slipGrad);

        answer.assemble(rStressLoad, sigmaS_loc);
    } else if ( type == InternalForcesVector ) {
        FloatMatrix Ksig, Ksig_el;
        FloatArray fe_sig, fe_u;
        FloatArray sigmaS, u_s, u_c;
        IntArray loc, loc_s, loc_c, masterDofIDs;

        // Fetch the current values of the Lagrange multiplier;
        lmSigmaSHom->giveUnknownVector(sigmaS, lmSigmaSIds, mode, tStep);

        // CONTRIBUTION FROM REINFORCEMENT
        FloatMatrix C;
        double gammaBoxInt = computeInterfaceLength(rebarSets);
        this->computeWeightMatrix(C, rebarSets);
//        C.resizeWithData(4,4);
        Ksig.clear();
        Ksig_el.clear();
        fe_sig.clear();
        fe_u.clear();

        for ( int i = 0; i < rebarSets.giveSize(); ++i ) {
            Set *steelSet = this->giveDomain()->giveSet(rebarSets.at(i+1));

            for (int pos = 1; pos <= steelSet->giveElementList().giveSize(); ++pos) {
                Element *es = this->giveDomain()->giveElement( steelSet->giveElementList().at(pos) );

                // Fetch the element information;
                es->giveLocationArray(loc_s, s, &masterDofIDs);

                // Fetch the nodal displacements
                // since computeVectorOf gives the unknown vector in local coordinate system, it needs to be rotated to global coordinates
                FloatMatrix G2L;
                es->computeGtoLRotationMatrix(G2L);
                es->computeVectorOf(mode, tStep, u_s);
                u_s.rotatedWith(G2L, 't');

                //Compute the stiffness matrix expansion
                this->integrateTangentRStressSteel(Ksig_el, es, rebarSets.at(i+1));
                Ksig.beTProductOf(C, Ksig_el);
                Ksig.times(1/gammaBoxInt);

                //Compute the contribution to internal force vector
                fe_u.beTProductOf(Ksig, sigmaS);
                fe_sig.beProductOf(Ksig, u_s);

                //Assemble
                answer.assemble(fe_u, loc_s);
                answer.assemble(fe_sig, sigmaS_loc);
            }
        }

        Ksig.clear();
        fe_sig.clear();
        fe_u.clear();

        // CONTRIBUTION FROM CONCRETE
        double omegaBox = PrescribedDispSlipHomogenization::domainSize(this->giveDomain(), this->giveSetNumber());
        Set *concreteBoundSet = this->giveDomain()->giveSet(this->set);
        const IntArray &boundaries = concreteBoundSet->giveBoundaryList();

        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *ec = this->giveDomain()->giveElement( boundaries.at(pos*2-1) );
            int boundary = boundaries.at(pos*2);

            // Fetch the element information;
            ec->giveLocationArray(loc_c, s, &masterDofIDs);

            // Fetch the nodal displacements
            ec->computeVectorOf(mode, tStep, u_c);

            //Compute the stiffness matrix expansion
            this->integrateTangentRStressConcrete(Ksig, ec, boundary);
            Ksig.times(1/omegaBox);

            //Compute the contribution to internal force vector
            fe_u.beTProductOf(Ksig, sigmaS);
            fe_sig.beProductOf(Ksig, u_c);

            // Note: The terms appear negative in the equations:
            fe_u.negated();
            fe_sig.negated();

            //Assemble
            answer.assemble(fe_u, loc_c);
            answer.assemble(fe_sig, sigmaS_loc);
        }
    }
}

void PrescribedDispSlipBCNeumannRC :: assemble(SparseMtrx &answer, TimeStep *tStep,
                                             CharType type, 
                                             const UnknownNumberingScheme &r_s, 
                                             const UnknownNumberingScheme &c_s, 
                                             double scale,
                                             void* lock)
{
    if ( type == TangentStiffnessMatrix || type == SecantStiffnessMatrix || type == ElasticStiffnessMatrix ) {
        if ( dispGradON ) {
            this->assembleOnStress(answer, r_s, c_s, scale );
        }

        if ( slipON ) {
            this->assembleOnTransferStress(answer, r_s, c_s, scale);
        }

        if ( slipGradON ) {
            this->assembleOnReinfStress(answer, r_s, c_s, scale);
        }
    } else   {
        OOFEM_LOG_DEBUG("Skipping assembly in PrescribedDispSlipBCNeumann::assemble().");
    }
}



void PrescribedDispSlipBCNeumannRC :: giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                                                       const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    int i = 0;
    if ( dispGradON ) {
        IntArray loc_r, loc_c, sigma_loc_r, sigma_loc_c;

        // Fetch the columns/rows for the stress contributions;
        mpSigmaHom->giveLocationArray( mSigmaIds, sigma_loc_r, r_s );
        mpSigmaHom->giveLocationArray( mSigmaIds, sigma_loc_c, c_s );

        Set *set                   = this->giveDomain()->giveSet( this->set );
        const IntArray &boundaries = set->giveBoundaryList();

        rows.resize( boundaries.giveSize() );
        cols.resize( boundaries.giveSize() );
        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *e = this->giveDomain()->giveElement( boundaries.at( pos * 2 - 1 ) );

            // Here, we could use only the nodes actually located on the boundary, but we don't.
            // Instead, we use all nodes belonging to the element, which is allowed because the
            // basis functions related to the interior nodes will be zero on the boundary.
            // Obviously, this is less efficient, so why do we want to do it this way?
            // Because it is easier when XFEM enrichments are present. /ES
            e->giveLocationArray( loc_r, r_s );
            e->giveLocationArray( loc_c, c_s );

            // For most uses, loc_r == loc_c, and sigma_loc_r == sigma_loc_c.
            rows[i] = loc_r;
            cols[i] = sigma_loc_c;
            i++;
            // and the symmetric part (usually the transpose of above)
            rows[i] = sigma_loc_r;
            cols[i] = loc_c;
            i++;
        }
    }

    if ( slipON ) {
        IntArray loc_r, loc_c, tau_loc_r, tau_loc_c;

        // Fetch the columns/rows for the transfer stress contributions;
        lmTauHom->giveLocationArray(lmTauIds, tau_loc_r, r_s);
        lmTauHom->giveLocationArray(lmTauIds, tau_loc_c, c_s);

        Set *concreteSet = this->giveDomain()->giveSet(this->concreteVolSet);
        const IntArray &conEl = concreteSet->giveElementList();
        IntArray steelEl;

        for ( int j = 1; j <= this->rebarSets.giveSize(); ++j ) {
            steelEl.followedBy(this->giveDomain()->giveSet(rebarSets.at(j))->giveElementList());
        }

        rows.resize( 2*steelEl.giveSize() + 2*conEl.giveSize() + i);
        cols.resize( 2*steelEl.giveSize() + 2*conEl.giveSize() + i);
        //Contribution from steel
        for ( int pos = 1; pos <= steelEl.giveSize(); ++pos ) {
            Element *e = this->giveDomain()->giveElement(steelEl.at(pos));

            e->giveLocationArray(loc_r, r_s);
            e->giveLocationArray(loc_c, c_s);

            // For most uses, loc_r == loc_c, and tau_loc_r == tau_loc_c.
            rows [ i ] = loc_r;
            cols [ i ] = tau_loc_c;
            i++;
            // and the symmetric part (usually the transpose of above)
            rows [ i ] = tau_loc_r;
            cols [ i ] = loc_c;
            i++;
        }

        //Contribution from concrete
        for ( int pos = 1; pos <= conEl.giveSize(); ++pos ) {
            Element *e = this->giveDomain()->giveElement(conEl.at(pos));

            e->giveLocationArray(loc_r, r_s);
            e->giveLocationArray(loc_c, c_s);

            // For most uses, loc_r == loc_c, and tau_loc_r == tau_loc_c.
            rows [ i ] = loc_r;
            cols [ i ] = tau_loc_c;
            i++;
            // and the symmetric part (usually the transpose of above)
            rows [ i ] = tau_loc_r;
            cols [ i ] = loc_c;
            i++;
        }
    }

    if ( slipGradON ) {
        IntArray loc_r, loc_c, sigmaS_loc_r, sigmaS_loc_c;

        // Fetch the columns/rows for the reinforcement stress contributions;
        lmSigmaSHom->giveLocationArray(lmSigmaSIds, sigmaS_loc_r, r_s);
        lmSigmaSHom->giveLocationArray(lmSigmaSIds, sigmaS_loc_c, c_s);

        Set *concreteSet = this->giveDomain()->giveSet(this->set);
        const IntArray &boundaries = concreteSet->giveBoundaryList();
        IntArray steelEl;

        for ( int i = 1; i <= this->rebarSets.giveSize(); ++i ) {
            steelEl.followedBy(this->giveDomain()->giveSet(rebarSets.at(i))->giveElementList());
        }

        rows.resize( 2*steelEl.giveSize() + boundaries.giveSize() + i );
        cols.resize( 2*steelEl.giveSize() + boundaries.giveSize() + i );
        //Contribution from steel
        for ( int pos = 1; pos <= steelEl.giveSize(); ++pos ) {
            Element *e = this->giveDomain()->giveElement(steelEl.at(pos));

            e->giveLocationArray(loc_r, r_s);
            e->giveLocationArray(loc_c, c_s);

            // For most uses, loc_r == loc_c, and sigmaS_loc_r == sigmaS_loc_c.
            rows [ i ] = loc_r;
            cols [ i ] = sigmaS_loc_c;
            i++;
            // and the symmetric part (usually the transpose of above)
            rows [ i ] = sigmaS_loc_r;
            cols [ i ] = loc_c;
            i++;
        }

        //Contribution from concrete
        for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
            Element *e = this->giveDomain()->giveElement( boundaries.at(pos*2-1) );

            e->giveLocationArray(loc_r, r_s);
            e->giveLocationArray(loc_c, c_s);

            // For most uses, loc_r == loc_c, and sigmaS_loc_r == sigmaS_loc_c.
            rows [ i ] = loc_r;
            cols [ i ] = sigmaS_loc_c;
            i++;
            // and the symmetric part (usually the transpose of above)
            rows [ i ] = sigmaS_loc_r;
            cols [ i ] = loc_c;
            i++;
        }
    }
}

void PrescribedDispSlipBCNeumannRC :: computeStress(FloatArray &sigma, TimeStep *tStep)
{
    if ( dispGradON ) {
        double volRVE  = this->domainSize( this->giveDomain(), this->giveSetNumber() );
        double areaRVE = PrescribedDispSlipHomogenization::domainSize( this->giveDomain(), this->giveSetNumber() );
        double thick   = volRVE / areaRVE;
        mpSigmaHom->giveUnknownVector( sigma, mSigmaIds, VM_Total, tStep );
        sigma.times( 1 / thick );
    }
}


void PrescribedDispSlipBCNeumannRC::computeTransferStress( FloatArray &bStress, TimeStep *tStep )
{
    if ( slipON ) {
        double volRVE = this->domainSize( this->giveDomain(), this->giveSetNumber() );
        lmTauHom->giveUnknownVector( bStress, lmTauIds, VM_Total, tStep );
        bStress.times( -1 / volRVE );
    }
}


void PrescribedDispSlipBCNeumannRC::computeReinfStress( FloatArray &rStress, TimeStep *tStep )
{
    if ( slipGradON ) {
        double volRVE = this->domainSize( this->giveDomain(), this->giveSetNumber() );
        lmSigmaSHom->giveUnknownVector( rStress, lmSigmaSIds, VM_Total, tStep );
        rStress.resizeWithValues(4);
        rStress.times( -1 / volRVE );
    }
}


void PrescribedDispSlipBCNeumannRC :: computeTangent(FloatMatrix &tangent, TimeStep *tStep)
{
    OOFEM_ERROR("Tangent not implemented")
}


void PrescribedDispSlipBCNeumannRC :: integrateTangentStress(FloatMatrix &oTangent, Element *e, int iBndIndex)
{
    FloatArray normal, n;
    FloatMatrix nMatrix, E_n;
    FloatMatrix contrib;

    FEInterpolation *interp = e->giveInterpolation(); // Geometry interpolation

    int nsd = e->giveDomain()->giveNumberOfSpatialDimensions();

    // Interpolation order
    int order = interp->giveInterpolationOrder();
    std :: unique_ptr< IntegrationRule > ir;

    ir = interp->giveBoundaryIntegrationRule(order, iBndIndex, e->giveGeometryType());

    oTangent.clear();

    for ( auto &gp: *ir ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        FEIElementGeometryWrapper cellgeo(e);

        // Evaluate the normal;
        double detJ = interp->boundaryEvalNormal(normal, iBndIndex, lcoords, cellgeo);

        // Compute global coordinates of Gauss point
        FloatArray globalCoord;
        interp->boundaryLocal2Global(globalCoord, iBndIndex, lcoords, cellgeo);

        // Compute local coordinates on the element
        FloatArray bulkElLocCoords;
        e->computeLocalCoordinates(bulkElLocCoords, globalCoord);

        // Evaluate the velocity/displacement coefficients
        interp->evalN(n, bulkElLocCoords, cellgeo);
        nMatrix.beNMatrixOf(n, nsd);

        E_n.resize(nsd*nsd, nsd);
        E_n.zero();

        if ( nsd == 3 ) {
            E_n.at(1, 1) = normal.at(1);
            E_n.at(2, 2) = normal.at(2);
            E_n.at(3, 3) = normal.at(3);
            E_n.at(4, 1) = normal.at(2);
            E_n.at(5, 1) = normal.at(3);
            E_n.at(6, 2) = normal.at(3);
            E_n.at(7, 2) = normal.at(1);
            E_n.at(8, 3) = normal.at(1);
            E_n.at(9, 3) = normal.at(2);
        } else if ( nsd == 2 ) {
            E_n.at(1, 1) = normal.at(1);
            E_n.at(2, 2) = normal.at(2);
            E_n.at(3, 1) = normal.at(2);
            E_n.at(4, 2) = normal.at(1);
        } else {
            E_n.at(1, 1) = normal.at(1);
        }

        contrib.beProductOf(E_n, nMatrix);

        oTangent.add(detJ * gp->giveWeight(), contrib);
    }
}


void PrescribedDispSlipBCNeumannRC::integrateTangentBStressSteel( FloatMatrix &oTangent, Element *e, const int &rebarSet )
{
    FloatMatrix Nstau, Ns, contrib;
    double Ss;
    IntArray DofMask;

    //Get number of DOFs per DofManager
    e->giveDofManDofIDMask(1, DofMask);
    int ndof = DofMask.giveSize();

    FEInterpolation *interp = e->giveInterpolation();
    int order = interp->giveInterpolationOrder();
    std :: unique_ptr< IntegrationRule > ir;
    ir = interp->giveIntegrationRule(order, e->giveGeometryType());

    oTangent.clear();

    for ( auto &gp : *ir ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        FEIElementGeometryWrapper cellgeo(e);

        //Compute shape functions
        FloatArray nu;
        interp->evalN(nu, lcoords, cellgeo);
        Ns.beNMatrixOf(nu, ndof);

        //Constant traction interpolation
        this->computeRebarDyad(Nstau, rebarSet);
        Nstau.resizeWithData(3,2);

        //Compute rebar perimeter
        CrossSection *cs = e->giveCrossSection();
        Ss = sqrt(4 * cs->give(CS_Area, gp) * M_PI);

        //Compute the integral
        contrib.beTProductOf(Nstau, Ns);
        contrib.times(Ss);
        double detJ = interp->giveTransformationJacobian(lcoords, cellgeo);
        oTangent.add(detJ * gp->giveWeight(), contrib);
    }
}


void PrescribedDispSlipBCNeumannRC::integrateTangentBStressConcrete( FloatMatrix &oTangent, Element *e )
{
    FloatMatrix Nctau, Nc, contrib;
    IntArray DofMask;

    //Get number of DOFs per DofManager
    e->giveDofManDofIDMask(1, DofMask);
    int ndof = DofMask.giveSize();

    FEInterpolation *interp = e->giveInterpolation();
    int order = interp->giveInterpolationOrder();
    std :: unique_ptr< IntegrationRule > ir;
    ir = interp->giveIntegrationRule(order, e->giveGeometryType());

    oTangent.clear();
    Nctau.resize(ndof,ndof);

    for ( auto &gp : *ir ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        FEIElementGeometryWrapper cellgeo(e);

        //Compute shape functions
        FloatArray nu;
        interp->evalN(nu, lcoords, cellgeo);
        Nc.beNMatrixOf(nu, ndof);

        //Constant traction interpolation
        Nctau.beUnitMatrix();

        //Compute the integral
        contrib.beTProductOf(Nctau, Nc);
        double detJ = interp->giveTransformationJacobian(lcoords, cellgeo);
        oTangent.add(detJ * gp->giveWeight(), contrib);
    }
}


void PrescribedDispSlipBCNeumannRC::integrateTangentRStressSteel( FloatMatrix &oTangent, Element *e, const int &rebarSet )
{
    FloatMatrix Nssig, Bs, contrib;
    double Ss;
    IntArray DofMask;

    //Get number of DOFs per DofManager
    e->giveDofManDofIDMask(1, DofMask);

    FEInterpolation *interp = e->giveInterpolation();
    int order = interp->giveInterpolationOrder();
    std :: unique_ptr< IntegrationRule > ir;
    ir = interp->giveIntegrationRule(order, e->giveGeometryType());
    //Cast into StructuralElement
    StructuralElement *se = dynamic_cast<StructuralElement*>(e);

    oTangent.clear();

    for ( auto &gp : *ir ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        FEIElementGeometryWrapper cellgeo(e);

        //Compute shape function derivatives
        FloatMatrix b;
        se->computeBmatrixAt(gp, b);
        //TODO: make it element independent; now hardcoded for libeam2d elements
        Bs.resize(2,6);
        Bs.at(1,1) = b.at(1,1);
        Bs.at(2,2) = b.at(1,1);
        Bs.at(1,4) = b.at(1,4);
        Bs.at(2,5) = b.at(1,4);

        //Constant traction interpolation
        this->computeRebarDyad(Nssig, rebarSet);
//        Nssig.resizeWithData(2,4);

        //Compute rebar perimeter
        CrossSection *cs = e->giveCrossSection();
        Ss = sqrt(4 * cs->give(CS_Area, gp) * M_PI);

        //Compute the integral
        contrib.beTProductOf(Nssig, Bs);
        contrib.times(Ss);
        double detJ = interp->giveTransformationJacobian(lcoords, cellgeo);
        oTangent.add(detJ * gp->giveWeight(), contrib);
    }
}


void PrescribedDispSlipBCNeumannRC::integrateTangentRStressConcrete( FloatMatrix &oTangent, Element *e, int iBndIndex )
{
    FloatArray normal, n;
    FloatMatrix nMatrix, E_n;
    FloatMatrix contrib;

    FEInterpolation *interp = e->giveInterpolation();

    int nsd = e->giveDomain()->giveNumberOfSpatialDimensions();

    //Interpolation order
    int order = interp->giveInterpolationOrder();
    std :: unique_ptr< IntegrationRule > ir = interp->giveBoundaryEdgeIntegrationRule(order, iBndIndex, e->giveGeometryType());

    oTangent.clear();

    for ( auto &gp : *ir ) {
        const FloatArray &lcoords = gp->giveNaturalCoordinates();
        FEIElementGeometryWrapper cellgeo(e);

        //Evaluate the normal
        double detJ = interp->boundaryEvalNormal(normal, iBndIndex, lcoords, cellgeo);

        //Compute global coordinates of Gauss point
        FloatArray globalCoord;
        interp->boundaryLocal2Global(globalCoord, iBndIndex, lcoords, cellgeo);

        //Compute local coordinates on the element
        FloatArray bulkElLocCoords;
        e->computeLocalCoordinates(bulkElLocCoords, globalCoord);

        //Evaluate the shape functions
        interp->evalN(n, bulkElLocCoords, cellgeo);
        nMatrix.beNMatrixOf(n, nsd);

        E_n.resize(nsd, nsd);
        E_n.zero();
        if ( nsd == 3 ) {
            E_n.at(1, 1) = normal.at(1);
            E_n.at(2, 2) = normal.at(2);
            E_n.at(3, 3) = normal.at(3);
        } else if ( nsd == 2 ) {
            E_n.at(1, 1) = normal.at(1);
            E_n.at(2,2) = normal.at(2);
        } else {
            E_n.at(1, 1) = normal.at(1);
        }

        contrib.beProductOf(E_n, nMatrix);

        oTangent.add(detJ * gp->giveWeight(), contrib);
    }

}


void PrescribedDispSlipBCNeumannRC::assembleOnStress( SparseMtrx &answer, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale )
{
    FloatMatrix Ke, KeT;
    IntArray loc_r, loc_c, sigma_loc_r, sigma_loc_c;
    Set *set = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = set->giveBoundaryList();

    // Fetch the columns/rows for the stress contributions;
    mpSigmaHom->giveLocationArray(mSigmaIds, sigma_loc_r, r_s);
    mpSigmaHom->giveLocationArray(mSigmaIds, sigma_loc_c, c_s);


    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *e = this->giveDomain()->giveElement( boundaries.at(pos * 2 - 1) );
        int boundary = boundaries.at(pos * 2);

        e->giveLocationArray(loc_r, r_s);
        e->giveLocationArray(loc_c, c_s);

        this->integrateTangentStress(Ke, e, boundary);

        Ke.negated();
        Ke.times(scale);
        KeT.beTranspositionOf(Ke);

        answer.assemble(sigma_loc_r, loc_c, Ke); // Contribution to delta_s_i equations
        answer.assemble(loc_r, sigma_loc_c, KeT); // Contributions to delta_v equations
    }
}


void PrescribedDispSlipBCNeumannRC::assembleOnTransferStress( SparseMtrx &answer, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale )
{
    FloatMatrix Ke, KeT;
    FloatMatrix Ktau;
    IntArray loc_r, loc_c, tau_loc_r, tau_loc_c;

    // Fetch the columns/rows for transferr stress contributions;
    lmTauHom->giveLocationArray(lmTauIds, tau_loc_r, r_s);
    lmTauHom->giveLocationArray(lmTauIds, tau_loc_c, c_s);

    //Contribution from reinforcement
    FloatMatrix C;
    double gammaBoxInt = computeInterfaceLength(rebarSets);
    this->computeWeightMatrix(C, rebarSets);

    for ( int i = 0; i < rebarSets.giveSize(); ++i ) {
        Set *steelSet = this->giveDomain()->giveSet(rebarSets.at(i+1));

        for ( int pos = 1; pos <= steelSet->giveElementList().giveSize(); ++pos ) {
            Element *es = this->giveDomain()->giveElement( steelSet->giveElementList().at(pos) );

            es->giveLocationArray(loc_r, r_s);
            es->giveLocationArray(loc_c, c_s);

            this->integrateTangentBStressSteel(Ktau, es, rebarSets.at(i+1));
            Ke.beTProductOf(C, Ktau);
            Ke.times(1/gammaBoxInt);
            KeT.beTranspositionOf(Ke);

            answer.assemble(tau_loc_r, loc_c, Ke);
            answer.assemble(loc_r, tau_loc_c, KeT);
        }
    }

    Ktau.clear();
    Ke.clear();
    KeT.clear();

    //Contribution from concrete
    Set *concreteSet = this->giveDomain()->giveSet(concreteVolSet);
    double omegaBox = PrescribedDispSlipHomogenization::domainSize(this->giveDomain(), this->giveSetNumber());

    for ( int pos = 1; pos <= concreteSet->giveElementList().giveSize(); ++pos ) {
        Element *ec = this->giveDomain()->giveElement( concreteSet->giveElementList().at(pos) );

        ec->giveLocationArray(loc_r, r_s);
        ec->giveLocationArray(loc_c, c_s);

        this->integrateTangentBStressConcrete(Ke, ec);
        Ke.times(1/omegaBox);
        Ke.negated();
        KeT.beTranspositionOf(Ke);

        answer.assemble(tau_loc_r, loc_c, Ke);
        answer.assemble(loc_r, tau_loc_c, KeT);
    }
}


void PrescribedDispSlipBCNeumannRC::assembleOnReinfStress( SparseMtrx &answer, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale )
{
    FloatMatrix Ke, KeT;
    FloatMatrix Ksig;
    IntArray loc_r, loc_c, sigmaS_loc_r, sigmaS_loc_c;

    // Fetch the columns/rows for the stress contributions;
    lmSigmaSHom->giveLocationArray(lmSigmaSIds, sigmaS_loc_r, r_s);
    lmSigmaSHom->giveLocationArray(lmSigmaSIds, sigmaS_loc_c, c_s);

    //Contribution from reinforcement
    FloatMatrix C;
    double gammaBoxInt = computeInterfaceLength(rebarSets);
    this->computeWeightMatrix(C, rebarSets);
//    C.resizeWithData(4,4);

    for ( int i = 0; i < rebarSets.giveSize(); ++i ) {
        Set *steelSet = this->giveDomain()->giveSet(rebarSets.at(i+1));

        for ( int pos = 1; pos <= steelSet->giveElementList().giveSize(); ++pos ) {
            Element *es = this->giveDomain()->giveElement( steelSet->giveElementList().at(pos) );

            es->giveLocationArray(loc_r, r_s);
            es->giveLocationArray(loc_c, c_s);

            this->integrateTangentRStressSteel(Ksig, es, rebarSets.at(i+1));
            Ke.beTProductOf(C, Ksig);
            Ke.times(1/gammaBoxInt);
            KeT.beTranspositionOf(Ke);

            answer.assemble(sigmaS_loc_r, loc_c, Ke);
            answer.assemble(loc_r, sigmaS_loc_c, KeT);
        }
    }

    Ksig.clear();
    Ke.clear();
    KeT.clear();

    //Contribution from concrete
    Set *concreteBoundSet = this->giveDomain()->giveSet(this->set);
    const IntArray &boundaries = concreteBoundSet->giveBoundaryList();
    double omegaBox = PrescribedDispSlipHomogenization::domainSize(this->giveDomain(), this->giveSetNumber());

    for ( int pos = 1; pos <= boundaries.giveSize() / 2; ++pos ) {
        Element *ec = this->giveDomain()->giveElement( boundaries.at(pos*2-1) );
        int boundary = boundaries.at(pos*2);

        ec->giveLocationArray(loc_r, r_s);
        ec->giveLocationArray(loc_c, c_s);

        this->integrateTangentRStressConcrete(Ke, ec, boundary);

        Ke.times(1/omegaBox);
        Ke.negated();
        KeT.beTranspositionOf(Ke);

        answer.assemble(sigmaS_loc_r, loc_c, Ke);
        answer.assemble(loc_r, sigmaS_loc_c, KeT);
    }
}


void PrescribedDispSlipBCNeumannRC::computeWeightMatrix( FloatMatrix &C, const IntArray &reinfSets )
{
    ///Assumption: the rebars cross the entire RVE, i.e. they do not end within the RVE
    /// - the rebars have circular cross section
    double gammaBoxInt = computeInterfaceLength(reinfSets);
    FloatMatrix dyadMatrix;
    for ( int i = 0; i < reinfSets.giveSize(); ++i ) {
        FloatMatrix dya;
        double perimeterCoeff = 0;

        Set *set = this->giveDomain()->giveSet(reinfSets.at(i+1));
        for ( int j = 0; j < set->giveElementList().giveSize(); ++j ) {
            Element *e = this->giveDomain()->giveElement(set->giveElementList().at(j+1));
            FEInterpolation *interp = e->giveInterpolation();
            int order = interp->giveInterpolationOrder();
            std :: unique_ptr< IntegrationRule > ir;
            ir = interp->giveIntegrationRule(order, e->giveGeometryType());

            for ( auto &gp : *ir ) {
                const FloatArray &lcoords = gp->giveNaturalCoordinates();
                FEIElementGeometryWrapper cellgeo(e);
                double detJ = interp->giveTransformationJacobian(lcoords, cellgeo);
                CrossSection *cs = e->giveCrossSection();
                double perim = sqrt(4 * cs->give(CS_Area, gp) * M_PI);
                perimeterCoeff = perimeterCoeff + perim * detJ * gp->giveWeight();
            }
        }

        this->computeRebarDyad(dya, reinfSets.at(i+1));
        dyadMatrix.add(perimeterCoeff, dya);
    }

    C.beInverseOf(dyadMatrix);
    C.times(gammaBoxInt);
}


void PrescribedDispSlipBCNeumannRC::computeRebarDyad( FloatMatrix &dyad, const int &reinfSet )
{
    ///Assuming that the reinforcement bar remains straight throughout the RVE
    ///it is then enough to compute the tangent vector for the first element

    FloatArray e_l;
    Set *set = this->giveDomain()->giveSet(reinfSet);
    IntArray ellist = set->giveElementList();
    Element *element = this->giveDomain()->giveElement(ellist.at(1));
    FEIElementGeometryWrapper cellgeo(element);

    e_l.resize(2);
    e_l.at(1) = cellgeo.giveVertexCoordinates(2).at(1) - cellgeo.giveVertexCoordinates(1).at(1);
    e_l.at(2) = cellgeo.giveVertexCoordinates(2).at(2) - cellgeo.giveVertexCoordinates(1).at(2);
    e_l.normalize();

    dyad.beDyadicProductOf(e_l,e_l);
}


double PrescribedDispSlipBCNeumannRC::computeInterfaceLength( const IntArray &reinfSets )
{
    double gamma = 0.;
    for ( int i=0; i < reinfSets.giveSize(); ++i ) {
        Set *set = this->giveDomain()->giveSet(reinfSets.at(i+1));
        for ( int j=0; j < set->giveElementList().giveSize(); ++j ) {
            Element *e = this->giveDomain()->giveElement(set->giveElementList().at(j+1));
            FEInterpolation *interp = e->giveInterpolation();
            int order = interp->giveInterpolationOrder();
            std :: unique_ptr< IntegrationRule > ir;
            ir = interp->giveIntegrationRule(order, e->giveGeometryType());

            for ( auto &gp : *ir ) {
                const FloatArray &lcoords = gp->giveNaturalCoordinates();
                FEIElementGeometryWrapper cellgeo(e);
                double detJ = interp->giveTransformationJacobian(lcoords, cellgeo);
                gamma = gamma + detJ * gp->giveWeight();
            }
        }
    }

    return gamma;
}


double PrescribedDispSlipBCNeumannRC::domainSize( Domain *d, int set )
{
    double omegaBox = PrescribedDispSlipHomogenization::domainSize(d, this->giveSetNumber());

    if ( this->giveDomain()->giveNumberOfSpatialDimensions() == 2 ) {
        //assuming that the RVE thickness is constant in 2D
        Element *e = this->giveDomain()->giveElement( this->giveDomain()->giveSet( this->giveSetNumber() )->giveBoundaryList().at(1) );
        std::unique_ptr<IntegrationRule> ir = e->giveInterpolation()->giveIntegrationRule( e->giveInterpolation()->giveInterpolationOrder(), e->giveGeometryType() );
        CrossSection *cs = e->giveCrossSection();
        GaussPoint *gp = ir->getIntegrationPoint(0);
        double thickness = cs->give(CS_Thickness, gp);
        return omegaBox * thickness;
    } else {
        return omegaBox;
    }
}

} /* namespace oofem */
