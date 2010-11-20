/* $Header: /home/cvs/bp/oofem/tm/src/transportelement.C,v 1.3.4.1 2004/04/05 15:19:53 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

#include "supgelement2.h"
#include "domain.h"
#include "timestep.h"
#include "node.h"
#include "dof.h"
#include "load.h"
#include "boundaryload.h"
#include "gausspnt.h"
#include "gaussintegrationrule.h"
#include "intarray.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "debug.h"
#include "verbose.h"

#include "elementside.h"
#include "mathfem.h"
#ifndef __MAKEDEPEND
#include <stdlib.h>
#include <stdio.h>
#endif

#ifdef __OOFEG
#include "oofeggraphiccontext.h"
#include "conTable.h"
#endif

namespace oofem {

SUPGElement2 :: SUPGElement2(int n, Domain *aDomain) :
    SUPGElement(n, aDomain)
    // Constructor. Creates an element with number n, belonging to aDomain.
{ 
  rotationMatrix     = NULL;
  rotationMatrixDefined = 0;

}


SUPGElement2 :: ~SUPGElement2()
// Destructor.
{ 

if ( rotationMatrix ) {
        delete rotationMatrix;
    }
}

IRResultType
SUPGElement2 :: initializeFrom(InputRecord *ir)
{
    //const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    //IRResultType result;                               // Required by IR_GIVE_FIELD macro

    //FMElement :: initializeFrom (ir);
    SUPGElement :: initializeFrom(ir);
    this->computeGaussPoints();
    return IRRT_OK;
}



void
SUPGElement2 ::  giveCharacteristicMatrix(FloatMatrix &answer,
                                          CharType mtrx, TimeStep *tStep)
//
// returns characteristics matrix of receiver according to mtrx
//
{
    if ( mtrx == AccelerationTerm_MB ) {
        this->computeAccelerationTerm_MB(answer, tStep);
    } else if ( mtrx == AdvectionDerivativeTerm_MB )  {
        this->computeAdvectionDerivativeTerm_MB(answer, tStep);
    } else if ( mtrx == DiffusionDerivativeTerm_MB )  {
        this->computeDiffusionDerivativeTerm_MB(answer, TangentStiffness, tStep);
    } else if ( mtrx == SecantDiffusionDerivativeTerm_MB )  {
        this->computeDiffusionDerivativeTerm_MB(answer, SecantStiffness, tStep);
    } else if ( mtrx == TangentDiffusionDerivativeTerm_MB )  {
        this->computeDiffusionDerivativeTerm_MB(answer, TangentStiffness, tStep);
    } else if ( mtrx == InitialDiffusionDerivativeTerm_MB )  {
        this->computeDiffusionDerivativeTerm_MB(answer, ElasticStiffness, tStep);
    } else if ( mtrx == PressureTerm_MB )  {
        this->computePressureTerm_MB(answer, tStep);
    } else if ( mtrx == LinearAdvectionTerm_MC )  {
        this->computeLinearAdvectionTerm_MC(answer, tStep);
    } else if ( mtrx == AdvectionDerivativeTerm_MC )  {
        this->computeAdvectionDerivativeTerm_MC(answer, tStep);
    } else if ( mtrx == DiffusionDerivativeTerm_MC )  {
        this->computeDiffusionDerivativeTerm_MC(answer, tStep);
    } else if ( mtrx == AccelerationTerm_MC )  {
        this->computeAccelerationTerm_MC(answer, tStep);
    } else if ( mtrx == PressureTerm_MC )  {
        this->computePressureTerm_MC(answer, tStep);
    } else if ( mtrx == LSICStabilizationTerm_MB )  {
        this->computeLSICStabilizationTerm_MB(answer, tStep);
    } else {
        _error("giveCharacteristicMatrix: Unknown Type of characteristic mtrx.");
    }

    return;
}


void
SUPGElement2 ::  giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,
                                          TimeStep *tStep)
//
// returns characteristics vector of receiver according to requested type
//
{
    if ( mtrx == AdvectionTerm_MB ) {
        this->computeAdvectionTerm_MB(answer, tStep);
    } else if ( mtrx == DiffusionTerm_MB )  {
        this->computeDiffusionTerm_MB(answer, tStep);
    } else if ( mtrx == AdvectionTerm_MC )  {
        this->computeAdvectionTerm_MC(answer, tStep);
    } else if ( mtrx == BCRhsTerm_MB )  {
        this->computeBCRhsTerm_MB(answer, tStep);
    } else if ( mtrx == BCRhsTerm_MC )  {
        this->computeBCRhsTerm_MC(answer, tStep);
    } else if ( mtrx == DiffusionTerm_MC )  {
        this->computeDiffusionTerm_MC(answer, tStep);
    } else {
        _error("giveCharacteristicVector: Unknown Type of characteristic mtrx.");
    }

    return;
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
SUPGElement2 :: updateInternalState(TimeStep *stepN)
{
    int i, j;
    IntegrationRule *iRule;
    FloatArray stress;

    // force updating strains & stresses
    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        iRule = integrationRulesArray [ i ];
        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
            computeDeviatoricStress(stress, iRule->getIntegrationPoint(j), stepN);
        }
    }
}

void
SUPGElement2 :: printOutputAt(FILE *file, TimeStep *stepN)
// Performs end-of-step operations.
{
    int i;

#ifdef __PARALLEL_MODE
    fprintf( file, "element %d [%8d] :\n", this->giveNumber(), this->giveGlobalNumber() );
#else
    fprintf(file, "element %d :\n", number);
#endif

    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        integrationRulesArray [ i ]->printOutputAt(file, stepN);
    }
}


#ifdef __OOFEG
int
SUPGElement2 :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                        int node, TimeStep *atTime)
{
    int indx = 1;
    Node *n = this->giveNode(node);

    if ( type == IST_Velocity ) {
        answer.resize( this->giveSpatialDimension() );
        int dofindx;
        if ( ( dofindx = n->findDofWithDofId(V_u) ) ) {
            answer.at(indx++) = n->giveDof(dofindx)->giveUnknown(EID_MomentumBalance, VM_Total, atTime);
        }

        if ( ( dofindx = n->findDofWithDofId(V_v) ) ) {
            answer.at(indx++) = n->giveDof(dofindx)->giveUnknown(EID_MomentumBalance, VM_Total, atTime);
        }

        if ( ( dofindx = n->findDofWithDofId(V_w) ) ) {
            answer.at(indx++) = n->giveDof(dofindx)->giveUnknown(EID_MomentumBalance, VM_Total, atTime);
        }

        return 1;
    } else if ( type == IST_Pressure ) {
        int dofindx;
        if ( ( dofindx = n->findDofWithDofId(P_f) ) ) {
            answer.resize(1);
            answer.at(1) = n->giveDof(dofindx)->giveUnknown(EID_ConservationEquation, VM_Total, atTime);
            return 1;
        } else {
            return 0;
        }
    } else {
        return Element :: giveInternalStateAtNode(answer, type, mode, node, atTime);
    }
}

#endif

int
SUPGElement2 :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type)
{
    if ( ( type == IST_Velocity ) ) {
        IntArray mask;
        int indx = 1;
        answer.resize(3);
        this->giveElementDofIDMask(EID_MomentumBalance, mask);
        if ( mask.findFirstIndexOf(V_u) ) {
            answer.at(1) = indx++;
        }

        if ( mask.findFirstIndexOf(V_v) ) {
            answer.at(2) = indx++;
        }

        if ( mask.findFirstIndexOf(V_w) ) {
            answer.at(3) = indx++;
        }

        return 1;
    } else if ( type == IST_Pressure ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else {
        return Element :: giveIntVarCompFullIndx(answer, type);
    }
}

/*
 * void
 * SUPGElement2::computeVectorOfPrescribed (EquationID ut, ValueModeType type, TimeStep* stepN, FloatArray& answer)
 * {
 * double scale;
 *
 * Element::computeVectorOfPrescribed (ut, type, stepN, answer);
 *
 * if (domain->giveEngngModel()->giveEquationScalingFlag()) {
 *  if (ut == EID_MomentumBalance) {
 *    scale = domain->giveEngngModel()->giveVariableScale(VST_Velocity);
 *  } else if (ut == EID_ConservationEquation) {
 *    scale = domain->giveEngngModel()->giveVariableScale(VST_Pressure);
 *  } else scale = 1.0;
 *  answer.times (1.0/scale);
 * }
 * }
 */

void
SUPGElement2 :: computeAccelerationTerm_MB(FloatMatrix &answer, TimeStep *atTime)
{
    FloatMatrix n, b;
    double dV, rho;
    int k, undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    GaussPoint *gp;

    answer.resize(undofs, undofs);
    answer.zero();
    IntegrationRule *iRule = this->integrationRulesArray [ 2 ];
    /* consistent part + supg stabilization term */
    for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        this->computeNuMatrix(n, gp);
        this->computeUDotGradUMatrix( b, gp, atTime->givePreviousStep() );
        dV  = this->computeVolumeAround(gp);
        rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, gp, atTime);
        /* consistent part */
        answer.plusProductUnsym(n, n, rho * dV);
        /* supg stabilization */
        answer.plusProductUnsym(b, n, rho * t_supg * dV);
    }

    if ( this->updateRotationMatrix() ) {
      //answer.rotatedWith(* this->rotationMatrix);
    
      FloatMatrix *Trans = this->rotationMatrix->GiveTransposition();
      answer.rotatedWith(* Trans);
      delete Trans;
      
    }

}

void
SUPGElement2 :: computeAdvectionTerm_MB(FloatArray &answer, TimeStep *atTime)
{
  FloatMatrix n, b, bn;
    FloatArray u, v(3);
    double dV, rho, coeff, sum;
    int i, k, undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    int isd, nsd = this->giveNumberOfSpatialDimensions();
    GaussPoint *gp;

    answer.resize(undofs);
    answer.zero();

    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);
    
    
    if ( this->updateRotationMatrix() ) {
      u.rotatedWith(this->rotationMatrix, 't');
    }
    
    IntegrationRule *iRule = this->integrationRulesArray [ 2 ];
    /* consistent part + supg stabilization term */
    for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
      gp = iRule->getIntegrationPoint(k);
      this->computeNuMatrix(n, gp);
      this->computeUDotGradUMatrix( bn, gp, atTime->givePreviousStep() );
      this->computeUDotGradUMatrix( b, gp, atTime);
      v.beProductOf(b, u);
      dV  = this->computeVolumeAround(gp);
      rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, gp, atTime);
      /* consistent part */
      coeff = rho * dV;
      for ( i = 1; i <= undofs; i++ ) {
	for ( sum = 0.0, isd = 1; isd <= nsd; isd++ ) {
	  sum += n.at(isd, i) * v.at(isd);
	}
	
	answer.at(i) += coeff * sum;
      }
      
      /* supg stabilization */
      coeff = t_supg * rho * dV;
      for ( i = 1; i <= undofs; i++ ) {
	for ( sum = 0, isd = 1; isd <= nsd; isd++ ) {
	  sum += bn.at(isd, i) * v.at(isd);
	}
	
	answer.at(i) += coeff * sum;
      }
      
      
      
    }
    
    if ( this->updateRotationMatrix() ) {
      answer.rotatedWith(* this->rotationMatrix, 'n');
    }
}

void
SUPGElement2 :: computeAdvectionDerivativeTerm_MB(FloatMatrix &answer, TimeStep *atTime)
{
  FloatMatrix n, b, bn, grad_u, grad_uN, N;
  FloatArray u; 
  double dV, rho;
  int k, undofs = this->computeNumberOfDofs(EID_MomentumBalance);
  GaussPoint *gp;
  
  answer.resize(undofs, undofs);
  answer.zero();
  //this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);
  IntegrationRule *iRule = this->integrationRulesArray [ 2 ];
  /* consistent part + supg stabilization term */
    for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
      gp = iRule->getIntegrationPoint(k);
      this->computeNuMatrix(n, gp);
      this->computeUDotGradUMatrix( bn, gp, atTime->givePreviousStep() );
      this->computeUDotGradUMatrix( b, gp, atTime);
      dV  = this->computeVolumeAround(gp);
      rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, gp, atTime);
      
      this->computeGradUMatrix(grad_u, gp, atTime);
     
      
      /* consistent part */
      answer.plusProductUnsym(n, b, rho * dV);
      
      grad_uN.beProductOf(grad_u, n);
      answer.plusProductUnsym(n, grad_uN, rho * dV);
      /* supg stabilization */
      answer.plusProductUnsym(bn, b, t_supg * rho * dV);
      answer.plusProductUnsym(bn, grad_uN, t_supg * rho * dV);
    }

    if ( this->updateRotationMatrix() ) {
      //answer.rotatedWith(* this->rotationMatrix);
    
      FloatMatrix *Trans = this->rotationMatrix->GiveTransposition();
      answer.rotatedWith(* Trans);
      delete Trans;



    }
}

void
SUPGElement2 :: computeDiffusionTerm_MB(FloatArray &answer, TimeStep *atTime)
{
    int k, undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    FloatArray u, eps, stress, bs, dDB_u;
    FloatMatrix b, un_gu, dDB;
    GaussPoint *gp;
    double coeff,i, sum, isd, nsd,  dV, Re = domain->giveEngngModel()->giveUnknownComponent(ReynoldsNumber, VM_Unknown, atTime, domain, NULL);

    answer.resize(undofs);
    answer.zero();
    
    nsd = this->giveNumberOfSpatialDimensions();
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);
 
    if ( this->updateRotationMatrix() ) {
        u.rotatedWith(this->rotationMatrix, 't');
    }

    IntegrationRule *iRule = this->integrationRulesArray [ 1 ];
    for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        dV  = this->computeVolumeAround(gp);
        this->computeBMatrix(b, gp);
        this->computeDivTauMatrix(dDB, gp, atTime);
	this->computeUDotGradUMatrix( un_gu, gp, atTime->givePreviousStep() );
	eps.beProductOf(b, u);
        ( ( FluidDynamicMaterial * ) this->giveMaterial() )->computeDeviatoricStressVector(stress, gp, eps, atTime);
	dDB_u.beProductOf(dDB, u);
	/* consistent part */
        stress.times(dV / Re);
        bs.beTProductOf(b, stress);
        answer.add(bs);
    
	/* SUPG term */	
	//answer.plusProductUnsym(un_gu,dDB_u, t_supg * dV * (-1.0) * (1./Re));
	
	coeff = (-1.0) * t_supg * dV;
        for ( i = 1; i <= undofs; i++ ) {
            for ( sum = 0, isd = 1; isd <= nsd; isd++ ) {
                sum += un_gu.at(isd, i) * dDB_u.at(isd);
            }

            answer.at(i) += coeff * sum;
	}
    }
   
    if ( this->updateRotationMatrix() ) {
      answer.rotatedWith(* this->rotationMatrix, 'n');
    }
}

void
SUPGElement2 :: computeDiffusionDerivativeTerm_MB(FloatMatrix &answer, MatResponseMode mode, TimeStep *atTime)
{
    int k, undofs = this->computeNumberOfDofs(EID_MomentumBalance);

    answer.resize(undofs, undofs);
    answer.zero();
    FloatMatrix _db, _d, _b, dDB, un_gu;
    double dV, Re = domain->giveEngngModel()->giveUnknownComponent(ReynoldsNumber, VM_Unknown, atTime, domain, NULL);
    GaussPoint *gp;
    FloatArray dDB_u;
    IntegrationRule *iRule = this->integrationRulesArray [ 1 ];
    for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        dV  = this->computeVolumeAround(gp);
        this->computeBMatrix(_b, gp);
        ( ( FluidDynamicMaterial * ) this->giveMaterial() )->giveDeviatoricStiffnessMatrix(_d, mode,gp, atTime);

	this->computeDivTauMatrix(dDB, gp, atTime);
	this->computeUDotGradUMatrix( un_gu, gp, atTime->givePreviousStep() );
	/* standard term */
        _db.beProductOf(_d, _b);
        answer.plusProductUnsym(_b, _db, dV); //answer.plusProduct (_b,_db,area);
        //answer.symmetrized() ;
        answer.times(1. / Re);
    
	/* SUPG term */	

	answer.plusProductUnsym(un_gu, dDB, t_supg * dV * (-1.0) * (1./Re));

    }

    
    if ( this->updateRotationMatrix() ) {
      //answer.rotatedWith(* this->rotationMatrix);
    
      FloatMatrix * Trans = this->rotationMatrix->GiveTransposition();
      answer.rotatedWith(* Trans);
      delete Trans;
    }

}


void
SUPGElement2 :: computePressureTerm_MB(FloatMatrix &answer, TimeStep *atTime)
{
    int k;
    double dV;
    GaussPoint *gp;
    FloatMatrix gu, np, b;
    int undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    int pndofs = this->computeNumberOfDofs(EID_ConservationEquation);

    answer.resize(undofs, pndofs);
    answer.zero();
    IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
    for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        dV  = this->computeVolumeAround(gp);
        this->computeDivUMatrix(gu, gp);
        this->computeNpMatrix(np, gp);

        /* standard term */
        answer.plusProductUnsym(gu, np, ( -1.0 ) * dV);
    }

    iRule = this->integrationRulesArray [ 1 ];
    for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        dV  = this->computeVolumeAround(gp);
        this->computeUDotGradUMatrix( b, gp, atTime->givePreviousStep() );
        this->computeGradPMatrix(np, gp);

        // supg term
        answer.plusProductUnsym(b, np, t_supg * dV);
    }

    if ( this->updateRotationMatrix() ) {
      

      FloatMatrix tmp = answer;
      answer.beProductOf(*rotationMatrix, tmp);
      
    }
}
void
SUPGElement2 :: computeLSICStabilizationTerm_MB(FloatMatrix &answer, TimeStep *atTime)
{
    int k, undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    double dV, rho;
    GaussPoint *gp;
    FloatMatrix b;

    answer.resize(undofs, undofs);
    answer.zero();

    IntegrationRule *iRule = this->integrationRulesArray [ 0 ];
    for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        dV  = this->computeVolumeAround(gp);
        rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, gp, atTime);
        this->computeDivUMatrix(b, gp);

        answer.plusProductSymmUpper(b, b, dV * rho * t_lsic);
    }

    answer.symmetrized();
    
    if ( this->updateRotationMatrix() ) {
      //answer.rotatedWith(* this->rotationMatrix);
      
      FloatMatrix *Trans = this->rotationMatrix->GiveTransposition();
      answer.rotatedWith(* Trans);
      delete Trans;
      
    }
    
}

void
SUPGElement2 :: computeLinearAdvectionTerm_MC(FloatMatrix &answer, TimeStep *atTime)
{
    int k;
    double dV;
    GaussPoint *gp;
    FloatMatrix gu, np;
    int undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    int pndofs = this->computeNumberOfDofs(EID_ConservationEquation);

    answer.resize(pndofs, undofs);
    answer.zero();
    IntegrationRule *iRule = this->integrationRulesArray [ 1 ];
    for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        dV  = this->computeVolumeAround(gp);
        this->computeDivUMatrix(gu, gp);
        this->computeNpMatrix(np, gp);

        /* standard term */
        answer.plusProductUnsym(np, gu, dV);
    }
    if ( this->updateRotationMatrix() ) {
      
      FloatMatrix tmp = answer;
      FloatMatrix *Trans = this->rotationMatrix->GiveTransposition();
      answer.beProductOf(tmp, *Trans);
      delete Trans;
    }
    
}

void
SUPGElement2 :: computeAdvectionTerm_MC(FloatArray &answer, TimeStep *atTime)
{
    // N_epsilon (due to PSPG stabilization)
    FloatMatrix g, b;
    FloatArray u, v;
    GaussPoint *gp;
    double dV, coeff, sum;
    int i, k, pndofs = this->computeNumberOfDofs(EID_ConservationEquation);
    int isd, nsd = this->giveNumberOfSpatialDimensions();

    answer.resize(pndofs);
    answer.zero();

    IntegrationRule *iRule = this->integrationRulesArray [ 1 ];
    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);
    
     if ( this->updateRotationMatrix() ) {
      u.rotatedWith(this->rotationMatrix, 't');
    }

    /* pspg stabilization term */
    for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        this->computeGradPMatrix(g, gp);
        this->computeUDotGradUMatrix( b, gp, atTime);
        v.beProductOf(b, u);
        dV  = this->computeVolumeAround(gp);
        coeff = dV * t_pspg;
        for ( i = 1; i <= pndofs; i++ ) {
            for ( sum = 0.0, isd = 1; isd <= nsd; isd++ ) {
                sum += g.at(isd, i) * v.at(isd);
            }

            answer.at(i) += coeff * sum;
        }
    }
}


void
SUPGElement2 :: computeAdvectionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *atTime)
{
    int pndofs = this->computeNumberOfDofs(EID_ConservationEquation);
    int undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    GaussPoint *gp;
    FloatMatrix g, b;
    double dV, coeff;
    int k;

    answer.resize(pndofs, undofs);
    answer.zero();

    IntegrationRule *iRule = this->integrationRulesArray [ 1 ];
    /* pspg stabilization term */
    for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        this->computeGradPMatrix(g, gp);
        this->computeUDotGradUMatrix( b, gp, atTime);
        dV  = this->computeVolumeAround(gp);
        coeff = dV * t_pspg;
        answer.plusProductUnsym(g, b, coeff);
    }
    if ( this->updateRotationMatrix() ) {
     
      FloatMatrix tmp = answer;
      FloatMatrix *Trans = this->rotationMatrix->GiveTransposition();
      answer.beProductOf(tmp, *Trans);
      delete Trans;
      
      
    }
    
    
    
}

void
SUPGElement2 :: computeDiffusionDerivativeTerm_MC(FloatMatrix &answer, TimeStep *atTime)
{
    int pndofs = this->computeNumberOfDofs(EID_ConservationEquation);
    int undofs = this->computeNumberOfDofs(EID_MomentumBalance);

    answer.resize(pndofs, undofs);
    FloatMatrix dDB, _d, g;
    double dV, coeff, rho;
    GaussPoint *gp;
    IntegrationRule *iRule = this->integrationRulesArray [ 1 ];
    int k;

    for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        dV  = this->computeVolumeAround(gp);
	rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, gp, atTime); 

	  coeff = (-1.0) * dV * t_pspg / rho;
	  //( ( FluidDynamicMaterial * ) this->giveMaterial() )->giveDeviatoricStiffnessMatrix(_d, TangentStiffness,gp, atTime);

	this->computeDivTauMatrix(dDB, gp, atTime);
	this->computeGradPMatrix(g, gp);

	answer.plusProductUnsym(g, dDB, coeff);
    }

    
    answer.zero();
}

void
SUPGElement2 :: computeDiffusionTerm_MC(FloatArray &answer, TimeStep *atTime)
{
    int pndofs = this->computeNumberOfDofs(EID_ConservationEquation);
    answer.resize(pndofs);
   
    /*
    FloatMatrix dDB, _d, g;
    double sum, dV, coeff, rho, isd, nsd = this->giveNumberOfSpatialDimensions();
    GaussPoint *gp;
    IntegrationRule *iRule = this->integrationRulesArray [ 1 ];
    int i,k;
    FloatArray u, dDB_u;

    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);
    
    for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        dV  = this->computeVolumeAround(gp);
	rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, gp, atTime); 

	coeff = (-1.0) * dV * t_pspg / rho;
	
	this->computeGradPMatrix(g, gp);
	this->computeDivTauMatrix(dDB, gp, atTime);
	
	dDB_u.beProductOf(dDB, u);

	for ( i = 1; i <= pndofs; i++ ) {
            for ( sum = 0, isd = 1; isd <= nsd; isd++ ) {
                sum += g.at(isd, i) * dDB_u.at(isd);
            }

            answer.at(i) += coeff * sum;
	}



    }
    */
    
    answer.zero();
}


void
SUPGElement2 :: computeAccelerationTerm_MC(FloatMatrix &answer, TimeStep *atTime)
{
    int pndofs = this->computeNumberOfDofs(EID_ConservationEquation);
    int undofs = this->computeNumberOfDofs(EID_MomentumBalance);
    int k;
    double dV, coeff;
    FloatMatrix g, n;
    GaussPoint *gp;

    answer.resize(pndofs, undofs);
    answer.zero();
    // pspg stabilization term: M_\epsilon term
    IntegrationRule *iRule = this->integrationRulesArray [ 1 ];

    for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        this->computeGradPMatrix(g, gp);
        this->computeNuMatrix(n, gp);
        dV  = this->computeVolumeAround(gp);
        coeff = dV * t_pspg;
        answer.plusProductUnsym(g, n, coeff);
    }

    if ( this->updateRotationMatrix() ) {
      
      FloatMatrix tmp = answer;
      FloatMatrix *Trans = this->rotationMatrix->GiveTransposition();
      answer.beProductOf(tmp, *Trans);
      delete Trans;
    }
    
    




}

void
SUPGElement2 :: computePressureTerm_MC(FloatMatrix &answer, TimeStep *atTime)
{
    int k;
    int pndofs = this->computeNumberOfDofs(EID_ConservationEquation);
    double dV, rho, coeff;
    GaussPoint *gp;
    FloatMatrix g;

    answer.resize(pndofs, pndofs);
    answer.zero();

    IntegrationRule *iRule = this->integrationRulesArray [ 1 ];

    for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        this->computeGradPMatrix(g, gp);
        dV  = this->computeVolumeAround(gp);
        rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, gp, atTime);
        coeff = dV * t_pspg / rho;
        answer.plusProductSymmUpper(g, g, coeff);
    }

    answer.symmetrized();
}


void
SUPGElement2 :: computeBCRhsTerm_MB(FloatArray &answer, TimeStep *atTime)
{
    int undofs = this->computeNumberOfDofs(EID_MomentumBalance);

    answer.resize(undofs);
    answer.zero();

    int i, id, k, n, nLoads;
    double dV, rho;
    Load *load;
    bcGeomType ltype;
    IntegrationRule *iRule = this->integrationRulesArray [ 1 ];
    GaussPoint *gp;
    FloatArray un, gVector, s, helpLoadVector;
    FloatMatrix b, nu;

    // add body load (gravity) termms
    nLoads    = this->giveBodyLoadArray()->giveSize();
    for ( i = 1; i <= nLoads; i++ ) {
        load  = domain->giveLoad( bodyLoadArray.at(i) );
        ltype = load->giveBCGeoType();
        if ( ( ltype == BodyLoadBGT ) && ( load->giveBCValType() == ForceLoadBVT ) ) {
            load->computeComponentArrayAt(gVector, atTime, VM_Total);
            if ( gVector.giveSize() ) {
                for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
                    gp = iRule->getIntegrationPoint(k);
                    this->computeUDotGradUMatrix( b, gp, atTime->givePreviousStep() );
                    this->computeNuMatrix(nu, gp);
                    dV  = this->computeVolumeAround(gp);
                    rho = this->giveMaterial()->giveCharacteristicValue(MRM_Density, gp, atTime);
                    s.beTProductOf(b, gVector);
                    s.times(t_supg * rho * dV);
                    answer.add(s);
                    s.beTProductOf(nu, gVector);
                    s.times(rho * dV);
                    answer.add(s);
                }
            }
        }
    }

    // integrate tractions
    // if no traction bc applied but side marked as with traction load
    // then zero traction is assumed !!!

    // loop over boundary load array
    nLoads    = this->giveBoundaryLoadArray()->giveSize() / 2;
    for ( i = 1; i <= nLoads; i++ ) {
        n     = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
        id    = boundaryLoadArray.at(i * 2);
        load  = dynamic_cast< Load * >( domain->giveLoad(n) );
        ltype = load->giveBCGeoType();
        if ( ltype == EdgeLoadBGT ) {
            this->computeEdgeLoadVector_MB(helpLoadVector, ( Load * ) load, id, atTime);
            if ( helpLoadVector.giveSize() ) {
                answer.add(helpLoadVector);
            }
        } else if ( ltype == SurfaceLoadBGT ) {
            this->computeSurfaceLoadVector_MB(helpLoadVector, ( Load * ) load, id, atTime);
            if ( helpLoadVector.giveSize() ) {
                answer.add(helpLoadVector);
            }
        } else {
            _error("computeForceLoadVector : unsupported load type class");
        }
    }
    if ( this->updateRotationMatrix() ) {
      
      //FloatMatrix *T;
      //int p, i, j, k;
      //double coeff; 
      //T =  this->rotationMatrix;
      
      answer.rotatedWith(this->rotationMatrix, 'n');
      
    } 
}



void
SUPGElement2 :: computeBCRhsTerm_MC(FloatArray &answer, TimeStep *atTime)
{
    int i, k, n, id, nLoads;
    double dV;
    Load *load;
    bcGeomType ltype;
    FloatArray s, gVector, helpLoadVector;
    FloatMatrix g;
    int pndofs = this->computeNumberOfDofs(EID_ConservationEquation);
    IntegrationRule *iRule = this->integrationRulesArray [ 1 ];
    GaussPoint *gp;

    answer.resize(pndofs);
    answer.zero();
    nLoads    = this->giveBodyLoadArray()->giveSize();
    for ( i = 1; i <= nLoads; i++ ) {
        load  = domain->giveLoad( bodyLoadArray.at(i) );
        ltype = load->giveBCGeoType();
        if ( ( ltype == BodyLoadBGT ) && ( load->giveBCValType() == ForceLoadBVT ) ) {
            load->computeComponentArrayAt(gVector, atTime, VM_Total);
            if ( gVector.giveSize() ) {
                for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
                    gp = iRule->getIntegrationPoint(k);
                    this->computeGradPMatrix(g, gp);
                    dV  = this->computeVolumeAround(gp);
                    s.beTProductOf(g, gVector);
                    s.times(t_pspg * dV);
                    answer.add(s);
                }
            }
        }
    }

    // integrate tractions
    // if no traction bc applied but side marked as with traction load
    // then zero traction is assumed !!!

    // loop over boundary load array
    nLoads    = this->giveBoundaryLoadArray()->giveSize() / 2;
    for ( i = 1; i <= nLoads; i++ ) {
        n     = boundaryLoadArray.at(1 + ( i - 1 ) * 2);
        id    = boundaryLoadArray.at(i * 2);
        load  = dynamic_cast< Load * >( domain->giveLoad(n) );
        ltype = load->giveBCGeoType();
        if ( ltype == EdgeLoadBGT ) {
            this->computeEdgeLoadVector_MC(helpLoadVector, ( Load * ) load, id, atTime);
            if ( helpLoadVector.giveSize() ) {
                answer.add(helpLoadVector);
            }
        } else if ( ltype == SurfaceLoadBGT ) {
            this->computeSurfaceLoadVector_MC(helpLoadVector, ( Load * ) load, id, atTime);
            if ( helpLoadVector.giveSize() ) {
                answer.add(helpLoadVector);
            }
        } else {
            _error("computeForceLoadVector : unsupported load type class");
        }
    }
}

void
SUPGElement2 :: updateStabilizationCoeffs(TimeStep *atTime)
{
    //TR1_2D_SUPG :: updateStabilizationCoeffs (atTime);
    /* UGN-Based Stabilization */
    double h_ugn, sum = 0.0, vnorm, t_sugn1, t_sugn2, t_sugn3, u_1, u_2, u_3, z, Re_ugn;
    double dscale, uscale, lscale, tscale, dt;
    //bool zeroFlag = false;
    int i, k, im1;
    FloatArray u, divu;
    FloatMatrix du;

    uscale = domain->giveEngngModel()->giveVariableScale(VST_Velocity);
    lscale = domain->giveEngngModel()->giveVariableScale(VST_Length);
    tscale = domain->giveEngngModel()->giveVariableScale(VST_Time);
    dscale = domain->giveEngngModel()->giveVariableScale(VST_Density);

    this->computeVectorOf(EID_MomentumBalance, VM_Total, atTime, u);
    u.times(uscale);
    double nu;

    // compute averaged viscosity based on rule of mixture
    GaussPoint *gp;

    dt = atTime->giveTimeIncrement() * tscale;

    IntegrationRule *iRule = this->integrationRulesArray [ 1 ];
    gp = iRule->getIntegrationPoint(0);
    nu = this->giveMaterial()->giveCharacteristicValue(MRM_Viscosity, gp, atTime);
    nu *= domain->giveEngngModel()->giveVariableScale(VST_Viscosity);

    for ( k = 0; k < iRule->getNumberOfIntegrationPoints(); k++ ) {
        gp = iRule->getIntegrationPoint(k);
        this->computeDivUMatrix(du, gp);
        divu.beProductOf(du, u);
        sum += divu.at(1);
    }

    sum *= ( 1. / lscale / iRule->getNumberOfIntegrationPoints() );

    /*
     * for (i=1; i<=3;i++) {
     * im1=i-1;
     * sum+= fabs(u.at((im1)*2+1)*b[im1]/lscale + u.at(im1*2+2)*c[im1]/lscale);
     * }
     */
    vnorm = 0.;
    int nsd = this->giveNumberOfSpatialDimensions();
    for ( i = 1; i <= numberOfDofMans; i++ ) {
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

    //this->t_lsic=0.0;
    //this->t_pspg=0.0;
}

void
SUPGElement2 :: computeEdgeLoadVector_MB(FloatArray &answer, Load *load, int id, TimeStep *stepN)
{
    _error("computeEdgeLoadVectorAt_MB: not implemented");
}

void
SUPGElement2 :: computeSurfaceLoadVector_MB(FloatArray &answer, Load *load, int id, TimeStep *stepN)
{
    _error("computeSurfaceLoadVectorAt_MB: not implemented");
}

void
SUPGElement2 :: computeEdgeLoadVector_MC(FloatArray &answer, Load *load, int id, TimeStep *stepN)
{
    _error("computeEdgeLoadVectorAt_MC: not implemented");
}

void
SUPGElement2 :: computeSurfaceLoadVector_MC(FloatArray &answer, Load *load, int id, TimeStep *stepN)
{
    _error("computeEdgeLoadVectorAt_MC: not implemented");
}

void
SUPGElement2 :: computeDeviatoricStress(FloatArray &answer, GaussPoint *gp, TimeStep *tStep)
{
    /* one should call material driver instead */
    FloatArray u, eps;
    FloatMatrix b;
    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, u);
    this->computeBMatrix(b, gp);
    eps.beProductOf(b, u);
    ( ( FluidDynamicMaterial * ) this->giveMaterial() )->computeDeviatoricStressVector(answer, gp, eps, tStep);
}





int
SUPGElement2 :: updateRotationMatrix()
{
    /* returns a tranformation matrix between local coordinate system
     * and global coordinate system, taking into account possible local
     * coordinate system in nodes.
     * if no transformation necessary - returns NULL
     */
    int isT_GtoL, isT_NtoG;
    FloatMatrix T_GtoL, T_NtoG;

    if ( rotationMatrixDefined ) {
        return ( rotationMatrix != NULL );
    }

    rotationMatrixDefined = 1;
    isT_GtoL = this->computeGtoLRotationMatrix(T_GtoL);
    isT_NtoG = this->computeGNDofRotationMatrix(T_NtoG, _toNodalCS);

#ifdef DEBUG
    if ( isT_GtoL ) {
        if ( ( !T_GtoL.isSquare() ) ||
            ( T_GtoL.giveNumberOfRows() != this->computeNumberOfDofs(EID_MomentumBalance) ) ) {
            _error("StructuralElement :: updateRotationMatrix - T_GtoL transformation matrix size mismatch");
        }
    }

    if ( isT_NtoG ) {
        if ( ( !T_NtoG.isSquare() ) ||
            ( T_NtoG.giveNumberOfRows() != this->computeNumberOfDofs(EID_MomentumBalance) ) ) {
            _error("StructuralElement :: updateRotationMatrix - T_NtoG transformation matrix size mismatch");
        }
    }

#endif

    if ( isT_GtoL && T_NtoG.isNotEmpty() ) {
        rotationMatrix = T_GtoL.Times(& T_NtoG);
    } else if ( isT_GtoL ) {
        rotationMatrix = T_GtoL.GiveCopy();
    } else if ( T_NtoG.isNotEmpty() ) {
        rotationMatrix = T_NtoG.GiveCopy();
    } else {
        rotationMatrix = NULL;
    }

    //delete T_GtoL;
    //delete T_GtoNTransp;
    return ( rotationMatrix != NULL );
}



int
SUPGElement2 :: computeGNDofRotationMatrix(FloatMatrix &answer, DofManTransfType mode)
{
    int i, j, k, lastRowPos = 0, lastColPos = 0, flag = 0;

    // test if transformation is necessary
    for ( i = 1; i <= numberOfDofMans; i++ ) {
        flag += this->giveDofManager(i)->requiresTransformation();
    }

    if ( flag == 0 ) {
        answer.beEmptyMtrx();
        return 0;
    }

    // initialize answer
    int gsize = this->computeGlobalNumberOfDofs(EID_MomentumBalance);
    if ( mode == _toGlobalCS ) {
        answer.resize(this->computeNumberOfDofs(EID_MomentumBalance), gsize);
    } else if ( mode == _toNodalCS ) {
        answer.resize( gsize, this->computeNumberOfDofs(EID_MomentumBalance) );
    } else {
        _error("computeGNDofRotationMatrix: unsupported DofManTrasfType value");
    }

    answer.zero();

    FloatMatrix dofManT;
    IntArray dofIDmask;
    int nr, nc;
    // loop over nodes
    for ( i = 1; i <= numberOfDofMans; i++ ) {
        this->giveDofManDofIDMask(i, EID_MomentumBalance, dofIDmask);
        this->giveDofManager(i)->computeDofTransformation(dofManT, & dofIDmask, mode);
        nc = dofManT.giveNumberOfColumns();
        nr = dofManT.giveNumberOfRows();
        for ( j = 1; j <= nr; j++ ) {
            for ( k = 1; k <= nc; k++ ) {
                // localize node contributions
                answer.at(lastRowPos + j, lastColPos + k) = dofManT.at(j, k);
            }
        }

        lastRowPos += nr;
        lastColPos += nc;
    }

    return 1;
}


int
SUPGElement2 :: computeGNLoadRotationMatrix(FloatMatrix &answer, DofManTransfType mode)
{
    int i, j, k, lastRowPos = 0, lastColPos = 0, flag = 0;

    // test if transformation is necessary
    for ( i = 1; i <= numberOfDofMans; i++ ) {
        flag += this->giveDofManager(i)->requiresTransformation();
    }

    if ( flag == 0 ) {
        answer.beEmptyMtrx();
        return 0;
    }

    // initialize answer
    int gsize = this->computeGlobalNumberOfDofs(EID_MomentumBalance);
    if ( mode == _toGlobalCS ) {
        answer.resize(this->computeNumberOfDofs(EID_MomentumBalance), gsize);
    } else if ( mode == _toNodalCS ) {
        answer.resize( gsize, this->computeNumberOfDofs(EID_MomentumBalance) );
    } else {
        _error("computeGNDofRotationMatrix: unsupported DofManTrasfType value");
    }

    answer.zero();

    FloatMatrix dofManT;
    IntArray dofIDmask;
    int nr, nc;
    // loop over nodes
    for ( i = 1; i <= numberOfDofMans; i++ ) {
        this->giveDofManDofIDMask(i, EID_MomentumBalance, dofIDmask);
        this->giveDofManager(i)->computeLoadTransformation(dofManT, & dofIDmask, mode);
        nc = dofManT.giveNumberOfColumns();
        nr = dofManT.giveNumberOfRows();
        for ( j = 1; j <= nr; j++ ) {
            for ( k = 1; k <= nc; k++ ) {
                // localize node contributions
                answer.at(lastRowPos + j, lastColPos + k) = dofManT.at(j, k);
            }
        }

        lastRowPos += nr;
        lastColPos += nc;
    }

    return 1;
}











} // end namespace oofem
