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

#include "shell7basexfem.h"
#include "shell7base.h"
#include "enrichmentitem.h"
#include "xfemmanager.h"
#include "constantpressureload.h"

namespace oofem {
IntArray Shell7BaseXFEM :: dofId_Midplane(3);
IntArray Shell7BaseXFEM :: dofId_Director(3);
IntArray Shell7BaseXFEM :: dofId_InhomStrain(1); 
bool Shell7BaseXFEM :: __initializedFieldDofId = Shell7BaseXFEM :: initDofId();

Shell7BaseXFEM :: Shell7BaseXFEM(int n, Domain *aDomain) : Shell7Base(n, aDomain), XfemElementInterface(this) 
{

}

int 
Shell7BaseXFEM :: checkConsistency()
{
    Shell7Base :: checkConsistency();
    this->xMan =  this->giveDomain()->giveXfemManager(1);
    return 1;
}


IRResultType Shell7BaseXFEM :: initializeFrom(InputRecord *ir)
{
    Shell7Base :: initializeFrom(ir);
    return IRRT_OK; 
}


Interface
*Shell7BaseXFEM :: giveInterface(InterfaceType it)
{
    if ( it != XfemElementInterfaceType ) {
        return Shell7Base :: giveInterface(it);
    } else if ( it == XfemElementInterfaceType ) {
        return ( XfemElementInterface * ) this;
    } else {
        return Shell7Base :: giveInterface(it);
    }
}


double 
Shell7BaseXFEM :: giveGlobalZcoord(GaussPoint *gp) 
{

    this->setupDelaminationXiCoordList();
    this->setupGPDelaminationGroupList();

    double xiRef = gp->giveCoordinate(3);
    int dGroup   = this->giveDelaminationGroupAt( xiRef );
    double xiMid = this->giveDelaminationGroupMidXi(dGroup);
    
    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * > (this->element->giveCrossSection());
    return (xiRef - xiMid)*layeredCS->computeIntegralThick()*0.5; // new xi-coord measured from dGroup c.s. 
    
}


void
Shell7BaseXFEM :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    // Returns the total id mask of the dof manager - regular id's + enriched id's

    // Continuous part
    Shell7Base ::giveDofManDofIDMask(inode, ut, answer);

    // Discontinuous part
    DofManager *dMan = this->giveDofManager(inode);
    for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { // Only one is supported at the moment
        EnrichmentItem *ei = this->xMan->giveEnrichmentItem(i);
        for ( int j = 1; j <= ei->giveNumberOfEnrichmentDomains(); j++ ) {
            if ( ei->isDofManEnrichedByEnrichmentDomain(dMan,j) ) {
                IntArray eiDofIdArray;
                ei->giveEIDofIdArray(eiDofIdArray, j);
                answer.followedBy(eiDofIdArray);
            }
        }
    }
}


void
Shell7BaseXFEM :: discGiveDofManDofIDMask(int inode,  int enrichmentdomainNumber, IntArray &answer) const
{
    // Returns the id mask corresponding to a given enrichment domain 
    ///@todo: input should also be *ei and then check *ed
    if ( enrichmentdomainNumber == 0 ) { // return mask corresponding to the regular id's
        Shell7Base ::giveDofManDofIDMask(inode, EID_MomentumBalance, answer);
    } else {
        // Discontinuous part
        DofManager *dMan = this->giveDofManager(inode);
        for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { // Only one is supported at the moment
            EnrichmentItem *ei = this->xMan->giveEnrichmentItem(i);
            if ( ei->isDofManEnrichedByEnrichmentDomain(dMan,enrichmentdomainNumber) ) {
                
                    ei->giveEIDofIdArray(answer, enrichmentdomainNumber);
            }
        }
    }
}


void
Shell7BaseXFEM :: evalCovarBaseVectorsAt(GaussPoint *gp, FloatMatrix &gcov, FloatArray &genEpsC)
{
    // Continuous part
    Shell7Base :: evalCovarBaseVectorsAt(gp, gcov, genEpsC);

    // Discontinuous part - ///@todo bad implementation regarding enr. functions - should be changed
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { // Only one is supported at the moment
        Delamination *dei =  dynamic_cast< Delamination * >( this->xMan->giveEnrichmentItem(i) ); // should check success
        //EnrichmentFunction *ef = dei->giveEnrichmentFunction(1);

        for ( int j = 1; j <= dei->giveNumberOfEnrichmentDomains(); j++ ) {
            if ( dei->isElementEnrichedByEnrichmentDomain(this, j) ) {
                double xi0 = dei->enrichmentDomainXiCoords.at(j);
                double H   = dei->heaviside(gp->giveCoordinate(3), xi0);        
                if ( H > 0.1 ) {
                    FloatArray dGenEps;
                    computeDiscGeneralizedStrainVector(dGenEps, gp, dei, j, tStep);
                    FloatMatrix gcovd; 
                    Shell7Base :: evalCovarBaseVectorsAt(gp, gcovd, dGenEps);
                    gcov.add(H,gcovd); 
                }

            }
        }
    }
}


void
Shell7BaseXFEM :: computeDiscGeneralizedStrainVector(FloatArray &answer, GaussPoint *gp, EnrichmentItem *ei, int enrichmentDomainNumber, TimeStep *tStep)
{
    FloatArray solVecD;
    IntArray eiDofIdArray;
    ei->giveEIDofIdArray(eiDofIdArray, enrichmentDomainNumber);
    this->giveSolutionVector(solVecD, eiDofIdArray, tStep);  
    FloatMatrix B11, B22, B32, B43, B53;
    Shell7Base :: computeBmatricesAt(gp, B11, B22, B32, B43, B53);
    Shell7Base :: computeGeneralizedStrainVector(answer, solVecD, B11, B22, B32, B43, B53);
}


int
Shell7BaseXFEM :: giveNumberOfDofs()
{
    // Continuous part
    int nDofs = Shell7Base ::giveNumberOfDofs();

    // Discontinuous part
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        DofManager *dMan = this->giveDofManager(i);
        for ( int j = 1; j <= this->xMan->giveNumberOfEnrichmentItems(); j++ ) {
            EnrichmentItem *ei = this->xMan->giveEnrichmentItem(j);
            for ( int k = 1; k <= ei->giveNumberOfEnrichmentDomains(); k++ ) {
                if ( ei->isDofManEnrichedByEnrichmentDomain(dMan,k) ) {
                    IntArray eiDofIdArray;
                    ei->giveEIDofIdArray(eiDofIdArray, k);
                    nDofs += eiDofIdArray.giveSize();
                }
            }
        }
    }
    return nDofs;
}


void 
Shell7BaseXFEM :: computeOrderingArray( IntArray &orderingArray, IntArray &activeDofsArray,  int enrichmentDomainNumber, SolutionField field)
{
    // Routine to extract vector given an array of dofid items
    // If a certain dofId does not exist a zero is used as value

    const IntArray &ordering_cont = this->giveOrdering(field);
    const IntArray &fieldDofId    = this->giveFieldDofId(field);

    IntArray ordering_temp, activeDofsArrayTemp;
    ordering_temp.resize(ordering_cont.giveSize());
    activeDofsArrayTemp.resize(ordering_cont.giveSize());


    int activeDofPos = 0, activeDofIndex = 0, orderingDofIndex = 0;
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        DofManager *dMan = this->giveDofManager(i);
        IntArray dofManDofIdMask, dofManDofIdMaskAll; 

        if ( enrichmentDomainNumber == 0 ) { // return mask corresponding to the regular id's
           Shell7Base ::giveDofManDofIDMask(i, EID_MomentumBalance, dofManDofIdMask);
        } else {
            EnrichmentItem *ei = this->xMan->giveEnrichmentItem(1); ///@todo: *ei should be input
            if ( ei->isDofManEnrichedByEnrichmentDomain(dMan,enrichmentDomainNumber) ) {
                ei->giveEIDofIdArray(dofManDofIdMask, enrichmentDomainNumber);
            }
        }
    

        for (int j = 1; j <= dofManDofIdMask.giveSize(); j++ ) {
            int pos = dMan->findDofWithDofId( (DofIDItem) dofManDofIdMask.at(j) );
            activeDofPos++;
            ordering_temp      .at(activeDofPos) = orderingDofIndex + pos;
            activeDofsArrayTemp.at(activeDofPos) = activeDofIndex   + j;
        }
        this->giveDofManDofIDMask(i, EID_MomentumBalance, dofManDofIdMaskAll);
        orderingDofIndex += dofManDofIdMaskAll.giveSize();
        activeDofIndex   += fieldDofId.giveSize();

    }
    
    // Reduce arrays to actual size 
    ///@todo will not work if there are several ei
    int numActiveDofs = activeDofPos;
    orderingArray.resize(numActiveDofs), activeDofsArray.resize(numActiveDofs);
    
    for ( int i = 1; i <= numActiveDofs; i++ ) {
        orderingArray.at(i) = ordering_temp.at(i); 
        activeDofsArray.at(i) = activeDofsArrayTemp.at(i);
    }
}

#if 0
void 
Shell7BaseXFEM :: edgeComputeOrderingArray( IntArray &orderingArray, IntArray &activeDofsArray, int iEdge, int enrichmentDomainNumber, SolutionField field)
{
    // Routine to extract vector given an array of dofid items
    // If a certain dofId does not exist a zero is used as value

    //const IntArray &ordering_cont = this->giveOrdering(field);
    
    IntArray ordering_cont; 
    this->giveEdgeDofMapping(ordering_cont, iEdge);
    ordering_cont.printYourself();
    const IntArray &fieldDofId    = this->giveFieldDofId(field);

    IntArray ordering_temp, activeDofsArrayTemp;
    ordering_temp.resize(ordering_cont.giveSize());
    activeDofsArrayTemp.resize(ordering_cont.giveSize());

    //new
    IntArray edgeNodes;
    this->fei->computeLocalEdgeMapping(edgeNodes, iEdge);

    int activeDofPos = 0, activeDofIndex = 0, orderingDofIndex = 0;
    //for ( int i = 1; i <= numberOfDofMans; i++ ) {
    for ( int i = 1; i <= edgeNodes.giveSize(); i++ ) {
        DofManager *dMan = this->giveDofManager(edgeNodes.at(i));
        IntArray dofManDofIdMask, dofManDofIdMaskAll; 

        if ( enrichmentDomainNumber == 0 ) { // return mask corresponding to the regular id's
           Shell7Base ::giveDofManDofIDMask(i, EID_MomentumBalance, dofManDofIdMask);
        } else {
            EnrichmentItem *ei = this->xMan->giveEnrichmentItem(1); ///@todo: *ei should be input
            if ( ei->isDofManEnrichedByEnrichmentDomain(dMan,enrichmentDomainNumber) ) {
                ei->giveEIDofIdArray(dofManDofIdMask, enrichmentDomainNumber);
            }
        }
    

        for (int j = 1; j <= dofManDofIdMask.giveSize(); j++ ) {
            int pos = dMan->findDofWithDofId( (DofIDItem) dofManDofIdMask.at(j) );
            activeDofPos++;
            ordering_temp      .at(activeDofPos) = orderingDofIndex + pos;
            activeDofsArrayTemp.at(activeDofPos) = activeDofIndex   + j;
        }
        this->giveDofManDofIDMask(i, EID_MomentumBalance, dofManDofIdMaskAll);
        orderingDofIndex += dofManDofIdMaskAll.giveSize();
        activeDofIndex   += fieldDofId.giveSize();

    }
    
    // Reduce arrays to actual size 
    ///@todo will not work if there are several ei
    int numActiveDofs = activeDofPos;
    orderingArray.resize(numActiveDofs), activeDofsArray.resize(numActiveDofs);
    
    for ( int i = 1; i <= numActiveDofs; i++ ) {
        orderingArray.at(i) = ordering_temp.at(i); 
        activeDofsArray.at(i) = activeDofsArrayTemp.at(i);
    }
}
#endif

void
Shell7BaseXFEM :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
//
// Computes internal forces as the array of: [f_c, f_d1, ..., f_dn]
{
    answer.resize( this->giveNumberOfDofs() );
    answer.zero();
    FloatArray solVec, temp;

    // Continuous part
    this->giveUpdatedSolutionVector(solVec, tStep);
    this->computeSectionalForces(temp, tStep, solVec, useUpdatedGpRecord);
    IntArray activeDofs, ordering; 
    this->computeOrderingArray(ordering, activeDofs, 0, All);
    answer.assemble(temp, ordering);

    // Disccontinuous part
    for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { // Only one is supported at the moment
        Delamination *dei =  dynamic_cast< Delamination * >( this->xMan->giveEnrichmentItem(i) ); // should check success

        for ( int j = 1; j <= dei->giveNumberOfEnrichmentDomains(); j++ ) {
            if ( dei->isElementEnrichedByEnrichmentDomain(this, j) ) {
                IntArray eiDofIdArray;
                dei->giveEIDofIdArray(eiDofIdArray, j);
                FloatArray solVecD;
                this->giveSolutionVector(solVecD, eiDofIdArray, tStep);  
                this->discComputeSectionalForces(temp, tStep, solVec, solVecD, useUpdatedGpRecord, dei, j);                      

                // Assemble
                this->computeOrderingArray(ordering, activeDofs,  j, All);
                FloatArray tempRed;
                tempRed.beSubArrayOf(temp, activeDofs);
                answer.assemble(tempRed, ordering);
            }
        }
    }
    
}

void
Shell7BaseXFEM :: discComputeSectionalForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, FloatArray &solVecD, int useUpdatedGpRecord, 
                     EnrichmentItem *ei, int enrichmentDomainNumber)
//
{
    Delamination *dei =  dynamic_cast< Delamination * >( ei );
    double xi0 = 0.0;
    if( enrichmentDomainNumber == 0 ) {
        xi0 = -1.0e6;
    } else {
        xi0 = dei->enrichmentDomainXiCoords.at(enrichmentDomainNumber);
    }

    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     // conversion of types
    FloatArray f1(18), f2(18), f3(6);
    f1.zero(); f2.zero(); f3.zero();

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRuleL = integrationRulesArray [ layer - 1 ];
        Material *mat = domain->giveMaterial( this->layeredCS->giveLayerMaterial(layer) );

        for ( int j = 1; j <= iRuleL->getNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRuleL->getIntegrationPoint(j - 1);


            if ( gp->giveCoordinate(3) > xi0 ) { // Should be enriched ///@todo not general!

                FloatMatrix B11, B22, B32, B43, B53;
                this->computeBmatricesAt(gp, B11, B22, B32, B43, B53);
                FloatArray genEps, genEpsD;
                this->computeGeneralizedStrainVector(genEps,  solVec,  B11, B22, B32, B43, B53);
                this->computeGeneralizedStrainVector(genEpsD, solVecD, B11, B22, B32, B43, B53);

                double zeta = giveGlobalZcoord(gp);
                FloatArray N, M, T, Ms;
                double Ts = 0.;
                this->computeSectionalForcesAt(N, M, T, Ms, Ts, gp, mat, tStep, genEps, genEpsD, zeta);

                // Computation of nodal forces: f = B^t*[N M T Ms Ts]^t
                FloatArray f1temp(18), f2temp(18), f3temp(6), temp;
                // f1 = BT11*N
                f1temp.beTProductOf(B11, N);

                // f2 = BT22*M + BT32*T
                f2temp.beTProductOf(B22, M);
                temp.beTProductOf(B32, T);
                f2temp.add(temp);

                // f3 = BT43*Ms + BT53*Ts
                f3temp.beTProductOf(B43, Ms);
                for ( int i = 1; i <= 6; i++ ) {
                    f3temp.at(i) += B53.at(1, i) * Ts;
                }

                double dV = this->computeVolumeAroundLayer(gp, layer);
                f1.add(dV, f1temp);
                f2.add(dV, f2temp);
                f3.add(dV, f3temp);
            
            }

        }
    }

    const IntArray &ordering_phibar = this->giveOrdering(Midplane);
    const IntArray &ordering_m      = this->giveOrdering(Director);
    const IntArray &ordering_gam    = this->giveOrdering(InhomStrain);
    
    answer.resize( Shell7Base :: giveNumberOfDofs() );
    answer.zero();
    answer.assemble(f1, ordering_phibar);
    answer.assemble(f2, ordering_m);
    answer.assemble(f3, ordering_gam);

}


void
Shell7BaseXFEM :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{

    int ndofs = this->giveNumberOfDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();
    
    FloatMatrix temp;
    FloatArray solVec;
    
    // Continuous part
    this->giveUpdatedSolutionVector(solVec, tStep);
    this->new_computeBulkTangentMatrix(temp, solVec, solVec, solVec, rMode, tStep);
    IntArray ordering, activeDofs;
    this->computeOrderingArray(ordering, activeDofs, 0, All);
    answer.assemble(temp, ordering, ordering);


    // Disccontinuous part
    for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { // Only one is supported at the moment
        Delamination *dei =  dynamic_cast< Delamination * >( this->xMan->giveEnrichmentItem(i) ); // should check success

    #if 1
        for ( int j = 1; j <= dei->giveNumberOfEnrichmentDomains(); j++ ) {
            if ( dei->isElementEnrichedByEnrichmentDomain(this, j) ) {
                IntArray eiDofIdArray;
                FloatArray solVecJ;
                dei->giveEIDofIdArray(eiDofIdArray, j);
                this->giveSolutionVector(solVecJ, eiDofIdArray, tStep);  

                for ( int k = 1; k <= dei->giveNumberOfEnrichmentDomains(); k++ ) {
                    if ( dei->isElementEnrichedByEnrichmentDomain(this, k) ) {
                        FloatArray solVecK;
                        dei->giveEIDofIdArray(eiDofIdArray, k);
                        this->giveSolutionVector(solVecK, eiDofIdArray, tStep);  
                        discComputeBulkTangentMatrix(temp, solVec, solVecJ, solVecK, rMode, tStep, dei, j, k);


                        // Assemble part correpsonding to active dofs
                        IntArray orderingJ, orderingK, activeDofsJ, activeDofsK;
                        computeOrderingArray(orderingJ, activeDofsJ, j, All);
                        computeOrderingArray(orderingK, activeDofsK, k, All);
                        FloatMatrix tempRed;
                        tempRed.beSubMatrixOf(temp, activeDofsJ, activeDofsK);
                        answer.assemble(tempRed, orderingJ, orderingK);
                        //orderingJ.printYourself();
                        //orderingK.printYourself();
                        //printf("Kdd \n");
                        //tempRed.printYourself();
                    }
                }
                

            }
        }


        // First column

        for ( int j = 1; j <= dei->giveNumberOfEnrichmentDomains(); j++ ) {
            if ( dei->isElementEnrichedByEnrichmentDomain(this, j) ) {
                IntArray eiDofIdArray;
                FloatArray solVecJ;
                dei->giveEIDofIdArray(eiDofIdArray, j);
                this->giveSolutionVector(solVecJ, eiDofIdArray, tStep);  
                discComputeBulkTangentMatrix(temp, solVec, solVecJ, solVec, rMode, tStep, dei, j, 0);


                // Assemble part correpsonding to active dofs
                IntArray orderingJ, activeDofsJ;
                computeOrderingArray(orderingJ, activeDofsJ, j, All);

                // First column K_{di,c}
                FloatMatrix tempRed;
                tempRed.beSubMatrixOf(temp, activeDofsJ, activeDofs);
                answer.assemble(tempRed, orderingJ, ordering);

                // First row  K_{c,di} = K_{di,c}^T
                FloatMatrix tempRedT;
                tempRedT.beTranspositionOf(tempRed);
                answer.assemble(tempRedT, ordering, orderingJ);
                

                //FloatMatrix test;
                //test.beSubMatrixOf(tempRed,1,7,15,21);
                //printf("Kdc \n");
                //test.printYourself();
                
                
            }
            
        }

         // First row
        #if 0
        for ( int k = 1; k <= dei->giveNumberOfEnrichmentDomains(); k++ ) {
            if ( dei->isElementEnrichedByEnrichmentDomain(this, k) ) {
                IntArray eiDofIdArray;
                FloatArray solVecK;
                dei->giveEIDofIdArray(eiDofIdArray, k);
                this->giveSolutionVector(solVecK, eiDofIdArray, tStep);  
                discComputeBulkTangentMatrix(temp, solVec, solVec, solVecK, rMode, tStep, dei, 0, k);


                // Assemble part correpsonding to active dofs
                IntArray orderingK, activeDofsK;
                computeOrderingArray(orderingK, activeDofsK, k, All);
                
                FloatMatrix tempRed;
                tempRed.beSubMatrixOf(temp, activeDofs, activeDofsK);
                answer.assemble(tempRed, ordering, orderingK);
            }
        }
        #endif

    #endif

    }
    //answer.symmetrized();
}


void
Shell7BaseXFEM :: discComputeBulkTangentMatrix(FloatMatrix &answer, FloatArray &solVec,  FloatArray &solVecI, FloatArray &solVecJ, MatResponseMode rMode, TimeStep *tStep,
                    EnrichmentItem *ei, int enrichmentDomainNumberI, int enrichmentDomainNumberJ)
{
    FloatMatrix A [ 3 ] [ 3 ], lambdaI [ 3 ], lambdaJ [ 3 ];
    FloatMatrix L(18,18);
    FloatMatrix B11, B22, B32, B43, B53, B;
    FloatArray S1g(3), S2g(3), S3g(3);
    FloatMatrix K(42,42), tempAnswer(42,42);
    K.zero(); tempAnswer.zero();
    double xi0I = 0.0; 
    double xi0J = 0.0;
    Delamination *dei =  dynamic_cast< Delamination * >( ei );
    if( enrichmentDomainNumberI == 0 ) {
        xi0I = -1.0e6;
    } else {
        xi0I = dei->enrichmentDomainXiCoords.at(enrichmentDomainNumberI);
    }
    if( enrichmentDomainNumberJ == 0 ) {
        xi0J = -1.0e6;
    } else {
        xi0J = dei->enrichmentDomainXiCoords.at(enrichmentDomainNumberJ);
    }


    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRule = integrationRulesArray [ layer - 1 ];
        Material *mat = domain->giveMaterial( this->layeredCS->giveLayerMaterial(layer) );

        for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(i);
            double zeta = giveGlobalZcoord(gp);

            if ( gp->giveCoordinate(3) > xi0I  &&   gp->giveCoordinate(3) > xi0J  ) { // Should be enriched ///@todo not general!

                this->computeBmatricesAt(gp, B11, B22, B32, B43, B53);
                FloatArray genEpsI, genEpsJ, genEps;
                this->computeGeneralizedStrainVector(genEpsI, solVecI, B11, B22, B32, B43, B53);
                this->computeGeneralizedStrainVector(genEpsJ, solVecJ, B11, B22, B32, B43, B53);
                this->computeGeneralizedStrainVector(genEps , solVec , B11, B22, B32, B43, B53);

                // Material stiffness
                Shell7Base :: computeLinearizedStiffness(gp, mat, tStep, S1g, S2g, S3g, A, genEps);


                this->computeLambdaMatrices(lambdaI, genEpsI, zeta);
                this->computeLambdaMatrices(lambdaJ, genEpsJ, zeta);

                this->computeBmatrixAt(gp, B, 0, 0);
                // L = sum_{i,j} (lambdaI_i)^T * A^ij * lambdaJ_j
                // Naive implementation - should be optimized 
                // note: L will only be symmetric if lambdaI = lambdaJ
                FloatMatrix temp;
                L.zero();
                for ( int j = 0; j < 3; j++ ) {
                    for ( int k = 0; k < 3; k++ ) {
                        this->computeTripleProduct(temp, lambdaI [ j ], A [ j ][ k ], lambdaJ [ k ]);
                        L.add(temp);
                    }
                }
     
                FloatMatrix Ktemp, K;
                Ktemp.beProductOf(L, B);
                double dV = this->computeVolumeAroundLayer(gp, layer);
                K.beTProductOf(B,Ktemp);
                tempAnswer.add(dV, K);
            }
        }
    }


    int ndofs = Shell7Base :: giveNumberOfDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();
    const IntArray &ordering = this->giveOrdering(All);
    answer.assemble(tempAnswer, ordering, ordering);

}


void
Shell7BaseXFEM :: computeMassMatrixNum(FloatMatrix &answer, TimeStep *tStep)
{
    // Num refers in this case to  numerical integration in both in-plane and through the thickness.
    // For analytically integrated throught he thickness, see computeMassMatrix


    FloatMatrix mass, temp;
    FloatArray solVec;
    this->giveUpdatedSolutionVector(solVec, tStep);
    int ndofs = this->giveNumberOfDofs();
    temp.resize(ndofs, ndofs);
    temp.zero();


    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >( this->giveCrossSection() );
    int numberOfLayers = layeredCS->giveNumberOfLayers();     // conversion of data

    FloatMatrix M11(18, 18), M12(18, 18), M13(18, 6), M22(18, 18), M23(18, 6), M33(6, 6);
    M11.zero();
    M12.zero();
    M13.zero();
    M22.zero();
    M23.zero();
    M33.zero();

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRuleL = integrationRulesArray [ layer - 1 ];
        Material *mat = domain->giveMaterial( layeredCS->giveLayerMaterial(layer) );

        for ( int j = 1; j <= iRuleL->getNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRuleL->getIntegrationPoint(j - 1);

            FloatMatrix N11, N22, N33;
            this->computeNmatricesAt(gp, N11, N22, N33);
            FloatArray xbar, m;
            double gam = 0.;
            //this->computeSolutionFields(xbar, m, gam, solVec, N11, N22, N33);
            FloatArray localCoords = * gp->giveCoordinates();
            this->giveUnknownsAt(localCoords, solVec, xbar, m, gam, tStep);
            //this->computeNmatrixAt(gp, N);
            //unknowns.beProductOf(N,a); // [xbar, m, gam]^T
            //m.setValues(3, unknowns.at(4), unknowns.at(5), unknowns.at(6) );
            //double gam = unknowns.at(7);


            /* Consistent mass matrix M = int{N^t*mass*N}
             *
             *         3    3    1
             *         3 [a*I  b*I   c*m      [A  B  C
             * mass =   3       d*I   e*m    =     D  E
             *         1  sym       f*m.m]     sym   F]
             */


            double zeta = giveGlobalZcoord(gp);
            double fac1 = 4;
            double fac2 = 2.0 * zeta * ( 2.0 + gam * zeta );
            double fac3 = 2.0 * zeta * zeta;
            double fac4 = zeta * zeta * ( 2.0 + gam * zeta ) * ( 2.0 + gam * zeta );
            double fac5 = zeta * zeta * zeta * ( 2.0 + gam * zeta );
            double fac6 = zeta * zeta * zeta * zeta;
            FloatMatrix mass11(3, 3), mass12(3, 3), mass13(3, 1), mass21(3, 3), mass22(3, 3), mass23(3, 1), mass31(1, 3), mass32(1, 3), mass33(1, 1);
            mass11.zero();
            mass12.zero();
            mass13.zero();
            mass21.zero();
            mass22.zero();
            mass23.zero();
            mass31.zero();
            mass32.zero();
            mass33.zero();
            mass.resize(7, 7);
            mass11.at(1, 1) = mass11.at(2, 2) = mass11.at(3, 3) = fac1;          // A
            mass12.at(1, 1) = mass12.at(2, 2) = mass12.at(3, 3) = fac2;          // B
            mass13.at(1, 1) = fac3 * m.at(1);
            mass13.at(2, 1) = fac3 * m.at(2);
            mass13.at(3, 1) = fac3 * m.at(3);            // C
            mass22.at(1, 1) = mass22.at(2, 2) = mass22.at(3, 3) = fac4;          // D
            mass23.at(1, 1) = fac5 * m.at(1);
            mass23.at(2, 1) = fac5 * m.at(2);
            mass23.at(3, 1) = fac5 * m.at(3);            // E
            mass33.at(1, 1) = fac6 * m.dotProduct(m);            // F
            mass21.beTranspositionOf(mass12);
            mass31.beTranspositionOf(mass13);
            mass32.beTranspositionOf(mass23);
            //mass.symmetrized();

            double dV = this->computeVolumeAroundLayer(gp, layer);
            double rho = mat->give('d', gp);

            FloatMatrix M11temp, M12temp, M13temp, M22temp, M23temp, M33temp;
            this->computeTripleProduct(M11temp, N11, mass11, N11);
            this->computeTripleProduct(M12temp, N11, mass12, N22);
            this->computeTripleProduct(M13temp, N11, mass13, N33);
            this->computeTripleProduct(M22temp, N22, mass22, N22);
            this->computeTripleProduct(M23temp, N22, mass23, N33);
            this->computeTripleProduct(M33temp, N33, mass33, N33);
            M11.add(0.25 * rho * dV, M11temp);
            M12.add(0.25 * rho * dV, M12temp);
            M13.add(0.25 * rho * dV, M13temp);
            M22.add(0.25 * rho * dV, M22temp);
            M23.add(0.25 * rho * dV, M23temp);
            M33.add(0.25 * rho * dV, M33temp);
        }

        
        answer.resize(ndofs, ndofs);
        answer.zero();

        const IntArray &ordering_phibar = this->giveOrdering(Midplane);
        const IntArray &ordering_m = this->giveOrdering(Director);
        const IntArray &ordering_gam = this->giveOrdering(InhomStrain);
        answer.assemble(M11, ordering_phibar, ordering_phibar);
        answer.assemble(M12, ordering_phibar, ordering_m);
        answer.assemble(M13, ordering_phibar, ordering_gam);
        answer.assemble(M22, ordering_m,      ordering_m);
        answer.assemble(M23, ordering_m,      ordering_gam);
        answer.assemble(M33, ordering_gam,    ordering_gam);

        FloatMatrix M21, M31, M32;
        M21.beTranspositionOf(M12);
        M31.beTranspositionOf(M13);
        M32.beTranspositionOf(M23);
        answer.assemble(M21, ordering_m,      ordering_phibar);
        answer.assemble(M31, ordering_gam,    ordering_phibar);
        answer.assemble(M32, ordering_gam,    ordering_m);
        answer.symmetrized();

    }
}


IntArray
Shell7BaseXFEM :: giveFieldDofId(SolutionField fieldType) const
{
    if ( fieldType == Midplane ) {
        return this->dofId_Midplane;
    } else if ( fieldType == Director  )   {
        return this->dofId_Director;
    } else if ( fieldType == InhomStrain  )   {
        return this->dofId_InhomStrain;
    } else if ( fieldType == All  )   {
        IntArray dofId;
        Shell7Base::giveDofManDofIDMask(0, EID_MomentumBalance, dofId);
        return dofId;
    } else {
        _error("giveFieldDofId: unknown fieldType");
        return 0;
    }
}


void
Shell7BaseXFEM :: computeEdgeLoadVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *tStep, ValueModeType mode)
{
    BoundaryLoad *edgeLoad = dynamic_cast< BoundaryLoad * >( load );
    if ( edgeLoad ) {
        answer.resize( this->computeNumberOfDofs(EID_MomentumBalance) );
        answer.zero();

        // Continuous part
        FloatArray fT;
        this->computeTractionForce(fT, iEdge, edgeLoad, tStep);

        IntArray activeDofs, ordering; 
        this->computeOrderingArray(ordering, activeDofs, 0, All);
        answer.assemble(fT, ordering);        

        
        //@todo add check to see if size components = 8
        FloatArray componentsTemp, coordsTemp(1);
        coordsTemp.at(1) = 0.0;
        edgeLoad->computeValueAt(componentsTemp, tStep, coordsTemp, VM_Total);
        double xi = componentsTemp.at(8); // use the 8th component to store the-xi coord where the load acts

        // Disccontinuous part
        FloatArray temp;
        for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { // Only one is supported at the moment
            Delamination *dei =  dynamic_cast< Delamination * >( this->xMan->giveEnrichmentItem(i) ); // should check success

            for ( int j = 1; j <= dei->giveNumberOfEnrichmentDomains(); j++ ) {
                if ( dei->isElementEnrichedByEnrichmentDomain(this, j) ) {
                    double xi0 = dei->enrichmentDomainXiCoords.at(j);
                    if ( xi > xi0 ) {
                        this->computeTractionForce(temp, iEdge, edgeLoad, tStep);
                        // Assemble
                        this->computeOrderingArray(ordering, activeDofs,  j, All);
                        FloatArray tempRed;
                        tempRed.beSubArrayOf(temp, activeDofs);
                        answer.assemble(tempRed, ordering);
                    }
                }
            }
        }

        return;
    } else {
        _error("Shell7Base :: computeEdgeLoadVectorAt: load type not supported");
        return;
    }
}








// Delamination specific

#if 1

void
Shell7BaseXFEM :: setupDelaminationXiCoordList()
{
    if ( this->delaminationXiCoordList.size()==0 ) {
        // Stores a paired list with the EnrichmentDomain# and the corresponding xi-coord of the delamination.
        int numEI = xMan->giveNumberOfEnrichmentItems();
        for ( int i = 1; i <= numEI; i++ ) {
            Delamination *dei =  dynamic_cast< Delamination * >( xMan->giveEnrichmentItem(i) );
            if ( dei ) {
                int numED = dei->giveNumberOfEnrichmentDomains(); // numEnrDomains max possible number
                int pos = 1;
                for ( int j = 1; j <= numED; j++ ) {
                    if( dei->isElementEnrichedByEnrichmentDomain(this->element, j) ) { 
                        std::pair<int, double> pid;
                        pid.first  = pos;
                        pid.second = dei->enrichmentDomainXiCoords.at(j); 
                        this->delaminationXiCoordList.push_back(pid); 
                        pos++;
                    }
                }
            } 
        }

        // Sort xi-coords in acending order.
        this->delaminationXiCoordList.sort(sortFunc);

    }
}

void 
Shell7BaseXFEM :: setupGPDelaminationGroupList() 
{
    // Creates a list wich stores the gp# and the dGroup# for quick access later.
    if ( this->gpDelaminationGroupList.size()==0 ) {
        int numberOfLayers = this->layeredCS->giveNumberOfLayers();  
        for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRule = integrationRulesArray [ layer - 1 ];

            for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
                GaussPoint *gp = iRule->getIntegrationPoint(i);
                std::pair<int, int> pid;
                pid.first  = gp->giveNumber();
                pid.second = giveDelaminationGroupAt( gp->giveCoordinate(3));
                this->gpDelaminationGroupList.push_back(pid); 
            }
        }
    }
}


int
Shell7BaseXFEM :: giveDelaminationGroupAt(double xi) 
{
    // Starts ordering from 0
    std::list< std::pair<int, double> >::const_iterator iter;
    iter = this->delaminationXiCoordList.begin();

    size_t nDelam = this->delaminationXiCoordList.size();
    for ( int j = 1; j <= (int)nDelam; j++ ) {
        double xiDelam = (*iter).second;
        if ( xi < xiDelam ) { //belong to the delamination group just below delamination #j. How to deal with points that lie on the boundary?
            return j-1;
        }
        iter++;
    }
    return (int) nDelam;

}


void 
Shell7BaseXFEM :: giveDelaminationGroupXiLimits(int &dGroup, double &xiTop, double &xiBottom)
{
    size_t nDelam = this->delaminationXiCoordList.size();
    if ( nDelam == 0 ) {
        xiBottom = xiTop = 0.0;
    } else {
        std::list< std::pair<int, double> >::const_iterator iter;
        iter = this->delaminationXiCoordList.begin(); 

        if ( dGroup == 0 ) {
            xiBottom = - this->layeredCS->giveMidSurfaceXiCoordFromBottom();
            xiTop = (*iter).second;
        } else if (dGroup == (int)nDelam) {
            std::advance(iter, dGroup-1);
            xiBottom = (*iter).second;
            xiTop    = -this->layeredCS->giveMidSurfaceXiCoordFromBottom() + 2.0;
        } else {
            std::advance(iter, dGroup-1);
            xiBottom = (*iter).second;
            iter++;
            xiTop    = (*iter).second;
        }

    #if DEBUG
        if ( xiBottom > xiTop ) {
            OOFEM_ERROR2("giveDelaminationGroupZLimits: Bottom xi-coord is larger than top xi-coord in dGroup. (%i)", dGroup);
        }
    #endif
    }
}


double 
Shell7BaseXFEM :: giveDelaminationGroupMidXi(int dGroup)
{
    double xiTop=0., xiBottom=0.;
    this->giveDelaminationGroupXiLimits(dGroup, xiTop, xiBottom);
    return 0.5 * ( xiTop + xiBottom );
}

#endif


void
Shell7BaseXFEM :: vtkEvalUpdatedGlobalCoordinateAt(FloatArray &localCoords, int layer, FloatArray &globalCoords, TimeStep *tStep)
{
    double zeta = localCoords.at(3)*this->layeredCS->giveLayerThickness(layer)*0.5 + this->layeredCS->giveLayerMidZ(layer);

    // Continuous part
    FloatArray solVec;
    this->giveUpdatedSolutionVector(solVec, tStep);
    FloatArray xc, mc; double gamc=0;
    this->giveUnknownsAt(localCoords, solVec, xc, mc, gamc, tStep); 
    double fac = ( zeta + 0.5 * gamc * zeta * zeta );
    
    globalCoords = xc;    
    globalCoords.add(fac,mc);

    // Discontinuous part
    for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { // Only one is supported at the moment
        Delamination *dei =  dynamic_cast< Delamination * >( this->xMan->giveEnrichmentItem(i) ); // should check success
        //EnrichmentFunction *ef = dei->giveEnrichmentFunction(1);

        for ( int j = 1; j <= dei->giveNumberOfEnrichmentDomains(); j++ ) {
            if ( dei->isElementEnrichedByEnrichmentDomain(this, j) ) {
                double zeta0 = dei->enrichmentDomainXiCoords.at(j)*this->layeredCS->computeIntegralThick()*0.5;
                double H = dei->heaviside(zeta, zeta0);

                if ( H > 0.1 ) {
                    FloatArray solVecD;
                    IntArray eiDofIdArray;
                    dei->giveEIDofIdArray(eiDofIdArray, j);
                    this->giveSolutionVector(solVecD, eiDofIdArray, tStep); 
                    FloatArray xd, md; double gamd=0;
                    this->giveUnknownsAt(localCoords, solVecD, xd, md, gamd, tStep);
                    double fac = ( zeta + 0.5 * gamd * zeta * zeta );
                    FloatArray xtemp(3);
                    xtemp = xd;
                    xtemp.add(fac,md);
                    globalCoords.add(xtemp); 
                }

            }
        }
    }

}

} // end namespace oofem
