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
};



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
    
    return (xiRef - xiMid)*layeredCS->computeIntegralThick(); // new xi-coord measured from dGroup c.s. 
    
}




void
Shell7BaseXFEM :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
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
    if ( enrichmentdomainNumber == 0 ) {
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
Shell7BaseXFEM :: evalCovarBaseVectorsAt(GaussPoint *gp, FloatArray &g1, FloatArray &g2, FloatArray &g3, FloatArray &genEpsC)
{
    // Continuous part
    FloatArray g1c, g2c, g3c;
    Shell7Base :: evalCovarBaseVectorsAt(gp, g1, g2, g3, genEpsC);

    // Discontinuous part - ///@todo naive imlementation should be changed
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { // Only one is supported at the moment
        Delamination *dei =  dynamic_cast< Delamination * >( this->xMan->giveEnrichmentItem(i) ); // should check success
        EnrichmentFunction *ef = dei->giveEnrichmentFunction(1);

        for ( int j = 1; j <= dei->giveNumberOfEnrichmentDomains(); j++ ) {
            if ( dei->isElementEnrichedByEnrichmentDomain(this, j) ) {
                double xi0 = dei->enrichmentDomainXiCoords.at(j);
                double H = dei->heaviside(gp->giveCoordinate(3), xi0);        
                
                FloatArray g1d_temp, g2d_temp, g3d_temp, dGenEps;
                computeDiscGeneralizedStrainComponents(dGenEps, gp, dei, j, tStep);
                Shell7Base :: evalCovarBaseVectorsAt(gp, g1d_temp, g2d_temp, g3d_temp, dGenEps);

                g1.add(H,g1d_temp); g2.add(H,g2d_temp); g3.add(H,g3d_temp);
            }
        }
    }
}


void
Shell7BaseXFEM :: computeDiscGeneralizedStrainComponents(FloatArray &answer, GaussPoint *gp, EnrichmentItem *ei, int enrichmentDomainNumber, TimeStep *tStep)
{
    FloatArray dSolVec;
    IntArray eiDofIdArray;
    ei->giveEIDofIdArray(eiDofIdArray, enrichmentDomainNumber);
    this->discGiveUpdatedSolutionVector(dSolVec, eiDofIdArray, tStep);
    FloatMatrix B11, B22, B32, B43, B53;
    Shell7Base :: computeBmatricesAt(gp, B11, B22, B32, B43, B53);
    Shell7Base :: computeGeneralizedStrainVector(answer, dSolVec, B11, B22, B32, B43, B53);
}



void
Shell7BaseXFEM :: discGiveUpdatedSolutionVector(FloatArray &answer, IntArray &eiDofIdArray, TimeStep *tStep)
{
    // Returns the solution vector of discontinuous dofs dx_d & dm_d
    FloatArray temp;
    temp_computeVectorOf(eiDofIdArray, VM_Total, tStep, temp);
    answer.resize( Shell7Base :: giveNumberOfDofs() );
    answer.assemble( temp, this->giveOrdering(AllInv) );
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
Shell7BaseXFEM :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
//
// Computes internal forces as a summation of: sectional forces + convective mass force
{

    answer.resize( this->giveNumberOfDofs() );
    answer.zero();
    FloatArray solVec;
    FloatArray temp, tempRed;

    // Continuous part
    this->giveUpdatedSolutionVector(solVec, tStep);
    this->computeSectionalForces(temp, tStep, solVec, useUpdatedGpRecord);
    IntArray activeDofs, ordering; 
    this->computeOrderingArray(ordering, activeDofs, 0, All);


    // Disccontinuous part
    for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { // Only one is supported at the moment
        Delamination *dei =  dynamic_cast< Delamination * >( this->xMan->giveEnrichmentItem(i) ); // should check success

        for ( int j = 1; j <= dei->giveNumberOfEnrichmentDomains(); j++ ) {
            if ( dei->isElementEnrichedByEnrichmentDomain(this, j) ) {
                IntArray eiDofIdArray;
                dei->giveEIDofIdArray(eiDofIdArray, j);
                this->discGiveUpdatedSolutionVector(solVec, eiDofIdArray, tStep);
                double xi0 = dei->enrichmentDomainXiCoords.at(j);
                
                this->discComputeSectionalForces(temp, tStep, solVec, useUpdatedGpRecord, xi0, dei, j);                      

                // Assemble
                this->computeOrderingArray(ordering, activeDofs,  j, All);
                tempRed.beSubArrayOf(temp, activeDofs);
                answer.assemble(tempRed, ordering);

            }
        }
    }
}


void 
Shell7BaseXFEM :: computeOrderingArray( IntArray &orderingArray, IntArray &activeDofsArray,  int enrichmentDomainNumber, SolutionField field)
{
    // Routine to extract vector given an array of dofid items
    // If a certain dofId does not exist a zero is used as value
    //IntArray eiDofIdArray;
    //this->xMan->giveEnrichmentItem(1)->giveEIDofIdArray(eiDofIdArray, enrichmentDomainNumber);
    //ei->giveEIDofIdArray(eiDofIdArray, enrichmentDomainNumber);

    IntArray ordering_cont = this->giveOrdering(field);
    IntArray fieldDofId    = this->giveFieldDofId(field);

    IntArray ordering_temp, activeDofsArrayTemp;
    ordering_temp.resize(ordering_cont.giveSize());
    activeDofsArrayTemp.resize(ordering_cont.giveSize());

    int activeDofPos = 0, activeDofIndex = 0, orderingDofIndex = 0;
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        DofManager *dMan = this->giveDofManager(i);
        IntArray dofManDofIdMask, dofManDofIdMaskAll; 
        this->discGiveDofManDofIDMask(i, enrichmentDomainNumber, dofManDofIdMask);

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
    
    // Reduce arrays to actual size ///@todo will not work if there are several ei
    int numActiveDofs = activeDofPos;
    IntArray ordering; orderingArray.resize(numActiveDofs), activeDofsArray.resize(numActiveDofs);
    
    for ( int i = 1; i <= numActiveDofs; i++ ) {
        orderingArray.at(i) = ordering_temp.at(i); 
        activeDofsArray.at(i) = activeDofsArrayTemp.at(i);
    }

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

        for ( int j = 1; j <= dei->giveNumberOfEnrichmentDomains(); j++ ) {
            if ( dei->isElementEnrichedByEnrichmentDomain(this, j) ) {
                IntArray eiDofIdArray;
                FloatArray solVecJ;
                dei->giveEIDofIdArray(eiDofIdArray, j);
                this->discGiveUpdatedSolutionVector(solVecJ, eiDofIdArray, tStep);
                for ( int k = j; k <= dei->giveNumberOfEnrichmentDomains(); k++ ) {
                    if ( dei->isElementEnrichedByEnrichmentDomain(this, k) ) {
                        FloatArray solVecK;
                        dei->giveEIDofIdArray(eiDofIdArray, k);
                        this->discGiveUpdatedSolutionVector(solVecK, eiDofIdArray, tStep);

                        discComputeBulkTangentMatrix(temp, solVec, solVecJ, solVecK, rMode, tStep, dei, j, k);


                        // Assemble part correpsonding to active dofs
                        IntArray orderingJ, orderingK, activeDofsJ, activeDofsK;
                        computeOrderingArray(orderingJ, activeDofsJ, j, All);
                        computeOrderingArray(orderingK, activeDofsK, k, All);

                        FloatMatrix tempRed;
                        tempRed.beSubMatrixOf(temp, activeDofsJ, activeDofsK);
                        answer.assemble(tempRed, orderingJ, orderingK);
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
                this->discGiveUpdatedSolutionVector(solVecJ, eiDofIdArray, tStep);                
                discComputeBulkTangentMatrix(temp, solVec, solVecJ, solVec, rMode, tStep, dei, j, 0);


                // Assemble part correpsonding to active dofs
                IntArray orderingJ, activeDofsJ;
                computeOrderingArray(orderingJ, activeDofsJ, j, All);

                FloatMatrix tempRed;
                tempRed.beSubMatrixOf(temp, activeDofsJ, activeDofs);
                answer.assemble(tempRed, orderingJ, ordering);
            }
            
        }

         // First row

        for ( int k = 1; k <= dei->giveNumberOfEnrichmentDomains(); k++ ) {
            if ( dei->isElementEnrichedByEnrichmentDomain(this, k) ) {
                IntArray eiDofIdArray;
                FloatArray solVecK;
                dei->giveEIDofIdArray(eiDofIdArray, k);
                this->discGiveUpdatedSolutionVector(solVecK, eiDofIdArray, tStep);                
                discComputeBulkTangentMatrix(temp, solVec, solVec, solVecK, rMode, tStep, dei, 0, k);


                // Assemble part correpsonding to active dofs
                IntArray orderingK, activeDofsK;
                computeOrderingArray(orderingK, activeDofsK, k, All);

                FloatMatrix tempRed;
                tempRed.beSubMatrixOf(temp, activeDofs, activeDofsK);
                answer.assemble(tempRed, ordering, orderingK);
            }
            
        }



    }
}


void
Shell7BaseXFEM :: discComputeBulkTangentMatrix(FloatMatrix &answer, FloatArray &solVec,  FloatArray &solVecI, FloatArray &solVecJ, MatResponseMode rMode, TimeStep *tStep,
                    EnrichmentItem *ei, int enrichmentDomainNumberI, int enrichmentDomainNumberJ)
{
 
    FloatMatrix A [ 3 ] [ 3 ], lambdaI [ 3 ], lambdaJ [ 3 ];
    FloatMatrix L(18,18);
    FloatMatrix B11, B22, B32, B43, B53, B;
    FloatArray S1g(3), S2g(3), S3g(3);
    FloatMatrix K(42,42);
    K.zero();

    //int ndofs = this->giveNumberOfDofs();
    int ndofs = Shell7Base :: giveNumberOfDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();

    Delamination *dei =  dynamic_cast< Delamination * >( ei );
    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRule = layerIntegrationRulesArray [ layer - 1 ];
        Material *mat = domain->giveMaterial( this->layeredCS->giveLayerMaterial(layer) );

        for ( int i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
            GaussPoint *gp = iRule->getIntegrationPoint(i);
            double zeta = giveGlobalZcoord(gp);
            double xi0I = 0.0; 
            double xi0J = 0.0;
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
                answer.add(dV, K);
            }
        }
    }



    IntArray ordering_phibar, ordering_m, ordering_gam;
    IntArray activeDofs_phibar, activeDofs_m, activeDofs_gam;

    //computeOrderingArray(ordering_phibar, activeDofs_phibar, ei, enrichmentDomainNumber, Midplane);
    //computeOrderingArray(ordering_m, activeDofs_m, ei, enrichmentDomainNumber, Director);
    //computeOrderingArray(ordering_gam, activeDofs_gam, ei, enrichmentDomainNumber, InhomStrain);


    /*
    FloatMatrix K11Ass, K12Ass, K13Ass, K22Ass, K23Ass, K33Ass, mat1, mat2(3,3);

    answer.assemble(K11, ordering_phibar, ordering_phibar);
    answer.assemble(K12, ordering_phibar, ordering_m);
    answer.assemble(K13, ordering_phibar, ordering_gam);
    answer.assemble(K22, ordering_m,      ordering_m);
    answer.assemble(K23, ordering_m,      ordering_gam);
    answer.assemble(K33, ordering_gam,    ordering_gam);

    FloatMatrix K21, K31, K32;
    K21.beTranspositionOf(K12);
    K31.beTranspositionOf(K13);
    K32.beTranspositionOf(K23);
    answer.assemble(K21, ordering_m,      ordering_phibar);
    answer.assemble(K31, ordering_gam,    ordering_phibar);
    answer.assemble(K32, ordering_gam,    ordering_m);
    */
}



void
Shell7BaseXFEM :: discComputeSectionalForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, int useUpdatedGpRecord, 
                    double xi0, EnrichmentItem *ei, int enrichmentDomainNumber)
//
{
    FloatMatrix B;
    FloatArray BtF, genEps;
    answer.resize( Shell7Base :: giveNumberOfDofs() );
    answer.zero();

    
    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >( this->giveCrossSection() );
    int numberOfLayers = layeredCS->giveNumberOfLayers();     // conversion of types
    FloatArray f1(18), f2(18), f3(6);
    f1.zero();
    f2.zero();
    f3.zero();
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRuleL = layerIntegrationRulesArray [ layer - 1 ];
        Material *mat = domain->giveMaterial( layeredCS->giveLayerMaterial(layer) );

        for ( int j = 1; j <= iRuleL->getNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRuleL->getIntegrationPoint(j - 1);

            if ( gp->giveCoordinate(3) > xi0 ) { // Should be enriched ///@todo not general!

                FloatMatrix B11, B22, B32, B43, B53;
                this->computeBmatricesAt(gp, B11, B22, B32, B43, B53);
                this->computeGeneralizedStrainVector(genEps, solVec, B11, B22, B32, B43, B53);

                double zeta = giveGlobalZcoord(gp);
                FloatArray N, M, T, Ms;
                double Ts = 0.;
                this->computeSectionalForcesAt(N, M, T, Ms, Ts, gp, mat, tStep, genEps, zeta);

                // Computation of sectional forces: f = B^t*[N M T Ms Ts]^t
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


    // Should assemble to xfem dofs
    //IntArray ordering_phibar, ordering_m, ordering_gam, ordering_all;
    //IntArray activeDofs_phibar, activeDofs_m, activeDofs_gam, activeDofs_all;



    //computeOrderingArray(ordering_phibar, activeDofs_phibar, ei, enrichmentDomainNumber, Midplane);
    //computeOrderingArray(ordering_m, activeDofs_m, ei, enrichmentDomainNumber, Director);
    //computeOrderingArray(ordering_gam, activeDofs_gam, ei, enrichmentDomainNumber, InhomStrain);
    //computeOrderingArray(ordering_all, activeDofs_all, ei, enrichmentDomainNumber, All);
/*
    ordering_phibar.printYourself();
    ordering_m.printYourself();
    ordering_gam.printYourself();

    activeDofs_phibar.printYourself();
    activeDofs_m.printYourself();
    activeDofs_gam.printYourself();
    */

    /*
    FloatArray f1Ass, f2Ass, f3Ass;

    f1Ass.beSubArrayOf(f1,activeDofs_phibar);
    f2Ass.beSubArrayOf(f2,activeDofs_m);
    f3Ass.beSubArrayOf(f3,activeDofs_gam);
    answer.assemble(f1Ass, ordering_phibar);
    answer.assemble(f2Ass, ordering_m);
    answer.assemble(f3Ass, ordering_gam);
    */

    //FloatArray f(42);
    answer.addSubVector(f1,1);
    answer.addSubVector(f2,19);
    answer.addSubVector(f3,37);
    //FloatArray fAss;
    //fAss.beSubArrayOf(f,activeDofs_all);
    //answer.assemble(fAss, ordering_all);
}






void
Shell7BaseXFEM :: computeMassMatrixNum(FloatMatrix &answer, TimeStep *tStep) {
    // Num refers in this case to  numerical integration in both in-plane and through the thickness.
    // For analytically integrated throught he thickness, see computeMassMatrix


    FloatMatrix N, Nt, Ntm, NtmN, mass, temp;
    FloatArray solVec, unknowns;
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
        IntegrationRule *iRuleL = layerIntegrationRulesArray [ layer - 1 ];
        Material *mat = domain->giveMaterial( layeredCS->giveLayerMaterial(layer) );

        for ( int j = 1; j <= iRuleL->getNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRuleL->getIntegrationPoint(j - 1);

            FloatMatrix N11, N22, N33;
            this->computeNmatricesAt(gp, N11, N22, N33);
            FloatArray xbar, m;
            double gam = 0.;
            this->computeSolutionFields(xbar, m, gam, solVec, N11, N22, N33);
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

        IntArray ordering_phibar = giveOrdering(Midplane);
        IntArray ordering_m = giveOrdering(Director);
        IntArray ordering_gam = giveOrdering(InhomStrain);
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






# if 0
void
Shell7BaseXFEM :: discEvalCovarBaseVectorsAt(GaussPoint *gp, FloatArray &g1d, FloatArray &g2d, FloatArray &g3d, FloatArray &dGenEps)
{

    FloatArray lcoords = * gp->giveCoordinates();
    double zeta = giveGlobalZcoord(gp);

    FloatArray dxdxi1, dxdxi2, dmdxi1, dmdxi2, m;
    this->discGiveGeneralizedStrainComponents(dGenEps, dxdxi1, dxdxi2, dmdxi1, dmdxi2, m);

    g1d = dxdxi1;
    g1d.add(zeta, dmdxi1);
    g2d = dxdxi2;
    g2d.add(zeta, dmdxi2);
    g3d = m;
    
}



void
Shell7BaseXFEM :: discComputeGeneralizedStrainVector(FloatArray &answer, const FloatArray &solVec, const FloatMatrix &B11,
                                                     const FloatMatrix &B22, const FloatMatrix &B32) {
    answer.resize(15);
    answer.zero();
    int ndofs_xm  = this->giveNumberOfFieldDofs(Midplane);
    int ndofs_gam = this->giveNumberOfFieldDofs(InhomStrain);
    for ( int i = 1; i <= ndofs_xm; i++ ) {
        for ( int j = 1; j <= 6; j++ ) {
            answer.at(j)  += B11.at(j, i) * solVec.at(i);            // dx/dxi
            answer.at(6 + j)  += B22.at(j, i) * solVec.at(i + ndofs_xm);      // dm/dxi
        }
    }

    for ( int i = 1; i <= ndofs_xm; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            answer.at(12 + j) += B32.at(j, i) * solVec.at(i + ndofs_xm);      // m
        }
    }

}


void
Shell7BaseXFEM :: discGiveGeneralizedStrainComponents(FloatArray &genEps, FloatArray &dphidxi1, FloatArray &dphidxi2, FloatArray &dmdxi1, 
         FloatArray &dmdxi2, FloatArray &m) {
    // generealized strain vector for discontinuous part  [dxdxi, dmdxi, m]^T
    dphidxi1.setValues( 3, genEps.at(1), genEps.at(2), genEps.at(3) );
    dphidxi2.setValues( 3, genEps.at(4), genEps.at(5), genEps.at(6) );
    dmdxi1.setValues( 3, genEps.at(7), genEps.at(8), genEps.at(9) );
    dmdxi2.setValues( 3, genEps.at(10), genEps.at(11), genEps.at(12) );
    m.setValues( 3, genEps.at(13), genEps.at(14), genEps.at(15) );
}
#endif









IntArray
Shell7BaseXFEM :: giveFieldDofId(SolutionField fieldType) const {
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
        _error("giveOrdering: unknown fieldType");
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
{   // Creates a list wich stores the gp# and the dGroup# for quick access later.
  if ( this->gpDelaminationGroupList.size()==0 ) {
    int numberOfLayers = this->layeredCS->giveNumberOfLayers();  
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRule = layerIntegrationRulesArray [ layer - 1 ];

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

    int nDelam = this->delaminationXiCoordList.size();
    for ( int j = 1; j <= nDelam; j++ ) {
        double xiDelam = (*iter).second;
        if ( xi < xiDelam ) { //belong to the delamination group just below delamination #j. How to deal with poins that lie onthe boundary?
            return j-1;            
        }
        iter++;
    }
    return nDelam;
            
}


void 
Shell7BaseXFEM :: giveDelaminationGroupXiLimits(int &dGroup, double &xiTop, double &xiBottom)
{
    int nDelam = this->delaminationXiCoordList.size();
    std::list< std::pair<int, double> >::const_iterator iter;
    iter = this->delaminationXiCoordList.begin(); 

    if ( dGroup == 0 ) {
        xiBottom = - this->layeredCS->giveMidSurfaceXiCoordFromBottom();
        xiTop = (*iter).second;
    } else if (dGroup == nDelam) {
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


double 
Shell7BaseXFEM :: giveDelaminationGroupMidXi(int dGroup)
{
    double xiTop=0., xiBottom=0.;
    this->giveDelaminationGroupXiLimits(dGroup, xiTop, xiBottom);
    return 0.5 * ( xiTop + xiBottom );
}

#endif














void 
Shell7BaseXFEM :: setupDelaminationXiCoordsAtGP() 
{
    std::pair<int, double> pid;
    std::list<std::pair<int, double> > *delaminationXiCoordList;

    xMan = this->giveDomain()->giveXfemManager(1);
    int numEI = xMan->giveNumberOfEnrichmentItems();
    for ( int i = 1; i <= numEI; i++ ) {
        Delamination *dei =  dynamic_cast< Delamination * >( xMan->giveEnrichmentItem(i) );
        if ( dei ) {
            if ( dei->isElementEnriched(this) ) {


                int nDelam = dei->giveNumberOfEnrichmentDomains(); // numEnrDomains max possible number
                int pos = 1;
                for ( int j = 1; j <= nDelam; j++ ) {
                    //if( this->isElementEnriched(element) ) {
                    //pid.first  = pos;
                    //pid.second = this->delaminationZCoords.at(i); 
                    //xiCoordList->push_back(pid); 
                    //pos++;
                }
            }
        } 
    }
}








} // end namespace oofem

