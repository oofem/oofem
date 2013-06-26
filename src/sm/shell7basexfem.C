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
#include "dofmanager.h"
#include "constantpressureload.h"
#include "simpleinterfacemat.h"
namespace oofem {
IntArray Shell7BaseXFEM :: dofId_Midplane(3);
IntArray Shell7BaseXFEM :: dofId_Director(3);
IntArray Shell7BaseXFEM :: dofId_InhomStrain(1); 
bool Shell7BaseXFEM :: __initializedFieldDofId = Shell7BaseXFEM :: initDofId();

Shell7BaseXFEM :: Shell7BaseXFEM(int n, Domain *aDomain) : Shell7Base(n, aDomain), XfemElementInterface(this) 
{
    czMat = NULL;
}

int 
Shell7BaseXFEM :: checkConsistency()
{
    Shell7Base :: checkConsistency();

    // check if defined xi-coords in delamination EI corresponds to actual layer boundaries
    //@todo should be improved
#if 1
    int numLayers = this->layeredCS->giveNumberOfLayers();
    FloatArray interfaceXi(numLayers-1), delamXiCoords;
    this->layeredCS->giveInterfaceXiCoords(interfaceXi);
    IntArray DelaminatedInterfaceList;
    this->computeDelaminatedInterfaceList(DelaminatedInterfaceList);
    for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { 
        EnrichmentItem *ei = this->xMan->giveEnrichmentItem(i);
        if ( Delamination *dei = dynamic_cast< Delamination * >( ei ) ) {
            dei->giveActiveDelaminationXiCoords(delamXiCoords, this);

            if( delamXiCoords.giveSize() && DelaminatedInterfaceList.giveSize() == 0  ) {
                OOFEM_ERROR("Shell7BaseXFEM :: checkConsistency - Delamination xi-coords in input does not correspond to layer interfaces");
            }

            //for ( int j = 1; j <= delamXiCoords.giveSize(); j++ ) {
            //bool flag = false;
            //    for ( int k = 1; k <= interfaceXi.giveSize(); k++ ) {
            //        if ( abs(delamXiCoords.at(j)-interfaceXi.at(k)) < 1.0e-3 ) {
            //            flag = true;
            //            //this->DelaminatedInterfaceList.followedBy(k); //check for duplicates? 
            //        }
            //    }
            //    if ( !flag ) {
            //        OOFEM_ERROR2("checkConsistency: Delamination xi-coord (%e) does not correspond to a layer interface (check input).", delamXiCoords.at(j));
            //    }
            //    
            //}
        }
    }
#endif
    return 1;
}

void
Shell7BaseXFEM :: postInitialize()
{
    Shell7Base :: postInitialize();

    this->xMan =  this->giveDomain()->giveXfemManager();
    if ( this->czMatNum > 0 ) {
        this->czMat = this->giveDomain()->giveMaterial(this->czMatNum);
    }

}

void
Shell7BaseXFEM :: updateYourself(TimeStep *tStep)
{
    StructuralElement :: updateYourself(tStep);
}


void
    Shell7BaseXFEM :: computeDelaminatedInterfaceList(IntArray &list) 
{
    int numLayers = this->layeredCS->giveNumberOfLayers();
    FloatArray interfaceXi(numLayers-1), delamXiCoords;
    this->layeredCS->giveInterfaceXiCoords(interfaceXi); // coords defining each layer interface in the cross section
    
    list.resize(0);
    for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { 
        EnrichmentItem *ei = this->xMan->giveEnrichmentItem(i);
        if ( Delamination *dei = dynamic_cast< Delamination * >( ei ) ) {
            dei->giveActiveDelaminationXiCoords(delamXiCoords, this);
            for ( int j = 1; j <= delamXiCoords.giveSize(); j++ ) {
                for ( int k = 1; k <= interfaceXi.giveSize(); k++ ) {
                    if ( abs(delamXiCoords.at(j)-interfaceXi.at(k)) < 1.0e-3 ) {
                        list.followedBy(k); //check for duplicates? 
                    }                    
                }
            }
        }
    }
}

void
Shell7BaseXFEM :: computeFailureCriteriaQuantities(FailureCriteria *fc, TimeStep *tStep) 
{
    // Compute necessary quantities for evaluation of failure criterias
    switch ( fc->giveType() ) {

    case FC_MaxShearStress:
        IntArray delaminatedInterfaceList;
        this->computeDelaminatedInterfaceList(delaminatedInterfaceList);

        // Stress ordering (1, 5, 9, 6, 3, 2) = (xx, yy, zz, yz, xz, xy)
        int numInterfaces = this->layeredCS->giveNumberOfLayers() - 1;
        std::vector < FloatArray > interLamStresses;
        fc->quantities.resize(numInterfaces); // will overwrite this every time
        for (int intface = 1; intface <= numInterfaces; intface++ ) { 
            
            if ( ! delaminatedInterfaceList.contains(intface) ) { // interface is whole
                this->computeInterLaminarStressesAt(intface, tStep, interLamStresses); // all 6 components in each evaluation point (ip)
                int numEvalPoints = interLamStresses.size();
                fc->quantities[intface-1].resize(numEvalPoints); // most often = numIP
            
                for ( int evalPoint = 1; evalPoint <= numEvalPoints; evalPoint++) {
                    FloatArray &values = fc->quantities[intface-1][evalPoint-1];  // one resulting shear stress
                    FloatArray &vS = interLamStresses[evalPoint-1];               // Stress in eval point
                    values.resize(1);                                             // scalar measure in this case
                    values.at(1) = sqrt( vS.at(2)*vS.at(2) + vS.at(3)*vS.at(3) ); // components can't be right here? shouldn't it be a traction vector?

                }
            }
        }

    };
}


IRResultType Shell7BaseXFEM :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro
    
    int material = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, material, _IFT_Shell7BaseXFEM_CohesiveZoneMaterial);
    this->czMatNum = material;

    Shell7Base :: initializeFrom(ir);
    return IRRT_OK; 
}

bool 
Shell7BaseXFEM :: hasCohesiveZone()
{
    if ( this->czMat !=NULL ) {
        return true;
    } else {
        return false;
    }
}

Interface
*Shell7BaseXFEM :: giveInterface(InterfaceType it)
{
    if ( it != XfemElementInterfaceType ) {
        return Shell7Base :: giveInterface(it);
    } else if ( it == XfemElementInterfaceType ) {
        return ( XfemElementInterface * ) this;
    } else {
        return Shell7Base :: giveInterface(it); //@todo remove
    }
}


double 
Shell7BaseXFEM :: giveGlobalZcoord(double xi) 
{
    
    this->setupDelaminationXiCoordList();
    //this->setupGPDelaminationGroupList();

    //double xiRef = gp->giveCoordinate(3);
    double xiRef = xi;
    int dGroup   = this->giveDelaminationGroupAt( xiRef );
    double xiMid = this->giveDelaminationGroupMidXi(dGroup);
    
    //return (xiRef - xiMid)*layeredCS->computeIntegralThick()*0.5; // new xi-coord measured from dGroup c.s. 
    return (xiRef+0. )*this->layeredCS->computeIntegralThick()*0.5; // new xi-coord measured from dGroup c.s. 
    
}


void
Shell7BaseXFEM :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    // Returns the total id mask of the dof manager - regular id's + enriched id's

    // Continuous part
    Shell7Base ::giveDofManDofIDMask(inode, ut, answer);
    XfemManager *xMan = this->giveDomain()->giveXfemManager();
    // Discontinuous part
    DofManager *dMan = this->giveDofManager(inode);
    //for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { // Only one is supported at the moment
    for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) { // Only one is supported at the moment
        //EnrichmentItem *ei = this->xMan->giveEnrichmentItem(i);
        EnrichmentItem *ei = xMan->giveEnrichmentItem(i);
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
Shell7BaseXFEM :: evalCovarBaseVectorsAt(FloatArray &lCoords, FloatMatrix &gcov, FloatArray &genEpsC)
{
    // Continuous part
    Shell7Base :: evalCovarBaseVectorsAt(lCoords, gcov, genEpsC);

    // Discontinuous part - ///@todo bad implementation regarding enr. functions - should be changed
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { // Only one is supported at the moment
        Delamination *dei =  dynamic_cast< Delamination * >( this->xMan->giveEnrichmentItem(i) ); // should check success
        //EnrichmentFunction *ef = dei->giveEnrichmentFunction(1);

        for ( int j = 1; j <= dei->giveNumberOfEnrichmentDomains(); j++ ) {
            if ( dei->isElementEnrichedByEnrichmentDomain(this, j) ) {
                double xi0 = dei->enrichmentDomainXiCoords.at(j);
                double H   = dei->heaviside(lCoords.at(3), xi0);        
                if ( H > 0.1 ) {
                    FloatArray dGenEps;
                    computeDiscGeneralizedStrainVector(dGenEps, lCoords, dei, j, tStep);
                    FloatMatrix gcovd; 
                    Shell7Base :: evalCovarBaseVectorsAt(lCoords, gcovd, dGenEps);
                    gcov.add(H,gcovd); 
                }

            }
        }
    }
}


void
Shell7BaseXFEM :: computeDiscGeneralizedStrainVector(FloatArray &answer, FloatArray &lCoords, EnrichmentItem *ei, int enrichmentDomainNumber, TimeStep *tStep)
{
    FloatArray solVecD;
    IntArray eiDofIdArray;
    ei->giveEIDofIdArray(eiDofIdArray, enrichmentDomainNumber);
    this->giveSolutionVector(solVecD, eiDofIdArray, tStep);  
    FloatMatrix B;
    Shell7Base :: computeBmatrixAt(lCoords, B, 0, 0);
    Shell7Base :: computeGeneralizedStrainVectorNew(answer, solVecD, B);
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

                // Cohesive zone model
                if ( this->hasCohesiveZone() ) {
                    FloatArray fCZ;
                    this->computeCohesiveForces( fCZ, tStep, solVec, solVecD, useUpdatedGpRecord, dei, j);
                    tempRed.beSubArrayOf(fCZ, activeDofs);
                    answer.assemble(tempRed, ordering);
                }
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

    int ndofs = Shell7Base :: giveNumberOfDofs();
    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     // conversion of types
    FloatArray f(ndofs);
    FloatArray genEps, genEpsD;
    FloatMatrix B;
    FloatArray ftemp, lCoords;

    f.zero();
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRuleL = integrationRulesArray [ layer - 1 ];
        Material *mat = domain->giveMaterial( this->layeredCS->giveLayerMaterial(layer) );

        for ( int j = 1; j <= iRuleL->giveNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRuleL->getIntegrationPoint(j - 1);
            lCoords = *gp->giveCoordinates();

            if ( gp->giveCoordinate(3) > xi0 ) { // Should be enriched ///@todo not general!

                this->computeBmatrixAt(lCoords, B);
                this->computeGeneralizedStrainVectorNew(genEps,  solVec,  B);
                this->computeGeneralizedStrainVectorNew(genEpsD, solVecD, B);

                double zeta = giveGlobalZcoord(gp->giveCoordinate(3));

                FloatArray sectionalForces;
                this->computeSectionalForcesAt(sectionalForces, gp, mat, tStep, genEps, genEpsD, zeta);

                // Computation of nodal forces: f = B^t*[N M T Ms Ts]^t
                ftemp.beTProductOf(B,sectionalForces);
                double dV = this->computeVolumeAroundLayer(gp, layer);
                f.add(dV, ftemp);
            
            }

        }
    }

    answer.resize(ndofs);  answer.zero();
    const IntArray &ordering_all = this->giveOrdering(All);
    answer.assemble(f, ordering_all);

}



void
Shell7BaseXFEM :: computeCohesiveForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, FloatArray &solVecD, int useUpdatedGpRecord, 
                     Delamination *dei, int enrichmentDomainNumber)
{
    //Computes the cohesive nodal forces for a given interface
    FloatArray answerTemp, lCoords;
    answerTemp.resize(Shell7Base :: giveNumberOfDofs() ); 
    answerTemp.zero();

    FloatMatrix N, B;  

    IntegrationRule *iRuleL = czIntegrationRulesArray [ enrichmentDomainNumber - 1 ];
    StructuralMaterial *mat = static_cast < StructuralMaterial * > (this->czMat);
    FloatMatrix lambdaN;
    FloatArray Fp, cTraction;
    for ( int i = 1; i <= iRuleL->giveNumberOfIntegrationPoints(); i++ ) {
        IntegrationPoint *ip = iRuleL->getIntegrationPoint(i - 1);
        lCoords = *ip->giveCoordinates();
        this->computeBmatrixAt(lCoords, B);
        this->computeNmatrixAt(lCoords, N);

        // Lambda matrix
        FloatArray genEpsD;
        genEpsD.beProductOf(B, solVecD);
        double xi = dei->enrichmentDomainXiCoords.at(enrichmentDomainNumber);
        double zeta = xi * this->layeredCS->computeIntegralThick() * 0.5;
        FloatMatrix lambda;
        this->computeLambdaNMatrix(lambda, genEpsD, zeta);
        
        // Compute jump vector
        FloatArray xd, unknowns;
        unknowns.beProductOf(N, solVecD);
        xd.beProductOf(lambda,unknowns); // spatial jump
       
        // Compute cohesive traction based on jump
        mat->giveRealStressVector(cTraction, FullForm, ip, xd, tStep);
        lambdaN.beProductOf(lambda,N);
        Fp.beTProductOf(lambdaN, cTraction);
        double dA = this->computeAreaAround(ip,xi);
        answerTemp.add(dA,Fp);
    }
    int ndofs = Shell7Base ::giveNumberOfDofs();
    answer.resize(ndofs);
    answer.zero();
    const IntArray &ordering = this->giveOrdering(All);
    answer.assemble(answerTemp, ordering);

}



void
Shell7BaseXFEM :: computeCohesiveTangent(FloatMatrix &answer, TimeStep *tStep)
{
    //Computes the cohesive tangent forces for a given interface
    FloatArray solVecD;
    FloatMatrix temp;
    int ndofs = this->giveNumberOfDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();

    IntArray eiDofIdArray, orderingJ, activeDofsJ;
    FloatMatrix tempRed;
    // Disccontinuous part (continuous part does not contribute)
    for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { // Only one is supported at the moment
        if ( Delamination *dei =  dynamic_cast< Delamination * >( this->xMan->giveEnrichmentItem(i) ) ) {
            for ( int j = 1; j <= dei->giveNumberOfEnrichmentDomains(); j++ ) {
                if ( dei->isElementEnrichedByEnrichmentDomain(this, j) ) {
                    dei->giveEIDofIdArray(eiDofIdArray, j);
                    this->giveSolutionVector(solVecD, eiDofIdArray, tStep);  
                    this->computeCohesiveTangentAt(temp, tStep, solVecD, dei, j);
                    // Assemble part correpsonding to active dofs
                    computeOrderingArray(orderingJ, activeDofsJ, j, All);
                    tempRed.beSubMatrixOf(temp, activeDofsJ, activeDofsJ);
                    answer.assemble(tempRed, orderingJ, orderingJ);
                }
            }
        }
    }
}



void
Shell7BaseXFEM :: computeCohesiveTangentAt(FloatMatrix &answer, TimeStep *tStep, FloatArray &solVecD, 
    Delamination *dei, int enrichmentDomainNumber)
{
    //Computes the cohesive tangent for a given interface
    FloatArray lCoords;
    FloatMatrix answerTemp, N, B, lambda, K, lambdaN, temp, tangent;
    int nDofs = Shell7Base :: giveNumberOfDofs();
    answerTemp.resize(nDofs, nDofs); 
    answerTemp.zero();

    IntegrationRule *iRuleL = czIntegrationRulesArray [ enrichmentDomainNumber - 1 ];
    StructuralMaterial *mat = static_cast < StructuralMaterial * > (this->czMat); // CZ-material
    double xi = dei->enrichmentDomainXiCoords.at(enrichmentDomainNumber);
    double zeta = this->giveGlobalZcoord(xi);
    this->computeLambdaNMatrixDis(lambda, zeta);

    for ( int i = 1; i <= iRuleL->giveNumberOfIntegrationPoints(); i++ ) {
        IntegrationPoint *ip = iRuleL->getIntegrationPoint(i - 1);
        lCoords = *ip->giveCoordinates();
        this->computeBmatrixAt(lCoords, B);
        this->computeNmatrixAt(lCoords, N);
        
        mat->giveCharacteristicMatrix(K, FullForm, TangentStiffness, ip, tStep);
        
        this->computeTripleProduct(temp, lambda, K, lambda);
        this->computeTripleProduct(tangent, N, temp, N);
        double dA = this->computeAreaAround(ip,xi);
        answerTemp.add(dA,tangent); 
    }

    answer.resize(nDofs, nDofs); 
    answer.zero();
    const IntArray &ordering = this->giveOrdering(All);
    answer.assemble(answerTemp, ordering, ordering);
}



//---------------------------------------------
// Optimized version
void
Shell7BaseXFEM :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{

    int ndofs = this->giveNumberOfDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();
    
    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     
    FloatMatrix tempRed, tempRedT;
    FloatMatrix KCC, KCD, KDD;
    IntArray orderingC, activeDofsC;
    this->computeOrderingArray(orderingC, activeDofsC, 0, All);
    std::vector<IntArray> orderingArrays;
    std::vector<IntArray> activeDofsArrays;

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRule = integrationRulesArray [ layer - 1 ];
        Material *mat = domain->giveMaterial( this->layeredCS->giveLayerMaterial(layer) );

        for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
            IntegrationPoint *ip = iRule->getIntegrationPoint(i);
            this->discComputeBulkTangentMatrix(KCC, KCD, KDD, ip, mat, layer, tStep);
            double xi = ip->giveCoordinate(3);

            // Continuous part K_{c,c}
            answer.assemble(KCC, orderingC, orderingC);

            // Discontinuous part
            for ( int m = 1; m <= this->xMan->giveNumberOfEnrichmentItems(); m++ ) { // Only one is supported at the moment
                if ( Delamination *dei = dynamic_cast< Delamination * >( this->xMan->giveEnrichmentItem(m) ) ) {
                    int numED = dei->giveNumberOfEnrichmentDomains();
                    if ( orderingArrays.size() == 0 ) {
                        orderingArrays.resize(numED);
                        activeDofsArrays.resize(numED);
                    }

                    for ( int j = 1; j <= numED; j++ ) {
                        if ( dei->isElementEnrichedByEnrichmentDomain(this, j) ) {
                            double xi0J = dei->enrichmentDomainXiCoords.at(j);
                            if( orderingArrays[j-1].giveSize() == 0 ) {
                                this->computeOrderingArray( orderingArrays[j-1], activeDofsArrays[j-1], j, All);
                            }
                            
                            // K_{c,dk} & K_{dk,c}
                            if ( xi > xi0J ) {
                                tempRed.beSubMatrixOf(KCD, activeDofsC, activeDofsArrays[j-1]);
                                answer.assemble(tempRed, orderingC, orderingArrays[j-1]);
                                tempRedT.beTranspositionOf(tempRed);
                                answer.assemble(tempRedT, orderingArrays[j-1], orderingC);
                            }

                            // K_{dk,dl}
                            for ( int k = 1; k <= numED; k++ ) {
                                if ( dei->isElementEnrichedByEnrichmentDomain(this, k) ) {
                                    double xi0K = dei->enrichmentDomainXiCoords.at(k);
                                    if ( xi > xi0J  &&  xi > xi0K  ) {
                                        tempRed.beSubMatrixOf(KDD, activeDofsArrays[j-1], activeDofsArrays[k-1]);
                                        answer.assemble(tempRed, orderingArrays[j-1], orderingArrays[k-1]);
                                    }
                                }
                            }


                        }

                    }
                }
            }

        }
    }


    // Cohesive zones
#if 1
    FloatMatrix Kcoh;
    if ( this->hasCohesiveZone() ) {
        this->computeCohesiveTangent(Kcoh, tStep);
        answer.add(Kcoh);
    }
#endif


    // Add contribution due to pressure load
#if 1

    int nLoads = this->boundaryLoadArray.giveSize() / 2;

    for ( int k = 1; k <= nLoads; k++ ) {     // For each pressure load that is applied
        int load_number = this->boundaryLoadArray.at(2 * k - 1);
        int iSurf = this->boundaryLoadArray.at(2 * k);         // load_id
        Load *load = this->domain->giveLoad(load_number);

        if ( ConstantPressureLoad* pLoad = dynamic_cast< ConstantPressureLoad * >( load ) ) {
            FloatMatrix K_pressure;

            IntegrationRule *iRule = specialIntegrationRulesArray [ 1 ];


            for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
                IntegrationPoint *ip = iRule->getIntegrationPoint(i);
                //this->discComputeBulkTangentMatrix(KCC, KCD, KDD, ip, mat, layer, tStep);
                this->computePressureTangentMatrixDis(KCC, KCD, KDD, ip, load, iSurf, tStep);
                // Continuous part
                answer.assemble(KCC, orderingC, orderingC);

                // Discontinuous part
                for ( int m = 1; m <= this->xMan->giveNumberOfEnrichmentItems(); m++ ) { // Only one is supported at the moment
                    if ( Delamination *dei = dynamic_cast< Delamination * >( this->xMan->giveEnrichmentItem(m) ) ) {
                        int numED = dei->giveNumberOfEnrichmentDomains();

                        for ( int j = 1; j <= numED; j++ ) {
                            if ( dei->isElementEnrichedByEnrichmentDomain(this, j) ) {
                                double xi0J = dei->enrichmentDomainXiCoords.at(j);
                                IntArray &orderingJ = orderingArrays[j-1];
                                IntArray &activeDofsJ = activeDofsArrays[j-1];

                                // con-dis & dis-con
                                if ( pLoad->giveLoadOffset() > xi0J ) {
                                    tempRed.beSubMatrixOf(KCD, activeDofsC, activeDofsJ);
                                    answer.assemble(tempRed, orderingC, orderingJ);
                                    tempRedT.beTranspositionOf(tempRed);
                                    answer.assemble(tempRedT, orderingJ, orderingC);
                                }

                                // dis-dis
                                for ( int k = 1; k <= dei->giveNumberOfEnrichmentDomains(); k++ ) {
                                    if ( dei->isElementEnrichedByEnrichmentDomain(this, k) ) {
                                        double xi0K = dei->enrichmentDomainXiCoords.at(k);
                                        if ( pLoad->giveLoadOffset() > xi0J  &&   pLoad->giveLoadOffset() > xi0K  ) {
                                            IntArray &orderingK = orderingArrays[k-1];
                                            IntArray &activeDofsK = activeDofsArrays[k-1];
                                            tempRed.beSubMatrixOf(KDD, activeDofsJ, activeDofsK);
                                            answer.assemble(tempRed, orderingJ, orderingK);
                                        }
                                    }
                                }


                            }

                        }
                    }
                }

            }
            


        }
    }
#endif
    
}

void
Shell7BaseXFEM :: discComputeBulkTangentMatrix(FloatMatrix &KCC, FloatMatrix &KCD, FloatMatrix &KDD, IntegrationPoint *ip, Material *mat, int layer, TimeStep *tStep)
{
    FloatArray genEpsC;
    FloatMatrix temp, B;
    FloatMatrix A [ 3 ] [ 3 ], lambdaC [ 3 ], lambdaD [ 3 ];
    
    FloatArray solVecC, lCoords;
    this->giveUpdatedSolutionVector(solVecC, tStep);            
    // continuous part
    lCoords = *ip->giveCoordinates();
    this->computeBmatrixAt(lCoords, B, 0, 0);
    this->computeGeneralizedStrainVectorNew(genEpsC, solVecC , B);
    Shell7Base :: computeLinearizedStiffness(ip, mat, tStep, A, genEpsC);
            
    double zeta = giveGlobalZcoord(ip->giveCoordinate(3));
    this->computeLambdaGMatrices(lambdaC, genEpsC, zeta);
    this->computeLambdaGMatricesDis(lambdaD, zeta);
    double dV = this->computeVolumeAroundLayer(ip, layer);

    FloatMatrix LCAC(18,18), LCAD(18,18), LDAD(18,18); 
    LCAC.zero(); LCAD.zero(); LDAD.zero();
    for ( int j = 0; j < 3; j++ ) {
        for ( int k = 0; k < 3; k++ ) {
            this->computeTripleProduct(temp, lambdaC [ j ], A [ j ][ k ], lambdaC [ k ] );
            LCAC.add(dV,temp);
            this->computeTripleProduct(temp, lambdaC [ j ], A [ j ][ k ], lambdaD [ k ] );
            LCAD.add(dV,temp);
            this->computeTripleProduct(temp, lambdaD [ j ], A [ j ][ k ], lambdaD [ k ] );
            LDAD.add(dV,temp);
        }
    }


    FloatMatrix KCCtemp, KCDtemp, KDDtemp;
    FloatMatrix LCACB, LDADB, LCADB;

    this->computeTripleProduct(KCCtemp, B, LCAC, B );    
    this->computeTripleProduct(KCDtemp, B, LCAD, B );    
    this->computeTripleProduct(KDDtemp, B, LDAD, B );    
    
    int ndofs = Shell7Base :: giveNumberOfDofs();
    KCC.resize(ndofs,ndofs); KCD.resize(ndofs,ndofs); KDD.resize(ndofs,ndofs); 
    KCC.zero(); KCD.zero(); KDD.zero();
    const IntArray &ordering = this->giveOrdering(All);
    KCC.assemble(KCCtemp, ordering, ordering);
    KCD.assemble(KCDtemp, ordering, ordering);
    KDD.assemble(KDDtemp, ordering, ordering);

}


void
Shell7BaseXFEM :: computeLambdaGMatricesDis(FloatMatrix lambda [ 3 ], double zeta)
{
    // computes the lambda^g matrices associated with the variation and linearization 
    // of the discontinuous base vectors gd_i.

    // thickness coefficients
    double a = zeta ;
    double c = 1.0 ;

    // lambda1 =  ( I,   0,  a*I,   0 ,  0,  0 ,  0,  0 )
    // lambda2 =  ( 0,   I,   0 ,  a*I,  0,  0 ,  0,  0 )
    FloatMatrix eye(3,3), aEye(3,3);
    eye.beUnitMatrix();
    aEye=eye;  aEye.times(a);
    lambda[ 0 ].resize(3,18);   lambda[ 0 ].zero();  
    lambda[ 1 ].resize(3,18);   lambda[ 1 ].zero();
    lambda[ 0 ].setSubMatrix(eye,1,1);  lambda[ 0 ].setSubMatrix(aEye,1,7);
    lambda[ 1 ].setSubMatrix(eye,1,4);  lambda[ 1 ].setSubMatrix(aEye,1,10);

    // lambda3 =  ( 0,  0,  0 ,  0 , c*I , 0 , 0 , 0 )
    lambda[ 2 ].resize(3,18);   lambda[ 2 ].zero();
    lambda[ 2 ].at(1,13) = lambda[ 2 ].at(2,14) = lambda[ 2 ].at(3,15) = c;

}

void
Shell7BaseXFEM :: computeLambdaNMatrixDis(FloatMatrix &lambda_xd, double zeta)
{
    // computes the lambda_x matrix associated with the variation and linearization of the 
    // discontinuous position vector x_di. (\delta x_{di} = lambda_{xd}*\delta n_x)
    
    // lambda_x =  ( I, a*I, 0 )
    lambda_xd.resize(3,7);
    lambda_xd.zero();
    lambda_xd.at(1,1) = lambda_xd.at(2,2) = lambda_xd.at(3,3) = 1.0;
    lambda_xd.at(1,4) = lambda_xd.at(2,5) = lambda_xd.at(3,6) = zeta;
}



void
Shell7BaseXFEM :: computePressureTangentMatrixDis(FloatMatrix &KCC, FloatMatrix &KCD, FloatMatrix &KDD, IntegrationPoint *ip, Load *load, const int iSurf, TimeStep *tStep)
{
    // Computes tangent matrix associated with the linearization of pressure loading. Assumes constant pressure.
    ConstantPressureLoad* pLoad = dynamic_cast< ConstantPressureLoad * >( load );
    
    FloatMatrix N, B, NLB, L(7, 18), gcov, W1, W2;
    FloatArray lcoords(3), solVec, pressure;
    FloatArray g1, g2, genEps, genEpsJ, genEpsK;
    FloatMatrix lambdaGC [ 3 ], lambdaNC, lambdaGD [ 3 ], lambdaND;
    double xi   = pLoad->giveLoadOffset();
    double zeta = this->giveGlobalZcoord(xi);
    this->giveUpdatedSolutionVector(solVec, tStep);
    // compute w1,w2, KC
    int ndof = Shell7Base :: giveNumberOfDofs();


    lcoords.at(1) = ip->giveCoordinate(1);
    lcoords.at(2) = ip->giveCoordinate(2);
    lcoords.at(3) = xi;     // local coord where load is applied

    this->computeNmatrixAt(lcoords, N);
    this->computeBmatrixAt(lcoords, B);
    genEps.beProductOf(B, solVec);  

    
    FloatMatrix LCC(18,18), LCD(18,18), LDD(18,18); 
    LCC.zero(); LCD.zero(); LDD.zero();
        
    //(xc+xd)*(g1xg2)=xc*g1xg2 + xd*g1xg2 -> xc*(W2*Dg1 - W1*Dg2) + xd*(W2*Dg1 - W1*Dg2)
    // Traction tangent, L =  lambdaN * ( W2*lambdaG_1 - W1*lambdaG_2  ) 
    load->computeValueAt(pressure, tStep, * ( ip->giveCoordinates() ), VM_Total);        // pressure component   
    this->evalCovarBaseVectorsAt(lcoords, gcov, genEps);
    g1.beColumnOf(gcov,1);
    g2.beColumnOf(gcov,2);
    W1 = this->giveAxialMatrix(g1);
    W2 = this->giveAxialMatrix(g2);
        
    this->computeLambdaGMatrices(lambdaGC, genEps, zeta);
    this->computeLambdaNMatrix(lambdaNC, genEps, zeta);
    this->computeLambdaGMatricesDis(lambdaGD, zeta);
    this->computeLambdaNMatrixDis(lambdaND, zeta);

    FloatMatrix W2L, W1L;
    W2L.beProductOf(W2,lambdaGC[0]);
    W1L.beProductOf(W1,lambdaGC[1]);
    W2L.subtract(W1L);
    LCC.beTProductOf(lambdaNC, W2L);
    LCC.times( -pressure.at(1) );

    W2L.beProductOf(W2,lambdaGD[0]);
    W1L.beProductOf(W1,lambdaGD[1]);
    W2L.subtract(W1L);
    LCD.beTProductOf(lambdaNC, W2L);
    LCD.times( -pressure.at(1) );

    W2L.beProductOf(W2,lambdaGD[0]);
    W1L.beProductOf(W1,lambdaGD[1]);
    W2L.subtract(W1L);
    LDD.beTProductOf(lambdaND, W2L);
    LDD.times( -pressure.at(1) );


    FloatMatrix KCCtemp, KCDtemp, KDDtemp;
    this->computeTripleProduct(KCCtemp, N, LCC, B );    
    this->computeTripleProduct(KCDtemp, N, LCD, B );    
    this->computeTripleProduct(KDDtemp, N, LDD, B );    
    
    int ndofs = Shell7Base :: giveNumberOfDofs();
    KCC.resize(ndofs,ndofs); KCD.resize(ndofs,ndofs); KDD.resize(ndofs,ndofs); 
    KCC.zero(); KCD.zero(); KDD.zero();
    const IntArray &ordering = this->giveOrdering(All);
    KCC.assemble(KCCtemp, ordering, ordering);
    KCD.assemble(KCDtemp, ordering, ordering);
    KDD.assemble(KDDtemp, ordering, ordering);
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

        for ( int j = 1; j <= iRuleL->giveNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRuleL->getIntegrationPoint(j - 1);

            FloatMatrix N11, N22, N33, N;
            //this->computeNmatricesAt(gp, N11, N22, N33);
            this->computeNmatrixAt(gp, N);
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


            double zeta = giveGlobalZcoord(gp->giveCoordinate(3));
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

        
        //@todo assumes no variation i xi1-dir
        FloatArray componentsTemp, coordsTemp(1);
        coordsTemp.at(1) = 0.0; // 
        edgeLoad->computeValueAt(componentsTemp, tStep, coordsTemp, VM_Total);
        double xi = 0.0; // defaults to geometric midplane
        if ( componentsTemp.giveSize() == 8 ) {
            xi = componentsTemp.at(8);   // use the 8th component to store the-xi coord where the load acts
        }

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
        _error("Shell7BaseXFEM :: computeEdgeLoadVectorAt: load type not supported");
        return;
    }
}


// Surface
void
Shell7BaseXFEM :: computeSurfaceLoadVectorAt(FloatArray &answer, Load *load,
                                         int iSurf, TimeStep *tStep, ValueModeType mode)
{
    BoundaryLoad *surfLoad = dynamic_cast< BoundaryLoad * >( load );

    if ( surfLoad ) {
        answer.resize( this->giveNumberOfDofs() );
        answer.zero();

        // Continuous part
        FloatArray solVec, force;
        this->giveUpdatedSolutionVector(solVec, tStep);
        this->computePressureForce(force, solVec, iSurf, surfLoad, tStep);

        IntArray activeDofs, ordering, eiDofIdArray;
        this->computeOrderingArray(ordering, activeDofs, 0, All);
        answer.assemble(force, ordering);     
       
        // Disccontinuous part
#if 1
        FloatArray componentsTemp, solVecD;

        double xi = 0.0; // defaults to geometric midplane
        if ( ConstantPressureLoad* pLoad = dynamic_cast< ConstantPressureLoad * >( load ) ) {
            xi = pLoad->giveLoadOffset();
        }
        
        FloatArray temp;
        for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { // Only one is supported at the moment
            Delamination *dei =  dynamic_cast< Delamination * >( this->xMan->giveEnrichmentItem(i) ); // should check success

            for ( int j = 1; j <= dei->giveNumberOfEnrichmentDomains(); j++ ) {
                if ( dei->isElementEnrichedByEnrichmentDomain(this, j) ) {
                    double xi0 = dei->enrichmentDomainXiCoords.at(j);
                    if ( xi > xi0 ) {
                        dei->giveEIDofIdArray(eiDofIdArray, j);
                        this->giveSolutionVector(solVecD, eiDofIdArray, tStep);  
                        
                        this->computePressureForce(temp, solVecD, iSurf, surfLoad, tStep);
                        // Assemble
                        this->computeOrderingArray(ordering, activeDofs,  j, All);
                        FloatArray tempRed;
                        tempRed.beSubArrayOf(temp, activeDofs);
                        answer.assemble(tempRed, ordering);
                    }
                }
            }
        }
#endif
        return;
    } else {
        _error("Shell7Base :: computeSurfaceLoadVectorAt: load type not supported");
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

            for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
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

    #ifdef DEBUG
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
    //double zeta = localCoords.at(3)*this->layeredCS->giveLayerThickness(layer)*0.5 + this->layeredCS->giveLayerMidZ(layer);
    double zeta = this->giveGlobalZcoordInLayer(localCoords.at(3), layer); //@todo should probably be giveGlobalZetaCoord
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
