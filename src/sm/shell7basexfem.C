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

#include "shell7basexfem.h"
#include "shell7base.h"
#include "enrichmentitem.h"
#include "xfemmanager.h"
#include "dofmanager.h"
#include "constantpressureload.h"
#include "simpleinterfacemat.h"
#include "connectivitytable.h"
#include "bilinearczmaterialFagerstrom.h"
#include "mathfem.h"
#include "node.h"

#include "geometry.h"

namespace oofem {


const double DISC_DOF_SCALE_FAC = 1.0e-3;

Shell7BaseXFEM :: Shell7BaseXFEM(int n, Domain *aDomain) : Shell7Base(n, aDomain), XfemElementInterface(this) 
{
    czMat = NULL;
}

int 
Shell7BaseXFEM :: checkConsistency()
{
    Shell7Base :: checkConsistency();
    return 1;
}

void
Shell7BaseXFEM :: postInitialize()
{
    Shell7Base :: postInitialize();

    this->xMan =  this->giveDomain()->giveXfemManager();

}

void
Shell7BaseXFEM :: computeFailureCriteriaQuantities(FailureCriteriaStatus *fcStatus, TimeStep *tStep) 
{
    
    // Compute necessary quantities for evaluation of failure criterias
#if 1
    if ( DamagedNeighborLayeredStatus *status = dynamic_cast<DamagedNeighborLayeredStatus *>(fcStatus) ) {
        /*
        Go through the neighbors of the element and check for each layer if the 
        corresponding cz is damaged (if applicable)

        */
        IntArray neighbors;
        IntArray elements(1);
        ConnectivityTable *conTable = this->giveDomain()->giveConnectivityTable();
        elements.at(1) = fcStatus->el->giveNumber();
        conTable->giveElementNeighbourList(neighbors, elements); 
        FloatArray damageArray(this->layeredCS->giveNumberOfLayers() - 1), damageArrayNeigh;
        fcStatus->quantities.resize( neighbors.giveSize() );

        for ( int i = 1; i <= neighbors.giveSize(); i++ ) {

            Shell7BaseXFEM *neighbor = 
                dynamic_cast< Shell7BaseXFEM * > (this->giveDomain()->giveElement( neighbors.at(i) ));
            if ( neighbor ) {
                neighbor->giveMaxCZDamages(damageArrayNeigh, tStep); // damage parameter for each interface

                // store maximum damage from neighbors
                for ( int j = 1; j <= damageArray.giveSize(); j++ ) {
                    if (damageArrayNeigh.at(j) > damageArray.at(j) ) {
                        damageArray.at(j) = damageArrayNeigh.at(j);
                    }
                }
            }

        }


        // debugging
        for ( int j = 1; j <= damageArray.giveSize(); j++ ) {
            damageArray.at(j) = 1.0;                    
        }

        status->layerDamageValues = damageArray;
    }
    //case FC_DamagedNeighborCZ:
    //std::vector < FloatArray > interLamStresses;
    //int numInterfaces = this->layeredCS->giveNumberOfLayers() - 1;
    //IntArray delaminatedInterfaceList;
    //switch ( fc->giveType() ) {
    //case FC_MaxShearStress:
    //    
    //    this->computeDelaminatedInterfaceList(delaminatedInterfaceList); ///@todo remove the need for this method - JB
    //    
    //    fc->quantities.resize(numInterfaces); // will overwrite this every time
    //    //fc->elQuantities[ this->giveNumber() - 1 ].
    //    for (int intface = 1; intface <= numInterfaces; intface++ ) { 
    //        
    //        if ( ! delaminatedInterfaceList.contains(intface) ) { // interface is whole
    //            this->computeInterLaminarStressesAt(intface, tStep, interLamStresses); // all 6 components in each evaluation point (ip)
    //            int numEvalPoints = interLamStresses.size();
    //            fc->quantities[intface-1].resize(numEvalPoints); // most often = numIP
    //        
    //            for ( int evalPoint = 1; evalPoint <= numEvalPoints; evalPoint++) {
    //                FloatArray &values = fc->quantities[intface-1][evalPoint-1];  // one resulting shear stress
    //                FloatArray &vS = interLamStresses[evalPoint-1];               // Stress in eval point
    //                values.resize(1);                                             // scalar measure in this case
    //                values.at(1) = sqrt( vS.at(2)*vS.at(2) + vS.at(3)*vS.at(3) ); // components can't be right here? shouldn't it be a traction vector?
    //            }
    //        }
    //    }
    //    break;
    
#endif  
}





IRResultType Shell7BaseXFEM :: initializeFrom(InputRecord *ir)
{
    //const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    //IRResultType result;                   // Required by IR_GIVE_FIELD macro
    
    if ( ir->hasField(_IFT_Shell7BaseXFEM_CohesiveZoneMaterial) ) {
        OOFEM_ERROR("'czmaterial' this keyword is not in use anymore! Instead define cz material for each interface in the cross secton, ex: interfacematerials 3 x x x ");
    }

    Shell7Base :: initializeFrom(ir);
    return IRRT_OK; 
}

bool 
Shell7BaseXFEM :: hasCohesiveZone(int interfaceNum)
{
    return this->layeredCS->giveInterfaceMaterialNum(interfaceNum) > 0;
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


void
Shell7BaseXFEM :: giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const
{
    // Returns the total id mask of the dof manager - regular id's + enriched id's

    // Continuous part
    Shell7Base ::giveDofManDofIDMask(inode, ut, answer);
    XfemManager *xMan = this->giveDomain()->giveXfemManager(); // xman not initialized on el. level??

    // Discontinuous part
    DofManager *dMan = this->giveDofManager(inode);
    for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) { 
        EnrichmentItem *ei = xMan->giveEnrichmentItem(i);
        if ( ei->isDofManEnriched(*dMan) ) {
            IntArray eiDofIdArray;
            ei->giveEIDofIdArray(eiDofIdArray); 
            answer.followedBy(eiDofIdArray);
        }
     }
}


void
Shell7BaseXFEM :: evalCovarBaseVectorsAt(FloatArray &lCoords, FloatMatrix &gcov, FloatArray &genEpsC)
{
    // Continuous part
    Shell7Base :: evalCovarBaseVectorsAt(lCoords, gcov, genEpsC);

    // Discontinuous part
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    FloatArray dGenEps;
    FloatMatrix gcovd; 
    std :: vector< double > ef;
    for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { 
        EnrichmentItem *ei = this->xMan->giveEnrichmentItem(i);

        if ( ei->isElementEnriched(this) ) {
            computeDiscGeneralizedStrainVector(dGenEps, lCoords, ei, tStep); 
            Shell7Base :: evalCovarBaseVectorsAt(lCoords, gcovd, dGenEps);
            gcov.add(gcovd); 
        }
    }
}


void
Shell7BaseXFEM :: computeDiscGeneralizedStrainVector(FloatArray &answer, FloatArray &lCoords, EnrichmentItem *ei, TimeStep *tStep)
{
    FloatArray solVecD;
    IntArray eiDofIdArray;
    ei->giveEIDofIdArray(eiDofIdArray); 
    this->giveDisSolutionVector(solVecD, eiDofIdArray, tStep);  
    FloatMatrix B;
    this->computeEnrichedBmatrixAt(lCoords, B, ei);
    Shell7Base :: computeGeneralizedStrainVectorNew(answer, solVecD, B);

    //answer.times(DISC_DOF_SCALE_FAC);
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
            if ( ei->isDofManEnriched(*dMan) ) {
                IntArray eiDofIdArray;
                ei->giveEIDofIdArray(eiDofIdArray); 
                nDofs += eiDofIdArray.giveSize();
            }
        }
    }
    return nDofs;
}


void
Shell7BaseXFEM :: giveDisSolutionVector(FloatArray &answer, const IntArray &dofIdArray, TimeStep *tStep)
{
    // Scales the discontinuous dofs with a factor
    Shell7Base :: giveSolutionVector(answer, dofIdArray, tStep);
    //answer.times(DISC_DOF_SCALE_FAC);
}

void
Shell7BaseXFEM :: edgeGiveUpdatedSolutionVector(FloatArray &answer, const int iedge, TimeStep *tStep)
{
    Shell7Base :: edgeGiveUpdatedSolutionVector(answer, iedge, tStep);
    answer.times(DISC_DOF_SCALE_FAC);
}


void 
Shell7BaseXFEM :: computeOrderingArray( IntArray &orderingArray, IntArray &activeDofsArray,  EnrichmentItem *ei)
{
    // Routine to extract vector given an array of dofid items
    // If a certain dofId does not exist a zero is used as value

    const IntArray &ordering_cont = this->giveOrdering(All);
    IntArray fieldDofId; 
    Shell7Base::giveDofManDofIDMask(0, EID_MomentumBalance, fieldDofId);


    IntArray ordering_temp, activeDofsArrayTemp;
    ordering_temp.resize(ordering_cont.giveSize());
    activeDofsArrayTemp.resize(ordering_cont.giveSize());


    int activeDofPos = 0, activeDofIndex = 0, orderingDofIndex = 0;
    
    IntArray dofManDofIdMask, dofManDofIdMaskAll;
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        DofManager *dMan = this->giveDofManager(i);
        

        if ( ei == NULL) { // return mask corresponding to the regular id's
           Shell7Base ::giveDofManDofIDMask(i, EID_MomentumBalance, dofManDofIdMask);
        } else {
            if ( ei->isDofManEnriched(*dMan) ) {
                ei->giveEIDofIdArray(dofManDofIdMask); 
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

        dofManDofIdMask.resize(0);
    }
    
    // Reduce arrays to actual size 
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
    this->computeOrderingArray(ordering, activeDofs, NULL); 
    answer.assemble(temp, ordering);

    // Disccontinuous part
    IntArray eiDofIdArray;
    FloatArray solVecD, tempRed, fCZ;
    for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { 
        
        EnrichmentItem *ei = this->xMan->giveEnrichmentItem(i);

        if ( ei->isElementEnriched(this) ) {

            ei->giveEIDofIdArray(eiDofIdArray); 
            this->giveDisSolutionVector(solVecD, eiDofIdArray, tStep);  
            this->discComputeSectionalForces(temp, tStep, solVec, solVecD, ei);

            this->computeOrderingArray(ordering, activeDofs, ei);
            tempRed.beSubArrayOf(temp, activeDofs);
            answer.assemble(tempRed, ordering);

            // Cohesive zone model
            if ( this->hasCohesiveZone(i) ) {
                Delamination *dei =  dynamic_cast< Delamination * >( this->xMan->giveEnrichmentItem(i) ); // should check success
                this->computeCohesiveForces( fCZ, tStep, solVec, solVecD, dei); 
                tempRed.beSubArrayOf(fCZ, activeDofs);
                answer.assemble(tempRed, ordering);
            }
        }
    }

}


void
Shell7BaseXFEM :: discComputeSectionalForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, FloatArray &solVecD, EnrichmentItem *ei)
//
{
    int ndofs = Shell7Base :: giveNumberOfDofs();
    int numberOfLayers = this->layeredCS->giveNumberOfLayers(); 
    FloatArray f(ndofs), genEps, genEpsD, ftemp, lCoords, sectionalForces, N;
    FloatMatrix B, BEnr;
    f.zero();

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRuleL = integrationRulesArray [ layer - 1 ];
        Material *mat = domain->giveMaterial( this->layeredCS->giveLayerMaterial(layer) );

        for ( int j = 1; j <= iRuleL->giveNumberOfIntegrationPoints(); j++ ) {
            GaussPoint *gp = iRuleL->getIntegrationPoint(j - 1);
            lCoords = *gp->giveCoordinates();           

            this->computeEnrichedBmatrixAt(lCoords, B, NULL);
            this->computeEnrichedBmatrixAt(lCoords, BEnr, ei);
            this->computeGeneralizedStrainVectorNew(genEps,  solVec,  B);
            this->computeGeneralizedStrainVectorNew(genEpsD, solVecD, BEnr);

            //double zeta = giveGlobalZcoord(gp->giveCoordinate(3));
            double zeta = giveGlobalZcoord( gp->giveCoordinate( 3 ), *gp->giveCoordinates( ) );
            this->computeSectionalForcesAt(sectionalForces, gp, mat, tStep, genEps, genEpsD, zeta);

            // Computation of nodal forces: f = B^t*[N M T Ms Ts]^t
            ftemp.beTProductOf( BEnr, sectionalForces );
            double dV = this->computeVolumeAroundLayer(gp, layer);
            f.add( dV, ftemp );            


        }
    }

    answer.resize(ndofs); answer.zero();
    const IntArray &ordering_all = this->giveOrdering(All);
    answer.assemble(f, ordering_all);

}


double 
Shell7BaseXFEM :: evaluateLevelSet(const FloatArray &lCoords, EnrichmentItem *ei)
{
    // Evaluates the corresponding level set depending on the type of enrichment
    double levelSet = 0.0;
    if ( dynamic_cast< Delamination * >( ei ) ) {
        levelSet = lCoords.at(3) - dynamic_cast< Delamination * >( ei )->giveDelamXiCoord();                   
    } else if ( dynamic_cast< Crack * >( ei ) ) {
        FloatArray N;
        const IntArray &elNodes = this->giveDofManArray();
        this->fei->evalN( N, lCoords, FEIElementGeometryWrapper(this) );
        ei->interpLevelSet(levelSet, N, elNodes);
    } else {
        OOFEM_ERROR1("Shell7BaseXFEm error in evaluation of levelset");
    }          
    return levelSet;
}

double 
Shell7BaseXFEM :: edgeEvaluateLevelSet(const FloatArray &lCoords, EnrichmentItem *ei)
{
    // Evaluates the corresponding level set depending on the type of enrichment
    double levelSet = 0.0;
    if ( dynamic_cast< Delamination * >( ei ) ) {
        double xiLoad = 0.0; ///@todo need info about load position
        levelSet = xiLoad - dynamic_cast< Delamination * >( ei )->giveDelamXiCoord();                   
    } else if ( dynamic_cast< Crack * >( ei ) ) {
        FloatArray N;
        const IntArray &elNodes = this->giveDofManArray();
        this->fei->edgeEvalN( N, 1, lCoords, FEIElementGeometryWrapper(this) );
        ei->interpLevelSet(levelSet, N, elNodes);
    } else {
        OOFEM_ERROR1("Shell7BaseXFEm error in evaluation of levelset");
    }          
    return levelSet;
}

void
Shell7BaseXFEM :: giveMaxCZDamages(FloatArray &answer, TimeStep *tStep)
{
    // for each interface get maximum interface damage of all the ip's
    int numZones = this->layeredCS->giveNumberOfLayers() - 1;
    answer.resize(numZones);
    StructuralInterfaceMaterial *mat;
    FloatArray ipValues;
    for ( int i = 0; i < numZones; i++ ) {
        if ( hasCohesiveZone(i+1) ) {
            IntegrationRule *iRule = czIntegrationRulesArray [ i ];
            double max = 0.0;    
            for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
                IntegrationPoint *ip = iRule->getIntegrationPoint(j);

                mat = static_cast < StructuralInterfaceMaterial * > (this->layeredCS->giveInterfaceMaterial(i+1) );
                //mat = static_cast < StructuralInterfaceMaterial * > ( this->czMat );
                mat->giveIPValue(ipValues, ip, IST_DamageScalar, tStep);

                double val = ipValues.at(1);
                if ( val > max ) {
                    max = val;
                }

            }
                answer(i) = max;
        } else {
            answer(i) = 0.0;
        }
    }

    

}


void
Shell7BaseXFEM :: computeCohesiveForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVecC, FloatArray &solVecD, Delamination *dei)
{
    //Computes the cohesive nodal forces for a given interface
    FloatArray answerTemp, lCoords(3);
    answerTemp.resize(Shell7Base :: giveNumberOfDofs() ); 
    answerTemp.zero();

    FloatMatrix N, B, F;  
    int delamNum = dei->giveNumber();
    IntegrationRule *iRuleL = czIntegrationRulesArray [ delamNum - 1 ]; ///@ todo does this work with giveNumber?

    StructuralInterfaceMaterial *intMat = static_cast < StructuralInterfaceMaterial * > 
        (this->layeredCS->giveInterfaceMaterial(delamNum) );

    FloatMatrix lambda, lambdaN, Q;
    FloatArray Fp, T, interfaceXiCoords, nCov, xd, unknowns, genEpsC, genEpsD;
    this->layeredCS->giveInterfaceXiCoords(interfaceXiCoords);

    for ( int i = 1; i <= iRuleL->giveNumberOfIntegrationPoints(); i++ ) {
        IntegrationPoint *ip = iRuleL->getIntegrationPoint(i - 1);
        lCoords.at(1) = ip->giveCoordinate(1);
        lCoords.at(2) = ip->giveCoordinate(2);
        lCoords.at(3) = dei->giveDelamXiCoord();

        this->computeBmatrixAt(lCoords, B);
        this->computeNmatrixAt(lCoords, N);

        // Lambda matrix
        genEpsD.beProductOf(B, solVecD);
        double xi = dei->giveDelamXiCoord();
        double zeta = xi * this->layeredCS->computeIntegralThick() * 0.5;
        this->computeLambdaNMatrix(lambda, genEpsD, zeta);
        
        // Compute jump vector
        unknowns.beProductOf(N, solVecD);
        xd.beProductOf(lambda,unknowns); // spatial jump
		genEpsC.beProductOf(B, solVecC);
		this->computeFAt(lCoords, F, genEpsC);

		// Transform xd and F to a local coord system
        this->evalInitialCovarNormalAt(nCov, lCoords);
        Q.beLocalCoordSys(nCov);
        xd.rotatedWith(Q,'n');
        F.rotatedWith(Q,'n');

        // Compute cohesive traction based on jump
        intMat->giveFirstPKTraction_3d(T, ip, xd, F, tStep);
        lambdaN.beProductOf(lambda,N);
        T.rotatedWith(Q,'t'); // transform back to global coord system

        Fp.beTProductOf(lambdaN, T);
        double dA = this->computeAreaAround(ip,xi);
        answerTemp.add(dA*DISC_DOF_SCALE_FAC,Fp);
    }

    int ndofs = Shell7Base :: giveNumberOfDofs();
    answer.resize(ndofs); answer.zero();
    const IntArray &ordering = this->giveOrdering(All);
    answer.assemble(answerTemp, ordering);
}

void
Shell7BaseXFEM :: updateYourself(TimeStep *tStep)
// Updates the receiver at end of step.
{
    Shell7Base :: updateYourself(tStep);

    for ( int i = 0; i < this->layeredCS->giveNumberOfLayers() - 1; i++ ) {
        if ( this->hasCohesiveZone(i+1) ) {
            czIntegrationRulesArray [ i ]->updateYourself(tStep);
        }
    }
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

    // Disccontinuous part (continuous part does not contribute)
    IntArray eiDofIdArray, orderingJ, activeDofsJ;
    FloatMatrix tempRed;
    for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) {

        Delamination *dei =  dynamic_cast< Delamination * >( this->xMan->giveEnrichmentItem(i) );
        if ( dei != NULL && dei->isElementEnriched(this) && this->hasCohesiveZone(i) ) {
            dei->giveEIDofIdArray(eiDofIdArray); 
            this->giveDisSolutionVector(solVecD, eiDofIdArray, tStep);  
            this->computeCohesiveTangentAt(temp, tStep, solVecD, dei);
            // Assemble part correpsonding to active dofs
            this->computeOrderingArray(orderingJ, activeDofsJ, dei); //
            tempRed.beSubMatrixOf(temp, activeDofsJ, activeDofsJ);
            answer.assemble(tempRed, orderingJ, orderingJ);
        }
    }
}



void
Shell7BaseXFEM :: computeCohesiveTangentAt(FloatMatrix &answer, TimeStep *tStep, FloatArray &solVecD, 
    Delamination *dei)
{
    //Computes the cohesive tangent for a given interface
    FloatArray lCoords(3);
    FloatMatrix answerTemp, N, lambda, K, temp, tangent;
    int nDofs = Shell7Base :: giveNumberOfDofs();
    answerTemp.resize(nDofs, nDofs); 
    answerTemp.zero();

    int delamNum = dei->giveNumber();
    IntegrationRule *iRuleL = czIntegrationRulesArray [ dei->giveNumber() - 1 ];
    StructuralInterfaceMaterial *intMat = static_cast < StructuralInterfaceMaterial * > 
        (this->layeredCS->giveInterfaceMaterial(delamNum) );

    double xi = dei->giveDelamXiCoord();
    //double zeta = this->giveGlobalZcoord(xi);
    //this->computeLambdaNMatrixDis(lambda, zeta);

    FloatMatrix Q;
    FloatArray nCov;
    FloatArray interfaceXiCoords;
    this->layeredCS->giveInterfaceXiCoords(interfaceXiCoords);
    for ( int i = 1; i <= iRuleL->giveNumberOfIntegrationPoints(); i++ ) {
        IntegrationPoint *ip = iRuleL->getIntegrationPoint(i - 1);
        lCoords.at(1) = ip->giveCoordinate(1);
        lCoords.at(2) = ip->giveCoordinate(2);
        lCoords.at(3) = dei->giveDelamXiCoord();
        double zeta = giveGlobalZcoord( xi, lCoords);
        this->computeLambdaNMatrixDis( lambda, zeta );
        this->computeNmatrixAt(lCoords, N);
                
        intMat->give3dStiffnessMatrix_dTdj(K, TangentStiffness, ip, tStep);
        this->evalInitialCovarNormalAt(nCov, lCoords);
        Q.beLocalCoordSys(nCov);
        K.rotatedWith(Q,'t');   // rotate back to global coord system

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



void 
Shell7BaseXFEM :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{

    int ndofs = this->giveNumberOfDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();
    
    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     
    FloatMatrix tempRed, tempRedT;
    FloatMatrix KCC, KCD, KDD, Bc;
    IntArray orderingC, activeDofsC;
    this->computeOrderingArray(orderingC, activeDofsC, NULL); //
    std::vector<IntArray> orderingArrays;
    std::vector<IntArray> activeDofsArrays;

    FloatMatrix A [ 3 ] [ 3 ];
    FloatArray lCoords, N; FloatArray genEpsC, solVecC;
    this->giveUpdatedSolutionVector(solVecC, tStep);

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRule = integrationRulesArray [ layer - 1 ];
        StructuralMaterial *material = static_cast< StructuralMaterial* >( domain->giveMaterial( this->layeredCS->giveLayerMaterial(layer) ) );

        for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
            IntegrationPoint *ip = iRule->getIntegrationPoint(i);
            lCoords = *ip->giveCoordinates();
            this->computeEnrichedBmatrixAt(lCoords, Bc, NULL);
            this->computeGeneralizedStrainVectorNew(genEpsC, solVecC , Bc);
            
            Shell7Base :: computeLinearizedStiffness(ip, material, tStep, A, genEpsC);
            
            // Continuous part K_{c,c}
            this->discComputeBulkTangentMatrix(KCC, ip, NULL, NULL, layer, A, tStep);
            answer.assemble(KCC, orderingC, orderingC);

            // Discontinuous part
            int numEI = this->xMan->giveNumberOfEnrichmentItems();
            for ( int m = 1; m <= numEI; m++ ) { 
                EnrichmentItem *eiM = this->xMan->giveEnrichmentItem(m);
                if ( orderingArrays.size() == 0 ) {
                    orderingArrays.resize(numEI);
                    activeDofsArrays.resize(numEI);
                }

                if ( eiM->isElementEnriched(this) ) {
                    
                    if( orderingArrays[m-1].giveSize() == 0 ) {
                        this->computeOrderingArray( orderingArrays[m-1], activeDofsArrays[m-1], eiM); 
                    }

                    this->discComputeBulkTangentMatrix(KCD, ip, NULL, eiM, layer, A, tStep);
                    tempRed.beSubMatrixOf(KCD, activeDofsC, activeDofsArrays[m-1]);
                    answer.assemble(tempRed, orderingC, orderingArrays[m-1]);
                    tempRedT.beTranspositionOf(tempRed);
                    answer.assemble(tempRedT, orderingArrays[m-1], orderingC);

                    // K_{dk,dl}
                    for ( int k = 1; k <= numEI; k++ ) {
                        EnrichmentItem *eiK = this->xMan->giveEnrichmentItem(k);
                        
                        if ( eiK->isElementEnriched(this) ) {
                            this->discComputeBulkTangentMatrix(KDD, ip, eiM, eiK, layer, A, tStep);
                            tempRed.beSubMatrixOf(KDD, activeDofsArrays[m-1], activeDofsArrays[k-1]);
                            answer.assemble(tempRed, orderingArrays[m-1], orderingArrays[k-1]);
                        }

                    }
                }
            }

        }
    }


    // Cohesive zones
#if 0
    FloatMatrix Kcoh;
    //if ( this->hasCohesiveZone() ) {
        this->computeCohesiveTangent(Kcoh, tStep);
        Kcoh.times(DISC_DOF_SCALE_FAC*DISC_DOF_SCALE_FAC);
        answer.add(Kcoh);
    //}
#endif


    // Add contribution due to pressure load
#if 0

    int nLoads = this->boundaryLoadArray.giveSize() / 2;

    for ( int k = 1; k <= nLoads; k++ ) {     // For each pressure load that is applied
        int load_number = this->boundaryLoadArray.at(2 * k - 1);
        int iSurf = this->boundaryLoadArray.at(2 * k);         // load_id
        Load *load = this->domain->giveLoad(load_number);
        std :: vector< double > efM, efK;

        if ( ConstantPressureLoad* pLoad = dynamic_cast< ConstantPressureLoad * >( load ) ) {

            IntegrationRule *iRule = specialIntegrationRulesArray [ 1 ];

            for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
                IntegrationPoint *ip = iRule->getIntegrationPoint(i);
                this->computePressureTangentMatrixDis(KCC, KCD, KDD, ip, load, iSurf, tStep);
                KCD.times(DISC_DOF_SCALE_FAC);
                KDD.times(DISC_DOF_SCALE_FAC*DISC_DOF_SCALE_FAC);

                // Continuous part
                answer.assemble(KCC, orderingC, orderingC);

                // Discontinuous part
                int numEI = this->xMan->giveNumberOfEnrichmentItems();
                for ( int m = 1; m <= numEI; m++ ) { 
                    Delamination *deiM = dynamic_cast< Delamination * >( this->xMan->giveEnrichmentItem(m) );

                    if ( deiM !=NULL && deiM->isElementEnriched(this) ) {
                        double levelSetM = pLoad->giveLoadOffset() - deiM->giveDelamXiCoord();
                        deiM->evaluateEnrFuncAt(efM, lCoords, levelSetM); 
                        
                        IntArray &orderingJ = orderingArrays[m-1];
                        IntArray &activeDofsJ = activeDofsArrays[m-1];

                        // con-dis & dis-con
                        if ( efM[0] > 0.1 ) {                            
                            tempRed.beSubMatrixOf(KCD, activeDofsC, activeDofsJ);
                            answer.assemble(tempRed, orderingC, orderingJ);
                            tempRedT.beTranspositionOf(tempRed);
                            answer.assemble(tempRedT, orderingJ, orderingC);
                        }

                        // dis-dis
                        for ( int k = 1; k <= numEI; k++ ) {
                            Delamination *deiK = dynamic_cast< Delamination * >( this->xMan->giveEnrichmentItem(k) );
                            if ( deiK != NULL && deiK->isElementEnriched(this) ) {
                                double levelSetK = pLoad->giveLoadOffset() - deiK->giveDelamXiCoord();
                                deiK->evaluateEnrFuncAt(efK, lCoords, levelSetK);   
                                if ( efM[0] > 0.1 && efK[0] > 0.1 ) {
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
#endif
    
    /*FloatMatrix test;
    test= answer;
    
    NEW_computeStiffnessMatrix(answer, rMode, tStep);
    test.subtract(answer);
    test.printYourself();*/
}


void
Shell7BaseXFEM :: discComputeBulkTangentMatrix(FloatMatrix &KdIJ, IntegrationPoint *ip, EnrichmentItem *ei1, EnrichmentItem *ei2, int layer, FloatMatrix A [ 3 ] [ 3 ],  TimeStep *tStep)
{

    FloatMatrix temp, B1, B2;
    FloatMatrix lambda1 [ 3 ], lambda2 [ 3 ];

    FloatArray lCoords;
    lCoords = *ip->giveCoordinates();
    this->computeEnrichedBmatrixAt(lCoords, B1, ei1);
    this->computeEnrichedBmatrixAt(lCoords, B2, ei2);

    //double zeta = giveGlobalZcoord(ip->giveCoordinate(3));
    double zeta = giveGlobalZcoord( ip->giveCoordinate( 3 ), *ip->giveCoordinates( ) );
    if ( ei1 ) {
        this->computeLambdaGMatricesDis(lambda1, zeta);
    } else {
        FloatArray solVecC, genEpsC;
        this->giveUpdatedSolutionVector(solVecC, tStep);  
        this->computeGeneralizedStrainVectorNew(genEpsC, solVecC , B1);
        this->computeLambdaGMatrices(lambda1, genEpsC, zeta);
    }
    
    if ( ei2 ) {
        this->computeLambdaGMatricesDis(lambda2, zeta);
    } else {
        FloatArray solVecC, genEpsC;
        this->giveUpdatedSolutionVector(solVecC, tStep);  
        this->computeGeneralizedStrainVectorNew(genEpsC, solVecC , B2);
        this->computeLambdaGMatrices(lambda2, genEpsC, zeta);
    }

    double dV = this->computeVolumeAroundLayer(ip, layer);

    FloatMatrix LDAD(18,18); 
    LDAD.zero();
    for ( int j = 0; j < 3; j++ ) {
        for ( int k = 0; k < 3; k++ ) {
            this->computeTripleProduct(temp, lambda1 [ j ], A [ j ][ k ], lambda2 [ k ] );
            LDAD.add(dV,temp);
        }
    }

    FloatMatrix KDDtemp;
    this->computeTripleProduct(KDDtemp, B1, LDAD, B2 );    
    
    int ndofs = Shell7Base :: giveNumberOfDofs();
    KdIJ.resize(ndofs,ndofs);
    KdIJ.zero();
    const IntArray &ordering = this->giveOrdering(All);
    
    KdIJ.assemble(KDDtemp, ordering, ordering);

}


void
Shell7BaseXFEM :: discComputeBulkTangentMatrix(FloatMatrix &KCC, FloatMatrix &KCD, FloatMatrix &KDD, IntegrationPoint *ip, Material *mat, int layer, TimeStep *tStep)
{
// old remove
    OOFEM_ERROR1("discComputeBulkTangentMatrix should not be used!!")
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
    
    FloatMatrix N, B, L(7, 18), gcov, W1, W2;
    FloatArray lcoords(3), solVec, pressure;
    FloatArray g1, g2, genEps;
    FloatMatrix lambdaGC [ 3 ], lambdaNC, lambdaGD [ 3 ], lambdaND;
    double xi   = pLoad->giveLoadOffset();
    //double zeta = this->giveGlobalZcoord(xi);
    this->giveUpdatedSolutionVector(solVec, tStep);
    // compute w1,w2, KC
    lcoords.at(1) = ip->giveCoordinate(1);
    lcoords.at(2) = ip->giveCoordinate(2);
    lcoords.at(3) = xi;     // local coord where load is applied

    double zeta = giveGlobalZcoord( xi, lcoords );

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
    ///@todo broken!

#if 0
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
#endif
}



void
Shell7BaseXFEM :: computeEdgeLoadVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *tStep, ValueModeType mode)
{
    BoundaryLoad *edgeLoad = dynamic_cast< BoundaryLoad * >( load );
    if ( edgeLoad ) {
        answer.resize( this->computeNumberOfDofs() );
        answer.zero();

        // Continuous part
        FloatArray fT;
        Shell7Base :: computeTractionForce(fT, iEdge, edgeLoad, tStep, mode);

        IntArray activeDofs, ordering; 
        this->computeOrderingArray(ordering, activeDofs, NULL); 
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
        FloatArray temp, tempRed;
        for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { // Only one is supported at the moment

            EnrichmentItem *ei = this->xMan->giveEnrichmentItem(i);
            if ( ei->isElementEnriched(this) ) {
                this->computeTractionForce(temp, iEdge, edgeLoad, tStep, mode, ei);
                this->computeOrderingArray(ordering, activeDofs, ei); 
                tempRed.beSubArrayOf(temp, activeDofs);
                answer.assemble(tempRed, ordering);                    
            }
        }
        return;
    } else {
        _error("Shell7BaseXFEM :: computeEdgeLoadVectorAt: load type not supported");
        return;
    }
}


void
Shell7BaseXFEM :: computeTractionForce(FloatArray &answer, const int iEdge, BoundaryLoad *edgeLoad, TimeStep *tStep, ValueModeType mode, EnrichmentItem *ei)
{
    // 
    IntegrationRule *iRule = specialIntegrationRulesArray [ 2 ];   // rule #3 for edge integration of distributed loads given in [*/m]
    GaussPoint *gp;

    FloatMatrix N, Q;
    FloatArray fT(7), components, lCoords;
    
    //BoundaryLoad :: CoordSystType coordSystType = edgeLoad->giveCoordSystMode();
    Load :: CoordSystType coordSystType = edgeLoad->giveCoordSystMode();

    FloatArray Nftemp(21), Nf(21);
    Nf.zero();
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        lCoords =*gp->giveCoordinates();

        edgeLoad->computeValueAt(components, tStep, lCoords, mode);
        this->edgeComputeEnrichedNmatrixAt(lCoords, N, ei);

        //if ( coordSystType ==  BoundaryLoad :: BL_UpdatedGlobalMode ) {
        if ( coordSystType ==  Load :: CST_UpdatedGlobal ) {
            
            // Updated global coord system
            FloatMatrix gcov;
            FloatArray lCoords = *gp->giveCoordinates();
            this->edgeEvalEnrCovarBaseVectorsAt(lCoords, iEdge, gcov, tStep, ei); 
            Q.beTranspositionOf(gcov);

            FloatArray distrForces(3), distrMoments(3), t1, t2;
            distrForces .setValues(3, components.at(1), components.at(2), components.at(3) );
            distrMoments.setValues(3, components.at(4), components.at(5), components.at(6) );
            t1.beTProductOf(Q, distrForces);
            t2.beTProductOf(Q, distrMoments);
            fT.addSubVector(t1,1);
            fT.addSubVector(t2,4);
            fT.at(7) = components.at(7); // don't do anything with the 'gamma'-load

        //} else if( coordSystType == BoundaryLoad :: BL_GlobalMode ) { 
        } else if( coordSystType == Load :: CST_Global ) { 
            // Undeformed global coord system
            for ( int i = 1; i <= 7; i++) {
                fT.at(i) = components.at(i);
            }
        } else {
            OOFEM_ERROR("Shell7Base :: computeTractionForce - does not support local coordinate system");
        }

        double dL = this->edgeComputeLengthAround(gp, iEdge);        
        
        Nftemp.beTProductOf(N, fT*dL);
        Nf.add(Nftemp);
    }

    IntArray mask;
    this->giveEdgeDofMapping(mask, iEdge);
    answer.resize( Shell7Base :: giveNumberOfDofs()  );
    answer.zero();
    answer.assemble(Nf, mask);

}

void
Shell7BaseXFEM :: edgeEvalEnrCovarBaseVectorsAt(FloatArray &lcoords, const int iedge, FloatMatrix &gcov, TimeStep *tStep, EnrichmentItem *ei)
{
    // Evaluates the covariant base vectors in the current configuration for an edge
    double zeta = lcoords.at(3);

    FloatArray solVecEdge;
    FloatMatrix B;
    IntArray edgeNodes;
    this->fei->computeLocalEdgeMapping(edgeNodes, iedge);
    this->edgeComputeEnrichedBmatrixAt(lcoords, B, ei);
    this->edgeGiveUpdatedSolutionVector(solVecEdge, iedge, tStep);

    FloatArray genEpsEdge;                 // generalized strain
    genEpsEdge.beProductOf(B, solVecEdge); // [dxdxi, dmdxi, m, dgamdxi, gam]^T

    FloatArray dxdxi, m, dmdxi;
    dxdxi.setValues( 3, genEpsEdge.at(1), genEpsEdge.at(2), genEpsEdge.at(3) );
    dmdxi.setValues( 3, genEpsEdge.at(4), genEpsEdge.at(5), genEpsEdge.at(6) );
        m.setValues( 3, genEpsEdge.at(7), genEpsEdge.at(8), genEpsEdge.at(9) );
    double dgamdxi = genEpsEdge.at(10);
    double gam     = genEpsEdge.at(11);

    double fac1 = ( zeta + 0.5 * gam * zeta * zeta );
    double fac2 = ( 0.5 * zeta * zeta );
    double fac3 = ( 1.0 + zeta * gam );
    
    FloatArray g1, g2, g3;
    g2 = dxdxi + fac1*dmdxi + fac2*dgamdxi*m; // base vector along the edge
    g3 = fac3*m;                              // director field

    g2.normalize();
    g3.normalize();
    g1.beVectorProductOf(g2, g3);
    g1.normalize();
    gcov.resize(3,3);
    gcov.setColumn(g1,1);
    gcov.setColumn(g2,2);
    gcov.setColumn(g3,3);
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
        this->computePressureForce(force, solVec, iSurf, surfLoad, tStep, mode);

        IntArray activeDofs, ordering, eiDofIdArray;
        this->computeOrderingArray(ordering, activeDofs, NULL); 
        answer.assemble(force, ordering);     
       
        // Disccontinuous part
#if 1
        FloatArray solVecD;
        double xi = 0.0; // defaults to geometric midplane
        if ( ConstantPressureLoad* pLoad = dynamic_cast< ConstantPressureLoad * >( load ) ) {
            xi = pLoad->giveLoadOffset();
        }
        std :: vector< double > ef;
        FloatArray temp, tempRed, lCoords(2);
        lCoords.zero();
        for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) {
            Delamination *dei =  dynamic_cast< Delamination * >( this->xMan->giveEnrichmentItem(i) ); 
            if ( dei != NULL && dei->isElementEnriched(this) ) {
                double levelSet = xi - dei->giveDelamXiCoord();
                dei->evaluateEnrFuncAt(ef, lCoords, levelSet); 
                if ( ef[0] > 0.1 ) {
                   dei->giveEIDofIdArray(eiDofIdArray); 
                   this->giveDisSolutionVector(solVecD, eiDofIdArray, tStep);  
                   this->computePressureForce(temp, solVecD, iSurf, surfLoad, tStep, mode);
                   temp.times(DISC_DOF_SCALE_FAC);
                   // Assemble
                   this->computeOrderingArray(ordering, activeDofs, dei);
                   tempRed.beSubArrayOf(temp, activeDofs);
                   answer.assemble(tempRed, ordering);
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



// Shifted N and B matrices
void
Shell7BaseXFEM :: computeEnrichedBmatrixAt(FloatArray &lcoords, FloatMatrix &answer, EnrichmentItem *ei)
{
    // Returns the enriched and shifted {B} matrix of the receiver, evaluated at gp. Such that
    // B*a = [dxbar_dxi, dwdxi, w, dgamdxi, gam]^T, where a is the vector of unknowns
 
    if ( ei && dynamic_cast< Crack*>(ei) ) { 

        int ndofs = Shell7Base :: giveNumberOfDofs();
        int ndofs_xm  = this->giveNumberOfFieldDofs(Midplane);
        answer.resize(18, ndofs);
        answer.zero();
        FloatArray N;
        FloatMatrix dNdxi;
        this->fei->evalN( N, lcoords, FEIElementGeometryWrapper(this) );
        this->fei->evaldNdxi( dNdxi, lcoords, FEIElementGeometryWrapper(this) );

        /*    18   18   6
         * 6 [B_u   0   0
         * 6   0   B_w  0
         * 3   0   N_w  0
         * 2   0    0  B_gam
         * 1   0    0  N_gam]
         */

        // Evaluate enrichment function at point given by lcoords
        std :: vector< double > efGP;
        double levelSetGP = this->evaluateLevelSet(lcoords, ei);
        ei->evaluateEnrFuncAt(efGP, lcoords, levelSetGP);

        int ndofman = this->giveNumberOfDofManagers();

        // First column
        for ( int i = 1, j = 0; i <= ndofman; i++, j += 3 ) {
            double factor = efGP [ 0 ] - EvaluateEnrFuncInDofMan(i, ei);

            answer.at(1, 1 + j) = dNdxi.at(i, 1) * factor;
            answer.at(2, 2 + j) = dNdxi.at(i, 1) * factor;
            answer.at(3, 3 + j) = dNdxi.at(i, 1) * factor;
            answer.at(4, 1 + j) = dNdxi.at(i, 2) * factor;
            answer.at(5, 2 + j) = dNdxi.at(i, 2) * factor;
            answer.at(6, 3 + j) = dNdxi.at(i, 2) * factor;
        }

        // Second column
        for ( int i = 1, j = 0; i <= ndofman; i++, j += 3 ) {
            double factor = efGP [ 0 ] - EvaluateEnrFuncInDofMan(i, ei);

            answer.at(7, ndofs_xm + 1 + j) = dNdxi.at(i, 1) * factor;
            answer.at(8, ndofs_xm + 2 + j) = dNdxi.at(i, 1) * factor;
            answer.at(9, ndofs_xm + 3 + j) = dNdxi.at(i, 1) * factor;
            answer.at(10, ndofs_xm + 1 + j) = dNdxi.at(i, 2) * factor;
            answer.at(11, ndofs_xm + 2 + j) = dNdxi.at(i, 2) * factor;
            answer.at(12, ndofs_xm + 3 + j) = dNdxi.at(i, 2) * factor;
            answer.at(13, ndofs_xm + 1 + j) = N.at(i) * factor;
            answer.at(14, ndofs_xm + 2 + j) = N.at(i) * factor;
            answer.at(15, ndofs_xm + 3 + j) = N.at(i) * factor;
        }

        // Third column
        for ( int i = 1, j = 0; i <= ndofman; i++, j += 1 ) {
            double factor = efGP [ 0 ] - EvaluateEnrFuncInDofMan(i, ei);

            answer.at(16, ndofs_xm * 2 + 1 + j) = dNdxi.at(i, 1) * factor;
            answer.at(17, ndofs_xm * 2 + 1 + j) = dNdxi.at(i, 2) * factor;
            answer.at(18, ndofs_xm * 2 + 1 + j) = N.at(i) * factor;
        }

        answer.times(DISC_DOF_SCALE_FAC);

    } else if ( ei && dynamic_cast< Delamination*>(ei) ){
        Shell7Base :: computeBmatrixAt(lcoords, answer);
        std :: vector< double > efGP;
        double levelSetGP = this->evaluateLevelSet(lcoords, ei);
        ei->evaluateEnrFuncAt(efGP, lcoords, levelSetGP);
        if ( efGP[0] > 0.1 ) {
            answer.times( efGP[0]*DISC_DOF_SCALE_FAC );
        } else {
            answer.times(0.0);
        }
    } else {
        Shell7Base :: computeBmatrixAt(lcoords, answer);
    }
}

double
Shell7BaseXFEM :: EvaluateEnrFuncInDofMan(int dofManNum, EnrichmentItem *ei)
{
        DofManager *dMan = this->giveDofManager(dofManNum);
        int globalNodeInd = dMan->giveGlobalNumber(); // global number in order to pick levelset value in that node
        double levelSetNode  = 0.0;
        ei->evalLevelSetNormalInNode( levelSetNode, globalNodeInd );
        std :: vector< double >efNode;
        const FloatArray &nodePos = * ( dMan->giveCoordinates() );
        ei->evaluateEnrFuncAt(efNode, nodePos, levelSetNode, globalNodeInd);
        if( efNode.size() ) {
            return efNode [ 0 ];
        } else {
            return 0.0;
        }
        //return efNode [ 0 ];
        //return 0.;
}


void
Shell7BaseXFEM :: computeEnrichedNmatrixAt(const FloatArray &lcoords, FloatMatrix &answer, EnrichmentItem *ei)
{
    // Returns the displacement interpolation matrix {N} of the receiver,
    // evaluated at aGaussPoint.

    int ndofs = Shell7Base :: giveNumberOfDofs();
    int ndofs_xm  = this->giveNumberOfFieldDofs(Midplane);
    answer.resize(7, ndofs);
    answer.zero();
    FloatArray N;
    this->fei->evalN( N, lcoords, FEIElementGeometryWrapper(this) );

    /*   nno*3 nno*3 nno
     * 3 [N_x   0    0
     * 3   0   N_m   0
     * 1   0    0  N_gmm ]
     */

    if ( ei && dynamic_cast< Crack*>(ei) ) { 
        std :: vector< double > efGP;
        double levelSetGP = this->evaluateLevelSet(lcoords, ei);
        ei->evaluateEnrFuncAt(efGP, lcoords, levelSetGP);
        
        for ( int i = 1, j = 0; i <= this->giveNumberOfDofManagers(); i++, j += 3 ) {
            double factor = efGP [ 0 ] - EvaluateEnrFuncInDofMan(i, ei);

            answer.at(1, 1 + j) = N.at(i) * factor;
            answer.at(2, 2 + j) = N.at(i) * factor;
            answer.at(3, 3 + j) = N.at(i) * factor;
            answer.at(4, ndofs_xm + 1 + j) = N.at(i) * factor;
            answer.at(5, ndofs_xm + 2 + j) = N.at(i) * factor;
            answer.at(6, ndofs_xm + 3 + j) = N.at(i) * factor;
            answer.at(7, ndofs_xm * 2 + i) = N.at(i) * factor;
        }
        answer.times(DISC_DOF_SCALE_FAC);
    } else if ( ei && dynamic_cast< Delamination*>(ei) ) {
        Shell7Base :: computeNmatrixAt(lcoords, answer);
        std :: vector< double > efGP;
        double levelSetGP = this->evaluateLevelSet(lcoords, ei);
        ei->evaluateEnrFuncAt(efGP, lcoords, levelSetGP);
        if ( efGP[0] > 0.1 ) {
            answer.times( efGP[0] * DISC_DOF_SCALE_FAC );
        } else {
            answer.times(0.0);
        }
    } else {
        Shell7Base :: computeNmatrixAt(lcoords, answer);
    }
}



void
Shell7BaseXFEM :: edgeComputeEnrichedNmatrixAt(FloatArray &lcoords, FloatMatrix &answer, EnrichmentItem *ei)
{
// Returns the displacement interpolation matrix {N} of the receiver 
// evaluated at gaussPoint along one edge.

    answer.resize( 7, this->giveNumberOfEdgeDofs() );
    answer.zero();

    FloatArray N;
    this->fei->edgeEvalN( N, 1, lcoords, FEIElementGeometryWrapper(this) );

    /*    9   9    3
     * 3 [N_x   0    0
     * 3   0   N_m   0
     * 1   0    0  N_gmm ]
     */
    int ndofs_xm = this->giveNumberOfEdgeDofs() / 7 * 3;   // numEdgeNodes * 3 dofs

    if ( ei && dynamic_cast< Crack*>(ei) ) { 
        std :: vector< double > efGP;
        double levelSetGP = this->evaluateLevelSet(lcoords, ei);
        ei->evaluateEnrFuncAt(efGP, lcoords, levelSetGP);

        for ( int i = 1, j = 0; i <= this->giveNumberOfEdgeDofManagers(); i++, j += 3 ) {
            double factor = efGP [ 0 ] - EvaluateEnrFuncInDofMan(i, ei);
            answer.at(1, 1 + j)   = N.at(i) * factor;
            answer.at(2, 2 + j)   = N.at(i) * factor;
            answer.at(3, 3 + j)   = N.at(i) * factor;
            answer.at(4, ndofs_xm + 1 + j) = N.at(i) * factor;
            answer.at(5, ndofs_xm + 2 + j) = N.at(i) * factor;
            answer.at(6, ndofs_xm + 3 + j) = N.at(i) * factor;
            answer.at(7, ndofs_xm * 2 + i) = N.at(i) * factor;
        }
        answer.times(DISC_DOF_SCALE_FAC);

    } else if ( ei && dynamic_cast< Delamination*>(ei) ) {
        Shell7Base :: edgeComputeNmatrixAt(lcoords, answer);
        std :: vector< double > efGP;
        double levelSetGP = this->edgeEvaluateLevelSet(lcoords, ei);
        ei->evaluateEnrFuncAt(efGP, lcoords, levelSetGP);
        if ( efGP[0] > 0.1 ) {
            answer.times( efGP[0] * DISC_DOF_SCALE_FAC );
        } else {
            answer.times(0.0);
        }
    } else {
        Shell7Base :: edgeComputeNmatrixAt(lcoords, answer);
    }

}


void
Shell7BaseXFEM :: edgeComputeEnrichedBmatrixAt(FloatArray &lcoords, FloatMatrix &answer, EnrichmentItem *ei)
{
/* Returns the  matrix {B} of the receiver, evaluated at aGaussPoint. Such that
 * B*a = [dxbar_dxi, dwdxi, w, dgamdxi, gam]^T, where a is the vector of unknowns
 */

    answer.resize( 11, this->giveNumberOfEdgeDofs() );
    answer.zero();
    FloatArray N, dNdxi;



    if ( ei && dynamic_cast< Crack*>(ei) ) { 

        this->fei->edgeEvalN( N, 1, lcoords, FEIElementGeometryWrapper(this) );
        int iedge = 0;
        this->fei->edgeEvaldNdxi( dNdxi, iedge, lcoords, FEIElementGeometryWrapper(this) );

        /*
         * 3 [B_u   0    0
         * 3   0   B_w   0
         * 3   0   N_w   0
         * 1   0    0  B_gam
         * 1   0    0  N_gam]
         */

        // Evaluate enrichment function at point given by lcoords
        std :: vector< double > efGP;
        double levelSetGP = this->evaluateLevelSet(lcoords, ei);
        ei->evaluateEnrFuncAt(efGP, lcoords, levelSetGP);


        int ndofs_xm = this->giveNumberOfEdgeDofs() / 7 * 3;   // numEdgeNodes * 3 dofs
        int ndofman = this->giveNumberOfEdgeDofManagers();
        // First row
        for ( int i = 1, j = 0; i <= ndofman; i++, j += 3  ) {
            double factor = efGP [ 0 ] - EvaluateEnrFuncInDofMan(i, ei);
            
            answer.at(1, 1 + j) = dNdxi.at(i) * factor;
            answer.at(2, 2 + j) = dNdxi.at(i) * factor;
            answer.at(3, 3 + j) = dNdxi.at(i) * factor;
        }

        // Second row
        for ( int i = 1, j = 0; i <= ndofman; i++, j += 3  ) {
            double factor = efGP [ 0 ] - EvaluateEnrFuncInDofMan(i, ei);

            answer.at(4, ndofs_xm + 1 + j) = dNdxi.at(i) * factor;
            answer.at(5, ndofs_xm + 2 + j) = dNdxi.at(i) * factor;
            answer.at(6, ndofs_xm + 3 + j) = dNdxi.at(i) * factor;
            answer.at(7, ndofs_xm + 1 + j) = N.at(i) * factor;
            answer.at(8, ndofs_xm + 2 + j) = N.at(i) * factor;
            answer.at(9, ndofs_xm + 3 + j) = N.at(i) * factor;
        }

        // Third row
        for ( int i = 1, j = 0; i <= ndofman; i++, j += 1  ) {
            double factor = efGP [ 0 ] - EvaluateEnrFuncInDofMan(i, ei);

            answer.at(10, ndofs_xm * 2 + 1 + j) = dNdxi.at(i) * factor;
            answer.at(11, ndofs_xm * 2 + 1 + j) = N.at(i) * factor;
        }


        answer.times(DISC_DOF_SCALE_FAC);

    } else if ( ei && dynamic_cast< Delamination*>(ei) ){
        Shell7Base :: edgeComputeBmatrixAt(lcoords, answer);
        std :: vector< double > efGP;
        double levelSetGP = this->evaluateLevelSet(lcoords, ei);
        ei->evaluateEnrFuncAt(efGP, lcoords, levelSetGP);
        if ( efGP[0] > 0.1 ) {
            answer.times( efGP[0]*DISC_DOF_SCALE_FAC );
        } else {
            answer.times(0.0);
        }
    } else {
        Shell7Base :: edgeComputeBmatrixAt(lcoords, answer);
    }

}






// Delamination specific


void
Shell7BaseXFEM :: vtkEvalUpdatedGlobalCoordinateAt(FloatArray &localCoords, int layer, FloatArray &globalCoords, TimeStep *tStep)
{
    //double zeta = this->giveGlobalZcoordInLayer(localCoords.at(3), layer); 
    double zeta = this->giveGlobalZcoord(localCoords.at(3), localCoords); 

    // Continuous part
    FloatArray solVec;
    this->giveUpdatedSolutionVector(solVec, tStep);
    FloatArray xc, mc; double gamc=0;
    Shell7Base :: giveUnknownsAt(localCoords, solVec, xc, mc, gamc, tStep); 
    double fac = ( zeta + 0.5 * gamc * zeta * zeta );
    globalCoords = xc;    
    globalCoords.add(fac,mc);

#if 1
    // Discontinuous part
    std :: vector< double > ef;
    FloatArray solVecD, xd, md, xtemp(3), N; double gamd=0;
    for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { 
        EnrichmentItem *ei = this->xMan->giveEnrichmentItem(i);

        if ( ei->isElementEnriched(this) ) {
            IntArray eiDofIdArray;
            ei->giveEIDofIdArray(eiDofIdArray); 
            this->giveDisSolutionVector(solVecD, eiDofIdArray, tStep); 
            this->giveDisUnknownsAt(localCoords, ei, solVecD, xd, md, gamd, tStep);
            double fac = ( zeta + 0.5 * gamd * zeta * zeta );
            xtemp = xd;
            xtemp.add(fac,md);
            globalCoords.add(xtemp); 
  

        }
        
    }
#endif
}


void
Shell7BaseXFEM :: giveDisUnknownsAt(FloatArray &lcoords, EnrichmentItem *ei, FloatArray &solVec, FloatArray &x, FloatArray &m, double gam, TimeStep *tStep)
{
    // returns the unknowns evaluated at a point (xi1, xi2, xi3)
    FloatArray vec;
    FloatMatrix NEnr;
    this->computeEnrichedNmatrixAt(lcoords, NEnr, ei);
    vec.beProductOf(NEnr, solVec);
    x.resize(3); x.zero();
    m.resize(3); m.zero();
    x.at(1) = vec.at(1);
    x.at(2) = vec.at(2);
    x.at(3) = vec.at(3);
    m.at(1) = vec.at(4);
    m.at(2) = vec.at(5);
    m.at(3) = vec.at(6);
    gam  = vec.at(7);

}



void
Shell7BaseXFEM :: giveCompositeExportData(VTKPiece &vtkPiece, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep )
{

    int numSubCells = 1;
    if ( this->allTri.size() ) {
        numSubCells = (int)this->allTri.size();
    }
    
    int numLayers = this->layeredCS->giveNumberOfLayers();
    int numCells = numLayers * numSubCells;
    
    const int numCellNodes  = 15; // quadratic wedge
    int numTotalNodes = numCellNodes*numCells;

    vtkPiece.setNumberOfCells(numCells);
    vtkPiece.setNumberOfNodes(numTotalNodes);

    std::vector <FloatArray> nodeCoords;
    int val    = 1;
    int offset = 0;
    int currentCell = 1;
    IntArray nodes(numCellNodes);

    // Compute fictious node coords
    int nodeNum = 1;
    for ( int layer = 1; layer <= numLayers; layer++ ) {
        
        for ( int subCell = 1; subCell <= numSubCells; subCell++ ) {

            // Node coordinates
            if ( numSubCells == 1 ) {
                Shell7Base :: giveFictiousNodeCoordsForExport(nodeCoords, layer); 
            } else {
                this->giveFictiousNodeCoordsForExport(nodeCoords, layer, subCell);       
            }
        
            for ( int node = 1; node <= numCellNodes; node++ ) {    
                vtkPiece.setNodeCoords(nodeNum, nodeCoords[node-1] );
                nodeNum += 1;
            }

            // Connectivity       
            for ( int i = 1; i <= numCellNodes; i++ ) {            
                nodes.at(i) = val++;
            }
            vtkPiece.setConnectivity(currentCell, nodes);
        
            // Offset
            offset += numCellNodes;
            vtkPiece.setOffset(currentCell, offset);

            // Cell types
            vtkPiece.setCellType(currentCell, 26); // Quadratic wedge

            currentCell++;
        }
    }




    // Export nodal variables from primary fields        
    vtkPiece.setNumberOfPrimaryVarsToExport(primaryVarsToExport.giveSize(), numTotalNodes);

    std::vector<FloatArray> updatedNodeCoords;
    FloatArray u(3);
    std::vector<FloatArray> values;
    for ( int fieldNum = 1; fieldNum <= primaryVarsToExport.giveSize(); fieldNum++ ) {
        UnknownType type = ( UnknownType ) primaryVarsToExport.at(fieldNum);
        nodeNum = 1;
        int currentCell = 1;
        for ( int layer = 1; layer <= numLayers; layer++ ) {            
            
            for ( int subCell = 1; subCell <= numSubCells; subCell++ ) {

                if ( type == DisplacementVector ) { // compute displacement as u = x - X
                    if ( numSubCells == 1 ) {
                        Shell7Base :: giveFictiousNodeCoordsForExport(nodeCoords, layer);
                        this->giveFictiousUpdatedNodeCoordsForExport(updatedNodeCoords, layer, tStep, 0);
                    } else {
                        this->giveFictiousNodeCoordsForExport(nodeCoords, layer, subCell);
                        this->giveFictiousUpdatedNodeCoordsForExport(updatedNodeCoords, layer, tStep, subCell);
                    }
                    for ( int j = 1; j <= numCellNodes; j++ ) {
                        u = updatedNodeCoords[j-1];
                        u.subtract(nodeCoords[j-1]);
                        vtkPiece.setPrimaryVarInNode(fieldNum, nodeNum, u);
                        nodeNum += 1;        
                    }

                } else {
                    ZZNodalRecoveryMI_recoverValues(values, layer, ( InternalStateType ) 1, tStep); // does not work well - fix
                    for ( int j = 1; j <= numCellNodes; j++ ) {
                        vtkPiece.setPrimaryVarInNode(fieldNum, nodeNum, values[j-1]);
                        nodeNum += 1;
                    }
                }

                currentCell++;
            }
        }
    }

#if 0

    // Export nodal variables from internal fields
    
    vtkPiece.setNumberOfInternalVarsToExport( internalVarsToExport.giveSize(), numTotalNodes );
    for ( int fieldNum = 1; fieldNum <= internalVarsToExport.giveSize(); fieldNum++ ) {
        InternalStateType type = ( InternalStateType ) internalVarsToExport.at(fieldNum);
        nodeNum = 1;
        //this->recoverShearStress(tStep);
        for ( int layer = 1; layer <= numCells; layer++ ) {            
            recoverValuesFromIP(values, layer, type, tStep);        
            for ( int j = 1; j <= numCellNodes; j++ ) {
                vtkPiece.setInternalVarInNode( fieldNum, nodeNum, values[j-1] );
                //ZZNodalRecoveryMI_recoverValues(el.nodeVars[fieldNum], layer, type, tStep);          
                nodeNum += 1;        
            }                                
        }  
    }


    // Export cell variables
    FloatArray average;
    vtkPiece.setNumberOfCellVarsToExport(cellVarsToExport.giveSize(), numCells);
    for ( int i = 1; i <= cellVarsToExport.giveSize(); i++ ) {
        InternalStateType type = ( InternalStateType ) cellVarsToExport.at(i);;
      
        for ( int layer = 1; layer <= numCells; layer++ ) {     
            IntegrationRule *iRuleL = integrationRulesArray [ layer - 1 ];
            VTKXMLExportModule::computeIPAverage(average, iRuleL, this, type, tStep);
            vtkPiece.setCellVar(i, layer, convV6ToV9Stress(average) );  
        }

    }

#endif 


#if 0


    // Export of XFEM related quantities


    //int numCells = this->layeredCS->giveNumberOfLayers();
    //const int numCellNodes  = 15; // quadratic wedge
    //int numTotalNodes = numCellNodes*numCells;

    XfemManager *xFemMan =  this->xMan;
    int nEnrIt = xFemMan->giveNumberOfEnrichmentItems();


    IntArray wedgeToTriMap;
    wedgeToTriMap.setValues(15, 1, 2, 3, 1, 2, 3, 4, 5, 6, 4, 5, 6, 1, 2, 3 );


    vtkPiece.setNumberOfInternalXFEMVarsToExport(xFemMan->vtkExportFields.giveSize(), nEnrIt, numTotalNodes);
    for ( int field = 1; field <= xFemMan->vtkExportFields.giveSize(); field++ ) {
        XFEMStateType xfemstype = ( XFEMStateType ) xFemMan->vtkExportFields [ field - 1 ];
        
        for ( int enrItIndex = 1; enrItIndex <= nEnrIt; enrItIndex++ ) {
            EnrichmentItem *ei = xFemMan->giveEnrichmentItem(enrItIndex);
            int nodeNum = 1;
            for ( int layer = 1; layer <= numCells; layer++ ) {


                for ( int subCell = 1; subCell <= numSubCells; subCell++ ) {

                    for ( int nodeIndx = 1; nodeIndx <= numCellNodes; nodeIndx++ ) {

                        Node *node = this->giveNode( wedgeToTriMap.at(nodeIndx) );
                        FloatArray valueArray;
                        const FloatArray *val = NULL;
                        if ( xfemstype == XFEMST_LevelSetPhi ) {
                            valueArray.resize(1);
                            val = & valueArray;
                            ei->evalLevelSetNormalInNode( valueArray.at(1), node->giveNumber() );
                        } else if ( xfemstype == XFEMST_LevelSetGamma ) {
                            valueArray.resize(1);
                            val = & valueArray;
                            ei->evalLevelSetTangInNode( valueArray.at(1), node->giveNumber() );
                        } else if ( xfemstype == XFEMST_NodeEnrMarker ) {
                            valueArray.resize(1);
                            val = & valueArray;
                            ei->evalNodeEnrMarkerInNode( valueArray.at(1), node->giveNumber() );
                        } else {
                            //OOFEM_WARNING2("VTKXMLExportModule::getNodalVariableFromXFEMST: invalid data in node %d", inode);
                        }

                        vtkPiece.setInternalXFEMVarInNode(field, enrItIndex, nodeNum, valueArray);
                        nodeNum += 1;
                    }
                }
            }
        }
    }
#endif    


}


void 
Shell7BaseXFEM :: giveFictiousNodeCoordsForExport(std::vector<FloatArray> &nodes, int layer, int subCell)
{

    // compute fictious node coords
    FloatArray nodeLocalXi1Coords, nodeLocalXi2Coords, nodeLocalXi3Coords;

    // need to return local coordinates corresponding to the nodes of the sub triangles
    giveLocalNodeCoordsForExport(nodeLocalXi1Coords, nodeLocalXi2Coords, nodeLocalXi3Coords, subCell);


    nodes.resize(15);
    for ( int i = 1; i <= 15; i++ ){
        FloatArray coords, localCoords(3);
        localCoords.at(1) = nodeLocalXi1Coords.at(i);
        localCoords.at(2) = nodeLocalXi2Coords.at(i);
        localCoords.at(3) = nodeLocalXi3Coords.at(i);
        
        this->vtkEvalInitialGlobalCoordinateAt(localCoords, layer, coords);
        nodes[i-1].resize(3); 
        nodes[i-1] = coords;
    }

}


void 
Shell7BaseXFEM :: giveFictiousUpdatedNodeCoordsForExport(std::vector<FloatArray> &nodes, int layer, TimeStep *tStep, int subCell)
{
    // compute fictious node coords
    FloatArray nodeLocalXi1Coords, nodeLocalXi2Coords, nodeLocalXi3Coords;
    if ( subCell == 0) { 
        giveLocalNodeCoordsForExport(nodeLocalXi1Coords, nodeLocalXi2Coords, nodeLocalXi3Coords);

        // must get local z-coord in terms of the total thickness not layerwise
        

    } else {
        giveLocalNodeCoordsForExport(nodeLocalXi1Coords, nodeLocalXi2Coords, nodeLocalXi3Coords, subCell);
    }
    nodes.resize(15);
    for ( int i = 1; i <= 15; i++ ){
        FloatArray coords, localCoords(3);
        localCoords.at(1) = nodeLocalXi1Coords.at(i);
        localCoords.at(2) = nodeLocalXi2Coords.at(i);


        // Map local layer cs to local shell cs
        double scaleFactor = 0.9999; // Will be numerically unstable with xfem if the endpoints lie at +-1
        double totalThickness = this->layeredCS->computeIntegralThick();
        double zMid_i = this->layeredCS->giveLayerMidZ(layer); // global z-coord
        double xiMid_i = 1.0 - 2.0 * ( totalThickness - this->layeredCS->giveMidSurfaceZcoordFromBottom() - zMid_i ) / totalThickness; // local z-coord
        double deltaxi = nodeLocalXi3Coords.at(i) * this->layeredCS->giveLayerThickness(layer) / totalThickness; // distance from layer mid
        nodeLocalXi3Coords.at(i) = xiMid_i + deltaxi * scaleFactor;

        localCoords.at(3) = nodeLocalXi3Coords.at(i);
        this->vtkEvalUpdatedGlobalCoordinateAt(localCoords, layer, coords, tStep);
        nodes[i-1].resize(3); 
        nodes[i-1] = coords;
    }
}


void
Shell7BaseXFEM :: giveLocalNodeCoordsForExport(FloatArray &nodeLocalXi1Coords, FloatArray &nodeLocalXi2Coords, FloatArray &nodeLocalXi3Coords) {
    // Local coords for a quadratic wedge element (VTK cell type 26)
    double z = 0.999;
    nodeLocalXi1Coords.setValues(15, 1., 0., 0., 1., 0., 0., .5, 0., .5, .5, 0., .5, 1., 0., 0.);      
    nodeLocalXi2Coords.setValues(15, 0., 1., 0., 0., 1., 0., .5, .5, 0., .5, .5, 0., 0., 1., 0.);
    nodeLocalXi3Coords.setValues(15, -z, -z, -z,  z,  z,  z, -z, -z, -z,  z,  z,  z, 0., 0., 0.);
}


void
Shell7BaseXFEM :: giveLocalNodeCoordsForExport(FloatArray &nodeLocalXi1Coords, FloatArray &nodeLocalXi2Coords, FloatArray &nodeLocalXi3Coords, int subCell) {
    // Local coords for a quadratic wedge element (VTK cell type 26)
    double scale = 0.999;
    double z = 1.0*scale;
    //nodeLocalXi1Coords.setValues(15, 1., 0., 0., 1., 0., 0., .5, 0., .5, .5, 0., .5, 1., 0., 0.);      
    //nodeLocalXi2Coords.setValues(15, 0., 1., 0., 0., 1., 0., .5, .5, 0., .5, .5, 0., 0., 1., 0.);
    //nodeLocalXi3Coords.setValues(15, -z, -z, -z,  z,  z,  z, -z, -z, -z,  z,  z,  z, 0., 0., 0.);

    FloatArray g1, g2, g3;
    g1 = this->allTri[subCell-1].giveVertex(1);
    g2 = this->allTri[subCell-1].giveVertex(2);
    g3 = this->allTri[subCell-1].giveVertex(3);

    FloatArray gs1, gs2, gs3;
    
    double alpha1 = scale; double alpha2 = (1.0-alpha1)*0.5; double alpha3 = alpha2;
    gs1 = alpha1*g1 + alpha2*g2 + alpha3*g3;
    gs2 = alpha2*g1 + alpha1*g2 + alpha3*g3;
    gs3 = alpha2*g1 + alpha3*g2 + alpha1*g3;

    FloatArray loc1, loc2, loc3;
    this->computeLocalCoordinates(loc1, gs1);
    this->computeLocalCoordinates(loc2, gs2);
    this->computeLocalCoordinates(loc3, gs3);

    FloatArray loc12, loc23, loc31;
    loc12 = 0.5 * (loc1 + loc2);
    loc23 = 0.5 * (loc2 + loc3);
    loc31 = 0.5 * (loc3 + loc1);
    double a = loc1.at(1);
    double b = loc2.at(1);
    double c = loc3.at(1);
    double d = loc12.at(1);
    double e = loc23.at(1);
    double f = loc31.at(1);
    nodeLocalXi1Coords.setValues(15, a, b, c, a, b, c, d, e, f, d, e, f, a, b, c);      
    
    a = loc1.at(2);
    b = loc2.at(2);
    c = loc3.at(2);
    d = loc12.at(2);
    e = loc23.at(2);
    f = loc31.at(2);
    nodeLocalXi2Coords.setValues(15, a, b, c, a, b, c, d, e, f, d, e, f, a, b, c);

    nodeLocalXi3Coords.setValues(15, -z, -z, -z,  z,  z,  z, -z, -z, -z,  z,  z,  z, 0., 0., 0.);
}


} // end namespace oofem
