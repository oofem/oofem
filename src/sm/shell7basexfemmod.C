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

#include "shell7basexfemmod.h"
#include "modshell7base.h"
#include "xfem/enrichmentitem.h"
#include "xfem/xfemmanager.h"
#include "dofmanager.h"
#include "constantpressureload.h"
//#include "InterfaceMaterials\Old\simpleinterfacemat.h"
#include "connectivitytable.h"
#include "InterfaceMaterials/intmatbilinczfagerstrom.h"
#include "mathfem.h"
#include "node.h"
#include "geometry.h"
#include "cltypes.h"
#include "xfem/enrichmentitems/crack.h"
#include "xfem/enrichmentitems/shellcrack.h"
namespace oofem {

/* Scale factor fo the discontinuous dofs. Implies that the corresponding 
   dofs must be scaled with 1/factor in the input file
*/
const double DISC_DOF_SCALE_FAC = 1.0e-3; 

ModShell7BaseXFEM :: ModShell7BaseXFEM(int n, Domain *aDomain) : ModShell7Base(n, aDomain), XfemElementInterface(this) 
{
}

ModShell7BaseXFEM :: ~ModShell7BaseXFEM()
{
    for ( auto &iRule: czIntegrationRulesArray ) {
        delete iRule;
    }
}

int 
ModShell7BaseXFEM :: checkConsistency()
{
    ModShell7Base :: checkConsistency();
    return 1;
}

void
ModShell7BaseXFEM :: postInitialize()
{
    ModShell7Base :: postInitialize();
    this->xMan =  this->giveDomain()->giveXfemManager();
    

    // Set up ordering arrays and arrays with the actived dofs
    int numEI = xMan->giveNumberOfEnrichmentItems();
    this->orderingArrays.resize(numEI);
    this->activeDofsArrays.resize(numEI);
    for ( int i = 1; i <= numEI; i++ ) {
	EnrichmentItem *ei = this->xMan->giveEnrichmentItem(i);
	if ( ei->isElementEnriched(this) ) {
	    this->computeOrderingArray( this->orderingArrays[i-1], activeDofsArrays[i-1], ei); 
	}
    }
}

void
ModShell7BaseXFEM :: computeFailureCriteriaQuantities(FailureCriteriaStatus *fcStatus, TimeStep *tStep) 
{
    
    // Compute necessary quantities for evaluation of failure criterias
    ///@todo Ugly code and not tested in a while so probably broken

    if ( DamagedNeighborLayeredStatus * status = dynamic_cast< DamagedNeighborLayeredStatus * >(fcStatus) ) {
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

            ModShell7BaseXFEM *neighbor = 
                dynamic_cast< ModShell7BaseXFEM * > (this->giveDomain()->giveElement( neighbors.at(i) ));
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
    
}





IRResultType ModShell7BaseXFEM :: initializeFrom(InputRecord *ir)
{
    //IRResultType result;                   // Required by IR_GIVE_FIELD macro
    
    if ( ir->hasField(_IFT_ModShell7BaseXFEM_CohesiveZoneMaterial) ) {
        OOFEM_ERROR("'czmaterial' this keyword is not in use anymore! Instead define cz material for each interface in the cross secton, ex: interfacematerials 3 x x x ");
    }

    ModShell7Base :: initializeFrom(ir);
    return IRRT_OK; 
}

bool 
ModShell7BaseXFEM :: hasCohesiveZone(int interfaceNum)
{
    if ( interfaceNum <= this->layeredCS->giveNumberOfLayers()-1  ) { 
        return this->layeredCS->giveInterfaceMaterialNum(interfaceNum) > 0;
    } else {
        return 0;
    }
}

Interface
*ModShell7BaseXFEM :: giveInterface(InterfaceType it)
{
    if ( it != XfemElementInterfaceType ) {
        return ModShell7Base :: giveInterface(it);
    } else if ( it == XfemElementInterfaceType ) {
        return static_cast< XfemElementInterface * >(this);
    } else {
        return ModShell7Base :: giveInterface(it); ///@todo remove
    }
}


void
ModShell7BaseXFEM :: giveDofManDofIDMask(int inode, IntArray &answer) const
{
    // Returns the total id mask of the dof manager - regular id's + enriched id's

    // Continuous part
    ModShell7Base :: giveDofManDofIDMask(inode, answer);
    XfemManager *xMan = this->giveDomain()->giveXfemManager(); // xman not initialized on el. level when this is first called

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
ModShell7BaseXFEM :: evalCovarBaseVectorsAt(FloatArray &lCoords, FloatMatrix &gcov, FloatArray &genEpsC)
{
    // Continuous part g_c
    ModShell7Base :: evalCovarBaseVectorsAt(lCoords, gcov, genEpsC);

    // Discontinuous part g_d
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    FloatArray dGenEps;
    FloatMatrix gcovd; 
    std :: vector< double > ef;
    for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { 
        EnrichmentItem *ei = this->xMan->giveEnrichmentItem(i);

        if ( ei->isElementEnriched(this) ) {
            computeDiscGeneralizedStrainVector(dGenEps, lCoords, ei, tStep); 
            ModShell7Base :: evalCovarBaseVectorsAt(lCoords, gcovd, dGenEps);
            gcov.add(gcovd); 
        }
    }
}


void
ModShell7BaseXFEM :: computeDiscGeneralizedStrainVector(FloatArray &answer, FloatArray &lCoords, EnrichmentItem *ei, TimeStep *tStep)
{
    FloatArray solVecD;
    IntArray eiDofIdArray;
    ei->giveEIDofIdArray(eiDofIdArray); 
    this->computeDiscSolutionVector(eiDofIdArray, tStep, solVecD);
    FloatMatrix B;
    this->computeEnrichedBmatrixAt(lCoords, B, ei);
    answer.beProductOf(B, solVecD);
}


void 
ModShell7BaseXFEM :: computeDiscSolutionVector(IntArray &dofIdArray , TimeStep *tStep, FloatArray &answer)
{
    FloatArray temp;
    IntArray newIdArray;
    newIdArray = dofIdArray;
    newIdArray.followedBy(0); // Extend it by one such that the length is 7 -> answer = size 42
    this->computeVectorOf(newIdArray, VM_Total, tStep, temp, true);
    answer.resize(temp.giveSize());
    answer.zero();
    answer.assemble( temp, this->giveOrderingNodes() );
}

int
ModShell7BaseXFEM :: giveNumberOfDofs()
{
    // Continuous part
    int nDofs = ModShell7Base :: giveNumberOfDofs();

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
ModShell7BaseXFEM :: edgeGiveUpdatedSolutionVector(FloatArray &answer, const int iedge, TimeStep *tStep)
{
    ModShell7Base :: edgeGiveUpdatedSolutionVector(answer, iedge, tStep);
    answer.times(DISC_DOF_SCALE_FAC); ///@todo should this be here?
}


void 
ModShell7BaseXFEM :: computeOrderingArray( IntArray &orderingArray, IntArray &activeDofsArray,  EnrichmentItem *ei)
{
    // Routine to extract vector given an array of dofid items
    // If a certain dofId does not exist a zero is used as value

    const IntArray &ordering_cont = this->giveOrderingDofTypes();
    IntArray fieldDofId; 
    ModShell7Base :: giveDofManDofIDMask(0, fieldDofId);


    IntArray ordering_temp, activeDofsArrayTemp;
    ordering_temp.resize(ordering_cont.giveSize());
    activeDofsArrayTemp.resize(ordering_cont.giveSize());


    int activeDofPos = 0, activeDofIndex = 0, orderingDofIndex = 0;
    
    IntArray dofManDofIdMask, dofManDofIdMaskAll;
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        DofManager *dMan = this->giveDofManager(i);
        

        if ( ei == NULL) { // return mask corresponding to the regular id's
            ModShell7Base :: giveDofManDofIDMask(i, dofManDofIdMask);
        } else {
            if ( ei->isDofManEnriched(*dMan) ) {
                ei->giveEIDofIdArray(dofManDofIdMask); 
            }
        }
    
        for (int j = 1; j <= dofManDofIdMask.giveSize(); j++ ) {
            int pos = dMan->findDofWithDofId( ( DofIDItem ) dofManDofIdMask.at(j) ) - dMan->begin() + 1;
            activeDofPos++;
            ordering_temp      .at(activeDofPos) = orderingDofIndex + pos;
            activeDofsArrayTemp.at(activeDofPos) = activeDofIndex   + j;
        }
        this->giveDofManDofIDMask(i, dofManDofIdMaskAll);
        orderingDofIndex += dofManDofIdMaskAll.giveSize();
        activeDofIndex   += fieldDofId.giveSize();

        dofManDofIdMask.clear();
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
ModShell7BaseXFEM :: giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord)
{
    // Computes internal forces as the array of: [f_c, f_d1, ..., f_dn]
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
	  this->computeDiscSolutionVector(eiDofIdArray, tStep, solVecD);
	  this->discComputeSectionalForces(temp, tStep, solVec, solVecD, ei);

	  tempRed.beSubArrayOf(temp, this->activeDofsArrays[i-1]);
	  answer.assemble(tempRed, this->orderingArrays[i-1]);

	  // Cohesive zone model
	  if ( this->hasCohesiveZone(i) ) {
	      this->computeCohesiveForces( fCZ, tStep, solVec, solVecD, this->xMan->giveEnrichmentItem(i)); 
	      tempRed.beSubArrayOf(fCZ, this->activeDofsArrays[i-1]);
	      answer.assemble(tempRed, this->orderingArrays[i-1]);
	  }
      }
    }

}


void
ModShell7BaseXFEM :: discComputeSectionalForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, FloatArray &solVecD, EnrichmentItem *ei)
//
{
    int ndofs = ModShell7Base :: giveNumberOfDofs();
    int numberOfLayers = this->layeredCS->giveNumberOfLayers(); 
    FloatArray f(ndofs), genEps, genEpsD, ftemp, lCoords, Nd;
    FloatMatrix B, BEnr;
    f.zero();

    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRuleL = integrationRulesArray [ layer - 1 ];
        Material *mat = domain->giveMaterial( this->layeredCS->giveLayerMaterial(layer) );

        for ( GaussPoint *gp: *iRuleL ) {
            lCoords = * gp->giveNaturalCoordinates();

            this->computeEnrichedBmatrixAt(lCoords, B, NULL);
            this->computeEnrichedBmatrixAt(lCoords, BEnr, ei);
	    genEps.beProductOf(B, solVec);
	    genEpsD.beProductOf(BEnr, solVecD);
            double zeta = giveGlobalZcoord(lCoords);
            this->computeSectionalForcesAt(Nd, gp, mat, tStep, genEps, genEpsD, zeta);

            // Computation of nodal forces: f = B^t*[N M T Ms Ts]^t
            ftemp.beTProductOf( BEnr, Nd );
            double dV = this->computeVolumeAroundLayer(gp, layer);
            f.add( dV, ftemp );            


        }
    }

    answer.resize(ndofs); answer.zero();
    const IntArray &ordering = this->giveOrderingDofTypes();
    answer.assemble(f, ordering);

}


void
ModShell7BaseXFEM :: computeSectionalForcesAt(FloatArray &sectionalForces, IntegrationPoint *ip, Material *mat, TimeStep *tStep, FloatArray &genEps, FloatArray &genEpsD, double zeta)
{
    // New, in terms of PK1 stress
    // \Lambda_i * P * G^I
    FloatArray PG1(3), PG2(3), PG3(3);
    FloatMatrix lambda[3], Gcon, P, PG;
    FloatArray PVector, temp;
    this->computeStressMatrix(P, genEps, ip, mat, tStep);

    FloatArray lCoords = *ip->giveNaturalCoordinates();
    this->evalInitialContravarBaseVectorsAt(lCoords, Gcon);
    PG.beProductOf(P,Gcon);
    PG1.beColumnOf(PG, 1);
    PG2.beColumnOf(PG, 2);
    PG3.beColumnOf(PG, 3);
    //this->computeLambdaGMatrices(lambda, genEpsD, zeta); // associated with the variation of the test functions   
    this->computeLambdaGMatricesDis(lambda, zeta);
    // f = lambda_1^T * P*G^1 + lambda_2^T * P*G^2 + lambda_3^T * P*G^3
    sectionalForces.resize(18);
    sectionalForces.zero();
    temp.beTProductOf(lambda[0], PG1);
    sectionalForces.add(temp);
    temp.beTProductOf(lambda[1], PG2);
    sectionalForces.add(temp);
    temp.beTProductOf(lambda[2], PG3);
    sectionalForces.add(temp);

}













double 
ModShell7BaseXFEM :: evaluateLevelSet(const FloatArray &lCoords, EnrichmentItem *ei)
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
        OOFEM_ERROR("error in evaluation of levelset");
    }          
    return levelSet;
}

double 
ModShell7BaseXFEM :: edgeEvaluateLevelSet(const FloatArray &lCoords, EnrichmentItem *ei)
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
        OOFEM_ERROR("error in evaluation of levelset");
    }          
    return levelSet;
}



double 
ModShell7BaseXFEM :: evaluateHeavisideGamma(double xi, EnrichmentItem *ei)
{
    // Evaluates the Heaviside product H(gamma_1) * H(gamma_2)

    double xiBottom = dynamic_cast< ShellCrack * >( ei )->xiBottom;
    double xiTop    = dynamic_cast< ShellCrack * >( ei )->xiTop;

    if ( ( xiBottom <= xi ) && ( xi <= xiTop ) ) {
            return 1.0;
    } else {
            return 0.0;
    }          

}


void
ModShell7BaseXFEM :: giveMaxCZDamages(FloatArray &answer, TimeStep *tStep)
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
            for ( GaussPoint *gp: *iRule ) {

                mat = static_cast < StructuralInterfaceMaterial * > (this->layeredCS->giveInterfaceMaterial(i+1) );
                mat->giveIPValue(ipValues, gp, IST_DamageScalar, tStep);

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
ModShell7BaseXFEM :: computeCohesiveForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVecC, FloatArray &solVecD, EnrichmentItem *ei)
{
    //Computes the cohesive nodal forces for a given interface

    Delamination *dei =  dynamic_cast< Delamination * >( ei ); 
    if ( dei ) { // For now only add cz for delaminations
        FloatArray answerTemp, lCoords(3);
        answerTemp.resize(ModShell7Base :: giveNumberOfDofs() ); 
        answerTemp.zero();

        FloatMatrix NEnr, BEnr, F;  
        int delamNum = dei->giveNumber();
        IntegrationRule *iRuleL = czIntegrationRulesArray [ delamNum - 1 ]; ///@ todo does this work with giveNumber? Probably not, will have to save num delam.

        StructuralInterfaceMaterial *intMat = static_cast < StructuralInterfaceMaterial * > 
            (this->layeredCS->giveInterfaceMaterial(delamNum) );

        FloatMatrix lambda, lambdaN, Q;
        FloatArray Fp, T, nCov, xd, unknowns, genEpsC, genEpsD;

    for ( GaussPoint *gp: *iRuleL ) {
        lCoords.at(1) = gp->giveNaturalCoordinate(1);
        lCoords.at(2) = gp->giveNaturalCoordinate(2);
        lCoords.at(3) = dei->giveDelamXiCoord();

        this->computeEnrichedBmatrixAt(lCoords, BEnr, dei);
        this->computeEnrichedNmatrixAt(lCoords, NEnr, dei);

        // Lambda matrix
        genEpsD.beProductOf(BEnr, solVecD);
        double xi = dei->giveDelamXiCoord();
        double zeta = xi * this->layeredCS->computeIntegralThick() * 0.5;
        this->computeLambdaNMatrix(lambda, genEpsD, zeta);
        
        // Compute jump vector
        unknowns.beProductOf(NEnr, solVecD);
        xd.beProductOf(lambda,unknowns); // spatial jump
	genEpsC.beProductOf(BEnr, solVecC);
	this->computeFAt(lCoords, F, genEpsC);

	// Transform xd and F to a local coord system
        this->evalInitialCovarNormalAt(nCov, lCoords);
        Q.beLocalCoordSys(nCov);
        xd.rotatedWith(Q,'n');
        F.rotatedWith(Q,'n');

        // Compute cohesive traction based on jump
        intMat->giveFirstPKTraction_3d(T, gp, xd, F, tStep);
	lambdaN.beProductOf(lambda,NEnr);
	//Q.printYourself();	
	//T.printYourself();
	T.rotatedWith(Q,'t'); // transform back to global coord system
//lambdaN.printYourself();

	Fp.beTProductOf(lambdaN, T);
	//Fp.printYourself();
        double dA = this->computeAreaAround(gp,xi);
        answerTemp.add(dA*DISC_DOF_SCALE_FAC,Fp);
        }

        int ndofs = ModShell7Base :: giveNumberOfDofs();
        answer.resize(ndofs); answer.zero();
        const IntArray &ordering = this->giveOrderingDofTypes();
        answer.assemble(answerTemp, ordering);
    }
}

void
ModShell7BaseXFEM :: updateYourself(TimeStep *tStep)
// Updates the receiver at end of step.
{
    ModShell7Base :: updateYourself(tStep);

    for ( int i = 0; i < this->layeredCS->giveNumberOfLayers() - 1; i++ ) {
        if ( this->hasCohesiveZone(i+1) ) {
            czIntegrationRulesArray [ i ]->updateYourself(tStep);
        }
    }
}

void
ModShell7BaseXFEM :: computeCohesiveTangent(FloatMatrix &answer, TimeStep *tStep)
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
	    this->computeDiscSolutionVector(eiDofIdArray, tStep, solVecD);
            this->computeCohesiveTangentAt(temp, tStep, solVecD, dei);
            // Assemble part correpsonding to active dofs
            this->computeOrderingArray(orderingJ, activeDofsJ, dei); //
            tempRed.beSubMatrixOf(temp, activeDofsJ, activeDofsJ);
            answer.assemble(tempRed, orderingJ, orderingJ);
        }
    }
}



void
ModShell7BaseXFEM :: computeCohesiveTangentAt(FloatMatrix &answer, TimeStep *tStep, FloatArray &solVecD, 
    Delamination *dei)
{
    //Computes the cohesive tangent for a given interface K_cz = N^t*lambda*t * dT/dj * lambda * N
    FloatArray lCoords(3);
    FloatMatrix answerTemp, N, lambda, K, temp, tangent;
    int nDofs = ModShell7Base :: giveNumberOfDofs();
    answerTemp.resize(nDofs, nDofs); 
    answerTemp.zero();

    int delamNum = dei->giveNumber();
    IntegrationRule *iRuleL = czIntegrationRulesArray [ dei->giveNumber() - 1 ];
    StructuralInterfaceMaterial *intMat = static_cast < StructuralInterfaceMaterial * > 
        (this->layeredCS->giveInterfaceMaterial(delamNum) );

    double xi = dei->giveDelamXiCoord();

    FloatMatrix Q;
    FloatArray nCov;
    FloatArray interfaceXiCoords;
    this->layeredCS->giveInterfaceXiCoords(interfaceXiCoords);
    for ( GaussPoint *gp: *iRuleL ) {
        lCoords.at(1) = gp->giveNaturalCoordinate(1);
        lCoords.at(2) = gp->giveNaturalCoordinate(2);
        lCoords.at(3) = dei->giveDelamXiCoord();

        double zeta = giveGlobalZcoord( lCoords);
        this->computeLambdaNMatrixDis( lambda, zeta );
        this->computeEnrichedNmatrixAt(lCoords, N, dei);
                
        intMat->give3dStiffnessMatrix_dTdj(K, TangentStiffness, gp, tStep);
        this->evalInitialCovarNormalAt(nCov, lCoords);
        Q.beLocalCoordSys(nCov);
        K.rotatedWith(Q,'t');   // rotate back to global coord system

        this->computeTripleProduct(temp, lambda, K, lambda); ///@todo Fix JB 090814
        this->computeTripleProduct(tangent, N, temp, N);
        double dA = this->computeAreaAround(gp, xi);
        answerTemp.add(dA,tangent); 
    }

    answer.resize(nDofs, nDofs); 
    answer.zero();
    const IntArray &ordering = this->giveOrderingDofTypes();
    answer.assemble(answerTemp, ordering, ordering);
}



void 
ModShell7BaseXFEM :: computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep)
{

    int ndofs = this->giveNumberOfDofs();
    answer.resize(ndofs, ndofs);
    answer.zero();
    
    int numberOfLayers = this->layeredCS->giveNumberOfLayers();     
    FloatMatrix tempRed, tempRedT;
    FloatMatrix KCC, KCD, KDD, Bc;
    IntArray orderingC, activeDofsC;
    this->computeOrderingArray(orderingC, activeDofsC, NULL); 

    FloatMatrix A [ 3 ] [ 3 ];
    FloatArray lCoords; 
    
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        IntegrationRule *iRule = integrationRulesArray [ layer - 1 ];
        StructuralMaterial *layerMaterial = static_cast< StructuralMaterial* >( domain->giveMaterial( this->layeredCS->giveLayerMaterial(layer) ) );

        for ( GaussPoint *gp: *iRule ) {
            lCoords = * gp->giveNaturalCoordinates();
            this->computeEnrichedBmatrixAt(lCoords, Bc, NULL);
            
            ModShell7Base :: computeLinearizedStiffness(gp, layerMaterial, tStep, A);
            
            // Continuous part K_{c,c}
            this->discComputeBulkTangentMatrix(KCC, gp, NULL, NULL, layer, A, tStep);
            answer.assemble(KCC, orderingC, orderingC);

            // Discontinuous part
            int numEI = this->xMan->giveNumberOfEnrichmentItems();
            for ( int m = 1; m <= numEI; m++ ) { 
                EnrichmentItem *eiM = this->xMan->giveEnrichmentItem(m);

                if ( eiM->isElementEnriched(this) ) {
                          
                    this->discComputeBulkTangentMatrix(KCD, gp, NULL, eiM, layer, A, tStep);
                    tempRed.beSubMatrixOf(KCD, activeDofsC, this->activeDofsArrays[m-1]);
                    answer.assemble(tempRed, orderingC, this->orderingArrays[m-1]);
                    tempRedT.beTranspositionOf(tempRed);
                    answer.assemble(tempRedT, this->orderingArrays[m-1], orderingC);

                    // K_{dk,dl}
                    for ( int k = 1; k <= numEI; k++ ) {
                        EnrichmentItem *eiK = this->xMan->giveEnrichmentItem(k);
                        
                        if ( eiK->isElementEnriched(this) ) {
                            this->discComputeBulkTangentMatrix(KDD, gp, eiM, eiK, layer, A, tStep);
                            if ( this->activeDofsArrays[m-1].giveSize() != 0 && this->activeDofsArrays[k-1].giveSize() != 0 ) {
                                tempRed.beSubMatrixOf(KDD, this->activeDofsArrays[m-1], this->activeDofsArrays[k-1]);
                                answer.assemble(tempRed, this->orderingArrays[m-1], this->orderingArrays[k-1]);
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
    this->computeCohesiveTangent(Kcoh, tStep);
    //Kcoh.times(DISC_DOF_SCALE_FAC*DISC_DOF_SCALE_FAC);
    answer.add(Kcoh);

#endif


    // Add contribution due to pressure load
#if 0

    int nLoads = this->boundaryLoadArray.giveSize() / 2;

    for ( int k = 1; k <= nLoads; k++ ) {     // For each pressure load that is applied
        int load_number = this->boundaryLoadArray.at(2 * k - 1);
        int iSurf = this->boundaryLoadArray.at(2 * k);         // load_id
        Load *load = this->domain->giveLoad(load_number);
        std :: vector< double > efM, efK;

        if ( ConstantPressureLoad * pLoad = dynamic_cast< ConstantPressureLoad * >(load) ) {

            IntegrationRule *iRule = specialIntegrationRulesArray [ 1 ];

            for ( GaussPoint *gp: *iRule ) {
                this->computePressureTangentMatrixDis(KCC, KCD, KDD, gp, load, iSurf, tStep);
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
    
}


void
ModShell7BaseXFEM :: discComputeBulkTangentMatrix(FloatMatrix &KdIJ, IntegrationPoint *ip, EnrichmentItem *ei1, EnrichmentItem *ei2, int layer, FloatMatrix A [ 3 ] [ 3 ],  TimeStep *tStep)
{

    FloatMatrix temp, B1, B2;
    FloatMatrix lambda1 [ 3 ], lambda2 [ 3 ];

    FloatArray lCoords = * ip->giveNaturalCoordinates();
    this->computeEnrichedBmatrixAt(lCoords, B1, ei1);
    this->computeEnrichedBmatrixAt(lCoords, B2, ei2);

    double zeta = giveGlobalZcoord( lCoords );
    if ( ei1 ) {
        this->computeLambdaGMatricesDis(lambda1, zeta);
    } else {
        FloatArray solVecC, genEpsC;
        this->giveUpdatedSolutionVector(solVecC, tStep);  
	genEpsC.beProductOf(B1,solVecC);
        this->computeLambdaGMatrices(lambda1, genEpsC, zeta);
    }
    
    if ( ei2 ) {
        this->computeLambdaGMatricesDis(lambda2, zeta);
    } else {
        FloatArray solVecC, genEpsC;
        this->giveUpdatedSolutionVector(solVecC, tStep);  
	genEpsC.beProductOf(B2,solVecC);
        this->computeLambdaGMatrices(lambda2, genEpsC, zeta);
    }

    double dV = this->computeVolumeAroundLayer(ip, layer);

    FloatMatrix LAL(18,18); 
    LAL.zero();
    for ( int j = 0; j < 3; j++ ) {
        for ( int k = 0; k < 3; k++ ) {
            this->computeTripleProduct(temp, lambda1 [ j ], A [ j ][ k ], lambda2 [ k ] ); ///@todo Fix 090814
            LAL.add(dV,temp);
        }
    }

    FloatMatrix KDDtemp;
    this->computeTripleProduct(KDDtemp, B1, LAL, B2 );    ///@todo Fix 090814
    
    int ndofs = ModShell7Base :: giveNumberOfDofs();
    KdIJ.resize(ndofs,ndofs);
    KdIJ.zero();
    const IntArray &ordering = this->giveOrderingDofTypes();
    
    KdIJ.assemble(KDDtemp, ordering, ordering);

}


void
ModShell7BaseXFEM :: computeLambdaGMatricesDis(FloatMatrix lambda [ 3 ], double zeta)
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
ModShell7BaseXFEM :: computeLambdaNMatrixDis(FloatMatrix &lambda_xd, double zeta)
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
ModShell7BaseXFEM :: computePressureTangentMatrixDis(FloatMatrix &KCC, FloatMatrix &KCD, FloatMatrix &KDD, IntegrationPoint *ip, Load *load, const int iSurf, TimeStep *tStep)
{
    // Computes tangent matrix associated with the linearization of pressure loading. Assumes constant pressure.
    ConstantPressureLoad* pLoad = dynamic_cast< ConstantPressureLoad * >( load );
    
    FloatMatrix N, B, L(7, 18), gcov, W1, W2;
    FloatArray lCoords(3), solVec, pressure;
    FloatArray g1, g2, genEps;
    FloatMatrix lambdaGC [ 3 ], lambdaNC, lambdaGD [ 3 ], lambdaND;
    double xi   = pLoad->giveLoadOffset();
    lCoords = *ip->giveNaturalCoordinates();
    lCoords.at(3) = xi;     // local coord where load is applied
    double zeta = this->giveGlobalZcoord( lCoords );
    this->giveUpdatedSolutionVector(solVec, tStep);

    // compute w1,w2, KC
    this->computeNmatrixAt(lCoords, N);
    this->computeBmatrixAt(lCoords, B);
    genEps.beProductOf(B, solVec);  

    
    FloatMatrix LCC(18,18), LCD(18,18), LDD(18,18); 
    LCC.zero(); LCD.zero(); LDD.zero();
        
    //(xc+xd)*(g1xg2)=xc*g1xg2 + xd*g1xg2 -> xc*(W2*Dg1 - W1*Dg2) + xd*(W2*Dg1 - W1*Dg2)
    // Traction tangent, L =  lambdaN * ( W2*lambdaG_1 - W1*lambdaG_2  ) 
    load->computeValueAt(pressure, tStep, * ( ip->giveNaturalCoordinates() ), VM_Total);        // pressure component
    this->evalCovarBaseVectorsAt(lCoords, gcov, genEps);
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
    this->computeTripleProduct(KCCtemp, N, LCC, B );    ///@todo Fix 090814
    this->computeTripleProduct(KCDtemp, N, LCD, B );    
    this->computeTripleProduct(KDDtemp, N, LDD, B );    
    
    int ndofs = ModShell7Base :: giveNumberOfDofs();
    KCC.resize(ndofs,ndofs); KCD.resize(ndofs,ndofs); KDD.resize(ndofs,ndofs); 
    KCC.zero(); KCD.zero(); KDD.zero();
    const IntArray &ordering = this->giveOrderingDofTypes();
    KCC.assemble(KCCtemp, ordering, ordering);
    KCD.assemble(KCDtemp, ordering, ordering);
    KDD.assemble(KDDtemp, ordering, ordering);
}

void
ModShell7BaseXFEM :: computeMassMatrixNum(FloatMatrix &answer, TimeStep *tStep)
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

        for ( GaussPoint *gp: *iRule ) {

            FloatMatrix N11, N22, N33, N;
            //this->computeNmatricesAt(gp, N11, N22, N33);
            this->computeNmatrixAt(gp, N);
            FloatArray xbar, m;
            double gam = 0.;
            //this->computeSolutionFields(xbar, m, gam, solVec, N11, N22, N33);
            FloatArray localCoords = * gp->giveNaturalCoordinates();
            this->giveUnknownsAt(localCoords, solVec, xbar, m, gam, tStep);
            //this->computeNmatrixAt(gp, N);
            //unknowns.beProductOf(N,a); // [xbar, m, gam]^T
            //m = {unknowns.at(4), unknowns.at(5), unknowns.at(6) };
            //double gam = unknowns.at(7);


            /* Consistent mass matrix M = int{N^t*mass*N}
             *
             *         3    3    1
             *         3 [a*I  b*I   c*m      [A  B  C
             * mass =   3       d*I   e*m    =     D  E
             *         1  sym       f*m.m]     sym   F]
             */


            double zeta = giveGlobalZcoord( gp->giveNaturalCoordinate(3) );
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
ModShell7BaseXFEM :: computeEdgeLoadVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *tStep, ValueModeType mode)
{
    BoundaryLoad *edgeLoad = dynamic_cast< BoundaryLoad * >(load);
    if ( edgeLoad ) {
        answer.resize( this->computeNumberOfDofs() );
        answer.zero();

        // Continuous part
        FloatArray fT;
        ModShell7Base :: computeTractionForce(fT, iEdge, edgeLoad, tStep, mode);

        IntArray activeDofs, ordering; 
        this->computeOrderingArray(ordering, activeDofs, NULL); 
        answer.assemble(fT, ordering);        

        FloatArray componentsTemp, coordsTemp(1);
        coordsTemp.at(1) = 0.0; // 
        edgeLoad->computeValueAt(componentsTemp, tStep, coordsTemp, VM_Total);
        ///@todo Add support for load defined over certain area
	//double xi = 0.0; // defaults to geometric midplane
        //if ( componentsTemp.giveSize() == 8 ) {
        //    xi = componentsTemp.at(8);   // use the 8th component to store the-xi coord where the load acts
        //}

        // Disccontinuous part
        FloatArray temp, tempRed;
        for ( int i = 1; i <= this->xMan->giveNumberOfEnrichmentItems(); i++ ) { // Only one is supported at the moment

            EnrichmentItem *ei = this->xMan->giveEnrichmentItem(i);
            if ( ei->isElementEnriched(this) ) {
                this->computeEnrTractionForce(temp, iEdge, edgeLoad, tStep, mode, ei);
                this->computeOrderingArray(ordering, activeDofs, ei); 
                tempRed.beSubArrayOf(temp, activeDofs);
                answer.assemble(tempRed, ordering);                    
            }
        }
        return;
    } else {
        OOFEM_ERROR("load type not supported");
        return;
    }
}


void
ModShell7BaseXFEM :: computeEnrTractionForce(FloatArray &answer, const int iEdge, BoundaryLoad *edgeLoad, TimeStep *tStep, ValueModeType mode, EnrichmentItem *ei)
{
    // 
    IntegrationRule *iRule = specialIntegrationRulesArray [ 2 ];   // rule #3 for edge integration of distributed loads given in [*/m]
    GaussPoint *gp;

    FloatMatrix N, Q;
    FloatArray fT(7), components, lCoords;
    Load :: CoordSystType coordSystType = edgeLoad->giveCoordSystMode();

    FloatArray Nftemp(21), Nf(21);
    Nf.zero();
    for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
        gp = iRule->getIntegrationPoint(i);
        FloatArray lCoords = *gp->giveNaturalCoordinates();

        edgeLoad->computeValueAt(components, tStep, lCoords, mode);
        this->edgeComputeEnrichedNmatrixAt(lCoords, N, ei);

        if ( coordSystType ==  Load :: CST_UpdatedGlobal ) {
            
            // Updated global coord system
            FloatMatrix gcov;
            this->edgeEvalEnrCovarBaseVectorsAt(lCoords, iEdge, gcov, tStep, ei); 
            Q.beTranspositionOf(gcov);

            FloatArray distrForces(3), distrMoments(3), t1, t2;
            distrForces = { components.at(1), components.at(2), components.at(3) };
            distrMoments = { components.at(4), components.at(5), components.at(6) };
            t1.beTProductOf(Q, distrForces);
            t2.beTProductOf(Q, distrMoments);
            fT.addSubVector(t1,1);
            fT.addSubVector(t2,4);
            fT.at(7) = components.at(7); // don't do anything with the 'gamma'-load

        } else if( coordSystType == Load :: CST_Global ) { 
            // Undeformed global coord system
            for ( int i = 1; i <= 7; i++) {
                fT.at(i) = components.at(i);
            }
        } else {
            OOFEM_ERROR("ModShell7Base :: computeTractionForce - does not support local coordinate system");
        }

        double dL = this->edgeComputeLengthAround(gp, iEdge);        
        
        Nftemp.beTProductOf(N, fT*dL);
        Nf.add(Nftemp);
    }

    IntArray mask;
    this->giveEdgeDofMapping(mask, iEdge);
    answer.resize( ModShell7Base :: giveNumberOfDofs()  );
    answer.zero();
    answer.assemble(Nf, mask);

}

void
ModShell7BaseXFEM :: edgeEvalEnrCovarBaseVectorsAt(FloatArray &lcoords, const int iedge, FloatMatrix &gcov, TimeStep *tStep, EnrichmentItem *ei)
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
    dxdxi = { genEpsEdge.at(1), genEpsEdge.at(2), genEpsEdge.at(3) };
    dmdxi = { genEpsEdge.at(4), genEpsEdge.at(5), genEpsEdge.at(6) };
    m = { genEpsEdge.at(7), genEpsEdge.at(8), genEpsEdge.at(9) };
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
ModShell7BaseXFEM :: computeSurfaceLoadVectorAt(FloatArray &answer, Load *load,
                                         int iSurf, TimeStep *tStep, ValueModeType mode)
{
    BoundaryLoad *surfLoad = dynamic_cast< BoundaryLoad * >(load);

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
        if ( ConstantPressureLoad * pLoad = dynamic_cast< ConstantPressureLoad * >(load) ) {
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
		   this->computeDiscSolutionVector(eiDofIdArray, tStep, solVecD);
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
        OOFEM_ERROR("load type not supported");
        return;
    }
}



// Shifted N and B matrices
void
ModShell7BaseXFEM :: computeEnrichedBmatrixAt(FloatArray &lcoords, FloatMatrix &answer, EnrichmentItem *ei)
{
    // Returns the enriched and shifted {B} matrix of the receiver, evaluated at gp. Such that
    // B*a = [dxbar_dxi, dwdxi, w, dgamdxi, gam]^T, where a is the vector of unknowns
 
    if ( ei && dynamic_cast< Crack*>(ei) ) { 

        int ndofs = ModShell7Base :: giveNumberOfDofs();
        int ndofs_xm  = 3 * this->giveNumberOfDofManagers();
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
        //printf("enr func in gp = %e \n", efGP[0]);
        int ndofman = this->giveNumberOfDofManagers();

        // First column
        for ( int i = 1, j = 0; i <= ndofman; i++, j += 3 ) {
            double factor = efGP [ 0 ] - EvaluateEnrFuncInDofMan(i, ei);
            //printf("enr func in dofman %d = %e \n", i, EvaluateEnrFuncInDofMan(i, ei));

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
        

        answer.times( this->evaluateHeavisideGamma(lcoords.at(3), ei) ); 
        answer.times(DISC_DOF_SCALE_FAC);

    } else if ( ei && dynamic_cast< Delamination*>(ei) ){
        ModShell7Base :: computeBmatrixAt(lcoords, answer);
        std :: vector< double > efGP;
        double levelSetGP = this->evaluateLevelSet(lcoords, ei);
        ei->evaluateEnrFuncAt(efGP, lcoords, levelSetGP);
        if ( efGP[0] > 0.1 ) {
            answer.times( efGP[0]*DISC_DOF_SCALE_FAC );
        } else {
            answer.times(0.0);
        }
    } else {
        ModShell7Base :: computeBmatrixAt(lcoords, answer);
    }
}

double
ModShell7BaseXFEM :: EvaluateEnrFuncInDofMan(int dofManNum, EnrichmentItem *ei)
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
        
}


void
ModShell7BaseXFEM :: computeEnrichedNmatrixAt(const FloatArray &lcoords, FloatMatrix &answer, EnrichmentItem *ei)
{
    // Returns the displacement interpolation matrix {N} of the receiver,
    // evaluated at aGaussPoint.

    int ndofs = ModShell7Base :: giveNumberOfDofs();
    int ndofs_xm  = 3 * this->giveNumberOfDofManagers();
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
        
        answer.times( this->evaluateHeavisideGamma(lcoords.at(3), ei) ); 
        answer.times(DISC_DOF_SCALE_FAC);

    } else if ( ei && dynamic_cast< Delamination*>(ei) ) {
        ModShell7Base :: computeNmatrixAt(lcoords, answer);
        std :: vector< double > efGP;
        double levelSetGP = this->evaluateLevelSet(lcoords, ei);
        ei->evaluateEnrFuncAt(efGP, lcoords, levelSetGP);
        if ( efGP[0] > 0.1 ) {
            answer.times( efGP[0] * DISC_DOF_SCALE_FAC );
        } else {
            answer.times(0.0);
        }
    } else {
        ModShell7Base :: computeNmatrixAt(lcoords, answer);
    }
}



void
ModShell7BaseXFEM :: edgeComputeEnrichedNmatrixAt(FloatArray &lcoords, FloatMatrix &answer, EnrichmentItem *ei)
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
        answer.times( this->evaluateHeavisideGamma(lcoords.at(3), ei) ); 
        answer.times(DISC_DOF_SCALE_FAC);

    } else if ( ei && dynamic_cast< Delamination*>(ei) ) {
        ModShell7Base :: edgeComputeNmatrixAt(lcoords, answer);
        std :: vector< double > efGP;
        double levelSetGP = this->edgeEvaluateLevelSet(lcoords, ei);
        ei->evaluateEnrFuncAt(efGP, lcoords, levelSetGP);
        if ( efGP[0] > 0.1 ) {
            answer.times( efGP[0] * DISC_DOF_SCALE_FAC );
        } else {
            answer.times(0.0);
        }
    } else {
        ModShell7Base :: edgeComputeNmatrixAt(lcoords, answer);
    }

}


void
ModShell7BaseXFEM :: edgeComputeEnrichedBmatrixAt(FloatArray &lcoords, FloatMatrix &answer, EnrichmentItem *ei)
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

        answer.times( this->evaluateHeavisideGamma(lcoords.at(3), ei) ); 
        answer.times(DISC_DOF_SCALE_FAC);

    } else if ( ei && dynamic_cast< Delamination*>(ei) ){
        ModShell7Base :: edgeComputeBmatrixAt(lcoords, answer);
        std :: vector< double > efGP;
        double levelSetGP = this->evaluateLevelSet(lcoords, ei);
        ei->evaluateEnrFuncAt(efGP, lcoords, levelSetGP);
        if ( efGP[0] > 0.1 ) {
            answer.times( efGP[0]*DISC_DOF_SCALE_FAC );
        } else {
            answer.times(0.0);
        }
    } else {
        ModShell7Base :: edgeComputeBmatrixAt(lcoords, answer);
    }

}






// Delamination specific


void
ModShell7BaseXFEM :: vtkEvalUpdatedGlobalCoordinateAt(FloatArray &localCoords, int layer, FloatArray &globalCoords, TimeStep *tStep)
{
    double zeta = this->giveGlobalZcoord( localCoords ); 

    // Continuous part
    FloatArray solVec;
    this->giveUpdatedSolutionVector(solVec, tStep);
    FloatArray xc, mc; double gamc=0;
    ModShell7Base :: giveUnknownsAt(localCoords, solVec, xc, mc, gamc, tStep); 
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
	    this->computeDiscSolutionVector(eiDofIdArray, tStep, solVecD);
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
ModShell7BaseXFEM :: giveDisUnknownsAt(FloatArray &lcoords, EnrichmentItem *ei, FloatArray &solVec, FloatArray &x, FloatArray &m, double gam, TimeStep *tStep)
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
ModShell7BaseXFEM :: giveCompositeExportData(std::vector< VTKPiece > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep )
{
    vtkPieces.resize(1);
    this->giveShellExportData(vtkPieces[0], primaryVarsToExport, internalVarsToExport, cellVarsToExport, tStep );
    //this->giveCZExportData(vtkPieces[1], primaryVarsToExport, internalVarsToExport, cellVarsToExport, tStep );
    
}



void
ModShell7BaseXFEM :: giveCompositeExportData(VTKPiece &vtkPiece, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep )
{

    this->giveCZExportData(vtkPiece, primaryVarsToExport, internalVarsToExport, cellVarsToExport, tStep );
    //this->giveShellExportData(vtkPiece, primaryVarsToExport, internalVarsToExport, cellVarsToExport, tStep );
    
}

void
ModShell7BaseXFEM :: giveShellExportData(VTKPiece &vtkPiece, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep )
{
    int numSubCells = 1;
    //if ( this->crackSubdivisions[0].size() ) {
    //    numSubCells = (int)this->crackSubdivisions[0].size();
    //}
    
    int numLayers = this->layeredCS->giveNumberOfLayers();
    //int numCells = numLayers * numSubCells;
    // Go through each layer and count the number of subcells
    int numCells = 0;
    for ( int i = 0; i < numLayers; i++ ) {
        numCells += (int)(this->crackSubdivisions[i].size());
    }

    int numCellNodes  = 15; // quadratic wedge

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

        numSubCells = (int)this->crackSubdivisions[layer - 1].size();
        for ( int subCell = 1; subCell <= numSubCells; subCell++ ) {

            // Node coordinates
            if ( numSubCells == 1 ) {
                ModShell7Base :: giveFictiousNodeCoordsForExport(nodeCoords, layer); 
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

            numSubCells = (int)this->crackSubdivisions[layer - 1].size();
            for ( int subCell = 1; subCell <= numSubCells; subCell++ ) {

                if ( type == DisplacementVector ) { // compute displacement as u = x - X
                    if ( numSubCells == 1 ) {
                        ModShell7Base :: giveFictiousNodeCoordsForExport(nodeCoords, layer);
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
                    NodalRecoveryMI_recoverValues(values, layer, ( InternalStateType ) 1, tStep); // does not work well - fix
                    for ( int j = 1; j <= numCellNodes; j++ ) {
                        vtkPiece.setPrimaryVarInNode(fieldNum, nodeNum, values[j-1]);
                        nodeNum += 1;
                    }
                }

                currentCell++;
            }
        }
    }



    // Export nodal variables from internal fields
    
    vtkPiece.setNumberOfInternalVarsToExport( internalVarsToExport.giveSize(), numTotalNodes );
    for ( int fieldNum = 1; fieldNum <= internalVarsToExport.giveSize(); fieldNum++ ) {
        InternalStateType type = ( InternalStateType ) internalVarsToExport.at(fieldNum);
        nodeNum = 1;
        //this->recoverShearStress(tStep);
        //int currentCell = 1;
        for ( int layer = 1; layer <= numLayers; layer++ ) {   
            numSubCells = (int)this->crackSubdivisions[layer - 1].size();

            for ( int subCell = 1; subCell <= numSubCells; subCell++ ) {
                recoverValuesFromIP(values, layer, type, tStep);        
                
                for ( int j = 1; j <= numCellNodes; j++ ) {
                    vtkPiece.setInternalVarInNode( fieldNum, nodeNum, values[j-1] );
                    //values[j-1].printYourself();
                    //ZZNodalRecoveryMI_recoverValues(el.nodeVars[fieldNum], layer, type, tStep);          
                    nodeNum += 1;        
                }                                
            }
        }  
    }


    // Export cell variables
    FloatArray average;
    vtkPiece.setNumberOfCellVarsToExport(cellVarsToExport.giveSize(), numCells);
    for ( int i = 1; i <= cellVarsToExport.giveSize(); i++ ) {
        InternalStateType type = ( InternalStateType ) cellVarsToExport.at(i);;
        InternalStateValueType valueType =  giveInternalStateValueType(type);
        int currentCell = 1;
        for ( int layer = 1; layer <= numLayers; layer++ ) {  
            numSubCells = (int)this->crackSubdivisions[layer - 1].size();
            for ( int subCell = 1; subCell <= numSubCells; subCell++ ) {
                IntegrationRule *iRuleL = integrationRulesArray [ layer - 1 ];
                VTKXMLExportModule::computeIPAverage(average, iRuleL, this, type, tStep);
                if ( valueType == ISVT_TENSOR_S3 ) {
                    vtkPiece.setCellVar(i, currentCell, convV6ToV9Stress(average) );  
                } else {
                    vtkPiece.setCellVar(i, currentCell, average);  
                }
                currentCell += 1;
            }
        }

    }




#if 1


    // Export of XFEM related quantities
    XfemManager *xFemMan =  this->xMan;
    int nEnrIt = xFemMan->giveNumberOfEnrichmentItems();

    IntArray wedgeToTriMap;
    // For each node in the wedge take the value from the triangle node given by the map below
    wedgeToTriMap = { 1, 2, 3, 1, 2, 3, 4, 5, 6, 4, 5, 6, 1, 2, 3 };

    vtkPiece.setNumberOfInternalXFEMVarsToExport(xFemMan->vtkExportFields.giveSize(), nEnrIt, numTotalNodes);
    for ( int field = 1; field <= xFemMan->vtkExportFields.giveSize(); field++ ) {
        XFEMStateType xfemstype = ( XFEMStateType ) xFemMan->vtkExportFields [ field - 1 ];
        
        for ( int enrItIndex = 1; enrItIndex <= nEnrIt; enrItIndex++ ) {
            EnrichmentItem *ei = xFemMan->giveEnrichmentItem(enrItIndex);
            int nodeNum = 1;
            for ( int layer = 1; layer <= numLayers; layer++ ) {
                FloatMatrix localNodeCoords;

                numSubCells = (int)this->crackSubdivisions[layer - 1].size();
                for ( int subCell = 1; subCell <= numSubCells; subCell++ ) {
                    FloatArray nodeLocalXi1Coords, nodeLocalXi2Coords, nodeLocalXi3Coords, nodeLocalXi3CoordsMapped;
                    if ( numSubCells == 1) {
                        this->interpolationForExport.giveLocalNodeCoords(localNodeCoords);
                    } else {
                        giveLocalNodeCoordsForExport(nodeLocalXi1Coords, nodeLocalXi2Coords, nodeLocalXi3Coords, subCell, layer, localNodeCoords);
                    }
                    mapXi3FromLocalToShell(nodeLocalXi3CoordsMapped, nodeLocalXi3Coords, layer);
                    for ( int nodeIndx = 1; nodeIndx <= numCellNodes; nodeIndx++ ) {

                        FloatArray lCoords;
                        //lCoords.setValues(3,  nodeLocalXi1Coords.at(nodeIndx), nodeLocalXi2Coords.at(nodeIndx), nodeLocalXi3CoordsMapped.at(nodeIndx));
                        lCoords.beColumnOf(localNodeCoords, nodeIndx);
                        Node *node = this->giveNode( wedgeToTriMap.at(nodeIndx) );
                        FloatArray valueArray;
                        //const FloatArray *val = NULL;
                        if ( xfemstype == XFEMST_LevelSetPhi ) {
                            valueArray.resize(1);
                            //val = & valueArray;
                            //ei->evalLevelSetNormalInNode( valueArray.at(1), node->giveNumber() );
                            valueArray.at(1) = this->evaluateLevelSet(lCoords, ei);
                        } else if ( xfemstype == XFEMST_LevelSetGamma ) {
                            valueArray.resize(1);
                            //val = & valueArray;
                            ei->evalLevelSetTangInNode( valueArray.at(1), node->giveNumber() );
                        } else if ( xfemstype == XFEMST_NodeEnrMarker ) {
                            valueArray.resize(1);
                            //val = & valueArray;
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
ModShell7BaseXFEM :: giveFictiousNodeCoordsForExport(std::vector<FloatArray> &nodes, int layer, int subCell)
{

    // compute fictious node coords
    FloatArray nodeLocalXi1Coords, nodeLocalXi2Coords, nodeLocalXi3Coords;
    FloatMatrix localNodeCoords;
    // need to return local coordinates corresponding to the nodes of the sub triangles
    giveLocalNodeCoordsForExport(nodeLocalXi1Coords, nodeLocalXi2Coords, nodeLocalXi3Coords, subCell, layer, localNodeCoords);

    //this->interpolationForExport.giveLocalNodeCoords(localNodeCoords);

    nodes.resize(localNodeCoords.giveNumberOfColumns());
    for ( int i = 1; i <= localNodeCoords.giveNumberOfColumns(); i++ ){
        FloatArray coords, localCoords(3);
        localCoords.at(1) = nodeLocalXi1Coords.at(i);
        localCoords.at(2) = nodeLocalXi2Coords.at(i);
        localCoords.at(3) = nodeLocalXi3Coords.at(i);
        localCoords.beColumnOf(localNodeCoords,i);

        this->vtkEvalInitialGlobalCoordinateAt(localCoords, layer, coords);
        nodes[i-1].resize(3); 
        nodes[i-1] = coords;
    }

}


void 
ModShell7BaseXFEM :: giveFictiousCZNodeCoordsForExport(std::vector<FloatArray> &nodes, int layer, int subCell)
{

    // compute fictious node coords
    FloatArray nodeLocalXi1Coords, nodeLocalXi2Coords, nodeLocalXi3Coords;
    FloatMatrix localNodeCoords;
    // need to return local coordinates corresponding to the nodes of the sub triangles
    giveLocalCZNodeCoordsForExport(nodeLocalXi1Coords, nodeLocalXi2Coords, nodeLocalXi3Coords, subCell, localNodeCoords);

    nodes.resize(localNodeCoords.giveNumberOfColumns());
    for ( int i = 1; i <= localNodeCoords.giveNumberOfColumns(); i++ ){
        FloatArray coords, localCoords(3);
        localCoords.beColumnOf(localNodeCoords,i);

        this->vtkEvalInitialGlobalCZCoordinateAt(localCoords, layer, coords);
        nodes[i-1].resize(3); 
        nodes[i-1] = coords;
    }

}

void 
ModShell7BaseXFEM :: giveFictiousUpdatedNodeCoordsForExport(std::vector<FloatArray> &nodes, int layer, TimeStep *tStep, int subCell)
{
    // compute fictious node coords
    FloatArray nodeLocalXi1Coords, nodeLocalXi2Coords, nodeLocalXi3Coords;

    FloatMatrix localNodeCoords;


    if ( subCell == 0) { 
        this->interpolationForExport.giveLocalNodeCoords(localNodeCoords);
        // must get local z-coord in terms of the total thickness not layerwise
        
    } else {
        giveLocalNodeCoordsForExport(nodeLocalXi1Coords, nodeLocalXi2Coords, nodeLocalXi3Coords, subCell, layer, localNodeCoords);
    }
    nodes.resize(localNodeCoords.giveNumberOfColumns());
    for ( int i = 1; i <= localNodeCoords.giveNumberOfColumns(); i++ ){
        FloatArray coords, localCoords(3);
        localCoords.beColumnOf(localNodeCoords, i);

        // Map local layer cs to local shell cs
        double scaleFactor = 0.9999; // Will be numerically unstable with xfem if the endpoints lie at +-1 - or will it?
        double totalThickness = this->layeredCS->computeIntegralThick();
        double zMid_i = this->layeredCS->giveLayerMidZ(layer); // global z-coord
        double xiMid_i = 1.0 - 2.0 * ( totalThickness - this->layeredCS->giveMidSurfaceZcoordFromBottom() - zMid_i ) / totalThickness; // local z-coord
        double deltaxi = localCoords.at(3) * this->layeredCS->giveLayerThickness(layer) / totalThickness; // distance from layer mid
        localCoords.at(3) = xiMid_i + deltaxi * scaleFactor;

        this->vtkEvalUpdatedGlobalCoordinateAt(localCoords, layer, coords, tStep);
        nodes[i-1].resize(3); 
        nodes[i-1] = coords;
    }
}

void 
ModShell7BaseXFEM :: giveFictiousUpdatedCZNodeCoordsForExport(std::vector<FloatArray> &nodes, int interface, TimeStep *tStep, int subCell)
{
    // compute fictious node coords
    FloatArray nodeLocalXi1Coords, nodeLocalXi2Coords, nodeLocalXi3Coords;

    FloatMatrix localNodeCoords;


    if ( subCell == 0) { 
        this->interpolationForCZExport.giveLocalNodeCoords(localNodeCoords);
        // must get local z-coord in terms of the total thickness not layerwise
        
    } else {
        giveLocalCZNodeCoordsForExport(nodeLocalXi1Coords, nodeLocalXi2Coords, nodeLocalXi3Coords, subCell, localNodeCoords);
    }
    nodes.resize(localNodeCoords.giveNumberOfColumns());
    for ( int i = 1; i <= localNodeCoords.giveNumberOfColumns(); i++ ){
        FloatArray coords, localCoords(3);
        localCoords.beColumnOf(localNodeCoords, i);
        localCoords.at(3) = 1.0;
        // Map local layer cs to local shell cs
        double scaleFactor = 0.9999; // Will be numerically unstable with xfem if the endpoints lie at +-1 - or will it?
        double totalThickness = this->layeredCS->computeIntegralThick();
        double zMid_i = this->layeredCS->giveLayerMidZ(interface); // global z-coord
        double xiMid_i = 1.0 - 2.0 * ( totalThickness - this->layeredCS->giveMidSurfaceZcoordFromBottom() - zMid_i ) / totalThickness; // local z-coord
        double deltaxi = localCoords.at(3) * this->layeredCS->giveLayerThickness(interface) / totalThickness; // distance from layer mid
        localCoords.at(3) = xiMid_i + deltaxi * scaleFactor;

        this->vtkEvalUpdatedGlobalCoordinateAt(localCoords, interface, coords, tStep);
        nodes[i-1].resize(3); 
        nodes[i-1] = coords;
    }
}


void
ModShell7BaseXFEM :: mapXi3FromLocalToShell(FloatArray &answer, FloatArray &local, int layer)
{
    // Map local layer cs to local shell cs
    answer.resize(15);
    for ( int i = 1; i <= 15; i++ ){
        double scaleFactor = 0.9999; // Will be numerically unstable with xfem if the endpoints lie at +-1
        double totalThickness = this->layeredCS->computeIntegralThick();
        double zMid_i = this->layeredCS->giveLayerMidZ(layer); // global z-coord
        double xiMid_i = 1.0 - 2.0 * ( totalThickness - this->layeredCS->giveMidSurfaceZcoordFromBottom() - zMid_i ) / totalThickness; // local z-coord
        double deltaxi = local.at(i) * this->layeredCS->giveLayerThickness(layer) / totalThickness; // distance from layer mid
        local.at(i) = xiMid_i + deltaxi * scaleFactor;

        answer.at(i) = local.at(i);
    }
}



void
ModShell7BaseXFEM :: giveLocalNodeCoordsForExport(FloatArray &nodeLocalXi1Coords, FloatArray &nodeLocalXi2Coords, FloatArray &nodeLocalXi3Coords, int subCell, int layer, FloatMatrix &localNodeCoords) 
{
    // Local coords for a quadratic wedge element - coords for subtriangles
    double scale = 0.9999;
    double z = 1.0*scale;

    // Triangle coordinates
    FloatArray g1, g2, g3;
    //g1 = this->allTri[subCell-1].giveVertex(1);
    //g2 = this->allTri[subCell-1].giveVertex(2);
    //g3 = this->allTri[subCell-1].giveVertex(3);
    g1 = this->crackSubdivisions[layer-1][subCell-1].giveVertex(1);
    g2 = this->crackSubdivisions[layer-1][subCell-1].giveVertex(2);
    g3 = this->crackSubdivisions[layer-1][subCell-1].giveVertex(3);

    FloatArray gs1, gs2, gs3;
    // Move the triangle nodes slightly towards the center to avoid numerical problems - controlled by 'scale' 
    double alpha1 = scale; double alpha2 = (1.0-alpha1)*0.5; double alpha3 = alpha2;
    gs1 = alpha1*g1 + alpha2*g2 + alpha3*g3;
    gs2 = alpha2*g1 + alpha1*g2 + alpha3*g3;
    gs3 = alpha2*g1 + alpha3*g2 + alpha1*g3;

    // Local coordinates for the (scaled) triangle coordinates
    FloatArray loc1, loc2, loc3;
    this->computeLocalCoordinates(loc1, gs1);
    this->computeLocalCoordinates(loc2, gs2);
    this->computeLocalCoordinates(loc3, gs3);

    // Compute coordinates for the three mid nodes 
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
    nodeLocalXi1Coords = { a, b, c, a, b, c, d, e, f, d, e, f, a, b, c };
    
    a = loc1.at(2);
    b = loc2.at(2);
    c = loc3.at(2);
    d = loc12.at(2);
    e = loc23.at(2);
    f = loc31.at(2);
    nodeLocalXi2Coords = { a, b, c, a, b, c, d, e, f, d, e, f, a, b, c };

    nodeLocalXi3Coords = { -z, -z, -z, z, z, z, -z, -z, -z, z, z, z, 0., 0., 0. };
    
    FloatMatrix localNodeCoordsT;
    localNodeCoordsT.resize(15,3);
    localNodeCoordsT.setColumn(nodeLocalXi1Coords,1);
    localNodeCoordsT.setColumn(nodeLocalXi2Coords,2);
    localNodeCoordsT.setColumn(nodeLocalXi3Coords,3);
    localNodeCoords.beTranspositionOf(localNodeCoordsT);

}


void
ModShell7BaseXFEM :: giveLocalCZNodeCoordsForExport(FloatArray &nodeLocalXi1Coords, FloatArray &nodeLocalXi2Coords, FloatArray &nodeLocalXi3Coords, int subCell, FloatMatrix &localNodeCoords) 
{
    // Local coords for a quadratic triangle element - coords for subtriangles
    double scale = 0.999;
    //double z = 1.0*scale;

    // Triangle coordinates
    FloatArray g1, g2, g3;
    g1 = this->allTri[subCell-1].giveVertex(1);
    g2 = this->allTri[subCell-1].giveVertex(2);
    g3 = this->allTri[subCell-1].giveVertex(3);

    FloatArray gs1, gs2, gs3;
    // Move the triangle nodes slightly towards the center to avoid numerical problems - controlled by 'scale' 
    double alpha1 = scale; double alpha2 = (1.0-alpha1)*0.5; double alpha3 = alpha2;
    gs1 = alpha1*g1 + alpha2*g2 + alpha3*g3;
    gs2 = alpha2*g1 + alpha1*g2 + alpha3*g3;
    gs3 = alpha2*g1 + alpha3*g2 + alpha1*g3;

    // Local coordinates for the (scaled) triangle coordinates
    FloatArray loc1, loc2, loc3;
    this->computeLocalCoordinates(loc1, gs1);
    this->computeLocalCoordinates(loc2, gs2);
    this->computeLocalCoordinates(loc3, gs3);

    // Compute coordinates for the three mid nodes 
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
    nodeLocalXi1Coords = { a, b, c, d, e, f };
    
    a = loc1.at(2);
    b = loc2.at(2);
    c = loc3.at(2);
    d = loc12.at(2);
    e = loc23.at(2);
    f = loc31.at(2);
    nodeLocalXi2Coords = { a, b, c, d, e, f };

    nodeLocalXi3Coords = { 0., 0., 0., 0., 0., 0. };
    
    FloatMatrix localNodeCoordsT;
    localNodeCoordsT.resize(6,3);
    localNodeCoordsT.setColumn(nodeLocalXi1Coords,1);
    localNodeCoordsT.setColumn(nodeLocalXi2Coords,2);
    localNodeCoordsT.setColumn(nodeLocalXi3Coords,3);
    localNodeCoords.beTranspositionOf(localNodeCoordsT);

}

void
ModShell7BaseXFEM :: giveCZExportData(VTKPiece &vtkPiece, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep )
{

    int numSubCells = 1;
    if ( this->allTri.size() ) {
        numSubCells = (int)this->allTri.size();
    }
    
    int numInterfaces = this->layeredCS->giveNumberOfLayers()-1;
    int numCells = numInterfaces * numSubCells;
    
    int numCellNodes = 6;   // quadratic triangle

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
    for ( int layer = 1; layer <= numInterfaces; layer++ ) {
        
        for ( int subCell = 1; subCell <= numSubCells; subCell++ ) {

            // Node coordinates
            if ( numSubCells == 1 ) {
                ModShell7Base :: giveFictiousCZNodeCoordsForExport(nodeCoords, layer); 
            } else {
                this->giveFictiousCZNodeCoordsForExport(nodeCoords, layer, subCell);       
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
            vtkPiece.setCellType(currentCell, 22); // Quadratic triangle

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
        for ( int layer = 1; layer <= numInterfaces; layer++ ) {            
            
            for ( int subCell = 1; subCell <= numSubCells; subCell++ ) {

                if ( type == DisplacementVector ) { // compute displacement as u = x - X
                    if ( numSubCells == 1 ) {
                        ModShell7Base :: giveFictiousCZNodeCoordsForExport(nodeCoords, layer);
                        this->giveFictiousUpdatedCZNodeCoordsForExport(updatedNodeCoords, layer, tStep, 0);
                    } else {
                        this->giveFictiousCZNodeCoordsForExport(nodeCoords, layer, subCell);
                        this->giveFictiousUpdatedCZNodeCoordsForExport(updatedNodeCoords, layer, tStep, subCell);
                    }
                    for ( int j = 1; j <= numCellNodes; j++ ) {
                        u = updatedNodeCoords[j-1];
                        u.subtract(nodeCoords[j-1]);
                        vtkPiece.setPrimaryVarInNode(fieldNum, nodeNum, u);
                        nodeNum += 1;        
                    }

                } else {
                    NodalRecoveryMI_recoverValues(values, layer, ( InternalStateType ) 1, tStep); // does not work well - fix
                    for ( int j = 1; j <= numCellNodes; j++ ) {
                        vtkPiece.setPrimaryVarInNode(fieldNum, nodeNum, values[j-1]);
                        nodeNum += 1;
                    }
                }

                currentCell++;
            }
        }
    }



    // Export nodal variables from internal fields
    
    vtkPiece.setNumberOfInternalVarsToExport( internalVarsToExport.giveSize(), numTotalNodes );
    for ( int fieldNum = 1; fieldNum <= internalVarsToExport.giveSize(); fieldNum++ ) {
        InternalStateType type = ( InternalStateType ) internalVarsToExport.at(fieldNum);
        nodeNum = 1;
        //this->recoverShearStress(tStep);
        //int currentCell = 1;
        for ( int layer = 1; layer <= numInterfaces; layer++ ) {            
            for ( int subCell = 1; subCell <= numSubCells; subCell++ ) {
                this->recoverValuesFromCZIP(values, layer, type, tStep);        
                
                for ( int j = 1; j <= numCellNodes; j++ ) {
                    vtkPiece.setInternalVarInNode( fieldNum, nodeNum, values[j-1] );
                    //ZZNodalRecoveryMI_recoverValues(el.nodeVars[fieldNum], layer, type, tStep);          
                    nodeNum += 1;        
                }                                
            }
        }  
    }


    // Export cell variables
    FloatArray average;
    vtkPiece.setNumberOfCellVarsToExport(cellVarsToExport.giveSize(), numCells);
    for ( int i = 1; i <= cellVarsToExport.giveSize(); i++ ) {
        InternalStateType type = ( InternalStateType ) cellVarsToExport.at(i);;
        InternalStateValueType valueType =  giveInternalStateValueType(type);
        int currentCell = 1;
        for ( int layer = 1; layer <= numInterfaces; layer++ ) {  
            for ( int subCell = 1; subCell <= numSubCells; subCell++ ) {
                if ( type == IST_CrossSectionNumber ) {
                    average = FloatArray{ -double(layer) }; // Set a negative number for interfaces
                } else {
                    IntegrationRule *iRuleL = integrationRulesArray [ layer - 1 ];
                    VTKXMLExportModule::computeIPAverage(average, iRuleL, this, type, tStep);
                }
                if ( valueType == ISVT_TENSOR_S3 ) {
                    vtkPiece.setCellVar(i, currentCell, convV6ToV9Stress(average) );  
                } else {
                    vtkPiece.setCellVar(i, currentCell, average);  
                }
                currentCell += 1;
            }
        }

    }




#if 0


    // Export of XFEM related quantities
    XfemManager *xFemMan =  this->xMan;
    int nEnrIt = xFemMan->giveNumberOfEnrichmentItems();

    IntArray wedgeToTriMap;
    // For each node in the wedge take the value from the triangle node given by the map below
    wedgeToTriMap.setValues(15, 1, 2, 3, 1, 2, 3, 4, 5, 6, 4, 5, 6, 1, 2, 3 );

    vtkPiece.setNumberOfInternalXFEMVarsToExport(xFemMan->vtkExportFields.giveSize(), nEnrIt, numTotalNodes);
    for ( int field = 1; field <= xFemMan->vtkExportFields.giveSize(); field++ ) {
        XFEMStateType xfemstype = ( XFEMStateType ) xFemMan->vtkExportFields [ field - 1 ];
        
        for ( int enrItIndex = 1; enrItIndex <= nEnrIt; enrItIndex++ ) {
            EnrichmentItem *ei = xFemMan->giveEnrichmentItem(enrItIndex);
            int nodeNum = 1;
            for ( int layer = 1; layer <= numInterfaces; layer++ ) {
                FloatMatrix localNodeCoords;
                
                for ( int subCell = 1; subCell <= numSubCells; subCell++ ) {
                    FloatArray nodeLocalXi1Coords, nodeLocalXi2Coords, nodeLocalXi3Coords, nodeLocalXi3CoordsMapped;
                    if ( numSubCells == 1) {
                        this->interpolationForExport.giveLocalNodeCoords(localNodeCoords);
                    } else {
                        giveLocalNodeCoordsForExport(nodeLocalXi1Coords, nodeLocalXi2Coords, nodeLocalXi3Coords, subCell, localNodeCoords);
                    }
                    mapXi3FromLocalToShell(nodeLocalXi3CoordsMapped, nodeLocalXi3Coords, layer);
                    for ( int nodeIndx = 1; nodeIndx <= numCellNodes; nodeIndx++ ) {

                        FloatArray lCoords;
                        //lCoords.setValues(3,  nodeLocalXi1Coords.at(nodeIndx), nodeLocalXi2Coords.at(nodeIndx), nodeLocalXi3CoordsMapped.at(nodeIndx));
                        lCoords.beColumnOf(localNodeCoords, nodeIndx);
                        //Node *node = this->giveNode( wedgeToTriMap.at(nodeIndx) );
                        Node *node = this->giveNode( nodeIndx );
                        FloatArray valueArray;
                        const FloatArray *val = NULL;
                        if ( xfemstype == XFEMST_LevelSetPhi ) {
                            valueArray.resize(1);
                            val = & valueArray;
                            valueArray.at(1) = this->evaluateLevelSet(lCoords, ei);
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
ModShell7BaseXFEM :: recoverValuesFromCZIP(std::vector<FloatArray> &recoveredValues, int interfce, InternalStateType type, TimeStep *tStep)
{
    // recover nodal values by coosing the ip closest to the node

    // composite element interpolator
    FloatMatrix localNodeCoords;
    //this->interpolationForExport.giveLocalNodeCoords(localNodeCoords);
    this->interpolationForCZExport.giveLocalNodeCoords(localNodeCoords);
    int numNodes = localNodeCoords.giveNumberOfColumns();
    recoveredValues.resize(numNodes);
    
    IntegrationRule *iRule = this->czIntegrationRulesArray [ interfce - 1 ];
    IntegrationPoint *ip;

    // Find closest ip to the nodes
    IntArray closestIPArray(numNodes);
    FloatArray nodeCoords, ipCoords, ipValues;

    for ( int i = 1; i <= numNodes; i++ ) {
        nodeCoords.beColumnOf(localNodeCoords, i);
        double distOld = 3.0; // should not be larger
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            ip = iRule->getIntegrationPoint(j);
            ipCoords = *ip->giveNaturalCoordinates();
            double dist = nodeCoords.distance(ipCoords);
            if ( dist < distOld ) {
                closestIPArray.at(i) = j;
                distOld = dist;
            }
        }
    }

   InternalStateValueType valueType =  giveInternalStateValueType(type);

    // recover ip values
    for ( int i = 1; i <= numNodes; i++ ) {
        
        if ( this->layeredCS->giveInterfaceMaterial(interfce) ) {
            ip = iRule->getIntegrationPoint( closestIPArray.at(i) );
            this->layeredCS->giveInterfaceMaterial(interfce)->giveIPValue(ipValues, ip, type, tStep);
        } else {
            ipValues.resize(0);
        }

        if ( ipValues.giveSize() == 0 && type == IST_AbaqusStateVector) {
            recoveredValues[i-1].resize(23);
            recoveredValues[i-1].zero();        
        } else if ( ipValues.giveSize() == 0 ) {
            recoveredValues[i-1].resize(giveInternalStateTypeSize(valueType));
            recoveredValues[i-1].zero();
        
        } else if ( valueType == ISVT_TENSOR_S3 ) {
            recoveredValues[i-1].resize(9);
            recoveredValues[i-1] = convV6ToV9Stress(ipValues);

        } else {
            recoveredValues[i-1] = ipValues;
        }
    }

}



void
ModShell7BaseXFEM :: computeTripleProduct(FloatMatrix &answer, const FloatMatrix &a, const FloatMatrix &b, const FloatMatrix &c)
{
    // Computes the product a^T*b*c
    FloatMatrix temp;
    temp.beTProductOf(a, b);
    answer.beProductOf(temp, c);
}


} // end namespace oofem
