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

namespace oofem {
    IntArray Shell7BaseXFEM :: dofId_Midplane(3);
	IntArray Shell7BaseXFEM :: dofId_Director(3);
    IntArray Shell7BaseXFEM :: dofId_InhomStrain(1); 
	bool Shell7BaseXFEM :: __initializedFieldDofId = Shell7BaseXFEM :: initDofId();
Shell7BaseXFEM :: Shell7BaseXFEM(int n, Domain *aDomain) : Shell7Base(n, aDomain), XfemElementInterface(this) 
{
    //this->xMan =  this->giveDomain()->giveXfemManager(1);
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
    XfemManager *xMan =  this->giveDomain()->giveXfemManager(1);
    for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) { // Only one is supported at the moment
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
Shell7BaseXFEM :: evalCovarBaseVectorsAt(GaussPoint *gp, FloatArray &g1, FloatArray &g2, FloatArray &g3, FloatArray &genEpsC)
{
    // Continuous part
    FloatArray g1c, g2c, g3c;
    Shell7Base :: evalCovarBaseVectorsAt(gp, g1, g2, g3, genEpsC);

    // Discontinuous part - ///@todo naive imlementation should be changed
    TimeStep *tStep = this->giveDomain()->giveEngngModel()->giveCurrentStep();
    xMan =  this->giveDomain()->giveXfemManager(1);
    for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) { // Only one is supported at the moment
        Delamination *dei =  dynamic_cast< Delamination * >( xMan->giveEnrichmentItem(i) ); // should check success
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
    temp_computeVectorOf(eiDofIdArray, VM_Total, tStep, answer);
    //eiDofIdArray.printYourself();
    //answer.printYourself();
}



int
Shell7BaseXFEM :: giveNumberOfDofs()
{
    // Continuous part
    int nDofs = Shell7Base ::giveNumberOfDofs();

    // Discontinuous part
    XfemManager *xMan =  this->giveDomain()->giveXfemManager(1);

    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        DofManager *dMan = this->giveDofManager(i);
        for ( int j = 1; j <= xMan->giveNumberOfEnrichmentItems(); j++ ) {
            EnrichmentItem *ei = xMan->giveEnrichmentItem(j);
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

    // compute number of xDofs
    xMan =  this->giveDomain()->giveXfemManager(1);
//    EnrichmentItem *ei = xMan->giveEnrichmentItem(1);
//    int nXdofs = ei->giveNumberOfEnrichmentDomains() * ei->giveEnrichesDofsWithIdArray()->giveSize() * this->giveNumberOfDofManagers(); 
    

    answer.resize( this->giveNumberOfDofs() );
    answer.zero();
    
    FloatArray solVec;

    // Continuous part
    this->giveUpdatedSolutionVector(solVec, tStep);
    this->computeSectionalForces(answer, tStep, solVec, useUpdatedGpRecord);

    // Disccontinuous part
    
    for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) { // Only one is supported at the moment
        Delamination *dei =  dynamic_cast< Delamination * >( xMan->giveEnrichmentItem(i) ); // should check success

        for ( int j = 1; j <= dei->giveNumberOfEnrichmentDomains(); j++ ) {
            if ( dei->isElementEnrichedByEnrichmentDomain(this, j) ) {
                IntArray eiDofIdArray;
                dei->giveEIDofIdArray(eiDofIdArray, j);
                this->discGiveUpdatedSolutionVector(solVec, eiDofIdArray, tStep);
                double xi0 = dei->enrichmentDomainXiCoords.at(j);
                this->discComputeSectionalForces(answer, tStep, solVec, useUpdatedGpRecord, xi0, dei, j);                      

            }
        }
    }
 
}


void 
Shell7BaseXFEM :: computeOrderingArray( IntArray &orderingArray, IntArray &activeDofsArray,  EnrichmentItem *ei, int enrichmentDomainNumber, SolutionField field)
{
    // Routine to extract vector given an array of dofid items
    // If a certain dofId does not exist a zero is used as value
    IntArray eiDofIdArray;
    ei->giveEIDofIdArray(eiDofIdArray, enrichmentDomainNumber);

    IntArray ordering_cont = this->giveOrdering(field);
    IntArray fieldDofId    = this->giveFieldDofId(field);

    IntArray ordering_temp, temp;
    ordering_temp.resize(ordering_cont.giveSize());
    temp.resize(ordering_cont.giveSize());

    int k = 0;
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        DofManager *dMan = this->giveDofManager(i);
        int pos = (i-1)*fieldDofId.giveSize();
        int pos2 = (i-1)*eiDofIdArray.giveSize();
        if ( ei->isDofManEnrichedByEnrichmentDomain(dMan,enrichmentDomainNumber) ){

            for (int j = 1; j <= fieldDofId.giveSize(); j++ ) {
                if ( dMan->hasDofID( (DofIDItem) fieldDofId.at(j) ) ) {
                    k++;
                    ordering_temp.at(k) = ordering_cont.at(k);
                    temp.at(k) = j + pos;
                }
            }
        }
    }
      
    IntArray ordering; orderingArray.resize(k), activeDofsArray.resize(k);
    ///@todo will not work if there are several ei
    int shift = Shell7Base :: giveNumberOfDofs(); 
    for ( int i = 1; i <= k; i++ ) {
        orderingArray.at(i) = ordering_temp.at(i) + shift;
        activeDofsArray.at(i) = temp.at(i);
    }

}


void
Shell7BaseXFEM :: discComputeSectionalForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, int useUpdatedGpRecord, 
                    double xi0, EnrichmentItem *ei, int enrichmentDomainNumber)
//
{
    FloatMatrix B;
    FloatArray BtF, f, genEps;
    //answer.resize( this->giveNumberOfDofs() );
    //answer.zero();

    
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
    IntArray ordering_phibar, ordering_m, ordering_gam;
    IntArray activeDofs_phibar, activeDofs_m, activeDofs_gam;

    computeOrderingArray(ordering_phibar, activeDofs_phibar, ei, enrichmentDomainNumber, Midplane);
    computeOrderingArray(ordering_m, activeDofs_m, ei, enrichmentDomainNumber, Director);
    computeOrderingArray(ordering_gam, activeDofs_gam, ei, enrichmentDomainNumber, InhomStrain);

/*
    ordering_phibar.printYourself();
    ordering_m.printYourself();
    ordering_gam.printYourself();

    activeDofs_phibar.printYourself();
    activeDofs_m.printYourself();
    activeDofs_gam.printYourself();
    */


    FloatArray f1Ass, f2Ass, f3Ass;
    f1Ass.beSubArrayOf(f1,activeDofs_phibar);
    f2Ass.beSubArrayOf(f2,activeDofs_m);
    f3Ass.beSubArrayOf(f3,activeDofs_gam);
    answer.assemble(f1Ass, ordering_phibar);
    answer.assemble(f2Ass, ordering_m);
    answer.assemble(f3Ass, ordering_gam);

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
    xMan =  this->giveDomain()->giveXfemManager(1); // When does one need several xfemman?
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
    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * >( this->giveCrossSection() );
    int numberOfLayers = layeredCS->giveNumberOfLayers();  

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
    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * > (this->element->giveCrossSection());
    
    std::list< std::pair<int, double> >::const_iterator iter;
    iter = this->delaminationXiCoordList.begin(); 

    if ( dGroup == 0 ) {
        xiBottom = - layeredCS->giveMidSurfaceXiCoordFromBottom();
        xiTop = (*iter).second;
    } else if (dGroup == nDelam) {
        std::advance(iter, dGroup-1);
        xiBottom = (*iter).second;
        xiTop    = -layeredCS->giveMidSurfaceXiCoordFromBottom() + 2.0;
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

