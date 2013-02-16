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

namespace oofem {
    

Shell7BaseXFEM :: Shell7BaseXFEM(int n, Domain *aDomain) : Shell7Base(n, aDomain), XfemElementInterface(this) 
{
    //xMan =  this->giveDomain()->giveEngngModel()->giveXfemManager(1);
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


void 
Shell7BaseXFEM :: giveNumberOFEnrichmentDofsForDofMan( int answer)
{
    XfemManager *xMan = this->giveDomain()->giveXfemManager(1);
    for ( int i =1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) {
        EnrichmentItem *ei = xMan->giveEnrichmentItem(i);
        for ( int j = 1; j <= ei->giveNumberOfEnrichmentDomains(); j++ ) {
            bool val = ei->isElementEnrichedByEnrichmentDomain(this, j);
        }
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
Shell7BaseXFEM :: evalCovarBaseVectorsAt(GaussPoint *gp, FloatArray &g1, FloatArray &g2, FloatArray &g3, FloatArray &genEpsC)
{
    // Continuous part
    FloatArray g1c, g2c, g3c;
    Shell7Base :: evalCovarBaseVectorsAt(gp, g1c, g2c, g3c, genEpsC);

    // Discontinuous part
    xMan =  this->giveDomain()->giveEngngModel()->giveXfemManager(1);
    for ( int i = 1; i <= xMan->giveNumberOfEnrichmentItems(); i++ ) {
        Delamination *dei =  dynamic_cast< Delamination * >( xMan->giveEnrichmentItem(i) ); // should check success
        for ( int j = 1; j <= dei->giveNumberOfEnrichmentDomains(); j++ ) {
            if ( dei->isElementEnrichedByEnrichmentDomain(this, j) ) {
            // FloatArray g1d, g2d, g3d, dGenEps;
            // discComputeGeneralizedStrainVector(dGenEps, dSolVec, B11, B22, B32);
            // discGiveGeneralizedStrainComponents(dGenEps, dxdxi1, dxdxi2, dmdxi1, dmdxi2, m);
            // discEvalCovarBaseVectorsAt(gp, g1d, g2d, g3d, dGenEps);
            // add contribution g1d += g1d_temp  
            }
        }
    }
}



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
Shell7BaseXFEM :: discGiveUpdatedSolutionVector(FloatArray &answer, TimeStep *tStep)
{
    // Returns the solution vector of discontinuous dofs dx_d & dm_d
    //this->giveInitialSolutionVector(answer);
    FloatArray temp;
    this->computeVectorOf(EID_MomentumBalance, VM_Total, tStep, temp);
    answer.assemble( temp, this->giveOrdering(AllInv) );
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

    xMan = this->giveDomain()->giveEngngModel()->giveXfemManager(1);
    int numEI = xMan->giveNumberOfEnrichmentItems();
    Element *e = this->element;
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

