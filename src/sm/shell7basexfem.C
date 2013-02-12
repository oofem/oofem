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
    
    this->delaminationXiCoordList.clear();
    this->gpDelaminationGroupList.clear();
}


void
Shell7BaseXFEM :: setupDelaminationXiCoordList()
{
    if ( this->delaminationXiCoordList.size()==0 ) {
    // Stores a paired list with the EnrichmentDomain# and the corresponding xi-coord of the delamination.
    xMan =  this->giveDomain()->giveEngngModel()->giveXfemManager(1);
    int numEI = xMan->giveNumberOfEnrichmentItems();
    for ( int i = 1; i <= numEI; i++ ) {
        Delamination *dei =  dynamic_cast< Delamination * >( xMan->giveEnrichmentItem(i) );
        if ( dei ) {
            int numED = dei->giveNumberOfEnrichmentDomains(); // numEnrDomains max possible number
            int pos = 1;
            for ( int j = 1; j <= numED; j++ ) {
                if( dei->isElementEnriched(this->element) ) { // should check particular ED
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

