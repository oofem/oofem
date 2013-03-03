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

#include "xfemmanager.h"
#include "flotmtrx.h"
#include "enrichmentitem.h"
#include "element.h"
#include "enrichmentfunction.h"
#include "cltypes.h"
#include "conTable.h"
#include "oofem_limits.h"
#include "usrdefsub.h"

namespace oofem {
EnrichmentItem :: EnrichmentItem(int n, XfemManager *xm, Domain *aDomain) : FEMComponent(n, aDomain)
{
    xmanager = xm;
    enrichmentFunction = 0;

    //new JB
    this->enrichmentFunctionList = new AList< EnrichmentFunction >(0);
    this->enrichementDomainList = new AList< BasicGeometry >(0); // Should be an enrichmentdomain type or something similar
    this->enrDomainList = new AList< EnrichmentDomain >(0); 
    numberOfEnrichmentFunctions = 0;    
    numberOfEnrichmentDomains = 0;  
    this->enrichesDofsWithIdArray = new IntArray;
}

EnrichmentItem :: ~EnrichmentItem()
{
    delete this->enrichesDofsWithIdArray;
}


// remove - should ask for specific geom. can be multiple?
BasicGeometry *EnrichmentItem :: giveGeometry()
{
    return this->enrichementDomainList->at(1);
}

BasicGeometry *EnrichmentItem :: giveGeometry(int n)
// Returns the n-th geometry.
{
    if ( this->enrichementDomainList->includes(n) ) {
        return this->enrichementDomainList->at(n);
    } else {
        OOFEM_ERROR2("giveGeometry: undefined geometry (%d)", n);
    }

    return NULL;
}



EnrichmentFunction *EnrichmentItem :: giveEnrichmentFunction(int n)
// Returns the n-th geometry.
{
    if ( enrichmentFunctionList->includes(n) ) {
        return enrichmentFunctionList->at(n);
    } else {
        OOFEM_ERROR2("giveEnrichmentFunction: undefined enrichment function (%d)", n);
    }

    return NULL;
}



bool EnrichmentItem :: isElementEnriched(const Element *element) 
{
    for ( int i = 1; i <= element->giveNumberOfDofManagers(); i++ ) {
        if ( this->isDofManEnriched( element->giveDofManager(i) ) ) {
            return true;
        }
    }
    return false;
}

bool EnrichmentItem :: isElementEnrichedByEnrichmentDomain(const Element *element, int edNumber) 
{
    for ( int i = 1; i <= element->giveNumberOfDofManagers(); i++ ) {
        DofManager *dMan = element->giveDofManager(i);
        if ( isDofManEnrichedByEnrichmentDomain(dMan, edNumber) ){
            return true;
        }
    }
    return false;
}


bool EnrichmentItem :: isDofManEnriched(const DofManager *dMan)
{
    for ( int i = 1; i <= this->enrichementDomainList->giveSize() ; i++ ) {
        if ( isDofManEnrichedByEnrichmentDomain(dMan, i) ){
            return true;
        }
    }
    return false;
}

bool EnrichmentItem :: isDofManEnrichedByEnrichmentDomain(const DofManager *dMan, int edNumber)
{
    EnrichmentDomain *ed = this->enrDomainList->at(edNumber);
    return ed->isDofManagerEnriched(dMan);
}





// Spatial query metods
#if 1
bool EnrichmentItem :: isOutside(BasicGeometry *bg)
{
    return this->giveGeometry()->isOutside(bg);
}

void EnrichmentItem :: computeIntersectionPoints(AList< FloatArray > *intersectionPoints, Element *element)
{
    this->giveGeometry()->computeIntersectionPoints(element, intersectionPoints);
}

int EnrichmentItem :: computeNumberOfIntersectionPoints(Element *element)
{
    return this->giveGeometry()->computeNumberOfIntersectionPoints(element);
}
#endif





IRResultType EnrichmentItem :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro

    this->geometry = 0;
    this->enrichmentFunction = 0;
    
    IR_GIVE_FIELD(ir, this->enrichmentDomainNumbers, IFT_EnrichmentItem_enrichmentdomains, "enrichmentdomains"); // Macro
    this->numberOfEnrichmentDomains = this->enrichmentDomainNumbers.giveSize();
    //int test = this->enrichmentDomainNumbers.maximum();


    IR_GIVE_OPTIONAL_FIELD(ir, enrichmentFunction, IFT_EnrichmentItem_enrichmentFunctionNr, "enrichmentfunction"); // Macro
    
    return IRRT_OK;
}


IRResultType Inclusion :: initializeFrom(InputRecord *ir)
{
    EnrichmentItem :: initializeFrom(ir);
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro
    int material = 0;
    IR_GIVE_FIELD(ir, material, IFT_EnrichmentItem_materialNr, "material"); // Macro
    this->mat = this->giveDomain()->giveMaterial(material);
    //this->numberOfEnrichmentFunctions = 1;
    //IR_GIVE_OPTIONAL_FIELD(ir, numberOfEnrichmentFunctions, IFT_XfemManager_numberOfEnrichmentFunctions, "numberofenrichmentfunctions"); // Macro
    return IRRT_OK;
}



int EnrichmentItem :: instanciateYourself(DataReader *dr)
{
    const char *__proc = "instanciateYourself"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro
    int i;
    std :: string name;
    EnrichmentItem *ei;
    EnrichmentFunction *ef;
    EnrichmentDomain *ed;
    InputRecord *mir;


    this->enrichmentFunctionList->growTo(numberOfEnrichmentFunctions);
    for ( i = 1; i <= this->numberOfEnrichmentFunctions; i++ ) {
        mir = dr->giveInputRecord(DataReader :: IR_enrichFuncRec, i);
        result = mir->giveRecordKeywordField(name);

        if ( result != IRRT_OK ) {
            IR_IOERR(giveClassName(), __proc, IFT_RecordIDField, "", mir, result);
        }
        
        ef = CreateUsrDefEnrichmentFunction( name.c_str(), i, this->xmanager->giveDomain() );
        if ( ef == NULL ) {
            OOFEM_ERROR2( "EnrichmentItem::instanciateYourself: unknown enrichment function (%s)", name.c_str() );
        }

        enrichmentFunctionList->put(i, ef);
        ef->initializeFrom(mir);
        
    }

    
    enrichementDomainList->growTo(numberOfEnrichmentDomains);
    enrDomainList->growTo(numberOfEnrichmentDomains);

    for ( i = 1; i <= numberOfEnrichmentDomains; i++ ) {
        mir = dr->giveInputRecord(DataReader :: IR_geoRec, i);
        result = mir->giveRecordKeywordField(name);
        if ( result != IRRT_OK ) {
            IR_IOERR(giveClassName(), __proc, IFT_RecordIDField, "", mir, result);
        }
#if 0
        ge = CreateUsrDefGeometry( name.c_str() );
        enrichementDomainList->put(i, ge);
        ge->initializeFrom(mir);
        
        if ( ge == NULL ) {
            OOFEM_ERROR2( "EnrichmentItem::instanciateYourself: unknown geometry (%s)", name.c_str() );
        }
#endif

        // new

        ed = CreateUsrDefEnrichmentDomain( name.c_str() ); 
        this->enrDomainList->put(i, ed);
        ed->initializeFrom(mir);
    }

#ifdef VERBOSE
    VERBOSE_PRINT0("Instanciated enrichment items ", nnode)
#endif
    return 1;
}


void
EnrichmentItem :: giveEIDofIdArray(IntArray &answer, int enrichmentDomainNumber)
{
    // Computes an array containing the dofId's that should be created as new dofs.
    IntArray *enrichesDofsWithIdArray = this->giveEnrichesDofsWithIdArray();
    
    int eiEnrSize = enrichesDofsWithIdArray->giveSize();
    answer.resize(eiEnrSize);
    int xDofAllocSize = eiEnrSize * this->giveNumberOfEnrDofs(); // number of new dof id's the ei will allocate
    for ( int i = 1; i <= eiEnrSize; i++ ) {
        answer.at(i) = this->giveStartOfDofIdPool() + (enrichmentDomainNumber-1)*xDofAllocSize + i; 
    }
   
}


void
EnrichmentItem :: computeDofIdArray(IntArray &answer, DofManager *dMan, int enrichmentDomainNumber)
{
    // Computes an array containing the dofId's that should be created as new dofs.
    IntArray *enrichesDofsWithIdArray = this->giveEnrichesDofsWithIdArray();
    
    int eiEnrSize = enrichesDofsWithIdArray->giveSize();
    
    // Go through the list of dofs that the EI supports and compare with the available dofs in the dofMan. 
    // Store matches in dofMask
    IntArray dofMask(eiEnrSize); dofMask.zero();
    int count = 0; 
    for ( int i = 1; i <= eiEnrSize; i++ ) {
        if ( dMan->hasDofID( (DofIDItem) enrichesDofsWithIdArray->at(i) ) ) {
            count++;
            dofMask.at(count) = i;
        }
    }

    answer.resize(count);
    int xDofAllocSize = eiEnrSize * this->giveNumberOfEnrDofs(); // number of new dof id's the ei will allocate
    for ( int i = 1; i <= count; i++ ) {
        answer.at(i) = this->giveStartOfDofIdPool() + (enrichmentDomainNumber-1)*xDofAllocSize + dofMask.at(i) ;
    }

}

int 
EnrichmentItem :: giveNumberOfEnrDofs() 
{ 
    // returns the array of dofs a particular EI s
    int temp=0;
    for ( int i = 1; i <= this->giveNumberOfEnrichmentfunctions(); i++ ) { 
        EnrichmentFunction *ef = this->giveEnrichmentFunction(i);
        // This is per regular dof
        temp += ef->giveNumberOfDofs(); // = number of functions associated with a particular enrichment function, e.g. 4 for branch function.

    }
    return temp; 
} 
    


/**
 * DELAMINATION
 */
Delamination :: Delamination(int n, XfemManager *xm, Domain *aDomain) : EnrichmentItem(n, xm, aDomain)
{ 
    //this->enrichesDofsWithIdArray->setValues(7, D_u, D_v, D_w, W_u, W_v, W_w, Gamma);
    this->enrichesDofsWithIdArray->setValues(6, D_u, D_v, D_w, W_u, W_v, W_w);
}


IRResultType Delamination :: initializeFrom(InputRecord *ir)
{
    this->numberOfEnrichmentFunctions = 1; // must be set before EnrichmentItem :: initializeFrom(ir) is called
    EnrichmentItem :: initializeFrom(ir);
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro
    int material = 0;
    

    IR_GIVE_FIELD(ir, this->enrichmentDomainXiCoords, IFT_EnrichmentItem_Delamination_delaminationXiCoords, "delaminationxicoords"); // Macro
    if ( this->numberOfEnrichmentDomains != this->enrichmentDomainXiCoords.giveSize() ) {
        OOFEM_ERROR3( "EnrichmentItem :: initializeFrom: size of enrichmentDomainXiCoords (%i) differs from numberOfEnrichmentDomains (%i)", 
           this->enrichmentDomainXiCoords.giveSize(), this->numberOfEnrichmentDomains );
    }

    

    //write an instanciate method
    
    return IRRT_OK;
}




    


double 
Delamination :: giveDelaminationZCoord(int n, Element *element) 
{ 
    AList<double> *xiCoordList;
    int nDelam = this->giveNumberOfEnrichmentDomains(); // max possible number
    int pos = 1;
    for ( int i = 1; i <= nDelam; i++ ) {
        if( this->isElementEnriched(element) ) {
            //xiCoordList.put( pos, this->delaminationZCoords.at(i) );
            pos++;
        } 
    }
    return 0.;;//
    

    
}; 

// Remove!
int
Delamination :: giveDelaminationGroupAt(double zeta) 
{
    //double zRef = Shell7Base :: giveLocalZetaCoord(gp);
    int nDelam = this->giveNumberOfEnrichmentDomains();
    for ( int j = 1; j <= nDelam; j++ ) {
        //double zDelam = this->giveDelaminationZCoord(j);
        double zDelam = 0.;
        if ( zeta  < zDelam ) { //belong to the delamination group just below delamination #j. How to deal with poins that lie onthe boundary?
            return j;            
        }

    }
    return nDelam + 1;
            
}

double 
Delamination :: heaviside(double xi, double xi0)
{
    if( xi < xi0 ) {
        return 0.0;
    } else {
        return 1.0;
    }
}


double 
Delamination :: giveDelaminationGroupMidZ(int dGroup, Element *e)
{
    double zTop=0., zBottom=0.;
    this->giveDelaminationGroupZLimits(dGroup, zTop, zBottom, e);
    return 0.5 * ( zTop + zBottom );
}


double 
Delamination :: giveDelaminationGroupThickness(int dGroup, Element *e)
{
    double zTop, zBottom;
    this->giveDelaminationGroupZLimits(dGroup, zTop, zBottom, e);
    return zTop - zBottom;
}

void 
Delamination :: giveDelaminationGroupZLimits(int &dGroup, double &zTop, double &zBottom, Element *e)
{
    int nDelam = this->giveNumberOfEnrichmentDomains();
    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * > (e->giveCrossSection());
    
    if ( dGroup == 1 ) {
        zBottom = - layeredCS->giveMidSurfaceZcoordFromBottom();
        zTop    =  0.;//this->giveDelaminationZCoord(dGroup);
    } else if (dGroup == nDelam + 1) {
        zBottom =  0.;//this->giveDelaminationZCoord(dGroup-1);
        zTop    = -layeredCS->giveMidSurfaceZcoordFromBottom() + layeredCS->computeIntegralThick();
    } else {
        zBottom =  0.;//this->giveDelaminationZCoord(dGroup-1);
        zTop    =  0.;//this->giveDelaminationZCoord(dGroup);
    }
#if DEBUG
    if ( zBottom > zTop ) {
        OOFEM_ERROR2("giveDelaminationGroupZLimits: Bottom z-coord is larger than top z-coord in dGroup. (%i)", dGroup);
    }
#endif
}


} // end namespace oofem
