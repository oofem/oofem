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

#include "xfemmanager.h"
#include "flotmtrx.h"
#include "enrichmentitem.h"
#include "geometry.h"
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
    geometry = 0;
    enrichmentFunction = 0;

    //new JB
    this->enrichmentFunctionList = new AList< EnrichmentFunction >(0);
    this->geometryList = new AList< BasicGeometry >(0);
    numberOfEnrichmentFunctions = 0;    
    numberOfGeometryItems = 0;  
}
// remove - should ask for specific geom. can be multiple?
BasicGeometry *EnrichmentItem :: giveGeometry()
{
    //return xmanager->giveGeometry(this->geometry);
    //return this->giveGeometry( this->geometry );
    return this->geometryList->at(1);
}

BasicGeometry *EnrichmentItem :: giveGeometry(int n)
// Returns the n-th geometry.
{
    if ( this->geometryList->includes(n) ) {
        return this->geometryList->at(n);
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



bool EnrichmentItem :: isElementEnriched(Element *element) 
{
    //return this->giveGeometry()->intersects(element);
    return this->giveGeometry()->isInside(element);
}


bool EnrichmentItem :: isDofManEnriched(int nodeNumber)
{
    bool ret = false;
    // gets neighbouring elements of a node
    const IntArray *neighbours = domain->giveConnectivityTable()->giveDofManConnectivityArray(nodeNumber);
    for ( int i = 1; i <= neighbours->giveSize(); i++ ) {
        // for each of the neighbouring elements finds out whether it interacts with this EnrichmentItem
        if ( this->isElementEnriched( domain->giveElement( neighbours->at(i) ) ) ) {
            ret = true;
            break;
        }
    }

    return ret;
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
//    IR_GIVE_FIELD(ir, numberOfEnrichmentFunctions, IFT_XfemManager_numberOfEnrichmentFunctions, "numberofenrichmentfunctions"); // Macro// Macro
    IR_GIVE_FIELD(ir, numberOfGeometryItems, IFT_XfemManager_numberOfGeometryItems, "numberofgeometryitems");
    

    //IR_GIVE_FIELD(ir, geometry, IFT_EnrichmentItem_geometryItemNr, "geometryitem"); // Macro
    IR_GIVE_FIELD(ir, enrichmentDomainNumbers, IFT_EnrichmentItem_geometryItemNr, "geometryitem"); // Macro
    int numEnrichementDomains = this->enrichmentDomainNumbers.giveSize();
   // this->geometryList->growTo( numEnrichementDomains  );
       // for ( int i = 1; i <= numEnrichementDomains; i++ ) {
       //     this->geometryList->at(i) = 
       // }


    
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
    numberOfEnrichmentFunctions = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfEnrichmentFunctions, IFT_XfemManager_numberOfEnrichmentFunctions, "numberofenrichmentfunctions"); // Macro
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
    BasicGeometry *ge;
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

    
    geometryList->growTo(numberOfGeometryItems);
    for ( i = 1; i <= numberOfGeometryItems; i++ ) {
        mir = dr->giveInputRecord(DataReader :: IR_geoRec, i);
        result = mir->giveRecordKeywordField(name);
        if ( result != IRRT_OK ) {
            IR_IOERR(giveClassName(), __proc, IFT_RecordIDField, "", mir, result);
        }

        ge = CreateUsrDefGeometry( name.c_str() );

        if ( ge == NULL ) {
            OOFEM_ERROR2( "EnrichmentItem::instanciateYourself: unknown geometry (%s)", name.c_str() );
        }

        geometryList->put(i, ge);
        ge->initializeFrom(mir);
    }

#ifdef VERBOSE
    VERBOSE_PRINT0("Instanciated enrichment items ", nnode)
#endif
    return 1;
}


// DELAMINATION

IRResultType Delamination :: initializeFrom(InputRecord *ir)
{
    EnrichmentItem :: initializeFrom(ir);
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro
    int material = 0;
    numberOfEnrichmentFunctions = 1;

    IR_GIVE_FIELD(ir, this->delaminationZCoords, IFT_EnrichmentItem_Delamination_delaminationZCoords, "delaminationzcoords"); // Macro
    this->numberOfDelaminations = this->delaminationZCoords.giveSize();

    //write an instanciate method
    
    return IRRT_OK;
}

int
Delamination :: giveDelaminationGroupAt(double zeta) 
{
    //double zRef = Shell7Base :: giveLocalZetaCoord(gp);
    int nDelam = this->giveNumberOfDelaminations();
    for ( int j = 1; j <= nDelam; j++ ) {
        double zDelam = this->giveDelaminationZCoord(j);
        if ( zeta  < zDelam ) { //belong to the delamination group just below delamination #j. How to deal with poins that lie onthe boundary?
            return j;            
        }

    }
    return nDelam + 1;
            
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
    int nDelam = this->giveNumberOfDelaminations();
    //CrossSection *cs = e->giveCrossSection();
    LayeredCrossSection *layeredCS = dynamic_cast< LayeredCrossSection * > (e->giveCrossSection());
    
    if ( dGroup == 1 ) {
        zBottom = - layeredCS->giveMidSurfaceZcoordFromBottom();
        zTop    =  this->giveDelaminationZCoord(dGroup);
    } else if (dGroup == nDelam + 1) {
        zBottom =  this->giveDelaminationZCoord(dGroup-1);
        zTop    = -layeredCS->giveMidSurfaceZcoordFromBottom() + layeredCS->computeIntegralThick();
    } else {
        zBottom =  this->giveDelaminationZCoord(dGroup-1);
        zTop    =  this->giveDelaminationZCoord(dGroup);
    }
#if DEBUG
    if ( zBottom > zTop ) {
        OOFEM_ERROR2("giveDelaminationGroupZLimits: Bottom z-coord is larger than top z-coord in dGroup. (%i)", dGroup);
    }
#endif
}


} // end namespace oofem
