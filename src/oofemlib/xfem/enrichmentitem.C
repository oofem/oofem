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
#include "floatmatrix.h"
#include "enrichmentitem.h"
#include "element.h"
#include "enrichmentfunction.h"
#include "enrichmentdomain.h"
#include "cltypes.h"
#include "connectivitytable.h"
#include "oofem_limits.h"
#include "classfactory.h"
#include "fracturemanager.h"
namespace oofem {

REGISTER_EnrichmentItem( CrackTip )
REGISTER_EnrichmentItem( CrackInterior )
REGISTER_EnrichmentItem( Inclusion )
REGISTER_EnrichmentItem( Delamination )

EnrichmentItem :: EnrichmentItem(int n, XfemManager *xMan, Domain *aDomain) : FEMComponent(n, aDomain)
{
    this->xMan = xMan;
    this->enrichmentFunctionList = new AList< EnrichmentFunction >(0);
    this->enrichmentDomainList = new AList< EnrichmentDomain >(0); 
    this->numberOfEnrichmentFunctions = 1;
    this->numberOfEnrichmentDomains = 1;
    this->startOfDofIdPool = -1;
    this->enrichesDofsWithIdArray = new IntArray;
}

EnrichmentItem :: ~EnrichmentItem()
{
    delete this->enrichesDofsWithIdArray;
}


EnrichmentFunction *EnrichmentItem :: giveEnrichmentFunction(int n)
{
    // Returns the n-th geometry.
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
    // Checks if any of the dofmanagers are enriched
    for ( int i = 1; i <= element->giveNumberOfDofManagers(); i++ ) {
        DofManager *dMan = element->giveDofManager(i);
        if ( isDofManEnrichedByEnrichmentDomain(dMan, edNumber) ){
            return true;
        }
    }
    return false;
}

bool EnrichmentItem :: isElementFullyEnrichedByEnrichmentDomain(const Element *element, int edNumber) 
{
    // Checks if all of the dofmanagers are enriched
    int count = 0;
    for ( int i = 1; i <= element->giveNumberOfDofManagers(); i++ ) {
        DofManager *dMan = element->giveDofManager(i);
        if ( isDofManEnrichedByEnrichmentDomain(dMan, edNumber) ){
            count++;
        }
    }
    return count == element->giveNumberOfDofManagers();
}

bool EnrichmentItem :: isDofManEnriched(DofManager *dMan)
{
    for ( int i = 1; i <= this->enrichmentDomainList->giveSize() ; i++ ) {
        if ( isDofManEnrichedByEnrichmentDomain(dMan, i) ){
            return true;
        }
    }
    return false;
}

bool EnrichmentItem :: isDofManEnrichedByEnrichmentDomain(DofManager *dMan, int edNumber)
{
    EnrichmentDomain *ed = this->enrichmentDomainList->at(edNumber);
    return ed->isDofManagerEnriched(dMan);
}


IRResultType EnrichmentItem :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, this->enrichmentDomainNumbers, _IFT_EnrichmentItem_domains);
    this->numberOfEnrichmentDomains = this->enrichmentDomainNumbers.giveSize();

    IR_GIVE_OPTIONAL_FIELD(ir, enrichmentFunction, _IFT_EnrichmentItem_function);
    
    return IRRT_OK;
}


IRResultType Inclusion :: initializeFrom(InputRecord *ir)
{
    EnrichmentItem :: initializeFrom(ir);
    const char *__proc = "initializeFrom"; 
    IRResultType result; 
    int material = 0;
    IR_GIVE_FIELD(ir, material, _IFT_Inclusion_material); 
    this->mat = this->giveDomain()->giveMaterial(material);
    this->numberOfEnrichmentFunctions = 1;
    // Not sure this should be an input paramter but should instead be determined by the ei which describes the physical model /JB
    //IR_GIVE_OPTIONAL_FIELD(ir, numberOfEnrichmentFunctions, _IFT_XfemManager_numberOfEnrichmentFunctions, "numberofenrichmentfunctions");
    return IRRT_OK;
}



int EnrichmentItem :: instanciateYourself(DataReader *dr)
{
    const char *__proc = "instanciateYourself"; // Required by IR_GIVE_FIELD macro
    IRResultType result; // Required by IR_GIVE_FIELD macro
    std :: string name;

    // Instanciate enrichment functions
    this->enrichmentFunctionList->growTo(numberOfEnrichmentFunctions);
    for ( int i = 1; i <= this->numberOfEnrichmentFunctions; i++ ) {
        InputRecord *mir = dr->giveInputRecord(DataReader :: IR_enrichFuncRec, i);
        result = mir->giveRecordKeywordField(name);

        if ( result != IRRT_OK ) {
            IR_IOERR(giveClassName(), __proc, "", mir, result);
        }
        
        EnrichmentFunction *ef = classFactory.createEnrichmentFunction( name.c_str(), i, this->xMan->giveDomain() );
        if ( ef == NULL ) {
            OOFEM_ERROR2( "EnrichmentItem::instanciateYourself: unknown enrichment function (%s)", name.c_str() );
        }
        enrichmentFunctionList->put(i, ef);
        ef->initializeFrom(mir);
    }

    // Instanciate enrichment domains
    enrichmentDomainList->growTo(numberOfEnrichmentDomains);

    for ( int i = 1; i <= numberOfEnrichmentDomains; i++ ) {
        InputRecord *mir = dr->giveInputRecord(DataReader :: IR_geoRec, i);
        result = mir->giveRecordKeywordField(name);
        if ( result != IRRT_OK ) {
            IR_IOERR(giveClassName(), __proc, "", mir, result);
        }

        EnrichmentDomain *ed = classFactory.createEnrichmentDomain( name.c_str() ); 
        if ( ed == NULL ) {
            OOFEM_ERROR2( "EnrichmentItem::instanciateYourself: unknown enrichment domain (%s)", name.c_str() );
        }
        ed->setNumber(i);
        this->enrichmentDomainList->put(i, ed);
        ed->initializeFrom(mir);

    }

    // Set start of the enrichment dof pool for the given EI
//    int xDofPoolAllocSize = this->giveEnrichesDofsWithIdArray()->giveSize() * this->giveNumberOfEnrDofs() * this->giveNumberOfEnrichmentDomains(); 
    int xDofPoolAllocSize = this->giveEnrichesDofsWithIdArray()->giveSize() * this->giveNumberOfEnrDofs() * 100; //@todo should allocate max number of dofs but this requires knowing the max number of layers in my case /JB
    this->startOfDofIdPool = this->giveDomain()->giveNextFreeDofID(xDofPoolAllocSize);

    return 1;
}


void 
EnrichmentItem :: addEnrichmentDomain( EnrichmentDomain *ed )
{
    // Appends the enrichment domain ed to the list
    // Does not check if there is a duplicate in the list
    this->numberOfEnrichmentDomains++;
    enrichmentDomainList->growTo(numberOfEnrichmentDomains);

    ed->setNumber(this->numberOfEnrichmentDomains);
    this->enrichmentDomainList->put(this->numberOfEnrichmentDomains, ed);
}

void
EnrichmentItem :: giveEIDofIdArray(IntArray &answer, int enrichmentDomainNumber)
{
    // Returns an array containing the dof Id's of the new enrichment dofs pertinent to the ei. 
    // Note: the dof managers may not support all/any of these new potential dof id's. Eg. a 
    // beam will not be able to account for a pressure dof. 
    IntArray *enrichesDofsWithIdArray = this->giveEnrichesDofsWithIdArray();
    int eiEnrSize = enrichesDofsWithIdArray->giveSize();

    answer.resize(eiEnrSize);
    int xDofAllocSize = eiEnrSize * this->giveNumberOfEnrDofs(); // number of new dof id's the ei will allocate
    for ( int i = 1; i <= eiEnrSize; i++ ) {
        answer.at(i) = this->giveStartOfDofIdPool() + (enrichmentDomainNumber-1)*xDofAllocSize + (i-1); 
    }
}

void
EnrichmentItem :: computeDofManDofIdArray(IntArray &answer, DofManager *dMan, int enrichmentDomainNumber)
{
    // Gives an array containing the dofId's that should be created as new dofs (which dofs to enrich).
    IntArray *enrichesDofsWithIdArray = this->giveEnrichesDofsWithIdArray();
    int eiEnrSize = enrichesDofsWithIdArray->giveSize();
    
    // Go through the list of dofs that the EI supports and compare with the available dofs in the dofMan.
    // If the dofMan has support for the particular dof add it to the list.
    // Store matches in dofMask
    IntArray dofMask(eiEnrSize); dofMask.zero();
    int count = 0; 
    for ( int i = 1; i <= eiEnrSize; i++ ) {
        if ( dMan->hasDofID( (DofIDItem) enrichesDofsWithIdArray->at(i) ) ) {
            count++;
            dofMask.at(count) = dMan->giveDofWithID( enrichesDofsWithIdArray->at(i) )->giveNumber();
        }
    }

    answer.resize(count);
    int xDofAllocSize = eiEnrSize * this->giveNumberOfEnrDofs(); // number of new dof id's the ei will allocate
    for ( int i = 1; i <= count; i++ ) {
        answer.at(i) = this->giveStartOfDofIdPool() + (enrichmentDomainNumber-1)*xDofAllocSize + dofMask.at(i)-1 ;
    }

}

int 
EnrichmentItem :: giveNumberOfEnrDofs() 
{ 
    // returns the number of new xfem dofs a particular EI wants to add per regular dof
    int temp=0;
    for ( int i = 1; i <= this->giveNumberOfEnrichmentfunctions(); i++ ) { 
        EnrichmentFunction *ef = this->giveEnrichmentFunction(i);
        temp += ef->giveNumberOfDofs(); // = number of functions associated with a particular enrichment function, e.g. 4 for the common branch functions.
    }
    return temp; 
} 

void 
EnrichmentItem :: updateGeometry(TimeStep *tStep, FractureManager *fMan)
{
    //this->needsUpdate = false;
    Domain *domain= this->giveDomain();   

    for ( int i = 1; i <= domain->giveNumberOfElements(); i++ ) { 
        printf( "\n -------------------------------\n");
        Element *el = domain->giveElement(i);

        for ( int j = 1; j <= fMan->failureCriterias->giveSize(); j++ ) {
            FailureCriteria *fc = fMan->failureCriterias->at(j);
            fMan->evaluateFailureCriteria(fc, el, tStep);

            if ( Delamination *dei = dynamic_cast< Delamination * > (this) )  {
                dei->updateGeometry(tStep, fMan, el, fc); //not an overloaded function, change the name
            }
        }
    }

}



Inclusion :: Inclusion(int n, XfemManager *xm, Domain *aDomain) : EnrichmentItem(n, xm, aDomain)
{ 
    this->enrichesDofsWithIdArray->setValues(3, D_u, D_v, D_w);
}


//------------------
// DELAMINATION
//------------------

void 
Delamination :: updateGeometry(TimeStep *tStep, FractureManager *fMan, Element *el, FailureCriteria *fc)
{
    for ( int i = 1; i <= fc->quantities.size(); i++ ) {
        if ( fc->hasFailed(i) ) { // interface has failed
            printf( "Element %d fails in interface %d \n", el->giveNumber(), i );
            fMan->setUpdateFlag(true);

            // should create a new *ed at the level given by interface i iff it does not already exist
            // add coord
            FloatArray xiCoords;
            dynamic_cast <LayeredCrossSection * > ( el->giveCrossSection() )->giveInterfaceXiCoords(xiCoords);
            double xi = xiCoords.at(i); // current xi-coord
            // find if xi is in this->enrichmentDomainXiCoords
            bool flag=false;
            int num = 0;
            for ( int j = 1; j <= this->enrichmentDomainXiCoords.giveSize(); j++ ) {
                if ( abs(xi-this->enrichmentDomainXiCoords.at(j)) < 1.0e-6 ) {
                    flag = true;
                    num = j;
                }
            }


            IntArray dofManNumbers, elDofMans;
            elDofMans = el->giveDofManArray();
            for ( int i = 1; i <= el->giveNumberOfDofManagers(); i++ ) {  
                // ugly piece of code that will skip enrichment of dofmans that have any bc's
                // which is not generally what you want
                #if 1
                bool hasBc= false;
                for ( int j = 1; j <= el->giveDofManager(i)->giveNumberOfDofs(); j++ ) {
                    if ( el->giveDofManager(i)->giveDof(j)->hasBc(tStep) ) {
                        hasBc = true;
                        continue;
                    }
                }
                #endif
                if ( !hasBc) {
                    dofManNumbers.followedBy(elDofMans.at(i));
                }
            }
            

            //dofManNumbers.printYourself();

            if ( flag ) { //in list only add dofmans
                dynamic_cast< DofManList * > ( this->giveEnrichmentDomain(num) )->addDofManagers( dofManNumbers );
                //dynamic_cast< DofManList * > ( this->giveEnrichmentDomain(num) )->updateEnrichmentDomain(dofManNumbers);
            } else { //create ed
                int numED = this->giveNumberOfEnrichmentDomains();
                EnrichmentDomain *ed = classFactory.createEnrichmentDomain( "DofManList" ); 
                DofManList *dml = dynamic_cast< DofManList * > ( ed );
                dml->addDofManagers( dofManNumbers ); // add the dofmans of the el to the list
                this->addEnrichmentDomain(ed);
                this->enrichmentDomainXiCoords.resizeWithValues(numED+1);
                this->enrichmentDomainXiCoords.at(numED+1) = xiCoords.at(i);

            }                        

        }
    }


}

Delamination :: Delamination(int n, XfemManager *xm, Domain *aDomain) : EnrichmentItem(n, xm, aDomain)
{ 
    this->enrichesDofsWithIdArray->setValues(6, D_u, D_v, D_w, W_u, W_v, W_w);
}


IRResultType Delamination :: initializeFrom(InputRecord *ir)
{
    this->numberOfEnrichmentFunctions = 1; // must be set before EnrichmentItem :: initializeFrom(ir) is called
    EnrichmentItem :: initializeFrom(ir);
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, this->enrichmentDomainXiCoords, _IFT_Delamination_xiCoords);
    if ( this->numberOfEnrichmentDomains != this->enrichmentDomainXiCoords.giveSize() ) { // move to checkConsistency
        OOFEM_ERROR3( "EnrichmentItem :: initializeFrom: size of enrichmentDomainXiCoords (%i) differs from numberOfEnrichmentDomains (%i)", 
                       this->enrichmentDomainXiCoords.giveSize(), this->numberOfEnrichmentDomains );
    }
    int material = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, material, _IFT_Delamination_CohesiveZoneMaterial);
    if ( material > 0 ) {
        this->mat = this->giveDomain()->giveMaterial(material);
    }
    //write an instanciate method

    //enrichmentDomainInterfaceList

    return IRRT_OK;
}


void 
Delamination :: giveActiveDelaminationXiCoords(FloatArray &xiCoords, Element *element) 
{
    // Goes through the list of delaminations and checks which are active for a given element
    int nDelam = this->giveNumberOfEnrichmentDomains(); // max possible number
    int pos = 1;
    xiCoords.resize(0);
    for ( int i = 1; i <= nDelam; i++ ) {
        if( this->isElementFullyEnrichedByEnrichmentDomain(element, i) ) {
            xiCoords.resizeWithValues(pos);
            xiCoords.at(pos) = this->giveDelaminationXiCoord(i);
            pos++;
        } 
    }
};


// remove???
double 
Delamination :: giveDelaminationZCoord(int n, Element *element) 
{
    //AList<double> *xiCoordList;
    int nDelam = this->giveNumberOfEnrichmentDomains(); // max possible number
    int pos = 1;
    for ( int i = 1; i <= nDelam; i++ ) {
        if( this->isElementEnriched(element) ) {
            //xiCoordList.put( pos, this->delaminationZCoords.at(i) );
            pos++;
        } 
    }
    return 0.;
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

// unneccesary??
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
#ifdef DEBUG
    if ( zBottom > zTop ) {
        OOFEM_ERROR2("giveDelaminationGroupZLimits: Bottom z-coord is larger than top z-coord in dGroup. (%i)", dGroup);
    }
#endif
}





} // end namespace oofem
