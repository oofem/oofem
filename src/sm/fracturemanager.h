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

#ifndef fracturemanager_h
#define fracturemanager_h

#include "alist.h"
#include "datareader.h"
#include "inputrecord.h"
#include "classtype.h"
#include "contextioresulttype.h"
#include "contextmode.h"
#include "timestep.h"
#include "interface.h"
/*
///@name Input fields for XfemManager
//@{
#define _IFT_XfemManager_numberOfGeometryItems "numberofgeometryitems"  // -> numberOfEnrichmentDomains
#define _IFT_XfemManager_numberOfEnrichmentItems "numberofenrichmentitems"
#define _IFT_XfemManager_numberOfEnrichmentFunctions "numberofenrichmentfunctions"
#define _IFT_XfemManager_name "xfemmanagername" ///< @todo Should this exist? / Mikael - No /JB
//@}
*/
namespace oofem {
class Domain;
//class EnrichmentItem;
class IntArray;
class Element;
//class DataStream;

/**
 * This class manages the fracture mechanics part
 *
 * @author Jim Brouzoulis
 */
#include "enumitem.h"
#define FailureCriteria_DEF \
    ENUM_ITEM_WITH_VALUE(FC_Undefined, 0) \
    ENUM_ITEM_WITH_VALUE(FC_MaxShearStress, 1) 

enum FailureCriteriaType {
    FailureCriteria_DEF
};

class Crack
{
protected:
    int numFronts; //will be tips in 2D, curves in 3D
    //PropagationLaw propLaw //Associated propagation law - should maybe be several
    //geometryDescription //Explicit -> list of nodes, EnrichmentDomain, …
    //updateGeometryDescription(){};
};


class FailureCriteria
{
private:    
    

    FailureCriteriaType type; 
    bool failedFlag;
public:
    FailureCriteria(FailureCriteriaType type){ this->type = type; };
    ~FailureCriteria(){}; // must destroy object correctly


    // list of all the quantities for each layer - quantities[layer][ip].arrayOfValues
    std::vector < std::vector < FloatArray > > quantities;
    FloatArray thresholds;


    FailureCriteriaType giveType() { return this->type; }
    bool evaluateFailureCriteria();
    bool evaluateFCQuantities() { return false; };

    bool hasFailed() { return failedFlag; }

/*
class LocalFailureCriteria : public FalureCriteria
{
maxStress, effectivePlasticStrain, tsaiHill, etc.
}

class NonLocalFailureCriteria : public FalureCriteria
{
J-integral, G, K
}
*/
};

class PropagationLaw
{
    //evaluatePropLaw(TimeStep *tStep){}; //Should update geometry, give direction and rate
};


class FailureModuleElementInterface : public Interface
{
public:
    FailureModuleElementInterface() : Interface() {}
    virtual const char *giveClassName() const { return "FailureModuleElementInterface"; }
    virtual void evaluateFailureCriteriaQuantities(FailureCriteria *fc, TimeStep *tStep) {};
};


class FractureManager
{
public:
    Domain *domain;

    
    AList < Crack           > *crackList;        // Keep track of all cracks - each crack may have several fronts/tips
    AList < FailureCriteria > *failureCriterias; // All failure criterias to evaluate
    AList < PropagationLaw  > *propagationLaws;
    
    bool needsUpdate;
    //createCrack(){}; //Should be able to create a new crack based on failure criterias?
    //removeCrack(int number){};
    
    void evaluateFailureCriterias(TimeStep *tStep); //Loop through all elements and evaluate criteria (if supported)
    virtual void evaluateFailureCriteria(FailureCriteria *fc, Element *el, TimeStep *tStep);
    
    void update(TimeStep *tStep);
    
    //evaluatePropagationLaws(){};
    /// Constructor.
    FractureManager(Domain *domain);
    /// Destructor.
    ~FractureManager();
    
    IRResultType initializeFrom(InputRecord *ir);
    int instanciateYourself(DataReader *dr);
    virtual classType giveClassID() const { return FractureManagerClass; }
    const char *giveClassName() const { return "FractureManager"; }
    const char *giveInputRecordName() const { return "FractureManager"; }
    void clear();
    Domain *giveDomain() { return this->domain; }
    
};





} // end namespace oofem
#endif // fracturemanager_h
