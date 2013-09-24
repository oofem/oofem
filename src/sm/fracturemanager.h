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
class IntArray;
class Element;

class FractureManager;

/**
 * This class manages the fracture mechanics part
 *
 * @author Jim Brouzoulis
 */
#include "enumitem.h"

///@name Input fields for FractureManager
//@{
#define _IFT_FracManager_Name "fracmanager"
#define _IFT_fracManager_criteriaList "criterialist"  
//@}

#define FailureCriteria_DEF \
    ENUM_ITEM_WITH_VALUE(FC_Undefined, 0) \
    ENUM_ITEM_WITH_VALUE(Local, 1) \
    ENUM_ITEM_WITH_VALUE(Nonlocal, 2) 

enum FailureCriteriaType {
    FailureCriteria_DEF
};





class FailureCriteria
{
    // abstract class from all the different failure criterias should be derived
private:    
    
    FractureManager *fMan;          // pointer to its corresponding manager
    FailureCriteriaType type;       // local, nonlocal
    //FailureCriteriaName name;     // max strain, von Mises, effectivePlasticStrain, tsaiHill, J-integral, G, K etc.
    bool failedFlag;                // is the criteria fulfilled?

     
    

public:
    FailureCriteria(FailureCriteriaType type, FractureManager *fMan)
    { 
        this->type = type;
        //this->failedFlag = false;
        this->fMan = fMan;
    };
    FailureCriteria(){};
    ~FailureCriteria(){}; // must destroy object correctly
    Element *el;

    virtual int instanciateYourself(DataReader *dr){ return 1;};
    bool evaluateFCQuantities(Element *el, TimeStep *tStep);
    
    std::vector < std::vector < FloatArray > > quantities;
    FloatArray thresholds;
    
    //New
    std :: vector< bool > failedFlags;

    FailureCriteriaType giveType() { return this->type; }
    virtual const char *giveClassName() const = 0; //{ return FractureManagerClass; }

    virtual bool computeFailureCriteriaQuantities(TimeStep *tStep);
    virtual bool evaluateFailureCriteria() = 0;

    bool hasFailed( int i) { return failedFlags.at(i-1); }


};



class DamagedNeighborLayered : public FailureCriteria
{

public:
    DamagedNeighborLayered(FailureCriteriaType type, FractureManager *fMan) 
    : FailureCriteria(type,  fMan) {}

    virtual bool evaluateFailureCriteria();
    virtual const char *giveClassName() const { return "DamagedNeighborLayered"; }
    FloatArray layerDamageValues;
};


class FailureCriteriaManager
{
    // stores all the data for a given type of failure criteria

private:    
    FailureCriteriaType type;       // local, nonlocal 
    FractureManager *fMan;          // pointer to its corresponding manager

public:
    FailureCriteriaManager(FailureCriteriaType type, FractureManager *fMan)
        {
            this->type = type;
            this->fMan = fMan;
        };
    ~FailureCriteriaManager(){}; // must destroy object correctly

    std :: vector< FailureCriteria *> list;
    FailureCriteriaType giveType() { return this->type; }
    //FractureManager *giveFractureManager() { return this->fMan; }
    //FailureCriteriaType setType(FailureCriteriaType type) { return this->type = type; }
};


class FailureModuleElementInterface : public Interface
{
public:
    FailureModuleElementInterface() : Interface() {}
    virtual const char *giveClassName() const { return "FailureModuleElementInterface"; }
    virtual void computeFailureCriteriaQuantities(FailureCriteria *fc, TimeStep *tStep) {};  
};




class FractureManager
{
private:
    bool updateFlag;
    Domain *domain;

public:

    void setUpdateFlag(bool flag) { this->updateFlag = flag; };
    bool giveUpdateFlag() { return this->updateFlag; };
      
    void evaluateFailureCriterias(TimeStep *tStep); //Loop through all elements and evaluate criteria (if supported)
    
    
    void evaluateYourself(TimeStep *tStep);
    void updateXFEM(TimeStep *tStep);
    void updateXFEM(FailureCriteria *fc, TimeStep *tStep);

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


    std :: vector< FailureCriteriaManager* > criteriaManagers;
  //std :: vector< CrackManager*           > crackManagers;   // Keep track of all cracks - each crack may have several fronts/tips
  //std :: vector< PropagationLawManager*  > propagationLawManagers;

    
};





} // end namespace oofem
#endif // fracturemanager_h
