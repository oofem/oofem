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
#define FailureCriteria_DEF \
    ENUM_ITEM_WITH_VALUE(FC_Undefined, 0) \
    ENUM_ITEM_WITH_VALUE(FC_MaxShearStress, 1) \
    ENUM_ITEM_WITH_VALUE(FC_DamagedNeighborCZ, 2) \
    ENUM_ITEM_WITH_VALUE(Local, 3) \
    ENUM_ITEM_WITH_VALUE(Nonlocal, 4) 

enum FailureCriteriaType {
    FailureCriteria_DEF
};

//enum QuantityType {
//    FailureCriteria_DEF
//};




class FailureCriteriaQuantity
{
private:    
    
    //QuantityType type;       // local, nonlocal

    std::vector < std::vector < FloatArray > > quantities;

public:
    //FailureCriteriaQuantity( QuantityType type)
    //{ 
      //  this->type = type;
    //};
    ~FailureCriteriaQuantity(){}; 
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
        this->failedFlag = false;
        this->fMan = fMan;
    };
    FailureCriteria(){};
    ~FailureCriteria(){}; // must destroy object correctly
    Element *el;
    int instanciateYourself(DataReader *dr);

    // list of all the quantities for each layer - quantities[layer][ip].arrayOfValues
    std::vector < std::vector < FloatArray > > quantities;
    FloatArray thresholds;
    
    //New
    std::vector < FailureCriteriaQuantity > elQuantities;
    //bool evaluateFCQuantities(FailureCriteriaQuantity *quantity, TimeStep *tStep){}; 
    std :: vector< bool > failedFlags;

    FailureCriteriaType giveType() { return this->type; }
    virtual bool evaluateFailureCriteria(TimeStep *tStep);
    virtual bool computeFailureCriteriaQuantities(TimeStep *tStep)
        { return evaluateFailureCriteria(tStep); };
    virtual bool evaluateFailureCriteria();
    bool evaluateFCQuantities(Element *el, TimeStep *tStep); 

    bool hasFailed() { return failedFlag; }
    bool hasFailed( int i) { return failedFlags.at(i-1); }


};



class DamagedNeighborLayered : public FailureCriteria
{

public:
    DamagedNeighborLayered(FailureCriteriaType type, FractureManager *fMan);
    virtual bool evaluateFailureCriteria();
    //std::vector < std::vector < FloatArray > > layerDamageValues;
    FloatArray layerDamageValues;
    //{ 
        //this->failedFlag = false;
        //elQuantities.resize( fMan->giveDomain()->giveNumberOfElements() );
    //};

};

class FailureCriteriaManager
{
    // stores all the data for a given type o failure criteria
private:    
    
    //FractureManager *fMan;          // pointer to its corresponding manager

public:
    FailureCriteriaManager()
    { 
    };
    ~FailureCriteriaManager(){}; // must destroy object correctly

    std :: vector< FailureCriteria *> list;
};


class FailureModuleElementInterface : public Interface
{
public:
    FailureModuleElementInterface() : Interface() {}
    virtual const char *giveClassName() const { return "FailureModuleElementInterface"; }
    virtual void computeFailureCriteriaQuantities(FailureCriteria *fc, TimeStep *tStep) {};
    virtual void computeFailureCriteriaQuantities(FailureCriteria *fc, FailureCriteriaQuantity quantities, FailureCriteriaType type,  TimeStep *tStep){}; 
};


class FractureManager
{
private:
    bool updateFlag;
    Domain *domain;

public:

    //AList < Crack           > *crackList;        // Keep track of all cracks - each crack may have several fronts/tips
    AList < FailureCriteria > *failureCriterias; // All failure criterias to evaluate
    //AList < PropagationLaw  > *propagationLaws;
    void setUpdateFlag(bool flag) { this->updateFlag = flag; };
    bool giveUpdateFlag() { return this->updateFlag; };
      
    void evaluateFailureCriterias(TimeStep *tStep); //Loop through all elements and evaluate criteria (if supported)
    
    
    void update(TimeStep *tStep);
    void updateXFEM(TimeStep *tStep);
    void updateXFEM(FailureCriteria *fc, TimeStep *tStep);

    void updateXFEM(FractureManager *fMan, TimeStep *tStep);//new

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


    
};





} // end namespace oofem
#endif // fracturemanager_h
