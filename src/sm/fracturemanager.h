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
class FailureCriteria;
class FailureCriteriaStatus;
/**
 * This class manages the fracture mechanics part
 *
 * @author Jim Brouzoulis
 */
#include "enumitem.h"

///@name Input fields for FractureManager
//@{
#define _IFT_FracManager_Name "fracmanager"
#define _IFT_FracManager_numcriterias "numcriterias"
#define _IFT_FracManager_verbose "verbose"

// Failure criterias
#define _IFT_DamagedNeighborLayered_Name "damagedneighborlayered"
#define _IFT_DamagedNeighborLayered_DamageThreshold "damagethreshold"
//@}

// IPLocal - each integration point, ELLocal - one per element
#define FailureCriteria_DEF \
    ENUM_ITEM_WITH_VALUE(FC_Undefined, 0) \
    ENUM_ITEM_WITH_VALUE(IPLocal, 1) \
    ENUM_ITEM_WITH_VALUE(ELLocal, 2) \
    ENUM_ITEM_WITH_VALUE(Nonlocal, 3) 

enum FailureCriteriaType {
    FailureCriteria_DEF
};





class FailureCriteriaStatus
{
    // abstract class from all the different failure criterias should be derived
private:    
    
    //FractureManager *fMan;         
    FailureCriteria *failCrit;      // pointer to the corresponding failure criteria
    FailureCriteriaType type;       // local, nonlocal
    //FailureCriteriaName name;     // max strain, von Mises, effectivePlasticStrain, tsaiHill, J-integral, G, K etc.
    bool failedFlag;                // is the criteria fulfilled?
    int number;
     
    

public:
    FailureCriteriaStatus(int number, FailureCriteria *failCrit)
    { 
        this->number = number;
        this->failCrit = failCrit;
    };

    FailureCriteriaStatus(Element *el, FailureCriteria *failCrit)
    { 
        this->el = el;
        this->failCrit = failCrit;
    };

    FailureCriteriaStatus(){};
    ~FailureCriteriaStatus(){}; // must destroy object correctly
    Element *el;

    std::vector < std::vector < FloatArray > > quantities;
    FloatArray thresholds;
    std :: vector< bool > failedFlags;


    //FailureCriteriaType giveType() { return this->giveFailureCriteria().giveType(); };
    FailureCriteria *giveFailureCriteria() { return this->failCrit; };


    bool hasFailed( int i) { return failedFlags.at(i-1); }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int instanciateYourself(DataReader *dr){ return 1;};
    virtual const char *giveClassName() const { return "FailureCriteriaStatus"; }

};





class FailureCriteria
{
private:    
    FailureCriteriaType type;       // local, nonlocal 
    FractureManager *fMan;          // pointer to its corresponding manager
    int number;
public:

    FailureCriteria(int number, FractureManager *fMan)
        {
            this->number = number;
            this->fMan = fMan;
        };
    ~FailureCriteria(){}; // must destroy object correctly

    std :: vector< FailureCriteriaStatus *> list;

    FailureCriteriaType giveType() { return this->type; }
    FractureManager *giveFractureManager() { return this->fMan; }
    void setType(FailureCriteriaType type) { this->type = type; };

    virtual IRResultType initializeFrom(InputRecord *ir);
    int instanciateYourself(DataReader *dr);
    virtual const char *giveClassName() const { return "FailureCriteria"; }


    virtual FailureCriteriaStatus *CreateStatus(Element *el, FailureCriteria *failCrit) const = 0;
    virtual bool computeFailureCriteriaQuantities(FailureCriteriaStatus *fcStatus, TimeStep *tStep);
    virtual bool evaluateFCQuantities(Element *el, TimeStep *tStep) { return false; }; // defaults to no implementation
    virtual bool evaluateFailureCriteria(FailureCriteriaStatus *fcStatus) = 0;
};




class DamagedNeighborLayeredStatus : public FailureCriteriaStatus
{

public:
    DamagedNeighborLayeredStatus(Element *el, FailureCriteria *failCrit) 
    : FailureCriteriaStatus(el, failCrit) {}

    FloatArray layerDamageValues; 
};


class DamagedNeighborLayered : public FailureCriteria
{
private:
    double DamageThreshold;

public:
    DamagedNeighborLayered(int number, FractureManager *fracMan) 
    : FailureCriteria(number,  fracMan) {}

    virtual bool evaluateFailureCriteria(FailureCriteriaStatus *fcStatus);
    virtual const char *giveClassName() const { return "DamagedNeighborLayered"; }
    virtual const char *giveInputRecordName() const { return _IFT_DamagedNeighborLayered_Name; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual FailureCriteriaStatus *CreateStatus(Element *el, FailureCriteria *failCrit) const
        { return new DamagedNeighborLayeredStatus(el, failCrit); };


};



class FailureModuleElementInterface : public Interface
{
public:
    FailureModuleElementInterface() : Interface() {}
    virtual const char *giveClassName() const { return "FailureModuleElementInterface"; }
    virtual void computeFailureCriteriaQuantities(FailureCriteriaStatus *fc, TimeStep *tStep) {};  
};




class FractureManager
{
private:
    bool updateFlag;
    Domain *domain;

public:
    /// Constructor.
    FractureManager(Domain *domain);
    /// Destructor.
    ~FractureManager();

    void setUpdateFlag(bool flag) { this->updateFlag = flag; };
    bool giveUpdateFlag() { return this->updateFlag; };
      
    void evaluateFailureCriterias(TimeStep *tStep); //Loop through all elements and evaluate criteria (if supported)
    
    
    void evaluateYourself(TimeStep *tStep);
    void updateXFEM(TimeStep *tStep);
    void updateXFEM(FailureCriteriaStatus *fc, TimeStep *tStep);


    
    IRResultType initializeFrom(InputRecord *ir);
    int instanciateYourself(DataReader *dr);
    virtual classType giveClassID() const { return FractureManagerClass; }
    const char *giveClassName() const { return "FractureManager"; }
    const char *giveInputRecordName() const { return "FractureManager"; }
    void clear();
    Domain *giveDomain() { return this->domain; }

    

    std :: vector< FailureCriteria* > criteriaList;
  //std :: vector< CrackManager*           > crackManagers;   // Keep track of all cracks - each crack may have several fronts/tips
  //std :: vector< PropagationLawManager*  > propagationLawManagers;

    
};





} // end namespace oofem
#endif // fracturemanager_h
