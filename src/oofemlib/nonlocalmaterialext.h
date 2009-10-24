/* $Header: /home/cvs/bp/oofem/oofemlib/src/nonlocalmaterialext.h,v 1.14.4.1 2004/04/05 15:19:43 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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


//
// class NonlocalMaterialExtension
//

#ifndef nonlocalmaterialext_h
#define nonlocalmaterialext_h

//#include "material.h"

#include "matstatus.h"
#include "dynalist.h"
#include "interface.h"
//
// local integration record - stores pointer to gp and its integration weight
// Also remote integration record should be defined in case of parallel scheme
// using domain decomposition.
//
/**
 * Structure containing reference to integration point and its corresponding nonlocal integration weight.
 * Used by nonlocal constitutive models based on integral averaging procedure, where in each integration
 * point the corresponding list of influencing integration points is kept, together with their weights.
 * This structure encapsulates the reference to influencing integration point and its corresponding weight.
 */
struct localIntegrationRecord {
    /// Refence to influencing integration point.
    GaussPoint *nearGp;
    /// Corresponding integration weight.
    double weight;
};

//template <class localIntegrationRecord> class dynaList;

/**
 * Abstract base class for all nonlocal constitutive model statuses. Introduces the list of
 * localIntegrationRecords stored in each integration point, where references to all influencing
 * integration points with their integration weights are kept. Also the total volume associated to
 * corresponding integration point is kept. Services for accessing of integration list as well as
 * services for manipulating and requesting integration volume are provided.
 * Generally speaking, the nonlocal weight function with "bounded" or limited support is assumed.
 * When nonlocal weight function unbounded support is used, then keeping the list of
 * influencing integration points has no sence (because all integration points in domain participate)
 * and should be appropriate not to use afore mentioned list and redefine the
 * MaterialStatus::buildNonlocalPointTable sevice to void service.
 *
 * @see localIntegrationRecord structure.
 * @see MaterialStatus::buildNonlocalPointTable service.
 */
class NonlocalMaterialStatusExtensionInterface : public Interface
{
protected:
    /** List containing localIntegrationRecord values.
     */
    dynaList< localIntegrationRecord >integrationDomainList;
    /// Nonlocal volume of corresponding integration point.
    double integrationScale;

public:
    /**
     * Constructor.
     */
    NonlocalMaterialStatusExtensionInterface();

    /// Destructor.
    ~NonlocalMaterialStatusExtensionInterface();

    /**
     * Returns integration list of receiver. Contains localIntegrationRecord structures, containing
     * references to integration points and their weights that influence to nonlocal average in
     * receiver's associated integration point.
     */
    dynaList< localIntegrationRecord > *giveIntegrationDomainList() { return & integrationDomainList; }
    /// Returns associated integration scale.
    double giveIntegrationScale() { return integrationScale; }
    /// Sets associated integration scale.
    void   setIntegrationScale(double val) { integrationScale = val; }
    /// clears the integration list of receiver
    void clear() {integrationDomainList.clear();}
};




/**
 * Abstract base class for all nonlocal materials. Nonlocal in sence, that response in particular
 * point depends not only on state in that point, but also takes into account state of surrounding
 * points. Response typically depends on some nonlocal quantity obtained as nonlocal average over
 * some characteristic volume.
 * General services for updating domain before nonlocal average, building table of influencing
 * integration points for given integration point and computing nonlocal weight function are
 * declared. The general service for building table of influencing
 * integration points for given integration point is also implemented.
 * The use of multiple inheritance is assumed. Typically, the class representing nonlocal
 * constitutive model is derived both from class representing local model and from this class or
 * from one of its derived classes (which declare services and variables corresponding to specific
 * analysis type).
 */
class NonlocalMaterialExtensionInterface : public Interface
{
protected:
    /*
     * It is necessary, mainly due to resulting efficiency, to compute variable(s)
     * which are nonlocally averaged in advance, before average process begins.
     * The loop over all  integration points is typically made to compute these variables.
     * To prevent doing this multiple times at the same solution state,
     * the modification time mark is kept.
     * In the present implementation, there is one comon stateCounter for the whole domain.
     * This implies, that all variables (which undergo averaging) for all material models of the domain
     * should be prepared at the same time in updateDomainBeforeNonlocAverage method.
     * It is believed that this is general enough and can somehow handle even the case of multiple models
     * with different parameters beeing averaged (but is this realistic?).
     * If this scheme will not be enough general, then the state counter
     * can be kept as attribute of NonlocalMaterialExtensionInterface, so independently for
     * each material model. Each model will be then updated in separate call. But in the case of
     * several material models of the same type (with diferent params) this will lead to
     * multiple update, which can not be avoided, althoug it is renundant.
     *
     * StateCounterType lastUpdatedStateCounter;
     */
    Domain *domain;
    /// map indicating regions to skip (region - cross section model)
    IntArray regionMap;
    /// flag indicating whether to keep nonlocal interaction tables of integration points cached
    bool permanentNonlocTableFlag;
public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param n material number
     * @param d domain to which new material will belong
     */
    NonlocalMaterialExtensionInterface(Domain *d);
    /// Destructor.
    ~NonlocalMaterialExtensionInterface()                { }


    /**
     * Updates data in all integration points before nonlocal average takes place.
     * It is necessary, mainly due to resulting efficiency, to compute variable(s)
     * which are nonlocally averaged in advance, before average process begins.
     * These variables must be stored in integration point's associated statuses.
     * This function updates the whole problem domain, by updating all integration points
     * values, which take part in nonlocal average process. All elements are updated using
     * Element::updateBeforeNonlocalAverage abstract service, which in turn updates all
     * integration points associated with particular element.
     * The service used to update element integration point depends on analysis type and
     * is specified by element-specific type (like StructuralElement) corresponding to anlysis type.
     * This service can be invoked multiple times, but update for specific material is done only once, because
     * last modification time mark is kept.
     * @see Element::updateBeforeNonlocalAverage.
     */
    void updateDomainBeforeNonlocAverage(TimeStep *atTime);

    /**
     * Builds list of integration points which take part in nonlocal average in given integration point.
     * This list is stored in integration point corresponding nonlocal status.
     * Generally speaking, the nonlocal weight function with "bounded" or limited support is assumed.
     * When nonlocal weight function unbounded support is used, then keeping the list of
     * influencing integration points would be wasting of space and should be cleared after averaging has
     * been finished in integration point. The endIPNonlocalAverage method will ensure this.
     */
    void buildNonlocalPointTable(GaussPoint *gp);

    /**
     * Rebuild list of integration points which take part
     * in nonlocal average in given integration point.
     * if contributingElems param is not NULL, then it is assumed that it contains
     * a COMPLETE list of elements contributing to receiver. Ifis equal to NULL
     * existing list is cleared and  buildNonlocalPointTable service is invoked.
     */
    void rebuildNonlocalPointTable(GaussPoint *gp, IntArray *contributingElems);
    /**
     * Returns integration list corresponding to given integration point.
     * Contains localIntegrationRecord structures, containing
     * references to integration points and their weights that influence to nonlocal average in
     * receiver's associated integration point.
     * Rebuilds the IP list by calling  buildNonlocalPointTable(GaussPoint *gp) if not available
     */
    dynaList< localIntegrationRecord > *giveIPIntegrationList(GaussPoint *gp);
    /**
     * Computes the value of nonlocal weight function in given point.
     * @param src coordinates of source point.
     * @param coord coordinates of point, where nonlocal weight function is evaluated.
     * @return value of weight function.
     */
    virtual double computeWeightFunction(const FloatArray &src, const FloatArray &coord) = 0;

    /**
     * Determines the number of material regions of domain.
     * In the current implementation the region is associated with cross section model.
     */
    int giveNumberOfRegions();
    /**
     * Returns the region id of given element
     * @param element pointer to element which region id is requsted
     * @return region id (number) for this element
     */
    //int giveElementRegion (Element* element);
    /**
     * Determines, whether receiver has bounded weighting function (limited support)
     * @return true if weighting function bounded, zero otherwise
     */
    virtual int hasBoundedSupport() {return 1;}
    /**
     * Determines the width (radius) of limited support of weighting function
     */
    virtual void giveSupportRadius(double &radius) { radius = 0.0; }

    /// returns reference to domain
    Domain *giveDomain() { return this->domain; }

    IRResultType initializeFrom(InputRecord *ir);
    /** Setups the input record string of receiver
     * @param str string to be filled by input record
     * @param keyword print record keyword (default true)
     */
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);
    /*
     * Creates new copy of associated status and inserts it into given integration point.
     * @param gp Integration point where newly created status will be stored.
     * @return reference to new status.
     */
    //MaterialStatus* CreateStatus (GaussPoint* gp)
    // {return  new NonlocalMaterialStatus (1,this->giveDomain(), gp);;}
    /**
     * Notifies the receiver, that the nonlocal averaging has been finished for given ip.
     * It deletes IP nonlocal table if permanentNonlocTableFlag is flase.
     * This can save significat memory, since nonlocal tables are not stored, but every time computed when needed,
     * but on the other hand computational time may significantly grow.
     */
    void endIPNonlocalAverage (GaussPoint*gp);

protected:
    /**
     * Returns true if the barrier is activated
     * by interaction of two given points. In this case the nonlocal influence
     * is not considered. Otherwise returns false.
     * @param c1 coordinates of first point
     * @param c2 coordinates of second point
     * @return true if barrier is activated, false otherwise
     */
    //bool isBarrierActivated (const FloatArray& c1, const FloatArray& c2) const;
    void applyBarrierConstraints(const FloatArray &gpCoords, const FloatArray &jGpCoords, double &weight);
};

#endif // nonlocalmaterialext_h


