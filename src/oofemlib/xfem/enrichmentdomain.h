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

#ifndef enrichmentdomain_h
#define enrichmentdomain_h

#include "domain.h"
#include "floatarray.h"
#include "node.h"
#include "contextioresulttype.h"
#include "contextmode.h"
#include "geometry.h"

namespace oofem {

///@name Input fields for Enrichment domains
//@{
#define _IFT_DofManList_list "list"
#define _IFT_DofManList_Name "dofmanlist"
#define _IFT_WholeDomain_Name "wholedomain"
#define _IFT_EDBGCircle_Name "circle"
//#define _IFT_BasicGeometryDomain<Line>_Name "line" // Odd one out, how should we treat these?
//@}

/**
 * Abstract representation of enrichment domain - the geometry description of the particular 
 * enrichment item. Includes BasicGeometry as one type of description, list of enriched dofmanagers etc.
 * Should be extended to handle implicit geometry descriptions like e.g. level-sets. 
 * @author Jim Brouzoulis
 */
class EnrichmentDomain 
{
public:
    EnrichmentDomain(){};
    virtual ~EnrichmentDomain(){};
    virtual IRResultType initializeFrom(InputRecord *ir){ return IRRT_OK; }; 
    virtual const char *giveClassName() const { return NULL; }
    virtual classType giveClassID() const { return EnrichmentDomainClass; }

    virtual bool isDofManagerEnriched(DofManager *dMan) = 0;
    // Default is to loop through the dofman and check if any of them are enriched
    virtual bool isElementEnriched(Element *element); 

    // Spatial search metohds
    virtual void computeIntersectionPoints(AList< FloatArray > *intersectionPoints, Element *element){};
    virtual int computeNumberOfIntersectionPoints(Element *element){return 0;};
};


/**
 * Base class for EnrichmentDomains that derive from BasicGeometry
 * ///@todo: Add additional basic geometry descriptions like polygon
 */
class EnrichmentDomain_BG : public EnrichmentDomain
{
public:
    BasicGeometry *bg;
    EnrichmentDomain_BG(){}; 
    virtual ~EnrichmentDomain_BG() { }
    virtual IRResultType initializeFrom(InputRecord *ir) { return this->bg->initializeFrom(ir); };
    virtual bool isDofManagerEnriched(DofManager *dMan){ return false; };

    virtual void computeIntersectionPoints(AList< FloatArray > *intersectionPoints, Element *element) { bg->computeIntersectionPoints(element, intersectionPoints); }
    virtual int computeNumberOfIntersectionPoints(Element *element) { return bg->computeNumberOfIntersectionPoints(element); };

};

class EDBGCircle : public EnrichmentDomain_BG
{
public:
    EDBGCircle(){ bg = new Circle; }; 
    virtual ~EDBGCircle() { }
    virtual IRResultType initializeFrom(InputRecord *ir) { return bg->initializeFrom(ir);  };
    virtual bool isDofManagerEnriched(DofManager *dMan);
    virtual bool isElementEnriched(Element *element);
    virtual void computeIntersectionPoints(AList< FloatArray > *intersectionPoints, Element *element) { bg->computeIntersectionPoints(element, intersectionPoints); }
    virtual int computeNumberOfIntersectionPoints(Element *element) { return static_cast<Circle *>(bg)->computeNumberOfIntersectionPoints(element); };
};

/**
 * List of DofManagers 
 * ///@todo: Add additional basic geometry descriptions like polygon
 */
class DofManList : public EnrichmentDomain
{
protected:
    std::list< int > dofManList;
public:
    DofManList(){ }
    virtual ~DofManList(){};
    virtual IRResultType initializeFrom(InputRecord *ir) ;
    virtual bool isDofManagerEnriched(DofManager *dMan);
};



/**
 * The whole computational domain is enriched which thus is a global enrichment
 * Mostly intended for debugging but may easily lead to a singular problem if the
 * solution is enriched with strong discontinuities.
 */
class WholeDomain : public EnrichmentDomain
{
public:
    WholeDomain(){ }
    virtual ~WholeDomain(){};
    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; } ;
    virtual bool isDofManagerEnriched(DofManager *dMan) { return true; };
    virtual bool isElementEnriched(Element *element) { return true; };
};

} // end namespace oofem
#endif  



