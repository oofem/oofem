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
#include "flotarry.h"
#include "node.h"
#include "contextioresulttype.h"
#include "contextmode.h"


#include "geometry.h"

namespace oofem {
/**
 * Abstract representation of enrichment domain
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


    // What does it need to be able to answer?
    // should be pure virtual and thus be supported by all representations
    //virtual bool isDofManagerEnriched(const DofManager *dMan) { return false; };
    virtual bool isDofManagerEnriched(DofManager *dMan) = 0 {};
    virtual bool isElementEnriched(Element *element) { return false; };

};

class EnrichmentDomain_BG : public EnrichmentDomain
{
private:
    
public:
BasicGeometry *bg;
    EnrichmentDomain_BG(){}; 
    virtual ~EnrichmentDomain_BG() { }
    virtual IRResultType initializeFrom(InputRecord *ir) { return this->bg->initializeFrom(ir); };

    //void setGeometry(BasicGeometry *geom) { this->bg = geom;}
    virtual bool isDofManagerEnriched(DofManager *dMan){ return false; };
   
};

class EDBGCircle : public EnrichmentDomain_BG
{
public:
    EDBGCircle(){ bg = new Circle; }; 
    virtual ~EDBGCircle() { }
    virtual IRResultType initializeFrom(InputRecord *ir) { return this->bg->initializeFrom(ir);  };
    virtual bool isDofManagerEnriched(DofManager *dMan);
    virtual bool isElementEnriched(const Element *element);
   
};


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



// The whole computational domain is enriched.
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



