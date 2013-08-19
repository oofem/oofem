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
#include "tipinfo.h"

namespace oofem {
class EnrichmentItem;

///@name Input fields for Enrichment domains
//@{
#define _IFT_DofManList_list "list"
#define _IFT_DofManList_Name "dofmanlist"
#define _IFT_WholeDomain_Name "wholedomain"
#define _IFT_EDBGCircle_Name "circle"

#define _IFT_EDCrack_Name "polygoncrack"

//#define _IFT_BasicGeometryDomain<Line>_Name "line" // Odd one out, how should we treat these?
//@}

/**
 * Abstract representation of enrichment domain - the geometry description of the particular
 * enrichment item. Includes BasicGeometry as one type of description, list of enriched dofmanagers etc.
 * Should be extended to handle implicit geometry descriptions like e.g. level-sets.
 * @author Jim Brouzoulis
 * @author Erik Svenning
 */
class EnrichmentDomain
{
public:
    EnrichmentDomain();
    virtual ~EnrichmentDomain() { }
    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }

    // Spatial search methods
    virtual void computeIntersectionPoints(std :: vector< FloatArray > &oIntersectionPoints, Element *element) { }
    virtual int computeNumberOfIntersectionPoints(Element *element) { return 0; }

    virtual const char *giveInputRecordName() const = 0;
    virtual const char *giveClassName() const = 0;


    /// Functions for computing signed distance in normal and tangential direction.
    /// Used by XFEM level set functions.
    virtual void computeNormalSignDist(double &oDist, const FloatArray &iPoint) const = 0;
    virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint) const = 0;

    // Use double dispatch to call the correct version of CallNodeEnrMarkerUpdate.
    virtual void CallNodeEnrMarkerUpdate(EnrichmentItem &iEnrItem, XfemManager &ixFemMan) {}

    virtual bool GiveClosestTipPosition(FloatArray &oCoords, const FloatArray &iCoords) const {oCoords.setValues(2, 0.0, 0.0); return false;}

    virtual bool GiveClosestTipInfo(const FloatArray &iCoords, TipInfo &oInfo) const {return false;}
};


/**
 * Base class for EnrichmentDomains that derive from BasicGeometry
 * @todo: Add additional basic geometry descriptions like polygon
 */
class EnrichmentDomain_BG : public EnrichmentDomain
{
public:
    BasicGeometry *bg;
    EnrichmentDomain_BG() { }
    virtual ~EnrichmentDomain_BG() { }
    virtual IRResultType initializeFrom(InputRecord *ir) { return this->bg->initializeFrom(ir); }

    virtual void computeIntersectionPoints(std :: vector< FloatArray > &oIntersectionPoints, Element *element) { printf("Entering BasicGeometry::computeIntersectionPoints()\n"); bg->computeIntersectionPoints(element, oIntersectionPoints); }
    virtual int computeNumberOfIntersectionPoints(Element *element) { return bg->computeNumberOfIntersectionPoints(element); }


    /// Functions for computing signed distance in normal and tangential direction.
    /// Used by XFEM level set functions.
    virtual void computeNormalSignDist(double &oDist, const FloatArray &iPoint) const { bg->computeNormalSignDist(oDist, iPoint); };
    virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint) const { bg->computeTangentialSignDist(oDist, iPoint); };

    // Use double dispatch to call the correct version of CallNodeEnrMarkerUpdate.
    virtual void CallNodeEnrMarkerUpdate(EnrichmentItem &iEnrItem, XfemManager &ixFemMan);
};

class EDBGCircle : public EnrichmentDomain_BG
{
public:
    EDBGCircle() { bg = new Circle; };
    virtual ~EDBGCircle() { delete bg; }
    virtual IRResultType initializeFrom(InputRecord *ir) { return bg->initializeFrom(ir); }

    virtual void computeIntersectionPoints(std :: vector< FloatArray > &oIntersectionPoints, Element *element) { bg->computeIntersectionPoints(element, oIntersectionPoints); }
    virtual int computeNumberOfIntersectionPoints(Element *element) { return static_cast< Circle * >( bg )->computeNumberOfIntersectionPoints(element); }

    virtual const char *giveInputRecordName() const { return _IFT_EDBGCircle_Name; }
    virtual const char *giveClassName() const { return "EDBGCircle"; }
};

class EDCrack : public EnrichmentDomain_BG
{
public:
    EDCrack() { bg = new PolygonLine; }
    virtual ~EDCrack() { delete bg; }
    virtual IRResultType initializeFrom(InputRecord *ir) { return bg->initializeFrom(ir); }

    virtual const char *giveInputRecordName() const { return _IFT_EDCrack_Name; }
    virtual const char *giveClassName() const { return "EDCrack"; }

    virtual bool GiveClosestTipPosition(FloatArray &oCoords, const FloatArray &iCoords) const;
    virtual bool GiveClosestTipInfo(const FloatArray &iCoords, TipInfo &oInfo) const;

};


/**
 * List of DofManagers
 * ///@todo: Add additional basic geometry descriptions like polygon
 */
class DofManList : public EnrichmentDomain
{
protected:
    std::vector< int > dofManList;
public:
    DofManList() { }
    virtual ~DofManList() { }

    const std :: vector< int > &giveDofManList() const { return dofManList; }

    virtual void computeNormalSignDist(double &oDist, const FloatArray &iPoint) const { OOFEM_ERROR("DofManList::computeNormalSignDist -- not implemented"); };
    virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint) const { OOFEM_ERROR("DofManList::computeTangentialSignDist -- not implemented"); };

    // Use double dispatch to call the correct version of CallNodeEnrMarkerUpdate.
    virtual void CallNodeEnrMarkerUpdate(EnrichmentItem &iEnrItem, XfemManager &ixFemMan);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual const char *giveInputRecordName() const { return _IFT_DofManList_Name; }
    virtual const char *giveClassName() const { return "DofManList"; }
};

/**
 * The whole computational domain is enriched which thus is a global enrichment
 * Mostly intended for debugging but may easily lead to a singular problem if the
 * solution is enriched with strong discontinuities.
 */
class WholeDomain : public EnrichmentDomain
{
public:
    WholeDomain() { }
    virtual ~WholeDomain() { }

    virtual void computeNormalSignDist(double &oDist, const FloatArray &iPoint) const { OOFEM_ERROR("WholeDomain::computeNormalSignDist -- not implemented"); };
    virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint) const { OOFEM_ERROR("WholeDomain::computeTangentialSignDist -- not implemented"); };

    // Use double dispatch to call the correct version of CallNodeEnrMarkerUpdate.
    virtual void CallNodeEnrMarkerUpdate(EnrichmentItem &iEnrItem, XfemManager &ixFemMan);

    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }

    virtual const char *giveInputRecordName() const { return _IFT_WholeDomain_Name; }
    virtual const char *giveClassName() const { return "WholeDomain"; }
};
} // end namespace oofem
#endif
