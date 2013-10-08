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
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
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
#define _IFT_DofManList_DelaminationLevel "xi"
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
    virtual void giveInputRecord(DynamicInputRecord &input) = 0;
    int number; // remove - JB
    int giveNumber() { return number; }; // remove - JB
    void setNumber(int i) { this->number = i; }; // remove - JB

    // Update of description
    virtual void updateEnrichmentDomain(){};

    virtual const char *giveInputRecordName() const = 0;
    virtual const char *giveClassName() const = 0;


    /// Functions for computing signed distance in normal and tangential direction.
    /// Used by XFEM level set functions.
    virtual void computeNormalSignDist(double &oDist, const FloatArray &iPoint) const = 0;
    virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint) const = 0;
    virtual void computeSurfaceNormalSignDist(double &oDist, const FloatArray &iPoint) const = 0;

    // Use double dispatch to call the correct version of CallNodeEnrMarkerUpdate.
    virtual void CallNodeEnrMarkerUpdate(EnrichmentItem &iEnrItem, XfemManager &ixFemMan) const {}


    virtual bool giveClosestTipInfo(const FloatArray &iCoords, TipInfo &oInfo) const {return false;}

    /// Return array with info about all tips
    virtual bool giveTipInfos(std::vector<TipInfo> &oInfo) const {return false;}

    /// Propagate tips
    virtual bool propagateTips(const std::vector<TipPropagation> &iTipProp) {return false;}
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
    virtual void giveInputRecord(DynamicInputRecord &input);

    /**
     * Functions for computing signed distance in normal and tangential direction.
     * Used by XFEM level set functions.
     */
    virtual void computeNormalSignDist(double &oDist, const FloatArray &iPoint) const { bg->computeNormalSignDist(oDist, iPoint); };
    virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint) const { bg->computeTangentialSignDist(oDist, iPoint); };
    virtual void computeSurfaceNormalSignDist(double &oDist, const FloatArray &iPoint) const { OOFEM_ERROR("EnrichmentDomain_BG::computeNormalSignDist -- not implemented"); };
    // Use double dispatch to call the correct version of CallNodeEnrMarkerUpdate.
    virtual void CallNodeEnrMarkerUpdate(EnrichmentItem &iEnrItem, XfemManager &ixFemMan) const;
};


class EDBGCircle : public EnrichmentDomain_BG
{
public:
    EDBGCircle() { bg = new Circle; };
    virtual ~EDBGCircle() { delete bg; }

    virtual IRResultType initializeFrom(InputRecord *ir) { return bg->initializeFrom(ir); }

    virtual const char *giveInputRecordName() const { return _IFT_EDBGCircle_Name; }
    virtual const char *giveClassName() const { return "EDBGCircle"; }
};

/**
 * EDCrack: Enrichment geometry described by a piecewise linear polygon.
 */
class EDCrack : public EnrichmentDomain_BG
{
public:
    EDCrack() { bg = new PolygonLine; }
    virtual ~EDCrack() { delete bg; }

    virtual IRResultType initializeFrom(InputRecord *ir) { return bg->initializeFrom(ir); }

    virtual const char *giveInputRecordName() const { return _IFT_EDCrack_Name; }
    virtual const char *giveClassName() const { return "EDCrack"; }

    virtual bool giveClosestTipInfo(const FloatArray &iCoords, TipInfo &oInfo) const;
    virtual bool giveTipInfos(std::vector<TipInfo> &oInfo) const;
    virtual bool propagateTips(const std::vector<TipPropagation> &iTipProp);

};


/**
 * List of DofManagers
 * ///@todo: Add additional basic geometry descriptions like polygon
 */
class DofManList : public EnrichmentDomain
{
protected:
    std::vector< int > dofManList;
    double xi;
public:
    DofManList() { }
    virtual ~DofManList() { }

    const std :: vector< int > &giveDofManList() const { return dofManList; }

    //virtual void computeNormalSignDist(double &oDist, const FloatArray &iPoint) const { OOFEM_ERROR("DofManList::computeNormalSignDist -- not implemented"); };
    virtual void computeNormalSignDist(double &oDist, const FloatArray &iPoint) const { oDist = 0.0; };
    //virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint) const { OOFEM_ERROR("DofManList::computeTangentialSignDist -- not implemented"); };
    virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint) const { oDist = 0.0;};
    virtual void computeSurfaceNormalSignDist(double &oDist, const FloatArray &iPoint) const; // new /JB
    // Use double dispatch to call the correct version of CallNodeEnrMarkerUpdate.
    virtual void CallNodeEnrMarkerUpdate(EnrichmentItem &iEnrItem, XfemManager &ixFemMan) const;

    virtual IRResultType initializeFrom(InputRecord *ir);
    void addDofManagers(IntArray &dofManNumbers);
    virtual void giveInputRecord(DynamicInputRecord &input);
    virtual void updateEnrichmentDomain(IntArray &dofManNumbers);

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
    virtual void computeSurfaceNormalSignDist(double &oDist, const FloatArray &iPoint) const { OOFEM_ERROR("WholeDomain::computeNormalSignDist -- not implemented"); };
    // Use double dispatch to call the correct version of CallNodeEnrMarkerUpdate.
    virtual void CallNodeEnrMarkerUpdate(EnrichmentItem &iEnrItem, XfemManager &ixFemMan) const;

    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual const char *giveInputRecordName() const { return _IFT_WholeDomain_Name; }
    virtual const char *giveClassName() const { return "WholeDomain"; }
};
} // end namespace oofem
#endif
