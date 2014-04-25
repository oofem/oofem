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

#include "oofemcfg.h"
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
//@}

/**
 * Abstract representation of enrichment domain - the geometry description of the particular
 * enrichment item. Includes BasicGeometry as one type of description, list of enriched dofmanagers etc.
 * Should be extended to handle implicit geometry descriptions like e.g. level-sets.
 * @author Jim Brouzoulis
 * @author Erik Svenning
 */
class OOFEM_EXPORT EnrichmentDomain
{
public:
    EnrichmentDomain();
    virtual ~EnrichmentDomain() { }
    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }
    virtual void giveInputRecord(DynamicInputRecord &input) = 0;

    virtual const char *giveInputRecordName() const = 0;
    virtual const char *giveClassName() const = 0;


    /// Functions for computing signed distance in normal and tangential direction.
    /// Used by XFEM level set functions.
    virtual void computeNormalSignDist(double &oDist, const FloatArray &iPoint) const = 0;
    virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint, double &oMinDistArcPos) const = 0;
    virtual void computeSurfaceNormalSignDist(double &oDist, const FloatArray &iPoint) const = 0;
    virtual void giveSubPolygon(std :: vector< FloatArray > &oPoints, const double &iXiStart, const double &iXiEnd) const { OOFEM_ERROR("not implemented."); }

    // Use double dispatch to call the correct version of CallNodeEnrMarkerUpdate.
    virtual void CallNodeEnrMarkerUpdate(EnrichmentItem &iEnrItem, XfemManager &ixFemMan) const { }

    // Return the center and radius of a sphere containing all nodes
    // enriched by the enrichment item.
    virtual void giveBoundingSphere(FloatArray &oCenter, double &oRadius) = 0;

    virtual bool giveClosestTipInfo(const FloatArray &iCoords, TipInfo &oInfo) const { return false; }

    /// Return array with info about all tips
    virtual bool giveTipInfos(TipInfo &oStartTipInfo, TipInfo &oEndTipInfo) const { return false; }

    /// Propagate tips
    virtual bool propagateTips(const std :: vector< TipPropagation > &iTipProp) { return false; }

    void setVtkDebug(bool iDebugVTK) { mDebugVTK = iDebugVTK; }
    bool getVtkDebug() const { return mDebugVTK; }

protected:
    bool mDebugVTK;
};


/**
 * Base class for EnrichmentDomains that derive from BasicGeometry
 * @todo: Add additional basic geometry descriptions like polygon
 */
class OOFEM_EXPORT EnrichmentDomain_BG : public EnrichmentDomain
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
    virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint, double &oMinDistArcPos) const { bg->computeTangentialSignDist(oDist, iPoint, oMinDistArcPos); };
    virtual void computeSurfaceNormalSignDist(double &oDist, const FloatArray &iPoint) const { OOFEM_ERROR("not implemented"); };
    virtual void giveSubPolygon(std :: vector< FloatArray > &oPoints, const double &iXiStart, const double &iXiEnd) const { bg->giveSubPolygon(oPoints, iXiStart, iXiEnd); }

    // Use double dispatch to call the correct version of CallNodeEnrMarkerUpdate.
    virtual void CallNodeEnrMarkerUpdate(EnrichmentItem &iEnrItem, XfemManager &ixFemMan) const;

    virtual void giveBoundingSphere(FloatArray &oCenter, double &oRadius);
};

class OOFEM_EXPORT EDBGCircle : public EnrichmentDomain_BG
{
public:
    EDBGCircle() {
        bg = new Circle;
    };
    virtual ~EDBGCircle() {
        delete bg;
    }

    virtual IRResultType initializeFrom(InputRecord *ir) { return bg->initializeFrom(ir); }

    virtual const char *giveInputRecordName() const { return _IFT_EDBGCircle_Name; }
    virtual const char *giveClassName() const { return "EDBGCircle"; }
    std :: string errorInfo(const char *func) const { return std :: string( giveClassName() ) + func; }

    virtual void giveBoundingSphere(FloatArray &oCenter, double &oRadius);
};

/**
 * EDCrack: Enrichment geometry described by a piecewise linear polygon.
 */
class OOFEM_EXPORT EDCrack : public EnrichmentDomain_BG
{
public:
    EDCrack() {
        bg = new PolygonLine;
    }
    virtual ~EDCrack() {
        delete bg;
    }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual const char *giveInputRecordName() const { return _IFT_EDCrack_Name; }
    virtual const char *giveClassName() const { return "EDCrack"; }

    virtual bool giveClosestTipInfo(const FloatArray &iCoords, TipInfo &oInfo) const;
    virtual bool giveTipInfos(TipInfo &oStartTipInfo, TipInfo &oEndTipInfo) const;
    virtual bool propagateTips(const std :: vector< TipPropagation > &iTipProp);
};


/**
 * List of DofManagers
 */
class OOFEM_EXPORT DofManList : public EnrichmentDomain
{
protected:
    std :: vector< int >dofManList;
    double xi;
public:
    DofManList() : xi(0.0) { }
    virtual ~DofManList() { }

    const std :: vector< int > &giveDofManList() const { return dofManList; }

    virtual void computeNormalSignDist(double &oDist, const FloatArray &iPoint) const { oDist = 0.0; };
    virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint, double &oMinDistArcPos) const { oDist = 0.0; };
    virtual void computeSurfaceNormalSignDist(double &oDist, const FloatArray &iPoint) const; // new /JB
    // Use double dispatch to call the correct version of CallNodeEnrMarkerUpdate.
    virtual void CallNodeEnrMarkerUpdate(EnrichmentItem &iEnrItem, XfemManager &ixFemMan) const;

    virtual IRResultType initializeFrom(InputRecord *ir);
    void addDofManagers(IntArray &dofManNumbers);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual const char *giveInputRecordName() const { return _IFT_DofManList_Name; }
    virtual const char *giveClassName() const { return "DofManList"; }

    virtual void giveBoundingSphere(FloatArray &oCenter, double &oRadius);
};

/**
 * The whole computational domain is enriched which thus is a global enrichment
 * Mostly intended for debugging but may easily lead to a singular problem if the
 * solution is enriched with strong discontinuities.
 */
class OOFEM_EXPORT WholeDomain : public EnrichmentDomain
{
public:
    WholeDomain() { }
    virtual ~WholeDomain() { }

    virtual void computeNormalSignDist(double &oDist, const FloatArray &iPoint) const { OOFEM_ERROR("not implemented"); };
    virtual void computeTangentialSignDist(double &oDist, const FloatArray &iPoint, double &oMinDistArcPos) const { OOFEM_ERROR("not implemented"); };
    virtual void computeSurfaceNormalSignDist(double &oDist, const FloatArray &iPoint) const { OOFEM_ERROR("not implemented"); };
    // Use double dispatch to call the correct version of CallNodeEnrMarkerUpdate.
    virtual void CallNodeEnrMarkerUpdate(EnrichmentItem &iEnrItem, XfemManager &ixFemMan) const;

    virtual IRResultType initializeFrom(InputRecord *ir) { return IRRT_OK; }
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual const char *giveInputRecordName() const { return _IFT_WholeDomain_Name; }
    virtual const char *giveClassName() const { return "WholeDomain"; }

    virtual void giveBoundingSphere(FloatArray &oCenter, double &oRadius);
};
} // end namespace oofem
#endif
