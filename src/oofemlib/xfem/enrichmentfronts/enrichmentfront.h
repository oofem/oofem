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


#ifndef ENRICHMENTFRONT_H_
#define ENRICHMENTFRONT_H_
#include "oofemcfg.h"
#include <vector>
#include "inputrecord.h"
#include "xfem/tipinfo.h"

#include <unordered_map>

namespace oofem {
class XfemManager;
class DofManager;
class FloatArray;
class InputRecord;
class DynamicInputRecord;
class GaussPoint;
enum NodeEnrichmentType : int;

struct EfInput
{
    EfInput() {}
    EfInput(const FloatArray &iPos, const double &iLevelSet, int iNodeInd, const FloatArray &iClosestPointOnCrack, const double &iArcPos, const FloatArray &iLocalTangDir) :
        mPos(iPos),
        mLevelSet(iLevelSet),
        mNodeInd(iNodeInd),
        mClosestPointOnCrack(iClosestPointOnCrack),
        mArcPos(iArcPos),
        mLocalTangDir(iLocalTangDir)
    {}

    ~EfInput() {}

    FloatArray mPos;
    double mLevelSet;
    int mNodeInd;
    FloatArray mClosestPointOnCrack;
    double mArcPos;
    FloatArray mLocalTangDir;
};


/**
 * Class EnrichmentFront: describes the edge or tip of an XFEM enrichment.
 * The purpose is to add a different treatment of the front than the "interior"
 * enrichments. We may, e.g.
 *  - Apply branch functions at a crack tip for the element containing the crack tip.
 *  - Apply branch functions on nodes within a certain radius from the crack tip.
 *  - Exclude nodes touched by the front.
 *
 *  The desired behavior is obtained by choosing a suitable EnrichmentFront.
 *
 * @author Erik Svenning
 * @date Feb 14, 2014
 */
class OOFEM_EXPORT EnrichmentFront
{
public:
    EnrichmentFront() { }
    virtual ~EnrichmentFront() { }

    /**
     *  MarkNodesAsFront:
     *  Intput:
     *  -ioNodeEnrMarker:   A vector with the same size as the number of nodes in the mesh
     *                      where the nodes corresponding to interior XFEM enrichments are
     *                      marked with 1, other entries are zero.
     *
     *  Output:
     *  -ioNodeEnrMarker:   Modifies the vector by marking tip nodes as 2, meaning that they
     *                      should get special treatment. May also modify the set of nodes
     *                      enriched by the interior enrichment.
     */
    virtual void MarkNodesAsFront(std :: unordered_map< int, NodeEnrichmentType > &ioNodeEnrMarkerMap, XfemManager &ixFemMan, const std :: unordered_map< int, double > &iLevelSetNormalDirMap, const std :: unordered_map< int, double > &iLevelSetTangDirMap, const TipInfo &iTipInfo) = 0;

    // The number of enrichment functions applied to tip nodes.
    virtual int  giveNumEnrichments(const DofManager &iDMan) const = 0;
    virtual int  giveMaxNumEnrichments() const = 0;


    // Evaluate the enrichment function and its derivative in front nodes.
    virtual void evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const EfInput &iEfInput) const = 0;
    virtual void evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const EfInput &iEfInput, const FloatArray &iGradLevelSet) const = 0;
    virtual void evaluateEnrFuncJumps(std :: vector< double > &oEnrFuncJumps, GaussPoint &iGP, int iNodeInd, bool iGPLivesOnCurrentCrack, const double &iNormalSignDist) const = 0;

    std :: string errorInfo(const char *func) const { return std :: string( giveClassName() ) + func; }

    virtual const char *giveClassName() const = 0;
    virtual const char *giveInputRecordName() const = 0;

    virtual IRResultType initializeFrom(InputRecord *ir) = 0;
    virtual void giveInputRecord(DynamicInputRecord &input) = 0;

    virtual double giveSupportRadius() const = 0;

    virtual bool propagationIsAllowed() const { return true; }

    const TipInfo &giveTipInfo() const { return mTipInfo; }

    void computeCrackTangent(FloatArray &oTangent, FloatArray &oNormal, bool &oFlipTangent, const EfInput &iEfInput) const;

protected:
    TipInfo mTipInfo;

    /**
     * Several enrichment fronts enrich all nodes in the tip element.
     * This help function accomplishes that.
     */
    void MarkTipElementNodesAsFront(std :: unordered_map< int, NodeEnrichmentType > &ioNodeEnrMarkerMap, XfemManager &ixFemMan,  const std :: unordered_map< int, double > &iLevelSetNormalDirMap, const std :: unordered_map< int, double > &iLevelSetTangDirMap, const TipInfo &iTipInfo);
};
} // end namespace oofem


#endif /* ENRICHMENTFRONT_H_ */
