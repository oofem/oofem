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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "xfem/enrichmentfronts/enrichmentfront.h"
#include "floatarray.h"

#ifndef ENRICHMENTFRONTINTERSECTION_H_
 #define ENRICHMENTFRONTINTERSECTION_H_

 #define _IFT_EnrFrontIntersection_Name "enrfrontintersection"
 #define _IFT_EnrFrontIntersection_Tangent "tangent"

namespace oofem {
class XfemManager;
class DofManager;
class FloatArray;
class InputRecord;
class DynamicInputRecord;
class LinElBranchFunction;

/**
 * EnrFrontIntersection
 *
 * An enrichment front capable of handling crack intersections.
 *
 * @author Erik Svenning
 * @date Apr 23, 2014
 */
class EnrFrontIntersection : public EnrichmentFront
{
public:
    EnrFrontIntersection();
    virtual ~EnrFrontIntersection();

    void MarkNodesAsFront(std :: unordered_map< int, NodeEnrichmentType > &ioNodeEnrMarkerMap, XfemManager &ixFemMan,  const std :: unordered_map< int, double > &iLevelSetNormalDirMap, const std :: unordered_map< int, double > &iLevelSetTangDirMap, const TipInfo &iTipInfo) override;

    int giveNumEnrichments(const DofManager &iDMan) const override;
    int giveMaxNumEnrichments() const override { return 1; }

    // Evaluate the enrichment function and its derivative in front nodes.
    void evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const EfInput &iEfInput) const override;
    void evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const EfInput &iEfInput, const FloatArray &iGradLevelSet) const override;
    void evaluateEnrFuncJumps(std :: vector< double > &oEnrFuncJumps, GaussPoint &iGP, int iNodeInd, bool iGPLivesOnCurrentCrack, const double &iNormalSignDist) const override;

    const char *giveClassName() const override { return "EnrFrontIntersection"; }
    const char *giveInputRecordName() const override { return _IFT_EnrFrontIntersection_Name; }

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    double giveSupportRadius() const override { return 0.0; }
    bool propagationIsAllowed() const override { return false; }

    void setTangent(FloatArray iTangent) { mTangent = std :: move(iTangent); }

protected:
    FloatArray mTangent;
};
} /* namespace oofem */

#endif /* ENRICHMENTFRONTINTERSECTION_H_ */
