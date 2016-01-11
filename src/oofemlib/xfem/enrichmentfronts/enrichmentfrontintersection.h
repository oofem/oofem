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

    virtual void MarkNodesAsFront(std :: unordered_map< int, NodeEnrichmentType > &ioNodeEnrMarkerMap, XfemManager &ixFemMan,  const std :: unordered_map< int, double > &iLevelSetNormalDirMap, const std :: unordered_map< int, double > &iLevelSetTangDirMap, const TipInfo &iTipInfo);

    virtual int  giveNumEnrichments(const DofManager &iDMan) const;
    virtual int  giveMaxNumEnrichments() const { return 1; }

    // Evaluate the enrichment function and its derivative in front nodes.
    virtual void evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const EfInput &iEfInput) const;
    virtual void evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const EfInput &iEfInput, const FloatArray &iGradLevelSet) const;
    virtual void evaluateEnrFuncJumps(std :: vector< double > &oEnrFuncJumps, GaussPoint &iGP, int iNodeInd, bool iGPLivesOnCurrentCrack, const double &iNormalSignDist) const;

    virtual const char *giveClassName() const { return "EnrFrontIntersection"; }
    virtual const char *giveInputRecordName() const { return _IFT_EnrFrontIntersection_Name; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual double giveSupportRadius() const { return 0.0; }
    virtual bool propagationIsAllowed() const { return false; }

    void setTangent(FloatArray iTangent) { mTangent = std :: move(iTangent); }

protected:
    FloatArray mTangent;
};
} /* namespace oofem */

#endif /* ENRICHMENTFRONTINTERSECTION_H_ */
