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

namespace oofem {

class XfemManager;
class TipInfo;
class DofManager;
class FloatArray;
class InputRecord;
class DynamicInputRecord;
/*
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
    EnrichmentFront() { };
    virtual ~EnrichmentFront() { };

    /*
     *  MarkNodesAsFront:
     *  Intput:
     *  -ioNodeEnrMarker:   A vector with the same size as the number of nodes in the mesh
     *                      where the nodes corresponding to interior XFEM enrichments are
     *                      marked with 1, other entries are zero.
     *
     *  Output:
     *  -ioNodeEnrMarker:	Modifies the vector by marking tip nodes as 2, meaning that they
     *                      should get special treatment. May also modify the set of nodes
     *                      enriched by the interior enrichment.
     */
    virtual void MarkNodesAsFront(std :: vector< int > &ioNodeEnrMarker, XfemManager &ixFemMan, const std :: vector< double > &iLevelSetNormalDir, const std :: vector< double > &iLevelSetTangDir, const std :: vector< TipInfo > &iTipInfo) = 0;

    // The number of enrichment functions applied to tip nodes.
    virtual int  giveNumEnrichments(const DofManager &iDMan) const = 0;
    virtual int  giveMaxNumEnrichments() const = 0;


    // Evaluate the enrichment function and its derivative in front nodes.
    virtual void evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iPos, const double &iLevelSet, int iNodeInd) const = 0;
    virtual void evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const FloatArray &iPos, const double &iLevelSet, const FloatArray &iGradLevelSet, int iNodeInd) const = 0;
    virtual void evaluateEnrFuncJumps(std :: vector< double > &oEnrFuncJumps) const = 0;

    virtual const char *giveClassName() const = 0;
    virtual const char *giveInputRecordName() const = 0;

    virtual IRResultType initializeFrom(InputRecord *ir) = 0;
    virtual void giveInputRecord(DynamicInputRecord &input) = 0;

    virtual bool giveElementTipCoord(FloatArray &oCoord, double &oArcPos, int iElIndex) const;

protected:
    std :: vector< TipInfo >mTipInfo;

    /**
     * Keep record of the tips associated with an enriched node:
     * pair.first -> node index
     * pair.second-> tip indices
     */
    std :: vector< std :: pair< int, std :: vector< int > > >mNodeTipIndices;

    void addTipIndexToNode(int iNodeInd, int iTipInd); // Help function for updating mNodeTipIndices
    void giveNodeTipIndices(int iNodeInd, std :: vector< int > &oTipIndices) const;
};

} // end namespace oofem


#endif /* ENRICHMENTFRONT_H_ */
