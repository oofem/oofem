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

#ifndef xfemelementinterface_h
#define xfemelementinterface_h

#include "interface.h"
#include "xfemmanager.h"
#include "geometry.h"
#include "matresponsemode.h"
#include "materialmode.h"

#define _IFT_XfemElementInterface_CohesiveZoneMaterial "czmaterial"
#define _IFT_XfemElementInterface_NumIntPointsCZ "nipcz"
#define _IFT_XfemElementInterface_PlaneStrain "useplanestrain"

namespace oofem {
class FloatArray;
class FloatMatrix;
class Triangle;
class Element;
class GaussPoint;
class Element;
class XfemManager;
class IntegrationRule;
class MaterialStatus;

/**
 * Provides Xfem interface for an element.
 * @author Erik Svenning
 */
class OOFEM_EXPORT XfemElementInterface : public Interface
{
public:
    Element *element;


    /// Index of enrichment items associated with cohesive zones
    std :: vector< int >mCZEnrItemIndices; // TODO: Not nice. /ES

    /**
     * Indices of enrichment items that give cohesive zone contributions
     * to a given GP, even though the GP is not located on any
     * of these enrichment items.
     */
    std :: vector< std :: vector< int > >mCZTouchingEnrItemIndices;

    /// Flag that tells if plane stress or plane strain is assumed
    bool mUsePlaneStrain;

    virtual const char *giveClassName() const { return "XfemElementInterface"; }
    std :: string errorInfo(const char *func) const { return std :: string( giveClassName() ) + func; }

public:
    /// Constructor.
    XfemElementInterface(Element *e);

    virtual ~XfemElementInterface();
	XfemElementInterface (const XfemElementInterface& src) = delete;
	XfemElementInterface &operator = (const XfemElementInterface &src) = delete;

    /// Creates enriched B-matrix.
    void XfemElementInterface_createEnrBmatrixAt(FloatMatrix &oAnswer, GaussPoint &iGP, Element &iEl);

    /// Creates enriched BH-matrix.
    void XfemElementInterface_createEnrBHmatrixAt(FloatMatrix &oAnswer, GaussPoint &iGP, Element &iEl);

    /**
     * Help function for computation of B and BH.
     * Avoid duplication of code.
     */
    void ComputeBOrBHMatrix(FloatMatrix &oAnswer, GaussPoint &iGP, Element &iEl, bool iComputeBH);

    /// Creates enriched N-matrix.
    void XfemElementInterface_createEnrNmatrixAt(FloatMatrix &oAnswer, const FloatArray &iLocCoord, Element &iEl, bool iSetDiscontContribToZero);

    /**
     * Creates enriched N-matrix for a chosen subset of element nodes.
     */
    void XfemElementInterface_createEnrNmatrixAt(FloatMatrix &oAnswer, const FloatArray &iLocCoord, Element &iEl, const std :: vector< int > &iLocNodeInd, bool iSetDiscontContribToZero);

    /**
     * Computes total number of enrichments in a node.
     */
    int XfemElementInterface_giveNumDofManEnrichments(const DofManager &iDMan, XfemManager &iXMan) const;

    /// Partitions the element into patches by a triangulation.
    virtual void XfemElementInterface_partitionElement(std :: vector< Triangle > &oTriangles, const std :: vector< FloatArray > &iPoints);
    /// Updates integration rule based on the triangulation.
    virtual bool XfemElementInterface_updateIntegrationRule();

    /// Returns an array of array of points. Each array of points defines the points of a subregion of the element.
    virtual void XfemElementInterface_prepareNodesForDelaunay(std :: vector< std :: vector< FloatArray > > &oPointPartitions, double &oCrackStartXi, double &oCrackEndXi, int iEnrItemIndex, bool &oIntersection);
    virtual void XfemElementInterface_prepareNodesForDelaunay(std :: vector< std :: vector< FloatArray > > &oPointPartitions, double &oCrackStartXi, double &oCrackEndXi, const Triangle &iTri, int iEnrItemIndex, bool &oIntersection);

    // Help functions for partitioning
    void putPointsInCorrectPartition(std :: vector< std :: vector< FloatArray > > &oPointPartitions, const std :: vector< FloatArray > &iIntersecPoints, const std :: vector< const FloatArray * > &iNodeCoord) const;

    /**
     * Partition a boundary segment to account for cracks cutting
     * the boundary. This is a necessary step to evaluate integrals
     * along an edge cut by one or several cracks.
     */
    void partitionEdgeSegment(int iBndIndex, std :: vector< Line > &oSegments, std :: vector< FloatArray > &oIntersectionPoints, const double &iTangDistPadding = 0.0);

    // TODO: Move to XfemStructuralElementInterface
    std :: vector< std :: unique_ptr< IntegrationRule > >mpCZIntegrationRules;

    MaterialMode giveMaterialMode();

    void updateYourselfCZ(TimeStep *tStep);

    void computeDisplacementJump(GaussPoint &iGP, FloatArray &oJump, const FloatArray &iSolVec, const FloatMatrix &iNMatrix);

    /**
     * Compute N-matrix for cohesive zone.
     */
    void computeNCohesive(FloatMatrix &oN, GaussPoint &iGP, int iEnrItemIndex, const std :: vector< int > &iTouchingEnrItemIndices);
};
} // end namespace oofem
#endif // xfemelementinterface_h
