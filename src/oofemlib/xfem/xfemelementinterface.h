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
#include "alist.h"
#include "xfemmanager.h"

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
class StructuralInterfaceMaterial;

//#define XFEM_DEBUG_VTK 1

/**
 * Provides Xfem interface for an element.
 * @author Erik Svenning
 */
class OOFEM_EXPORT XfemElementInterface : public Interface
{
public:
    Element *element;

    // Cohesive Zone variables
    StructuralInterfaceMaterial *mpCZMat;
    int mCZMaterialNum;
    int mCSNumGaussPoints;
    std::vector<IntegrationRule*> mpCZIntegrationRules;

    /// Index of enrichment items associated with cohesive zones
    std::vector<int> mCZEnrItemIndices; // TODO: Not nice. /ES

    /// Flag that tells if plane stress or plane strain is assumed
    bool mUsePlaneStrain;

public:
    /// Constructor.
    XfemElementInterface(Element *e);

    virtual ~XfemElementInterface();

    /// Creates enriched B-matrix.
    void XfemElementInterface_createEnrBmatrixAt(FloatMatrix &oAnswer, GaussPoint &iGP, Element &iEl);

    /// Creates enriched N-matrix.
    void XfemElementInterface_createEnrNmatrixAt(FloatMatrix &oAnswer, const FloatArray &iLocCoord, Element &iEl);

    /// Partitions the element into patches by a triangulation.
    virtual void XfemElementInterface_partitionElement(std :: vector< Triangle > &oTriangles, const std :: vector< FloatArray > &iPoints);
    /// Updates integration rule based on the triangulation.
    virtual bool XfemElementInterface_updateIntegrationRule();

    /// Returns an array of array of points. Each array of points defines the points of a subregion of the element.
    virtual void XfemElementInterface_prepareNodesForDelaunay(std :: vector< std :: vector< FloatArray > > &oPointPartitions, double &oCrackStartXi, double &oCrackEndXi, int iEnrItemIndex, bool &oIntersection);
    virtual void XfemElementInterface_prepareNodesForDelaunay(std :: vector< std :: vector< FloatArray > > &oPointPartitions, double &oCrackStartXi, double &oCrackEndXi, const Triangle &iTri, int iEnrItemIndex, bool &oIntersection);

    // Help functions for partitioning
    void putPointsInCorrectPartition(std :: vector< std :: vector< FloatArray > > &oPointPartitions, const std :: vector< FloatArray > &iIntersecPoints, const std::vector<const FloatArray*> &iNodeCoord) const;

    /**
     * XfemElementInterface_computeConstitutiveMatrixAt.
     * The reason for having a special implementation of this function is
     * that the enrichment item may have a different material than
     * the bulk material inside.
     */
    virtual void XfemElementInterface_computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *, TimeStep *tStep);
    virtual void XfemElementInterface_computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *stepN);

    /**
     * Cohesive Zone functions
     */
    bool hasCohesiveZone() const { return ( mpCZMat != NULL && mpCZIntegrationRules.size() > 0 ); }

    void computeCohesiveForces(FloatArray &answer, TimeStep *tStep);
    void computeGlobalCohesiveTractionVector(FloatArray &oT, const FloatArray &iJump, const FloatArray &iCrackNormal, const FloatMatrix &iNMatrix, GaussPoint &iGP, TimeStep *tStep);

    void computeCohesiveTangent(FloatMatrix &answer, TimeStep *tStep);
    void computeCohesiveTangentAt(FloatMatrix &answer, TimeStep *tStep);

    void XfemElementInterface_computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass, const double *ipDensity = NULL);

    virtual IRResultType initializeCZFrom(InputRecord *ir);
    MaterialMode giveMaterialMode();
    virtual void giveCZInputRecord(DynamicInputRecord &input);

    void initializeCZMaterial();

    void updateYourselfCZ(TimeStep *tStep);

    void computeDisplacementJump(GaussPoint &iGP, FloatArray &oJump, const FloatArray &iSolVec, const FloatMatrix &iNMatrix);

    /**
     * Compute N-matrix for cohesive zone.
     */
    void computeNCohesive(FloatMatrix &oN, GaussPoint &iGP, int iEnrItemIndex);

};
} // end namespace oofem
#endif // xfemelementinterface_h
