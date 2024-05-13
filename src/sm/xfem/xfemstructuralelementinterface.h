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

#ifndef XFEMSTRUCTURALELEMENTINTERFACE_H_
#define XFEMSTRUCTURALELEMENTINTERFACE_H_

#include "xfem/xfemelementinterface.h"
#include "internalstatetype.h"
namespace oofem {
class Material;
class IntegrationRule;
class ExportRegion;
class StructuralFE2MaterialStatus;
/**
 * Provides Xfem interface for a structural element.
 * @author Erik Svenning
 * @date Feb 14, 2014
 */
class OOFEM_EXPORT XfemStructuralElementInterface : public XfemElementInterface
{
public:
    XfemStructuralElementInterface(Element *e);

    /// Updates integration rule based on the triangulation.
    bool XfemElementInterface_updateIntegrationRule() override;

    MaterialStatus *giveClosestGP_MatStat(double &oClosestDist, std :: vector< std :: unique_ptr< IntegrationRule > > &iRules, const FloatArray &iCoord);

    double computeEffectiveSveSize(StructuralFE2MaterialStatus *iFe2Ms);

    virtual void XfemElementInterface_computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *, TimeStep *tStep);
    virtual void XfemElementInterface_computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);

    virtual bool hasCohesiveZone() const { return ( mpCZMat != nullptr && mpCZIntegrationRules.size() > 0 ); }

    virtual void computeCohesiveForces(FloatArray &answer, TimeStep *tStep);
    virtual void computeGlobalCohesiveTractionVector(FloatArray &oT, const FloatArray &iJump, const FloatArray &iCrackNormal, const FloatMatrix &iNMatrix, GaussPoint &iGP, TimeStep *tStep);

    virtual void computeCohesiveTangent(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeCohesiveTangentAt(FloatMatrix &answer, TimeStep *tStep);

    virtual void XfemElementInterface_computeConsistentMassMatrix(FloatMatrix &answer, TimeStep *tStep, double &mass, const double *ipDensity = NULL);

    virtual void initializeCZFrom(InputRecord &ir);
    virtual void giveCZInputRecord(DynamicInputRecord &input);

    virtual void initializeCZMaterial();

    bool useNonStdCz();

    virtual void XfemElementInterface_computeDeformationGradientVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);


    /**
     * Identify enrichment items with an intersection enrichment
     * front that touches a given enrichment item.
     */
    void giveIntersectionsTouchingCrack(std :: vector< int > &oTouchingEnrItemIndices, const std :: vector< int > &iCandidateIndices, int iEnrItemIndex, XfemManager &iXMan);

    // Cohesive Zone variables
    Material *mpCZMat = nullptr;
    int mCZMaterialNum = -1;
    int mCSNumGaussPoints = 4;

    // Options for nonstandard cohesive zone model
    bool mIncludeBulkJump = true;
    bool mIncludeBulkCorr = true;

    // Store element subdivision for postprocessing
    std :: vector< Triangle >mSubTri;

    /// VTK Interface
    void giveSubtriangulationCompositeExportData(std :: vector< ExportRegion > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep);

    /// Help functions for VTK export.
    void computeIPAverageInTriangle(FloatArray &answer, IntegrationRule *iRule, Element *elem, InternalStateType isType, TimeStep *tStep, const Triangle &iTri);

};
} /* namespace oofem */

#endif /* XFEMSTRUCTURALELEMENTINTERFACE_H_ */
