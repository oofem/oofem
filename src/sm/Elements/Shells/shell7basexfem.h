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

#ifndef Shell7BaseXFEM_h
#define Shell7BaseXFEM_h

#include "sm/Elements/Shells/shell7base.h"
#include "sm/xfem/enrichmentitems/delamination.h"
#include "xfem/xfemelementinterface.h"
#include "fei3dtrquad.h"

///@name Input fields for el
//@{
#define _IFT_Shell7BaseXFEM_CohesiveZoneMaterial "czmaterial"
//@}


namespace oofem {
class FEI3dTrQuad;
class BoundaryLoad;
class EnrichmentItem;
class ShellCrack;
/**
 * This class represent a 7 parameter shell element.
 * Each node has 7 degrees of freedom (displ. vec., director vec., inhomogeneous thickness strain ).
 * Add ref. to paper!
 * @author Jim Brouzoulis
 * @date 2012-11-01
 */
#define _ExportCZ
class Shell7BaseXFEM : public Shell7Base, public XfemElementInterface
{
protected:
    XfemManager *xMan;
    void updateYourself(TimeStep *tStep) override;
    void postInitialize() override;
    void computeOrderingArray(IntArray &orderingArray, IntArray &activeDofsArray,  EnrichmentItem *ei);

    void evalCovarBaseVectorsAt(const FloatArray &lCoords, FloatMatrix &gcon, FloatArray &solVec, TimeStep *tStep) override;
    void discGiveInitialSolutionVector(FloatArray &answer, IntArray &eiDofIdArray); // should be replaced with general function
    void computeDiscGeneralizedStrainVector(FloatArray &dGenEps, const FloatArray &lCoords, EnrichmentItem *ei, TimeStep *tStep);
    void computeDiscSolutionVector(IntArray &dofIdArray , TimeStep *tStep, FloatArray &solVecD);
    void computeInterfaceJumpAt(int interf, FloatArray &lCoords, TimeStep *tStep, FloatArray &answer);

    // Internal forces
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) override;
    void discComputeSectionalForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, FloatArray &solVecD, EnrichmentItem *ei);
    void computeSectionalForcesAt(FloatArray &sectionalForces, IntegrationPoint *ip, Material *mat, TimeStep *tStep, FloatArray &genEps, double zeta);

    double evaluateLevelSet(const FloatArray &lCoords, EnrichmentItem *ei);
//    double edgeEvaluateLevelSet(const FloatArray &lCoords, EnrichmentItem *ei, const int edge);
    double evaluateHeavisideXi(double xi, ShellCrack *ei);
    double evaluateHeavisideXi(double xi, Delamination *ei);
    double evaluateCutHeaviside(const double xi, const double xiBottom, const double xiTop) const;
    void computeCohesiveForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, FloatArray &solVecD, EnrichmentItem *ei, EnrichmentItem *coupledToEi);

    // Tangent matrices
    void computeLambdaGMatricesDis(FloatMatrix lambdaD [ 3 ], double zeta);
    void computeLambdaNMatrixDis(FloatMatrix &lambda_xd, double zeta);  
    virtual void OLDcomputeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;

    virtual void discComputeBulkTangentMatrix(FloatMatrix &KdIJ, IntegrationPoint *ip, EnrichmentItem *eiI, EnrichmentItem *eiJ, int layer, FloatMatrix A [ 3 ] [ 3 ], TimeStep *tStep);
    virtual void discComputeStiffness(FloatMatrix &LCC, FloatMatrix &LDD, FloatMatrix &LDC, IntegrationPoint *ip, int layer, FloatMatrix A [ 3 ] [ 3 ], TimeStep *tStep);

    double EvaluateEnrFuncInDofMan(int dofManNum, EnrichmentItem *ei);
    void computeEnrichedBmatrixAt(const FloatArray &lCoords, FloatMatrix &answer, EnrichmentItem *ei);
    void computeEnrichedNmatrixAt(const FloatArray &lCoords, FloatMatrix &answer, EnrichmentItem *ei);
    void computeCohesiveNmatrixAt(const FloatArray &lCoords, FloatMatrix &answer, EnrichmentItem *ei);

    void edgeGiveUpdatedSolutionVector(FloatArray &answer, const int iedge, TimeStep *tStep) override;

    void edgeEvalEnrCovarBaseVectorsAt(const FloatArray &lCoords, const int iedge, FloatMatrix &gcov, TimeStep *tStep, EnrichmentItem *ei);
    void computeCohesiveTangent(FloatMatrix &answer, TimeStep *tStep);
    void computeCohesiveTangentAt(FloatMatrix &answer, TimeStep *tStep, Delamination *dei, EnrichmentItem *ei_j, EnrichmentItem *ei_k);

    void edgeComputeEnrichedNmatrixAt(const FloatArray &lCoords, FloatMatrix &answer, EnrichmentItem *ei, const int edge);
    void edgeComputeEnrichedBmatrixAt(const FloatArray &lCoords, FloatMatrix &answer, EnrichmentItem *ei, const int edge);
    
    void computePressureTangentMatrixDis(FloatMatrix &KCC, FloatMatrix &KCD, FloatMatrix &KDD, IntegrationPoint *ip, Load *load, const int iSurf, TimeStep *tStep);

    // External loads
    // Overloaded, as the element is using enhanced approximation 
    void computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep, bool global) override;
    // overloaded, the computeBoundaryEdgeLoadVector returns full element DOFs
    void giveBoundaryLocationArray(IntArray &locationArray, const IntArray &bNodes, const UnknownNumberingScheme &s, IntArray *dofIdArray) override
    { this->giveLocationArray (locationArray, s, dofIdArray);}

    void computeEnrTractionForce(FloatArray &answer, const int iedge, BoundaryLoad *edgeLoad, TimeStep *tStep, 
        ValueModeType mode, EnrichmentItem *ei);

    // Mass matrices
    void computeMassMatrixNum(FloatMatrix &answer, TimeStep *tStep) override;

    // VTK
    void giveCompositeExportData(std::vector < VTKPiece > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep ) override;
    void giveCompositeExportData(VTKPiece &vtkPiece, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep ) override;
    void giveShellExportData(VTKPiece &vtkPiece, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep ) override;
    void giveCZExportData(VTKPiece &vtkPiece, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep );

    void vtkEvalUpdatedGlobalCoordinateAt(const FloatArray &localCoords, int layer, FloatArray &globalCoords, TimeStep *tStep) override;
    void giveDisUnknownsAt(const FloatArray &lCoords, EnrichmentItem *ei, FloatArray &solVec, FloatArray &x, FloatArray &m, double &gam, TimeStep *tStep);
    IntArray DelaminatedInterfaceList;
    void computeFailureCriteriaQuantities(FailureCriteriaStatus *fc, TimeStep *tStep) override;


    // Subdivision
    std :: vector< Triangle > allTri;

    std :: vector < std :: vector< Triangle > > crackSubdivisions;
    IntArray numSubDivisionsArray;

    void giveFictiousNodeCoordsForExport(std::vector<FloatArray> &nodes, int layer, int subCell);
    void giveFictiousCZNodeCoordsForExport(std::vector<FloatArray> &nodes, int layer, int subCell);
    void giveFictiousUpdatedNodeCoordsForExport(std::vector<FloatArray> &nodes, int layer, TimeStep *tStep, int subCell);
    void giveFictiousUpdatedCZNodeCoordsForExport(std::vector<FloatArray> &nodes, int layer, TimeStep *tStep, int subCell);
    void giveLocalNodeCoordsForExport(FloatArray &nodeLocalXi1Coords, FloatArray &nodeLocalXi2Coords, FloatArray &nodeLocalXi3Coords, int subCell, int layer, FloatMatrix &localNodeCoords);
    void giveLocalCZNodeCoordsForExport(FloatArray &nodeLocalXi1Coords, FloatArray &nodeLocalXi2Coords, FloatArray &nodeLocalXi3Coords, int subCell, FloatMatrix &localNodeCoords);
    void mapXi3FromLocalToShell(FloatArray &answer, FloatArray &local, int layer);
    void recoverValuesFromCZIP(std::vector<FloatArray> &recoveredValues, int interfce, InternalStateType type, TimeStep *tStep);

    FEI3dTrQuad interpolationForCZExport;
    FEI3dWedgeQuad interpolationForExport;

    std::vector< IntArray > orderingArrays;
    std::vector< IntArray > activeDofsArrays;

    void computeTripleProduct(FloatMatrix &answer, const FloatMatrix &a, const FloatMatrix &b, const FloatMatrix &c);

    // Recovery of through thickness stresses by momentum balance
    void recoverShearStress(TimeStep *tStep) override;

public:
    Shell7BaseXFEM(int n, Domain * d);
    virtual ~Shell7BaseXFEM();
    int checkConsistency() override;

    void giveMaxCZDamages(FloatArray &answer, TimeStep *tStep);
    const char *giveClassName() const override { return "Shell7BaseXFEM"; }
    std :: string errorInfo(const char *func) const { return std :: string(giveClassName()) + func; }
    Interface *giveInterface(InterfaceType it) override;

    IRResultType initializeFrom(InputRecord *ir) override;
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    int giveNumberOfDofs() override;

    // Recovery of through thickness stresses by momentum balance
    void giveFailedInterfaceNumber(IntArray &failedInterfaces, FloatArray &initiationFactor, TimeStep *tStep, bool recoverStresses = true);
    void giveAverageTransverseInterfaceStress(std::vector<FloatMatrix> &transverseStress, TimeStep *tStep);
    void giveRecoveredTransverseInterfaceStress(std::vector<FloatMatrix> &transverseStress, TimeStep *tStep) override;

    bool hasCohesiveZone(int interfaceNum);
    std :: vector< std :: unique_ptr< IntegrationRule > > czIntegrationRulesArray;
private:
    void jump(FloatMatrix lambda, FloatArray deltaUnknowns);
};
} // end namespace oofem
#endif
