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

#ifndef Shell7Base_h
#define Shell7Base_h

#include "sm/CrossSections/layeredcrosssection.h"
#include "sm/Elements/nlstructuralelement.h"
#include "eleminterpmapperinterface.h"
#include "nodalaveragingrecoverymodel.h"
#include "vtkxmlexportmodule.h"
#include "zznodalrecoverymodel.h"
#include "fei3dwedgequad.h"
#include "fei3dtrquad.h"
#include "fracturemanager.h"
#include "cltypes.h"
#include <vector>
#include "sm/EngineeringModels/structengngmodel.h"

///@name Input fields for shell7base
//@{
#define _IFT_Shell7base_recoverStress "recoverstress"
//@}

namespace oofem {
class BoundaryLoad;

#define _ExportCZ

/**
 * This class represent a 7 parameter shell element.
 * Each node has 7 degrees of freedom (displ. vec., director vec., inhomogeneous thickness strain ).
 * @todo Add ref. to paper!
 * @author Jim Brouzoulis
 * @date 2012-11-01
 */
class Shell7Base : public NLStructuralElement, public NodalAveragingRecoveryModelInterface, public LayeredCrossSectionInterface,
    public VTKXMLExportModuleElementInterface, public ZZNodalRecoveryModelInterface, public FailureModuleElementInterface
{
public:
    Shell7Base(int n, Domain *d); // constructor

    void giveDofManDofIDMask(int inode, IntArray &) const override;
    int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords) override;
    virtual int computeGlobalCoordinatesOnEdge(FloatArray &answer, const FloatArray &lcoords, const int iEdge);
    int computeNumberOfDofs() override { return this->giveNumberOfDofs(); }
    int checkConsistency() override;
    void postInitialize() override;
    void printOutputAt(FILE *file, TimeStep *tStep) override;

    // Definition & identification
    const char *giveClassName() const override { return "Shell7Base"; }
    MaterialMode giveMaterialMode() override { return _3dMat; }


    // Element specific
    virtual int giveNumberOfInPlaneIP() { return numInPlaneIP; };
    int giveNumberOfDofs() override;
    virtual int giveNumberOfEdgeDofs() = 0;
    virtual int giveNumberOfEdgeDofManagers() = 0;
    FloatMatrixF< 3, 3 >evalInitialCovarBaseVectorsAt(const FloatArrayF< 3 > &lCoords);

    static void giveGeneralizedStrainComponents(FloatArray genEps, FloatArrayF< 3 > &dphidxi1, FloatArrayF< 3 > &dphidxi2, FloatArrayF< 3 > &dmdxi1,
                                                FloatArrayF< 3 > &dmdxi2, FloatArrayF< 3 > &m, double &dgamdxi1, double &dgamdxi2, double &gam);
    static FloatMatrixF< 3, 3 >giveDualBase(FloatMatrixF< 3, 3 > &base1);
    LayeredCrossSection *giveLayeredCS() { return this->layeredCS; }

    // Overloaded, as the element is using enhanced approximation
    void computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep, bool global) override;

    void computeConstitutiveMatrix_dPdF_At(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override
    { OOFEM_ERROR("calling of this funciton is not allowed"); }

protected:
    // Recover transverse stresses using momentum balance, cf. Främby, Fagerström & Bouzoulis, 'Adaptive modelling of delamination initiation and propagation using an equivalent single-layer shell approach', IJNME, 2016
    bool recoverStress;

    Interface *giveInterface(InterfaceType it) override;
    LayeredCrossSection *layeredCS;

    static FEI3dTrQuad interpolationForCZExport;
    static FEI3dWedgeQuad interpolationForExport;

    FEInterpolation3d *fei;

    virtual const IntArray &giveOrderingDofTypes() const = 0;
    virtual const IntArray &giveOrderingNodes() const = 0;
    virtual const IntArray &giveOrderingEdgeNodes() const = 0;

    std::vector< FloatArrayF< 3 > >initialNodeDirectors;
    FloatArrayF< 3 > &giveInitialNodeDirector(int i) {
        return this->initialNodeDirectors [ i - 1 ];
    }

    FloatArray initialSolutionVector;
    FloatArray &giveInitialSolutionVector() {
        return this->initialSolutionVector;
    }

    std::vector< FloatArray >initialEdgeSolutionVectors;
    FloatArray &giveInitialEdgeSolutionVector(int i) {
        return this->initialEdgeSolutionVectors [ i - 1 ];
    }

    // Element specific methods
    int numInPlaneIP;
    void computeGaussPoints() override = 0;
    virtual double computeVolumeAroundLayer(GaussPoint *mastergp, int layer) = 0;
    virtual double computeAreaAround(GaussPoint *gp, double xi) = 0;
    void giveSurfaceDofMapping(IntArray &answer, int iSurf) const override = 0;
    void giveEdgeDofMapping(IntArray &answer, int iEdge) const override = 0;

    void initializeFrom(InputRecord &ir) override;

    // Integration
    virtual double edgeComputeLengthAround(GaussPoint *gp, const int iedge);


    // Base vectors and directors
    virtual void setupInitialNodeDirectors();
    FloatArrayF< 3 >evalInitialDirectorAt(const FloatArrayF< 3 > &lCoords);



    FloatMatrixF< 3, 3 >evalInitialContravarBaseVectorsAt(const FloatArrayF< 3 > &lCoords);

    virtual FloatMatrixF< 3, 3 >evalCovarBaseVectorsAt(const FloatArrayF< 3 > &lCoords, FloatArray &genEps, TimeStep *tStep);

    virtual FloatArrayF< 3 >evalCovarNormalAt(const FloatArrayF< 3 > &lCoords, FloatArray &genEpsC, TimeStep *tStep);
    virtual FloatArrayF< 3 >evalInitialCovarNormalAt(const FloatArrayF< 3 > &lCoords);

    FloatArrayF< 3 >edgeEvalInitialDirectorAt(const FloatArrayF< 1 > &lCoords, const int iEdge);

    std::pair< FloatArrayF< 3 >, FloatArrayF< 3 > >edgeEvalInitialCovarBaseVectorsAt(const FloatArrayF< 1 > &lCoords, const int iedge);

    FloatMatrixF< 3, 3 >edgeEvalCovarBaseVectorsAt(const FloatArrayF< 3 > &lCoords, const int iedge, TimeStep *tStep);

    virtual double giveGlobalZcoord(const FloatArrayF< 3 > &lCoords);
    virtual double giveGlobalZcoordInLayer(double xi, int layer);

    FloatMatrixF< 3, 3 >giveAxialMatrix(const FloatArrayF< 3 > &vec);

    // Stress and strain
    FloatMatrixF< 3, 3 >computeFAt(const FloatArrayF< 3 > &lCoords, FloatArray &genEps, TimeStep *tStep);
    FloatMatrixF< 3, 3 >computeStressMatrix(FloatArray &genEps, GaussPoint *gp, Material *mat, TimeStep *tStep);

    virtual void computeCauchyStressVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    // Mass matrices
    void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep) override;
    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep) override;    // analytically integrated through the thickness
    virtual void computeMassMatrixNum(FloatMatrix &answer, TimeStep *tStep); // numerical integration in B_X
    virtual void giveMassFactorsAt(GaussPoint *gp, FloatArray &answer, double &gam);
    void computeConvectiveMassForce(FloatArray &answer, TimeStep *tStep);
    void computeThicknessMappingCoeff(GaussPoint *gp, FloatArray &answer); // for analytically integrated mass matrix


    // Tangent matrices
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    virtual void computeBulkTangentMatrix(FloatMatrix &answer, FloatArray &solVec, TimeStep *tStep);
    void computeLinearizedStiffness(GaussPoint *gp,  StructuralMaterial *mat, TimeStep *tStep, FloatMatrix A[ 3 ] [ 3 ]);
    void computePressureTangentMatrix(FloatMatrix &answer, Load *load, const int iSurf, TimeStep *tStep);
    std::array< FloatMatrixF< 3, 18 >, 3 >computeLambdaGMatrices(FloatArray &genEps, double zeta);
    FloatMatrixF< 3, 7 >computeLambdaNMatrix(FloatArray &genEps, double zeta);

    // Internal forces
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0) override;
    void computeSectionalForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, int useUpdatedGpRecord = 0);
    void computeSectionalForcesAt(FloatArray &sectionalForces, IntegrationPoint *ip, Material *mat, TimeStep *tStep, FloatArray &genEpsC, double zeta);

    // External forces
    void computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *tStep, ValueModeType mode) override;
    void computePressureForce(FloatArray &answer, FloatArray solVec, const int iSurf, BoundaryLoad *surfLoad, TimeStep *tStep, ValueModeType mode);
    void computePressureForceAt(GaussPoint *gp, FloatArray &answer, const int iSurf, FloatArray genEps, BoundaryLoad *surfLoad, TimeStep *tStep, ValueModeType mode);
    virtual void computeTractionForce(FloatArray &answer, const int iedge, BoundaryLoad *edgeLoad, TimeStep *tStep, ValueModeType mode, bool map2elementDOFs = false);


    // Solution vectors
    void computeVectorOfDofIDs(const IntArray &dofIdArray, ValueModeType u, TimeStep *tStep, FloatArray &answer);
    void temp_computeBoundaryVectorOf(IntArray &dofIdArray, int boundary, ValueModeType u, TimeStep *tStep, FloatArray &answer);

    virtual void edgeGiveUpdatedSolutionVector(FloatArray &answer, const int iedge, TimeStep *tStep);

    void setupInitialSolutionVector();
    void setupInitialEdgeSolutionVector();

    //void giveInitialSolutionVector(FloatArray &answer);
    void giveUpdatedSolutionVector(FloatArray &answer, TimeStep *tStep);
    void giveUnknownsAt(const FloatArrayF< 3 > &lcoords, const FloatArray &solVec, FloatArrayF< 3 > &x, FloatArrayF< 3 > &m, double &gam, TimeStep *tStep);

    // Nodal averaging interface:
    void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep) override;

    // ZZ recovery
    void ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatArray &answer, GaussPoint *gp, InternalStateType type);
    void NodalRecoveryMI_computeNValProduct(FloatMatrix &answer, int layer, InternalStateType type, TimeStep *tStep);
    void NodalRecoveryMI_computeNNMatrix(FloatArray &answer, int layer, InternalStateType type);
    void NodalRecoveryMI_recoverValues(std::vector< FloatArray > &recoveredValues, int layer, InternalStateType type, TimeStep *tStep);

    // VTK interface
    virtual FloatArrayF< 3 >vtkEvalInitialGlobalCoordinateAt(const FloatArrayF< 3 > &localCoords, int layer);
    virtual FloatArrayF< 3 >vtkEvalUpdatedGlobalCoordinateAt(const FloatArrayF< 3 > &localCoords, int layer, TimeStep *tStep);
    virtual FloatArrayF< 3 >vtkEvalInitialGlobalCZCoordinateAt(const FloatArrayF< 3 > &localCoords, int interface);

    void giveCompositeExportData(std::vector< ExportRegion > &vtkPieces, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep ) override;
    virtual void giveShellExportData(ExportRegion &vtkPiece, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep );

    std::vector< FloatArray >giveFictiousNodeCoordsForExport(int layer);
    std::vector< FloatArray >giveFictiousCZNodeCoordsForExport(int interface);
    std::vector< FloatArray >giveFictiousUpdatedNodeCoordsForExport(int layer, TimeStep *tStep);
    //void giveLocalNodeCoordsForExport(FloatArray &nodeLocalXi1Coords, FloatArray &nodeLocalXi2Coords, FloatArray &nodeLocalXi3Coords);

    // Recovery of through thickness stresses by momentum balance
    enum stressRecoveryType { copyIPvalue, LSfit, L2fit };
    virtual void giveRecoveredTransverseInterfaceStress(std::vector< FloatMatrix > &transverseStress, TimeStep *tStep);
    void giveTractionBC(FloatMatrix &tractionTop, FloatMatrix &tractionBtm, TimeStep *tStep);
    void recoverValuesFromIP(std::vector< FloatArray > &nodes, int layer, InternalStateType type, TimeStep *tStep, stressRecoveryType SRtype = copyIPvalue);
    void CopyIPvaluesToNodes(std::vector< FloatArray > &recoveredValues, int layer, InternalStateType type, TimeStep *tStep);
    void nodalLeastSquareFitFromIP(std::vector< FloatArray > &recoveredValues, int layer, InternalStateType type, TimeStep *tStep);
    virtual void recoverShearStress(TimeStep *tStep);
    void giveLayerContributionToSR(FloatMatrix &dSmat, FloatMatrix &dSmatLayerIP, int layer, double zeroThicknessLevel, TimeStep *tStep);
    void fitRecoveredStress2BC(std::vector< FloatMatrix > &answer1, std::vector< FloatMatrix > &answer2, std::vector< FloatMatrix > &dSmat, std::vector< FloatMatrix > &dSmatIP, FloatMatrix &SmatOld, FloatMatrix &tractionBtm, FloatMatrix &tractionTop, double zeroThicknessLevel, FloatArray fulfillBC, int startLayer, int endLayer);
    void updateLayerTransvStressesSR(FloatMatrix &dSmatLayerIP, int layer);
    void updateLayerTransvShearStressesSR(FloatMatrix &dSmatLayerIP, FloatMatrix &SmatOld, int layer);
    void updateLayerTransvNormalStressSR(FloatMatrix &dSzzMatLayerIP, FloatArray &SzzMatOld, int layer);
    //    void computeBmatrixForStressRecAt(FloatArray &lcoords, FloatMatrix &answer, int layer, bool intSzz = false);
    //     void givePolynomial2GradientForStressRecAt(FloatArray &answer, FloatArray &coords);
    void giveZintegratedPolynomialGradientForStressRecAt(FloatArray &answer, FloatArray &coords);
    //    void giveZintegratedPolynomial2GradientForStressRecAt(FloatArray &answer, FloatArray &coords);
    void giveZ2integratedPolynomial2GradientForStressRecAt(FloatArray &answer, FloatArray &coords);
    void giveL2contribution(FloatMatrix &ipValues, FloatMatrix &Nbar, int layer, InternalStateType type, TimeStep *tStep);
    void giveSPRcontribution(FloatMatrix &eltIPvalues, FloatMatrix &eltPolynomialValues, int layer, InternalStateType type, TimeStep *tStep);

    // N and B matrices
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li = 1, int ui = ALL_STRAINS) override { answer.clear(); }

    virtual void computeBmatrixAt(const FloatArray &lCoords, FloatMatrix &answer, int li = 1, int ui = ALL_STRAINS);
    void computeNmatrixAt(const FloatArray &iLocCoords, FloatMatrix &answer) override;
    virtual void edgeComputeNmatrixAt(const FloatArray &lCoords, FloatMatrix &answer);

    void computeStrainVectorInLayer(FloatArray &answer, const FloatArray &masterGpStrain, GaussPoint *masterGp, GaussPoint *slaveGp, TimeStep *tStep) override
    {
        OOFEM_ERROR("ComputeStrainVectorInLayer - Should not be called! Not meaningful for this element.");
    }
    virtual void edgeComputeBmatrixAt(const FloatArray &lCoords, FloatMatrix &answer, int li = 1, int ui = ALL_STRAINS);

    static FloatArrayF< 9 >convV6ToV9Stress(const FloatArrayF< 6 > &V6);

    int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep) override;
    void computeInterLaminarStressesAt(int interfaceNum, TimeStep *tStep, std::vector< FloatArray > &interLamStresses);
    virtual void evaluateFailureCriteriaQuantities(FailureCriteriaStatus *fc, TimeStep *tStep);
    int giveSymVoigtIndex(int ind1, int ind2);
    int giveVoigtIndex(int ind1, int ind2);
    std::vector< std::vector< int > >voigtIndices;
};
} // end namespace oofem
#endif
