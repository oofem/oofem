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
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef Shell7Base_h
#define Shell7Base_h


#include "eleminterpmapperinterface.h"
#include "nodalaveragingrecoverymodel.h"
#include "layeredcrosssection.h"

#include "nlstructuralelement.h"
#include <vector>
namespace oofem {
class BoundaryLoad;

/**
 * This class represent a 7 parameter shell element.
 * Each node has 7 degrees of freedom (displ. vec., director vec., inhomogeneous thickness strain ).
 * Add ref. to paper!
 * @author Jim Brouzoulis
 * @date 2012-11-01
 */
class Shell7Base : public NLStructuralElement, public NodalAveragingRecoveryModelInterface, public LayeredCrossSectionInterface
{
public:
    Shell7Base(int n, Domain *d); // constructor
    virtual ~Shell7Base() {}    // destructor -> declaring as virtual will make each subclass call their respective destr.
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;
    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    virtual int computeNumberOfDofs(EquationID ut) { return this->giveNumberOfDofs(); }


    // Definition & identification
    virtual const char *giveClassName() const {     return "Shell7Base"; }
    virtual classType giveClassID() const { return Shell7BaseClass; }
    virtual MaterialMode giveMaterialMode() { return _3dMat; }


    // Element specific
    virtual int giveNumberOfDofs() = 0;
    virtual int giveNumberOfEdgeDofs() = 0;
    virtual int giveNumberOfEdgeDofManagers() = 0;
    virtual Element_Geometry_Type giveGeometryType() const = 0;
    virtual FEInterpolation *giveInterpolation() = 0;
    virtual integrationDomain  giveIntegrationDomain() const = 0;


protected:
    virtual Interface *giveInterface(InterfaceType it);
    int numberOfGaussPoints;
    IntegrationRule **layerIntegrationRulesArray;

    enum SolutionField {
        Midplane,       ///< phi_bar 7 x_bar (3 dofs)
        Director,       ///< m (3 dofs)
        InhomStrain,    ///< gamma (1 dofs) - inhomogenious thickness strain
        All,
        AllInv,
        EdgeInv,
    };

    virtual IntArray giveOrdering(SolutionField fieldType) const = 0;

    std :: vector< FloatArray >initialNodeDirectors;

    FloatArray &giveInitialNodeDirector(int i) {
        return this->initialNodeDirectors [ i - 1 ];
    };

    // Element specific methods
    virtual void computeGaussPoints() = 0;
    virtual void giveLocalNodeCoords(FloatArray &nodeLocalXiCoords, FloatArray &nodeLocalEtaCoords) = 0;
    virtual double computeVolumeAroundLayer(GaussPoint *mastergp, int layer) = 0;
    virtual double computeAreaAround(GaussPoint *gp) = 0;
    virtual void giveSurfaceDofMapping(IntArray &answer, int iSurf) const = 0;
    virtual void giveEdgeDofMapping(IntArray &answer, int iEdge) const = 0;


    virtual IRResultType initializeFrom(InputRecord *ir);



    // Integration
    virtual double edgeComputeLengthAround(GaussPoint *gp, const int iedge);


    // Base vectors and directors
    virtual void setupInitialNodeDirectors();
    void evalInitialDirectorAt(GaussPoint *gp, FloatArray &answer);
    void evalInitialCovarBaseVectorsAt(GaussPoint *gp, FloatArray &G1, FloatArray &G2, FloatArray &G3);
    void evalInitialContravarBaseVectorsAt(GaussPoint *gp, FloatArray &G1, FloatArray &G2, FloatArray &G3);
    void giveDualBase(const FloatArray &G1, const FloatArray &G2, const FloatArray &G3, FloatArray &g1, FloatArray &g2, FloatArray &g3);
    virtual void evalCovarBaseVectorsAt(GaussPoint *gp, FloatArray &g1, FloatArray &g2, FloatArray &g3, FloatArray &solVec);
    void evalContravarBaseVectorsAt(GaussPoint *gp, FloatArray &g1, FloatArray &g2, FloatArray &g3, FloatArray &solVec);
    void edgeEvalInitialDirectorAt(GaussPoint *gp, FloatArray &answer, const int iEdge);
    void edgeEvalInitialCovarBaseVectorsAt(GaussPoint *gp, const int iedge, FloatArray &G1, FloatArray &G3);
    void edgeEvalCovarBaseVectorsAt(GaussPoint *gp, const int iedge, FloatArray &g1, FloatArray &g3, TimeStep *tStep);
virtual double giveGlobalZcoord(GaussPoint *gp);
//double giveLayerZetaCoord(GaussPoint *gp, int layer);


    // Stress and strain
    void computeFAt(GaussPoint *gp, FloatMatrix &answer, FloatArray &genEps);
    void computeCovarStressAt(GaussPoint *gp, FloatArray &answer);
    void giveGeneralizedStrainComponents(FloatArray genEps, FloatArray &dphidxi1, FloatArray &dphidxi2, FloatArray &dmdxi1,
                                         FloatArray &dmdxi2, FloatArray &m, double &dgamdxi1, double &dgamdxi2, double &gam);
    void computeStressResultantsAt(GaussPoint *gp, FloatArray &Svec, FloatArray &S1g, FloatArray &S2g, FloatArray &S3g, FloatArray &solVec);
    void computeStressVector(FloatArray &answer, FloatArray &genEps, GaussPoint *gp, Material *mat, TimeStep *stepN);
    virtual void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, FloatArray &genEps);


    // Mass matrices
    virtual void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep);    // analytically integrated through the thickness
    virtual void computeMassMatrixNum(FloatMatrix &answer, TimeStep *tStep); // numerical integration in B_X
    virtual void giveMassFactorsAt(GaussPoint *gp, FloatArray &answer, double &gam);
    void computeConvectiveMassForce(FloatArray &answer, TimeStep *tStep);
    void computeThicknessMappingCoeff(GaussPoint *gp, FloatArray &answer); // for analytically integrated mass matrix


    // Tangent matrices
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual void computeBulkTangentMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    void computeLinearizedStiffness(GaussPoint * gp,  Material * mat, TimeStep * tStep,
                                    FloatArray & S1g, FloatArray & S2g, FloatArray & S3g, FloatMatrix A [ 3 ] [ 3 ], FloatArray & solVec);
    void computePressureTangentMatrix(FloatMatrix &answer, Load *load, const int iSurf, TimeStep *tStep);


    // Internal forces
    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);
    void computeSectionalForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, int useUpdatedGpRecord = 0);
    void computeSectionalForcesOld(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, int useUpdatedGpRecord = 0);
    void computeSectionalForcesAt(FloatArray &N, FloatArray  &M, FloatArray &T, FloatArray  &Ms, double &Ts, GaussPoint *gp, Material *mat, TimeStep *tStep, FloatArray &genEps, double zeta);

    // External forces
    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *stepN, ValueModeType mode);
    void computeEdgeLoadVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *tStep, ValueModeType mode);
    void computeSurfaceLoadVectorAt(FloatArray &answer, Load *load, int iSurf, TimeStep *tStep, ValueModeType mode);
    void computePressureForce(FloatArray &answer, FloatArray solVec, const int iSurf, BoundaryLoad *surfLoad, TimeStep *tStep);
    void computePressureForceAt(GaussPoint *gp, FloatArray &answer, const int iSurf, FloatArray genEps, BoundaryLoad *surfLoad, TimeStep *tStep);
    virtual void computeTractionForce(FloatArray &answer, const int iedge, BoundaryLoad *edgeLoad, TimeStep *tStep);


    // Transformation
    void giveCoordTransMatrix(FloatMatrix &answer, FloatArray &g1, FloatArray &g2, FloatArray &g3,
                              FloatArray &G1, FloatArray &G2, FloatArray &G3);

    void giveCoordTransMatrix(FloatMatrix &answer, FloatArray &g1, FloatArray &g2, FloatArray &g3);

    void giveBondTransMatrix(FloatMatrix &answer, FloatMatrix &Q);
    void transInitialCartesianToInitialContravar(GaussPoint *gp, const FloatArray &VoightMatrix, FloatArray &answer);
    void transInitialCartesianToInitialContravar(GaussPoint *gp, const FloatMatrix &Stiffness, FloatMatrix &answer);


    // Solution vectors
    void computeGeneralizedStrainVector(FloatArray &answer, const FloatArray &solVec, const FloatMatrix &B11,
                                        const FloatMatrix &B22, const FloatMatrix &B32, const FloatMatrix &B43, const FloatMatrix  &B53);
    void computeSolutionFields(FloatArray &xbar, FloatArray &m, double &gam, const FloatArray &solVec, const FloatMatrix &N11, const FloatMatrix &N22, const FloatMatrix &N33);

    void edgeGiveUpdatedSolutionVector(FloatArray &answer, const int iedge, TimeStep *tStep);
    void edgeGiveInitialSolutionVector(FloatArray &answer, const int iedge);

    void giveInitialSolutionVector(FloatArray &answer);
    void giveUpdatedSolutionVector(FloatArray &answer, TimeStep *tStep);



    // Nodal averaging interface:
    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep);
    virtual void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side, InternalStateType type, TimeStep *tStep);
    virtual int  NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type);





    // N and B matrices

    int giveFieldSize(SolutionField fieldType);
    int giveNumberOfFieldDofs(SolutionField fieldType);
    virtual void computeFieldBmatrix(FloatMatrix & answer, FloatMatrix & dNdxi, SolutionField);
    virtual void computeFieldNmatrix(FloatMatrix & answer, FloatArray & N, SolutionField);
    void computeBmatricesAt(GaussPoint *gp, FloatMatrix &B11, FloatMatrix &B22, FloatMatrix &B32, FloatMatrix &B43, FloatMatrix &B53);
    void computeNmatricesAt(GaussPoint *gp, FloatMatrix &N11, FloatMatrix &N22, FloatMatrix &N33);
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li = 1, int ui = ALL_STRAINS);
    virtual void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void edgeComputeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void edgeComputeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li = 1, int ui = ALL_STRAINS);


    // Misc
    void computeTripleProduct(FloatMatrix &answer, const FloatMatrix &a, const FloatMatrix &b, const FloatMatrix &c);
    int giveVoigtIndex(const int ind1, const int ind2);
    void giveTensorForm(const FloatMatrix &matrix, FloatArray &tensor);
    void compareMatrices(const FloatMatrix &matrix1, const FloatMatrix &matrix2, FloatMatrix &answer);
};
} // end namespace oofem
#endif
