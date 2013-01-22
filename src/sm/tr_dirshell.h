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

#ifndef trdirshell_h
#define trdirshell_h


#include "eleminterpmapperinterface.h"
#include "nodalaveragingrecoverymodel.h"
#include "layeredcrosssection.h"

#include "nlstructuralelement.h"
namespace oofem {


class FEI2dTrQuad;
class FEI3dTrQuad;
class BoundaryLoad;

/**
 * This class represent a 7 parameter shell element. 
 * Each node has 7 degrees of freedom (displ. vec., director vec., inhomogeneous thickness strain ).
 *
 * @author Jim Brouzoulis
 * @date 2012-11-01
 */
class TrDirShell : public NLStructuralElement, public NodalAveragingRecoveryModelInterface, public LayeredCrossSectionInterface
{
protected:
    int numberOfGaussPoints;
	IntegrationRule **layerIntegrationRulesArray;
	FloatMatrix testmatrix;
    static FEI2dTrQuad interpolation_quad;
    static FEI3dTrQuad interpolation;
    static bool __initialized;
	static IntArray ordering_x;
	static IntArray ordering_m;
	static IntArray ordering_gam;
    static IntArray ordering_all;

    static bool initOrdering() {
        ordering_x.setValues(18, 1, 2, 3, 8, 9, 10, 15, 16, 17, 22, 23, 24, 29, 30 ,31, 36, 37, 38);
		ordering_m.setValues(18, 4, 5, 6, 11, 12, 13, 18, 19, 20, 25, 26, 27, 32, 33 ,34, 39, 40, 41);
		ordering_gam.setValues(6, 7, 14, 21, 28, 35, 42);
        ordering_all.setValues(42, 1, 2, 3, 8, 9, 10, 15, 16, 17, 22, 23, 24, 29, 30 ,31, 36, 37, 38,
                                   4, 5, 6, 11, 12, 13, 18, 19, 20, 25, 26, 27, 32, 33 ,34, 39, 40, 41,
                                   7, 14, 21, 28, 35, 42);
        return true;
    }

	FloatArray initialNodeDirectors[6];
	void setupInitialNodeDirectors();
	//FloatArray &giveInitialNodeDirector(int i){return this->initialNodeDirectors[i-1];};
    FloatArray giveInitialNodeDirector(int i);


    virtual void computeGaussPoints();

	// Interpolation functions
	virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li = 1, int ui = ALL_STRAINS);
	virtual void computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);

    // -Edge
	void edgeEvalN(FloatArray &answer, const FloatArray &lcoords);
    void edgeEvaldNdxi(FloatMatrix &answer, const FloatArray &lcoords);
    void edgeComputeNmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    void edgeComputeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li = 1, int ui = ALL_STRAINS);
    void edgeEvalCovarBaseVectorsAt(GaussPoint *gp, const int iedge, FloatArray &g1, FloatArray &g3, TimeStep *tStep);
    void edgeGiveUpdatedSolutionVector(FloatArray &answer,const int iedge, TimeStep *tStep);
    virtual double edgeComputeAreaAround(GaussPoint *gp, const int iedge);
    void edgeEvalInitialDirectorAt(GaussPoint *gp, FloatArray &answer, const int iEdge);
    void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    void edgeEvalInitialCovarBaseVectorsAt(GaussPoint *gp, const int iedge, FloatArray &G1, FloatArray &G3);

    virtual double computeVolumeAround(GaussPoint *gp);
    virtual double computeVolumeAroundLayer(GaussPoint *mastergp, int layer);
    virtual double computeThicknessAroundLayer(GaussPoint *gp, int layer);
    virtual double computeAreaAround(GaussPoint *gp);

    void giveSurfaceDofMapping(IntArray &answer, int iSurf) const;

   	// Base vectors and directors
	void evalInitialDirectorAt(GaussPoint *gp, FloatArray &answer);
	
	void evalInitialCovarBaseVectorsAt(GaussPoint *gp, FloatArray &G1, FloatArray &G2, FloatArray &G3);
	void evalInitialContravarBaseVectorsAt(GaussPoint *gp, FloatArray &G1, FloatArray &G2, FloatArray &G3);


	void giveDualBase(const FloatArray &G1, const FloatArray &G2, const FloatArray &G3, FloatArray &g1, FloatArray &g2, FloatArray &g3 );
	void evalCovarBaseVectorsAt(GaussPoint *gp, FloatArray &g1, FloatArray &g2, FloatArray &g3, TimeStep *tStep, FloatArray &solVec);
	void evalContravarBaseVectorsAt(GaussPoint *gp, FloatArray &g1, FloatArray &g2, FloatArray &g3, TimeStep *tStep, FloatArray &solVec);
    double giveLocalZetaCoord(GaussPoint *gp);
    double giveLayerZetaCoord(GaussPoint *gp, int layer);

	
    void giveUpdatedSolutionVector(FloatArray &answer, TimeStep *tStep);
    void computeThicknessMappingCoeff(GaussPoint *gp, FloatArray &answer); // for analytically integrated mass matrix

    // Loads
    void computeEdgeLoadVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *tStep, ValueModeType mode);
    void computeSurfaceLoadVectorAt(FloatArray &answer, Load *load, int iSurf, TimeStep *tStep, ValueModeType mode);
    void computeConvectiveMassForce(FloatArray &answer, TimeStep *tStep);
    void computeSectionalForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, int useUpdatedGpRecord = 0);
    void computePressureForce(FloatArray &answer, FloatArray solVec, const int iSurf, BoundaryLoad *surfLoad, TimeStep *tStep);
    void computePressureForceAt(GaussPoint *gp, FloatArray &answer, const int iSurf, FloatArray genEps, BoundaryLoad *surfLoad, TimeStep *tStep);

    void computeTractionForce(FloatArray &answer, const int iedge, BoundaryLoad *edgeLoad, TimeStep *tStep);
    void computePressureTangentMatrix(FloatMatrix &answer, Load *load, const int iSurf, TimeStep *tStep);

	// Stress and strain
	void computeFAt(GaussPoint *gp, FloatMatrix &answer, TimeStep *stepN, FloatArray &genEps);
	void computeCovarStressAt(GaussPoint *gp, FloatArray &answer);
	void transInitialCartesianToInitialContravar(GaussPoint *gp, const FloatArray &VoightMatrix, FloatArray &answer);
    void transInitialCartesianToInitialContravar(GaussPoint *gp, const FloatMatrix &Stiffness, FloatMatrix &answer);
    void giveGeneralizedStrainComponents(FloatArray genEps, FloatArray &dphidxi1, FloatArray &dphidxi2, FloatArray &dmdxi1, 
         FloatArray &dmdxi2, FloatArray &m, double &dgamdxi1, double &dgamdxi2, double &gam);
    void computeStressResultantsAt(GaussPoint *gp, FloatArray &Svec, FloatArray &S1g, FloatArray &S2g, FloatArray &S3g, TimeStep *tStep, FloatArray &solVec);
    void computeStressVector(FloatArray &answer, FloatArray &genEps, GaussPoint *gp, Material *mat, TimeStep *stepN );



    void computeSectionalForcesAt(FloatArray &answer, GaussPoint *gp, Material *mat, TimeStep *tStep, FloatArray &genEps, double zeta);


public:
    TrDirShell(int n, Domain *d);	// constructor
    virtual ~TrDirShell() { }		// destructor -> declaring as virtual will make each subclass call their respective destr.
	virtual int computeNumberOfDofs(EquationID ut) { return 42; } // 6*7dofs
	virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;

	
	virtual IRResultType initializeFrom(InputRecord *ir);

	virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
	

    // definition & identification
    virtual const char *giveClassName() const { return "TrDirShell"; }
    virtual classType giveClassID() const { return TrDirShellClass; }
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_triangle_2; }
	virtual MaterialMode giveMaterialMode() { return _3dMat; }
    virtual FEInterpolation *giveInterpolation();

	virtual integrationDomain  giveIntegrationDomain() { return _Triangle; } // write new wedge-like type

	virtual void printOutputAt(FILE *file, TimeStep *tStep);



	virtual void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, FloatArray &genEps);

	
	virtual void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
	virtual void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep);    // analytically integrated through the thickness
    virtual void computeMassMatrixNum(FloatMatrix &answer, TimeStep *tStep); // numerical integration in B_X


    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual void computeBulkTangentMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);

    virtual void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord = 0);
    

    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *forLoad, TimeStep *stepN, ValueModeType mode);

    void giveCoordTransMatrix(FloatMatrix &answer, FloatArray &g1, FloatArray &g2, FloatArray &g3,
                                                 FloatArray &G1, FloatArray &G2, FloatArray &G3);

    void giveCoordTransMatrix(FloatMatrix &answer, FloatArray &g1, FloatArray &g2, FloatArray &g3);

    void giveBondTransMatrix(FloatMatrix &answer, FloatMatrix &Q);

    void computeLinearizedStiffness(GaussPoint *gp,  Material *mat, TimeStep *tStep,
                FloatArray &S1g, FloatArray &S2g, FloatArray &S3g, FloatMatrix A[3][3], FloatArray &solVec);
    

     int giveVoigtIndex(const int ind1, const int ind2);

     void giveTensorForm(const FloatMatrix &matrix, FloatArray &tensor);

     void compareMatrices(const FloatMatrix &matrix1, const FloatMatrix &matrix2, FloatMatrix &answer);

     virtual Interface *giveInterface(InterfaceType it);
    // Nodal averaging interface:
    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep);
    virtual void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side, InternalStateType type, TimeStep *tStep);
    virtual int NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type);


    // layered cross section
        virtual void computeStrainVectorInLayer(FloatArray &answer, GaussPoint *masterGp,
            GaussPoint *slaveGp, TimeStep *tStep) {};

};



} // end namespace oofem
#endif 
