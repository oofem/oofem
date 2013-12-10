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

#include "shell7base.h"
#include "xfemelementinterface.h"


///@name Input fields for el
//@{
#define _IFT_Shell7BaseXFEM_CohesiveZoneMaterial "czmaterial"
//@}


namespace oofem {

class FEI3dTrQuad;
class BoundaryLoad;
class EnrichmentItem;

/**
 * This class represent a 7 parameter shell element. 
 * Each node has 7 degrees of freedom (displ. vec., director vec., inhomogeneous thickness strain ).
 * Add ref. to paper!
 * @author Jim Brouzoulis
 * @date 2012-11-01
 */

class Shell7BaseXFEM : public Shell7Base, public XfemElementInterface
{
protected:
    Material *czMat; // cohesive zone material
    //StructuralInterfaceMaterial *mat; // should replace with this
    int czMatNum;
    XfemManager *xMan;

    virtual void updateYourself(TimeStep *tStep);
    virtual void postInitialize();
    void computeOrderingArray( IntArray &orderingArray, IntArray &activeDofsArray,  EnrichmentItem *ei);

    virtual void evalCovarBaseVectorsAt(FloatArray &lCoords, FloatMatrix &gcon, FloatArray &solVec);
    void discGiveInitialSolutionVector(FloatArray &answer, IntArray &eiDofIdArray); // should be replaced with general function
    void computeDiscGeneralizedStrainVector(FloatArray &dGenEps, FloatArray &lCoords, EnrichmentItem *ei, TimeStep *tStep);
    void giveDisSolutionVector(FloatArray &answer, const IntArray &dofIdArray, TimeStep *tStep);

    // Internal forces
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    //void discComputeSectionalForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, FloatArray &solVecD, Delamination *dei);
    void discComputeSectionalForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, FloatArray &solVecD, EnrichmentItem *ei);
    double evaluateLevelSet(const FloatArray &lCoords, EnrichmentItem *ei);
    void computeCohesiveForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, FloatArray &solVecD, Delamination *dei);

    // Tangent matrices
    void computeLambdaGMatricesDis(FloatMatrix lambdaD [ 3 ], double zeta);
    void computeLambdaNMatrixDis(FloatMatrix &lambda_xd, double zeta);
    virtual void computeStiffnessMatrixOLD(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual void discComputeBulkTangentMatrix(FloatMatrix &KCC, FloatMatrix &KCD, FloatMatrix &KDD, IntegrationPoint *ip, Material *mat, int layer, TimeStep *tStep);

    virtual void computeShiftedBmatrixAt(FloatArray &lCoords, FloatMatrix &answer, int li = 1, int ui = ALL_STRAINS);
    virtual void computeShiftedNmatrixAt(const FloatArray &iLocCoords, FloatMatrix &answer);

    void computeCohesiveTangent(FloatMatrix &answer, TimeStep *tStep);
    void computeCohesiveTangentAt(FloatMatrix &answer, TimeStep *tStep, FloatArray &solVecD, Delamination *dei);

    void computePressureTangentMatrixDis(FloatMatrix &KCC, FloatMatrix &KCD, FloatMatrix &KDD, IntegrationPoint *ip, Load *load, const int iSurf, TimeStep *tStep);

    // External loads
    virtual void computeEdgeLoadVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *tStep, ValueModeType mode);
    virtual void computeSurfaceLoadVectorAt(FloatArray &answer, Load *load,
                                         int iSurf, TimeStep *tStep, ValueModeType mode);

    // Mass matrices
    void computeMassMatrixNum(FloatMatrix &answer, TimeStep *tStep); 

    // VTK
    virtual void giveCompositeExportData(VTKPiece &vtkPiece, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep );
    virtual void vtkEvalUpdatedGlobalCoordinateAt(FloatArray &localCoords, int layer, FloatArray &globalCoords, TimeStep *tStep);

    IntArray DelaminatedInterfaceList;
    void computeFailureCriteriaQuantities(FailureCriteriaStatus *fc, TimeStep *tStep);

public:
    Shell7BaseXFEM(int n, Domain *d);   
    virtual ~Shell7BaseXFEM() {};		
    virtual int checkConsistency();
    
    void giveMaxCZDamages(FloatArray &answer, TimeStep *tStep);
    virtual const char *giveClassName()  const { return "Shell7BaseXFEM"; }
    virtual Interface *giveInterface(InterfaceType it);
    
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;
    virtual int giveNumberOfDofs();

    bool hasCohesiveZone(int interfaceNum);
    IntegrationRule **czIntegrationRulesArray;
    IntegrationRule giveCZIntegrationRulesArray() { return **czIntegrationRulesArray; };
};



} // end namespace oofem
#endif 
