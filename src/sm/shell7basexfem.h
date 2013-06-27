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

#ifndef Shell7BaseXFEM_h
#define Shell7BaseXFEM_h

//#include "eleminterpmapperinterface.h"
//#include "nodalaveragingrecoverymodel.h"
//#include "layeredcrosssection.h"
//#include "nlstructuralelement.h"
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
    int czMatNum;
    virtual double giveGlobalZcoord(GaussPoint *gp);
    std::list< std::pair<int, double> > delaminationXiCoordList;
    void setupDelaminationXiCoordList();
    void setupDelaminationXiCoordsAtGP();

    std::list<std::pair<int, int> > gpDelaminationGroupList; // (gp#, dGroup#)
    void setupGPDelaminationGroupList();
    int giveDelaminationGroupAt(double zeta); 
    void giveDelaminationGroupXiLimits(int &dGroup, double &zTop, double &zBottom);
    double giveDelaminationGroupMidXi(int dGroup);

    XfemManager *xMan;
    
    static bool __initializedFieldDofId;
    static IntArray dofId_Midplane;
    static IntArray dofId_Director;
    static IntArray dofId_InhomStrain;
    static bool initDofId() {
        dofId_Midplane.setValues(3, D_u, D_v, D_w);
        dofId_Director.setValues(3,W_u, W_v, W_w);
        dofId_InhomStrain.setValues(1, Gamma);
        return true;
    }

    virtual IntArray giveFieldDofId(SolutionField fieldType) const;

    static bool sortFunc(std::pair<int, double> a, std::pair<int, double> b) {
        return a.second < b.second;
    }
    IntegrationRule **czIntegrationRulesArray;

    virtual void evalCovarBaseVectorsAt(GaussPoint *gp, FloatMatrix &gcon, FloatArray &solVec);
    void discGiveInitialSolutionVector(FloatArray &answer, IntArray &eiDofIdArray);
    void computeDiscGeneralizedStrainVector(FloatArray &dGenEps, GaussPoint *gp, EnrichmentItem *ei, int enrichmentDomainNumber, TimeStep *tStep);

    // compute solution vector
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    void discComputeSectionalForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, FloatArray &solVecD, int useUpdatedGpRecord,  
          EnrichmentItem *ei, int enrichmentDomainNumber);
    void computeCohesiveForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, FloatArray &solVecD, int useUpdatedGpRecord,  
          Delamination *dei, int enrichmentDomainNumber);
    void computeCohesiveTangent(FloatMatrix &answer, TimeStep *tStep);
    void computeCohesiveTangentAt(FloatMatrix &answer, TimeStep *tStep, FloatArray &solVec, FloatArray &solVecD, Delamination *dei, int enrichmentDomainNumber);

    void computeOrderingArray(IntArray &orderingArray, IntArray &activeDofsArray, int enrichmentDomainNumber, SolutionField field);
    //void edgeComputeOrderingArray( IntArray &orderingArray, IntArray &activeDofsArray, int iEdge, int enrichmentDomainNumber, SolutionField field);

    // Tangent matrices
    virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    virtual void discComputeBulkTangentMatrix(FloatMatrix &answer, FloatArray &solVec, FloatArray &solVecI, FloatArray &solVecJ, MatResponseMode rMode, TimeStep *tStep,
         EnrichmentItem *ei, int enrichmentDomainNumberI, int enrichmentDomainNumberJ);

    // Loads
    virtual void computeEdgeLoadVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *tStep, ValueModeType mode);
    virtual void computeSurfaceLoadVectorAt(FloatArray &answer, Load *load,
                                         int iSurf, TimeStep *tStep, ValueModeType mode);

    // Mass matrices
    void computeMassMatrixNum(FloatMatrix &answer, TimeStep *tStep); 

    // VTK
    virtual void vtkEvalUpdatedGlobalCoordinateAt(FloatArray &localCoords, int layer, FloatArray &globalCoords, TimeStep *tStep);

	virtual void updateYourself(TimeStep *tStep);


public:
    // constructor
    Shell7BaseXFEM(int n, Domain *d);   // : Shell7Base(n, d),  XfemElementInterface(this);
    virtual ~Shell7BaseXFEM() {};		// destructor -> declaring as virtual will make each subclass call their respective destr.
    virtual int checkConsistency();

    virtual const char *giveClassName()  const { return "Shell7BaseXFEM"; }
    virtual classType giveClassID()      const { return Shell7BaseXFEMClass; }
    virtual Interface *giveInterface(InterfaceType it);
    
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;
    void discGiveDofManDofIDMask(int inode,  int enrichmentdomainNumber, IntArray &answer) const;
    virtual int giveNumberOfDofs();
    bool hasCohesiveZone();
};



} // end namespace oofem
#endif 
