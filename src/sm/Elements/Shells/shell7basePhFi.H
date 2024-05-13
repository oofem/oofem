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

#ifndef Shell7BasePhFi_h
#define Shell7BasePhFi_h

#include "eleminterpmapperinterface.h"
#include "nodalaveragingrecoverymodel.h"
#include "layeredcrosssection.h"
#include "nlstructuralelement.h"
#include "vtkxmlexportmodule.h"
#include "zznodalrecoverymodel.h"
#include "fei3dwedgequad.h"
#include "fei3dtrquad.h"
#include "fracturemanager.h"
#include "cltypes.h"
#include <vector>
#include "phasefieldelement.h"
#include "shell7base.h"


namespace oofem {
class BoundaryLoad;

#define _ExportCZ

/**
 * This class represent a 7 parameter shell element with isotropic damage in each layer modelled by phase field theory.
 * Each node has 7 degrees of freedom (displ. vec., director vec., inhomogeneous thickness strain ).
 * @todo Add ref. to paper!
 * @author Martin Fagerstr�m
 * @date 2014-06-03
 */
class Shell7BasePhFi : public Shell7Base, public PhaseFieldElement
{
public:
    int startIDdamage, endIDdamage;

    Shell7BasePhFi(int n, Domain *d); // constructor
    virtual ~Shell7BasePhFi() {}
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;
    virtual void giveDofManDofIDMask_u(IntArray &answer) const;
    virtual void giveDofManDofIDMask_d(IntArray &answer) const;
    //virtual void giveUpdatedSolutionVector(FloatArray &answer, TimeStep *tStep);


    double computeGInLayer(int layer, GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN);
    double computeGInLayer_dist(int layer, int indx, GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN);
    double computeGprimInLayer(int layer, GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN);
    double computeGprimInLayer_dist(int layer, int indx, GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN);
    double computeGbisInLayer(int layer, GaussPoint *gp, ValueModeType valueMode, TimeStep *stepN);
    double computeDamageInLayerAt(int layer, GaussPoint *gp, ValueModeType valueMode,  TimeStep *stepN);
    double computeOldDamageInLayerAt(int layer, GaussPoint *gp, ValueModeType valueMode,  TimeStep *stepN);
    double computeDamageInLayerAt_dist(int layer, int indx, GaussPoint *gp, ValueModeType valueMode,  TimeStep *stepN);
    IntArray computeDamageIndexArray(int layer);

    int computeNumberOfDofs() override { return this->giveNumberOfDofs(); }
    //int checkConsistency() override;
    void postInitialize() override;

    // Definition & identification
    const char *giveClassName() const override { return "Shell7BasePhFi"; }
   // virtual classType giveClassID() const { return Shell7BasePhFiClass; }


    // Element specific
    int giveNumberOfDofs() override;
    virtual int giveNumberOfuDofs();
    virtual int giveNumberOfdDofs();
    //virtual int giveNumberOfEdgeDofs() = 0; //todo: used for PhFi???
    //virtual int giveNumberOfEdgeDofManagers() = 0;

    double neg_MaCauley(double par);
    double neg_MaCauleyPrime(double par);

protected:

    enum SolutionFieldPhFi {
        Midplane,       ///< phi_bar 7 x_bar (3 dofs)
        Director,       ///< m (3 dofs)
        InhomStrain,    ///< gamma (1 dofs) - inhomogenious thickness strain
        Displacement,   ///< phi_bar, m, gamma
        Damage,         ///< d (n dofs where n is the number of layers in the layeredCS
        All,
        AllInv,
        EdgeInv,
    };

    int numberOfLayers;

    virtual const IntArray &giveOrdering(SolutionField fieldType) const 
    { 
        OOFEM_ERROR("Shell7BasePhFi :: giveOrdering not implemented: Use Shell7BasePhFi :: giveOrderingPhFi instead");
        return 0; 
    };

    virtual const IntArray &giveOrderingPhFi(SolutionFieldPhFi fieldType) const = 0;
    virtual const IntArray &giveOrdering_Disp() const = 0;
    virtual const IntArray &giveOrdering_Damage() const = 0;
    void initializeFrom(InputRecord &ir) override;

    // Tangent matrices
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    virtual void new_computeBulkTangentMatrix(FloatMatrix &answer, FloatArray &solVec, FloatArray &solVecI, FloatArray &solVecJ, MatResponseMode rMode, TimeStep *tStep);

    void computeStiffnessMatrix_uu(FloatMatrix &answer, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_ud(FloatMatrix &answer, MatResponseMode, TimeStep *);
    void computeStiffnessMatrixNum_ud(FloatMatrix &answer, MatResponseMode, TimeStep *);
    void computeStiffnessMatrixNum_dd(FloatMatrix &answer, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_dd(FloatMatrix &answer, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_du(FloatMatrix &answer, MatResponseMode, TimeStep *);

    // Internal forces
    virtual void computeSectionalForces(FloatArray &answer, TimeStep *tStep, FloatArray &solVec, int useUpdatedGpRecord = 0);
    virtual void computeSectionalForces_dist(FloatArray &answer, int indx, TimeStep *tStep, FloatArray &solVec, int useUpdatedGpRecord = 0);
    virtual void computeSectionalForcesAt_d(double &ans_scalar, FloatArray &answer, IntegrationPoint *ip, Material *mat, TimeStep *tStep,  double zeta, int layer, FloatArray &Ddam_Dxi);
    virtual void computeSectionalForcesAt_d_dist(double &ans_scalar, FloatArray &answer, int indx, IntegrationPoint *ip, Material *mat, TimeStep *tStep,  double zeta, int layer, FloatArray &Ddam_Dxi);


    // Solution vectors
    //void giveSolutionVector_u(FloatArray &answer, const IntArray &dofIdArray, TimeStep *tStep);
    void giveSolutionVector_d(FloatArray &answer, const IntArray &dofIdArray, TimeStep *tStep);
    void computeVectorOfDofIDs(const IntArray &dofIdArray, ValueModeType u, TimeStep *stepN, FloatArray &answer);
    void giveUpdatedSolutionVector_d(FloatArray &answer, TimeStep *tStep);

    // VTK interface
    virtual void giveShellExportData(ExportRegion &vtkPiece, IntArray &primaryVarsToExport, IntArray &internalVarsToExport, IntArray cellVarsToExport, TimeStep *tStep );


    // N and B matrices
    virtual void computeBdmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li = 1, int ui = ALL_STRAINS)
    { OOFEM_ERROR1("Shell7BasePhFi :: computeBdmatrixAt should not be called");    };

    virtual void computeBdmatrixAt(FloatArray &lCoords, FloatMatrix &answer);

    virtual void computeNdvectorAt(const FloatArray &iLocCoords, FloatArray &answer);
    virtual void computeNdMatrixAt(const FloatArray &iLocCoords, FloatMatrix &answer); // same as vector but saved as a FloatMatrix
};
} // end namespace oofem
#endif
