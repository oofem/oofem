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

#ifndef huertaerrorestimator_h
#define huertaerrorestimator_h

#include "../sm/refinedelement.h"
#include "../sm/refinedmesh.h"
#include "errorestimator.h"
#include "interface.h"
#include "floatarray.h"
#include "statecountertype.h"
#include "remeshingcrit.h"

#include <vector>

///@name Input fields for HuertaErrorEstimator
//@{
#define _IFT_HuertaErrorEstimator_normtype "normtype"
#define _IFT_HuertaErrorEstimator_refinelevel "refinelevel"
#define _IFT_HuertaErrorEstimator_requirederror "requirederror"
#define _IFT_HuertaErrorEstimator_skipsteps "skipsteps"
#define _IFT_HuertaErrorEstimator_initialskipsteps "initialskipsteps"
#define _IFT_HuertaErrorEstimator_werror "werror"
#define _IFT_HuertaErrorEstimator_perfectCSect "perCSect"
#define _IFT_HuertaErrorEstimator_impCSect "impCSect"
#define _IFT_HuertaErrorEstimator_impPos "imppos"
#define _IFT_HuertaErrorEstimator_exact "exact"
//@}

///@name Input fields for HuertaRemeshingCriteria
//@{
#define _IFT_HuertaRemeshingCriteria_requirederror "requirederror"
#define _IFT_HuertaRemeshingCriteria_minelemsize "minelemsize"
#define _IFT_HuertaRemeshingCriteria_noremesh "noremesh"
#define _IFT_HuertaRemeshingCriteria_werror "werror"
#define _IFT_HuertaRemeshingCriteria_refinecoeff "refinecoeff"
//@}

///@name Input fields for HuertaErrorEstimatorInterface
//@{
#define _IFT_HuertaErrorEstimatorInterface_coords "coords"
//@}

namespace oofem {
class Element;
class GaussPoint;

/**
 * The implementation of Zienkiewicz Zhu Error Estimator.
 * The basic task is to evaluate the stress error on associated domain.
 * The algorithm is written in general way, so it is possible to to evaluate
 * different errors (for example temperature error). Then corresponding
 * member attribute identifying the type of quantity used should be declared and initialized
 * (for example in instanciateYourself() service). Then the modification is required
 * only when requesting element contributions.
 *
 * This task requires the special element algorithms, which are supported at element level
 * using interface concept.
 * This estimator also provides the compatible Remeshing Criteria, which
 * based on error measure will evaluate the required mesh density of a new domain.
 */
class HuertaErrorEstimator : public ErrorEstimator
{
public:
    /// Type of norm used.
    enum NormType { L2Norm, EnergyNorm };
    /// Mode of analysis.
    enum AnalysisMode { HEE_linear, HEE_nlinear };

protected:
    /// Global error norm.
    double globalENorm;
    /// Global weighted error norm.
    double globalWENorm;
    /// Global norm of primary unknown.
    double globalUNorm;
    /// Global error estimate (relative)
    double globalErrorEstimate;
    /// Cache storing element norms.
    FloatArray eNorms;
    /// Type of norm used.
    NormType normType;
    /// Actual state counter.
    StateCounterType stateCounter;
    /// Primary unknown nodal error.
    FloatArray primaryUnknownError;
    /// Refinement level.
    int refineLevel;
    /// Fine mesh.
    std :: vector< RefinedElement >refinedElementList;
    /// Mesh refinement.
    RefinedMesh refinedMesh;
    /// Linear analysis flag.
    AnalysisMode mode;
    /// Required error to obtain.
    double requiredError;
    /// Weighted error flag.
    bool wError;

    double lastError;
    int stepsToSkip, skippedSteps, maxSkipSteps, initialSkipSteps;

public:
    /// Constructor
    HuertaErrorEstimator(int n, Domain * d) : ErrorEstimator(n, d), eNorms(0), primaryUnknownError(0),
        refinedElementList(), refinedMesh()
    {
        eeType = EET_HEE;
        stateCounter = 0;
        normType = EnergyNorm;
        refineLevel = 1;
        mode = HEE_linear;
        wError = false;
        lastError = -1.0;
        stepsToSkip = skippedSteps = initialSkipSteps = 0;
    }

    /// Destructor
    virtual ~HuertaErrorEstimator() { }

    /**
     * Returns refinement level
     */
    int giveRefinementLevel() { return this->refineLevel; }

    virtual double giveElementError(EE_ErrorType type, Element *elem, TimeStep *tStep);

    virtual double giveValue(EE_ValueType type, TimeStep *tStep);

    virtual int estimateError(EE_ErrorMode err_mode, TimeStep *tStep);

    virtual RemeshingCriteria *giveRemeshingCrit();

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual const char *giveInputRecordName() const { return NULL; }
    virtual const char *giveClassName() const { return "HuertaErrorEstimator"; }

    AnalysisMode giveAnalysisMode() { return mode; }

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

private:
    /**
     * Builds refined mesh
     */
    void buildRefinedMesh();

    /**
     * Solves the refined element problem.
     * @param elemId Element id.
     * @param localNodeIdArray Array of local problem node ids.
     * @param globalNodeIdArray Array of global problem node ids.
     * @param tStep Time step.
     */
    void solveRefinedElementProblem(int elemId, IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                    TimeStep *tStep);
    /**
     * Solves the refined patch problem.
     * @param nodeId Node id.
     * @param localNodeIdArray Array of local problem node ids.
     * @param globalNodeIdArray Array of global problem node ids.
     * @param tStep Time step.
     */
    void solveRefinedPatchProblem(int nodeId, IntArray &localNodeIdArray,
                                  IntArray &globalNodeIdArray, TimeStep *tStep);
    /**
     * Solves the refined whole problem.
     * @param localNodeIdArray Array of local problem node ids.
     * @param globalNodeIdArray Array of global problem node ids.
     * @param tStep Time step.
     */
    void solveRefinedWholeProblem(IntArray &localNodeIdArray, IntArray &globalNodeIdArray, TimeStep *tStep);
    /**
     * Extracts nodal vector from global vector for each dof of all element nodes.
     * @param element Element.
     * @param vector Global vector.
     * @param answer Element nodal vector.
     * @param dofs Number of dofs at each node.
     * @param tStep Time step.
     */
    void extractVectorFrom(Element *element, FloatArray &vector, FloatArray &answer, int dofs, TimeStep *tStep);

    void setupRefinedProblemProlog(const char *problemName, int problemId, IntArray &localNodeIdArray,
                                   int nodes, int elems, int csects, int mats, int loads, int funcs,
                                   IntArray &controlNode, IntArray &controlDof, TimeStep *tStep);
    void setupRefinedProblemEpilog1(int csects, int mats, int loads, int nlbarriers);
    void setupRefinedProblemEpilog2(int tfuncs);
};


/**
 * The element interface corresponding to HuertaErrorEstimator.
 * It declares necessary services provided by element to be compatible with HuertaErrorEstimator.
 */
class HuertaErrorEstimatorInterface : public Interface
{
public:
    /// Mode for problem setup.
    enum SetupMode { CountMode = 0, NodeMode = 1, ElemMode = 2, BCMode = 3 };

public:
    /// Constructor
    HuertaErrorEstimatorInterface() { }

    virtual void HuertaErrorEstimatorI_setupRefinedElementProblem(RefinedElement *refinedElement, int level, int nodeId,
                                                                  IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                                  HuertaErrorEstimatorInterface :: SetupMode mode, TimeStep *tStep,
                                                                  int &localNodeId, int &localElemId, int &localBcId,
                                                                  IntArray &controlNode, IntArray &controlDof,
                                                                  HuertaErrorEstimator :: AnalysisMode aMode) = 0;

    virtual void HuertaErrorEstimatorI_computeNmatrixAt(GaussPoint *gp, FloatMatrix &answer) = 0;

protected:
    void setupRefinedElementProblem1D(Element *element, RefinedElement *refinedElement,
                                      int level, int nodeId, IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                      HuertaErrorEstimatorInterface :: SetupMode mode, TimeStep *tStep, int nodes,
                                      FloatArray **corner, FloatArray &midNode,
                                      int &localNodeId, int &localElemId, int &localBcId,
                                      IntArray &controlNode, IntArray &controlDof,
                                      HuertaErrorEstimator :: AnalysisMode aMode, const char *edgetype);

    void setupRefinedElementProblem2D(Element *element, RefinedElement *refinedElement,
                                      int level, int nodeId, IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                      HuertaErrorEstimatorInterface :: SetupMode mode, TimeStep *tStep, int nodes,
                                      FloatArray **corner, FloatArray *midSide, FloatArray &midNode,
                                      int &localNodeId, int &localElemId, int &localBcId,
                                      IntArray &controlNode, IntArray &controlDof,
                                      HuertaErrorEstimator :: AnalysisMode aMode, const char *quadtype);

    void setupRefinedElementProblem3D(Element *element, RefinedElement *refinedElement,
                                      int level, int nodeId, IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                      HuertaErrorEstimatorInterface :: SetupMode mode, TimeStep *tStep, int nodes,
                                      FloatArray **corner, FloatArray *midSide, FloatArray *midFace, FloatArray &midNode,
                                      int &localNodeId, int &localElemId, int &localBcId,
                                      int hexaSideNode [ 1 ] [ 3 ], int hexaFaceNode [ 1 ] [ 3 ],
                                      IntArray &controlNode, IntArray &controlDof,
                                      HuertaErrorEstimator :: AnalysisMode aMode, const char *hexatype);
};


/**
 * The class representing Huerta remeshing criteria.
 * The basic task is to evaluate the required mesh density (at nodes) on given domain,
 * based on information provided by the compatible error estimator.
 *
 * The remeshing criteria is maintained by the corresponding error estimator. This is mainly due to fact, that is
 * necessary for given EE to create compatible RC. In our concept, the EE is responsible.
 */
class HuertaRemeshingCriteria : public RemeshingCriteria
{
public:
    /// Mode of receiver, allows to use it in more general situations.
    enum HuertaRemeshingCriteriaModeType { primaryUnknownBased };

protected:
    /// Array of nodal mesh densities.
    FloatArray nodalDensities;
    /// Remeshing strategy proposed.
    RemeshingStrategy remeshingStrategy;
    /// Actual values (densities) state counter.
    StateCounterType stateCounter;
    /// Mode of receiver.
    HuertaRemeshingCriteriaModeType mode;
    /// Required error to obtain.
    double requiredError;
    /// Minimum element size alloved.
    double minElemSize;
    /// Refinement coefficient.
    double refineCoeff;
    /// Remeshing flag.
    bool noRemesh;
    /// Weighted error flag.
    bool wError;

public:
    /// Constructor
    HuertaRemeshingCriteria(int n, ErrorEstimator * e);
    /// Destructor
    virtual ~HuertaRemeshingCriteria() { }

    virtual double giveRequiredDofManDensity(int num, TimeStep *tStep, int relative = 0);
    virtual double giveDofManDensity(int num);
    virtual RemeshingStrategy giveRemeshingStrategy(TimeStep *tStep);
    virtual int estimateMeshDensities(TimeStep *tStep);
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual const char *giveInputRecordName() const { return NULL; }
    virtual const char *giveClassName() const { return "HuertaErrorEstimator"; }
};

} // end namespace oofem
#endif // huertaerrorestimator_h
