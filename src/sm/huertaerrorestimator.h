/* $Header: /home/cvs/bp/oofem/sm/src/huertaerrorestimator.h,v 1.2.4.1 2004/04/05 15:19:46 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   *************************************************************************
//   *** CLASS ERROR ESTIMATOR (INDICATOR) ACORDING TO ZIENKIEWICZ AND ZHU ***
//   *************************************************************************

#ifndef huertaerrorestimator_h
#define huertaerrorestimator_h

#include "errorestimator.h"
#include "interface.h"
#include "refinedelement.h"
#include "refinedmesh.h"
#include "alist.h"

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
    // type of norm used
    enum NormType { L2Norm, EnergyNorm };
    // mode of analysis
    enum AnalysisMode { HEE_linear, HEE_nlinear };

protected:
    /// global error norm
    double globalENorm;
    /// global weighted error norm
    double globalWENorm;
    /// global norm of primary unknown
    double globalUNorm;
    /// cache storing element norms
    FloatArray eNorms;
    /// type of norm used
    NormType normType;
    /// actual state counter.
    StateCounterType stateCounter;
    /// primary unknown nodal error
    FloatArray primaryUnknownError;
    /// refinement level
    int refineLevel;
    /// fine mesh
    AList< RefinedElement >refinedElementList;
    /// mesh refinement
    RefinedMesh refinedMesh;
    /// linear analysis flag
    AnalysisMode mode;
    /// required error to obtain
    double requiredError;
    /// weighted error flag
    bool wError;

    double lastError;
    int stepsToSkip, skippedSteps, maxSkipSteps, initialSkipSteps;

public:
    /// Constructor
    HuertaErrorEstimator(int n, Domain *d) : ErrorEstimator(n, d), eNorms(0), primaryUnknownError(0),
        refinedElementList(0), refinedMesh()
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
    ~HuertaErrorEstimator() { }
    /** Returns refinement level
     */
    int giveRefinementLevel(void) { return this->refineLevel; }
    /** Returns the element error of requested type. The estimateError service should be called before.
     * @param type error type
     * @param elem element for which error requested
     * @param tStep time step
     */
    virtual double giveElementError(EE_ErrorType type, Element *elem, TimeStep *tStep);
    /** Returns the characteristic value of given type.
     * The estimateError service should be called before. Intended to be used by remeshingCriterias to query
     * various values provided by specific error estimator.
     * @param type value type
     * @param tStep time step
     */
    virtual double giveValue(EE_ValueType type, TimeStep *tStep);
    /**
     * Estimates the error on associated domain at given timeStep.
     * @param tStep time step
     */
    virtual int estimateError(EE_ErrorMode mode, TimeStep *tStep);
    /** Returns reference to associated remeshing criteria.
     */
    virtual RemeshingCriteria *giveRemeshingCrit();
    /** Initializes receiver according to object description stored in input record.
     * This function is called immediately after creating object using
     * constructor. InitString can be imagined as data record in component database
     * belonging to receiver. Receiver may use value-name extracting functions
     * to extract particular field from record.
     * @see readInteger, readDouble and similar functions */
    virtual IRResultType initializeFrom(InputRecord *ir);
    /// Returns class name of the receiver.
    const char *giveClassName() const { return "HuertaErrorEstimator"; }
    /** Returns classType id of receiver.
     * @see FEMComponent::giveClassID
     */
    classType                giveClassID() const { return HuertaErrorEstimatorClass; }

    AnalysisMode giveAnalysisMode() { return mode; }

    /**
     * Stores context of receiver into given stream.
     * Only non-temp internal history variables are stored.
     * @param stream stream where to write data
     * @param obj pointer to integration point, which invokes this method
     * @return contextIOResultType.
     */
    contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    /**
     * Restores context of receiver from given stream.
     * @param stream stream where to read data
     * @param obj pointer to integration point, which invokes this method
     * @return contextIOResultType.
     */
    contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);

private:
    /** Builds refined mesh
     */
    void buildRefinedMesh(void);

    /** Solves the refined element problem.
     * @param elemId element id
     * @param localNodeIdArray array of local problem node ids
     * @param globalNodeIdArray array of global problem node ids
     * @param tStep time step
     */
    void solveRefinedElementProblem(int elemId, IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                    TimeStep *tStep);
    /** Solves the refined patch problem.
     * @param nodeId node id
     * @param localNodeIdArray array of local problem node ids
     * @param globalNodeIdArray array of global problem node ids
     * @param tStep time step
     */
    void solveRefinedPatchProblem(int nodeId, IntArray &localNodeIdArray,
                                  IntArray &globalNodeIdArray, TimeStep *tStep);
    /** Solves the refined whole problem.
     * @param localNodeIdArray array of local problem node ids
     * @param globalNodeIdArray array of global problem node ids
     * @param tStep time step
     */
    void solveRefinedWholeProblem(IntArray &localNodeIdArray, IntArray &globalNodeIdArray, TimeStep *tStep);
    /** Extracts nodal vector from global vector for each dof of all element nodes
     * @param element element
     * @param vector global vector
     * @param answer element nodal vector
     * @param dofs number of dofs at each node
     * @param tStep time step
     */
    void extractVectorFrom(Element *element, FloatArray &vector, FloatArray &answer, int dofs, TimeStep *tStep);

    void setupRefinedProblemProlog(const char *problemName, int problemId, IntArray &localNodeIdArray,
                                   int nodes, int elems, int csects, int mats, int loads, int ltfuncs,
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
    // mode for problem setup
    enum SetupMode { CountMode = 0, NodeMode = 1, ElemMode = 2, BCMode = 3 };

public:
    /// Constructor
    HuertaErrorEstimatorInterface() { }

    /// Returns reference to corresponding element
    virtual Element *HuertaErrorEstimatorI_giveElement() = 0;

    virtual void HuertaErrorEstimatorI_setupRefinedElementProblem(RefinedElement *refinedElement, int level, int nodeId,
                                                                  IntArray &localNodeIdArray, IntArray &globalNodeIdArray,
                                                                  HuertaErrorEstimatorInterface :: SetupMode mode, TimeStep *tStep,
                                                                  int &localNodeId, int &localElemId, int &localBcId,
                                                                  IntArray &controlNode, IntArray &controlDof,
                                                                  HuertaErrorEstimator :: AnalysisMode aMode) = 0;

    virtual void HuertaErrorEstimatorI_computeLocalCoords(FloatArray &answer, const FloatArray &coords) = 0;
    virtual void HuertaErrorEstimatorI_computeNmatrixAt(GaussPoint *aGaussPoint, FloatMatrix &answer) = 0;

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

    void setupRefinedElementProblem3D(Element * element, RefinedElement * refinedElement,
                                      int level, int nodeId, IntArray & localNodeIdArray, IntArray & globalNodeIdArray,
                                      HuertaErrorEstimatorInterface :: SetupMode mode, TimeStep * tStep, int nodes,
                                      FloatArray * * corner, FloatArray * midSide, FloatArray * midFace, FloatArray & midNode,
                                      int & localNodeId, int & localElemId, int & localBcId,
                                      int hexaSideNode [ 1 ] [ 3 ], int hexaFaceNode [ 1 ] [ 3 ],
                                      IntArray & controlNode, IntArray & controlDof,
                                      HuertaErrorEstimator :: AnalysisMode aMode, const char *hexatype);
};


/**
 * The class representing Huerta remeshing criteria.
 * The basic task is to evaluate the required mesh density (at nodes) on given domain,
 * based on informations provided by the compatible error ertimator.
 *
 * The remeshing criteria is maintained by the corresponding error estimator. This is mainly due to fact, that is
 * necessary for given EE to create compatible RC. In our concept, the EE is responsible.
 *
 */
class HuertaRemeshingCriteria : public RemeshingCriteria
{
public:
    /// mode of receiver, allows to use it in more general situations
    enum HuertaRemeshingCriteriaModeType { primaryUnknownBased };

protected:
    /// Array of nodal mesh densities
    FloatArray nodalDensities;
    /// Remeshing strategy proposed
    RemeshingStrategy remeshingStrategy;
    /// actual values (densities) state counter.
    StateCounterType stateCounter;
    /// mode of receiver
    HuertaRemeshingCriteriaModeType mode;
    /// required error to obtain
    double requiredError;
    /// minimum element size alloved
    double minElemSize;
    /// refinement coefficient
    double refineCoeff;
    /// remeshing flag
    bool noRemesh;
    /// weighted error flag
    bool wError;

public:
    /// Constructor
    HuertaRemeshingCriteria(int n, ErrorEstimator *e);
    /// Destructor
    ~HuertaRemeshingCriteria() { }

    /** Returns the required mesh size n given dof manager.
     * The mesh density is defined as a required element size
     * (in 1D the element length, in 2D the square from element area).
     * @param num dofman  number
     * @param tStep time step
     * @param relative if zero, then actual density is returned, otherwise the relative density to current is returned.
     */
    virtual double giveRequiredDofManDensity(int num, TimeStep *tStep, int relative = 0);
    /**
     * Returns existing mesh size for given dof manager.
     * @param num dofMan number
     */
    virtual double giveDofManDensity(int num);
    /**
     * Determines, if the remeshing is needed, and if nedded, the type of strategy used
     */
    virtual RemeshingStrategy giveRemeshingStrategy(TimeStep *tStep);
    /**
     * Estimates the nodal densities.
     * @param tStep time step
     */
    virtual int estimateMeshDensities(TimeStep *tStep);
    /** Initializes receiver acording to object description stored in input record.
     * This function is called immediately after creating object using
     * constructor. InitString can be imagined as data record in component database
     * belonging to receiver. Receiver use value-name extracting functions
     * to extract particular field from record.
     * @see readInteger, readDouble and similar functions */
    virtual IRResultType initializeFrom(InputRecord *ir);

    /// Returns "HuertaErrorEstimator" - class name of the receiver.
    const char *giveClassName() const { return "HuertaErrorEstimator"; }
    /** Returns HuertaRemeshingCriteriaClass - classType id of receiver.
     * @see FEMComponent::giveClassID
     */
    classType                giveClassID() const { return HuertaRemeshingCriteriaClass; }

protected:
};

/**
 * The corresponding element interface to HuertaRemeshingCriteria class.
 * Declares the necessary services, which have to be provided by particular elements.
 */

class HuertaRemeshingCriteriaInterface : public Interface
{
public:
    /// Constructor
    HuertaRemeshingCriteriaInterface() : Interface() { }
    /**
     * Determines the characteristic size of element. This quantity is defined as follows:
     * For 1D it is the element length, for 2D it is the square root of element area.
     */
    virtual double HuertaRemeshingCriteriaI_giveCharacteristicSize() = 0;
    /**
     * Returns the polynomial order of receiver trial functions.
     */
    virtual int HuertaRemeshingCriteriaI_givePolynOrder() = 0;
};
} // end namespace oofem
#endif // huertaerrorestimator_h
