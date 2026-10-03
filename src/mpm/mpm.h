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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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
#ifndef mpm_h
#define mpm_h

/**
 * Multiphysics module 
 * Classes:
 * - ElementBase(Element) defining geometry
 * - Variable class representing unknown field (or test feld) in a weak psolution. The variable has its interpolation, type (scalar, vector), size.
    When test field, it keeps reference to its primary (unknown) variable. The history parameter dermines how many time steps to remember. 
 * - Term class represnting a term to evaluate on element. Paramaters element(geometry), variables
 * - Element - responsible for defining and performing integration (of terms), assembly of term contributions. 
 */

#include "element.h"
#include "dofmanager.h"
#include "gausspoint.h"
#include "feinterpol.h"
#include "intarray.h"
#include "masterdof.h"
#include "integrationrule.h"
#include "gaussintegrationrule.h"
#include "classfactory.h"
#include "enum.h"
#include "dofiditem.h"
#include "stateoperator.h"
#include "stateoperator.h"

#include <map>


namespace oofem {

class MPElement;
class EngngModel;

/* 
 * Note: someone should be able to return for given cell and variable vector of unknowns.
 * this depends on interpolation (constant, linear, etc), cell type and variable (defined the physical meaning of unknown(s))
 * Interpolation should identify (or even introduce) cell nodes needed (quadratic element, linear interpolation), variable should assign to these nodes DOFs.
 * interpolation.getCellNodes(cell)
 */

/**
 * @brief Class representing unknown field (or test field) in a weak solution.
 * The variable has its interpolation, type (scalar, vector), size.
 * When test (dual) field, it keeps reference to its primary (unknown) variable.
 * @todo The history parameter determines how many time steps to remember. 
 */

#define ENUM_TYPE VariableType
#define ENUM_DEF ENUM_ITEM(scalar) ENUM_ITEM(vector)
#define ENUM_CLASS
#include "enum-impl.h"

class OOFEM_EXPORT Variable {
    public:
    typedef oofem::VariableType VariableType;

    std::string name;
    const FEInterpolation* interpolation;
    /**
     * When the receiver is a test (weighting) function, the unknown field it weights; null when
     * the receiver is itself an unknown.
     *
     * Declared in the input as `dualto "<name>"` on the test variable and resolved by
     * EngngModel::instanciateMPM. This is the only thing that distinguishes the two roles: a deck
     * declares them as two Variable records that differ just by name, and may give them different
     * interpolations (a non-symmetric, Petrov-Galerkin formulation). Assembly needs the
     * distinction because nodal unknowns have to be read through the unknown field's
     * interpolation, not its weighting function's.
     */
    Variable* dualVar;
    /// Name given by `dualto`, resolved into @ref dualVar after all variables have been read.
    std::string dualVarName;
    VariableType type;
    /// What the field physically is. Core-level, so a material can name it in its state layout.
    FieldType q;
    int size;
    IntArray dofIDs;

    Variable () : name(), interpolation(nullptr), dualVar(NULL), type(VariableType::scalar), q(FT_Displacements), size(0) {}
    Variable (const FEInterpolation* i, FieldType q, Variable::VariableType t, int size, Variable* dual = NULL, std :: initializer_list< int > dofIDs={}, std::string name="") : 
        interpolation(i), 
        dualVar(dual), 
        q(q), 
        dofIDs(dofIDs) {
        this->type = t;
        this->size = size;
    }
    Variable (const FEInterpolation* i, FieldType q, Variable::VariableType t, int size, IntArray& dofIDs, Variable* dual = NULL, std::string name="") : 
        name(name),
        interpolation(i), 
        dualVar(dual), 
        q(q), 
        dofIDs(dofIDs) {
        this->type = t;
        this->size = size;
    }


    /**
     * Completes the set-up of the receiver once every variable of the problem has been read.
     *
     * Resolves the `dualto` name into @ref dualVar, so the declaration order of the Variable
     * records does not matter. Only relevant to variables built from an input record; those
     * constructed directly (from python, or the hardwired fields of the classic up/tm elements)
     * pass their dual to the constructor.
     */
    void postInitialize(EngngModel *problem);

    /// True if the receiver is a test (weighting) function of some unknown field.
    bool isTestField() const { return this->dualVar != nullptr; }

    /// Returns DodIF mask in node; need generalization (which dofMan)
    const IntArray& getDofManDofIDs () const {return this->dofIDs;}

    void initializeFrom(const std::shared_ptr<InputRecord> &ir); // enable instantiation from input record
};


/**
 * @brief Class representing a weak form expression to be evaluated (integrated).
 * It defines two key methods:
 * - evaluate method to evaluate the term value, all unknowns known
 * - evaluate_lin to evaluate the consistent linearization of the term, so if Term is T(u), depending on unknown u,
 *   this term evaluates dT/du, which typically contributes to the LHS. 
 */
class OOFEM_EXPORT Term {
    public:
    const Variable* field;
    const Variable* testField;
    MaterialMode mode;

    public:
    Term () : field(nullptr), testField(nullptr), mode(MaterialMode::_Unknown) {}
    Term (const Variable* testField, const Variable* unknownField, MaterialMode m=MaterialMode::_Unknown) : field(unknownField), testField(testField) {mode=m;}
    
    // evaluate linearized term contribution to weak form on given cell at given point 
    virtual void evaluate_lin (FloatMatrix& , MPElement& cell, GaussPoint* gp, TimeStep* tStep) const =0;
    // evaluate contribution (all vars known) on given cell
    virtual void evaluate (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tStep) const =0;
    virtual void getDimensions(Element& cell) const =0;
    virtual void initializeCell(Element& cell) const =0;
    virtual IntegrationRule* giveElementIntegrationRule(Element* e) const {return NULL;};

    virtual void initializeFrom(const std::shared_ptr<InputRecord> &ir, EngngModel* problem); // enable instantiation from input record

};


/**
 * MPMSymbolic terms extend standard Terms to allow for cell initialization.
 * The symbolic terms are assumed to be evaluated on generic cells (and not on problem-specific elements).
 * Therefore the need to ensure that proper DOFs and integration rules are set-up.
 */
class OOFEM_EXPORT GenericCellTerm : public Term {
    protected:
    int nip=0; // assumed order of interpolation for the term
    public:
    GenericCellTerm() : Term() {}
    GenericCellTerm (const Variable *testField, const Variable* unknownField, MaterialMode m)  : Term(testField, unknownField, m) {};
    void initializeFrom(const std::shared_ptr<InputRecord> &ir, EngngModel* problem) override {
        Term::initializeFrom(ir, problem);
        IR_GIVE_OPTIONAL_FIELD(ir, nip, "nip");
    }
    void initializeCell(Element& cell) const override {
        // initialize cell for interpolation use
        // @TODO: prevent multiple initialization for same interpolation
        this->field->interpolation->initializeCell(&cell);
        this->testField->interpolation->initializeCell(&cell);

        // allocate necessary DOFs
        IntArray enodes, einteranlnodes, dofIDs;
        // process term field
        dofIDs = this->field->getDofManDofIDs();
        this->field->interpolation->giveCellDofMans(enodes, einteranlnodes, &cell);
        for (auto i: enodes) {
            DofManager* dman =  cell.giveDofManager(i);
            for (auto d: dofIDs) {
                if (!dman->hasDofID((DofIDItem) d)) {
                    // create a DOF
                    MasterDof* dof = new MasterDof(dman, (DofIDItem)d);
                    dman->appendDof(dof);
                }
            }
        }
        for (auto i: einteranlnodes) {
            DofManager* dman =  cell.giveInternalDofManager(i);
            for (auto d: dofIDs) {
                if (!dman->hasDofID((DofIDItem) d)) {
                    // create a DOF
                    MasterDof* dof = new MasterDof(dman, (DofIDItem)d);
                    dman->appendDof(dof);
                }
            }
        }
        // process testField
        dofIDs = this->testField->getDofManDofIDs();
        this->testField->interpolation->giveCellDofMans(enodes, einteranlnodes, &cell);
        for (auto i: enodes) {
            DofManager* dman =  cell.giveDofManager(i);
            for (auto d: dofIDs) {
                if (!dman->hasDofID((DofIDItem)d)) {
                    // create a DOF
                    MasterDof* dof = new MasterDof(dman, (DofIDItem)d);
                    dman->appendDof(dof);
                }
            }
        }
        for (auto i: einteranlnodes) {
            DofManager* dman =  cell.giveInternalDofManager(i);
            for (auto d: dofIDs) {
                if (!dman->hasDofID((DofIDItem)d)) {
                    // create a DOF
                    MasterDof* dof = new MasterDof(dman, (DofIDItem)d);
                    dman->appendDof(dof);
                }
            }
        }
        // set up the integration rule on cell
        // get required number of IPs
        int myorder = this->field->interpolation->giveInterpolationOrder() *  this->testField->interpolation->giveInterpolationOrder(); 
        GaussIntegrationRule ir(0, &cell);
        int nip = ir.getRequiredNumberOfIntegrationPoints(cell.giveIntegrationDomain(), myorder);
        if (this->nip>0) {
            nip = this->nip;
        }
        // create nd insert it toelement if not exist yet.
        std::vector< std :: unique_ptr< IntegrationRule > > &irvec = cell.giveIntegrationRulesArray();
        bool found = false;
        int size = irvec.size();
        for (int i = 0; i< size; i++) {
            if (irvec[i].get()->giveNumberOfIntegrationPoints() == nip) {
                found = true;
                break;
            }
        }
        if (!found) {
            // need to insert right one
            irvec.resize( size +1);
            irvec [ size] = std::make_unique<GaussIntegrationRule>(size, &cell);
            //irvec [ size ]->SetUpPointsOnSquare(nip, this->mode);
            irvec[size]->setUpIntegrationPoints(cell.giveIntegrationDomain(), nip, this->mode);
            OOFEM_LOG_INFO("Integration rule with %d nip created for cell %d\n",nip,cell.giveNumber());
        }
    }
    IntegrationRule* giveElementIntegrationRule(Element* e) const override {
        int myorder = this->field->interpolation->giveInterpolationOrder() *  this->testField->interpolation->giveInterpolationOrder(); 
        GaussIntegrationRule ir(0, e);
        int nip = ir.getRequiredNumberOfIntegrationPoints(e->giveIntegrationDomain(), myorder);
        if (this->nip>0) {
            nip = this->nip;
        }
        std::vector< std :: unique_ptr< IntegrationRule > > &irvec = e->giveIntegrationRulesArray();
        int size = irvec.size();
        for (int i = 0; i< size; i++) {
            if (irvec[i].get()->giveNumberOfIntegrationPoints() == nip) {
                return irvec[i].get();
            }
        }
        return NULL;
    }
};



/**
* Element code sample:
* term=Poison(Variable(interpolation, Temperature, 1), Variable(Interpolation, Temperature,1));
* this-> Assemble(integrationRule, term, destination); // where to integrate (volume, surface,edge?)
* 
*/

/**
 * @brief Base class for elements based on mp (multi-physics) concept
 * 
 */
class OOFEM_EXPORT MPElement : public Element {
   public:

    MPElement (int n, Domain * aDomain) : 
        Element(n, aDomain)
    {}

    void initialize () {
        // loop over variables and allocate nodal dofs (for unknownFields)
    }

    void integrateTerm_dw (FloatMatrix& answer, const Term& term, IntegrationRule* iRule, TimeStep* tstep) {
        // need for integration domain and rule.
        // who should determine integration domain? Element or term? Term is just integrand, not integral
        // so integral type (surface, volume, etc) defined by element ---
        FloatMatrix dw;
        for ( GaussPoint *igp : * iRule ) {
            term.evaluate_lin(dw, *this, igp, tstep);
            dw.times(this->computeVolumeAround(igp));
            answer.add(dw);
        }
    } 

    void integrateTerm_c (FloatArray& answer, const Term& term, IntegrationRule* iRule, TimeStep* tstep) {
        // need for integration domain and rule.
        // who should determine integration domain? Element or term? Term is just integrand, not integral
        // so integral type (surface, volume, etc) defined by element ---
        FloatArray dw;
        for ( GaussPoint *igp : * iRule ) {
            term.evaluate(dw, *this, igp, tstep);
            dw.times(this->computeVolumeAround(igp));
            answer.add(dw);
        }
    } 

    void integrateSurfaceTerm_dw (FloatMatrix& answer, const Term& term, IntegrationRule* iRule, int isurf, TimeStep* tstep) {
            // need for integration domain and rule.
            // who should determine integration domain? Element or term? Term is just integrand, not integral
            // so integral type (surface, volume, etc) defined by element ---
            FloatMatrix dw;
            for ( GaussPoint *igp : * iRule ) {
                term.evaluate_lin(dw, *this, igp, tstep);
                dw.times(this->computeSurfaceVolumeAround(igp, isurf));
                answer.add(dw);
            }
        } 

    void integrateSurfaceTerm_c (FloatArray& answer, const Term& term, IntegrationRule* iRule, int isurf, TimeStep* tstep) {
            // need for integration domain and rule.
            // who should determine integration domain? Element or term? Term is just integrand, not integral
            // so integral type (surface, volume, etc) defined by element ---
            FloatArray dw;
            for ( GaussPoint *igp : * iRule ) {
                term.evaluate(dw, *this, igp, tstep);
                dw.times(this->computeSurfaceVolumeAround(igp,isurf));
                answer.add(dw);
            }
        } 

    void integrateEdgeTerm_dw (FloatMatrix& answer, const Term& term, IntegrationRule* iRule, int iedge, TimeStep* tstep) {
            // need for integration domain and rule.
            // who should determine integration domain? Element or term? Term is just integrand, not integral
            // so integral type (surface, volume, etc) defined by element ---
            FloatMatrix dw;
            for ( GaussPoint *igp : * iRule ) {
                term.evaluate_lin(dw, *this, igp, tstep);
                dw.times(this->computeEdgeVolumeAround(igp, iedge));
                answer.add(dw);
            }
        } 

    void integrateEdgeTerm_c (FloatArray& answer, const Term& term, IntegrationRule* iRule, int iedge, TimeStep* tstep) {
            // need for integration domain and rule.
            // who should determine integration domain? Element or term? Term is just integrand, not integral
            // so integral type (surface, volume, etc) defined by element ---
            FloatArray dw;
            for ( GaussPoint *igp : * iRule ) {
                term.evaluate(dw, *this, igp, tstep);
                dw.times(this->computeEdgeVolumeAround(igp, iedge));
                answer.add(dw);
            }
        } 

    /**
     *  Returns local code numbers corresponding to specific variable identified by variablequantity.
     *  This essentialy allows to assemble partial element contributions separately.
     *  It is assumed, that partial element contributions are obtained using standard 
     *  getCharacteristicMatrix and getCharacteristicVector methods. 
     * @param answer code numbers corrresponding to given variable
     * @param q variable type 
     */
  virtual void getDofManLocalCodeNumbers(IntArray& answer, const FieldType q, int n) const = 0;
  virtual void getInternalDofManLocalCodeNumbers(IntArray& answer, const FieldType q, int n) const = 0;
  virtual void getLocalCodeNumbers (IntArray& answer, const FieldType q ) const {
    IntArray dl;
    answer.resize(0);
    
    for (int i=1; i<= this->giveNumberOfDofManagers(); i++) {
      this->getDofManLocalCodeNumbers(dl, q, i);
      answer.followedBy(dl);
    }
    for (int i=1; i<= this->giveNumberOfInternalDofManagers(); i++) {
      this->getInternalDofManLocalCodeNumbers(dl, q, i);
      answer.followedBy(dl);
    }

  }
  virtual int getNumberOfSurfaceDOFs() const =0;
  virtual int getNumberOfEdgeDOFs() const =0;


  /**
     Returns mapping from quantity dofs to local surface dofs
  */
  virtual void getSurfaceLocalCodeNumbers (IntArray& answer, const FieldType q) const =0;
  virtual void getEdgeLocalCodeNumbers (IntArray& answer, const FieldType q) const =0;
  /** @brief Returns element code numbers of the unknowns associated with given boundary entity. 
   * @param answer 
   * @param q 
   * @param isurf
   */ 
  virtual void getSurfaceElementCodeNumbers (IntArray& answer, const FieldType q, int isurf ) const {
    IntArray dl, sn = this->getGeometryInterpolation()->boundarySurfaceGiveNodes(isurf, this->giveGeometryType());
    answer.resize(0);
    for (int i : sn) {
      this->getDofManLocalCodeNumbers(dl, q, i);
      answer.followedBy(dl);
    }
  }
  virtual void getEdgeElementCodeNumbers (IntArray& answer, const FieldType q, int isurf ) const {
    IntArray dl, sn = this->getGeometryInterpolation()->boundaryEdgeGiveNodes(isurf, this->giveGeometryType());
    answer.resize(0);
    for (int i : sn) {
      this->getDofManLocalCodeNumbers(dl, q, i);
      answer.followedBy(dl);
    }
  }

    /** 
     * Returns boundary entity unknown vector
     * @param ibc boundary entity ID
     * @param bt boundary type ('s' for surface, 'e' for edge)
    */
    virtual void getBoundaryUnknownVector(FloatArray& answer, const Variable* field, ValueModeType mode, int ibc, char bt, TimeStep* tStep) {
        FloatArray uloc;
        IntArray bNodes, dofs=field->getDofManDofIDs();
        answer.clear();
        if (bt == 's') {
            bNodes = this->giveBoundarySurfaceNodes(ibc);
        } else {
            bNodes = this->giveBoundaryEdgeNodes(ibc);
        }
        std::list<double> ans;
        for (int i : bNodes) {
            this->giveDofManager(i)->giveUnknownVector(uloc, dofs, mode, tStep);
            for(const double& u: uloc){ ans.push_back(u); }
        }
        answer=FloatArray::fromList(ans);
    }
  
    /// @brief  Assembles the partial element contribution into local element matrix
    /// @param answer 
    /// @param contrib 
    /// @param t 
    void assembleTermContribution (FloatMatrix& answer, FloatMatrix& contrib, const Term& t) {
        IntArray uloc, tloc;
        this->getLocalCodeNumbers(uloc, t.field->q);
        this->getLocalCodeNumbers(tloc, t.testField->q);
        answer.assemble(contrib, tloc, uloc);
    }

    void assembleTermContributionT (FloatMatrix& answer, FloatMatrix& contrib, const Term& t) {
        IntArray uloc, tloc;
        this->getLocalCodeNumbers(uloc, t.field->q);
        this->getLocalCodeNumbers(tloc, t.testField->q);
        answer.assembleT(contrib, tloc, uloc);
    }
    
    void assembleTermContribution (FloatArray& answer, FloatArray& contrib, const Term& t) {
        IntArray loc;
        this->getLocalCodeNumbers(loc, t.testField->q);
        answer.assemble(contrib, loc);
    }

    /**
     * @brief Returns vector of nodal unknowns for given Variable
     * 
     * @param answer 
     * @param field 
     * @param tstep 
     */
    virtual const void getUnknownVector(FloatArray& answer, const Variable* field, ValueModeType mode, TimeStep* tstep) {
        IntArray nodes, internalNodes, dofs;
        field->interpolation->giveCellDofMans(nodes, internalNodes, this);
        std::vector<double> ans;
        for (int i : nodes) {
            dofs=field->getDofManDofIDs();
            FloatArray uloc;
            this->giveDofManager(i)->giveUnknownVector(uloc, dofs, mode, tstep);
            ans.insert(ans.end(),uloc.begin(),uloc.end());
        }
        for (int i : internalNodes) {
            dofs=field->getDofManDofIDs();
            FloatArray uloc;
            this->giveInternalDofManager(i)->giveUnknownVector(uloc, dofs, mode, tstep);
            ans.insert(ans.end(),uloc.begin(),uloc.end());
        }
        answer = FloatArray::fromVector(ans);
    }

    /**
     * @name Interpolation operators for a Variable at an integration point.
     *
     * These are the shared implementations of the operators the symbolic terms expose as Grad_s,
     * Grad and N. They live here rather than in the symbolic layer so that both the symbolic
     * functors and assembleStateVector use one definition.
     */
    //@{
    /// Symmetric gradient (strain) operator of a vector field; shape depends on the material mode.
    void computeGradSymMatrixAt(FloatMatrix &answer, const Variable *v, GaussPoint *gp) const;
    /// Gradient operator of a scalar field; nsd rows.
    void computeGradMatrixAt(FloatMatrix &answer, const Variable *v, GaussPoint *gp) const;
    /// Interpolation (shape function) operator of a field.
    void computeNMatrixAt(FloatMatrix &answer, const Variable *v, GaussPoint *gp) const;
    /// Divergence operator of a vector field.
    void computeDivMatrixAt(FloatMatrix &answer, const Variable *v, GaussPoint *gp) const;
    /// Dispatches to the operator named by @p op.
    void computeStateOperatorAt(FloatMatrix &answer, StateOperator op, const Variable *v, GaussPoint *gp) const;
    //@}

    /// Maps a field, by what it physically is, onto the Variable supplying it.
    typedef std::map<int, const Variable *> StateVariableMap;

    /**
     * Records the given primary field as a source of state on the receiver, under what it
     * physically is (Variable::q).
     *
     * Called from Integral::initialize for the unknown field of each term acting on the receiver,
     * alongside the creation of the dofs and integration rules that term needs. Resolution is
     * therefore per cell, which is what lets a domain carry several materials with different state
     * layouts, each fed by its own fields.
     *
     * Test (weighting) fields must not be registered: they are declared with the same dof ids as
     * the field they weight, and may carry a different interpolation. Variable::isTestField says
     * which is which.
     */
    void registerStateVariable(const Variable *v);

    /// Returns the primary field registered for the given physical field, or nullptr.
    const Variable *giveStateVariableSource(FieldType field) const {
        auto it = this->stateVariables.find( (int) field );
        return ( it == this->stateVariables.end() ) ? nullptr : it->second;
    }

    /**
     * Assembles the generalized state vector described by @p layout at the given point, from the
     * fields registered on the receiver.
     *
     * This is the C++ counterpart of the packing input decks used to do by hand, in expressions of
     * the form "flux = vcat(eps, pw, pa)". Each entry names a field and the operator to apply to
     * it.
     *
     * @param layout Layout advertised by the material through giveStateVariableIDs.
     */
    void assembleStateVector(FloatArray &answer, const StateVariableLayout &layout, GaussPoint *gp, TimeStep *tStep);

    /**
     * Pushes the current state to the material at every integration point of the receiver.
     *
     * Materials advertising an empty layout are skipped -- they do not participate in the
     * push/pull protocol and keep using their physics-specific entry points.
     */
    void updateTempState(TimeStep *tStep);

protected:
    /// Which primary field supplies each state quantity on the receiver; see registerStateVariable.
    StateVariableMap stateVariables;

public:

  virtual double computeSurfaceVolumeAround(GaussPoint* igp, int iSurf)
  {return igp->giveWeight()*this->getGeometryInterpolation()->boundarySurfaceGiveTransformationJacobian(iSurf, igp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this));}
  virtual double computeEdgeVolumeAround(GaussPoint* igp, int iEdge) 
  {return igp->giveWeight()*this->getGeometryInterpolation()->boundaryEdgeGiveTransformationJacobian(iEdge, igp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this));}
  virtual double computeVolumeAround(GaussPoint* igp) override
  {return igp->giveWeight()*this->getGeometryInterpolation()->giveTransformationJacobian(igp->giveNaturalCoordinates(), FEIElementGeometryWrapper(this));}
  
   //const FEInterpolation& getGeometryInterpolation () const override = 0;
   FEInterpolation *giveInterpolation() const override { return const_cast<FEInterpolation*>( this->getGeometryInterpolation()); }

    /**
     * Returns transformation matrix from local boundary (edge/surface) c.s  to element local coordinate system
     * of load/flux vector components. If no transformation is necessary, answer is empty matrix (default);
     * @param answer Computed rotation matrix.
     * @param iSurf Surface/edge number.
     * @param lc local coordinates at which to compute transformation
     * @param dofIDs dofID ask identifying relevant DOFs
     * @param btype boundary type ('s'for surface, 'e' for edge)
     * 
     * @return Nonzero if transformation matrix is not empty matrix, zero otherwise.
     * 
     */
    virtual int computeFluxLBToLRotationMatrix(FloatMatrix &answer, int iSurf, const FloatArray& lc, const FieldType q, char btype) {
        answer.clear(); 
        return 0;
    }
    IntArray giveBoundarySurfaceNodes(int boundary, bool includeHierarchical=false) const override {
        return this->getGeometryInterpolation()->boundarySurfaceGiveNodes(boundary, this->giveGeometryType(), includeHierarchical);
    }
    IntArray giveBoundaryEdgeNodes(int boundary, bool includeHierarchical=false) const override {
        return this->getGeometryInterpolation()->boundaryEdgeGiveNodes(boundary, this->giveGeometryType(), includeHierarchical);
    }
    virtual void giveCharacteristicMatrixFromBC(FloatMatrix &answer, CharType type, TimeStep *tStep, GeneralBoundaryCondition *bc, int boundaryID) {
        answer.clear();
    }
    virtual void giveCharacteristicVectorFromBC(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep, GeneralBoundaryCondition *bc, int boundaryID) {
        answer.clear();
    }

};




} // end namespace oofem
#endif // mpm_h
