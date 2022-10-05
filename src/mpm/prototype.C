/**
 * Multiphysics element template 
 * Classes:
 * - ElementBase(Element) defining geometry
 * - Variable class representing unknown field (or test feld) in a weak psolution. The variable has its interpolation, type (scalar, vector), size.
    When test field, it keeps reference to its primary (unknown) variable. The history parameter dermines how many time steps to remember. 
 * - Term class represnting a term to evaluate on element. Paramaters element(geometry), variables
 * - Element - responsible for defining and performing integration (of terms), assembly of term contributions. 
 */

#include "element.h"
#include "gausspoint.h"
#include "feinterpol.h"
#include "intarray.h"
#include "classfactory.h"

// for demo
#include "fei2dtrlin.h"
#include "gaussintegrationrule.h"


namespace oofem {

class Variable {
    public:
    enum VariableType {
        scalar,
        vector
    };

    enum VariableQuantity {
        Displacement,
        Temperature
    };

    FEInterpolation& interpolation;  
    Variable* dualVar; //? or just bool?
    VariableType type;
    VariableQuantity q;
    int size;
    IntArray dofIDs;

    Variable (FEInterpolation& i, Variable::VariableQuantity q, Variable::VariableType t, int size, Variable* dual = NULL) : 
        interpolation(i), 
        dualVar(dual), 
        q(q) {
        this->type = t;
        this->size = size;
    }

    /// Returns DodIF mask in node
    const IntArray& getDofManDofIDs () {return this->dofIDs;}
};

class Term {
    public:
    Variable& field;
    Variable& testField;

    public:
    Term (Variable& unknownField, Variable &testField) : field(unknownField), testField(testField) {}
    
    // evaluate term contribution to weak form on given cell at given point 
    virtual void evaluate_dw (FloatMatrix& , Element& cell, const FloatArray& coodrs) =0;
    // evaluate contribution (all vars known) on given cell
    virtual void evaluate_c (FloatArray&, Element& cell, const FloatArray& coords)=0;
    virtual void getDimensions_dw(Element& cell) =0;
    virtual void initializeCell(Element& cell) =0;
};

class PoissonTerm : public Term {
    protected:
        double c;
    public:
    PoissonTerm (Variable& unknownField, Variable &testField, double c) : Term(unknownField, testField) {
        this->c = c;
    }

    void grad(FloatMatrix& answer, Variable &v, FEInterpolation& interpol, Element& cell, const FloatArray& coords) {
        interpol.evaldNdx(answer, coords, FEIElementGeometryWrapper(&cell));
    }


    void evaluate_dw (FloatMatrix& answer, Element& e, const FloatArray& coords) override {
        FEInterpolation & si = field.interpolation;
        FEInterpolation &ti = testField.interpolation;
        FloatMatrix bs, bt;
        this->grad(bs, this->field,si,e,coords);
        this->grad(bt, this->testField,ti,e,coords);

        FloatMatrix gc, c(2,2);
        c.at(1,1) = this->c;
        c.at(2,2) = this->c;

        gc.beProductOf(bs, c);
        answer.beProductTOf(gc, bt);
    }

    void evaluate_c (FloatArray&, Element& cell, const FloatArray& coords) override {}
    void getDimensions_dw(Element& cell) override {}
    void initializeCell(Element& cell) override {}
};


/**
* Element code sample:
* term=Poison(Variable(interpolation, Temperature, 1), Variable(Interpolation, Temperature,1));
* this-> Assemble(integrationRule, term, destination); // where to integrate (volume, surface,edge?)
* 
*/
class MPElement : public Element {
   public:

    MPElement (int n, Domain * aDomain) : 
        Element(n, aDomain)
    {}

    void initialize () {
        // loop over variables and allocate nodal dofs (for unknownFields)
    }

    void integrateTerm_dw (FloatMatrix& answer, Term& term, IntegrationRule* iRule) {
        // need for integration domain and rule.
        // who should determine integration domain? Element or term? Term is just integrand, not integral
        // so integral type (surface, volume, etc) defined by element ---
        FloatMatrix dw;
        for ( GaussPoint *igp : * iRule ) {
            term.evaluate_dw(dw, *this, igp->giveNaturalCoordinates());
            answer.add(dw);
        }
    } 

    /// @brief  returns local code number on element to assemble variable/term contribution
    /// @param answer 
    /// @param v 
    virtual void getLocalCodeNumbers (IntArray& answer, Variable& v) =0;
    

    void assembleTermContribution (FloatMatrix& answer, FloatMatrix& contrib, Term& t) {
        IntArray uloc, tloc;
        this->getLocalCodeNumbers(uloc, t.field);
        this->getLocalCodeNumbers(tloc, t.testField);
        answer.assemble(contrib, uloc, tloc);
    }
};

#define _IFT_PoissonElement_Name "pe"


class PoissonElement : public MPElement {
     protected:
        FEI2dTrLin interpol;
        Variable t;
        Variable dt;
        PoissonTerm p;
        GaussIntegrationRule ir;
    public:
    PoissonElement(int n, Domain* d): 
        MPElement(n,d), 
        interpol(1,2), 
        t(interpol, Variable::VariableQuantity::Temperature, Variable::VariableType::scalar, 3), 
        dt(interpol, Variable::VariableQuantity::Temperature, Variable::VariableType::scalar, 3, &t),
        p(t,dt,1.0),
        ir(1, this)
    {
        numberOfDofMans  = 3;
        numberOfGaussPoints = 1;
        ir.SetUpPointsOnTriangle(numberOfGaussPoints, _2dHeat);
    }


    void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep) override {
        if (type == ConductivityMatrix) {
            FloatMatrix term;
            answer.resize(3,3);
            this->integrateTerm_dw (term, this->p, &this->ir) ;
            this->assembleTermContribution(answer, term, this->p);
            answer.printYourself("Conductivity");

        }
    }
    void getLocalCodeNumbers (IntArray& answer, Variable& v) override {
        answer.enumerate(this->giveNumberOfDofManagers());
    }
    
    const char *giveInputRecordName() const override {return "pe";}
};

REGISTER_Element(PoissonElement)


} // end namespace oofem