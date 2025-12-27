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
#ifndef termlibrary_h
#define termlibrary_h

#include "mpm.h"
#include "boundaryload.h"
#include "bodyload.h"

namespace oofem {

#define _IFT_BTamNTerm_Name "BTamNTerm"
#define _IFT_NTamTBTerm_Name "NTamTBTerm"
#define _IFT_NTcN_Name "NTcN"


/**
 * @brief A Linear momentum balance equation term ($B^T\sigma(u)$)
 * 
 */
class BTSigTerm : public Term {
    protected:
    public:
    BTSigTerm (const Variable *testField, const Variable* unknownField) ;

    /**
     * @brief Evaluates the linearization of $B^T\sigma(u)$, i.e. $B^TDBu$
     * 
     * @param answer 
     * @param e 
     * @param coords 
     */
    void evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const override;
    /**
     * @brief Evaluates Internal forces vector, i.e. $b^T\sigma(u)$
     * 
     * @param cell 
     * @param coords 
     */
    void evaluate (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const override;
    void getDimensions(Element& cell) const override;
    void initializeCell(Element& cell) const override;

    protected:
    /**
     * @brief Evaluates B matrix; i.e. $LN$ where $L$ is operator matrix and $N$ is interpolation matrix of unknowns
     * 
     * @param answer B matrix
     * @param v 
     * @param interpol 
     * @param cell 
     * @param coords 
     */
    void grad(FloatMatrix& answer, const Variable *v, const FEInterpolation* interpol, const Element& cell, const FloatArray& coords, const MaterialMode mmode) const  ;
    
};

/**
 * @brief A continuity equation term ($Hp=(\grad N_p)^T f(p)$, where $f=\bf{k}/\mu \grad p$. 
 */
class gNTfTerm : public Term {
    protected:
    MatResponseMode lhsType, rhsType;
    public:
    gNTfTerm (const Variable *testField, const Variable* unknownField, MatResponseMode lhsType, MatResponseMode rhsType) ;

    /**
     * @brief Evaluates $\bf{H}$ matrix, the linearization of $w^T(\grad N)^T f(p)$, i.e. $(\grad N)^T \bf{k}/\mu \grad p = \bf{H}$
     * 
     * @param answer 
     * @param e 
     * @param coords 
     */
    void evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const override;
    /**
     * @brief Evaluates Internal forces vector, i.e. $w^T(\grad N)^T f(p)$
     * 
     * @param cell 
     * @param coords 
     */
    void evaluate (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const  override;
    void getDimensions(Element& cell) const  override;
    void initializeCell(Element& cell) const override;

    protected:
    /**
     * @brief Evaluates B matrix; i.e. $\grad N$ where $N$ is interpolation matrix of unknown (p)
     * 
     * @param answer B matrix
     * @param v 
     * @param interpol 
     * @param cell 
     * @param coords 
     */
    void grad(FloatMatrix& answer, const Variable *v, const FEInterpolation* interpol, const Element& cell, const FloatArray& coords) const ;
    
};

/**
 * @brief A continuity equation term $Qp=(B)^T \alpha\bf{m}N_p$. 
 */
class BTamNTerm : public MPMSymbolicTerm {
    protected:
        MatResponseMode aType = BiotConstant;
    public:
    BTamNTerm() : MPMSymbolicTerm() {}
    BTamNTerm (const Variable *testField, const Variable* unknownField, MatResponseMode at = BiotConstant) ;

    /**
     * @brief Evaluates the linearization of receiver, i.e. the LHS term
     * 
     * @param answer 
     * @param e 
     * @param coords 
     */
    void evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const override;
    /**
     * @brief Evaluates Internal forces vector, i.e. $w^T(\grad N)^T f(p)$
     * 
     * @param cell 
     * @param coords 
     */
    void evaluate (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const override;
    void getDimensions(Element& cell) const override;
    void initializeFrom(InputRecord &ir, EngngModel* problem) override {
        MPMSymbolicTerm::initializeFrom(ir, problem);
        int value = 0;
        IR_GIVE_FIELD(ir, value, "atype");
        aType = static_cast<MatResponseMode>(value);
    }


    protected:
    /**
     * @brief Evaluates B matrix; i.e. $\grad N$ where $N$ is interpolation matrix of unknown (p)
     * 
     * @param answer B matrix
     * @param v 
     * @param interpol 
     * @param cell 
     * @param coords 
     */
    void grad(FloatMatrix& answer, const Variable *v, const FEInterpolation* interpol, const Element& cell, const FloatArray& coords, const MaterialMode mmode) const ;
    
};

/**
 * @brief A continuity equation term $Q^T(du\over dt)=(N)^T \alpha\bf{m}^TB du\over dt$. 
 */
class NTamTBTerm : public MPMSymbolicTerm {
    protected:
        MatResponseMode aType = BiotConstant;
        ValueModeType unknownFieldVMT = VM_Velocity;
    public:
    NTamTBTerm () : MPMSymbolicTerm() {}
    NTamTBTerm (const Variable *testField, const Variable* unknownField, MatResponseMode at = BiotConstant, ValueModeType unknownFieldVMT = VM_Velocity) ;

    /**
     * @brief Evaluates the linearization of receiver, i.e. the LHS term
     * 
     * @param answer 
     * @param e 
     * @param coords 
     */
    void evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const override;
    /**
     * @brief Evaluates Internal forces vector, i.e. $w^T(\grad N)^T f(p)$
     * 
     * @param cell 
     * @param coords 
     */
    void evaluate (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const override;
    void getDimensions(Element& cell) const override;
    void initializeFrom(InputRecord &ir, EngngModel* problem) override {
        MPMSymbolicTerm::initializeFrom(ir, problem);
        int value = 0;
        IR_GIVE_FIELD(ir, value, "atype");
        aType = static_cast<MatResponseMode>(value);

        value = unknownFieldVMT; // VM_Velocity
        IR_GIVE_OPTIONAL_FIELD(ir, value, "uvmt" );
        unknownFieldVMT = static_cast<ValueModeType>(value);
    }
    protected:
    /**
     * @brief Evaluates B matrix; i.e. $\grad N$ where $N$ is interpolation matrix of unknown (p)
     * 
     * @param answer B matrix
     * @param v 
     * @param interpol 
     * @param cell 
     * @param coords 
     */
    void grad(FloatMatrix& answer, const Variable *v, const FEInterpolation* interpol, const Element& cell, const FloatArray& coords, const MaterialMode mmode) const ;
    
};


/**
 * @brief A continuity equation compressibility matrix $S=(N_p)^T c\ N_p$, where $c=({\alpha-n}\over{K_s}+{n}\over{K_w})$. 
 */
class NTcN : public MPMSymbolicTerm {
    protected:
        MatResponseMode ctype;
        ValueModeType unknownFieldVMT = VM_Velocity;
    public:
    NTcN () : MPMSymbolicTerm() {}
    NTcN (const Variable *testField, const Variable* unknownField, MatResponseMode ctype, ValueModeType uFieldVMT = VM_Velocity) ;

    /**
     * @brief Evaluates the linearization of term (the lhs contribution)
     * 
     * @param answer 
     * @param e 
     * @param coords 
     */
    void evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const override;
    /**
     * @brief Evaluates Internal forces vector, i.e. $w^T(\grad N)^T f(p)$
     * 
     * @param cell 
     * @param coords 
     */
    void evaluate (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const override;
    void getDimensions(Element& cell) const override;
    void initializeFrom(InputRecord &ir, EngngModel* problem) override {
        MPMSymbolicTerm::initializeFrom(ir, problem);
        int value = 0;
        IR_GIVE_FIELD(ir, value, "ctype");
        ctype = static_cast<MatResponseMode>(value);

        value = unknownFieldVMT; // VM_Velocity
        IR_GIVE_OPTIONAL_FIELD(ir, value, "uvmt" );
        unknownFieldVMT = static_cast<ValueModeType>(value);
    }
    
};

/**
 * @brief An external flux functor
 * 
 */
class NTfFunctor {
    public:
  virtual void evaluate(FloatArray& answer, const FloatArray& coords, MPElement& cell, const Variable *testField, TimeStep* tStep) const = 0;
};


class BoundaryFluxFunctor: public NTfFunctor {
    protected:
        BoundaryLoad *load;
        IntArray dofIDs;
        int isurf;
        char type; // indicates boundary type: 'e' for edge, 's' for surface
    public:
    BoundaryFluxFunctor(BoundaryLoad *load, int surf, const IntArray& dofIDs, char btype) : load(load), dofIDs(dofIDs), isurf(surf), type(btype) {}

    void evaluate(FloatArray& answer, const FloatArray& lcoords, MPElement& cell, const Variable *testField, TimeStep* tStep) const override {

        ValueModeType mode = VM_Total;
        if ( load->giveFormulationType() == Load :: FT_Entity ) {
            load->computeValues(answer, tStep, lcoords, dofIDs, mode);
        } else {
            FloatArray globalIPcoords;
            testField->interpolation->local2global(globalIPcoords, lcoords, FEIElementGeometryWrapper(&cell) );
            load->computeValues(answer, tStep, globalIPcoords, dofIDs, mode);
        }

        ///@todo Make sure this part is correct.
        // We always want the global values in the end, so we might as well compute them here directly:
        // transform force
        if ( load->giveCoordSystMode() == Load :: CST_Global ) {
            // then just keep it in global c.s
        } else {
            FloatMatrix T;
            // then to global c.s
            if ( cell.computeFluxLBToLRotationMatrix(T, isurf, lcoords, testField->q, type )) {
                answer.rotatedWith(T, 'n');
            }
        }
    }
};

class BodyFluxFunctor: public NTfFunctor {
    protected:
        BodyLoad *load;
        IntArray dofIDs;
    public:
    BodyFluxFunctor(BodyLoad *load, const IntArray& dofIDs) : load(load), dofIDs(dofIDs) {}

    void evaluate(FloatArray& answer, const FloatArray& lcoords, MPElement& cell, const Variable *testField, TimeStep* tStep) const override {

        ValueModeType mode = VM_Total;
        if ( load->giveFormulationType() == Load :: FT_Entity ) {
            load->computeValues(answer, tStep, lcoords, dofIDs, mode);
        } else {
            FloatArray globalIPcoords;
            testField->interpolation->local2global(globalIPcoords, lcoords, FEIElementGeometryWrapper(&cell) );
            load->computeValues(answer, tStep, globalIPcoords, dofIDs, mode);
        }

        ///@todo Make sure this part is correct.
        // We always want the global values in the end, so we might as well compute them here directly:
        // transform force
        if ( load->giveCoordSystMode() == Load :: CST_Global ) {
            // then just keep it in global c.s
        } else {
            OOFEM_ERROR("Body load in local coordinate system not supported yet.")
        }
    }
};


/**
 * @brief A external flux term $S=(N)^T f$, where $f$ is functor evaluating the flux. 
 */
class NTf_Surface : public Term {
    protected:
        const NTfFunctor& f;
        int isurf;
    public:
    NTf_Surface (const Variable *testField, const NTfFunctor& f, int surf) ;

    /**
     * @brief Evaluates the linearization of term (the lhs contribution)
     * 
     * @param answer 
     * @param e 
     * @param coords 
     */
    void evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const override {}
    /**
     * @brief Evaluates Internal forces vector, i.e. $w^T(\grad N)^T f(p)$
     * 
     * @param cell 
     * @param coords 
     */
    void evaluate (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const override ;
    void getDimensions(Element& cell) const override {}
    void initializeCell(Element& cell) const override {}
    
};

/**
 * @brief A external flux term $S=(N)^T f$, where $f$ is functor evaluating the flux. 
 */
class NTf_Edge : public Term {
    protected:
        const NTfFunctor& f;
        int isurf;
    public:
   NTf_Edge (const Variable *testField, const NTfFunctor& f, int surf) ;

    /**
     * @brief Evaluates the linearization of term (the lhs contribution)
     * 
     * @param answer 
     * @param e 
     * @param coords 
     */
    void evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const override {}
    /**
     * @brief Evaluates Internal forces vector, i.e. $w^T(\grad N)^T f(p)$
     * 
     * @param cell 
     * @param coords 
     */
    void evaluate (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const override ;
    void getDimensions(Element& cell) const override {}
    void initializeCell(Element& cell) const override {}
    
};

/**
 * @brief A external flux term $S=(N)^T f$, where $f$ is functor evaluating the flux. 
 */
class NTf_Body : public Term {
    protected:
        const NTfFunctor& f;
    public:
    NTf_Body (const Variable *testField, const NTfFunctor& f) ;

    /**
     * @brief Evaluates the linearization of term (the lhs contribution)
     * 
     * @param answer 
     * @param e 
     * @param coords 
     */
    void evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const override {}
    /**
     * @brief Evaluates Internal forces vector, i.e. $w^T(\grad N)^T f(p)$
     * 
     * @param cell 
     * @param coords 
     */
    void evaluate (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const override ;
    void getDimensions(Element& cell) const override {}
    void initializeCell(Element& cell) const override {}
    
};

} // end namespace oofem
#endif // termlibrary_h
