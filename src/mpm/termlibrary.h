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
#ifndef termlibrary_h
#define termlibrary_h

#include "mpm.h"

namespace oofem {
/**
 * @brief A Linear momentum balance equation term ($B^T\sigma(u)$)
 * 
 */
class BTSigTerm : public Term {
    protected:
    public:
    BTSigTerm (const Variable& unknownField, const Variable &testField) ;

    /**
     * @brief Evaluates the linearization of $B^T\sigma(u)$, i.e. $B^TDBu$
     * 
     * @param answer 
     * @param e 
     * @param coords 
     */
    void evaluate_dw (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const override;
    /**
     * @brief Evaluates Internal forces vector, i.e. $b^T\sigma(u)$
     * 
     * @param cell 
     * @param coords 
     */
    void evaluate_c (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const override;
    void getDimensions_dw(Element& cell) const override;
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
    void grad(FloatMatrix& answer, const Variable &v, const FEInterpolation& interpol, const Element& cell, const FloatArray& coords) const  ;
    
};

/**
 * @brief A continuity equation term ($Hp=(\grad N_p)^T f(p)$, where $f=\bf{k}/\mu \grad p$. 
 */
class gNTfTerm : public Term {
    protected:
    public:
    gNTfTerm (const Variable& unknownField, const Variable &testField) ;

    /**
     * @brief Evaluates $\bf{H}$ matrix, the linearization of $w^T(\grad N)^T f(p)$, i.e. $(\grad N)^T \bf{k}/\mu \grad p = \bf{H}$
     * 
     * @param answer 
     * @param e 
     * @param coords 
     */
    void evaluate_dw (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const override;
    /**
     * @brief Evaluates Internal forces vector, i.e. $w^T(\grad N)^T f(p)$
     * 
     * @param cell 
     * @param coords 
     */
    void evaluate_c (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const  override;
    void getDimensions_dw(Element& cell) const  override;
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
    void grad(FloatMatrix& answer, const Variable &v, const FEInterpolation& interpol, const Element& cell, const FloatArray& coords) const ;
    
};

/**
 * @brief A continuity equation term $Qp=(B)^T \alpha\bf{m}N_p$. 
 */
class BTamNTerm : public Term {
    protected:
    public:
    BTamNTerm (const Variable& unknownField, const Variable &testField) ;

    /**
     * @brief Evaluates the linearization of receiver, i.e. the LHS term
     * 
     * @param answer 
     * @param e 
     * @param coords 
     */
    void evaluate_dw (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const override;
    /**
     * @brief Evaluates Internal forces vector, i.e. $w^T(\grad N)^T f(p)$
     * 
     * @param cell 
     * @param coords 
     */
    void evaluate_c (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const override;
    void getDimensions_dw(Element& cell) const override;
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
    void grad(FloatMatrix& answer, const Variable &v, const FEInterpolation& interpol, const Element& cell, const FloatArray& coords) const ;
    
};

/**
 * @brief A continuity equation term $Q^T(du\over dt)=(N)^T \alpha\bf{m}^TB du\over dt$. 
 */
class NTamTBTerm : public Term {
    protected:
    public:
    NTamTBTerm (const Variable& unknownField, const Variable &testField) ;

    /**
     * @brief Evaluates the linearization of receiver, i.e. the LHS term
     * 
     * @param answer 
     * @param e 
     * @param coords 
     */
    void evaluate_dw (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const override;
    /**
     * @brief Evaluates Internal forces vector, i.e. $w^T(\grad N)^T f(p)$
     * 
     * @param cell 
     * @param coords 
     */
    void evaluate_c (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const override;
    void getDimensions_dw(Element& cell) const override;
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
    void grad(FloatMatrix& answer, const Variable &v, const FEInterpolation& interpol, const Element& cell, const FloatArray& coords) const ;
    
};


/**
 * @brief A continuity equation compressibility matrix $S=(N_p)^T c\ N_p$, where $c=({\alpha-n}\over{K_s}+{n}\over{K_w})$. 
 */
class NTcN : public Term {
    protected:
    public:
    NTcN (const Variable& unknownField, const Variable &testField) ;

    /**
     * @brief Evaluates the linearization of term (the lhs contribution)
     * 
     * @param answer 
     * @param e 
     * @param coords 
     */
    void evaluate_dw (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const override;
    /**
     * @brief Evaluates Internal forces vector, i.e. $w^T(\grad N)^T f(p)$
     * 
     * @param cell 
     * @param coords 
     */
    void evaluate_c (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const override;
    void getDimensions_dw(Element& cell) const override;
    void initializeCell(Element& cell) const override;
    
};



} // end namespace oofem
#endif // termlibrary_h