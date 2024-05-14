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
#ifndef termlibrary2_h
#define termlibrary2_h

#include "mpm.h"
#include "boundaryload.h"

// file containing various term definitions 

namespace oofem {

/**
 * @brief Evaluates $\delta B = B_{div}$ matrix; 
 * 
 * @param answer deltaB matrix
 * @param v 
 * @param interpol 
 * @param cell 
 * @param coords 
 */
void deltaB(FloatMatrix& answer, const Variable &v, const FEInterpolation& interpol, const Element& cell, const FloatArray& coords, const MaterialMode mmode) ;
/**
 * @brief Evaluates $B$ = matrix; 
 * 
 * @param answer B matrix
 * @param v 
 * @param interpol 
 * @param cell 
 * @param coords 
 */
void evalB(FloatMatrix& answer, const Variable &v, const FEInterpolation& interpol, const Element& cell, const FloatArray& coords, const MaterialMode mmode) ;

/** evaluetes volume fraction by interpolating the values in gp
 **/
  double evalVolumeFraction(const Variable&vf, MPElement& e, const FloatArray& coords, TimeStep* tstep);


  
/**
 * @brief Fluid mass balance term ($N_p u_{i,i} = N_p^T\delta\B$)
 * 
 */
class NTBdivTerm : public Term {
    protected:
    ValueModeType m;
    public:
    NTBdivTerm (const Variable &testField, const Variable& unknownField, ValueModeType m) ;

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
    
};


/// @brief evaluates ∫_Ω▒〖P 〖δv〗_(i,i) ϕdΩ〗
class deltaBTfiNpTerm : public Term {
    protected:
        const Variable& volumeFraction;
    public:
    deltaBTfiNpTerm (const Variable &testField, const Variable& unknownField, const Variable& volumeFraction) ;

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

};


/// @brief evalueates ∫_Ω▒〖N_d^T ∅_(,i) N_p  dΩ〖 r〗_p
class NdTdvfNpTerm : public Term {
    protected:
        const Variable& volumeFraction;
    public:
    /// @brief Constructor
    /// @param unknownField 
    /// @param testField 
    /// @param volumeFraction 
    NdTdvfNpTerm (const Variable &testField, const Variable& unknownField, const Variable& volumeFraction) ;

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

};

/// @brief evalueates ∫_Ω▒〖μ(v_(i,j)+v_(j,i) ) 〖δv〗_(i,j) dΩ〗
class BTmuBTerm : public Term {
    protected:
    public:
    /// @brief Constructor
    /// @param unknownField 
    /// @param testField 
    BTmuBTerm (const Variable &testField, const Variable& unknownField) ;

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

};

/// @brief evalueates ∫_Ω▒〖μ∅(〗 u ̇_(i,j)+u ̇_(j,i))〖δv〗_(i,j) dΩ
class BTmuVfBTerm : public Term {
    protected:
        const Variable& volumeFraction; 
    public:
    /// @brief Constructor
    /// @param unknownField 
    /// @param testField 
    /// @param volumeFraction 
    BTmuVfBTerm (const Variable &testField, const Variable& unknownField, const Variable& volumeFraction) ;

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

};

/// @brief evalueates ∫_Ω▒〖〖δv〗_i∅μS_ij^(-1) v_j dΩ+  〗
class NTmuVfSNTerm : public Term {
    protected:
        const Variable& volumeFraction; 
    public:
    /// @brief Constructor
    /// @param unknownField 
    /// @param testField 
    /// @param volumeFraction 
    NTmuVfSNTerm (const Variable &testField, const Variable& unknownField, const Variable& volumeFraction) ;

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

};

/// @brief evaluates ∫_Ω▒(δB)^T N_p  dΩ
class deltaBTNpTerm : public Term {
    protected:
    public:
    deltaBTNpTerm (const Variable &testField, const Variable& unknownField) ;

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

};

} // end namespace oofem
#endif // termlibrary_h
