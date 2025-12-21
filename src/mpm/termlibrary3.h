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
#ifndef termlibrary3_h
#define termlibrary3_h

#include "mpm.h"
#include "boundaryload.h"
#include "termlibrary.h"

// file containing various term definitions 

namespace oofem {

 /**
 * @brief A Linear momentum balance equation term ($B^T\sigma(u)$)
 * 
 */
class TMBTSigTerm : public BTSigTerm {
    protected:
    const Variable* temperatureField;
    public:
    TMBTSigTerm (const Variable *testField, const Variable* unknownField, const Variable* temperatureField) ;
    /**
     * @brief Evaluates Internal forces vector, i.e. $b^T\sigma(u)$
     * 
     * @param cell 
     * @param coords 
     */
    void evaluate (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const override;    
    void computeTMgeneralizedStrain (FloatArray& answer, FloatMatrix&B, MPElement& cell, const FloatArray& lcoords, MaterialMode mmode, TimeStep* tstep) const;
};
/**
 * @brief momentum balance term (lhs only) (\int (\partial N)^T D (-\alpha\Pi))
 * 
 */
class BDalphaPiTerm : public Term {
    protected:
    ValueModeType m;
    public:
    BDalphaPiTerm (const Variable *testField, const Variable* unknownField, ValueModeType m) ;

    /**
     * @brief Evaluates the linearization of $B^T\sigma(u)$, i.e. $B^TDBu$
     * 
     * @param answer 
     * @param e 
     * @param coords 
     */
    void evaluate_lin (FloatMatrix& answer, MPElement& e, GaussPoint* gp, TimeStep* tstep) const override;
    /**
     * @brief Empty, this should be evaluated by BTSigTerm term$
     * 
     * @param cell 
     * @param coords 
     */
    void evaluate (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const override;
    void getDimensions(Element& cell) const override;
    void initializeCell(Element& cell) const override;

    protected:
    
};

class TMgNTfTerm : public gNTfTerm {
    protected:
    public:
    TMgNTfTerm (const Variable *testField, const Variable* unknownField, MatResponseMode lhsType, MatResponseMode rhsType) ;

    
    /**
     * @brief Evaluates Internal forces vector, i.e. $w^T(\grad N)^T f(p)$
     * 
     * @param cell 
     * @param coords 
     */
    void evaluate (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const  override;
    
};


/// @brief evaluates lhs of ∫(𝐵𝑢)𝑇(𝑑𝜎(𝜀,𝑇)/𝑑𝑇)
class BTdSigmadT : public Term {
    protected:
    public:
    BTdSigmadT (const Variable *testField, const Variable* unknownField) ;

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


/// @brief evaluates flux due to convection BC ∫(Nt)a(T-Te)ds, Te being external temperature and a being convective heat transfer coefficient 
class NTaTmTe : public Term {
    protected:
        BoundaryLoad* bl;
        int boundaryID;
        char boundaryType;
    public:
    NTaTmTe (const Variable *testField, const Variable* unknownField, BoundaryLoad* bl, int boundaryID, char boundaryType) ;

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


class InternalTMFluxSourceTerm : public TMBTSigTerm {
    protected:
    public:
    InternalTMFluxSourceTerm (const Variable *testField, const Variable* unknownField, const Variable* temperatureField) ;
    /**
     * @brief Evaluates Internal forces vector, i.e. $b^T\sigma(u)$
     * 
     * @param cell 
     * @param coords 
     */
    void evaluate (FloatArray&, MPElement& cell, GaussPoint* gp, TimeStep* tstep) const override;    
};



} // end namespace oofem
#endif // termlibrary_h
