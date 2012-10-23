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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef mixedgradientpressurecneumann_h
#define mixedgradientpressurecneumann_h

#include "mixedgradientpressurebc.h"
#include "boundary.h"
#include "dof.h"
#include "bctype.h"
#include "valuemodetype.h"
#include "classtype.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "dynalist.h"

namespace oofem {
class MasterDof;
class Node;
class IntegrationRule;
class SparseMtrx;
class SparseLinearSystemNM;

/**
 * Applies a mean deviatoric shear rate and pressure (Neumann boundary condition).
 * Introduces unknowns for the deviatoric stress, and this boundary condition adds contributions to both the left and right hand side of the system.
 * The approach expresses the deviatoric stresses/strain(rates) in a special (nonsymmetric) base;
 * @f[
 \boldsymbol A_{\mathrm{dev}} = \sum_{i=1}^{n} A_{\mathrm{dev},i} \boldsymbol E_i
 * @f]
 * where the bases @f$ \boldsymbol E_i @f$ can be chosen as 
 * @f[
 \boldsymbol E_1 = \frac{1}{\sqrt{2}}\begin{pmatrix}1 & 0 \\ 0 & -1\end{pmatrix}, \quad
 \boldsymbol E_2 = \begin{pmatrix}0 & 1 \\ 0 & 0\end{pmatrix}, \quad
 \boldsymbol E_3 = \begin{pmatrix}0 & 0 \\ 1 & 0\end{pmatrix}.
 * @f]
 * in 2D, and
 * @f[
 \boldsymbol E_1 = \frac{1}{\sqrt{6}}\begin{pmatrix}2 & 0 & 0\\ 0 & -1 & 0 \\ 0 & 0 & -1\end{pmatrix}, \quad
 \boldsymbol E_2 = \frac{1}{\sqrt{2}}\begin{pmatrix}0 & 0 & 0\\ 0 &  1 & 0 \\ 0 & 0 & -1\end{pmatrix}, \quad
 \boldsymbol E_3 = \begin{pmatrix}0 & 1 & 0\\ 0 & 0 & 0 \\ 0 & 0 & 0\end{pmatrix}, \ldots\quad
 \boldsymbol E_8 = \begin{pmatrix}0 & 0 & 0\\ 0 & 0 & 0 \\ 0 & 1 & 0\end{pmatrix}
 * @f]
 * in 3D (in 1D no deviatoric component exists).
 * 
 * Either bulk elements + side number can be used, or a boundary element (line in 2D, surface in 3D) can be used directly.
 * 
 * @note The 2D case assumes plane strain(rate), as is the case in _2dFlow.
 * @note This is only applicable to momentum balance equation. Both solid or fluids, incompressible or compressible, should work.
 * @note Should typically be prescribed on the entire external boundary of an representative volume element.
 * @note Should be applied to element boundaries, not DOFs.
 * @note The implementation doesn't assume that the stress is symmetric, so rigid body rotations are automatically removed.
 * @note Rigid body translations must be controlled separately.
 *
 * @see MixedGradientPressureNeumann
 * 
 * @author Mikael Öhman
 */
class MixedGradientPressureNeumann : public MixedGradientPressureBC
{
protected:
    /// Prescribed gradient @f$ d_{\mathrm{dev},ij} @f$ in Voigt form.
    FloatArray devGradient;
    /**
     * The volumetric part of what was sent in (needed to return the difference).
     * If caller takes care and sends in a deviatoric gradient, then this will be zero and the return value for the volumetric part will be the true volumetric change.
     */
    double volGradient;

    /// Prescribed pressure.
    double pressure;

    /// DOF-manager containing the unknown deviatoric stress.
    Node *sigmaDev;
    
    /// Element boundaries to integrate over. Boundary number 0 indicates that the element is a boundary element itself.
    dynaList< std::pair<int,int> > boundaries;

public:
    /**
     * Creates boundary condition with given number, belonging to given domain.
     * @param n Boundary condition number.
     * @param d Domain to which new object will belongs.
     */
    MixedGradientPressureNeumann(int n, Domain *d);

    /// Destructor
    virtual ~MixedGradientPressureNeumann();

    /**
     * Returns the number of internal DOF managers (=2).
     * This boundary condition stores its own DOF managers, one for @f$ d_{\mathrm{dev},ij} @f$ in which the DOFs are prescribed
     * and one for @f$ d_{\mathrm{vol}} @f$ for single free volumetric strain rate.
     */
    virtual int giveNumberOfInternalDofManagers();
    /**
     * Returns the volumetric DOF manager for i == 1, and the deviatoric manager for i == 2.
     */
    virtual DofManager *giveInternalDofManager(int i);

    /// Not relevant for this boundary condition.
    virtual bcType giveType() const { return UnknownBT; }

    /**
     * Initializes receiver according to object description stored in input record.
     * The input record contains two fields;
     * - devGradient \#columns { d_11 d_22 ... d_21 ... } (required)
     * - pressure p (required)
     * - elementSides List of element numbers and sides (interleaved) to apply boundary condition to.
     * - elements List of boundary elements to apply boundary condition to.
     * The gradient should be in Voigt notation (only the deviatoric part will be used)
     */
    virtual IRResultType initializeFrom(InputRecord *ir);
    
    virtual void addElementSide(int elem, int side);
    virtual void addElement(int elem);
    void clearElements();

    virtual int giveInputRecordString(std :: string &str, bool keyword = true);

    virtual void scale(double s);

    /**
     * Computes the homogenized fields through sensitivity analysis.
     * @param[out] stressDev Computes the homogenized deviatoric stress.
     * @param[out] vol Computes the homogenized volumetric gradient.
     * @param eid Equation ID that fields belong to.
     * @param tStep Time step for which field to obtain.
     */
    virtual void computeFields(FloatArray &sigmaDev, double &vol, EquationID eid, TimeStep *tStep);

    /**
     * Computes the macroscopic tangents through sensitivity analysis.
     * @param[out] Ed Tangent @f$ \frac{\partial \sigma_{\mathrm{dev}}}{\partial d_{\mathrm{dev}}} @f$.
     * @param[out] Ep Tangent @f$ \frac{\partial \sigma_{\mathrm{dev}}}{\partial p} @f$.
     * @param[out] Cd Tangent @f$ \frac{\partial d_{\mathrm{vol}}}{\partial d_{\mathrm{dev}}} @f$.
     * @param[out] Cp Tangent @f$ \frac{\partial d_{\mathrm{vol}}}{\partial p} @f$.
     * @param eid Equation ID for which the tangents belong.
     * @param tStep Time step for the tangents.
     */
    virtual void computeTangents(FloatMatrix &Ed, FloatArray &Ep, FloatArray &Cd, double &Cp, EquationID eid, TimeStep *tStep);

    /**
     * Set prescribed pressure.
     * @param p New prescribed pressure.
     */
    virtual void setPrescribedPressure(double p) { pressure = p; }

    /**
     * Sets the prescribed tensor from the matrix from given Voigt notation.
     * Assumes use of double values (gamma) for off-diagonal, usually the way for strain in Voigt form.
     * @param ddev Vector in Voigt format.
     */
    virtual void setPrescribedDeviatoricGradientFromVoigt(const FloatArray &ddev);

    virtual double assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
                                  CharType type, ValueModeType mode,
                                  const UnknownNumberingScheme &s, Domain *domain);
    
    virtual void assemble(SparseMtrx *answer, TimeStep *tStep, EquationID eid,
                          CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, Domain *domain);
    
    virtual void giveLocationArrays(AList<IntArray> &rows, AList<IntArray> &cols, EquationID eid, CharType type,
                                    const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, Domain *domain);

    virtual const char *giveClassName() const { return "MixedGradientPressureNeumann"; }
    virtual classType giveClassID() const { return MixedGradientPressureNeumannClass; }
    
protected:
    /// Helper function that creates suitable integration rule
    IntegrationRule *CreateIntegrationRule(Element *e, int order);
    
    /// Helper function that integrates the deviatoric tangent contribution from a single element boundary.
    void integrateDevTangent(FloatMatrix &answer, Element *e, int boundary);
    /// Helper function that integrates the volumetric tangent contribution from a single element boundary.
    void integrateVolTangent(FloatArray &answer, Element *e, int boundary);

    /// Converts from deviatoric to (normal) cartesian base (arrays are second order 2D tensors in Voigt notation)
    void fromDeviatoricBase2D(FloatArray &cartesian, FloatArray &deviatoric);
    /// Converts from deviatoric to (normal) cartesian base (arrays are second order 3D tensors in Voigt notation)
    void fromDeviatoricBase3D(FloatArray &cartesian, FloatArray &deviatoric);
    /// Converts from deviatoric to (normal) cartesian base (arrays are fourth order 2D tensors in Voigt notation)
    void fromDeviatoricBase2D(FloatMatrix &cartesian, FloatMatrix &deviatoric);
    /// Converts from deviatoric to (normal) cartesian base (arrays are fourth order 3D tensors in Voigt notation)
    void fromDeviatoricBase3D(FloatMatrix &cartesian, FloatMatrix &deviatoric);
};
} // end namespace oofem

#endif // mixedgradientpressurecneumann_h

