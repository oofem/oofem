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

#ifndef mixedgradientpressuredirichlet_h
#define mixedgradientpressuredirichlet_h

#include "mixedgradientpressurebc.h"
#include "boundary.h"
#include "dof.h"
#include "bctype.h"
#include "valuemodetype.h"
#include "classtype.h"
#include "flotarry.h"
#include "flotmtrx.h"

namespace oofem {
class MasterDof;
class Node;
class SparseMtrx;
class SparseLinearSystemNM;

/**
 * Prescribes @f$ v_i = d_{ij}(x_j-\bar{x}_j) = (d_{\mathrm{dev},ij}+ \frac13 d_{\mathrm{vol}} I_{ij})(x_j+\bar{x}_j)} @f$
 * where @f$ d_{\mathrm{vol}} @f$ is unknown, loaded by a pressure.
 * This is mixed Dirichlet boundary condition in multiscale analysis where @f$ d = \partial_x s@f$
 * would a macroscopic gradient at the integration point, i.e. this is a boundary condition for prolongation.
 * It is also convenient to use when one wants to test a arbitrary specimen for shear.
 *
 * @note The 2D case assumes plane strain(rate), as is the case in _2dFlow.
 * @note It is possible to use a non-deviatoric strain rate, in which case, @f$ d_{\mathrm{vol}} @f$ obtains the value of the actual volumetric strain rate minus the volumetric part of @f$ d_{\mathrm{dev} @f$.
 * @note This is only applicable to momentum balance equation. Both solid or fluids, incompressible or compressible, should work.
 * @note Can only be prescribed to active dofs, so they must be denoted as such in the in the input file.
 *
 * @see MixedGradientPressureNeumann
 * 
 * @author Mikael Öhman
 */
class MixedGradientPressureDirichlet : public MixedGradientPressureBC
{
protected:
    /// Prescribed gradient @f$ d_{\mathrm{dev},ij} @f$ in Voigt form.
    FloatArray devGradient;

    /// Center coordinate @f$ \bar{x}_i @f$
    FloatArray centerCoord;

    /// Prescribed pressure.
    double pressure;

    /// DOF-manager containing the unknown volumetric strain(rate).
    Node *voldman;

    /// DOF-manager containing the known deviatoric strain(rate).
    Node *devdman;

    Dof *giveVolDof();

public:
    /**
     * Creates boundary condition with given number, belonging to given domain.
     * @param n Boundary condition number.
     * @param d Domain to which new object will belongs.
     */
    MixedGradientPressureDirichlet(int n, Domain *d);

    /// Destructor
    virtual ~MixedGradientPressureDirichlet();

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
     * - cCoords \#columns x_1 x_2 ... (optional, default 0)
     * The gradient should be deviatoric and in Voigt notation.
     * The prescribed tensor's columns must be equal to the size of the center coordinates.
     * The size of the center coordinates must be equal to the size of the coordinates in the applied nodes.
     */
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int giveInputRecordString(std :: string &str, bool keyword = true);

    virtual void scale(double s) { devGradient.times(s); pressure *= s; }

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

    /**
     * Set the center coordinate.
     * @param x New center coordinate.
     */
    virtual void setCenterCoordinate(const FloatArray &x) { centerCoord = x; }
    /// Returns the center coordinate
    virtual FloatArray &giveCenterCoordinate() { return centerCoord; }

    virtual double assembleVector(FloatArray &answer, TimeStep *tStep, EquationID eid,
                                  CharType type, ValueModeType mode,
                                  const UnknownNumberingScheme &s, Domain *domain);

    virtual bool requiresActiveDofs() { return true; }
    virtual bool isPrimaryDof(ActiveDof *dof);

    virtual double giveBcValue(ActiveDof *dof, ValueModeType mode, TimeStep *tStep);
    virtual bool hasBc(ActiveDof *dof, TimeStep *tStep);

    /// Returns true is DOF represents one of the deviatoric parts.
    bool isDevDof(ActiveDof *dof);

    virtual int giveNumberOfMasterDofs(ActiveDof *dof);
    virtual Dof *giveMasterDof(ActiveDof *dof, int mdof);
    virtual void computeDofTransformation(ActiveDof *dof, FloatArray &masterContribs);

    double giveUnknown(double vol, const FloatArray &dev, ValueModeType mode, TimeStep *tStep, ActiveDof *dof);
    virtual double giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *tStep, ActiveDof *dof);
    virtual double giveUnknown(EquationID type, ValueModeType mode, TimeStep *tStep, ActiveDof *dof);

    virtual const char *giveClassName() const { return "MixedGradientPressureDirichlet"; }
    virtual classType giveClassID() const { return MixedGradientPressureDirichletClass; }
};
} // end namespace oofem

#endif // mixedgradientpressuredirichlet_h

