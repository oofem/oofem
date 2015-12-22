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

#ifndef mixedgradientpressuredirichlet_h
#define mixedgradientpressuredirichlet_h

#include "mixedgradientpressurebc.h"
#include "boundarycondition.h"
#include "dof.h"
#include "bctype.h"
#include "valuemodetype.h"
#include "floatarray.h"
#include "floatmatrix.h"

#include <memory>

#define _IFT_MixedGradientPressureDirichlet_Name "mixedgradientpressuredirichlet"

namespace oofem {
class MasterDof;
class Node;
class SparseMtrx;
class SparseLinearSystemNM;

/**
 * Prescribes @f$ v_i = d_{ij}(x_j-\bar{x}_j) = (d_{\mathrm{dev},ij}+ \frac13 d_{\mathrm{vol}} I_{ij})(x_j+\bar{x}_j)) @f$
 * where @f$ d_{\mathrm{vol}} @f$ is unknown, loaded by a pressure.
 * This is mixed Dirichlet boundary condition in multiscale analysis where @f$ d = \partial_x s@f$
 * would a macroscopic gradient at the integration point, i.e. this is a boundary condition for prolongation.
 * It is also convenient to use when one wants to test a arbitrary specimen for shear.
 *
 * @note The 2D case assumes plane strain(rate), as is the case in _2dFlow.
 * @note It is possible to use a non-deviatoric strain rate, in which case, @f$ d_{\mathrm{vol}} @f$ obtains the value of the actual volumetric strain rate minus the volumetric part of @f$ d_{\mathrm{dev}} @f$.
 * @note This is only applicable to momentum balance equation. Both solid or fluids, incompressible or compressible, should work.
 * @note Can only be prescribed to active dofs, so they must be denoted as such in the in the input file.
 *
 * @see MixedGradientPressureNeumann
 *
 * @author Mikael Öhman
 */
class OOFEM_EXPORT MixedGradientPressureDirichlet : public MixedGradientPressureBC
{
protected:
    /// Prescribed gradient @f$ d_{\mathrm{dev},ij} @f$ in Voigt form.
    FloatArray devGradient;

    /// Center coordinate @f$ \bar{x}_i @f$
    FloatArray centerCoord;

    /// Prescribed pressure.
    double pressure;

    /// DOF-manager containing the unknown volumetric strain(rate).
    std :: unique_ptr< Node > voldman;
    int vol_id;

    /// DOF-manager containing the known deviatoric strain(rate).
    std :: unique_ptr< Node > devdman;
    IntArray dev_id;

    Dof *giveVolDof();

public:
    /**
     * Creates boundary condition with given number, belonging to given domain.
     * @param n Boundary condition number.
     * @param d Domain to which new object will belongs.
     */
    MixedGradientPressureDirichlet(int n, Domain * d);

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
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual void scale(double s) {
        devGradient.times(s);
        pressure *= s;
    }

    virtual void computeFields(FloatArray &stressDev, double &vol, TimeStep *tStep);
    virtual void computeTangents(FloatMatrix &Ed, FloatArray &Ep, FloatArray &Cd, double &Cp, TimeStep *tStep);

    virtual void setPrescribedPressure(double p) { pressure = p; }
    virtual void setPrescribedDeviatoricGradientFromVoigt(const FloatArray &ddev);

    /**
     * Set the center coordinate.
     * @param x New center coordinate.
     */
    virtual void setCenterCoordinate(const FloatArray &x) { centerCoord = x; }
    /// Returns the center coordinate
    virtual FloatArray &giveCenterCoordinate() { return centerCoord; }

    virtual void assembleVector(FloatArray &answer, TimeStep *tStep,
                                CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, FloatArray *eNorms = NULL);

    virtual bool requiresActiveDofs() { return true; }
    virtual bool isPrimaryDof(ActiveDof *dof);

    virtual double giveBcValue(Dof *dof, ValueModeType mode, TimeStep *tStep);
    virtual bool hasBc(Dof *dof, TimeStep *tStep);

    /// Returns true is DOF represents one of the deviatoric parts.
    bool isDevDof(Dof *dof);

    virtual int giveNumberOfMasterDofs(ActiveDof *dof);
    virtual Dof *giveMasterDof(ActiveDof *dof, int mdof);
    virtual void computeDofTransformation(ActiveDof *dof, FloatArray &masterContribs);

    double giveUnknown(double vol, const FloatArray &dev, ValueModeType mode, TimeStep *tStep, ActiveDof *dof);
    virtual double giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *tStep, ActiveDof *dof);
    virtual double giveUnknown(ValueModeType mode, TimeStep *tStep, ActiveDof *dof);

    virtual const char *giveClassName() const { return "MixedGradientPressureDirichlet"; }
    virtual const char *giveInputRecordName() const { return _IFT_MixedGradientPressureDirichlet_Name; }
};
} // end namespace oofem

#endif // mixedgradientpressuredirichlet_h
