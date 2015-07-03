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

#ifndef PRESCRIBEDGRADIENTBCNEUMANN_H_
#define PRESCRIBEDGRADIENTBCNEUMANN_H_

#include "prescribedgradienthomogenization.h"
#include "activebc.h"

#include <memory>

#define _IFT_PrescribedGradientBCNeumann_Name   "prescribedgradientbcneumann"

namespace oofem {
class Node;
class Element;
/**
 * Imposes a prescribed gradient weakly on the boundary
 * with a Neumann boundary condition.
 *
 * @author Erik Svenning
 * @author Mikael Öhman
 */
class OOFEM_EXPORT PrescribedGradientBCNeumann : public ActiveBoundaryCondition, public PrescribedGradientHomogenization
{
public:
    PrescribedGradientBCNeumann(int n, Domain *d);
    virtual ~PrescribedGradientBCNeumann();

    virtual int giveNumberOfInternalDofManagers() { return 1; }
    virtual DofManager *giveInternalDofManager(int i);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual bcType giveType() const { return UnknownBT; }

    virtual void scale(double s);

    virtual void assembleVector(FloatArray &answer, TimeStep *tStep,
                                CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, FloatArray *eNorm = NULL);

    virtual void assemble(SparseMtrx &answer, TimeStep *tStep,
                          CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);

    virtual void giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                                    const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s);

    virtual const char *giveClassName() const { return "PrescribedGradientBCNeumann"; }
    virtual const char *giveInputRecordName() const { return _IFT_PrescribedGradientBCNeumann_Name; }

    virtual void computeField(FloatArray &sigma, TimeStep *tStep);
    virtual void computeTangent(FloatMatrix &tangent, TimeStep *tStep);

    void giveStressLocationArray(IntArray &oCols, const UnknownNumberingScheme &r_s);

protected:
    /// DOF-manager containing the unknown homogenized stress.
    std :: unique_ptr< Node > mpSigmaHom;
    IntArray mSigmaIds;

    /// Help function that integrates the tangent contribution from a single element boundary.
    void integrateTangent(FloatMatrix &oTangent, Element *e, int iBndIndex);
};
} /* namespace oofem */

#endif /* PRESCRIBEDGRADIENTBCNEUMANN_H_ */
