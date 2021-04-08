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
 *               Copyright (C) 1993 - 2021   Borek Patzak
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

#ifndef PRESCRIBEDGRADIENTBCNEUMANNRC_H_
#define PRESCRIBEDGRADIENTBCNEUMANNRC_H_

#include "prescribedgradientbcneumann.h"

#define _IFT_PrescribedGradientBCNeumannRC_Name   "prescribedgradientbcneumannrc"

namespace oofem {
class Node;
class Element;
/**
 * Prescribes a displacement gradient with Neumann boundary condition.
 * Used in multiscale analyses of reinforced concrete structures.
 * In the RVE, applied to element boundary of the concrete solid.
 * Currently, works only in 2D plane stress mode.
 *
 *
 * @author Adam Sciegaj
 */
class OOFEM_EXPORT PrescribedGradientBCNeumannRC : public PrescribedGradientBCNeumann
{
public:
    PrescribedGradientBCNeumannRC(int n, Domain *d);
    virtual ~PrescribedGradientBCNeumannRC();

    void computeField(FloatArray &sigma, TimeStep *tStep) override;

    void assembleVector(FloatArray &answer, TimeStep *tStep,
                        CharType type, ValueModeType mode,
                        const UnknownNumberingScheme &s, FloatArray *eNorm=nullptr, void* lock=nullptr) override;

    void assemble(SparseMtrx &answer, TimeStep *tStep,
                  CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s,
                  double scale = 1.0, void* lock=nullptr) override;

protected:
    virtual double domainSize(Domain *d, int set) override;
};
} /* namespace oofem */

#endif /* PRESCRIBEDGRADIENTBCNEUMANNRC_H_ */
