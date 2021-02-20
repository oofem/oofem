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

#ifndef PRESCRIBEDDISPSLIPBCNEUMANN_H_
#define PRESCRIBEDDISPSLIPBCNEUMANN_H_

#include "prescribeddispsliphomogenization.h"
#include "activebc.h"

#include <memory>

#define _IFT_PrescribedDispSlipBCNeumannRC_Name   "prescribeddispslipbcneumannrc"

namespace oofem {
class Node;
class Element;
/**
 * Prescribes macroscopic displacement gradient, reinforcement slip field and slip gradient
 * on the reinforced concrete RVE using Neumann boundary conditions, cf.
 * Sciegaj, A., Larsson, F., Lundgren, K., Nilenius, F., & Runesson, K. (2019). A multiscale model for reinforced concrete with macroscopic variation of reinforcement slip. Computational Mechanics, 63(2), 139–158. https://doi.org/10.1007/s00466-018-1588-3
 * Macroscopic slip and slip gradient fields are optional (see PrescribedDispSlipHomogenization).
 * If not specified, this BC will prescribe only the displacement gradient, i.e., it will work the same was as a PrescribedGradient bc.
 * Works with 2D RVEs comprising solid elements (concrete), reinforcement (beam/truss elements) and interface elements in between.
 * Used in multiscale analyses of reinforced concrete structures. Currently, only orthogonal reinforcement in X and Y direction is supported.
 *
 * @author Adam Sciegaj
 *
 */
class OOFEM_EXPORT PrescribedDispSlipBCNeumannRC : public ActiveBoundaryCondition, public PrescribedDispSlipHomogenization
{
public:
    PrescribedDispSlipBCNeumannRC(int n, Domain *d);
    virtual ~PrescribedDispSlipBCNeumannRC();

    int giveNumberOfInternalDofManagers() override { return 1; }
    DofManager *giveInternalDofManager(int i) override;

    void initializeFrom(InputRecord &ir) override;
    void giveInputRecord(DynamicInputRecord &input) override;

    bcType giveType() const override { return UnknownBT; }

    void scale(double s) override;

    void assembleVector(FloatArray &answer, TimeStep *tStep,
                        CharType type, ValueModeType mode,
                        const UnknownNumberingScheme &s, FloatArray *eNorm=nullptr, void* lock=nullptr) override;

    void assemble(SparseMtrx &answer, TimeStep *tStep,
                  CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, 
                  double scale = 1.0, void* lock=nullptr) override;

    void giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type,
                            const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s) override;

    const char *giveClassName() const override { return "PrescribedDispSlipBCNeumannRC"; }
    const char *giveInputRecordName() const override { return _IFT_PrescribedDispSlipBCNeumannRC_Name; }

    void computeStress(FloatArray &sigma, TimeStep *tStep) override;
    void computeTransferStress(FloatArray &bStress, TimeStep *tStep) override {};
    void computeReinfStress(FloatArray &rStress, TimeStep *tStep) override {};

    void computeTangent(FloatMatrix &tangent, TimeStep *tStep) override;

    void giveStressLocationArray(IntArray &oCols, const UnknownNumberingScheme &r_s);

protected:
    std :: unique_ptr< Node > mpSigmaHom;
    IntArray mSigmaIds;

    void integrateTangent(FloatMatrix &oTangent, Element *e, int iBndIndex);
    virtual double domainSize(Domain *d, int set) override;
};
} /* namespace oofem */

#endif /* PRESCRIBEDDISPSLIPBCNEUMANN_H_ */
