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

#ifndef PRESCRIBEDDISPSLIPBCNEUMANN_H_
#define PRESCRIBEDDISPSLIPBCNEUMANN_H_

#include "prescribeddispsliphomogenization.h"
#include "activebc.h"
#include "dynamicinputrecord.h"

#include <memory>

#define _IFT_PrescribedDispSlipBCNeumannRC_Name   "prescribeddispslipbcneumannrc"
#define _IFT_PrescribedDispSlipBCNeumannRC_ConcreteVolSet "concretevolset"
#define _IFT_PrescribedDispSlipBCNeumannRC_RebarSets "rebarsets"

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
 * For the optional slip and slip gradient, input fields concretevolset and rebarsets must be specified. concretevolset is the element set containing
 * all concrete elements. rebarsets is an array of sets, each containing all reinforcement bars in one direction.
 *
 *
 * @author Adam Sciegaj
 *
 */
class OOFEM_EXPORT PrescribedDispSlipBCNeumannRC : public ActiveBoundaryCondition, public PrescribedDispSlipHomogenization
{
public:
    PrescribedDispSlipBCNeumannRC(int n, Domain *d);
    virtual ~PrescribedDispSlipBCNeumannRC();

    int giveNumberOfInternalDofManagers() override;
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
    void computeTransferStress(FloatArray &bStress, TimeStep *tStep) override;
    void computeReinfStress(FloatArray &rStress, TimeStep *tStep) override;

    void computeTangent(FloatMatrix &tangent, TimeStep *tStep) override;

protected:
    bool dispGradON = false; /// on/off flag specifying whether the displacement gradient should be applied via Neumann BCs
    bool slipON = false; /// on/off flag specifying whether the slip field should be applied via Neumann BCs
    bool slipGradON = false; /// on/off flag specifying whether the slip gradient should be applied via Neumann BCs

    /// DofManager for effective stress
    std :: unique_ptr< Node > mpSigmaHom;
    IntArray mSigmaIds;

    /// DofManager for effective transfer stress
    std :: unique_ptr< Node > lmTauHom;
    IntArray lmTauIds;

    /// DofManager for effective reinforcement membrane stress
    std :: unique_ptr< Node > lmSigmaSHom;
    IntArray lmSigmaSIds;

    /// Element set containing concrete elements (in the volume of RVE)
    int concreteVolSet=0;
    /// IntArray containing sets of individual reinforcement bars
    IntArray rebarSets;

    void integrateTangentStress(FloatMatrix &oTangent, Element *e, int iBndIndex);
    void assembleVectorStress(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode, const UnknownNumberingScheme &s, FloatArray *eNorm);
    void assembleOnStress(SparseMtrx &answer, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale);

    void integrateTangentBStressSteel(FloatMatrix &oTangent, Element *e, const int &rebarSet);
    void integrateTangentBStressConcrete(FloatMatrix &oTangent, Element *e);
    void assembleVectorBStress(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode, const UnknownNumberingScheme &s);
    void assembleOnTransferStress(SparseMtrx &answer, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale);

    void integrateTangentRStressSteel(FloatMatrix &oTangent, Element *e, const int &rebarSet);
    void integrateTangentRStressConcrete(FloatMatrix &oTangent, Element *e, int iBndIndex);
    void assembleVectorRStress(FloatArray &answer, TimeStep *tStep, CharType type, ValueModeType mode, const UnknownNumberingScheme &s);
    void assembleOnReinfStress(SparseMtrx &answer, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s, double scale);

    void computeWeightMatrix(FloatMatrix &C, const IntArray &reinfSets);
    void computeRebarDyad(FloatMatrix &dyad, const int &reinfSet);
    double computeInterfaceLength(const IntArray &reinfSets);

    virtual double domainSize(Domain *d, int set) override;
};
} /* namespace oofem */

#endif /* PRESCRIBEDDISPSLIPBCNEUMANN_H_ */
