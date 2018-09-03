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

#ifndef sprnodalrecoverymodel_h
#define sprnodalrecoverymodel_h

#include "nodalrecoverymodel.h"
#include "interface.h"

#define _IFT_SPRNodalRecoveryModel_Name "spr"

namespace oofem {
class GaussPoint;
class SPRNodalRecoveryModelInterface;
class ProcessCommunicator;

enum SPRPatchType {
    SPRPatchType_none = 0,
    SPRPatchType_2dxy = 1,
    SPRPatchType_3dBiLin,
    SPRPatchType_2dquadratic,
    SPRPatchType_3dBiQuadratic
};

/**
 * The Superconvergent Patch Recovery (SPR) nodal recovery model is based on paper of Zienkiewicz and Zhu
 * "The Superconvergent Patch recovery and Posteriori Error Estimates. Part 1: The Recovery Technique",
 * Int. Journal for Num. Meth in Engng, vol. 33, 1331-1364, 1992.
 * The recovery uses local discrete least square smoothing over an element patch surrounding the particular
 * node considered.
 */
class OOFEM_EXPORT SPRNodalRecoveryModel : public NodalRecoveryModel
{
protected:
    /**
     * Helper structure to pass required arguments to packing/unpacking functions
     * needed in parallel mode
     */
    struct parallelStruct {
        FloatArray *dofManValues;
        IntArray *dofManPatchCount;
        IntArray *regionNodalNumbers;
        int regionValSize;
        parallelStruct(FloatArray *a, IntArray *b, IntArray *c, int d) :
            dofManValues(a), dofManPatchCount(b), regionNodalNumbers(c), regionValSize(d) { }
    };

public:
    /// Constructor.
    SPRNodalRecoveryModel(Domain * d);
    /// Destructor.
    virtual ~SPRNodalRecoveryModel();

    int recoverValues(Set elementSet, InternalStateType type, TimeStep *tStep) override;

    const char *giveClassName() const override { return "SPRNodalRecoveryModel"; }

private:
    /**
     * Initializes the region table indicating regions to skip.
     * @param regionMap Region table, the nonzero entry for region indicates region to skip due to
     * unsupported elements or incompatible value size
     * @param regionTypes SPRPatchType of each region.
     * @param type Determines the type of internal variable to be recovered
     */
    void initRegionMap(IntArray &regionMap, IntArray &regionTypes, InternalStateType type);

    void determinePatchAssemblyPoints(IntArray &pap, SPRPatchType regType, Set &elemset);
    void initPatch(IntArray &patchElems, IntArray &dofManToDetermine, IntArray &pap, int papNumber, Set &elementList);
    void computePatch(FloatMatrix &a, IntArray &patchElems, int &regionValSize,
                      SPRPatchType regType, InternalStateType type, TimeStep *tStep);
    void determineValuesFromPatch(FloatArray &dofManValues, IntArray &dofManCount,
                                  IntArray &regionNodalNumbers, IntArray &dofManToDetermine,
                                  FloatMatrix &a, SPRPatchType type);
    void computePolynomialTerms(FloatArray &P, const FloatArray &coords, SPRPatchType type);
    int  giveNumberOfUnknownPolynomialCoefficients(SPRPatchType regType);
    SPRPatchType determinePatchType(Set &elementList);

#ifdef __PARALLEL_MODE
    void initCommMaps();
    void exchangeDofManValues(FloatArray &dofManValues,
                              IntArray &dofManPatchCount, IntArray &regionNodalNumbers,
                              int regionValSize);
    int packSharedDofManData(parallelStruct *s, ProcessCommunicator &processComm);
    int unpackSharedDofManData(parallelStruct *s, ProcessCommunicator &processComm);
#endif
};

/**
 * The element interface required by ZZNodalRecoveryModel.
 */
class OOFEM_EXPORT SPRNodalRecoveryModelInterface : public Interface
{
public:
    /// Constructor.
    SPRNodalRecoveryModelInterface() { }

    /// @name The element interface required by SPRNodalRecoveryModelInterface
    //@{
    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap) = 0;
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap) = 0;
    virtual int SPRNodalRecoveryMI_giveNumberOfIP() = 0;
    virtual SPRPatchType SPRNodalRecoveryMI_givePatchType() = 0;
    //@}
};
} // end namespace oofem
#endif // sprnodalrecoverymodel_h
