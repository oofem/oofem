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

#ifndef zznodalrecoverymodel_h
#define zznodalrecoverymodel_h

#include "compiler.h"

#include "nodalrecoverymodel.h"
#include "interface.h"

namespace oofem {
class GaussPoint;
class ZZNodalRecoveryModelInterface;

/**
 * The nodal recovery model based on paper of Zienkiewicz and Zhu "A Simple Estimator and Adaptive
 * Procedure for Practical Engineering Analysis". The recovery is based
 * on nodal averaging or projection process in which it is assumed that the stress @f$ \sigma_\star @f$ is
 * interpolated by the same function as the displacement.
 */
class ZZNodalRecoveryModel : public NodalRecoveryModel
{
protected:
    /**
     * Helper structure to pass required arguments to packing/unpacking functions
     * needed in parallel mode.
     */
    struct parallelStruct {
        FloatArray *lhs;
        FloatMatrix *rhs;
        IntArray *regionNodalNumbers;
        parallelStruct(FloatArray *a, FloatMatrix *b, IntArray *c):
            lhs(a), rhs(b), regionNodalNumbers(c) { }
    };

public:
    /// Constructor.
    ZZNodalRecoveryModel(Domain *d);
    /// Destructor.
    virtual ~ZZNodalRecoveryModel();

    virtual int recoverValues(InternalStateType type, TimeStep *tStep);

private:
    /**
     * Initializes the region table indicating regions to skip.
     * @param regionMap Region table, the nonzero entry for region indicates region to skip due to
     * unsupported elements or incompatible value size.
     * @param regionValSize Contains the record size for each region.
     * @param type Determines the type of internal variable to be recovered.
     */
    void initRegionMap(IntArray &regionMap, IntArray &regionValSize, InternalStateType type);

#ifdef __PARALLEL_MODE
    void initCommMaps();
    void exchangeDofManValues(int ireg, FloatArray &lhs, FloatMatrix &rhs, IntArray &rn);
    int packSharedDofManData(parallelStruct *s, ProcessCommunicator &processComm);
    int unpackSharedDofManData(parallelStruct *s, ProcessCommunicator &processComm);
#endif
};

/**
 * The element interface required by ZZNodalRecoveryModel.
 */
class ZZNodalRecoveryModelInterface : public Interface
{
public:
    /// Constructor
    ZZNodalRecoveryModelInterface() { }

    /// @name The element interface required by ZZNodalRecoveryModel
    //@{
    /**
     * Computes the element contribution to @f$ \int_\Omega N^{\mathrm{T}}\alpha\;\mathrm{d}\Omega @f$,
     * where @f$ \alpha @f$ is quantity to be recovered (for example stress or strain vector).
     * The size of answer should be recordSize*numberofDofManagers.
     * @param answer Contains the result.
     * @param type Determines the type of internal variable to be recovered.
     * @param tStep Time step.
     */
    virtual void ZZNodalRecoveryMI_computeNValProduct(FloatMatrix &answer, InternalStateType type, TimeStep *tStep);
    /**
     * Computes the element contribution to @f$\int_\Omega N^{\mathrm{T}} \cdot N\;\mathrm{d}\Omega @f$ term.
     * The size of answer should be [recordSize*numberofDofManagers].
     * @param answer Contain diagonalized result.
     * @param type Determines the type of internal variable to be recovered.
     */
    virtual void ZZNodalRecoveryMI_computeNNMatrix(FloatArray &answer, InternalStateType type);
    /**
     * Returns the size of DofManger record required to hold recovered values for given mode.
     * Default implementation uses element giveIPValueSize method.
     * @param type Determines the type of internal variable to be recovered.
     * @return Size of DofManger record required to hold recovered values.
     */
    virtual int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);
    /**
     * Returns the corresponding element to interface.
     */
    virtual Element *ZZNodalRecoveryMI_giveElement() = 0;
    /**
     * Evaluates N matrix (interpolation estimated value matrix).
     * Default implementation requires element to provide valid interpolation via giveInterpolation method.
     */
    virtual void ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatArray &answer, GaussPoint *aGaussPoint,
                                                                     InternalStateType type);
    //@}
};
} // end namespace oofem
#endif // zznodalrecoverymodel_h
