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

#ifndef assemblercallback_h
#define assemblercallback_h

#include "oofemenv.h"
#include "valuemodetype.h" ///@todo We shouldn't have this for assembling vectors or matrices(!) / Mikael
#include "matresponsemode.h"
#include "chartype.h"

namespace oofem {
class IntArray;
class FloatArray;
class FloatMatrix;
class Element;
class DofManager;
class TimeStep;
class NodalLoad;
class BodyLoad;
class BoundaryLoad;
class SurfaceLoad;
class EdgeLoad;
class UnknownNumberingScheme;
class SparseMtrx;
class ActiveBoundaryCondition;

/**
 * Callback class for assembling specific types of vectors.
 * Default implementations are that no contributions are considered (empty vectors on output).
 * @author Mikael Öhman
 */
class OOFEM_EXPORT VectorAssembler
{
public:
    virtual void vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const;
    virtual void vectorFromLoad(FloatArray &vec, Element &element, BodyLoad *load, TimeStep *tStep, ValueModeType mode) const;
    virtual void vectorFromSurfaceLoad(FloatArray &vec, Element &element, SurfaceLoad *load, int boundary, TimeStep *tStep, ValueModeType mode) const;
    virtual void vectorFromEdgeLoad(FloatArray &vec, Element &element, EdgeLoad *load, int edge, TimeStep *tStep, ValueModeType mode) const;
    virtual void vectorFromNodeLoad(FloatArray &vec, DofManager &dman, NodalLoad *load, TimeStep *tStep, ValueModeType mode) const;
    virtual void assembleFromActiveBC(FloatArray &answer, ActiveBoundaryCondition &bc, TimeStep* tStep, ValueModeType mode, const UnknownNumberingScheme &s, FloatArray *eNorms, void* lock=nullptr) const;

    /// Default implementation takes all the DOF IDs
    virtual void locationFromElement(IntArray &loc, Element &element, const UnknownNumberingScheme &s, IntArray *dofIds = nullptr) const;
    /// Default implementation takes all the DOF IDs
    virtual void locationFromElementNodes(IntArray &loc, Element &element, const IntArray &bNodes, const UnknownNumberingScheme &s, IntArray *dofIds = nullptr) const;
};

/**
 * Callback class for assembling specific types of matrices
 * @author Mikael Öhman
 */
class OOFEM_EXPORT MatrixAssembler
{
public:
    virtual void matrixFromElement(FloatMatrix &mat, Element &element, TimeStep *tStep) const;
    virtual void matrixFromLoad(FloatMatrix &mat, Element &element, BodyLoad *load, TimeStep *tStep) const;
    virtual void matrixFromSurfaceLoad(FloatMatrix &mat, Element &element, SurfaceLoad *load, int boundary, TimeStep *tStep) const;
    virtual void matrixFromEdgeLoad(FloatMatrix &mat, Element &element, EdgeLoad *load, int edge, TimeStep *tStep) const;
    virtual void assembleFromActiveBC(SparseMtrx &k, ActiveBoundaryCondition &bc, TimeStep* tStep, const UnknownNumberingScheme &s_r, const UnknownNumberingScheme &s_c, void* lock=nullptr) const;

    virtual void locationFromElement(IntArray &loc, Element &element, const UnknownNumberingScheme &s, IntArray *dofIds = nullptr) const;
    virtual void locationFromElementNodes(IntArray &loc, Element &element, const IntArray &bNodes, const UnknownNumberingScheme &s, IntArray *dofIds = nullptr) const;
};


/**
 * Implementation for assembling internal forces vectors in standard monolithic, nonlinear FE-problems
 * @author Mikael Öhman
 */
class OOFEM_EXPORT InternalForceAssembler : public VectorAssembler
{
public:
    void vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const override;
    void vectorFromLoad(FloatArray &vec, Element &element, BodyLoad *load, TimeStep *tStep, ValueModeType mode) const override;
    void vectorFromSurfaceLoad(FloatArray &vec, Element &element, SurfaceLoad *load, int boundary, TimeStep *tStep, ValueModeType mode) const override;
    void vectorFromEdgeLoad(FloatArray &vec, Element &element, EdgeLoad *load, int edge, TimeStep *tStep, ValueModeType mode) const override;
    void assembleFromActiveBC(FloatArray &answer, ActiveBoundaryCondition &bc, TimeStep* tStep, ValueModeType mode, const UnknownNumberingScheme &s, FloatArray *eNorms, void* lock=nullptr) const override;
};

/**
 * Implementation for assembling external forces vectors in standard monolithic FE-problems
 * @author Mikael Öhman
 */
class OOFEM_EXPORT ExternalForceAssembler : public VectorAssembler
{
public:
    void vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const override; ///@todo Temporary: Remove when switch to sets is complete
    void vectorFromLoad(FloatArray &vec, Element &element, BodyLoad *load, TimeStep *tStep, ValueModeType mode) const override;
    void vectorFromSurfaceLoad(FloatArray &vec, Element &element, SurfaceLoad *load, int boundary, TimeStep *tStep, ValueModeType mode) const override;
    void vectorFromEdgeLoad(FloatArray &vec, Element &element, EdgeLoad *load, int edge, TimeStep *tStep, ValueModeType mode) const override;
    void vectorFromNodeLoad(FloatArray &vec, DofManager &dman, NodalLoad *load, TimeStep *tStep, ValueModeType mode) const override;
    void assembleFromActiveBC(FloatArray &answer, ActiveBoundaryCondition &bc, TimeStep* tStep, ValueModeType mode, const UnknownNumberingScheme &s, FloatArray *eNorms, void*lock=nullptr) const override;
};

/**
 * Implementation for assembling reference (external) forces vectors
 * @author Mikael Öhman
 */
class OOFEM_EXPORT ReferenceForceAssembler : public VectorAssembler
{
public:
    void vectorFromLoad(FloatArray &vec, Element &element, BodyLoad *load, TimeStep *tStep, ValueModeType mode) const override;
    void vectorFromSurfaceLoad(FloatArray &vec, Element &element, SurfaceLoad *load, int boundary, TimeStep *tStep, ValueModeType mode) const override;
    void vectorFromEdgeLoad(FloatArray &vec, Element &element, EdgeLoad *load, int edge, TimeStep *tStep, ValueModeType mode) const override;
    void vectorFromNodeLoad(FloatArray &vec, DofManager &dman, NodalLoad *load, TimeStep *tStep, ValueModeType mode) const override;
};

/**
 * Implementation for assembling lumped mass matrix (diagonal components) in vector form.
 * @author Mikael Öhman
 */
class OOFEM_EXPORT LumpedMassVectorAssembler : public VectorAssembler
{
public:
    void vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const override;
};

/**
 * Implementation for assembling the intertia forces vector (i.e. C * dT/dt or M * a)
 * @author Mikael Öhman
 */
class OOFEM_EXPORT InertiaForceAssembler : public VectorAssembler
{
public:
    void vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const override;
};


/**
 * Implementation for assembling forces computed by multiplication with a matrix.
 * This is useful for computing; f = K * u for extrapolated forces, without constructing the K-matrix.
 * @author Mikael Öhman
 */
class OOFEM_EXPORT MatrixProductAssembler : public VectorAssembler
{
protected:
    const MatrixAssembler &mAssem;

public:
    MatrixProductAssembler(const MatrixAssembler &m): VectorAssembler(), mAssem(m) {}

    void vectorFromElement(FloatArray &vec, Element &element, TimeStep *tStep, ValueModeType mode) const override;
    void vectorFromLoad(FloatArray &vec, Element &element, BodyLoad *load, TimeStep *tStep, ValueModeType mode) const override;
    void vectorFromSurfaceLoad(FloatArray &vec, Element &element, SurfaceLoad *load, int boundary, TimeStep *tStep, ValueModeType mode) const override;
    void vectorFromEdgeLoad(FloatArray &vec, Element &element, EdgeLoad *load, int edge, TimeStep *tStep, ValueModeType mode) const override;
};


/**
 * Implementation for assembling tangent matrices in standard monolithic FE-problems
 * @author Mikael Öhman
 */
class OOFEM_EXPORT TangentAssembler : public MatrixAssembler
{
protected:
    ///@todo This is more general than just material responses; we should make a "TangentType"
    MatResponseMode rmode;

public:
    TangentAssembler(MatResponseMode m = TangentStiffness): MatrixAssembler(), rmode(m) {}

    void matrixFromElement(FloatMatrix &mat, Element &element, TimeStep *tStep) const override;
    void matrixFromLoad(FloatMatrix &mat, Element &element, BodyLoad *load, TimeStep *tStep) const override;
    void matrixFromSurfaceLoad(FloatMatrix &mat, Element &element, SurfaceLoad *load, int boundary, TimeStep *tStep) const override;
    void matrixFromEdgeLoad(FloatMatrix &mat, Element &element, EdgeLoad *load, int edge, TimeStep *tStep) const override;
    void assembleFromActiveBC(SparseMtrx &k, ActiveBoundaryCondition &bc, TimeStep* tStep, const UnknownNumberingScheme &s_r, const UnknownNumberingScheme &s_c, void*lock=nullptr) const override;
};


/**
 * Implementation for assembling the consistent mass matrix
 * @author Mikael Öhman
 */
class OOFEM_EXPORT MassMatrixAssembler : public MatrixAssembler
{
public:
    void matrixFromElement(FloatMatrix &mat, Element &element, TimeStep *tStep) const override;
};


/**
 * Callback class for assembling effective tangents composed of stiffness and mass matrix.
 * @author Mikael Öhman
 */
class OOFEM_EXPORT EffectiveTangentAssembler : public MatrixAssembler
{
protected:
    MatResponseMode rmode;
    double lumped;
    double k, m;

public:
    EffectiveTangentAssembler(MatResponseMode mode, bool lumped, double k, double m);
    void matrixFromElement(FloatMatrix &mat, Element &element, TimeStep *tStep) const override;
    void matrixFromLoad(FloatMatrix &mat, Element &element, BodyLoad *load, TimeStep *tStep) const override;
    void matrixFromSurfaceLoad(FloatMatrix &mat, Element &element, SurfaceLoad *load, int boundary, TimeStep *tStep) const override;
    void matrixFromEdgeLoad(FloatMatrix &mat, Element &element, EdgeLoad *load, int edge, TimeStep *tStep) const override;
    void assembleFromActiveBC(SparseMtrx &k, ActiveBoundaryCondition &bc, TimeStep* tStep, const UnknownNumberingScheme &s_r, const UnknownNumberingScheme &s_c, void *lock=nullptr) const override;
};

}
#endif // assemblercallback_h
