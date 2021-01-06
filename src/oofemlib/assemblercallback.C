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

#include "assemblercallback.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "element.h"
#include "dofmanager.h"
#include "activebc.h"

#include "nodalload.h"
#include "bodyload.h"
#include "boundaryload.h"

namespace oofem {

void VectorAssembler :: vectorFromElement(FloatArray& vec, Element& element, TimeStep* tStep, ValueModeType mode) const { vec.clear(); }

void VectorAssembler :: vectorFromLoad(FloatArray& vec, Element& element, BodyLoad* load, TimeStep* tStep, ValueModeType mode) const { vec.clear(); }

void VectorAssembler :: vectorFromSurfaceLoad(FloatArray& vec, Element& element, SurfaceLoad* load, int boundary, TimeStep* tStep, ValueModeType mode) const { vec.clear(); }

void VectorAssembler :: vectorFromEdgeLoad(FloatArray& vec, Element& element, EdgeLoad* load, int edge, TimeStep* tStep, ValueModeType mode) const { vec.clear(); }

void VectorAssembler :: vectorFromNodeLoad(FloatArray& vec, DofManager& dman, NodalLoad* load, TimeStep* tStep, ValueModeType mode) const { vec.clear(); }

void VectorAssembler :: assembleFromActiveBC(FloatArray &answer, ActiveBoundaryCondition &bc, TimeStep* tStep, ValueModeType mode, const UnknownNumberingScheme &s, FloatArray *eNorms, void* lock) const { }

void VectorAssembler :: locationFromElement(IntArray& loc, Element& element, const UnknownNumberingScheme& s, IntArray* dofIds) const
{
    element.giveLocationArray(loc, s, dofIds);
}

void VectorAssembler :: locationFromElementNodes(IntArray& loc, Element& element, const IntArray& bNodes, const UnknownNumberingScheme& s, IntArray* dofIds) const
{
    element.giveBoundaryLocationArray(loc, bNodes, s, dofIds);
}

void MatrixAssembler :: matrixFromElement(FloatMatrix& mat, Element& element, TimeStep* tStep) const { mat.clear(); }

void MatrixAssembler :: matrixFromLoad(FloatMatrix& mat, Element& element, BodyLoad* load, TimeStep* tStep) const { mat.clear(); }

void MatrixAssembler :: matrixFromSurfaceLoad(FloatMatrix& mat, Element& element, SurfaceLoad* load, int boundary, TimeStep* tStep) const { mat.clear(); }

void MatrixAssembler :: matrixFromEdgeLoad(FloatMatrix& mat, Element& element, EdgeLoad* load, int edge, TimeStep* tStep) const { mat.clear(); }

void MatrixAssembler :: assembleFromActiveBC(SparseMtrx &k, ActiveBoundaryCondition &bc, TimeStep* tStep, const UnknownNumberingScheme &s_r, const UnknownNumberingScheme &s_c, void *lock) const {}

void MatrixAssembler :: locationFromElement(IntArray& loc, Element& element, const UnknownNumberingScheme& s, IntArray* dofIds) const
{
    element.giveLocationArray(loc, s, dofIds);
}

void MatrixAssembler :: locationFromElementNodes(IntArray& loc, Element& element, const IntArray& bNodes, const UnknownNumberingScheme& s, IntArray* dofIds) const
{
    element.giveBoundaryLocationArray(loc, bNodes, s, dofIds);
}


void MatrixProductAssembler :: vectorFromElement(FloatArray& vec, Element& element, TimeStep* tStep, ValueModeType mode) const
{
    FloatMatrix mat;
    this->mAssem.matrixFromElement(mat, element, tStep);
    vec.beProductOf(mat, this->vec);
}

void MatrixProductAssembler :: vectorFromLoad(FloatArray& vec, Element& element, BodyLoad* load, TimeStep* tStep, ValueModeType mode) const
{
    FloatMatrix mat;
    this->mAssem.matrixFromLoad(mat, element, load, tStep);
    vec.beProductOf(mat, this->vec);
}

void MatrixProductAssembler :: vectorFromSurfaceLoad(FloatArray& vec, Element& element, SurfaceLoad* load, int boundary, TimeStep* tStep, ValueModeType mode) const
{
    FloatMatrix mat;
    this->mAssem.matrixFromSurfaceLoad(mat, element, load, boundary, tStep);
    vec.beProductOf(mat, this->vec);
}

void MatrixProductAssembler :: vectorFromEdgeLoad(FloatArray& vec, Element& element, EdgeLoad* load, int edge, TimeStep* tStep, ValueModeType mode) const
{
    FloatMatrix mat;
    this->mAssem.matrixFromEdgeLoad(mat, element, load, edge, tStep); // ??????
    vec.beProductOf(mat, this->vec);
}


void InternalForceAssembler :: vectorFromElement(FloatArray& vec, Element& element, TimeStep* tStep, ValueModeType mode) const
{
    element.giveCharacteristicVector(vec, InternalForcesVector, mode, tStep);
    //element.computeInternalForces(vec, tStep);
}

void InternalForceAssembler :: vectorFromLoad(FloatArray& vec, Element& element, BodyLoad* load, TimeStep* tStep, ValueModeType mode) const
{
    element.computeLoadVector(vec, load, InternalForcesVector, mode, tStep);
    //element.computeInternalForcesFromLoad(vec, load, tStep);
}

void InternalForceAssembler :: vectorFromSurfaceLoad(FloatArray& vec, Element& element, SurfaceLoad* load, int boundary, TimeStep* tStep, ValueModeType mode) const
{
    element.computeBoundarySurfaceLoadVector(vec, load, boundary, InternalForcesVector, mode, tStep);
    //element.computeInternalForcesFromBoundaryLoad(vec, load, boundary, tStep);
}

void InternalForceAssembler :: vectorFromEdgeLoad(FloatArray& vec, Element& element, EdgeLoad* load, int edge, TimeStep* tStep, ValueModeType mode) const
{
    element.computeBoundaryEdgeLoadVector(vec, load, edge, InternalForcesVector, mode, tStep);
    //element.computeInternalForcesFromEdgeLoad(vec, load, edge, tStep);
}

void InternalForceAssembler :: assembleFromActiveBC(FloatArray &answer, ActiveBoundaryCondition &bc, TimeStep* tStep, ValueModeType mode, const UnknownNumberingScheme &s, FloatArray *eNorms, void* lock) const
{
    bc.assembleVector(answer, tStep, InternalForcesVector, mode, s, eNorms, lock);
    //bc.assembleInternalForces(answer, tStep, s, eNorms);
}


void ExternalForceAssembler :: vectorFromElement(FloatArray& vec, Element& element, TimeStep* tStep, ValueModeType mode) const
{
    ///@todo To be removed when sets are used for loads.
    element.giveCharacteristicVector(vec, ExternalForcesVector, mode, tStep);
}

void ExternalForceAssembler :: vectorFromLoad(FloatArray& vec, Element& element, BodyLoad* load, TimeStep* tStep, ValueModeType mode) const
{
    vec.clear();
    if ( ! load->reference )
        element.computeLoadVector(vec, load, ExternalForcesVector, mode, tStep);
        //element.computeExternalForcesFromLoad(vec, load, tStep);
}

void ExternalForceAssembler :: vectorFromSurfaceLoad(FloatArray& vec, Element& element, SurfaceLoad* load, int boundary, TimeStep* tStep, ValueModeType mode) const
{
    vec.clear();
    if ( ! load->reference )
        element.computeBoundarySurfaceLoadVector(vec, load, boundary, ExternalForcesVector, mode, tStep);
        //element.computeExternalForcesFromSurfaceLoad(vec, load, boundary, tStep);
}

void ExternalForceAssembler :: vectorFromEdgeLoad(FloatArray& vec, Element& element, EdgeLoad* load, int edge, TimeStep* tStep, ValueModeType mode) const
{
    vec.clear();
    if ( ! load->reference )
        element.computeBoundaryEdgeLoadVector(vec, load, edge, ExternalForcesVector, mode, tStep);
        //element.computeExternalForcesFromEdgeLoad(vec, load, edge, tStep);
}

void ExternalForceAssembler :: vectorFromNodeLoad(FloatArray& vec, DofManager& dman, NodalLoad* load, TimeStep* tStep, ValueModeType mode) const
{
    vec.clear();
    if ( ! load->reference )
        dman.computeLoadVector(vec, load, ExternalForcesVector, tStep, mode);
        //dman.computeExternalForcesFromLoad(vec, load, tStep);
}

void ExternalForceAssembler :: assembleFromActiveBC(FloatArray &answer, ActiveBoundaryCondition &bc, TimeStep* tStep, ValueModeType mode, const UnknownNumberingScheme &s, FloatArray *eNorms, void* lock) const
{
    bc.assembleVector(answer, tStep, ExternalForcesVector, mode, s, eNorms, lock);
    //bc.assembleExternalForces(answer, tStep, s, eNorms);
}


void ReferenceForceAssembler :: vectorFromLoad(FloatArray& vec, Element& element, BodyLoad* load, TimeStep* tStep, ValueModeType mode) const
{
    vec.clear();
    if ( load->reference )
        element.computeLoadVector(vec, load, ExternalForcesVector, mode, tStep);
    //element.computeExternalForcesFromLoad(vec, load, tStep);
}


void ReferenceForceAssembler :: vectorFromSurfaceLoad(FloatArray& vec, Element& element, SurfaceLoad* load, int boundary, TimeStep* tStep, ValueModeType mode) const
{
    vec.clear();
    if ( load->reference )
        element.computeBoundarySurfaceLoadVector(vec, load, boundary, ExternalForcesVector, mode, tStep);
    //element.computeExternalForcesFromBoundaryLoad(vec, load, boundary, tStep);
}

void ReferenceForceAssembler :: vectorFromEdgeLoad(FloatArray& vec, Element& element, EdgeLoad* load, int edge, TimeStep* tStep, ValueModeType mode) const
{
    vec.clear();
    if ( load->reference )
        element.computeBoundaryEdgeLoadVector(vec, load, edge, ExternalForcesVector, mode, tStep);
    //element.computeExternalForcesFromEdgeLoad(vec, load, edge, tStep);
}

void ReferenceForceAssembler :: vectorFromNodeLoad(FloatArray& vec, DofManager& dman, NodalLoad* load, TimeStep* tStep, ValueModeType mode) const
{
    vec.clear();
    if ( load->reference )
        dman.computeLoadVector(vec, load, ExternalForcesVector, tStep, mode);
        //dman.computeExternalForcesFromLoad(vec, load, tStep);
}


void LumpedMassVectorAssembler :: vectorFromElement(FloatArray& vec, Element& element, TimeStep* tStep, ValueModeType mode) const
{
    element.giveCharacteristicVector(vec, LumpedMassMatrix, mode, tStep);
    //element.computeLumpedMassMatrix(vec, tStep);
}


void InertiaForceAssembler :: vectorFromElement(FloatArray& vec, Element& element, TimeStep* tStep, ValueModeType mode) const
{
    element.giveCharacteristicVector(vec, InertiaForcesVector, mode, tStep);
    //element.computeInertiaForces(vec, tStep);
}


void TangentAssembler :: matrixFromElement(FloatMatrix& mat, Element& element, TimeStep* tStep) const
{
    if ( this->rmode == TangentStiffness ) {
        element.giveCharacteristicMatrix(mat, TangentStiffnessMatrix, tStep);
    } else if ( this->rmode == ElasticStiffness ) {
        element.giveCharacteristicMatrix(mat, ElasticStiffnessMatrix, tStep);
    } else if ( this->rmode == SecantStiffness ) {
        element.giveCharacteristicMatrix(mat, SecantStiffnessMatrix, tStep);
    }
    //element.computeTangentMatrix(mat, this->rmode, tStep);
}

void TangentAssembler :: matrixFromLoad(FloatMatrix& mat, Element& element, BodyLoad* load, TimeStep* tStep) const
{
    mat.clear();
    //element.computeTangentFromLoad(mat, load, this->rmode, tStep);
}

void TangentAssembler :: matrixFromSurfaceLoad(FloatMatrix& mat, Element& element, SurfaceLoad* load, int boundary, TimeStep* tStep) const
{
    element.computeTangentFromSurfaceLoad(mat, load, boundary, this->rmode, tStep);
}

void TangentAssembler :: matrixFromEdgeLoad(FloatMatrix& mat, Element& element, EdgeLoad* load, int edge, TimeStep* tStep) const
{
  //mat.clear();
  element.computeTangentFromEdgeLoad(mat, load, edge, this->rmode, tStep);
}

void TangentAssembler :: assembleFromActiveBC(SparseMtrx &k, ActiveBoundaryCondition &bc, TimeStep* tStep, const UnknownNumberingScheme &s_r, const UnknownNumberingScheme &s_c, void *lock) const
{
    bc.assemble(k, tStep, TangentStiffnessMatrix, s_r, s_c, 1.0, lock);
}



void MassMatrixAssembler :: matrixFromElement(FloatMatrix& mat, Element& element, TimeStep* tStep) const
{
    element.giveCharacteristicMatrix(mat, MassMatrix, tStep);
    //element.computeMassMatrix(mat, tStep);
}



EffectiveTangentAssembler :: EffectiveTangentAssembler(MatResponseMode mode, bool lumped, double k, double m) :
    MatrixAssembler(), rmode(mode), lumped(lumped), k(k), m(m)
{}

void EffectiveTangentAssembler :: matrixFromElement(FloatMatrix &answer, Element &el, TimeStep *tStep) const
{
    FloatMatrix massMatrix;

    if ( this->rmode == TangentStiffness ) {
        el.giveCharacteristicMatrix(answer, TangentStiffnessMatrix, tStep);
    } else if ( this->rmode == ElasticStiffness ) {
        el.giveCharacteristicMatrix(answer, ElasticStiffnessMatrix, tStep);
    } else if ( this->rmode == SecantStiffness ) {
        el.giveCharacteristicMatrix(answer, SecantStiffnessMatrix, tStep);
    }
    //element.computeTangentMatrix(answer, this->rmode, tStep);
    answer.times(this->k);

    el.giveCharacteristicMatrix(massMatrix, this->lumped ? LumpedMassMatrix : MassMatrix, tStep);
    answer.add(this->m, massMatrix);
}

void EffectiveTangentAssembler :: matrixFromLoad(FloatMatrix& mat, Element& element, BodyLoad* load, TimeStep* tStep) const
{
    mat.clear();
    //element.computeTangentFromLoad(mat, load, this->rmode, tStep);
    //mat.times(this->k);
}

void EffectiveTangentAssembler :: matrixFromSurfaceLoad(FloatMatrix& mat, Element& element, SurfaceLoad* load, int boundary, TimeStep* tStep) const
{
    element.computeTangentFromSurfaceLoad(mat, load, boundary, this->rmode, tStep);
    mat.times(this->k);
}

void EffectiveTangentAssembler :: matrixFromEdgeLoad(FloatMatrix& mat, Element& element, EdgeLoad* load, int edge, TimeStep* tStep) const
{
    mat.clear();
    //element.computeTangentFromEdgeLoad(mat, load, edge, this->rmode, tStep);
    //mat.times(this->k);
}

void EffectiveTangentAssembler :: assembleFromActiveBC(SparseMtrx &k, ActiveBoundaryCondition &bc, TimeStep* tStep, const UnknownNumberingScheme &s_r, const UnknownNumberingScheme &s_c, void* lock) const
{
    // TODO: Crucial part to add: We have to support a scaling factor for this method to support effective tangents.
    bc.assemble(k, tStep, TangentStiffnessMatrix, s_r, s_c, this->k, lock);
}


}
