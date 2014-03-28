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
 *             OOFEM : Object Oriented Finite ElementGeometry Code
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

#include "elementevaluator.h"
#include "crosssection.h"
#include "integrationrule.h"
#include "errorestimator.h"
#include "intarray.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "primaryfield.h"
#include "verbose.h"
#include "entityrenumberingscheme.h"
#include "error.h"
#include "classfactory.h"
#include "datastream.h"
#include "materialmapperinterface.h"
#include "contextioerr.h"
#include "mathfem.h"
#include "feinterpol.h"
#include "feinterpol1d.h"
#include "feinterpol2d.h"
#include "feinterpol3d.h"
#include "function.h"
#include "dofmanager.h"
#include "node.h"
#include "dynamicinputrecord.h"
#include "matstatmapperint.h"

#include <cstdio>

namespace oofem {

ElementEvaluator::ElementEvaluator() : bodyLoadArray(), boundaryLoadArray()
{
		
}
	


void
ElementEvaluator :: giveCharacteristicMatrix(FloatMatrix &answer,
                                    CharType mtrx, TimeStep *tStep)
//
// returns characteristics matrix of receiver according to mtrx
//
{
    OOFEM_SIMPLE_ERROR("giveCharacteristicMatrix: Unknown Type of characteristic mtrx.");
}


void
ElementEvaluator::giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep)
//
// returns characteristics vector of receiver according to mtrx
//
{
    OOFEM_SIMPLE_ERROR("giveCharacteristicVector: Unknown Type of characteristic mtrx.");
}


void
ElementEvaluator::computeLoadVector(FloatArray &answer, Load *load, CharType type, ValueModeType mode, TimeStep *tStep)
{
    OOFEM_SIMPLE_ERROR("computeLoadVector: Unknown load type.");
}


void
ElementEvaluator::computeBoundaryLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep)
{
    OOFEM_SIMPLE_ERROR("computeBoundaryLoadVector: Unknown load type.");
}


void
ElementEvaluator::computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int edge, CharType type, ValueModeType mode, TimeStep *tStep)
{
    ///@todo Change the load type to "BoundaryEdgeLoad" maybe?
    OOFEM_SIMPLE_ERROR("computeBoundaryEdgeLoadVector: Unknown load type.");
}


double
ElementEvaluator::giveCharacteristicValue(CharType mtrx, TimeStep *tStep)
//
// returns characteristics value of receiver according to CharType
//
{
    OOFEM_SIMPLE_ERROR("giveCharacteristicValue: Unknown Type of characteristic mtrx.");
    return 0.;
}

void ElementEvaluator::computeVectorOf(EquationID type, ValueModeType u, TimeStep *tStep, FloatArray &answer, ElementGeometry* elementGeometry)
// Forms the vector containing the values of the unknown 'u' (e.g., the
// Total value) of the dofs of the callers local cs.
{
	int k, nDofs;
	IntArray dofIDMask;
	FloatMatrix G2L;
	FloatArray vec;
	answer.resize(this->computeNumberOfGlobalDofs());
	k = 0;
	for (int i = 1; i <= elementGeometry->giveNumberOfDofManagers(); i++) {
		elementGeometry->giveDofManDofIDMask(i, type, dofIDMask);
		elementGeometry->giveDofManager(i)->giveUnknownVector(vec, dofIDMask, u, tStep);
		nDofs = vec.giveSize();
		for (int j = 1; j <= nDofs; j++) {
			answer.at(++k) = vec.at(j);
		}
	}
	
	for (int i = 1; i <= elementGeometry->giveNumberOfInternalDofManagers(); i++) {
		elementGeometry->giveInternalDofManDofIDMask(i, type, dofIDMask);
		elementGeometry->giveInternalDofManager(i)->giveUnknownVector(vec, dofIDMask, u, tStep);
		nDofs = vec.giveSize();
		for (int j = 1; j <= nDofs; j++) {
			answer.at(++k) = vec.at(j);
		}
	}
	answer.resizeWithValues(k);
	
	if (this->computeGtoLRotationMatrix(G2L)) {
		answer.rotatedWith(G2L, 'n');
	}
}

void ElementEvaluator::computeBoundaryVectorOf(const IntArray &bNodes, EquationID type, ValueModeType u, TimeStep *tStep, FloatArray &answer, ElementGeometry* elementGeometry)
// Forms the vector containing the values of the unknown 'u' (e.g., the
// Total value) of the dofs of the callers local cs.
{
	int k;
	IntArray dofIDMask;
	FloatMatrix G2L;
	FloatArray vec;
	
	k = 0;
	for (int i = 1; i <= bNodes.giveSize(); i++) {
		elementGeometry->giveDofManDofIDMask(bNodes.at(i), type, dofIDMask);
		k += dofIDMask.giveSize();
	}
	answer.resize(k);

	k = 0;
	for (int i = 1; i <= bNodes.giveSize(); i++) {
		elementGeometry->giveDofManDofIDMask(bNodes.at(i), type, dofIDMask);
		elementGeometry->giveDofManager(bNodes.at(i))->giveUnknownVector(vec, dofIDMask, u, tStep);
		for (int j = 1; j <= vec.giveSize(); j++) {
			answer.at(++k) = vec.at(j);
		}
	}
	
	if (this->computeGtoLRotationMatrix(G2L)) {
		OOFEM_SIMPLE_ERROR("ElementEvaluator :: computeBoundaryVector - Local coordinate system is not implemented yet");
	}
}


void ElementEvaluator::computeVectorOf(PrimaryField &field, ValueModeType u, TimeStep *tStep, FloatArray &answer, ElementGeometry* elementGeometry)
// Forms the vector containing the values of the unknown 'u' (e.g., the
// Total value) of the dofs of the receiver's nodes (in nodal cs).
// Dofs containing expected unknowns (of expected type) are determined
// using this->GiveNodeDofIDMask function
{
	int k, nDofs;
	IntArray dofIDMask;
	FloatMatrix G2L;
	FloatArray vec;
	answer.resize(this->computeNumberOfGlobalDofs());

	k = 0;
	for (int i = 1; i <= elementGeometry->giveNumberOfDofManagers(); i++) {
		elementGeometry->giveDofManDofIDMask(i, field.giveEquationID(), dofIDMask);
		elementGeometry->giveDofManager(i)->giveUnknownVector(vec, dofIDMask, field, u, tStep);
		nDofs = vec.giveSize();
		for (int j = 1; j <= nDofs; j++) {
			answer.at(++k) = vec.at(j);
		}
	}

	for (int i = 1; i <= elementGeometry->giveNumberOfInternalDofManagers(); i++) {
		elementGeometry->giveInternalDofManDofIDMask(i, field.giveEquationID(), dofIDMask);
		elementGeometry->giveInternalDofManager(i)->giveUnknownVector(vec, dofIDMask, field, u, tStep);
		nDofs = vec.giveSize();
		for (int j = 1; j <= nDofs; j++) {
			answer.at(++k) = vec.at(j);
		}
	}
	answer.resizeWithValues(k);

	if (this->computeGtoLRotationMatrix(G2L)) {
		answer.rotatedWith(G2L, 'n');
	}
}


void ElementEvaluator::computeVectorOfPrescribed(EquationID ut, ValueModeType mode, TimeStep *tStep, FloatArray &answer, ElementGeometry* elementGeometry)
// Forms the vector containing the prescribed values of the unknown 'u'
// (e.g., the prescribed displacement) of the dofs of the receiver's
// nodes. Puts 0 at each free dof.
{
	int k, nDofs;
	IntArray dofIDMask, dofMask;
	FloatMatrix G2L;
	FloatArray vec;

	answer.resize(this->computeNumberOfGlobalDofs());
	
	k = 0;
	for (int i = 1; i <= elementGeometry->giveNumberOfDofManagers(); i++) {
		elementGeometry->giveDofManDofIDMask(i, ut, dofIDMask);
		elementGeometry->giveDofManager(i)->givePrescribedUnknownVector(vec, dofIDMask, mode, tStep);
		nDofs = vec.giveSize();
		for (int j = 1; j <= nDofs; j++) {
			answer.at(++k) = vec.at(j);
		}
	}

	for (int i = 1; i <= elementGeometry->giveNumberOfInternalDofManagers(); i++) {
		elementGeometry->giveInternalDofManDofIDMask(i, ut, dofIDMask);
		elementGeometry->giveInternalDofManager(i)->givePrescribedUnknownVector(vec, dofIDMask, mode, tStep);
		nDofs = vec.giveSize();
		for (int j = 1; j <= nDofs; j++) {
			answer.at(++k) = vec.at(j);
		}
	}
	answer.resizeWithValues(k);

	if (this->computeGtoLRotationMatrix(G2L)) {
		answer.rotatedWith(G2L, 'n');
	}
}


int	ElementEvaluator::computeNumberOfGlobalDofs()
{
	return this->computeNumberOfDofs();
}


int	ElementEvaluator::computeNumberOfPrimaryMasterDofs(EquationID ut, ElementGeometry* elementGeometry)
{
	int answer = 0;
	IntArray nodeDofIDMask, dofMask;

	for (int i = 1; i <= elementGeometry->giveNumberOfDofManagers(); i++) {
		elementGeometry->giveDofManDofIDMask(i, ut, nodeDofIDMask);
		elementGeometry->giveDofManager(i)->giveDofArray(nodeDofIDMask, dofMask);
		answer += elementGeometry->giveDofManager(i)->giveNumberOfPrimaryMasterDofs(dofMask);
	}

	for (int i = 1; i <= elementGeometry->giveNumberOfInternalDofManagers(); i++) {
		elementGeometry->giveInternalDofManDofIDMask(i, ut, nodeDofIDMask);
		elementGeometry->giveInternalDofManager(i)->giveDofArray(nodeDofIDMask, dofMask);
		answer += elementGeometry->giveInternalDofManager(i)->giveNumberOfPrimaryMasterDofs(dofMask);
	}
	return answer;
}


bool
ElementEvaluator::computeGtoLRotationMatrix(FloatMatrix &answer)
{
	answer.clear();
	return false;
}

bool ElementEvaluator::giveRotationMatrix(FloatMatrix &answer, EquationID eid, ElementGeometry* elementGeometry)
{
	bool is_GtoL, is_NtoG;
	FloatMatrix GtoL, NtoG;
	IntArray nodes;
	nodes.enumerate(elementGeometry->giveNumberOfDofManagers());

	is_GtoL = this->computeGtoLRotationMatrix(GtoL);
	is_NtoG = this->computeDofTransformationMatrix(NtoG, nodes, true, eid, elementGeometry);

#ifdef DEBUG
	if (is_GtoL) {
		if (GtoL.giveNumberOfColumns() != this->computeNumberOfGlobalDofs()) {
			OOFEM_SIMPLE_ERROR("ElementEvaluator :: updateRotationMatrix - GtoL transformation matrix size mismatch in columns");
		}
		if (GtoL.giveNumberOfRows() != this->computeNumberOfDofs()) {
			OOFEM_SIMPLE_ERROR("ElementEvaluator :: updateRotationMatrix - GtoL transformation matrix size mismatch in rows");
		}
	}
	if (is_NtoG) {
		if (NtoG.giveNumberOfColumns() != this->computeNumberOfPrimaryMasterDofs(eid, elementGeometry)) {
			OOFEM_SIMPLE_ERROR("ElementEvaluator :: updateRotationMatrix - NtoG transformation matrix size mismatch in columns");
		}
		if (NtoG.giveNumberOfRows() != this->computeNumberOfGlobalDofs()) {
			OOFEM_SIMPLE_ERROR("ElementEvaluator :: updateRotationMatrix - NtoG transformation matrix size mismatch in rows");
		}
	}
#endif

	if (is_GtoL && NtoG.isNotEmpty()) {
		answer.beProductOf(GtoL, NtoG);
	} else if (is_GtoL) {
		answer = GtoL;
	} else if (is_NtoG) {
		answer = NtoG;
	} else {
		answer.clear();
		return false;
	}
	return true;
}


bool ElementEvaluator::computeDofTransformationMatrix(FloatMatrix &answer, const IntArray &nodes, bool includeInternal, EquationID eid, ElementGeometry* elementGeometry)
{
	bool flag = false;
	int numberOfDofMans = nodes.giveSize();

	// test if transformation is necessary
	for (int i = 1; i <= numberOfDofMans; i++) {
		flag = flag || elementGeometry->giveDofManager(nodes.at(i))->requiresTransformation();
	}
	
	if (!flag) {
		answer.clear();
		return false;
	}
	// initialize answer
	int gsize = this->computeNumberOfPrimaryMasterDofs(eid, elementGeometry);
	answer.resize(this->computeNumberOfGlobalDofs(), gsize);
	answer.zero();

	FloatMatrix dofManT;
	IntArray dofIDmask;
	int nr, nc, lastRowPos = 0, lastColPos = 0;
	// loop over nodes
	for (int i = 1; i <= numberOfDofMans; i++) {
		elementGeometry->giveDofManDofIDMask(nodes.at(i), eid, dofIDmask);
		if (!elementGeometry->giveDofManager(nodes.at(i))->computeM2GTransformation(dofManT, dofIDmask)) {
			dofManT.resize(dofIDmask.giveSize(), dofIDmask.giveSize());
			dofManT.zero();
			dofManT.beUnitMatrix();
		}
		nc = dofManT.giveNumberOfColumns();
		nr = dofManT.giveNumberOfRows();
		for (int j = 1; j <= nr; j++) {
			for (int k = 1; k <= nc; k++) {
				// localize node contributions
				answer.at(lastRowPos + j, lastColPos + k) = dofManT.at(j, k);
			}
		}
		lastRowPos += nr;
		lastColPos += nc;
	}
	if (includeInternal) {
		for (int i = 1; i <= elementGeometry->giveNumberOfInternalDofManagers(); i++) {
			elementGeometry->giveInternalDofManDofIDMask(i, eid, dofIDmask);
			if (!elementGeometry->giveInternalDofManager(nodes.at(i))->computeM2GTransformation(dofManT, dofIDmask)) {
				dofManT.resize(dofIDmask.giveSize(), dofIDmask.giveSize());
				dofManT.zero();
				dofManT.beUnitMatrix();
			}
			nc = dofManT.giveNumberOfColumns();
			nr = dofManT.giveNumberOfRows();
			for (int j = 1; j <= nr; j++) {
				for (int k = 1; k <= nc; k++) {
					// localize node contributions
					answer.at(lastRowPos + j, lastColPos + k) = dofManT.at(j, k);
				}
			}

			lastRowPos += nr;
			lastColPos += nc;
		}
	}
	answer.resizeWithData(answer.giveNumberOfRows(), lastColPos);
	return true;
}


IntArray *ElementEvaluator::giveBodyLoadArray()
// Returns the array which contains the number of every body load that act
// on the receiver.
{
	return &bodyLoadArray;
}


IntArray *ElementEvaluator::giveBoundaryLoadArray()
// Returns the array which contains the number of every body load that act
// on the receiver.
{
	return &boundaryLoadArray;
}


void ElementEvaluator::giveLocationArray(IntArray &locationArray, EquationID eid, const UnknownNumberingScheme &s, ElementGeometry* elementGeometry, IntArray *dofIdArray) const
// Returns the location array of the receiver. This array is obtained by
// simply appending the location array of every node of the receiver.
{
	IntArray masterDofIDs, nodalArray, dofIDMask;
	locationArray.resize(0);
	if (dofIdArray) {
		dofIdArray->resize(0);
	}
	for (int i = 1; i <= elementGeometry->giveNumberOfDofManagers(); i++) {
		elementGeometry->giveDofManDofIDMask(i, eid, dofIDMask);
		elementGeometry->giveDofManager(i)->giveLocationArray(dofIDMask, nodalArray, s);
		locationArray.followedBy(nodalArray);
		if (dofIdArray) {
			elementGeometry->giveDofManager(i)->giveMasterDofIDArray(dofIDMask, masterDofIDs);
			dofIdArray->followedBy(masterDofIDs);
		}
	}
	for (int i = 1; i <= elementGeometry->giveNumberOfInternalDofManagers(); i++) {
		elementGeometry->giveInternalDofManDofIDMask(i, eid, dofIDMask);
		elementGeometry->giveInternalDofManager(i)->giveLocationArray(dofIDMask, nodalArray, s);
		locationArray.followedBy(nodalArray);
		if (dofIdArray) {
			elementGeometry->giveInternalDofManager(i)->giveMasterDofIDArray(dofIDMask, masterDofIDs);
			dofIdArray->followedBy(masterDofIDs);
		}
	}
}

void ElementEvaluator::giveLocationArray(IntArray &locationArray, const IntArray &dofIDMask, const UnknownNumberingScheme &s, ElementGeometry* elementGeometry, IntArray *dofIdArray) const
{
	IntArray masterDofIDs, nodalArray, ids = dofIDMask;
	locationArray.resize(0);
	if (dofIdArray) {
		dofIdArray->resize(0);
	}
	for (int i = 1; i <= elementGeometry->giveNumberOfDofManagers(); i++) {
		if (dofIDMask.giveSize() == 0) {
			elementGeometry->giveDefaultDofManDofIDMask(i, ids);
		}
		elementGeometry->giveDofManager(i)->giveLocationArray(ids, nodalArray, s);
		locationArray.followedBy(nodalArray);
		if (dofIdArray) {
			elementGeometry->giveDofManager(i)->giveMasterDofIDArray(ids, masterDofIDs);
			dofIdArray->followedBy(masterDofIDs);
		}
	}
	for (int i = 1; i <= elementGeometry->giveNumberOfInternalDofManagers(); i++) {
		if (dofIDMask.giveSize() == 0) {
			elementGeometry->giveDefaultInternalDofManDofIDMask(i, ids);
		}
		elementGeometry->giveInternalDofManager(i)->giveLocationArray(ids, nodalArray, s);
		locationArray.followedBy(nodalArray);
		if (dofIdArray) {
			elementGeometry->giveInternalDofManager(i)->giveMasterDofIDArray(ids, masterDofIDs);
			dofIdArray->followedBy(masterDofIDs);
		}
	}
}


void ElementEvaluator::giveBoundaryLocationArray(IntArray &locationArray, const IntArray &bNodes, EquationID eid, const UnknownNumberingScheme &s, ElementGeometry* elementGeometry, IntArray *dofIdArray)
{
	IntArray masterDofIDs, nodalArray, dofIDMask;
	locationArray.resize(0);
	if (dofIdArray) {
		dofIdArray->resize(0);
	}
	
	for (int i = 1; i <= bNodes.giveSize(); i++) {
		elementGeometry->giveDofManDofIDMask(bNodes.at(i), eid, dofIDMask);
		elementGeometry->giveDofManager(bNodes.at(i))->giveLocationArray(dofIDMask, nodalArray, s);
		locationArray.followedBy(nodalArray);
		if (dofIdArray) {
			elementGeometry->giveDofManager(bNodes.at(i))->giveMasterDofIDArray(dofIDMask, masterDofIDs);
			dofIdArray->followedBy(masterDofIDs);
		}
	}
}


void ElementEvaluator::giveBoundaryLocationArray(IntArray &locationArray, const IntArray &bNodes, const IntArray &dofIDMask, const UnknownNumberingScheme &s, ElementGeometry* elementGeometry, IntArray *dofIdArray)
{
	IntArray masterDofIDs, nodalArray, ids = dofIDMask;
	locationArray.resize(0);
	if (dofIdArray) {
		dofIdArray->resize(0);
	}
	
	for (int i = 1; i <= bNodes.giveSize(); i++) {
		if (dofIDMask.giveSize() == 0) {
			elementGeometry->giveDefaultDofManDofIDMask(bNodes.at(i), ids);
		}
		elementGeometry->giveDofManager(bNodes.at(i))->giveLocationArray(ids, nodalArray, s);
		locationArray.followedBy(nodalArray);
		if (dofIdArray) {
			elementGeometry->giveDofManager(bNodes.at(i))->giveMasterDofIDArray(ids, masterDofIDs);
			dofIdArray->followedBy(masterDofIDs);
		}
	}
}

IRResultType ElementEvaluator::initializeFrom(InputRecord *ir)
{
	IRResultType result;                          // Required by IR_GIVE_FIELD macro
#  ifdef VERBOSE
			// VERBOSE_PRINT1("Instanciating element ",number);
#  endif			
	bodyLoadArray.resize(0);
	IR_GIVE_OPTIONAL_FIELD(ir, bodyLoadArray, _IFT_ElementEvaluator_bodyload);

	boundaryLoadArray.resize(0);
	IR_GIVE_OPTIONAL_FIELD(ir, boundaryLoadArray, _IFT_ElementEvaluator_boundaryload);
	
	return IRRT_OK;
}


void ElementEvaluator::giveInputRecord(DynamicInputRecord &input)
{
	if (bodyLoadArray.giveSize() > 0) {
		input.setField(bodyLoadArray, _IFT_ElementEvaluator_bodyload);
	}
	if (boundaryLoadArray.giveSize() > 0) {
		input.setField(boundaryLoadArray, _IFT_ElementEvaluator_boundaryload);
	}
	

}



contextIOResultType ElementEvaluator::saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full element context (saves state variables, that completely describe
// current state)
//
{
	contextIOResultType iores;
	if ((iores = bodyLoadArray.storeYourself(stream, mode)) != CIO_OK) {
		THROW_CIOERR(iores);
	}
	if ((iores = boundaryLoadArray.storeYourself(stream, mode)) != CIO_OK) {
		THROW_CIOERR(iores);
	}

	return CIO_OK;
}


contextIOResultType ElementEvaluator::restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full element context (saves state variables, that completely describe
// current state)
//
{
	contextIOResultType iores;
	
	if (mode & CM_Definition) {
		if ((iores = bodyLoadArray.restoreYourself(stream, mode)) != CIO_OK) {
			THROW_CIOERR(iores);
		}
		if ((iores = boundaryLoadArray.restoreYourself(stream, mode)) != CIO_OK) {
			THROW_CIOERR(iores);
		}
	}

	return CIO_OK;
}


} // end namespace oofem

