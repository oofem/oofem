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

#include "element.h"
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
#include "loadtimefunction.h"
#include "dofmanager.h"
#include "node.h"
#include "dynamicinputrecord.h"
#include "matstatmapperint.h"

#include <cstdio>

namespace oofem {
Element :: Element(int n, Domain *aDomain) :
    FEMComponent(n, aDomain), dofManArray(), bodyLoadArray(), boundaryLoadArray()
{
    material           = 0;
    numberOfDofMans    = 0;
    numberOfIntegrationRules = 0;
    activityLtf = 0;
    integrationRulesArray  = NULL;
}


Element :: ~Element()
{
    if ( integrationRulesArray ) {
        for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
            delete integrationRulesArray [ i ];
        }

        delete[] integrationRulesArray;
    }
}

#if 0
void
Element :: computeVectorOf(const IntArray &dofIDMask, ValueModeType u, TimeStep *stepN, FloatArray &answer)
{
    int k, nDofs;
    FloatMatrix G2L;
    FloatArray vec;
    answer.resize( dofIDMask.giveSize() * ( this->giveNumberOfDofManagers() + this->giveNumberOfInternalDofManagers() ) );

    k = 0;
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        this->giveDofManager(i)->giveUnknownVector(vec, dofIDMask, u, stepN);
        nDofs = vec.giveSize();
        for ( int j = 1; j <= nDofs; j++ ) {
            answer.at(++k) = vec.at(j);
        }
    }

    for ( int i = 1; i <= giveNumberOfInternalDofManagers(); i++ ) {
        this->giveInternalDofManager(i)->giveUnknownVector(vec, dofIDMask, u, stepN);
        nDofs = vec.giveSize();
        for ( int j = 1; j <= nDofs; j++ ) {
            answer.at(++k) = vec.at(j);
        }
    }

    if (this->computeGtoLRotationMatrix(G2L)) {
        answer.rotatedWith(G2L, 'n');
    }
}
#endif

void
Element :: computeVectorOf(EquationID type, ValueModeType u, TimeStep *stepN, FloatArray &answer)
// Forms the vector containing the values of the unknown 'u' (e.g., the
// Total value) of the dofs of the callers local cs.
{
    int k, nDofs;
    IntArray dofIDMask;
    FloatMatrix G2L;
    FloatArray vec;
    answer.resize( this->computeNumberOfGlobalDofs() );

    k = 0;
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        this->giveDofManDofIDMask(i, type, dofIDMask);
        this->giveDofManager(i)->giveUnknownVector(vec, dofIDMask, u, stepN);
        nDofs = vec.giveSize();
        for ( int j = 1; j <= nDofs; j++ ) {
            answer.at(++k) = vec.at(j);
        }
    }

    for ( int i = 1; i <= giveNumberOfInternalDofManagers(); i++ ) {
        this->giveInternalDofManDofIDMask(i, type, dofIDMask);
        this->giveInternalDofManager(i)->giveUnknownVector(vec, dofIDMask, u, stepN);
        nDofs = vec.giveSize();
        for ( int j = 1; j <= nDofs; j++ ) {
            answer.at(++k) = vec.at(j);
        }
    }
    answer.resizeWithValues(k);

    if (this->computeGtoLRotationMatrix(G2L)) {
        answer.rotatedWith(G2L, 'n');
    }
}


void
Element :: computeBoundaryVectorOf(const IntArray &bNodes, EquationID type, ValueModeType u, TimeStep *stepN, FloatArray &answer)
// Forms the vector containing the values of the unknown 'u' (e.g., the
// Total value) of the dofs of the callers local cs.
{
    int k;
    IntArray dofIDMask;
    FloatMatrix G2L;
    FloatArray vec;

    k = 0;
    for ( int i = 1; i <= bNodes.giveSize(); i++ ) {
        this->giveDofManDofIDMask(bNodes.at(i), type, dofIDMask);
        k += dofIDMask.giveSize();
    }
    answer.resize(k);

    k = 0;
    for ( int i = 1; i <= bNodes.giveSize(); i++ ) {
        this->giveDofManDofIDMask(bNodes.at(i), type, dofIDMask);
        this->giveDofManager(bNodes.at(i))->giveUnknownVector(vec, dofIDMask, u, stepN);
        for ( int j = 1; j <= vec.giveSize(); j++ ) {
            answer.at(++k) = vec.at(j);
        }
    }

    if (this->computeGtoLRotationMatrix(G2L)) {
        OOFEM_ERROR("Element :: computeBoundaryVector - Local coordinate system is not implemented yet");
    }
}


void
Element :: computeVectorOf(PrimaryField &field, ValueModeType u, TimeStep *stepN, FloatArray &answer)
// Forms the vector containing the values of the unknown 'u' (e.g., the
// Total value) of the dofs of the receiver's nodes (in nodal cs).
// Dofs containing expected unknowns (of expected type) are determined
// using this->GiveNodeDofIDMask function
{
    int k, nDofs;
    IntArray dofIDMask;
    FloatMatrix G2L;
    FloatArray vec;
    answer.resize( this->computeNumberOfGlobalDofs() );

    k = 0;
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        this->giveDofManDofIDMask(i, field.giveEquationID(), dofIDMask);
        this->giveDofManager(i)->giveUnknownVector(vec, dofIDMask, field, u, stepN);
        nDofs = vec.giveSize();
        for ( int j = 1; j <= nDofs; j++ ) {
            answer.at(++k) = vec.at(j);
        }
    }

    for ( int i = 1; i <= giveNumberOfInternalDofManagers(); i++ ) {
        this->giveInternalDofManDofIDMask(i, field.giveEquationID(), dofIDMask);
        this->giveInternalDofManager(i)->giveUnknownVector(vec, dofIDMask, field, u, stepN);
        nDofs = vec.giveSize();
        for ( int j = 1; j <= nDofs; j++ ) {
            answer.at(++k) = vec.at(j);
        }
    }
    answer.resizeWithValues(k);

    if (this->computeGtoLRotationMatrix(G2L)) {
        answer.rotatedWith(G2L, 'n');
    }
}


void
Element :: computeVectorOfPrescribed(EquationID ut, ValueModeType mode, TimeStep *stepN, FloatArray &answer)
// Forms the vector containing the prescribed values of the unknown 'u'
// (e.g., the prescribed displacement) of the dofs of the receiver's
// nodes. Puts 0 at each free dof.
{
    int k, nDofs;
    IntArray dofIDMask, dofMask;
    FloatMatrix G2L;
    FloatArray vec;

    answer.resize( this->computeNumberOfGlobalDofs() );

    k = 0;
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        this->giveDofManDofIDMask(i, ut, dofIDMask);
        this->giveDofManager(i)->givePrescribedUnknownVector(vec, dofIDMask, mode, stepN);
        nDofs = vec.giveSize();
        for ( int j = 1; j <= nDofs; j++ ) {
            answer.at(++k) = vec.at(j);
        }
    }

    for ( int i = 1; i <= giveNumberOfInternalDofManagers(); i++ ) {
        this->giveInternalDofManDofIDMask(i, ut, dofIDMask);
        this->giveInternalDofManager(i)->givePrescribedUnknownVector(vec, dofIDMask, mode, stepN);
        nDofs = vec.giveSize();
        for ( int j = 1; j <= nDofs; j++ ) {
            answer.at(++k) = vec.at(j);
        }
    }
    answer.resizeWithValues(k);

    if (this->computeGtoLRotationMatrix(G2L)) {
        answer.rotatedWith(G2L, 'n');
    }
}


int
Element :: computeNumberOfGlobalDofs()
{
    return this->computeNumberOfDofs();
}


int
Element :: computeNumberOfPrimaryMasterDofs(EquationID ut)
{
    int answer = 0;
    IntArray nodeDofIDMask, dofMask;

    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        this->giveDofManDofIDMask(i, ut, nodeDofIDMask);
        this->giveDofManager(i)->giveDofArray(nodeDofIDMask, dofMask);
        answer += this->giveDofManager(i)->giveNumberOfPrimaryMasterDofs(dofMask);
    }

    for ( int i = 1; i <= giveNumberOfInternalDofManagers(); i++ ) {
        this->giveInternalDofManDofIDMask(i, ut, nodeDofIDMask);
        this->giveInternalDofManager(i)->giveDofArray(nodeDofIDMask, dofMask);
        answer += this->giveInternalDofManager(i)->giveNumberOfPrimaryMasterDofs(dofMask);
    }
    return answer;
}


bool
Element :: giveRotationMatrix(FloatMatrix &answer, EquationID eid)
{
    bool is_GtoL, is_NtoG;
    FloatMatrix GtoL, NtoG;
    IntArray nodes;
    nodes.enumerate(this->giveNumberOfDofManagers());

    is_GtoL = this->computeGtoLRotationMatrix(GtoL);
    is_NtoG = this->computeDofTransformationMatrix(NtoG, nodes, true, eid);

#ifdef DEBUG
    if ( is_GtoL ) {
        if ( GtoL.giveNumberOfColumns() != this->computeNumberOfGlobalDofs() ) {
            _error("Element :: updateRotationMatrix - GtoL transformation matrix size mismatch in columns");
        }
        if ( GtoL.giveNumberOfRows() != this->computeNumberOfDofs() ) {
            _error("Element :: updateRotationMatrix - GtoL transformation matrix size mismatch in rows");
        }
    }
    if ( is_NtoG ) {
        if ( NtoG.giveNumberOfColumns() != this->computeNumberOfPrimaryMasterDofs(eid) ) {
            _error("Element :: updateRotationMatrix - NtoG transformation matrix size mismatch in columns");
        }
        if ( NtoG.giveNumberOfRows() != this->computeNumberOfGlobalDofs() ) {
            _error("Element :: updateRotationMatrix - NtoG transformation matrix size mismatch in rows");
        }
    }
#endif

    if ( is_GtoL && NtoG.isNotEmpty() ) {
        answer.beProductOf(GtoL, NtoG);
    } else if ( is_GtoL ) {
        answer = GtoL;
    } else if ( is_NtoG ) {
        answer = NtoG;
    } else {
        answer.beEmptyMtrx();
        return false;
    }
    return true;
}


bool
Element :: computeDofTransformationMatrix(FloatMatrix &answer, const IntArray &nodes, bool includeInternal, EquationID eid)
{
    bool flag = false;
    int numberOfDofMans = nodes.giveSize();

    // test if transformation is necessary
    for (int i = 1; i <= numberOfDofMans; i++ ) {
        flag = flag || this->giveDofManager(nodes.at(i))->requiresTransformation();
    }

    if ( !flag ) {
        answer.beEmptyMtrx();
        return false;
    }

    // initialize answer
    int gsize = this->computeNumberOfPrimaryMasterDofs(eid);
    answer.resize(this->computeNumberOfGlobalDofs(), gsize);
    answer.zero();

    FloatMatrix dofManT;
    IntArray dofIDmask;
    int nr, nc, lastRowPos = 0, lastColPos = 0;
    // loop over nodes
    for (int i = 1; i <= numberOfDofMans; i++ ) {
        this->giveDofManDofIDMask(nodes.at(i), eid, dofIDmask);
        if (!this->giveDofManager(nodes.at(i))->computeM2GTransformation(dofManT, dofIDmask)) {
            dofManT.resize(dofIDmask.giveSize(), dofIDmask.giveSize());
            dofManT.zero();
            dofManT.beUnitMatrix();
        }
        nc = dofManT.giveNumberOfColumns();
        nr = dofManT.giveNumberOfRows();
        for (int j = 1; j <= nr; j++ ) {
            for (int k = 1; k <= nc; k++ ) {
                // localize node contributions
                answer.at(lastRowPos + j, lastColPos + k) = dofManT.at(j, k);
            }
        }

        lastRowPos += nr;
        lastColPos += nc;
    }
    if ( includeInternal ) {
        for (int i = 1; i <= this->giveNumberOfInternalDofManagers(); i++ ) {
            this->giveInternalDofManDofIDMask(i, eid, dofIDmask);
            if (!this->giveInternalDofManager(nodes.at(i))->computeM2GTransformation(dofManT, dofIDmask)) {
                dofManT.resize(dofIDmask.giveSize(), dofIDmask.giveSize());
                dofManT.zero();
                dofManT.beUnitMatrix();
            }
            nc = dofManT.giveNumberOfColumns();
            nr = dofManT.giveNumberOfRows();
            for (int j = 1; j <= nr; j++ ) {
                for (int k = 1; k <= nc; k++ ) {
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


IntArray *
Element :: giveBodyLoadArray()
// Returns the array which contains the number of every body load that act
// on the receiver.
{
    return & bodyLoadArray;
}


IntArray *
Element :: giveBoundaryLoadArray()
// Returns the array which contains the number of every body load that act
// on the receiver.
{
    return & boundaryLoadArray;
}


void
Element :: giveLocationArray(IntArray &locationArray, EquationID eid, const UnknownNumberingScheme &s, IntArray *dofIdArray) const
// Returns the location array of the receiver. This array is obtained by
// simply appending the location array of every node of the receiver.
{
    IntArray masterDofIDs, nodalArray, dofIDMask;
    locationArray.resize(0);
    if (dofIdArray) dofIdArray->resize(0);
    for ( int i = 1; i <= this->numberOfDofMans; i++ ) {
        this->giveDofManDofIDMask(i, eid, dofIDMask);
        this->giveDofManager(i)->giveLocationArray(dofIDMask, nodalArray, s);
        locationArray.followedBy(nodalArray);
        if (dofIdArray) {
            this->giveDofManager(i)->giveMasterDofIDArray(dofIDMask, masterDofIDs);
            dofIdArray->followedBy(masterDofIDs);
        }
    }
    for ( int i = 1; i <= this->giveNumberOfInternalDofManagers(); i++ ) {
        this->giveInternalDofManDofIDMask(i, eid, dofIDMask);
        this->giveInternalDofManager(i)->giveLocationArray(dofIDMask, nodalArray, s);
        locationArray.followedBy(nodalArray);
        if (dofIdArray) {
            this->giveInternalDofManager(i)->giveMasterDofIDArray(dofIDMask, masterDofIDs);
            dofIdArray->followedBy(masterDofIDs);
        }
    }
}


void
Element :: giveLocationArray(IntArray &locationArray, const IntArray &dofIDMask, const UnknownNumberingScheme &s, IntArray *dofIdArray) const
{
    IntArray masterDofIDs, nodalArray, ids = dofIDMask;
    locationArray.resize(0);
    if (dofIdArray) dofIdArray->resize(0);
    for ( int i = 1; i <= this->numberOfDofMans; i++ ) {
        if ( dofIDMask.giveSize() == 0 ) {
            this->giveDefaultDofManDofIDMask(i, ids);
        }
        this->giveDofManager(i)->giveLocationArray(ids, nodalArray, s);
        locationArray.followedBy(nodalArray);
        if (dofIdArray) {
            this->giveDofManager(i)->giveMasterDofIDArray(ids, masterDofIDs);
            dofIdArray->followedBy(masterDofIDs);
        }
    }
    for ( int i = 1; i <= this->giveNumberOfInternalDofManagers(); i++ ) {
        if ( dofIDMask.giveSize() == 0 ) {
            this->giveDefaultInternalDofManDofIDMask(i, ids);
        }
        this->giveInternalDofManager(i)->giveLocationArray(ids, nodalArray, s);
        locationArray.followedBy(nodalArray);
        if (dofIdArray) {
            this->giveInternalDofManager(i)->giveMasterDofIDArray(ids, masterDofIDs);
            dofIdArray->followedBy(masterDofIDs);
        }
    }
}


void
Element :: giveBoundaryLocationArray(IntArray &locationArray, const IntArray &bNodes, EquationID eid, const UnknownNumberingScheme &s, IntArray *dofIdArray)
{
    IntArray masterDofIDs, nodalArray, dofIDMask;
    locationArray.resize(0);
    if (dofIdArray) dofIdArray->resize(0);
    for ( int i = 1; i <= bNodes.giveSize(); i++ ) {
        this->giveDofManDofIDMask(bNodes.at(i), eid, dofIDMask);
        this->giveDofManager(bNodes.at(i))->giveLocationArray(dofIDMask, nodalArray, s);
        locationArray.followedBy(nodalArray);
        if (dofIdArray) {
            this->giveDofManager(bNodes.at(i))->giveMasterDofIDArray(dofIDMask, masterDofIDs);
            dofIdArray->followedBy(masterDofIDs);
        }
    }
}


void
Element :: giveBoundaryLocationArray(IntArray &locationArray, const IntArray &bNodes, const IntArray &dofIDMask, const UnknownNumberingScheme &s, IntArray *dofIdArray)
{
    IntArray masterDofIDs, nodalArray, ids = dofIDMask;
    locationArray.resize(0);
    if (dofIdArray) dofIdArray->resize(0);
    for ( int i = 1; i <= bNodes.giveSize(); i++ ) {
        if ( dofIDMask.giveSize() == 0 ) {
            this->giveDefaultDofManDofIDMask(bNodes.at(i), ids);
        }
        this->giveDofManager(bNodes.at(i))->giveLocationArray(ids, nodalArray, s);
        locationArray.followedBy(nodalArray);
        if (dofIdArray) {
            this->giveDofManager(bNodes.at(i))->giveMasterDofIDArray(ids, masterDofIDs);
            dofIdArray->followedBy(masterDofIDs);
        }
    }
}


Material *Element :: giveMaterial()
// Returns the material of the receiver.
{
#ifdef DEBUG
    if ( !material ) {
        // material = this -> readInteger("mat") ;
        _error("giveMaterial: material not defined");
    }
#endif
    return domain->giveMaterial(material);
}


CrossSection *Element :: giveCrossSection()
// Returns the crossSection of the receiver.
{
#ifdef DEBUG
    if ( !crossSection ) {
        _error("giveCrossSection: crossSection not defined");
    }
#endif
    return domain->giveCrossSection(crossSection);
}


int
Element :: giveRegionNumber()
{
    return this->giveCrossSection()->giveNumber();
}


DofManager *
Element :: giveDofManager(int i) const
// Returns the i-th node of the receiver.
{
#ifdef DEBUG
    if ( ( i <= 0 ) || ( i > dofManArray.giveSize() ) ) {
        OOFEM_ERROR2("giveNode: Node %i is not defined", i);
    }
#endif
    return domain->giveDofManager(dofManArray.at(i));
}


Node *
Element :: giveNode(int i) const
// Returns the i-th node of the receiver.
{
#ifdef DEBUG
    if ( ( i <= 0 ) || ( i > dofManArray.giveSize() ) ) {
        _error("giveNode: Node is not defined");
    }
#endif
    return domain->giveNode(dofManArray.at(i));
}


ElementSide *
Element :: giveSide(int i) const
// Returns the i-th side of the receiver.
{
#ifdef DEBUG
    if ( ( i <= 0 ) || ( i > dofManArray.giveSize() ) ) {
        _error("giveNode: Side is not defined");
    }
#endif
    return domain->giveSide(dofManArray.at(i));
}


void
Element :: setDofManagers(const IntArray &_dmans)
{
    this->dofManArray = _dmans;
}

void
Element :: setBodyLoads(const IntArray &_bodyLoads)
{
    this->bodyLoadArray = _bodyLoads;
}

void
Element :: setIntegrationRules(AList< IntegrationRule > *irlist)
{
    if ( integrationRulesArray ) {
        for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
            delete integrationRulesArray [ i ];
        }

        delete[] integrationRulesArray;
    }

    numberOfIntegrationRules = irlist->giveSize();
    integrationRulesArray = new IntegrationRule * [ irlist->giveSize() ];

    for ( int j = 0; j < irlist->giveSize(); j++ ) {
        integrationRulesArray [ j ] =  irlist->at(j + 1);
        irlist->unlink(j + 1);
    }
}


void
Element :: giveCharacteristicMatrix(FloatMatrix &answer,
                                     CharType mtrx, TimeStep *tStep)
//
// returns characteristics matrix of receiver according to mtrx
//
{
    _error("giveCharacteristicMatrix: Unknown Type of characteristic mtrx.");
}


void
Element :: giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep)
//
// returns characteristics vector of receiver according to mtrx
//
{
    _error("giveCharacteristicVector: Unknown Type of characteristic mtrx.");
}


void
Element :: computeLoadVector(FloatArray &answer, Load *load, CharType type, ValueModeType mode, TimeStep *tStep)
{
    _error("computeLoadVector: Unknown load type.");
}


void
Element :: computeBoundaryLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep)
{
    _error("computeBoundaryLoadVector: Unknown load type.");
}


void
Element :: computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int edge, CharType type, ValueModeType mode, TimeStep *tStep)
{
    ///@todo Change the load type to "BoundaryEdgeLoad" maybe?
    _error("computeBoundaryEdgeLoadVector: Unknown load type.");
}


double
Element :: giveCharacteristicValue(CharType mtrx, TimeStep *tStep)
//
// returns characteristics value of receiver according to CharType
//
{
    _error("giveCharacteristicValue: Unknown Type of characteristic mtrx.");
    return 0.;
}


IRResultType
Element :: initializeFrom(InputRecord *ir)
{
    const char *__proc = "initializeFrom"; // Required by IR_GIVE_FIELD macro
    IRResultType result;                          // Required by IR_GIVE_FIELD macro

#  ifdef VERBOSE
    // VERBOSE_PRINT1("Instanciating element ",number);
#  endif
    //IR_GIVE_FIELD(ir, material, _IFT_Element_mat);
    IR_GIVE_OPTIONAL_FIELD(ir, material, _IFT_Element_mat);

    //IR_GIVE_FIELD(ir, crossSection, _IFT_Element_crosssect);
    IR_GIVE_OPTIONAL_FIELD(ir, crossSection, _IFT_Element_crosssect);

    IR_GIVE_FIELD(ir, dofManArray, _IFT_Element_nodes);

    bodyLoadArray.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, bodyLoadArray, _IFT_Element_bodyload);

    boundaryLoadArray.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, boundaryLoadArray, _IFT_Element_boundaryload);

    elemLocalCS.resize(0, 0);

    if ( ir->hasField(_IFT_Element_lcs) ) { //local coordinate system
        double n1 = 0.0, n2 = 0.0;
        FloatArray triplets;
        triplets.resize(0);
        IR_GIVE_OPTIONAL_FIELD(ir, triplets, _IFT_Element_lcs);
        elemLocalCS.resize(3, 3);
        for ( int j = 1; j <= 3; j++ ) {
            elemLocalCS.at(j, 1) = triplets.at(j);
            n1 += triplets.at(j) * triplets.at(j);
            elemLocalCS.at(j, 2) = triplets.at(j + 3);
            n2 += triplets.at(j + 3) * triplets.at(j + 3);
        }

        n1 = sqrt(n1);
        n2 = sqrt(n2);
        for ( int j = 1; j <= 3; j++ ) { // normalize e1' e2'
            elemLocalCS.at(j, 1) /= n1;
            elemLocalCS.at(j, 2) /= n2;
        }

        // vector e3' computed from vector product of e1', e2'
        elemLocalCS.at(1, 3) = ( elemLocalCS.at(2, 1) * elemLocalCS.at(3, 2) - elemLocalCS.at(3, 1) * elemLocalCS.at(2, 2) );
        elemLocalCS.at(2, 3) = ( elemLocalCS.at(3, 1) * elemLocalCS.at(1, 2) - elemLocalCS.at(1, 1) * elemLocalCS.at(3, 2) );
        elemLocalCS.at(3, 3) = ( elemLocalCS.at(1, 1) * elemLocalCS.at(2, 2) - elemLocalCS.at(2, 1) * elemLocalCS.at(1, 2) );
        //elemLocalCS.printYourself();
    }

#ifdef __PARALLEL_MODE
    partitions.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, partitions, _IFT_Element_partitions);
    if ( ir->hasField(_IFT_Element_remote) ) {
        parallel_mode = Element_remote;
    } else {
        parallel_mode = Element_local;
    }

#endif

    IR_GIVE_OPTIONAL_FIELD(ir, activityLtf, _IFT_Element_activityltf);

    IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, _IFT_Element_nip);

    return IRRT_OK;
}


void
Element :: giveInputRecord(DynamicInputRecord &input)
{
    FEMComponent :: giveInputRecord(input);

    input.setField(material, _IFT_Element_mat);

    input.setField(crossSection, _IFT_Element_crosssect);

    input.setField(dofManArray, _IFT_Element_nodes);

    if(bodyLoadArray.giveSize() > 0) {
        input.setField(bodyLoadArray, _IFT_Element_bodyload);
    }


    if(boundaryLoadArray.giveSize() > 0) {
        input.setField(boundaryLoadArray, _IFT_Element_boundaryload);
    }


    if( elemLocalCS.giveNumberOfRows() > 0 ) {
        FloatArray triplets(6);
        for ( int j = 1; j <= 3; j++ ) {
            triplets.at(j) = elemLocalCS.at(j, 1);
            triplets.at(j + 3) = elemLocalCS.at(j, 2);
        }
        input.setField(triplets, _IFT_Element_lcs);
    }


#ifdef __PARALLEL_MODE
    if ( this->giveDomain()->giveEngngModel()->isParallel() && partitions.giveSize() > 0 ) {
        input.setField(this->partitions, _IFT_Element_partitions);
        if ( this->parallel_mode == Element_remote ) {
            input.setField(_IFT_Element_remote);
        }
    }
#endif

    if( activityLtf > 0 ) {
        input.setField(activityLtf, _IFT_Element_activityltf);
    }

    input.setField(numberOfGaussPoints, _IFT_Element_nip);
}


void
Element :: postInitialize()
{
    this->computeGaussPoints();
}


void
Element :: printOutputAt(FILE *file, TimeStep *stepN)
// Performs end-of-step operations.
{
    fprintf( file, "element %d (%8d) :\n", this->giveLabel(), this->giveNumber() );

    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        integrationRulesArray [ i ]->printOutputAt(file, stepN);
    }
}


void
Element :: updateYourself(TimeStep *tStep)
// Updates the receiver at end of step.
{
#  ifdef VERBOSE
    // VERBOSE_PRINT1("Updating Element ",number)
#  endif

    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        integrationRulesArray [ i ]->updateYourself(tStep);
    }

}


bool 
Element :: isActivated(TimeStep *tStep)
{
    if ( activityLtf ) {
        if ( tStep ) {
            return ( domain->giveLoadTimeFunction(activityLtf)->evaluate(tStep, VM_Total) > 1.e-3 );
        } else {
            return false;
        }
    } else {
        return true;
    }
}

void
Element :: initForNewStep()
// initializes receiver to new time step or can be used
// if current time step must be restarted
{
    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        integrationRulesArray [ i ]->initForNewStep();
    }
}


contextIOResultType Element :: saveContext(DataStream *stream, ContextMode mode, void *obj)
//
// saves full element context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;
    int _val;

    if ( ( iores = FEMComponent :: saveContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( ( mode & CM_Definition ) ) {
        if ( !stream->write(& numberOfDofMans, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream->write(& material, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream->write(& crossSection, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

#ifdef __PARALLEL_MODE
        if ( mode & CM_DefinitionGlobal ) {
            // send global numbers instead of local ones
            int s = dofManArray.giveSize();
            IntArray globDN(s);
            for ( int i = 1; i <= s; i++ ) {
                globDN.at(i) = this->giveDofManager(i)->giveGlobalNumber();
            }

            if ( ( iores = globDN.storeYourself(stream, mode) ) != CIO_OK ) {
                THROW_CIOERR(iores);
            }
        } else {
            if ( ( iores = dofManArray.storeYourself(stream, mode) ) != CIO_OK ) {
                THROW_CIOERR(iores);
            }
        }

#else
        if ( ( iores = dofManArray.storeYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

#endif
        if ( ( iores = bodyLoadArray.storeYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( ( iores = boundaryLoadArray.storeYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( !stream->write(& numberOfIntegrationRules, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
            _val = integrationRulesArray [ i ]->giveIntegrationRuleType();
            if ( !stream->write(& _val, 1) ) {
                THROW_CIOERR(CIO_IOERR);
            }
        }

#ifdef __PARALLEL_MODE
        int _mode;
        if ( !stream->write(& globalNumber, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        _mode = parallel_mode;
        if ( !stream->write(& _mode, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( ( iores = partitions.storeYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

#endif
    }

    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        if ( ( iores = integrationRulesArray [ i ]->saveContext(stream, mode, obj) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}


contextIOResultType Element :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
//
// restores full element context (saves state variables, that completely describe
// current state)
//
{
    contextIOResultType iores;
    int _nrules;

    if ( ( iores = FEMComponent :: restoreContext(stream, mode, obj) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    if ( mode & CM_Definition ) {
        if ( !stream->read(& numberOfDofMans, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream->read(& material, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream->read(& crossSection, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( ( iores = dofManArray.restoreYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( ( iores = bodyLoadArray.restoreYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( ( iores = boundaryLoadArray.restoreYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( !stream->read(& _nrules, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        // restore integration rules
        IntArray dtypes(_nrules);
        for ( int i = 1; i <= _nrules; i++ ) {
            if ( !stream->read(& dtypes.at(i), 1) ) {
                THROW_CIOERR(CIO_IOERR);
            }
        }

        if ( _nrules != numberOfIntegrationRules ) {
            // delete old int rule array
            if ( integrationRulesArray ) {
                for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
                    delete integrationRulesArray [ i ];
                }

                delete[] integrationRulesArray;
            }

            // AND ALLOCATE NEW ONE
            integrationRulesArray = new IntegrationRule * [ _nrules ];
            for ( int i = 0; i < _nrules; i++ ) {
                integrationRulesArray [ i ] = classFactory.createIRule( ( IntegrationRuleType ) dtypes(i), i + 1, this );
            }

            numberOfIntegrationRules = _nrules;
        } else {
            for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
                if ( integrationRulesArray [ i ]->giveIntegrationRuleType() != dtypes(i) ) {
                    delete integrationRulesArray [ i ];
                    integrationRulesArray [ i ] = classFactory.createIRule( ( IntegrationRuleType ) dtypes(i), i + 1, this );
                }
            }
        }

#ifdef __PARALLEL_MODE
        int _mode;
        if ( !stream->read(& globalNumber, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream->read(& _mode, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        parallel_mode = ( elementParallelMode ) _mode;
        if ( ( iores = partitions.restoreYourself(stream, mode) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

#endif
    }


    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        if ( ( iores = integrationRulesArray [ i ]->restoreContext(stream, mode, this) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    return CIO_OK;
}


double
Element :: computeVolumeAreaOrLength()
// the element computes its volume, area or length
// (depending on the spatial dimension of that element)
{
    GaussPoint *gp;
    double answer = 0.;
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    if ( iRule ) {
        for ( int i = 0; i < iRule->giveNumberOfIntegrationPoints(); i++ ) {
            gp  = iRule->getIntegrationPoint(i);
            answer += this->computeVolumeAround(gp);
        }

        return answer;
    }

    return -1.; // means "cannot be evaluated"
}


double
Element :: computeMeanSize()
// Computes the size of the element defined as its length,
// square root of area or cube root of volume (depending on spatial dimension)
{
    // if the method "giveArea" is properly implemented
    // for the particular type of element, the mean size is the square root of area
    // 8 July 2010 - does not seem to work any more (Element does not inherit from ElementGeometry)
    // double area = this->giveArea();
    // if (area>0.)
    //  return sqrt(area);

    // if "giveArea" is not implemented (default value 0.),
    // then the contributing areas or volumes are collected from Gauss points
    double volume = this->computeVolumeAreaOrLength();
    if ( volume < 0. ) {
        return -1.; // means "cannot be evaluated"
    }

    int dim = this->giveSpatialDimension();
    switch ( dim ) {
    case 1: return volume;

    case 2: return sqrt(volume);

    case 3: return pow(volume, 1./3.);
    }

    return -1.; // means "cannot be evaluated"
}


double
Element :: computeVolume()
{
    FEInterpolation3d *fei = dynamic_cast<FEInterpolation3d*>(this->giveInterpolation());
#ifdef DEBUG
    if (!fei) {
        OOFEM_ERROR("Element :: computeVolume - Function not overloaded and necessary interpolator isn't available");
        return 0.0;
    }
#endif
    return fei->giveVolume(FEIElementGeometryWrapper(this));
}


double
Element :: computeArea()
{
    FEInterpolation2d *fei = dynamic_cast<FEInterpolation2d*>(this->giveInterpolation());
#ifdef DEBUG
    if (!fei) {
        OOFEM_ERROR("Element :: computeArea - Function not overloaded and necessary interpolator isn't available");
        return 0.0;
    }
#endif
    return fei->giveArea(FEIElementGeometryWrapper(this));
}


double
Element :: computeLength()
{
    FEInterpolation1d *fei = dynamic_cast<FEInterpolation1d*>(this->giveInterpolation());
#ifdef DEBUG
    if (!fei) {
        OOFEM_ERROR("Element :: computeLength - Function not overloaded and necessary interpolator isn't available");
        return 0.0;
    }
#endif
    return fei->giveLength(FEIElementGeometryWrapper(this));
}


double
Element :: giveLenghtInDir(const FloatArray &normalToCrackPlane)
//
// returns receivers projection length (for some material models)
// to direction given by normalToCrackPlane;
//
{
    FloatArray *coords;
    double maxDis, minDis, dis;
    int nnode = giveNumberOfNodes();

    coords = this->giveNode(1)->giveCoordinates();
    minDis = maxDis = normalToCrackPlane.dotProduct(*coords, coords->giveSize());

    for (int i = 2; i <= nnode; i++ ) {
        coords = this->giveNode(i)->giveCoordinates();
        dis = normalToCrackPlane.dotProduct(*coords, coords->giveSize());
        if ( dis > maxDis ) {
            maxDis = dis;
        } else if ( dis < minDis ) {
            minDis = dis;
        }
    }

    return maxDis - minDis;
}


int
Element :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    FEInterpolation *fei = this->giveInterpolation();
#ifdef DEBUG
    if (!fei) {
        answer.resize(0);
        return false;
    }
#endif
    fei->local2global(answer, lcoords, FEIElementGeometryWrapper(this));
    return true;
}


bool
Element :: computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords)
{
    FEInterpolation *fei = this->giveInterpolation();
    if (fei) {
        return fei->global2local(answer, gcoords, FEIElementGeometryWrapper(this));
    } else {
        return false;
    }
}


int
Element :: giveLocalCoordinateSystem(FloatMatrix &answer)
{
    if ( elemLocalCS.isNotEmpty() ) {
        answer = elemLocalCS;
        return 1;
    } else   {
        answer.beEmptyMtrx();
    }

    return 0;
}


void
Element :: computeMidPlaneNormal(FloatArray &answer, const GaussPoint *)
// valid only for plane elements (shells, plates, ....)
// computes mid-plane normal at gaussPoint - for materials with orthotrophy
{
    _error("Unable to compute mid-plane normal, not supported");
}


int
Element :: giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime)
{
    if ( type == IST_ErrorIndicatorLevel ) {
        ErrorEstimator *ee = this->giveDomain()->giveErrorEstimator();
        if ( ee ) {
            answer.resize(1);
            answer.at(1) = ee->giveElementError(indicatorET, this, atTime);
        } else {
            answer.resize(0);
            return 0;
        }

        return 1;
    } else if ( type == IST_InternalStressError ) {
        ErrorEstimator *ee = this->giveDomain()->giveErrorEstimator();
        if ( ee ) {
            answer.resize(1);
            answer.at(1) = ee->giveElementError(internalStressET, this, atTime);
        } else {
            answer.resize(0);
            return 0;
        }

        return 1;
    } else if ( type == IST_PrimaryUnknownError ) {
        ErrorEstimator *ee = this->giveDomain()->giveErrorEstimator();
        if ( ee ) {
            answer.resize(1);
            answer.at(1) = ee->giveElementError(primaryUnknownET, this, atTime);
        } else {
            answer.resize(0);
            return 0;
        }

        return 1;
    } else {
        return this->giveCrossSection()->giveIPValue(answer, aGaussPoint, type, atTime);
    }
}

int
Element :: giveSpatialDimension()
{
    ///@todo Just ask the interpolator instead?
    switch ( this->giveGeometryType() ) {
    case EGT_point:
        return 0;

    case EGT_line_1:
    case EGT_line_2:
        return 1;

    case EGT_triangle_1:
    case EGT_triangle_2:
    case EGT_quad_1:
    case EGT_quad_2:
    case EGT_quad9_2:
    case EGT_quad_1_interface:
    case EGT_quad_21_interface:
        return 2;

    case EGT_tetra_1:
    case EGT_tetra_2:
    case EGT_hexa_1:
    case EGT_hexa_2:
    case EGT_hexa_27:
    case EGT_wedge_1:
    case EGT_wedge_2:
        return 3;

    case EGT_Composite:
    case EGT_unknown:
        break;
    }

    _error("giveSpatialDimension: failure (maybe new element type was registered)");
    return 0; //to make compiler happy
}


int
Element :: giveNumberOfBoundarySides()
{
    ///@todo Just ask the interpolator instead?
    switch ( this->giveGeometryType() ) {
    case EGT_point:
        return 0;

    case EGT_line_1:
    case EGT_line_2:
    case EGT_quad_1_interface:
    case EGT_quad_21_interface:
        return 2;

    case EGT_triangle_1:
    case EGT_triangle_2:
        return 3;

    case EGT_quad_1:
    case EGT_quad_2:
    case EGT_quad9_2:
        return 4;

    case EGT_tetra_1:
    case EGT_tetra_2:
        return 4;

    case EGT_wedge_1:
    case EGT_wedge_2:
        return 5;

    case EGT_hexa_1:
    case EGT_hexa_2:
    case EGT_hexa_27:
        return 6;

    case EGT_Composite:
    case EGT_unknown:
        break;
    }

    _error2( "giveSpatialDimension: failure, unsupported geometry type (%s)",
            __Element_Geometry_TypeToString( this->giveGeometryType() ) );
    return 0; // to make compiler happy
}


int
Element :: adaptiveMap(Domain *oldd, TimeStep *tStep)
{
    int result = 1;
    IntegrationRule *iRule;
    MaterialModelMapperInterface *interface = static_cast< MaterialModelMapperInterface * >
                                              ( this->giveMaterial()->giveInterface(MaterialModelMapperInterfaceType) );

    if ( !interface ) {
        return 0;
    }

    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        iRule = integrationRulesArray [ i ];
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            result &= interface->MMI_map(iRule->getIntegrationPoint(j), oldd, tStep);
        }
    }

    return result;
}

int
Element :: mapStateVariables(const Domain &iOldDom, const TimeStep &iTStep)
{
    int result = 1;

    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
    	IntegrationRule *iRule = integrationRulesArray [ i ];
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {


        	GaussPoint &gp = *(iRule->getIntegrationPoint(j));

        	MaterialStatus *ms = dynamic_cast<MaterialStatus*>(gp.giveMaterialStatus() );
        	if(ms == NULL) {
        		OOFEM_ERROR("In Element :: mapStateVariables(): failed to fetch MaterialStatus.\n");
        	}

            MaterialStatusMapperInterface *interface = dynamic_cast< MaterialStatusMapperInterface * > ( ms );
            if ( interface == NULL ) {
        		OOFEM_ERROR("In Element :: mapStateVariables(): Failed to fetch MaterialStatusMapperInterface.\n");
            }

            result &= interface->MSMI_map(gp, iOldDom, iTStep, *( ms ) );
        }
    }

    return result;
}


int
Element :: adaptiveFinish(TimeStep *tStep)
{
    int result = 1;
    IntegrationRule *iRule;
    MaterialModelMapperInterface *interface = static_cast< MaterialModelMapperInterface * >
                                              ( this->giveMaterial()->giveInterface(MaterialModelMapperInterfaceType) );

    if ( !interface ) {
        return 0;
    }

    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        iRule = integrationRulesArray [ i ];
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            result &= interface->MMI_finish(tStep);
        }
    }

    return result;
}


void
Element :: updateLocalNumbering(EntityRenumberingFunctor &f)
{
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
        dofManArray.at(i) = f(dofManArray.at(i), ERS_DofManager);
    }
}


integrationDomain
Element :: giveIntegrationDomain() const
{
    FEInterpolation *fei = this->giveInterpolation();
    return fei ? fei->giveIntegrationDomain() : _Unknown_integrationDomain;
}


Element_Geometry_Type
Element :: giveGeometryType() const
{
    FEInterpolation *fei = this->giveInterpolation();
    return fei ? fei->giveGeometryType() : EGT_unknown;
}


bool
Element :: computeGtoLRotationMatrix(FloatMatrix &answer)
{
    answer.beEmptyMtrx();
    return false;
}


#ifdef __PARALLEL_MODE
int
Element :: packUnknowns(CommunicationBuffer &buff, TimeStep *stepN)
{
    int result = 1;
    IntegrationRule *iRule;

    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        iRule = integrationRulesArray [ i ];
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            result &= this->giveCrossSection()->packUnknowns( buff, stepN, iRule->getIntegrationPoint(j) );
        }
    }

    return result;
}


int
Element :: unpackAndUpdateUnknowns(CommunicationBuffer &buff, TimeStep *stepN)
{
    int result = 1;
    IntegrationRule *iRule;

    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        iRule = integrationRulesArray [ i ];
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            result &= this->giveCrossSection()->unpackAndUpdateUnknowns( buff, stepN, iRule->getIntegrationPoint(j) );
        }
    }

    return result;
}


int
Element :: estimatePackSize(CommunicationBuffer &buff)
{
    int result = 0;
    IntegrationRule *iRule;

    for ( int i = 0; i < numberOfIntegrationRules; i++ ) {
        iRule = integrationRulesArray [ i ];
        for ( int j = 0; j < iRule->giveNumberOfIntegrationPoints(); j++ ) {
            result += this->giveCrossSection()->estimatePackSize( buff, iRule->getIntegrationPoint(j) );
        }
    }

    return result;
}


double
Element :: predictRelativeComputationalCost()
{
    int nip;
    double wgt = 0;
    IntegrationRule *iRule = this->giveDefaultIntegrationRulePtr();
    nip = iRule->giveNumberOfIntegrationPoints();
    for ( int j = 0; j < nip; j++ ) {
        wgt += this->giveCrossSection()->predictRelativeComputationalCost( iRule->getIntegrationPoint(j) );
    }

    return ( this->giveRelativeSelfComputationalCost() * wgt );
}
#endif


#ifdef __OOFEG
void
Element :: drawYourself(oofegGraphicContext &gc)
{
    OGC_PlotModeType mode = gc.giveIntVarPlotMode();

    if ( mode == OGC_rawGeometry ) {
        this->drawRawGeometry(gc);
    } else if ( mode == OGC_elementAnnotation ) {
        this->drawAnnotation(gc);
    } else if ( mode == OGC_deformedGeometry ) {
        this->drawDeformedGeometry(gc, DisplacementVector);
    } else if ( mode == OGC_eigenVectorGeometry ) {
        this->drawDeformedGeometry(gc, EigenVector);
    } else if ( mode == OGC_scalarPlot ) {
        this->drawScalar(gc);
    } else if ( mode == OGC_elemSpecial ) {
        this->drawSpecial(gc);
    } else {
        _error("drawYourself : unsupported mode");
    }
}


void
Element :: drawAnnotation(oofegGraphicContext &gc)
{
    int i, count = 0;
    Node *node;
    WCRec p [ 1 ]; /* point */
    GraphicObj *go;
    char num [ 30 ];

    p [ 0 ].x = p [ 0 ].y = p [ 0 ].z = 0.0;
    // compute element center
    for ( i = 1; i <= numberOfDofMans; i++ ) {
        if ( ( node = this->giveNode(i) ) ) {
            p [ 0 ].x += node->giveCoordinate(1);
            p [ 0 ].y += node->giveCoordinate(2);
            p [ 0 ].z += node->giveCoordinate(3);
            count++;
        }
    }

    p [ 0 ].x /= count;
    p [ 0 ].y /= count;
    p [ 0 ].z /= count;

    EASValsSetLayer(OOFEG_ELEMENT_ANNOTATION_LAYER);
    EASValsSetColor( gc.getElementColor() );
 #ifdef __PARALLEL_MODE
    sprintf( num, "%d(%d)", this->giveNumber(), this->giveGlobalNumber() );
 #else
    sprintf( num, "%d", this->giveNumber() );
 #endif
    go = CreateAnnText3D(p, num);
    EGWithMaskChangeAttributes(COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}


int
Element :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                   int node, TimeStep *atTime)
{
    if ( type == IST_RelMeshDensity ) {
        ErrorEstimator *ee = this->giveDomain()->giveErrorEstimator();
        if ( ee ) {
            answer.resize(1);
            answer.at(1) = this->giveDomain()->giveErrorEstimator()->giveRemeshingCrit()->
                           giveRequiredDofManDensity(this->giveNode(node)->giveNumber(), atTime, 1);
            return 1;
        } else {
            answer.resize(0);
            return 0;
        }
    } else {
        if ( mode == ISM_recovered ) {
            const FloatArray *nodval;
            NodalRecoveryModel* smoother = this->giveDomain()->giveSmoother();
            int result = smoother->giveNodalVector( nodval, this->giveNode(node)->giveNumber(),
                                    smoother->giveElementVirtualRegionNumber(this->number) );
            if ( nodval ) {
                answer = * nodval;
            } else {
                answer.resize(0);
            }

            return result;
        } else {
            return 0;
        }
    }
}


#endif
} // end namespace oofem
