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
#include "function.h"
#include "dofmanager.h"
#include "node.h"
#include "gausspoint.h"
#include "unknownnumberingscheme.h"
#include "dynamicinputrecord.h"
#include "matstatmapperint.h"
#include "cltypes.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
#endif

#include <cstdio>

namespace oofem {
Element :: Element(int n, Domain *aDomain) :
    FEMComponent(n, aDomain), dofManArray(), crossSection(0), bodyLoadArray(), boundaryLoadArray(), integrationRulesArray()
{
    material           = 0;
    numberOfDofMans    = 0;
    activityTimeFunction = 0;
}


Element :: ~Element()
{
}


void
Element :: computeVectorOf(ValueModeType u, TimeStep *tStep, FloatArray &answer)
{
    IntArray dofIDMask;
    FloatMatrix G2L;
    FloatArray vec;

    answer.reserve( this->computeNumberOfGlobalDofs() );

    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        this->giveDofManDofIDMask(i, dofIDMask);
        this->giveDofManager(i)->giveUnknownVector(vec, dofIDMask, u, tStep, true);
        answer.append(vec);
    }

    for ( int i = 1; i <= giveNumberOfInternalDofManagers(); i++ ) {
        this->giveInternalDofManDofIDMask(i, dofIDMask);
        this->giveInternalDofManager(i)->giveUnknownVector(vec, dofIDMask, u, tStep, true);
        answer.append(vec);
    }

    if ( this->computeGtoLRotationMatrix(G2L) ) {
        answer.rotatedWith(G2L, 'n');
    }
}


void
Element :: computeVectorOf(const IntArray &dofIDMask, ValueModeType u, TimeStep *tStep, FloatArray &answer, bool padding)
{
    FloatMatrix G2L;
    FloatArray vec;

    answer.reserve( dofIDMask.giveSize() * ( this->giveNumberOfDofManagers() + this->giveNumberOfInternalDofManagers() ) );

    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        this->giveDofManager(i)->giveUnknownVector(vec, dofIDMask, u, tStep, padding);
        answer.append(vec);
    }

    for ( int i = 1; i <= giveNumberOfInternalDofManagers(); i++ ) {
        this->giveInternalDofManager(i)->giveUnknownVector(vec, dofIDMask, u, tStep, padding);
        answer.append(vec);
    }

    ///@todo This rotation matrix needs to have the dofidmask/eid/primary field or something passed to it, otherwise it won't work generally.
    if ( this->computeGtoLRotationMatrix(G2L) ) {
        OOFEM_WARNING("The transformation matrix from global -> element local c.s. is not fully supported for this function (yet)");
        answer.rotatedWith(G2L, 'n');
    }
}


void
Element :: computeBoundaryVectorOf(const IntArray &bNodes, const IntArray &dofIDMask, ValueModeType u, TimeStep *tStep, FloatArray &answer, bool padding)
{
    FloatMatrix G2L;
    FloatArray vec;

    answer.reserve( dofIDMask.giveSize() * bNodes.giveSize() );

    for ( int bNode: bNodes ) {
        this->giveDofManager( bNode )->giveUnknownVector(vec, dofIDMask, u, tStep, padding);
        answer.append(vec);
    }

    if ( this->computeGtoLRotationMatrix(G2L) ) {
        OOFEM_ERROR("Local coordinate system is not implemented yet");
    }
}


void
Element :: computeVectorOf(PrimaryField &field, const IntArray &dofIDMask, ValueModeType u, TimeStep *tStep, FloatArray &answer, bool padding)
{
    FloatMatrix G2L;
    FloatArray vec;
    answer.reserve( this->computeNumberOfGlobalDofs() );

    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        this->giveDofManager(i)->giveUnknownVector(vec, dofIDMask, field, u, tStep);
        answer.append(vec);
    }

    for ( int i = 1; i <= giveNumberOfInternalDofManagers(); i++ ) {
        this->giveInternalDofManager(i)->giveUnknownVector(vec, dofIDMask, field, u, tStep);
        answer.append(vec);
    }

    if ( this->computeGtoLRotationMatrix(G2L) ) {
        answer.rotatedWith(G2L, 'n');
    }
}


void
Element :: computeVectorOfPrescribed(ValueModeType u, TimeStep *tStep, FloatArray &answer)
{
    IntArray dofIDMask;
    FloatMatrix G2L;
    FloatArray vec;

    answer.reserve( this->computeNumberOfGlobalDofs() );

    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        this->giveDofManDofIDMask(i, dofIDMask);
        this->giveDofManager(i)->givePrescribedUnknownVector(vec, dofIDMask, u, tStep);
        answer.append(vec);
    }

    for ( int i = 1; i <= giveNumberOfInternalDofManagers(); i++ ) {
        this->giveInternalDofManDofIDMask(i, dofIDMask);
	this->giveInternalDofManager(i)->givePrescribedUnknownVector(vec, dofIDMask, u, tStep);
        answer.append(vec);
    }

    if ( this->computeGtoLRotationMatrix(G2L) ) {
        answer.rotatedWith(G2L, 'n');
    }

}

  
void
Element :: computeVectorOfPrescribed(const IntArray &dofIDMask, ValueModeType mode, TimeStep *tStep, FloatArray &answer)
{
    FloatMatrix G2L;
    FloatArray vec;

    answer.reserve( dofIDMask.giveSize() * this->computeNumberOfGlobalDofs() );

    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        this->giveDofManager(i)->givePrescribedUnknownVector(vec, dofIDMask, mode, tStep);
        answer.append(vec);
    }

    for ( int i = 1; i <= giveNumberOfInternalDofManagers(); i++ ) {
        this->giveInternalDofManager(i)->givePrescribedUnknownVector(vec, dofIDMask, mode, tStep);
        answer.append(vec);
    }

    if ( this->computeGtoLRotationMatrix(G2L) ) {
        answer.rotatedWith(G2L, 'n');
    }
}


int
Element :: computeNumberOfGlobalDofs()
{
    return this->computeNumberOfDofs();
}


int
Element :: computeNumberOfPrimaryMasterDofs()
{
    int answer = 0;
    IntArray nodeDofIDMask;

    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        this->giveDofManDofIDMask(i, nodeDofIDMask);
        answer += this->giveDofManager(i)->giveNumberOfPrimaryMasterDofs(nodeDofIDMask);
    }

    for ( int i = 1; i <= giveNumberOfInternalDofManagers(); i++ ) {
        this->giveInternalDofManDofIDMask(i, nodeDofIDMask);
        answer += this->giveInternalDofManager(i)->giveNumberOfPrimaryMasterDofs(nodeDofIDMask);
    }
    return answer;
}


bool
Element :: giveRotationMatrix(FloatMatrix &answer)
{
    bool is_GtoL, is_NtoG;
    FloatMatrix GtoL, NtoG;
    IntArray nodes;
    nodes.enumerate( this->giveNumberOfDofManagers() );

    is_GtoL = this->computeGtoLRotationMatrix(GtoL);
    is_NtoG = this->computeDofTransformationMatrix(NtoG, nodes, true);

#ifdef DEBUG
    if ( is_GtoL ) {
        if ( GtoL.giveNumberOfColumns() != this->computeNumberOfGlobalDofs() ) {
            OOFEM_ERROR("GtoL transformation matrix size mismatch in columns");
        }
        if ( GtoL.giveNumberOfRows() != this->computeNumberOfDofs() ) {
            OOFEM_ERROR("GtoL transformation matrix size mismatch in rows");
        }
    }
    if ( is_NtoG ) {
        if ( NtoG.giveNumberOfColumns() != this->computeNumberOfPrimaryMasterDofs() ) {
            OOFEM_ERROR("NtoG transformation matrix size mismatch in columns");
        }
        if ( NtoG.giveNumberOfRows() != this->computeNumberOfGlobalDofs() ) {
            OOFEM_ERROR("NtoG transformation matrix size mismatch in rows");
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
        answer.clear();
        return false;
    }
    return true;
}


bool
Element :: computeDofTransformationMatrix(FloatMatrix &answer, const IntArray &nodes, bool includeInternal)
{
    bool flag = false;

    // test if transformation is necessary
    for ( int n : nodes ) {
        flag = flag || this->giveDofManager( n )->requiresTransformation();
    }

    if ( !flag ) {
        answer.clear();
        return false;
    }

    // initialize answer
    int gsize = this->computeNumberOfPrimaryMasterDofs();
    answer.resize(this->computeNumberOfGlobalDofs(), gsize);
    answer.zero();

    FloatMatrix dofManT;
    IntArray dofIDmask;
    int nr, nc, lastRowPos = 0, lastColPos = 0;
    // loop over nodes
    for ( int n : nodes ) {
        this->giveDofManDofIDMask(n, dofIDmask);
        if ( !this->giveDofManager( n )->computeM2GTransformation(dofManT, dofIDmask) ) {
            dofManT.resize( dofIDmask.giveSize(), dofIDmask.giveSize() );
            dofManT.zero();
            dofManT.beUnitMatrix();
        }
        nc = dofManT.giveNumberOfColumns();
        nr = dofManT.giveNumberOfRows();
        for ( int j = 1; j <= nr; j++ ) {
            for ( int k = 1; k <= nc; k++ ) {
                // localize node contributions
                answer.at(lastRowPos + j, lastColPos + k) = dofManT.at(j, k);
            }
        }

        lastRowPos += nr;
        lastColPos += nc;
    }
    if ( includeInternal ) {
        for ( int i = 1; i <= this->giveNumberOfInternalDofManagers(); i++ ) {
            this->giveInternalDofManDofIDMask(i, dofIDmask);
            if ( !this->giveInternalDofManager( nodes.at(i) )->computeM2GTransformation(dofManT, dofIDmask) ) {
                dofManT.resize( dofIDmask.giveSize(), dofIDmask.giveSize() );
                dofManT.zero();
                dofManT.beUnitMatrix();
            }
            nc = dofManT.giveNumberOfColumns();
            nr = dofManT.giveNumberOfRows();
            for ( int j = 1; j <= nr; j++ ) {
                for ( int k = 1; k <= nc; k++ ) {
                    // localize node contributions
                    answer.at(lastRowPos + j, lastColPos + k) = dofManT.at(j, k);
                }
            }

            lastRowPos += nr;
            lastColPos += nc;
        }
    }
    answer.resizeWithData(lastRowPos, lastColPos);
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
Element :: giveLocationArray(IntArray &locationArray, const UnknownNumberingScheme &s, IntArray *dofIdArray) const
{
    IntArray masterDofIDs, nodalArray, ids;
    locationArray.clear();
    if ( dofIdArray ) {
        dofIdArray->clear();
    }
    for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
        this->giveDofManDofIDMask(i, ids);
        this->giveDofManager(i)->giveLocationArray(ids, nodalArray, s);
        locationArray.followedBy(nodalArray);
        if ( dofIdArray ) {
            this->giveDofManager(i)->giveMasterDofIDArray(ids, masterDofIDs);
            dofIdArray->followedBy(masterDofIDs);
        }
    }
    for ( int i = 1; i <= this->giveNumberOfInternalDofManagers(); i++ ) {
        this->giveInternalDofManDofIDMask(i, ids);
        this->giveInternalDofManager(i)->giveLocationArray(ids, nodalArray, s);
        locationArray.followedBy(nodalArray);
        if ( dofIdArray ) {
            this->giveInternalDofManager(i)->giveMasterDofIDArray(ids, masterDofIDs);
            dofIdArray->followedBy(masterDofIDs);
        }
    }
}


void
Element :: giveLocationArray(IntArray &locationArray, const IntArray &dofIDMask, const UnknownNumberingScheme &s, IntArray *dofIdArray) const
{
    IntArray masterDofIDs, nodalArray;
    locationArray.clear();
    if ( dofIdArray ) {
        dofIdArray->clear();
    }
    for ( int i = 1; i <=this->giveNumberOfDofManagers(); i++ ) {
        this->giveDofManager(i)->giveLocationArray(dofIDMask, nodalArray, s);
        locationArray.followedBy(nodalArray);
        if ( dofIdArray ) {
            this->giveDofManager(i)->giveMasterDofIDArray(dofIDMask, masterDofIDs);
            dofIdArray->followedBy(masterDofIDs);
        }
    }
    for ( int i = 1; i <= this->giveNumberOfInternalDofManagers(); i++ ) {
        this->giveInternalDofManager(i)->giveLocationArray(dofIDMask, nodalArray, s);
        locationArray.followedBy(nodalArray);
        if ( dofIdArray ) {
            this->giveInternalDofManager(i)->giveMasterDofIDArray(dofIDMask, masterDofIDs);
            dofIdArray->followedBy(masterDofIDs);
        }
    }
}


void
Element :: giveBoundaryLocationArray(IntArray &locationArray, const IntArray &bNodes, const UnknownNumberingScheme &s, IntArray *dofIdArray)
{
    IntArray masterDofIDs, nodalArray, dofIDMask;
    locationArray.clear();
    if ( dofIdArray ) {
        dofIdArray->clear();
    }
    for ( int i = 1; i <= bNodes.giveSize(); i++ ) {
        this->giveDofManDofIDMask(bNodes.at(i), dofIDMask);
        this->giveDofManager( bNodes.at(i) )->giveLocationArray(dofIDMask, nodalArray, s);
        locationArray.followedBy(nodalArray);
        if ( dofIdArray ) {
            this->giveDofManager( bNodes.at(i) )->giveMasterDofIDArray(dofIDMask, masterDofIDs);
            dofIdArray->followedBy(masterDofIDs);
        }
    }
}


void
Element :: giveBoundaryLocationArray(IntArray &locationArray, const IntArray &bNodes, const IntArray &dofIDMask, const UnknownNumberingScheme &s, IntArray *dofIdArray)
{
    IntArray masterDofIDs, nodalArray;
    locationArray.clear();
    if ( dofIdArray ) {
        dofIdArray->clear();
    }
    for ( int i = 1; i <= bNodes.giveSize(); i++ ) {
        this->giveDofManager( bNodes.at(i) )->giveLocationArray(dofIDMask, nodalArray, s);
        locationArray.followedBy(nodalArray);
        if ( dofIdArray ) {
            this->giveDofManager( bNodes.at(i) )->giveMasterDofIDArray(dofIDMask, masterDofIDs);
            dofIdArray->followedBy(masterDofIDs);
        }
    }
}


Material *Element :: giveMaterial()
{
#ifdef DEBUG
    if ( !material ) {
        OOFEM_ERROR("material not defined");
    }
#endif
    return domain->giveMaterial(material);
}


CrossSection *Element :: giveCrossSection()
{
#ifdef DEBUG
    if ( !crossSection ) {
        OOFEM_ERROR("crossSection not defined");
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
{
#ifdef DEBUG
    if ( ( i <= 0 ) || ( i > dofManArray.giveSize() ) ) {
        OOFEM_ERROR("Node %i is not defined", i);
    }
#endif
    return domain->giveDofManager( dofManArray.at(i) );
}

void
Element :: addDofManager(DofManager *dMan)
{
#ifdef DEBUG
    if ( dMan == NULL ) {
        OOFEM_ERROR("Element :: addDofManager - dMan is a null pointer");
    }
#endif
    int size =  dofManArray.giveSize();
    this->dofManArray.resizeWithValues( size + 1 );
    this->dofManArray.at(size + 1) = dMan->giveGlobalNumber();
}

ElementSide *
Element :: giveSide(int i) const
{
#ifdef DEBUG
    if ( ( i <= 0 ) || ( i > dofManArray.giveSize() ) ) {
        OOFEM_ERROR("Side is not defined");
    }
#endif
    return domain->giveSide( dofManArray.at(i) );
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
Element :: setIntegrationRules(std :: vector< std :: unique_ptr< IntegrationRule > > irlist)
{
    integrationRulesArray = std :: move(irlist);
}


void
Element :: giveCharacteristicMatrix(FloatMatrix &answer,
                                    CharType mtrx, TimeStep *tStep)
//
// returns characteristics matrix of receiver according to mtrx
//
{
    OOFEM_ERROR("Unknown Type of characteristic mtrx.");
}


void
Element :: giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep)
//
// returns characteristics vector of receiver according to mtrx
//
{
    OOFEM_ERROR("Unknown Type of characteristic mtrx.");
}


void
Element :: computeLoadVector(FloatArray &answer, BodyLoad *load, CharType type, ValueModeType mode, TimeStep *tStep)
{
    answer.clear();
    OOFEM_ERROR("Unknown load type.");
}


void
Element :: computeBoundarySurfaceLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep, bool global)
{
    answer.clear();
    OOFEM_ERROR("Unknown load type.");
}


void
Element :: computeTangentFromSurfaceLoad(FloatMatrix &answer, SurfaceLoad *load, int boundary, MatResponseMode rmode, TimeStep *tStep)
{
    answer.clear();
}

  void
Element :: computeTangentFromEdgeLoad(FloatMatrix &answer, EdgeLoad *load, int boundary, MatResponseMode rmode, TimeStep *tStep)
{
    answer.clear();
}

void
Element :: computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int edge, CharType type, ValueModeType mode, TimeStep *tStep, bool global)
{
    ///@todo Change the load type to "BoundaryEdgeLoad" maybe?
    answer.clear();
    OOFEM_ERROR("Unknown load type.");
}


double
Element :: giveCharacteristicValue(CharType mtrx, TimeStep *tStep)
//
// returns characteristics value of receiver according to CharType
//
{
    OOFEM_ERROR("Unknown Type of characteristic mtrx.");
    return 0.;
}


void
Element :: initializeFrom(InputRecord &ir)
{
#  ifdef VERBOSE
    // VERBOSE_PRINT1("Instanciating element ",number);
#  endif
    //IR_GIVE_FIELD(ir, material, _IFT_Element_mat);
    material = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, material, _IFT_Element_mat);

    //IR_GIVE_FIELD(ir, crossSection, _IFT_Element_crosssect);
    crossSection = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, crossSection, _IFT_Element_crosssect);

    IR_GIVE_FIELD(ir, dofManArray, _IFT_Element_nodes);

    bodyLoadArray.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, bodyLoadArray, _IFT_Element_bodyload);

    boundaryLoadArray.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, boundaryLoadArray, _IFT_Element_boundaryload);

    elemLocalCS.clear();

    if ( ir.hasField(_IFT_Element_lcs) ) { //local coordinate system
        double n1 = 0.0, n2 = 0.0;
        FloatArray triplets;
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
    }

    partitions.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, partitions, _IFT_Element_partitions);
    if ( ir.hasField(_IFT_Element_remote) ) {
        parallel_mode = Element_remote;
    } else {
        parallel_mode = Element_local;
    }

    activityTimeFunction = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, activityTimeFunction, _IFT_Element_activityTimeFunction);

    IR_GIVE_OPTIONAL_FIELD(ir, numberOfGaussPoints, _IFT_Element_nip);
}


void
Element :: giveInputRecord(DynamicInputRecord &input)
{
    FEMComponent :: giveInputRecord(input);

    input.setField(material, _IFT_Element_mat);

    input.setField(crossSection, _IFT_Element_crosssect);

    input.setField(dofManArray, _IFT_Element_nodes);

    if ( bodyLoadArray.giveSize() > 0 ) {
        input.setField(bodyLoadArray, _IFT_Element_bodyload);
    }


    if ( boundaryLoadArray.giveSize() > 0 ) {
        input.setField(boundaryLoadArray, _IFT_Element_boundaryload);
    }


    if ( elemLocalCS.giveNumberOfRows() > 0 ) {
        FloatArray triplets(6);
        for ( int j = 1; j <= 3; j++ ) {
            triplets.at(j) = elemLocalCS.at(j, 1);
            triplets.at(j + 3) = elemLocalCS.at(j, 2);
        }
        input.setField(triplets, _IFT_Element_lcs);
    }


    if ( partitions.giveSize() > 0 ) {
        input.setField(this->partitions, _IFT_Element_partitions);
        if ( this->parallel_mode == Element_remote ) {
            input.setField(_IFT_Element_remote);
        }
    }

    if ( activityTimeFunction > 0 ) {
        input.setField(activityTimeFunction, _IFT_Element_activityTimeFunction);
    }

    input.setField(numberOfGaussPoints, _IFT_Element_nip);
}


void
Element :: postInitialize()
{
    this->computeGaussPoints();
}


void
Element :: printOutputAt(FILE *file, TimeStep *tStep)
{
    fprintf( file, "element %d (%8d) :\n", this->giveLabel(), this->giveNumber() );

    for ( int i = 1; i <= this->giveNumberOfInternalDofManagers(); ++i ) {
        DofManager *dman = this->giveInternalDofManager(i);
        dman->printOutputAt(file, tStep);
    }

    if (this->isActivated(tStep) ) {
        for ( auto &iRule: integrationRulesArray ) {
            iRule->printOutputAt(file, tStep);
        }
    } else {
        fprintf(file, "is not active in current time step\n");
    }
}


void
Element :: updateYourself(TimeStep *tStep)
// Updates the receiver at end of step.
{
#  ifdef VERBOSE
    // VERBOSE_PRINT1("Updating Element ",number)
#  endif

    for ( auto &iRule: integrationRulesArray ) {
        iRule->updateYourself(tStep);
    }
}


bool
Element :: isActivated(TimeStep *tStep)
{
    if ( activityTimeFunction ) {
        if ( tStep ) {
            return ( domain->giveFunction(activityTimeFunction)->evaluateAtTime( tStep->giveIntrinsicTime() ) > 1.e-3 );
        } else {
            return false;
        }
    } else {
        return true;
    }
}


bool
Element :: isCast(TimeStep *tStep)
{
    // this approach used to work when material was assigned to element
    //    if ( tStep->giveIntrinsicTime() >= this->giveMaterial()->giveCastingTime() )  
    if ( tStep ) {

        double castingTime;
        double tNow = tStep->giveIntrinsicTime();

        for ( auto &iRule : integrationRulesArray ) {
            for ( auto &gp: *iRule ) {
                castingTime = this->giveCrossSection()->giveMaterial(gp)->giveCastingTime();

                if ( tNow < castingTime ) {
                    return false;
                }
            }
        }
        return true;

    } else {
        return false;
    }
}

void
Element :: initForNewStep()
// initializes receiver to new time step or can be used
// if current time step must be restarted
{
    for ( auto &iRule: integrationRulesArray ) {
        for ( auto &gp: *iRule ) {
            this->giveCrossSection()->giveMaterial(gp)->initTempStatus(gp);
        }
    }
}


IntArray
Element::giveBoundaryEdgeNodes(int boundary) const
{
    return this->giveInterpolation()->boundaryEdgeGiveNodes(boundary);
}

IntArray
Element::giveBoundarySurfaceNodes(int boundary) const
{
    return this->giveInterpolation()->boundarySurfaceGiveNodes(boundary);
}

std::unique_ptr<IntegrationRule>
Element::giveBoundaryEdgeIntegrationRule(int order, int boundary)
{
    return this->giveInterpolation()->giveBoundaryEdgeIntegrationRule(order, boundary);
}

std::unique_ptr<IntegrationRule>
Element::giveBoundarySurfaceIntegrationRule(int order, int boundary)
{
    return this->giveInterpolation()->giveBoundarySurfaceIntegrationRule(order, boundary);
}


void Element :: saveContext(DataStream &stream, ContextMode mode)
{
    FEMComponent :: saveContext(stream, mode);

    if ( ( mode & CM_Definition ) ) {
        contextIOResultType iores;
        if ( !stream.write(numberOfDofMans) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream.write(material) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream.write(crossSection) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( mode & CM_DefinitionGlobal ) {
            // send global numbers instead of local ones
            int s = dofManArray.giveSize();
            IntArray globDN(s);
            for ( int i = 1; i <= s; i++ ) {
                globDN.at(i) = this->giveDofManager(i)->giveGlobalNumber();
            }

            if ( ( iores = globDN.storeYourself(stream) ) != CIO_OK ) {
                THROW_CIOERR(iores);
            }
        } else {
            if ( ( iores = dofManArray.storeYourself(stream) ) != CIO_OK ) {
                THROW_CIOERR(iores);
            }
        }

        if ( ( iores = bodyLoadArray.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( ( iores = boundaryLoadArray.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        int numberOfIntegrationRules = (int)integrationRulesArray.size();
        if ( !stream.write(numberOfIntegrationRules) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        for ( auto &iRule: integrationRulesArray ) {
            int _val = iRule->giveIntegrationRuleType();
            if ( !stream.write(_val) ) {
                THROW_CIOERR(CIO_IOERR);
            }
        }

        int _mode;
        if ( !stream.write(globalNumber) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        _mode = parallel_mode;
        if ( !stream.write(_mode) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( ( iores = partitions.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    for ( auto &iRule: integrationRulesArray ) {
        iRule->saveContext(stream, mode);
    }
}


void Element :: restoreContext(DataStream &stream, ContextMode mode)
{
    contextIOResultType iores;
    int _nrules;

    FEMComponent :: restoreContext(stream, mode);

    if ( mode & CM_Definition ) {
        if ( !stream.read(numberOfDofMans) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream.read(material) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream.read(crossSection) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( ( iores = dofManArray.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( ( iores = bodyLoadArray.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( ( iores = boundaryLoadArray.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }

        if ( !stream.read(_nrules) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        // restore integration rules
        IntArray dtypes(_nrules);
        for ( int i = 1; i <= _nrules; i++ ) {
            if ( !stream.read(dtypes.at(i)) ) {
                THROW_CIOERR(CIO_IOERR);
            }
        }

        if ( _nrules != (int)integrationRulesArray.size() ) {

            // AND ALLOCATE NEW ONE
            integrationRulesArray.resize( _nrules );
            for ( int i = 0; i < _nrules; i++ ) {
                integrationRulesArray [ i ] = classFactory.createIRule( ( IntegrationRuleType ) dtypes[i], i + 1, this );
            }
        } else {
            for ( int i = 0; i < _nrules; i++ ) {
                if ( integrationRulesArray [ i ]->giveIntegrationRuleType() != dtypes[i] ) {
                    integrationRulesArray [ i ] = classFactory.createIRule( ( IntegrationRuleType ) dtypes[i], i + 1, this );
                }
            }
        }

        int _mode;
        if ( !stream.read(globalNumber) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        if ( !stream.read(_mode) ) {
            THROW_CIOERR(CIO_IOERR);
        }

        parallel_mode = ( elementParallelMode ) _mode;
        if ( ( iores = partitions.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }


    for ( auto &iRule: integrationRulesArray ) {
        iRule->restoreContext(stream, mode);
    }
}


double
Element :: computeVolumeAreaOrLength()
// the element computes its volume, area or length
// (depending on the spatial dimension of that element)
{
    double answer = 0.;
    IntegrationRule *iRule = this->giveDefaultIntegrationRulePtr();
    if ( iRule ) {
        for ( auto &gp: *iRule ) {
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
    double volume = this->computeVolumeAreaOrLength();
    if ( volume < 0. ) {
        return -1.; // means "cannot be evaluated"
    }

    int dim = this->giveSpatialDimension();
    switch ( dim ) {
    case 1: return volume;

    case 2: return sqrt(volume);

    case 3: return cbrt(volume);
    }

    return -1.; // means "cannot be evaluated"
}


double
Element :: computeVolume()
{
    FEInterpolation3d *fei = dynamic_cast< FEInterpolation3d * >( this->giveInterpolation() );
#ifdef DEBUG
    if ( !fei ) {
        OOFEM_ERROR("Function not overloaded and necessary interpolator isn't available");
        return 0.0;
    }
#endif
    return fei->giveVolume( FEIElementGeometryWrapper(this) );
}


double
Element :: computeArea()
{
    FEInterpolation2d *fei = dynamic_cast< FEInterpolation2d * >( this->giveInterpolation() );
#ifdef DEBUG
    if ( !fei ) {
        OOFEM_ERROR("Function not overloaded and necessary interpolator isn't available");
        return 0.0;
    }
#endif
    return fei->giveArea( FEIElementGeometryWrapper(this) );
}


double
Element :: computeLength()
{
    FEInterpolation1d *fei = dynamic_cast< FEInterpolation1d * >( this->giveInterpolation() );
#ifdef DEBUG
    if ( !fei ) {
        OOFEM_ERROR("Function not overloaded and necessary interpolator isn't available");
        return 0.0;
    }
#endif
    return fei->giveLength( FEIElementGeometryWrapper(this) );
}


double
Element :: giveLengthInDir(const FloatArray &normalToCrackPlane)
//
// returns receiver's size projected onto the direction given by normalToCrackPlane
// (exploited by the crack band approach)
//
{
    double maxDis, minDis;
    int nnode = giveNumberOfNodes();

    const auto &coords = this->giveNode(1)->giveCoordinates();
    minDis = maxDis = normalToCrackPlane.dotProduct( coords, coords.giveSize() );

    for ( int i = 2; i <= nnode; i++ ) {
        const auto &coords = this->giveNode(i)->giveCoordinates();
        double dis = normalToCrackPlane.dotProduct( coords, coords.giveSize() );
        if ( dis > maxDis ) {
            maxDis = dis;
        } else if ( dis < minDis ) {
            minDis = dis;
        }
    }

    return maxDis - minDis;
}

double 
Element :: giveCharacteristicLengthForPlaneElements(const FloatArray &normalToCrackPlane) 
//
// returns receiver's size projected onto the direction given by normalToCrackPlane
// or, if that direction is not in-plane, the square root of element area
// (this can happen if the crack normal is set to the maximum principal stress direction
//  and the in-plane principal stresses are negative)
//
{
    if ( normalToCrackPlane.at(3) < 0.5 ) { // check whether the projection direction is in-plane 
        return this->giveLengthInDir(normalToCrackPlane);
    } else { // if not, compute the size from element area
        return this->computeMeanSize();
    }
}

double 
Element :: giveCharacteristicLengthForAxisymmElements(const FloatArray &normalToCrackPlane) 
//
// returns receiver's size projected onto the direction given by normalToCrackPlane
// or, if that direction is not in-plane, the distance from the axis of symmetry
// multiplied by pi, assuming that two symmetrically located radial cracks will form
// (this can happen if the cracking is caused by the hoop stress)
//
{
    if ( normalToCrackPlane.at(3) < 0.5 ) { // check whether the projection direction is in-plane 
        return this->giveLengthInDir(normalToCrackPlane);
    } else { // if not, take the average distance from axis of symmetry multiplied by pi
        double r = 0.;
        for ( int i = 1; i <= this->giveNumberOfDofManagers(); i++ ) {
            r += this->giveNode(i)->giveCoordinate(1);
        }
        r = r * M_PI / ( ( double ) this->giveNumberOfDofManagers() );
        return r;
    }
}

int
Element :: computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords)
{
    FEInterpolation *fei = this->giveInterpolation();
#ifdef DEBUG
    if ( !fei ) {
        answer.clear();
        return false;
    }
#endif
    fei->local2global( answer, lcoords, FEIElementGeometryWrapper(this) );
    return true;
}


bool
Element :: computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords)
{
    FEInterpolation *fei = this->giveInterpolation();
    if ( fei ) {
        return fei->global2local( answer, gcoords, FEIElementGeometryWrapper(this) );
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
    } else {
        answer.clear();
    }

    return 0;
}

void
Element :: giveLocalCoordinateSystemVector(InternalStateType isttype, FloatArray &answer)
{
    int col = 0;
    FloatMatrix rotMat;
    
    if ( isttype == IST_X_LCS ) {
        col = 1;
    } else if ( isttype == IST_Y_LCS ) {
        col = 2;
    } else if ( isttype == IST_Z_LCS ) {
        col = 3;
    } else {
        OOFEM_ERROR("Only IST_X_LCS, IST_Y_LCS, IST_Z_LCS options are permitted, you provided %s", __InternalStateTypeToString(isttype));
    }
    
    if ( !this->giveLocalCoordinateSystem(rotMat) ) {
        rotMat.resize(3, 3);
        rotMat.beUnitMatrix();
    }
    answer.beRowOf(rotMat, col);
}


void
Element :: computeMidPlaneNormal(FloatArray &answer, const GaussPoint *)
// valid only for plane elements (shells, plates, ....)
// computes mid-plane normal at gaussPoint - for materials with orthotrophy
{
    OOFEM_ERROR("Unable to compute mid-plane normal, not supported");
}


int
Element :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( type == IST_ErrorIndicatorLevel ) {
        ErrorEstimator *ee = this->giveDomain()->giveErrorEstimator();
        if ( ee ) {
            answer.resize(1);
            ///@todo Which "error type" should be used? Why are there several? I don't see the point of this enum when 
            /// there could be different function calls just as well (and different IST values)
            answer.at(1) = ee->giveElementError(internalStressET, this, tStep);
        } else {
            answer.clear();
            return 0;
        }

        return 1;
    } else if ( type == IST_InternalStressError ) {
        ErrorEstimator *ee = this->giveDomain()->giveErrorEstimator();
        if ( ee ) {
            answer.resize(1);
            answer.at(1) = ee->giveElementError(internalStressET, this, tStep);
        } else {
            answer.clear();
            return 0;
        }
        return 1;
    } else if ( type == IST_PrimaryUnknownError ) {
        ErrorEstimator *ee = this->giveDomain()->giveErrorEstimator();
        if ( ee ) {
            answer.resize(1);
            answer.at(1) = ee->giveElementError(primaryUnknownET, this, tStep);
        } else {
            answer.clear();
            return 0;
        }
        return 1;
    } else if ( type == IST_CrossSectionNumber ) {
        answer.resize(1);
        answer.at(1) = gp->giveCrossSection()->giveNumber();
        return 1;
    } else if ( type == IST_ElementNumber ) {
        answer.resize(1);
        answer.at(1) = this->giveNumber();
        return 1;
    } else if ( type == IST_X_LCS || type == IST_Y_LCS || type == IST_Z_LCS ) {
        this->giveLocalCoordinateSystemVector(type, answer);
        return 1;
    } else {
        return this->giveCrossSection()->giveIPValue(answer, gp, type, tStep);
    }
}

int
Element :: giveGlobalIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    InternalStateValueType valtype = giveInternalStateValueType(type);
    if ( elemLocalCS.isNotEmpty() && valtype != ISVT_SCALAR ) {
        FloatArray ans;
        FloatMatrix full;
        int ret = this->giveIPValue(ans, gp, type, tStep);
        if ( ret == 0 ) return 0;
        if ( valtype == ISVT_VECTOR ) {
            ///@todo Check transpose here
            answer.beProductOf(elemLocalCS, ans);
            return 1;
        }

        // Tensors are more complicated to transform, easiest to write them in matrix form first
        if ( valtype == ISVT_TENSOR_S3E ) {
            ans.at(4) *= 0.5;
            ans.at(5) *= 0.5;
            ans.at(6) *= 0.5;
        }  else if ( valtype != ISVT_TENSOR_G && valtype != ISVT_TENSOR_S3 ) {
            OOFEM_ERROR("Unsupported internal state value type for computing global IP value");
        }
        ///@todo Check transpose here
        full.beMatrixForm(ans);
        full.rotatedWith(elemLocalCS, 'n');
        answer.beVectorForm(full);
        if ( valtype == ISVT_TENSOR_S3 ) {
            answer.resizeWithValues(6);
        } else if ( valtype == ISVT_TENSOR_S3E ) {
            answer.at(4) += answer.at(7);
            answer.at(5) += answer.at(8);
            answer.at(6) += answer.at(9);
            answer.resizeWithValues(6);
        }
        return 1;
    } else {
        return this->giveIPValue(answer, gp, type, tStep);
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

    OOFEM_ERROR("failure (maybe new element type was registered)");
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

    OOFEM_ERROR("failure, unsupported geometry type (%s)",
            __Element_Geometry_TypeToString( this->giveGeometryType() ));
    return 0; // to make compiler happy
}


int
Element :: adaptiveMap(Domain *oldd, TimeStep *tStep)
{
    int result = 1;
    CrossSection *cs = this->giveCrossSection();

    for ( auto &iRule: integrationRulesArray ) {
        for ( auto &gp: *iRule ) {
            MaterialModelMapperInterface *interface = static_cast< MaterialModelMapperInterface * >
                ( cs->giveMaterial(gp)->giveInterface(MaterialModelMapperInterfaceType) );
            if ( interface ) {
                result &= interface->MMI_map(gp, oldd, tStep);
            } else {
                result = 0;
            }
        }
    }

    return result;
}

int
Element :: mapStateVariables(Domain &iOldDom, const TimeStep &iTStep)
{
    int result = 1;

    // create source set (this is quite inefficient here as done for each element.
    // the alternative MaterialModelMapperInterface approach allows to cache sets on material model
    Set sourceElemSet = Set(0, & iOldDom);
    int materialNum = this->giveMaterial()->giveNumber();
    sourceElemSet.setElementList( iOldDom.giveElementsWithMaterialNum(materialNum) );

    for ( auto &iRule: integrationRulesArray ) {
        for ( GaussPoint *gp: *iRule ) {

            MaterialStatus *ms = dynamic_cast< MaterialStatus * >( gp->giveMaterialStatus() );
            if ( ms == NULL ) {
                OOFEM_ERROR("failed to fetch MaterialStatus.");
            }

            MaterialStatusMapperInterface *interface = dynamic_cast< MaterialStatusMapperInterface * >(ms);
            if ( interface == NULL ) {
                OOFEM_ERROR("Failed to fetch MaterialStatusMapperInterface.");
            }

            result &= interface->MSMI_map( *gp, iOldDom, sourceElemSet, iTStep, * ( ms ) );
        }
    }

    return result;
}


int
Element :: adaptiveFinish(TimeStep *tStep)
{
    MaterialModelMapperInterface *interface = static_cast< MaterialModelMapperInterface * >
                                              ( this->giveMaterial()->giveInterface(MaterialModelMapperInterfaceType) );

    if ( !interface ) {
        return 0;
    }
#if 0
    int result = 1;
    for ( auto &iRule: integrationRulesArray ) {
        for ( GaussPoint *gp: *iRule ) {
            result &= interface->MMI_finish(tStep);
        }
    }

    return result;
#else
    return interface->MMI_finish(tStep);
#endif
}


void
Element :: updateLocalNumbering(EntityRenumberingFunctor &f)
{
    for ( auto &dnum : dofManArray ) {
        dnum = f(dnum, ERS_DofManager);
    }
    
}


integrationDomain
Element :: giveIntegrationDomain() const
{
    FEInterpolation *fei = this->giveInterpolation();
    return fei ? fei->giveIntegrationDomain() : _UnknownIntegrationDomain;
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
    answer.clear();
    return false;
}


int
Element :: packUnknowns(DataStream &buff, TimeStep *tStep)
{
    int result = 1;

    for ( auto &iRule: integrationRulesArray ) {
        for ( GaussPoint *gp: *iRule ) {
            result &= this->giveCrossSection()->packUnknowns( buff, tStep, gp );
        }
    }

    return result;
}


int
Element :: unpackAndUpdateUnknowns(DataStream &buff, TimeStep *tStep)
{
    int result = 1;

    for ( auto &iRule: integrationRulesArray ) {
        for ( GaussPoint *gp: *iRule ) {
            result &= this->giveCrossSection()->unpackAndUpdateUnknowns( buff, tStep, gp );
        }
    }

    return result;
}


int
Element :: estimatePackSize(DataStream &buff)
{
    int result = 0;

    for ( auto &iRule: integrationRulesArray ) {
        for ( GaussPoint *gp: *iRule ) {
            result += this->giveCrossSection()->estimatePackSize( buff, gp );
        }
    }

    return result;
}


double
Element :: predictRelativeComputationalCost()
{
    double wgt = 0;
    IntegrationRule *iRule = this->giveDefaultIntegrationRulePtr();
    for ( GaussPoint *gp: *iRule ) {
        wgt += this->giveCrossSection()->predictRelativeComputationalCost( gp );
    }

    return ( this->giveRelativeSelfComputationalCost() * wgt );
}


#ifdef __OOFEG
void
Element :: drawYourself(oofegGraphicContext &gc, TimeStep *tStep)
{
    OGC_PlotModeType mode = gc.giveIntVarPlotMode();

    if ( mode == OGC_rawGeometry ) {
        this->drawRawGeometry(gc, tStep);
    } else if ( mode == OGC_elementAnnotation ) {
        this->drawAnnotation(gc, tStep);
    } else if ( mode == OGC_deformedGeometry ) {
        this->drawDeformedGeometry(gc, tStep, DisplacementVector);
    } else if ( mode == OGC_eigenVectorGeometry ) {
        this->drawDeformedGeometry(gc, tStep, EigenVector);
    } else if ( mode == OGC_scalarPlot ) {
        this->drawScalar(gc, tStep);
    } else if ( mode == OGC_elemSpecial ) {
        this->drawSpecial(gc, tStep);
    } else {
        OOFEM_ERROR("unsupported mode");
    }
}


void
Element :: drawAnnotation(oofegGraphicContext &gc, TimeStep *tStep)
{
    int count = 0;
    Node *node;
    WCRec p [ 1 ]; /* point */
    GraphicObj *go;
    char num [ 30 ];

    p [ 0 ].x = p [ 0 ].y = p [ 0 ].z = 0.0;
    // compute element center
    for ( int i = 1; i <= numberOfDofMans; i++ ) {
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
    sprintf( num, "%d(%d)", this->giveNumber(), this->giveGlobalNumber() );

    go = CreateAnnText3D(p, num);
    EGWithMaskChangeAttributes(COLOR_MASK | LAYER_MASK, go);
    EMAddGraphicsToModel(ESIModel(), go);
}


int
Element :: giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                   int node, TimeStep *tStep)
{
    if ( type == IST_RelMeshDensity ) {
        ErrorEstimator *ee = this->giveDomain()->giveErrorEstimator();
        if ( ee ) {
            answer.resize(1);
            answer.at(1) = this->giveDomain()->giveErrorEstimator()->giveRemeshingCrit()->
            giveRequiredDofManDensity(this->giveNode(node)->giveNumber(), tStep, 1);
            return 1;
        } else {
            answer.clear();
            return 0;
        }
    } else {
        if ( mode == ISM_recovered ) {
            const FloatArray *nodval;
            NodalRecoveryModel *smoother = this->giveDomain()->giveSmoother();
            int result = smoother->giveNodalVector( nodval, this->giveNode(node)->giveNumber() );
            if ( nodval ) {
                answer = * nodval;
            } else {
                answer.clear();
            }

            return result;
        } else {
            return 0;
        }
    }
}


#endif
} // end namespace oofem
