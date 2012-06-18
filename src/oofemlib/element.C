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

#include "element.h"
#include "crosssection.h"
#include "integrationrule.h"
#include "errorestimator.h"
#include "intarray.h"
#include "flotarry.h"
#include "flotmtrx.h"
#include "primaryfield.h"
#include "verbose.h"
#include "entityrenumberingscheme.h"
#include "error.h"
#include "usrdefsub.h"
#include "datastream.h"
#include "materialmapperinterface.h"
#include "contextioerr.h"
#include "mathfem.h"
#include "feinterpol.h"
#include "feinterpol1d.h"
#include "feinterpol2d.h"
#include "feinterpol3d.h"

#ifndef __MAKEDEPEND
 #include <cstdio>
#endif

namespace oofem {
Element :: Element(int n, Domain *aDomain) :
    FEMComponent(n, aDomain), dofManArray(), bodyLoadArray(), boundaryLoadArray()
{
    material           = 0;
    numberOfDofMans    = 0;
    numberOfIntegrationRules = 0;
    locationArray      = NULL;
    integrationRulesArray  = NULL;
}


Element :: ~Element()
{
    int i;

    delete locationArray;
    if ( integrationRulesArray ) {
        for ( i = 0; i < numberOfIntegrationRules; i++ ) {
            delete integrationRulesArray [ i ];
        }

        delete[] integrationRulesArray;
    }
}


void
Element :: computeVectorOf(EquationID type, ValueModeType u, TimeStep *stepN, FloatArray &answer)
// Forms the vector containing the values of the unknown 'u' (e.g., the
// Total value) of the dofs of the callers local cs.
{
    int i, j, k, nDofs, size;
    IntArray elementNodeMask;
    FloatMatrix G2L;
    FloatArray vec;
    answer.resize( size = this->computeNumberOfGlobalDofs(type) );

    k = 0;
    for ( i = 1; i <= numberOfDofMans; i++ ) {
        this->giveDofManDofIDMask(i, type, elementNodeMask);
        this->giveDofManager(i)->giveUnknownVector(vec, elementNodeMask, type, u, stepN);
        nDofs = vec.giveSize();
        for ( j = 1; j <= nDofs; j++ ) {
            answer.at(++k) = vec.at(j);
        }
    }

    for ( i = 1; i <= giveNumberOfInternalDofManagers(); i++ ) {
        this->giveInternalDofManDofIDMask(i, type, elementNodeMask);
        this->giveInternalDofManager(i)->giveUnknownVector(vec, elementNodeMask, type, u, stepN);
        nDofs = vec.giveSize();
        for ( j = 1; j <= nDofs; j++ ) {
            answer.at(++k) = vec.at(j);
        }
    }

    if (this->computeGtoLRotationMatrix(G2L)) {
        answer.rotatedWith(G2L, 'n');
    }
}


void
Element :: computeVectorOf(PrimaryField &field, ValueModeType u, TimeStep *stepN, FloatArray &answer)
// Forms the vector containing the values of the unknown 'u' (e.g., the
// Total value) of the dofs of the receiver's nodes (in nodal cs).
// Dofs cointaining expected unknowns (of expected type) are determined
// using this->GiveNodeDofIDMask function
{
    int i, j, k, nDofs, size;
    IntArray elementNodeMask;
    FloatMatrix G2L;
    FloatArray vec;
    answer.resize( size = this->computeNumberOfGlobalDofs( field.giveEquationID() ) );

    k = 0;
    for ( i = 1; i <= numberOfDofMans; i++ ) {
        this->giveDofManDofIDMask(i, field.giveEquationID(), elementNodeMask);
        this->giveDofManager(i)->giveUnknownVector(vec, elementNodeMask, field, u, stepN);
        nDofs = vec.giveSize();
        for ( j = 1; j <= nDofs; j++ ) {
            answer.at(++k) = vec.at(j);
        }
    }

    for ( i = 1; i <= giveNumberOfInternalDofManagers(); i++ ) {
        this->giveInternalDofManDofIDMask(i, field.giveEquationID(), elementNodeMask);
        this->giveInternalDofManager(i)->giveUnknownVector(vec, elementNodeMask, field, u, stepN);
        nDofs = vec.giveSize();
        for ( j = 1; j <= nDofs; j++ ) {
            answer.at(++k) = vec.at(j);
        }
    }

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
    int i, j, k, size, nDofs;
    IntArray elementNodeMask, dofMask;
    FloatMatrix G2L;
    FloatArray vec;

    answer.resize( size = this->computeNumberOfGlobalDofs(ut) );

    k = 0;
    for ( i = 1; i <= numberOfDofMans; i++ ) {
        this->giveDofManDofIDMask(i, ut, elementNodeMask);
        this->giveDofManager(i)->givePrescribedUnknownVector(vec, elementNodeMask, mode, stepN);
        nDofs = vec.giveSize();
        for ( j = 1; j <= nDofs; j++ ) {
            answer.at(++k) = vec.at(j);
        }
    }

    for ( i = 1; i <= giveNumberOfInternalDofManagers(); i++ ) {
        this->giveInternalDofManDofIDMask(i, ut, elementNodeMask);
        this->giveInternalDofManager(i)->givePrescribedUnknownVector(vec, elementNodeMask, mode, stepN);
        nDofs = vec.giveSize();
        for ( j = 1; j <= nDofs; j++ ) {
            answer.at(++k) = vec.at(j);
        }
    }

    if (this->computeGtoLRotationMatrix(G2L)) {
        answer.rotatedWith(G2L, 'n');
    }
}


int
Element :: computeNumberOfGlobalDofs(EquationID eid)
{
    return this->computeNumberOfDofs(eid);
}


int
Element :: computeNumberOfPrimaryMasterDofs(EquationID ut)
{
    if ( this->locationArray ) {
        return this->locationArray->giveSize();
    } else {
        int i, answer = 0;
        IntArray nodeDofIDMask, dofMask;

        for ( i = 1; i <= numberOfDofMans; i++ ) {
            this->giveDofManDofIDMask(i, ut, nodeDofIDMask);
            this->giveDofManager(i)->giveDofArray(nodeDofIDMask, dofMask);
            answer += this->giveDofManager(i)->giveNumberOfPrimaryMasterDofs(dofMask);
        }

        for ( i = 1; i <= giveNumberOfInternalDofManagers(); i++ ) {
            this->giveInternalDofManDofIDMask(i, ut, nodeDofIDMask);
            this->giveInternalDofManager(i)->giveDofArray(nodeDofIDMask, dofMask);
            answer += this->giveDofManager(i)->giveNumberOfPrimaryMasterDofs(dofMask);
        }
        return answer;
    }
}


bool
Element :: giveRotationMatrix(FloatMatrix &answer, EquationID eid)
{
    bool is_GtoL, is_NtoG;
    FloatMatrix GtoL, NtoG;

    is_GtoL = this->computeGtoLRotationMatrix(GtoL);
    is_NtoG = this->computeDofTransformationMatrix(NtoG, eid);

#ifdef DEBUG
    if ( is_GtoL ) {
        if ( GtoL.giveNumberOfColumns() != this->computeNumberOfGlobalDofs(eid) ) {
            _error("StructuralElement :: updateRotationMatrix - GtoL transformation matrix size mismatch in columns");
        }
        if ( GtoL.giveNumberOfRows() != this->computeNumberOfDofs(eid) ) {
            _error("StructuralElement :: updateRotationMatrix - GtoL transformation matrix size mismatch in rows");
        }
    }
    if ( is_NtoG ) {
        if ( NtoG.giveNumberOfColumns() != this->computeNumberOfPrimaryMasterDofs(eid) ) {
            _error("StructuralElement :: updateRotationMatrix - NtoG transformation matrix size mismatch in columns");
        }
        if ( NtoG.giveNumberOfRows() != this->computeNumberOfGlobalDofs(eid) ) {
            _error("StructuralElement :: updateRotationMatrix - NtoG transformation matrix size mismatch in rows");
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
Element :: computeDofTransformationMatrix(FloatMatrix &answer, EquationID eid)
{
    bool flag = false;
    int numberOfDofMans = this->giveNumberOfDofManagers();

    // test if transformation is necessary
    for (int i = 1; i <= numberOfDofMans; i++ ) {
        flag = flag || this->giveDofManager(i)->requiresTransformation();
    }

    if ( !flag ) {
        answer.beEmptyMtrx();
        return false;
    }

    // initialize answer
    int gsize = this->computeNumberOfPrimaryMasterDofs(eid);
    answer.resize(this->computeNumberOfGlobalDofs(eid), gsize);
    answer.zero();

    FloatMatrix dofManT;
    IntArray dofIDmask;
    int nr, nc, lastRowPos = 0, lastColPos = 0;
    // loop over nodes
    for (int i = 1; i <= numberOfDofMans; i++ ) {
        this->giveDofManDofIDMask(i, eid, dofIDmask);
        if (!this->giveDofManager(i)->computeM2GTransformation(dofManT, dofIDmask)) {
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
Element :: giveLocationArray(IntArray &locationArray, EquationID ut, const UnknownNumberingScheme &s) const
// Returns the location array of the receiver. This array is obtained by
// simply appending the location array of every node of the receiver.
{
    IntArray nodeDofIDMask;
    IntArray nodalArray;
    int i;

    if ( s.isDefault() && this->locationArray ) {
        locationArray = * this->locationArray;
        return;
    } else {
        locationArray.resize(0);
        for ( i = 1; i <= numberOfDofMans; i++ ) {
            this->giveDofManDofIDMask(i, ut, nodeDofIDMask);
            this->giveDofManager(i)->giveLocationArray(nodeDofIDMask, nodalArray, s);
            locationArray.followedBy(nodalArray);
        }
    }
}


void
Element :: invalidateLocationArray()
{
    // invalitaes current location array in receiver
    // next call of giveLocationArray() will asemble
    // new location array
    // used mainly for model supporting dynamic changes of
    // static system


    // force assembling
    // of new location array
    if ( locationArray ) {
        delete locationArray;
    }

    locationArray = NULL;
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
    int n;
#ifdef DEBUG
    if ( ( i <= 0 ) || ( i > dofManArray.giveSize() ) ) {
        OOFEM_ERROR2("giveNode: Node %i is not defined", i);
    }
#endif
    n = dofManArray.at(i);
    return domain->giveDofManager(n);
}


Node *
Element :: giveNode(int i) const
// Returns the i-th node of the receiver.
{
    int n;
#ifdef DEBUG
    if ( ( i <= 0 ) || ( i > dofManArray.giveSize() ) ) {
        _error("giveNode: Node is not defined");
    }
#endif
    n = dofManArray.at(i);
    return domain->giveNode(n);
}


ElementSide *
Element :: giveSide(int i) const
// Returns the i-th side of the receiver.
{
    int n;
#ifdef DEBUG
    if ( ( i <= 0 ) || ( i > dofManArray.giveSize() ) ) {
        _error("giveNode: Side is not defined");
    }
#endif
    n = dofManArray.at(i);
    return domain->giveSide(n);
}


void
Element :: setDofManagers(const IntArray &_dmans)
{
    this->dofManArray = _dmans;
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
    IR_GIVE_FIELD(ir, material, IFT_Element_mat, "mat"); // Macro

    IR_GIVE_FIELD(ir, crossSection, IFT_Element_crosssect, "crosssect"); // Macro

    IR_GIVE_FIELD(ir, dofManArray, IFT_Element_nodes, "nodes"); // Macro

    //sideArray.resize(0);
    //IR_GIVE_OPTIONAL_FIELD (ir, sideArray, IFT_Element_sides, "sides"); // Macro

    bodyLoadArray.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, bodyLoadArray, IFT_Element_bodyload, "bodyloads"); // Macro

    boundaryLoadArray.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, boundaryLoadArray, IFT_Element_boundaryload, "boundaryloads"); // Macro

    elemLocalCS.resize(0, 0);

    if ( ir->hasField(IFT_Element_lcs, "lcs") ) { //local coordinate system
        double n1 = 0.0, n2 = 0.0;
        int j;
        FloatArray triplets;
        triplets.resize(0);
        IR_GIVE_OPTIONAL_FIELD(ir, triplets, IFT_Element_lcs, "lcs");
        elemLocalCS.resize(3, 3);
        for ( j = 1; j <= 3; j++ ) {
            elemLocalCS.at(j, 1) = triplets.at(j);
            n1 += triplets.at(j) * triplets.at(j);
            elemLocalCS.at(j, 2) = triplets.at(j + 3);
            n2 += triplets.at(j + 3) * triplets.at(j + 3);
        }

        n1 = sqrt(n1);
        n2 = sqrt(n2);
        for ( j = 1; j <= 3; j++ ) { // normalize e1' e2'
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
 #ifndef __ENABLE_COMPONENT_LABELS
    globalNumber = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, globalNumber, IFT_Element_globnum, "globnum"); // Macro
 #endif
    partitions.resize(0);
    IR_GIVE_OPTIONAL_FIELD(ir, partitions, IFT_Element_partitions, "partitions"); // Macro
    // if (hasString (initString, "shared")) parallel_mode = Element_shared;
    if ( ir->hasField(IFT_Element_remote, "remote") ) {
        parallel_mode = Element_remote;
    } else {
        parallel_mode = Element_local;
    }

#endif

    return IRRT_OK;
}


void
Element :: postInitialize()
{
    this->computeGaussPoints();
}


void Element :: printOutputAt(FILE *file, TimeStep *stepN)
// Performs end-of-step operations.
{
    int i;

#if defined ( __PARALLEL_MODE ) || defined ( __ENABLE_COMPONENT_LABELS )
    fprintf( file, "element %d (%8d) :\n", this->giveLabel(), this->giveNumber() );
#else
    fprintf(file, "element %d :\n", number);
#endif

    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        integrationRulesArray [ i ]->printOutputAt(file, stepN);
    }
}


void
Element :: updateYourself(TimeStep *tStep)
// Updates the receiver at end of step.
{
    int i;

#  ifdef VERBOSE
    // VERBOSE_PRINT1("Updating Element ",number)
#  endif
    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        integrationRulesArray [ i ]->updateYourself(tStep);
    }
}


void
Element :: initForNewStep()
// initializes receiver to new time step or can be used
// if current time step must be restarted
//
// call material->initGpForNewStep() for all GPs.
//
{
    int i;

    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
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
    int i, _val;

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
            for ( i = 1; i <= s; i++ ) {
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

        for ( i = 0; i < numberOfIntegrationRules; i++ ) {
            _val = integrationRulesArray [ i ]->giveClassID();
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

    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
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
    int i, _nrules;

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
        for ( i = 1; i <= _nrules; i++ ) {
            if ( !stream->read(& dtypes.at(i), 1) ) {
                THROW_CIOERR(CIO_IOERR);
            }
        }

        if ( _nrules != numberOfIntegrationRules ) {
            // delete old int rule array
            if ( integrationRulesArray ) {
                for ( i = 0; i < numberOfIntegrationRules; i++ ) {
                    delete integrationRulesArray [ i ];
                }

                delete[] integrationRulesArray;
            }

            // AND ALLOCATE NEW ONE
            integrationRulesArray = new IntegrationRule * [ _nrules ];
            for ( i = 0; i < _nrules; i++ ) {
                integrationRulesArray [ i ] = CreateUsrDefIRuleOfType( ( classType ) dtypes(i), i + 1, this );
            }

            numberOfIntegrationRules = _nrules;
        } else {
            for ( i = 0; i < numberOfIntegrationRules; i++ ) {
                if ( integrationRulesArray [ i ]->giveClassID() != dtypes(i) ) {
                    delete integrationRulesArray [ i ];
                    integrationRulesArray [ i ] = CreateUsrDefIRuleOfType( ( classType ) dtypes(i), i + 1, this );
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


    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
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
    int i;
    GaussPoint *gp;
    double answer = 0.;
    IntegrationRule *iRule = integrationRulesArray [ giveDefaultIntegrationRule() ];
    if ( iRule ) {
        for ( i = 0; i < iRule->getNumberOfIntegrationPoints(); i++ ) {
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


int
Element :: computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords)
{
    FEInterpolation *fei = this->giveInterpolation();
    if (fei) {
        return fei->global2local(answer, gcoords, FEIElementGeometryWrapper(this));
    } else {
        answer.resize(0);
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
Element :: giveIntVarCompFullIndx(IntArray &answer, InternalStateType type)
{
    if ( ( type == IST_ErrorIndicatorLevel ) || ( type == IST_RelMeshDensity ) ||
        ( type == IST_InternalStressError ) || ( type == IST_PrimaryUnknownError ) ) {
        answer.resize(1);
        answer.at(1) = 1;
        return 1;
    } else {
        return this->giveCrossSection()->giveIntVarCompFullIndx( answer, type,
                                                                this->giveDefaultIntegrationRulePtr()->
                                                                getIntegrationPoint(0)->giveMaterialMode(), this->giveMaterial() );
    }
}


InternalStateValueType
Element :: giveIPValueType(InternalStateType type)
{
    if ( ( type == IST_ErrorIndicatorLevel ) || ( type == IST_RelMeshDensity ) ||
        ( type == IST_InternalStressError ) || ( type == IST_PrimaryUnknownError ) ) {
        return ISVT_SCALAR;
    } else {
        return this->giveCrossSection()->giveIPValueType( type, this->giveMaterial() );
    }
}


int
Element :: giveIPValueSize(InternalStateType type, GaussPoint *gp)
{
    if ( ( type == IST_ErrorIndicatorLevel ) || ( type == IST_RelMeshDensity ) ||
        ( type == IST_InternalStressError ) || ( type == IST_PrimaryUnknownError ) ) {
        return 1;
    } else {
        return this->giveCrossSection()->giveIPValueSize(type, gp);
    }
}


int
Element :: giveSpatialDimension(void) const
{
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
        return 2;

    case EGT_tetra_1:
    case EGT_tetra_2:
    case EGT_hexa_1:
    case EGT_hexa_2:
        return 3;

    case EGT_Composite:
    case EGT_unknown:
        break;
    }

    _error("giveSpatialDimension: failure (maybe new element type was registered)");
    return 0; //to make compiler happy
}


int
Element :: giveNumberOfBoundarySides(void) const
{
    switch ( this->giveGeometryType() ) {
    case EGT_point:
        return 0;

    case EGT_line_1:
    case EGT_line_2:
        return 2;

    case EGT_triangle_1:
    case EGT_triangle_2:
        return 3;

    case EGT_quad_1:
    case EGT_quad_2:
        return 4;

    case EGT_tetra_1:
    case EGT_tetra_2:
        return 4;

    case EGT_hexa_1:
    case EGT_hexa_2:
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
    int i, j, result = 1;
    IntegrationRule *iRule;
    MaterialModelMapperInterface *interface = ( MaterialModelMapperInterface * )
                                              this->giveMaterial()->giveInterface(MaterialModelMapperInterfaceType);

    if ( !interface ) {
        return 0;
    }

    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        iRule = integrationRulesArray [ i ];
        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
            result &= interface->MMI_map(iRule->getIntegrationPoint(j), oldd, tStep);
        }
    }

    return result;
}


int
Element :: adaptiveFinish(TimeStep *tStep)
{
    int i, j, result = 1;
    IntegrationRule *iRule;
    MaterialModelMapperInterface *interface = ( MaterialModelMapperInterface * )
                                              this->giveMaterial()->giveInterface(MaterialModelMapperInterfaceType);

    if ( !interface ) {
        return 0;
    }

    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        iRule = integrationRulesArray [ i ];
        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
            result &= interface->MMI_finish(tStep);
        }
    }

    return result;
}


void
Element :: updateLocalNumbering(EntityRenumberingFunctor &f)
{
    int i;
    for ( i = 1; i <= numberOfDofMans; i++ ) {
        dofManArray.at(i) = f(dofManArray.at(i), ERS_DofManager);
    }
}


#ifdef __PARALLEL_MODE
int
Element :: packUnknowns(CommunicationBuffer &buff, TimeStep *stepN)
{
    int i, j, result = 1;
    IntegrationRule *iRule;

    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        iRule = integrationRulesArray [ i ];
        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
            result &= this->giveCrossSection()->packUnknowns( buff, stepN, iRule->getIntegrationPoint(j) );
        }
    }

    return result;
}


int
Element :: unpackAndUpdateUnknowns(CommunicationBuffer &buff, TimeStep *stepN)
{
    int i, j, result = 1;
    IntegrationRule *iRule;

    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        iRule = integrationRulesArray [ i ];
        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
            result &= this->giveCrossSection()->unpackAndUpdateUnknowns( buff, stepN, iRule->getIntegrationPoint(j) );
        }
    }

    return result;
}


int
Element :: estimatePackSize(CommunicationBuffer &buff)
{
    int i, j, result = 0;
    IntegrationRule *iRule;

    for ( i = 0; i < numberOfIntegrationRules; i++ ) {
        iRule = integrationRulesArray [ i ];
        for ( j = 0; j < iRule->getNumberOfIntegrationPoints(); j++ ) {
            result += this->giveCrossSection()->estimatePackSize( buff, iRule->getIntegrationPoint(j) );
        }
    }

    return result;
}


double
Element :: predictRelativeComputationalCost()
{
    int j, nip;
    double wgt = 0;
    IntegrationRule *iRule = this->giveDefaultIntegrationRulePtr();
    nip = iRule->getNumberOfIntegrationPoints();
    for ( j = 0; j < nip; j++ ) {
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
