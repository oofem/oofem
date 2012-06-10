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

#include "dof.h"
#include "dofmanager.h"
#include "timestep.h"
#include "boundary.h"
#include "initial.h"
#include "datastream.h"
#include "oofem_limits.h"
#include "contextioerr.h"
#include "unknownnumberingscheme.h"

#ifndef __MAKEDEPEND
 #include <cstring>
 #include <cstdlib>
 #include <cstdarg>
#endif

namespace oofem {
Dof :: Dof(int i, DofManager *aNode, DofIDItem id)
// Constructor. Creates a new d.o.f., with number i, belonging
// to aNode
{
    number         = i;
    dofManager     = aNode;
    dofID          = id;
}

int Dof :: giveEquationNumber(const UnknownNumberingScheme &s) {
    return s.giveDofEquationNumber(this);
}

void Dof :: giveEquationNumbers(IntArray &masterEqNumbers, const UnknownNumberingScheme &s)
{
    masterEqNumbers.resize(1);
    masterEqNumbers.at(1) = s.giveDofEquationNumber(this);
}


int
Dof :: giveDofManNumber() const { return this->dofManager->giveNumber(); } // termitovo

#ifdef __PARALLEL_MODE
int
Dof :: giveDofManGlobalNumber() const { return this->dofManager->giveGlobalNumber(); }
#endif

void Dof :: printSingleOutputAt(FILE *File, TimeStep *stepN, char ch,
                                EquationID type, ValueModeType mode, double scale)
// Prints in the data file the unknown 'u' (for example, the displacement
// 'd') of the receiver, at stepN.
{
    double x = scale * this->giveUnknown(type, mode, stepN);
    fprintf(File, "  dof %d   %c % .8e\n", number, ch, x);
}



void Dof :: printMultipleOutputAt(FILE *File, TimeStep *stepN, char *ch,
                                  EquationID type, ValueModeType *mode, int nite)
// Prints in the data file the unknown 'u' (for example, the displacement
// 'd') of the receiver, at stepN.
{
    double x;

    fprintf(File, "  dof %d", number);
    for (int i = 1; i <= nite; i++ ) {
        x = this->giveUnknown(type, mode [ i - 1 ], stepN);
        fprintf(File, "   %c % .8e", ch [ i - 1 ], x);
    }

    fprintf(File, "\n");
}


void Dof :: printYourself()
// Prints the receiver on screen.
{
    printf( "dof %d  of %s %d :\n", number, dofManager->giveClassName(), dofManager->giveNumber() );
    // printf ("equation %d    bc %d \n",equationNumber,bc) ;

    // printOutputAt (node->giveDomain()->giveEngngModel()->giveCurrentStep());
}


void Dof :: error(const char *file, int line, const char *format, ...) const
{
    char buffer [ MAX_ERROR_MSG_LENGTH ];
    va_list args;

    va_start(args, format);
    vsprintf(buffer, format, args);
    va_end(args);

    __OOFEM_ERROR5(file, line, "Class: %s, number %d, of DofManager %d\n%s",
                   this->giveClassName(), number, dofManager->giveNumber(), buffer);
}

char *
Dof :: giveDofIDName(char *s)
{
    // returns character string representing receiver's dofID
    switch ( this->giveDofID() ) {
    case D_u:
        return strcpy(s, "D_u");

    case D_v:
        return strcpy(s, "D_v");

    case D_w:
        return strcpy(s, "D_w");

    case R_u:
        return strcpy(s, "R_u");

    case R_v:
        return strcpy(s, "R_v");

    case R_w:
        return strcpy(s, "R_w");

    case T_f:
        return strcpy(s, "T_f");

    case G_0:
        return strcpy(s, "G_0");

    case G_1:
        return strcpy(s, "G_1");

    default:
        sprintf(s, "%3d", this->number);
        return s;
    }

    // return s;
}

double
Dof :: giveBcValue(ValueModeType mode, TimeStep *tStep)
{
    if ( this->hasBc(tStep) ) {
        double rel = 0.0;
        if ( mode == VM_Incremental && tStep->isTheFirstStep() && hasIcOn(VM_Total) ) {
            rel = giveIc()->give(VM_Total);
        }

        return this->giveBc()->give(this, mode, tStep) - rel;
    } else {
        return 0.0;
    }
}

contextIOResultType
Dof :: saveContext(DataStream *stream, ContextMode mode, void *obj)
{
    if ( stream == NULL ) {
        THROW_CIOERR(CIO_IOERR);
    }

    // store dofid
    int _val = dofID;
    if ( !stream->write(& _val, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }


    if ( mode & CM_Definition ) {
        if ( !stream->write(& number, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }

    return CIO_OK;
}


contextIOResultType
Dof :: restoreContext(DataStream *stream, ContextMode mode, void *obj)
{
    // restore dofid
    int _val;
    if ( !stream->read(& _val, 1) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    dofID = ( DofIDItem ) _val;


    if ( mode & CM_Definition ) {
        if ( !stream->read(& number, 1) ) {
            THROW_CIOERR(CIO_IOERR);
        }
    }

    return CIO_OK;
}

} // end namespace oofem
