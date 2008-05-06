/* $Header: /home/cvs/bp/oofem/oofemlib/src/rigidarmslavedof.C,v 1.9 2003/04/06 14:08:25 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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


//   file SLAVEDOF.CC

#include "rigidarmslavedof.h"
#include "rigidarmnode.h"
#include "dofmanager.h"
#include "domain.h"
#include "timestep.h"
#include "boundary.h"
#include "initial.h"

#include "flotarry.h"
#include "dictionr.h"

#include "debug.h"

#ifndef __MAKEDEPEND
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#endif


RigidArmSlaveDof :: RigidArmSlaveDof(int i, DofManager *aNode, DofID id) : Dof(i, aNode, id)
{ }

/*
 * BoundaryCondition*  RigidArmSlaveDof :: giveBc ()
 * // Returns the boundary condition the receiver is subjected to.
 * {
 * to do
 * //return  this->giveMasterDof()->giveBc();
 * }
 *
 * InitialCondition*  RigidArmSlaveDof :: giveIc ()
 * // Returns the initial condition on the receiver. Not used.
 * {
 * to do
 * //return this->giveMasterDof()->giveIc ();
 * }
 */

double RigidArmSlaveDof :: giveUnknown(EquationID type, ValueModeType mode, TimeStep *stepN)
// The key method of class Dof. Returns the value of the unknown 'u'
// (e.g., the displacement) of the receiver, at stepN. This value may,
// or may not, be already available. It may depend on a boundary (if it
// is not a predicted unknown) or initial condition. stepN is not the
// current time step n, it is assumed to be the previous one (n-1).
{
    FloatMatrix t;
    DofManager *master = ( ( RigidArmNode * ) ( dofManager ) )->giveMasterDofMngr();
    FloatArray answer, masterUnknwns( master->giveNumberOfDofs() );

    return master->giveDof(this->number)->giveUnknown(type, mode, stepN);

    /*
     * this->computeTransformation (t, _toSlave);
     * for (i=1; i<=master->giveNumberOfDofs(); i++) {
     * if (this->giveUnknownType () == master->giveDof(i)->giveUnknownType())
     * masterUnknwns.at(i) = master->giveDof(i)->giveUnknown (this->giveUnknownType(), mode, stepN);
     * }
     *
     * answer.beProductOf (t,masterUnknwns);
     * return answer.at(1);
     */
}


double RigidArmSlaveDof :: giveUnknown(PrimaryField &field, ValueModeType mode, TimeStep *stepN)
{
    FloatMatrix t;
    DofManager *master = ( ( RigidArmNode * ) ( dofManager ) )->giveMasterDofMngr();
    FloatArray answer, masterUnknwns( master->giveNumberOfDofs() );

    return master->giveDof(this->number)->giveUnknown(field, mode, stepN);
}


double RigidArmSlaveDof :: giveLocalUnknown(EquationID type, ValueModeType mode, TimeStep *stepN)
// The key method of class Dof. Returns the value of the unknown 'u'
// (e.g., the displacement) of the receiver, at stepN. This value may,
// or may not, be already available. It may depend on a boundary (if it
// is not a predicted unknown) or initial condition. stepN is not the
// current time step n, it is assumed to be the previous one (n-1).
{
    int i;
    FloatMatrix t;
    Node *master = ( ( RigidArmNode * ) ( dofManager ) )->giveMasterDofMngr();
    FloatArray answer, masterUnknwns( master->giveNumberOfDofs() );

    for ( i = 1; i <= master->giveNumberOfDofs(); i++ ) {
        //if (this->giveUnknownType () == master->giveDof(i)->giveUnknownType())
        // masterUnknwns.at(i) = master->giveDof(i)->giveUnknown (this->giveUnknownType(), mode, stepN);
        masterUnknwns.at(i) = master->giveDof(i)->giveUnknown(type, mode, stepN);
    }

    master->computeGNDofTransformation(t, NULL);
    masterUnknwns.rotatedWith(t, 't'); // transform master unknowns to global c.s
    this->computeDofTransformation(t, _toSlave); // express unknowns in receiver space
    answer.beProductOf(t, masterUnknwns);

    return answer.at(1);
}



int RigidArmSlaveDof :: hasBc(TimeStep *tStep)
// Returns True if the receiver is subjected to a boundary condition, else
// returns False. If necessary, reads the answer in the data file.
{
    // RigidArmSlaveDof can not be subjected to bc
    // it is only mapping to master
    return 0;
}


int RigidArmSlaveDof :: hasIc()
// Returns True if the receiver is subjected to an initial condition,
// else returns False.
{
    // RigidArmSlaveDof can not be subjected to ic
    // it is only mapping to master
    return 0;
}


int RigidArmSlaveDof :: hasIcOn(ValueModeType u)
// Returns True if the unknown 'u' (e.g., the displacement 'd') of the
// receiver is subjected to an initial condition, else returns False.
{
    // RigidArmSlaveDof can not be subjected to ic
    // it is only mapping to master
    return 0;
}

int RigidArmSlaveDof :: giveBcIdValue()
{
    return 0;
}

contextIOResultType RigidArmSlaveDof :: saveContext(FILE *stream, void *obj)
//
// saves full node context (saves state variables, that completely describe
// current state)
//
{
    return CIO_OK;
}


contextIOResultType RigidArmSlaveDof :: restoreContext(FILE *stream, void *obj)
//
// restores full node context (saves state variables, that completely describe
// current state)
//
{
    return CIO_OK;
}


void
RigidArmSlaveDof :: computeDofTransformation(FloatMatrix &answer, RigidArmSlaveDofTransfType ttype)
{
    // computes trasformation of local slave dof to ALL master dofs (in Global cs).
    // the non-relevant are included, because of a lot of decisons to be done,
    // which will probably lead to reduced efficiency.
    Node *from, *to;
    RigidArmNode *owner  = ( RigidArmNode * ) dofManager;
    DofManager *master   = owner->giveMasterDofMngr();
    Node *masterNode = dynamic_cast< Node * >( master );
    int indx, masterdofs;

    if ( masterNode == NULL ) {
        _error("computeTransformation: master must be Node");
    }

    if ( ttype == _toMaster ) {
        to    =  masterNode;
        from  =  ( Node * ) dofManager;
    } else {
        to    = ( Node * ) dofManager;
        from  = masterNode;
    }

    masterdofs = master->giveNumberOfDofs();
    answer.resize(1, masterdofs);
    answer.zero();

    // compute transformation due to rigid arm between master and receiver in Global CS.
    if ( ( this->dofID == D_u ) || ( this->dofID == D_v ) || ( this->dofID == D_w ) ) {
        // special class to represent rigidly connected dofMAnager ??
        if ( ( indx = master->findDofWithDofId(this->dofID) ) ) {
            answer.at(1, indx) = 1.0;
        }

        if ( this->dofID == D_u ) {
            if ( ( indx = master->findDofWithDofId(R_v) ) ) {
                if ( owner->giveSlaveDofMask(indx) ) {
                    answer.at(1, indx) =  ( to->giveCoordinate(3) - from->giveCoordinate(3) );
                }
            }

            if ( ( indx = master->findDofWithDofId(R_w) ) ) {
                if ( owner->giveSlaveDofMask(indx) ) {
                    answer.at(1, indx) = -( to->giveCoordinate(2) - from->giveCoordinate(2) );
                }
            }
        } else if ( this->dofID == D_v ) {
            if ( ( indx = master->findDofWithDofId(R_u) ) ) {
                if ( owner->giveSlaveDofMask(indx) ) {
                    answer.at(1, indx) = -( to->giveCoordinate(3) - from->giveCoordinate(3) );
                }
            }

            if ( ( indx = master->findDofWithDofId(R_w) ) ) {
                if ( owner->giveSlaveDofMask(indx) ) {
                    answer.at(1, indx) =  ( to->giveCoordinate(1) - from->giveCoordinate(1) );
                }
            }
        } else { // if (this->dofID == D_w)
            if ( ( indx = master->findDofWithDofId(R_u) ) ) {
                if ( owner->giveSlaveDofMask(indx) ) {
                    answer.at(1, indx) =  ( to->giveCoordinate(2) - from->giveCoordinate(2) );
                }
            }

            if ( ( indx = master->findDofWithDofId(R_v) ) ) {
                if ( owner->giveSlaveDofMask(indx) ) {
                    answer.at(1, indx) = -( to->giveCoordinate(1) - from->giveCoordinate(1) );
                }
            }
        }

        // other components are zero
    } else  if ( ( this->dofID == R_u ) || ( this->dofID == R_v ) || ( this->dofID == R_w ) ) {
        // special class to represent rigidly connected dofMAnager ??
        if ( ( indx = master->findDofWithDofId(this->dofID) ) ) {
            answer.at(1, indx) = 1.0;
        }

        /*
         * if (this->dofID == R_u) {
         * if (indx = master->findDofWithDofId (D_v)) answer.at(1,indx) =  (to->giveCoordinate(3)-from->giveCoordinate(3));
         * if (indx = master->findDofWithDofId (D_w)) answer.at(1,indx) = -(to->giveCoordinate(2)-from->giveCoordinate(2));
         * } else if (this->dofID == R_v) {
         * if (indx = master->findDofWithDofId (D_u)) answer.at(1,indx) = -(to->giveCoordinate(3)-from->giveCoordinate(3));
         * if (indx = master->findDofWithDofId (D_w)) answer.at(1,indx) =  (to->giveCoordinate(1)-from->giveCoordinate(1));
         * } else { // if (this->dofID == D_w)
         * if (indx = master->findDofWithDofId (D_u)) answer.at(1,indx) =  (to->giveCoordinate(2)-from->giveCoordinate(2));
         * if (indx = master->findDofWithDofId (D_v)) answer.at(1,indx) = -(to->giveCoordinate(1)-from->giveCoordinate(1));
         * }
         */
        // other components are zero
    } else {
        if ( ( indx = master->findDofWithDofId(this->dofID) ) ) {
            answer.at(1, indx) = 1.0;
        } else {
            _error("computeTransformation: incoppatible dofs of slave-master");
        }
    }
}




void
RigidArmSlaveDof :: computeLoadTransformation(FloatMatrix &answer, RigidArmSlaveDofTransfType ttype)
{
    // computes trasformation of local slave dof to ALL master dofs (in Global cs).
    // the non-relevant are included, because of a lot of decisons to be done,
    // which will probably lead to reduced efficiency.
    Node *from, *to;
    RigidArmNode *owner  = ( RigidArmNode * ) dofManager;
    DofManager *master   = owner->giveMasterDofMngr();
    Node *masterNode = dynamic_cast< Node * >( master );
    int indx, masterdofs;

    if ( masterNode == NULL ) {
        _error("computeTransformation: master must be Node");
    }

    if ( ttype == _toMaster ) {
        to    =  masterNode;
        from  =  ( Node * ) dofManager;
    } else {
        to    = ( Node * ) dofManager;
        from  = masterNode;
    }

    masterdofs = master->giveNumberOfDofs();
    answer.resize(1, masterdofs);
    answer.zero();

    // compute transformation due to rigid arm between master and receiver in Global CS.
    if ( ( this->dofID == D_u ) || ( this->dofID == D_v ) || ( this->dofID == D_w ) ) {
        // special class to represent rigidly connected dofMAnager ??
        if ( ( indx = master->findDofWithDofId(this->dofID) ) ) {
            answer.at(1, indx) = 1.0;
        }

        /*
         * if (this->dofID == D_u) {
         * if (indx = master->findDofWithDofId (R_v)) answer.at(1,indx) = -(to->giveCoordinate(3)-from->giveCoordinate(3));
         * if (indx = master->findDofWithDofId (R_w)) answer.at(1,indx) =  (to->giveCoordinate(2)-from->giveCoordinate(2));
         * } else if (this->dofID == D_v) {
         * if (indx = master->findDofWithDofId (R_u)) answer.at(1,indx) =  (to->giveCoordinate(3)-from->giveCoordinate(3));
         * if (indx = master->findDofWithDofId (R_w)) answer.at(1,indx) = -(to->giveCoordinate(1)-from->giveCoordinate(1));
         * } else { // if (this->dofID == D_w)
         * if (indx = master->findDofWithDofId (R_u)) answer.at(1,indx) = -(to->giveCoordinate(2)-from->giveCoordinate(2));
         * if (indx = master->findDofWithDofId (R_v)) answer.at(1,indx) =  (to->giveCoordinate(1)-from->giveCoordinate(1));
         * }
         */

        // other components are zero
    } else  if ( ( this->dofID == R_u ) || ( this->dofID == R_v ) || ( this->dofID == R_w ) ) {
        // special class to represent rigidly connected dofMAnager ??
        if ( ( indx = master->findDofWithDofId(this->dofID) ) ) {
            answer.at(1, indx) = 1.0;
        }

        if ( this->dofID == R_u ) {
            if ( ( indx = master->findDofWithDofId(D_v) ) ) {
                if ( owner->giveSlaveDofMask(indx) ) {
                    answer.at(1, indx) =  ( to->giveCoordinate(3) - from->giveCoordinate(3) );
                }
            }

            if ( ( indx = master->findDofWithDofId(D_w) ) ) {
                if ( owner->giveSlaveDofMask(indx) ) {
                    answer.at(1, indx) = -( to->giveCoordinate(2) - from->giveCoordinate(2) );
                }
            }
        } else if ( this->dofID == R_v ) {
            if ( ( indx = master->findDofWithDofId(D_u) ) ) {
                if ( owner->giveSlaveDofMask(indx) ) {
                    answer.at(1, indx) = -( to->giveCoordinate(3) - from->giveCoordinate(3) );
                }
            }

            if ( ( indx = master->findDofWithDofId(D_w) ) ) {
                if ( owner->giveSlaveDofMask(indx) ) {
                    answer.at(1, indx) =  ( to->giveCoordinate(1) - from->giveCoordinate(1) );
                }
            }
        } else { // if (this->dofID == D_w)
            if ( ( indx = master->findDofWithDofId(D_u) ) ) {
                if ( owner->giveSlaveDofMask(indx) ) {
                    answer.at(1, indx) =  ( to->giveCoordinate(2) - from->giveCoordinate(2) );
                }
            }

            if ( ( indx = master->findDofWithDofId(D_v) ) ) {
                if ( owner->giveSlaveDofMask(indx) ) {
                    answer.at(1, indx) = -( to->giveCoordinate(1) - from->giveCoordinate(1) );
                }
            }
        }

        // other components are zero
    } else {
        if ( ( indx = master->findDofWithDofId(this->dofID) ) ) {
            answer.at(1, indx) = 1.0;
        } else {
            _error("computeTransformation: incoppatible dofs of slave-master");
        }
    }
}


void RigidArmSlaveDof :: printSingleOutputAt(FILE *File, TimeStep *stepN, char ch,
                                             EquationID type, ValueModeType mode, double scale)
// Prints in the data file the unknown 'u' (for example, the displacement
// 'd') of the receiver, at stepN.
{
    double x;

    x    = scale * this->giveLocalUnknown(type, mode, stepN);
    fprintf(File, "  dof %d   %c % .8e\n", number, ch, x);
}



void RigidArmSlaveDof :: printMultipleOutputAt(FILE *File, TimeStep *stepN, char *ch,
                                               EquationID type, IntArray &mode, int nite)
{
    int i;
    double x;

    fprintf(File, "  dof %d", number);
    for ( i = 1; i <= nite; i++ ) {
        x    = this->giveLocalUnknown(type, ( ValueModeType ) mode.at(i), stepN);
        fprintf(File, "   %c % .8e", ch [ i - 1 ], x);
    }

    fprintf(File, "\n");
}
