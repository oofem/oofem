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
 *               Copyright (C) 1993 - 2025   Borek Patzak
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

#include "seek.h"
#include <math.h>

/*
 * following functions work only for forward direction.
 * It is no possible to seek in backward direction.
 * (because we assume also stdin redirection.)
 * The user must ensure the proper order of requests.
 */


static stateType currState = {
    -1., -1., -1., -1, -1, -1, -1, -1, -1, -1, er_undefined
};

void resetStateExceptCurrStep()
{
    /* currState.currEigval = -1;*/
    currState.currNodeSide = -1;
    currState.currElement = -1;
    currState.currDof = -1;
    currState.currGP = -1;
    currState.currReactionRecFound = -1;
    currState.currQuasiReactionRecFound = -1;
    currState.currLoaLevelRecFound = -1;
    currState.currElementRec = er_undefined;
}


void getCurrState(stateType *cs)
{
    cs->currStep = currState.currStep;
    cs->currEigval = currState.currEigval;
    cs->currLoaLevel = currState.currLoaLevel;
    cs->currNodeSide = currState.currNodeSide;
    cs->currElement = currState.currElement;
    cs->currDof = currState.currDof;
    cs->currGP = currState.currGP;
    cs->currReactionRecFound = currState.currReactionRecFound;
    cs->currLoaLevelRecFound = currState.currLoaLevelRecFound;
    cs->currQuasiReactionRecFound = currState.currQuasiReactionRecFound;
}

int seekSolutionStep(Tokenizer *t, double stepVal)
{
    /* seeks records in stream connected with tokezer t,
     * until record with requested stepVal is found.
     * if succesfull nonzero is returned, otherwise
     * zero return value indicates no such record.
     *
     * the stepVal == -1 means any next time step
     * it value can be extracted from currState.
     */
    double val;

    // check reached state
    if ( ( stepVal != -1 ) && ( fabs(currState.currStep - stepVal) < SEEK_TOL ) ) {
        return 1;
    }

    while ( !t->isEOF() ) {
        if ( t->giveNumberOfTokens() > 3 ) {
            if ( !strncmp(t->giveToken(1), "Output", 6) && !strncmp(t->giveToken(3), "time", 4) ) {
                sscanf(t->giveToken(4), "%lf", & val);
                if ( ( stepVal == -1 ) || ( fabs(val - stepVal) < SEEK_TOL ) ) {
                    resetStateExceptCurrStep();
                    currState.currStep = val;
                    return 1;
                }

                //   } else if (!strncmp(t->giveToken(1),"Reached",6) && !strncmp(t->giveToken(3),"level",5)) {
                //    sscanf(t->giveToken(5),"%lf",&val);
                //    if (fabs(val-stepVal) < 1.e-2) return 1;
            }
        }

        t->giveLineFromInput();
    }

    return 0;
}


int seekEigvalSolutionStep(Tokenizer *t, double stepVal, double &eigVal)
{
    /* seeks Eigenvalue record  in stream connected with tokezer t,
     * until record with requested stepVal Eigenvalue is found.
     * if succesfull nonzero is returned, otherwise
     * zero return value indicates no such record.
     * This function also seeks proper solution step.
     *
     * the stepVal == -1 means any next time step
     * it value can be extracted from currState.
     */
    double val;

    while ( !t->isEOF() ) {
        if ( ( stepVal != -1 ) && ( fabs(currState.currEigval - stepVal) < SEEK_TOL ) ) {
            return 1;
        }

        if ( t->giveNumberOfTokens() > 3 ) {
            if ( !strncmp(t->giveToken(1), "Output", 6) && !strncmp(t->giveToken(3), "eigen", 5) ) {
                sscanf(t->giveToken(6), "%lf", & val);
                if ( ( stepVal == -1 ) || ( fabs(val - stepVal) < 1.e-2 ) ) {
                    // move to next record
                    t->giveLineFromInput();
                    if ( t->giveNumberOfTokens() < 10 ) {
                        return 0;
                    }

                    if ( !strncmp(t->giveToken(7), "eigen", 5) && !strncmp(t->giveToken(8), "value", 5) ) {
                        sscanf(t->giveToken(10), "%lf", & eigVal);
                        resetStateExceptCurrStep();
                        currState.currEigval = val;
                        return 1;
                    } else {
                        return 0; // no Eigval record
                    }
                }
            }
        }

        t->giveLineFromInput();
    }

    return 0;
}



int seekNLSolutionStep(Tokenizer *t, double stepVal, double &loadlevel, double &nite)
{
    /* seeks Nonlinear solution  record  in stream connected with tokezer t,
     * until record with requested stepVal  is found.
     * if succesfull nonzero is returned, otherwise
     * zero return value indicates no such record.
     * On success return value val contains reached load level.
     * The proper solution step must be seeked
     */

    if ( currState.currLoaLevelRecFound != 1 ) {
        while ( !t->isEOF() ) {
            if ( t->giveNumberOfTokens() >= 5 ) {
                if ( !strncmp(t->giveToken(1), "Reached", 7) && !strncmp(t->giveToken(3), "level", 5) ) {
                    currState.currLoaLevelRecFound = 1;
                    sscanf(t->giveToken(5), "%lf", & loadlevel);
                    //resetStateExceptCurrStep();
                    currState.currLoaLevel = loadlevel;
                    sscanf(t->giveToken(7), "%lf", & nite);
                    return 1;
                }
            }

            t->giveLineFromInput();
        }
    } else {
        if ( !strncmp(t->giveToken(1), "Reached", 7) && !strncmp(t->giveToken(3), "level", 5) ) {
            sscanf(t->giveToken(5), "%lf", & loadlevel);
            //resetStateExceptCurrStep();
            currState.currLoaLevel = loadlevel;
            sscanf(t->giveToken(7), "%lf", & nite);
            return 1;
        } else {
            return 0;
        }
    }
    return 0;
}


int seekFirstQuasiReactionRecord(Tokenizer *t)
{
    // seeks first quasi-reaction record in stream connected with tokenizer t,
    //  until first quasi-rection record is found.
    //
    // Before this function call a proper solution step must be seeked.
    // If successfull nonzero is returned, zero otherwise.
    //
    if ( currState.currQuasiReactionRecFound == 1 ) {
        return 1;
    }

    while ( !t->isEOF() ) {
        if ( t->giveNumberOfTokens() == 3 ) {
            if ( !strncmp(t->giveToken(1), "Quasi", 5) && !strncmp(t->giveToken(3), "table:", 6) ) { // header found
                // skip 3 lines and check for separator
                for ( int i = 0; i < 3; t->giveLineFromInput(), i++ ) {
                    ;
                }

                // check separator
                if ( !strncmp(t->giveToken(1), "=", 1) ) {
                    currState.currQuasiReactionRecFound = 1;
                    t->giveLineFromInput();
                    return 1;
                }

                return 0;
            }
        }

        t->giveLineFromInput();
    }

    return 0;
}

int seekQuasiReactionRecord(Tokenizer *t, int node, int dof, double &val)
{
    // seeks Node quasi-reaction  record in stream connected with tokezer t,
    // until record for requested node and dof  is found
    //
    // Before this function call , a proper solution step and first reaction record must be seeked.
    // if succesfull nonzero is returned, otherwise
    // zero return value indicates no such record.

    int n, d;

    while ( !t->isEOF() ) {
        if ( t->giveNumberOfTokens() != 3 ) {
            return 0;
        }

        sscanf(t->giveToken(1), "%d", & n);
        sscanf(t->giveToken(2), "%d", & d);
        if ( ( n == node ) && ( d == dof ) ) {
            sscanf(t->giveToken(3), "%lf", & val);
            return 1;
        }

        t->giveLineFromInput();
    }

    return 0;
}




int seekNodeRecord(Tokenizer *t, int node)
{
    /* seeks Node record in stream connected with tokezer t,
     * until record with requested Node is found
     *
     * Before this function call , a proper solution step must be seeked.
     * if succesfull nonzero is returned, otherwise
     * zero return value indicates no such record.
     */
    int n;

    // check current state
    if ( currState.currNodeSide == node ) {
        return 1;
    }

    while ( !t->isEOF() ) {
        if ( t->giveNumberOfTokens() > 1 ) {
            if ( !strncmp(t->giveToken(1), "Node", 4) || !strncmp(t->giveToken(1), "Side", 4) ||
                !strncmp(t->giveToken(1), "RigidArmNode", 12) || !strncmp(t->giveToken(1), "HangingNode", 11) ) {
                sscanf(t->giveToken(2), "%d", & n);
                if ( node == n ) {
                    resetStateExceptCurrStep();
                    ;
                    currState.currNodeSide = node;
                    t->giveLineFromInput(); // skip one line to move to dof reccords
                    return 1; // node record found
                }
            }
        }

        t->giveLineFromInput();
    }

    return 0;
}

int seekDofRecord(Tokenizer *t, int dof, char u, double &val)
{
    /* seeks Dof record in stream connected with tokezer t,
     * until record with requested Dof is found
     *
     * Before this function call , a proper node record must be seeked.
     * if succesfull nonzero is returned, otherwise
     * zero return value indicates no such record.
     * return parameter val is (if success) corresponding unknown value.
     */
    int d;
    int pos;

    if ( currState.currDof == dof ) {
        return readDofUnknown(t, dof, u, val);
    }

    while ( !t->isEOF() ) {
        if ( t->giveNumberOfTokens() > 3 ) {
            if ( !strncmp(t->giveToken(1), "Node", 4) ) {
                //    pos = 3;
                //    // check if we met node record again - then error
                //    if (init == 0) return 0; else init  = 0;
                return 0;
            } else {
                pos = 1;
            }

            if ( !strncmp(t->giveToken(pos), "dof", 3) ) {
                sscanf(t->giveToken(pos + 1), "%d", & d);
                if ( d == dof ) { // dof record found - now let us find unknown
                    currState.currDof = dof; // do not reset - currNode must be properly set
                    return readDofUnknown(t, dof, u, val);
                }
            } else {
                // parsed record does not belog to node record set
                // => no such dof in node exists
                return 0;
            }
        }

        t->giveLineFromInput();
    }

    return 0;
}


int readDofUnknown(Tokenizer *t, int dof, char u, double &val)
{
    /* reads Dof unknown from currLine of stream connected with tokezer t,
     *
     * Before this function call , a proper node and dof records must be seeked.
     * if succesfull nonzero is returned, otherwise
     * zero return value indicates no such record.
     * return parameter val is (if success) corresponding unknown value.
     */
    int d;
    int pos;

    if ( t->giveNumberOfTokens() > 3 ) {
        if ( !strncmp(t->giveToken(1), "Node", 4) ) {
            pos = 3;
            // check if we met node record again - then error
        } else {
            pos = 1;
        }

        if ( !strncmp(t->giveToken(pos), "dof", 3) ) {
            sscanf(t->giveToken(pos + 1), "%d", & d);
            if ( d == dof ) { // dof record found - now let us find unknown
                int pp = pos + 2;
                while ( pp < t->giveNumberOfTokens() ) {
                    if ( t->giveToken(pp) [ 0 ] == u ) { // first character  - unknown type
                        sscanf(t->giveToken(pp + 1), "%lf", & val);
                        return 1;
                    } else {
                        pp += 2;
                    }
                }
            } else {
                // dof record exist, but no such unknown
                return 0;
            }
        }

        // parsed record does not belog to node record set
        // => no such dof in node exists
        return 0;
    }

    return 0;
}


int seekElementRecord(Tokenizer *t, int elem)
{
    /* seeks Element record in stream connected with tokezer t,
     * until record with requested Element  is encountered.
     *
     * Before this function call , a proper solution step must be seeked.
     * if succesfull nonzero is returned, otherwise
     * zero return value indicates no such record.
     */
    int e;

    if ( currState.currElement == elem ) {
        return 1;
    }

    while ( !t->isEOF() ) {
        if ( t->giveNumberOfTokens() >= 2 ) {
            if ( !strncmp(t->giveToken(1), "element", 7) ) {
                sscanf(t->giveToken(2), "%d", & e);
                if ( elem == e ) { // element record found
                    resetStateExceptCurrStep();
                    currState.currElement = elem;
                    t->giveLineFromInput(); // skip one line to move to GP reccords
                    return 1;
                }
            }
        }

        t->giveLineFromInput();
    }

    return 0;
}

int seekGPRecord(Tokenizer *t, int gp)
{
    /* seeks GP record in stream connected with tokezer t,
     * until record with requested  gp is encountered.
     *
     * Before this function call , a proper solution step and element must be seeked.
     * if succesfull nonzero is returned, otherwise
     * zero return value indicates no such record.
     */
    int g;

    if ( currState.currGP == gp ) {
        return 1;
    }

    while ( !t->isEOF() ) {
        if ( !strncmp(t->giveToken(1), "GP", 2) ) {
            sscanf(t->giveToken(2), "%d", & g);
            if ( g == gp ) { // gp record found
                currState.currGP = gp;
                currState.currElementRec = er_strain;
                return 1;
            }
        } else if ( !strncmp(t->giveToken(1), "element", 7) ) {
            return 0;                                        // next element rec. found
        }

        t->giveLineFromInput();
    }

    return 0;
}

int seekStressStrainGPRecord(Tokenizer *t, int ifstress, int ifstrain,
                             int stressstraincomp, double &val)
{
    /* seeks GP record in stream connected with tokezer t,
     * until record with requested data  is encountered.
     *
     * Before this function call , a proper GP and element  must be seeked.
     * if succesfull nonzero is returned, otherwise
     * zero return value indicates no such record.
     * return parameter value val contains unknown value of unknown.
     */
    int pos = 1;

    if ( ifstress && ( currState.currElementRec == er_strain ) ) {
        t->giveLineFromInput(); // skip one line
        currState.currElementRec = er_stress;
    }

    // if (ifstatus) {t->giveLineFromInput(); t->giveLineFromInput();} // skip two lines

    if ( ifstrain ) {
        // check proper record
        if ( !( pos = t->hasToken("strains") ) ) {
            return 0;
        }
    } else if ( ifstress ) {
        if ( !( pos = t->hasToken("stresses") ) ) {
            return 0;
        }
    } /* else {
       * if (!strncmp(t->giveToken(1),"status",6)) return 0;
       * } */

    if ( ifstrain || ifstress ) { // read directly value
        if ( ( pos + stressstraincomp ) > t->giveNumberOfTokens() ) {
            return 0;                                               // no such component exists
        }

        sscanf(t->giveToken(pos + stressstraincomp), "%lf", & val);
        return 1;
    } /*else {
       * // status record
       *
       * }*/

    return 0;
}


int seekAndParseMaterialStatusRecord(Tokenizer *t, char *keyword, int valIndex, double &val)
{
    /* seeks GP status record in stream connected with tokezer t,
     * until record with requested data  is encountered.
     *
     * Before this function call , a proper GP and element  must be seeked.
     * if succesfull nonzero is returned, otherwise
     * zero return value indicates no such record.
     * return parameter value val contains unknown value of unknown.
     */
    int result = 0;

    if ( currState.currElementRec != er_status ) {
        while ( !t->isEOF() ) {
            if ( t->giveNumberOfTokens() >= 2 ) {
                if ( !strncmp(t->giveToken(1), "status", 6) ) {
                    // status found
                    currState.currElementRec = er_status;
                    break;
                }
            }

            t->giveLineFromInput();
        }
    }

    // seek keyword
    for ( int i = 1; i <= t->giveNumberOfTokens(); i++ ) {
        //if ( !strncmp(t->giveToken(i), keyword, strlen(keyword) - 1) ) {
        if ( !strncmp( t->giveToken(i), keyword, strlen(keyword) ) ) {
            // keyword found
            if ( i + valIndex > t->giveNumberOfTokens() ) {
                return 0;
            }

            sscanf(t->giveToken(i + valIndex), "%lf", & val);
            result = 1;
            break;
        }
    }

    // if status entry not found do not generate error and set value to zero
    if ( result == 0 ) {
        result = 1;
        val = 0.0;
    }

    return result;
}


int seekBeamRecord(Tokenizer *t, int ifforce, int ifdispl,
                   int stressstraincomp, double &val)
{
    /* seeks GP record in stream connected with tokezer t,
     * until record with requested data  is encountered.
     *
     * Before this function call , a proper GP and element  must be seeked.
     * if succesfull nonzero is returned, otherwise
     * zero return value indicates no such record.
     * return parameter value val contains unknown value of unknown.
     */
    int pos = 1;

    if ( ifforce && ( currState.currElementRec != er_stress ) ) {
        // if (ifforce && (currState.currStress != 1)) {
        t->giveLineFromInput(); // skip one line
        currState.currElementRec = er_stress;
    }

    // if (ifstatus) {t->giveLineFromInput(); t->giveLineFromInput();} // skip two lines

    if ( ifdispl ) {
        // check proper record
        if ( strncmp(t->giveToken(2), "displacements", 13) ) {
            return 0;
        }

        pos = 2;
    } else if ( ifforce ) {
        if ( strncmp(t->giveToken(3), "forces", 6) ) {
            return 0;
        }

        pos = 3;
    } /* else {
       * if (!strncmp(t->giveToken(1),"status",6)) return 0;
       * } */

    if ( ifforce || ifdispl ) { // read directly value
        if ( ( pos + stressstraincomp ) > t->giveNumberOfTokens() ) {
            return 0;                                               // no such component exists
        }

        sscanf(t->giveToken(pos + stressstraincomp), "%lf", & val);
        return 1;
    } /*else {
       * // status record
       *
       * }*/

    return 0;
}


int seekFirstReactionsRecord(Tokenizer *t)
{
    /* seeks first reaction  record in stream connected with tokezer t,
     * until record with requested first reactions record is found
     *
     * Before this function call , a proper solution step must be seeked.
     * if succesfull nonzero is returned, otherwise
     * zero return value indicates no such record.
     */
    if ( currState.currReactionRecFound == 1 ) {
        return 1;
    }

    while ( !t->isEOF() ) {
        if ( t->giveNumberOfTokens() > 3 ) {
            if ( !strncmp(t->giveToken(1), "R", 1) && !strncmp(t->giveToken(2), "E", 1) &&
                !strncmp(t->giveToken(3), "A", 1) ) { // header found
                // skip 3 lines
                for ( int i = 0; i < 4; t->giveLineFromInput(), i++ ) {
                    ;
                }

                // check
                if ( !strncmp(t->giveToken(5), "reaction", 8) ) {
                    currState.currReactionRecFound = 1;
                    return 1;
                }

                return 0;
            }
        }

        t->giveLineFromInput();
    }

    return 0;
}


int seekReactionRecord(Tokenizer *t, int node, int dof, double &val)
{
    /* seeks Node reaction  record in stream connected with tokezer t,
     * until record for requested node and dof  is found
     *
     * Before this function call , a proper solution step and first reaction record must be seeked.
     * if succesfull nonzero is returned, otherwise
     * zero return value indicates no such record.
     */
    int n, d;

    while ( !t->isEOF() ) {
        if ( t->giveNumberOfTokens() < 6 ) {
            return 0;
        }

        sscanf(t->giveToken(2), "%d", & n);
        sscanf(t->giveToken(4), "%d", & d);
        if ( ( n == node ) && ( d == dof ) ) {
            sscanf(t->giveToken(6), "%lf", & val);
            return 1;
        }

        t->giveLineFromInput();
    }

    return 0;
}


int seekStepUserTime(Tokenizer *t, double &utime)
{
    //if (currState.currQuasiReactionRecFound == 1) return 1;
    while ( !t->isEOF() ) {
        if ( t->giveNumberOfTokens() == 9 ) {
            if ( !strncmp(t->giveToken(1), "User", 4) && !strncmp(t->giveToken(6), "step", 4) ) { // header found
                sscanf(t->giveToken(8), "%lf", & utime);
                return 1;
            }
        }

        t->giveLineFromInput();
    }

    return 0;
}















