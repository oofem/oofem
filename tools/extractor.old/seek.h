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
#ifndef seek_h
#define seek_h

#include "tokenizer.h"

#define SEEK_TOL  1.e-4

enum elemRec { er_undefined, er_strain, er_stress, er_status };
struct stateType {
    double currStep;
    double currEigval;
    double currLoaLevel;
    int currNodeSide, currElement, currDof, currGP, currReactionRecFound;
    int currQuasiReactionRecFound;
    int currLoaLevelRecFound;
    elemRec currElementRec;
};

void resetStateExceptCurrStep();
void getCurrState(stateType *cs);
int seekSolutionStep(Tokenizer *t, double stepVal);
int seekNodeRecord(Tokenizer *t, int node);
int seekDofRecord(Tokenizer *t, int dof, char u, double &val);
int readDofUnknown(Tokenizer *t, int dof, char u, double &val);
int seekElementRecord(Tokenizer *t, int elem);
int seekGPRecord(Tokenizer *t, int gp);
int seekStressStrainGPRecord(Tokenizer *t, int ifstress, int ifsrain,
                             int stressstraincomp, double &val);
int seekBeamRecord(Tokenizer *t, int ifforce, int ifdispl,
                   int stressstraincomp, double &val);
int seekFirstReactionsRecord(Tokenizer *t);
int seekReactionRecord(Tokenizer *t, int node, int dof, double &val);
int seekEigvalSolutionStep(Tokenizer *t, double stepVal, double &eigVal);
int seekNLSolutionStep(Tokenizer *t, double stepVal, double &loadlevel, double &nite);
int seekAndParseMaterialStatusRecord(Tokenizer *t, char *keyword, int valIndex, double &val);
int seekFirstQuasiReactionRecord(Tokenizer *t);
int seekQuasiReactionRecord(Tokenizer *t, int node, int dof, double &val);
int seekStepUserTime(Tokenizer *t, double &utime);

#endif // seek_h
