/* $Header: /home/cvs/bp/oofem/sm/src/lspace.h,v 1.6 2003/05/19 13:04:00 bp Exp $ */
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


//   ***************************************
//   *** Macro Linear Hexahedral element ***
//   ***************************************

#ifndef macrolspace_h
#define macrolspace_h

#include "lspace.h"
#include "sparsemtrx.h"
#include "engngm.h"
//#include "nlstructuralelement.h"
//#include "fei3dhexalin.h"

class MacroLSpace : public LSpace
{
/*
This class implements a macroelement. It is derived from eight-node brick element. The stiffness matrix is computed from underlying RVE and is condensed to 24 DoFs to corner nodes.
*/

protected:

public:
 MacroLSpace (int,Domain*) ;                      // constructor
 ~MacroLSpace ()  {} ;                            // destructor

 const char* giveClassName () const { return "MacroLSpace" ;}

 virtual void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
 //friend class EngngModel;//access to assemble()
 



protected:
  SparseMtrx *stiffnessMatrixMicro;
  
};




#endif //macrolspace_h
