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

#include "microplane.h"

namespace oofem {
Microplane :: Microplane(IntegrationRule *ir, int n, MaterialMode mode) :
    GaussPoint(ir, n, FloatArray(), 0., mode)

    // Constructor. Creates a Microplane belonging to element e, with number
    // n, with coordinates a, with weight w.
{ }


Microplane :: ~Microplane()
// Destructor.
{ }



void
Microplane :: printOutputAt(FILE *File, TimeStep *tStep)
// Prints the strains and stresses on the data file.
{
    /*
     * int  i ;
     * // int n;
     *
     *
     * fprintf (File,"  GP %d :",number) ;
     * if (matStatus) matStatus -> printOutputAt (File,tStep) ;
     *
     * if ( numberOfGp != 0)  // layered material
     * {
     *  fprintf (File,"Layers report \n{\n");
     *  for ( i = 0; i< numberOfGp ; i++)
     *   {
     *    gaussPoints[i]->printOutputAt (File,tStep);
     *   }
     *  fprintf (File,"} end layers report\n");
     * }
     */
}
} // end namespace oofem
