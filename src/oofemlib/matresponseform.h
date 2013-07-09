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

#ifndef matesponseform_h
#define matesponseform_h

namespace oofem {
/**
 * Type representing the form of returned characteristic value (for cross section and material models).
 * The response can be returned in so called full or reduced form.
 * Generally, the full form contain all components, even if they are generally always zero (based on MaterialMode of
 * given integration point). On the other hand, the reduced form contain only generally nonzero components.
 * For example the "full-like" strain vector contains six components. For integration point in plane stress mode
 * the "reduced-like" strain vector contains only 3 generally nonzero components. The contents of full and reduced forms
 * is defined by corresponding (material or cross section level) base classes.
 */
enum MatResponseForm {
    ReducedForm, ///< Only stiffness for necessary stresses are given.
    FullForm,    ///< all component of 3d stresses are available, even if they equal 0.
};
} // end namespace oofem
#endif // matesponseform_h

