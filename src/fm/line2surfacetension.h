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
 *               Copyright (C) 1993 - 2010   Borek Patzak
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

#ifndef line2surfacetension_h
#define line2surfacetension_h

#include "linesurfacetension.h"

namespace oofem {
/**
 * 3 node line elements for surface tension.
 * @see SurfaceTension2D
 * @author Mikael Öhman
 */
class Line2SurfaceTension : public LineSurfaceTension
{
    //protected:
    //static FEI2dLineQuad interpolation; // TODO: Implement this

public:
    /**
     * Constructor. Creates an element with number n belonging to domain aDomain.
     * @param n Element's number
     * @param aDomain Pointer to the domain to which element belongs.
     */
    Line2SurfaceTension(int, Domain *);
    /// Destructor.
    ~Line2SurfaceTension();

    /**
     * Gives the load from surface tension.
     * If element has no capability to compute requested type of characteristic vector
     * error function is invoked.
     * @param answer requested characteristic vector
     * @param type    only supports loadvector type
     * @param mode    ignored.
     * @param tStep   ignored (independent).
     */
    //void giveCharacteristicVector(FloatArray & answer, CharType, ValueModeType, TimeStep *);

    /**
     * Gives the tangent from surface tension.
     */
    void computeTangent(FloatMatrix &answer, TimeStep *tStep);

    void computeLoadVector(FloatArray &answer, ValueModeType mode, TimeStep *tStep);

    /**
     * Returns the type of geometry.
     * @return Quadratic line element type.
     */
    Element_Geometry_Type giveGeometryType() const { return EGT_line_2; }

    /**
     * Always 6 DOFs (V_u, V_v)*3.
     * @param[in] ut ignored.
     * @return 6 degrees of freedom.
     */
    virtual int computeNumberOfDofs(EquationID ut) { return 6; }
};
}

#endif // line2surfacetension_h
