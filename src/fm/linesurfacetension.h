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

#ifndef linesurfacetension_h
#define linesurfacetension_h

#include "fmelement.h"
#include "structuralelement.h"
#include "domain.h"

namespace oofem {

/**
 * Implements the load and tangent for surfacen tension boundary potential.
 * This class is a base class for higher order and specifically introduces simplifications for linear elements, and components are exact.
 * The derivations are to long to write here, see documentation for further detail.
 *
 * This class could serve as a reference for the math if these parts were to ever be moved towards the boundary conditions.
 * I'm not convinced that it belongs to a boundary condition class.
 * The load is dependent on the solution, thus there is a "stiffness"/tangent matrix as well.
 * @author Mikael Öhman
 */
class LineSurfaceTension : public FMElement//, public StructuralElement <-- I'd like to be both?
{
protected:
	/// Material parameter for the surface energy. Defines the traction as ${\bf t}_s = 2\kappa\gamma_s{\bf n}$.
	double gamma_s;
	/// Flags denoting if a a node is on the boundary.
	bool bflag1, bflag2;
	/// True if used as FMElement
	bool fmtype;

public:
    /**
     * Constructor. Creates an element with number n belonging to domain aDomain.
     * @param n Element's number
     * @param aDomain Pointer to the domain to which element belongs.
     */
	LineSurfaceTension (int, Domain *);
    /// Destructor.
	~LineSurfaceTension ();

	/**
	 * Initializes the element.
	 * Reads
	 * - gamma_s (required)
	 */
	IRResultType initializeFrom(InputRecord *);

    /**
     * Only computes internal load from surface tension.
     * @see Element :: giveCharacteristicVector
     */
	virtual void giveCharacteristicVector(FloatArray & answer, CharType, ValueModeType, TimeStep *tStep);

    /**
     * The element has only a load, the characteristic matrix is always zero.
     * @param answer 4 by 4 matrix with zeros.
     */
	virtual void giveCharacteristicMatrix(FloatMatrix & answer, CharType, TimeStep *tStep);

	/**
	 * Computes the load vector.
	 */
	virtual void computeLoadVector(FloatArray & answer, ValueModeType, TimeStep *tStep);

	/**
	 * Computes tangent to load vector.
	 */
	virtual void computeTangent(FloatMatrix & answer, TimeStep *tStep);

	/**
	 * Returns the type of geometry.
	 * @return Linear line element type.
	 */
	virtual Element_Geometry_Type giveGeometryType() const { return EGT_line_1; }

	/**
	 * (V_u, V_v)x2
	 * @return 4 degrees of freedom.
	 */
	virtual int computeNumberOfDofs(EquationID ut) { return 4;}

	/**
	 * Returns the displacement or velocities, depending on the domain type.
	 */
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;

};

} // end namespace oofem

#endif /* SURFACETENSION2D_H_ */
