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
 *               Copyright (C) 1993 - 2014   Borek Patzak
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

#ifndef baseelement_h
#define baseelement_h


#define _IFT_BaseElement_partitions "partitions"
#define _IFT_BaseElement_remote "remote"

namespace oofem {
class ElementEvaluator;
class ElementGeometry;

/**
 * Abstract temporary class providing proper functioning of the old Element class within the new design 
 * of the code. The design is now based on decoupling of orginal Element class into
 * ElementGeometry and ElementEvaluator classes which provides more flexibility for
 * multiphysics problems or problems with several different degrees of freedom, e.g. 
 * gradient damage models
 */
class OOFEM_EXPORT BaseElement 
{
public:
    /**
     * Constructor. Creates an abstract base element class.     
     */
	BaseElement(){;}
    /// Virtual destructor.
	virtual ~BaseElement(){;}
	 /// @return ElementGeometry which receiver belongs to.
	virtual ElementGeometry* giveElementGeometry() = 0;
	 /// @return ElementEvaluator which receiver belongs to.
	virtual ElementEvaluator* giveElementEvaluator() = 0;   
};
} // end namespace oofem
#endif //baseelement_h
