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

#ifndef contactelement_h
#define contactelement_h

#include "domain.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "feinterpol.h"
#include "element.h"
#include "aabb.h"

namespace oofem {

/**
 * Abstract base class for all contact finite elements. Derived classes should be  base
 * classes for specific analysis type (for example base class for structural analysis,
 * thermal analysis or magnetostatics one). These derived classes then declare
 * analysis-specific part of interface and they provide default implementation
 * for these methods.
 * This abstract class declares (and possibly implements) general data and methods
 * common to all element types. General methods for obtaining characteristic vectors,
 * matrices and values are introduced and should be used instead of calling directly
 * specific member functions (these must be overloaded by derived analysis-specific
 * classes in order to invoke proper method according to type of component requested).
 */
  class OOFEM_EXPORT ContactElement : public Element
{
protected:
 public:
    /**
     * Constructor. Creates an element with number n belonging to domain aDomain.
     * @param n Element's number
     * @param aDomain Pointer to the domain to which element belongs.
     */
     ContactElement(int n, Domain * aDomain);
    /// Virtual destructor.
    virtual ~ContactElement();
    /**@name General element functions */
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) = 0;
    virtual FloatArray computeNormalVectorAt(const FloatArray &lCoords) = 0;
    void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep) override
      
    {answer.clear();}
    void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep) override {answer.clear();} 

    void updateYourself(TimeStep *tStep) override{;}
    void printOutputAt(FILE *file, TimeStep *tStep) override{;}

    
    virtual AABB computeAABB();
};



} // end namespace oofem
#endif //contactelement_h
