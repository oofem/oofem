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

#ifndef planestresselement_h
#define planestresselement_h

#include "../sm/Elements/nlstructuralelement.h"

///@name Input fields for PlaneStressElement

namespace oofem {
class TimeStep;
class Node;
class Material;
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;

/**
 * Base class for plane stress elements with geometrical nonlinearities.
 *
 * @author Jim Brouzoulis
 */
class PlaneStressElement : public NLStructuralElement
{

public:
    /**
     * Constructor. Creates element with given number, belonging to given domain.
     * @param n Element number.
     * @param d Domain to which new material will belong.
     */
    PlaneStressElement(int n, Domain * d);
    /// Destructor.
    virtual ~PlaneStressElement() { }

    virtual int computeNumberOfDofs();
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;
    virtual double computeVolumeAround(GaussPoint *gp);
    
    virtual MaterialMode giveMaterialMode() { return _PlaneStress; }
    
    virtual IRResultType initializeFrom(InputRecord *ir);
    
    
protected:
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) ;
    virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeGaussPoints();
    
    
    
    // Edge support
    virtual void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp);
    virtual void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp);
    
};

} // end namespace oofem
#endif // planestresselement_h
