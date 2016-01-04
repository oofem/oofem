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

#ifndef structural2delement_h
#define structural2delement_h

#include "Elements/nlstructuralelement.h"
#include "feinterpol2d.h"

#define _IFT_Structural2DElement_materialCoordinateSystem "matcs" ///< [optional] Support for material directions based on element orientation.

namespace oofem {

/**
 * Base class for planar 2D elements.
 *
 * @author Jim Brouzoulis
 */
class Structural2DElement : public NLStructuralElement
{

protected:
    /** 
     * To facilitate the transformation of 2d elements into 3d, the complexity of transformation from 3d to
     *  local 2d system can be efficiently hidden in custom FEICellGeometry wrapper, that performs the transformation
     * into local system. This way, the existing 2d intrpolation classes can be used.
     * The element maintain its FEICellGeometry, which is accesible through the giveCellGeometryWrapper. 
     * Generalization to 3d then would require only substitution of the geometry warpper and definition of 
     * element transformation matrix.
     */
    FEICellGeometry* cellGeometryWrapper;

    bool matRotation;

public:
    /**
     * Constructor. Creates element with given number, belonging to given domain.
     * @param n Element number.
     * @param d Domain to which new material will belong.
     */
    Structural2DElement(int n, Domain * d);
    /// Destructor.
    virtual ~Structural2DElement();
    virtual void postInitialize();
    virtual int giveNumberOfNodes() const;
    /**
     * Returns the Cell Geometry Wrapper. Default inplementation creates FEIElementGeometryWrapper.
     */
    virtual FEICellGeometry* giveCellGeometryWrapper();
    
    virtual int computeNumberOfDofs();
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;
    virtual double computeVolumeAround(GaussPoint *gp);
    
    virtual IRResultType initializeFrom(InputRecord *ir);
    
    virtual double giveCharacteristicLength(const FloatArray &normalToCrackPlane);
    
protected:
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) = 0 ;
    virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer) = 0;
    virtual void computeGaussPoints();

    void giveMaterialOrientationAt( FloatArray &x, FloatArray &y, const FloatArray &lcoords);
    
    // Edge support
    virtual void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp);
    virtual void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp);
    virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 ); }
};



class PlaneStressElement : public Structural2DElement
{
public:
    PlaneStressElement(int n, Domain * d);
    virtual ~PlaneStressElement() { }
    virtual MaterialMode giveMaterialMode() { return _PlaneStress; }
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

protected:
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) ;
    virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer);
};


class PlaneStrainElement : public Structural2DElement
{
public:
    PlaneStrainElement(int n, Domain * d);
    virtual ~PlaneStrainElement() { }
    virtual MaterialMode giveMaterialMode() { return _PlaneStrain; }
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

protected:
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) ;
    virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer);
};


class AxisymElement : public Structural2DElement
{
public:
    AxisymElement(int n, Domain * d);
    virtual ~AxisymElement() { }
    virtual MaterialMode giveMaterialMode() { return _3dMat; }
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

    virtual double giveCharacteristicLength(const FloatArray &crackToNormalPlane);
    virtual double computeVolumeAround(GaussPoint *gp);

protected:
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) ;
    virtual void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer);
    virtual void computeGaussPoints();
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
};


} // end namespace oofem
#endif // structural2delement_h
