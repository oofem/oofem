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

#ifndef structural2delement_h
#define structural2delement_h

#include "sm/Elements/nlstructuralelement.h"
#include "feinterpol2d.h"

namespace oofem {
class ParamKey;
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
    FEICellGeometry *cellGeometryWrapper;

    bool matRotation;

    static ParamKey IPK_Structural2DElement_materialCoordinateSystem; ///< [optional] Material coordinate system (local) for the element.
public:
    /**
     * Constructor. Creates element with given number, belonging to given domain.
     * @param n Element number.
     * @param d Domain to which new material will belong.
     */
    Structural2DElement(int n, Domain *d);
    /// Destructor.
    virtual ~Structural2DElement();
    void postInitialize() override;
    int giveNumberOfNodes() const override;
    /**
     * Returns the Cell Geometry Wrapper. Default inplementation creates FEIElementGeometryWrapper.
     */
    virtual FEICellGeometry *giveCellGeometryWrapper();

    int computeNumberOfDofs() override;
    void giveDofManDofIDMask(int inode, IntArray &answer) const override;
    double computeVolumeAround(GaussPoint *gp) override;

    void initializeFrom(InputRecord &ir, int priority) override;

    double giveCharacteristicLength(const FloatArray &normalToCrackPlane) override;

protected:
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) override = 0;
    void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer) override = 0;
    void computeGaussPoints() override;

    void giveMaterialOrientationAt(FloatArray &x, FloatArray &y, const FloatArray &lcoords);

    // Edge support
    void giveEdgeDofMapping(IntArray &answer, int iEdge) const override;
    double computeEdgeVolumeAround(GaussPoint *gp, int iEdge) override;
    int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp) override;
    int testElementExtension(ElementExtension ext) override { return ( ( ext == Element_EdgeLoadSupport ) ? 1 : 0 ); }
    
};


class PlaneStressElement : public Structural2DElement
{
public:
    PlaneStressElement(int n, Domain *d);
    virtual ~PlaneStressElement() { }
    MaterialMode giveMaterialMode() override { return _PlaneStress; }
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrix_dPdF_At(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    
protected:
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) override;
    void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer) override;
};


class PlaneStrainElement : public Structural2DElement
{
public:
    PlaneStrainElement(int n, Domain *d);
    virtual ~PlaneStrainElement() { }
    MaterialMode giveMaterialMode() override { return _PlaneStrain; }
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrix_dPdF_At(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;

protected:
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) override;
    void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer) override;
};


class AxisymElement : public Structural2DElement
{
public:
    AxisymElement(int n, Domain *d);
    virtual ~AxisymElement() { }
    MaterialMode giveMaterialMode() override { return _3dMat; }
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;
    void computeConstitutiveMatrix_dPdF_At(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override;

    double giveCharacteristicLength(const FloatArray &crackToNormalPlane) override;
    double computeVolumeAround(GaussPoint *gp) override;

protected:
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int lowerIndx = 1, int upperIndx = ALL_STRAINS) override;
    void computeBHmatrixAt(GaussPoint *gp, FloatMatrix &answer) override;
    void computeGaussPoints() override;
    double computeEdgeVolumeAround(GaussPoint *gp, int iEdge) override;
};
} // end namespace oofem
#endif // structural2delement_h
