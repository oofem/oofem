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

#ifndef tr_warp_h
#define tr_warp_h

#include "structuralelement.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"
#include "fei2dtrlin.h"

#define _IFT_Tr_Warp_Name "trwarp"


namespace oofem {
/**
 * Triangle (2d) element with linear approximation for free warping analysis.
 */
class Tr_Warp : public StructuralElement, public SpatialLocalizerInterface, public ZZNodalRecoveryModelInterface
{
protected:
    static FEI2dTrLin interp;


public:
    Tr_Warp(int n, Domain *d);
    virtual ~Tr_Warp();
    virtual void computeFirstMomentOfArea(FloatArray &answer);
    virtual double computeVolumeAround(GaussPoint *gp);
    void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
    virtual void computeLocalForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode);
    virtual void computeEdgeLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode);



    // definition
    virtual const char *giveInputRecordName() const { return _IFT_Tr_Warp_Name; }
    virtual const char *giveClassName() const { return "Tr_WarpElement"; }

    virtual int computeNumberOfDofs() { return 4; }
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual MaterialMode giveMaterialMode() { return _Warping; }
    virtual double giveThicknessAt(const FloatArray &gcoords);
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer);
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int li, int ui);
    void giveDofManDofIDMask(int inode, IntArray &answer) const;

    virtual Interface *giveInterface(InterfaceType t);

    virtual void giveInternalDofManDofIDMask(int inode, IntArray &answer) const;
    virtual DofManager *giveInternalDofManager(int i) const;



    virtual int giveNumberOfInternalDofManagers() const { return 1; }



    virtual Element *ZZNodalRecoveryMI_giveElement() { return this; }
    void ZZNodalRecoveryMI_computeNNMatrix(FloatArray &answer, InternalStateType type);
    bool ZZNodalRecoveryMI_computeNValProduct(FloatMatrix &answer, InternalStateType type,
                                              TimeStep *tStep);
    virtual Element *SpatialLocalizerI_giveElement() { return this; }
    virtual int SpatialLocalizerI_containsPoint(const FloatArray &coords);
    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords);

    virtual FEInterpolation *giveInterpolation() const { return & this->interp; }
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_triangle_1; }
#ifdef __OOFEG
    // Graphics output
    //void drawYourself(oofegGraphicContext&);
    //virtual void drawRawGeometry(oofegGraphicContext&) {}
    //virtual void drawDeformedGeometry(oofegGraphicContext&, UnknownType) {}
#endif

protected:
    virtual void computeGaussPoints();
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    void transformCoordinates(FloatArray &answer, FloatArray &c, const int CGnumber);
    virtual void  postInitialize();
};
} // end namespace oofem
#endif // tr_warp_h
