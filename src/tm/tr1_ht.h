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
 *               Copyright (C) 1993 - 2012   Borek Patzak
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

#ifndef tr1_ht_h
#define tr1_ht_h

#include "transportelement.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"
#include "fei2dtrlin.h"

namespace oofem {
/**
 * Triangle (2d) element with linear approximation for heat transfer.
 * @todo Use the interpolation classes.
 */
class Tr1_ht : public TransportElement, public SpatialLocalizerInterface, public ZZNodalRecoveryModelInterface
{
protected:
    int numberOfGaussPoints;
    static FEI2dTrLin interp;

public:
    Tr1_ht(int n, Domain *d);
    virtual ~Tr1_ht();

    virtual void computeInternalSourceRhsVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode);
    virtual double computeVolumeAround(GaussPoint *gp);

    virtual int computeGlobalCoordinates(FloatArray &answer, const FloatArray &lcoords);
    virtual int computeLocalCoordinates(FloatArray &answer, const FloatArray &gcoords);
    // definition
    virtual const char *giveClassName() const { return "Tr1_htElement"; }
    virtual classType giveClassID() const { return Tr1_htClass; }

    virtual int computeNumberOfDofs(EquationID ut) { return ( emode == HeatTransferEM ) ? 3 : 6; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_triangle_1; }

    virtual Interface *giveInterface(InterfaceType t);

    virtual int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);
    virtual Element *ZZNodalRecoveryMI_giveElement() { return this; }
    virtual void ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatMatrix &answer, GaussPoint *aGaussPoint, InternalStateType type);

    virtual Element *SpatialLocalizerI_giveElement() { return this; }
    virtual int SpatialLocalizerI_containsPoint(const FloatArray &coords);
    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords);

    virtual FEInterpolation *giveInterpolation() { return &this->interp; }
    virtual FEInterpolation *giveInterpolation(DofIDItem id);

#ifdef __OOFEG
    // Graphics output
    //void drawYourself(oofegGraphicContext&);
    //virtual void drawRawGeometry(oofegGraphicContext&) {}
    //virtual void drawDeformedGeometry(oofegGraphicContext&, UnknownType) {}
#endif

protected:
    virtual void computeGaussPoints();

    virtual void computeGradientMatrixAt(FloatMatrix &answer, GaussPoint *gp);
    virtual void computeNmatrixAt(FloatMatrix &n, FloatArray *lcoords);
    virtual void computeNSubMatrixAt(FloatMatrix &n, FloatArray *lcoords);

    virtual void computeEgdeNMatrixAt(FloatMatrix &n, GaussPoint *gp);
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
    virtual void giveEdgeDofMapping(IntArray &mask, int iEdge);
    virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);

    virtual int giveApproxOrder(int unknownIndx) { return 1; }
};

/**
 * Class for heat and mass transfer.
 */
class Tr1_hmt: public Tr1_ht
{
public:
    Tr1_hmt(int n, Domain *d);
    virtual const char *giveClassName() const { return "Tr1_hmt"; }
    virtual classType giveClassID() const { return Tr1_hmtClass; }
};

} // end namespace oofem
#endif // tr1_ht_h
