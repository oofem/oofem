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

#ifndef brick1_ht_h
#define brick1_ht_h

#include "transportelement.h"
#include "spatiallocalizer.h"
#include "zznodalrecoverymodel.h"
#include "fei3dhexalin.h"

namespace oofem {
/**
 * Brick (3d) elements with linear approximation for heat and mass transfer.
 */
class Brick1_ht : public TransportElement, public SpatialLocalizerInterface, public ZZNodalRecoveryModelInterface
{
protected:
    static FEI3dHexaLin interpolation;
    int numberOfGaussPoints;

public:
    Brick1_ht(int n, Domain *d);
    virtual ~Brick1_ht();

    virtual void computeInternalSourceRhsVectorAt(FloatArray &answer, TimeStep *, ValueModeType mode);
    virtual double computeVolumeAround(GaussPoint *gp);

    // definition
    virtual const char *giveClassName() const { return "Brick1_ht"; }
    virtual classType giveClassID() const { return Brick1_htClass; }

    virtual int computeNumberOfDofs(EquationID ut) { return ( emode == HeatTransferEM ) ? 8 : 16; }
    virtual void giveDofManDofIDMask(int inode, EquationID, IntArray &) const;
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_hexa_1; }

    virtual Interface *giveInterface(InterfaceType t);
    virtual int testElementExtension(ElementExtension ext)
    { return ( ext == Element_EdgeLoadSupport ) || ( ext == Element_SurfaceLoadSupport ); }

    virtual int ZZNodalRecoveryMI_giveDofManRecordSize(InternalStateType type);
    virtual Element *ZZNodalRecoveryMI_giveElement() { return this; }
    virtual void ZZNodalRecoveryMI_ComputeEstimatedInterpolationMtrx(FloatMatrix &answer, GaussPoint *aGaussPoint, InternalStateType type);

    virtual Element *SpatialLocalizerI_giveElement() { return this; }
    virtual int SpatialLocalizerI_containsPoint(const FloatArray &coords);
    virtual double SpatialLocalizerI_giveDistanceFromParametricCenter(const FloatArray &coords);

#ifdef __OOFEG
    // Graphics output
    void drawRawGeometry(oofegGraphicContext &);
    virtual void drawScalar(oofegGraphicContext &context);
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

    virtual IntegrationRule *GetSurfaceIntegrationRule(int approxOrder);
    virtual void computeSurfaceNMatrixAt(FloatMatrix &n, GaussPoint *gp);
    virtual double computeSurfaceVolumeAround(GaussPoint *gp, int iEdge);
    virtual void giveSurfaceDofMapping(IntArray &mask, int iEdge);
    virtual void computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iSurf);

    virtual int giveApproxOrder(int unknownIndx) { return 1; }
};

class Brick1_hmt : public Brick1_ht
{
public:
    Brick1_hmt(int n, Domain *d);

    virtual const char *giveClassName() const { return "Brick1_hmt"; }
    virtual classType giveClassID() const { return Brick1_hmtClass; }
};

} // end namespace oofem
#endif // brick1_ht_h
