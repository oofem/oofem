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

#ifndef qwedge_h
#define qwedge_h

#include "../sm/Elements/nlstructuralelement.h"
#include "../sm/ErrorEstimators/huertaerrorestimator.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "eleminterpmapperinterface.h"
#include "sprnodalrecoverymodel.h"

#define _IFT_QWedge_Name "qwedge"
#define _IFT_QWedge_materialCoordinateSystem "matcs" ///< [optional] Experimental support for material directions

namespace oofem {
class FEI3dWedgeQuad;

/**
 * This class implements an Quadratic 3d  15 - node structural finite element.
 * Each node has 3 degrees of freedom.
 *
 * @author Mikael Öhman (among others)
 */
class QWedge : public NLStructuralElement, public SPRNodalRecoveryModelInterface, public ZZNodalRecoveryModelInterface, public NodalAveragingRecoveryModelInterface

{
protected:
    static FEI3dWedgeQuad interpolation;

    bool matRotation;

public:
    QWedge(int, Domain *);
    virtual ~QWedge() { }

    virtual FEInterpolation *giveInterpolation() const;

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;
    virtual double computeVolumeAround(GaussPoint *);
    virtual int giveNumberOfIPForMassMtrxIntegration() { return 9; }

    virtual Interface *giveInterface(InterfaceType);
    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_SurfaceLoadSupport ) ? 1 : 0 ); }

    /**
     * Function determines the 3 basis vectors for the local coordinate system that should be used by materials.
     * @param x First base vector of the local c.s.
     * @param y Second base vector
     * @param z Third base vector
     * @param lcoords Local coordinates to evaluate the local c.s. for.
     */
    void giveMaterialOrientationAt(FloatArray &x, FloatArray &y, FloatArray &z, const FloatArray &lcoords);

    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    virtual int SPRNodalRecoveryMI_giveNumberOfIP();
    virtual void SPRNodalRecoveryMI_computeIPGlobalCoordinates(FloatArray &coords, GaussPoint *gp);
    virtual SPRPatchType SPRNodalRecoveryMI_givePatchType();

    virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep);

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_QWedge_Name; }
    virtual const char *giveClassName() const { return "QWedge"; }
    virtual int computeNumberOfDofs() { return 45; }
    virtual MaterialMode giveMaterialMode();

protected:
    virtual void computeGaussPoints();

    virtual void computeBmatrixAt(GaussPoint *, FloatMatrix &, int = 1, int = ALL_STRAINS);
    //virtual void computeBFmatrixAt(GaussPoint *, FloatMatrix &);
    virtual void computeBHmatrixAt(GaussPoint *, FloatMatrix &);
};
} // end namespace oofem
#endif
