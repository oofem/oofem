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

#ifndef tria1platesubsoil_H
#define tria1platesubsoil_H

#include "../sm/Elements/structuralelement.h"
#include "zznodalrecoverymodel.h"
#include "sprnodalrecoverymodel.h"

#define _IFT_Tria1PlateSubSoil_Name "tria1platesubsoil"

namespace oofem {
class FEI2dTrLin;

/**
 * This class implements an triangular four-node plate subsoil element with linear interpolation in xy plane.
 * Each node has 1 degree of freedom (out-of-plane displacement).
 *
 * Loading types supported;
 *
 * Reference:
 * Bittnar, Sejnoha: Numerical Methods in Structural Mechanics, Thomas Telford, Jan 1, 1996, ISBN-13: 978-0784401705
 * @author Borek Patzak
 */
class Tria1PlateSubSoil : public StructuralElement,
public ZZNodalRecoveryModelInterface,
public SPRNodalRecoveryModelInterface
{
protected:
    static FEI2dTrLin interp_lin;

public:
    Tria1PlateSubSoil(int n, Domain * d);
    virtual ~Tria1PlateSubSoil() { }

    virtual FEInterpolation *giveInterpolation() const;
    virtual FEInterpolation *giveInterpolation(DofIDItem id) const;

    virtual MaterialMode giveMaterialMode()  { return _2dPlateSubSoil; }
    virtual int testElementExtension(ElementExtension ext) { return ( ( ext == Element_SurfaceLoadSupport ) ? 1 : 0 ); }

    // definition & identification
    virtual const char *giveInputRecordName() const { return _IFT_Tria1PlateSubSoil_Name; }
    virtual const char *giveClassName() const { return "Tria1PlateSubSoil"; }
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual int computeNumberOfDofs() { return 3; }
    virtual void giveDofManDofIDMask(int inode, IntArray &) const;

    virtual void computeMidPlaneNormal(FloatArray &answer, const GaussPoint *gp);

    virtual double giveCharacteristicLength(const FloatArray &normalToCrackPlane);
    virtual double computeVolumeAround(GaussPoint *gp);

    virtual void computeLumpedMassMatrix(FloatMatrix &answer, TimeStep *tStep);
    virtual void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep)
    { computeLumpedMassMatrix(answer, tStep); }

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    virtual Interface *giveInterface(InterfaceType it);

protected:
    virtual void computeGaussPoints();
    virtual void computeBodyLoadVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode);
    virtual void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int = 1, int = ALL_STRAINS);

    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap);
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap);
    virtual int SPRNodalRecoveryMI_giveNumberOfIP() { return this->numberOfGaussPoints; }
    virtual SPRPatchType SPRNodalRecoveryMI_givePatchType() { return SPRPatchType_2dxy; }
    /**
     * @name Surface load support
     */
    //@{
    virtual void computeSurfaceNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *gp);
    virtual void giveSurfaceDofMapping(IntArray &answer, int iSurf) const;
    virtual IntegrationRule *GetSurfaceIntegrationRule(int iSurf);
    virtual double computeSurfaceVolumeAround(GaussPoint *gp, int iSurf);
    virtual void computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iSurf);
    virtual int computeLoadLSToLRotationMatrix(FloatMatrix &answer, int iSurf, GaussPoint *gp);
    //@}
};
} // end namespace oofem
#endif // tria1platesubsoil_H
