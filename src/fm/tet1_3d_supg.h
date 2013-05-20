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

#ifndef tet1_3d_supg_h
#define tet1_3d_supg_h

#include "supgelement2.h"
#include "levelsetpcs.h"
#include "fei3dtetlin.h"

#define _IFT_Tet1_3D_SUPG_Name "tet1supg"

namespace oofem {

/**
 * Class representing 3d linear tetrahedral element
 * for solving incompressible fluid with SUPG solver
 */
class Tet1_3D_SUPG : public SUPGElement2, public LevelSetPCSElementInterface
{
protected:
    static FEI3dTetLin interpolation;

public:
    Tet1_3D_SUPG(int n, Domain *d);
    virtual ~Tet1_3D_SUPG();

    // definition
    virtual const char *giveClassName() const { return "Tet1_3D_SUPG"; }
    virtual const char *giveInputRecordName() const { return _IFT_Tet1_3D_SUPG_Name; }
    virtual classType giveClassID() const { return Tet1_3D_SUPGClass; }
    virtual Element_Geometry_Type giveGeometryType() const { return EGT_tetra_1; }
    virtual MaterialMode giveMaterialMode() { return _3dFlow; }
    virtual void giveElementDofIDMask(EquationID, IntArray & answer) const;
    virtual void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;
    virtual int computeNumberOfDofs(EquationID ut);
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual Interface *giveInterface(InterfaceType t);

    virtual double computeCriticalTimeStep(TimeStep *tStep);
    virtual double computeVolumeAround(GaussPoint *gp);

    virtual void updateStabilizationCoeffs(TimeStep *tStep);

    virtual double LS_PCS_computeF(LevelSetPCS *, TimeStep *);
    virtual void LS_PCS_computedN(FloatMatrix &answer);
    virtual double LS_PCS_computeVolume();
    virtual double LS_PCS_computeS(LevelSetPCS *, TimeStep *);
    virtual void LS_PCS_computeVOFFractions(FloatArray &answer, FloatArray &fi);

#ifdef __OOFEG
    void drawRawGeometry(oofegGraphicContext &);
#endif

protected:
    virtual void computeGaussPoints();
    virtual void computeNuMatrix(FloatMatrix &answer, GaussPoint *gp);
    virtual void computeUDotGradUMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *atTime);
    virtual void computeBMatrix(FloatMatrix &anwer, GaussPoint *gp);
    virtual void computeDivUMatrix(FloatMatrix &answer, GaussPoint *gp);
    virtual void computeNpMatrix(FloatMatrix &answer, GaussPoint *gp);
    virtual void computeGradPMatrix(FloatMatrix &answer, GaussPoint *gp);
    virtual void computeDivTauMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *atTime);
    virtual void computeGradUMatrix(FloatMatrix &answer, GaussPoint *gp, TimeStep *atTime);
    virtual int  giveNumberOfSpatialDimensions();

    virtual int giveTermIntergationRuleIndex(CharType termType);
};
} // end namespace oofem
#endif // tet1_3d_supg_h
