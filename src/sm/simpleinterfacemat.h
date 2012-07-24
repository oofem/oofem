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

//   **********************************************************************
//   *** CLASS SIMPLE INTERFACE MATERIAL FOR INTERFACE ELEMENTS   *********
//   **********************************************************************

#ifndef simpleinterfacemat_h
#define simpleinterfacemat_h

#include "material.h"
#include "linearelasticmaterial.h"
#include "structuralmaterial.h"
#include "structuralms.h"

namespace oofem {

/**
 * This class implements associated Material Status to SimpleInterfaceMaterial.
 */
class SimpleInterfaceMaterialStatus : public StructuralMaterialStatus
{
protected:
    FloatArray shearStressShift, tempShearStressShift;

public:
    /// Constructor
    SimpleInterfaceMaterialStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor
    virtual ~SimpleInterfaceMaterialStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    // definition
    virtual const char *giveClassName() const { return "SimpleInterfaceMaterialStatus"; }
    virtual classType giveClassID() const { return SimpleInterfaceMaterialStatusClass; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    FloatArray giveShearStressShift();
    FloatArray giveTempShearStressShift();
    void setTempShearStressShift(FloatArray newShearStressShift) { tempShearStressShift = newShearStressShift; };

    virtual contextIOResultType saveContext(DataStream *stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream *stream, ContextMode mode, void *obj = NULL);
};


/**
 * Base class representing general isotropic damage model.
 * It is based on isotropic damage concept, assuming that damage evolution law
 * is postulated in explicit form, relatin damage parameter (omega) to scalar measure
 * of the largest strain level ever reached in material (kappa).
 */
class SimpleInterfaceMaterial : public StructuralMaterial
{
protected:
    double kn;
    double stiffCoeff;
    double frictCoeff;
    /// Normal distance which needs to be closed when interface element should act in compression (distance is 0 by default).
    double normalClearance;

public:
    /// Constructor
    SimpleInterfaceMaterial(int n, Domain *d);
    /// Destructor
    virtual ~SimpleInterfaceMaterial();

    virtual int hasNonLinearBehaviour() { return 1; }
    virtual int hasMaterialModeCapability(MaterialMode mode);
    virtual const char *giveClassName() const { return "SimpleInterfaceMaterial"; }
    virtual classType giveClassID() const { return SimpleInterfaceMaterialClass; }

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseForm from,
                                               MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *atTime);

    virtual void giveRealStressVector(FloatArray & answer,  MatResponseForm, GaussPoint *,
                              const FloatArray &, TimeStep *);

    virtual void  giveCharacteristicMatrix(FloatMatrix &answer,
                                           MatResponseForm form,
                                           MatResponseMode mode,
                                           GaussPoint *gp,
                                           TimeStep *atTime);

    virtual int giveStressStrainComponentIndOf(MatResponseForm, MaterialMode mmode, int);
    virtual void giveStressStrainMask(IntArray & answer, MatResponseForm, MaterialMode mmode) const;
    virtual int giveSizeOfReducedStressStrainVector(MaterialMode);
    void giveReducedCharacteristicVector(FloatArray &answer, GaussPoint *gp,
                                         const FloatArray &charVector3d);
    void giveFullCharacteristicVector(FloatArray &answer,  GaussPoint *gp,
                                      const FloatArray &strainVector);

#ifdef __OOFEG
#endif

    virtual int giveIPValue(FloatArray &answer, GaussPoint *aGaussPoint, InternalStateType type, TimeStep *atTime);
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type, MaterialMode mmode);
    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    virtual int giveIPValueSize(InternalStateType type, GaussPoint *aGaussPoint);
    virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual int giveInputRecordString(std :: string &str, bool keyword = true);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new SimpleInterfaceMaterialStatus(1, domain, gp); }
};
} // end namespace oofem
#endif // simpleinterfacemat_h
