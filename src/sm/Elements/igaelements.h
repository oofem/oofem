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

#ifndef igaelements_h
#define igaelements_h

#include "iga/iga.h"
#include "iga/feibspline.h"
#include "iga/feinurbs.h"
#include "iga/feitspline.h"
#include "sm/Elements/PlaneStress/planestresselementevaluator.h"
#include "sm/Elements/3D/space3delementevaluator.h"
#include "floatarray.h"
#include "floatmatrix.h"
#include "matresponsemode.h"
#include "mathfem.h"

#define _IFT_BsplinePlaneStressElement_Name "bsplineplanestresselement"
#define _IFT_NURBSPlaneStressElement_Name "nurbsplanestresselement"
#define _IFT_TSplinePlaneStressElement_Name "tsplineplanestresselement"
#define _IFT_NURBSSpace3dElement_Name "nurbs3delement"

namespace oofem {
class BsplinePlaneStressElement : public IGAElement, public PlaneStressStructuralElementEvaluator
{
protected:
    BSplineInterpolation interpolation;

public:
    BsplinePlaneStressElement(int n, Domain * aDomain);

    void initializeFrom(InputRecord &ir, int priority) override;
    int checkConsistency() override;

    void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep) override {
        PlaneStressStructuralElementEvaluator :: giveCharacteristicMatrix(answer, mtrx, tStep);
    }
    void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *t) override {
        PlaneStressStructuralElementEvaluator :: giveCharacteristicVector(answer, type, mode, t);
    }

    FEInterpolation *giveInterpolation() const override { return const_cast< BSplineInterpolation * >(& this->interpolation); }
    Element *giveElement() override { return this; }
    void giveDofManDofIDMask(int inode, IntArray &answer) const override {
        PlaneStressStructuralElementEvaluator :: giveDofManDofIDMask(inode, answer);
    }
    int computeNumberOfDofs() override { return numberOfDofMans * 2; }
    void updateInternalState(TimeStep *tStep) override { PlaneStressStructuralElementEvaluator :: updateInternalState(tStep); }
    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_BsplinePlaneStressElement_Name; }
    const char *giveClassName() const override { return "BsplinePlaneStressElement"; }
    Element_Geometry_Type giveGeometryType() const override {return EGT_unknown;}


#ifdef __OOFEG
    // Graphics output
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType ut) override {
        drawIGAPatchDeformedGeometry(this, this, gc, tStep, ut);
    }
#endif

protected:
    int giveNsd(const Element_Geometry_Type) override { return 2; }
};


class NURBSPlaneStressElement : public IGAElement, public PlaneStressStructuralElementEvaluator
{
protected:
    NURBSInterpolation interpolation;

public:
    NURBSPlaneStressElement(int n, Domain * aDomain);

    void initializeFrom(InputRecord &ir, int priority) override;
    int checkConsistency() override;

    void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep) override {
        PlaneStressStructuralElementEvaluator :: giveCharacteristicMatrix(answer, mtrx, tStep);
    }
    void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *t) override {
        PlaneStressStructuralElementEvaluator :: giveCharacteristicVector(answer, type, mode, t);
    }

    FEInterpolation *giveInterpolation() const override { return const_cast< NURBSInterpolation * >(& this->interpolation); }
    Element *giveElement() override { return this; }
    void giveDofManDofIDMask(int inode, IntArray &answer) const override {
        PlaneStressStructuralElementEvaluator :: giveDofManDofIDMask(inode, answer);
    }
    int computeNumberOfDofs() override { return numberOfDofMans * 2; }
    void updateInternalState(TimeStep *tStep) override { PlaneStressStructuralElementEvaluator :: updateInternalState(tStep); }
    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_NURBSPlaneStressElement_Name; }
    const char *giveClassName() const override { return "NURBSPlaneStressElement"; }
    Element_Geometry_Type giveGeometryType() const override {return EGT_unknown;}
#ifdef __OOFEG
    //
    // Graphics output
    //
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType ut) override {
        drawIGAPatchDeformedGeometry(this, this, gc, tStep, ut);
    }

#endif

protected:
    int giveNsd(const Element_Geometry_Type) override { return 2; }
};



class TSplinePlaneStressElement : public IGATSplineElement, public PlaneStressStructuralElementEvaluator
{
protected:
    TSplineInterpolation interpolation;

public:
    TSplinePlaneStressElement(int n, Domain * aDomain);

    void initializeFrom(InputRecord &ir, int priority) override {
        IGATSplineElement :: initializeFrom(ir, priority);
        //PlaneStressStructuralElementEvaluator::initializeFrom(ir);
    }

    void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep) override {
        PlaneStressStructuralElementEvaluator :: giveCharacteristicMatrix(answer, mtrx, tStep);
    }
    void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *t) override {
        PlaneStressStructuralElementEvaluator :: giveCharacteristicVector(answer, type, mode, t);
    }

    FEInterpolation *giveInterpolation() const override { return const_cast< TSplineInterpolation * >(& this->interpolation); }
    Element *giveElement() override { return this; }
    void giveDofManDofIDMask(int inode, IntArray &answer) const override {
        PlaneStressStructuralElementEvaluator :: giveDofManDofIDMask(inode, answer);
    }
    int computeNumberOfDofs() override { return numberOfDofMans * 2; }
    void updateInternalState(TimeStep *tStep) override { PlaneStressStructuralElementEvaluator :: updateInternalState(tStep); }
    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_TSplinePlaneStressElement_Name; }
    const char *giveClassName() const override { return "TSplinePlaneStressElement"; }
    Element_Geometry_Type giveGeometryType() const override {return EGT_unknown;}

#ifdef __OOFEG
    // Graphics output
    void  drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
#endif

protected:
    int giveNsd(const Element_Geometry_Type) override { return 2; }
};


class NURBSSpace3dElement : public IGAElement, public Space3dStructuralElementEvaluator
{
protected:
    NURBSInterpolation interpolation;

public:
    NURBSSpace3dElement(int n, Domain * aDomain);

    void initializeFrom(InputRecord &ir, int priority) override;
    int checkConsistency() override;

    void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep) override {
        Space3dStructuralElementEvaluator :: giveCharacteristicMatrix(answer, mtrx, tStep);
    }
    void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *t) override {
        Space3dStructuralElementEvaluator :: giveCharacteristicVector(answer, type, mode, t);
    }

    FEInterpolation *giveInterpolation() const override { return const_cast< NURBSInterpolation * >(& this->interpolation); }
    Element *giveElement() override { return this; }

    void giveDofManDofIDMask(int inode, IntArray &answer) const override {
        Space3dStructuralElementEvaluator :: giveDofManDofIDMask(inode, answer);
    }
    int computeNumberOfDofs() override { return numberOfDofMans * 3; }
    void updateInternalState(TimeStep *tStep) override { Space3dStructuralElementEvaluator :: updateInternalState(tStep); }
    // definition & identification
    const char *giveInputRecordName() const override { return _IFT_NURBSSpace3dElement_Name; }
    const char *giveClassName() const override { return "NURBSSpace3dElement"; }
    Element_Geometry_Type giveGeometryType() const override {return EGT_unknown;}

#ifdef __OOFEG
    // Graphics output
    void drawScalar(oofegGraphicContext &gc, TimeStep *tStep) override;
    void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType ut) override {
        drawIGAPatchDeformedGeometry(this, this, gc, tStep, ut);
    }
#endif

protected:
    int giveNsd(const Element_Geometry_Type) override { return 3; }
};
} // end namespace oofem
#endif //igaelements_h
