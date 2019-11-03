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

#ifndef htselement_h
#define htselement_h

#include "sm/Elements/structuralelement.h"
#include "gaussintegrationrule.h"

#define _IFT_HTSelement_Name "htselement"

namespace oofem {
  /**
 * Implements a Hybrid-Trefftz element
 * See http://en.wikipedia.org/wiki/Trefftz_method for description.
 * @author Jan Novak (among others?) 
 */
class HTSelement : public StructuralElement
{
protected:
    int numberOfEdges;
    //debug
    double lambda, mu;
    double cgX, cgY;
    int numberOfStressDofs;
    int numberOfDofs;

public:
    HTSelement(int n, Domain * d);
    virtual ~HTSelement() { }

    void initializeFrom(InputRecord &ir) override;

    const char *giveInputRecordName() const override { return _IFT_HTSelement_Name; }
    const char *giveClassName() const override { return "HTSelement"; }

protected:
    void computeBmatrixAt(GaussPoint *gp, FloatMatrix &answer, int, int) override { }
    void computeNmatrixAt(const FloatArray &iLocCoord, FloatMatrix &answer) override { }
    void computeStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep) override;
    double computeVolumeAroundSide(GaussPoint *gp, int elemSideNumber);
    Node *giveSideNode(int elementSideNumber, int nodeNumber);
    double giveSideLength(int sideNumber);
    int computeNumberOfDofs() override { return 4 * numberOfEdges; }
    void computeGaussPoints() override;
    void giveDofManDofIDMask(int inode, IntArray &) const override;
    virtual StructuralElement *giveStructuralElement() { return this; }
    //jak se pocita deformace???
    void computeStrainVector(FloatArray &answer, GaussPoint *gp, TimeStep *tStep) override { answer.resize(numberOfStressDofs); }
    //dodelat vypocet napeti!!!
    void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep) override { answer.resize(numberOfStressDofs); }
    void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) override {
        OOFEM_ERROR("Function not defined for this element and should never be called. This is a bug.");
    }
    //dodelat internal forces, budou potreba pro nelinearni vypocet
    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord) override { answer.resize(numberOfDofs); }

    void computePuVectorAt(FloatArray &answer, FloatMatrix N, FloatArray u, GaussPoint *gp, int sideNumber);
    void computePsVectorAt(FloatArray &answer, FloatArray t, GaussPoint *gp);
    void computePrescribedDisplacementLoadVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode);
    int testElementExtension(ElementExtension ext) override { return ( ext == Element_EdgeLoadSupport ); }

    void computeFMatrixAt(FloatMatrix &answer, FloatMatrix N, GaussPoint *gp, int sideNumber);
    void computeAMatrixAt(FloatMatrix &answer, FloatMatrix N, GaussPoint *gp, int sideNumber);
    void computeUvMatrixAt(FloatMatrix &answer, GaussPoint *gp, int sideNubmer);
    void computeSvMatrixAt(FloatMatrix &answer, GaussPoint *gp, int sideNumber);
    void computeUgammaMatrixAt(FloatMatrix &answer, GaussPoint *gp);
    void computeOutwardNormalMatrix(FloatMatrix &answer, int sideNumber);

    void computeCenterOfGravity();
    int giveNumberOfNodes() const override { return numberOfEdges; }
    //uv functions
    void uv1(FloatArray &answer, double x, double y);
    void uv2(FloatArray &answer, double x, double y);
    void uv3(FloatArray &answer, double x, double y);
    void uv4(FloatArray &answer, double x, double y);
    void uv5(FloatArray &answer, double x, double y);
    void uv6(FloatArray &answer, double x, double y);
    void uv7(FloatArray &answer, double x, double y);
    void uv8(FloatArray &answer, double x, double y);
    void uv9(FloatArray &answer, double x, double y);
    void uv10(FloatArray &answer, double x, double y);
    void uv11(FloatArray &answer, double x, double y);
    void uv12(FloatArray &answer, double x, double y);
    void uv25_4(FloatArray &answer, double x, double y);

    //sv functions
    void sv1(FloatArray &answer, double x, double y);
    void sv2(FloatArray &answer, double x, double y);
    void sv3(FloatArray &answer, double x, double y);
    void sv4(FloatArray &answer, double x, double y);
    void sv5(FloatArray &answer, double x, double y);
    void sv6(FloatArray &answer, double x, double y);
    void sv7(FloatArray &answer, double x, double y);
    void sv8(FloatArray &answer, double x, double y);
    void sv9(FloatArray &answer, double x, double y);
    void sv10(FloatArray &answer, double x, double y);
    void sv11(FloatArray &answer, double x, double y);
    void sv12(FloatArray &answer, double x, double y);
    void sv25_4(FloatArray &answer, double x, double y);

    //u_gamma functions
    double u_gammaConst(GaussPoint *gp);
    double u_gammaLin(GaussPoint *gp);
};
} // end namespace oofem
#endif // htselement_h
