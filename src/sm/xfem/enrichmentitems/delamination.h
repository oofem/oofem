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

#ifndef DELAMINATION_H_
#define DELAMINATION_H_

#include "xfem/listbasedei.h"

///@name Input fields for Delamination
//@{
#define _IFT_Delamination_Name "delamination"
#define _IFT_Delamination_xiCoord "delaminationxicoord"
#define _IFT_Delamination_interfacenum "interfacenum"
#define _IFT_Delamination_csnum "csnum"
#define _IFT_Delamination_CohesiveZoneMaterial "czmaterial"
#define _IFT_Delamination_initiationFactor "initiationfactor"
#define _IFT_Delamination_initiationRadius "initiationradius"
#define _IFT_Delamination_averageStresses "averageStresses"
//@}

namespace oofem {
class Material;

/**
 * Delamination.
 * */
class OOFEM_EXPORT Delamination : public ListBasedEI
{
protected:
    Material *mat;  // Material for cohesive zone model
    IntArray interfaceNum; // starting and ending interfaceNum for the delamination
    IntArray crossSectionNum;     // use this to keep track of the interface coordinates
    int matNum; // still used??
    double delamXiCoord;    // defines at what local xi-coord the delamination is defined

    // New 110814 JB defines between what local xi-coords the delamination is defined
    double xiBottom;
    double xiTop;
    
    // Adaptive enrichment, 
    // cf. Främby, Fagerström & Bouzoulis, 'Adaptive modelling of delamination initiation and propagation using an equivalent single-layer shell approach', IJNME, 2016 
    double initiationFactor;   // knock-down factor on initiation values, (0,1] //JF
    double initiationRadius;   // radius around around newlye initiated element nodes to be included //JF
    bool recoverStresses;      // recover tranverse stresses using momentum balance (default). //JF

public:
    Delamination(int n, XfemManager *xm, Domain *aDomain);
    virtual ~Delamination() { }

    void initializeFrom(InputRecord &ir) override;
    int instanciateYourself(DataReader &dr) override;
    void postInitialize() override;
    void appendInputRecords(DynamicDataReader &oDR) override;

    const char *giveClassName() const override { return "Delamination"; }
    const char *giveInputRecordName() const override { return _IFT_Delamination_Name; }

    double giveDelamXiCoord() const { return xiBottom; }     // coord where the delamination is defined
    double giveBoundingDelamXiCoord() const { return xiTop; } // coord where the delamination enrichment should stop, default is the shell surface
    int giveDelamInterfaceNum() const { return interfaceNum.at(1); }
    IntArray giveDelamCrossSectionNum() const { return crossSectionNum; }
    double giveInitiationFactor() const { return initiationFactor; }
    //void updateGeometry(FailureCriteriaStatus *fc, TimeStep *tStep) override;
    bool hasPropagatingFronts() const override { return true; }
    bool hasInitiationCriteria() override;

    void propagateFronts(bool &oFrontsHavePropagated) override;
    void findInitiationFronts(bool &failureChecked, const IntArray &CSnumbers, std :: vector< IntArray > &CSinterfaceNumbers, std :: vector< IntArray > &CSDofManNumbers, std :: vector< FloatArray > &initiationFactors, TimeStep *tStep);

    void evaluateEnrFuncInNode(std :: vector< double > &oEnrFunc, const Node &iNode) const override { OOFEM_ERROR("Not implemented.") }

    void evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iGlobalCoord, const FloatArray &iLocalCoord, int iNodeInd, const Element &iEl) const override;
    void evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iGlobalCoord, const FloatArray &iLocalCoord, int iNodeInd, const Element &iEl, const FloatArray &iN, const IntArray &iElNodes) const override;

    void evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const FloatArray &iGlobalCoord, const FloatArray &iLocalCoord, int iNodeInd, const Element &iEl) const override { OOFEM_ERROR("Not implemented."); }
    void evaluateEnrFuncDerivAt(std :: vector< FloatArray > &oEnrFuncDeriv, const FloatArray &iGlobalCoord, const FloatArray &iLocalCoord, int iNodeInd, const Element &iEl, const FloatArray &iN, const FloatMatrix &idNdX, const IntArray &iElNodes) const override { OOFEM_ERROR("Not implemented."); }

    void evalLevelSetNormal(double &oLevelSet, const FloatArray &iGlobalCoord, const FloatArray &iN, const IntArray &iNodeInd) const override;
    void evalLevelSetTangential(double &oLevelSet, const FloatArray &iGlobalCoord, const FloatArray &iN, const IntArray &iNodeInd) const override { OOFEM_ERROR("Not implemented."); }
    void evalGradLevelSetNormal(FloatArray &oGradLevelSet, const FloatArray &iGlobalCoord, const FloatMatrix &idNdX, const IntArray &iNodeInd) const override { OOFEM_ERROR("Not implemented."); }

    void evaluateEnrFuncAt(std :: vector< double > &oEnrFunc, const FloatArray &iPos, const double &iLevelSet) const;
};
} // end namespace oofem


#endif /* DELAMINATION_H_ */
