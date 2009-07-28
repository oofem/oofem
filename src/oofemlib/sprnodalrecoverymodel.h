/* $Header: /home/cvs/bp/oofem/oofemlib/src/sprnodalrecoverymodel.h,v 1.2 2003/04/06 14:08:26 bp Exp $ */
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
 *               Copyright (C) 1993 - 2008   Borek Patzak
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

//   **************************************
//   *** CLASS SPR NODAL RECOVERY MODEL ***
//   **************************************

#ifndef sprnodalrecoverymodel_h
#define sprnodalrecoverymodel_h

#include "compiler.h"

#include "nodalrecoverymodel.h"
#include "interface.h"

class GaussPoint;
class SPRNodalRecoveryModelInterface;

enum SPRPatchType { SPRPatchType_2dxy, SPRPatchType_3dBiLin, SPRPatchType_2dquadratic };

/**
 * The Superconvergent Patch Recovery (SPR) nodal recovery model is based on paper of Zienkiewicz and Zhu
 * "The Superconvergent Patch recovery and Posteriori Error Estimates. Part 1: The Recovery Technique",
 * Int. Journal for Num. Meth in Engng, vol. 33, 1331-1364, 1992.
 * The recovery uses local discrete least square smoothing over an element patch surrounding the particular
 * node considered.
 */
class SPRNodalRecoveryModel : public NodalRecoveryModel
{
protected:
  /** Helper structure to pass required argumets to packing/unpacking functions
      needed in parallel mode */
  struct parallelStruct {
    FloatArray* dofManValues;
    IntArray* dofManPatchCount;
    IntArray* regionNodalNumbers;
    int regionValSize;
    parallelStruct (FloatArray* a, IntArray* b, IntArray* c, int d) {
      dofManValues=a; dofManPatchCount=b; regionNodalNumbers=c; regionValSize=d;}
  };

public:
    /// Constructor
    SPRNodalRecoveryModel(Domain *d);
    /// Destructor
    ~SPRNodalRecoveryModel();
    /** Recovers the nodal values for all regions of given Domain.
     * @param d domain of interest
     * @param type determines the type of internal variable to be recovered
     * @param tStep time step
     */
    int recoverValues(InternalStateType type, TimeStep *tStep);
private:
    /**
     * Initializes the region table indicating regions to skip.
     * @param regionMap region tabl, the nonzero entry for region indicates region to skip due to
     * unsupported elements or incompatible value size
     * @param regionValSize contains the record size for each region
     * @param type determines the type of internal variable to be recovered
     */
    void initRegionMap(IntArray &regionMap, IntArray &regionValSize, IntArray &regionTypes, InternalStateType type);


    void determinePatchAssemblyPoints(IntArray &pap, int ireg, SPRPatchType regType);
    void initPatch(IntArray &patchElems, IntArray &dofManToDetermine, IntArray &pap, int papNumber, int ireg);
    void computePatch(FloatMatrix &a, IntArray &patchElems, int papNumber, int regionValSize,
                      SPRPatchType regType, InternalStateType type, TimeStep *tStep);
    void determineValuesFromPatch(FloatArray &dofManValues, IntArray &dofManCount,
                                  IntArray &regionNodalNumbers, IntArray &dofManToDetermine,
                                  FloatMatrix &a, int papNumber, int regionValSize,
                                  SPRPatchType type);
    void computePolynomialTerms(FloatArray &P, FloatArray &coords, SPRPatchType type);
    int  giveNumberOfUnknownPolynomialCoefficients(SPRPatchType regType);

#ifdef __PARALLEL_MODE
    void initCommMaps ();
    void exchangeDofManValues   (int ireg, FloatArray& dofManValues, 
				 IntArray& dofManPatchCount, IntArray& regionNodalNumbers,
				 int regionValSize);
    int  packSharedDofManData   (parallelStruct* s, ProcessCommunicator &processComm);
    int  unpackSharedDofManData (parallelStruct* s, ProcessCommunicator &processComm);
#endif
};

/**
 * The element interface required by ZZNodalRecoveryModel.
 */
class SPRNodalRecoveryModelInterface : public Interface
{
public:
    /// Constructor
    SPRNodalRecoveryModelInterface() { }


    virtual void SPRNodalRecoveryMI_giveSPRAssemblyPoints(IntArray &pap) = 0;
    virtual void SPRNodalRecoveryMI_giveDofMansDeterminedByPatch(IntArray &answer, int pap) = 0;
    virtual int SPRNodalRecoveryMI_giveDofManRecordSize(InternalStateType type) = 0;
    virtual int SPRNodalRecoveryMI_giveNumberOfIP() = 0;
    //virtual void SPRNodalRecoveryMI_giveIPValue (FloatArray& answer, int ipNum, InternalStateType type) = 0;
    virtual void SPRNodalRecoveryMI_computeIPGlobalCoordinates(FloatArray &coords, GaussPoint *gp) = 0;
    virtual SPRPatchType SPRNodalRecoveryMI_givePatchType() = 0;
};


#endif // sprnodalrecoverymodel_h






