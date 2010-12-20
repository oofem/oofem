/* $Header: /home/cvs/bp/oofem/tm/src/transportelement.h,v 1.3 2003/04/23 14:22:15 bp Exp $ */
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

//   **********************************************************
//   *** CLASS GENERAL ELEMENT CLASS FOR TRANSPORT PROBLEMS ***
//   **********************************************************


#ifndef transportelement_h
#define transportelement_h


#include "element.h"
#include "femcmpnn.h"
#include "domain.h"
#include "flotmtrx.h"

#include "primaryfield.h"
#include "matresponsemode.h"

namespace oofem {
class TimeStep;
class Node;
class Material;
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;

/**
 * This abstract class represent a general base element class for
 * transport problems.
 * In actual implementation, the same approximation order of all unknowns is assumed,
 * but this can be easily implemented.
 */
class TransportElement : public Element, public EIPrimaryFieldInterface
{
public:
    enum ElementMode { HeatTransferEM, HeatMass1TransferEM };
protected:
    ElementMode emode;

public:
    // constructor
    TransportElement(int, Domain *, ElementMode em = HeatTransferEM);
    ~TransportElement();                        // destructor

    // characteristic  matrix
    void giveCharacteristicMatrix(FloatMatrix & answer, CharType, TimeStep *);
    void giveCharacteristicVector(FloatArray & answer, CharType, ValueModeType, TimeStep *);

    /** Computes the capacity matrix of the receiver */
    virtual void computeCapacityMatrix(FloatMatrix &answer, TimeStep *);
    /** Computes the conductivity matrix of the receiver */
    virtual void computeConductivityMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    /** Computes the RHS contribution to balance equation(s) due to boundary conditions */
    virtual void computeBCVectorAt(FloatArray &answer, TimeStep *, ValueModeType mode);
    /* Computes the LHS contribution to balance equation(s) due to boundary conditions */
    virtual void computeBCMtrxAt(FloatMatrix &answer, TimeStep *, ValueModeType mode);
    /** Computes the contribution to balance equation(s) due to internal sources */
    virtual void computeInternalSourceRhsVectorAt(FloatArray &answer, TimeStep *, ValueModeType mode) = 0;
    /** Computes the LHS contribution to balance equation(s) due to material internal source */
    virtual void computeIntSourceLHSMatrix(FloatMatrix &answer, TimeStep *tStep);
    /** Computes the part of internal source LHS contribution corresponding to unknown identified by rmode parameter */
    virtual void computeIntSourceLHSSubMatrix(FloatMatrix &answer, MatResponseMode rmode, int iri, TimeStep *tStep);
    /**
     * Computes a flow vector in an integration point
     * @param answer flow vector
     * @param gp integration point
     * @param stepN time step
     */
    virtual void computeFlow(FloatArray &answer, GaussPoint *gp, TimeStep *stepN);

    // time step termination
    /**
     * Updates a state vector in each integration point of element
     * @param tStep finished time step
     */
    void                  updateInternalState(TimeStep *);
    void                  printOutputAt(FILE *, TimeStep *);
    virtual int    checkConsistency();

    // definition
    const char *giveClassName() const { return "TransportElement"; }
    classType                giveClassID() const { return TransportElementClass; }

    virtual void giveElementDofIDMask(EquationID, IntArray & answer) const;

    /**
     * Evaluates the value of field at given point of interest (should be located inside receiver's volume) using
     * element interpolation.
     */
    virtual int  EIPrimaryFieldI_evaluateFieldVectorAt(FloatArray &answer, PrimaryField &pf,
                                                       FloatArray &coords, IntArray &dofId, ValueModeType mode,
                                                       TimeStep *atTime);

    /**
     * Returns the mask of reduced indexes of Internal Variable component .
     * @param answer mask of Full VectorSize, with components beeing the indexes to reduced form vectors.
     * @param type determines the internal variable requested (physical meaning)
     * @returns nonzero if ok or error is generated for unknown mat mode.
     */
    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type);

#ifdef __OOFEG
    int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                int node, TimeStep *atTime);
    //
    // Graphics output
    //
    //void          drawYourself (oofegGraphicContext&);
    //virtual void  drawRawGeometry (oofegGraphicContext&) {}
    //virtual void  drawDeformedGeometry(oofegGraphicContext&, UnknownType) {}
#endif

protected:
    /**
     * Computes constitutive matrix of receiver.
     * @param answer computed answer
     * @param rMode material response mode of answer
     * @param gp integration point for which constitutive matrix is computed
     * @param tStep time step
     */
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer,
                                             MatResponseMode rMode, GaussPoint *,
                                             TimeStep *tStep);
    /** computes the matrix $P=\int_\Omega N^T N\;d\Omega$
     * the result should be multiplied by corresponding coefficint and localized
     * @param answer contains the result
     * @param rmode determines the capacity coeff to be used.
     * @param iri index of integration rule to use
     * @param tStep time step
     */
    virtual void  computeCapacitySubMatrix(FloatMatrix &answer, MatResponseMode rmode, int iri, TimeStep *tStep);
    /** computes the matrix $C=\int_\Omega \sum_i \left( dN\over x_i\right)^T dN\over x_i\right d\Omega$
     * the result should be multiplied by corresponding coefficint and localized
     * @param answer contains the result
     * @param nsd number of spatial dimension
     * @param iri index of integration rule to use
     * @param rmode determines the constitutive submatrix to be used.
     * @param tStep solution step
     */
    virtual void  computeConductivitySubMatrix(FloatMatrix &answer, int nsd, int iri, MatResponseMode rmode, TimeStep *tStep);
    /** computes the part of RHS due to applied BCs. The part corresponds to unknown identified by indx param.
     * the result should be localized. The part corresponding to Dirichlet BC is not taken into account.
     * @param answer contains the result
     * @param tStep solution step
     * @param mode CharTypeMode
     * @param indx unknown index
     */
    void computeBCSubVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode, int indx);
    /** computes the part of LHS due to applied BCs.
     * the result should be localized.
     * @param answer contains the result
     * @param tStep solution step
     * @param mode CharTypeMode
     * @param indx unknown index
     */
    void computeBCSubMtrxAt(FloatMatrix &answer, TimeStep *, ValueModeType mode, int indx);
    /** computes the part of RHS due to applied BCs on particular edge.
     * The part corresponds to unknown identified by indx param.
     * the result should be localized. The part corresponding to Dirichlet BC is not taken into account.
     * @param answer contains the result
     * @param load edge load object
     * @param iEdge element edge number subjected to bc
     * @param tStep solution step
     * @param mode CharTypeMode
     * @param indx unknown index
     */
    void computeEdgeBCSubVectorAt(FloatArray &answer, Load *load, int iSurf,
                                  TimeStep *tStep, ValueModeType mode, int indx);

    /** computes the part of RHS due to applied BCs on particular surface.
     * The part corresponds to unknown identified by indx param.
     * the result should be localized. The part corresponding to Dirichlet BC is not taken into account.
     * @param answer contains the result
     * @param load surface load object
     * @param iSurf element surface number subjected to bc
     * @param tStep solution step
     * @param mode CharTypeMode
     * @param indx unknown index
     */
    void computeSurfaceBCSubVectorAt(FloatArray &answer, Load *load, int iSurf,
                                     TimeStep *tStep, ValueModeType mode, int indx);
    /*
     * Computes the load vector due to the Dirichlet boundary conditions acting on the
     * receiver's nodes, at stepN.
     */
    //void  computeDirichletBcRhsVectorAt (FloatArray& answer, TimeStep* stepN, CharTypeMode mode);
    virtual void  computeGradientMatrixAt(FloatMatrix &answer, GaussPoint *)  = 0;
    virtual void  computeNmatrixAt(FloatMatrix &n, FloatArray *)  = 0;
    /* computes the submatrix of interpolation matrix cooresponding to single unknown.
     * currently, the same approximation order is assumed, but it can be extended.
     */
    virtual void  computeNSubMatrixAt(FloatMatrix &n, FloatArray *)  = 0;
    /** Computes the contribution to balance equation(s) due to internal sources */
    virtual void          computeInternalSourceRhsSubVectorAt(FloatArray &answer, TimeStep *, ValueModeType mode, int indx);

    virtual void computeEgdeNMatrixAt(FloatMatrix &n, GaussPoint *gp) = 0;
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge) = 0;
    virtual void giveEdgeDofMapping(IntArray &mask, int iEdge) = 0;
    virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge) = 0;

    virtual IntegrationRule *GetSurfaceIntegrationRule(int approxOrder) { return NULL; }
    virtual void computeSurfaceNMatrixAt(FloatMatrix &answer, GaussPoint *gp) { answer.resize(0, 0); }
    virtual double computeSurfaceVolumeAround(GaussPoint *gp, int iEdge) { return 0.; }
    virtual void giveSurfaceDofMapping(IntArray &mask, int iEdge)  { mask.resize(0); }
    virtual void computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iSurf)  { answer.resize(0); }

    virtual int giveApproxOrder(int unknownIndx) = 0;

    /**
     * Assembles the given source matrix of size (ndofs, ndofs) into target matrix answer.
     * The coefficients in src matrix are assumed to represent local sub-matrix for specific DOF.
     * The answer matrix is element global one, and it typically contatins the values for different DOFs.
     * In this routine is assumed, that the number of dofs per node is the same.
     * The DOFs in answer are ordered in a same way as the DOFS at the element level are
     * (all Dofs corresponding to given node subsequently, followed by next node DOFS).
     * @param answer receiver matrix
     * @param src source mtrx containing local submatrix corresponding to single DOF
     * @param ndofs number of DOFs per node (assumed same for all nodes)
     * @param rdof rows of src are localized into rows of answer corresponding to rdof-th dof
     * @param cdof columns of src are localized into columns of answer corresponding to cdof-th dof
     * @param coeff all coefficients of src are multiplied by coeff before localized
     */
    void assembleLocalContribution(FloatMatrix &answer, FloatMatrix &src, int ndofs, int rdof, int cdof, double coeff);
    /**
     * Assembles the given source vector of size (ndofs) into target answer.
     * The coefficients in src vector are assumed to represent local sub-vector for specific DOF.
     * The answer vector is element global one, and it typically contatins the values for different DOFs.
     * In this routine is assumed, that the number of dofs per node is the same.
     * The DOFs in answer are ordered in a same way as the DOFS at the element level are
     * (all Dofs corresponding to given node subsequently, followed by next node DOFS).
     * @param answer receiver vector
     * @param src source vector containing local subvector corresponding to single DOF
     * @param ndofs number of DOFs per node (assumed same for all nodes)
     * @param rdof rows of src are localized into rows of answer corresponding to rdof-th dof
     * @param coeff all coefficients of src are multiplied by coeff before localized
     */
    void assembleLocalContribution(FloatArray &answer, FloatArray &src, int ndofs, int rdof, double coeff);
};
} // end namespace oofem
#endif // transportelement_h
