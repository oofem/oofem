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
 *               Copyright (C) 1993 - 2011   Borek Patzak
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

#ifndef transportelement_h
#define transportelement_h

#include "element.h"
#include "flotmtrx.h"

#include "primaryfield.h"
#include "matresponsemode.h"

namespace oofem {
/**
 * This abstract class represent a general base element class for transport problems.
 * In actual implementation, the same approximation order of all unknowns is assumed,
 * but this can be easily implemented.
 * @todo Alot of missing documentation.
 * @todo Nonlinear problems.
 */
class TransportElement : public Element, public EIPrimaryFieldInterface
{
public:
    enum ElementMode { HeatTransferEM, HeatMass1TransferEM };

protected:
    ElementMode emode;

public:
    TransportElement(int n, Domain *d, ElementMode em = HeatTransferEM);
    virtual ~TransportElement();

    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep);
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep);

    /** Computes the capacity matrix of the receiver */
    virtual void computeCapacityMatrix(FloatMatrix &answer, TimeStep *tStep);
    /** Computes the conductivity matrix of the receiver */
    virtual void computeConductivityMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    /** Computes the RHS contribution to balance equation(s) due to boundary conditions */
    virtual void computeBCVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode);
    /** Computes the LHS contribution to balance equation(s) due to boundary conditions */
    virtual void computeBCMtrxAt(FloatMatrix &answer, TimeStep *tStep, ValueModeType mode);
    /** Computes the contribution to balance equation(s) due to internal sources */
    virtual void computeInternalSourceRhsVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode) = 0;
    /** Computes the LHS contribution to balance equation(s) due to material internal source */
    virtual void computeIntSourceLHSMatrix(FloatMatrix &answer, TimeStep *tStep);
    /** Computes the part of internal source LHS contribution corresponding to unknown identified by rmode parameter */
    virtual void computeIntSourceLHSSubMatrix(FloatMatrix &answer, MatResponseMode rmode, int iri, TimeStep *tStep);
    /**
     * Computes a flow vector in an integration point.
     * @param answer Flow vector.
     * @param gp Integration point.
     * @param tStep Time step.
     */
    virtual void computeFlow(FloatArray &answer, GaussPoint *gp, TimeStep *tStep);

    // time step termination
    virtual void updateInternalState(TimeStep *tStep);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    virtual int checkConsistency();

    virtual const char *giveClassName() const { return "TransportElement"; }
    virtual classType giveClassID() const { return TransportElementClass; }

    virtual void giveElementDofIDMask(EquationID, IntArray & answer) const;

    virtual int EIPrimaryFieldI_evaluateFieldVectorAt(FloatArray &answer, PrimaryField &pf,
                                                      FloatArray &coords, IntArray &dofId, ValueModeType mode,
                                                      TimeStep *atTime);

    virtual int giveIntVarCompFullIndx(IntArray &answer, InternalStateType type);

#ifdef __OOFEG
    int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                int node, TimeStep *atTime);
    // Graphics output
    //void drawYourself (oofegGraphicContext&);
    //virtual void drawRawGeometry(oofegGraphicContext&) {}
    //virtual void drawDeformedGeometry(oofegGraphicContext&, UnknownType) {}
#endif

protected:
    /**
     * Computes constitutive matrix of receiver.
     * @param answer Computed answer.
     * @param rMode Material response mode of answer.
     * @param gp Integration point for which constitutive matrix is computed.
     * @param tStep Time step.
     */
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer,
                                             MatResponseMode rMode, GaussPoint *gp,
                                             TimeStep *tStep);
    /**
     * Computes the matrix @f$ P=\int_\Omega N^{\mathrm{T}} N\;\mathrm{d}\Omega @f$.
     * The result should be multiplied by corresponding coefficient and localized.
     * @param answer Contains the result.
     * @param rmode Determines the capacity coefficient to be used.
     * @param iri Index of integration rule to use.
     * @param tStep Time step.
     */
    virtual void computeCapacitySubMatrix(FloatMatrix &answer, MatResponseMode rmode, int iri, TimeStep *tStep);
    /**
     * Computes the matrix @f$ C=\int_\Omega \sum_i \left( \frac{\partial N}{\partial x_i}\right)^{\mathrm{T}} \frac{\partial N}{\partial x_i} \;\mathrm{d}\Omega @f$.
     * The result should be multiplied by corresponding coefficient and localized.
     * @todo Is this expression really correct?
     * @param answer Contains the result.
     * @param nsd Number of spatial dimension.
     * @param iri Index of integration rule to use.
     * @param rmode Determines the constitutive sub matrix to be used.
     * @param tStep Solution step.
     */
    virtual void computeConductivitySubMatrix(FloatMatrix &answer, int nsd, int iri, MatResponseMode rmode, TimeStep *tStep);
    /**
     * Computes the part of RHS due to applied BCs. The part corresponds to unknown identified by indx param.
     * The result should be localized. The part corresponding to Dirichlet BC is not taken into account.
     * @param answer Contains the result.
     * @param tStep Solution step.
     * @param mode Mode of value (incremental, total, ...).
     * @param indx Unknown index.
     */
    void computeBCSubVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode, int indx);
    /**
     * Computes the part of LHS due to applied BCs.
     * the result should be localized.
     * @param answer Contains the result.
     * @param tStep Solution step.
     * @param mode Mode of value (incremental, total, ...).
     * @param indx Unknown index.
     */
    void computeBCSubMtrxAt(FloatMatrix &answer, TimeStep *tStep, ValueModeType mode, int indx);
    /**
     * Computes the part of RHS due to applied BCs on particular edge.
     * The part corresponds to unknown identified by indx param.
     * The result should be localized. The part corresponding to Dirichlet BC is not taken into account.
     * @param answer Contains the result.
     * @param load Edge load object.
     * @param iEdge Element edge number subjected to bc.
     * @param tStep Solution step.
     * @param mode Mode of value (incremental, total, ...).
     * @param indx Unknown index.
     */
    void computeEdgeBCSubVectorAt(FloatArray &answer, Load *load, int iEdge,
                                  TimeStep *tStep, ValueModeType mode, int indx);

    /**
     * Computes the part of RHS due to applied BCs on particular surface.
     * The part corresponds to unknown identified by indx param.
     * The result should be localized. The part corresponding to Dirichlet BC is not taken into account.
     * @param answer Contains the result.
     * @param load Surface load object.
     * @param iSurf Element surface number subjected to bc.
     * @param tStep Solution step.
     * @param mode Mode of value (incremental, total, ...).
     * @param indx Unknown index.
     */
    void computeSurfaceBCSubVectorAt(FloatArray &answer, Load *load, int iSurf,
                                     TimeStep *tStep, ValueModeType mode, int indx);
    /*
     * Computes the load vector due to the Dirichlet boundary conditions acting on the
     * receiver's nodes, at stepN.
     */
    //void  computeDirichletBcRhsVectorAt (FloatArray& answer, TimeStep* stepN, CharTypeMode mode);
    virtual void computeGradientMatrixAt(FloatMatrix &answer, GaussPoint *gp)  = 0;
    /**
     * Computes the sub-matrix of interpolation matrix corresponding to single unknown.
     * Currently, the same approximation order is assumed, but it can be extended.
     */
    virtual void computeNSubMatrixAt(FloatMatrix &n, FloatArray *lcoords)  = 0;
    virtual void computeNmatrixAt(FloatMatrix &n, FloatArray *lcoords)  = 0;
    /**
     * Computes the contribution to balance equation(s) due to internal sources
     */
    virtual void computeInternalSourceRhsSubVectorAt(FloatArray &answer, TimeStep *, ValueModeType mode, int indx);

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
     * The answer matrix is element global one, and it typically contains the values for different DOFs.
     * In this routine is assumed, that the number of dofs per node is the same.
     * The DOFs in answer are ordered in a same way as the DOFS at the element level are
     * (all Dofs corresponding to given node subsequently, followed by next node DOFS).
     * @param answer Receiver matrix.
     * @param src Source matrix containing local sub matrix corresponding to single DOF.
     * @param ndofs Number of DOFs per node (assumed same for all nodes).
     * @param rdof Rows of src are localized into rows of answer corresponding to rdof-th dof.
     * @param cdof Columns of src are localized into columns of answer corresponding to cdof-th dof.
     * @param coeff All coefficients of src are multiplied by coeff before localized.
     */
    void assembleLocalContribution(FloatMatrix &answer, FloatMatrix &src, int ndofs, int rdof, int cdof, double coeff);
    /**
     * Assembles the given source vector of size (ndofs) into target answer.
     * The coefficients in src vector are assumed to represent local sub-vector for specific DOF.
     * The answer vector is element global one, and it typically contains the values for different DOFs.
     * In this routine is assumed, that the number of dofs per node is the same.
     * The DOFs in answer are ordered in a same way as the DOFS at the element level are
     * (all Dofs corresponding to given node subsequently, followed by next node DOFS).
     * @param answer Receiver vector.
     * @param src Source vector containing local subvector corresponding to single DOF.
     * @param ndofs Number of DOFs per node (assumed same for all nodes).
     * @param rdof Rows of src are localized into rows of answer corresponding to rdof-th dof.
     * @param coeff All coefficients of src are multiplied by coeff before localized.
     */
    void assembleLocalContribution(FloatArray &answer, FloatArray &src, int ndofs, int rdof, double coeff);
};
} // end namespace oofem
#endif // transportelement_h
