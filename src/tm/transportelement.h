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

#ifndef transportelement_h
#define transportelement_h

#include "element.h"
#include "floatmatrix.h"
#include "primaryfield.h"
#include "matresponsemode.h"

namespace oofem {
class TransportCrossSection;

/**
 * This abstract class represent a general base element class for transport problems.
 * In actual implementation, the same approximation order of all unknowns is assumed,
 * but this can be easily implemented.
 */
class TransportElement : public Element, public EIPrimaryFieldInterface
{
public:
    enum ElementMode { HeatTransferEM, HeatMass1TransferEM, Mass1TransferEM };

protected:
    ElementMode emode;
    /// Stefan–Boltzmann constant W/m2/K4
    static const double stefanBoltzmann;

public:
    TransportElement(int n, Domain * d, ElementMode em = HeatTransferEM);
    virtual ~TransportElement();

    virtual void giveCharacteristicMatrix(FloatMatrix &answer, CharType type, TimeStep *tStep);
    virtual void giveCharacteristicVector(FloatArray &answer, CharType type, ValueModeType mode, TimeStep *tStep);

    virtual void giveDofManDofIDMask(int inode, IntArray &answer) const;

    virtual void computeInternalForcesVector(FloatArray &answer, TimeStep *tStep);
    virtual void computeExternalForcesVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode);
    virtual void computeInertiaForcesVector(FloatArray &answer, TimeStep *tStep);
    virtual void computeLumpedCapacityVector(FloatArray &answer, TimeStep *tStep);

    virtual void computeLoadVector(FloatArray &answer, Load *load, CharType type, ValueModeType mode, TimeStep *tStep);
    virtual void computeBoundaryLoadVector(FloatArray &answer, BoundaryLoad *load, int boundary, CharType type, ValueModeType mode, TimeStep *tStep);
    virtual void computeTangentFromBoundaryLoad(FloatMatrix &answer, BoundaryLoad *load, int boundary, MatResponseMode rmode, TimeStep *tStep);
    virtual void computeBoundaryEdgeLoadVector(FloatArray &answer, BoundaryLoad *load, int edge, CharType type, ValueModeType mode, TimeStep *tStep);

    //Compute radiative heat tranfer coefficient as 4*eps*sigma*T_infty^3
    virtual double getRadiativeHeatTranferCoef(BoundaryLoad *bLoad, TimeStep *tStep);
    
    virtual void computeField(ValueModeType mode, TimeStep *tStep, const FloatArray &lcoords, FloatArray &answer);

    TransportCrossSection * giveTransportCrossSection();
    virtual Material * giveMaterial();

    /**
     * Gives the thickness at some global coordinate.
     * For solid elements, the value returned is 1.0 (which is the default implementation).
     * @param gcoords Global coordinates.
     * @return Thickness of element at given coordinate.
     * @todo Move this into the base element?
     */
    virtual double giveThicknessAt(const FloatArray &gcoords) { return 1.0; }

    /** Computes the capacity matrix of the receiver */
    virtual void computeCapacityMatrix(FloatMatrix &answer, TimeStep *tStep);
    /** Computes the conductivity matrix of the receiver */
    virtual void computeConductivityMatrix(FloatMatrix &answer, MatResponseMode rMode, TimeStep *tStep);
    /** Computes the RHS contribution to balance equation(s) due to boundary conditions */
    virtual void computeBCVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode);
    /** Computes the LHS contribution to balance equation(s) due to boundary conditions */
    virtual void computeBCMtrxAt(FloatMatrix &answer, TimeStep *tStep, ValueModeType mode);
    /** Computes the contribution to balance equation(s) due to internal sources */
    virtual void computeInternalSourceRhsVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode);
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
    virtual int checkConsistency();

    virtual int EIPrimaryFieldI_evaluateFieldVectorAt(FloatArray &answer, PrimaryField &pf,
                                                      FloatArray &coords, IntArray &dofId, ValueModeType mode,
                                                      TimeStep *tStep);

#ifdef __OOFEG
    int giveInternalStateAtNode(FloatArray &answer, InternalStateType type, InternalStateMode mode,
                                int node, TimeStep *tStep);
    // Graphics output
    //virtual void drawRawGeometry(oofegGraphicContext &gc, TimeStep *tStep) {}
    //virtual void drawDeformedGeometry(oofegGraphicContext &gc, TimeStep *tStep, UnknownType) {}
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
     * @param iri Index of integration rule to use.
     * @param rmode Determines the constitutive sub matrix to be used.
     * @param tStep Solution step.
     */
    virtual void computeConductivitySubMatrix(FloatMatrix &answer, int iri, MatResponseMode rmode, TimeStep *tStep);
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
    void computeBodyBCSubVectorAt(FloatArray &answer, Load *load, TimeStep *tStep, ValueModeType mode, int indx);
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

    /**
     * Computes the basis functions.
     * @param answer The basis functions evaluated at lcoord.
     * @param lcoord The local coordinate.
     */
    virtual void computeNAt(FloatArray &answer, const FloatArray &lcoord);
    /**
     * Computes the interpolation matrix corresponding to all unknowns.
     * In the default implementation the same approximation order is assumed, but it can be extended.
     */
    virtual void computeNmatrixAt(FloatMatrix &answer, const FloatArray &lcoords);
    virtual void computeBmatrixAt(FloatMatrix &answer, const FloatArray &lcoords);
    /**
     * Computes the gradient matrix corresponding to one unknown.
     */
    virtual void computeGradientMatrixAt(FloatMatrix &answer, const FloatArray &lcoords);
    /**
     * Computes the contribution to balance equation(s) due to internal sources
     */
    virtual void computeInternalSourceRhsSubVectorAt(FloatArray &answer, TimeStep *, ValueModeType mode, int indx);

    /**
     * Computes the basis functions at the edge for one unknown.
     */
    virtual void computeEgdeNAt(FloatArray &answer, int iEdge, const FloatArray &lcoord);
    /**
     * Gives the node indexes for given edge.
     */
    virtual void giveEdgeDofMapping(IntArray &mask, int iEdge);
    /**
     * Computes the length around a integration point on a edge.
     */
    virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge) = 0;
    virtual void computeEdgeIpGlobalCoords(FloatArray &answer, const FloatArray &lcoord, int iEdge);

    virtual IntegrationRule *GetSurfaceIntegrationRule(int approxOrder) { return NULL; }
    virtual void computeSurfaceNAt(FloatArray &answer, int iSurf, const FloatArray &lcoord);
    virtual double computeSurfaceVolumeAround(GaussPoint *gp, int iSurf) { return 0.; }
    virtual void giveSurfaceDofMapping(IntArray &mask, int iSurf);
    virtual void computeSurfIpGlobalCoords(FloatArray &answer, const FloatArray &lcoord, int iSurf);

    virtual int giveApproxOrder(int unknownIndx);

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
     */
    void assembleLocalContribution(FloatMatrix &answer, FloatMatrix &src, int ndofs, int rdof, int cdof);
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
     */
    void assembleLocalContribution(FloatArray &answer, FloatArray &src, int ndofs, int rdof);
};
} // end namespace oofem
#endif // transportelement_h
