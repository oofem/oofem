/*
 * Q1TrBase.h
 *
 *  Created on: Feb 25, 2010
 *      Author: carl
 */

#ifndef Q1TRBASE_H_
#define Q1TRBASE_H_
#include "transportelement.h"
#include "fei2dtrlin.h"
#include "domain.h"
#include "nodalaveragingrecoverymodel.h"
#include "transportmaterial.h"

namespace oofem {

/** Element class for the Darcyflow engineering model. Linear description of seepage */
class Tr1Darcy : public TransportElement, public NodalAveragingRecoveryModelInterface
{
protected:
/*    double b [ 3 ];
    double c [ 3 ];*/
	double area;
	int numberOfGaussPoints;
	static FEI2dTrLin interpolation_lin;

public:

	virtual void  computeGradientMatrixAt(FloatMatrix &answer, GaussPoint *); // Empty function. Temporary workaround.
    virtual void computeInternalSourceRhsVectorAt(FloatArray &answer, TimeStep *, ValueModeType mode) { return; };	// Empty function. Temporary workaround.
    virtual void  computeNmatrixAt(FloatMatrix &n, FloatArray *)  { return; }; // Empty function. Temporary workaround.
    virtual void  computeNSubMatrixAt(FloatMatrix &n, FloatArray *) { return; }; // Empty function. Temporary workaround.
    virtual void  computeEgdeNMatrixAt(FloatMatrix &n, GaussPoint *gp) { return; }; // Empty function. Temporary workaround.
    virtual void  giveEdgeDofMapping(IntArray &mask, int iEdge) { return; }; // Empty function. Temporary workaround.
    virtual void  computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge) { return; }; // Empty function. Temporary workaround.
    virtual int   giveApproxOrder(int unknownIndx) { return 0; }; // Empty function. Temporary workaround.	virtual void computeInternalSourceRhsVectorAt(FloatArray &answer, TimeStep *tStep, ValueModeType mode) {};

	Tr1Darcy(int, Domain *);
	virtual ~Tr1Darcy ();
	IRResultType initializeFrom(InputRecord *ir);
	void initGeometry();

	void giveDofManDofIDMask(int inode, EquationID ut, IntArray &answer) const;
	void giveCharacteristicVector(FloatArray &answer, CharType mtrx, ValueModeType mode,TimeStep *tStep);
	void giveCharacteristicMatrix(FloatMatrix &answer, CharType mtrx, TimeStep *tStep);
	void computeStiffnessMatrix(FloatMatrix &answer, TimeStep *atTime);
	void computeLoadVector(FloatArray &answer, TimeStep *atTime);
	void computeGaussPoints();
	int computeNumberOfDofs(EquationID ut);
	double computeVolume() { return this->area; }

	void computeEdgeBCSubVectorAt(FloatArray &answer, Load *load, int iEdge, TimeStep *tStep);
	void computeInternalForcesVector (FloatArray &answer, TimeStep *atTime);

	double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
	void giveEdgeDofMappingV(IntArray &answer, int iEdge);
	void giveNEdge_xieta(FloatMatrix &answer, GaussPoint *aGaussPoint);

	Element_Geometry_Type giveGeometryType() const { return EGT_triangle_1; }

	// From NodalAveragingRecoveryModelInterface
	const char *giveClassName() const { return "Tr1Darcy"; };
	virtual void NodalAveragingRecoveryMI_computeNodalValue(FloatArray &answer, int node, InternalStateType type, TimeStep *tStep);
	virtual void NodalAveragingRecoveryMI_computeSideValue(FloatArray &answer, int side, InternalStateType type, TimeStep *tStep);
	virtual int NodalAveragingRecoveryMI_giveDofManRecordSize(InternalStateType type) { return 2; }
	Interface * giveInterface(InterfaceType interface);
};

}

#endif /* Q1TRBASE_H_ */
