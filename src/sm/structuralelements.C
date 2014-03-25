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

#include "structuralelements.h"
#include "feinterpol.h"
#include "crosssection.h"
#include "node.h"

#include "feinterpol2d.h"
#include "feinterpol3d.h"

namespace oofem {

	// 1D Element
	OneDimensionalStructuralElement :: OneDimensionalStructuralElement(int n, Domain *aDomain) : ElementGeometry(n, aDomain), BaseElement(), OneDimensionalStructuralElementEvaluator()
	{}
	
	double OneDimensionalStructuralElement :: computeVolumeAround(GaussPoint *gp) 
	{
		double detJ = abs( this->giveInterpolation()->giveTransformationJacobian( * gp->giveCoordinates(), FEIElementGeometryWrapper(this) ) );
		double weight  = gp->giveWeight();
		return detJ *weight *this->giveCrossSection()->give(CS_Area, gp);
	} 
	
//////////////////////////////////////// PlaneStress Element
	PlaneStressElement :: PlaneStressElement(int n, Domain *aDomain) : ElementGeometry(n, aDomain), BaseElement()
	{}
	

	double PlaneStressElement :: computeVolumeAround(GaussPoint *gp) 
	{
		double determinant, weight, thickness;
		determinant = abs( this->giveInterpolation(displacementInterpolationNumber)->giveTransformationJacobian(* gp->giveCoordinates(),FEIElementGeometryWrapper(this)) );
		weight      = gp->giveWeight();
		thickness   = this->giveCrossSection()->give(CS_Thickness,gp);
		return determinant * weight * thickness;		
	}

	void PlaneStressElement :: computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
	{
		/*
		 *
		 * computes interpolation matrix for element edge.
		 * we assemble locally this matrix for only nonzero
		 * shape functions.
		 * (for example only two nonzero shape functions for 2 dofs are
		 * necessary for linear plane stress quad edge).
		 * These nonzero shape functions are then mapped to
		 * global element functions.
		 *
		 * Using mapping technique will allow to assemble shape functions
		 * without regarding particular side
		 */
				
		FEInterpolation2d *interpolation =  static_cast< FEInterpolation2d *>(this->giveInterpolation());
		int interpOrder = interpolation->giveInterpolationOrder();
		FloatArray n(interpOrder+1);


		interpolation->edgeEvalN( n, iedge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
		IntArray bNodes;
		interpolation->boundaryEdgeGiveNodes(bNodes,iedge);
		int numberBoundaryNodes = bNodes.giveSize();

		answer.resize(2, 2*numberBoundaryNodes);
		answer.zero();
		for (int i = 1; i <= numberBoundaryNodes; i++)
		{
			answer.at(1, 2 * i - 1 ) = n.at(i);
			answer.at(2, 2 * i - 0 ) = n.at(i);

		}
		
	}



double PlaneStressElement :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
	FEInterpolation2d *interpolation =  static_cast< FEInterpolation2d *>(this->giveInterpolation());


	FloatArray c(2);
    this->computeEdgeIpGlobalCoords(c, gp, iEdge);
    double result = interpolation->edgeGiveTransformationJacobian( iEdge, * gp->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this) );
	return c.at(1) * result * gp->giveWeight();

}


void
PlaneStressElement :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
	FEInterpolation2d *interpolation = static_cast< FEInterpolation2d *>(this->giveInterpolation());
    interpolation->edgeLocal2global( answer, iEdge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
}

/////////////////////////////////////////////// PlaneStrain Element

	PlaneStrainElement :: PlaneStrainElement(int n, Domain *aDomain) : ElementGeometry(n, aDomain), BaseElement(), PlaneStrainElementEvaluator()
	{}


	double PlaneStrainElement :: computeVolumeAround(GaussPoint *gp) 
	{
		double determinant, weight, thickness;
		determinant = abs( this->giveInterpolation(displacementInterpolationNumber)->giveTransformationJacobian(* gp->giveCoordinates(),FEIElementGeometryWrapper(this)) );
		weight      = gp->giveWeight();
		thickness   = this->giveCrossSection()->give(CS_Thickness,gp);
		return determinant * weight * thickness;
	}

		void PlaneStrainElement :: computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
	{
		/*
		 *
		 * computes interpolation matrix for element edge.
		 * we assemble locally this matrix for only nonzero
		 * shape functions.
		 * (for example only two nonzero shape functions for 2 dofs are
		 * necessary for linear plane stress quad edge).
		 * These nonzero shape functions are then mapped to
		 * global element functions.
		 *
		 * Using mapping technique will allow to assemble shape functions
		 * without regarding particular side
		 */
				
		FEInterpolation2d *interpolation = static_cast< FEInterpolation2d *>(this->giveInterpolation());
		int interpOrder = interpolation->giveInterpolationOrder();
		FloatArray n(interpOrder+1);


		interpolation->edgeEvalN( n, iedge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
		IntArray bNodes;
		interpolation->boundaryEdgeGiveNodes(bNodes,iedge);
		int numberBoundaryNodes = bNodes.giveSize();

		answer.resize(2, 2*numberBoundaryNodes);
		answer.zero();
		for (int i = 1; i <= numberBoundaryNodes; i++)
		{
			answer.at(1, 2 * i - 1 ) = n.at(i);
			answer.at(2, 2 * i - 0 ) = n.at(i);

		}
		
	}



double PlaneStrainElement :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
	FEInterpolation2d *interpolation = static_cast< FEInterpolation2d *>(this->giveInterpolation());


	FloatArray c(2);
    this->computeEdgeIpGlobalCoords(c, gp, iEdge);
    double result = interpolation->edgeGiveTransformationJacobian( iEdge, * gp->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this) );
	return c.at(1) * result * gp->giveWeight();

}


void
PlaneStrainElement :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
	FEInterpolation2d *interpolation = static_cast< FEInterpolation2d *>(this->giveInterpolation());
    interpolation->edgeLocal2global( answer, iEdge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
}



/////////////////////////////////////////////// Axisymmetric Element
	AxisymmetricStructuralElement :: AxisymmetricStructuralElement(int n, Domain *aDomain) : ElementGeometry(n, aDomain), BaseElement(), AxisymmetricStructuralElementEvaluator()
	{}

	double AxisymmetricStructuralElement :: computeVolumeAround(GaussPoint *gp) 
	{
		double determinant, weight, r, x;
		FloatArray n(4);
		this->giveInterpolation()->evalN(n, *gp->giveCoordinates(), FEIElementGeometryWrapper(this));

		r = 0.;
		for ( int i = 1; i <= numberOfDofMans; i++ ) {
			x  = this->giveNode(i)->giveCoordinate(1);
			r += x * n.at(i);
		}
		
		determinant = abs(this->giveInterpolation()->giveTransformationJacobian(*gp->giveCoordinates(), FEIElementGeometryWrapper(this)));
		weight = gp->giveWeight();
		return determinant * weight * r;
	}

	void AxisymmetricStructuralElement :: computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
	{
		/*
		 *
		 * computes interpolation matrix for element edge.
		 * we assemble locally this matrix for only nonzero
		 * shape functions.
		 * (for example only two nonzero shape functions for 2 dofs are
		 * necessary for linear plane stress quad edge).
		 * These nonzero shape functions are then mapped to
		 * global element functions.
		 *
		 * Using mapping technique will allow to assemble shape functions
		 * without regarding particular side
		 */
				
		FEInterpolation2d *interpolation =  static_cast< FEInterpolation2d *>(this->giveInterpolation());
		int interpOrder = interpolation->giveInterpolationOrder();
		FloatArray n(interpOrder+1);


		interpolation->edgeEvalN( n, iedge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
		IntArray bNodes;
		interpolation->boundaryEdgeGiveNodes(bNodes,iedge);
		int numberBoundaryNodes = bNodes.giveSize();

		answer.resize(2, 2*numberBoundaryNodes);
		answer.zero();
		for (int i = 1; i <= numberBoundaryNodes; i++)
		{
			answer.at(1, 2 * i - 1 ) = n.at(i);
			answer.at(2, 2 * i - 0 ) = n.at(i);

		}
		
	}



double AxisymmetricStructuralElement :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
	FEInterpolation2d *interpolation =  static_cast< FEInterpolation2d *>(this->giveInterpolation());


	FloatArray c(2);
    this->computeEdgeIpGlobalCoords(c, gp, iEdge);
    double result = interpolation->edgeGiveTransformationJacobian( iEdge, * gp->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this) );
	return c.at(1) * result * gp->giveWeight();

}


void
AxisymmetricStructuralElement :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
	FEInterpolation2d *interpolation =  static_cast< FEInterpolation2d *>(this->giveInterpolation());
    interpolation->edgeLocal2global( answer, iEdge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
}


///////////////////////////////////////////// 3D Element/////////////////////////////////////////////////////////////////////////////////
ThreeDimensionalStructuralElement :: ThreeDimensionalStructuralElement(int n, Domain *aDomain) : ElementGeometry(n, aDomain), BaseElement(), ThreeDimensionalStructuralElementEvaluator()
{}

double ThreeDimensionalStructuralElement :: computeVolumeAround(GaussPoint *gp) 
{
		double determinant, weight;
		determinant = abs( this->giveInterpolation(displacementInterpolationNumber)->giveTransformationJacobian(* gp->giveCoordinates(),FEIElementGeometryWrapper(this)) );
		weight      = gp->giveWeight();
		return determinant * weight;
}

void ThreeDimensionalStructuralElement :: computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp)
	{
		/*
		 *
		 * computes interpolation matrix for element edge.
		 * we assemble locally this matrix for only nonzero
		 * shape functions.
		 * (for example only two nonzero shape functions for 2 dofs are
		 * necessary for linear plane stress quad edge).
		 * These nonzero shape functions are then mapped to
		 * global element functions.
		 *
		 * Using mapping technique will allow to assemble shape functions
		 * without regarding particular side
		 */
				
		FEInterpolation3d *interpolation =  static_cast< FEInterpolation3d *>(this->giveInterpolation());
		int interpOrder = interpolation->giveInterpolationOrder();
		FloatArray n(interpOrder+1);


		interpolation->edgeEvalN( n, iedge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
		IntArray bNodes;
		interpolation->boundaryEdgeGiveNodes(bNodes,iedge);
		int numberBoundaryNodes = bNodes.giveSize();

		answer.resize(3, 3*numberBoundaryNodes);
		answer.zero();
		for (int i = 1; i <= numberBoundaryNodes; i++)
		{
			answer.at(1, 3 * i - 2 ) = n.at(i);
			answer.at(2, 3 * i - 1 ) = n.at(i);
			answer.at(3, 3 * i - 0 ) = n.at(i);
		}
		
	}



double ThreeDimensionalStructuralElement :: computeEdgeVolumeAround(GaussPoint *gp, int iEdge)
{
	FEInterpolation3d *interpolation =  static_cast< FEInterpolation3d *>(this->giveInterpolation());

    double result = interpolation->edgeGiveTransformationJacobian( iEdge, * gp->giveCoordinates(),
                                                                       FEIElementGeometryWrapper(this) );
    return result *gp->giveWeight();
}


void
ThreeDimensionalStructuralElement :: computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge)
{
	FEInterpolation3d *interpolation =  static_cast< FEInterpolation3d *>(this->giveInterpolation());
    interpolation->edgeLocal2global( answer, iEdge, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
}


void
ThreeDimensionalStructuralElement :: computeSurfaceNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *sgp)
{
	FEInterpolation3d *interpolation =  static_cast< FEInterpolation3d *>(this->giveInterpolation());
	
	//@todo replace with interpolation->giveNumberOfSurfaceNodes??
	IntArray nodes;
	interpolation->computeLocalSurfaceMapping(nodes, iSurf);
	int nNodesPerSurface = nodes.giveSize();

	FloatArray n(nNodesPerSurface);
    interpolation->surfaceEvalN( n, iSurf, * sgp->giveCoordinates(), FEIElementGeometryWrapper(this) );

	answer.resize(3, 3 * nNodesPerSurface);
    answer.zero();

	for (int i = 1; i <= nNodesPerSurface; i++)
	{
		answer.at(1, 3*i-2) = n.at(i);
		answer.at(2, 3*i-1) = n.at(i);
		answer.at(3, 3*i-0) = n.at(i);
		

	}

}


void
ThreeDimensionalStructuralElement :: giveSurfaceDofMapping(IntArray &answer, int iSurf) const
{
	//@todo check if is it working for all elements, i.e. edge ...
    IntArray nodes;
    const int ndofsn = 3;

	FEInterpolation3d *interpolation =  static_cast< FEInterpolation3d *>(this->giveInterpolation());
		
    interpolation->computeLocalSurfaceMapping(nodes, iSurf);
	int nNodesPerSurface = nodes.giveSize();
	answer.resize(ndofsn*nNodesPerSurface);

	for ( int i = 1; i <= nNodesPerSurface; i++ ) {
        answer.at(i * ndofsn - 2) = nodes.at(i) * ndofsn - 2;
        answer.at(i * ndofsn - 1) = nodes.at(i) * ndofsn - 1;
        answer.at(i * ndofsn) = nodes.at(i) * ndofsn;
    }

}




double
ThreeDimensionalStructuralElement :: computeSurfaceVolumeAround(GaussPoint *gp, int iSurf)
{
	FEInterpolation3d *interpolation =  static_cast< FEInterpolation3d *>(this->giveInterpolation());

    double determinant, weight, volume;
    determinant = abs( interpolation->surfaceGiveTransformationJacobian( iSurf, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) ) );

    weight = gp->giveWeight();
    volume = determinant * weight;

    return volume;
}


void
ThreeDimensionalStructuralElement :: computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int isurf)
{
	FEInterpolation3d *interpolation =  static_cast< FEInterpolation3d *>(this->giveInterpolation());
	interpolation->surfaceLocal2global( answer, isurf, * gp->giveCoordinates(), FEIElementGeometryWrapper(this) );
}




} // end namespace oofem
