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
 *               Copyright (C) 1993 - 2014   Borek Patzak
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

#ifndef structuralelements_h
#define structuralelements_h

#include "structuralelementevaluator2.h"
#include "elementgeometry.h"
#include "baseelement.h"
#include "gausspoint.h"

#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "spatiallocalizer.h"
#include "directerrorindicatorrc.h"
#include "eleminterpmapperinterface.h"
#include "zzerrorestimator.h"
#include "mmashapefunctprojection.h"
#include "huertaerrorestimator.h"


namespace oofem {
class Domain;
class InputRecord;

class OneDimensionalStructuralElement : public OneDimensionalStructuralElementEvaluator, public ElementGeometry, public BaseElement	
{
	public:
		OneDimensionalStructuralElement(int n, Domain *aDomain);
		~OneDimensionalStructuralElement(){;}
		
		virtual double computeVolumeAround(GaussPoint *gp);

		virtual int giveApproxOrder(){return 1;}//this->giveInterpolation()->giveInterpolationOrder();}
//		virtual IRResultType initializeFrom(InputRecord *ir){ElementGeometry :: initializeFrom(ir);}
                                                                      


    

};

class PlaneStressElement : public PlaneStressElementEvaluator2, public ElementGeometry, public BaseElement
{
	public:
		PlaneStressElement(int n, Domain *aDomain);
		~PlaneStressElement(){;}
	
		virtual double computeVolumeAround(GaussPoint *gp);
		virtual int giveApproxOrder(){return 1;}//this->giveInterpolation()->giveInterpolationOrder();}
//		virtual IRResultType initializeFrom(InputRecord *ir){ElementGeometry :: initializeFrom(ir);}

			/**
		 * @name Edge load support
		 */
		//@{
		virtual void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp);
		//virtual void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
		virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
		virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
		//virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp);
		//@}

};

class PlaneStrainElement : public PlaneStrainElementEvaluator, public ElementGeometry, public BaseElement
{
	public:
		PlaneStrainElement(int n, Domain *aDomain);
		~PlaneStrainElement(){;}
	
		virtual double computeVolumeAround(GaussPoint *gp);
		virtual int giveApproxOrder(){return 1;}//this->giveInterpolation()->giveInterpolationOrder();}
//		virtual IRResultType initializeFrom(InputRecord *ir){ElementGeometry :: initializeFrom(ir);}

			/**
		 * @name Edge load support
		 */
		//@{
		virtual void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp);
		//virtual void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
		virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
		virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
		//virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp);
		//@}

};

class AxisymmetricStructuralElement : public AxisymmetricStructuralElementEvaluator, public ElementGeometry, public BaseElement
{
	public:
		AxisymmetricStructuralElement(int n, Domain *aDomain);
		~AxisymmetricStructuralElement(){;}
	
		virtual double computeVolumeAround(GaussPoint *gp);
		virtual int giveApproxOrder(){return 1;}//this->giveInterpolation()->giveInterpolationOrder();}
//		virtual IRResultType initializeFrom(InputRecord *ir){ElementGeometry :: initializeFrom(ir);}

			/**
		 * @name Edge load support
		 */
		//@{
		virtual void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp);
		//virtual void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
		virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
		virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
		//virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp);
		//@}


};

class ThreeDimensionalStructuralElement : public ThreeDimensionalStructuralElementEvaluator, public ElementGeometry, public BaseElement
{
	public:
		ThreeDimensionalStructuralElement(int n, Domain *aDomain);
		~ThreeDimensionalStructuralElement(){;}
		
		virtual double computeVolumeAround(GaussPoint *gp);


		/**
		 * @name Edge load support
		 */
		//@{
		virtual void computeEgdeNMatrixAt(FloatMatrix &answer, int iedge, GaussPoint *gp);
		//virtual void giveEdgeDofMapping(IntArray &answer, int iEdge) const;
		virtual double computeEdgeVolumeAround(GaussPoint *gp, int iEdge);
		virtual void computeEdgeIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iEdge);
		//virtual int computeLoadLEToLRotationMatrix(FloatMatrix &answer, int iEdge, GaussPoint *gp);
		//@}

		/**
		 * @name Surface load support
		 */
		//@{
		virtual void computeSurfaceNMatrixAt(FloatMatrix &answer, int iSurf, GaussPoint *gp);
		virtual void giveSurfaceDofMapping(IntArray &answer, int iSurf) const;
		//virtual IntegrationRule *GetSurfaceIntegrationRule(int iSurf);
		virtual double computeSurfaceVolumeAround(GaussPoint *gp, int iSurf);
		virtual void computeSurfIpGlobalCoords(FloatArray &answer, GaussPoint *gp, int iSurf);
		//virtual int computeLoadLSToLRotationMatrix(FloatMatrix &answer, int iSurf, GaussPoint *gp);
		//@}


		virtual int giveApproxOrder(){return 1;}//this->giveInterpolation()->giveInterpolationOrder();}
		//virtual IRResultType initializeFrom(InputRecord *ir){ElementGeometry :: initializeFrom(ir);}

};

} // end namespace oofem
#endif // structuralelements_h
