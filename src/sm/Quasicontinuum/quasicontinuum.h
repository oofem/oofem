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

#ifndef quasicontinuum_h
#define quasicontinuum_h

#include "floatarray.h"
#include "element.h"
#include "node.h"
#include "qcnode.h"


namespace oofem {
/**
 * General simplification for Quasicontinuum simulation.
 * 
 */
class Quasicontinuum 

{
protected:
    std::vector<IntArray> interpolationMeshNodes;
    IntArray interpolationElementNumbers;
    IntArray interpolationElementIndices;
    std::vector <IntArray*> connectivityTable;
    int nDimensions;
    IntArray elemList;
    IntArray nodeList;


public:
    Quasicontinuum();
    virtual ~Quasicontinuum();

    void setNoDimensions(Domain *d);
    void setupInterpolationMesh(Domain *d, int generateInterpolationElements, int interpolationElementsMaterialNumber, std::vector<IntArray> *newMeshNodes);
    void createInterpolationElements(Domain *d);
    void addCrosssectionToInterpolationElements(Domain *d);

    void applyApproach1(Domain *d);
    void applyApproach2(Domain *d, int homMtrxType, double volumeOfInterpolationMesh);
    void applyApproach3(Domain *d, int homMtrxType);
    
    void homogenizationOfStiffMatrix(double &homogenizedE, double &homogenizedNu, FloatMatrix *Diso );
    void createGlobalStiffnesMatrix(FloatMatrix *Diso, double &S0, Domain *d, int homMtrxType,  double volumeOfInterpolationMesh);

    void computeStiffnessTensorOf1Link(FloatMatrix &D1, double &S0, Element *e, Domain *d);
    bool stiffnessAssignment( std::vector<FloatMatrix*> &individualStiffnessTensors, FloatArray &individialS0, Domain *d, Element *e, qcNode *qn1, qcNode *qn2 );

    void computeIntersectionsOfLinkWithInterpElements( IntArray &intersected, FloatArray &lengths, Domain *d, Element *e, qcNode *qn1, qcNode *qn2);
    bool computeIntersectionsOfLinkWith2DTringleElements( IntArray &intersected, FloatArray &lengths, Domain *d, Element *e, qcNode *qn1, qcNode *qn2);
    bool computeIntersectionsOfLinkWith3DTetrahedraElements( IntArray &intersected, FloatArray &lengths, Domain *d, Element *e, qcNode *qn1, qcNode *qn2);

    void initializeConnectivityTableForInterpolationElements( Domain *d );

    bool intersectionTestSegmentTrianglePlucker3D (FloatArray  &intersectCoords, FloatArray *A, FloatArray *B, FloatArray *C, FloatArray *X1, FloatArray *X2 );

    int intersectionTestSegmentTetrahedra3D (FloatArray  &intersectCoordsX, FloatArray  &intersectCoordsY, FloatArray  &intersectCoordsZ, FloatArray *A, FloatArray *B, FloatArray *C, FloatArray *D, FloatArray *X1, FloatArray *X2 );

    bool intersectionTestSegmentSegment2D (FloatArray  &intersectCoords, FloatArray *A1, FloatArray *A2, FloatArray *B1, FloatArray *B2 );
	
    int intersectionTestSegmentTriangle2D (FloatArray  &intersectCoordsX, FloatArray  &intersectCoordsY, FloatArray *A, FloatArray *B, FloatArray *C, FloatArray *U1, FloatArray *U2 );

    void transformStiffnessTensorToMatrix (FloatMatrix  *matrix, FloatMatrix *tensor );
    
};
} // end namespace oofem
#endif // quasicontinuum_h
