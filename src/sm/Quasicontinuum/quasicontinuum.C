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
 *               Copyright (C) 1993 - 2025   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; eitherc
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

#include "sm/Quasicontinuum/quasicontinuum.h"
#include "sm/EngineeringModels/qclinearstatic.h"
#include "qcnode.h"
#include "element.h"
#include "mathfem.h"

#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "sm/CrossSections/simplecrosssection.h"
#include "sm/Materials/isolinearelasticmaterial.h"
#include "sm/Materials/anisolinearelasticmaterial.h"
#include "interfacetype.h"
#include "sm/Materials/qcmaterialextensioninterface.h"


namespace oofem {
//REGISTER_Quasicontinuum(QCFullsolveddomain);

Quasicontinuum :: Quasicontinuum()
{}


Quasicontinuum :: ~Quasicontinuum()
{}


void
Quasicontinuum :: setNoDimensions(Domain *d)
// rewrite existing interpolation mesh
{
    int dim = d->giveNumberOfSpatialDimensions();
    if ( dim == 2 || dim == 3 ) {
        this->nDimensions = dim;
    } else {
        OOFEM_ERROR("Invalid number of dimensions. Only 2d and 3d domains are supported in QC simulation. \n");
    }
}


void
Quasicontinuum :: setupInterpolationMesh(Domain *d, int generateInterpolationElements, int interpolationElementsMaterialNumber, std :: vector< IntArray > &newMeshNodes)
// rewrite existing interpolation mesh
{
    int noNewEl = 0;
    if ( generateInterpolationElements == 0 ) {
        // nodes are not needed
        this->interpolationMeshNodes.clear();
        // elements are defined in imput file by material number
        this->interpolationElementNumbers = d->giveElementsWithMaterialNum(interpolationElementsMaterialNumber);
    } else if ( ( generateInterpolationElements == 1 ) || ( generateInterpolationElements == 2 ) ) {
        // nodes were (generated and ) loaded from t3d output file
        this->interpolationMeshNodes = newMeshNodes;
        // elements will be placed in the end of element list
        this->interpolationElementNumbers.clear();
        int noEl = d->giveNumberOfElements();
        noNewEl = newMeshNodes.size();
        this->interpolationElementNumbers.resize(noNewEl);
        for ( int i = 1; i <= noNewEl; i++ ) {
            this->interpolationElementNumbers.at(i) = noEl + i;
        }
    } else {
        OOFEM_ERROR("Invalid input field \"GenInterpElem\". Value: %d is not supported", generateInterpolationElements);
    }

    interpolationElementIndices.resize(d->giveNumberOfElements() + noNewEl);
    for ( int i = 1; i <= interpolationElementNumbers.giveSize(); i++ ) {
        interpolationElementIndices.at( interpolationElementNumbers.at(i) ) = i;
    }
}


void
Quasicontinuum :: createInterpolationElements(Domain *d)
// create new elements. Material will be defined after post initialization
{
    /*
     * // new crosssection
     * int ncrosssect  = d->giveNumberOfCrossSectionModels() + 1;
     * d->resizeCrossSectionModels(ncrosssect);
     * CrossSection *crossSection;
     * crossSection = classFactory.createCrossSection("SimpleCS", ncrosssect, d);
     * DynamicInputRecord irCS;
     * irCS.setField(1.0, _IFT_SimpleCrossSection_thick);
     * irCS.setField(1.0, _IFT_SimpleCrossSection_width);
     * crossSection->initializeFrom(&irCS);
     * d->setCrossSection(ncrosssect, crossSection);
     */

    int ninterpelem = interpolationMeshNodes.size();

    /*
     * #if DEBUG
     * if  (ninterpelem==0)  {
     * OOFEM_WARNING("No interpolation elements generated");
     * }
     * #endif
     */

    // add interpolation elements (with new CS and no mat)
    if ( interpolationMeshNodes.size() != 0 ) {
        // for 2d
        const char *elemType;
        if ( nDimensions == 2 ) {
            //elemType  = "TrPlaneStrain";
            elemType  = "TrPlaneStress2d";
            //for 3d
        } else {
            elemType  = "LTRSpace";
        }

        int nelem = d->giveNumberOfElements();
        DynamicInputRecord irEl;
        // over interpolation elements
        for ( int i = 1; i <= ninterpelem; i++ ) {
            int elemNumber = nelem + i;
            d->resizeElements(elemNumber);
            auto elem = classFactory.createElement(elemType, elemNumber, d);
            irEl.setField(interpolationMeshNodes [ i - 1 ], Element::IPK_Element_nodes.getNameCStr());
            //irEl.setField( ncrosssect, _IFT_Element_crosssect);
            //irEl.setField( nmat, _IFT_Element_mat);
            elem->initializeFrom(irEl, 1);
            elem->setGlobalNumber(elemNumber);
            d->setElement(elemNumber, std::move(elem));
        }
    }
}


void
Quasicontinuum :: addCrosssectionToInterpolationElements(Domain *d)
{
    // new crosssection
    int ncrosssect = d->giveNumberOfCrossSectionModels() + 1;
    d->resizeCrossSectionModels(ncrosssect);
    auto crossSection = classFactory.createCrossSection("SimpleCS", ncrosssect, d);
    DynamicInputRecord irCS;
    irCS.setField(1.0, _IFT_SimpleCrossSection_thick);
    irCS.setField(1.0, _IFT_SimpleCrossSection_width);
    crossSection->initializeFrom(irCS);
    d->setCrossSection(ncrosssect, std::move(crossSection));

    for ( int i = 1; i <= interpolationElementNumbers.giveSize(); i++ ) {
        int elNum = interpolationElementNumbers.at(i);
        d->giveElement(elNum)->setCrossSection(ncrosssect);
    }
}


void
Quasicontinuum :: applyApproach1(Domain *d)
// Approach 1 --- Hangingnodes (interpolation elementss with zero stiffnes)
{
    // new material with zero stiffness
    int nmat = d->giveNumberOfMaterialModels() + 1;
    d->resizeMaterials(nmat);
    auto mat = classFactory.createMaterial("IsoLE", nmat, d);
    DynamicInputRecord irMat;
    irMat.setField(1.e-20, _IFT_IsotropicLinearElasticMaterial_e);
    irMat.setField(0.0, _IFT_IsotropicLinearElasticMaterial_n);
    irMat.setField(0.0, _IFT_IsotropicLinearElasticMaterial_talpha);
    irMat.setField(0.0, _IFT_Material_density);
    mat->initializeFrom(irMat);
    d->setMaterial(nmat, std::move(mat));


    // add material and CS to interpolation elements
    int ninterpelem = interpolationElementNumbers.giveSize();
    int matNum = nmat;
    // over interpolation elements
    for ( int i = 1; i <= ninterpelem; i++ ) {
        int elNum = interpolationElementNumbers.at(i);
        d->giveElement(elNum)->setMaterial(matNum);
    }
}


void
Quasicontinuum :: applyApproach2(Domain *d, int homMtrxType, double volumeOfInterpolationMesh)
// Approach 2 --- Global homogenization
{
    // decision which links will solved explicitly/homogenized
    elemList.resize( d->giveNumberOfElements() );
    elemList.zero();
    nodeList.resize( d->giveNumberOfDofManagers() );
    nodeList.zero();

    // ovel all elements // TO DO improve by using list of links only
    for ( int i = 1; i <= d->giveNumberOfElements(); i++ ) {
        if ( d->giveElement(i)->giveNumberOfNodes() > 2 ) {
            elemList.at(i) = 1; // all interpolation elements will be active
            continue; // skip interpolation elements
        }

        Element *e = d->giveElement(i);
        int n1 = e->giveNode(1)->giveQcNodeType();
        int n2 = e->giveNode(2)->giveQcNodeType();
        // repnode--repnode - save link to elemList
        if ( ( n1 == 1 ) &&  ( n2 == 1 ) ) {
            elemList.at(i) = 1;
            nodeList.at( e->giveNode(1)->giveNumber() ) = 1;
            nodeList.at( e->giveNode(2)->giveNumber() ) = 1;

            // hanging--hanging
        } else if ( ( n1 == 2 ) && ( n2 == 2 ) ) {
            // nothing to store, link will be homogenized-removed
            // eventualy, link with end/ends outside intero element can be solved explicitly

            // repnode--hanging
        } else if ( ( n1 == 1 ) &&  ( n2 == 2 ) ) {
            int mn = e->giveNode(1)->giveNumber();
            int hn = e->giveNode(2)->giveNumber();
            qcNode *qn = static_cast< qcNode * >( e->giveNode(2) );
            int masterElem = qn->giveMasterElementNumber();
            IntArray nodesOfMasterElem = d->giveElement(masterElem)->giveDofManArray();
            // is masterNode one of nodes of master element?
            if ( !nodesOfMasterElem.contains(mn) ) { // master node is not vertex of hn master elements
                elemList.at(i) = 1;
                nodeList.at(mn) = 1;
                nodeList.at(hn) = 1;
            } else {
                // nothing to store, link will be homogenized-removed
            }

            // hanging--repnode
        } else if ( ( n1 == 2 ) &&  ( n2 == 1 ) ) {
            int mn = e->giveNode(2)->giveNumber();
            int hn = e->giveNode(1)->giveNumber();
            qcNode *n = static_cast< qcNode * >( e->giveNode(1) );
            int masterElem = n->giveMasterElementNumber();
            IntArray nodesOfMasterElem = d->giveElement(masterElem)->giveDofManArray();
            // is masterNode one of nodes of master element?
            if ( !nodesOfMasterElem.contains(mn) ) { // master node is not vertex of hn master elements
                elemList.at(i) = 1;
                nodeList.at(mn) = 1;
                nodeList.at(hn) = 1;
            } else {
                // nothing to store, link will be homogenized-removed
            }


            // no QcNode
        } else {
            OOFEM_WARNING("Element %d with non-qcNode can not be homogenized", i);
            elemList.at(i) = 1;
            nodeList.at( e->giveNode(1)->giveNumber() ) = 1;
            nodeList.at( e->giveNode(2)->giveNumber() ) = 1;
        }
    } // ovel all links


    /* // qcRenumbering ir realized instead of this...
     * // deactivate homogenized nodes
     * for (int i=1; i<=d->giveNumberOfDofManagers(); i++) {
     * if ( !nodeList(i) ) {
     * /// d->giveNode(i)->deactivateYourself();
     * }
     * }
     *
     * // deactivate homogenized elements
     * for (int i=1; i<=d->giveNumberOfElements(); i++) {
     * if ( !elemList(i) ) {
     * /// d->giveElement(i)->deactivateYourself();
     * }
     * }
     */

    // set up lists of active nodes and elements
    QClinearStatic *em = dynamic_cast< QClinearStatic * >( d->giveEngngModel() );
    if ( em ) {
        em->setActivatedNodeList(nodeList, d);
        em->setActivatedElementList(elemList);
    } else {
        OOFEM_ERROR("Quasicontinuum can be applied only in QClinearStatic engngm");
    }


    // global homogenization
    double S0;
    FloatMatrix Diso(6, 6);
    createGlobalStiffnesMatrix(Diso, S0, d, homMtrxType, volumeOfInterpolationMesh);

    double homogenizedE, homogenizedNu;
    homogenizationOfStiffMatrix(homogenizedE, homogenizedNu, Diso);


    // new material with homogenized stiffness
    int nmat = d->giveNumberOfMaterialModels() + 1;
    d->resizeMaterials(nmat);
    std::unique_ptr<Material> mat;
    DynamicInputRecord irMat;
    // isotropic
    if ( homMtrxType == 1 ) {
        mat = classFactory.createMaterial("IsoLE", nmat, d);
        irMat.setField(homogenizedE, _IFT_IsotropicLinearElasticMaterial_e);
        irMat.setField(homogenizedNu, _IFT_IsotropicLinearElasticMaterial_n);
        irMat.setField(0.0, _IFT_IsotropicLinearElasticMaterial_talpha);
        irMat.setField(0.0, _IFT_Material_density);
        // anisotropic
    } else if ( homMtrxType == 2 ) {
        //OOFEM_ERROR("anisotropic homog. is not inmplemented yet");
        mat = classFactory.createMaterial("AnisoLE", nmat, d);
        FloatArray stiff;
        //FloatArray alpha(6);
        //alpha.zero();
        if ( nDimensions == 2 ) {
            stiff = VecX({
                Diso.at(1, 1), Diso.at(1, 2), 0, 0, 0, Diso.at(1, 6), \
                Diso.at(2, 2), 0, 0, 0, Diso.at(2, 6), 33, 0, 0, 0, 44, \
                0, 0, 55, 0, Diso.at(6, 6)
            });
        } else if ( nDimensions == 3 ) {
            stiff = VecX({
                Diso.at(1, 1), Diso.at(1, 2), Diso.at(1, 3), Diso.at(1, 4), Diso.at(1, 5), Diso.at(1, 6), \
                Diso.at(2, 2), Diso.at(2, 3), Diso.at(2, 4), Diso.at(2, 5), Diso.at(2, 6), \
                Diso.at(3, 3), Diso.at(3, 4), Diso.at(3, 5), Diso.at(3, 6), \
                Diso.at(4, 4), Diso.at(4, 5), Diso.at(4, 6), \
                Diso.at(5, 5), Diso.at(5, 6), Diso.at(6, 6)
            });
        } else {
            OOFEM_ERROR("Invalid number of dimensions. Only 2d and 3d domains are supported in QC simulation. \n");
        }
        irMat.setField(stiff, _IFT_AnisotropicLinearElasticMaterial_stiff);
        //irMat.setField(alpha, _IFT_AnisotropicLinearElasticMaterial_talpha);
        irMat.setField(0.0, _IFT_Material_density);
    } else {
        OOFEM_ERROR("Invalid homMtrxType");
    }
    mat->initializeFrom(irMat);
    d->setMaterial(nmat, std::move(mat));


    // add material to interpolation elements
    int ninterpelem = interpolationElementNumbers.giveSize();
    int matNum = nmat;
    // over interpolation elements
    for ( int i = 1; i <= ninterpelem; i++ ) {
        int elNum = interpolationElementNumbers.at(i);
        d->giveElement(elNum)->setMaterial(matNum);
    }
}


void
Quasicontinuum :: applyApproach3(Domain *d, int homMtrxType)
// Approach 3 --- Local homogenization (for each element individually)
{
    initializeConnectivityTableForInterpolationElements(d);

    // decision which links will solved explicitly/homogenized
    elemList.resize( d->giveNumberOfElements() );
    elemList.zero();
    nodeList.resize( d->giveNumberOfDofManagers() );
    nodeList.zero();

    bool stiffnesAssigned = false;
    FloatMatrix stfTensorOf1Link(9, 9);
    double s0Of1Link;

    int noIntEl = interpolationElementNumbers.giveSize();
    std :: vector< FloatMatrix > individualStiffnessMatrices;
    std :: vector< FloatMatrix > individualStiffnessTensors;
    individualStiffnessMatrices.resize(noIntEl, FloatMatrix(6, 6));
    individualStiffnessTensors.resize(noIntEl, FloatMatrix(9, 9));

    FloatArray individualS0;
    individualS0.resize(noIntEl);
    individualS0.zero();


    // ovel all elements // TO DO improve by using list of links only
    for ( int i = 1; i <= d->giveNumberOfElements(); i++ ) {
        if ( d->giveElement(i)->giveNumberOfNodes() > 2 ) {
            elemList.at(i) = 1; // all interpolation elements will be active
            continue; // skip interpolation elements
        }

        Element *e = d->giveElement(i);
        int n1 = e->giveNode(1)->giveQcNodeType();
        int n2 = e->giveNode(2)->giveQcNodeType();
        // repnode--repnode - save link to elemList
        if ( ( n1 == 1 ) &&  ( n2 == 1 ) ) {
            elemList.at(i) = 1;
            nodeList.at( e->giveNode(1)->giveNumber() ) = 1;
            nodeList.at( e->giveNode(2)->giveNumber() ) = 1;

            // hanging--hanging
        } else if ( ( n1 == 2 ) && ( n2 == 2 ) ) {
            // one or both ends ale not located outside interpolation elements
            // link will be solved explicitly
            qcNode *qn1 = static_cast< qcNode * >( e->giveNode(1) );
            qcNode *qn2 = static_cast< qcNode * >( e->giveNode(2) );
            int mn = e->giveNode(1)->giveNumber();
            int hn = e->giveNode(2)->giveNumber();
            if ( qn1->giveMasterElementNumber() < 0 || qn2->giveMasterElementNumber() < 0  ) {
                elemList.at(i) = 1;
                nodeList.at(mn) = 1;
                nodeList.at(hn) = 1;
            }
            // both end are inside interp. elements
            // link will be homogenized
            else {
                stiffnesAssigned = false;
                stiffnesAssigned = stiffnessAssignment(individualStiffnessTensors, individualS0, d, e, qn1, qn2);
                //	 stiffnesAssigned = stiffnessAssignment( individualStiffnessTensors, individualS0 );
                if ( !stiffnesAssigned ) {
                    // stiffness of this link was no successfully asigned to interpolation elements - link will be solved explicitly
                    elemList.at(i) = 1;
                    nodeList.at(mn) = 1;
                    nodeList.at(hn) = 1;
                }
            }

            // repnode--hanging
        } else if ( ( n1 == 1 ) &&  ( n2 == 2 ) ) {
            int mn = e->giveNode(1)->giveNumber();
            int hn = e->giveNode(2)->giveNumber();
            qcNode *qn = static_cast< qcNode * >( e->giveNode(2) );
            int masterElem = qn->giveMasterElementNumber();
            IntArray nodesOfMasterElem = d->giveElement(masterElem)->giveDofManArray();
            // is masterNode one of nodes of master element?
            if ( !nodesOfMasterElem.contains(mn) ) { // master node is not vertex of hn master elements
                elemList.at(i) = 1;
                nodeList.at(mn) = 1;
                nodeList.at(hn) = 1;
            } else {
                // stiffness of link is assigned to masterElement
                computeStiffnessTensorOf1Link(stfTensorOf1Link, s0Of1Link, e, d);
                individualStiffnessTensors [ interpolationElementIndices.at(masterElem) - 1 ].add(stfTensorOf1Link);
                individualS0 [ interpolationElementIndices.at(masterElem) - 1 ] += s0Of1Link;
            }

            // hanging--repnode
        } else if ( ( n1 == 2 ) &&  ( n2 == 1 ) ) {
            int mn = e->giveNode(2)->giveNumber();
            int hn = e->giveNode(1)->giveNumber();
            qcNode *n = static_cast< qcNode * >( e->giveNode(1) );
            int masterElem = n->giveMasterElementNumber();
            IntArray nodesOfMasterElem = d->giveElement(masterElem)->giveDofManArray();
            // is masterNode one of nodes of master element?
            if ( !nodesOfMasterElem.contains(mn) ) { // master node is not vertex of hn master elements
                elemList.at(i) = 1;
                nodeList.at(hn) = 1;
                nodeList.at(mn) = 1;
            } else {
                // stiffness of link is assigned to masterElement
                computeStiffnessTensorOf1Link(stfTensorOf1Link, s0Of1Link, e, d);
                individualStiffnessTensors [ interpolationElementIndices.at(masterElem) - 1 ].add(stfTensorOf1Link);
                individualS0 [ interpolationElementIndices.at(masterElem) - 1 ] += s0Of1Link;
            }


            // no QcNode
        } else {
            OOFEM_WARNING("Element %d with non-qcNode can not be homogenized", i);
            elemList.at(i) = 1;
            nodeList.at( e->giveNode(1)->giveNumber() ) = 1;
            nodeList.at( e->giveNode(2)->giveNumber() ) = 1;
        }
    } // ovel all links


    // convert stiffness tensors to stiff matrices and divide by volume of element
    double volumeofEl, th;
    for ( int i = 1; i <= noIntEl; i++ ) {
        transformStiffnessTensorToMatrix(individualStiffnessMatrices [ i - 1 ], individualStiffnessTensors [ i - 1 ]);
        volumeofEl = d->giveElement( interpolationElementNumbers.at(i) )->computeVolumeAreaOrLength();
        if ( d->giveNumberOfSpatialDimensions() == 2 ) {
            th = d->giveElement( interpolationElementNumbers.at(i) )->giveCrossSection()->give(CS_Thickness, NULL);
            volumeofEl = volumeofEl / th;
        }
        individualStiffnessMatrices [ i - 1 ].times(1 / volumeofEl);
    }
    // ...


    // deactivate homogenized nodes
    for ( int i = 1; i <= d->giveNumberOfDofManagers(); i++ ) {
        if ( !nodeList.at(i) ) {
            // d->giveNode(i)->deactivateYourself();
        }
    }

    // deactivate homogenized elements
    for ( int i = 1; i <= d->giveNumberOfElements(); i++ ) {
        if ( !elemList.at(i) ) {
            // d->giveElement(i)->deactivateYourself();
        }
    }

    // set up lists of active nodes and elements
    QClinearStatic *em = dynamic_cast< QClinearStatic * >( d->giveEngngModel() );
    if ( em ) {
        em->setActivatedNodeList(nodeList, d);
        em->setActivatedElementList(elemList);
    } else {
        OOFEM_ERROR("Quasicontinuum can be applied only in QClinearStatic engngm");
    }


    // new materials with homogenized stiffness
    int originalNumberOfMaterials = d->giveNumberOfMaterialModels();
    d->resizeMaterials(originalNumberOfMaterials + noIntEl);
    int nmat = originalNumberOfMaterials;
    DynamicInputRecord irMat;

    // isotropic
    if ( homMtrxType == 1 ) {
        double homogenizedE, homogenizedNu;
        for ( int i = 1; i <= noIntEl; i++ ) {
            nmat++;
            homogenizationOfStiffMatrix(homogenizedE, homogenizedNu, individualStiffnessMatrices [ i - 1 ]);
            auto mat = classFactory.createMaterial("IsoLE", nmat, d);
            irMat.setField(homogenizedE, _IFT_IsotropicLinearElasticMaterial_e);
            irMat.setField(homogenizedNu, _IFT_IsotropicLinearElasticMaterial_n);
            irMat.setField(0.0, _IFT_IsotropicLinearElasticMaterial_talpha);
            irMat.setField(0.0, _IFT_Material_density);

            mat->initializeFrom(irMat);
            d->setMaterial(nmat, std::move(mat));
        }

        // anisotropic
    } else if ( homMtrxType == 2 ) {
        for ( int i = 1; i <= noIntEl; i++ ) {
            FloatMatrix &Da = individualStiffnessMatrices [ i - 1 ];
            nmat++;
            auto mat = classFactory.createMaterial("AnisoLE", nmat, d);
            FloatArray stiff;
            // FloatArray alpha(6);
            if ( nDimensions == 2 ) {
                stiff = VecX({
                    Da.at(1, 1), Da.at(1, 2), 0, 0, 0, Da.at(1, 6), \
                    Da.at(2, 2), 0, 0, 0, Da.at(2, 6), 33, 0, 0, 0, 44, \
                    0, 0, 55, 0, Da.at(6, 6)
                });
            } else if ( nDimensions == 3 ) {
                stiff = VecX({
                    Da.at(1, 1), Da.at(1, 2), Da.at(1, 3), Da.at(1, 4), Da.at(1, 5), Da.at(1, 6), \
                    Da.at(2, 2), Da.at(2, 3), Da.at(2, 4), Da.at(2, 5), Da.at(2, 6), \
                    Da.at(3, 3), Da.at(3, 4), Da.at(3, 5), Da.at(3, 6), \
                    Da.at(4, 4), Da.at(4, 5), Da.at(4, 6), \
                    Da.at(5, 5), Da.at(5, 6), Da.at(6, 6)
                });
            } else {
                OOFEM_ERROR("Invalid number of dimensions. Only 2d and 3d domains are supported in QC simulation. \n");
            }
            irMat.setField(stiff, _IFT_AnisotropicLinearElasticMaterial_stiff);
            // irMat.setField(alpha, _IFT_AnisotropicLinearElasticMaterial_talpha);
            irMat.setField(0.0, _IFT_Material_density);

            mat->initializeFrom(irMat);
            d->setMaterial(nmat, std::move(mat));
        }
    } else {
        OOFEM_ERROR("Invalid homMtrxType");
    }


    // add material to interpolation elements
    int elNum;
    int matNum = originalNumberOfMaterials;
    // over interpolation elements
    for ( int i = 1; i <= noIntEl; i++ ) {
        matNum++;
        elNum = interpolationElementNumbers.at(i);
        d->giveElement(elNum)->setMaterial(matNum);
    }
}


void
Quasicontinuum :: homogenizationOfStiffMatrix(double &homogenizedE, double &homogenizedNu, const FloatMatrix &Diso)
// identification of E and Nu from isotropic stiffness matrix
{
    // for 2d (planme stress only)
    if ( this->nDimensions == 2 ) {
        double a = Diso.at(1, 1);
        double b = Diso.at(1, 2);
        double d = Diso.at(2, 2);
        double f = Diso.at(6, 6); // xy stiff is in last column
        // Norma A -  using eigen vectors - not implemented in 3D
        if ( 33 * a + 2 * b + 33 * d + 4 * f == 0 ) { // zero matrix = empty interpolation element -> zero stiffness
            homogenizedE  = 1.0e-20;
            homogenizedNu = 0;
        } else {
            homogenizedE  = 4 * ( a + 2 * b + d ) * ( 4 * a - 8 * b + 4 * d + f ) / ( 33 * a + 2 * b + 33 * d + 4 * f );
            homogenizedNu = ( a + 66 * b + d - 4 * f ) / ( 33 * a + 2 * b + 33 * d + 4 * f );
        }
        // Norma B -  sum of squares of differences (with weight functions)
        /*
         * if (3*a+12*b+3*d==0) {// zero matrix = empty interpolation element -> zero stiffness
         *  homogenizedE  = 1.0e-20;
         *  homogenizedNu = 0;
         * } else {
         *  homogenizedE  = (a+2*b+d)*(a+14*b+d)/2/(3*a+12*b+3*d);
         *  homogenizedNu = 2*(a-b+d)/(3*a+12*b+3*d);
         * }
         */

        //for 3d
    } else if ( this->nDimensions == 3 ) {
        // Norma A -  using eigen vectors
        //- not implemented in 3D

        // Norma B -  sum of squares of differences (with weight functions)
        double a = Diso.at(1, 1) + Diso.at(2, 2) + Diso.at(3, 3);
        double b = Diso.at(4, 4) + Diso.at(5, 5) + Diso.at(6, 6);

        if ( 5 * a + 13 * b == 0 ) { // zero matrix = empty interpolation element -> zero stiffness
            homogenizedE  = 1.0e-20;
            homogenizedNu = 0;
        } else {
            homogenizedE  = ( a + 2 * b ) * ( a + 11 * b ) / ( 15 * a + 39 * b );
            homogenizedNu = ( 2 * a + b ) / ( 5 * a + 13 * b );
        }
    } else {
        OOFEM_ERROR("Invalid number of dimensions. Only 2d and 3d domains are supported in QC simulation. \n");
    }
}


void
Quasicontinuum :: computeStiffnessTensorOf1Link(FloatMatrix &D1, double &S0, Element *e, Domain *d)
{
    // stiffness tensor 3^4 of 1 link represented by 9x9 mtrx
    D1.resize(9, 9);
    D1.zero();

    S0 = 0.;

    QCMaterialExtensionInterface *qcmei =  static_cast< QCMaterialExtensionInterface * >( e->giveMaterial()->giveInterface(QCMaterialExtensionInterfaceType) );
    if ( !qcmei ) {
        OOFEM_ERROR("Material doesn't implement the required QC interface!");
    }
    double Etruss = qcmei->giveQcElasticParamneter();
    double Ry     = qcmei->giveQcPlasticParamneter();

    double area = 0.;
    SimpleCrossSection *cs = dynamic_cast< SimpleCrossSection * >( e->giveCrossSection() );
    if ( cs ) {
        area = cs->give(CS_Area, NULL);
    } else {
        OOFEM_ERROR( "Invalid CrossSection of link %d. simpleCS is only supported CS of links in QC simulation. \n", e->giveGlobalNumber() );
    }

    double x1 = e->giveDofManager(1)->giveCoordinate(1);
    double y1 = e->giveDofManager(1)->giveCoordinate(2);
    double z1 = e->giveDofManager(1)->giveCoordinate(3);
    double x2 = e->giveDofManager(2)->giveCoordinate(1);
    double y2 = e->giveDofManager(2)->giveCoordinate(2);
    double z2 = e->giveDofManager(2)->giveCoordinate(3);
    double TrussLength = sqrt( ( x2 - x1 ) * ( x2 - x1 ) + ( y2 - y1 ) * ( y2 - y1 ) + ( z2 - z1 ) * ( z2 - z1 ) );
    double n [ 3 ];
    n [ 0 ] = ( x2 - x1 ) / TrussLength;
    n [ 1 ] = ( y2 - y1 ) / TrussLength;
    n [ 2 ] = ( z2 - z1 ) / TrussLength;
    // create stiffness tensor
    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            for ( int k = 1; k <= 3; k++ ) {
                for ( int l = 1; l <= 3; l++ ) {
                    D1.at(3 * ( i - 1 ) + k, 3 * ( j - 1 ) + l) = TrussLength * Etruss * n [ i - 1 ] * n [ j - 1 ] * n [ k - 1 ] * n [ l - 1 ];
                }
            }
        }
    }

    // homogenized yeld limit
    if ( !std :: isinf(Ry) ) {
        S0 += TrussLength * area * Ry / 3;
    }
}


void
Quasicontinuum :: createGlobalStiffnesMatrix(FloatMatrix &D, double &S0, Domain *d, int homMtrxType,  double volumeOfInterpolationMesh)
{
    int noe = d->giveNumberOfElements();

    S0 = 0;
    //double StiffessTensor[81]; // stiffness tensor 3^4 represented by 9x9 mtrx
    FloatMatrix StiffnessTensor(9, 9);

    //double D1[81]; // stiffness tensor 3^4 of 1 link represented by 9x9 mtrx
    FloatMatrix D1(9, 9);

    double n [ 3 ];
    // over all elements
    for ( int t = 1; t <= noe; t++ ) {
        // ??km?? TO DO: use list of interpolation elements instead of skipping...
        if ( d->giveElement(t)->giveNumberOfDofManagers() > 2 ) {
            continue; // skip interpolation elements
        }

        // TO DO: function "computeStiffnessMatrixOf1Link" can be used here
        QCMaterialExtensionInterface *qcmei =  static_cast< QCMaterialExtensionInterface * >( d->giveElement(t)->giveMaterial()->giveInterface(QCMaterialExtensionInterfaceType) );
        if ( !qcmei ) {
            OOFEM_ERROR("Material doesn't implement the required QC interface!");
        }
        double Etruss = qcmei->giveQcElasticParamneter();
        double Ry     = qcmei->giveQcPlasticParamneter();

        double area = 0.;
        SimpleCrossSection *cs = dynamic_cast< SimpleCrossSection * >( d->giveElement(t)->giveCrossSection() );
        if ( cs ) {
            area = cs->give(CS_Area, NULL);
        } else {
            OOFEM_ERROR( "Invalid CrossSection of link %d. simpleCS is only supported CS of links in QC simulation. \n", d->giveElement(t)->giveGlobalNumber() );
        }

        double x1 = d->giveElement(t)->giveDofManager(1)->giveCoordinate(1);
        double y1 = d->giveElement(t)->giveDofManager(1)->giveCoordinate(2);
        double z1 = d->giveElement(t)->giveDofManager(1)->giveCoordinate(3);
        double x2 = d->giveElement(t)->giveDofManager(2)->giveCoordinate(1);
        double y2 = d->giveElement(t)->giveDofManager(2)->giveCoordinate(2);
        double z2 = d->giveElement(t)->giveDofManager(2)->giveCoordinate(3);
        double TrussLength = sqrt( ( x2 - x1 ) * ( x2 - x1 ) + ( y2 - y1 ) * ( y2 - y1 ) + ( z2 - z1 ) * ( z2 - z1 ) );
        n [ 0 ] = ( x2 - x1 ) / TrussLength;
        n [ 1 ] = ( y2 - y1 ) / TrussLength;
        n [ 2 ] = ( z2 - z1 ) / TrussLength;
        //  stiffness tensor
        for ( int i = 1; i <= 3; i++ ) {
            for ( int j = 1; j <= 3; j++ ) {
                for ( int k = 1; k <= 3; k++ ) {
                    for ( int l = 1; l <= 3; l++ ) {
                        D1.at(3 * ( i - 1 ) + k, 3 * ( j - 1 ) + l) = TrussLength * Etruss * n [ i - 1 ] * n [ j - 1 ] * n [ k - 1 ] * n [ l - 1 ];
                    }
                }
            }
        }

        StiffnessTensor.add(D1); // add stfMtrX of 1 link to global stfMtrx

        if ( !std :: isinf(Ry) ) {
            S0 += TrussLength * area * Ry / 3;
        }
    } // for t=1:note

    // divide by volume of element
    StiffnessTensor.times(1 / volumeOfInterpolationMesh);
    S0 = S0 / volumeOfInterpolationMesh;

    transformStiffnessTensorToMatrix(D, StiffnessTensor);
}


bool
Quasicontinuum :: stiffnessAssignment(std :: vector< FloatMatrix > &individualStiffnessTensors, FloatArray &individualS0, Domain *d, Element *e, qcNode *qn1, qcNode *qn2)
{
    FloatMatrix stfTensorOf1Link(9, 9);
    double s0Of1Link;

    // numbers of elements containing ending nodes
    int end1 = qn1->giveMasterElementNumber();
    int end2 = qn2->giveMasterElementNumber();

    /*
     * if ( (end1==0) && (end2==0) ) { // both ends are located outside
     * // this case is solved explicitly before StiffnessAssignment
     * } elseif ( (end1==0) || (end2==0) ) { // both ends are located outside
     * // this case is solved explicitly before StiffnessAssignment
     * }
     */

    if ( end1 == end2 ) {  // both ends are in the same el.
        // stiffness of link is assigned to this el.
        computeStiffnessTensorOf1Link(stfTensorOf1Link, s0Of1Link, e, d);
        int nel = end1; // number of interp element
        int indx = interpolationElementIndices.at(nel); // number in list of interp elem
        individualStiffnessTensors [ indx - 1 ].add(stfTensorOf1Link);
        individualS0 [ indx - 1 ] += s0Of1Link;

        return true;
    } else {   // ends are located in different elements -> all intersected el needs to be found
        IntArray intersected; // to store numbers of intersected elem
        std::vector<double> lengths_; // to store lengths ofintersections

        computeIntersectionsOfLinkWithInterpElements(intersected, lengths_, d, e, qn1, qn2);
        FloatArray lengths=FloatArray::fromVector(lengths_);

        // check: is there any intersectted elements
        int noIntersected = lengths.giveSize();
        if ( noIntersected == 0 ) {
            return false;
        }

        // check: is stiff. assigned to all parts of the link?
        // TO DO: the same check is done in "computeIntersectionsOfLinkWith...".
        //        Use only return value here
        double totalLength = distance(qn1->giveCoordinates(), qn2->giveCoordinates());
        lengths.times(1 / totalLength); // it is better to use relative lengths
        double sumLength = lengths.sum();
        if ( (sumLength < 0.9999) || (sumLength > 1.0001) ) {
            return false;
        }

        // if check OK ( sumLength == 1 )
        computeStiffnessTensorOf1Link(stfTensorOf1Link, s0Of1Link, e, d);

        for ( int i = 1; i <= noIntersected; i++ ) { // add stiff to all intersected elem
            int nel = intersected [ i - 1 ]; // number of element
            int indx = interpolationElementIndices.at(nel); // number in list of interp elem
            individualStiffnessTensors [ indx - 1 ].add(lengths.at(i), stfTensorOf1Link);
            individualS0 [ indx - 1 ] += lengths.at(i) * s0Of1Link;
        }
        return true;
    }  // ends are located in different elements

    return false;
}


void
Quasicontinuum :: computeIntersectionsOfLinkWithInterpElements(IntArray &intersected, std::vector<double> &lengths, Domain *d, Element *e, qcNode *qn1, qcNode *qn2)
{
    int dim = d->giveNumberOfSpatialDimensions();
    if ( dim == 2 ) {
        computeIntersectionsOfLinkWith2DTringleElements(intersected, lengths, d, e, qn1, qn2);
    } else if ( dim == 3 ) {
        computeIntersectionsOfLinkWith3DTetrahedraElements(intersected, lengths, d, e, qn1, qn2);
    } else {
        OOFEM_ERROR("unsupported number of dimensions");
    }
}


bool
Quasicontinuum :: computeIntersectionsOfLinkWith2DTringleElements(IntArray &intersected, std::vector<double> &lengths, Domain *d, Element *e, qcNode *qn1, qcNode *qn2)
{
    intersected.clear();
    lengths.clear();

    // numbers of elements containing ending nodes
    int end1 = qn1->giveMasterElementNumber();
    int end2 = qn2->giveMasterElementNumber();

    if ( end2 < end1 ) { // swap?
        std::swap(end1, end2);
        std::swap(qn1, qn2);
    }

    // coordinates of ending nodes
    const auto &X1 = qn1->giveCoordinates();
    const auto &X2 = qn2->giveCoordinates();

    double TotalLength = distance(X2, X1);

    int numberOfIntersected = 0; // number of intersected elements

    FloatArray vector(2);
    std::vector<FloatArray> intersectCoords;
    double iLen, sumLength = 0.0;

    // ftiffness assignment for starting element (end1)
    int iel = end1;
    FloatArray iLens;

    const auto &A = d->giveElement(iel)->giveNode(1)->giveCoordinates();
    const auto &B = d->giveElement(iel)->giveNode(2)->giveCoordinates();
    const auto &C = d->giveElement(iel)->giveNode(3)->giveCoordinates();
    numberOfIntersected = intersectionTestSegmentTriangle2D(intersectCoords, A, B, C, X1, X2);

    switch ( numberOfIntersected ) {
    case 0:     // no intersection with this element
        lengths.push_back(0);
        intersected.followedBy(iel);     // cislo elementu
        break;
    case 1:     // link starts in this element
        iLen = distance(X1, intersectCoords[0]);
        if ( iLen / TotalLength >= 1.e-6 ) { // intersection is not negligible
            lengths.push_back(iLen);
            intersected.followedBy(iel);
        } else {     // intersection is negligible
            lengths.push_back(0);
            intersected.followedBy(iel);     // cislo elementu
        }
        break;
    default:     // 2 or more intersections

        //lengths of all possible cobmination of intersection points
        iLens.resize(numberOfIntersected);
        for ( int nn = 1; nn <= numberOfIntersected; nn++ ) {
            iLens.at(nn) = distance(X1, intersectCoords[nn - 1]);
        }
        iLen = iLens.at( iLens.giveIndexMaxElem() );     // real length is maximal one
        if ( iLen / TotalLength >= 1.e-6 ) { // intersection is not negligible
            lengths.push_back(iLen);
            intersected.followedBy(iel);
        } else {     // intersection is negligible
            lengths.push_back(0);
            intersected.followedBy(iel);     // cislo elementu
        }
        break;
    }


    // ftiffness assignment for the rest of the link
    IntArray neighboursList, alreadyTestedElemList;
    alreadyTestedElemList.clear();

    int nome = interpolationElementNumbers.giveSize();

    //bool secondEndReached = false;
    int nextElement = end1;

    for ( int ii = 1; ii <= nome; ii++ ) { // searching for nex intersected element
        int nel = nextElement; // previous element number
        nextElement = 0; // next element number (to be found)

        // list of connected elements
        neighboursList = connectivityTable [ interpolationElementIndices.at(nel) - 1 ];

        // delete already intersected elements
        for ( int k = 1; k <= intersected.giveSize(); k++ ) {
            int index = neighboursList.findFirstIndexOf( intersected.at(k) );
            if ( index != 0 ) {
                neighboursList.erase(index);
            }
        }

        // delete already tested elements
        for ( int k = 1; k <= alreadyTestedElemList.giveSize(); k++ ) {
            int index = neighboursList.findFirstIndexOf( alreadyTestedElemList.at(k) );
            if ( index != 0 ) {
                neighboursList.erase(index);
            }
        }


        if ( neighboursList.giveSize() == 0 ) { // no elements to continue
            //secondEndReached = false;
            break; // ends searching for nex intersected element
        }

        for ( int k = 1; k <= neighboursList.giveSize(); k++ ) { // testing of neighbouring elements
            iel = neighboursList.at(k);

            const auto &A = d->giveElement(iel)->giveNode(1)->giveCoordinates();
            const auto &B = d->giveElement(iel)->giveNode(2)->giveCoordinates();
            const auto &C = d->giveElement(iel)->giveNode(3)->giveCoordinates();
            numberOfIntersected = intersectionTestSegmentTriangle2D(intersectCoords, A, B, C, X1, X2);

            bool breakFor = false; // go to next element?
            switch ( numberOfIntersected ) {
            case 0: // no intersection with (iel) element -> go to next element
                alreadyTestedElemList.followedBy(iel);
                continue;
            case 1: // link ends here or there in only fake intersection
                if ( iel == end2 ) { // end of the link in in this element
                    iLen = distance(X2, intersectCoords[0]);
                    if ( iLen / TotalLength >= 1.e-6 ) { // intersection is not negligible
                        lengths.push_back(iLen);
                        intersected.followedBy(iel);
                        // no break here, some elements still needs to be tested
                    } else { // intersection is negligible
                        lengths.push_back(0);
                        intersected.followedBy(iel); // cislo elementu
                    }
                } else { // fake intersection (only touch with edge)
                    lengths.push_back(0);
                    intersected.followedBy(iel);
                }
                break;
            default: // 2 or more intersections
                if ( iel == end2 ) { // end of the link is inthis element
                    //lengths of all possible cobmination of intersection points
                    iLens.resize(numberOfIntersected);
                    for ( int nn = 1; nn <= numberOfIntersected; nn++ ) {
                        iLens.at(nn) = distance(X2, intersectCoords[nn - 1]);
                    }
                    iLen = iLens.at( iLens.giveIndexMaxElem() ); // real length is maximal
                    if ( iLen / TotalLength >= 1.e-6 ) { // intersection is not negligible
                        lengths.push_back(iLen);
                        intersected.followedBy(iel);
                        nextElement = iel;
                        breakFor = true;
                        break; // go to next element
                    } else { // intersection is negligible
                        lengths.push_back(0);
                        intersected.followedBy(iel);
                    }
                } else { // end of the link is NOT inthis elem
                    iLens.resize(3); // combinations: 12 13 23
                    iLens.zero();
                    int m = 0;
                    for ( int nn = 1; nn < numberOfIntersected; nn++ ) {
                        for ( int ii = nn + 1; ii <= numberOfIntersected; ii++ ) {
                            m += 1;
                            iLens.at(m) = distance(intersectCoords[ii - 1], intersectCoords[nn - 1]);
                        }
                    }
                    iLen = iLens.at( iLens.giveIndexMaxElem() ); // real length is maximum
                    if ( iLen / TotalLength >= 1.e-6 ) { // intersection is not negligible
                        lengths.push_back(iLen);
                        intersected.followedBy(iel);
                        nextElement = iel;
                        breakFor = true;
                        break; // go to next element
                    } else { // intersection is negligible
                        lengths.push_back(0);
                        intersected.followedBy(iel);
                    }
                } // end of the link is inthis element
            } // switch

            if ( breakFor ) {
                break;
            }                        // go to next interpolation element
        } // over neigbouring elements


        // check: all length of link is assembles
        sumLength = std::accumulate(lengths.begin(), lengths.end(), 0.);
        if ( ( fabs(sumLength / TotalLength - 1) ) <= ( 0.0001 ) ) {
            return true;
        }

        // check: there is no neighnouring element to continue
        if ( nextElement == 0 ) {
            return false;
        }
    } // searching for nex intersected element

    return false;
}


bool
Quasicontinuum :: computeIntersectionsOfLinkWith3DTetrahedraElements(IntArray &intersected, std::vector<double> &lengths, Domain *d, Element *e, qcNode *qn1, qcNode *qn2)
{
    intersected.clear();
    lengths.clear();

    // numbers of elements containing ending nodes
    int end1 = qn1->giveMasterElementNumber();
    int end2 = qn2->giveMasterElementNumber();

    if ( end2 < end1 ) { // swap?
        std::swap(end1, end2);
        std::swap(qn1, qn2);
    }

    // coordinates of ending nodes
    const auto &X1 = qn1->giveCoordinates();
    const auto &X2 = qn2->giveCoordinates();

    double TotalLength = distance(X2, X1);

    int numberOfIntersected = 0; // number of intersected elements

    std::vector<FloatArray> intersectCoords;
    double iLen, sumLength = 0.0;

    // ftiffness assignment for starting element (end1)

    int iel = end1;
    FloatArray iLens;

    const auto &A = d->giveElement(iel)->giveNode(1)->giveCoordinates();
    const auto &B = d->giveElement(iel)->giveNode(2)->giveCoordinates();
    const auto &C = d->giveElement(iel)->giveNode(3)->giveCoordinates();
    const auto &D = d->giveElement(iel)->giveNode(4)->giveCoordinates();
    numberOfIntersected = intersectionTestSegmentTetrahedra3D(intersectCoords, A, B, C, D, X1, X2);

    switch ( numberOfIntersected ) {
    case 0:     // no intersection with this element
        lengths.push_back(0);
        intersected.followedBy(iel);     // cislo elementu
        break;
    case 1:     // link starts in this element
        iLen = distance(X1, intersectCoords[0]);
        if ( iLen / TotalLength >= 1.e-6 ) { // intersection is not negligible
            lengths.push_back(iLen);
            intersected.followedBy(iel);
        } else {     // intersection is negligible
            lengths.push_back(0);
            intersected.followedBy(iel);
        }
        break;
    default:     // 2 or more intersections

        //lengths of all possible cobmination of intersection points
        iLens.resize(numberOfIntersected);
        for ( int nn = 1; nn <= numberOfIntersected; nn++ ) {
            iLens.at(nn) = distance(X1, intersectCoords[nn - 1]);
        }
        iLen = iLens.at( iLens.giveIndexMaxElem() );     // real length is maximal one
        if ( iLen / TotalLength >= 1.e-6 ) { // intersection is not negligible
            lengths.push_back(iLen);
            intersected.followedBy(iel);
        } else {     // intersection is negligible
            lengths.push_back(0);
            intersected.followedBy(iel);
        }
        break;
    }

    // ftiffness assignment for the rest of the link

    bool breakFor = false;
    IntArray neighboursList, alreadyTestedElemList;
    alreadyTestedElemList.clear();

    int nome = interpolationElementNumbers.giveSize();

    //bool secondEndReached = false;
    int nextElement = end1;

    for ( int ii = 1; ii <= nome; ii++ ) { // searching for nex intersected element
        int nel = nextElement; // previous element number
        nextElement = 0; // next element number (to be found)

        // list of connected elements
        neighboursList = connectivityTable [ interpolationElementIndices.at(nel) - 1 ];

        // delete already intersected elements
        for ( int k = 1; k <= intersected.giveSize(); k++ ) {
            int index = neighboursList.findFirstIndexOf( intersected.at(k) );
            if ( index != 0 ) {
                neighboursList.erase(index);
            }
        }

        // delete already tested elements
        for ( int k = 1; k <= alreadyTestedElemList.giveSize(); k++ ) {
            int index = neighboursList.findFirstIndexOf( alreadyTestedElemList.at(k) );
            if ( index != 0 ) {
                neighboursList.erase(index);
            }
        }

        if ( neighboursList.giveSize() == 0 ) { // no elements to continue
            //secondEndReached = false;
            break; // ends searching for nex intersected element
        }

        for ( int k = 1; k <= neighboursList.giveSize(); k++ ) { // testing of neighbouring elements
            int iel = neighboursList.at(k);

            const auto &A = d->giveElement(iel)->giveNode(1)->giveCoordinates();
            const auto &B = d->giveElement(iel)->giveNode(2)->giveCoordinates();
            const auto &C = d->giveElement(iel)->giveNode(3)->giveCoordinates();
            const auto &D = d->giveElement(iel)->giveNode(4)->giveCoordinates();
            int numberOfIntersected = intersectionTestSegmentTetrahedra3D(intersectCoords, A, B, C, D, X1, X2);


            breakFor = false; // go to next interp elem
            switch ( numberOfIntersected ) {
            case 0: // no intersection with (iel) element -> go to next element
                alreadyTestedElemList.followedBy(iel);
                continue; // zadna tuhost se prvku neprirazuje
            case 1: // link ends here or there in only fake intersection
                if ( iel == end2 ) { // end of the link in in this element
                    iLen = distance(X2, intersectCoords[0]);
                    if ( iLen / TotalLength >= 1.e-6 ) { // intersection is not negligible
                        lengths.push_back(iLen);
                        intersected.followedBy(iel);
                        // no break here, some elements still needs to be tested
                    } else { // intersection is negligible
                        lengths.push_back(0);
                        intersected.followedBy(iel); // cislo elementu
                    }
                } else { // fake intersection (only touch with edge)
                    lengths.push_back(0);
                    intersected.followedBy(iel);
                }
                break;
            default: // 2 or more intersections
                if ( iel == end2 ) { // end of the link is inthis element
                    //lengths of all possible cobmination of intersection points
                    iLens.resize(numberOfIntersected);
                    for ( int nn = 1; nn <= numberOfIntersected; nn++ ) {
                        iLens.at(nn) = distance(X2, intersectCoords[nn - 1]);
                    }
                    iLen = iLens.at( iLens.giveIndexMaxElem() ); // real length is maximal
                    if ( iLen / TotalLength >= 1.e-6 ) { // intersection is not negligible
                        lengths.push_back(iLen);
                        intersected.followedBy(iel);
                        nextElement = iel;
                        breakFor = true;
                        break; // go to next element
                    } else { // intersection is negligible
                        lengths.push_back(0);
                        intersected.followedBy(iel);
                    }
                } else { // end of the link is NOT inthis elem
                    iLens.resize(6); // combinations: 12 13 14 23 24 34
                    iLens.zero();
                    int m = 0;
                    for ( int nn = 1; nn <= numberOfIntersected; nn++ ) {
                        for ( int ii = nn + 1; ii <= numberOfIntersected; ii++ ) {
                            m += 1;
                            iLens.at(m) = distance(intersectCoords[ii - 1], intersectCoords[nn - 1]);
                        }
                    }
                    iLen = iLens.at( iLens.giveIndexMaxElem() ); // real length is maximum
                    if ( iLen / TotalLength >= 1.e-6 ) { // intersection is not negligible
                        lengths.push_back(iLen);
                        intersected.followedBy(iel);
                        nextElement = iel;
                        breakFor = true;
                        break; // go to next element
                    } else { // intersection is negligible
                        lengths.push_back(0);
                        intersected.followedBy(iel);
                    }
                } // end of the link is inthis element
            } // switch

            if ( breakFor ) {
                break;
            }                        // go to next interp elem
        } // over neigbouring elements


        // check: all length of link is assembles
        sumLength = std::accumulate(lengths.begin(), lengths.end(), 0.);
        if ( ( fabs(sumLength / TotalLength - 1) ) <= ( 0.0001 ) ) {
            return true;
        }

        // check: there is no neighnouring element to continue
        if ( nextElement == 0 ) {
            return false;
        }
    } // searching for nex intersected element

    return false;
}


void
Quasicontinuum :: initializeConnectivityTableForInterpolationElements(Domain *d)
{
    // for each element save elements with shared node (edge, face,,,)
    int nome = interpolationElementNumbers.giveSize();
    int noVertices1, noVertices2;
    IntArray vertices1;
    connectivityTable.clear();
    connectivityTable.resize(nome);

    for ( int j = 1; j <= nome; j++ ) { // all interpolation elements
        int elNum1 = interpolationElementNumbers.at(j);
        Element *e1 = d->giveElement(elNum1);

        // vector of vertices numbers of el j
        noVertices1 = e1->giveNumberOfNodes();
        vertices1.resize(noVertices1);
        for ( int k = 1; k <= noVertices1; k++ ) { // fill in vector
            vertices1.at(k) = e1->giveNode(k)->giveNumber();
        }
        for ( int i = 1; i <= nome; i++ ) {
            if ( i == j ) {
                continue;
            }                       // is not itself neighbour
            // vector of vertices numbers of el j
            int elNum2 = interpolationElementNumbers.at(i);
            Element *e2 = d->giveElement(elNum2);
            noVertices2 = e2->giveNumberOfNodes();
            for ( int l = 1; l <= noVertices2; l++ ) { // fill in vector
                if ( vertices1.contains( e2->giveNode(l)->giveNumber() ) ) {
                    connectivityTable [ j - 1 ].followedBy(elNum2); // add number of neighbour
                    break;
                }
            }
        }
    }
}


bool
Quasicontinuum :: intersectionTestSegmentTrianglePlucker3D(FloatArray  &intersectCoords, const FloatArray &A, const FloatArray &B, const FloatArray &C, const FloatArray &X1, const FloatArray &X2)
{
    // test if line (X2-X1) intersects tringle ABC (for 3D!)
    intersectCoords.clear();

    // directional vektor
    double L = distance(X2, X1);
    double Ur1 = ( X2.at(1) - X1.at(1) ) / L;
    double Ur2 = ( X2.at(2) - X1.at(2) ) / L;
    double Ur3 = ( X2.at(3) - X1.at(3) ) / L;
    // cross product (needed for Plucker. coords)
    double Vr1 = Ur2 * X1.at(3) - Ur3 *X1.at(2);
    double Vr2 = Ur3 * X1.at(1) - Ur1 *X1.at(3);
    double Vr3 = Ur1 * X1.at(2) - Ur2 *X1.at(1);
    // Plucker. coords of line {U,UxX1}
    // Pr = [Ur, Vr];

    //  vectors of edges: (order is important!)
    FloatArray U1(3), U2(3), U3(3);
    for ( int i = 1; i <= 3; i++ ) {
        U1.at(i) = B.at(i) - A.at(i);
        U2.at(i) = C.at(i) - B.at(i);
        U3.at(i) = A.at(i) - C.at(i);
    }

    // cross product (needed for Plucker. coords)
    FloatArray V1(3), V2(3), V3(3);
    // cross(U1,A);
    V1.at(1) = U1.at(2) * A.at(3) - U1.at(3) * A.at(2);
    V1.at(2) = U1.at(3) * A.at(1) - U1.at(1) * A.at(3);
    V1.at(3) = U1.at(1) * A.at(2) - U1.at(2) * A.at(1);
    // cross(U2,B);
    V2.at(1) = U2.at(2) * B.at(3) - U2.at(3) * B.at(2);
    V2.at(2) = U2.at(3) * B.at(1) - U2.at(1) * B.at(3);
    V2.at(3) = U2.at(1) * B.at(2) - U2.at(2) * B.at(1);
    // cross(U3,C);
    V3.at(1) = U3.at(2) * C.at(3) - U3.at(3) * C.at(2);
    V3.at(2) = U3.at(3) * C.at(1) - U3.at(1) * C.at(3);
    V3.at(3) = U3.at(1) * C.at(2) - U3.at(2) * C.at(1);


    // permuted inner products Pr(.)P1 = Ur.V1 + U1.Vr
    double S1 = Ur1 * V1.at(1) + Ur2 *V1.at(2) + Ur3 *V1.at(3) + U1.at(1) * Vr1 + U1.at(2) * Vr2 + U1.at(3) * Vr3; // dot(Ur,V1) + dot(U1,Vr);
    double S2 = Ur1 * V2.at(1) + Ur2 *V2.at(2) + Ur3 *V2.at(3) + U2.at(1) * Vr1 + U2.at(2) * Vr2 + U2.at(3) * Vr3; // dot(Ur,V2) + dot(U2,Vr);
    double S3 = Ur1 * V3.at(1) + Ur2 *V3.at(2) + Ur3 *V3.at(3) + U3.at(1) * Vr1 + U3.at(2) * Vr2 + U3.at(3) * Vr3; // dot(Ur,V3) + dot(U3,Vr);

    // tolerance ??
    double TOL = 1e-11;
    if ( fabs(S1) <= TOL ) {
        S1 = 0;
    }
    if ( fabs(S2) <= TOL ) {
        S2 = 0;
    }
    if ( fabs(S3) <= TOL ) {
        S3 = 0;
    }

    // in planar test
    if ( ( S1 == 0 ) && ( S2 == 0 ) && ( S3 == 0 ) ) {
        // segment is in the same plane as triangle
        // in this case it means no ontersection
        // for tetrahedra elements there will be another 2 intersections with remaining faces
        intersectCoords.clear();
        return false; // no intersection
    }


    // same orientation -> intersection with line (not necessarily with segment)
    if ( ( ( S1 >= 0 ) && ( S2 >= 0 ) && ( S3 >= 0 ) )  ||  ( ( S1 <= 0 ) && ( S2 <= 0 ) && ( S3 <= 0 ) ) ) {
        // calculation of intersection coords
        double S = S1 + S2 + S3;
        // from barycentric coords
        double u1 = S1 / S;
        double u2 = S2 / S;
        double u3 = S3 / S;
        // to cartesian
        // Px = [x y z] = u1*C + u2*A +u3*B; // multiple by opposite vertex
        double x = u1 * C.at(1) + u2 *A.at(1) + u3 *B.at(1);
        double y = u1 * C.at(2) + u2 *A.at(2) + u3 *B.at(2);
        double z = u1 * C.at(3) + u2 *A.at(3) + u3 *B.at(3);

        // is interection on line segment

        // intersection parameter
        double t;
        // Pk = [x y z] = X1 + t*(X2-X1)
        if ( ( X2.at(1) - X1.at(1) ) != 0 ) { //  x-coord will be used
            t = ( x - X1.at(1) ) / ( X2.at(1) - X1.at(1) ); // parameter on segment
        } else {
            if ( ( X2.at(2) - X1.at(2) ) != 0 ) { // y-coord will be used
                t = ( y - X1.at(2) ) / ( X2.at(2) - X1.at(2) ); // parameter on segment
            } else {   // z-coord will be used
                t = ( z - X1.at(3) ) / ( X2.at(3) - X1.at(3) ); // parameter on segment
            }
        }

        // is interection on line segment? (0<t<1)
        double EPS = 1.0e-12; // numerical tolerance
        if ( ( 0 < t + EPS ) && ( t - EPS < 1 ) ) {
            intersectCoords.resize(3);
            intersectCoords.at(1) = x;
            intersectCoords.at(2) = y;
            intersectCoords.at(3) = z;
            return true;
        } else {
            intersectCoords.clear();
            return false; // no intersection
        }
    } else {   //different orientation -> no intersection
        intersectCoords.clear();
        return false;  // no intersection
    } // if same orientation
}


int
Quasicontinuum :: intersectionTestSegmentTetrahedra3D(std::vector<FloatArray> &intersectCoords, 
                                                      const FloatArray &A, const FloatArray &B, const FloatArray &C, const FloatArray &D, const FloatArray &X1, const FloatArray &X2)
{
    // test if line (X2-X1) intersects tetrahedra ABCD (for 3D!)
    intersectCoords.clear();
    FloatArray intersectCoord;

    // face 1
    if ( intersectionTestSegmentTrianglePlucker3D(intersectCoord, A, B, C, X1, X2) ) {
        intersectCoords.push_back(intersectCoord);
    }

    // face 2
    if ( intersectionTestSegmentTrianglePlucker3D(intersectCoord, A, D, B, X1, X2) ) {
        intersectCoords.push_back(intersectCoord);
    }

    // face 3
    if ( intersectionTestSegmentTrianglePlucker3D(intersectCoord, B, D, C, X1, X2) ) {
        intersectCoords.push_back(intersectCoord);
    }

    // face 4
    if ( intersectionTestSegmentTrianglePlucker3D(intersectCoord, A, C, D, X1, X2) ) {
        intersectCoords.push_back(intersectCoord);
    }

    return (int)intersectCoords.size();
}


int
Quasicontinuum :: intersectionTestSegmentTriangle2D(std::vector<FloatArray> &intersectCoords, 
                                                    const FloatArray &A, const FloatArray &B, const FloatArray &C, const FloatArray &U1, const FloatArray &U2)
{
    // test if line (X2-X1) intersects tringle ABC (for 2D!)
    intersectCoords.clear();
    FloatArray intersectCoord;

    // edge 1
    if ( intersectionTestSegmentSegment2D(intersectCoord, A, B, U1, U2) ) {
        intersectCoords.push_back(intersectCoord);
    }

    // edge 2
    if ( intersectionTestSegmentSegment2D(intersectCoord, B, C, U1, U2) ) {
        intersectCoords.push_back(intersectCoord);
    }

    // edge 3
    if ( intersectionTestSegmentSegment2D(intersectCoord, A, C, U1, U2) ) {
        intersectCoords.push_back(intersectCoord);
    }

    return (int)intersectCoords.size();
}


bool
Quasicontinuum :: intersectionTestSegmentSegment2D(FloatArray &intersectCoords, const FloatArray &A1, const FloatArray &A2, const FloatArray &B1, const FloatArray &B2)
{
    //intersectionTest of 2 segments A1-A2 B1-B2 (2D only)

    double a1x = A1.at(1);
    double a1y = A1.at(2);
    double a2x = A2.at(1);
    double a2y = A2.at(2);
    double b1x = B1.at(1);
    double b1y = B1.at(2);
    double b2x = B2.at(1);
    double b2y = B2.at(2);

    //  y=f*x+g
    double f1 = std :: numeric_limits< float > :: infinity();
    double f2 = std :: numeric_limits< float > :: infinity();
    if ( a1x != a2x ) {
        f1 = -( a2y - a1y ) / ( a1x - a2x );
    }
    if ( b1x != b2x ) {
        f2 = -( b2y - b1y ) / ( b1x - b2x );
    }

    double g1 = a1y - f1 * a1x;
    double g2 = b1y - f2 * b1x;

    // % parallels
    if ( f1 == f2 ) {
        intersectCoords.clear();
        return false;
    }

    // intersection...
    double x, y;
    if ( a1x == a2x ) {
        x = a1x;
        y = f2 * a1x + g2;
    } else if ( b1x == b2x ) {
        x = b1x;
        y = f1 * b1x + g1;
    } else {
        x = ( g2 - g1 ) / ( f1 - f2 );
        y = f1 * x + g1;
    }

    // is intersection on both segments? (it is not enough to test x-coords only)
    double x1, x2, x3, x4;
    if ( a1x <= a2x ) {
        x1 = a1x;
        x2 = a2x;
    } else {
        x1 = a2x;
        x2 = a1x;
    }
    if ( b1x <= b2x ) {
        x3 = b1x;
        x4 = b2x;
    } else {
        x3 = b2x;
        x4 = b1x;
    }
    double y1, y2, y3, y4;
    if ( a1y <= a2y ) {
        y1 = a1y;
        y2 = a2y;
    } else {
        y1 = a2y;
        y2 = a1y;
    }
    if ( b1y <= b2y ) {
        y3 = b1y;
        y4 = b2y;
    } else {
        y3 = b2y;
        y4 = b1y;
    }

    // test
    double EPS = 1e-10;
    if ( ( x1 <= x + EPS ) && ( x <= x2 + EPS ) && ( x3 <= x + EPS ) && ( x <= x4 + EPS )  &&  ( y1 <= y + EPS ) && ( y <= y2 + EPS ) && ( y3 <= y + EPS ) && ( y <= y4 + EPS ) ) {
        intersectCoords = Vec2(x, y);
        return true;
    } else {
        intersectCoords.clear();
        return false;
    }
}


void
Quasicontinuum :: transformStiffnessTensorToMatrix(FloatMatrix &matrix, const FloatMatrix &tensor)
{
    // transform stiff. tensor 9x9 to stiff. matrix 6x6 (universal for 3D and 2D)
    IntArray indicesI = {1, 2, 3, 2, 1, 1};
    IntArray indicesJ = {1, 2, 3, 3, 3, 2};
    IntArray indicesK = {1, 2, 3, 2, 1, 1};
    IntArray indicesL = {1, 2, 3, 3, 3, 2};

    for ( int ii = 1; ii <= 6; ii++ ) {
        for ( int jj = 1; jj <= 6; jj++ ) {
            int i = indicesI.at(ii);
            int j = indicesJ.at(ii);

            int k = indicesK.at(jj);
            int l = indicesL.at(jj);

            matrix.at(ii, jj) = tensor.at(3 * ( i - 1 ) + k, 3 * ( j - 1 ) + l);
        }
    }
}
} // end namespace oofem
