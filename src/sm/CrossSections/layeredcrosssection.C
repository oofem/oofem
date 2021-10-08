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

#include "sm/CrossSections/layeredcrosssection.h"
#include "sm/Elements/structuralelement.h"
#include "sm/Materials/structuralmaterial.h"
#include "sm/Materials/structuralms.h"
#include "gausspoint.h"
#include "material.h"
#include "floatarray.h"
#include "floatarrayf.h"
#include "floatmatrixf.h"
#include "contextioerr.h"
#include "gaussintegrationrule.h"
#include "mathfem.h"
#include "classfactory.h"
#include "lobattoir.h"
#include "dynamicinputrecord.h"
#include "cltypes.h"
#include "simplecrosssection.h"

namespace oofem {
REGISTER_CrossSection(LayeredCrossSection);


FloatArrayF<6>
LayeredCrossSection :: giveRealStress_3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    if ( gp->giveIntegrationRule()->giveIntegrationDomain() == _Cube || gp->giveIntegrationRule()->giveIntegrationDomain() == _Wedge ) {
        // Determine which layer the gp belongs to. This code assumes that the gauss point are created consistently (through CrossSection::setupIntegrationPoints)
        int ngps = gp->giveIntegrationRule()->giveNumberOfIntegrationPoints();
        int gpnum = gp->giveNumber();
        int gpsperlayer = ngps / this->numberOfLayers;
        int layer = ( gpnum - 1 ) / gpsperlayer + 1;
        auto layerMat = this->domain->giveMaterial( this->giveLayerMaterial(layer) );
        if ( this->layerRots.at(layer) != 0. ) {
            double rot = this->layerRots.at(layer);
            double c = cos(rot * M_PI / 180.);
            double s = sin(rot * M_PI / 180.);

            FloatArrayF<6> rotStrain = {
                c * c * strain.at(1) - c * s * strain.at(6) + s * s * strain.at(2),
                c * c * strain.at(2) + c * s * strain.at(6) + s * s * strain.at(1),
                strain.at(3),
                c * strain.at(4) + s * strain.at(5),
                c * strain.at(5) - s * strain.at(4),
                ( c * c - s * s ) * strain.at(6) + 2 * c * s * ( strain.at(1) - strain.at(2) ),
            };

            auto rotStress = static_cast< StructuralMaterial * >(layerMat)->giveRealStressVector_3d(rotStrain, gp, tStep);

            return {
                c * c * rotStress.at(1) + 2 * c * s * rotStress.at(6) + s * s * rotStress.at(2),
                c * c * rotStress.at(2) - 2 * c * s * rotStress.at(6) + s * s * rotStress.at(1),
                rotStress.at(3),
                c * rotStress.at(4) - s * rotStress.at(5),
                c * rotStress.at(5) + s * rotStress.at(4),
                ( c * c - s * s ) * rotStress.at(6) - c * s * ( rotStress.at(1) - rotStress.at(2) ),
            };
        } else {
            return static_cast< StructuralMaterial * >(layerMat)->giveRealStressVector_3d(strain, gp, tStep);
        }
    } else {
        OOFEM_ERROR("Only cubes and wedges are meaningful for layered cross-sections");
        return zeros<6>();
    }
}


FloatArrayF<6>
LayeredCrossSection :: giveRealStress_3dDegeneratedShell(const FloatArrayF<6> &reducedStrain, GaussPoint *gp, TimeStep *tStep) const
{
    ///@todo - check-V
    return zeros<6>();
}


FloatArrayF<4>
LayeredCrossSection :: giveRealStress_PlaneStrain(const FloatArrayF<4> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not supported");
    return zeros<4>();
}


FloatArrayF<3>
LayeredCrossSection :: giveRealStress_PlaneStress(const FloatArrayF<3> &strain, GaussPoint *masterGp, TimeStep *tStep) const
{
    //strain eps_x, eps_y, gamma_xy
    //stress sig_x, sig_y, tau_xy
    //answer n_x, n_y, n_xy
    
    FloatArray layerStrain;
   
    //double bottom = this->give(CS_BottomZCoord, masterGp);
    //double top = this->give(CS_TopZCoord, masterGp);
    
    auto element = dynamic_cast< StructuralElement * >( masterGp->giveElement() );
    double totThick = 0.0;
    
    FloatArrayF<3> answer;
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
      for (int igp = 0; igp < layerIntegrationPoints.at(layer); igp++) {
        auto layerGp = this->giveSlaveGaussPoint(masterGp, layer - 1, igp);
        auto layerMat = this->domain->giveMaterial( this->giveLayerMaterial(layer) );
        auto interface = static_cast< LayeredCrossSectionInterface * >( element->giveInterface(LayeredCrossSectionInterfaceType) );
        auto lgpw = layerGp->giveWeight();
        
        // resolve current layer z-coordinate
        double layerThick = this->layerThicks.at(layer);
        totThick += layerThick * lgpw;
        //double layerZeta = layerGp->giveNaturalCoordinate(3);
        //double layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );

        // Compute the layer stress
        interface->computeStrainVectorInLayer(layerStrain, strain, masterGp, layerGp, tStep);
        auto reducedLayerStress = dynamic_cast< StructuralMaterial * >(layerMat)->giveRealStressVector_PlaneStress(layerStrain, layerGp, tStep);
        answer.at(1) += reducedLayerStress.at(1) * layerThick* lgpw;
        answer.at(2) += reducedLayerStress.at(2) * layerThick* lgpw;
        answer.at(3) += reducedLayerStress.at(3) * layerThick* lgpw; // * ( 5. / 6. );
      }
    }
    answer*=(1./totThick);
    auto status = static_cast< StructuralMaterialStatus * >( domain->giveMaterial( layerMaterials.at(1) )->giveStatus(masterGp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
    
    return answer;
}


FloatArrayF<1>
LayeredCrossSection :: giveRealStress_1d(const FloatArrayF<1> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not supported");
    return zeros<1>();
}


FloatArrayF<2>
LayeredCrossSection :: giveRealStress_Warping(const FloatArrayF<2> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not supported");
    return zeros<2>();
}


FloatMatrixF<6,6>
LayeredCrossSection :: giveStiffnessMatrix_3d(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    if ( gp->giveIntegrationRule()->giveIntegrationDomain() != _Cube && gp->giveIntegrationRule()->giveIntegrationDomain() != _Wedge ) {
        OOFEM_ERROR("Only cubes and wedges are meaningful for layered cross-sections");
    }
    // Determine which layer the gp belongs to. This code assumes that the gauss point are created consistently (through CrossSection::setupIntegrationPoints)
    int ngps = gp->giveIntegrationRule()->giveNumberOfIntegrationPoints();
    int gpnum = gp->giveNumber();
    int gpsperlayer = ngps / this->numberOfLayers;
    int layer = ( gpnum - 1 ) / gpsperlayer + 1;
    auto layerMat = this->domain->giveMaterial( this->giveLayerMaterial(layer) );
    auto tangent = static_cast< StructuralMaterial * >(layerMat)->give3dMaterialStiffnessMatrix(rMode, gp, tStep);

    if ( this->layerRots.at(layer) != 0. ) {
        double rot = this->layerRots.at(layer);
        double c = cos(rot * M_PI / 180.);
        double s = sin(rot * M_PI / 180.);

        FloatMatrixF<6,6> rotTangent = {
                 c *c,       s *s,  0.,  0.,  0.,         -c *s,
                 s *s,       c *c,  0.,  0.,  0.,          c *s,
                   0.,         0.,  1.,  0.,  0.,            0.,
                   0.,         0.,  0.,   c,   s,            0.,
                   0.,         0.,  0.,  -s,   c,            0.,
            2 * c * s, -2 * c * s,  0.,  0.,  0., c * c - s * s,
        };

        return unrotate(tangent, rotTangent);
    } else {
        return tangent;
    }
}


FloatMatrixF<3,3>
LayeredCrossSection :: giveStiffnessMatrix_PlaneStress(MatResponseMode rMode, GaussPoint *masterGp, TimeStep *tStep) const
{
    FloatMatrixF<3,3> answer;
    double totThick = 0.;
    
    //Average stiffness over all layers
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
      for (int igp=0; igp< layerIntegrationPoints.at(layer); igp++) {
        auto slaveGP = this->giveSlaveGaussPoint(masterGp, layer - 1, igp);
        auto layerMat = this->domain->giveMaterial( this->giveLayerMaterial(layer) );
        auto sgpw = slaveGP->giveWeight();
        double layerThick = this->layerThicks.at(layer);
        totThick += layerThick * sgpw;
        auto subAnswer = dynamic_cast< StructuralMaterial * >(layerMat)->givePlaneStressStiffMtrx(rMode, slaveGP, tStep);
        answer += layerThick * sgpw * subAnswer;
      }
    }
    //answer.at(3,3) *= (5./6.);
    return answer * (1./totThick);
}


FloatMatrixF<4,4>
LayeredCrossSection :: giveStiffnessMatrix_PlaneStrain(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not supported");
    return FloatMatrixF<4,4>();
}


FloatMatrixF<1,1>
LayeredCrossSection :: giveStiffnessMatrix_1d(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not supported");
    return FloatMatrixF<1,1>();
}


FloatArrayF<3>
LayeredCrossSection :: giveGeneralizedStress_Beam2d(const FloatArrayF<3> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    FloatArray layerStrain;
    auto element = static_cast< StructuralElement * >( gp->giveElement() );
    auto interface = static_cast< LayeredCrossSectionInterface * >( element->giveInterface(LayeredCrossSectionInterfaceType) );


    // perform integration over layers
    double bottom = this->give(CS_BottomZCoord, gp);
    double top = this->give(CS_TopZCoord, gp);

    if ( interface == nullptr ) {
        OOFEM_ERROR("element with no layer support encountered");
    }

    FloatArrayF<3> answer;
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
      for (int igp=0; igp< layerIntegrationPoints.at(layer); igp++) {
        auto layerGp = this->giveSlaveGaussPoint(gp, layer - 1, igp);
        auto layerMat = static_cast< StructuralMaterial * >( domain->giveMaterial( layerMaterials.at(layer) ) );
        auto lgpw = layerGp->giveWeight();

        // resolve current layer z-coordinate
        double layerThick = this->layerThicks.at(layer);
        double layerWidth = this->layerWidths.at(layer);
        double layerZeta = layerGp->giveNaturalCoordinate(3);
        double layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );

        // Compute the layer stress
        interface->computeStrainVectorInLayer(layerStrain, strain, gp, layerGp, tStep);

        FloatArrayF<2> reducedLayerStress;
        if ( this->layerRots.at(layer) != 0. ) {
            OOFEM_ERROR("Rotation not supported for beams");
        } else {
            reducedLayerStress = layerMat->giveRealStressVector_2dBeamLayer(layerStrain, layerGp, tStep);
        }

        answer.at(1) += reducedLayerStress.at(1) * layerWidth * layerThick * lgpw; //Nx
        answer.at(2) += reducedLayerStress.at(1) * layerWidth * layerThick * lgpw * layerZCoord;//My
        answer.at(3) += reducedLayerStress.at(2) * layerWidth * layerThick * lgpw * beamShearCoeffxz; //Vz
      }
    }

    // Create material status according to the first layer material
    ///@todo This should be replaced with a general "CrossSectionStatus"
    //CrossSectionStatus *status = new CrossSectionStatus(gp);
    //gp->setMaterialStatus(status);
    auto status = static_cast< StructuralMaterialStatus * >( domain->giveMaterial( layerMaterials.at(1) )->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);

    return answer;
}


FloatArrayF<6>
LayeredCrossSection :: giveGeneralizedStress_Beam3d(const FloatArrayF<6> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not supported");
    return zeros<6>();
}


FloatArrayF<5>
LayeredCrossSection :: giveGeneralizedStress_Plate(const FloatArrayF<5> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    FloatArray layerStrain;
    auto element = static_cast< StructuralElement * >( gp->giveElement() );
    auto interface = static_cast< LayeredCrossSectionInterface * >( element->giveInterface(LayeredCrossSectionInterfaceType) );

    // perform integration over layers
    double bottom = this->give(CS_BottomZCoord, gp);
    double top = this->give(CS_TopZCoord, gp);

    if ( interface == nullptr ) {
        OOFEM_ERROR("element with no layer support encountered");
    }

    FloatArrayF<5> answer;
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
      for (int igp = 0; igp < layerIntegrationPoints.at(layer); igp++ ) {
        auto layerGp = this->giveSlaveGaussPoint(gp, layer - 1, igp);
        auto layerMat = static_cast< StructuralMaterial * >( domain->giveMaterial( layerMaterials.at(layer) ) );
        auto lgpw = layerGp->giveWeight();

        // resolve current layer z-coordinate
        double layerThick = this->layerThicks.at(layer);
        double layerWidth = this->layerWidths.at(layer);
        double layerZeta = layerGp->giveNaturalCoordinate(3);
        double layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );

        // Compute the layer stress
        interface->computeStrainVectorInLayer(layerStrain, strain, gp, layerGp, tStep);

        FloatArrayF<5> reducedLayerStress;
        if ( this->layerRots.at(layer) != 0. ) {
            double rot = this->layerRots.at(layer);
            double c = cos(rot * M_PI / 180.);
            double s = sin(rot * M_PI / 180.);

            FloatArrayF<5> rotStrain = {
                c *c * layerStrain.at(1) - c * s * layerStrain.at(5) + s * s * layerStrain.at(2),
                c * c * layerStrain.at(2) + c * s * layerStrain.at(5) + s * s * layerStrain.at(1),
                c * layerStrain.at(3) + s * layerStrain.at(4),
                c * layerStrain.at(4) - s * layerStrain.at(3),
                ( c * c - s * s ) * layerStrain.at(5) + c * s * ( layerStrain.at(1) - layerStrain.at(2) ),
            }; 

            auto rotStress = layerMat->giveRealStressVector_PlateLayer(rotStrain, layerGp, tStep);

            reducedLayerStress = {
                c *c * rotStress.at(1) + 2 * c * s * rotStress.at(5) + s * s * rotStress.at(2),
                c * c * rotStress.at(2) - 2 * c * s * rotStress.at(5) + s * s * rotStress.at(1),
                c * rotStress.at(3) - s * rotStress.at(4),
                c * rotStress.at(4) + s * rotStress.at(3),
                ( c * c - s * s ) * rotStress.at(5) - c * s * ( rotStress.at(1) - rotStress.at(2) ),
            };
        } else {
            reducedLayerStress = layerMat->giveRealStressVector_PlateLayer(layerStrain, layerGp, tStep);
        }

        answer.at(1) += reducedLayerStress.at(1) * layerWidth * layerThick * lgpw * layerZCoord;
        answer.at(2) += reducedLayerStress.at(2) * layerWidth * layerThick * lgpw * layerZCoord;
        answer.at(3) += reducedLayerStress.at(5) * layerWidth * layerThick * lgpw * layerZCoord;
        answer.at(4) += reducedLayerStress.at(4) * layerWidth * layerThick * lgpw * (5./6.);
        answer.at(5) += reducedLayerStress.at(3) * layerWidth * layerThick * lgpw * (5./6.);
      }
    }

    // now we must update master gp
    // Create material status according to the first layer material
    ///@todo This should be replaced with a general "CrossSectionStatus"
    //CrossSectionStatus *status = new CrossSectionStatus(gp);
    //gp->setMaterialStatus(status);
    auto status = static_cast< StructuralMaterialStatus * >( domain->giveMaterial( layerMaterials.at(1) )->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
    return answer;
}


FloatArrayF<8>
LayeredCrossSection :: giveGeneralizedStress_Shell(const FloatArrayF<8> &strain, GaussPoint *gp, TimeStep *tStep) const
{
    FloatArray layerStrain;
    auto element = static_cast< StructuralElement * >( gp->giveElement() );
    auto interface = static_cast< LayeredCrossSectionInterface * >( element->giveInterface(LayeredCrossSectionInterfaceType) );

    // perform integration over layers
    double bottom = this->give(CS_BottomZCoord, gp);
    double top = this->give(CS_TopZCoord, gp);

    if ( interface == nullptr ) {
        OOFEM_ERROR("element with no layer support encountered");
    }

    FloatArrayF<8> answer;
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
      for (int igp=0; igp< layerIntegrationPoints.at(layer); igp++) {
        auto layerGp = this->giveSlaveGaussPoint(gp, layer - 1, igp);
        auto layerMat = static_cast< StructuralMaterial * >( domain->giveMaterial( layerMaterials.at(layer) ) );
        auto lgpw = layerGp->giveWeight();

        // resolve current layer z-coordinate
        double layerThick = this->layerThicks.at(layer);
        double layerWidth = this->layerWidths.at(layer);
        double layerZeta = layerGp->giveNaturalCoordinate(3);
        double layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );

        // Compute the layer stress
        interface->computeStrainVectorInLayer(layerStrain, strain, gp, layerGp, tStep); // FIXME convert to return value fixed size array.

        FloatArrayF<5> reducedLayerStress;
        if ( this->layerRots.at(layer) != 0. ) {
            double rot = this->layerRots.at(layer);
            double c = cos(rot * M_PI / 180.);
            double s = sin(rot * M_PI / 180.);

            FloatArrayF<5> rotStrain = {
                c *c * layerStrain.at(1) - c * s * layerStrain.at(5) + s * s * layerStrain.at(2),
                c * c * layerStrain.at(2) + c * s * layerStrain.at(5) + s * s * layerStrain.at(1),
                c * layerStrain.at(3) + s * layerStrain.at(4),
                c * layerStrain.at(4) - s * layerStrain.at(3),
                ( c * c - s * s ) * layerStrain.at(5) + c * s * ( layerStrain.at(1) - layerStrain.at(2) ),
            };

            auto rotStress = layerMat->giveRealStressVector_PlateLayer(rotStrain, layerGp, tStep);

            reducedLayerStress = {
                c *c * rotStress.at(1) + 2 * c * s * rotStress.at(5) + s * s * rotStress.at(2),
                c * c * rotStress.at(2) - 2 * c * s * rotStress.at(5) + s * s * rotStress.at(1),
                c * rotStress.at(3) - s * rotStress.at(4),
                c * rotStress.at(4) + s * rotStress.at(3),
                ( c * c - s * s ) * rotStress.at(5) - c * s * ( rotStress.at(1) - rotStress.at(2) ),
            };
        } else {
            reducedLayerStress = layerMat->giveRealStressVector_PlateLayer(layerStrain, layerGp, tStep);
        }

        // 1) membrane terms sx, sy, sxy
        answer.at(1) += reducedLayerStress.at(1) * layerWidth * layerThick *lgpw;
        answer.at(2) += reducedLayerStress.at(2) * layerWidth * layerThick *lgpw;
        answer.at(3) += reducedLayerStress.at(5) * layerWidth * layerThick * lgpw;
        // 2) bending terms mx, my, mxy
        answer.at(4) += reducedLayerStress.at(1) * layerWidth * layerThick * layerZCoord *lgpw;
        answer.at(5) += reducedLayerStress.at(2) * layerWidth * layerThick * layerZCoord * lgpw;
        answer.at(6) += reducedLayerStress.at(5) * layerWidth * layerThick * layerZCoord * lgpw;
        // 3) shear terms qx, qy
        answer.at(7) += reducedLayerStress.at(4) * layerWidth * layerThick *lgpw * (5./6.);
        answer.at(8) += reducedLayerStress.at(3) * layerWidth * layerThick *lgpw * (5./6.);
      }
    }

   
    // now we must update master gp
    ///@todo This should be replaced with a general "CrossSectionStatus"
    //CrossSectionStatus *status = new CrossSectionStatus(gp);
    //gp->setMaterialStatus(status);
    // Create material status according to the first layer material
    auto status = static_cast< StructuralMaterialStatus * >( domain->giveMaterial( layerMaterials.at(1) )->giveStatus(gp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
    
    return answer;
}

FloatArrayF<9>
LayeredCrossSection :: giveGeneralizedStress_ShellRot(const FloatArrayF<9> &strain, GaussPoint *gp, TimeStep *tStep) const
{


  FloatArrayF<9> answer;
  FloatArrayF<8> rstrain;
  for (int i=1; i<=8; i++) {
    rstrain.at(i)=strain.at(i);
  }
  FloatArray ra = this->giveGeneralizedStress_Shell(rstrain, gp, tStep);
  for (int i=1; i<=8; i++) {
    answer.at(i)=ra.at(i);
  }
  answer.at(9) = this->give(CS_DrillingStiffness, gp)*strain.at(9);

  
  ///@todo This should be replaced with a general "CrossSectionStatus"
  //CrossSectionStatus *status = new CrossSectionStatus(gp);
  //gp->setMaterialStatus(status);
  // Create material status according to the first layer material
  auto status = static_cast< StructuralMaterialStatus * >( domain->giveMaterial( layerMaterials.at(1) )->giveStatus(gp) );
  status->letTempStrainVectorBe(strain);
  status->letTempStressVectorBe(answer);

  return answer;
}



FloatArrayF<4>
LayeredCrossSection :: giveGeneralizedStress_MembraneRot(const FloatArrayF<4> &strain, GaussPoint *masterGp, TimeStep *tStep) const
{
    //strain eps_x, eps_y, gamma_xy
    //stress sig_x, sig_y, tau_xy
    //answer n_x, n_y, n_xy
    
    FloatArray layerStrain;
   
    //double bottom = this->give(CS_BottomZCoord, masterGp);
    //double top = this->give(CS_TopZCoord, masterGp);
    
    auto element = dynamic_cast< StructuralElement * >( masterGp->giveElement() );
    double totThick = 0.0;
    
    FloatArrayF<4> answer;
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
      for (int igp=0; igp<layerIntegrationPoints.at(layer); igp++) {
        auto layerGp = this->giveSlaveGaussPoint(masterGp, layer - 1, igp);
        auto layerMat = this->domain->giveMaterial( this->giveLayerMaterial(layer) );
        auto interface = static_cast< LayeredCrossSectionInterface * >( element->giveInterface(LayeredCrossSectionInterfaceType) );
        auto lgpw = layerGp->giveWeight();
        
        // resolve current layer z-coordinate
        double layerThick = this->layerThicks.at(layer);
        totThick += layerThick * lgpw;
        //double layerZeta = layerGp->giveNaturalCoordinate(3);
        //double layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );

        // Compute the layer stress
        interface->computeStrainVectorInLayer(layerStrain, strain, masterGp, layerGp, tStep);
	// extract membrane part only
	FloatArrayF<3> layerStrainMembrane(layerStrain.at(1), layerStrain.at(2), layerStrain.at(3));
        auto reducedLayerStress = dynamic_cast< StructuralMaterial * >(layerMat)->giveRealStressVector_PlaneStress(layerStrainMembrane, layerGp, tStep);
        answer.at(1) += reducedLayerStress.at(1) * layerThick *lgpw;
        answer.at(2) += reducedLayerStress.at(2) * layerThick *lgpw;
        answer.at(3) += reducedLayerStress.at(3) * layerThick *lgpw;
      }
    }
    
    // assume rotation term elastic response
    auto de= this->giveMembraneRotStiffMtrx(ElasticStiffness, masterGp, tStep);
    answer.at(4)=strain.at(4)*de.at(4,4) * totThick;
    answer*=(1./totThick);

    auto status = static_cast< StructuralMaterialStatus * >( domain->giveMaterial( layerMaterials.at(1) )->giveStatus(masterGp) );
    status->letTempStrainVectorBe(strain);
    status->letTempStressVectorBe(answer);
    
    return answer;
}

FloatArrayF<3>
LayeredCrossSection :: giveGeneralizedStress_PlateSubSoil(const FloatArrayF<3> &generalizedStrain, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not supported in given cross-section (yet).");
    return zeros<3>();
}

void
LayeredCrossSection :: giveCharMaterialStiffnessMatrix(FloatMatrix &answer,
                                                       MatResponseMode rMode,
                                                       GaussPoint *gp,
                                                       TimeStep *tStep)
//
// only interface to material class, forcing returned matrix to be in reduced form.
//
{
    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _2dBeam ) {
        answer = this->give2dBeamStiffMtrx(rMode, gp, tStep);
    } else if ( mode == _3dBeam ) {
        answer = this->give3dBeamStiffMtrx(rMode, gp, tStep);
    } else if ( mode == _2dPlate ) {
        answer = this->give2dPlateStiffMtrx(rMode, gp, tStep);
    } else if ( mode == _3dShell ) {
        answer = this->give3dShellStiffMtrx(rMode, gp, tStep);
    } else {
        int ngps = gp->giveIntegrationRule()->giveNumberOfIntegrationPoints();
        int gpnum = gp->giveNumber();
        int gpsperlayer = ngps / this->numberOfLayers;
        int layer = ( gpnum - 1 ) / gpsperlayer + 1;
        auto mat = static_cast< StructuralMaterial * >( domain->giveMaterial( this->giveLayerMaterial(layer) ) );
        if ( mat->hasMaterialModeCapability( gp->giveMaterialMode() ) ) {
            mat->giveStiffnessMatrix(answer, rMode, gp, tStep);
        } else {
            OOFEM_ERROR("unsupported StressStrainMode");
        }
    }
}


FloatMatrixF<5,5>
LayeredCrossSection :: give2dPlateStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const

//
// assumption sigma_z = 0.
//
// General strain layer vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
// returned strain or stress vector has the form:
// 2) strainVectorShell {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy}
//
{
    // perform integration over layers
    double bottom = this->give(CS_BottomZCoord, gp);
    double top = this->give(CS_TopZCoord, gp);

    FloatMatrixF<5,5> answer;
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
      for (int igp=0; igp< layerIntegrationPoints.at(layer); igp++) {
        auto layerGp = giveSlaveGaussPoint(gp, layer - 1, igp);
        auto lgpw = layerGp->giveWeight();

        ///@todo Just using the gp number doesn't nicely support more than 1 gp per layer. Must rethink.
        auto mat = static_cast< StructuralMaterial * >( domain->giveMaterial( this->giveLayerMaterial(layer) ) );
        auto layerMatrix = mat->givePlateLayerStiffMtrx(rMode, layerGp, tStep);
        if ( this->layerRots.at(layer) != 0. ) {
            double rot = this->layerRots.at(layer);
            double c = cos(rot * M_PI / 180.);
            double s = sin(rot * M_PI / 180.);

            FloatMatrixF<5,5> rotTangent = {
                    c * c,      s * s, 0., 0.,         -c *s,
                    s * s,      c * c, 0., 0.,          c *s,
                       0.,         0.,  c,  s,            0.,
                       0.,         0., -s,  c,            0.,
                2 * c * s, -2 * c * s, 0., 0., c * c - s * s,
            };
            layerMatrix = unrotate(layerMatrix, rotTangent);
        }

        //
        // resolve current layer z-coordinate
        //
        double layerThick = this->layerThicks.at(layer);
        double layerWidth  = this->layerWidths.at(layer);
        double layerZeta   = layerGp->giveNaturalCoordinate(3);
        double layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );
        double layerZCoord2 = layerZCoord * layerZCoord;
        //
        // perform integration
        //
        // 1) bending terms mx, my, mxy
        answer.at(1, 1) += layerMatrix.at(1, 1) *lgpw* layerWidth * layerThick * layerZCoord2;
        answer.at(1, 2) += layerMatrix.at(1, 2) *lgpw* layerWidth * layerThick * layerZCoord2;
        answer.at(1, 3) += layerMatrix.at(1, 5) *lgpw* layerWidth * layerThick * layerZCoord2;

        answer.at(2, 1) += layerMatrix.at(2, 1) *lgpw* layerWidth * layerThick * layerZCoord2;
        answer.at(2, 2) += layerMatrix.at(2, 2) *lgpw* layerWidth * layerThick * layerZCoord2;
        answer.at(2, 3) += layerMatrix.at(2, 5) *lgpw* layerWidth * layerThick * layerZCoord2;

        answer.at(3, 1) += layerMatrix.at(5, 1) *lgpw* layerWidth * layerThick * layerZCoord2;
        answer.at(3, 2) += layerMatrix.at(5, 2) *lgpw* layerWidth * layerThick * layerZCoord2;
        answer.at(3, 3) += layerMatrix.at(5, 5) *lgpw* layerWidth * layerThick * layerZCoord2;

        // 2) shear terms qx = qxz, qy = qyz
        answer.at(4, 4) += layerMatrix.at(4, 4) *lgpw* layerWidth * layerThick * (5./6.);
        answer.at(4, 5) += layerMatrix.at(4, 3) *lgpw* layerWidth * layerThick * (5./6.);
        answer.at(5, 4) += layerMatrix.at(3, 4) *lgpw* layerWidth * layerThick * (5./6.);
        answer.at(5, 5) += layerMatrix.at(3, 3) *lgpw* layerWidth * layerThick * (5./6.);
      }
    }
    return answer;
}


FloatMatrixF<8,8>
LayeredCrossSection :: give3dShellStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
//
// assumption sigma_z = 0.
//
// General strain layer vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
// returned strain or stress vector has the form:
// 2) strainVectorShell {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy}
//
{
    // perform integration over layers
    double bottom = this->give(CS_BottomZCoord, gp);
    double top = this->give(CS_TopZCoord, gp);
    double shearcoeff = 5./6.;
    
    FloatMatrixF<8,8> answer;
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
      for (int igp=0; igp<layerIntegrationPoints.at(layer); igp++) {
        auto layerGp = giveSlaveGaussPoint(gp, layer - 1, igp);
        auto lgpw = layerGp->giveWeight();

        ///@todo The logic in this whole class is pretty messy to support both slave-gp's and normal gps. Rethinking the approach is necessary.
        /// Just using the gp number doesn't nicely support more than 1 gp per layer. Must rethink.
        auto mat = static_cast< StructuralMaterial * >( domain->giveMaterial( this->giveLayerMaterial(layer) ) );
        auto layerMatrix = mat->givePlateLayerStiffMtrx(rMode, layerGp, tStep);
        if ( this->layerRots.at(layer) != 0. ) {
            double rot = this->layerRots.at(layer);
            double c = cos(rot);
            double s = sin(rot);

            FloatMatrixF<5,5> rotTangent = {
                    c * c,      s * s, 0., 0.,         -c *s,
                    s * s,      c * c, 0., 0.,          c *s,
                       0.,         0.,  c,  s,            0.,
                       0.,         0., -s,  c,            0.,
                2 * c * s, -2 * c * s, 0., 0., c * c - s * s,
            };
            layerMatrix = unrotate(layerMatrix, rotTangent);
        }

        //
        // resolve current layer z-coordinate
        //
        double layerThick = this->layerThicks.at(layer);
        double layerWidth  = this->layerWidths.at(layer);
        double layerZeta   = layerGp->giveNaturalCoordinate(3);
        double layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );
        double layerZCoord2 = layerZCoord * layerZCoord;
        //
        // perform integration
        //
        // 1) membrane terms sx, sy, sxy
        answer.at(1, 1) += layerMatrix.at(1, 1) * lgpw * layerWidth * layerThick;
        answer.at(1, 2) += layerMatrix.at(1, 2) * lgpw * layerWidth * layerThick;
        answer.at(1, 3) += layerMatrix.at(1, 5) * lgpw * layerWidth * layerThick;

        answer.at(2, 1) += layerMatrix.at(2, 1) * lgpw * layerWidth * layerThick;
        answer.at(2, 2) += layerMatrix.at(2, 2) * lgpw * layerWidth * layerThick;
        answer.at(2, 3) += layerMatrix.at(2, 5) * lgpw * layerWidth * layerThick;

        answer.at(3, 1) += layerMatrix.at(5, 1) * lgpw * layerWidth * layerThick;
        answer.at(3, 2) += layerMatrix.at(5, 2) * lgpw * layerWidth * layerThick;
        answer.at(3, 3) += layerMatrix.at(5, 5) * lgpw * layerWidth * layerThick;

        // 2) bending terms mx, my, mxy

        answer.at(4, 4) += layerMatrix.at(1, 1) * lgpw * layerWidth * layerThick * layerZCoord2;
        answer.at(4, 5) += layerMatrix.at(1, 2) * lgpw * layerWidth * layerThick * layerZCoord2;
        answer.at(4, 6) += layerMatrix.at(1, 5) * lgpw * layerWidth * layerThick * layerZCoord2;

        answer.at(5, 4) += layerMatrix.at(2, 1) * lgpw * layerWidth * layerThick * layerZCoord2;
        answer.at(5, 5) += layerMatrix.at(2, 2) * lgpw * layerWidth * layerThick * layerZCoord2;
        answer.at(5, 6) += layerMatrix.at(2, 5) * lgpw * layerWidth * layerThick * layerZCoord2;

        answer.at(6, 4) += layerMatrix.at(5, 1) * lgpw * layerWidth * layerThick * layerZCoord2;
        answer.at(6, 5) += layerMatrix.at(5, 2) * lgpw * layerWidth * layerThick * layerZCoord2;
        answer.at(6, 6) += layerMatrix.at(5, 5) * lgpw * layerWidth * layerThick * layerZCoord2;

        // 3) shear terms qx, qy
        answer.at(7, 7) += layerMatrix.at(4, 4) * lgpw * layerWidth * layerThick * shearcoeff;
        answer.at(7, 8) += layerMatrix.at(4, 3) * lgpw * layerWidth * layerThick * shearcoeff;
        answer.at(8, 7) += layerMatrix.at(3, 4) * lgpw * layerWidth * layerThick * shearcoeff;
        answer.at(8, 8) += layerMatrix.at(3, 3) * lgpw * layerWidth * layerThick * shearcoeff;
      }
    }
    return answer;
}

FloatMatrixF<9,9>
LayeredCrossSection :: give3dShellRotStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
//
// assumption sigma_z = 0.
//
// General strain layer vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
// returned strain or stress vector has the form:
// 2) strainVectorShell {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy, eps_normalRotation}
//
{
  FloatMatrix d, answer;
  d = this->give3dShellStiffMtrx(rMode, gp, tStep);
  answer.resize(9,9);
  answer.zero();
  answer.assemble(d, {1,2,3,4,5,6,7,8});
  answer.at(9,9) = this->give(CS_DrillingStiffness, gp);
  //answer.printYourself("De");

  return answer;

}

FloatMatrixF<6,6>
LayeredCrossSection :: give3dDegeneratedShellStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    ///@todo - check-V
    return FloatMatrixF<6,6>();
}


FloatMatrixF<3,3>
LayeredCrossSection :: give2dBeamStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
//
// assumption sigma_z = 0.
//
// General strain layer vector has one of the following forms:
// 1) strainVector3d {eps_x,eps_y,eps_z,gamma_yz,gamma_zx,gamma_xy}
//
// returned strain or stress vector has the form:
// 2) strainVectorShell {eps_x,eps_y,gamma_xy, kappa_x, kappa_y, kappa_xy, gamma_zx, gamma_zy}
//
{
    // perform integration over layers
    double bottom = this->give(CS_BottomZCoord, gp);
    double top = this->give(CS_TopZCoord, gp);

    FloatMatrixF<3,3> answer;
    for ( int i = 1; i <= numberOfLayers; i++ ) {
      for (int igp=0; igp<layerIntegrationPoints.at(i); igp++) {
        auto layerGp = giveSlaveGaussPoint(gp, i - 1, igp);
        auto lgpw = layerGp->giveWeight();

        ///@todo The logic in this whole class is pretty messy to support both slave-gp's and normal gps. Rethinking the approach is necessary.
        /// Just using the gp number doesn't nicely support more than 1 gp per layer. Must rethink.
        auto mat = static_cast< StructuralMaterial * >( domain->giveMaterial( this->giveLayerMaterial(i) ) );
        auto layerMatrix = mat->give2dBeamLayerStiffMtrx(rMode, layerGp, tStep);
        if ( this->layerRots.at(i) != 0. ) {
            OOFEM_ERROR("Doesn't support layer rotations.");
        }

        //
        // resolve current layer z-coordinate
        //
        double layerThick = this->layerThicks.at(i);
        double layerWidth  = this->layerWidths.at(i);
        double layerZeta   = layerGp->giveNaturalCoordinate(3);
        double layerZCoord = 0.5 * ( ( 1. - layerZeta ) * bottom + ( 1. + layerZeta ) * top );
        double layerZCoord2 = layerZCoord * layerZCoord;
        //
        // perform integration
        //
        // 1) membrane terms sx
        answer.at(1, 1) += layerMatrix.at(1, 1) * lgpw * layerWidth * layerThick;
        answer.at(1, 3) += layerMatrix.at(1, 2) * lgpw * layerWidth * layerThick;
        // 2) bending terms my
        answer.at(2, 2) += layerMatrix.at(1, 1) * lgpw * layerWidth * layerThick * layerZCoord2;
        answer.at(2, 3) += layerMatrix.at(1, 2) * lgpw * layerWidth * layerThick * layerZCoord2;
        // 3) shear terms qx
        answer.at(3, 1) += layerMatrix.at(2, 1) * lgpw * layerWidth * layerThick * beamShearCoeffxz;
        answer.at(3, 3) += layerMatrix.at(2, 2) * lgpw * layerWidth * layerThick * beamShearCoeffxz;
      }
    }
    return answer;
}


FloatMatrixF<6,6>
LayeredCrossSection :: give3dBeamStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not implemented");
    return FloatMatrixF<6,6>();
}


FloatMatrixF<4,4>
LayeredCrossSection :: giveMembraneRotStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
  auto d = this->giveStiffnessMatrix_PlaneStress(rMode, gp, tStep);
  auto de= this->giveStiffnessMatrix_PlaneStress(ElasticStiffness, gp, tStep);
  auto ds = assemble<4,4>(d, {0, 1, 2}, {0, 1, 2});
  ds.at(4, 4) = 2.0*de.at(3,3); //this->give(CS_DrillingStiffness, gp);
  return ds;
}

FloatMatrixF<3,3>
LayeredCrossSection :: give2dPlateSubSoilStiffMtrx(MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep) const
{
    OOFEM_ERROR("Not implemented");
    return FloatMatrixF<3,3>();
}



FloatArray *
LayeredCrossSection :: imposeStressConstrainsOnGradient(GaussPoint *gp, FloatArray *gradientStressVector3d)
//
// returns modified gradient of stress vector, which is used to
// bring stresses back to yield surface.
//
// imposes zeros on places, where zero stress occurs. if energetically connected
// strain is zero, we do not impose zero there, because stress exist and
// must be taken into account when computing yeld function. In such case
// a problem is assumed to be full 3d with some explicit strain equal to 0.
//
// On the other hand, if some stress is imposed to be zero, we understand
// such case as subspace of 3d case (like a classical plane stess problem, with no
// tracing of ez, sigma_z)
//
{
    MaterialMode mode = gp->giveMaterialMode();
    int size = gradientStressVector3d->giveSize();
    if ( size != 6 ) {
        OOFEM_ERROR("size mismatch");
    }

    switch ( mode ) {
    case _PlateLayer:
        gradientStressVector3d->at(3) = 0.;
        break;
    case _2dBeamLayer:
        for ( int i = 2; i <= 5; i++ ) {
            gradientStressVector3d->at(i) = 0.;
        }

        break;
    default:
        StructuralCrossSection :: imposeStressConstrainsOnGradient(gp, gradientStressVector3d);
    }

    return gradientStressVector3d;
}


FloatArray *
LayeredCrossSection :: imposeStrainConstrainsOnGradient(GaussPoint *gp, FloatArray *gradientStrainVector3d)
//
// returns modified gradient of strain vector, which is used to
// compute plastic strain increment.
//
// imposes zeros on places, where zero strain occurs or energetically connected stress
// is prescribed to be zero.
//
{
    MaterialMode mode = gp->giveMaterialMode();
    if ( gradientStrainVector3d->giveSize() != 6 ) {
        OOFEM_ERROR("size mismatch");
    }

    switch ( mode ) {
    case _PlateLayer:
        gradientStrainVector3d->at(3) = 0.;
        break;
    case _2dBeamLayer:
        for ( int i = 2; i <= 5; i++ ) {
            gradientStrainVector3d->at(i) = 0.;
        }

        break;
    default:
        StructuralCrossSection :: imposeStrainConstrainsOnGradient(gp, gradientStrainVector3d);
    }

    return gradientStrainVector3d;
}


void
LayeredCrossSection :: initializeFrom(InputRecord &ir)
{
    CrossSection :: initializeFrom(ir);

    IR_GIVE_FIELD(ir, numberOfLayers, _IFT_LayeredCrossSection_nlayers);
    if ( numberOfLayers <= 0 ) {
        throw ValueInputException(ir, _IFT_LayeredCrossSection_nlayers, "numberOfLayers <= 0 is not allowed");
    }


    IR_GIVE_FIELD(ir, layerMaterials, _IFT_LayeredCrossSection_layermaterials);
    if ( numberOfLayers != layerMaterials.giveSize() ) {  
        if ( layerMaterials.giveSize() == 1 ) {
            OOFEM_WARNING("Assuming same material in all layers");
            double temp = layerMaterials.at(1);
            layerMaterials.resize(numberOfLayers); layerMaterials.zero();
            layerMaterials.add(temp);
        } else {
            throw ValueInputException(ir, _IFT_LayeredCrossSection_layermaterials, "numberOfLayers does not equal given number of materials. ");
        }
    }

    IR_GIVE_FIELD(ir, layerThicks, _IFT_LayeredCrossSection_thicks);
    if ( numberOfLayers != layerThicks.giveSize() ) {  
        if ( layerThicks.giveSize() == 1 ) {
            OOFEM_WARNING("Assuming same thickness in all layers");
            double temp = layerThicks.at(1);
            layerThicks.resize(numberOfLayers); layerThicks.zero();
            layerThicks.add(temp);
        } else {
            throw ValueInputException(ir, _IFT_LayeredCrossSection_thicks, "numberOfLayers does not equal given number of thicknesses. ");
        }
    }

    layerWidths.resize(numberOfLayers);
    layerWidths.zero();
    IR_GIVE_OPTIONAL_FIELD(ir, layerWidths, _IFT_LayeredCrossSection_widths);
    layerRots.resize(numberOfLayers);
    layerRots.zero();
    IR_GIVE_OPTIONAL_FIELD(ir, layerRots, _IFT_LayeredCrossSection_layerRotations);

    if ( numberOfLayers != layerRots.giveSize() ) {  //|| ( numberOfLayers != layerWidths.giveSize() ) ) || numberOfLayers != layerThicks.giveSize() || numberOfLayers != layerMaterials.giveSize()
        throw ValueInputException(ir, _IFT_LayeredCrossSection_layerRotations, "numberOfLayers does not equal given number of layer rotations. ");
    }

    // Interface materials // add check if correct numbers
    interfacerMaterials.resize(numberOfLayers - 1);
    interfacerMaterials.zero();
    IR_GIVE_OPTIONAL_FIELD(ir, interfacerMaterials, _IFT_LayeredCrossSection_interfacematerials);

    numberOfIntegrationPoints = 1;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfIntegrationPoints, _IFT_LayeredCrossSection_nintegrationpoints);
    this->layerIntegrationPoints.resize(numberOfLayers);
    for (int i=1; i<=numberOfLayers; i++) this->layerIntegrationPoints.at(i)=numberOfIntegrationPoints;
    
    IR_GIVE_OPTIONAL_FIELD(ir, layerIntegrationPoints, _IFT_LayeredCrossSection_nlayerintegrationpoints);
    if (layerIntegrationPoints.giveSize()!=numberOfLayers) {
      throw ValueInputException(ir, _IFT_LayeredCrossSection_nlayerintegrationpoints, "size of layerIntegrationPoints does not equal given number of layers. ");
    }


    this->totalThick = layerThicks.sum();
    // read z-coordinate of mid-surface measured from bottom layer
    midSurfaceZcoordFromBottom = 0.5 * this->computeIntegralThick();  // Default: geometric midplane
    midSurfaceXiCoordFromBottom = 1.0; // add to IR
    IR_GIVE_OPTIONAL_FIELD(ir, midSurfaceZcoordFromBottom, _IFT_LayeredCrossSection_midsurf);

    this->setupLayerMidPlanes();
    
    this->area = this->layerThicks.dotProduct(this->layerWidths);
    IR_GIVE_OPTIONAL_FIELD(ir, beamShearCoeffxz, _IFT_LayeredCrossSection_shearcoeff_xz);

    double value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_SimpleCrossSection_drillStiffness);
    propertyDictionary.add(CS_DrillingStiffness, value);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_SimpleCrossSection_relDrillStiffness);
    propertyDictionary.add(CS_RelDrillingStiffness, value);

    value = 0.0;
    IR_GIVE_OPTIONAL_FIELD(ir, value, _IFT_SimpleCrossSection_drillType);
    propertyDictionary.add(CS_DrillingType, value);

}

void LayeredCrossSection :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralCrossSection :: giveInputRecord(input);

    input.setField(this->numberOfLayers, _IFT_LayeredCrossSection_nlayers);
    input.setField(this->layerMaterials, _IFT_LayeredCrossSection_layermaterials);
    input.setField(this->layerThicks, _IFT_LayeredCrossSection_thicks);
    input.setField(this->layerWidths, _IFT_LayeredCrossSection_widths);
    input.setField(this->layerRots, _IFT_LayeredCrossSection_layerRotations);
    input.setField(this->interfacerMaterials, _IFT_LayeredCrossSection_interfacematerials);
    input.setField(this->layerIntegrationPoints, _IFT_LayeredCrossSection_nlayerintegrationpoints);
    input.setField(this->midSurfaceZcoordFromBottom, _IFT_LayeredCrossSection_midsurf);
}

void LayeredCrossSection :: createMaterialStatus(GaussPoint &iGP)
{
    for ( int i = 1; i <= numberOfLayers; i++ ) {
      for (int k=0; k<layerIntegrationPoints.at(i); k++) {
        GaussPoint *layerGp = giveSlaveGaussPoint(& iGP, i - 1, k);
        StructuralMaterial *mat = static_cast< StructuralMaterial * >( domain->giveMaterial( this->giveLayerMaterial(i) ) );
        MaterialStatus *matStat = mat->CreateStatus(layerGp);
        layerGp->setMaterialStatus(matStat);
      }
    }
}


void
LayeredCrossSection :: setupLayerMidPlanes()
{
    // z-coord of each layer midplane measured from the global cross section z-coord
    this->layerMidZ.resize(this->numberOfLayers);
    double layerBottomZ = -midSurfaceZcoordFromBottom; // initialize to the bottom coord
    for ( int j = 1; j <= numberOfLayers; j++ ) {
        double thickness = this->layerThicks.at(j);
        this->layerMidZ.at(j) = layerBottomZ + thickness * 0.5;
        layerBottomZ += thickness;
    }
}


Material *
LayeredCrossSection :: giveMaterial(IntegrationPoint *ip) const
{
    ///@todo We should keep track in integration point (integration rule) what material from layer is assigned. Otherwise difficulties due to different elements and IP numbering.
    if ( ip->giveIntegrationRule()->giveIntegrationDomain() == _Cube ||
        ip->giveIntegrationRule()->giveIntegrationDomain() == _Wedge
    ) {
        return domain->giveMaterial( layerMaterials.at(1) );
        //return this->domain->giveMaterial( this->giveLayerMaterial(ip->giveNumber()) );
    }
    
    if ( ip->hasSlaveGaussPoint() ) {
        return domain->giveMaterial( layerMaterials.at(1) );//virtual master, has no material assigned in input file
    } else {
        return domain->giveMaterial( layerMaterials.at(1) );//virtual master, has no material assigned in input file
        //OOFEM_ERROR("Not implemented.")
    }
    return nullptr;
}


int
LayeredCrossSection :: setupIntegrationPoints(IntegrationRule &irule, int npoints, Element *element)
{
    ///@todo We must send arrays for integration points instead of just a single scalar.
    if ( element->giveIntegrationDomain() == _Cube ) {
#if 0
        ///@todo "npoints" should be an intarray
        return irule.SetUpPointsOnCubeLayers(npoints.at(1), npoints.at(2), this->numberOfIntegrationPoints,
                                             element->giveMaterialMode(), this->layerThicks);

#else
        int points1 = ( int ) floor(cbrt( double ( npoints ) ) + 0.5);
        // If numberOfIntegrationPoints > 0 then use that, otherwise use the element's default.
        return irule.SetUpPointsOnCubeLayers(points1, points1, this->numberOfIntegrationPoints ? numberOfIntegrationPoints : points1,
                                             element->giveMaterialMode(), this->layerThicks);

#endif
    } else if ( element->giveIntegrationDomain() == _Wedge ) {
#if 0
        ///@todo "npoints" should be an intarray
        return irule.SetUpPointsOnWedgeLayers(npoints.at(1), this->numberOfIntegrationPoints,
                                              element->giveMaterialMode(), this->layerThicks);

#else
        if ( npoints == 2 ) {
            return irule.SetUpPointsOnWedgeLayers(1, this->numberOfIntegrationPoints,
                                                  element->giveMaterialMode(), this->layerThicks);
        } else {
            return irule.SetUpPointsOnWedgeLayers(3, this->numberOfIntegrationPoints,
                                                  element->giveMaterialMode(), this->layerThicks);
        }
#endif
    } else {
        return irule.setUpIntegrationPoints( element->giveIntegrationDomain(), npoints, element->giveMaterialMode() );
    }
}

int
LayeredCrossSection :: setupIntegrationPoints(IntegrationRule &irule, int nPointsXY, int nPointsZ, Element *element)
{
    switch ( element->giveIntegrationDomain() ) {
    case _3dDegShell:
        return irule.SetUpPointsOn3dDegShellLayers(nPointsXY, nPointsZ, element->giveMaterialMode(), this->layerThicks);
    default:
        OOFEM_ERROR("Unknown mode (%d)", element->giveIntegrationDomain());
    }
    return 0;
}



int
LayeredCrossSection :: giveSlaveGPIndex (int ilayer, int igp) const
{
  // slave gps stored at master in sequence
  // need to take into account variable number of IPs per layer
  int indx=0;
  // count number of gps in preceeding layers
  for (int i=0; i<ilayer; i++) indx+=layerIntegrationPoints(i);
  indx+=igp;
  return indx;
}


GaussPoint *
LayeredCrossSection :: giveSlaveGaussPoint(GaussPoint *masterGp, int ilayer, int igp) const
//
// return the i-th slave gauss point of master gp
// if slave gp don't exists - create them
// ilayer and igp numbered from 0
//
{
  auto slave = masterGp->giveSlaveGaussPoint(giveSlaveGPIndex(ilayer,igp));
    if ( slave == nullptr ) {
        // check for proper dimensions - slave can be NULL if index too high or if not
        // slaves previously defined
        if ( ilayer > this->numberOfLayers ) {
            OOFEM_ERROR("no such layer defined");
        }

        // create new slave record in masterGp
        // (requires that this is friend of gp)
        const auto &masterCoords = masterGp->giveNaturalCoordinates();
        // resolve slave material mode
        auto masterMode = masterGp->giveMaterialMode();
        auto slaveMode = this->giveCorrespondingSlaveMaterialMode(masterMode);

        double bottom = this->give(CS_BottomZCoord, masterGp);
        double top = this->give(CS_TopZCoord, masterGp);

        ///@todo Generalize to multiple integration points per layer
        int numberOfSlaves = 0;
        for (int i=0; i<numberOfLayers; i++) numberOfSlaves+= layerIntegrationPoints(i);
        masterGp->gaussPoints.resize( numberOfSlaves);
        // helper 1d rule
        double currentZTopCoord = -midSurfaceZcoordFromBottom;
        for ( int j = 0; j < numberOfLayers; j++ ) {
          FloatArray sgpc(layerIntegrationPoints.at(j+1));
          FloatArray sgpw(layerIntegrationPoints.at(j+1));
          GaussIntegrationRule::giveLineCoordsAndWeights(layerIntegrationPoints.at(j+1), sgpc, sgpw);

          currentZTopCoord += this->layerThicks.at(j + 1);
          for (int k = 0; k < layerIntegrationPoints.at(j+1); k++) { 
            FloatArray zCoord(3);

            double currentZCoord = currentZTopCoord - this->layerThicks.at(j + 1) * (1+sgpc(k)) / 2.0; // z-coord of layer gp
            if ( masterCoords.giveSize() > 0 ) {
                zCoord.at(1) = masterCoords.at(1); // gp x-coord of mid surface
            }

            if ( masterCoords.giveSize() > 1 ) {
                zCoord.at(2) = masterCoords.at(2); // gp y-coord of mid surface
            }

            zCoord.at(3) = ( 2.0 * currentZCoord - top - bottom ) / ( top - bottom );
            //printf("SGP %d: currentZTopCoord %e, currentZCoord %e\n", j*numberOfIntegrationPoints+k, currentZTopCoord, currentZCoord);
            // in gp - is stored isoparametric coordinate (-1,1) of z-coordinate
            masterGp->gaussPoints [giveSlaveGPIndex(j,k)] = new GaussPoint(masterGp->giveIntegrationRule(), j + 1, zCoord, sgpw(k)/2.0, slaveMode);

            // test - remove!
//             masterGp->gaussPoints [ j ] = new GaussPoint(masterGp->giveIntegrationRule(), j + 1, zCoord, 1.0, slaveMode);
          }
        }

        slave = masterGp->gaussPoints [giveSlaveGPIndex(ilayer, igp) ];
    }

    return slave;
}

double
LayeredCrossSection :: computeIntegralThick() const
{
    return totalThick;
}


void
LayeredCrossSection :: printYourself()
// Prints the receiver on screen.
{
    printf("Cross Section with properties: \n");
    propertyDictionary.printYourself();
    printf("Layer Materials: \n");
    layerMaterials.printYourself();
    printf("Thickness of each layer: \n");
    layerThicks.printYourself();
    if ( layerWidths.giveSize() ) {
        printf("Width of each layer: \n");
        layerWidths.printYourself();
    }
    printf("Number of integration points per layer: %i \n", this->numberOfIntegrationPoints);
    printf("MidSurfaceZCoordinate from bottom: %f \n", midSurfaceZcoordFromBottom);
}


void
LayeredCrossSection :: saveIPContext(DataStream &stream, ContextMode mode, GaussPoint *masterGp)
{
    CrossSection :: saveIPContext(stream, mode, masterGp);

    // saved master gp record;
    // and now save slave gp of master:
    for ( int i = 1; i <= numberOfLayers; i++ ) {
      for (int k=0; k<layerIntegrationPoints.at(i); k++) {
        GaussPoint *slaveGP = this->giveSlaveGaussPoint(masterGp, i - 1, k);
        StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( domain->giveMaterial( layerMaterials.at(i) ) );
        mat->saveIPContext(stream, mode, slaveGP);
      }
    }
}


void
LayeredCrossSection :: restoreIPContext(DataStream &stream, ContextMode mode, GaussPoint *masterGp)
{
    CrossSection :: restoreIPContext(stream, mode, masterGp);

    // and now save slave gp of master:
    for ( int i = 1; i <= numberOfLayers; i++ ) {
      for (int k=0; k<layerIntegrationPoints.at(i); k++) {
        // creates also slaves if they don't exists
        GaussPoint *slaveGP = this->giveSlaveGaussPoint(masterGp, i - 1, k);
        StructuralMaterial *mat = dynamic_cast< StructuralMaterial * >( domain->giveMaterial( layerMaterials.at(i) ) );
        mat->restoreIPContext(stream, mode, slaveGP);
      }
    }
}


MaterialMode
LayeredCrossSection :: giveCorrespondingSlaveMaterialMode(MaterialMode masterMode)
//
// returns corresponding slave material mode to master mode
//
{
    if ( masterMode == _2dPlate ) {
        return _PlateLayer;
    } else if ( masterMode == _2dBeam ) {
        return _2dBeamLayer;
    } else if (( masterMode == _PlaneStress ) || ( masterMode == _PlaneStressRot )) {
        return _PlaneStress;    
    } else if (( masterMode == _3dShell ) || (masterMode == _3dShellRot)) {
        return _PlateLayer;
    } else if ( masterMode == _3dDegeneratedShell ) {
        return _3dDegeneratedShell;
    } else if ( masterMode == _3dMat ) {
        return _3dMat;
    } else {
        throw std::runtime_error("unsupported material mode");
    }

    return _Unknown;
}


double
LayeredCrossSection :: give(CrossSectionProperty aProperty, GaussPoint *gp) const
{
    if ( aProperty == CS_Thickness ) {
        return this->computeIntegralThick();
    } else if ( aProperty == CS_TopZCoord ) {
        this->computeIntegralThick();
        return totalThick - midSurfaceZcoordFromBottom;
    } else if ( aProperty == CS_BottomZCoord ) {
        return -midSurfaceZcoordFromBottom;
    } else if ( aProperty == CS_Area ) {
        return this->giveArea();
    } else if ( aProperty == CS_NumLayers ) {
        return this->numberOfLayers;
    //} else if (aProperty == CS_Layer ) {
    //    return this->giveLayer(gp);
    }

    return CrossSection :: give(aProperty, gp);
}

int 
LayeredCrossSection :: giveLayer(GaussPoint *gp) const	//@todo: works only for equal thickness of each layer
{
    FloatArray lCoords;
    int noLayers = this->giveNumberOfLayers();
    double dh = 2.0/noLayers;
    lCoords = gp->giveNaturalCoordinates();
    double lowXi = -1.0;

    for (int i = 1; i <= noLayers; i++)
    {
        if (lCoords.at(3) > lowXi && lCoords.at(3) < lowXi+dh)
        {
            return i;
        }
        lowXi+=dh;
    }
    OOFEM_ERROR("LayeredCrossSection :: giveLayer - the actual integration point can not be associated with a layer in the cross section");
}

double
LayeredCrossSection :: give(CrossSectionProperty aProperty, const FloatArray &coords, Element *elem, bool local) const
{
    if ( aProperty == CS_Thickness ) {
        return this->computeIntegralThick();
    } else if ( aProperty == CS_TopZCoord ) {
        this->computeIntegralThick();
        return totalThick - midSurfaceZcoordFromBottom;
    } else if ( aProperty == CS_BottomZCoord ) {
        return -midSurfaceZcoordFromBottom;
    } else if ( aProperty == CS_Area ) {
        return this->giveArea();
    } else if ( aProperty == CS_NumLayers ) {
        return this->numberOfLayers;
    }

    return CrossSection :: give(aProperty, coords, elem, local);
}


int
LayeredCrossSection :: giveNumberOfLayers() const
{
    return this->numberOfLayers;
}

double
LayeredCrossSection :: giveArea() const
{
    return area;
}


bool LayeredCrossSection :: isCharacteristicMtrxSymmetric(MatResponseMode rMode) const
{
    for ( int i = 1; i <= this->numberOfLayers; i++ ) {
        if ( !this->domain->giveMaterial( this->giveLayerMaterial(i) )->isCharacteristicMtrxSymmetric(rMode) ) {
            return false;
        }
    }
    return true;
}


void
LayeredCrossSection :: giveInterfaceXiCoords(FloatArray &answer) const
{
    // returns an array with the xi-coords corresponding to the boundaries where
    // the layers meet (size = number of layers -1)

    int numInterfaces = this->giveNumberOfLayers() - 1;
    answer.resize(numInterfaces);
    double totalThickness = this->computeIntegralThick();
    for ( int i = 1; i <= numInterfaces; i++  ) {
        double midZ = this->giveLayerMidZ(i);
        double interfaceZ  = midZ + this->giveLayerThickness(i) * 0.5;
        answer.at(i) = interfaceZ * ( 2.0 / totalThickness );
    }
}

void
LayeredCrossSection :: setupLayeredIntegrationRule(std :: vector< std :: unique_ptr< IntegrationRule > > &integrationRulesArray, Element *el, int numInPlanePoints)
{
    // Loop over each layer and set up an integration rule as if each layer was an independent element
    // @todo - only works for wedge integration at the moment
    int numberOfLayers     = this->giveNumberOfLayers();
    int numPointsThickness = this->giveNumIntegrationPointsInLayer();

    integrationRulesArray.clear();
    integrationRulesArray.reserve( numberOfLayers );
    for ( int i = 0; i < numberOfLayers; i++ ) {
        integrationRulesArray.emplace_back( new LayeredIntegrationRule(i + 1, el) );
        integrationRulesArray.back()->SetUpPointsOnWedge(numInPlanePoints, numPointsThickness, _3dMat);
    }
    this->mapLayerGpCoordsToShellCoords(integrationRulesArray);
}


void
LayeredCrossSection :: mapLayerGpCoordsToShellCoords(std :: vector< std :: unique_ptr< IntegrationRule > > &layerIntegrationRulesArray)
/*
 * Maps the local xi-coord (z-coord) in each layer [-1,1] to the corresponding
 * xi-coord in the cross section coordinate system.
 * Also renames the gp numbering from layerwise to global (1,2,1,2 -> 1,2,3,4)
 *  xi
 * --------  1               --------  1
 |           |                         |
 |           |                         |
 | -------- -1       =>      --------  x
 | --------  1               --------  x
 |           |                         |
 | -------- -1               -------- -1
 */
{
    double scaleFactor = 0.999; // Will be numerically unstable with xfem if the endpoints lie at +-1
    double totalThickness = this->computeIntegralThick();
    int number = 1;
    for ( int layer = 1; layer <= numberOfLayers; layer++ ) {
        for ( GaussPoint *gp: *layerIntegrationRulesArray [ layer - 1 ] ) {

            // Map local layer cs to local shell cs
            double zMid_i = this->giveLayerMidZ(layer); // global z-coord
            double xiMid_i = 1.0 - 2.0 * ( totalThickness - this->midSurfaceZcoordFromBottom - zMid_i ) / totalThickness; // local z-coord
            double deltaxi = gp->giveNaturalCoordinates().at(3) * this->giveLayerThickness(layer) / totalThickness; // distance from layer mid
            double xinew = xiMid_i + deltaxi * scaleFactor;
            FloatArray lcoords = gp->giveNaturalCoordinates();
            lcoords.at(3) = xinew;
            gp->setNaturalCoordinates(lcoords);
            gp->number = number;   // fix gp ordering
            number++;
        }
    }
}


LayeredIntegrationRule :: LayeredIntegrationRule(int n, Element *e,
                                                 int startIndx, int endIndx, bool dynamic) :
    IntegrationRule(n, e, startIndx, endIndx, dynamic) { }

LayeredIntegrationRule :: LayeredIntegrationRule(int n, Element *e) :
    IntegrationRule(n, e) { }


int
LayeredIntegrationRule :: SetUpPointsOnWedge(int nPointsTri, int nPointsThickness, MaterialMode mode)
{
    // Set up integration rule for a specific layer

    int nPoints = nPointsTri * nPointsThickness;
    //@todo - is not really a Gauss point but rather a hybrid.
    this->gaussPoints.resize( nPoints );

    // uses Gauss integration in the plane and Lobatto in the thickness
    FloatArray coords_xi1, coords_xi2, coords_xi, weights_tri, weights_thickness;
    GaussIntegrationRule   :: giveTriCoordsAndWeights(nPointsTri, coords_xi1, coords_xi2, weights_tri);
    //LobattoIntegrationRule :: giveLineCoordsAndWeights(nPointsThickness, coords_xi, weights_thickness );
    GaussIntegrationRule :: giveLineCoordsAndWeights(nPointsThickness, coords_xi, weights_thickness);

    // Assumes that the integration rules of the layers are the same such that the ordering of the ip's are also
    // the same =>  upperInterfacePoints.at(i) of one layer is paired with lowerInterfacePoints.at(i) of the next.
    // This will be used to estimate interlaminar stresses, sice values in the two ip will generally be different
    // due to beam/plate/shell theory assumptions.
    if ( nPointsThickness != 1 ) { // otherwise there are no points on the interface
        this->lowerInterfacePoints.resize(nPointsTri);
        this->upperInterfacePoints.resize(nPointsTri);
    }
    for ( int i = 1, ind = 0; i <= nPointsThickness; i++ ) {
        for ( int j = 1; j <= nPointsTri; j++ ) {
            this->gaussPoints [ ind ] =
                new GaussPoint(this, 1, {coords_xi1.at(j), coords_xi2.at(j), coords_xi.at(i)},
                               weights_tri.at ( j ) *weights_thickness.at ( i ), mode);

            // store interface points
            if ( i == 1 && nPointsThickness > 1 ) { //then lower surface
                this->lowerInterfacePoints.at(j) = ind;
            } else if ( i == nPointsThickness && nPointsThickness > 1 ) {  //then upper surface
                this->upperInterfacePoints.at(j) = ind;
            }
            ind++;
        }
    }
    return this->giveNumberOfIntegrationPoints();
}


int
LayeredCrossSection :: checkConsistency()
//
// check internal consistency
// mainly tests, whether material and crossSection data
// are safe for conversion to "Structural" versions
//
{
    int result = 1;
    for ( int i = 1; this->giveNumberOfLayers(); i++ ) {
        Material *mat = this->giveDomain()->giveMaterial( this->giveLayerMaterial(i) );
        if ( !dynamic_cast< StructuralMaterial * >(mat) ) {
            OOFEM_WARNING("material %s without structural support", mat->giveClassName() );
            result = 0;
            continue;
        }
    }
    return result;
}


int
LayeredCrossSection :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    if ( gp->giveIntegrationRule()->giveIntegrationDomain() == _Cube || gp->giveIntegrationRule()->giveIntegrationDomain() == _Wedge ) {
        // Determine which layer the gp belongs to. This code assumes that the gauss point are created consistently (through CrossSection::setupIntegrationPoints)
        int ngps = gp->giveIntegrationRule()->giveNumberOfIntegrationPoints();
        int gpnum = gp->giveNumber();
        int gpsperlayer = ngps / this->numberOfLayers;
        int layer = ( gpnum - 1 ) / gpsperlayer + 1;
        Material *layerMat = this->domain->giveMaterial( this->giveLayerMaterial(layer) );
        if ( this->layerRots.at(layer) != 0. ) {
            FloatArray rotVal; // the requested value in the material c.s.
            InternalStateValueType valType = giveInternalStateValueType(type);
            double rot = this->layerRots.at(layer);
            double c = cos(rot * M_PI / 180.);
            double s = sin(rot * M_PI / 180.);

            int ret = layerMat->giveIPValue(rotVal, gp, type, tStep);
            if ( ret == 0 ) {
                return 0;
            }

            // Determine how to rotate it according to the value type
            if ( valType == ISVT_TENSOR_S3 ) {
                answer = {
                    c *c * rotVal.at(1) + 2 * c * s * rotVal.at(6) + s * s * rotVal.at(2),
                    c * c * rotVal.at(2) - 2 * c * s * rotVal.at(6) + s * s * rotVal.at(1),
                    rotVal.at(3),
                    c * rotVal.at(4) - s * rotVal.at(5),
                    c * rotVal.at(5) + s * rotVal.at(4),
                    ( c * c - s * s ) * rotVal.at(6) - c * s * ( rotVal.at(1) - rotVal.at(2) ),
                };
            } else if ( valType == ISVT_TENSOR_S3E ) {
                answer = {
                    c *c * rotVal.at(1) + c * s * rotVal.at(6) + s * s * rotVal.at(2),
                    c * c * rotVal.at(2) - c * s * rotVal.at(6) + s * s * rotVal.at(1),
                    rotVal.at(3),
                    c * rotVal.at(4) - s * rotVal.at(5),
                    c * rotVal.at(5) + s * rotVal.at(4),
                    ( c * c - s * s ) * rotVal.at(6) - 2 * c * s * ( rotVal.at(1) - rotVal.at(2) ),
                };
            } else if ( valType == ISVT_VECTOR ) {
                answer = {
                    c *rotVal.at(1) - s * rotVal.at(2), s * rotVal.at(1) +c * rotVal.at(2), rotVal.at(3)
                };
            } else if ( valType == ISVT_SCALAR ) {
                answer = rotVal;
            } else {
                return 0;
            }
            return 1;
        } else {
            return layerMat->giveIPValue(answer, gp, type, tStep);
        }
    } else {
        //return CrossSection :: giveIPValue(answer, gp, type, tStep);

        ///@todo so far this only works for elements where each layer has its own integration rule
        int layer = gp->giveIntegrationRule()->giveNumber();
        return this->giveDomain()->giveMaterial( this->giveLayerMaterial(layer) )->giveIPValue(answer, gp, type, tStep);
    }
}


double
LayeredCrossSection :: give(int aProperty, GaussPoint* gp) const
{
    double average = 0.;
    for ( int layer = 1; layer <= numberOfLayers; ++layer ) {
        Material *mat = this->giveDomain()->giveMaterial( giveLayerMaterial(layer) );
        average += mat->give(aProperty, gp) * giveLayerThickness(layer);
    }
    return average / this->totalThick;
}

} // end namespace oofem
