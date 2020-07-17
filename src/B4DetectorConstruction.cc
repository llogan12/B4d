//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B4DetectorConstruction.hh"
#include "G4AssemblyVolume.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Tubs.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSCylinderSurfaceFlux.hh"
#include "G4SDManager.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//static unsigned int layers = 1;
B4DetectorConstruction::B4DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::~B4DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
 

  // Envelope parameters
  //
  G4double env_sizeXY = 2.5*cm, env_sizeZ = 20*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_Galactic");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
/*  G4double world_sizeXY = 3*cm;
  G4double world_sizeZ  = 20.5*cm;
 
   G4Element *element_C  = new G4Element("Carbon",   "C",  6,  12.011*g/mole);
  G4Element *element_H  = new G4Element("Hydrogen", "H",  1,  1.00794*g/mole);
  G4Element *element_O  = new G4Element("Oxyzen",   "O",  8,  15.999*g/mole);
 G4Material* world_mat = new G4Material("PMMA", 1.18*g/cm3, 3);
 world_mat -> AddElement(element_C, 5);
 world_mat -> AddElement(element_H, 8);
 world_mat -> AddElement(element_O, 2);
*/
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        env_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Envelope
  //  
  
/*

 G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  */     
  // Shape 1
  //  
//   G4Material* shape0_mat = nist->FindOrBuildMaterial("G4_Galactic");
//   G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_Galactic");
   G4Material* shape2_mat = nist->FindOrBuildMaterial("G4_Be");
//   G4Material* shape3_mat = nist->FindOrBuildMaterial("G4_Galactic");
//   G4Material* shape4_mat = nist->FindOrBuildMaterial("G4_GRAPHITE");

/*  G4Box* solidShape0 =
    new G4Box ("Shape0",                                //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 10*mm);         //its size

  G4LogicalVolume* logicShape0  =
    new G4LogicalVolume(solidShape0,         //its solid
                        shape0_mat,          //its material
                        "Shape0");           //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0, 0, -9*cm),//at (0,0,0)
                    logicShape0,             //its logical volume
                    "Shape0",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    false);          //overlaps checking

    
    G4Box* solidShape1 =    
    new G4Box ("Shape1",                                //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.1*mm);         //its size
      
  G4LogicalVolume* logicShape1  =                         
    new G4LogicalVolume(solidShape1,         //its solid
                        shape1_mat,          //its material
                        "Shape1");           //its name
               
    new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0, 0, 9.9*mm),//at (0,0,0)
                    logicShape1,             //its logical volume
                    "Shape1",                //its name
                    logicShape0,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    false);          //overlaps checking
*/ 
  G4Box* solidShape2 =    
    new G4Box ("Shape2",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*mm);         //its size
      
  G4LogicalVolume* logicShape2  =                         
    new G4LogicalVolume(solidShape2,         //its solid
                        shape2_mat,            //its material
                        "Shape2");           //its name
  
/*  G4Box* solidShape3 =
    new G4Box ("Shape3",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY,0.1*cm);         //its size

  G4LogicalVolume* logicShape3  =
    new G4LogicalVolume(solidShape3,         //its solid
                        shape3_mat,            //its material
                        "Shape3");           //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0, -3.8*cm),   //at (0,0,0)
                    logicShape3,             //its logical volume
                    "Shape3",                //its name
                    logicWorld,             //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    false);          //overlaps checking
*/  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0, -4.2*cm),   //at (0,0,0)
                    logicShape2,             //its logical volume
                    "Shape2",                //its name
                    logicWorld,             //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    false);          //overlaps checking
/*
 G4Box* solidShape4 =
    new G4Box ("Shape4",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.6*cm);         //its size

  G4LogicalVolume* logicShape4  =
    new G4LogicalVolume(solidShape4,         //its solid
                        shape4_mat,            //its material
                        "Shape4");           //its name

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,-3.4*cm),   //at (0,0,0)
                    logicShape4,             //its logical volume
                    "Shape4",                //its name
                    logicWorld,             //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    false);          //overlaps checking
*/
  fScoringVolume = logicShape2;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*void RE06DetectorConstruction::SetupDetectors()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);
  G4String filterName, particleName;

  G4SDParticleFilter* gammaFilter 
    = new G4SDParticleFilter(filterName="neutronFilter",particleName="neutron");  
  G4SDParticleFilter* positronFilter 
    = new G4SDParticleFilter(filterName="positronFilter",particleName="e+");
//////////////////////////////////////////////////////   
 //   G4VPrimitiveScorer* primitive;
 //   primitive = new G4PSNofSecondary("nNeutron",j);
 //   primitive->SetFilter(gammaFilter);
 //   det->RegisterPrimitive(primitive);
 //   primitive = new G4PSNofSecondary("nPositron",j);
 //   primitive->SetFilter(positronFilter);
 //   det->RegisterPrimitive(primitive);






}*/

