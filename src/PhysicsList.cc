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
/// \file eventgenerator/exgps/src/PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"

// particles
//
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ChargeExchangePhysics.hh"

#include "G4ChargeExchangeProcess.hh"
#include "G4ChargeExchange.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
#include <iomanip>   

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4DataQuestionaire.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronPhysicsQGSP_BIC_AllHP.hh"
#include "G4HadronElasticPhysicsPHP.hh"
#include "G4SystemOfUnits.hh"
#include "G4IonPhysicsPHP.hh"
#include "G4WarnPLStatus.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4Neutron.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4_DECLARE_PHYSCONSTR_FACTORY(G4ChargeExchangePhysics);
PhysicsList::PhysicsList() 
: G4VModularPhysicsList(){
  SetVerboseLevel(1);

  defaultCutValue = 7.*mm;

  // EM Physics
  RegisterPhysics( new G4EmStandardPhysics());

  // Synchroton Radiation & GN Physics
  RegisterPhysics( new G4EmExtraPhysics());

  // Decays
  RegisterPhysics( new G4DecayPhysics() );

   // Hadron Elastic scattering
//  RegisterPhysics( new G4HadronElasticPhysics());

   // Hadron Physics
  RegisterPhysics( new G4HadronElasticPhysicsPHP());
  RegisterPhysics(  new G4HadronPhysicsQGSP_BIC_AllHP());
//  RegisterPhysics( new G4HadronPhysicsQGSP_BIC_HP);
  // Stopping Physics
  RegisterPhysics( new G4StoppingPhysics());
//  RegisterPhysics( new G4StepLimiterPhysics());

  // Ion Physics
  RegisterPhysics( new G4IonPhysicsPHP());
 // RegisterPhysics( new G4ChargeExchangePhysics());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::SetCuts()
{
 G4VUserPhysicsList::SetCuts(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
