#include <stdexcept>
#include <iostream>
#include "globals.hh"
#include "CLHEP/Random/Random.h"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "OutputManager.hh"
#include "G4UIExecutive.hh"
#include "G4VisExecutive.hh"

//---------------------------------------------------------------------------

int main(int argc, char** argv)
{
  CLHEP::HepRandom::createInstance();

  unsigned int seed = time(0);
  unsigned int devrandseed = 0;
  FILE *fdrand = fopen("/dev/urandom", "r");
  if( fdrand ){
    fread(&devrandseed, sizeof(int), 1, fdrand);
    seed += devrandseed;
    fclose(fdrand);
  }

  CLHEP::HepRandom::setTheSeed(seed);
  
  //detect interactive mode (if no arguments) and define UI session
  G4UIExecutive* ui = nullptr;
  if (argc == 1) ui = new G4UIExecutive(argc,argv);

  G4RunManager* runManager = new G4RunManager;
  PhysicsList*  phys       = new PhysicsList();
  runManager->SetUserInitialization(phys);

  DetectorConstruction*   detCon     = new DetectorConstruction();
  PrimaryGeneratorAction* pga        = new PrimaryGeneratorAction();
  OutputManager*          outManager = new OutputManager();
  EventAction*            event      = new EventAction( outManager, pga );

  runManager->SetUserInitialization(detCon);
  runManager->SetUserAction(pga);
  runManager->SetUserAction(new SteppingAction(detCon));
  runManager->SetUserAction(event);

  //initialize visualization
  G4VisManager* visManager = nullptr;

  //get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (ui)  {
   //interactive mode
   visManager = new G4VisExecutive(argc,argv);
   visManager->Initialize();
   UImanager->ExecuteMacroFile("macros/vis_beam.mac");
   ui->SessionStart();
   delete ui;
   delete visManager;
  }
  else  {
   //batch mode
  	G4String command = "/control/execute ";
        G4String fileName = argv[1];
	UImanager->ApplyCommand(command+fileName);
	G4String commandr = "/run/beamOn ";
	G4int nev = pga->GetNEvents();
	char snev[50];
	sprintf( snev, "%d",nev );
	UImanager->ApplyCommand(commandr+snev);
  }

  delete outManager;

  return 0;
}

//---------------------------------------------------------------------------
