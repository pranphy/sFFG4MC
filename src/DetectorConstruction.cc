#include "DetectorConstruction.hh"
#include "VirtualDetectorSD.hh"
#include "RealDetectorSD.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SDManager.hh"
#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"

#include "G4VisAttributes.hh"
#include "G4String.hh"
#include "globals.hh"

#include <fstream>

using namespace CLHEP;
using namespace std;

//---------------------------------------------------------------------------

DetectorConstruction::DetectorConstruction()
{
  fNistManager  = G4NistManager::Instance();
  fDetMessenger = new DetectorMessenger(this);

  fNPSAngle        = 15.5 *deg; 
  fNPSDist         = 301.0 *cm;
  fNPSShieldThick  = 1.0 *cm;
  fHCALAngle       = 42.5 *deg;
  fHCALDist        = 427.0 *cm;
  fHCALShieldThick = 5.0 *cm;
  fSCWinThick      = 0.050 *cm;
  fTarLength       = 10.0 *cm;
  fBeamline        = 0;
}

//---------------------------------------------------------------------------

DetectorConstruction::~DetectorConstruction() 
{
  delete fDetMessenger;
}

//---------------------------------------------------------------------------
G4VPhysicalVolume* DetectorConstruction::Construct()
{ 

  G4int SDcount = 1;

  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  
  //---------------------------------------------------------------------------
  // Create experimental hall
  //---------------------------------------------------------------------------
  
  G4Box* expHall_box           = new G4Box("expHall_box",
					   10.75 *m, 10.75 *m, 10.75 *m );
  

  fexpHall_log = new G4LogicalVolume(expHall_box,
				     fNistManager->FindOrBuildMaterial("G4_AIR"),  
				     "expHall_log", 0, 0, 0);
  
  fExpHall                     = new G4PVPlacement(0, G4ThreeVector(),
						   fexpHall_log, "expHall", 0, false, 0);

  //---------------------------------------------------------------------------
  // Create scattering chamber, target and exit beamline
  // Modified from original NPS simulation:
  // https://github.com/gboon18/HallC_NPS
  //---------------------------------------------------------------------------
  
  if( fBeamline )
    BuildBeamline();
  else
    BuildTarget();
  
  //--------------------------------------------------------------------------- 
  // Create "NPS" electron arm
  //--------------------------------------------------------------------------- 
  
  G4double pbwo4_X = 20.*mm; 
  G4double pbwo4_Y = 20.*mm; 
  G4double pbwo4_Z = 200.*mm;
  
  G4Box* pbwo4_solid = new G4Box("pbwo4_solid", 0.5*pbwo4_X, 0.5*pbwo4_Y, 0.5*pbwo4_Z);
  
  G4LogicalVolume* pbwo4_log = new G4LogicalVolume(pbwo4_solid,
						   fNistManager->FindOrBuildMaterial("G4_PbWO4"),
						   "pbwo4_log");
  
  G4double NPS_x, NPS_y, NPS_z, NPS_th, NPS_ph;
  G4double NPS_xprime, NPS_yprime, NPS_zprime;
  char stmp[50];
  
  for( int ix=0 ; ix < fNPSNrow ; ix++ ) {
    for( int iy=0 ; iy < fNPSNcol ; iy++ ) {
    
      sprintf( stmp, "nps%d", SDcount );
      
      NPS_x  = 0.0;
      NPS_y  = -pbwo4_Y + (ix*pbwo4_Y);
      NPS_z  = fNPSDist;
    
      NPS_th = fNPSAngle;
      
      NPS_ph = 0 + ((360./fNPSNcol) * iy) *deg;

      NPS_yprime = NPS_y * std::cos(NPS_th) + NPS_z * std::sin(NPS_th); 
      NPS_zprime = -NPS_y * std::sin(NPS_th) + NPS_z * std::cos(NPS_th); 
      
      NPS_xprime = NPS_x * std::cos(NPS_ph) + NPS_yprime * std::sin(NPS_ph); 
      NPS_yprime = NPS_x * std::sin(NPS_ph) + NPS_yprime * std::cos(NPS_ph); 
      
      G4Transform3D NPS_t3d = G4Translate3D(G4ThreeVector(NPS_xprime, NPS_yprime, NPS_zprime)) 
	* G4RotateZ3D(NPS_ph).inverse() * G4RotateX3D(NPS_th).inverse();
      
      fDetVol[SDcount] = new G4PVPlacement(NPS_t3d, pbwo4_log, stmp, fexpHall_log, false, SDcount);
      
      SDcount++;
    }
  }


  //--------------------------------------------------------------------------- 
  // Create electron arm shields
  //--------------------------------------------------------------------------- 

  G4double NPSshield_X = pbwo4_X; 
  G4double NPSshield_Y = ((fNPSNrow) * pbwo4_Y) + 5.*mm;  
  G4double NPSshield_Z = fNPSShieldThick;
  
  G4Box* NPSshield_solid = new G4Box("NPSshield_solid", 0.5*NPSshield_X, 0.5*NPSshield_Y, 0.5*NPSshield_Z);
  
  G4LogicalVolume* NPSshield_log = new G4LogicalVolume( NPSshield_solid,
						        fNistManager->FindOrBuildMaterial("G4_Pb"),
						        "NPSshield_log");

  G4double NPSshield_x, NPSshield_y, NPSshield_z, NPSshield_th, NPSshield_ph;
  G4double NPSshield_xprime, NPSshield_yprime, NPSshield_zprime;
  
  for( int iy=0 ; iy < fNPSNcol ; iy++ ) {
    
    sprintf( stmp, "npsshield%d", iy );
    
    NPSshield_x  = 0.0;
    NPSshield_y  = -pbwo4_Y;
    NPSshield_z  = fNPSDist - 0.5*pbwo4_Z - NPSshield_Z - 20 *mm;
    
    NPSshield_th = fNPSAngle;
    
    NPSshield_ph = 0 + ((360./fNPSNcol) * iy) *deg;
    
    NPSshield_yprime = NPSshield_z * std::sin(NPSshield_th); 
    NPSshield_zprime = NPSshield_z * std::cos(NPSshield_th); 
    
    NPSshield_xprime = NPSshield_x * std::cos(NPSshield_ph) + NPSshield_yprime * std::sin(NPSshield_ph); 
    NPSshield_yprime = NPSshield_x * std::sin(NPSshield_ph) + NPSshield_yprime * std::cos(NPSshield_ph); 
    
    G4Transform3D NPSshield_t3d = G4Translate3D(G4ThreeVector(NPSshield_xprime, NPSshield_yprime, 0.0)) * G4RotateZ3D(NPSshield_ph).inverse() 
      * G4Translate3D(G4ThreeVector(0.0, 0.0, NPSshield_zprime)) * G4RotateX3D(NPSshield_th).inverse();
    
    new G4PVPlacement(NPSshield_t3d, NPSshield_log, stmp, fexpHall_log, false, 0);
    
  }

  
  //--------------------------------------------------------------------------- 
  // Create "HCAL" proton arm
  //--------------------------------------------------------------------------- 

  // 44 pairs of 10mm scint + 13mm Fe = 1012 mm
  G4int    HCALNpairs = 44;
  G4double feabs_Z    = 13.0*mm;
  G4double scintabs_X = 150.0*mm; 
  G4double scintabs_Y = 150.0*mm; 
  G4double scintabs_Z = 1012.0*mm;
  
  G4Box* scintabs_solid = new G4Box("scintabs_solid", 0.5*scintabs_X, 0.5*scintabs_Y, 0.5*scintabs_Z);
  
  G4Box* feabs_solid    = new G4Box("feabs_solid", 0.5*scintabs_X, 0.5*scintabs_Y, 0.5*feabs_Z);
  
  G4LogicalVolume* scintabs_log = new G4LogicalVolume(scintabs_solid,
						      fNistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"),
						      "scintabs_log");
  
  G4LogicalVolume* feabs_log = new G4LogicalVolume(feabs_solid,
						   fNistManager->FindOrBuildMaterial("G4_Fe"),
						   "feabs_log");
  
  for( int iz=0 ; iz < HCALNpairs ; iz++ ) {
    
    sprintf( stmp, "feabs%d", iz );
    new G4PVPlacement(0, G4ThreeVector( 0., 0., -scintabs_Z/2 + (iz+0.5)*(10*mm + feabs_Z)  ), 
		      feabs_log, stmp, scintabs_log, false, 9999 );   
  }
  
  G4double HCAL_x, HCAL_y, HCAL_z, HCAL_th, HCAL_ph;
  G4double HCAL_xprime, HCAL_yprime, HCAL_zprime;
  
  for( int ix=0 ; ix < fHCALNrow ; ix++ ) {
    for( int iy=0 ; iy < fHCALNcol ; iy++ ) {
      
      sprintf( stmp, "hcal%d", SDcount );
      
      HCAL_x  = 0.0;
      HCAL_y  = -scintabs_Y + (ix*scintabs_Y);
      HCAL_z  = fHCALDist;
      
      HCAL_th = fHCALAngle;
      
      HCAL_ph = 0 + ((360./fHCALNcol) * iy) *deg;
      
      HCAL_yprime = HCAL_y * std::cos(HCAL_th) + HCAL_z * std::sin(HCAL_th); 
      HCAL_zprime = -HCAL_y * std::sin(HCAL_th) + HCAL_z * std::cos(HCAL_th); 
      
      HCAL_xprime = HCAL_x * std::cos(HCAL_ph) + HCAL_yprime * std::sin(HCAL_ph); 
      HCAL_yprime = HCAL_x * std::sin(HCAL_ph) + HCAL_yprime * std::cos(HCAL_ph); 
      
      G4Transform3D HCAL_t3d = G4Translate3D(G4ThreeVector(HCAL_xprime, HCAL_yprime, HCAL_zprime)) 
	* G4RotateZ3D(HCAL_ph).inverse() * G4RotateX3D(HCAL_th).inverse();
      
      fDetVol[SDcount] = new G4PVPlacement(HCAL_t3d, scintabs_log, stmp, fexpHall_log, false, SDcount);
      
      SDcount++;
    }
  }


  //--------------------------------------------------------------------------- 
  // Create hadron arm hodoscope layer
  //--------------------------------------------------------------------------- 

  G4double hodoscint_X = 30.0 *mm;
  G4double hodoscint_Y = 30.0 *mm;
  G4double hodoscint_Z = 100.0*mm;
  
  G4Box* hodoscint_solid = new G4Box("hodoscint_solid", 0.5*hodoscint_X, 0.5*hodoscint_Y, 0.5*hodoscint_Z);
  
  G4LogicalVolume* hodoscint_log = new G4LogicalVolume(hodoscint_solid,
						       fNistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE"),
						       "hodoscint_log");
  
  G4double HODO_x, HODO_y, HODO_z, HODO_th, HODO_ph;
  G4double HODO_xprime, HODO_yprime, HODO_zprime;

  for( int ix=0 ; ix < fHodoNrow ; ix++ ) {
    for( int iy=0 ; iy < fHodoNcol ; iy++ ) {
      
      sprintf( stmp, "hodo%d", SDcount );
      
      HODO_x  = 0.0;
      HODO_y  = -(15*hodoscint_Y)/2 + (ix*hodoscint_Y);
      HODO_z  = fHCALDist - 0.5*scintabs_Z - 0.5*hodoscint_Z - 20 *mm;
      
      HODO_th = fHCALAngle;
      
      HODO_ph = 0 + ((360./fHodoNcol) * iy) *deg;
      
      HODO_yprime = HODO_y * std::cos(HODO_th) + HODO_z * std::sin(HODO_th); 
      HODO_zprime = -HODO_y * std::sin(HODO_th) + HODO_z * std::cos(HODO_th); 
      
      HODO_xprime = HODO_x * std::cos(HODO_ph) + HODO_yprime * std::sin(HODO_ph); 
      HODO_yprime = HODO_x * std::sin(HODO_ph) + HODO_yprime * std::cos(HODO_ph); 
      
      G4Transform3D HODO_t3d = G4Translate3D(G4ThreeVector(HODO_xprime, HODO_yprime, HODO_zprime))
	* G4RotateZ3D(HODO_ph).inverse() * G4RotateX3D(HODO_th).inverse();
      
      fDetVol[SDcount] = new G4PVPlacement(HODO_t3d, hodoscint_log, stmp, fexpHall_log, false, SDcount);
      
      SDcount++;
    }
  }

  //--------------------------------------------------------------------------- 
  // Create hadron arm shield
  //--------------------------------------------------------------------------- 

  G4double HCALshield_X = scintabs_X; 
  G4double HCALshield_Y = ((fHCALNrow) * scintabs_Y) ; 
  G4double HCALshield_Z = fNPSShieldThick;
  
  G4Box* HCALshield_solid = new G4Box("HCALshield_solid", 0.5*HCALshield_X, 0.5*HCALshield_Y, 0.5*HCALshield_Z);
  
  G4LogicalVolume* HCALshield_log = new G4LogicalVolume(HCALshield_solid,
							fNistManager->FindOrBuildMaterial("G4_Pb"),
							"HCALshield_log");
  
  G4double HCALshield_x, HCALshield_y, HCALshield_z, HCALshield_th, HCALshield_ph;
  G4double HCALshield_xprime, HCALshield_yprime, HCALshield_zprime;
  
  for( int iy=0 ; iy < fHCALNcol ; iy++ ) {
    
    sprintf( stmp, "hcalshield%d", iy );
    
    HCALshield_x  = 0.0;
    HCALshield_y  = -scintabs_Y;
    HCALshield_z  = fHCALDist - 0.5*scintabs_Z - hodoscint_Z - HCALshield_Z - 20 *mm;
    
    HCALshield_th = fHCALAngle;
    
    HCALshield_ph = 0 + ((360./fHCALNcol) * iy) *deg;
    
    HCALshield_yprime = HCALshield_z * std::sin(HCALshield_th); 
    HCALshield_zprime = HCALshield_z * std::cos(HCALshield_th); 
    
    HCALshield_xprime = HCALshield_x * std::cos(HCALshield_ph) + HCALshield_yprime * std::sin(HCALshield_ph); 
    HCALshield_yprime = HCALshield_x * std::sin(HCALshield_ph) + HCALshield_yprime * std::cos(HCALshield_ph); 
    
    G4Transform3D HCALshield_t3d = G4Translate3D(G4ThreeVector(HCALshield_xprime, HCALshield_yprime, HCALshield_zprime)) 
      * G4RotateZ3D(HCALshield_ph).inverse() * G4RotateX3D(HCALshield_th).inverse();
    
    new G4PVPlacement(HCALshield_t3d, HCALshield_log, stmp, fexpHall_log, false, 0);
    
  }
    
  //---------------------------------------------------------------------------
  // Set Logical Attributes
  //---------------------------------------------------------------------------

  // Senstive detector
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  
  fVirtualDetectorSD = new VirtualDetectorSD("VirtualDetectorSD", fNSD );
  SDman->AddNewDetector( fVirtualDetectorSD );
  pbwo4_log->SetSensitiveDetector( fVirtualDetectorSD );
  scintabs_log->SetSensitiveDetector( fVirtualDetectorSD );
  hodoscint_log->SetSensitiveDetector( fVirtualDetectorSD );

  fRealDetectorSD = new RealDetectorSD("RealDetectorSD", fNSD );
  SDman->AddNewDetector( fRealDetectorSD );
  pbwo4_log->SetSensitiveDetector( fRealDetectorSD );
  scintabs_log->SetSensitiveDetector( fRealDetectorSD );
  hodoscint_log->SetSensitiveDetector( fRealDetectorSD );

  // Visualisation
  fLogicTarget->SetVisAttributes(G4Colour::Blue());
  pbwo4_log->SetVisAttributes(G4Colour::Red());
  NPSshield_log->SetVisAttributes(G4Colour::Cyan());
  scintabs_log->SetVisAttributes(G4Colour::Yellow());
  hodoscint_log->SetVisAttributes(G4Colour::Green());
  HCALshield_log->SetVisAttributes(G4Colour::Cyan());

  //feabs_log->SetVisAttributes(G4VisAttributes::Invisible);
  //fexpHall_log->SetVisAttributes(G4VisAttributes::Invisible);

  //---------------------------------------------------------------------------

  return fExpHall;

}

//---------------------------------------------------------------------------

void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}

//---------------------------------------------------------------------------

void DetectorConstruction::BuildBeamline()
{
  
  G4Material* VacuumMaterial                = fNistManager->FindOrBuildMaterial("G4_Galactic");
  G4Material* ChamberMaterial               = fNistManager->FindOrBuildMaterial("G4_Al");
  G4Material* WindowMaterial                = fNistManager->FindOrBuildMaterial("G4_Al");
  G4Material* TargetMaterial                = fNistManager->FindOrBuildMaterial("G4_lH2");
  G4Material* TargetCellMaterial            = fNistManager->FindOrBuildMaterial("G4_Al");
  G4Material* BeampipeMaterial              = fNistManager->FindOrBuildMaterial("G4_Al");
  
  const G4double     inch = 2.54*cm;
  
  G4double ChamberOuterRadius   = 22.5*inch;
  G4double ChamberInnerRadius   = 20.5*inch;

  G4double ChamberHight         = 24.25*2*inch;//20180313 to be modified
  G4double WindowHight          = 19*inch;
  G4double WindowSubHight       = 15*inch;
  G4double WindowFrameHight     = 20*inch;//20130313 guessed
  G4double WindowClampHight     = 20*inch;

  G4double WindowFrameThickness = 1.25*inch;//20130313 guessed
  G4double WindowClampThickness = 0.750*inch;
  G4double WindowThickness      = fSCWinThick; //0.020*inch;

//   G4double ChamberHight         = 24.25*2*inch;//20180313 to be modified
//   G4double WindowHight          = 19*inch;
//   G4double WindowSubHight       = 15*inch;
//   G4double WindowFrameThickness = 1.25*inch;//20130313 guessed
//   G4double WindowFrameHight     = 20*inch;//20130313 guessed
//   G4double WindowClampThickness = 0.750*inch;
//   G4double WindowClampHight     = 20*inch;
//   G4double WindowThickness      = 0.020*inch;
  
  G4double TargetLength          = fTarLength; //100.*mm;//(129.759 - 29.759 + 50.)*mm; 
  G4double TargetRadius          = 0.5*50.*mm;//20.179*mm;
  G4double TargetCellLength      = (0.125 + 150.)*mm;//(5. + 129.759 - 29.759 + 50.)*mm; //(129.887 - 29.759 + 50.)*mm; 
  G4double TargetCellRadius      = (0.125 + 0.5*50.)*mm;//(20.179 + 5.)*mm;//20.32*mm;
  G4double TargetWindowThickness = 0.125*mm;//5.*mm;//0.128*mm;
  
  // G4double WindowInnerJoint1InnerRadius   = 0.5*0.846*inch;
  G4double WindowInnerJoint1OuterRadius   = 0.5*1.469*inch;
  G4double WindowInnerJoint1Thickness     = 0.5*inch;
  G4double WindowInnerJoint2InnerRadius   = 0.5*1.068*inch;
  G4double WindowInnerJoint2OuterRadius   = 0.5*1.50*inch;
  G4double WindowInnerJoint2Thickness     = 0.109*inch;
  // G4double WindowOuterJoint1InnerRadius   = 0.5*1.7*inch;//no need
  //G4double WindowOuterJoint1OuterRadius   = 0.5*2.38*inch;//no need
  //G4double WindowOuterJoint1Thickness     = 0.016*inch;//no need
  G4double WindowOuterJoint2InnerRadius   = 0.5*0.68*inch;
  G4double WindowOuterJoint2OuterRadius   = 0.5*1.062*inch;
  G4double WindowOuterJoint2Thickness     = (0.62 + 1.562 - 0.137 - 0.06)*inch;
  G4double WindowOuterJoint2_1InnerRadius = WindowOuterJoint2OuterRadius;
  G4double WindowOuterJoint2_1OuterRadius = 0.5*1.5*inch;
  G4double WindowOuterJoint2_1Thickness   = 0.19*inch;
  G4double WindowOuterJoint2_1Position    = 0.62*inch;
  //G4double WindowOuterJoint2_2InnerRadius = WindowOuterJoint2InnerRadius;
  G4double WindowOuterJoint2_2dx          = 1.12*inch;
  G4double WindowOuterJoint2_2dy          = 1.12*inch ;
  G4double WindowOuterJoint2_2dz          = 0.137*inch;
  G4double WindowOuterJoint2_3InnerRadius = 0.5*0.68*inch;//=WindowOuterJoint2InnerRadius;
  G4double WindowOuterJoint2_3OuterRadius = 0.5*0.738*inch;
  G4double WindowOuterJoint2_3Thickness   = 0.06*inch;
  
  G4double Beampipe1Innerdx1    = 18.915*mm;
  G4double Beampipe1Innerdx2    = 63.4*mm;
  G4double Beampipe1Innerdy1    = Beampipe1Innerdx1;
  G4double Beampipe1Innerdy2    = Beampipe1Innerdx2;
  G4double Beampipe1Outerdx1    = 25.265*mm;
  G4double Beampipe1Outerdx2    = 69.749*mm;
  G4double Beampipe1Outerdy1    = Beampipe1Outerdx1;
  G4double Beampipe1Outerdy2    = Beampipe1Outerdx2;
  G4double Beampipe1Length      = 1685.925*mm;
  G4double Beampipe2OuterRadius = 0.5*168.275*mm;
  G4double Beampipe2InnerRadius = 0.5*154.051*mm;
  G4double Beampipe2Length      = 129.633*inch;//20180417 changedrom 2620.138*mm;<-this is the actual length of the pipe.
  //however, there is a gap between beampipe1 and 2. In order toill the gap. The length of the beampipe2 currently is longer.
  G4double Beampipe2FrontCellThickness = Beampipe2OuterRadius -Beampipe2InnerRadius; //made up
  G4double Beampipe3OuterRadius = 0.5*273.05*mm;
  G4double Beampipe3InnerRadius = 0.5*254.508*mm;
  G4double Beampipe3Length      = 5102.225*mm;
  G4double Beampipe3FrontCellThickness = Beampipe3OuterRadius -Beampipe3InnerRadius; //made up
  
  //---------------------------------------------------------------------------
  //Scattering Chamber
  //---------------------------------------------------------------------------

  G4RotationMatrix *xChambRot = new G4RotationMatrix;  
  xChambRot->rotateX(90*degree);                     
  
  G4Tubs* sChamberOuter = new G4Tubs("ChamberOuter_sol",
				     ChamberInnerRadius,
				     ChamberOuterRadius,
				     0.5*ChamberHight,
				     0.,
				     twopi);
  
  G4double WindowStartTheta = (3.+10.+90.+180.)*pi/180.;//window beginning position
  G4double WindowDeltaTheta = 124.*pi/180.;//window size 127deg
  G4double WindowFrameStartTheta = WindowStartTheta-(4.5)*pi/180.;//to match the center of the window
  G4double WindowFrameDeltaTheta = WindowDeltaTheta+9*pi/180.;//total 133 deg.

  G4Tubs* sWindowSub = new G4Tubs("WindowSub_sol",
				  ChamberInnerRadius-1*cm,
				  ChamberOuterRadius+1*cm,
				  0.5*WindowSubHight,
				  WindowStartTheta,			    
				  WindowDeltaTheta);
  
  G4RotationMatrix *zWindowRot = new G4RotationMatrix;  
  zWindowRot->rotateZ(90*degree);                     
  G4SubtractionSolid* sChamber_sub_Window = new G4SubtractionSolid("Chamber_sub_Window",
								   sChamberOuter,
								   sWindowSub,
								   zWindowRot,
								   G4ThreeVector());

  G4LogicalVolume* LogicChamber = new G4LogicalVolume(sChamber_sub_Window, ChamberMaterial, "ChamberOuter_log");

  new G4PVPlacement(xChambRot, G4ThreeVector(), LogicChamber, "ChamberOuter_pos", fexpHall_log, false, 0);

  //---------------------------------------------------------------------------

  G4Tubs* sWindowFrame_before_sub = new G4Tubs("WindowFrame_before_sub_sol",
					       ChamberOuterRadius,
					       ChamberOuterRadius + WindowFrameThickness,//1.25inch is the thickness
					       0.5*WindowFrameHight,//20inch is the hight
					       WindowFrameStartTheta,//to match the center of the window
					       WindowFrameDeltaTheta);//total 133 deg.

  G4Tubs* sWindowFrame_sub = new G4Tubs("WindowFrame_sub_sol",
					ChamberOuterRadius - 1*cm,
					ChamberOuterRadius + WindowFrameThickness + 1*cm,//1.25inch is the thickness of the window frame
					0.5*WindowSubHight,
					WindowStartTheta,
					WindowDeltaTheta);
  
  G4RotationMatrix *pseudoWindowRot = new G4RotationMatrix;  
  pseudoWindowRot->rotateZ(0*degree);                     
  
  G4SubtractionSolid* sWindowFrame = new G4SubtractionSolid("WindowFrame_sol",
							    sWindowFrame_before_sub,
							    sWindowFrame_sub,
							    pseudoWindowRot,
							    G4ThreeVector());

  G4LogicalVolume* LogicWindowFrame = new G4LogicalVolume(sWindowFrame, ChamberMaterial, "WindowFrame_log");

  G4RotationMatrix *zxWindowRot = new G4RotationMatrix(0,-90*degree,-90*degree);

  new G4PVPlacement(zxWindowRot, G4ThreeVector(), LogicWindowFrame, "WindowFrame_pos", fexpHall_log, false, 0);

  //---------------------------------------------------------------------------

  G4Tubs* sWindowClamp_before_sub = new G4Tubs("WindowClamp_before_sub_sol",
					       ChamberOuterRadius + WindowFrameThickness + WindowThickness,//1.25inch is the thickness of the window frame
					       ChamberOuterRadius + WindowFrameThickness + WindowThickness + WindowClampThickness,//0.75inch is the thickness of the window clamp
					       0.5*WindowClampHight,//20inch is the hight
					       WindowFrameStartTheta,//to match the center of the window
					       WindowFrameDeltaTheta);//total 133 deg.
  
  G4Tubs* sWindowClamp_sub = new G4Tubs("WindowClamp_sub_sol",
					ChamberOuterRadius + WindowFrameThickness + WindowThickness -1*cm,
					ChamberOuterRadius + WindowFrameThickness + WindowThickness + WindowClampThickness +1*cm,//1.25inch is the thickness
					0.5*WindowSubHight,
					WindowStartTheta,
					WindowDeltaTheta);

  G4SubtractionSolid* sWindowClamp = new G4SubtractionSolid("WindowClamp_sol",
  							    sWindowClamp_before_sub,
  							    sWindowClamp_sub,
  							    pseudoWindowRot,
  							    G4ThreeVector());
  
  G4LogicalVolume* LogicWindowClamp = new G4LogicalVolume(sWindowClamp, ChamberMaterial, "WindowClamp_log");

  new G4PVPlacement(zxWindowRot, G4ThreeVector(), LogicWindowClamp, "WindowClamp_pos", fexpHall_log, false, 0 );
  
  //---------------------------------------------------------------------------
  // Vacuum inside the chamber
  //---------------------------------------------------------------------------

  G4Tubs* sChamberInner = new G4Tubs("ChamberInner_sol",
				     0.,
				     ChamberInnerRadius,
				     0.5*ChamberHight,
				     0.,
				     twopi);

  G4LogicalVolume* LogicInnerChamber = new G4LogicalVolume(sChamberInner, VacuumMaterial, "ChamberInner_log");

  new G4PVPlacement(xChambRot, G4ThreeVector(), LogicInnerChamber, "ChamberInner_pos", fexpHall_log, false, 0 );

  //---------------------------------------------------------------------------
  // Target Cell & Target
  //---------------------------------------------------------------------------

  //---------------------------------------------------------------------------
  // Target Cell & Target
  //---------------------------------------------------------------------------

  G4RotationMatrix *xTargetRot = new G4RotationMatrix;  
  xTargetRot->rotateX(-90*degree);                     

  G4Tubs* sTargetCell = new G4Tubs("TargetCell_sol",           
				   0., TargetCellRadius, (0.5*TargetLength)+TargetWindowThickness, 0.,twopi); 

  G4LogicalVolume*   LogicTargetCell = new G4LogicalVolume(sTargetCell, TargetCellMaterial, "TargetCell_log");   

  new G4PVPlacement(xTargetRot, G4ThreeVector(0., -43.9*cm, 0.), 
		    LogicTargetCell, "TargetCell_pos",  LogicInnerChamber, false, 0 );  

  //---------------------------------------------------------------------------

  G4Tubs* sTarget = new G4Tubs("Target_sol",         
			       0., TargetRadius, 0.5*TargetLength, 0.,twopi); 

  fLogicTarget = new G4LogicalVolume(sTarget, TargetMaterial, "Target_log");  
  
  new G4PVPlacement(0, G4ThreeVector(0., 0., 0.5*TargetWindowThickness), fLogicTarget, "Target_pos", LogicTargetCell, false, 0 ); 

  //---------------------------------------------------------------------------
  // beamline
  // window joint
  //---------------------------------------------------------------------------

//   G4Tubs* sWindowInnerJoint1 = new G4Tubs("WindowInnerJoint1_sol",
// 					  WindowOuterJoint2OuterRadius,
// 					  WindowInnerJoint1OuterRadius,
// 					  0.5*WindowInnerJoint1Thickness,
// 					  0.,
// 					  twopi);
  
//   G4LogicalVolume* LogicWindowInnerJoint1 = new G4LogicalVolume(sWindowInnerJoint1, BeampipeMaterial, "WindowInnerJoint1_log");

//   new G4PVPlacement(0,
//   		    G4ThreeVector(0,0,ChamberOuterRadius + WindowFrameThickness - WindowThickness -WindowInnerJoint2Thickness - 0.5*WindowInnerJoint1Thickness),
//   		    LogicWindowInnerJoint1,
//   		    "WindowInnerJoint1_pos",
//   		    fexpHall_log,
//   		    false,
//   		    0 );

//   //---------------------------------------------------------------------------

//   G4Tubs* sWindowInnerJoint2 = new G4Tubs("WindowInnerJoint2_sol",
// 					  WindowInnerJoint2InnerRadius,
// 					  WindowInnerJoint2OuterRadius,
// 					  0.5*WindowInnerJoint2Thickness,
// 					  0.,
// 					  twopi);

//   G4LogicalVolume* LogicWindowInnerJoint2 = new G4LogicalVolume(sWindowInnerJoint2, BeampipeMaterial, "WindowInnerJoint2_log");

//   new G4PVPlacement(0,
//   		    G4ThreeVector(0,0,ChamberOuterRadius + WindowFrameThickness - WindowThickness - 0.5*WindowInnerJoint2Thickness),
//   		    LogicWindowInnerJoint2,
//   		    "WindowInnerJoint2_pos",
//   		    fexpHall_log,
//   		    false,
//   		    0 );
  		    
//   //---------------------------------------------------------------------------

//   G4Tubs* sWindowOuterJoint2 = new G4Tubs("WindowOuterJoint2_sol",
// 					  0.,
// 					  WindowOuterJoint2OuterRadius,
// 					  0.5*WindowOuterJoint2Thickness,
// 					  0.,
// 					  twopi);

//   G4LogicalVolume* LogicWindowOuterJoint2 = new G4LogicalVolume(sWindowOuterJoint2,
// 								BeampipeMaterial,
// 								"WindowOuterJoint2_log");

//   new G4PVPlacement(0,
//   		    G4ThreeVector(0,0,ChamberOuterRadius + WindowFrameThickness +0.5*WindowOuterJoint2Thickness - WindowOuterJoint2_1Position),
//   		    LogicWindowOuterJoint2,
//   		    "WindowOuterJoint2_pos",
//   		    fexpHall_log,
//   		    false,
//   		    0 );
  		    

//   //---------------------------------------------------------------------------
//   //vacuum inside the OuterJoint2
//   //---------------------------------------------------------------------------

//   G4Tubs* sWindowOuterJoint2Vacuum = new G4Tubs("WindowOuterJoint2Vacuum_sol",
// 						0.,
// 						WindowOuterJoint2InnerRadius,
// 						0.5*WindowOuterJoint2Thickness,
// 						0.,
// 						twopi);
  
//   G4LogicalVolume* LogicWindowOuterJoint2Vacuum = new G4LogicalVolume(sWindowOuterJoint2Vacuum,
// 								       VacuumMaterial,
// 								       "WindowOuterJoint2Vacuum_log");

//   new G4PVPlacement(0,
//   		    G4ThreeVector(0,0,0),
//   		    LogicWindowOuterJoint2Vacuum,
//   		    "WindowOuterJoint2Vacuum_pos",
//   		    LogicWindowOuterJoint2,
//   		    false,
//   		    0 );
  		    

//   //---------------------------------------------------------------------------

//   G4Tubs* sWindowOuterJoint2_1 = new G4Tubs("WindowOuterJoint2_1_sol",
// 					    WindowOuterJoint2_1InnerRadius,
// 					    WindowOuterJoint2_1OuterRadius,
// 					    0.5*WindowOuterJoint2_1Thickness,
// 					    0.,
// 					    twopi);
  
//   G4LogicalVolume* LogicWindowOuterJoint2_1 = new G4LogicalVolume(sWindowOuterJoint2_1,
// 								  BeampipeMaterial,
// 								  "WindowOuterJoint2_1_log");

//   new G4PVPlacement(0,
//   		    G4ThreeVector(0,0,ChamberOuterRadius + WindowFrameThickness + 0.5*WindowOuterJoint2_1Thickness),
//   		    LogicWindowOuterJoint2_1,
//   		    "WindowOuterJoint2_1_pos",
//   		    fexpHall_log,
//   		    false,
//   		    0 );
  		    
//   //---------------------------------------------------------------------------

//   G4RotationMatrix *zPipeRot = new G4RotationMatrix;  
//   zPipeRot->rotateZ(45*degree);                     
  
//   G4Box* sWindowOuterJoint2_2 = new G4Box("WindowOuterJoint2_2_sol",
// 					  0.5*WindowOuterJoint2_2dx,
// 					  0.5*WindowOuterJoint2_2dy,
// 					  0.5*WindowOuterJoint2_2dz);

//   G4LogicalVolume* LogicWindowOuterJoint2_2 = new G4LogicalVolume(sWindowOuterJoint2_2,
// 								  BeampipeMaterial,
// 								  "WindowOuterJoint2_2_log");

//   new G4PVPlacement(zPipeRot,
//   		    G4ThreeVector(0,0,ChamberOuterRadius + WindowFrameThickness + WindowOuterJoint2Thickness - WindowOuterJoint2_1Position + 0.5*WindowOuterJoint2_2dz),
//   		    LogicWindowOuterJoint2_2,
//   		    "WindowOuterJoint2_2_pos",
//   		    fexpHall_log,
//   		    false,
//   		    0 );
  		    

//   //---------------------------------------------------------------------------
//   //vacuum inside the OuterJoint2_2
//   //---------------------------------------------------------------------------

//   G4Tubs* sWindowOuterJoint2_2Vacuum = new G4Tubs("WindowOuterJoint2_2Vacuum_sol",
// 						  0.,
// 						  WindowOuterJoint2InnerRadius,
// 						  0.5*WindowOuterJoint2_2dz,
// 						  0.,
// 						  twopi);
  
//   G4LogicalVolume* LogicWindowOuterJoint2_2Vacuum = new G4LogicalVolume(sWindowOuterJoint2_2Vacuum,
// 									 VacuumMaterial,
// 									 "WindowOuterJoint2_2Vacuum_log");

//   new G4PVPlacement(0,
//   		    G4ThreeVector(0,0,0),
//   		    LogicWindowOuterJoint2_2Vacuum,
//   		    "WindowOuterJoint2_2Vacuum_pos",
//   		    LogicWindowOuterJoint2_2,
//   		    false,
//   		    0 );
  		    
//   //---------------------------------------------------------------------------

//   G4Tubs* sWindowOuterJoint2_3 = new G4Tubs("WindowOuterJoint2_3_sol",
// 					    0.,
// 					    WindowOuterJoint2_3OuterRadius,
// 					    0.5*WindowOuterJoint2_3Thickness,
// 					    0.,
// 					    twopi);

//   G4LogicalVolume* LogicWindowOuterJoint2_3 = new G4LogicalVolume(sWindowOuterJoint2_3,
// 								   BeampipeMaterial,
// 								   "WindowOuterJoint2_3_log");

//   new G4PVPlacement(0,
//   		    G4ThreeVector(0,0,ChamberOuterRadius + WindowFrameThickness + WindowOuterJoint2Thickness - WindowOuterJoint2_1Position 
// 				  + WindowOuterJoint2_2dz + 0.5*WindowOuterJoint2_3Thickness),
//   		    LogicWindowOuterJoint2_3,
//   		    "WindowOuterJoint2_3_pos",
//   		    fexpHall_log,
//   		    false,
//   		    0 );
  		    

//   //---------------------------------------------------------------------------
//   // vacuum inside the OuterJoint2_3
//   //---------------------------------------------------------------------------

//   G4Tubs* sWindowOuterJoint2_3Vacuum = new G4Tubs("WindowOuterJoint2_3Vacuum_sol",
// 						  0.,
// 						  WindowOuterJoint2_3InnerRadius,
// 						  0.5*WindowOuterJoint2_3Thickness,
// 						  0.,
// 						  twopi);

//   G4LogicalVolume* LogicWindowOuterJoint2_3Vacuum = new G4LogicalVolume(sWindowOuterJoint2_3Vacuum,
// 									 VacuumMaterial,
// 									 "WindowOuterJoint2_3Vacuum_log");
  
//   new G4PVPlacement(0,
//   		    G4ThreeVector(0,0,0),
//   		    LogicWindowOuterJoint2_3Vacuum,
//   		    "WindowOuterJoint2_3Vacuum_pos",
//   		    LogicWindowOuterJoint2_3,
//   		    false,
//   		    0 );
  		    
//   //---------------------------------------------------------------------------
  
//   G4Trd* sBeampipe1 = new G4Trd("Beampipe1_sol",
// 				0.5*Beampipe1Outerdx1,
// 				0.5*Beampipe1Outerdx2,
// 				0.5*Beampipe1Outerdy1,
// 				0.5*Beampipe1Outerdy2,
// 				0.5*Beampipe1Length);
  
//   G4LogicalVolume* LogicBeampipe1 = new G4LogicalVolume(sBeampipe1,
// 							BeampipeMaterial,
// 							"Beampipe1_log");
  
//   new G4PVPlacement(zPipeRot,
//   		    G4ThreeVector(0,0,ChamberOuterRadius + WindowFrameThickness + WindowOuterJoint2Thickness - WindowOuterJoint2_1Position 
// 				  + WindowOuterJoint2_2dz  + 0.5*Beampipe1Length),
//   		    LogicBeampipe1,
//   		    "Beampipe1_pos",
//   		    fexpHall_log,
//   		    false,
//   		    0 );
  		    
  
//   //---------------------------------------------------------------------------
  
//   G4Trd* sBeampipe1Vacuum = new G4Trd("Beampipe1Vacuum_sol",
// 				      0.5*Beampipe1Innerdx1,
// 				      0.5*Beampipe1Innerdx2,
// 				      0.5*Beampipe1Innerdy1,
// 				      0.5*Beampipe1Innerdy2,
// 				      0.5*Beampipe1Length);
  
//   G4LogicalVolume* LogicBeampipe1Vacuum = new G4LogicalVolume(sBeampipe1Vacuum,
// 							       VacuumMaterial,
// 							       "Beampipe1Vacuum_log");
  
//   new G4PVPlacement(0,
//   		    G4ThreeVector(0,0,0),
//   		    LogicBeampipe1Vacuum,
//   		    "Beampipe1Vacuum_pos",
//   		    LogicBeampipe1,
//   		    false,
//   		    0 );
  		    

//   //---------------------------------------------------------------------------
//   //Beampipe2
//   //---------------------------------------------------------------------------

//   G4Tubs* sBeampipe2 = new G4Tubs("Beampipe_sol",
// 				  0.,
// 				  Beampipe2OuterRadius,
// 				  0.5*Beampipe2Length,
// 				  0.,
// 				  twopi);

//   G4LogicalVolume* LogicBeampipe2 = new G4LogicalVolume(sBeampipe2,
// 							BeampipeMaterial,
// 							"Beampipe2_log");

//   new G4PVPlacement(0,
//   		    G4ThreeVector(0,0,ChamberOuterRadius + WindowFrameThickness + WindowOuterJoint2Thickness - WindowOuterJoint2_1Position 
// 				  + WindowOuterJoint2_2dz + Beampipe1Length + 0.5*Beampipe2Length),
//   		    LogicBeampipe2,
//   		    "Beampipe2_pos",
//   		    fexpHall_log,
//   		    false,
//   		    0 );
  		    
//   //---------------------------------------------------------------------------
//   // vacuum inside the Beampipe2
//   //---------------------------------------------------------------------------

//   G4Tubs* sBeampipe2Vacuum = new G4Tubs("Beampipe2Vacuum_sol",
// 					0.,
// 					Beampipe2InnerRadius,
// 					0.5*Beampipe2Length,
// 					0.,
// 					twopi);

//   G4LogicalVolume* LogicBeampipe2Vacuum = new G4LogicalVolume(sBeampipe2Vacuum,
// 							       VacuumMaterial,
// 							       "Beampipe2Vacuum_log");

//   new G4PVPlacement(0,
//   		    G4ThreeVector(0,0,0),
//   		    LogicBeampipe2Vacuum,
//   		    "Beampipe2Vacuum_pos",
//   		    LogicBeampipe2,
//   		    false,
//   		    0 );
  		    

//   //---------------------------------------------------------------------------
//   // Beampipe2 front cover
//   //---------------------------------------------------------------------------

//   G4Tubs* sBeampipe2FrontCell_before_sub = new G4Tubs("Beampipe2FrontCell_before_sub_sol",
// 						      0.,
// 						      Beampipe2OuterRadius,
// 						      0.5*Beampipe2FrontCellThickness,
// 						      0.,
// 						      twopi);
  
//   G4SubtractionSolid* sBeampipe2FrontCell = new G4SubtractionSolid("Beampipe2FrontCell_sol",
// 								   sBeampipe2FrontCell_before_sub,
// 								   sBeampipe1,
// 								   zPipeRot,
// 								   G4ThreeVector(0, 0, -0.5*Beampipe1Length + 0.5*Beampipe2FrontCellThickness));

//   G4LogicalVolume* LogicBeampipe2FrontCell = new G4LogicalVolume(sBeampipe2FrontCell,
// 								 BeampipeMaterial,
// 								 "Beampipe2FrontCell_log");

//   new G4PVPlacement(0,
//   		    G4ThreeVector(0,0,ChamberOuterRadius + WindowFrameThickness + WindowOuterJoint2Thickness - WindowOuterJoint2_1Position 
// 				  + WindowOuterJoint2_2dz + Beampipe1Length - 0.5*Beampipe2FrontCellThickness),
//   		    LogicBeampipe2FrontCell,
//   		    "Beampipe2FrontCell_pos",
//   		    fexpHall_log,
//   		    false,
//   		    0 );
  		    
//   //---------------------------------------------------------------------------
//   // Beampipe3
//   //---------------------------------------------------------------------------

//   G4Tubs* sBeampipe3 = new G4Tubs("Beampipe_sol",
// 				  0.,
// 				  Beampipe3OuterRadius,
// 				  0.5*Beampipe3Length,
// 				  0.,
// 				  twopi);
  
//   G4LogicalVolume* LogicBeampipe3 = new G4LogicalVolume(sBeampipe3,
// 							 BeampipeMaterial,
// 							 "Beampipe3_log");

//   new G4PVPlacement(0,
//   		    G4ThreeVector(0,0,ChamberOuterRadius + WindowFrameThickness + WindowOuterJoint2Thickness - WindowOuterJoint2_1Position 
// 				  + WindowOuterJoint2_2dz  + Beampipe1Length + Beampipe2Length + 0.5*Beampipe3Length),
//   		    LogicBeampipe3,
//   		    "Beampipe3_pos",
//   		    fexpHall_log,
//   		    false,
//   		    0 );
  		    
//   //---------------------------------------------------------------------------
//   //vacuum inside the Beampipe3
//   //---------------------------------------------------------------------------

//   G4Tubs* sBeampipe3Vacuum = new G4Tubs("Beampipe3Vacuum_sol",
// 					0.,
// 					Beampipe3InnerRadius,
// 					0.5*Beampipe3Length,
// 					0.,
// 					twopi);

//   G4LogicalVolume* LogicBeampipe3Vacuum = new G4LogicalVolume(sBeampipe3Vacuum,
// 							       VacuumMaterial,
// 							       "Beampipe3Vacuum_log");

//   new G4PVPlacement(0,
//   		    G4ThreeVector(0,0,0),
//   		    LogicBeampipe3Vacuum,
//   		    "Beampipe3Vacuum_pos",
//   		    LogicBeampipe3,
//   		    false,
//   		    0 );
  		    
//   //---------------------------------------------------------------------------
  
//   G4Tubs* sBeampipe3FrontCell = new G4Tubs("Beampipe3FrontCell_before_sub_sol",
// 					   Beampipe2OuterRadius,
// 					   Beampipe3OuterRadius,
// 					   0.5*Beampipe3FrontCellThickness,
// 					   0.,
// 					   twopi);
  
//   G4LogicalVolume* LogicBeampipe3FrontCell = new G4LogicalVolume(sBeampipe3FrontCell,
// 								  BeampipeMaterial,
// 								  "Beampipe3FrontCell_log");
  
//   new G4PVPlacement(0,
//   		    G4ThreeVector(0,0,ChamberOuterRadius + WindowFrameThickness + WindowOuterJoint2Thickness - WindowOuterJoint2_1Position 
// 				  + WindowOuterJoint2_2dz + Beampipe1Length + Beampipe2Length - 0.5*Beampipe3FrontCellThickness),
//   		    LogicBeampipe3FrontCell,
//   		    "Beampipe3FrontCell_pos",
//   		    fexpHall_log,
//   		    false,
//   		    0 );
  		    
  //---------------------------------------------------------------------------
  // Al window for the chamber
  //---------------------------------------------------------------------------

//   G4Tubs* sWindow = new G4Tubs("Window_sol",
// 			       ChamberOuterRadius + WindowFrameThickness,//1.25inch is the thickness of the window frame
// 			       ChamberOuterRadius + WindowFrameThickness + WindowThickness,
// 			       0.5*WindowHight,//actual size is 19*inch however, they are covered by the window clamp
// 			       WindowFrameStartTheta,//to match the center of the window
// 			       WindowFrameDeltaTheta);//total 133 deg.

//   G4RotationMatrix *yRot = new G4RotationMatrix; 
//   yRot->rotateY(-90*degree);                    

//   G4SubtractionSolid* sWindow_sub_InnerJoint2 = new G4SubtractionSolid("Window_sub_InnerJoint2",
//   								       sWindow,
//   								       sWindowInnerJoint2,
//   								       yRot,
//   								       G4ThreeVector(ChamberOuterRadius + WindowFrameThickness 
// 										     - WindowThickness - 0.5*WindowInnerJoint2Thickness, 0, 0));

//   G4SubtractionSolid* sWindow_sub_InnerJoint2_sub_OuterJoint2 = new G4SubtractionSolid("Window_sub_InnerJoint2_sub_OutherJoint2",
//   										       sWindow_sub_InnerJoint2,
//   										       sWindowOuterJoint2,
//   										       yRot,
//   										       G4ThreeVector(ChamberOuterRadius + WindowFrameThickness 
// 												     +0.5*WindowOuterJoint2Thickness - WindowOuterJoint2_1Position, 0, 0));

//   G4SubtractionSolid* sWindow_sub_InnerJoint2_sub_OuterJoint2_sub_OuterJoint2_1 = new G4SubtractionSolid("Window_sub_InnerJoint2_sub_OutherJoint2",
// 													 sWindow_sub_InnerJoint2_sub_OuterJoint2,
// 													 sWindowOuterJoint2_1,
// 													 yRot,
// 													 G4ThreeVector(ChamberOuterRadius + WindowFrameThickness 
// 														       + 0.5*WindowOuterJoint2_1Thickness, 0, 0));

//   G4LogicalVolume*  LogicChamberWindow = new G4LogicalVolume(sWindow_sub_InnerJoint2_sub_OuterJoint2_sub_OuterJoint2_1,
// 							     WindowMaterial,
// 							     "Window_log");
  
//   new G4PVPlacement(zxWindowRot,
//   		    G4ThreeVector(),
//   		    LogicChamberWindow,
//   		    "Window_pos",
//   		    fexpHall_log,
//   		    false,
//   		    0 );
  		    

//   //---------------------------------------------------------------------------
//   //extra detailed vaccum in the chamber
//   //---------------------------------------------------------------------------

//   G4Tubs* sChamberWindowVacuum = new G4Tubs("ChamberWindowVacuum_sol",
// 					    ChamberInnerRadius,
// 					    ChamberOuterRadius+WindowFrameThickness,
// 					    0.5*WindowSubHight,
// 					    WindowStartTheta,
// 					    WindowDeltaTheta);
  
//   G4SubtractionSolid* sChamberWindowVacuum_sub_InnerJoint2 = new G4SubtractionSolid("ChamberWindowVacuum_sub_InnerJoint2",
// 										    sChamberWindowVacuum,
// 										    sWindowInnerJoint2,
// 										    yRot,
// 										    G4ThreeVector(ChamberOuterRadius + WindowFrameThickness - WindowThickness 
// 												  - 0.5*WindowInnerJoint2Thickness, 0, 0));

//   G4SubtractionSolid* sChamberWindowVacuum_sub_InnerJoint2_sub_OuterJoint2 = new G4SubtractionSolid("ChamberWindowVacuum_sub_InnerJoint2_sub_OutherJoint2",
// 												    sChamberWindowVacuum_sub_InnerJoint2,
// 												    sWindowOuterJoint2,
// 												    yRot,
// 												    G4ThreeVector(ChamberOuterRadius + WindowFrameThickness 
// 														  + 0.5*WindowOuterJoint2Thickness - WindowOuterJoint2_1Position, 0, 0));

//   G4SubtractionSolid* sChamberWindowVacuum_sub_InnerJoint2_sub_OuterJoint2_sub_InnerJoint1 = new G4SubtractionSolid("ChamberWindowVacuum_sub_InnerJoint2_sub_OutherJoint2_sub_InnerJoint1",
// 														    sChamberWindowVacuum_sub_InnerJoint2_sub_OuterJoint2,
// 														    sWindowInnerJoint1,
// 														    yRot,
// 														    G4ThreeVector(ChamberOuterRadius + WindowFrameThickness - WindowThickness 
// 																  - WindowInnerJoint2Thickness - 0.5*WindowInnerJoint1Thickness, 0, 0));



//   G4LogicalVolume*  LogicInnerChamber2 = new G4LogicalVolume(sChamberWindowVacuum_sub_InnerJoint2_sub_OuterJoint2_sub_InnerJoint1,//sChamberWindowVacuum,
// 							     VacuumMaterial,
// 							     "ChamberWindowVacuum_log");
  
//   new G4PVPlacement(zxWindowRot,
//   		    G4ThreeVector(),
//   		    LogicInnerChamber2,
//   		    "ChamberWindowVacuum_pos",
//   		    fexpHall_log,
//   		    false,
//   		    0 );
  
}

//---------------------------------------------------------------------------

void DetectorConstruction::BuildTarget()
{
  
  G4Material* VacuumMaterial                = fNistManager->FindOrBuildMaterial("G4_Galactic");
  G4Material* TargetMaterial                = fNistManager->FindOrBuildMaterial("G4_lH2");
  G4Material* TargetCellMaterial            = fNistManager->FindOrBuildMaterial("G4_Al");
  
  const G4double     inch = 2.54*cm;

  G4double ChamberInnerRadius   = 20.5*inch;
  G4double ChamberHight         = 3*24.25*2*inch;//20180313 to be modified

  G4double BeamlineInnerRadius   = 200. *mm;
  G4double BeamlineHight         = 1300. *cm;
  
  G4double TargetLength          = fTarLength; //100.*mm;//(129.759 - 29.759 + 50.)*mm; 
  G4double TargetRadius          = 0.5*50.*mm;//20.179*mm;
  G4double TargetCellLength      = (0.125 + 150.)*mm;//(5. + 129.759 - 29.759 + 50.)*mm; //(129.887 - 29.759 + 50.)*mm; 
  G4double TargetCellRadius      = (0.125 + 0.5*50.)*mm;//(20.179 + 5.)*mm;//20.32*mm;
  G4double TargetWindowThickness = 0.125*mm;//5.*mm;//0.128*mm;

  //---------------------------------------------------------------------------
  // Vacuum inside scattering chamber and exit beamline
  //---------------------------------------------------------------------------

  G4RotationMatrix *xChambRot = new G4RotationMatrix;  
  xChambRot->rotateX(90*degree);                     

  G4Tubs* sChamberInner = new G4Tubs("ChamberInner_sol",
				     0.,
				     ChamberInnerRadius,
				     0.5*ChamberHight,
				     0.,
				     twopi);

//   G4LogicalVolume* LogicInnerChamber = new G4LogicalVolume(sChamberInner, VacuumMaterial, "ChamberInner_log");
  
//   new G4PVPlacement(xChambRot, G4ThreeVector(), LogicInnerChamber, "ChamberInner_pos", fexpHall_log, false, 0 );

  G4Tubs* sBeamlineInner = new G4Tubs("BeamlineInner_sol",
				      0.,
				      BeamlineInnerRadius,
				      BeamlineHight,
				      0.,
				      twopi);

  G4UnionSolid* vacUnion = new G4UnionSolid("vacUnion", sChamberInner, sBeamlineInner, xChambRot, G4ThreeVector(0., -BeamlineHight, 0.));

  G4LogicalVolume* LogicInnerChamber = new G4LogicalVolume(vacUnion, VacuumMaterial, "ChamberInner_log");
  
  new G4PVPlacement(xChambRot, G4ThreeVector(), LogicInnerChamber, "ChamberInner_pos", fexpHall_log, false, 0 );

  //LogicInnerChamber->SetVisAttributes(G4VisAttributes::Invisible);

  //---------------------------------------------------------------------------
  

//   G4LogicalVolume* LogicInnerBeamline = new G4LogicalVolume(sBeamlineInner, VacuumMaterial, "BeamlineInner_log");
  
//   new G4PVPlacement(0, G4ThreeVector(), LogicInnerBeamline, "BeamlineInner_pos", fexpHall_log, false, 0 );
  
  //---------------------------------------------------------------------------
  // Target Cell & Target
  //---------------------------------------------------------------------------

  G4RotationMatrix *xTargetRot = new G4RotationMatrix;  
  xTargetRot->rotateX(-90*degree);                     

  G4Tubs* sTargetCell = new G4Tubs("TargetCell_sol",           
				   0., TargetCellRadius, (0.5*TargetLength)+TargetWindowThickness, 0.,twopi); 

  G4LogicalVolume*   LogicTargetCell = new G4LogicalVolume(sTargetCell, TargetCellMaterial, "TargetCell_log");   

  new G4PVPlacement(xTargetRot, G4ThreeVector(), LogicTargetCell, "TargetCell_pos",  LogicInnerChamber, false, 0 );  

  //---------------------------------------------------------------------------

  G4Tubs* sTarget = new G4Tubs("Target_sol",         
			       0., TargetRadius, 0.5*TargetLength, 0.,twopi); 

  fLogicTarget = new G4LogicalVolume(sTarget, TargetMaterial, "Target_log");  
  
  new G4PVPlacement(0, G4ThreeVector(0., 0., 0.5*TargetWindowThickness), fLogicTarget, "Target_pos", LogicTargetCell, false, 0 ); 

}

//---------------------------------------------------------------------------
