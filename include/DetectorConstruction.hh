#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "G4NistManager.hh"
#include "DetectorMessenger.hh"

//---------------------------------------------------------------------------

class G4VPhysicalVolume;
class VirtualDetectorSD;
class RealDetectorSD;

//---------------------------------------------------------------------------

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  DetectorConstruction();
  ~DetectorConstruction();

  G4VPhysicalVolume* Construct();

  void UpdateGeometry();
  void BuildBeamline();
  void BuildTarget();
  
  inline G4VPhysicalVolume* GetExpHall()           { return fExpHall;           };
  inline G4VPhysicalVolume* GetDetVol(G4int i)     { return fDetVol[i];         };
  inline VirtualDetectorSD* GetVirtualDetectorSD() { return fVirtualDetectorSD; };
  inline G4int              GetNoSD()              { return fNSD;               };

  void SetNPSAngle            ( G4double a  ) { fNPSAngle         = a;  }
  void SetNPSDistance         ( G4double d  ) { fNPSDist          = d;  }
  void SetNPSShieldThickness  ( G4double nt ) { fNPSShieldThick   = nt; }
  void SetHCALAngle           ( G4double an ) { fHCALAngle        = an; }
  void SetHCALDistance        ( G4double di ) { fHCALDist         = di; }
  void SetHCALShieldThickness ( G4double t  ) { fHCALShieldThick  = t;  }
  void SetWindowThickness     ( G4double th ) { fSCWinThick       = th; }
  void SetTargetLength        ( G4double l  ) { fTarLength        = l; }
  void SetBeamlineOn          ( G4int b  )    { fBeamline         = b; }

  private:

  G4NistManager*     fNistManager;
  DetectorMessenger* fDetMessenger;

  G4double           fNPSAngle;
  G4double           fNPSDist;
  G4double           fNPSShieldThick;
  G4double           fHCALAngle;
  G4double           fHCALDist;
  G4double           fHCALShieldThick;
  G4double           fSCWinThick;
  G4double           fTarLength;
  G4int              fBeamline;

  static const G4int fNPSNrow  = 5;
  static const G4int fNPSNcol  = 5;

  static const G4int fHCALNrow = 3;
  static const G4int fHCALNcol = 3;

  static const G4int fHodoNrow = 15;
  static const G4int fHodoNcol = 15;
  

  static const G4int fNSD = ( fNPSNrow*fNPSNcol + fHodoNrow*fHodoNcol + fHCALNrow*fHCALNcol + 1 ); 

  G4VPhysicalVolume* fExpHall;
  G4VPhysicalVolume* fDetVol[fNSD];
  G4LogicalVolume*   fexpHall_log;
  G4LogicalVolume*   fLogicTarget;

  VirtualDetectorSD* fVirtualDetectorSD;
  RealDetectorSD*    fRealDetectorSD;

};
#endif

//---------------------------------------------------------------------------

