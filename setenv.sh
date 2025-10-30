source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsenv

# source G4Beamline
source G4beamline-3.08/bin/g4bl-setup.sh
source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh

# ---- Geant4 data (CVMFS) ----
export G4ABLADATA=/cvmfs/geant4.cern.ch/share/data/G4ABLA3.3
export G4LEDATA=/cvmfs/geant4.cern.ch/share/data/G4EMLOW8.8
export G4ENSDFSTATEDATA=/cvmfs/geant4.cern.ch/share/data/G4ENSDFSTATE2.3
export G4PARTICLEXSDATA=/cvmfs/geant4.cern.ch/share/data/G4PARTICLEXS4.2
export G4LEVELGAMMADATA=/cvmfs/geant4.cern.ch/share/data/PhotonEvaporation5.7
export G4RADIOACTIVEDATA=/cvmfs/geant4.cern.ch/share/data/RadioactiveDecay5.6
export G4INCLDATA=/cvmfs/geant4.cern.ch/share/data/G4INCL1.3
export G4NEUTRONHPDATA=/cvmfs/geant4.cern.ch/share/data/G4NDL4.7.1
export G4PIIDATA=/cvmfs/geant4.cern.ch/share/data/G4PII1.3
export G4REALSURFACEDATA=/cvmfs/geant4.cern.ch/share/data/RealSurface2.2
export G4SAIDXSDATA=/cvmfs/geant4.cern.ch/share/data/G4SAIDDATA2.0

