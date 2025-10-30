# SPSH4BeamlineSimulation

This repository hosts a realistic G4Beamline model of the CERN SPS H4 beamline. Although the experimental area (PPE-164) setup is simulated for ECAL Test Beam studies, it can be modified to suit all detectors conducting test beam studies at the H4 Beamline, in order to accommodate the prevailing conditions at that time.

Load cmssw libraries:
> Use `.csh` if you're a **csh** user. Check with `echo $0` on your terminal (lxplus).

```
source /cvmfs/cms.cern.ch/cmsset_default.sh
```

Setup CMSSW & pull the framework:

```
cmsrel CMSSW_15_0_6
cd CMSSW_15_0_6/src
cmsenv
git cms-init
git clone https://github.com/asimsek/SPSH4BeamlineSimulation
cd SPSH4BeamlineSimulation
``` 


## nTuple Production

[!WARNING]
> Download **G4beamline v3.08** (Linux64 tar file) from [Muons INC website](https://www.muonsinc.com/Website1/tiki-index.php?page=G4beamline#Download) or copy from my EOS area: `/eos/cms/store/user/asimsek/G4beamline-3.08-Linux64.tgz` 




```
tar -xzf G4beamline-3.08-Linux64.tgz
source setenv.sh
```


```
g4bl g4bl_files/H4_Electron_v9_2021.in viewer=none InitEnergy=100 totNumEv=1000 jobID=1
```



## Plotting


#### Main Plotter:

```
python3 ./utils/plotter.py --initialEnergy <initialEnergy> --inputFile <inputRootFileName> --outputFolder <outputFolderName> --crystalball
python3 ./utils/plotter.py --initialEnergy 100 --inputFile /eos/cms/store/user/asimsek/SPSH4Results/v9/H4_Positron_v9_NormalBeam_Air_FullHodoscope_100GeV.root --outputFolder outPDFs/100GeV/ --crystalball
```

[!TIP]
> Use `--goodTree` argument along with the rest of them if you want the triggered events in the nTuple (Good Tree).
> Use `--crystalball` if you want to use Crystall Ball function instead of Gauss.

#### Particle production:

[!TIP]
> Script checks all the root files under the given directory/path!


```
g++ -O2 utils/dump_vd_counts.cc $(root-config --cflags --libs) -o utils/dump_vd_counts
./utils/dump_vd_counts --energy 100 --in ./
```

[!TIP]
> This provides a CSV file containing the number of particles per detector along the beamline. <br>
> The decision of the c++ usage is related with the performance. 
> You can plot this distribution by using the `plot_particle_production.py` script.


```
python3 ./utils/plot_particle_production.py <csv_file> <initial_energy>
```






------------

------------



## New Features!

**July 28, 2021**
###### - [v8 G4BL Script]
 - Geometry Improvements
 - Good Particles

**July 12, 2021**
###### - [v6 G4BL Script]
 - Full Hodoscope design, Preliminary ECAL TB Box design with crystals are added.
 - PPE144 (B21 Area) are re-designed.
 - World Material set as Air.
 - Needs to work on crystals.

**June 18, 2021**
###### - [v4 G4BL Script]
 - Scintillators, Delay Wire Chambers (DWCs) are included. Need to verify with in-situ visits.


**June 06, 2021**
###### - [v3 G4BL Script]
 - Beam Pipes, Al & Mylar Windows, Trims, Air Zones included.
 - First design of the hodoscope included. Need to check with ECAL people for the realistic values (geometry, DUT position etc).


**May 16, 2021**
###### - [v2 G4BL Script]
 - Cleaner & Organized Input File!
 - Synchrotron Radiation (SR) Effect is applied!?
 - Reference Momentum set to 100 GeV.
 - Beam Conditions: 90% e+ and 10% proton.
 - Sigma Calculation after SR Effect.
 - 90° view angle & wireframe for visualization.
 - killSecondaries setting is ON for checking the SR Effect without secondaries.


