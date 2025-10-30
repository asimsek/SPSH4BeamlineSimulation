# -*- coding: utf-8 -*-

from ROOT import *
from array import *
from random import gauss
from tqdm import tqdm
from prettytable import PrettyTable
from dateutil.relativedelta import relativedelta
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
import os, sys
import math
import array
import time

#### Disable All Printouts of Canvases! kPrint, kInfo, kWarning, kError, kBreak, kSysError, kFatal
gROOT.ProcessLine("gErrorIgnoreLevel = kWarning;")

start_time = datetime.now()

gROOT.SetBatch(True)
execfile('tdrStyle.py')


ener = "250"
mult1 = 0.97

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--initialEnergy", type=str, help="Initial Energy of the Beam [300, 250, 200, 150, 100]", default="%s" % (ener))
parser.add_argument("--inputFile", type=str, help="Input Root File created with G4BL [*.root]", default="/eos/cms/store/user/asimsek/SPSH4Results/v9/H4_Positron_v9_NormalBeam_Air_FullHodoscope_%sGeV.root" % (ener))
parser.add_argument("--outputFolder", type=str, help="Output Folder for PDFs and Output Roots [Results]", default="Results_%sGeV_H4_v9_gauss" % (ener))


args = parser.parse_args()
initialEnergy = args.initialEnergy
inputFile = args.inputFile
outputFolder = args.outputFolder


################# Variables ###############
#### Don't use this (includes all virtual detectors - not necessary!! Use full beamline one)
#virtualDetectors = ["BeamStartPoint","BeforeQuad1","AfterQuad3","BeforeColi1","AfterColi1","BeforeQuad4","BeforeB5a","BeforeColi2","AfterColi3","BeforeTrim01","AfterTrim01","AfterB5h","AfterQuad8a","BeforeTrim02","AfterTrim02","BeforeQuad8b","AfterQuad7b","BeforeAir01","BeforeQuad6b","BeforeB6a","AfterB7b","AfterB7c","BeforeTrim03","AfterTrim03","BeforeB8a","AafterB8b","BeforeB9a","AfterB9b","AfterCol10","AfterScint06","AfterTrim06","BeforeQ23","AfterQ24","NA64","BeforeXTDV22601","AfterXTDV22601","BeforeEcalAirZone","ECALFront","ECALBack"]


### Full Beamline
#virtualDetectors = ["BeamStartPoint","AfterQuad3","BeforeQuad4","BeforeB5a","BeforeColi2","AfterB5h","AfterQuad8a","BeforeTrim02","AfterTrim02","BeforeQuad8b","BeforeAir01","BeforeQuad6b","BeforeB6a","AfterB7b","AfterB7c","BeforeTrim03","AfterTrim03","BeforeB8a","AafterB8b","AfterB9b","AfterCol10","AfterScint06","AfterTrim06","BeforeQ23","AfterQ24","NA64","BeforeXTDV22601","AfterXTDV22601","BeforeEcalAirZone","ECALFront"]

### After Q24
#virtualDetectors = ["BeamStartPoint", "AfterQ24","NA64","BeforeXTDV22601","Det022640","ECALFront"]
virtualDetectors = ["ECALFront"]

#### Only Start & End
#virtualDetectors = ["BeamStartPoint","ECALFront"]


detectorListForPlotting = ["ECALFront"]


colors = [42,46,32,36,16,2]
sigma = 3.0 #### 
mean = 0.0
sigmaGauss = 0.0
tailEnergyCut = 0.0
tailInvestigation = False ## Make it true for tail investigation
###########################################


def percentage(part, whole):
  Percentage = 100 * float(part)/float(whole)
  #Percentage = "%0.3f" % Percentage
  return Percentage

#####################################################
######## Seconds to more understandable time format

def diffTime(t_a, t_b):
    t_diff = relativedelta(t_b, t_a)  # later/end time comes first!
    return '{h}h {m}m {s}s'.format(h=t_diff.hours, m=t_diff.minutes, s=t_diff.seconds)

###############################################
####### Read Mean & sigma values from txt file
with open('GaussTable.txt') as f:
	datafile = f.readlines()
	for line in datafile:
		if str(initialEnergy) in line:
			info = line.split("|")

			meanN = float(info[1])
			sigmaGauss = float(info[2])
			binning = float(info[3])
			print(" ->> Energy: \033[0;32m%s\033[0;0m - Mean: \033[0;32m%f\033[0;0m - Sigma: \033[0;32m%f\033[0;0m" %(str(info[0]), meanN, sigmaGauss) )
###############################################


###### Defining the x-axis momentum range
xmin_ = float(meanN)-((float(meanN)*sigma)/100)
xmax_ = float(meanN)+((float(meanN)*sigma)/100)
#########################################


#xmin_ = 0.0    #Uncomment this line to show tail starting from 0.0
ctr = -1
ctr3 = 0

if initialEnergy == "120" and xmin_ > 110:
	xmin_ = 116

#########################################
#### Definition of the Dictionaries
EventSize = {}
PositronSize = {}
ElectronSize = {}
ProtonSize = {}
PhotonSize = {}
MuonSize = {}
PionMinusSize = {}
PionPlusSize = {}

PositronRate = {}
ElectronRate = {}
ProtonRate = {}
#PhotonRate = {}
PhtnRt = {}
MuonRate = {}
PionMinusRate = {}
PionPlusRate = {}

#### Definition of the Arrays
TotalEventArray = array.array('i', [])

TotalEventArray = array.array('i', [])
PositronArray = array.array('i', [])
ElectronArray = array.array('i', [])
ProtonArray = array.array('i', [])
PhotonArray = array.array('i', [])
MuonArray = array.array('i', [])
PionMinusArray = array.array('i', [])
PionPlusArray = array.array('i', [])


PhotonEventSizeInitZ = {}
table1 = PrettyTable()

#########################################
#########################################
########### Root File Check ############
inputRoot = TFile.Open(inputFile, 'READ')
if inputRoot:
	print ("\033[1;31m ->> Root file has been opened! \033[0;0m")

	table1.field_names = ["Detector", "ALL Particles", "Positron", "Electron", "Photon", "Muon", "Proton"]

	########### Virtual Detector Loop ############
	#for i, Det in enumerate(virtualDetectors):
	for i in tqdm(range(len(virtualDetectors))):
		Det = virtualDetectors[i]
		######### Create an Output Folder ############
		#if Det=="ECALFront" or Det=="BeamStartPoint":
		if Det in detectorListForPlotting:
			import os ## I'm importing os here to pass some stupid python bugs!! (don't worry, it's on purpose!)
			os.system("mkdir -p %s/%s" % (outputFolder, Det))
			
		#print (" ->>\033[1;95m Output Folder:\033[0m \033[1;96m%s\033[0m" % (outputFolder))

		############################################################
		#################### Tail Investigation ####################
		############################################################

		if tailInvestigation == True:
			tailEnergyCut = float(meanN)-(float(meanN)*sigma/100)
			#print ("Tail Energy Cut: %f" % (tailEnergyCut))
		else:
			tailEnergyCut = float(initialEnergy)*2
			#print ("Tail Energy Cut: %f" % (tailEnergyCut))
		Tree1 = TTree()
		inputRoot.GetObject("%s" % (str(Det)), Tree1)
		tqdm.write("\033[;1m --> Collecting data from\033[0;0m \033[1;36m%s\033[0;0m" % (str(Det)))
		#print ("\033[;1m --> Collecting data from\033[0;0m \033[1;36m%s\033[0;0m" % (str(Det)) )
		EventSize[str(Det)] = Tree1.GetEntries()
		PositronSize[str(Det)] = Tree1.GetEntries("PDGid==-11 && (sqrt((Px*Px)+(Py*Py)+(Pz*Pz))/1000) <= %0.2f" % (tailEnergyCut) )
		ElectronSize[str(Det)] = Tree1.GetEntries("PDGid==11 && (sqrt((Px*Px)+(Py*Py)+(Pz*Pz))/1000) <= %0.2f" % (tailEnergyCut) )
		ProtonSize[str(Det)] = Tree1.GetEntries("PDGid==2212 && (sqrt((Px*Px)+(Py*Py)+(Pz*Pz))/1000) <= %0.2f" % (tailEnergyCut) )
		PhotonSize[str(Det)] = Tree1.GetEntries("PDGid==22 && (sqrt((Px*Px)+(Py*Py)+(Pz*Pz))/1000) <= %0.2f" % (tailEnergyCut) )
		MuonSize[str(Det)] = Tree1.GetEntries("PDGid==13 && (sqrt((Px*Px)+(Py*Py)+(Pz*Pz))/1000) <= %0.2f" % (tailEnergyCut) )
		PionMinusSize[str(Det)] = Tree1.GetEntries("PDGid==-211 && (sqrt((Px*Px)+(Py*Py)+(Pz*Pz))/1000) <= %0.2f" % (tailEnergyCut) )
		PionPlusSize[str(Det)] = Tree1.GetEntries("PDGid==211 && (sqrt((Px*Px)+(Py*Py)+(Pz*Pz))/1000) <= %0.2f" % (tailEnergyCut) )
		
		TotalEventArray.append(EventSize[str(Det)])
		PositronArray.append(PositronSize[str(Det)])
		ElectronArray.append(ElectronSize[str(Det)])
		ProtonArray.append(ProtonSize[str(Det)])
		PhotonArray.append(PhotonSize[str(Det)])
		MuonArray.append(MuonSize[str(Det)])
		PionMinusArray.append(PionMinusSize[str(Det)])
		PionPlusArray.append(PionPlusSize[str(Det)])

		PositronRate[str(Det)] = percentage(PositronSize[str(Det)], EventSize[str(Det)])
		ElectronRate[str(Det)] = percentage(ElectronSize[str(Det)], EventSize[str(Det)])
		ProtonRate[str(Det)] = percentage(ProtonSize[str(Det)], EventSize[str(Det)])
		PhtnRt[str(Det)] = percentage(PhotonSize[str(Det)], EventSize[str(Det)])
		MuonRate[str(Det)] = percentage(MuonSize[str(Det)], EventSize[str(Det)])
		PionMinusRate[str(Det)] = percentage(PionMinusSize[str(Det)], EventSize[str(Det)])
		PionPlusRate[str(Det)] = percentage(PionPlusSize[str(Det)], EventSize[str(Det)])

		ParticleRates = [PositronRate[str(Det)], ElectronRate[str(Det)], PhtnRt[str(Det)], ProtonRate[str(Det)]]

		if PhtnRt[str(Det)] <= 0.0001: PhtnRt[str(Det)] = 0.0001 

		#Total = float(PositronRate[str(Det)]) + float(ElectronRate[str(Det)]) + float(ProtonRate[str(Det)]) + float(PhotonRate[str(Det)]) + float(MuonRate[str(Det)]) + float(PionMinusRate[str(Det)]) + float(PionPlusRate[str(Det)])

		table1.add_row((str(Det), EventSize[str(Det)], PositronSize[str(Det)], ElectronSize[str(Det)], PhotonSize[str(Det)], MuonSize[str(Det)], ProtonSize[str(Det)]))
		
		del Tree1

		if tailInvestigation == True:
			continue

		
		###################################################################
		###################################################################
		###################################################################

		#if Det=="ECALFront":
		if Det in detectorListForPlotting:

			TreePt = inputRoot.Get("GoodParticle")
			#TreePt = inputRoot.Get("%s" % (str(Det)))
			ctr += 1
                        TreePt.Draw("sqrt((%s_Px*%s_Px)+(%s_Py*%s_Py)+(%s_Pz*%s_Pz))/1000>>histoPTotal(%f ,%f, %f)" % (Det, Det, Det, Det, Det, Det, binning, xmin_, xmax_) )
			#TreePt.Draw("sqrt((Px*Px)+(Py*Py)+(Pz*Pz))/1000>>histoPTotal(%f ,%f, %f)" % (binning, xmin_, xmax_) )

			c1 = TCanvas('c1_%s' % (Det), 'c1', 900, 1200)
			c1.cd()

			if Det == "ECALFront" or Det == "ECALBack":
				Detect = "ECAL"
			else:
				Detect = str(Det)
			histoPTotal.SetTitle("Total Momentum Spread at the %s - All Particles" % (str(Detect)))
			histoPTotal.GetYaxis().SetTitle("Number Of Events")
			histoPTotal.GetXaxis().SetTitle("P_{Total} [GeV]")
			histoPTotal.SetLineColor(colors[ctr])
			histoPTotal.SetFillColor(colors[ctr])
			histoPTotal.SetMarkerColor(kBlue-2)
			histoPTotal.SetFillStyle(3001)
			histoPTotal.GetXaxis().SetRangeUser(xmin_, xmax_)

			meanFit = meanN
			sigmaFit = sigmaGauss

			gaussFit = TF1("gaussfit", "[0]*exp(-0.5*((x-[1])/[2])^2)", float(meanFit)-((float(meanFit)*sigmaFit)/100), float(meanFit)+((float(meanFit)*sigmaFit)/100))
			#gaussFit.SetParameter(0, 1000)
			gaussFit.SetParameter(1, meanFit)
			gaussFit.SetParameter(2, sigmaFit)
			gaussFit.SetParLimits(1,float(meanFit)-((float(meanFit)*sigmaFit)/100),float(meanFit)+((float(meanFit)*sigmaFit)/100));
			gaussFit.SetParLimits(2, -1*sigmaFit, sigmaFit)
			histoPTotal.Fit(gaussFit, "") # "q"
			histoPTotal.Draw("")
			gaussFit.Draw("SAME")

			chi2 = gaussFit.GetChisquare()
			NDF = gaussFit.GetNDF()
			#meanFit = gaussFit.GetParameter(1)
			#sigmaFit = gaussFit.GetParameter(2)
			meanMC = histoPTotal.GetMean()
			sigmaMC = histoPTotal.GetStdDev()

			

			#print (" -> Mean=%0.3f  |  Sigma=%0.2f  |  Chi2=%0.2f  |  NDF=%0.2f" % (meanFit, sigmaFit, chi2, NDF) )
			leg = TLegend(0.67, 0.80, 0.99, 0.97)
			leg.SetBorderSize(5)
			leg.SetTextSize(0.022)
			leg.AddEntry(histoPTotal, "Monte Carlo [MC]", "f")
			leg.AddEntry(histoPTotal, "Mean: %0.2f" % (meanMC) , "")
			leg.AddEntry(histoPTotal, "Std. Dev.: %0.2f" % (sigmaMC) , "")
			leg.AddEntry(histoPTotal, "-----------------------", "")
			leg.AddEntry(gaussFit, "Th. #DeltaP/P due to SR", "l")
			leg.AddEntry(gaussFit, "Mean: %0.2f" % (meanFit) , "")
			leg.AddEntry(gaussFit, "Std. Dev.: %0.6f" % (sigmaFit) , "")
			leg.Draw()
			c1.SaveAs("%s/%s/PTotal_%s.pdf" % (outputFolder, Det, Det) )
			c1.Close()

			
			del c1
			del histoPTotal
			del TreePt

			###################################################################
			################# All vs Good Particles Only e^{+} ################
			###################################################################

		if Det=="ECALFront2":
			
			TreeGoodParticles = inputRoot.Get("NTuples/GoodParticle")
			TreeAllParticles = inputRoot.Get("%s" % (str(Det)))
			TreeGoodParticles.Draw(("sqrt((%s_Px*%s_Px)+(%s_Py*%s_Py)+(%s_Pz*%s_Pz))/1000>>histoPTotalGood(%f ,%f, %f)" % (str(Det),str(Det),str(Det),str(Det),str(Det),str(Det), xmax_, 0, xmax_)), ("%s_PDGid==-11" % (str(Det))) )
			TreeAllParticles.Draw( ("sqrt((Px*Px)+(Py*Py)+(Pz*Pz))/1000>>histoPTotal2(%f ,%f, %f)" % (xmax_, 0, xmax_)), "PDGid==-11" )

			c1Good = TCanvas('c1Good', 'c1Good', 900, 1200)
			c1Good.cd()
			c1Good.SetLogy(True)


			histoPTotal2.SetTitle("Triggered Event Effect (Good Particles) - Only e^{+}")
			histoPTotal2.GetYaxis().SetTitle("Number Of Events")
			histoPTotal2.GetXaxis().SetTitle("P_{Total} [GeV]")
			histoPTotal2.SetLineColor(colors[ctr])
			histoPTotal2.SetFillColor(colors[ctr])
			histoPTotal2.SetFillStyle(3001)
			histoPTotal2.SetNdivisions(-502)
			histoPTotal2.GetXaxis().SetRangeUser(0, xmax_)

			histoPTotalGood.SetLineColor(colors[ctr+1])
			histoPTotalGood.SetFillColor(colors[ctr+1])
			histoPTotalGood.SetFillStyle(3001)
			histoPTotalGood.SetNdivisions(-502)
			histoPTotalGood.GetXaxis().SetRangeUser(0, xmax_)

			meanMC2 = histoPTotal2.GetMean()
			sigmaMC2 = histoPTotal2.GetStdDev()

			meanMCGood = histoPTotalGood.GetMean()
			sigmaMCGood = histoPTotalGood.GetStdDev()


			histoPTotal2.Draw("")
			histoPTotalGood.Draw("SAME")


			leg = TLegend(0.67, 0.80, 0.99, 0.97)
			leg.SetBorderSize(5)
			leg.SetTextSize(0.022)
			leg.AddEntry(histoPTotal2, "All e^{+} [MC]", "f")
			leg.AddEntry(histoPTotal2, "Event Size: %d" % (histoPTotal2.GetEntries()) , "")
			leg.AddEntry(histoPTotal2, "-----------------------", "")
			leg.AddEntry(histoPTotalGood, "Good e^{+} [MC]", "f")
			leg.AddEntry(histoPTotalGood, "Event Size: %d" % (histoPTotalGood.GetEntries()) , "")
			leg.Draw()


			c1Good.SaveAs("%s/%s/PTotalGoodParticles_%s.pdf" % (outputFolder, Det, Det) )
			c1Good.Close()

			del histoPTotal2
			del histoPTotalGood
			del TreeGoodParticles
			del TreeAllParticles
			del c1Good

			

			###################################################################
			#################### PTotal, X & Y Stack Plot #####################
			###################################################################

		if Det in detectorListForPlotting:
			
			Tree2 = inputRoot.Get("%s" % (str(Det)))
			if Det == "ECALFront" or Det == "ECALBack":
				Detect = "ECAL"
			else:
				Detect = str(Det)
			ths1 = THStack("PtotalStack","Total Momentum Spectrum at %s;P_{Total} [GeV];Number Of Events" % (str(Detect)))
			thsX1 = THStack("XStack","X - Beam Position at %s;X-Position [mm];Number Of Events" % (str(Detect)))
			thsY1 = THStack("YStack","Y - Beam Position at %s;Y-Position [mm];Number Of Events" % (str(Detect)))
			thsInitZ = THStack("InitZStack","Beam Particle Production Areas - %s;Z-Position [m];Number Of Events" % (str(Detect)))

			leg6 = TLegend(0.42, 0.76, 0.74, 0.90)
			leg6.SetBorderSize(5)
			leg6.SetTextSize(0.022)

			leg7 = TLegend(0.73, 0.80, 0.99, 0.97)
			leg7.SetBorderSize(5)
			leg7.SetTextSize(0.019)

			leg8 = TLegend(0.73, 0.80, 0.99, 0.97)
			leg8.SetBorderSize(5)
			leg8.SetTextSize(0.019)


			leg9 = TLegend(0.73, 0.80, 0.99, 0.97)
			leg9.SetBorderSize(5)
			leg9.SetTextSize(0.019)

			histoPTot = []
			hx = []
			hy = []
			hInitz = []
			PDGids = [-11, 11, 22, 2212]
			ParticlesPDGid = ["Positron","Electron","Photon","Proton"]
			xminParticles = 0
			#xminParticles = xmin_
			ctr2 = -1
			#ix = detectorListForPlotting.index(Det)
			#print ("#%d - %s - ctr2:%d" % (ix, Det, ctr2) )
			for x in PDGids:
				ctr2 += 1
				ctr3 += 1
				histoPTot.append(TH1F("histoPTot2[%d]" % (ctr3), "Total Momentum Spectrum of %s" % (ParticlesPDGid[ctr2]), int(binning), int(xminParticles), int(xmax_+3)))
				hx.append(TH1F("hx[%d]" % (ctr3), "X - Beam Profile of %s" % (ParticlesPDGid[ctr2]), 200, -100, 100))
				hy.append(TH1F("hy[%d]" % (ctr3), "Y - Beam Profile of %s" % (ParticlesPDGid[ctr2]), 200, -100, 100))
				hInitz.append(TH1F("hInitz[%d]" % (ctr3), "Beam Particle Production Area of %s" % (ParticlesPDGid[ctr2]), 1010, -10, 1000))

				Tree2.Draw("sqrt((Px*Px)+(Py*Py)+(Pz*Pz))/1000>>histoPTot2[%d]" % (ctr3), "PDGid==%d" % (x))
				Tree2.Draw("x>>hx[%d]" % (ctr3), "PDGid==%d" % (x))
				Tree2.Draw("y>>hy[%d]" % (ctr3), "PDGid==%d" % (x))
				Tree2.Draw("InitZ/1000>>hInitz[%d]" % (ctr3), "PDGid==%d" % (x))


				histoPTot[ctr2].SetLineColor(colors[ctr2])
				histoPTot[ctr2].SetFillColor(colors[ctr2])
				histoPTot[ctr2].SetFillStyle(3001)
				histoPTot[ctr2].SetNdivisions(-502)

				hx[ctr2].SetLineColor(colors[ctr2])
				hx[ctr2].SetFillColor(colors[ctr2])
				hx[ctr2].SetFillStyle(3001)
				#hx[ctr2].SetNdivisions(-502)

				hy[ctr2].SetLineColor(colors[ctr2])
				hy[ctr2].SetFillColor(colors[ctr2])
				hy[ctr2].SetFillStyle(3001)

				hInitz[ctr2].SetLineColor(colors[ctr2])
				hInitz[ctr2].SetFillColor(colors[ctr2])
				hInitz[ctr2].SetFillStyle(3001)

				gROOT.cd()
				ths1.Add(histoPTot[ctr2])
				gROOT.cd()
				thsX1.Add(hx[ctr2])
				gROOT.cd()
				thsY1.Add(hy[ctr2])
				gROOT.cd()
				thsInitZ.Add(hInitz[ctr2])

				leg6.AddEntry(histoPTot[ctr2], "%s | %f%%" % ( ParticlesPDGid[ctr2], ParticleRates[ctr2] ) , "f")
				leg7.AddEntry(hx[ctr2], "%s | %f%%" % ( ParticlesPDGid[ctr2], ParticleRates[ctr2] ) , "f")
				leg8.AddEntry(hy[ctr2], "%s | %f%%" % ( ParticlesPDGid[ctr2], ParticleRates[ctr2] ) , "f")
				leg9.AddEntry(hInitz[ctr2], "%s | %f%%" % ( ParticlesPDGid[ctr2], ParticleRates[ctr2] ) , "f")

		
			#ths1.ls()

			cTstack1 = TCanvas('cTstack1', 'cTstack1', 900, 1200)
			cTstack1.cd()
			cTstack1.SetLogy(True)

			ths1.Draw("nostack")
			leg6.Draw()

			cTstack1.SaveAs("%s/%s/PTotalStack_%s.pdf" % (outputFolder, Det, Det) )
			cTstack1.Close()

			cTstackXPos = TCanvas('cTstackXPos', 'cTstackXPos', 900, 1200)
			cTstackXPos.cd()
			#cTstackXPos.SetLogy(True)

			thsX1.Draw("nostack")
			leg7.Draw()

			cTstackXPos.SaveAs("%s/%s/XPositionStack_%s.pdf" % (outputFolder, Det, Det) )
			cTstackXPos.Close()


			cTstackYPos = TCanvas('cTstackYPos', 'cTstackYPos', 900, 1200)
			cTstackYPos.cd()
			#cTstackYPos.SetLogy(True)

			thsY1.Draw("nostack")
			leg8.Draw()

			cTstackYPos.SaveAs("%s/%s/YPositionStack_%s.pdf" % (outputFolder, Det, Det) )
			cTstackYPos.Close()

			cTstackInitZ = TCanvas('cTstackInitZ', 'cTstackInitZ', 1200, 900)
			cTstackInitZ.cd()
			#cTstackInitZ.SetLogy(True)

			thsInitZ.Draw("nostack")
			leg9.Draw()

			cTstackInitZ.SaveAs("%s/%s/InitZStack_%s.pdf" % (outputFolder, Det, Det) )
			cTstackInitZ.Close()

			del histoPTot
			del hy
			del hx
			del hInitz
			del ths1
			del thsX1
			del thsY1
			del cTstack1
			del cTstackXPos
			del cTstackYPos
			del cTstackInitZ
			del Tree2


			###### Zoom in InitZ (After 600m)
			thsInitZA600 = THStack("InitZA600Stack","Beam Particle Production Areas - %s;Z-Position [m];Number Of Events" % (str(Detect)))
			Tree3 = inputRoot.Get("%s" % (str(Det)))
			hInitzA600 = []
			ctr4 = -1

			leg10 = TLegend(0.73, 0.80, 0.99, 0.97)
			leg10.SetBorderSize(5)
			leg10.SetTextSize(0.019)


			for x2 in PDGids:
				ctr4 += 1
				hInitzA600.append(TH1F("hInitz2A600[%d]" % (ctr3 + ctr4), "Beam Particle Production Area of %s - After 600m" % (ParticlesPDGid[ctr4]), 100, 600, 700))
				Tree3.Draw("InitZ/1000>>hInitz2A600[%d]" % (ctr3 + ctr4), "PDGid==%d" % (x2))

				gROOT.cd()
				thsInitZA600.Add(hInitzA600[ctr4])
				leg10.AddEntry(hInitzA600[ctr4], "All %s | %f%%" % ( ParticlesPDGid[ctr4], ParticleRates[ctr4] ) , "f")

				hInitzA600[ctr4].SetLineColor(colors[ctr4])
				hInitzA600[ctr4].SetFillColor(colors[ctr4])
				hInitzA600[ctr4].SetFillStyle(3001)

			cTstackInitZA600 = TCanvas('cTstackInitZA600', 'cTstackInitZA600', 1200, 900)
			cTstackInitZA600.cd()
			#cTstackInitZA600.SetLogy(True)

			thsInitZA600.Draw("nostack")
			leg10.Draw()

			cTstackInitZA600.SaveAs("%s/%s/InitZStackAfter600m_%s.pdf" % (outputFolder, Det, Det) )
			cTstackInitZA600.Close()

			
			del hInitzA600
			del thsInitZA600
			del cTstackInitZA600
			

			############################################################
			############################################################

		#if Det=="ECALFront" or Det=="BeamStartPoint":
		if Det in detectorListForPlotting: 

			
			Tree3 = inputRoot.Get("%s" % (str(Det)))
			Tree3.Draw("x>>histoX(%f ,%f, %f)" % (200, -100, 100) )

			histoX.SetTitle("Beam position %s - All Particles" % (str(Det)))
			histoX.GetYaxis().SetTitle("Number Of Events")
			histoX.GetXaxis().SetTitle("X Position (mm)")
			histoX.GetYaxis().SetTitleOffset(1.3)
			histoX.GetXaxis().SetTitleOffset(0.8)
			histoX.GetXaxis().SetNdivisions(10)
			histoX.GetYaxis().SetNdivisions(16)
			histoX.SetLineColor(colors[ctr+1])
			histoX.SetFillColor(colors[ctr+1])
			histoX.SetMarkerColor(kBlue-2)
			histoX.SetFillStyle(3001)

			c2 = TCanvas('c2', 'c2', 900, 1200)
			c2.cd()

			histoX.Draw("HIST")
			
			leg2 = TLegend(0.68, 0.85, 0.98, 0.93)
			leg2.SetBorderSize(4)
			leg2.SetTextSize(0.022)
			leg2.AddEntry(histoX, "Entries: %0.2f" % (histoX.GetEntries()), "f")
			leg2.AddEntry(histoX, "Mean: %0.3e" % (histoX.GetMean()), "")
			leg2.AddEntry(histoX, "Std. Dev.: %.3e" % (histoX.GetStdDev()), "")
			leg2.Draw()
			
			pave = TPaveText(0.18, 0.78, 0.48, 0.92, "blNDC")
			pave.SetFillColor(0)
			pave.SetBorderSize(0)
			pave.SetFillStyle(0)
			pave.SetTextAlign(12)
			pave.SetTextSize(0.025)
			t1 = pave.AddText("Particle Rates")
			t1.SetTextColor(kRed+3)
			t1.SetTextSize(0.03)
			t2 = pave.AddText(" ")
			t2.SetTextColor(kRed+3)
			t2.SetTextSize(0.03)
			pave.AddText("#gamma: %0.3f%%" % (ParticleRates[2]) )
			pave.AddText("e^{+}: %0.3f%%" % (PositronRate[str(Det)]) )
			pave.AddText("e^{-}: %0.3f%%" % (ElectronRate[str(Det)]) )
 			pave.AddText("#mu^{-}: %0.3f%%" % (MuonRate[str(Det)]) )
 			pave.AddText("Proton: %0.3f%%" % (ProtonRate[str(Det)]) )
 			pave.Draw()

			c2.SaveAs("%s/%s/XPos_%s.pdf" % (outputFolder, Det, Det))
			c2.Close()
			

			del histoX
			del c2
			del Tree3

			############################################################
			############################################################


			
			Tree33 = inputRoot.Get("%s" % (str(Det)))
			Tree33.Draw("y>>histoY(%f ,%f, %f)" % (200, -100, 100) )

			histoY.SetTitle("Beam position %s - All Particles" % (str(Det)))
			histoY.GetYaxis().SetTitle("Number Of Events")
			histoY.GetXaxis().SetTitle("Y Position (mm)")
			histoY.GetYaxis().SetTitleOffset(1.3)
			histoY.GetXaxis().SetTitleOffset(0.8)
			histoY.GetXaxis().SetNdivisions(10)
			histoY.GetYaxis().SetNdivisions(32)
			histoY.SetLineColor(colors[ctr+1])
			histoY.SetFillColor(colors[ctr+1])
			histoY.SetMarkerColor(kBlue-2)
			histoY.SetFillStyle(3001)


			c3 = TCanvas('c3', 'c3', 900, 1200)
			c3.cd()

			histoY.Draw("HIST")
			

			leg3 = TLegend(0.68, 0.85, 0.98, 0.93)
			leg3.SetBorderSize(4)
			leg3.SetTextSize(0.025);
			leg3.AddEntry(histoY, "Entries: %0.2f" % (histoY.GetEntries()), "f")
			leg3.AddEntry(histoY, "Mean: %0.2f" % (histoY.GetMean()), "")
			leg3.AddEntry(histoY, "Std. Dev.: %.3e" % (histoY.GetStdDev()), "")
			leg3.Draw()

			pave = TPaveText(0.18, 0.78, 0.48, 0.92, "blNDC")
			pave.SetFillColor(0)
			pave.SetBorderSize(0)
			pave.SetFillStyle(0)
			pave.SetTextAlign(12)
			pave.SetTextSize(0.025)
			t1 = pave.AddText("Particle Rates")
			t1.SetTextColor(kRed+3)
			t1.SetTextSize(0.03)
			t2 = pave.AddText(" ")
			t2.SetTextColor(kRed+3)
			t2.SetTextSize(0.03)
			pave.AddText("#gamma: %0.3f%%" % (ParticleRates[2]) )
			pave.AddText("e^{+}: %0.3f%%" % (PositronRate[str(Det)]) )
			pave.AddText("e^{-}: %0.3f%%" % (ElectronRate[str(Det)]) )
 			pave.AddText("#mu^{-}: %0.3f%%" % (MuonRate[str(Det)]) )
 			pave.AddText("Proton: %0.3f%%" % (ProtonRate[str(Det)]) )
 			pave.Draw()


			c3.SaveAs("%s/%s/YPos_%s.pdf" % (outputFolder, Det, Det))
			c3.Close()

			

			del histoY
			del c3
			del Tree33

		############################################################
		############################################################

		if Det=="BeamStartPoint":

			
			Tree4 = inputRoot.Get("%s" % (str(Det)))
			Tree4.Draw("x:atan(Px/Pz)>>histodPhiX(%f ,%f, %f, %f, %f, %f)" % (100, -0.002, 0.002, 100, -3, 3),"", "colz")

			histodPhiX.SetTitle("Injection Point Beam Profile - #phi vs X %s" % (str(Det)))
			histodPhiX.GetYaxis().SetTitle("X_{[mm]}")
			histodPhiX.GetXaxis().SetTitle("#phi_{[rad]}")
			histodPhiX.GetYaxis().SetTitleOffset(1.3)
			histodPhiX.GetXaxis().SetTitleOffset(0.8)
			histodPhiX.GetYaxis().SetLabelSize(0.04)
			histodPhiX.GetXaxis().SetLabelSize(0.04)
			histodPhiX.GetXaxis().SetNdivisions(4)
			histodPhiX.GetYaxis().SetNdivisions(16)

			c4 = TCanvas('c4', 'c4', 900, 1200)
			c4.cd()			
			gPad.SetRightMargin(0.18)

			histodPhiX.Draw("colz")

			leg4 = TLegend(0.65, 0.9, 0.91, 0.95)
			leg4.SetBorderSize(5)
			leg4.SetTextSize(0.022)
			leg4.AddEntry(histodPhiX, "Entries: %0.0f" % (histodPhiX.GetEntries()), "p")
			leg4.Draw()

			pave = TPaveText(0.18, 0.78, 0.48, 0.92, "blNDC")
			pave.SetFillColor(0)
			pave.SetBorderSize(0)
			pave.SetFillStyle(0)
			pave.SetTextAlign(12)
			pave.SetTextSize(0.025)
			t1 = pave.AddText("Particle Rates")
			t1.SetTextColor(kRed+3)
			t1.SetTextSize(0.03)
			t2 = pave.AddText(" ")
			t2.SetTextColor(kRed+3)
			t2.SetTextSize(0.03)
			pave.AddText("#gamma: %0.3f%%" % (ParticleRates[2]) )
			pave.AddText("e^{+}: %0.3f%%" % (PositronRate[str(Det)]) )
			pave.AddText("e^{-}: %0.3f%%" % (ElectronRate[str(Det)]) )
 			pave.AddText("#mu^{-}: %0.3f%%" % (MuonRate[str(Det)]) )
 			pave.AddText("Proton: %0.3f%%" % (ProtonRate[str(Det)]) )
 			pave.Draw()


			c4.SaveAs("%s/%s/dPhiXPos_%s.pdf" % (outputFolder, Det, Det))
			c4.Close()

			del Tree4
			del c4
			del histodPhiX

			############################################################
			############################################################


			
			Tree5 = inputRoot.Get("%s" % (str(Det)))
			Tree5.Draw("y:atan(Py/Pz)>>histodPhiY(%f ,%f, %f, %f, %f, %f)" % (100, -0.002, 0.002, 100, -3, 3),"", "colz")

			histodPhiY.SetTitle("Injection Point Beam Profile - #phi vs X %s" % (str(Det)))
			histodPhiY.GetYaxis().SetTitle("Y_{[mm]}")
			histodPhiY.GetXaxis().SetTitle("#phi_{[rad]}")
			histodPhiY.GetYaxis().SetTitleOffset(1.3)
			histodPhiY.GetXaxis().SetTitleOffset(0.8)
			histodPhiY.GetYaxis().SetLabelSize(0.04)
			histodPhiY.GetXaxis().SetLabelSize(0.04)
			histodPhiY.GetXaxis().SetNdivisions(4)
			histodPhiY.GetYaxis().SetNdivisions(16)



			c5 = TCanvas('c5', 'c5', 900, 1200)
			c5.cd()
			gPad.SetRightMargin(0.18)


			histodPhiY.Draw("colz")

			leg5 = TLegend(0.65, 0.9, 0.91, 0.95)
			leg5.SetBorderSize(5)
			leg5.SetTextSize(0.022)
			leg5.AddEntry(histodPhiY, "Entries: %0.0f" % (histodPhiY.GetEntries()), "p")
			leg5.Draw()

			pave = TPaveText(0.18, 0.78, 0.48, 0.92, "blNDC")
			pave.SetFillColor(0)
			pave.SetBorderSize(0)
			pave.SetFillStyle(0)
			pave.SetTextAlign(12)
			pave.SetTextSize(0.025)
			t1 = pave.AddText("Particle Rates")
			t1.SetTextColor(kRed+3)
			t1.SetTextSize(0.03)
			t2 = pave.AddText(" ")
			t2.SetTextColor(kRed+3)
			t2.SetTextSize(0.03)
			pave.AddText("#gamma: %0.3f%%" % (ParticleRates[2]) )
			pave.AddText("e^{+}: %0.3f%%" % (PositronRate[str(Det)]) )
			pave.AddText("e^{-}: %0.3f%%" % (ElectronRate[str(Det)]) )
 			pave.AddText("#mu^{-}: %0.3f%%" % (MuonRate[str(Det)]) )
 			pave.AddText("Proton: %0.3f%%" % (ProtonRate[str(Det)]) )
 			pave.Draw()


			c5.SaveAs("%s/%s/dPhiYPos_%s.pdf" % (outputFolder, Det, Det))
			c5.Close()

			del Tree5
			del c5
			del histodPhiY


		############################################################
		############################################################


		#if Det=="ECALFront" or Det=="BeamStartPoint":
		if Det in detectorListForPlotting:

			if Det=="ECALFront":
				binYAxis = 300
				BeginYAxis = -150
				EndYAxis = 150
				BeginPtot = -10
			if Det=="BeamStartPoint":
				binYAxis = 6
				BeginYAxis = -3
				EndYAxis = 3
				BeginPtot = float(initialEnergy) - 3

			
			Tree6 = TTree()
			Tree6 = inputRoot.Get("%s" % (str(Det)))
			Tree6.Draw("y:sqrt((Px*Px)+(Py*Py)+(Pz*Pz))/1000>>histoYvsPTotal(%f ,%f, %f, %f, %f, %f)" % (binning+10, -10, float(initialEnergy), binYAxis, BeginYAxis, EndYAxis), "", "colz")
			
			histoYvsPTotal.SetTitle("Beam Total Momentum vs Beam Y Position - %s" % (str(Det)))
			histoYvsPTotal.GetYaxis().SetTitle("Beam Y-Position [mm]")
			histoYvsPTotal.GetXaxis().SetTitle("P_{Total} [GeV]")
			histoYvsPTotal.GetYaxis().SetTitleOffset(1.3)
			histoYvsPTotal.GetXaxis().SetTitleOffset(0.8)
			histoYvsPTotal.GetYaxis().SetLabelSize(0.04)
			histoYvsPTotal.GetXaxis().SetLabelSize(0.04)
			histoYvsPTotal.GetXaxis().SetNdivisions(8)
			histoYvsPTotal.GetYaxis().SetNdivisions(32)


			c6 = TCanvas('c6', 'c6', 900, 1200)
			c6.cd()
			gPad.SetRightMargin(0.18)

			histoYvsPTotal.Draw("colz")



			leg6 = TLegend(0.65, 0.9, 0.91, 0.95)
			leg6.SetBorderSize(5)
			leg6.SetTextSize(0.022)
			leg6.AddEntry(histoYvsPTotal, "Entries: %0.2f" % (histoYvsPTotal.GetEntries()), "p")
			leg6.Draw()

			pave = TPaveText(0.18, 0.78, 0.48, 0.92, "blNDC")
			pave.SetFillColor(0)
			pave.SetBorderSize(0)
			pave.SetFillStyle(0)
			pave.SetTextAlign(12)
			pave.SetTextSize(0.025)
			t1 = pave.AddText("Particle Rates")
			t1.SetTextColor(kRed+3)
			t1.SetTextSize(0.03)
			t2 = pave.AddText(" ")
			t2.SetTextColor(kRed+3)
			t2.SetTextSize(0.03)
			pave.AddText("#gamma: %0.3f%%" % (ParticleRates[2]) )
			pave.AddText("e^{+}: %0.3f%%" % (PositronRate[str(Det)]) )
			pave.AddText("e^{-}: %0.3f%%" % (ElectronRate[str(Det)]) )
 			pave.AddText("#mu^{-}: %0.3f%%" % (MuonRate[str(Det)]) )
 			pave.AddText("Proton: %0.3f%%" % (ProtonRate[str(Det)]) )
 			pave.Draw()


			c6.SaveAs("%s/%s/YvsPTotal_%s.pdf" % (outputFolder, Det, Det))
			c6.Close()

			

			del histoYvsPTotal
			del c6
			del Tree6

	##### Print Particle Event Sizes (Pretty Table) #####
	print(table1.get_string(title="Particles"))

	print ("\033[1;31m --> Please wait, generating a histogram showing the amount of each particle by position!\033[0;0m")

	print ("\033[1;36m --> Done!\033[0;0m")

	

