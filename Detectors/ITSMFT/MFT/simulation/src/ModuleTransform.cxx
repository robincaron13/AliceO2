// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file ModuleTransform.cxx
/// \brief This macro find the transformation to apply to MFT detector elements using the initial misalignments
/// \author robin.caron@cern.ch (based on MyMuonModule macros)
/// \date 08/12/2020

/// This macro performs the misalignment on an existing MFT geometry
/// and  find the transformation to apply to MFT detector elements using the
/// initial misalignments based on the standard definition of the detector elements.
///
//-----------------------------------------------------------------------------

#include <fstream>
#include <Riostream.h>
#include <TGeoMatrix.h>
#include <TString.h>
#include <TFitter.h>
#include <TNtuple.h>
#include <TMath.h>
#include <TFile.h>
#include <iostream>
#include "TTree.h"
#include <vector>
#include <TRandom.h>

#include "TH1.h"
#include "TH2.h"
#include "TProfile2D.h"


#include "MFTSimulation/ModuleTransform.h"
#include "MFTSimulation/GeometryMisAligner.h"

#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"
#include "MFTBase/GeometryBuilder.h"

#include "MFTSimulation/Detector.h"

#include "DetectorsCommonDataFormats/DetID.h"
#include "DetectorsCommonDataFormats/AlignParam.h"

#include "Framework/Logger.h"

using namespace std;

using namespace o2::mft;
using namespace o2::detectors;

class TFile;

static void MyFcnDisk(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
static void MyFcnChip(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

ofstream myfilestored_alignparams;
std::vector<const o2::detectors::AlignParam> vecAlignParam;

ifstream myfile_idealpositions;
ifstream myfile_surveypositions;
ofstream mysurveyposition;
ifstream mysurveyposition2;
ofstream myidealposition;
ifstream myidealposition2;

ifstream myfile_idealpositionsdisk;
ifstream myfile_surveypositionsdisk;
ofstream mysurveypositiondisk;
ofstream myidealpositiondisk;

TFile *myhistos = new TFile("align_positions.root","recreate");
Int_t factor=10;
TH1D *histo_x = new TH1D("histo_x","",100,-0.01*factor,0.01*factor);
TH1D *histo_y = new TH1D("histo_y","",100,-0.05*factor,0.05*factor);
TH1D *histo_z = new TH1D("histo_z","",100,-0.1*factor,0.1*factor);
TH1D *histo_rx = new TH1D("histo_rx","",100,-2,2);
TH1D *histo_ry = new TH1D("histo_ry","",100,-2,2);
TH1D *histo_rz = new TH1D("histo_rz","",100,-2,2);
TH1D *histo_xsurvey = new TH1D("histo_xsurvey","",100,-0.01*factor,0.01*factor);
TH1D *histo_ysurvey = new TH1D("histo_ysurvey","",100,-0.05*factor,0.05*factor);
TH1D *histo_zsurvey = new TH1D("histo_zsurvey","",100,-0.1*factor,0.1*factor);
TH2D *histo2Dxy_ideal = new TH2D("histo2Dxy_ideal","",2000,-20,20,2000,-20,20);
TH2D *histo2Dxy_ladderIdeal = new TH2D("histo2Dxy_ladderIdeal","",2000,-16,-10,4000,0,12);
TH2D *histo2Dxy_ladderMisaligned = new TH2D("histo2Dxy_ladderMisaligned","",2000,-16,-10,4000,0,12);

ClassImp(o2::mft::ModuleTransform);

ModuleTransform::ModuleTransform()
  : TObject(),
    fModuleId(-1)
    , fFitter(0x0)
    , gMinuit(new TMinuit(6))
{
  // Default constructor
    /// Standard constructor
    for (Int_t i = 0; i < 3; i++) {
      for (Int_t j = 0; j < 4; j++) {
          fPadMisAlig[i][j] = 0.0;
      }
    }
    fHalfId=0;fDiskId=0;fLadderId=0;fSensorId=0;
    fsymName=nullptr;
    myfilestored_alignparams.open ("align_par.txt");
}
ModuleTransform::~ModuleTransform()
{
  // Destructor
    myfilestored_alignparams.close();
    
    myhistos->Write();
    myhistos->Close();
}

Double_t ModuleTransform::MyChi2Chip(Double_t* par)
{
  /// Returns the chisquare between local2global transform of local button targets and their surveyed position
  TGeoTranslation transTemp;
  TGeoRotation rotTemp;
  TGeoCombiTrans trfTemp;

  Double_t lChi2 = 0.;

  trfTemp.SetTranslation(transTemp);
  trfTemp.SetRotation(rotTemp);
  trfTemp.Clear();
    
  trfTemp.RotateZ(TMath::RadToDeg() * par[5]);
  trfTemp.RotateY(TMath::RadToDeg() * par[4]);
  trfTemp.RotateX(TMath::RadToDeg() * par[3]);
  trfTemp.SetTranslation(par[0], par[1], par[2]);
    
    std::vector<TGeoCombiTrans> mIdealVec;
    std::vector<TGeoCombiTrans> mAlignedVec;
    
    
      mIdealVec.clear();
      mAlignedVec.clear();
    
    TGeoCombiTrans kiChipTransformer ;
    TGeoCombiTrans kaChipTransformer ;
    
    TGeoCombiTrans kiPadTransformer ;
    TGeoCombiTrans kaPadTransformer ;
    
    Float_t x,y,z;
    
    Double_t xpadideal[4];
    Double_t ypadideal[4];
    Double_t zpadideal[4];
    Int_t nlines = 0;
   
    Int_t half,disk,ladder,sensor,pad;
    
    myidealposition2.open ("idealpositionsSensor.txt");
    while(1){
        myidealposition2 >> half >> disk >> ladder >> sensor >> pad >> x >> y >> z;
        if (!myidealposition2.good() ) break;
        if (nlines<4){
            xpadideal[nlines] = x;
            ypadideal[nlines] = y;
            zpadideal[nlines] = z;
            
            if (nlines==0) kiChipTransformer.SetTranslation(xpadideal[nlines], ypadideal[nlines], zpadideal[nlines]);

            kiPadTransformer.SetTranslation( xpadideal[nlines] - kiChipTransformer.GetTranslation()[0], ypadideal[nlines] - kiChipTransformer.GetTranslation()[1], zpadideal[nlines] - kiChipTransformer.GetTranslation()[2]);

            mIdealVec.push_back(kiPadTransformer);
        }
        nlines++;
    }
    
    Double_t xpad[4];
    Double_t ypad[4];
    Double_t zpad[4];
    nlines = 0;
    
    mysurveyposition2.open ("surveySensorPositions.txt");
    while(1){
        mysurveyposition2 >> half >> disk >> ladder >> sensor >> pad >> x >> y >> z;

        if (!mysurveyposition2.good() ) break;
        if (nlines<4){
            xpad[nlines] = x;
            ypad[nlines] = y;
            zpad[nlines] = z;

            if (nlines==0) kaChipTransformer.SetTranslation(xpadideal[nlines] +xpad[nlines], ypadideal[nlines] +ypad[nlines], zpadideal[nlines] +zpad[nlines]);
            
            kaPadTransformer.SetTranslation(xpad[nlines], ypad[nlines], zpad[nlines]);

            mAlignedVec.push_back(kaPadTransformer);
        }
        nlines++;
    }
    
    myidealposition2.close();
    mysurveyposition2.close();


    Double_t pl[3] = {0};
    Double_t pg[3] = {0};
    Double_t apg[3] = {0};

    Double_t ep[3] = {1, 1, 1};
    
    TGeoCombiTrans miPadTrf;
    TGeoCombiTrans maPadLocTrf;
    
    TGeoHMatrix miChipGloTrf = kiChipTransformer * trfTemp;
    TGeoCombiTrans miChipGlo(miChipGloTrf);
    
        
    for (int i=0; i<mIdealVec.size(); i++){
        
        miPadTrf = mIdealVec.at(i);
        
        if( mAlignedVec.size()==mIdealVec.size() ) maPadLocTrf = mAlignedVec.at(i);

        TGeoHMatrix maPadLoc2Trf = maPadLocTrf * miPadTrf ;
        TGeoCombiTrans maPadLoc(maPadLoc2Trf);
        
        TGeoHMatrix ma0ChipPadGlo = kiChipTransformer * maPadLoc ;
        TGeoCombiTrans maPadGlo(ma0ChipPadGlo);

        TGeoHMatrix miPadGloTrf = miChipGlo * miPadTrf ;
        TGeoCombiTrans miPadGlo(miPadGloTrf);

        // trfmiPadGlo.LocalToMaster(pl, pg) this could be also used
        // trfmiPadGlo.MasterToLocal(pl, pg) or this
        
        for(int k=0; k<3; k++){
            pg[k] = miPadGlo.GetTranslation()[k];
            apg[k] = maPadGlo.GetTranslation()[k];
        }

        // maPadGlo.LocalToMaster(pl, apg) this could be also used
        
        lChi2 += (pg[0] - apg[0]) * (pg[0] - apg[0]) / (ep[0] * ep[0]);
        lChi2 += (pg[1] - apg[1]) * (pg[1] - apg[1]) / (ep[1] * ep[1]);
        lChi2 += (pg[2] - apg[2]) * (pg[2] - apg[2]) / (ep[2] * ep[2]);
        
    }
  
    mIdealVec.clear();
    mAlignedVec.clear();
    
  return lChi2;
}


//_________________________________________________________________________
TGeoCombiTrans ModuleTransform::GetMatrixPadPosition() const
{
  /// Misalign given transformation and return the misaligned transformation.
  TGeoTranslation deltaTrans(fPadMisAlig[0][0], fPadMisAlig[1][0], fPadMisAlig[2][0]);
  TGeoRotation deltaRot;
  deltaRot.RotateX(0.0);
  deltaRot.RotateY(0.0);
  deltaRot.RotateZ(0.0);

  TGeoCombiTrans deltaTransf(deltaTrans, deltaRot);
  

  return TGeoCombiTrans(deltaTransf);
}


//_____________________________________________________________________________
static void MyFcnChip(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  /// Standard function as needed by Minuit-like minimization procedures.
  /// For the set of parameters par calculates and returns chi-squared.
    
  // smuggle a C++ object into a C function
  ModuleTransform* aMuonModule = (ModuleTransform*)gMinuit->GetObjectFit();

    
  f = aMuonModule->MyChi2Chip(par);
  if (iflag == 3) {
  }
  if (npar) {
  }

}


//_____________________________________________________________________________
static void MyFcnDisk(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  /// Standard function as needed by Minuit-like minimization procedures.
  /// For the set of parameters par calculates and returns chi-squared.

  // smuggle a C++ object into a C function
  ModuleTransform* aMuonModule = (ModuleTransform*)gMinuit->GetObjectFit();
    
  f = aMuonModule->MyChi2Disk(par);
  if (iflag == 3) {
  }
  if (npar) {
  }

}


//_____________________________________________________________________________
Int_t ModuleTransform::GetModuleMeanTransform(TGeoCombiTrans& modTransf, Double_t* parErr, Int_t level=4)
{
    /// Main function to obtain the misalignments from the surveyed position of the button targets;
    /// level to the last degree of freedom to align
    /// level 0 = MFT
    /// level 1 = half-MFT
    /// level 2 = half-disk
    /// level 3 = ladder
    /// level 4 = chip

    Double_t x,y,z;
    Double_t xp[1000];
    Double_t yp[1000];
    Double_t zp[1000];
    
    Double_t xpi[1000];
    Double_t ypi[1000];
    Double_t zpi[1000];
    
    Int_t nlines = 0;
    
    Double_t half,disk,sensor,pad;
    Double_t ladder=0.0;
    Double_t side=0;

    Int_t pflag=0;
    
    if(level == 4){
        myfile_idealpositions.open ("idealpositions.txt");
        while( 1){
            myfile_idealpositions >> half>> disk >> ladder >> sensor >> pad >> x >> y >> z;
            if (!myfile_idealpositions.good() ) break;
            if ( (Int_t(half)==fHalfId) && (Int_t(disk)==fDiskId) && (Int_t(ladder)==fLadderId) && (Int_t(sensor)==fSensorId) ) {
                xpi[pflag]=x;
                ypi[pflag]=y;
                zpi[pflag]=z;
                pflag++;
            }
            nlines++;
        }
        
        pflag=0;
        nlines = 0;
        
        myfile_surveypositions.open ("survey_true_padpositions.txt");
        while( 1){
            myfile_surveypositions >> half>> disk >> ladder >> sensor >> pad >> x >> y >> z;
            if (!myfile_surveypositions.good() ) break;
            if ( (Int_t(half)==fHalfId) && (Int_t(disk)==fDiskId) && (Int_t(ladder)==fLadderId) && (Int_t(sensor)==fSensorId) ) {
                xp[pflag]=x;
                yp[pflag]=y;
                zp[pflag]=z;
                pflag++;
            }
            nlines++;
        }
        
        mysurveyposition.open ("surveySensorPositions.txt");
        for(int k=0; k<4;k++){
            mysurveyposition << fHalfId <<" "<<fDiskId <<" " << fLadderId <<" " <<fSensorId<<" ";
            mysurveyposition<< k <<" "<<xp[k] <<" "<< yp[k] <<" "<< zp[k] ; //write to file
            mysurveyposition <<"\n";
            
            histo_xsurvey->Fill(xp[k]);
            histo_ysurvey->Fill(yp[k]);
            histo_zsurvey->Fill(zp[k]);
        }
        mysurveyposition.close();
        
        myidealposition.open ("idealpositionsSensor.txt");
        for(int k=0; k<4;k++){
            myidealposition << fHalfId <<" "<<fDiskId <<" " << fLadderId <<" " <<fSensorId<<" ";
            myidealposition<< k <<" "<< xpi[k] <<" "<< ypi[k] <<" "<< zpi[k] ; //write to file
            myidealposition <<"\n";
            
            histo2Dxy_ideal->Fill(xpi[k],ypi[k]);
            if ( (fHalfId==1) && (fDiskId==4) && (fLadderId>16) ){ histo2Dxy_ladderIdeal->Fill(xpi[k],ypi[k]);
            }
            Double_t xpad0=xpi[k]+xp[k];
            Double_t ypad0=ypi[k]+yp[k];
            Double_t zpad0=(zpi[k]+zp[k]);

            if ( (fHalfId==1) && (fDiskId==4) && (fLadderId>16) ){ histo2Dxy_ladderMisaligned->Fill(xpad0,ypad0);
            }
        }
        myidealposition.close();
        
        myfile_surveypositions.close();
        myfile_idealpositions.close();

    }
    
    if(level == 2){
        nlines=0;
        
        myfile_idealpositionsdisk.open ("idealpositionsdisk.txt");
        while( 1){
            myfile_idealpositionsdisk >> half>> disk >> side >> x >> y >> z;
            if (!myfile_idealpositionsdisk.good() ) break;
            if ( (Int_t(half)==fHalfId) && (Int_t(disk)==fDiskId) && (Int_t(side)==0) ) {
                xpi[0]=x;
                ypi[0]=y;
                zpi[0]=z;
            }
            if ( (Int_t(half)==fHalfId) && (Int_t(disk)==fDiskId) && (Int_t(side)==1) ) {
                xpi[1]=x;
                ypi[1]=y;
                zpi[1]=z;
            }
            nlines++;

        }
        nlines=0;
        myfile_surveypositionsdisk.open ("surveypositionsdisk_afterP2.txt");
        while( 1){
            myfile_surveypositionsdisk >> half>> disk >> side >> x >> y >> z;
            if (!myfile_surveypositionsdisk.good() ) break;
            if ( (Int_t(half)==fHalfId) && (Int_t(disk)==fDiskId) && (Int_t(side)==0) ) {
                xp[0]=x;
                yp[0]=y;
                zp[0]=z;
            }
            if ( (Int_t(half)==fHalfId) && (Int_t(disk)==fDiskId) && (Int_t(side)==1) ) {
                xp[1]=x;
                yp[1]=y;
                zp[1]=z;
            }
            nlines++;
        }
        
        myidealpositiondisk.open ("idealpositionsdiskonly.txt");
        for(int k=0; k<2;k++){
            myidealpositiondisk << fHalfId <<" "<<fDiskId <<" ";
            myidealpositiondisk<< k <<" "<<xpi[k] <<" "<< ypi[k] <<" "<< zpi[k] ; //write to file
            myidealpositiondisk <<"\n";
            histo2Dxy_ideal->Fill(xpi[k],ypi[k]);
        }
        mysurveypositiondisk.open ("surveypositionsdiskonly.txt");
        for(int k=0; k<2;k++){
            mysurveypositiondisk << fHalfId <<" "<<fDiskId <<" ";
            mysurveypositiondisk << k <<" "<<xp[k] <<" "<< yp[k] <<" "<< zp[k] ; //write to file
            mysurveypositiondisk <<"\n";
            histo_xsurvey->Fill(xp[k]);
            histo_ysurvey->Fill(yp[k]);
            histo_zsurvey->Fill(zp[k]);
            
            if ( (Int_t(half)==0) && (Int_t(disk)==4) && (Int_t(side)==0) ) histo2Dxy_ladderIdeal->Fill(xpi[k],ypi[k]);
            Double_t xpad0=xpi[k]+xp[k];
            Double_t ypad0=ypi[k]+yp[k];
            Double_t zpad0=(zpi[k]+zp[k]);

            if ( (Int_t(half)==0) && (Int_t(disk)==4) && (Int_t(side)==0) ) histo2Dxy_ladderMisaligned->Fill(xpad0,ypad0);

        }
        myidealpositiondisk.close();
        mysurveypositiondisk.close();

        myfile_idealpositionsdisk.close();
        myfile_surveypositionsdisk.close();

    }
    
    Int_t value = 1;
    Double_t partest[20] = {0};
    
    if(level == 4) gMinuit->SetFCN( MyFcnChip );
    if(level == 2) gMinuit->SetFCN( MyFcnDisk );

    // Vector of step, initial min and max value
    Double_t vstrt[6]={0.0};
    Double_t stp[6];
    Double_t bmin[6];
    Double_t bmax[6];
    
    //translation transformation
    vstrt[0] = modTransf.GetTranslation()[0];
    vstrt[1] = modTransf.GetTranslation()[1];
    vstrt[2] = modTransf.GetTranslation()[2];
    stp[0]= stp[1] =stp[2] =0.0001;
    bmin[0] = bmin[1] =bmin[2] =vstrt[0] -1.0;
    bmax[0] =bmax[1] =bmax[2] = vstrt[0] +1.0;
        
    //rotation transformation
    vstrt[3] = 0;
    vstrt[4] = 0;
    vstrt[5] = 0;
    stp[3] = stp[4] = stp[5]= 0.0001;
    bmin[3] = bmin[4]= bmin[5]= vstrt[3] - 0.5;
    bmax[3] = bmax[4]= bmax[5]= vstrt[3] + 0.5;

    Double_t arglist[10];
    Int_t ierflg = 0;
    
    gMinuit->SetPrintLevel(-1);
    
    // initialize the 6 minuit alignment parameters
    gMinuit->mnparm(0, "dx", vstrt[0], stp[0], bmin[0],bmax[0], ierflg);
    gMinuit->mnparm(1, "dy", vstrt[1], stp[1], bmin[1],bmax[1], ierflg);
    gMinuit->mnparm(2, "dz", vstrt[2], stp[2], bmin[2],bmax[2], ierflg);
    gMinuit->mnparm(3, "rx", vstrt[3], stp[3], bmin[3],bmax[3], ierflg);
    gMinuit->mnparm(4, "ry", vstrt[4], stp[4], bmin[4],bmax[4], ierflg);
    gMinuit->mnparm(5, "rz", vstrt[5], stp[5], bmin[5],bmax[5], ierflg);
    
    gMinuit->SetErrorDef(1);

    gMinuit->mnexcm("SET NOW", arglist ,1,ierflg);
    
    arglist[0] = 2;
  
    gMinuit->mnexcm("SET STR", arglist, 1,ierflg);
    
    // since only 2 pin/alignment markers per disk
    // produces no effects of a rotation on y
    if(level == 2) gMinuit->FixParameter(3);
    
    arglist[0] = 300;
    gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);
    //we could also call: gMinuit->mnexcm("MINIMIZE", arglist,1,ierflg)

    arglist[0] = 0;
    //we could also call: gMinuit->mnexcm("IMPROVE", arglist,1,ierflg)
    gMinuit->mnexcm("SCAN", arglist,1,ierflg);

    
    Double_t parlist[10]={0.0};
    Double_t errlist[10]={0.0};

    for (Int_t iPar = 0; iPar < 6; iPar++) {
      parErr[iPar] = errlist[iPar];
    }
    
  modTransf.Clear();
    
  modTransf.RotateZ( TMath::RadToDeg() * gMinuit->GetParameter(5,parlist[5],errlist[5]));
  modTransf.RotateY( TMath::RadToDeg() * gMinuit->GetParameter(4,parlist[4],errlist[4]));
  modTransf.RotateX( TMath::RadToDeg() * gMinuit->GetParameter(3,parlist[3],errlist[3]));
    
  modTransf.SetTranslation(gMinuit->GetParameter(0,parlist[0],errlist[0]),
                           gMinuit->GetParameter(1,parlist[1],errlist[1]),
                           gMinuit->GetParameter(2,parlist[2],errlist[2]));

    
    Int_t ipad=0;
    
    /*
    for (int ipad = 0; ipad < 4; ipad++){
        myfilestored_alignparams << fHalfId <<" "<<fDiskId <<" "<< fLadderId <<" " <<fSensorId<<" "<<ipad <<" ";

//        for (int j = 0; j < 3; j++){
//          //myfilestored_alignparams<< ipad <<" " << parlist[j] ; //write to file
//            myfilestored_alignparams <<vstrt[j]<<" " ; //write to file
//
//        }

        if(ipad==0) myfilestored_alignparams <<gRandom->Gaus(0.0,0.0001)<<" "<<gRandom->Gaus(0.0,0.0001)<<" "<<gRandom->Gaus(0.0,0.0001) <<" ";
        if(ipad==1) myfilestored_alignparams <<gRandom->Gaus(0.0,0.0001)<<" "<<gRandom->Gaus(0.0,0.0001)<<" "<<gRandom->Gaus(0.0,0.0001) <<" ";
        if(ipad==2) myfilestored_alignparams <<gRandom->Gaus(0.0,0.0001)<<" "<<gRandom->Gaus(1.5*0.9848- 1.5,0.0001) <<" "<<gRandom->Gaus(1.5*0.1736,0.0001) <<" ";
        if(ipad==3) myfilestored_alignparams <<gRandom->Gaus(0.0,0.0001)<<" "<<gRandom->Gaus(1.5*0.9848- 1.5,0.0001)<<" "<<gRandom->Gaus(1.5*0.1736,0.0001) <<" ";

        myfilestored_alignparams <<"\n";
    }
    */


    myfilestored_alignparams << fHalfId <<" "<<fDiskId <<" "<< fLadderId <<" " <<fSensorId<<" ";
    for (int j = 0; j < 3; j++){
      myfilestored_alignparams<< parlist[j]<< " " ; //write to file
    }
    for (int j = 3; j < 6; j++){
        myfilestored_alignparams<< TMath::RadToDeg() * parlist[j] << " "; //write to file
    }
    myfilestored_alignparams <<"\n";
    

    histo_x->Fill(parlist[0]);
    histo_y->Fill(parlist[1]);
    histo_z->Fill(parlist[2]);
    histo_rx->Fill(TMath::RadToDeg() * parlist[3]);
    histo_ry->Fill(TMath::RadToDeg() * parlist[4]);
    histo_rz->Fill(TMath::RadToDeg() * parlist[5]);
        
    
    TGeoCombiTrans mAlignParams;
    mAlignParams.RotateZ(TMath::RadToDeg() * parlist[5]);
    mAlignParams.RotateY(TMath::RadToDeg() * parlist[4]);
    mAlignParams.RotateX(TMath::RadToDeg() * parlist[3]);
    mAlignParams.SetTranslation(parlist[0], parlist[1], parlist[2]);
        
    
    o2::detectors::AlignParam lAlignParams;
    lAlignParams.setRotation(TMath::RadToDeg() * parlist[3], TMath::RadToDeg() * parlist[4], TMath::RadToDeg() * parlist[5]);
    lAlignParams.setTranslation(parlist[0], parlist[1], parlist[2]);
    lAlignParams.setSymName(fsymName);
    
    // Store the AlignParam object in vector (using the obtained alignment parameters from minimization)
    vecAlignParam.push_back(lAlignParams);
    
    //LOG(INFO) << " transf.par x = "<<parlist[0]<<" +/- "<<parErr[0]<< "   -> fill align par file " ;
    
  return 1;
}

void ModuleTransform::PrintAlignResults(){
    
    for (o2::detectors::AlignParam AlignPars : vecAlignParam) {
        AlignPars.Print();
    }
}




Double_t ModuleTransform::MyChi2Disk(Double_t* par)
{
  /// Returns the chisquare between local2global transform of local button targets and their surveyed position
  TGeoTranslation transTemp;
  TGeoRotation rotTemp;
  TGeoCombiTrans trfTemp;

  Double_t lChi2 = 0.;

  trfTemp.SetTranslation(transTemp);
  trfTemp.SetRotation(rotTemp);
  trfTemp.Clear();
    
  trfTemp.RotateZ(TMath::RadToDeg() * par[5]);
  trfTemp.RotateY(TMath::RadToDeg() * par[4]);
  trfTemp.RotateX(TMath::RadToDeg() * par[3]);
  trfTemp.SetTranslation(par[0], par[1], par[2]);
    
    std::vector<TGeoCombiTrans> mIdealVec;
    std::vector<TGeoCombiTrans> mAlignedVec;
    
    
      mIdealVec.clear();
      mAlignedVec.clear();
    
    TGeoCombiTrans kiChipTransformer ;
    TGeoCombiTrans kaChipTransformer ; // to be changed
    
    TGeoCombiTrans kiPadTransformer ;
    TGeoCombiTrans kaPadTransformer ; // to be changed
    
    Float_t x,y,z;
    
    Double_t xpadideal[4];
    Double_t ypadideal[4];
    Double_t zpadideal[4];
    Int_t nlines = 0;
   
    Int_t half,disk,ladder,sensor,pad, side;
    
    myidealposition2.open ("idealpositionsdiskonly.txt");
    while(1){
        myidealposition2 >> half >> disk >> side >> x >> y >> z;
        if (!myidealposition2.good() ) break;
        if (nlines<4){
            xpadideal[nlines] = x;
            ypadideal[nlines] = y;
            zpadideal[nlines] = z;
            
            if (nlines==0) kiChipTransformer.SetTranslation(0, 0, zpadideal[nlines]);

            kiPadTransformer.SetTranslation( xpadideal[nlines] - kiChipTransformer.GetTranslation()[0], ypadideal[nlines] - kiChipTransformer.GetTranslation()[1], zpadideal[nlines] - kiChipTransformer.GetTranslation()[2]);

            mIdealVec.push_back(kiPadTransformer);
        }
        nlines++;
    }
    
    Double_t xpad[4];
    Double_t ypad[4];
    Double_t zpad[4];
    nlines = 0;
    
    mysurveyposition2.open ("surveypositionsdiskonly.txt");
    while(1){
        mysurveyposition2 >> half >> disk >> side  >> x >> y >> z;
        if (!mysurveyposition2.good() ) break;
        if (nlines<4){
            xpad[nlines] = x;
            ypad[nlines] = y;
            zpad[nlines] = z;

            if (nlines==0) kaChipTransformer.SetTranslation(xpadideal[nlines] +xpad[nlines], ypadideal[nlines] +ypad[nlines], zpadideal[nlines] +zpad[nlines]);

            kaPadTransformer.SetTranslation(xpad[nlines], ypad[nlines], zpad[nlines]);

            mAlignedVec.push_back(kaPadTransformer);
        }
        nlines++;
    }
    
    myidealposition2.close();
    mysurveyposition2.close();


    Double_t pl[3] = {0};
    Double_t pg[3] = {0};
    Double_t apg[3] = {0};

    Double_t ep[3] = {1, 1, 1};
    
    TGeoCombiTrans miPadTrf;
    TGeoCombiTrans maPadLocTrf;
    
    TGeoHMatrix miChipGloTrf = kiChipTransformer * trfTemp;
    TGeoCombiTrans miChipGlo(miChipGloTrf);
    
        
    for (int i=0; i<mIdealVec.size(); i++){
        
        miPadTrf = mIdealVec.at(i);
        if( mAlignedVec.size()==mIdealVec.size() ) maPadLocTrf = mAlignedVec.at(i);

        TGeoHMatrix maPadLoc2Trf = maPadLocTrf * miPadTrf ;
        TGeoCombiTrans maPadLoc(maPadLoc2Trf);
        
        TGeoHMatrix ma0ChipPadGlo = kiChipTransformer * maPadLoc ;
        TGeoCombiTrans maPadGlo(ma0ChipPadGlo);

        TGeoHMatrix miPadGloTrf = miChipGlo * miPadTrf ;
        TGeoCombiTrans miPadGlo(miPadGloTrf);

        for(int k=0; k<3; k++){
            pg[k] = miPadGlo.GetTranslation()[k];
            apg[k] = maPadGlo.GetTranslation()[k];
        }
        
        lChi2 += (pg[0] - apg[0]) * (pg[0] - apg[0]) / (ep[0] * ep[0]);
        lChi2 += (pg[1] - apg[1]) * (pg[1] - apg[1]) / (ep[1] * ep[1]);
        lChi2 += (pg[2] - apg[2]) * (pg[2] - apg[2]) / (ep[2] * ep[2]);
            
    }
  
    mIdealVec.clear();
    mAlignedVec.clear();
    
  return lChi2;
}
