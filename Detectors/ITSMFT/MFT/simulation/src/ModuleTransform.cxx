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
#include <TMinuit.h>
#include <TMath.h>

//#include "AliLog.h"
//
//#include "AliMUONGeometryTransformer.h"
//#include "AliMUONGeometryModuleTransformer.h"
//#include "AliMUONGeometryDetElement.h"
//
//#include "AliMpExMap.h"
//#include "AliMpExMapIterator.h"

#include "MFTSimulation/ModuleTransform.h"
#include "MFTSimulation/GeometryMisAligner.h"

#include "Framework/Logger.h"

using namespace o2::mft;
//using namespace o2::detectors;
class TFitter;

void MyFcn(int& npar, double* g, double& f, double* par, int iflag);


ClassImp(o2::mft::ModuleTransform);

ModuleTransform::ModuleTransform()
  : TObject(),
    fModuleId(-1)
    // fIdealGeoTransformer(0x0),
    // fAlignGeoTransformer(0x0),
    // detElements(0x0),
    , fFitter(0x0)
{
  // Default constructor
}
ModuleTransform::~ModuleTransform()
{
  // Destructor
}

Double_t ModuleTransform::MyChi2(Double_t* par)
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

  //cout << "MyChi2 trfTemp: " << endl;
  //	trfTemp.Print();

  // module transformers
  //	Int_t oiMt = fModuleId%16;
  Int_t oiMt = fModuleId % 20;
  //	const AliMUONGeometryModuleTransformer *kiModuleTransformer = fIdealGeoTransformer->GetModuleTransformer(oiMt, true);
  //	const AliMUONGeometryModuleTransformer *kaModuleTransformer = fAlignGeoTransformer->GetModuleTransformer(oiMt, true);
  TGeoCombiTrans kiModuleTransformer ;
  TGeoCombiTrans kaModuleTransformer ; // to be changed

  TGeoCombiTrans iModTrf = kiModuleTransformer;
  //TGeoCombiTrans iModTrf = TGeoCombiTrans(*kiModuleTransformer->GetTransformation());
  //cout << "MyChi2 iModTrf: " << endl;
  //	iModTrf.Print();

  TGeoHMatrix matModGlo = iModTrf * trfTemp;
  TGeoCombiTrans trfModGlo(matModGlo);
  //cout << "MyChi2 : trfModGlo" << endl;

  Double_t pl[3] = {0};
  Double_t pg[3] = {0};
  Double_t apg[3] = {0};

  Double_t ep[3] = {0.015, 0.015, 0.015};

  //	AliMpExMap *detElements = kiModuleTransformer->GetDetElementStore();
  //	TIter next(detElements->CreateIterator());
  //	AliMUONGeometryDetElement* idetElement;
  //	AliMUONGeometryDetElement* adetElement;
  Int_t iDe(-1);
  TString lDetElemName;

  TGeoCombiTrans iDELocTrf ;
  //TGeoHMatrix iDELocTrf = lAP.createMatrix();  //create local matrix transformation

  TGeoHMatrix matDEGlo = trfModGlo * iDELocTrf;
  TGeoCombiTrans trfDEGlo(matDEGlo);

  trfDEGlo.LocalToMaster(pl, pg);
  //        kaModuleTransformer->Local2Global(idetElement->GetId(),pl[0], pl[1], pl[2], apg[0], apg[1], apg[2]);
  kaModuleTransformer.LocalToMaster(pl, apg);

  lChi2 += (pg[0] - apg[0]) * (pg[0] - apg[0]) / (ep[0] * ep[0]);
  lChi2 += (pg[1] - apg[1]) * (pg[1] - apg[1]) / (ep[1] * ep[1]);
  lChi2 += (pg[2] - apg[2]) * (pg[2] - apg[2]) / (ep[2] * ep[2]);

  //    o2::detectors::AlignParam lAP;
  //
  //    Int_t nHalf = mGeometryTGeo->getNumberOfHalfs();
  //    TString sname;
  //    for (Int_t hf = 0; hf < nHalf; hf++) {
  //
  //      Int_t nDisks = mGeometryTGeo->getNumberOfDisksPerHalf(hf);
  //
  //        for (Int_t dk = 0; dk < nDisks; dk++) {
  //
  //            Int_t nLadders = 0;
  //
  //            for (Int_t sensor = mGeometryTGeo->getMinSensorsPerLadder(); sensor < mGeometryTGeo->getMaxSensorsPerLadder() + 1; sensor++) {
  //              nLadders += mGeometryTGeo->getNumberOfLaddersPerDisk(hf, dk, sensor);
  //            }
  //
  //            for (Int_t lr = 0; lr < nLadders; lr++) {
  //
  //                Int_t nSensorsPerLadder = mGeometryTGeo->getNumberOfSensorsPerLadder(hf, dk, lr);
  //
  //                for (Int_t sr = 0; sr < nSensorsPerLadder; sr++) {
  //
  //                    sname = mGeometryTGeo->composeSymNameChip(hf, dk, lr, sr);
  //
  //                    lAP.setSymName(sname);
  //
  //                    lAP.setLocalParams(localDeltaTransform);
  //
  //                    //TGeoCombiTrans iDELocTrf = MisAlignSensor();
  //                    TGeoHMatrix iDELocTrf = lAP.createMatrix();  //create local matrix transformation
  //
  //                    TGeoHMatrix matDEGlo = trfModGlo*iDELocTrf;
  //                    TGeoCombiTrans trfDEGlo(matDEGlo);
  //
  //                    trfDEGlo.LocalToMaster(pl, pg);
  //            //        kaModuleTransformer->Local2Global(idetElement->GetId(),pl[0], pl[1], pl[2], apg[0], apg[1], apg[2]);
  //
  //
  //                    lChi2 += (pg[0]-apg[0])*(pg[0]-apg[0])/(ep[0]*ep[0]);
  //                    lChi2 += (pg[1]-apg[1])*(pg[1]-apg[1])/(ep[1]*ep[1]);
  //                    lChi2 += (pg[2]-apg[2])*(pg[2]-apg[2])/(ep[2]*ep[2]);
  //
  //                }
  //
  //            }
  //
  //
  //        }
  //
  //    }

  //
  //	while ( ( idetElement = static_cast<AliMUONGeometryDetElement*>(next()) ) )
  //	{
  //// 		if (fModuleId<4 && ((idetElement->GetId()%100)==1 || (idetElement->GetId()%100)==2)) {
  //// 			continue;
  //// 		}
  //// 		if (fModuleId>15 && ((idetElement->GetId()%100)==0 || (idetElement->GetId()%100)==3)) {
  //// 			continue;
  //// 		}
  //		++iDe;
  //
  //		// local transformation of this detection element.
  //		//TGeoCombiTrans iDELocTrf = TGeoCombiTrans(*idetElement->GetLocalTransformation());
  //        TGeoCombiTrans iDELocTrf = MisAlignSensor();
  //
  //		TGeoHMatrix matDEGlo = trfModGlo*iDELocTrf;
  //		TGeoCombiTrans trfDEGlo(matDEGlo);
  //
  //		cout << "MyChi2 : trfDEGlo" << endl;
  ////		trfDEGlo.Print();
  //
  //		trfDEGlo.LocalToMaster(pl, pg);
  ////		kaModuleTransformer->Local2Global(idetElement->GetId(),pl[0], pl[1], pl[2], apg[0], apg[1], apg[2]);
  //
  //
  //
  ////		printf("pl: %f,%f,%f\n",pl[0], pl[1], pl[2]);
  ////		printf("pg: %f,%f,%f\n",pg[0], pg[1], pg[2]);
  ////		printf("apg: %f,%f,%f\n",apg[0], apg[1], apg[2]);
  //
  //		lChi2 += (pg[0]-apg[0])*(pg[0]-apg[0])/(ep[0]*ep[0]);
  //		lChi2 += (pg[1]-apg[1])*(pg[1]-apg[1])/(ep[1]*ep[1]);
  //		lChi2 += (pg[2]-apg[2])*(pg[2]-apg[2])/(ep[2]*ep[2]);
  //
  ////		if (fModuleId>3 && fModuleId<16) {
  //		if (fModuleId>3 && fModuleId<20) {
  //			pl[0]=-iDELocTrf.GetTranslation()[0];
  //			pl[1]=-20;
  //		} else {
  //			pl[0]=100.;
  //			pl[1]=0.;
  //		}
  //		trfDEGlo.LocalToMaster(pl, pg);
  //		kaModuleTransformer->Local2Global(idetElement->GetId(),
  //																			pl[0], pl[1], pl[2],
  //																			apg[0], apg[1], apg[2]);
  //
  //		//		printf("pl: %f,%f,%f\n",pl[0], pl[1], pl[2]);
  //		//		printf("pg: %f,%f,%f\n",pg[0], pg[1], pg[2]);
  //		//		printf("apg: %f,%f,%f\n",apg[0], apg[1], apg[2]);
  //
  //		lChi2 += (pg[0]-apg[0])*(pg[0]-apg[0])/(ep[0]*ep[0]);
  //		lChi2 += (pg[1]-apg[1])*(pg[1]-apg[1])/(ep[1]*ep[1]);
  //		lChi2 += (pg[2]-apg[2])*(pg[2]-apg[2])/(ep[2]*ep[2]);
  //
  ////		if (fModuleId>3 && fModuleId<16) {
  //		if (fModuleId>3 && fModuleId<20) {
  //			pl[0]=+iDELocTrf.GetTranslation()[0];
  //			pl[1]=+20;
  //		} else {
  //			pl[0]=0.;
  //			pl[1]=100.;
  //		}
  //		trfDEGlo.LocalToMaster(pl, pg);
  //		kaModuleTransformer->Local2Global(idetElement->GetId(),
  //																			pl[0], pl[1], pl[2],
  //																			apg[0], apg[1], apg[2]);
  //
  //		//		printf("pl: %f,%f,%f\n",pl[0], pl[1], pl[2]);
  //		//		printf("pg: %f,%f,%f\n",pg[0], pg[1], pg[2]);
  //		//		printf("apg: %f,%f,%f\n",apg[0], apg[1], apg[2]);
  //
  //		lChi2 += (pg[0]-apg[0])*(pg[0]-apg[0])/(ep[0]*ep[0]);
  //		lChi2 += (pg[1]-apg[1])*(pg[1]-apg[1])/(ep[1]*ep[1]);
  //		lChi2 += (pg[2]-apg[2])*(pg[2]-apg[2])/(ep[2]*ep[2]);
  //  }

  //cout << "MyChi2 : lChi2: " << lChi2 << endl;

  return lChi2;
}

//_____________________________________________________________________________
void MyFcn(int& npar, double* g, double& f, double* par, int iflag)
{
  ///
  /// Standard function as needed by Minuit-like minimization procedures.
  /// For the set of parameters par calculates and returns chi-squared.
  ///
  TMinuit *gMinuit;  //initialize TMinuit

  // smuggle a C++ object into a C function
  ModuleTransform* aMuonModule = (ModuleTransform*)gMinuit->GetObjectFit();

  f = aMuonModule->MyChi2(par);
  if (iflag == 3) {
  }
  if (npar) {
  }
  if (g) {
  } // no warnings about unused stuff...
}

//_____________________________________________________________________________
Int_t ModuleTransform::GetModuleMeanTransform(TGeoCombiTrans& modTransf, Double_t* parErr)
{
  /// Main function to obtain the misalignments from the surveyed position of the button targets;

  //void (ModuleTransform::*func)(int& npar, double* g, double& f, double* par, int iflag);
  //func = &ModuleTransform::MyFcn;
    
//    LKFMinuit * ModuleTransform = new LKFMinuit();
    //double  LKFMinuit::minuitFunction( const double * x);
    //int ndim = 6;
    //ROOT::Math::Functor fcn(ndim, ModuleTransform, &LKFMinuit::minuitFunction);
    //ROOT::Fit::Fitter fitter;
    
  //TFitter fitter(100);
    TMinuit *gMinuit;
    
  gMinuit->SetObjectFit(this);
  fFitter->SetFCN( MyFcn );

  fFitter->SetParameter(0, "dx", 0, 0.01, -2, 2);
  fFitter->SetParameter(1, "dy", 0, 0.01, -2, 2);
  fFitter->SetParameter(2, "dz", 0, 0.01, -5, 5);
  fFitter->SetParameter(3, "rx", 0, 0.0001, -0.05, 0.05);
  fFitter->SetParameter(4, "ry", 0, 0.0001, -0.05, 0.05);

  //  fFitter->SetParameter(3,"rx",0,0.0001,0,0);
  //  fFitter->SetParameter(4,"ry",0,0.0001,0,0);
  //  fFitter->SetParameter(5,"rz",0,0.0001,0,0);
  fFitter->SetParameter(5, "rz", 0, 0.0001, -0.05, 0.05);

  double arglist[20];
  arglist[0] = 2;
  fFitter->ExecuteCommand("SET PRINT", arglist, 1);
  fFitter->ExecuteCommand("SET ERR", arglist, 1);
  arglist[0] = 0;
  //fFitter->ExecuteCommand("SIMPLEX", arglist, 1);
  //  fFitter->ExecuteCommand("MINIMIZE", arglist, 1);
  fFitter->ExecuteCommand("MIGRAD", arglist, 1);
  fFitter->ExecuteCommand("IMPROVE", arglist, 1);
  //  fFitter->ExecuteCommand("MINOS", arglist, 1);
  //  fFitter->ExecuteCommand("CALL 3", arglist,0);

  for (int j = 0; j < 3; j++)
    printf("%10.3f ", fFitter->GetParameter(j));
  for (int j = 3; j < 6; j++)
    printf("%10.3f ", 1000 * fFitter->GetParameter(j));
  printf("\n");
  for (int j = 0; j < 3; j++)
    printf("%10.3f ", fFitter->GetParError(j));
  for (int j = 3; j < 6; j++)
    printf("%10.3f ", 1000 * fFitter->GetParError(j));
  printf("\n");

  modTransf.Clear();
  modTransf.RotateZ(TMath::RadToDeg() * fFitter->GetParameter(5));
  modTransf.RotateY(TMath::RadToDeg() * fFitter->GetParameter(4));
  modTransf.RotateX(TMath::RadToDeg() * fFitter->GetParameter(3));
  modTransf.SetTranslation(fFitter->GetParameter(0), fFitter->GetParameter(1), fFitter->GetParameter(2));

  for (Int_t iPar = 0; iPar < 6; iPar++) {
    parErr[iPar] = fFitter->GetParError(iPar);
  }

  return 1;
}
