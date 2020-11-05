// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GeometryMisAligner.cxx
/// \brief This macro performs the misalignment on an existing MFT geometry
/// \author robin.caron@cern.ch (based on MUON/MCH AliRoot macros)
/// \date 01/07/2020

/// This macro performs the misalignment on an existing MFT geometry
/// based on the standard definition of the detector elements in
/// the MFTGeometryTransformer class.
///
/// It uses GeometryMisAligner :
/// - Creates a new MFTGeometryTransformer and GeometryMisAligner
/// - Loads the geometry from the specified geometry file (default is geometry.root)
/// - Creates a second MFTGeometryTransformer by misaligning the existing
///   one using GeometryMisAligner::MisAlign
/// - User has to specify the magnitude of the alignments, in the Cartesian
///   co-ordiantes (which are used to apply translation misalignments) and in the
///   spherical co-ordinates (which are used to apply angular displacements)
/// - User can also set misalignment ranges by hand using the methods :
///   SetMaxCartMisAlig, SetMaxAngMisAlig, SetXYAngMisAligFactor
///   (last method takes account of the fact that the misalingment is greatest in
///   the XY plane, since the detection elements are fixed to a support structure
///   in this plane. Misalignments in the XZ and YZ plane will be very small
///   compared to those in the XY plane, which are small already - of the order
///   of microns)
/// - Default behavior generates a "residual" misalignment using gaussian
///   distributions. Uniform distributions can still be used, see
///   GeometryMisAligner.
/// - User can also generate module misalignments using SetModuleCartMisAlig
///   and SetModuleAngMisAlig
///
/// Note: If the detection elements are allowed to be misaligned in all
/// directions, this has consequences for the alignment algorithm, which
/// needs to know the number of free parameters. Eric only allowed 3 :
/// x,y,theta_xy, but in principle z and the other two angles are alignable
/// as well.
///
//-----------------------------------------------------------------------------

#include "MFTSimulation/GeometryMisAligner.h"

#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"
#include "MFTBase/GeometryBuilder.h"

#include "MFTSimulation/Detector.h"

#include "DetectorsCommonDataFormats/DetID.h"
#include "DetectorsCommonDataFormats/AlignParam.h"

#include "MFTSimulation/ModuleTransform.h"

//#include "MFTGeometryTransformer.h"

//#include "MFTGeometryModuleTransformer.h"
//#include "MFTGeometryDetElement.h"

//#include "AliMpExMap.h"
//#include "AliMpExMapIterator.h"
//
//#include "AliAlignObjMatrix.h"
//#include "AliMathBase.h"
//#include "AliLog.h"

#include <TClonesArray.h>
#include <TGeoMatrix.h>
#include <TMatrixDSym.h>
#include <TMath.h>
#include <TRandom.h>
#include <Riostream.h>
//#include <TObject.h>
#include <vector>

#include "Framework/Logger.h"

using namespace o2::mft;
using namespace o2::detectors;

ClassImp(o2::mft::GeometryMisAligner);

//______________________________________________________________________________
GeometryMisAligner::GeometryMisAligner(Double_t cartXMisAligM, Double_t cartXMisAligW, Double_t cartYMisAligM, Double_t cartYMisAligW, Double_t angMisAligM, Double_t angMisAligW)
  : TObject(),
    fUseUni(kFALSE),
    fUseGaus(kTRUE),
    fXYAngMisAligFactor(0.0),
    fZCartMisAligFactor(0.0)
{
  /// Standard constructor
  for (Int_t i = 0; i < 6; i++) {
    for (Int_t j = 0; j < 2; j++) {
      fDetElemMisAlig[i][j] = 0.0;
      fModuleMisAlig[i][j] = 0.0;
      fSensorMisAlig[i][j] = 0.0;
      fHalfMisAlig[i][j] = 0.0;
    }
  }
  fDetElemMisAlig[0][0] = cartXMisAligM;
  fDetElemMisAlig[0][1] = cartXMisAligW;
  fDetElemMisAlig[1][0] = cartYMisAligM;
  fDetElemMisAlig[1][1] = cartYMisAligW;
  fDetElemMisAlig[5][0] = angMisAligM;
  fDetElemMisAlig[5][1] = angMisAligW;
}

//______________________________________________________________________________
GeometryMisAligner::GeometryMisAligner(Double_t cartMisAligM, Double_t cartMisAligW, Double_t angMisAligM, Double_t angMisAligW)
  : TObject(),
    fUseUni(kFALSE),
    fUseGaus(kTRUE),
    fXYAngMisAligFactor(0.0),
    fZCartMisAligFactor(0.0)
{
  /// Standard constructor
  for (Int_t i = 0; i < 6; i++) {
    for (Int_t j = 0; j < 2; j++) {
      fDetElemMisAlig[i][j] = 0.0;
      fModuleMisAlig[i][j] = 0.0;
      fSensorMisAlig[i][j] = 0.0;
      fHalfMisAlig[i][j] = 0.0;
    }
  }
  fDetElemMisAlig[0][0] = cartMisAligM;
  fDetElemMisAlig[0][1] = cartMisAligW;
  fDetElemMisAlig[1][0] = cartMisAligM;
  fDetElemMisAlig[1][1] = cartMisAligW;
  fDetElemMisAlig[5][0] = angMisAligM;
  fDetElemMisAlig[5][1] = angMisAligW;
}

//______________________________________________________________________________
GeometryMisAligner::GeometryMisAligner(Double_t cartMisAlig, Double_t angMisAlig)
  : TObject(),
    fUseUni(kTRUE),
    fUseGaus(kFALSE),
    fXYAngMisAligFactor(0.0),
    fZCartMisAligFactor(0.0)
{
  /// Standard constructor
  for (Int_t i = 0; i < 6; i++) {
    for (Int_t j = 0; j < 2; j++) {
      fDetElemMisAlig[i][j] = 0.0;
      fModuleMisAlig[i][j] = 0.0;
      fSensorMisAlig[i][j] = 0.0;
      fHalfMisAlig[i][j] = 0.0;
    }
  }
  fDetElemMisAlig[0][1] = cartMisAlig;
  fDetElemMisAlig[1][1] = cartMisAlig;
  fDetElemMisAlig[5][1] = angMisAlig;
}

//_____________________________________________________________________________
GeometryMisAligner::GeometryMisAligner()
  : TObject(),
    fUseUni(kTRUE),
    fUseGaus(kFALSE),
    fXYAngMisAligFactor(0.0),
    fZCartMisAligFactor(0.0)
{
  /// Default constructor
  for (Int_t i = 0; i < 6; i++) {
    for (Int_t j = 0; j < 2; j++) {
      fDetElemMisAlig[i][j] = 0.0;
      fModuleMisAlig[i][j] = 0.0;
      fSensorMisAlig[i][j] = 0.0;
      fHalfMisAlig[i][j] = 0.0;
    }
  }
}

//______________________________________________________________________________
GeometryMisAligner::~GeometryMisAligner()
{
  /// Destructor
}

//_________________________________________________________________________
void GeometryMisAligner::SetXYAngMisAligFactor(Double_t factor)
{
  /// Set XY angular misalign factor

  if (TMath::Abs(factor) > 1.0 && factor > 0.) {
    fXYAngMisAligFactor = factor;
    fDetElemMisAlig[3][0] = fDetElemMisAlig[5][0] * factor; // These lines were
    fDetElemMisAlig[3][1] = fDetElemMisAlig[5][1] * factor; // added to keep
    fDetElemMisAlig[4][0] = fDetElemMisAlig[5][0] * factor; // backward
    fDetElemMisAlig[4][1] = fDetElemMisAlig[5][1] * factor; // compatibility
  } else
    LOG(ERROR) << "Invalid XY angular misalign factor, " << factor;
}

//_________________________________________________________________________
void GeometryMisAligner::SetZCartMisAligFactor(Double_t factor)
{
  /// Set XY angular misalign factor
  if (TMath::Abs(factor) < 1.0 && factor > 0.) {
    fZCartMisAligFactor = factor;
    fDetElemMisAlig[2][0] = fDetElemMisAlig[0][0];          // These lines were added to
    fDetElemMisAlig[2][1] = fDetElemMisAlig[0][1] * factor; // keep backward compatibility
  } else
    LOG(ERROR) << "Invalid Z cartesian misalign factor, " << factor;
}

//_________________________________________________________________________
void GeometryMisAligner::GetUniMisAlign(double cartMisAlig[3], double angMisAlig[3], const double lParMisAlig[6][2]) const
{
  /// Misalign using uniform distribution
  /**
    misalign the centre of the local transformation
    rotation axes :
    fAngMisAlig[1,2,3] = [x,y,z]
    Assume that misalignment about the x and y axes (misalignment of z plane)
    is much smaller, since the entire detection plane has to be moved (the
    detection elements are on a support structure), while rotation of the x-y
    plane is more free.
  */
  cartMisAlig[0] = gRandom->Uniform(-lParMisAlig[0][1] + lParMisAlig[0][0], lParMisAlig[0][0] + lParMisAlig[0][1]);
  cartMisAlig[1] = gRandom->Uniform(-lParMisAlig[1][1] + lParMisAlig[1][0], lParMisAlig[1][0] + lParMisAlig[1][1]);
  cartMisAlig[2] = gRandom->Uniform(-lParMisAlig[2][1] + lParMisAlig[2][0], lParMisAlig[2][0] + lParMisAlig[2][1]);

  angMisAlig[0] = gRandom->Uniform(-lParMisAlig[3][1] + lParMisAlig[3][0], lParMisAlig[3][0] + lParMisAlig[3][1]);
  angMisAlig[1] = gRandom->Uniform(-lParMisAlig[4][1] + lParMisAlig[4][0], lParMisAlig[4][0] + lParMisAlig[4][1]);
  angMisAlig[2] = gRandom->Uniform(-lParMisAlig[5][1] + lParMisAlig[5][0], lParMisAlig[5][0] + lParMisAlig[5][1]); // degrees
}

//_________________________________________________________________________
void GeometryMisAligner::GetGausMisAlign(double cartMisAlig[3], double angMisAlig[3], const double lParMisAlig[6][2]) const
{
  /// Misalign using gaussian distribution
  /**
    misalign the centre of the local transformation
    rotation axes :
    fAngMisAlig[1,2,3] = [x,y,z]
    Assume that misalignment about the x and y axes (misalignment of z plane)
    is much smaller, since the entire detection plane has to be moved (the
    detection elements are on a support structure), while rotation of the x-y
    plane is more free.
  */
  cartMisAlig[0] = gRandom->Gaus(lParMisAlig[0][0], lParMisAlig[0][1]); //, 3. * lParMisAlig[0][1]);
  cartMisAlig[1] = gRandom->Gaus(lParMisAlig[1][0], lParMisAlig[1][1]); //, 3. * lParMisAlig[1][1]);
  cartMisAlig[2] = gRandom->Gaus(lParMisAlig[2][0], lParMisAlig[2][1]); //, 3. * lParMisAlig[2][1]);

  angMisAlig[0] = gRandom->Gaus(lParMisAlig[3][0], lParMisAlig[3][1]); //, 3. * lParMisAlig[3][1]);
  angMisAlig[1] = gRandom->Gaus(lParMisAlig[4][0], lParMisAlig[4][1]); //, 3. * lParMisAlig[4][1]);
  angMisAlig[2] = gRandom->Gaus(lParMisAlig[5][0], lParMisAlig[5][1]); //, 3. * lParMisAlig[5][1]); // degrees
}

//_________________________________________________________________________
TGeoCombiTrans GeometryMisAligner::MisAlignSensor() const
{
  /// Misalign given transformation and return the misaligned transformation.
  /// Use misalignment parameters for sensor on the ladder.
  /// Note that applied misalignments are small deltas with respect to the detection
  /// element own ideal local reference frame. Thus deltaTransf represents
  /// the transformation to go from the misaligned d.e. local coordinates to the
  /// ideal d.e. local coordinates.
  /// Also note that this -is not- what is in the ALICE alignment framework known
  /// as local nor global (see GeometryMisAligner::MisAlign)

  Double_t cartMisAlig[3] = {0, 0, 0};
  Double_t angMisAlig[3] = {0, 0, 0};

  if (fUseUni) {
    GetUniMisAlign(cartMisAlig, angMisAlig, fSensorMisAlig);
  } else {
    if (!fUseGaus) {
      LOG(INFO) << "Neither uniform nor gausian distribution is set! Will use gausian...";
    }
    GetGausMisAlign(cartMisAlig, angMisAlig, fSensorMisAlig);
  }

  TGeoTranslation deltaTrans(cartMisAlig[0], cartMisAlig[1], cartMisAlig[2]);
  TGeoRotation deltaRot;
  deltaRot.RotateX(angMisAlig[0]);
  deltaRot.RotateY(angMisAlig[1]);
  deltaRot.RotateZ(angMisAlig[2]);

  TGeoCombiTrans deltaTransf(deltaTrans, deltaRot);
  //TGeoHMatrix newTransfMat = transform * deltaTransf;

  return TGeoCombiTrans(deltaTransf);
}

//_________________________________________________________________________
TGeoCombiTrans GeometryMisAligner::MisAlignDetElem() const
{
  /// Misalign given transformation and return the misaligned transformation.
  /// Use misalignment parameters for detection elements.
  /// Note that applied misalignments are small deltas with respect to the detection
  /// element own ideal local reference frame. Thus deltaTransf represents
  /// the transformation to go from the misaligned d.e. local coordinates to the
  /// ideal d.e. local coordinates.
  /// Also note that this -is not- what is in the ALICE alignment framework known
  /// as local nor global (see GeometryMisAligner::MisAlign)

  Double_t cartMisAlig[3] = {0, 0, 0};
  Double_t angMisAlig[3] = {0, 0, 0};

  if (fUseUni) {
    GetUniMisAlign(cartMisAlig, angMisAlig, fDetElemMisAlig);
  } else {
    if (!fUseGaus) {
      LOG(INFO) << "Neither uniform nor gausian distribution is set! Will use gausian...";
    }
    GetGausMisAlign(cartMisAlig, angMisAlig, fDetElemMisAlig);
  }

  TGeoTranslation deltaTrans(cartMisAlig[0], cartMisAlig[1], cartMisAlig[2]);
  TGeoRotation deltaRot;
  deltaRot.RotateX(angMisAlig[0]);
  deltaRot.RotateY(angMisAlig[1]);
  deltaRot.RotateZ(angMisAlig[2]);

  TGeoCombiTrans deltaTransf(deltaTrans, deltaRot);
  //TGeoHMatrix newTransfMat = transform * deltaTransf;

  //LOG(INFO) << "Translation Module by (" << cartMisAlig[0] << "," << cartMisAlig[1] << "," << cartMisAlig[2] << ") on (X,Y,Z) axis.";
  //LOG(INFO) << "Rotated Module by (" << angMisAlig[0] << "," << angMisAlig[1] << "," << angMisAlig[2] << ") about (X,Y,Z) axis.";

  return TGeoCombiTrans(deltaTransf);
}

//_________________________________________________________________________
TGeoCombiTrans
  GeometryMisAligner::MisAlignHalf() const
{
  /// Misalign given transformation and return the misaligned transformation.
  /// Use misalignment parameters for half-MFT.
  /// Note that applied misalignments are small deltas with respect to the module
  /// own ideal local reference frame. Thus deltaTransf represents
  /// the transformation to go from the misaligned module local coordinates to the
  /// ideal module local coordinates.
  /// Also note that this -is not- what is in the ALICE alignment framework known
  /// as local nor global (see GeometryMisAligner::MisAlign)

  Double_t cartMisAlig[3] = {0, 0, 0};
  Double_t angMisAlig[3] = {0, 0, 0};

  if (fUseUni) {
    GetUniMisAlign(cartMisAlig, angMisAlig, fHalfMisAlig);
  } else {
    if (!fUseGaus) {
      LOG(INFO) << "Neither uniform nor gausian distribution is set! Will use gausian...";
    }
    GetGausMisAlign(cartMisAlig, angMisAlig, fHalfMisAlig);
  }

  TGeoTranslation deltaTrans(cartMisAlig[0], cartMisAlig[1], cartMisAlig[2]);
  TGeoRotation deltaRot;
  deltaRot.RotateX(angMisAlig[0]);
  deltaRot.RotateY(angMisAlig[1]);
  deltaRot.RotateZ(angMisAlig[2]);

  TGeoCombiTrans deltaTransf(deltaTrans, deltaRot);
  //TGeoHMatrix newTransfMat = transform * deltaTransf;

  return TGeoCombiTrans(deltaTransf);
}

//_________________________________________________________________________
TGeoCombiTrans
  GeometryMisAligner::MisAlignModule() const
{
  /// Misalign given transformation and return the misaligned transformation.
  /// Use misalignment parameters for modules.
  /// Note that applied misalignments are small deltas with respect to the module
  /// own ideal local reference frame. Thus deltaTransf represents
  /// the transformation to go from the misaligned module local coordinates to the
  /// ideal module local coordinates.
  /// Also note that this -is not- what is in the ALICE alignment framework known
  /// as local nor global (see GeometryMisAligner::MisAlign)

  Double_t cartMisAlig[3] = {0, 0, 0};
  Double_t angMisAlig[3] = {0, 0, 0};

  if (fUseUni) {
    GetUniMisAlign(cartMisAlig, angMisAlig, fModuleMisAlig);
  } else {
    if (!fUseGaus) {
      LOG(INFO) << "Neither uniform nor gausian distribution is set! Will use gausian...";
    }
    GetGausMisAlign(cartMisAlig, angMisAlig, fModuleMisAlig);
  }

  TGeoTranslation deltaTrans(cartMisAlig[0], cartMisAlig[1], cartMisAlig[2]);
  TGeoRotation deltaRot;
  deltaRot.RotateX(angMisAlig[0]);
  deltaRot.RotateY(angMisAlig[1]);
  deltaRot.RotateZ(angMisAlig[2]);

  TGeoCombiTrans deltaTransf(deltaTrans, deltaRot);
  //TGeoHMatrix newTransfMat = transform * deltaTransf;

  //LOG(INFO) << "Translation Module by (" << cartMisAlig[0] << "," << cartMisAlig[1] << "," << cartMisAlig[2] << ") on (X,Y,Z) axis.";
  //LOG(INFO) << "Rotated Module by (" << angMisAlig[0] << "," << angMisAlig[1] << "," << angMisAlig[2] << ") about (X,Y,Z) axis.";

  return TGeoCombiTrans(deltaTransf);
}

//______________________________________________________________________
bool GeometryMisAligner::matrixToAngles(const double* rot, double& psi, double& theta, double& phi)
{
  /// Calculates the Euler angles in "x y z" notation
  /// using the rotation matrix
  /// Returns false in case the rotation angles can not be
  /// extracted from the matrix
  //
  if (std::abs(rot[0]) < 1e-7 || std::abs(rot[8]) < 1e-7) {
    LOG(ERROR) << "Failed to extract roll-pitch-yall angles!";
    return false;
  }
  psi = std::atan2(-rot[5], rot[8]);
  theta = std::asin(rot[2]);
  phi = std::atan2(-rot[1], rot[0]);
  return true;
}

//______________________________________________________________________
void GeometryMisAligner::MisAlign(Bool_t verbose)
{
  /// Takes the internal geometry module transformers, copies them to
  /// new geometry module transformers.
  /// Calculates  module misalignment parameters and applies these
  /// to the new module transformer.
  /// Calculates the module misalignment delta transformation in the
  /// Alice Alignment Framework newTransf = delta * oldTransf.
  /// Add a module misalignment to the new geometry transformer.
  /// Gets the Detection Elements from the module transformer.
  /// Calculates misalignment parameters and applies these
  /// to the local transformation of the Detection Element.
  /// Obtains the new global transformation by multiplying the new
  /// module transformer transformation with the new local transformation.
  /// Applies the new global transform to a new detection element.
  /// Adds the new detection element to a new module transformer.
  /// Calculates the d.e. misalignment delta transformation in the
  /// Alice Alignment Framework (newGlobalTransf = delta * oldGlobalTransf).
  /// Add a d.e. misalignment to the new geometry transformer.
  /// Adds the new module transformer to a new geometry transformer.
  /// Returns the new geometry transformer.

  mGeometryTGeo = GeometryTGeo::Instance();

  o2::detectors::AlignParam lAP;

  std::vector<std::vector<o2::detectors::AlignParam>> lAPvec;

  std::vector<o2::detectors::AlignParam> lAPvecModule;  // Storage of all AlignParam for each Module
  std::vector<o2::detectors::AlignParam> lAPvecDetElem; // Storage of all AlignParam for each DetElement

  // The function to get misalignment transformations of DE from surveyed positions
  o2::mft::ModuleTransform mTransform;

  Int_t nAlignID = 0;
  //static TGeoHMatrix matIGTransf;
  double lPsi, lTheta, lPhi = 0.;

  Int_t nHalf = mGeometryTGeo->getNumberOfHalfs();

  for (Int_t hf = 0; hf < nHalf; hf++) {

    Int_t nDisks = mGeometryTGeo->getNumberOfDisksPerHalf(hf);

    // New module transformation
    TGeoCombiTrans localDeltaTransform = MisAlignHalf();
    // localDeltaTransform.Print();

    TString sname = mGeometryTGeo->composeSymNameHalf(hf);

    lAP.setSymName(sname);
    lAP.setAlignableID(nAlignID++);

    lAP.setLocalParams(localDeltaTransform);

    if (!matrixToAngles(localDeltaTransform.GetRotationMatrix(), lPsi, lTheta, lPhi)) {
      LOG(ERROR) << "Problem extracting angles! from Half";
    }
    LOG(DEBUG) << "** LocalDeltaTransform Half: " << fmt::format("{} : {} | X: {:+f} Y: {:+f} Z: {:+f} | pitch: {:+f} roll: {:+f} yaw: {:+f}\n", lAP.getSymName(), lAP.getAlignableID(), localDeltaTransform.GetTranslation()[0], localDeltaTransform.GetTranslation()[1], localDeltaTransform.GetTranslation()[2], localDeltaTransform.GetRotationMatrix()[0], localDeltaTransform.GetRotationMatrix()[1], localDeltaTransform.GetRotationMatrix()[2]);

    lAP.setLocalParams(localDeltaTransform);
    // Apply misalignment of the half to the ideal geometry
    lAP.applyToGeometry();

    for (Int_t dk = 0; dk < nDisks; dk++) {

      // New module transformation
      //TGeoCombiTrans localDeltaTransform = MisAlignModule();
      localDeltaTransform = MisAlignModule();

      // localDeltaTransform.Print();
      //TString sname = mGeometryTGeo->composeSymNameDisk(hf, dk);
      sname = mGeometryTGeo->composeSymNameDisk(hf, dk);

      lAP.setSymName(sname);
      lAP.setAlignableID(nAlignID++);

      if (!matrixToAngles(localDeltaTransform.GetRotationMatrix(), lPsi, lTheta, lPhi)) {
        LOG(ERROR) << "Problem extracting angles!";
      }
      LOG(DEBUG) << "**** LocalDeltaTransform Disk: " << fmt::format("{} : {} | X: {:+f} Y: {:+f} Z: {:+f} | pitch: {:+f} roll: {:+f} yaw: {:+f}\n", lAP.getSymName(), lAP.getAlignableID(), localDeltaTransform.GetTranslation()[0], localDeltaTransform.GetTranslation()[1], localDeltaTransform.GetTranslation()[2], localDeltaTransform.GetRotationMatrix()[0], localDeltaTransform.GetRotationMatrix()[1], localDeltaTransform.GetRotationMatrix()[2]);

      // Set the local delta transformation to the module
      lAP.setLocalParams(localDeltaTransform);
      // lAP.Print();

      LOG(DEBUG) << "**** AlignParam Disk: " << fmt::format("{} : {} | X: {:+f} Y: {:+f} Z: {:+f} | pitch: {:+f} roll: {:+f} yaw: {:+f}\n", lAP.getSymName(), lAP.getAlignableID(), lAP.getX(), lAP.getY(), lAP.getZ(), lAP.getPsi(), lAP.getTheta(), lAP.getPhi());

      Int_t nLadders = 0;

      for (Int_t sensor = mGeometryTGeo->getMinSensorsPerLadder(); sensor < mGeometryTGeo->getMaxSensorsPerLadder() + 1; sensor++) {
        nLadders += mGeometryTGeo->getNumberOfLaddersPerDisk(hf, dk, sensor);
      }

      if (verbose) {
        LOG(INFO) << "----- Half-Disk:" << dk << " ----- Number of Ladders: " << nLadders << " -----";
      }

      // Apply misalignment to the ideal geometry
      lAP.applyToGeometry();

      // Store AlignParam (misalignment parameters)
      lAPvecModule.push_back(lAP);

      for (Int_t lr = 0; lr < nLadders; lr++) {

        localDeltaTransform = MisAlignDetElem();

        sname = mGeometryTGeo->composeSymNameLadder(hf, dk, lr);

        Int_t nSensorsPerLadder = mGeometryTGeo->getNumberOfSensorsPerLadder(hf, dk, lr);
        TString path = "/cave_1/barrel_1/" + sname;

        lAP.setSymName(sname);

        lAP.setAlignableID(nAlignID++);

        //        if (!matrixToAngles(localDeltaTransform.GetRotationMatrix(), lPsi, lTheta, lPhi)) {
        //          LOG(ERROR) << "Problem extracting angles!";
        //        }

        LOG(DEBUG) << "LocalDeltaTransform Ladder: " << fmt::format("{} : {} | X: {:+f} Y: {:+f} Z: {:+f} | pitch: {:+f} roll: {:+f} yaw: {:+f}\n", lAP.getSymName(), lAP.getAlignableID(), localDeltaTransform.GetTranslation()[0], localDeltaTransform.GetTranslation()[1], localDeltaTransform.GetTranslation()[2], lPsi, lTheta, lPhi);

        //        if (!lAP.setLocalParams(localDeltaTransform)) {
        //          LOG(ERROR) << "  Could not set local params for " << sname;
        //        }

        // Set the local transformations
        lAP.setLocalParams(localDeltaTransform);

        LOG(DEBUG) << "AlignParam Ladder: " << fmt::format("  {} : {} | X: {:+f} Y: {:+f} Z: {:+f} | pitch: {:+f} roll: {:+f} yaw: {:+f}\n", lAP.getSymName(), lAP.getAlignableID(), lAP.getX(), lAP.getY(), lAP.getZ(), lAP.getPsi(), lAP.getTheta(), lAP.getPhi());

        if (verbose) {
          LOG(INFO) << "---------> MisAligned element : " << sname << "   Ladder: " << lr;
        }

        // Apply misaligned detection element to the geometry
        lAP.applyToGeometry();

        // Store AlignParam (misalignment parameters)
        lAPvecDetElem.push_back(lAP);

        for (Int_t sr = 0; sr < nSensorsPerLadder; sr++) {

          localDeltaTransform = MisAlignSensor();

          sname = mGeometryTGeo->composeSymNameChip(hf, dk, lr, sr);

          lAP.setSymName(sname);
          lAP.setAlignableID(nAlignID++);

          lAP.setLocalParams(localDeltaTransform);
          lAP.applyToGeometry();

          Double_t parErr[6] = {0.002, 0.001, 0.003, 0, 0, 0};
          Int_t resultTransformFcn = mTransform.GetModuleMeanTransform(localDeltaTransform, parErr);
          LOG(INFO) << " ---->  result: ModuleTransform " << resultTransformFcn;
        }
      }
    }
  }

  lAPvec.push_back(lAPvecModule);
  lAPvec.push_back(lAPvecDetElem);

  //  for (o2::detectors::AlignParam MDalignparam : lAPvecModule) {
  //    MDalignparam.Print();
  //    for (o2::detectors::AlignParam DEalignparam : lAPvecDetElem) {
  //      DEalignparam.Print();
  //    }
  //  }

  // return newGeometryTransformer;
}

void GeometryMisAligner::SetAlignmentResolution(const TClonesArray* misAlignArray, Int_t rChId, Double_t rChResX, Double_t rChResY, Double_t rDeResX, Double_t rDeResY)
{

  Int_t chIdMin = (rChId < 0) ? 0 : rChId;
  Int_t chIdMax = (rChId < 0) ? 9 : rChId;
  Double_t chResX = (rChResX < 0) ? fModuleMisAlig[0][1] : rChResX;
  Double_t chResY = (rChResY < 0) ? fModuleMisAlig[1][1] : rChResY;
  Double_t deResX = (rDeResX < 0) ? fDetElemMisAlig[0][1] : rDeResX;
  Double_t deResY = (rDeResY < 0) ? fDetElemMisAlig[1][1] : rDeResY;

  TMatrixDSym mChCorrMatrix(6);
  mChCorrMatrix[0][0] = chResX * chResX;
  mChCorrMatrix[1][1] = chResY * chResY;
  //  mChCorrMatrix.Print();

  TMatrixDSym mDECorrMatrix(6);
  mDECorrMatrix[0][0] = deResX * deResX;
  mDECorrMatrix[1][1] = deResY * deResY;
  //  mDECorrMatrix.Print();

  /*
  AliAlignObjMatrix *alignMat = 0x0;

  for(Int_t chId=chIdMin; chId<=chIdMax; chId++) {
    TString chName1;
    TString chName2;
    if (chId<4){
      chName1 = Form("GM%d",chId);
      chName2 = Form("GM%d",chId);
    } else {
      chName1 = Form("GM%d",4+(chId-4)*2);
      chName2 = Form("GM%d",4+(chId-4)*2+1);
    }
    
    for (int i=0; i<misAlignArray->GetEntries(); i++) {
      alignMat = (AliAlignObjMatrix*)misAlignArray->At(i);
      TString volName(alignMat->GetSymName());
      if((volName.Contains(chName1)&&
      ((volName.Last('/')==volName.Index(chName1)+chName1.Length())||
       (volName.Length()==volName.Index(chName1)+chName1.Length())))||
     (volName.Contains(chName2)&&
      ((volName.Last('/')==volName.Index(chName2)+chName2.Length())||
       (volName.Length()==volName.Index(chName2)+chName2.Length())))){
    volName.Remove(0,volName.Last('/')+1);
    if (volName.Contains("GM")) {
      //    alignMat->Print("NULL");
      alignMat->SetCorrMatrix(mChCorrMatrix);
    } else if (volName.Contains("DE")) {
      //    alignMat->Print("NULL");
      alignMat->SetCorrMatrix(mDECorrMatrix);
    }
      }
    }
  }
    
    */
}

