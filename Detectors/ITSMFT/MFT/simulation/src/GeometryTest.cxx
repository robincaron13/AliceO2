// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include "MFTSimulation/GeometryTest.h"

#include "DetectorsBase/GeometryManager.h"
#include "DetectorsBase/MaterialManager.h"
#include "Math/GenVector/Cartesian3D.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TH2F.h"
#include <iostream>

#include "MFTSimulation/GeometryMisAligner.h"

#include "MFTBase/Geometry.h"
#include "MFTBase/GeometryTGeo.h"
#include "MFTBase/GeometryBuilder.h"

#include "TPRegexp.h"
#include "TGLViewer.h"
#include "TGLRnrCtx.h"
#include "TVirtualPad.h"

#include "MFTSimulation/Detector.h"

#include "DetectorsCommonDataFormats/DetID.h"
#include "DetectorsCommonDataFormats/AlignParam.h"

namespace o2
{
namespace mft
{
namespace test
{

void createRegularGeometry()
{
  if (gGeoManager && gGeoManager->GetTopVolume()) {
    std::cerr << "Can only call this function with an empty geometry, i.e. gGeoManager==nullptr "
              << " or gGeoManager->GetTopVolume()==nullptr\n";
  }
  TGeoManager* g = new TGeoManager("MFT-BASICS", "ALICE MFT Regular Geometry");
  //o2::passive::Cave("CAVE", "Cave (for MFT Basics)").ConstructGeometry();
  // o2::passive::Dipole("DIPO", "Alice Dipole (for MCH Basics)").ConstructGeometry();
  //  o2::passive::Compensator("COMP", "Alice Compensator Dipole (for MCH Basics)").ConstructGeometry();
  //o2::passive::Pipe("PIPE", "Beam pipe (for MFT Basics)").ConstructGeometry();
  //o2::passive::Shil("SHIL", "Small angle beam shield (for MFT Basics)").ConstructGeometry();
  //  o2::passive::Absorber("ABSO", "Absorber (for MCH Basics)").ConstructGeometry();
  o2::mft::Detector(true).ConstructGeometry();
}

void addAlignableVolumes()
{
  if (!gGeoManager) {
    std::cerr << "gGeoManager == nullptr, must create a geometry first\n";
    return;
  }
  // If not closed, we need to close it
  if (!gGeoManager->IsClosed()) {
    gGeoManager->CloseGeometry();
  }
  // Then add the alignable volumes
  o2::mft::Detector(true).addAlignableVolumes();

  LOG(INFO) << "MFT  addAlignableVolumes ";
}

void MisalignGeometry()
{
  // create a regular geometry
  createRegularGeometry();

  LOG(INFO) << "Created MFT Regular Geometry ";

  if (!gGeoManager) {
    std::cerr << "gGeoManager == nullptr, must create a geometry first\n";
    return;
  }
  // If not closed, we need to close it
  if (!gGeoManager->IsClosed()) {
    gGeoManager->CloseGeometry();
  }

  // Then add the alignable volumes
  o2::mft::Detector(true).addAlignableVolumes();

  // The misaligner
  o2::mft::GeometryMisAligner aGMA;

  aGMA.SetModuleCartMisAlig(0.1, 0., 0.2, 0., 0.3, 0.);
  //aGMA.SetModuleAngMisAlig(0.1, 0., 0.2, 0., 0.3, 0.);

  aGMA.SetCartMisAlig(0.1, 0., 0.1, 0., 0.1, 0.);
  aGMA.SetAngMisAlig(0.1, 0.0);

  aGMA.MisAlign(true);

  LOG(INFO) << "Misaligned MFT Geometry ";
}

} // namespace test
} // namespace mft
} // namespace o2
