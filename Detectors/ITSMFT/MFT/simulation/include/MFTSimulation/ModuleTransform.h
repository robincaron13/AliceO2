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

//-----------------------------------------------------------------------------
#ifndef MODULETRANSFORM_H
#define MODULETRANSFORM_H

#include <TObject.h>
#include <TMinuit.h>

class TGeoCombiTrans;
class TFitter;
class TMinuit;
class AlignParam;

namespace o2
{
namespace mft
{
class GeometryTGeo;
}
} // namespace o2

namespace o2
{
namespace mft
{
class ModuleTransform : public TObject
{
 public:
  ModuleTransform();
  ~ModuleTransform() override;

  /// Not implemented
  ModuleTransform(const ModuleTransform& right);
  /// Not implemented
  ModuleTransform& operator=(const ModuleTransform& right);

  /// Set align param
  void SetMisAlignPadPosition(const Int_t iPad, Double_t xpos, Double_t ypos, Double_t zpos)
  {
    fPadMisAlig[0][iPad] = xpos;
    fPadMisAlig[1][iPad] = ypos;
    fPadMisAlig[2][iPad] = zpos;
  }

  Double_t MyChi2Disk(Double_t* par);
  Double_t MyChi2Chip(Double_t* par);

  Int_t GetModuleMeanTransform(TGeoCombiTrans& modTransf, Double_t* parErr, Int_t level);

  /// Set sensor id name
  void SetSensorId(Int_t half, Int_t disk, Int_t ladder, Int_t sensor)
  {
    fHalfId = half;
    fDiskId = disk;
    fLadderId = ladder;
    fSensorId = sensor;
  }

  /// Set sensor id name
  void SetSymName(const char* m) { fsymName = m; }

  void PrintAlignResults();

 private:
  Int_t fModuleId;

  TGeoCombiTrans GetMatrixPadPosition() const;

  GeometryTGeo* mGeometryTGeo; //! access to geometry details

  Int_t fHalfId;
  Int_t fDiskId;
  Int_t fLadderId;
  Int_t fSensorId;
  const char* fsymName;

  Double_t fPadMisAlig[3][4]; ///< Displacements of the pad along x,y,z (translations)

  TFitter* fFitter; // fitter object
  TMinuit* gMinuit; // minuit object

  ClassDefOverride(ModuleTransform, 0)
};

} // namespace mft
} // namespace o2

#endif
