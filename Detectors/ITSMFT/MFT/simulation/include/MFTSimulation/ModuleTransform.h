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

class TGeoCombiTrans;
class TFitter;

//class AliMUONGeometryTransformer;
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

  //void SetModuleId(Int_t mId) { fModuleId = mId; }
  //	void SetIdealGeoTransformer(AliMUONGeometryTransformer *iGeoTransformer) {fIdealGeoTransformer = iGeoTransformer;}
  //	void SetAlignGeoTransformer(AliMUONGeometryTransformer *aAlignTransformer) {fAlignGeoTransformer = aAlignTransformer;}
  //void SetIdealGeoTransformer();
  //void SetAlignGeoTransformer();

  Double_t MyChi2(Double_t* par);

  //void MyFcn(int& npar, double* g, double& f, double* par, int iflag);

  Int_t GetModuleMeanTransform(TGeoCombiTrans& modTransf, Double_t* parErr);

 private:
  Int_t fModuleId;

  //    AliMUONGeometryTransformer *fIdealGeoTransformer;
  //    AliMUONGeometryTransformer *fAlignGeoTransformer;

  //    AliMpExMap *detElements;

  TFitter* fFitter;
    
    
  ClassDefOverride(ModuleTransform, 0)
};

} // namespace mft
} // namespace o2

#endif
