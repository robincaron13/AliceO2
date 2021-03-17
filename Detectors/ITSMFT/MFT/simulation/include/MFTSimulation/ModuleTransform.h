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

  //    /// Set geometry instance
  //    void SetGeometry(GeometryTGeo* geo)
  //    {
  //        mGeometryTGeo = geo;
  //    }
  //

  /// Set align param
  void SetMisAlignPadPosition(const Int_t iPad, Double_t xpos, Double_t ypos, Double_t zpos)
  {
    fPadMisAlig[0][iPad] = xpos;
    fPadMisAlig[1][iPad] = ypos;
    fPadMisAlig[2][iPad] = zpos;
  }

  Double_t MyChi2Disk(Double_t* par);
  Double_t MyChi2Chip(Double_t* par);

  //function to insert in the tfitter for minuit minimization
  //void MyFcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

  Int_t GetModuleMeanTransform(TGeoCombiTrans& modTransf, Double_t* parErr, Int_t level);

  //Int_t SetMatrixTransform(TGeoCombiTrans& mDETransf, Double_t* parTrf);

  /// Set sensor id name
  void SetSensorId(Int_t half, Int_t disk, Int_t ladder, Int_t sensor)
  {
    fHalfId = half;
    fDiskId = disk;
    fLadderId = ladder;
    fSensorId = sensor;
  }
  //void InitMyChi2(Int_t &npad, Double_t &matrix[3][4]);

  /// Set sensor id name
  void SetSymName(const char* m) { fsymName = m; }

  void PrintAlignResults();

  //    void GetPadPositions(const Int_t iPad, Double_t matrixPadMisAlig[3][4])
  //    {
  //        matrixPadMisAlig[0][iPad] = fPadMisAlig[0][iPad];
  //        matrixPadMisAlig[1][iPad] = fPadMisAlig[1][iPad];
  //        matrixPadMisAlig[2][iPad] = fPadMisAlig[2][iPad];
  //    }

 private:
  Int_t fModuleId;

  //TGeoCombiTrans fiDETrf;
  TGeoCombiTrans GetMatrixPadPosition() const;

  GeometryTGeo* mGeometryTGeo; //! access to geometry details
                               //    o2::detectors::AlignParam* fAliPar;

  Int_t fHalfId;
  Int_t fDiskId;
  Int_t fLadderId;
  Int_t fSensorId;
  const char* fsymName;

  Double_t fPadMisAlig[3][4]; ///< Mean and width of the displacements of the modules along x,y,z (translations) and about x,y,z (rotations)

  //    AliMUONGeometryTransformer *fIdealGeoTransformer;
  //    AliMUONGeometryTransformer *fAlignGeoTransformer;
  //    AliMpExMap *detElements;

  TFitter* fFitter; // fitter object
  TMinuit* gMinuit; // minuit object

  ClassDefOverride(ModuleTransform, 0)
};

} // namespace mft
} // namespace o2

#endif
