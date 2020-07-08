// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file GeometryBuilder.h
/// \brief Class describing MFT Geometry Builder
/// \author Raphael Tieulent <raphael.tieulent@cern.ch>
/// \date 09/06/2015

#ifndef ALICEO2_MFT_GEOMETRYBUILDER_H_
#define ALICEO2_MFT_GEOMETRYBUILDER_H_

class TGeoHMatrix;

namespace o2
{
namespace mft
{

class GeometryBuilder
{

 public:
  GeometryBuilder() = default;
  ~GeometryBuilder() = default;

  // static methods
  static TGeoHMatrix Multiply(const TGeoHMatrix& m1, const TGeoHMatrix& m2);
  static TGeoHMatrix Multiply(const TGeoHMatrix& m1, const TGeoHMatrix& m2,
                              const TGeoHMatrix& m3);
  static TGeoHMatrix Multiply(const TGeoHMatrix& m1, const TGeoHMatrix& m2,
                              const TGeoHMatrix& m3, const TGeoHMatrix& m4);

  void buildGeometry();

 private:
  ClassDef(GeometryBuilder, 1);
};
} // namespace mft
} // namespace o2

#endif
