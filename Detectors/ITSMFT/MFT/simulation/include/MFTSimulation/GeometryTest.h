// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_MFT_GEOMETRYTEST_H
#define ALICEO2_MFT_GEOMETRYTEST_H

#include <iostream>

class TH2;

namespace o2
{
namespace mft
{
namespace test
{

/// creates MFT regular geometry
void createRegularGeometry();

/// add alignable volumes
/// useull for tests.
void addAlignableVolumes();

/// misalign the MFT regular geometry
void MisalignGeometry();

class Dummy
{
  // to force Root produce a dictionary for namespace test (seems it is doing it fully if there are only functions in the namespace)
};

} // namespace test
} // namespace mft
} // namespace o2

#endif //ALICEO2_MFT_GEOMETRYTEST_H
