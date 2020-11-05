// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifdef __CLING__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace o2;
#pragma link C++ namespace o2::mft;
//#pragma link C++ namespace o2::mft::test;

#pragma link C++ class o2::mft::Detector + ;
#pragma link C++ class o2::base::DetImpl < o2::mft::Detector> + ;
#pragma link C++ class o2::mft::DigitizerTask + ;
#pragma link C++ class o2::mft::GeometryMisAligner + ;
//#pragma link C++ class o2::mft::test::Dummy;
#pragma link C++ class o2::mft::ModuleTransform + ;

#pragma link C++ function o2::mft::test::addAlignableVolumes;
#pragma link C++ function o2::mft::test::createRegularGeometry;
#pragma link C++ function o2::mft::test::MisalignGeometry;

#endif
