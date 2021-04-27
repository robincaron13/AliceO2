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

#pragma link C++ class o2::zdc::Hit + ;
#pragma link C++ class std::vector < o2::zdc::Hit> + ;

#pragma link C++ class o2::zdc::MCLabel + ;
#include "SimulationDataFormat/MCTruthContainer.h"
#pragma link C++ class o2::dataformats::MCTruthContainer < o2::zdc::MCLabel> + ;

#pragma link C++ class o2::zdc::ChannelData + ;
#pragma link C++ class o2::zdc::BCData + ;
#pragma link C++ class o2::zdc::OrbitData + ;
#pragma link C++ class std::vector < o2::zdc::ChannelData> + ;
#pragma link C++ class std::vector < o2::zdc::BCData> + ;
#pragma link C++ class std::vector < o2::zdc::OrbitData> + ;

#pragma link C++ class std::vector < o2::zdc::RecEvent> + ;
#pragma link C++ class std::vector < o2::zdc::OrbitRawData> + ;
#pragma link C++ class std::vector < o2::zdc::OrbitRecData> + ;

#pragma link C++ struct o2::zdc::CTFHeader + ;
#pragma link C++ struct o2::zdc::CTF + ;
#pragma link C++ class o2::ctf::EncodedBlocks < o2::zdc::CTFHeader, 12, uint32_t> + ;

#endif
