// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef _HMPID_DIGITS_TO_RAW_SPEC_H_
#define _HMPID_DIGITS_TO_RAW_SPEC_H_

#include <vector>
#include "Framework/DataProcessorSpec.h"
#include "Framework/Task.h"
#include "DetectorsRaw/HBFUtils.h"
#include "DetectorsRaw/RawFileWriter.h"

#include "HMPIDBase/Common.h"
#include "HMPIDBase/Digit.h"
#include "HMPIDBase/Trigger.h"
#include "HMPIDSimulation/HmpidCoder2.h"

namespace o2
{
namespace hmpid
{

class DigitsToRawSpec : public framework::Task
{
 public:
  DigitsToRawSpec() = default;
  ~DigitsToRawSpec() override = default;
  void init(framework::InitContext& ic) final;
  void run(framework::ProcessingContext& pc) final;
  void endOfStream(framework::EndOfStreamContext& ec) override;

 private:
  void readRootFile();
  std::string mBaseFileName = "";
  std::string mDirectoryName = "";
  std::string mBaseRootFileName = "";
  bool mSkipEmpty = false;
  bool mDumpDigits = false;
  std::string mFileFor = "all";

  std::vector<o2::hmpid::Digit> mDigits;
  long mDigitsReceived;
  int mEventsReceived;
  HmpidCoder2* mCod;
  ExecutionTimer mExTimer;
  TTree* mDigTree;
};

o2::framework::DataProcessorSpec getDigitsToRawSpec();

} // end namespace hmpid
} // end namespace o2

#endif
