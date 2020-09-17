#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "MFTSimulation/GeometryTest.h"

#include "TGeoManager.h"
#include <iostream>


void MisalignMFTGeometry(const std::string inputGeom = ""){
    
    o2::base::GeometryManager::loadGeometry(inputGeom);

    auto gm = o2::mft::GeometryTGeo::Instance(); // geometry manager for mapping


    LOG(INFO) << "MisalignMFTGeometry ";

    o2::mft::test::MisalignGeometry();
    

}

#endif

