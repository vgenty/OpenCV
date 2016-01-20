//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class larlite::HitImageManager+;

#pragma link C++ class larlite::BaseAlgoImage+;
#pragma link C++ class larlite::AlgoImageMaker+;

#pragma link C++ class larlite::BaseAlgoCluster+;
#pragma link C++ class larlite::AlgoClusterHull+;
#pragma link C++ class larlite::AlgoClusterHoughConnect+;
#pragma link C++ class larlite::AlgoClusterHoughSimilar+;

#pragma link C++ class ProtoCluster+;
//ADD_NEW_CLASS ... do not change this line
#endif






