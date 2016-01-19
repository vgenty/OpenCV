//
// cint script to generate libraries
// Declaire namespace & classes you defined
// #pragma statement: order matters! Google it ;)
//

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ namespace imutil+;

#pragma link C++ class imutil::Coordinate_t+;
#pragma link C++ class imutil::Contour_t+;
#pragma link C++ class imutil::ContourArray_t+;
#pragma link C++ class imutil::ContourTreeMaker+;

//ADD_NEW_CLASS ... do not change this line
#endif

