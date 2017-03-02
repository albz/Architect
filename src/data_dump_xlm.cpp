

// #include <cstring>
// #include <boost/filesystem.hpp>
//
// extern "C" {
//
// #if defined (__xlC__)
//   void create_folder(char * folderName, size_t len) {
// #elif defined (_MSC_VER)
//   void CREATE_FOLDER(char * folderName, size_t len) {
// #else
//   void create_folder_(char * folderName, size_t len) {
// #endif
//     std::string fname(folderName, 0, len);
//     if (!boost::filesystem::exists(fname)) {
//       boost::filesystem::create_directories(fname);
//     }
    /*
    else {
      std::string oldFolder = fname+".old";
      boost::filesystem::rename(fname, oldFolder);
      boost::filesystem::create_directories(fname);
    }
    */
//   }
//
//
// }
