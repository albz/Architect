// #include <cstring>
// #include <iostream>
// #include <stdexcept>
//
//
// extern "C" {
//
//   #if defined (__xlC__)
//     void print_vtk(double num1, double num2) {
//   #elif defined (_MSC_VER)
//     void PRINT_VTK(double num1, double num2) {
//   #else
//     void print_vtk_(double num1, double num2) {
//   #endif
//     // std::cout << num1 << std::endl;
//     // std::cout << num2 << std::endl;
//     // std::cout << "end" << std::endl;
//     }
//   }
// //
// // #if defined (__xlC__)
// //   void create_folder(char * folderName, size_t len) {
// // #elif defined (_MSC_VER)
// //   void CREATE_FOLDER(char * folderName, size_t len) {
// // #else
// //   void create_folder_(char * folderName, size_t len) {
// // #endif
// //     std::string fname(folderName, 0, len);
// //     if (!boost::filesystem::exists(fname)) {
// //       boost::filesystem::create_directories(fname);
// //     }
//     /*
//     else {
//       std::string oldFolder = fname+".old";
//       boost::filesystem::rename(fname, oldFolder);
//       boost::filesystem::create_directories(fname);
//     }
//     */
// //   }
// //
// //
// // }
