 #include <iostream>
 // #include <string>
 // #include <sstream>
 #include <fstream>
 #include <cstring>
 #include <base64.h>
 // #include <climits>
 // #include <stdio.h>

 using namespace std;

extern "C" {
  // remember that in fortran values are passed by reference, so we need to receive them
  // using pointers in C/C++

  // Examples //
  // void print_int_(int * num1, int * num2) {
  //   std::cout << "num1:" << *num1 << " and num2:" << *num2 << std::endl;
  // }
  // void print_double_(double * num1, double * num2) {
  //   std::cout << "num1:" << *num1 << " and num2:" << *num2 << std::endl;
  // }


  void print_matrix_(double * meshZ, double * meshR, double * vector_fromF, int * dim1, int * dim2){
    double value;
    int* flen = new(int);
    char* s;
    float val = 0.0;

    //open file
    ofstream outfile;
    outfile.open("rho_tot.vtr");

    //header
    outfile << "<?xml version=\"1.0\"?>" << endl;
    outfile << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=";
    // if(isLittleEndian()){
      outfile << "\"LittleEndian\" header_type=\"UInt64\">\n";
    // } else {
    //   outfile << "\"BigEndian\">\n";
    // }

    //rectilinear grid
    int x_left=0; int y_left=0; int z_left=0;
    int x_right= *dim2;
    int y_right= 1;
    int z_right= *dim1;

    outfile << "<RectilinearGrid WholeExtent=\"";
    outfile << x_left << " ";
    outfile << x_right << " ";
    outfile << y_left << " ";
    outfile << y_right << " ";
    outfile << z_left << " ";
    outfile << z_right ;
    outfile << "\">" << endl;

    outfile << "<Piece Extent=\"";
    outfile << x_left << " ";
    outfile << x_right << " ";
    outfile << y_left << " ";
    outfile << y_right << " ";
    outfile << z_left << " ";
    outfile << z_right ;
    outfile << "\">" << endl;


    // This writes the point coordinates on the three axis
    outfile << "<Coordinates>" << endl;
    outfile << "<DataArray type=\"Float64\" Name=Z_axis";
    outfile << " format=\"binary\">\n";

    // long long byte_number = 8*((long long)(*dim1)+1);
    // char bins[8+8*((*dim1)+1)];
    // memcpy(bins, (char*)&byte_number, 8);
    // for(int i=0; i<*dim1; i++){
    //   value = (double)meshZ[i];
    //   memcpy(bins+(8+8*i), (char*)&value, 8);
    // }
    // value = value = (double)meshZ[*dim1];
    // memcpy(bins+(8+8*(*dim1)), (char*)&value, 8);
    // s=base64(bins, 8+8*(*dim1+1), flen);
    // outfile.write(s, *flen);
    // outfile << endl;
    // outfile << "</DataArray>" << endl;


    for(int i=0;i<*dim1;i++){
      for(int j=0;j<*dim2;j++){
        cout << vector_fromF[i*(*dim2)+j] << "\t"; // << *num2 << std::endl;
    }
    cout << "\n";
  }
  }

}
