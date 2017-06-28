 #include <iostream>
 // #include <string>
 // #include <sstream>
 #include <fstream>
 #include <cstring>
 #include "base64.h"
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
    int dimension1= *dim1;
    int dimension2= *dim2;

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
    int x_right= dimension1-1;
    int y_right= 1;
    int z_right= dimension2-1;

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
    outfile << "<DataArray type=\"Float64\" Name=\"Z_axis\"";
    outfile << " format=\"binary\">\n";

    // long long byte_number = 8*((long long)(dimension1)+1);
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
    long long byte_number = 8*((long long)dimension1);
    char bins[8+8*(dimension1)];
    memcpy(bins, (char*)&byte_number, 8);
    for(int i=0; i<dimension1; i++){
      value = meshZ[i];
      memcpy(bins+(8+8*i), (char*)&value, 8);
    }
    s=base64(bins, 8+8*(dimension1), flen);
    outfile.write(s, *flen);
    outfile << "</DataArray>" << endl;


    //---third dimension---//
    outfile << "<DataArray type=\"Float64\" Name=\"T_axis\"";
    outfile << " format=\"binary\">\n";
    byte_number = 8*((long long)2);
    char bins3[8+8*(2)];
    memcpy(bins3, (char*)&byte_number, 8);
    value = 0.0;
    memcpy(bins3+(8), (char*)&value, 8);
    value = 0.0;
    memcpy(bins3+(8+8*1), (char*)&value, 8);
    s=base64(bins3, 8+8*2, flen);
    outfile.write(s, *flen);
    outfile << endl;
    outfile << "</DataArray>" << endl;


    //---second dimension---//
    outfile << "<DataArray type=\"Float64\" Name=\"R_axis\"";
    outfile << " format=\"binary\">\n";
    byte_number = 8*((long long)dimension2);
    char bins2[8+8*(dimension2)];
    memcpy(bins2, (char*)&byte_number, 8);
    for(int i=0; i<dimension2; i++){
      value = (double)meshR[i];
      memcpy(bins2+(8+8*i), (char*)&value, 8);
    }
    s=base64(bins2, 8+8*(dimension2), flen);
    outfile.write(s, *flen);
    outfile << endl;
    outfile << "</DataArray>" << endl;

    outfile << "</Coordinates>" << endl;

    // outfile << "<CellData Scalars=\"" << "rho_tot" << "\">\n";
    // write_variable(outfile, data);


    //---scalar quantity---//
    outfile << "<CellData Scalars=\"" << "rho_tot" << "\">\n";
    byte_number = 8*((long long)(dimension1*2*dimension2));
    char binsS[8+8*(dimension1*2*dimension2)];
    memcpy(binsS, (char*)&byte_number, 8);
      outfile << "<DataArray type=\"Float64\" Name=";
      outfile << "\"" << "rho_tot" << "\"";
      outfile << " format=\"binary\">\n";
      int byte_count=0;
      for(int k=0; k<dimension2-1; k++){
        for(int j=0; j<2; j++){
          for(int i=0; i<dimension1-1; i++){
            value = (double)vector_fromF[i*dimension2+k];
            memcpy(binsS+(8+8*byte_count), (char*)&value, 8);
            byte_count++;
          }
        }
      }
      s=base64(binsS, 8+8*dimension1*2*dimension2, flen);
      outfile.write(s, *flen);
      outfile << endl;
      outfile << "</DataArray>" << endl;
      outfile << "</CellData>" << endl;

      outfile << "</Piece>" << endl;
      outfile << "</RectilinearGrid>" << endl;


    outfile << "</VTKFile>" << endl;




  //   for(int i=0;i<*dim1;i++){
  //     for(int j=0;j<*dim2;j++){
  //       cout << vector_fromF[i*(*dim2)+j] << "\t"; // << *num2 << std::endl;
  //   }
  //   cout << "\n";
  // }
   }

}
