// !*****************************************************************************************************!
// !             Copyright 2014-2017 Alberto Marocchino, Francesco Massimo                               !
// !*****************************************************************************************************!
//
// !*****************************************************************************************************!
// !  This file is part of architect.                                                                    !
// !                                                                                                     !
// !  Architect is free software: you can redistribute it and/or modify                                  !
// !  it under the terms of the GNU General Public License as published by                               !
// !  the Free Software Foundation, either version 3 of the License, or                                  !
// !  (at your option) any later version.                                                                !
// !                                                                                                     !
// !  Architect is distributed in the hope that it will be useful,                                       !
// !  but WITHOUT ANY WARRANTY; without even the implied warranty of                                     !
// !  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                      !
// !  GNU General Public License for more details.                                                       !
// !                                                                                                     !
// !  You should have received a copy of the GNU General Public License                                  !
// !  along with architect.  If not, see <http://www.gnu.org/licenses/>.                                 !
// !*****************************************************************************************************!
//
// by A. Marocchino with the valuable support of M. Cipriani, S. Sinigardi
//
//

 #include <iostream>
 #include <string>
 #include <fstream>
 #include <cstring>
 #include "base64.h"

 using namespace std;

extern "C" {
  // remember that in fortran values are passed by reference, so we need to receive them using pointers in C/C++

  void print_matrix_(int * idvariable, double * meshZ, double * meshR, double * vector_fromF, int * dim1, int * dim2, char * filename, int * filename_length ){
    double value;
    int* flen = new(int);
    char* s,file_name;
    float val = 0.0;
    int dimension1= *dim1;
    int dimension2= *dim2;

    //--- open file ---//
    ofstream outfile(filename);

    //--- header ---//
    outfile << "<?xml version=\"1.0\"?>" << endl;
    outfile << "<VTKFile type=\"RectilinearGrid\" version=\"0.1\" byte_order=";
    // if(isLittleEndian()){
    outfile << "\"LittleEndian\" header_type=\"UInt64\">\n";
    // } else {
    //   outfile << "\"BigEndian\">\n";
    // }

    //rectilinear grid
    int x_left=0; int y_left=0; int z_left=0;
    int z_right= dimension1-1;
    int y_right= dimension2-1;
    int x_right= 1;

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

    //--- second 2 points dimension ---//
    outfile << "<DataArray type=\"Float64\" Name=\"T_axis\"";
    outfile << " format=\"binary\">\n";
    long long byte_number = 8*((long long)2);
    char bins3[8+8*(2)];
    memcpy(bins3, (char*)&byte_number, 8);
    value = 0.0;
    memcpy(bins3+(8), (char*)&value, 8);
    value = 0.1;
    memcpy(bins3+(8+8*1), (char*)&value, 8);
    s=base64(bins3, 8+8*2, flen);
    outfile.write(s, *flen);
    outfile << endl;
    outfile << "</DataArray>" << endl;


    //--- R axis ---//
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

    //--- Z axis ---//
    outfile << "<DataArray type=\"Float64\" Name=\"Z_axis\"";
    outfile << " format=\"binary\">\n";
    byte_number = 8*((long long)dimension1);
    char bins[8+8*(dimension1)];
    memcpy(bins, (char*)&byte_number, 8);
    for(int i=0; i<dimension1; i++){
      value = meshZ[i];
      memcpy(bins+(8+8*i), (char*)&value, 8);
    }
    s=base64(bins, 8+8*(dimension1), flen);
    outfile.write(s, *flen);
    outfile << "</DataArray>" << endl;


    outfile << "</Coordinates>" << endl;

    //---scalar quantity---//
    if(*idvariable==1){outfile << "<CellData Scalars=\"" << "rho_bunch" << "\">\n";};
    if(*idvariable==2){outfile << "<CellData Scalars=\"" << "rho_bck" << "\">\n";};
    if(*idvariable==3){outfile << "<CellData Scalars=\"" << "rho_tot" << "\">\n";};
    if(*idvariable==4){outfile << "<CellData Scalars=\"" << "Er_bck" << "\">\n";};
    if(*idvariable==5){outfile << "<CellData Scalars=\"" << "Er_bunch" << "\">\n";};
    if(*idvariable==6){outfile << "<CellData Scalars=\"" << "Er_tot" << "\">\n";};
    if(*idvariable==7){outfile << "<CellData Scalars=\"" << "Ez_bck" << "\">\n";};
    if(*idvariable==8){outfile << "<CellData Scalars=\"" << "Ez_bunch" << "\">\n";};
    if(*idvariable==9){outfile << "<CellData Scalars=\"" << "Ez_tot" << "\">\n";};
    byte_number = 8*((long long)(dimension1*2*dimension2));
    char binsS[8+8*(dimension1*2*dimension2)];
    memcpy(binsS, (char*)&byte_number, 8);
    outfile << "<DataArray type=\"Float64\" Name=";
    if(*idvariable==1){outfile << "\"" << "rho_bunch" << "\"";};
    if(*idvariable==2){outfile << "\"" << "rho_bck" << "\"";};
    if(*idvariable==3){outfile << "\"" << "rho_tot" << "\"";};
    if(*idvariable==4){outfile << "\"" << "Er_bck" << "\"";};
    if(*idvariable==5){outfile << "\"" << "Er_bunch" << "\"";};
    if(*idvariable==6){outfile << "\"" << "Er_tot" << "\"";};
    if(*idvariable==7){outfile << "\"" << "Ez_bck" << "\"";};
    if(*idvariable==8){outfile << "\"" << "Ez_bunch" << "\"";};
    if(*idvariable==9){outfile << "\"" << "Ez_tot" << "\"";};
    outfile << " format=\"binary\">\n";
    int byte_count=0;
    for(int j=0; j<2; j++){
      for(int i=0; i<dimension1-1; i++){
        for(int k=0; k<dimension2-1; k++){
          value = (double)vector_fromF[k*(dimension1-1)+i];
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
  }

}
