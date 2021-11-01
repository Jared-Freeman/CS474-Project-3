#ifndef READIMAGE_H
#define READIMAGE_H


namespace std { //std namespace

#include <iostream>
#include <fstream>
#include <stdlib.h>

#include "image.h"

int readImage(char fname[], ImageType& image)
{
 int i, j;
 int N, M, Q;
 unsigned char *charImage;
 char header [100], *ptr;
 ifstream ifp;

 ifp.open(fname, ios::in | ios::binary);

 if (!ifp) {
   cout << "Can't read image: " << fname << endl;
   exit(1);
 }

 // read header
 //cout << "header\n";
 
 
 ifp.getline(header,100,'\n');
 if ( (header[0]!=80) ||    /* 'P' */
      (header[1]!=53) ) {   /* '5' */
      cout << "Image " << fname << " is not PGM" << endl;
      exit(1);
 }

ifp.getline(header,100,'\n');
 while(header[0]=='#')
   ifp.getline(header,100,'\n');

 M=strtol(header,&ptr,0);
 N=atoi(ptr);

 ifp.getline(header,100,'\n');
 Q=strtol(header,&ptr,0);

 charImage = (unsigned char *) new unsigned char [M*N];

 ifp.read( reinterpret_cast<char *>(charImage), (M*N)*sizeof(unsigned char));

 if (ifp.fail()) {
   cout << "Image " << fname << " has wrong size" << endl;
   exit(1);
 }

 ifp.close();

 //
 // Convert the unsigned characters to integers
 //
 
 //cout << "uc to i\n";

 int val;

 ImageType temp_img(N, M, Q);

 for(i=0; i<N; i++)
   for(j=0; j<M; j++) {
     val = (int)charImage[i*M+j];
     temp_img.setPixelVal(i, j, val);  
     
     //int b;
     //temp_img.getPixelVal(i, j, b);
     //cout << b << " ";   
   }
   //cout << endl;

 image.CopyImageData(temp_img);
 
 //cout << "deletion\n";
 
 
 delete [] charImage;

 return (1);

}

} //end namespace
#endif
