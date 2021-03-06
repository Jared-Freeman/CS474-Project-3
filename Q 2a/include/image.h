#ifndef IMAGE_H
#define IMAGE_H

// a simple example - you would need to add more funtions

#include <iostream>

class ImageType {
 public:
   ImageType();
   ImageType(int, int, int);
   ImageType(ImageType&);
   ~ImageType();
   void getImageInfo(int&, int&, int&);
   void setImageInfo(int, int, int);
   void setPixelVal(int, int, double);
   void getPixelVal(int, int, double&);

   void RemapPixelValues();
      
   void ClearImageData()
   {
			if(pixelValue != NULL)
			{
				//std::cout << pixelValue << std::endl;
			
				int i;

				for(i=0; i<N; i++)
				{
					delete [] pixelValue[i];
				}
				if(M > 0) 
					delete [] pixelValue;
			}
   }
   
   void CopyImageData(ImageType& oldImage)
   {
 		ClearImageData();  
   		 int i, j;

			 N = oldImage.N;
			 M = oldImage.M;
			 Q = oldImage.Q;

			 pixelValue = new double* [N];
			 for(i=0; i<N; i++) {
				 pixelValue[i] = new double[M];
				 for(j=0; j<M; j++)
					 pixelValue[i][j] = oldImage.pixelValue[i][j];
			 }/**/
   		
   }
 private:
 
   int N, M, Q;
   double **pixelValue;
};


#endif
