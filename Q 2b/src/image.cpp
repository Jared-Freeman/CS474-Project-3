#include <stdlib.h>

#include "image.h"

ImageType::ImageType()
{
 N = 0;
 M = 0;
 Q = 0;

 pixelValue = NULL;
}

ImageType::ImageType(int tmpN, int tmpM, int tmpQ)
{
 int i, j;

 N = tmpN;
 M = tmpM;
 Q = tmpQ;

 pixelValue = new double* [N];
 for(i=0; i<N; i++) {
   pixelValue[i] = new double[M];
   for(j=0; j<M; j++)
     pixelValue[i][j] = 0;
 }
}

ImageType::ImageType(ImageType& oldImage)
{
 int i, j;

 N = oldImage.N;
 M = oldImage.M;
 Q = oldImage.Q;

 pixelValue = new double* [N];
 for(i=0; i<N; i++) {
   pixelValue[i] = new double[M];
   for(j=0; j<M; j++)
     pixelValue[i][j] = oldImage.pixelValue[i][j];
 }
}

ImageType::~ImageType()
{
 int i;

 for(i=0; i<N; i++)
   delete [] pixelValue[i];
 delete [] pixelValue;
}


void ImageType::getImageInfo(int& rows, int& cols, int& levels)
{
 rows = N;
 cols = M;
 levels = Q;
} 

void ImageType::setImageInfo(int rows, int cols, int levels)
{
 N= rows;
 M= cols;
 Q= levels;
} 

void ImageType::setPixelVal(int i, int j, double val)
{
 pixelValue[i][j] = val;
}

void ImageType::getPixelVal(int i, int j, double& val)
{
 val = pixelValue[i][j];
}

void ImageType::RemapPixelValues()
{
    //get min and max from pxvals
    int max = -2147483647;
    int min = 2147483647;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            double val;
            getPixelVal(i, j, val);
            if (val > max) max = val;
            if (val < min) min = val;
        }
    }

    double slope = ((double)Q) / ((double)max - (double)min);

    // std::cout << "SLOPE: " << slope << std::endl;

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            double val;
            getPixelVal(i, j, val);
            int output = slope * (val - min);
            setPixelVal(i, j, output);
        }
    }
}
