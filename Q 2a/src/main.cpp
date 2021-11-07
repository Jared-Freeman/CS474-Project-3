#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <cstring>
#include <map>
#include <math.h>

#include "image.h"
#include "WriteImage.h"
#include "ReadImage.h"
#include "freeman_arg_parse.h" //A small utility I wrote for extracting command line args.

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define PI 3.14159265

bool FLAG_DEBUG = true;

std::string output_path = "../images";
// std::string image_input_path = "../images/";

int ClampPxVal(int val, int lo, int hi);
void PrintHistogram(std::map<int, float> hist);
void WriteImageToFile(std::string filename, ImageType& img);
void generateTestImage(int size, double** arr, int innerSize);

void fft(float data[], unsigned long nn, int isign);
void printArrayReal(float arr[], int SIZE);
void printArrayImag(float arr[], int SIZE);
void normalizeArray(float arr[], int valCount, int size);
void printMagnitude(float arr[], int size);
void DFT_WriteToCSV(float arr[], int SIZE, std::string filepath);

void fft2D(unsigned int N, unsigned int M, ImageType& i_real, ImageType& i_imag, int isign);


int main() 
{
  //create test image. create empty imaginary image object (init vals to 0)

  const int TEST_SIZE = 512;
  int whiteSquare = 32;
  double** testImage;
  testImage = new double* [TEST_SIZE];

  for (int i = 0; i < TEST_SIZE; i++)
    testImage[i] = new double[TEST_SIZE];

  generateTestImage(TEST_SIZE, testImage, whiteSquare);

  ImageType img_real(TEST_SIZE, TEST_SIZE, 255);
  for(int i=0; i<TEST_SIZE; i++)
  {
    for(int j=0; j<TEST_SIZE; j++)
    {
      img_real.setPixelVal(i,j,testImage[i][j]);

      double val;
      img_real.getPixelVal(i,j,val);
    }
  }

  ImageType img_imag;
  img_imag.CopyImageData(img_real);
  for(int i=0; i<TEST_SIZE; i++)
  {
    for(int j=0; j<TEST_SIZE; j++)
    {
      img_imag.setPixelVal(i,j,0);
    }
  }

  WriteImageToFile(output_path + "/test_image_raw.pgm", img_real);

  //2d ffts:
  //forward t
  fft2D(TEST_SIZE, TEST_SIZE, img_real, img_imag, -1);

  WriteImageToFile(output_path + "/test_image_real_frequency_domain.pgm", img_real);

  //backward t
  fft2D(TEST_SIZE, TEST_SIZE, img_real, img_imag, 1);


  WriteImageToFile(output_path + "/test_image_fwd_bck_transformed.pgm", img_real);


  // fft(testArr, realCount, -1); // forward fft
  // normalizeArray(testArr, realCount, SIZE);
  // //Save DFT Data (aka magnitude, real, imaginary, phase)
  // DFT_WriteToCSV(testArr, SIZE, output_path);
  // fft(testArr, realCount, 1); // inverse fft



  return 0;
}

void fft2D(unsigned int N, unsigned int M, ImageType& i_real, ImageType& i_imag, int isign)
{
  unsigned int SIZE = 2 * M + 1;
  float *arr = new float[SIZE];


  for(int k=0; k<SIZE; k++) arr[k] = 0.0f;


  //compute rows
  for(int i=0; i<N; i++)
  {
    // std::cout << i << std::endl;

    //load values into work array
    float* ptr_r = arr;
    ptr_r += 1;
    for(int j=0, k=1, l=2; j<M; j++, k+=2, l+=2)
    {
      double real, imag;
      i_real.getPixelVal(i,j,real); 
      i_imag.getPixelVal(i,j,imag);

      // *ptr_r = (float)real;
      // ptr_r += 2;
      arr[k] = real;
      arr[l] = imag;

      // if(i==255) std::cout << real << " ";
      // if(i == 255 && j==N-1)
      // {
      //   for(int b=0; b<50; b++) {std::cout << arr[b] << " ";}
      //   printArrayReal(arr, SIZE);
      //   printArrayImag(arr, SIZE);
      // }
    }

    //compute dft (row)
    fft(arr, M, isign);

    //copy values back into image row
    for(int j=0; j<M; j++)
    {
      double real, imag;
      real = arr[2 * j + 1];
      imag = arr[2 * j + 2];
      i_real.setPixelVal(i,j,real);
      i_imag.setPixelVal(i,j,imag);
    }


    //clear work array    
    for(int k=0; k<SIZE; k++) arr[k] = 0;
  }


  //compute columns
  for(int j=0; j<M; j++)
  {
    //load values into work array
    for(int i=0, k=1, l=2; i<N; i++, k+=2, l+=2)
    {
      double real, imag;
      i_real.getPixelVal(i,j,real);
      i_imag.getPixelVal(i,j,imag);
      arr[k] = real;
      arr[l] = imag;
    }

    //compute dft (column)
    fft(arr, N, isign);

    //copy values back into image col   
    for(int i=0; i<N; i++) 
    {
      double real, imag;
      real = arr[2 * i + 1];
      imag = arr[2 * i + 2];
      i_real.setPixelVal(i,j,real);
      i_imag.setPixelVal(i,j,imag);
    }
    //clear work array    
    for(int k=0; k<SIZE; k++) arr[k] = 0;
  }

  // return;

  //TODO: reenable
  if(isign < 0)
  {
    for(int i=0; i<N; i++) 
    {
      //clear work array
      for(int k=0; k<SIZE; k++) arr[k] = 0;

      //copy into work array
      for(int j=0; j<M; j++)
        {
        double real, imag;
        real = arr[2 * j + 1];
        imag = arr[2 * j + 2];
        i_real.setPixelVal(i,j,real);
        i_imag.setPixelVal(i,j,imag);        
      }

      //normalize
      normalizeArray(arr, N*M, SIZE);   

      //copy back into image storage
      for(int j=0; j<M; j++)
      {
      double real, imag;
      real = arr[2 * j + 1];
      imag = arr[2 * j + 2];
      i_real.setPixelVal(i,j,real);
      i_imag.setPixelVal(i,j,imag);
    }
    }
  }

  delete [] arr;
}

void generateTestImage(int size, double** arr, int innerSize)
{
  int leftBound = ((size/2) - (innerSize/2));
  int upperBound = ((size/2) - (innerSize/2));
  int rightBound = ((size/2) + (innerSize/2)) -1;
  int lowerBound = ((size/2) + (innerSize/2)) -1;

  // DEBUG ///////////////////////////////////////////
  // for(int i=0; i<size; i++)
  // {
  //   for (int j = 0; j < size; j++)
  //   {
  //     arr[i][j] = 255;
  //     // std::cout << i << ", " << j << ", val: " ;
  //     // std::cout << arr[i][j] << " | ";
  //   }
  // }

  for (int i = 0; i < size; i++)
  {
    for (int j = 0; j < size; j++)
    {
      if ( (i >= upperBound && i <= lowerBound) && (j >= leftBound && j <= rightBound) )
        arr[i][j] = 255;
      else
        arr[i][j] = 0;
    }
  }
}

void DFT_WriteToCSV(float arr[], int SIZE, std::string filepath)
{
	std::cout << "Writing DFT values as .csv files in directory: " << filepath << "\n";

  std::ofstream os;
  os.open (filepath + "/DFT_Real.csv");
  for (int i = 1; i < SIZE; i = i + 2)
    os << arr[i] << "\n";
  os.close();
  os.open (filepath + "/DFT_Imaginary.csv");
  for (int i = 2; i < SIZE; i = i + 2)
    os << arr[i] << "\n";
  os.close();
  os.open (filepath + "/DFT_Magnitude.csv");
  int i = 1;
  while (i < SIZE)
  {
    os << sqrt(pow(arr[i], 2) + pow(arr[i+1], 2)); 
    // |F(u)| = sqrt (R(u)^2 + I(u)^2)

    i += 2;
    os << "\n";
  }
  os.close();
  os.open (filepath + "/DFT_Phase.csv");
  i = 1;
  while (i < SIZE)
  {
    os << atan2(arr[i], arr[i+1]); 

    i += 2;
    os << "\n";
  }
  os.close();
}

void printArrayReal(float arr[], int SIZE)
{
  std::cout << "Array, real components: \n";
  for (int i = 1; i < SIZE; i = i + 2)
    std::cout << arr[i] << " ";
  std::cout << "\n\n";
}

void printArrayImag(float arr[], int SIZE)
{
  std::cout << "Array, imagninary components: \n";
  for (int i = 2; i < SIZE; i = i + 2)
    std::cout << arr[i] << " ";
  std::cout << "\n\n";
}

void fft(float data[], unsigned long nn, int isign)
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	float tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) {
		if (j > i) {
			SWAP(data[j],data[i]);
			SWAP(data[j+1],data[i+1]);
		}
		m=n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}


void normalizeArray(float arr[], int valCount, int SIZE)
{
  for (int i = 0; i < SIZE; i++)
    arr[i] = arr[i] * 1/valCount;
}

void printMagnitude(float arr[], int size)
{
  int i = 1;
  std::cout << "Magnitude of array: \n";
  while (i < size)
  {
    std::cout << sqrt(pow(arr[i], 2) + pow(arr[i+1], 2)); 
    // |F(u)| = sqrt (R(u)^2 + I(u)^2)

    i += 2;
    std::cout << "  ";
  }
  std::cout << "\n\n";
}

void PrintHistogram(std::map<int, float> hist)
{
  for(auto it : hist)
  {
    std::cout << it.first << " " << it.second << std::endl;
  }
  std::cout << std::endl;
}

void WriteImageToFile(std::string filename, ImageType& img)
{
  std::string out_file = filename;
  char *cstr = new char[out_file.length() + 1];
  strcpy(cstr, out_file.c_str());

  std::writeImage(cstr, img);

  std::cout << " * Saved image: " << out_file << "\n";

  delete [] cstr;
}

int ClampPxVal(int val, int lo, int hi)
{
  if(val < lo) return lo;
  else if (val > hi) return hi;
  else return val;
}

#undef SWAP
