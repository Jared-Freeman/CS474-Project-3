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
#include "Mask.h"
#include "freeman_arg_parse.h" //A small utility I wrote for extracting command line args.

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define PI 3.14159265

bool FLAG_DEBUG = true;

std::string output_path = "../images";
// std::string image_input_path = "../images/";

int ClampPxVal(int val, int lo, int hi);
void PrintHistogram(std::map<int, float> hist);
void WriteImageToFile(std::string filename, ImageType& img);

void fft(float data[], unsigned long nn, int isign);
void printArrayReal(float arr[], int SIZE);
void printArrayImag(float arr[], int SIZE);
void normalizeArray(float arr[], int valCount, int size);
void printMagnitude(float arr[], int size);
void DFT_WriteToCSV(float arr[], int SIZE, std::string filepath);

void fft2D(unsigned int N, unsigned int M, ImageType i_real, ImageType i_imag, int isign);


int main() 
{
  //create test image. create empty imaginary image object (init vals to 0)

  //2d fft





  // fft(testArr, realCount, -1); // forward fft
  // normalizeArray(testArr, realCount, SIZE);
  // //Save DFT Data (aka magnitude, real, imaginary, phase)
  // DFT_WriteToCSV(testArr, SIZE, output_path);
  // fft(testArr, realCount, 1); // inverse fft



  return 0;
}

void fft2D(unsigned int N, unsigned int M, ImageType i_real, ImageType i_imag, int isign)
{
  unsigned int SIZE = N * M + 1;
  float arr[SIZE];
  for(int k=0; k<SIZE; k++) arr[k] = 0;


  //compute rows
  for(int i=0; i<N; i++)
  {
    //load values into work array
    for(int j=0; j<M; j++)
    {
      double real, imag;
      i_real.getPixelVal(i,j,real);
      i_imag.getPixelVal(i,j,imag);
      arr[2 * j + 1] = real;
      arr[2 * j + 2] = imag;
    }

    //compute dft
    fft(arr, SIZE, isign);

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
    for(int i=0; i<N; i++)
    {
      double real, imag;
      i_real.getPixelVal(i,j,real);
      i_imag.getPixelVal(i,j,imag);
      arr[2 * j + 1] = real;
      arr[2 * j + 2] = imag;
    }

    //compute dft
    fft(arr, SIZE, isign);

    //copy values back into image col   
    for(int i=0; i<N; i++) 
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

  //normalize after operations? during?
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
      normalizeArray(arr, M, SIZE);   

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
