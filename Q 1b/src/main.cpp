#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define PI 3.14159265

std::string output_path = "../DFT_Data";

void fft(float data[], unsigned long nn, int isign);
void printArrayReal(float arr[], int SIZE);
void printArrayImag(float arr[], int SIZE);
void normalizeArray(float arr[], int valCount, int size);
void printMagnitude(float arr[], int size);
void DFT_WriteToCSV(float arr[], int SIZE, std::string filepath);

int main() 
{
  const unsigned int N = 128;

  const unsigned int SIZE = 2 * N + 1;
  int realCount = N; // quantity of real values in array
  float testArr[SIZE];
  for(int i=0; i < SIZE; i++)
  {
    testArr[i] = 0;
  }
  // Sample cos function cos(2 * pi * 8 * x)
  for(int i=0; i < N; i++)
  {
    testArr[2 * i + 1] = cos(16.0 * PI * ((double)i * .125 / N));
  }

  // printArrayReal(testArr, SIZE);
  // printArrayImag(testArr, SIZE);

  fft(testArr, realCount, -1); // forward fft
  // should be: [13/4, 1/4 (-2 + j), -1/4, -1/4 (2 + j)]
  // AKA:       [3.25,      imag,   -0.25,      imag]
  normalizeArray(testArr, realCount, SIZE);
  DFT_WriteToCSV(testArr, SIZE, output_path);

  // printArrayReal(testArr, SIZE);
  // printArrayImag(testArr, SIZE);
  // printMagnitude(testArr, SIZE);

  //Save DFT Data (aka magnitude, real, imaginary, phase later)

  fft(testArr, realCount, 1); // inverse fft
  
  printArrayReal(testArr, SIZE);


  return 0;
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



#undef SWAP
