#include <iostream>
#include <cmath>
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void fft(float data[], unsigned long nn, int isign);
void printArrayReal(float arr[], int SIZE);
void printArrayImag(float arr[], int SIZE);
void normalizeArray(float arr[], int valCount, int size);
void printMagnitude(float arr[], int size);

int main() 
{
  const unsigned int SIZE = 9;
  int realCount = 4; // quantity of real values in array
  float testArr[SIZE] = {0, 2, 0, 3, 0, 4, 0, 4, 0};

  printArrayReal(testArr, SIZE);

  fft(testArr, realCount, -1); // forward fft
  // should be: [13/4, 1/4 (-2 + j), -1/4, -1/4 (2 + j)]
  // AKA:       [3.25,      imag,   -0.25,      imag]
  normalizeArray(testArr, realCount, SIZE);

  printArrayReal(testArr, SIZE);
  printArrayImag(testArr, SIZE);
  printMagnitude(testArr, SIZE);

  fft(testArr, realCount, 1); // inverse fft
  
  printArrayReal(testArr, SIZE);

  return 0;
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
