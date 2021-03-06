#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <cstring>
#include <map>

#include "image.h"
#include "WriteImage.h"
#include "ReadImage.h"
#include "freeman_arg_parse.h" //A small utility I wrote for extracting command line args.

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define PI 3.14159265

bool FLAG_DEBUG = true;

std::string output_path = "../images";

int ClampPxVal(int val, int lo, int hi);
void WriteImageToFile(std::string filename, ImageType& img);
void ReadImageFromFile(std::string filename, ImageType& img);

int ProcessTestImages(int argc, char** argv);

void fft(double data[], unsigned long nn, int isign);
void normalizeArray(double arr[], int valCount, int size);
void DFT_WriteToCSV(double arr[], int SIZE, std::string filepath);

void fft2D(unsigned int N, unsigned int M, ImageType& i_real, ImageType& i_imag, int isign);

void zeroMagnitude(unsigned int N, unsigned int M, ImageType& i_real, ImageType& i_imag);

// NOTES: This program does NOT check if images have power-of-two dimensions. 
//        fft() will only work under these conditions.
//        Furthermore, as of 11/7/2021 this program has not been tested with N=/=M dimensioned images yet.
int main(int argc, char** argv)
{
    // Setup: Create image classes to hold real, imaginary, and magnitude data of processed image
    ImageType image;
    ReadImageFromFile(output_path + "/lenna.pgm", image);
    int N, M, Q;
    double q;

    image.getImageInfo(N, M, Q);

    ImageType img_real(N, M, 255);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            image.getPixelVal(i, j, q);
            img_real.setPixelVal(i, j, q);
        }
    }

    ImageType img_imag;
    img_imag.CopyImageData(img_real);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            img_imag.setPixelVal(i, j, 0);
        }
    }

    // Transformations
    //2d ffts:
    //forward t
    fft2D(N, M, img_real, img_imag, -1);

    zeroMagnitude(N, M, img_real, img_imag);

    //backward t
    fft2D(N, M, img_real, img_imag, 1);

    img_real.RemapPixelValues();

    WriteImageToFile(output_path + "/lenna_zeroMagnitude.pgm", img_real);

    // Test using input file: ////////////////////////////////////////////////
    ProcessTestImages(argc, argv);

    return 0;
}

//TODO: PROBABLY only works if N==M currently.
//Need to double check all ranges and uses of those vars to ensure they're correct
//Also the work array would need to be resized (SIZE == 2N + 1) for the columns fft
void fft2D(unsigned int N, unsigned int M, ImageType& i_real, ImageType& i_imag, int isign)
{
    unsigned int SIZE = 2 * M + 1;
    double* arr = new double[SIZE];

    for (int k = 0; k < SIZE; k++) arr[k] = 0.0f;

    //compute rows
    for (int i = 0; i < N; i++)
    {
        // std::cout << i << std::endl;

        //load values into work array
        double* ptr_r = arr;
        ptr_r += 1;
        for (int j = 0, k = 1, l = 2; j < M; j++, k += 2, l += 2)
        {
            double real, imag;
            i_real.getPixelVal(i, j, real);
            i_imag.getPixelVal(i, j, imag);

            arr[k] = real;
            arr[l] = imag;
        }

        //compute dft (row)
        fft(arr, M, isign);

        //copy values back into image row
        for (int j = 0; j < M; j++)
        {
            double real, imag;
            real = arr[2 * j + 1];
            imag = arr[2 * j + 2];

            if (isign < 0) real /= N;
            if (isign < 0) imag /= N;

            i_real.setPixelVal(i, j, real);
            i_imag.setPixelVal(i, j, imag);
        }


        //clear work array    
        for (int k = 0; k < SIZE; k++) arr[k] = 0;
    }


    //compute columns
    for (int j = 0; j < M; j++)
    {
        //load values into work array
        for (int i = 0, k = 1, l = 2; i < N; i++, k += 2, l += 2)
        {
            double real, imag;
            i_real.getPixelVal(i, j, real);
            i_imag.getPixelVal(i, j, imag);


            if (isign < 0) real /= N;
            if (isign < 0) imag /= N;

            arr[k] = real;
            arr[l] = imag;
        }

        //compute dft (column)
        fft(arr, N, isign);
        
      //copy values back into image col   
        for (int i = 0; i < N; i++)
        {
            double real, imag;
            real = arr[2 * i + 1];
            imag = arr[2 * i + 2];
            i_real.setPixelVal(i, j, real);
            i_imag.setPixelVal(i, j, imag);
        }
        //clear work array    
        for (int k = 0; k < SIZE; k++) arr[k] = 0;
    }

    return;
}

void DFT_WriteToCSV(double arr[], int SIZE, std::string filepath)
{
    std::cout << "Writing DFT values as .csv files in directory: " << filepath << "\n";

    std::ofstream os;
    os.open(filepath + "/DFT_Real.csv");
    for (int i = 1; i < SIZE; i = i + 2)
        os << arr[i] << "\n";
    os.close();
    os.open(filepath + "/DFT_Imaginary.csv");
    for (int i = 2; i < SIZE; i = i + 2)
        os << arr[i] << "\n";
    os.close();
    os.open(filepath + "/DFT_Magnitude.csv");
    int i = 1;
    while (i < SIZE)
    {
        os << sqrt(pow(arr[i], 2) + pow(arr[i + 1], 2));
        // |F(u)| = sqrt (R(u)^2 + I(u)^2)

        i += 2;
        os << "\n";
    }
    os.close();
    os.open(filepath + "/DFT_Phase.csv");
    i = 1;
    while (i < SIZE)
    {
        os << atan2(arr[i], arr[i + 1]);

        i += 2;
        os << "\n";
    }
    os.close();
}

void fft(double data[], unsigned long nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2) {
        if (j > i) {
            SWAP(data[j], data[i]);
            SWAP(data[j + 1], data[i + 1]);
        }
        m = n >> 1;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax = 2;
    while (n > mmax) {
        istep = mmax << 1;
        theta = isign * (6.28318530717959 / mmax);
        wtemp = sin(0.5 * theta);
        wpr = -2.0 * wtemp * wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m = 1; m < mmax; m += 2) {
            for (i = m; i <= n; i += istep) {
                j = i + mmax;
                tempr = wr * data[j] - wi * data[j + 1];
                tempi = wr * data[j + 1] + wi * data[j];
                data[j] = data[i] - tempr;
                data[j + 1] = data[i + 1] - tempi;
                data[i] += tempr;
                data[i + 1] += tempi;
            }
            wr = (wtemp = wr) * wpr - wi * wpi + wr;
            wi = wi * wpr + wtemp * wpi + wi;
        }
        mmax = istep;
    }
}


void normalizeArray(double arr[], int valCount, int SIZE)
{
    for (int i = 0; i < SIZE; i++)
        arr[i] = arr[i] * 1 / valCount;
}

void WriteImageToFile(std::string filename, ImageType& img)
{
    std::string out_file = filename;
    char* cstr = new char[out_file.length() + 1];
    strcpy(cstr, out_file.c_str());

    std::writeImage(cstr, img);

    std::cout << " * Saved image: " << out_file << "\n";

    delete[] cstr;
}

void ReadImageFromFile(std::string filename, ImageType& img)
{
    std::string in_file = filename;
    char* cstr = new char[in_file.length() + 1];
    strcpy(cstr, in_file.c_str());

    std::readImage(cstr, img);

    std::cout << " * Saved image: " << in_file << "\n";

    delete[] cstr;
}

int ClampPxVal(int val, int lo, int hi)
{
    if (val < lo) return lo;
    else if (val > hi) return hi;
    else return val;
}


int ProcessTestImages(int argc, char** argv)
{
    //TODO: Print description of process to console (per program)

    //Extract args into vector
    std::vector<std::string> args;
    for (int i = 1; i < argc; i++)
    {
        std::string next_element(argv[i]);
        args.push_back(next_element);
    }

    //Fill data structures using args
    std::vector<std::string> imagePaths = ExtractArgs("-in", args);
    std::vector<std::string> outputPaths = ExtractArgs("-out", args);
    if (outputPaths.size() > 0 && outputPaths[0].length() > 0)
    {
        if (outputPaths[0][outputPaths[0].length() - 1] != '/')
        {
            outputPaths[0] = outputPaths[0] + "/";
            //std::cout << outputPaths[0] << "\n";
        }
    }
    if (outputPaths.size() > 0)
    {
        std::cout << "Output paths specification is currently disabled for this program. Instead change the string \"output_path\" in main.cpp.\n";
    }


    //Process each image
    for (int i = 0; i < imagePaths.size(); i++)
    {
        std::cout << "_________________________\n";
        std::cout << "image " << i << ": \"" << imagePaths[i] << "\"\n";

        ImageType next_image;
        char* cstr = new char[imagePaths[i].length() + 1];
        strcpy(cstr, imagePaths[i].c_str());
        std::readImage(cstr, next_image);

        // DO STUFF
        int N, M, Q;
        next_image.getImageInfo(N, M, Q);

        ImageType img_imag;
        img_imag.CopyImageData(next_image);
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                img_imag.setPixelVal(i, j, 0);
            }
        }


        // GENERATE OUTPUT FILENAME CONVENTION
        //determine original file name from path string
        std::string original_filename = "";
        for (int l = imagePaths[i].length(); imagePaths[i][l] != '/' && l >= 0; l--)
        {
            if (imagePaths[i][l] != '/')
            {
                //std::cout << imagePaths[i][l] << std::endl;

                std::string temp;
                temp += imagePaths[i][l];
                original_filename.insert(0, temp);
            }
            //throw out extension
            if (imagePaths[i][l] == '.')
            {
                original_filename = "";
            }
        }
        std::string out_file = output_path + "/" + original_filename;


        WriteImageToFile(out_file + "_test_image_raw.pgm", next_image);
        //2d ffts:
        //forward t
        fft2D(N, M, next_image, img_imag, -1);
        WriteImageToFile(out_file + "_test_image_real_frequency_domain.pgm", next_image);
        //backward t
        fft2D(N, M, next_image, img_imag, 1);
        WriteImageToFile(out_file + "_test_image_fwd_bck_transformed.pgm", next_image);

        // END DO STUFF



        delete[] cstr;

        std::cout << "\n";
    }

    return 0;
}

void zeroMagnitude(unsigned int N, unsigned int M, ImageType& i_real, ImageType& i_imag)
{
    double Q_real, Q_imag, theta;

    for (int u = 0; u < N; u++)
    {
        for (int v = 0; v < M; v++)
        {
            i_real.getPixelVal(u, v, Q_real);
            i_imag.getPixelVal(u, v, Q_imag);
            theta = atan2(Q_imag, Q_real);

            i_real.setPixelVal(u, v, cos(theta));
            i_imag.setPixelVal(u, v, sin(theta));
        }
    }
}


#undef SWAP
