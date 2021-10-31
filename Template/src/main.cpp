 /*
  * @description: 
  * 
  *
  * @author: 
  * 
  *
  */
  
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cstring>
#include <map>
#include <math.h>

#include "image.h"
#include "WriteImage.h"
#include "ReadImage.h"
#include "Mask.h"

#include "freeman_arg_parse.h" //A small utility I wrote for extracting command line args.

bool FLAG_DEBUG = true;

int ClampPxVal(int val, int lo, int hi);
void PrintHistogram(std::map<int, float> hist);
void WriteImageToFile(std::string filename, ImageType& img);
/*
int main(int argc, char** argv)
{
  //TODO: Print description of process to console (per program)

  //Extract args into vector
  std::vector<std::string> args;
  for(int i=1; i < argc; i++)
  {
    std::string next_element(argv[i]);  
    args.push_back(next_element);
  }
  
	//Fill data structures using args
  std::vector<std::string> imagePaths = ExtractArgs("-in", args);
  std::vector<std::string> outputPaths = ExtractArgs("-out", args);
  if(outputPaths.size() > 0 && outputPaths[0].length() > 0)
  {
    if(outputPaths[0][outputPaths[0].length()-1] != '/')
    {
      outputPaths[0] = outputPaths[0] + "/";
      //std::cout << outputPaths[0] << "\n";
    }
  }
  if(outputPaths.size() > 1)
  {
    std::cout << "Can only specify one output path. Output will be saved in the first output directory.\n";
  }  
  if(outputPaths.size() < 1)
  {
    std::cout << "Please specify an output directory using the \"-out <output_path>\" args.\n";
    return 0;
  } 
  
    
  //Process each image
  for(int i=0; i < imagePaths.size(); i++)
  {
    std::cout << "_________________________\n";
    std::cout << "image " << i << ": \"" << imagePaths[i] << "\"\n";
    
    ImageType next_image;
    ImageType smoothed_image;
    ImageType gmsk_image;
    ImageType out_image;
    char *cstr = new char[imagePaths[i].length() + 1];
    strcpy(cstr, imagePaths[i].c_str());
    std::readImage(cstr, next_image);

    // DO STUFF
    // ...

    // GENERATE OUTPUT FILENAME CONVENTION
    //determine original file name from path string
    std::string original_filename = "";
    for(int l = imagePaths[i].length(); imagePaths[i][l] != '/' && l >= 0; l--)
    {
      if(imagePaths[i][l] != '/')
      {
        //std::cout << imagePaths[i][l] << std::endl;

        std::string temp;
        temp += imagePaths[i][l];
        original_filename.insert(0, temp);
      }
      //throw out extension
      if(imagePaths[i][l] == '.')
      {
        original_filename = "";
      }
    }      
    std::string out_file = outputPaths[0] + original_filename;


    // MASK GENERATION

    int** arr;
    int sz;
    int arr_i;


    // Gaussian 7x7 mask
    sz = 7;
    arr = new int*[sz];
    arr_i = 0;
    
    arr[arr_i++] = new int[sz]  {1,1,2,2,2,1,1};
    arr[arr_i++] = new int[sz]  {1,2,2,4,2,2,1};
    arr[arr_i++] = new int[sz]  {2,2,4,8,4,2,2};
    arr[arr_i++] = new int[sz]  {2,4,8,16,8,4,2};
    arr[arr_i++] = new int[sz]  {2,2,4,8,4,2,2};
    arr[arr_i++] = new int[sz]  {1,2,2,4,2,2,1};
    arr[arr_i++] = new int[sz]  {1,1,2,2,2,1,1};

    ImageMask mask_gaussian_7(sz, sz, arr);
    mask_gaussian_7.ApplyMask(next_image, smoothed_image);


    gmsk_image.CopyImageData(next_image);
    int N, M, Q;
    next_image.getImageInfo(N,M,Q);

    for(int i=0; i<N; i++)
    {
      for(int j=0; j<M; j++)
      {
        int v1,v2;
        next_image.getPixelVal(i,j,v1);
        smoothed_image.getPixelVal(i,j,v2);

        // std::cout << "[" << v1 - v2 << "]  ";
        gmsk_image.setPixelVal(i,j,v1-v2);
      }
    }

    out_image.CopyImageData(next_image);
    for(int i=0; i<N; i++)
    {
      for(int j=0; j<M; j++)
      {
        int v1,v2;
        next_image.getPixelVal(i,j,v1);
        gmsk_image.getPixelVal(i,j,v2);

        //std::cout << "[" << v1 << " " << v2 << "]  ";
        out_image.setPixelVal(i,j,ClampPxVal(v1+v2, 0, Q));
      }
    }
    WriteImageToFile(out_file + "_unsharp.pgm", out_image);

    for(int k=2; k <= 8; k *= 2)
    {
      out_image.CopyImageData(next_image);
      for(int i=0; i<N; i++)
      {
        for(int j=0; j<M; j++)
        {
          int v1,v2;
          next_image.getPixelVal(i,j,v1);
          gmsk_image.getPixelVal(i,j,v2);

          //std::cout << "[" << v1 << " " << v2 << "]  ";
          out_image.setPixelVal(i,j,ClampPxVal(v1+(k * v2), 0, Q));
        }
      }
      WriteImageToFile(out_file + "_hi_boost_" + std::to_string(k) + ".pgm", out_image);
    }

    for(int i=0; i<sz; i++)
      delete [] arr[i];
    delete [] arr;


    // END DO STUFF



    delete [] cstr;

    std::cout << "\n";
  }
  
  return 0;
}
*/
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