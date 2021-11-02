
#include "Mask.h"

#include <iostream>

ImageMask::ImageMask():
    m_rows(0)
,   m_cols(0)
,   m_intensity_values(nullptr)
{}

ImageMask::ImageMask(int rows, int cols, int** intensity_values):
    m_rows(rows)
,   m_cols(cols)
{
    AllocateMemory();
    for(int i=0; i<rows; i++)
    {
        for(int j=0; j<cols; j++)
        {
            m_intensity_values[i][j] = intensity_values[i][j];            
        }
    }
}

ImageMask::ImageMask(int rows, int cols, float** intensity_values):
    m_rows(rows)
,   m_cols(cols)
{
    AllocateMemory();
    for(int i=0; i<rows; i++)
    {
        for(int j=0; j<cols; j++)
        {
            m_intensity_values[i][j] = intensity_values[i][j];            
        }
    }
}

ImageMask::ImageMask(ImageType& mask_image)
{
    int Q;
    mask_image.getImageInfo(m_rows, m_cols, Q);
    AllocateMemory();

    int val = 0;
    for(int i=0; i < m_rows; i++)
    {
        for(int j=0; j < m_cols; j++)
        {
            mask_image.getPixelVal(i, j, val);
            m_intensity_values[i][j] = val;
        }
    }
    
}

ImageMask::~ImageMask()
{
    if(m_intensity_values != nullptr)
    {
        int i;
        for(i=0; i < m_rows ; i++)
            delete [] m_intensity_values[i];
        delete [] m_intensity_values;
    }
}

// should be called after rows, cols init
void ImageMask::AllocateMemory()
{
    m_intensity_values = new float* [m_rows];
    for(int i=0; i<m_rows; i++)
    {
        m_intensity_values[i] = new float[m_cols];
    }
}

void ImageMask::ApplyMask(ImageType& source_image, ImageType& output_image, bool flag_normalize)
{

    //init source image vals
    int N, M, Q;
    source_image.getImageInfo(N, M, Q);

    //resize output image to match source image
    output_image.CopyImageData(source_image);

    int weighted_sum = 0;
    //for each pixel value in source_image
    for(int i=0; i < N; i++)
    {
        for(int j=0; j < M; j++)
        {
            //Calculate weighted sum by sampling all values in window specified my mask dims
            weighted_sum = GetWeightedSum(N, M, i, j, source_image, flag_normalize);

            //Store calculated value in output_image
            output_image.setPixelVal(i, j, weighted_sum);
        }
    }
    
}


void ImageMask::ClampIndex(int& val, int lo, int hi)
{
    if(lo > hi) 
    {
        return;
    }
    else if(val < lo)
    {
        val = lo;
    }
    else if(val > hi)
    {
        val = hi;
    }
}

// desc: Returns a rounded, weighted sum of pixel values in the ref image as specified
// by this Mask
// note: If a Mask is even, we skew "left, up" meaning the window will be formed
// with more samples leftward and upward
int ImageMask::GetWeightedSum(int N, int M, int i, int j, ImageType& ref_image, bool flag_normalize)
{
    float partial_sum = 0;
    int ox_s, ox_e, oy_s, oy_e; //offsets (start, end indices for x and y)
    int cur_weight = 0;
    int xm_start = 0, ym_start = 0;

    ox_s = i - (int)(m_rows / 2);
    ox_e = i + (int)(m_rows / 2);
    //check for off-by-1 (because mask rows even)
    if(m_rows % 2 == 0) ox_e--;

    ClampIndex(ox_s, 0, N);
    xm_start = m_rows - (ox_e - ox_s) - 1;
    ClampIndex(ox_e, 0, N);

    oy_s = j - (int)(m_cols / 2);
    oy_e = j + (int)(m_cols / 2);
    //check for off-by-1 (because mask cols even)
    if(m_cols % 2 == 0) oy_e--;

    ClampIndex(oy_s, 0, M);
    ym_start = m_cols - (oy_e - oy_s) - 1;
    ClampIndex(oy_e, 0, M);

    // std::cout << "" << ox_s << "," << ox_e << "|" << oy_s << "," << oy_e << "|";
    // std::cout << xm_start << "," << ym_start << "  ";// << std::endl;

    // scan across window bounds, add weighted partial sums
    int px_val = 0;
    float normalization_value = 0; //we don't hold a constant for this because we sometimes "zero out" sections of mask vals

    for(int k = ox_s, x = xm_start; k < ox_e; k++, x++)
    {
        for(int l = oy_s, y = ym_start; l < oy_e; l++, y++)
        {
            ref_image.getPixelVal(k, l, px_val);
            partial_sum += px_val * m_intensity_values[x][y];
            normalization_value += m_intensity_values[x][y];
        }
    }

    if(flag_normalize) partial_sum /= normalization_value;

    return (int)partial_sum;
}