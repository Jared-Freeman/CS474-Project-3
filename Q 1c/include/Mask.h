
#include <array>
#include "image.h"

class ImageMask 
{
	public:
		//ctors
		ImageMask();
		ImageMask(int rows, int cols, int** intensity_values);
		ImageMask(int rows, int cols, float** intensity_values);
		ImageMask(ImageType& mask_image);
		
		//dtor
		~ImageMask();
		
		//methods
		void ApplyMask(ImageType& source_image, ImageType& output_image
			, bool flag_normalize = true);
		
	private:
		//members
		float** m_intensity_values;
		int m_rows;
		int m_cols;
		
		//private methods
		void AllocateMemory();
		void ClampIndex(int& val, int lo, int hi);
		int GetWeightedSum(
				int N
			, int M
			, int i
			, int j
			, ImageType& ref_image
			, bool flag_normalize = true
			);
		
};
