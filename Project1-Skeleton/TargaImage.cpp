///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.cpp                          Author:     Stephen Chenney
//                                              Modified:   Eric McDaniel
//                                              Date:       Fall 2004
//                                              Modified:   Feng Liu
//                                              Date:       Winter 2011
//                                              Why:        Change the library file 
//      Implementation of TargaImage methods.  You must implement the image
//  modification functions.
//
///////////////////////////////////////////////////////////////////////////////

#include "Globals.h"
#include "TargaImage.h"
#include "libtarga.h"
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <map>

using namespace std;

// constants
const int           RED             = 0;                // red channel
const int           GREEN           = 1;                // green channel
const int           BLUE            = 2;                // blue channel
const unsigned char BACKGROUND[3]   = { 0, 0, 0 };      // background color


// Computes n choose s, efficiently
double Binomial(int n, int s)
{
    double        res;

    res = 1;
    for (int i = 1 ; i <= s ; i++)
        res = (n - i + 1) * res / i ;

    return res;
}// Binomial


///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage() : width(0), height(0), data(NULL)
{}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h) : width(w), height(h)
{
   data = new unsigned char[width * height * 4];
   intensities = new float[width * height];
   for (int i = 0; i < width * height; i++) {
       intensities[i] = 0.0;
   }
   ClearToBlack();
}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables to values given.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h, unsigned char *d)
{
    int i;

    width = w;
    height = h;
    data = new unsigned char[width * height * 4];

    for (i = 0; i < width * height * 4; i++)
	    data[i] = d[i];

    intensities = new float [width * height];
    for (int i = 0; i < width * height; i++) {
            intensities[i] = 0.0;
    }
    
}// TargaImage

int compare(const void* a, const void* b)
{
    const float* x = (float*)a;
    const float* y = (float*)b;

    if (*x > *y)
        return 1;
    else if (*x < *y)
        return -1;

    return 0;
}


///////////////////////////////////////////////////////////////////////////////
//
//      Copy Constructor.  Initialize member to that of input
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(const TargaImage& image) 
{
   width = image.width;
   height = image.height;
   data = NULL; 
   if (image.data != NULL) {
      data = new unsigned char[width * height * 4];
      memcpy(data, image.data, sizeof(unsigned char) * width * height * 4);
   }
   
}


///////////////////////////////////////////////////////////////////////////////
//
//      Destructor.  Free image memory.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::~TargaImage()
{
    if (data)
        delete[] data;
    if (intensities)
        delete[] intensities;
}// ~TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Converts an image to RGB form, and returns the rgb pixel data - 24 
//  bits per pixel. The returned space should be deleted when no longer 
//  required.
//
///////////////////////////////////////////////////////////////////////////////
unsigned char* TargaImage::To_RGB(void)
{
    unsigned char   *rgb = new unsigned char[width * height * 3];
    int		    i, j;

    if (! data)
	    return NULL;

    // Divide out the alpha
    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = i * width * 4;
	    int out_offset = i * width * 3;

	    for (j = 0 ; j < width ; j++)
        {
	        RGBA_To_RGB(data + (in_offset + j*4), rgb + (out_offset + j*3));
	    }
    }

    return rgb;
}// TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Save the image to a targa file. Returns 1 on success, 0 on failure.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Save_Image(const char *filename)
{
    TargaImage	*out_image = Reverse_Rows();

    if (! out_image)
	    return false;

    if (!tga_write_raw(filename, width, height, out_image->data, TGA_TRUECOLOR_32))
    {
	    cout << "TGA Save Error: %s\n", tga_error_string(tga_get_last_error());
	    return false;
    }

    delete out_image;

    return true;
}// Save_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Load a targa image from a file.  Return a new TargaImage object which 
//  must be deleted by caller.  Return NULL on failure.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Load_Image(char *filename)
{
    unsigned char   *temp_data;
    TargaImage	    *temp_image;
    TargaImage	    *result;
    int		        width, height;

    if (!filename)
    {
        cout << "No filename given." << endl;
        return NULL;
    }// if

    temp_data = (unsigned char*)tga_load(filename, &width, &height, TGA_TRUECOLOR_32);
    if (!temp_data)
    {
        cout << "TGA Error: %s\n", tga_error_string(tga_get_last_error());
	    width = height = 0;
	    return NULL;
    }
    temp_image = new TargaImage(width, height, temp_data);
    free(temp_data);

    result = temp_image->Reverse_Rows();

    delete temp_image;

    return result;
}// Load_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::To_Grayscale()
{
    if (!data)
        return false;
    int k = 0;
    
    for (int i = 0; i < height; i++)
    {
        int in_offset = i * width * 4;
        for (int j = 0; j < width; j++)
        {
            unsigned char* channels = data + (in_offset + j * 4);

            unsigned char intensity = 0.299 * channels[0] + 0.587 * channels[1] + 0.114 * channels[2];
            for (int i = 0; i < 3; i++) {
                channels[i] = intensity;
            }
            intensities[k++] = intensity/256.0;
        }
    }

    return true;
}// To_Grayscale


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Uniform()
{
    if (!data)
        return false;

    const int initial_bits = 8;
    int levels[3] = { 8, 8, 4 };
    int scales[3];
    for (int i = 0; i < 3; i++) {
        int required_bits = log2(levels[i]);
        scales[i] = pow(2, (initial_bits - required_bits));
    }

    for (int i = 0; i < height; i++)
    {
        int in_offset = i * width * 4;
        for (int j = 0; j < width; j++)
        {
            unsigned char* channels = data + (in_offset + j * 4);
            for (int i = 0; i < 3; i++) {
                channels[i] = (unsigned char)((int)channels[i] * 1.0/scales[i]) * scales[i];
            }
        }
    }
    return true;
}// Quant_Uniform


///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Populosity()
{
    ClearToBlack();
    return false;
}// Quant_Populosity


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()
{
    if (!data)
        return false;

    for (int i = 0; i < height; i++)
    {
        int in_offset = i * width * 4;
        for (int j = 0; j < width; j++)
        {
            unsigned char* channels = data + (in_offset + j * 4);
            float intensity = 0.299 * channels[0] + 0.587 * channels[1] + 0.114 * channels[2];
            if (intensity/256 < 0.5 ) {
                channels[0] = 0;
                channels[1] = 0;
                channels[2] = 0;
            }
            else {
                channels[0] = 255;
                channels[1] = 255;
                channels[2] = 255;
            }
        }
    }

    return true;
}// Dither_Threshold


///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()
{
    if (!data)
        return false;

    To_Grayscale();
    //25 + (std::rand() % (63 - 25 + 1))
    for (int i = 0; i < height; i++)
    {
        int in_offset = i * width * 4;
        for (int j = 0; j < width; j++)
        {
            unsigned char* channels = data + (in_offset + j * 4);
            float intensity = channels[0]/256.0;
            intensity += 0.4 * rand() / RAND_MAX - 0.2;

            if (intensity < 0.5) {
                channels[0] = 0;
                channels[1] = 0;
                channels[2] = 0;
            }
            else {
                channels[0] = 255;
                channels[1] = 255;
                channels[2] = 255;
            }
        }
    }
    return true;
}// Dither_Random


///////////////////////////////////////////////////////////////////////////////
//
//      Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_FS()
{
    ClearToBlack();
    return false;
}// Dither_FS

//assumes that the image is grayscale
float calc_avg_intensity(TargaImage * img) {

    float intensity_sum = 0.0;
    for (int i = 0; i < img->height; i++)
    {
        int in_offset = i * img->width * 4;
        for (int j = 0; j < img->width; j++)
        {
            unsigned char* channels = img->data + (in_offset + j * 4);
            intensity_sum += channels[0]/256.0;
           
        }
    }
    return (intensity_sum / (img->height * img->width));
}


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Bright()
{
    if (!data)
        return false;

    To_Grayscale();
    float avg_intensity = calc_avg_intensity(this);
    int idx = (int)((1 - avg_intensity) * width * height) ;

    qsort(intensities, width * height, sizeof(float), compare);
    float threshold = intensities[idx];

    for (int i = 0; i < height; i++)
    {
        int in_offset = i * width * 4;
        for (int j = 0; j < width; j++)
        {
            unsigned char* channels = data + (in_offset + j * 4);
            float intensity = channels[0] / 256.0;
            if (intensity < threshold) {
                channels[0] = 0;
                channels[1] = 0;
                channels[2] = 0;
            }
            else {
                channels[0] = 255;
                channels[1] = 255;
                channels[2] = 255;
            }
        }
    }

    return true;
}// Dither_Bright


///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()
{
    float mask[4][4] = { { 0.7500,  0.3750,  0.6250,  0.2500}, {0.0625,  1.0000,  0.8750,  0.4375}, {0.5000,  0.8125,  0.9375,  0.1250}, {0.1875,  0.5625,  0.3125,  0.6875}};
    To_Grayscale();
    
    for (int i = 0; i < height; i++)
    {
        int in_offset = i * width * 4;
        for (int j = 0; j < width; j++)
        {
            unsigned char* channels = data + (in_offset + j * 4);
            float intensity = channels[0]/256.0;
            float threshold = mask[i % 4][j % 4];
            if (intensity < threshold) {
                channels[0] = 0;
                channels[1] = 0;
                channels[2] = 0;
            }
            else {
                channels[0] = 255;
                channels[1] = 255;
                channels[2] = 255;
            }
        }
    }
    return true;
        
}// Dither_Cluster


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Color()
{
    ClearToBlack();
    return false;
}// Dither_Color


///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Comp_Over(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout <<  "Comp_Over: Images not the same size\n";
        return false;
    }

    for (int i = 0; i < height; i++)
    {
        int in_offset = i * width * 4;
        for (int j = 0; j < width; j++)
        {
            unsigned char* channels = data + (in_offset + j * 4);
            unsigned char* pImgChannels = pImage->data + (in_offset + j * 4);
            float G = 1 - (channels[3] / 255.0);
            for (int k = 0; k < 4; k++) {
                channels[k] = (unsigned char)(channels[k] + G * pImgChannels[k]);
            }
        }
    }

    return true;

}// Comp_Over


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_In(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_In: Images not the same size\n";
        return false;
    }

    for (int i = 0; i < height; i++)
    {
        int in_offset = i * width * 4;
        for (int j = 0; j < width; j++)
        {
            unsigned char* channels = data + (in_offset + j * 4);
            unsigned char* pImgChannels = pImage->data + (in_offset + j * 4);
            float F = pImgChannels[3] / 255.0;
            float G = 0.0;
            for (int k = 0; k < 4; k++) {
                channels[k] = (unsigned char)(channels[k] + G * pImgChannels[k]);
            }
        }
    }
    return true;
}// Comp_In


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Out(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Out: Images not the same size\n";
        return false;
    }

    for (int i = 0; i < height; i++)
    {
        int in_offset = i * width * 4;
        for (int j = 0; j < width; j++)
        {
            unsigned char* channels = data + (in_offset + j * 4);
            unsigned char* pImgChannels = pImage->data + (in_offset + j * 4);
            float F = 1- (pImgChannels[3] / 255.0);
            float G = 0.0;
            for (int k = 0; k < 4; k++) {
                channels[k] = (unsigned char)(channels[k] + G * pImgChannels[k]);
            }
        }
    }
    return true;
}// Comp_Out


///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Atop(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Atop: Images not the same size\n";
        return false;
    }

    for (int i = 0; i < height; i++)
    {
        int in_offset = i * width * 4;
        for (int j = 0; j < width; j++)
        {
            unsigned char* channels = data + (in_offset + j * 4);
            unsigned char* pImgChannels = pImage->data + (in_offset + j * 4);
            float F = pImgChannels[3] / 255.0;
            float G = 1 - (channels[3] / 255.0);
            for (int k = 0; k < 4; k++) {
                channels[k] = (unsigned char)(channels[k] + G * pImgChannels[k]);
            }
        }
    }
    return true;
}// Comp_Atop


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Xor(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Xor: Images not the same size\n";
        return false;
    }
    for (int i = 0; i < height; i++)
    {
        int in_offset = i * width * 4;
        for (int j = 0; j < width; j++)
        {
            unsigned char* channels = data + (in_offset + j * 4);
            unsigned char* pImgChannels = pImage->data + (in_offset + j * 4);
            float F = 1- (pImgChannels[3] / 255.0);
            float G = 1 - (channels[3] / 255.0);
            for (int k = 0; k < 4; k++) {
                channels[k] = (unsigned char)(channels[k] + G * pImgChannels[k]);
            }
        }
    }
    return true;
    
}// Comp_Xor


///////////////////////////////////////////////////////////////////////////////
//
//      Calculate the difference bewteen this imag and the given one.  Image 
//  dimensions must be equal.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Difference(TargaImage* pImage)
{
    if (!pImage)
        return false;

    if (width != pImage->width || height != pImage->height)
    {
        cout << "Difference: Images not the same size\n";
        return false;
    }// if

    for (int i = 0 ; i < width * height * 4 ; i += 4)
    {
        unsigned char        rgb1[3];
        unsigned char        rgb2[3];

        RGBA_To_RGB(data + i, rgb1);
        RGBA_To_RGB(pImage->data + i, rgb2);

        data[i] = abs(rgb1[0] - rgb2[0]);
        data[i+1] = abs(rgb1[1] - rgb2[1]);
        data[i+2] = abs(rgb1[2] - rgb2[2]);
        data[i+3] = 255;
    }

    return true;
}// Difference


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Box()
{
    if (!data)
        return false;

    const int filter_size = 5;
    int box_filter[filter_size][filter_size] = { {1,1,1,1,1}, {1,1,1,1,1}, {1,1,1,1,1},  {1,1,1,1,1},  {1,1,1,1,1} };
    float denominator = 25.0;

    for (int i = 0; i < height-filter_size; i++)
    {
        int i_offset = i * width * 4;
        for (int j = 0; j < width-filter_size; j++)
        {
            int j_offset = i_offset + j * 4;

            int red = 0, green = 0, blue = 0;
            for (int h = 0; h < filter_size; h++) {
                int h_offset = j_offset + h* width * 4;
                for (int k = 0; k < filter_size; k++) {
                    unsigned char* channels = data + (h_offset + k * 4);
                    red += channels[0] * box_filter[h][k];
                    green += channels[1] * box_filter[h][k];
                    blue += channels[2] * box_filter[h][k];
                }
            }
            //convolve
            unsigned char* index = data + filter_size/2 *4 * width + j_offset + filter_size/2 * 4;
            index[0] = unsigned char(red / denominator);
            index[1] = unsigned char(green / denominator);
            index[2] = unsigned char(blue / denominator);
        }
    }

    return true;
}// Filter_Box



///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
    if (!data)
        return false;

    const int filter_size = 5;
    int barlett_filter[filter_size][filter_size] = { {1,2,3,2,1}, {2,4,6,4,2}, {3,6,9,6,3}, {2,4,6,4,2},  {1,2,3,2,1} };
    float denominator = 81.0;

    for (int i = 0; i < height - filter_size; i++)
    {
        int i_offset = i * width * 4;
        for (int j = 0; j < width - filter_size; j++)
        {
            int j_offset = i_offset + j * 4;

            int red = 0, green = 0, blue = 0;
            for (int h = 0; h < filter_size; h++) {
                int h_offset = j_offset + h * width * 4;
                for (int k = 0; k < filter_size; k++) {
                    unsigned char* channels = data + (h_offset + k * 4);
                    red += channels[0] * barlett_filter[h][k];
                    green += channels[1] * barlett_filter[h][k];
                    blue += channels[2] * barlett_filter[h][k];
                }
            }
            //convolve
            unsigned char* index = data + filter_size / 2 * 4 * width + j_offset + filter_size / 2 * 4;
            index[0] = unsigned char(red / denominator);
            index[1] = unsigned char(green / denominator);
            index[2] = unsigned char(blue / denominator);
        }
    }

    return true;
}// Filter_Bartlett


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
    if (!data)
        return false;

    const int filter_size = 5;
    int gauss_filter[filter_size][filter_size] = { {1,4,6,4,1}, {4,16,24,16,4}, {6,24,36,24,6}, {4,16,24,16,4}, {1,4,6,4,1} };
    float denominator = 256.0;

    for (int i = 0; i < height - filter_size; i++)
    {
        int i_offset = i * width * 4;
        for (int j = 0; j < width - filter_size; j++)
        {
            int j_offset = i_offset + j * 4;

            int red = 0, green = 0, blue = 0;
            for (int h = 0; h < filter_size; h++) {
                int h_offset = j_offset + h * width * 4;
                for (int k = 0; k < filter_size; k++) {
                    unsigned char* channels = data + (h_offset + k * 4);
                    red += channels[0] * gauss_filter[h][k];
                    green += channels[1] * gauss_filter[h][k];
                    blue += channels[2] * gauss_filter[h][k];
                }
            }
            //convolve
            unsigned char* index = data + filter_size / 2 * 4 * width + j_offset + filter_size / 2 * 4;
            index[0] = unsigned char(red / denominator);
            index[1] = unsigned char(green / denominator);
            index[2] = unsigned char(blue / denominator);
        }
    }

    return true;
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
void copy(int arr[], int brr[], int n) {
    for (int i = 0; i < n; i++) {
        brr[i] = arr[i];
    }

}
int* create_1D_mask(int n) {
    int * oneD = new int[n];
    int * temp= new int [n];
    for (int i = 0; i < n; i++) {
        oneD[i] = 0;
        temp[i] = 0;
    }
    oneD[0] = 1;
    temp[0] = 1;
    for (int i = 1; i < n; i++) {
        for (int j = 1; j <= i; j++) {
            int res = temp[j - 1] + temp[j];
            oneD[j] = res;
        }
        copy(oneD, temp, n);
    }
    return oneD;
}
int** create_2D_mask(int n) {
    int* oneD = create_1D_mask(n);
    int** twoD = new int* [n];
    for (int i = 0; i < n; i++) {
        twoD[i] = new int[n];
        for (int j = 0; j < n; j++) {
            twoD[i][j] = oneD[i] * oneD[j];
        }
    }
    return twoD;
}

float cal_denominator(int** mask, int n) {
    float sum = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            sum += mask[i][j];
        }
    }
    return sum;
}

bool TargaImage::Filter_Gaussian_N( unsigned int N )
{
    if (!data)
        return false;

    try {
        const int filter_size = N;
        int** gauss_n_mask = create_2D_mask(filter_size);
        float denominator = cal_denominator(gauss_n_mask, filter_size);

        for (int i = 0; i < height - filter_size; i++)
        {
            int i_offset = i * width * 4;
            for (int j = 0; j < width - filter_size; j++)
            {
                int j_offset = i_offset + j * 4;

                int red = 0, green = 0, blue = 0;
                for (int h = 0; h < filter_size; h++) {
                    int h_offset = j_offset + h * width * 4;
                    for (int k = 0; k < filter_size; k++) {
                        unsigned char* channels = data + (h_offset + k * 4);
                        red += channels[0] * gauss_n_mask[h][k];
                        green += channels[1] * gauss_n_mask[h][k];
                        blue += channels[2] * gauss_n_mask[h][k];
                    }
                }
                //convolve
                unsigned char* index = data + filter_size / 2 * 4 * width + j_offset + filter_size / 2 * 4;
                index[0] = unsigned char(red / denominator);
                index[1] = unsigned char(green / denominator);
                index[2] = unsigned char(blue / denominator);
            }
        }
        return true;
    }
    catch (...) {
        cout << "An exception occured while performing Gaussian-N filter on image." << endl;
    }
    return false;
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
    ClearToBlack();
    return false;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
    ClearToBlack();
    return false;
}// Filter_Enhance


///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
void paint_layer(TargaImage* canvas, TargaImage* referenceImage, int brushRadius) {
    for (int i = 0; i < referenceImage->height; i++) {
        for (int j = 0; j < canvas->height; j++) {

        }
    }
}

bool TargaImage::NPR_Paint()
{
    ClearToBlack();
    return false;

}



///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size()
{
    if (!data)
        return false;

    const int filter_size = 3;
    float bartlett_filter[filter_size][filter_size] = { {1.0/16, 1.0/8, 1.0/16}, {1.0/8, 1.0/4, 1.0/8}, {1.0/16, 1.0/8, 1.0/16} };
    unsigned char* half_data = new unsigned char[width/2 * height/2 * 4];

    int g = 0;
    for (int i = 0; i < height/2; i++)
    {
        for (int j = 0; j < width/2; j++)
        {
            int i_offset = i * 2 * width * 4;
            float val_r = 0.0;
            float val_g = 0.0;
            float val_b = 0.0;
            float val_a = 0.0;


            //apply 3x3 bartlett filter
            for (int k = 0; k < 3; k++) {
                int j_2 = j * 2;
                for (int h = 0; h < 3; h++) {
                    
                    unsigned char* channels = data + (i_offset + j_2 * 4);
                    val_r += bartlett_filter[k][h] * channels[0];
                    val_g += bartlett_filter[k][h] * channels[1];
                    val_b += bartlett_filter[k][h] * channels[2];
                    val_a += bartlett_filter[k][h] * channels[3];
                    j_2++;
                }
                i_offset += width * 4;
            }
            half_data[g++] = (unsigned char)val_r;
            half_data[g++] = (unsigned char)val_g;
            half_data[g++] = (unsigned char)val_b;
            half_data[g++] = (unsigned char)val_a;
        }
    }

    delete[] data;
    data = half_data;
    this->width = width / 2;
    this->height = height / 2;
    return true;

}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
    if (!data)
        return false;

    float filter1[3][3] = { {1.0 / 16, 1.0 / 8, 1.0 / 16}, {1.0 / 8, 1.0 / 4, 1.0 / 8}, {1.0 / 16, 1.0 / 8, 1.0 / 16} };
    float filter2[4][4] = { {1.0 / 64, 3.0 / 64, 3.0 / 64, 1.0 / 64}, {3.0 / 64, 9.0 / 64, 9.0 / 64, 3.0 / 64}, {3.0 / 64, 9.0 / 64, 9.0 / 64, 3.0 / 64},{1.0 / 64, 3.0 / 64, 3.0 / 64, 1.0 / 64} };
    float filter3[4][3] = { {1.0 / 32, 2.0 /32, 1.0 / 32}, {3.0 / 32, 6.0 / 32, 3.0 / 32},  {3.0 / 32, 6.0 / 32, 3.0 / 32} ,{1.0 / 32, 2.0 / 32, 1.0 / 32} };
    unsigned char* double_data = new unsigned char[width* 2 * height * 2 * 4];

    int g = 0;
    for (int i = 0; i < height * 2; i++)
    {
        for (int j = 0; j < width * 2; j++)
        {
            float val_r = 0.0;
            float val_g = 0.0;
            float val_b = 0.0;
            float val_a = 0.0;

            //choose a filter

            if (i % 2 == 0 && j % 2 == 0) {
                int i_offset = i/2 * width * 4;
                for (int k = 0; k < 3; k++) {
                    int j_2 = j/2;
                    for (int h = 0; h < 3; h++) {

                        unsigned char* channels = data + (i_offset + j_2 * 4);
                        val_r += filter1[k][h] * channels[0];
                        val_g += filter1[k][h] * channels[1];
                        val_b += filter1[k][h] * channels[2];
                        val_a += filter1[k][h] * channels[3];
                        j_2++;
                    }
                    i_offset += width * 4;
                }
                double_data[g++] = (unsigned char)val_r;
                double_data[g++] = (unsigned char)val_g;
                double_data[g++] = (unsigned char)val_b;
                double_data[g++] = (unsigned char)val_a;
            }
            else if (i % 2 == 1 && j % 2 == 1) {
                int i_2 = abs(i / 2 - 1);
                int k_limit = 4;
                int h_limit = 4;
                int i_offset = i_2 * width * 4;
                for (int k = 0; k < k_limit; k++) {
                    int j_2 = abs(j/2 -1);
                       
                    for (int h = 0; h < h_limit; h++) {

                        unsigned char* channels = data + (i_offset + j_2 * 4);
                        val_r += filter2[k][h] * channels[0];
                        val_g += filter2[k][h] * channels[1];
                        val_b += filter2[k][h] * channels[2];
                        val_a += filter2[k][h] * channels[3];
                        j_2++;
                    }
                    i_offset += width * 4;
                }
                double_data[g++] = (unsigned char)val_r;
                double_data[g++] = (unsigned char)val_g;
                double_data[g++] = (unsigned char)val_b;
                double_data[g++] = (unsigned char)val_a;
            }
            else if (i % 2 == 0 && j % 2 == 1) {
                int i_2 = abs(i / 2 - 1);
                int k_limit = 4;
                int h_limit = 3;
                
                int i_offset = i_2 * width * 4;
                for (int k = 0; k < k_limit; k++) {
                    int j_2 = abs(j/2-1);
                    
                    for (int h = 0; h < h_limit; h++) {

                        unsigned char* channels = data + (i_offset + j_2 * 4);
                        val_r += filter3[k][h] * channels[0];
                        val_g += filter3[k][h] * channels[1];
                        val_b += filter3[k][h] * channels[2];
                        val_a += filter3[k][h] * channels[3];
                        j_2++;
                    }
                    i_offset += width * 4;
                }
                double_data[g++] = (unsigned char)val_r;
                double_data[g++] = (unsigned char)val_g;
                double_data[g++] = (unsigned char)val_b;
                double_data[g++] = (unsigned char)val_a;
            }
            else
            {
                 int i_2 = i / 2 - 1;
                if (i_2 < 0)
                    i_2 = 0;
                int i_offset = i_2 * width * 4;
                for (int k = 0; k < 4; k++) {
                    int j_2 = j / 2 - 1;
                    if (j_2 < 0)
                        j_2 = 0;
                    for (int h = 0; h < 3; h++) {

                        unsigned char* channels = data + (i_offset + j_2 * 4);
                        val_r += filter3[k][h] * channels[0];
                        val_g += filter3[k][h] * channels[1];
                        val_b += filter3[k][h] * channels[2];
                        val_a += filter3[k][h] * channels[3];
                        j_2++;
                    }
                    i_offset += width * 4;
                }
                double_data[g++] = (unsigned char)val_r;
                double_data[g++] = (unsigned char)val_g;
                double_data[g++] = (unsigned char)val_b;
                double_data[g++] = (unsigned char)val_a;
            }
            
        }
    }

    delete[] data;
    data = double_data;
    this->width = width * 2;
    this->height = height * 2;
    return true;
}// Double_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
    ClearToBlack();
    return false;
}// Resize


//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
// 
// 


///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Rotate(float angleDegrees)
{
    ClearToBlack();
    return false;
}// Rotate


//////////////////////////////////////////////////////////////////////////////
//
//      Given a single RGBA pixel return, via the second argument, the RGB
//      equivalent composited with a black background.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb)
{
    const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

    unsigned char  alpha = rgba[3];

    if (alpha == 0)
    {
        rgb[0] = BACKGROUND[0];
        rgb[1] = BACKGROUND[1];
        rgb[2] = BACKGROUND[2];
    }
    else
    {
	    float	alpha_scale = (float)255 / (float)alpha;
	    int	val;
	    int	i;

	    for (i = 0 ; i < 3 ; i++)
	    {
	        val = (int)floor(rgba[i] * alpha_scale);
	        if (val < 0)
		    rgb[i] = 0;
	        else if (val > 255)
		    rgb[i] = 255;
	        else
		    rgb[i] = val;
	    }
    }
}// RGA_To_RGB


///////////////////////////////////////////////////////////////////////////////
//
//      Copy this into a new image, reversing the rows as it goes. A pointer
//  to the new image is returned.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Reverse_Rows(void)
{
    unsigned char   *dest = new unsigned char[width * height * 4];
    TargaImage	    *result;
    int 	        i, j;

    if (! data)
    	return NULL;

    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = (height - i - 1) * width * 4;
	    int out_offset = i * width * 4;

	    for (j = 0 ; j < width ; j++)
        {
	        dest[out_offset + j * 4] = data[in_offset + j * 4];
	        dest[out_offset + j * 4 + 1] = data[in_offset + j * 4 + 1];
	        dest[out_offset + j * 4 + 2] = data[in_offset + j * 4 + 2];
	        dest[out_offset + j * 4 + 3] = data[in_offset + j * 4 + 3];
        }
    }

    result = new TargaImage(width, height, dest);
    delete[] dest;
    return result;
}// Reverse_Rows


///////////////////////////////////////////////////////////////////////////////
//
//      Clear the image to all black.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::ClearToBlack()
{
    memset(data, 0, width * height * 4);
}// ClearToBlack


///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::Paint_Stroke(const Stroke& s) {
   int radius_squared = (int)s.radius * (int)s.radius;
   for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) {
      for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) {
         int x_loc = (int)s.x + x_off;
         int y_loc = (int)s.y + y_off;
         // are we inside the circle, and inside the image?
         if ((x_loc >= 0 && x_loc < width && y_loc >= 0 && y_loc < height)) {
            int dist_squared = x_off * x_off + y_off * y_off;
            if (dist_squared <= radius_squared) {
               data[(y_loc * width + x_loc) * 4 + 0] = s.r;
               data[(y_loc * width + x_loc) * 4 + 1] = s.g;
               data[(y_loc * width + x_loc) * 4 + 2] = s.b;
               data[(y_loc * width + x_loc) * 4 + 3] = s.a;
            } else if (dist_squared == radius_squared + 1) {
               data[(y_loc * width + x_loc) * 4 + 0] = 
                  (data[(y_loc * width + x_loc) * 4 + 0] + s.r) / 2;
               data[(y_loc * width + x_loc) * 4 + 1] = 
                  (data[(y_loc * width + x_loc) * 4 + 1] + s.g) / 2;
               data[(y_loc * width + x_loc) * 4 + 2] = 
                  (data[(y_loc * width + x_loc) * 4 + 2] + s.b) / 2;
               data[(y_loc * width + x_loc) * 4 + 3] = 
                  (data[(y_loc * width + x_loc) * 4 + 3] + s.a) / 2;
            }
         }
      }
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
               unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
   radius(iradius),x(ix),y(iy),r(ir),g(ig),b(ib),a(ia)
{
}

