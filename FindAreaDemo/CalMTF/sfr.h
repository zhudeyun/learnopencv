#define PI 3.14159265359
#define MaxNumRows 150 //Region of imterest should be smaller than 150x150
#define MaxNumCols 150
#define NumBin 4  //Number of Bins used for averaging, do not change unless you know what you are doing
#define RowMarginRoiSearch 5	//10  // Define expanded area used for ROI finetune
#define ColMarginRoiSearch 3
#define BackgroundCountThreshold 5
#define BackgroundCountCeiling 250
#define BackgroundCalcCeiling 6

#define MinNumOfPixels 3  //To check image quality
#define MaxBackgroundLevel 5

#define BitsShifted 10	//number of bits shifted when simulating floating point
#define DoubleBitsShifted 20	//number of bits shifted when simulating floating point
#define MtfDataBitsShifted 20	//number of bits shifted when simulating floating point

int ImageMTFLSF(unsigned char*TargetImage, int ImgXDim, int ImgYDim, int RoiStartPtX, int RoiStartPtY, int RoiWidth, int RoiHeight, int SearchRoiWidth, int BackgroundCalcMarginWidth, int PixelPitch, const int *SpatialFreqs, long *MTF, int numOfFreqs, int *LsfPeakValue, int *CentralPtLoc, long *LineSlope);
int ImageMtfLSF_rev57(unsigned char*TargetImage, int ImgXDim, int ImgYDim, int RoiStartPtX, int RoiStartPtY, int RoiWidth, int RoiHeight, int SearchRoiWidth, int BackgroundCalcMarginWidth, int PixelPitch, const int *SpatialFreqs, long *MTF, int numOfFreqs, int *LsfPeakValue, int *CentralPtLoc, long *LineSlope);
int ImageMtfLSF_H_rev57(unsigned char*TargetImage, int ImgXDim, int ImgYDim, int RoiStartPtX, int RoiStartPtY, int RoiWidth, int RoiHeight, int SearchRoiWidth, int BackgroundCalcMarginWidth, int PixelPitch, const int *SpatialFreqs, long *MTF, int numOfFreqs, int *LsfPeakValue, int *CentralPtLoc, long *LineSlope);
int XSearch(unsigned char*TargetImage, int ImgXDim, int ImgYDim, int threshold, int Hori_SearchStart, int Hori_SearchEnd, int *temp_lines);
int SearchEdge(unsigned char*TargetImage, int ImageX, int ImageY, int SearchRoiWidth, int Threshold, int *nX_new, int *temp_lines);