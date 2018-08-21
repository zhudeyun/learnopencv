#include "sfr.h"
#include <math.h>
#include <iostream>
using namespace std;

// macros
static long dmaxarg1, dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ? (dmaxarg1): (dmaxarg2))

long CosT[400];
void CalcCosTables()
{
	int i;
	for (i = 0; i<400; i++){
		CosT[i] = (long)((cos(i*PI / 360)) * 1024);
	}
}

#define LUTsize 512
#define LUTsizeInBits 9
#define LUThalfSize 256

long CosLUT[LUTsize], SinLUT[LUTsize], CosDeltaLUT[LUTsize], SinDeltaLUT[LUTsize];

void IniLUT()
{
	int i;
	for (i = 0; i<LUTsize; i++){
		CosLUT[i] = (long)(cos(i*PI / LUThalfSize) * 1024 * 1024);
		SinLUT[i] = (long)(sin(i*PI / LUThalfSize) * 1024 * 1024);
		CosDeltaLUT[i] = (long)((cos((i + 1)*PI / LUThalfSize) - cos(i*PI / LUThalfSize)) * 1024 * 1024);
		SinDeltaLUT[i] = (long)((sin((i + 1)*PI / LUThalfSize) - sin(i*PI / LUThalfSize)) * 1024 * 1024);
	}
}

void linear_regression_MTF(int num_points, long *y_points, long *b, long *a)
{
	long x_avg = 0;
	long y_avg = 0;
	long  numerator = 0;
	long denominator = 0;
	int i;
	float temp_b, temp_a;

	for (i = 0; i<num_points; i++) {
		x_avg += i;
		y_avg += *(y_points + i);
	}
	x_avg = ((long)(((int)x_avg + num_points / 2) / num_points));
	y_avg = ((long)(((int)y_avg + num_points / 2) / num_points));

	for (i = 0; i<num_points; i++) {
		numerator += (i - x_avg) * ((long)*(y_points + i) - y_avg);
		denominator += (i - x_avg) * (i - x_avg);
	}
	if (denominator == 0) denominator++;
	temp_b = (float)numerator / (float)(denominator);
	temp_a = (float)y_avg - (temp_b * (float)x_avg);

	*b = (long)(temp_b * 1024);
	*a = (long)(temp_a * 1024);
}


//-----------------------------
int ImageMTFLSF(unsigned char*TargetImage, int ImgXDim, int ImgYDim, int RoiStartPtX, int RoiStartPtY, int RoiWidth, int RoiHeight, int SearchRoiWidth, int BackgroundCalcMarginWidth, int PixelPitch, const int *SpatialFreqs, long *MTF, int numOfFreqs, int *LsfPeakValue, int *CentralPtLoc, long *LineSlope)
/*   Calcualte MTF based on a vertically oriented slit image

X means Column, Horizontal and Width
Y means Row, Veritcal and Height

TargetImage: pointer to the current frame (typedef  unsigned char byte).

static int ImgXDim; //ImgXDim: number of pixels in the current frame along X direction. =1280 for DR.
static int ImgYDim; //ImgYDim: number of pixels in the current frame along Y direction. =960 for DR.

RoiPtX: the x-coordinate of the upper-left corner of the rectangle.
Roi stands for Region of interest, which is also a rectangle.
RoiPtY: the y-coordinate of the upper-left corner of the rectangle.
RoiWidth:  the width of the rectangle.
RoiHeight:  the height of the rectangle.
PixelPitch: the pitch of imaging sensor, in nanometer.
SpatialFreqs: pointer to the frequencies at with MTF will be caluclated, in cyc/mm and shifted.
MTF: pointer to the MTF. shifted using number of bits defined by MtfDataBitsShifted
numOfFreqs: number of frequencies of interest.
LsfPeakValue: the peak value of LSF (line spread function), for monitoring purpose.
CentralPtLoc: the x-coordinate of the slit at y=ImgYDim/2, used to estimate the pointing and magnification
note this is not the true center of the slit.

Subroutine Return:
0 - normal execution.
1 - not vertically slanted image or slit not extent throughout the ROI.
2 - background too noisy. SNR too low.
3 - The slope of line too big or too small
4 - zero counts at LSF bin. The slope of line might be too small (aligned too close to y axis).
*/
//rev.5 replace all the float point with long, number of bits shifted is determined by BitsShifted
//rev.51 add calculation for central pt location, number of bits shifted is determined by BitsShifted
//rev.52 larger area is used to calculate the background dark level
//rev. 57 deal with dark images
{
	int ReturnStatus = 0;
	int i, j, CurX, CurY, CurPixelValue;
	long dt, dt1;
	long f_dt, f_dt1, f_CurPixelValue;

	int UpdatedRoiStartPtX; // X value of updated ROI
	long BackgroundLevel; //dark background, be delete when calculate LSF
	long Centroid[MaxNumRows], tmpCentroid;

	long line_slope, offset; // the line curve of centroids

	long tmp_peak_loc, wid, wid_left, wid_right; //parameters for Hamming window

	int numRows, tmp_numRows;

	//LSF alignment
	long binArray1[MaxNumCols*NumBin], binArray2[MaxNumCols*NumBin];
	long lsf[MaxNumCols*NumBin], tmp_lsf[MaxNumCols*NumBin];
	long new_offset, new_centroid;
	int new_start, newNumPixels, new_center, old_center, pixel_shift;
	int ling;
	int numZeroCount = 0;

	//FFT
	long a, b, MTF0;
	long g0, g1;
	double tmp_a, tmp_b;

	//Just need initial once in the programe
	CalcCosTables();
	IniLUT();

	//finetune ROI based on top rows and bottom rows
	tmpCentroid = 0;
	for (j = 0; j<ColMarginRoiSearch; j++){
		//top rows
		dt = 0;
		dt1 = 0;
		CurY = RoiStartPtY + j;
		for (i = -RowMarginRoiSearch; i<SearchRoiWidth + RowMarginRoiSearch; i++) {
			CurX = RoiStartPtX + i;
			CurPixelValue = *(TargetImage + CurY * ImgXDim + CurX);
			if ((CurPixelValue > BackgroundCountThreshold) && (CurPixelValue < BackgroundCountCeiling))
			{
				dt += CurPixelValue* i;
				dt1 += CurPixelValue;
			}
		}
		if (dt1<MinNumOfPixels)
			return 10;
		else
			tmpCentroid += (dt << BitsShifted) / dt1;

		// bot rows
		dt = 0;
		dt1 = 0;
		CurY = RoiStartPtY + RoiHeight - 1 - j;
		for (i = -RowMarginRoiSearch; i<SearchRoiWidth + RowMarginRoiSearch; i++) {
			CurX = RoiStartPtX + i;
			CurPixelValue = *(TargetImage + CurY * ImgXDim + CurX);
			if ((CurPixelValue > BackgroundCountThreshold) && (CurPixelValue < BackgroundCountCeiling))
			{
				dt += CurPixelValue* i;
				dt1 += CurPixelValue;
			}
		}
		if (dt1<MinNumOfPixels)
			return 11;
		else
			tmpCentroid += (dt << BitsShifted) / dt1;
	}
	tmpCentroid = (tmpCentroid / (2 * ColMarginRoiSearch)) >> BitsShifted;
	UpdatedRoiStartPtX = RoiStartPtX + tmpCentroid - (RoiWidth >> 1);

	//calculate background, based on the left and right columns just outside the ROI
	dt = 0;
	dt1 = 0;
	for (i = 0; i<BackgroundCalcMarginWidth; i++) {
		//left bank
		CurX = UpdatedRoiStartPtX + i;
		for (j = 0; j<RoiHeight; j++) {
			CurY = RoiStartPtY + j;
			if (*(TargetImage + CurY * ImgXDim + CurX)<BackgroundCalcCeiling){ //added in Rev 7
				dt += *(TargetImage + CurY * ImgXDim + CurX);
				dt1 += 1;
			}
		}
		//right bank
		CurX = UpdatedRoiStartPtX + RoiWidth - i;
		for (j = 0; j<RoiHeight; j++) {
			CurY = RoiStartPtY + j;
			if (*(TargetImage + CurY * ImgXDim + CurX)<BackgroundCalcCeiling){    //added in Rev 7
				dt += *(TargetImage + CurY * ImgXDim + CurX);
				dt1 += 1;
			}
		}
	}

	if (dt1<MinNumOfPixels)  
		return 12;

	BackgroundLevel = (dt << BitsShifted) / dt1;

	if (BackgroundLevel>(MaxBackgroundLevel << BitsShifted))
		ReturnStatus = 2;

	//locate centroids of each row
	for (j = 0; j<RoiHeight; j++){
		f_dt = 0;
		f_dt1 = 0;
		CurY = RoiStartPtY + j;
		for (i = 0; i<RoiWidth; i++) {
			CurX = UpdatedRoiStartPtX + i;
			f_CurPixelValue = (*(TargetImage + CurY * ImgXDim + CurX) << BitsShifted) - BackgroundLevel;
			f_dt += f_CurPixelValue* i;
			f_dt1 += f_CurPixelValue;
		}
		if (f_dt1 >> BitsShifted == 0)
		{
			Centroid[j] = 0;
		}
		else
		{
			Centroid[j] = f_dt / (f_dt1 >> BitsShifted);
		}
	}

	//fit
	linear_regression_MTF(RoiHeight, Centroid, &line_slope, &offset);

	if ((abs(line_slope)>200 * 1024) || (abs(line_slope)< 5 * 1024))  //line_slope is shifted by BitsShifted
		//return 3;

	// relocate centroids based on Hamming window
	//find the size of minimized
	for (j = 0; j<RoiHeight; j++){
		f_dt = 0;
		f_dt1 = 0;
		CurY = RoiStartPtY + j;
		//parameters for hamming_window
		tmp_peak_loc = (j*line_slope + offset) >> DoubleBitsShifted; // x=a*y+b
		wid_left = tmp_peak_loc;
		wid_right = RoiWidth - 1 - tmp_peak_loc;
		wid = DMAX(wid_left, wid_right);
		for (i = 0; i<RoiWidth; i++) {
			CurX = UpdatedRoiStartPtX + i;
			f_CurPixelValue = (((*(TargetImage + CurY * ImgXDim + CurX) << BitsShifted) - BackgroundLevel)*((566231 + 471 * CosT[abs((i - tmp_peak_loc) * 360 / wid)]) >> BitsShifted)) >> BitsShifted;
			f_dt += f_CurPixelValue* i;
			f_dt1 += f_CurPixelValue;
		}
		if (f_dt1 >> BitsShifted == 0)
		{
			Centroid[j] = 0;
		}
		else
		{
			Centroid[j] = f_dt / (f_dt1 >> BitsShifted);
		}		
	}

	//fit
	linear_regression_MTF(RoiHeight, Centroid, &line_slope, &offset);

	*CentralPtLoc = (UpdatedRoiStartPtX << BitsShifted) + (((RoiHeight >> 1)*line_slope + offset) >> BitsShifted);
	numRows = ((RoiHeight*abs(line_slope)) >> DoubleBitsShifted << DoubleBitsShifted) / abs(line_slope);

	if (numRows<10) numRows = RoiHeight;
	*LineSlope = line_slope;
	//project
	new_offset = (NumBin*(1 - numRows)*line_slope) >> DoubleBitsShifted;
	newNumPixels = RoiWidth*NumBin;

	for (j = 0; j<MaxNumCols*NumBin; j++)
	{
		binArray1[j] = 0;
		binArray2[j] = 0;
	}

	for (j = 0; j < numRows; j++)
	{
		for (i = 0; i < RoiWidth; i++)
		{
			ling = (((i << DoubleBitsShifted) - j*line_slope) >> (DoubleBitsShifted - 2)) + 51; //shift -2 becasue NumBin=4
			binArray1[ling] += 1;
			binArray2[ling] += (*(TargetImage + (RoiStartPtY + j) * ImgXDim + UpdatedRoiStartPtX + i) << BitsShifted) - BackgroundLevel;
		}
	}

	//check zero counts
	new_start = 52 + (new_offset >> 1);
	j = 0;
	f_dt = 0;
	f_dt1 = 0;

	for (i = new_start; i < new_start + newNumPixels; i++)
	{
		if (binArray1[i] == 0)
		{
			numZeroCount++;
			binArray1[i] = (binArray1[i - 1] + binArray1[i + 1]) >> 1;
		}
		if (binArray1[i] == 0)
		{
			lsf[j] = 0;
		}
		else
		{
			lsf[j] = binArray2[i] / binArray1[i];
		}
		
		//calculate centroid
		f_dt += lsf[j] * j;
		f_dt1 += lsf[j];
		j++;
	}

	new_centroid = f_dt / (f_dt1 >> BitsShifted);

	if (numZeroCount>0)
		ReturnStatus = 4;

	//apply hamming window & shift the peak to center
	new_centroid = new_centroid >> BitsShifted;
	wid_left = new_centroid;
	wid_right = newNumPixels - 1 - new_centroid;
	wid = DMAX(wid_left, wid_right);

	new_center = newNumPixels >> 2 - 1;
	old_center = new_centroid;
	pixel_shift = old_center - new_center;

	if (pixel_shift>0)
	{
		pixel_shift = -pixel_shift;
		for (i = 0; i < newNumPixels; i++)
		{
			tmp_lsf[newNumPixels - 1 - i] = lsf[i] * ((566231 + 471 * CosT[abs((i - new_centroid) * 360 / wid)]) >> BitsShifted) >> BitsShifted;
		}
	}
	else
	{
		for (i = 0; i < newNumPixels; i++)
		{
			tmp_lsf[i] = lsf[i] * ((566231 + 471 * CosT[abs((i - new_centroid) * 360 / wid)]) >> BitsShifted) >> BitsShifted;
		}
	}

	for (i = 0; i< newNumPixels + pixel_shift; i++)
	{
		lsf[i - pixel_shift] = tmp_lsf[i];
	}
	for (i = newNumPixels + pixel_shift; i<newNumPixels; i++)
	{
		lsf[newNumPixels - i - 1] = tmp_lsf[i];
	}
	//fft
	for (i = 0, MTF0 = 0, b = 0; i < newNumPixels; i++)
	{
		MTF0 += lsf[i];
		if (lsf[i]> b)	b = lsf[i];
	}
	*LsfPeakValue = b;

	for (j = 0; j <numOfFreqs; j++)
	{
		g0 = ((long)(PixelPitch*SpatialFreqs[j]) << (BitsShifted + LUTsizeInBits)) / (1000 * NumBin);
		for (i = 0, a = 0, b = 0; i < newNumPixels; i++)
		{
			g1 = g0*i;
			g1 = g1 - ((g1 >> (BitsShifted + LUTsizeInBits)) << (BitsShifted + LUTsizeInBits));
			a += lsf[i] * ((CosLUT[g1 >> BitsShifted] + (((g1 - ((g1 >> BitsShifted) << BitsShifted))*CosDeltaLUT[g1 >> BitsShifted]) >> BitsShifted)) >> BitsShifted);
			b += lsf[i] * ((SinLUT[g1 >> BitsShifted] + (((g1 - ((g1 >> BitsShifted) << BitsShifted))*SinDeltaLUT[g1 >> BitsShifted]) >> BitsShifted)) >> BitsShifted);
		}
		tmp_a = ((double)a) / MTF0;
		tmp_b = ((double)b) / MTF0;
		MTF[j] = ((long)sqrt(tmp_a*tmp_a + tmp_b*tmp_b)) << BitsShifted;
	}
	return  ReturnStatus;
}


int SearchEdge(unsigned char*TargetImage,int ImageX, int ImageY, int SearchRoiWidth, int Threshold, int *nX_new, int *temp_lines) {
	// This routine is used to locate the target dynamically
	// Inputs
	//		- P_SCANBUFFER,
	//		- ImageXDimension
	//		- ImageYDimension
	//		- SearchRoiWidth. Width of the search window. Normally 80 or 60.
	//		- Threshold. Threshold to smooth out the background noise. Recommend 10 though 8~15 is OK
	// Outputs
	//		- nX_new[], the X coordinate of the top-left edge of updated search window,
	//					size of the array is limited to 4.
	// Return
	//		- 0 normal
	//		- 1 Check threhold setting

	int nX, i, j;
	int Hori_SearchStart, Hori_SearchEnd;
	int temp[10];
	char ResultString[80];

	for (j = 0; j<4; j++) {
		switch (j) {
		case 0: Hori_SearchStart = 80; Hori_SearchEnd = 300; break;
		case 1: Hori_SearchStart = nX + 30; Hori_SearchEnd = nX + 300; break;
		case 2: Hori_SearchStart = nX + 30; Hori_SearchEnd = nX + 300; break;
		case 3: Hori_SearchStart = nX + 30; Hori_SearchEnd = nX + 300; break;
		}
		if (Hori_SearchEnd >= ImageX)
			Hori_SearchEnd = ImageX - 80;
		//rsh nX=XSearch(P_SCANBUFFER, ImageXDimension, ImageYDimension, Threshold,Hori_SearchStart,Hori_SearchEnd,temp);
		nX = XSearch(TargetImage,ImageX, ImageY, Threshold, Hori_SearchStart, Hori_SearchEnd, temp);
		//sprintf(ResultString, "nX%d: %d, start: %d, end: %d\r\n", j, nX, Hori_SearchStart, Hori_SearchEnd); xputstring(ResultString);
		for (i = 0; i<10; i++) {
			temp_lines[j * 10 + i] = temp[i];
		}

		if ((nX <= Hori_SearchStart) || (nX >= Hori_SearchEnd)) {
			j = 100;
			break;
		}
		else nX_new[j] = nX - SearchRoiWidth / 2;

	}

	if (j>4) {
		//nX_new[0]=ImageX/2 - 160;//480-SearchRoiWidth/2;
		//nX_new[1]=ImageX/2 - 45;//595-SearchRoiWidth/2;
		//nX_new[2]=ImageX/2 + 45;//685-SearchRoiWidth/2;
		//nX_new[3]=ImageX/2 + 160;//800-SearchRoiWidth/2;
		return 1;
	}
	else return 0;
}

int XSearch(unsigned char*TargetImage,int ImgXDim, int ImgYDim, int threshold, int Hori_SearchStart, int Hori_SearchEnd, int *temp_lines) {
	const int Hori_SearchStep = 1;
	//const int Vert_SearchStart=450, Vert_SearchEnd=510, Vert_SearchStep=10;
	const int Vert_SearchStart = ImgYDim / 2 - 30, Vert_SearchEnd = ImgYDim / 2 + 30, Vert_SearchStep = 10;
	const int numOfPt_Discarded = 2;
	int number_of_lines;  // must be > 5
	int i, j;
	int CurRowIndex;
	int CurPixelValue, PrevPixelValue;
	int PixelMinusTwo, PixelMinusThree;
	int nX_lines[20], temp_nX;
	int nX;

	number_of_lines = 0;
	for (j = Vert_SearchStart; j <= Vert_SearchEnd; j = j + Vert_SearchStep){
		nX_lines[number_of_lines] = 0;
		CurPixelValue = 0;
		PrevPixelValue = 0;
		CurRowIndex = j * ImgXDim;
		for (i = Hori_SearchStart; i <= Hori_SearchEnd; i = i + Hori_SearchStep) {
			CurPixelValue = *(TargetImage + CurRowIndex + i);
			if ((PrevPixelValue>threshold) && (CurPixelValue<PrevPixelValue)){
				PixelMinusTwo = *(TargetImage + CurRowIndex + i - 1 - Hori_SearchStep);
				PixelMinusThree = *(TargetImage + CurRowIndex + i - 2 - Hori_SearchStep);
				if ((PixelMinusTwo>threshold) && (PixelMinusThree>threshold)) {
					nX_lines[number_of_lines] = i;
					break;
				}
			}
			PrevPixelValue = CurPixelValue;
		}
		number_of_lines++;
	}

	//Sort and pick the middle ones
	for (i = 0; i<number_of_lines - 1; i++) {
		for (j = i + 1; j<number_of_lines; j++) {
			if (nX_lines[i] < nX_lines[j]) {
				temp_nX = nX_lines[i];  nX_lines[i] = nX_lines[j];  nX_lines[j] = temp_nX;
			}
		}
	}

	for (i = 0; i<10; i++) {
		temp_lines[i] = nX_lines[i];
	}

	nX = 0;
	if (number_of_lines > 2 * numOfPt_Discarded) {
		for (i = numOfPt_Discarded; i<number_of_lines - numOfPt_Discarded; i++) nX = nX + nX_lines[i];
		nX = nX / (number_of_lines - 2 * numOfPt_Discarded) - Hori_SearchStep;
		return nX;
	}
	else
		return 0;
}

//-----------------------------
int ImageMtfLSF_rev57(unsigned char*TargetImage, int ImgXDim, int ImgYDim, int RoiStartPtX, int RoiStartPtY, int RoiWidth, int RoiHeight, int SearchRoiWidth, int BackgroundCalcMarginWidth, int PixelPitch, const int *SpatialFreqs, long *MTF, int numOfFreqs, int *LsfPeakValue, int *CentralPtLoc, long *LineSlope)
/*   Calcualte MTF based on a vertically oriented slit image

X means Column, Horizontal and Width
Y means Row, Veritcal and Height

P_SCANBUFFER: pointer to the current frame (typedef  unsigned char byte).

static int ImgXDim; //ImgXDim: number of pixels in the current frame along X direction. =1280 for DR.
static int ImgYDim; //ImgYDim: number of pixels in the current frame along Y direction. =960 for DR.

RoiPtX: the x-coordinate of the upper-left corner of the rectangle.
Roi stands for Region of interest, which is also a rectangle.
RoiPtY: the y-coordinate of the upper-left corner of the rectangle.
RoiWidth:  the width of the rectangle.
RoiHeight:  the height of the rectangle.
PixelPitch: the pitch of imaging sensor, in nanometer.
SpatialFreqs: pointer to the frequencies at with MTF will be caluclated, in cyc/mm and shifted.
MTF: pointer to the MTF. shifted using number of bits defined by MtfDataBitsShifted
numOfFreqs: number of frequencies of interest.
LsfPeakValue: the peak value of LSF (line spread function), for monitoring purpose.
CentralPtLoc: the x-coordinate of the slit at y=ImgYDim/2, used to estimate the pointing and magnification
note this is not the true center of the slit.

Subroutine Return:
0 - normal execution.
1 - not vertically slanted image or slit not extent throughout the ROI.
2 - background too noisy. SNR too low.
3 - The slope of line too big or too small
4 - zero counts at LSF bin. The slope of line might be too small (aligned too close to y axis).
*/
//rev.5 replace all the float point with long, number of bits shifted is determined by BitsShifted
//rev.51 add calculation for central pt location, number of bits shifted is determined by BitsShifted
//rev.52 larger area is used to calculate the background dark level
//rev. 57 deal with dark images
{
	int ReturnStatus = 0;
	int i, j, CurX, CurY, CurPixelValue;
	long dt, dt1;
	long f_dt, f_dt1, f_CurPixelValue;

	int UpdatedRoiStartPtX; // X value of updated ROI
	long BackgroundLevel; //dark background, be delete when calculate LSF
	long Centroid[MaxNumRows], tmpCentroid;

	long line_slope, offset; // the line curve of centroids

	long tmp_peak_loc, wid, wid_left, wid_right; //parameters for Hamming window

	int numRows, tmp_numRows;

	//LSF alignment
	long binArray1[MaxNumCols*NumBin], binArray2[MaxNumCols*NumBin];
	long lsf[MaxNumCols*NumBin], tmp_lsf[MaxNumCols*NumBin];
	long new_offset, new_centroid;
	int new_start, newNumPixels, new_center, old_center, pixel_shift;
	int ling;
	int numZeroCount = 0;

	//FFT
	long a, b, MTF0;
	long g0, g1;
	double tmp_a, tmp_b;

	//Just need initial once in the programe
	CalcCosTables();
	IniLUT();

	//finetune ROI based on top rows and bottom rows
	tmpCentroid = 0;
	for (j = 0; j<ColMarginRoiSearch; j++){
		//top rows
		dt = 0;
		dt1 = 0;
		CurY = RoiStartPtY + j;
		for (i = -RowMarginRoiSearch; i<SearchRoiWidth + RowMarginRoiSearch; i++) {
			CurX = RoiStartPtX + i;
			CurPixelValue = *(TargetImage + CurY * ImgXDim + CurX);
			if (CurPixelValue<0)
				return 21;
			if ((CurPixelValue > BackgroundCountThreshold) && (CurPixelValue < BackgroundCountCeiling))
			{
				dt += CurPixelValue* i;
				dt1 += CurPixelValue;
			}
		}
		if (dt1<MinNumOfPixels)
			return 1;
		else
		{
			if (dt1 == 0)
				return 5;
			tmpCentroid += (dt << BitsShifted) / dt1;
		}
			

		// bot rows
		dt = 0;
		dt1 = 0;
		CurY = RoiStartPtY + RoiHeight - 1 - j;
		for (i = -RowMarginRoiSearch; i<SearchRoiWidth + RowMarginRoiSearch; i++) {
			CurX = RoiStartPtX + i;
			CurPixelValue = *(TargetImage + CurY * ImgXDim + CurX);
			if ((CurPixelValue > BackgroundCountThreshold) && (CurPixelValue < BackgroundCountCeiling))
			{
				dt += CurPixelValue* i;
				dt1 += CurPixelValue;
			}
		}
		if (dt1<MinNumOfPixels)
			return 1;
		else
		{
			if (dt1 == 0)
				return 6;
			tmpCentroid += (dt << BitsShifted) / dt1;
		}
			
	}
	tmpCentroid = (tmpCentroid / (2 * ColMarginRoiSearch)) >> BitsShifted;
	UpdatedRoiStartPtX = RoiStartPtX + tmpCentroid - (RoiWidth >> 1);

	//calculate background, based on the left and right columns just outside the ROI
	dt = 0;
	dt1 = 0;
	for (i = 0; i<BackgroundCalcMarginWidth; i++) {
		//left bank
		CurX = UpdatedRoiStartPtX + i;
		for (j = 0; j<RoiHeight; j++) {
			CurY = RoiStartPtY + j;
			if (*(TargetImage + CurY * ImgXDim + CurX)<BackgroundCalcCeiling){ //added in Rev 7
				dt += *(TargetImage + CurY * ImgXDim + CurX);
				dt1 += 1;
			}
		}
		//right bank
		CurX = UpdatedRoiStartPtX + RoiWidth - i;
		for (j = 0; j<RoiHeight; j++) {
			CurY = RoiStartPtY + j;
			if (*(TargetImage + CurY * ImgXDim + CurX)<BackgroundCalcCeiling){    //added in Rev 7
				dt += *(TargetImage + CurY * ImgXDim + CurX);
				dt1 += 1;
			}
		}
	}

	if (dt1<MinNumOfPixels)  return 1;
	if (dt1 == 0)
		return 10;
	BackgroundLevel = (dt << BitsShifted) / dt1;

	if (BackgroundLevel>(MaxBackgroundLevel << BitsShifted))
		ReturnStatus = 2;

	//locate centroids of each row
	for (j = 0; j<RoiHeight; j++){
		f_dt = 0;
		f_dt1 = 0;
		CurY = RoiStartPtY + j;
		for (i = 0; i<RoiWidth; i++) {
			CurX = UpdatedRoiStartPtX + i;
			f_CurPixelValue = (*(TargetImage + CurY * ImgXDim + CurX) << BitsShifted) - BackgroundLevel;
			f_dt += f_CurPixelValue* i;
			f_dt1 += f_CurPixelValue;
		}
		if ((f_dt1 >> BitsShifted) == 0)
			return 11;
		Centroid[j] = f_dt / (f_dt1 >> BitsShifted);
	}

	//fit
	linear_regression_MTF(RoiHeight, Centroid, &line_slope, &offset);

	if ((abs(line_slope)>200 * 1024) || (abs(line_slope)< 5 * 1024))  //line_slope is shifted by BitsShifted
		return 3;

	// relocate centroids based on Hamming window
	//find the size of minimized
	for (j = 0; j<RoiHeight; j++){
		f_dt = 0;
		f_dt1 = 0;
		CurY = RoiStartPtY + j;
		//parameters for hamming_window
		tmp_peak_loc = (j*line_slope + offset) >> DoubleBitsShifted; // x=a*y+b
		wid_left = tmp_peak_loc;
		wid_right = RoiWidth - 1 - tmp_peak_loc;
		wid = DMAX(wid_left, wid_right);
		for (i = 0; i<RoiWidth; i++) {
			CurX = UpdatedRoiStartPtX + i;
			f_CurPixelValue = (((*(TargetImage + CurY * ImgXDim + CurX) << BitsShifted) - BackgroundLevel)*((566231 + 471 * CosT[abs((i - tmp_peak_loc) * 360 / wid)]) >> BitsShifted)) >> BitsShifted;
			f_dt += f_CurPixelValue* i;
			f_dt1 += f_CurPixelValue;
		}
		if ((f_dt1 >> BitsShifted) == 0)
			return 12;
		Centroid[j] = f_dt / (f_dt1 >> BitsShifted);
	}

	//fit
	linear_regression_MTF(RoiHeight, Centroid, &line_slope, &offset);

	*CentralPtLoc = (UpdatedRoiStartPtX << BitsShifted) + (((RoiHeight >> 1)*line_slope + offset) >> BitsShifted);
	if (abs(line_slope) == 0)
		return 13;
	numRows = ((RoiHeight*abs(line_slope)) >> DoubleBitsShifted << DoubleBitsShifted) / abs(line_slope);

	if (numRows<10) numRows = RoiHeight;
	*LineSlope = line_slope;
	//project
	new_offset = (NumBin*(1 - numRows)*line_slope) >> DoubleBitsShifted;
	newNumPixels = RoiWidth*NumBin;

	for (j = 0; j<MaxNumCols*NumBin; j++)
	{
		binArray1[j] = 0;
		binArray2[j] = 0;
	}

	for (j = 0; j < numRows; j++)
	{
		for (i = 0; i < RoiWidth; i++)
		{
			ling = (((i << DoubleBitsShifted) - j*line_slope) >> (DoubleBitsShifted - 2)) + 51; //shift -2 becasue NumBin=4
			binArray1[ling] += 1;
			binArray2[ling] += (*(TargetImage + (RoiStartPtY + j) * ImgXDim + UpdatedRoiStartPtX + i) << BitsShifted) - BackgroundLevel;
		}
	}

	//check zero counts
	new_start = 52 + (new_offset >> 1);
	j = 0;
	f_dt = 0;
	f_dt1 = 0;

	for (i = new_start; i < new_start + newNumPixels; i++)
	{
		if (binArray1[i] == 0)
		{
			numZeroCount++;
			binArray1[i] = (binArray1[i - 1] + binArray1[i + 1]) >> 1;
		}
		if (binArray1[i] == 0)
			return 14;
		lsf[j] = binArray2[i] / binArray1[i];
		//calculate centroid
		f_dt += lsf[j] * j;
		f_dt1 += lsf[j];
		j++;
	}
	if ((f_dt1 >> BitsShifted) == 0)
		return 15;
	new_centroid = f_dt / (f_dt1 >> BitsShifted);

	if (numZeroCount>0)
		ReturnStatus = 4;

	//apply hamming window & shift the peak to center
	new_centroid = new_centroid >> BitsShifted;
	wid_left = new_centroid;
	wid_right = newNumPixels - 1 - new_centroid;
	wid = DMAX(wid_left, wid_right);

	new_center = newNumPixels >> 2 - 1;
	old_center = new_centroid;
	pixel_shift = old_center - new_center;

	if (pixel_shift>0)
	{
		pixel_shift = -pixel_shift;
		for (i = 0; i < newNumPixels; i++)
		{
			if (wid == 0)
				return 16;
			tmp_lsf[newNumPixels - 1 - i] = lsf[i] * ((566231 + 471 * CosT[abs((i - new_centroid) * 360 / wid)]) >> BitsShifted) >> BitsShifted;
		}
	}
	else
	{
		for (i = 0; i < newNumPixels; i++)
		{
			if (wid == 0)
				return 17;
			tmp_lsf[i] = lsf[i] * ((566231 + 471 * CosT[abs((i - new_centroid) * 360 / wid)]) >> BitsShifted) >> BitsShifted;
		}
	}

	for (i = 0; i< newNumPixels + pixel_shift; i++)
	{
		lsf[i - pixel_shift] = tmp_lsf[i];
	}
	for (i = newNumPixels + pixel_shift; i<newNumPixels; i++)
	{
		lsf[newNumPixels - i - 1] = tmp_lsf[i];
	}
	//fft
	for (i = 0, MTF0 = 0, b = 0; i < newNumPixels; i++)
	{
		MTF0 += lsf[i];
		if (lsf[i]> b)	b = lsf[i];
	}
	*LsfPeakValue = b;

	for (j = 0; j <numOfFreqs; j++)
	{
		g0 = ((long)(PixelPitch*SpatialFreqs[j]) << (BitsShifted + LUTsizeInBits)) / (1000 * NumBin);
		for (i = 0, a = 0, b = 0; i < newNumPixels; i++)
		{
			g1 = g0*i;
			g1 = g1 - ((g1 >> (BitsShifted + LUTsizeInBits)) << (BitsShifted + LUTsizeInBits));
			a += lsf[i] * ((CosLUT[g1 >> BitsShifted] + (((g1 - ((g1 >> BitsShifted) << BitsShifted))*CosDeltaLUT[g1 >> BitsShifted]) >> BitsShifted)) >> BitsShifted);
			b += lsf[i] * ((SinLUT[g1 >> BitsShifted] + (((g1 - ((g1 >> BitsShifted) << BitsShifted))*SinDeltaLUT[g1 >> BitsShifted]) >> BitsShifted)) >> BitsShifted);
		}
		if (MTF0 == 0) return 18;
		tmp_a = ((double)a) / MTF0;
		tmp_b = ((double)b) / MTF0;
		MTF[j] = ((long)sqrt(tmp_a*tmp_a + tmp_b*tmp_b)) << BitsShifted;
	}
	return  ReturnStatus;
 }

int ImageMtfLSF_H_rev57(unsigned char*TargetImage, int ImgXDim, int ImgYDim, int RoiStartPtX, int RoiStartPtY, int RoiWidth, int RoiHeight, int SearchRoiWidth, int BackgroundCalcMarginWidth, int PixelPitch, const int *SpatialFreqs, long *MTF, int numOfFreqs, int *LsfPeakValue, int *CentralPtLoc, long *LineSlope)
/*   Calcualte MTF based on a Horizontally oriented slit image

X means Column, Horizontal and Width
Y means Row, Veritcal and Height

P_SCANBUFFER: pointer to the current frame (typedef  unsigned char byte).

static int ImgXDim; //ImgXDim: number of pixels in the current frame along X direction. =1280 for DR.
static int ImgYDim; //ImgYDim: number of pixels in the current frame along Y direction. =960 for DR.

RoiPtX: the x-coordinate of the upper-left corner of the rectangle.
Roi stands for Region of interest, which is also a rectangle.
RoiPtY: the y-coordinate of the upper-left corner of the rectangle.
RoiWidth:  the width of the rectangle.
RoiHeight:  the height of the rectangle.
PixelPitch: the pitch of imaging sensor, in nanometer.
SpatialFreqs: pointer to the frequencies at with MTF will be caluclated, in cyc/mm and shifted.
MTF: pointer to the MTF. shifted using number of bits defined by MtfDataBitsShifted
numOfFreqs: number of frequencies of interest.
LsfPeakValue: the peak value of LSF (line spread function), for monitoring purpose.
CentralPtLoc: the x-coordinate of the slit at y=ImgYDim/2, used to estimate the pointing and magnification
note this is not the true center of the slit.

Subroutine Return:
0 - normal execution.
1 - not vertically slanted image or slit not extent throughout the ROI.
2 - background too noisy. SNR too low.
3 - The slope of line too big or too small
4 - zero counts at LSF bin. The slope of line might be too small (aligned too close to y axis).
*/
//rev.5 replace all the float point with long, number of bits shifted is determined by BitsShifted
//rev.51 add calculation for central pt location, number of bits shifted is determined by BitsShifted
{
	int ReturnStatus = 0;
	int i, j, CurX, CurY, CurPixelValue;
	long dt, dt1;
	long f_dt, f_dt1, f_CurPixelValue;

	int UpdatedRoiStartPtX; // X value of updated ROI
	long BackgroundLevel; //dark background, be delete when calculate LSF
	long Centroid[MaxNumRows], tmpCentroid;

	long line_slope, offset; // the line curve of centroids

	long tmp_peak_loc, wid, wid_left, wid_right; //parameters for Hamming window

	int numRows, tmp_numRows;

	//LSF alignment
	long binArray1[MaxNumCols*NumBin], binArray2[MaxNumCols*NumBin];
	long lsf[MaxNumCols*NumBin], tmp_lsf[MaxNumCols*NumBin];
	long new_offset, new_centroid;
	int new_start, newNumPixels, new_center, old_center, pixel_shift;
	int ling;
	int numZeroCount = 0;

	//FFT
	long a, b, MTF0;
	long g0, g1;
	double tmp_a, tmp_b;

	//Just need initial once in the programe
	CalcCosTables();
	IniLUT();

	//finetune ROI based on top rows and bottom rows
	tmpCentroid = 0;
	for (j = 0; j<ColMarginRoiSearch; j++){
		//top rows
		dt = 0;
		dt1 = 0;
		CurY = RoiStartPtY + j;
		for (i = -RowMarginRoiSearch; i<SearchRoiWidth + RowMarginRoiSearch; i++) {
			CurX = RoiStartPtX + i;
			if (CurX<0)
				return 20;
			CurPixelValue = *(TargetImage + CurX * ImgXDim + CurY); //Horizontal slit
			if ((CurPixelValue > BackgroundCountThreshold) && (CurPixelValue < BackgroundCountCeiling))
			{
				dt += CurPixelValue* i;
				dt1 += CurPixelValue;
			}
		}
		if (dt1<MinNumOfPixels)
			return 1;
		else
		{
			if (dt1 == 0)
				return 7;
			tmpCentroid += (dt << BitsShifted) / dt1;
		}
			

		// bot rows
		dt = 0;
		dt1 = 0;
		CurY = RoiStartPtY + RoiHeight - 1 - j;
		for (i = -RowMarginRoiSearch; i<SearchRoiWidth + RowMarginRoiSearch; i++) {
			CurX = RoiStartPtX + i;
			CurPixelValue = *(TargetImage + CurX * ImgXDim + CurY); //Horizontal slit
			if ((CurPixelValue > BackgroundCountThreshold) && (CurPixelValue < BackgroundCountCeiling))
			{
				dt += CurPixelValue* i;
				dt1 += CurPixelValue;
			}
		}
		if (dt1<MinNumOfPixels)
			return 1;
		else
		{
			if (dt1 == 0)
				return 8;
			tmpCentroid += (dt << BitsShifted) / dt1;
		}
			
	}
	tmpCentroid = (tmpCentroid / (2 * ColMarginRoiSearch)) >> BitsShifted;
	UpdatedRoiStartPtX = RoiStartPtX + tmpCentroid - (RoiWidth >> 1);

	//calculate background, based on the left and right columns just outside the ROI
	dt = 0;
	dt1 = 0;

	for (i = 0; i<BackgroundCalcMarginWidth; i++) {
		//left bank
		CurX = UpdatedRoiStartPtX + i;
		for (j = 0; j<RoiHeight; j++) {
			CurY = RoiStartPtY + j;
			if (*(TargetImage + CurX * ImgXDim + CurY)<BackgroundCalcCeiling){ //added in Rev 7
				dt += *(TargetImage + CurX * ImgXDim + CurY);
				dt1 += 1;
			}
		}
		//right bank
		CurX = UpdatedRoiStartPtX + RoiWidth - i;
		for (j = 0; j<RoiHeight; j++) {
			CurY = RoiStartPtY + j;
			if (*(TargetImage + CurX * ImgXDim + CurY)<BackgroundCalcCeiling){    //added in Rev 7
				dt += *(TargetImage + CurX * ImgXDim + CurY);
				dt1 += 1;
			}
		}
	}

	if (dt1<MinNumOfPixels)  return 1;
	if (dt1 == 0) return 19;
	BackgroundLevel = (dt << BitsShifted) / dt1;

	if (BackgroundLevel>(MaxBackgroundLevel << BitsShifted))
		ReturnStatus = 2;

	//locate centroids of each row
	for (j = 0; j<RoiHeight; j++){
		f_dt = 0;
		f_dt1 = 0;
		CurY = RoiStartPtY + j;
		for (i = 0; i<RoiWidth; i++) {
			CurX = UpdatedRoiStartPtX + i;
			f_CurPixelValue = (*(TargetImage + CurX * ImgXDim + CurY) << BitsShifted) - BackgroundLevel; //Horizontal slit
			f_dt += f_CurPixelValue* i;
			f_dt1 += f_CurPixelValue;
		}
		if ((f_dt1 >> BitsShifted) == 0) return 20;
		Centroid[j] = f_dt / (f_dt1 >> BitsShifted);
	}

	//fit
	linear_regression_MTF(RoiHeight, Centroid, &line_slope, &offset);

	if ((abs(line_slope)>200 * 1024) || (abs(line_slope)< 5 * 1024))  //line_slope is shifted by BitsShifted
		return 3;

	// relocate centroids based on Hamming window
	//find the size of minimized
	for (j = 0; j<RoiHeight; j++){
		f_dt = 0;
		f_dt1 = 0;
		CurY = RoiStartPtY + j;
		//parameters for hamming_window
		tmp_peak_loc = (j*line_slope + offset) >> DoubleBitsShifted; // x=a*y+b
		wid_left = tmp_peak_loc;
		wid_right = RoiWidth - 1 - tmp_peak_loc;
		wid = DMAX(wid_left, wid_right);
		for (i = 0; i<RoiWidth; i++) {
			CurX = UpdatedRoiStartPtX + i;
			f_CurPixelValue = (((*(TargetImage + CurX * ImgXDim + CurY) << BitsShifted) - BackgroundLevel)*((566231 + 471 * CosT[abs((i - tmp_peak_loc) * 360 / wid)]) >> BitsShifted)) >> BitsShifted; //Horizontal slit
			f_dt += f_CurPixelValue* i;
			f_dt1 += f_CurPixelValue;
		}
		if ((f_dt1 >> BitsShifted) == 0) return 21;
		Centroid[j] = f_dt / (f_dt1 >> BitsShifted);
	}

	//fit
	linear_regression_MTF(RoiHeight, Centroid, &line_slope, &offset);

	*CentralPtLoc = (UpdatedRoiStartPtX << BitsShifted) + (((RoiHeight >> 1)*line_slope + offset) >> BitsShifted);
	if (abs(line_slope) == 0) return 22;
	numRows = ((RoiHeight*abs(line_slope)) >> DoubleBitsShifted << DoubleBitsShifted) / abs(line_slope);

	if (numRows<10) numRows = RoiHeight;
	*LineSlope = line_slope;

	//project
	new_offset = (NumBin*(1 - numRows)*line_slope) >> DoubleBitsShifted;
	newNumPixels = RoiWidth*NumBin;

	for (j = 0; j<MaxNumCols*NumBin; j++)
	{
		binArray1[j] = 0;
		binArray2[j] = 0;
	}

	for (j = 0; j < numRows; j++)
	{
		for (i = 0; i < RoiWidth; i++)
		{
			ling = (((i << DoubleBitsShifted) - j*line_slope) >> (DoubleBitsShifted - 2)) + 51; //shift -2 becasue NumBin=4
			binArray1[ling] += 1;
			binArray2[ling] += (*(TargetImage + (UpdatedRoiStartPtX + i) * ImgXDim + RoiStartPtY + j) << BitsShifted) - BackgroundLevel; //Horizontal slit
		}
	}

	//check zero counts
	new_start = 52 + (new_offset >> 1);
	j = 0;
	f_dt = 0;
	f_dt1 = 0;

	for (i = new_start; i < new_start + newNumPixels; i++)
	{
		if (binArray1[i] == 0)
		{
			numZeroCount++;
			binArray1[i] = (binArray1[i - 1] + binArray1[i + 1]) >> 1;
		}
		if (binArray1[i] == 0) return 23;
		lsf[j] = binArray2[i] / binArray1[i];
		//calculate centroid
		f_dt += lsf[j] * j;
		f_dt1 += lsf[j];
		j++;
	}
	if ((f_dt1 >> BitsShifted) == 0) 
		return 24;
	new_centroid = f_dt / (f_dt1 >> BitsShifted);

	if (numZeroCount>0)
		ReturnStatus = 4;

	//apply hamming window & shift the peak to center
	new_centroid = new_centroid >> BitsShifted;
	wid_left = new_centroid;
	wid_right = newNumPixels - 1 - new_centroid;
	wid = DMAX(wid_left, wid_right);

	new_center = newNumPixels >> 2 - 1;
	old_center = new_centroid;
	pixel_shift = old_center - new_center;

	if (pixel_shift>0)
	{
		pixel_shift = -pixel_shift;
		for (i = 0; i < newNumPixels; i++)
		{
			if (wid == 0) return 25;
			tmp_lsf[newNumPixels - 1 - i] = lsf[i] * ((566231 + 471 * CosT[abs((i - new_centroid) * 360 / wid)]) >> BitsShifted) >> BitsShifted;
		}
	}
	else
	{
		for (i = 0; i < newNumPixels; i++)
		{
			if (wid == 0) return 26;
			tmp_lsf[i] = lsf[i] * ((566231 + 471 * CosT[abs((i - new_centroid) * 360 / wid)]) >> BitsShifted) >> BitsShifted;
		}
	}

	for (i = 0; i< newNumPixels + pixel_shift; i++)
	{
		lsf[i - pixel_shift] = tmp_lsf[i];
	}
	for (i = newNumPixels + pixel_shift; i<newNumPixels; i++)
	{
		lsf[newNumPixels - i - 1] = tmp_lsf[i];
	}
	//fft
	for (i = 0, MTF0 = 0, b = 0; i < newNumPixels; i++)
	{
		MTF0 += lsf[i];
		if (lsf[i]> b)	b = lsf[i];
	}
	*LsfPeakValue = b;

	for (j = 0; j <numOfFreqs; j++)
	{
		g0 = ((long)(PixelPitch*SpatialFreqs[j]) << (BitsShifted + LUTsizeInBits)) / (1000 * NumBin);
		for (i = 0, a = 0, b = 0; i < newNumPixels; i++)
		{
			g1 = g0*i;
			g1 = g1 - ((g1 >> (BitsShifted + LUTsizeInBits)) << (BitsShifted + LUTsizeInBits));
			a += lsf[i] * ((CosLUT[g1 >> BitsShifted] + (((g1 - ((g1 >> BitsShifted) << BitsShifted))*CosDeltaLUT[g1 >> BitsShifted]) >> BitsShifted)) >> BitsShifted);
			b += lsf[i] * ((SinLUT[g1 >> BitsShifted] + (((g1 - ((g1 >> BitsShifted) << BitsShifted))*SinDeltaLUT[g1 >> BitsShifted]) >> BitsShifted)) >> BitsShifted);
		}
		if (MTF0 == 0) return 27;
		tmp_a = ((double)a) / MTF0;
		tmp_b = ((double)b) / MTF0;
		MTF[j] = ((long)sqrt(tmp_a*tmp_a + tmp_b*tmp_b)) << BitsShifted;
	}
	return  ReturnStatus;
}