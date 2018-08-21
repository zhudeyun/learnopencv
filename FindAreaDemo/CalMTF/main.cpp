#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "nullptr_emulation.h"
#include "sfr.h"

using namespace std;
using namespace cv;

static int pixelPeriod = 4;
unsigned char* p = nullptr;
const int roiNum = 4;
float MTFS[4];

int CalSFR(Mat OriginImageX1, int Roi[][4], int H, float sfr[]);
int CalMTF(float sfr[roiNum]);


void main()
{
	float sfrs[roiNum] = {0};
	CalMTF(sfrs);
	getchar();
	return;
}

int CalMTF(float sfr[roiNum])
{
	Mat gray,dst,displayImg;
	Mat src = imread("lensfocus.bmp");
	if (!src.data)
	{
		cout<<"image not exit"<<endl;
		return 1;
	}	
	//init roi
	int rois[roiNum][4];
	for (int i = 0; i < roiNum; i++)
	{
		rois[i][1] = 260;	//y
		rois[i][2] = 80;	//w
		rois[i][3] = 120;	//h
	}
	rois[0][0] = 200;rois[1][0]=350;rois[2][0]=450;rois[3][0]=580;	//x

	cvtColor(src,gray,CV_RGB2GRAY);

	//draw roi
	displayImg = gray.clone();
	for (int i=0;i<roiNum;i++)
	{
		rectangle(displayImg,Point(rois[i][0],rois[i][1]),Point(rois[i][0]+rois[i][2],rois[i][1]+rois[i][3]),Scalar(255));
	}
	imshow("roi display",displayImg);
	waitKey(0);

	CalSFR(gray,rois,0,sfr);
	return 0;
}
int CalSFR(Mat OriginImageX1, int Roi[][4], int H, float sfr[])
{
	float mtfs[3] = { 0, 0, 0 };
	float mtfValue = 0;
	Mat grayX1, workImageX1, cedgeX1;

	grayX1 = OriginImageX1.clone();

	int imageTypeX1 = OriginImageX1.type();

	// if image is 8UC3, need to convert to 8 bit gray scale
	if (imageTypeX1 == CV_8UC3)
		cvtColor(OriginImageX1, grayX1, COLOR_BGR2GRAY);

	// adapt filter on the image
	//blur(gray, workImage, Size(3, 3));
	workImageX1 = grayX1.clone();

	grayX1.release();


	// if image is 16UC1, first downsample to 8 bit
	if (imageTypeX1 == CV_16U)
		workImageX1.convertTo(workImageX1, CV_8U, 0.00390625);
	
	char ResultString[80];

	int i, d_tmpNumROI;
	int nX, nY, nW, nH;
	int Freqs[5], tmpFreq;
	long MTFs[5];
	int NumOfFreqs;
	int PeakLSF;
	int CentralLoc;
	long LineSlope;

	//new.....
	float fCustomizedPitch = 6;	// 2.2;
	float fPeriod;
	int iPeriod;
	int startFreq = 30;
	int endFreq = 70;
	int stepFreq = 20;
	//.....

	PeakLSF = 0;
	CentralLoc = 0;
	LineSlope = 0;

	Freqs[0] = 44;
	NumOfFreqs = 1;

	int imageWidth = workImageX1.cols;
	int imageHeight = workImageX1.rows;
	
	unsigned char* bufferX1 = new unsigned char[imageHeight * imageWidth];	
	for (int i = 0; i < imageHeight; ++i) {
		p = workImageX1.ptr<uchar>(i);
		for (int j = 0; j < imageWidth; ++j) {
			bufferX1[i * imageWidth + j] = p[j];
		}
	}
	
	/*if (p != nullptr)
	{
		delete[] p;
	}*/
	

	workImageX1.release();
	cedgeX1.release();
	
	int nresult;
	
	for (d_tmpNumROI = 0; d_tmpNumROI<roiNum; d_tmpNumROI++) 
	{
		nX = Roi[d_tmpNumROI][0];
		nY = Roi[d_tmpNumROI][1];
		nW = Roi[d_tmpNumROI][2];
		nH = Roi[d_tmpNumROI][3];
		/*
		switch (d_tmpNumROI)
		{
		case 2:
		case 10:
		case 13:
		case 5:
			nresult = ImageMtfLSF_rev57(bufferX1, imageWidth, imageHeight, nX, nY, 40, nH, nW, 10, pixelPeriod, Freqs, MTFs, NumOfFreqs, &PeakLSF, &CentralLoc, &LineSlope);
			break;
		case 1:
		case 16:
		case 9:
		case 14:
		case 6:
			nresult = ImageMtfLSF_rev57(bufferX2, imageWidth, imageHeight, nX, nY, 40, nH, nW, 10, pixelPeriod, Freqs, MTFs, NumOfFreqs, &PeakLSF, &CentralLoc, &LineSlope);
			break;
		case 3:
		case 11:
		case 0:
		case 17:
		case 8:
			nresult = ImageMtfLSF_H_rev57(bufferY1, imageWidth, imageHeight, nY, nX, 40, nW, nH, 10, pixelPeriod, Freqs, MTFs, NumOfFreqs, &PeakLSF, &CentralLoc, &LineSlope);
			break;
		case 12:
		case 4:
		case 15:
		case 7:
			nresult = ImageMtfLSF_H_rev57(bufferY2, imageWidth, imageHeight, nY, nX, 40, nW, nH, 10, pixelPeriod, Freqs, MTFs, NumOfFreqs, &PeakLSF, &CentralLoc, &LineSlope);
			break;;
		default:
			break;
		}*/
		
		nresult = ImageMtfLSF_rev57(bufferX1, imageWidth, imageHeight, nX, nY, 40, nH, nW, 10, pixelPeriod, Freqs, MTFs, NumOfFreqs, &PeakLSF, &CentralLoc, &LineSlope);

		cout << nresult << endl;
		if (nresult != 0)
		{
			return nresult;
		}
			
		for (i = 0; i < NumOfFreqs; i++)
		{
			mtfValue = double(MTFs[i]) / 1024 / 1024;
			printf("MTF: %6.4f, ", mtfValue);
		}
		MTFS[d_tmpNumROI] = mtfValue;
		printf("PeakLSF: %7.3f, ", double(PeakLSF) / 1024);
		printf("CentralLoc: %7.3f, ", double(CentralLoc) / 1024);
		printf("LineSlope: %7.3f, ", double(LineSlope) / 1024 / 1024);
		printf("Res: %d,\r\n ", nresult);
	}
	for (int y = 0; y < roiNum; y++)
	{
		sfr[y] = MTFS[y];
	}
	if (bufferX1 != nullptr)
		delete[] bufferX1;
	return 0;
}