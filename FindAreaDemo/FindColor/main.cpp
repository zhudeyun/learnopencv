#include <stdio.h>
#include <iostream>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace std;
using namespace cv;
char* RedWindow = "Red";
char* GreenWindow = "Green";
char* YellowWindow = "Yellow";
char* WhiteWindow = "White";
//Red
int iLowRedH = 0;	//156;
int iHighRedH = 10;	//180;
int iLowRedS = 43; 
int iHighRedS = 255;
int iLowRedV = 46;
int iHighRedV = 255;
//Green
int iLowGreenH = 35;
int iHighGreenH = 77;
int iLowGreenS = 43; 
int iHighGreenS = 255;
int iLowGreenV = 46;
int iHighGreenV = 255;
//Yellow->Orange
int iLowYellowH = 11;	//26;
int iHighYellowH = 25;	//34;
int iLowYellowS = 43; 
int iHighYellowS = 255;
int iLowYellowV = 46;
int iHighYellowV = 255;
//White
int iLowWhiteH = 0;
int iHighWhiteH =180;
int iLowWhiteS = 0; 
int iHighWhiteS = 30;
int iLowWhiteV = 221;
int iHighWhiteV = 255;


Mat srcRed,hsvImgRed;
Mat srcGreen,hsvImgGreen;
Mat srcYellow,hsvImgYellow;
Mat srcWhite,hsvImgWhite;

enum COLOR_TYPE
{
	RED,
	GREEN,
	YELLOW,
	WHITE
};
COLOR_TYPE ct;
static void thresh_callback_red(int, void*);
static void thresh_callback_green(int, void*);
static void thresh_callback_yellow(int, void*);
static void thresh_callback_white(int, void*);
vector<MatND> getHSVHist(Mat &src);
void filteredRed(const Mat &inputImage, Mat &resultGray, Mat &resultColor);
void filteredColor(int colorType,const Mat &inputImage, Mat &resultGray, Mat &resultColor);
void TestRed();
void TestGreen();
void TestYellow();
void TestWhite();
void TestCalibration();
void TestColor(COLOR_TYPE color, char* filename);
void main()
{
	//TestCalibration();
	TestColor(RED,"red.bmp");
	waitKey(0);
	TestColor(GREEN,"green.bmp");
	waitKey(0);
	TestColor(YELLOW,"yellow.bmp");
	waitKey(0);
	TestColor(WHITE,"white.bmp");
	waitKey(0);
	
	return;
}
void TestCalibration()
{
	Mat imgCal = imread("red_standard.png",CV_LOAD_IMAGE_COLOR);
	getHSVHist(imgCal);
	getchar();
}
void TestRed()
{
	//red
	srcRed = imread("red.bmp",CV_LOAD_IMAGE_COLOR);
	namedWindow("red source",CV_WINDOW_NORMAL);
	imshow("red source",srcRed);
	//blur(srcRed,srcRed,Size(3,3)); //高斯模糊 
	//imshow("blured image",srcRed);
	//waitKey(0);
	cvtColor(srcRed,hsvImgRed,COLOR_BGR2HSV);
	//直方图
	/*vector<Mat> hsvSplitRed; 
	split(hsvImgRed, hsvSplitRed); 
	imshow("split 0",hsvSplitRed[0]);
	imshow("split 1",hsvSplitRed[1]);
	imshow("split 2",hsvSplitRed[2]);
	waitKey(0);
	equalizeHist(hsvSplitRed[2],hsvSplitRed[2]);  
	imshow("equalizeHist",hsvSplitRed[2]);

	merge(hsvSplitRed,hsvImgRed);  
	imshow("merged img",hsvImgRed);
	waitKey(0);*/

	namedWindow(RedWindow, CV_WINDOW_NORMAL);
	createTrackbar("iLowRedH：", RedWindow,&iLowRedH,180,thresh_callback_red );
	createTrackbar("iHighRedH：", RedWindow,&iHighRedH,180,thresh_callback_red );
	//createTrackbar("iLowRedS：", RedWindow,&iLowRedS,255,thresh_callback_red );
	//createTrackbar("iHighRedS：", RedWindow,&iHighRedS,255,thresh_callback_red );
	//createTrackbar("iLowRedV：", RedWindow,&iLowRedV,255,thresh_callback_red );
	//createTrackbar("iHighRedV：", RedWindow,&iHighRedV,255,thresh_callback_red );
	thresh_callback_red(0,0);

}
void TestGreen()
{
	srcGreen = imread("green.bmp",CV_LOAD_IMAGE_COLOR);
	namedWindow("green source",CV_WINDOW_NORMAL);
	imshow("green source",srcGreen);
	cvtColor(srcGreen,hsvImgGreen,COLOR_BGR2HSV);

	//直方图
	vector<Mat> hsvSplitGreen; 
	split(hsvImgGreen, hsvSplitGreen); 
	equalizeHist(hsvSplitGreen[2],hsvSplitGreen[2]);  
	merge(hsvSplitGreen,hsvImgGreen); 

	namedWindow(GreenWindow, CV_WINDOW_NORMAL);
	createTrackbar("iLowGreenH：", GreenWindow,&iLowGreenH,180,thresh_callback_green );
	createTrackbar("iHighGreenH：", GreenWindow,&iHighGreenH,180,thresh_callback_green );
	//createTrackbar("iLowGreenS：", GreenWindow,&iLowGreenS,255,thresh_callback_green );
	//createTrackbar("iHighGreenS：", GreenWindow,&iHighGreenS,255,thresh_callback_green );
	//createTrackbar("iLowGreenV：", GreenWindow,&iLowGreenV,255,thresh_callback_green );
	//createTrackbar("iHighGreenV：", GreenWindow,&iHighGreenV,255,thresh_callback_green );
	thresh_callback_green(0,0);
}
void TestYellow()
{
	srcYellow = imread("yellow.bmp",CV_LOAD_IMAGE_COLOR);
	namedWindow("yellow source",CV_WINDOW_NORMAL);
	imshow("yellow source",srcYellow);
	cvtColor(srcYellow,hsvImgYellow,COLOR_BGR2HSV);

	//直方图
	vector<Mat> hsvSplitYellow; 
	split(hsvImgYellow, hsvSplitYellow); 
	equalizeHist(hsvSplitYellow[2],hsvSplitYellow[2]);  
	merge(hsvSplitYellow,hsvImgYellow); 

	namedWindow(YellowWindow, CV_WINDOW_NORMAL);
	createTrackbar("iLowGreenH：", YellowWindow,&iLowYellowH,180,thresh_callback_yellow );
	createTrackbar("iHighGreenH：", YellowWindow,&iHighYellowH,180,thresh_callback_yellow );
	/*createTrackbar("iLowGreenS：", YellowWindow,&iLowYellowS,255,thresh_callback_yellow );
	createTrackbar("iHighGreenS：", YellowWindow,&iHighYellowS,255,thresh_callback_yellow );
	createTrackbar("iLowGreenV：", YellowWindow,&iLowYellowV,255,thresh_callback_yellow );
	createTrackbar("iHighGreenV：", YellowWindow,&iHighYellowV,255,thresh_callback_yellow );*/
	thresh_callback_yellow(0,0);
}
void TestWhite()
{
	srcWhite = imread("white.bmp",CV_LOAD_IMAGE_COLOR);
	namedWindow("white source",CV_WINDOW_NORMAL);
	imshow("white source",srcWhite);

	//blur(srcWhite,srcWhite,Size(3,3));

	cvtColor(srcWhite,hsvImgWhite,COLOR_BGR2HSV);
	/*
	vector<Mat> hsvSplitWhite;
	split(hsvImgWhite,hsvSplitWhite);
	equalizeHist(hsvSplitWhite[2],hsvSplitWhite[2]);
	merge(hsvSplitWhite,hsvImgWhite);*/

	namedWindow(WhiteWindow, CV_WINDOW_NORMAL);
	createTrackbar("iLowWhiteH：", WhiteWindow,&iLowWhiteH,180,thresh_callback_white );
	createTrackbar("iHighWhiteH：", WhiteWindow,&iHighWhiteH,180,thresh_callback_white );
	/*createTrackbar("iLowWhiteS：", WhiteWindow,&iLowWhiteS,255,thresh_callback_white );
	createTrackbar("iHighWhiteS：", WhiteWindow,&iHighWhiteS,255,thresh_callback_white );
	createTrackbar("iLowWhiteV：", WhiteWindow,&iLowWhiteV,255,thresh_callback_white );
	createTrackbar("iHighWhiteV：", WhiteWindow,&iHighWhiteV,255,thresh_callback_white );*/
	thresh_callback_white(0,0);	
}
static void thresh_callback_red(int, void*)
{
	// 查找指定范围内的颜色
	Mat imgThresholdedRed;
	inRange(hsvImgRed, Scalar(iLowRedH, iLowRedS, iLowRedV), Scalar(iHighRedH, iHighRedS, iHighRedV), imgThresholdedRed);

	Mat element = getStructuringElement(MORPH_RECT, Size(5, 5));
	//开操作 (去除一些噪点)
	morphologyEx(imgThresholdedRed, imgThresholdedRed, MORPH_OPEN, element);

	//闭操作 (连接一些连通域)
	morphologyEx(imgThresholdedRed, imgThresholdedRed, MORPH_CLOSE, element);	
	

	//Mat rec = imgThresholdedRed(Rect(100,35,100,80));
	//imshow("rect red",rec);
	/*Scalar scalar = mean(imgThresholdedRed);
	char text[100];
	sprintf(text,"APL:%f",scalar);
	putText(imgThresholdedRed,text,Point(10,150),FONT_HERSHEY_SIMPLEX,1,Scalar(255));*/
	imshow(RedWindow,imgThresholdedRed);
}
static void thresh_callback_green(int, void*)
{
	// 查找指定范围内的颜色
	Mat imgThresholdedGreen;
	inRange(hsvImgGreen, Scalar(iLowGreenH, iLowGreenS, iLowGreenV), Scalar(iHighGreenH, iHighGreenS, iHighGreenV), imgThresholdedGreen);
	//imshow("imgThresholdedGreen",imgThresholdedGreen);
	Mat element = getStructuringElement(MORPH_RECT, Size(10, 10));
	
	//imshow("close",imgThresholdedGreen);
	//开操作 (去除一些噪点)
	morphologyEx(imgThresholdedGreen, imgThresholdedGreen, MORPH_OPEN, element);

	//闭操作 (连接一些连通域)
	morphologyEx(imgThresholdedGreen, imgThresholdedGreen, MORPH_CLOSE, element);

	//imshow("open",imgThresholdedGreen);
	

	//Mat rec = imgThresholdedGreen(Rect(100,35,100,80));
	//imshow("rect green",rec);
	/*Scalar scalar = mean(imgThresholdedGreen);
	char text[100];
	sprintf(text,"APL:%f",scalar);
	putText(imgThresholdedGreen,text,Point(10,150),FONT_HERSHEY_SIMPLEX,1,Scalar(255));*/

	imshow(GreenWindow,imgThresholdedGreen);
}
static void thresh_callback_yellow(int, void*)
{
	// 查找指定范围内的颜色
	Mat imgThresholdedYellow;
	inRange(hsvImgYellow, Scalar(iLowYellowH, iLowYellowS, iLowYellowV), Scalar(iHighYellowH, iHighYellowS, iHighYellowV), imgThresholdedYellow);
	
	Mat element = getStructuringElement(MORPH_RECT, Size(5, 5));
	//开操作 (去除一些噪点)	
	morphologyEx(imgThresholdedYellow, imgThresholdedYellow, MORPH_OPEN, element);

	//闭操作 (连接一些连通域)
	morphologyEx(imgThresholdedYellow, imgThresholdedYellow, MORPH_CLOSE, element);
	
	

	//Mat rec = imgThresholdedYellow(Rect(100,35,100,80));
	//imshow("rect yellow",rec);
	/*Scalar scalar = mean(imgThresholdedYellow);
	char text[100];
	sprintf(text,"APL:%f",scalar);
	putText(imgThresholdedYellow,text,Point(10,150),FONT_HERSHEY_SIMPLEX,1,Scalar(255));*/


	imshow(YellowWindow,imgThresholdedYellow);
}

static void thresh_callback_white(int, void*)
{
	// 查找指定范围内的颜色
	Mat imgThresholdedWhite;
	inRange(hsvImgWhite, Scalar(iLowWhiteH, iLowWhiteS, iLowWhiteV), Scalar(iHighWhiteH, iHighWhiteS, iHighWhiteV), imgThresholdedWhite);
	imshow("white before open",imgThresholdedWhite);
	waitKey(0);

	Mat element = getStructuringElement(MORPH_RECT, Size(7, 7));
	//开操作 (去除一些噪点)	
	morphologyEx(imgThresholdedWhite, imgThresholdedWhite, MORPH_OPEN, element);
	imshow("white open",imgThresholdedWhite);
	waitKey(0);

	//闭操作 (连接一些连通域)
	morphologyEx(imgThresholdedWhite, imgThresholdedWhite, MORPH_CLOSE, element);
	imshow("white close",imgThresholdedWhite);
	waitKey(0);

	

	//Mat rec = imgThresholdedWhite(Rect(100,35,100,80));
	//imshow("rect white",rec);
	/*Scalar scalar = mean(imgThresholdedWhite);
	char text[100];
	sprintf(text,"APL:%f",scalar);
	putText(imgThresholdedWhite,text,Point(10,150),FONT_HERSHEY_SIMPLEX,1,Scalar(255));*/


	imshow(WhiteWindow,imgThresholdedWhite);
}

vector<MatND> getHSVHist(Mat &src){

	//输入图片得是三通道彩色图片
	assert (!src.empty() && src.channels() == 3);

	//rgb转hsv图像
	Mat hsv;
	cvtColor(src, hsv, CV_BGR2HSV);

	//h的范围是0~180，所以选取30个bin
	//s和v的范围都是0~255，那就选择51个bin
	int hbins = 30;
	int sbins = 51;
	int vbins = 51;
	int hHistSize[] = {hbins};
	int sHistSize[] = {sbins};
	int vHistSize[] = {vbins};

	float hranges[] = {0, 180};
	float sranges[] = {0, 255};
	float vranges[] = {0, 255};
	const float* hRanges[] = {hranges};
	const float* sRanges[] = {sranges};
	const float* vRanges[] = {vranges};
	vector<MatND> hist;

	int hChannels[] = {0};
	int sChannels[] = {1};
	int vChannels[] = {2};
	MatND hHist, sHist, vHist;
	calcHist(&hsv, 1, hChannels, Mat(), hHist, 1, hHistSize, hRanges);
	calcHist(&hsv, 1, sChannels, Mat(), sHist, 1, sHistSize, sRanges);
	calcHist(&hsv, 1, vChannels, Mat(), vHist, 1, vHistSize, vRanges);
	hist.push_back(hHist);
	hist.push_back(sHist);
	hist.push_back(vHist);
	normalize( hist[0], hist[0], 0, 1, NORM_MINMAX, -1, Mat() );
	normalize( hist[1], hist[1], 0, 1, NORM_MINMAX, -1, Mat() );
	normalize( hist[2], hist[2], 0, 1, NORM_MINMAX, -1, Mat() );

	int i;  
	int start = -1, end = -1;
	for(i = 0; i < 30; i++)  
	{  
		float value = hist[0].at<float>(i);
		if (value  > 0)
		{
			if (start == -1)
			{
				start = i;
				end = i;
			}
			else
				end = i;
			cout << "H Value" << i << ": " << value << endl;
		}
		else
		{
			if (start != -1)
				cout <<"H:" <<start*6 <<"~"<<(end+1)*6-1<<endl;
			start = end = -1;
		}
	}  
	if (start != -1)
		cout <<"H:" <<start*5 <<"~"<<(end+1)*5-1<<endl;

	start = -1, end = -1;
	for(i = 0; i < 51; i++)  
	{  
		float value = hist[1].at<float>(i);
		if (value  > 0)
		{
			if (start == -1)
			{
				start = i;
				end = i;
			}
			else
				end = i;
			cout << "S Value" << i << ": " << value << endl;
		}
		else
		{
			if (start != -1)
				cout <<"S:"<< start*5 <<"~"<<(end+1)*5-1<<endl;
			start = end = -1;
		}
	}  
	if (start != -1)
		cout <<"S:" <<start*5 <<"~"<<(end+1)*5-1<<endl;

	start = -1, end = -1;
	for(i = 0; i < 51; i++)  
	{  
		float value = hist[2].at<float>(i);
		if (value  > 0)
		{
			if (start == -1)
			{
				start = i;
				end = i;
			}
			else
				end = i;
			cout << "V Value" << i << ": " << value << endl;
		}
		else
		{
			if (start != -1)
				cout <<"V:" <<start*5 <<"~"<<(end+1)*5-1<<endl;
			start = end = -1;
		}
	}  
	if (start != -1)
		cout <<"V:" <<start*5 <<"~"<<(end+1)*5-1<<endl;

	return hist;
}

void filteredRed(const Mat &inputImage, Mat &resultGray, Mat &resultColor)
{
	Mat hsvImage;
	cvtColor(inputImage, hsvImage, CV_BGR2HSV);
	resultGray = Mat(hsvImage.rows, hsvImage.cols,CV_8U,cv::Scalar(0));  
	resultColor = Mat(hsvImage.rows, hsvImage.cols,CV_8UC3,cv::Scalar(0, 0, 0));
	double H=0.0,S=0.0,V=0.0;   
	for(int i=0;i<hsvImage.rows;i++)
	{
		for(int j=0;j<hsvImage.cols;j++)
		{
			H=hsvImage.at<Vec3b>(i,j)[0];
			S=hsvImage.at<Vec3b>(i,j)[1];
			V=hsvImage.at<Vec3b>(i,j)[2];

			if((H >= 0 && H<10) || (H>=125&&H<=180))
			{       
				if(S>=43 && S<=255)
				{
					if(V>=46 && V<=255)
					{
						resultGray.at<uchar>(i,j)=255;
						resultColor.at<Vec3b>(i, j)[0] = inputImage.at<Vec3b>(i, j)[0];
						resultColor.at<Vec3b>(i, j)[1] = inputImage.at<Vec3b>(i, j)[1];
						resultColor.at<Vec3b>(i, j)[2] = inputImage.at<Vec3b>(i, j)[2];
					}					
				}
			}
		}
	}
}
void filteredColor(int colorType,const Mat &inputImage, Mat &resultGray, Mat &resultColor)
{
	int Hmin,Hmax,Hmin2,Hmax2,Smin,Smax,Vmin,Vmax;
	switch(colorType)
	{
	case RED:
		Hmin=0;Hmax=10;Hmin2=125;Hmax2=180;Smin=43;Smax=255;Vmin=46;Vmax=255;
		break;
	case GREEN:
		Hmin=35;Hmax=77;Smin=43;Smax=255;Vmin=46;Vmax=255;
		break;
	case YELLOW:
		Hmin=11;Hmax=25;Smin=43;Smax=255;Vmin=46;Vmax=255;
		break;
	case WHITE:
		Hmin=0;Hmax=180;Smin=0;Smax=30;Vmin=221;Vmax=255;
		break;
	default:
		break;
	}
	Mat hsvImage;
	cvtColor(inputImage, hsvImage, CV_BGR2HSV);
	resultGray = Mat(hsvImage.rows, hsvImage.cols,CV_8U,cv::Scalar(0));  
	resultColor = Mat(hsvImage.rows, hsvImage.cols,CV_8UC3,cv::Scalar(0, 0, 0));
	double H=0.0,S=0.0,V=0.0;   
	for(int i=0;i<hsvImage.rows;i++)
	{
		for(int j=0;j<hsvImage.cols;j++)
		{
			H=hsvImage.at<Vec3b>(i,j)[0];
			S=hsvImage.at<Vec3b>(i,j)[1];
			V=hsvImage.at<Vec3b>(i,j)[2];

			if(S>=Smin && S<=Smax)
			{
				if(V>=Vmin && V<=Vmax)
				{
					if (colorType== RED)
					{
						if((H >= Hmin && H<Hmax) || (H>=Hmin2&&H<=Hmax2))
						{
							resultGray.at<uchar>(i,j)=255;
							resultColor.at<Vec3b>(i, j)[0] = inputImage.at<Vec3b>(i, j)[0];
							resultColor.at<Vec3b>(i, j)[1] = inputImage.at<Vec3b>(i, j)[1];
							resultColor.at<Vec3b>(i, j)[2] = inputImage.at<Vec3b>(i, j)[2];
						}
					}
					else
					{
						if(H >= Hmin && H<Hmax)
						{
							resultGray.at<uchar>(i,j)=255;
							resultColor.at<Vec3b>(i, j)[0] = inputImage.at<Vec3b>(i, j)[0];
							resultColor.at<Vec3b>(i, j)[1] = inputImage.at<Vec3b>(i, j)[1];
							resultColor.at<Vec3b>(i, j)[2] = inputImage.at<Vec3b>(i, j)[2];
						}
					}					
				}					
			}
		}
	}
}
void TestColor(COLOR_TYPE color, char* filename)
{
	char* window_name_source = "source";
	char* window_name_result_gray = "result gray";
	char* window_name_result_gray_after_open = "result";
	Mat inputImage, resultGray,resultColor;
	inputImage = imread(filename,CV_LOAD_IMAGE_COLOR);
	namedWindow(window_name_source,CV_WINDOW_NORMAL);
	imshow(window_name_source,inputImage);
	filteredColor(color,inputImage,resultGray,resultColor);
	namedWindow(window_name_result_gray,CV_WINDOW_NORMAL);
	imshow(window_name_result_gray,resultGray);
	//开操作 (去除一些噪点)
	Mat element = getStructuringElement(MORPH_RECT, Size(5, 5));
	morphologyEx(resultGray, resultGray, MORPH_OPEN, element);
	//闭操作 (连接一些连通域)
	morphologyEx(resultGray, resultGray, MORPH_CLOSE, element);
	namedWindow(window_name_result_gray_after_open,CV_WINDOW_NORMAL);
	imshow(window_name_result_gray_after_open,resultGray);
}