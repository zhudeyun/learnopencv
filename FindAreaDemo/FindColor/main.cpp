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
int iLowRedH = 156;
int iHighRedH = 180;
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
//Yellow
int iLowYellowH = 26;
int iHighYellowH = 34;
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

static void thresh_callback_red(int, void*);
static void thresh_callback_green(int, void*);
static void thresh_callback_yellow(int, void*);
static void thresh_callback_white(int, void*);

void main()
{
	//red
	srcRed = imread("red.bmp",CV_LOAD_IMAGE_COLOR);
	namedWindow("red source",CV_WINDOW_NORMAL);
	imshow("red source",srcRed);
	cvtColor(srcRed,hsvImgRed,COLOR_BGR2HSV);
	namedWindow(RedWindow, CV_WINDOW_AUTOSIZE);
	createTrackbar("iLowRedH：", RedWindow,&iLowRedH,180,thresh_callback_red );
	createTrackbar("iHighRedH：", RedWindow,&iHighRedH,180,thresh_callback_red );
	//createTrackbar("iLowRedS：", RedWindow,&iLowRedS,255,thresh_callback_red );
	//createTrackbar("iHighRedS：", RedWindow,&iHighRedS,255,thresh_callback_red );
	//createTrackbar("iLowRedV：", RedWindow,&iLowRedV,255,thresh_callback_red );
	//createTrackbar("iHighRedV：", RedWindow,&iHighRedV,255,thresh_callback_red );
	thresh_callback_red(0,0);

	//green
	srcGreen = imread("green.bmp",CV_LOAD_IMAGE_COLOR);
	namedWindow("green source",CV_WINDOW_NORMAL);
	imshow("green source",srcGreen);
	cvtColor(srcGreen,hsvImgGreen,COLOR_BGR2HSV);
	namedWindow(GreenWindow, CV_WINDOW_AUTOSIZE);
	createTrackbar("iLowGreenH：", GreenWindow,&iLowGreenH,180,thresh_callback_green );
	createTrackbar("iHighGreenH：", GreenWindow,&iHighGreenH,180,thresh_callback_green );
	//createTrackbar("iLowGreenS：", GreenWindow,&iLowGreenS,255,thresh_callback_green );
	//createTrackbar("iHighGreenS：", GreenWindow,&iHighGreenS,255,thresh_callback_green );
	//createTrackbar("iLowGreenV：", GreenWindow,&iLowGreenV,255,thresh_callback_green );
	//createTrackbar("iHighGreenV：", GreenWindow,&iHighGreenV,255,thresh_callback_green );
	thresh_callback_green(0,0);

	//yellow
	srcYellow = imread("yellow.bmp",CV_LOAD_IMAGE_COLOR);
	namedWindow("yellow source",CV_WINDOW_NORMAL);
	imshow("yellow source",srcYellow);
	cvtColor(srcYellow,hsvImgYellow,COLOR_BGR2HSV);
	namedWindow(YellowWindow, CV_WINDOW_AUTOSIZE);
	createTrackbar("iLowGreenH：", YellowWindow,&iLowYellowH,180,thresh_callback_yellow );
	createTrackbar("iHighGreenH：", YellowWindow,&iHighYellowH,180,thresh_callback_yellow );
	/*createTrackbar("iLowGreenS：", YellowWindow,&iLowYellowS,255,thresh_callback_yellow );
	createTrackbar("iHighGreenS：", YellowWindow,&iHighYellowS,255,thresh_callback_yellow );
	createTrackbar("iLowGreenV：", YellowWindow,&iLowYellowV,255,thresh_callback_yellow );
	createTrackbar("iHighGreenV：", YellowWindow,&iHighYellowV,255,thresh_callback_yellow );*/
	thresh_callback_yellow(0,0);

	//white
	srcWhite = imread("white.bmp",CV_LOAD_IMAGE_COLOR);
	namedWindow("white source",CV_WINDOW_NORMAL);
	imshow("white source",srcWhite);
	cvtColor(srcWhite,hsvImgWhite,COLOR_BGR2HSV);
	namedWindow(WhiteWindow, CV_WINDOW_AUTOSIZE);
	createTrackbar("iLowWhiteH：", WhiteWindow,&iLowWhiteH,180,thresh_callback_white );
	createTrackbar("iHighWhiteH：", WhiteWindow,&iHighWhiteH,180,thresh_callback_white );
	/*createTrackbar("iLowWhiteS：", WhiteWindow,&iLowWhiteS,255,thresh_callback_white );
	createTrackbar("iHighWhiteS：", WhiteWindow,&iHighWhiteS,255,thresh_callback_white );
	createTrackbar("iLowWhiteV：", WhiteWindow,&iLowWhiteV,255,thresh_callback_white );
	createTrackbar("iHighWhiteV：", WhiteWindow,&iHighWhiteV,255,thresh_callback_white );*/
	thresh_callback_white(0,0);
	
	waitKey(0);
	
	return;
}
static void thresh_callback_red(int, void*)
{
	// 查找指定范围内的颜色
	Mat imgThresholdedRed;
	inRange(hsvImgRed, Scalar(iLowRedH, iLowRedS, iLowRedV), Scalar(iHighRedH, iHighRedS, iHighRedV), imgThresholdedRed);
	Mat element = getStructuringElement(MORPH_RECT, Size(5, 5));
	//闭操作 (连接一些连通域)
	morphologyEx(imgThresholdedRed, imgThresholdedRed, MORPH_CLOSE, element);
	//开操作 (去除一些噪点)
	morphologyEx(imgThresholdedRed, imgThresholdedRed, MORPH_OPEN, element);	

	//Mat rec = imgThresholdedRed(Rect(100,35,100,80));
	//imshow("rect red",rec);
	Scalar scalar = mean(imgThresholdedRed);
	char text[100];
	sprintf(text,"APL:%f",scalar);
	putText(imgThresholdedRed,text,Point(10,150),FONT_HERSHEY_SIMPLEX,1,Scalar(255));
	imshow(RedWindow,imgThresholdedRed);
}
static void thresh_callback_green(int, void*)
{
	// 查找指定范围内的颜色
	Mat imgThresholdedGreen;
	inRange(hsvImgGreen, Scalar(iLowGreenH, iLowGreenS, iLowGreenV), Scalar(iHighGreenH, iHighGreenS, iHighGreenV), imgThresholdedGreen);
		
	Mat element = getStructuringElement(MORPH_RECT, Size(5, 5));
	//闭操作 (连接一些连通域)
	morphologyEx(imgThresholdedGreen, imgThresholdedGreen, MORPH_CLOSE, element);
	//开操作 (去除一些噪点)
	morphologyEx(imgThresholdedGreen, imgThresholdedGreen, MORPH_OPEN, element);

	//Mat rec = imgThresholdedGreen(Rect(100,35,100,80));
	//imshow("rect green",rec);
	Scalar scalar = mean(imgThresholdedGreen);
	char text[100];
	sprintf(text,"APL:%f",scalar);
	putText(imgThresholdedGreen,text,Point(10,150),FONT_HERSHEY_SIMPLEX,1,Scalar(255));

	imshow(GreenWindow,imgThresholdedGreen);
}
static void thresh_callback_yellow(int, void*)
{
	// 查找指定范围内的颜色
	Mat imgThresholdedYellow;
	inRange(hsvImgYellow, Scalar(iLowYellowH, iLowYellowS, iLowYellowV), Scalar(iHighYellowH, iHighYellowS, iHighYellowV), imgThresholdedYellow);
	
	Mat element = getStructuringElement(MORPH_RECT, Size(5, 5));
	//闭操作 (连接一些连通域)
	morphologyEx(imgThresholdedYellow, imgThresholdedYellow, MORPH_CLOSE, element);
	//开操作 (去除一些噪点)	
	morphologyEx(imgThresholdedYellow, imgThresholdedYellow, MORPH_OPEN, element);
	

	//Mat rec = imgThresholdedYellow(Rect(100,35,100,80));
	//imshow("rect yellow",rec);
	Scalar scalar = mean(imgThresholdedYellow);
	char text[100];
	sprintf(text,"APL:%f",scalar);
	putText(imgThresholdedYellow,text,Point(10,150),FONT_HERSHEY_SIMPLEX,1,Scalar(255));


	imshow(YellowWindow,imgThresholdedYellow);
}

static void thresh_callback_white(int, void*)
{
	// 查找指定范围内的颜色
	Mat imgThresholdedWhite;
	inRange(hsvImgWhite, Scalar(iLowWhiteH, iLowWhiteS, iLowWhiteV), Scalar(iHighWhiteH, iHighWhiteS, iHighWhiteV), imgThresholdedWhite);

	Mat element = getStructuringElement(MORPH_RECT, Size(5, 5));
	//闭操作 (连接一些连通域)
	morphologyEx(imgThresholdedWhite, imgThresholdedWhite, MORPH_CLOSE, element);
	//开操作 (去除一些噪点)	
	morphologyEx(imgThresholdedWhite, imgThresholdedWhite, MORPH_OPEN, element);
	

	//Mat rec = imgThresholdedWhite(Rect(100,35,100,80));
	//imshow("rect white",rec);
	Scalar scalar = mean(imgThresholdedWhite);
	char text[100];
	sprintf(text,"APL:%f",scalar);
	putText(imgThresholdedWhite,text,Point(10,150),FONT_HERSHEY_SIMPLEX,1,Scalar(255));


	imshow(WhiteWindow,imgThresholdedWhite);
}