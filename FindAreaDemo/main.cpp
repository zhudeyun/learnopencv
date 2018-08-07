#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace std;
using namespace cv;

int g_nThreshold1 = 100;
int g_nThreshold2 = 100;

static void on_Track();
void FindEdge();

void main()
{
	FindEdge();
	return;
}

static void On_Canny(int, void*)
{
	Mat gray;
	//载入原始图  
	Mat src = imread("../cat.jpg");  //工程目录下应该有一张名为1.jpg的素材图
	Mat src1=src.clone();

	//显示原始图 
	imshow("src", src); 

	//----------------------------------------------------------------------------------
	//	一、最简单的canny用法，拿到原图后直接用。
	//----------------------------------------------------------------------------------
	cvtColor(src,gray,CV_BGR2GRAY);
	Canny(gray,gray,g_nThreshold1,g_nThreshold2);
	imshow("Canny窗口",gray);
}
void FindEdge()
{
	
	//创建处理窗口
	namedWindow("Canny窗口", 2);
	//创建两条轨迹条
	//创建轨迹条
	createTrackbar("threshold1：", "Canny窗口",&g_nThreshold1,255,On_Canny );
	createTrackbar("threshold2：","Canny窗口",&g_nThreshold2,255*2,On_Canny );
	On_Canny(0,0);
	waitKey(0);
}