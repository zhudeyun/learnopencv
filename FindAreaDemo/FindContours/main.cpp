#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>

using namespace std;
using namespace cv;
Mat src,src_gray,blurimg,edge;
static void thresh_callback(int, void*);
void DrawRect(vector<vector<Point>> contours,Mat result);

int thresh=100;
RNG rng(12345);
char* source_window = "Demo窗口";

void main()
{
	
	src=imread("2.jpg");
	//imshow("src",src);
	cvtColor(src,src_gray,CV_BGR2GRAY);
	//imshow("gray",src_gray);
	blur(src_gray,blurimg,Size(3,3));
	//imshow("blur",blurimg);

	//创建处理窗口
	
	namedWindow(source_window, CV_WINDOW_AUTOSIZE);
	//创建轨迹条
	createTrackbar("threshold1：", source_window,&thresh,255,thresh_callback );
	thresh_callback(0,0);
	waitKey(0);
	return;
}

static void thresh_callback(int, void*)
{
	Mat canny_output;
	vector<vector<Point> > contours;
	vector<Vec4i> hierarchy;

	/// 使用Canndy检测边缘
	Canny(src_gray, canny_output, thresh, thresh * 2, 3);
	imshow(source_window,canny_output);
	/// 找到轮廓
	findContours(canny_output, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));

	/*/// 计算矩
	vector<Moments> mu(contours.size());
	for (int i = 0; i < contours.size(); i++)
	{
		mu[i] = moments(contours[i], false);
	}

	///  计算中心矩:
	vector<Point2f> mc(contours.size());
	for (int i = 0; i < contours.size(); i++)
	{
		mc[i] = Point2f(mu[i].m10 / mu[i].m00, mu[i].m01 / mu[i].m00);
	}*/

	/// 绘制轮廓
	Mat drawing = Mat::zeros(canny_output.size(), CV_8UC3);
	for (int i = 0; i< contours.size(); i++)
	{
		Scalar color = Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
		drawContours(drawing, contours, i, color, 2, 8, hierarchy, 0, Point());
		//circle(drawing, mc[i], 4, color, -1, 8, 0);
	}

	/// 显示到窗口中
	namedWindow("Contours", CV_WINDOW_AUTOSIZE);
	imshow("Contours", drawing);

	///// 通过m00计算轮廓面积并且和OpenCV函数比较
	//printf("\t Info: Area and Contour Length \n");
	//for (int i = 0; i< contours.size(); i++)
	//{
	//	//printf(" * Contour[%d] - Area (M_00) = %.2f - Area OpenCV: %.2f - Length: %.2f \n", i, mu[i].m00, contourArea(contours[i]), arcLength(contours[i], true));
	//	printf(" * Contour[%d] - Area OpenCV: %.2f - Length: %.2f \n", i,  contourArea(contours[i]), arcLength(contours[i], true));
	//	Scalar color = Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
	//	drawContours(drawing, contours, i, color, 2, 8, hierarchy, 0, Point());
	//	//circle(drawing, mc[i], 4, color, -1, 8, 0);
	//}
}

void DrawRect(vector<vector<Point>> contours,Mat result)
{
	// 轮廓表示为一个矩形
	for(int i=0;i<contours.size();i++)
	{
		Rect r = boundingRect(Mat(contours[i]));
		rectangle(result, r, Scalar(255), 2);
	}	
}
