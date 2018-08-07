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
	//����ԭʼͼ  
	Mat src = imread("../cat.jpg");  //����Ŀ¼��Ӧ����һ����Ϊ1.jpg���ز�ͼ
	Mat src1=src.clone();

	//��ʾԭʼͼ 
	imshow("src", src); 

	//----------------------------------------------------------------------------------
	//	һ����򵥵�canny�÷����õ�ԭͼ��ֱ���á�
	//----------------------------------------------------------------------------------
	cvtColor(src,gray,CV_BGR2GRAY);
	Canny(gray,gray,g_nThreshold1,g_nThreshold2);
	imshow("Canny����",gray);
}
void FindEdge()
{
	
	//����������
	namedWindow("Canny����", 2);
	//���������켣��
	//�����켣��
	createTrackbar("threshold1��", "Canny����",&g_nThreshold1,255,On_Canny );
	createTrackbar("threshold2��","Canny����",&g_nThreshold2,255*2,On_Canny );
	On_Canny(0,0);
	waitKey(0);
}