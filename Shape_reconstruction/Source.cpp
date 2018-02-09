//
//  cvFindContours()
//
// http://robocraft.ru
//

//#include <cv.h>
#include <opencv2/opencv.hpp>
#include <opencv2/highgui.hpp>
#include <stdlib.h>
#include <stdio.h>

IplImage* image0 = 0;
IplImage* image = 0;
IplImage* gray = 0;
IplImage* bin = 0;
IplImage* dst = 0;
IplImage* dst0 = 0;

int main34(int argc, char* argv[])
{
	char* filename = argc >= 2 ? argv[1] :// "C:\\Users\\Alexandra\\My\\Downloads\\Concha\\Concha\\NewSliceSeries\\Concha_020_x1y1x0y90z0x0y0z0copy\\17.jpg";
		"C:\\Users\\Alexandra\\My\\Downloads\\Concha\\Concha\\NewSliceSeries\\Concha_020_x1y1x0y90z0x0y0z0copy\\bin_contours\\15_lu_inv.jpg";//"C:\\Users\\Alexandra\\My\\Downloads\\Concha\\Concha\\NewSliceSeries\\Concha_020_x1y1x0y90z0x0y0z0copy\\17.jpg";//"C:\\Users\\Alexandra\\My\\Downloads\\Concha\\Concha\\NewSliceSeries\\Concha_020_x1y1x0y90z0x0y0z0copy\\17.jpg";
	image = cvLoadImage(filename, 1);
	
	printf("[i] image: %s\n", filename);
	assert(image != 0);

	//image0 = cvCloneImage(image);
	//cvSmooth(image, image0, CV_GAUSSIAN, 5, 5);

	gray = cvCreateImage(cvGetSize(image), IPL_DEPTH_8U, 1);
	bin = cvCreateImage(cvGetSize(image), IPL_DEPTH_8U, 1);

	dst = cvCloneImage(image);

	cvNamedWindow("original", CV_WINDOW_AUTOSIZE);
	cvNamedWindow("binary", CV_WINDOW_AUTOSIZE);
	cvNamedWindow("contours", CV_WINDOW_AUTOSIZE);

	cvCvtColor(image, gray, CV_RGB2GRAY);

	cvInRangeS(gray, cvScalar(40), cvScalar(150), bin); // atoi(argv[2])

	CvMemStorage* storage = cvCreateMemStorage(0);
	CvSeq* contours = 0;

	//cv::Mat bin0(*bin);
	//cv::bitwise_not(bin0, bin0); //cv::subtract(cv::Scalar:all(255),src,dst);

	int contoursCont = cvFindContours(bin, storage, &contours, sizeof(CvContour), CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE, cvPoint(0, 0));

	for (CvSeq* seq0 = contours; seq0 != 0; seq0 = seq0->h_next) {
		cvDrawContours(dst, seq0, CV_RGB(255, 216, 0), CV_RGB(0, 0, 250), 0, 1, 8); // 
	}

	cvShowImage("original", image);
	cvShowImage("binary", bin);
	cvShowImage("contours", dst);

	cvWaitKey(0);

	cvReleaseImage(&image);
	cvReleaseImage(&gray);
	cvReleaseImage(&bin);
	cvReleaseImage(&dst);
	cvDestroyAllWindows();
	return 0;
}