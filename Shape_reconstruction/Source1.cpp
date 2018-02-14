#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace cv;
using namespace std;

//namespace X {

	Mat src; Mat src_gray;
	int thresh = 100;
	int max_thresh = 255;
	RNG rng(12345);

	/// Function header
	void thresh_callback(int, void*);

	/** @function main */
	int main09(int argc, char** argv)
	{
		/// Load source image and convert it to gray
		char* filename = argc >= 2 ? argv[1] : "C:\\Users\\Alexandra\\My\\Downloads\\Concha\\Concha\\NewSliceSeries\\Concha_020_x1y1x0y90z0x0y0z0copy\\17.jpg";
		//"C:\\Users\\Alexandra\\My\\Downloads\\Concha\\Concha\\NewSliceSeries\\Concha_020_x1y1x0y90z0x0y0z0copy\\bin_contours\\15_lu_inv.jpg";//"C:\\Users\\Alexandra\\My\\Downloads\\Concha\\Concha\\NewSliceSeries\\Concha_020_x1y1x0y90z0x0y0z0copy\\17.jpg";//"C:\\Users\\Alexandra\\My\\Downloads\\Concha\\Concha\\NewSliceSeries\\Concha_020_x1y1x0y90z0x0y0z0copy\\17.jpg";
		//"C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\fig0.png";
		//"C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\figoutline.png";
		src = imread(filename, 1);

		/// Convert image to gray and blur it
		cvtColor(src, src_gray, CV_BGR2GRAY);
		blur(src_gray, src_gray, Size(3, 3));
		//src_gray = GaussianBlur(src_gray, Size(5, 5), Size(0));
		threshold(src_gray, src_gray, 20, 255, THRESH_BINARY + THRESH_OTSU);
		//threshold(src_gray, dst, threshold_value, max_BINARY_value, threshold_type);
		namedWindow("thresh", CV_WINDOW_AUTOSIZE);
		imshow("thresh", src_gray);

		/// Create Window
		char* source_window = "Source";
		namedWindow(source_window, CV_WINDOW_AUTOSIZE);
		imshow(source_window, src);

		createTrackbar(" Canny thresh:", "Source", &thresh, max_thresh, thresh_callback);
		thresh_callback(0, 0);

		waitKey(0);
		return(0);
	}

	/** @function thresh_callback */
	void thresh_callback(int, void*)
	{
		Mat canny_output;
		vector<vector<Point> > contours;
		vector<Vec4i> hierarchy;

		/// Detect edges using canny
		Canny(src_gray, canny_output, thresh, thresh * 2, 3);

		namedWindow("Canny", CV_WINDOW_AUTOSIZE);
		imshow("Canny", canny_output);

		/// Find contours
		findContours(canny_output, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));

		/// Draw contours
		Mat drawing = Mat::zeros(canny_output.size(), CV_8UC3);
		for (int i = 0; i < contours.size(); i++)
		{
			Scalar color = Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
			drawContours(drawing, contours, i, color, 2, 8, hierarchy, 0, Point());
		}

		/// Show in a window
		namedWindow("Contours", CV_WINDOW_AUTOSIZE);
		imshow("Contours", drawing);
	}
	/*
	vector<vector<Point>> selectContours(std::string& filename = std::string("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\figoutline.png") )
	{
		vector<vector<Point> > contours;
		vector<Vec4i> hierarchy;
		Mat src;
		src = imread(filename, 1);
		findContours(src, contours, hierarchy, CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));
		return contours;
	}*/
//}