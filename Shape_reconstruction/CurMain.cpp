
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#include "atlimage.h"
#include "Processing.h"
#include <fstream>
//#include "Reconstruct.cpp"
bool intersect(Point p1, Point q1, Point p2, Point q2, Point& res);
#include "paint_utils.h"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <assert.h>
#include "opencv2/core/core.hpp"
#include "opencv2/objdetect/objdetect.hpp"
#include <set>
#include <iostream>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/geometries/adapted/boost_polygon/polygon.hpp>

#include <boost/foreach.hpp>
//#include "atlimage.h"
#include <ctime>


namespace
{
	ofstream file3;
	bool	ISBlack(COLORREF color, bool inverted = false) {
		bool ans = GetBValue(color) < 128
			|| GetRValue(color) < 128
			|| GetGValue(color) < 128;
		if (inverted)
			return !ans;
		return ans;
	}
}


std::vector<Cell> main21(CImage& source1, CImage& source2, int firstH_, int secondH_, std::vector<Point> cPoint)
{
	CImage source, sourceboth, edge1, edge2, edge;
	int firstH = firstH_, secondH = secondH_;

	int width = min(source1.GetWidth(), source2.GetWidth()), height = min(source1.GetHeight(), source2.GetHeight());

	BitRaster* srcim1 = new BitRaster(width, height);
	BitRaster* srcim2 = new BitRaster(width, height);

	int t = clock();

	std::vector<std::vector<int>> imcolors(width);
	int sq1 = 0, sq2 = 0;

	for (int i = 0; i < width; ++i) {
		imcolors[i].resize(height);
		for (int j = 0; j < height; ++j) {
			srcim1->setBit(i, j, ISBlack(source1.GetPixel(i, j)));
			srcim2->setBit(i, j, ISBlack(source2.GetPixel(i, j)));
		}
	}
//	int t = clock();

	TPolFigure* fig1 = new TPolFigure(srcim1, 0);
	TPolFigure* fig2 = new TPolFigure(srcim2, 0);

	typedef boost::geometry::model::d2::point_xy<double> point2d;
	typedef boost::geometry::model::polygon<boost::geometry::model::d2::point_xy<double> > polygon;

	//polygon green;//, blue;
	boost::geometry::model::multi_polygon<polygon> multi, green, blue; //can be vector of polygons
	polygon multi1;

	//boost::geometry::read_wkt(
	//	"MULTIPOLYGON(((23 117,24 110,25 105,27 97,30 89,34 81,39 73,42 69,50 60,52 58,59 52,68 46,75 42,77 41,84 38,93 35,101 33,106 32,114 31,131 31,139 32,149 34,161 38,166 40,172 43,177 46,186 52,193 58,197 62,203 69,206 73,210 79,215 89,217 94,218 97,220 104,221 109,222 115,222 134,221 140,220 145,218 152,217 155,215 160,210 170,206 176,203 180,199 185,192 192,186 197,179 202,172 206,166 209,161 211,153 214,145 216,139 217,132 218,113 218,106 217,100 216,96 215,84 211,79 209,73 206,68 203,59 197,53 192,48 187,42 180,39 176,37 173,34 168,30 160,28 155,27 152,25 145,23 133, 23 117)))"
		//	"MULTIPOLYGON(((67 126,68 122,69 119,70 117,74 111,79 106,83 103,88 100,94 97,103 94,107 93,112 92,118 91,135 91,141 92,146 93,150 94,159 97,167 101,170 103,174 106,179 111,182 115,184 119,185 122,186 126,186 134,185 138,184 141,182 145,179 149,173 155,169 158,159 163,150 166,142 168,135 169,118 169,111 168,106 167,94 163,86 159,83 157,79 154,74 149,70 143,69 141,68 138,67 134, 67 126)))"
	//	, green);

//	boost::geometry::read_wkt(
		//	"MULTIPOLYGON(((27 125,28 119,29 116,32 110,34 107,42 99,45 97,51 94,54 93,60 92,68 92,73 93,77 94,85 98,95 108,99 116,100 119,101 124,101 134,100 139,99 142,95 150,85 160,77 164,74 165,69 166,59 166,54 165,51 164,43 160,33 150,29 142,28 139,27 133, 27 125)),((140 126,141 122,143 116,145 112,147 109,151 104,156 100,159 98,163 96,169 94,173 93,184 93,188 94,194 96,200 99,206 104,211 110,214 116,216 122,217 126,217 137,215 145,211 153,207 158,205 160,200 164,192 168,184 170,173 170,169 169,163 167,157 164,151 159,146 153,143 147,141 141,140 137, 140 126)))", blue);
//		"MULTIPOLYGON(((47 79,48 75,50 71,53 67,57 64,61 62,64 61,72 61,75 62,79 64,83 67,86 71,88 75,89 79,89 88,88 91,86 95,83 99,82 100,79 102,75 104,72 105,64 105,61 104,57 102,53 99,50 95,48 91,47 87, 47 79)),((129 89,131 85,135 81,139 79,145 79,148 80,152 83,154 86,155 89,155 95,154 98,152 101,151 102,148 104,145 105,139 105,136 104,133 102,130 98,129 95, 129 89)),((118 154,119 149,120 146,122 142,124 139,131 132,136 129,138 128,146 126,155 126,163 128,165 129,170 132,177 139,179 142,181 146,182 149,183 153,183 163,182 167,181 170,178 176,170 184,167 186,163 188,160 189,155 190,146 190,138 188,134 186,131 184,123 176,120 170,119 167,118 162, 118 154)),((57 146,58 143,59 141,64 136,66 135,70 134,76 134,80 135,82 136,87 141,88 143,89 146,89 151,88 154,87 156,82 161,80 162,77 163,69 163,66 162,64 161,59 156,58 154,57 151, 57 146)))"
//		, blue);


	{
		TConnected* Component1 = fig1->Components->first();            //if connected fig2
		blue.resize(fig1->Components->cardinal());
		int i = 0;
		while (Component1) {
			Point* itersecond = Component1->Border->ListPoints->first();
			//std::vector<point2d> curpoly;
			//blue.resize(i + 1);
			while (itersecond) {
				boost::geometry::append(blue[i].outer(), point2d(itersecond->X, itersecond->Y));
				//curpoly.push_back(point2d(itersecond->X, itersecond->Y) );
				//file3 << itersecond->X << "," << itersecond->Y << "," << secondH_ << std::endl;
				//file3 << itersecond->X << " " << itersecond->Y << ",";
				itersecond = itersecond->getNext();
			}
			boost::geometry::append(blue[i].outer(), point2d(Component1->Border->ListPoints->first()->X, Component1->Border->ListPoints->first()->Y));
			Component1 = Component1->getNext();
			// Create a polygon object and assign the points to it.
			//boost::geometry::model::polygon<point2d> polygon2;
			//boost::geometry::assign_points(polygon2, curpoly);
			++i;
			//blue.push_back(polygon2);
			//file3 << ";";
		}
		Component1 = fig2->Components->first();            //if connected fig2
		i = 0;
		green.resize(fig2->Components->cardinal());
		while (Component1) {
			Point* itersecond = Component1->Border->ListPoints->first();
			//green.resize(i+1);
			while (itersecond) {

				boost::geometry::append(green[i].outer(), point2d(itersecond->X, itersecond->Y));
				//curpoly.outer.push_back(point2d(itersecond->X, itersecond->Y));
				//file3 << itersecond->X << "," << itersecond->Y << "," << secondH_ << std::endl;
				//file3 << itersecond->X << " " << itersecond->Y << ",";
				itersecond = itersecond->getNext();
			}
			boost::geometry::append(green[i].outer(), point2d(Component1->Border->ListPoints->first()->X, Component1->Border->ListPoints->first()->Y));
			++i;
			Component1 = Component1->getNext();
			//green.push_back(curpoly);
			//file3 << ";";
		}
	}

	boost::geometry::sym_difference(green, blue, multi);
	//boost::geometry::difference(green, blue, multi);
	std::vector<boost::geometry::model::d2::point_xy<double>> output_points;
	boost::geometry::intersection(green, blue, output_points);
	//file3 << std::endl;
	
	//newimage.Save(_T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\WithPoints01.png"), Gdiplus::ImageFormatPNG);
	//newimage.Destroy();
	std::vector<std::vector<double> > imageF;
	imageF.resize(600);
	for (int i = 0; i < 600; ++i) {
		imageF[i].resize(600);
		std::fill(imageF[i].begin(), imageF[i].begin() + 256, 0);
	}

	double k = 1.;
	int z = 0;
	BOOST_FOREACH(polygon const& p, multi)
	{
		std::cout << "inners ";
		for (int i = 0; i < p.inners().size(); ++i) {
			for (int j = 0; j < p.inners()[i].size(); ++j)
				std::cout << p.inners()[i].at(j).x() << "-" << p.inners()[i].at(j).y() << " ";
		}
		std::cout << "outers ";
		for (int i = 0; i < p.outer().size(); ++i) {
			//std::cout << p.outer()[i].x() << "-" << p.outer()[i].y() << " ";
			PaintLine(imageF, &Point(p.outer()[i].x(), p.outer()[i].y()), &Point(p.outer()[(i + 1) % (p.outer().size())].x(), p.outer()[(i + 1) % (p.outer().size())].y()), 0.3*k);
		}
		//std::cout << z++ << ": " << boost::geometry::area(p) << std::endl;
		k += 1.;
	}
	PaintInFile(imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\bordsSource6_.png"));
	//return 0;
	std::set<Point> Points;
	for (int i = 0; i < output_points.size(); ++i)
	{
		Points.insert(Point(output_points[i].x(), output_points[i].y()));
	}
	TPolFigure* fig = new TPolFigure(multi, 0);
	fig->MakeTriangDel();
	fig->CutSkeleton(0);

	Reconstruct rec = Reconstruct();
	rec.skeleton = fig;



	TConnected* Component = rec.skeleton->Components->first();
	while (Component) {
	/*	int curcolor = -1;
		bool converted = false;
		Element* el = Component->Border->Elements[0];  //fter only vertices
		while (el && curcolor == -1)
		{
			if (!el->isVertex)
			{
				el = el->getNext();
				continue;
			}
			double x = ((Vertex*)el)->p->X;
			double y = ((Vertex*)el)->p->Y;
			if (Points.find(Point(x, y)) != Points.end())
			{
				converted = true;
			}
			bool inblue = boost::geometry::covered_by(boost::geometry::model::d2::point_xy<double>(x, y), blue);
			bool ingreen = boost::geometry::covered_by(boost::geometry::model::d2::point_xy<double>(x, y), green);
			if (inblue && !ingreen)
				curcolor = firstH;
			if (!inblue && ingreen)
				curcolor = secondH;
			el = el->getNext();
		}
		if (curcolor == -1)
			std::cout << "beda";
		if (converted)
			curcolor = firstH + secondH - curcolor;
		rec.SetHeightforBorders2(Component, Points, firstH, secondH, curcolor);*/
		double x = ((Vertex*)Component->Border->Elements[0])->p->X;
		double y = ((Vertex*)Component->Border->Elements[0])->p->Y;
		bool inblue = boost::geometry::covered_by(boost::geometry::model::d2::point_xy<double>(x, y), blue);
		bool ingreen = boost::geometry::covered_by(boost::geometry::model::d2::point_xy<double>(x, y), green);
		if (inblue )//&& !boost::geometry::covered_by(boost::geometry::model::d2::point_xy<double>(x, y), green) )//|| boost::geometry::within(boost::geometry::model::d2::point_xy<double>(x, y), blue))            //green big floor   within
			rec.SetHeightforBorders2(Component, Points, firstH, secondH, firstH);
		else if ( ingreen )//&& !boost::geometry::covered_by(boost::geometry::model::d2::point_xy<double>(x, y), blue) )//|| boost::geometry::within(boost::geometry::model::d2::point_xy<double>(x, y), green))
			rec.SetHeightforBorders2(Component, Points, firstH, secondH, secondH);
		else
		{
			double x = ((Vertex*)Component->Border->Elements[1])->p->X;
			double y = ((Vertex*)Component->Border->Elements[1])->p->Y;
			if (boost::geometry::covered_by(boost::geometry::model::d2::point_xy<double>(x, y), blue) )//&& !boost::geometry::covered_by(boost::geometry::model::d2::point_xy<double>(x, y), green))            //green big floor   within
				rec.SetHeightforBorders2(Component, Points, firstH, secondH, firstH);
			else if (boost::geometry::covered_by(boost::geometry::model::d2::point_xy<double>(x, y), green) )//&& !boost::geometry::covered_by(boost::geometry::model::d2::point_xy<double>(x, y), blue))
				rec.SetHeightforBorders2(Component, Points, firstH, secondH, secondH);
			else
			{
				double x = ((Vertex*)Component->Border->Elements[Component->Border->NumbElem-1])->p->X;
				double y = ((Vertex*)Component->Border->Elements[Component->Border->NumbElem - 1])->p->Y;
				if (boost::geometry::covered_by(boost::geometry::model::d2::point_xy<double>(x, y), blue) )//&& !boost::geometry::covered_by(boost::geometry::model::d2::point_xy<double>(x, y), green))            //green big floor   within
					rec.SetHeightforBorders2(Component, Points, firstH, secondH, firstH);
				else if (boost::geometry::covered_by(boost::geometry::model::d2::point_xy<double>(x, y), green) )//&& !boost::geometry::covered_by(boost::geometry::model::d2::point_xy<double>(x, y), blue))
					rec.SetHeightforBorders2(Component, Points, firstH, secondH, secondH);
				else				rec.SetHeightforBorders2(Component, Points, firstH, secondH, secondH); //std::cout << "beda";
			}
		}
		Component = Component->getNext();
	}

	rec.mainPart();
	std::cout << std::endl << "time" << (clock() - t)/1000.0 << std::endl;
	free(srcim1);
	free(srcim2);
	free(fig1);
	free(fig2);
	//PaintInFile(rec->imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\after_cur_out.png"));
	return rec.cells;
}

int main7(int argc, char *argv[]) {
	
	std::vector<std::pair<CImage, int> > sources;
	std::ifstream fin("C:/Users/Alexandra/My/Shape_reconstruction/data/slices.txt"); //yesno.txt   slices.txt
	std::string path;
	fin >> path;
	int size;
	fin >> size;
	sources.resize(size);
	for (int i = 0; i < size; ++i)
	{
		std::string name;
		int height;
		fin >> name >> height;
		sources[i].first.Load((path + name).c_str());
		//sources[i].first = cv::imread( (path + name).c_str());
		sources[i].second = height;
		//cv::Mat dst;
		//cv::resize(sources[i].first, dst, cv::Size(64, 128));
		//cv::imwrite(path + "res\\\\" +name, dst);
	}

	std::vector<std::vector<Cell>> cvec;
	std::vector<Point> cpoint;
	file3.open("C:/Users/Alexandra/My/Shape_reconstruction/data/Grid_more.txt", ios::out | ios::app);//*/
	for (int i = 0; i < sources.size() - 1; ++i)
	{
		//if success - better to divide
		cvec.push_back(main21(sources[i].first, sources[i + 1].first, sources[i].second, sources[i+1].second, cpoint)); //tmp
																														  //cvec.push_back(main12(sources[i+1].first, sources[i].first, sources[i].second, sources[i+1].second, cpoint)); //tmp

																														  //if (i > 0)
																														  //	contours_matching(cvec[i-1], cvec[i], sources[i].second);
																														  //cpoint = calc_upper_border(cvec[i], sources[i+1].second);
	}

	
	file3.close();

	main2(argc, argv, cvec);

	for (int i = 0; i < sources.size(); ++i)
	{
		sources[i].first.Destroy();
	}

	return 0;
}
//}