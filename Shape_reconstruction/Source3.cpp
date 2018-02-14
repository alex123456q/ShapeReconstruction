#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#include "atlimage.h"
#include "Processing.h"
#include <fstream>
//#include "Reconstruct.cpp"
bool intersect(Point p1, Point q1, Point p2, Point q2, Point& res);
#include "paint_utils.h"
//#define inverted true
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

namespace Y {
	bool	ISBlack(COLORREF color, bool inverted = false) {
		bool ans = GetBValue(color) < 128
			|| GetRValue(color) < 128
			|| GetGValue(color) < 128;
		if (inverted)
			return !ans;
		return ans;
	}
	bool check(Point* itersecond, Point* itersum)
	{
		return (itersecond->X == itersum->X && itersecond->Y == itersum->Y);
	}
	bool check_eps(Point* itersecond, Point* itersum)
	{
		return (abs(itersecond->X - itersum->X) < 2 && abs(itersecond->Y - itersum->Y) < 2);
	}

	double dist(Point* itersecond, Point* itersum)
	{
		return (pow(itersecond->X - itersum->X, 2) + pow(itersecond->Y - itersum->Y, 2));

	}


	std::vector<std::vector<cv::Point>> selectContours(std::string& filename = std::string("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\figoutline.png"))
	{
		std::vector<std::vector<cv::Point> > contours;
		std::vector<cv::Vec4i> hierarchy;
		cv::Mat src;
		src = cv::imread(filename, 1);
		cvtColor(src, src, CV_BGR2GRAY);
		cv::findContours(src, contours, hierarchy, CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE, cv::Point(0, 0));
		return contours;
	}
	std::vector<Cell> main12(CImage& source1, CImage& source2, int firstH, int secondH)
	{
		CImage source, sourceboth, edge1, edge2, edge;

		int width = min(source1.GetWidth(), source2.GetWidth()), height = min(source1.GetHeight(), source2.GetHeight());

		source.Create(width, height, 32);
		sourceboth.Create(width, height, 32);
		edge1.Create(width, height, 32);
		edge2.Create(width, height, 32);
		edge.Create(width, height, 32);

		BitRaster* srcim1 = new BitRaster(width, height);
		BitRaster* srcim2 = new BitRaster(width, height);
		BitRaster* srcim = new BitRaster(width, height);
		BitRaster* srcimboth = new BitRaster(width, height);
		BitRaster* srcimedge = new BitRaster(width, height);

		for (int i = 0; i < width; ++i)
			for (int j = 0; j < height; ++j) {
				srcim1->setBit(i, j, ISBlack(source1.GetPixel(i, j)));
				srcim2->setBit(i, j, ISBlack(source2.GetPixel(i, j)));
			}

		const TPolFigure const* fig1 = new TPolFigure(srcim1, 0);// AreaIgnore  /*площадь игнорируемых контуров*/
		const TPolFigure const* fig2 = new TPolFigure(srcim2, 0);

		std::vector<std::vector<double>> imageF;

		TConnected* Component1 = fig1->Components->first();
		TConnected* Component2 = fig2->Components->first();            //if connected fig2
		Point*  iterfirst = Component1->Border->ListPoints->first();//Componentsum->HoleList[0]->ListPoints->first();
		Point* itersecond = Component2->Border->ListPoints->first();
		//Component2->Border->ListPoints;
		Point* itermin = itersecond;
		const Point const *  iterfirst1 = Component1->Border->ListPoints->first();//Componentsum->HoleList[0]->ListPoints->first();
		/*const*/ Point const * itersecond2 = Component2->Border->ListPoints->first();

		double distt, mindist = 100000;
		while (Component2) { //only one
			itersecond = Component2->Border->ListPoints->first();
			while (itersecond) {
				distt = dist(itersecond, iterfirst);
				if (distt < mindist)
				{
					mindist = distt;
					itermin = itersecond;
				}
				itersecond = itersecond->getNext();
			}
			Component2 = Component2->getNext();
		}
		itersecond = itermin;
		Component2 = fig2->Components->first();            //if connected fig2
		vector<Cell> cells;
		const TPolFigure const* fig11 = new TPolFigure(srcim1, 0);// AreaIgnore  /*площадь игнорируемых контуров*/
		const TPolFigure const* fig22 = new TPolFigure(srcim2, 0);// AreaIgnore  /*площадь игнорируемых контуров*/

		//while (Component1) { //only one
		iterfirst = Component1->Border->ListPoints->first();
		while (iterfirst) {

			Point* itersecond = itermin;
			while (dist(itermin, iterfirst) < dist(itermin, iterfirst->getNextLooped()))
			{
				Cell newcell = Cell();

				newcell.borders.push_back(std::pair<Point, Point>(Point(iterfirst->X, iterfirst->Y), Point(itermin->X, itermin->Y)));
				newcell.borders_color.push_back(std::pair<double, double>(firstH, secondH));


				newcell.borders.push_back(std::pair<Point, Point>(Point(itermin->X, itermin->Y), Point(itermin->getNextLooped()->X, itermin->getNextLooped()->Y)));
				newcell.borders_color.push_back(std::pair<double, double>(secondH, secondH));

				/*if (!(itermin->getNext()))
					itermin = fig2->Components->first()->Border->ListPoints->first();//fig22->Components->first()->Border->ListPoints->first();
				else
					itermin = itermin->getNextLooped();*/

				cells.push_back(newcell);
				itermin = itermin->getNextLooped();
				if (itermin == itersecond)
					break;
			}


			/*while (iterfirst)
			{
				while (iterfirst && dist(itersecond, iterfirst) < dist(itersecond->getNextLooped(), iterfirst))
				{
					Cell newcell = Cell();

					newcell.borders.push_back(std::pair<Point, Point>(Point(iterfirst->X, iterfirst->Y), Point(itersecond->X, itersecond->Y)));
					newcell.borders_color.push_back(std::pair<double, double>(firstH, secondH));


					newcell.borders.push_back(std::pair<Point, Point>(Point(iterfirst->X, iterfirst->Y), Point(iterfirst->getNextLooped()->X, iterfirst->getNextLooped()->Y)));
					newcell.borders_color.push_back(std::pair<double, double>(firstH, firstH));

					cells.push_back(newcell);
					iterfirst = iterfirst->getNext();
				}
				itersecond = itersecond->getNextLooped();
				if (itersecond == itermin)
					break;
			}*/


			Cell newcell = Cell();
			//if (!(iterfirst->getNext()))
			//	iterfirst = fig11->Components->first()->Border->ListPoints->first();
			//else 
			//	iterfirst = iterfirst->getNext();

			//newcell.borders.push_back(std::pair<Point, Point>(*(itermin), *(itermin->getNextLooped())));
			//newcell.borders_color.push_back(std::pair<double, double>(secondH, secondH));

			newcell.borders.push_back(std::pair<Point, Point>(Point(iterfirst->X, iterfirst->Y), Point(iterfirst->getNextLooped()->X, iterfirst->getNextLooped()->Y)));
			newcell.borders_color.push_back(std::pair<double, double>(firstH, firstH));
			cells.push_back(newcell);
			//		if (iterfirst)
			iterfirst = iterfirst->getNext();


			if (!iterfirst)
				break;

			//Point* iterfirst1 = iterfirst;
			//Point* itermin = itersecond;
			/*while (dist(itermin, iterfirst) < dist(itermin->getNextLooped(), iterfirst))
				{
					Cell newcell = Cell();

					newcell.borders.push_back(std::pair<Point, Point>(Point(iterfirst->X, iterfirst->Y), Point(itermin->X, itermin->Y)));
					newcell.borders_color.push_back(std::pair<double, double>(firstH, secondH));


					newcell.borders.push_back(std::pair<Point, Point>(Point(iterfirst->X, iterfirst->Y), Point(iterfirst->getNextLooped()->X, iterfirst->getNextLooped()->Y)));
					newcell.borders_color.push_back(std::pair<double, double>(secondH, secondH));

					cells.push_back(newcell);
					iterfirst = iterfirst->getNext();
					if (!iterfirst)
						break;
				}*/
		}
		//	Component1 = Component1->getNext();
		//}

		source.Destroy();
		free(srcim1);
		free(srcim2);
		free(srcim);
		//free(fig1);
		//free(fig2);
		//	PaintInFile(rec->imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\after_cur_out.png"));
		return cells;
	}

	int main(int argc, char *argv[]) {
		//static HBITMAP bmpSource = NULL;
		std::vector<std::pair<CImage, int> > sources;
		std::ifstream fin("C:/Users/Alexandra/My/Shape_reconstruction/data/slices.txt");
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
			sources[i].second = height;
		}
		std::vector<std::vector<Cell>> cvec;

		for (int i = 0; i < sources.size() - 1; ++i)
		{
			//if success - better to divide
			cvec.push_back(main12(sources[i].first, sources[i + 1].first, sources[i].second, sources[i + 1].second)); //tmp
		}

		main2(argc, argv, cvec);

		for (int i = 0; i < sources.size(); ++i)
		{
			sources[i].first.Destroy();
		}

		return 0;
	}
}