
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
	//#define inverted true
//namespace old_main
//{

ofstream file3;

double dist(Point& itersecond, Point& itersum)
{
	return (pow(itersecond.X - itersum.X, 2) + pow(itersecond.Y - itersum.Y, 2));
}

double dist(Point* itersecond, Point* itersum)
{
	return (pow(itersecond->X - itersum->X, 2) + pow(itersecond->Y - itersum->Y, 2));
}

std::vector<Point> calc_upper_border(std::vector<Cell>& downer, int h1)
{
	std::vector<Point> downer_up;
	for (int i = 0; i < downer.size(); ++i)
		for (int j = 0; j < downer[i].borders.size(); ++j)
		{
			if (downer[i].borders_color[j].first == h1 && downer[i].borders_color[j].second == h1)
				if (downer_up.size() > 0 &&
					downer[i].borders[j].first.X == downer_up.back().X && downer[i].borders[j].first.Y == downer_up.back().Y)
					break;
				else
					downer_up.push_back(downer[i].borders[j].first);
		}
	return downer_up;
}
void contours_matching(std::vector<Cell>& upper, std::vector<Cell>& downer, int h1)  //downer <-> upper
{
	std::vector<Point> downer_up;
	for (int i = 0; i < downer.size(); ++i)
		for (int j = 0; j < downer[i].borders.size(); ++j)
		{
			if (downer[i].borders_color[j].first == h1 && downer[i].borders_color[j].second == h1)
				if (downer_up.size() > 0 &&
					downer[i].borders[j].first.X == downer_up.back().X && downer[i].borders[j].first.Y == downer_up.back().Y)
					break;
				else
				    downer_up.push_back(downer[i].borders[j].first);
		}
	//std::vector<Point> upper_down;
	double distt, mindist = 1000000;
	int zmin = 0;

	bool f = false;
	//Point p;
	/*for (int i = 0; i < upper.size(); ++i)
	{
		for (int j = 0; j < upper[i].borders.size(); ++j)
		{
			if (upper[i].borders_color[j].first == h1)
			{
				f = true;
				p = upper[i].borders[j].first;
				break;
			}
		}
		if (f)
			break;
	}

	for (int k = 0; k < downer_up.size(); ++k)
	{
		if (distt=dist(p, downer_up[k]) < mindist)
		{
			mindist = distt;
			zmin = k;
		}
	}*/

	for (int i = 0; i < upper.size(); ++i)
		for (int j = 0; j < upper[i].borders.size(); ++j)
		{
			if (upper[i].borders_color[j].first == h1)
			{
				/*while (dist(upper[i].borders[j].first, downer_up[zmin]) > dist(upper[i].borders[j].first, downer_up[(zmin + 1) % (downer_up.size())]))
					zmin = (zmin + 1) % (downer_up.size());
				if (dist(upper[i].borders[j].first, downer_up[zmin]) > dist(upper[i].borders[j].first, downer_up[(zmin + 2) % (downer_up.size())]))
					zmin = (zmin + 2) % (downer_up.size());

				upper[i].borders[j].first = downer_up[zmin];*/
				mindist = 1000000;
				for (int k = 0; k < downer_up.size(); ++k)
				{
					if ( ( distt = dist(upper[i].borders[j].first, downer_up[k]) ) < mindist)
					{
						mindist = distt;
						zmin = k;
					}
				}
				upper[i].borders[j].first = downer_up[zmin];

				if (upper[i].borders_color[j].second == h1)
				{
					int t = 0;

					while (dist(upper[i].borders[j].second, downer_up[zmin]) > dist(upper[i].borders[j].second, downer_up[(zmin + t) % (downer_up.size())]))
						t++;
					//if (dist(upper[i].borders[j].first, downer_up[zmin]) > dist(upper[i].borders[j].first, downer_up[(zmin + 2) % (downer_up.size())]))
					//	t = 2;
					upper[i].borders[j].second = downer_up[(zmin + t) % (downer_up.size())];
				}
				else
				{
					;
				}
				//upper_down.push_back(upper[i].borders[j].first);
				

				//distt = dist(downer_up[0], upper_down.back());

			}
		}
//	if (upper_down.size() != downer_up.size())
//	{
//		std::cout << "Nonoo";
//	}
	/*int j = zmin;
	for (int i = 0; i < upper_down.size(); ++i)
	{
		while()
	}*/
	
}


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

	CImage source0;
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
	std::vector<Cell> main12(CImage& source1, CImage& source2, int firstH_, int secondH_, std::vector<Point> cPoint)
	{
		CImage source, sourceboth, edge1, edge2, edge;
		int firstH = firstH_, secondH = secondH_;

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


		std::vector<std::vector<int>> imcolors(width);
		int sq1 = 0, sq2 = 0;

		for (int i = 0; i < width; ++i) {
			imcolors[i].resize(height);
			for (int j = 0; j < height; ++j) {
				imcolors[i][j] = ISBlack(source2.GetPixel(i, j))*secondH + ISBlack(source1.GetPixel(i, j))*firstH;   //0 .. 255 - getrvalue

				sourceboth.SetPixel(i, j, (RGB(255, 255, 255) - source1.GetPixel(i, j)) & (RGB(255, 255, 255) - source2.GetPixel(i, j)));
				if (i > 0 && j > 0 && i < width - 1 && j < height - 1) {
					if (!ISBlack(source1.GetPixel(i, j)) &&
						(ISBlack(source1.GetPixel(i - 1, j)) ||
							ISBlack(source1.GetPixel(i, j - 1)) ||
							ISBlack(source1.GetPixel(i + 1, j)) ||
							ISBlack(source1.GetPixel(i, j + 1)) ||
							ISBlack(source1.GetPixel(i - 1, j - 1)) ||
							ISBlack(source1.GetPixel(i + 1, j + 1))
							)
						)
						edge1.SetPixel(i, j, RGB(255, 255, 255));
					if (!ISBlack(source2.GetPixel(i, j)) &&
						(ISBlack(source2.GetPixel(i - 1, j)) ||
							ISBlack(source2.GetPixel(i, j - 1)) ||
							ISBlack(source2.GetPixel(i + 1, j)) ||
							ISBlack(source2.GetPixel(i, j + 1)) ||
							ISBlack(source2.GetPixel(i - 1, j - 1)) ||
							ISBlack(source2.GetPixel(i + 1, j + 1))
							)
						)
						edge2.SetPixel(i, j, RGB(255, 255, 255));
					edge.SetPixel(i, j, edge1.GetPixel(i, j) & edge2.GetPixel(i, j));
					//source.SetPixel(i, j, RGB(255, 255, 255));            //if innner???
					   //source2.SetPixel(i, j, source2.GetPixel(i, j) | source0.GetPixel(i, j));
					source.SetPixel(i, j, /*RGB(255, 255, 255)-*/(source1.GetPixel(i, j) ^ source2.GetPixel(i, j) | edge.GetPixel(i, j)));
					//source.SetPixel(i, j, (source1.GetPixel(i, j) & (RGB(255, 255, 255) - source2.GetPixel(i, j) ) | edge.GetPixel(i, j))); //diff
					
					//source0.SetPixel(i, j, source.GetPixel(i, j));

				}
				sq1 += !ISBlack(source1.GetPixel(i, j));
				sq2 += !ISBlack(source2.GetPixel(i, j));
				srcim1->setBit(i, j, ISBlack(source1.GetPixel(i, j)));
				srcim2->setBit(i, j, ISBlack(source2.GetPixel(i, j)));
				srcim->setBit(i, j, ISBlack(source.GetPixel(i, j)/* || edge.GetPixel(i, j)*/, true));
				srcimboth->setBit(i, j, ISBlack(sourceboth.GetPixel(i, j), true));
				srcimedge->setBit(i, j, ISBlack(edge.GetPixel(i, j), true));
			}
		}

		//if (sq1 < sq2)
		//{
		//	firstH = secondH_;
		//	secondH = firstH_;
		//}

		std::vector<std::pair<int, int>> ppps;
		for (int i = 1; i < width - 2; ++i) {
			for (int j = 1; j < height - 2; ++j)
			{
				std::set<int> cols;
				cols.insert(imcolors[i][j - 1]);
				cols.insert(imcolors[i][j + 1]);
				cols.insert(imcolors[i - 1][j]);
				cols.insert(imcolors[i + 1][j]);
				if (cols.size() == 4)
					ppps.push_back(std::pair<double, double>(i, j));
			}
		}

		TPolFigure* fig1 = new TPolFigure(srcim1, 0);// AreaIgnore  /*ïëîùàäü èãíîðèðóåìûõ êîíòóðîâ*/
		TPolFigure* fig2 = new TPolFigure(srcim2, 0);
		TPolFigure* fig = new TPolFigure(srcim, 0);
		TPolFigure* figboth = new TPolFigure(srcimboth, 0);
		TPolFigure* figedge = new TPolFigure(srcimedge, 0);

		std::vector<std::vector<double>> imageF;
		/*	imageF.resize(width);
			for (int i = 0; i < width; ++i) {
				imageF[i].resize(height);
				std::fill(imageF[i].begin(), imageF[i].begin() + height, -1);
			}*/
			//for (int i = 0; i < width; ++i) {
			//	for (int j = 0; j < height; ++j) {
			//		if (/*srcim*//*srcimedge*/srcim1->getBit(i, j))
			//			imageF[i][j] = 1;
			//	}
			//}
		//	PaintBorders(fig1, imageF);
		 //  PaintInFile(imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\fig1.png"));
		/*	    for (int i = 0; i < source1.GetWidth(); ++i){
					imageF[i].resize(source1.GetHeight());
					std::fill (imageF[i].begin(),imageF[i].begin()+source1.GetHeight(), -1);
				}*/
				//	for (int i = 0; i < width; ++i) {
				//		for (int j = 0; j < height; ++j) {
				//			if (/*srcim*//*srcimedge*/srcim2->getBit(i, j))
				//				imageF[i][j] = 1;
				//		}
				//	}
			//	   PaintBorders(fig2, imageF);
			//	   PaintInFile(imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\fig2.png"));
			/*	   for (int i = 0; i < source1.GetWidth(); ++i) {
					   imageF[i].resize(source1.GetHeight());
					   std::fill(imageF[i].begin(), imageF[i].begin() + source1.GetHeight(), -1);
				   }
					for (int i = 0; i < width; ++i) {
						for (int j = 0; j < height; ++j) {
							if (/*srcim*//*srcimedge*///srcim->getBit(i, j))
							/*	   				imageF[i][j] = 1;
										}
									}
								   PaintInFile(imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\figim.png"));*/
		imageF.resize(source1.GetWidth());
		for (int i = 0; i < source1.GetWidth(); ++i) {
			imageF[i].resize(source1.GetHeight());
			std::fill(imageF[i].begin(), imageF[i].begin() + source1.GetHeight(), -1);
		}
		for (int i = 0; i < width; ++i) {
			for (int j = 0; j < height; ++j) {
				if (srcim->getBit(i, j))
					imageF[i][j] = 1;
			}
		}
		//PaintBorders(fig, imageF);
	//	PaintInFile(imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\figulina.png"));
		CImage newimage;
		newimage.Create(imcolors.size(), imcolors[0].size(), 32);
		for (int i = 0; i < imcolors.size()/*image.GetWidth()*/; ++i)
			for (int j = 0; j < imcolors[0].size()/*image.GetHeight()*/; ++j) {
				/*if (imageF[i][j] < -0.5) {
					imageF[i][j] = 0;
				}*/
				int val = imcolors[i][j] / (secondH + firstH + 0.0) * 255; //srcim->getBit(i, j) ? 230 : 110;// imcolors[i][j]*255 / 23;
				//int val = imageF[i][j] * 255;
				newimage.SetPixel(i, j,
					RGB(val, val, val));
			}
		
		/*std::vector<cv::Point> cvpoints; //map with 
		std::vector<std::vector<cv::Point>> cont1 = selectContours(std::string("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\fig1.png"));
		std::vector<std::vector<cv::Point>> cont2 = selectContours(std::string("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\fig2.png"));
		for (int i = 0; i < cont1.size(); ++i)
			for (int j = 0; j < cont1[i].size(); ++j)
				for (int i2 = 0; i2 < cont2.size(); ++i2)
					for (int j2 = 0; j2 < cont2[i2].size(); ++j2)
						if (cont1[i][j].x == cont2[i2][j2].x && cont1[i][j].y == cont2[i2][j2].y)
							cvpoints.push_back(cont1[i][j]);
							*/


		{
			TConnected* Component1 = fig1->Components->first();            //if connected fig2
			while (Component1) {
				Point* itersecond = Component1->Border->ListPoints->first();
				while (itersecond) {
					file3 << itersecond->X << "," << itersecond->Y << "," << secondH_ << std::endl;
					itersecond = itersecond->getNext();
				}
				Component1 = Component1->getNext();
			}
		}
	//	return std::vector<Cell>();
		std::set<Point> sPoints;
		std::set<Point> sPoints1;
		std::set<Point> sPointsNo;

		std::vector<Point> sPointsEdge;
		TConnected* Componentsum = fig->Components->first();
		Point* itersum = NULL;
		if (Componentsum )
			itersum = Componentsum->Border->ListPoints->first();
		bool found = false;
		Point* itersecond = fig1->Components->first()->Border->ListPoints->first();
		int holen = 0;
		while (Componentsum) {
			Point* itersum = Componentsum->Border->ListPoints->first();//Componentsum->HoleList[0]->ListPoints->first();
			bool f = false;
			while (itersum)
			{
				found = false;
				TConnected* Component2 = fig2->Components->first();            //if connected fig2
				while (Component2) {
					itersecond = Component2->Border->ListPoints->first();
					while (itersecond) {
						if (check(itersecond, itersum)) {
							sPoints.insert(Point(itersum->X, itersum->Y));
							newimage.SetPixel(itersum->X, itersum->Y,
								RGB(255, 0, 0));

							/*newimage.SetPixel(itersum->X - 1, itersum->Y - 1,
								RGB(255, 0, 0));
							newimage.SetPixel(itersum->X - 1, itersum->Y,
								RGB(255, 0, 0));
							newimage.SetPixel(itersum->X, itersum->Y - 1,
								RGB(255, 0, 0));
							newimage.SetPixel(itersum->X, itersum->Y + 1,
								RGB(255, 0, 0));
							newimage.SetPixel(itersum->X + 1, itersum->Y,
								RGB(255, 0, 0));
							newimage.SetPixel(itersum->X + 1, itersum->Y + 1,
								RGB(255, 0, 0));*/

							found = true;
							break;
						}
						itersecond = itersecond->getNext();
					}
					if (found)
						break;
					Component2 = Component2->getNext();
				}

				if (!found)

				{
					TConnected* Component1 = fig1->Components->first();            //if connected fig2
					while (Component1) {
						itersecond = Component1->Border->ListPoints->first();
						while (itersecond) {
							if (check(itersecond, itersum)) {
								sPoints1.insert(Point(itersum->X, itersum->Y));
								newimage.SetPixel(itersum->X, itersum->Y,
									RGB(0, 255, 0));

								/*newimage.SetPixel(itersum->X - 1, itersum->Y - 1,
									RGB(0, 255, 0));
								newimage.SetPixel(itersum->X - 1, itersum->Y,
									RGB(0, 255, 0));
								newimage.SetPixel(itersum->X, itersum->Y - 1,
									RGB(0, 255, 0));
								newimage.SetPixel(itersum->X, itersum->Y + 1,
									RGB(0, 255, 0));
								newimage.SetPixel(itersum->X + 1, itersum->Y,
									RGB(0, 255, 0));
								newimage.SetPixel(itersum->X + 1, itersum->Y + 1,
									RGB(0, 255, 0));*/

								found = true;
								break;
							}
							itersecond = itersecond->getNext();
						}
						if (found)
							break;
						Component1 = Component1->getNext();
					}

				}
				if (!found)
				{

				/*	int val = -1;
					int x = itersum->X;
					int y = itersum->Y;
					//set+vector
					std::map<int, int> mvals;
					//for (int i = x - 1; i <= x + 1; ++i)
					//	for (int j = y - 1; j <= y + 1; ++j)
					//	{
					int i = x, j = y;
							mvals[imcolors[i][j]]++;
							mvals[imcolors[i-1][j]]++;
							mvals[imcolors[i+1][j]]++;
							mvals[imcolors[i][j-1]]++;
							mvals[imcolors[i][j+1]]++;


					//	}
					if (mvals.size() == 2) {
						//max
						if ((*mvals.begin()).first == 0)
							val = (*mvals.rbegin()).first;
						else
							val = ((*mvals.begin()).first == firstH ? secondH : firstH);
					}
					else
						val = -1;
					if (val == -1)
						newimage.SetPixel(itersum->X, itersum->Y,
							RGB(0, 0, 255));
					else
						if (val == secondH) {
							newimage.SetPixel(itersum->X, itersum->Y,
								RGB(255, 0, 0));
							sPoints.insert(Point(itersum->X, itersum->Y));
						}
						else {
							newimage.SetPixel(itersum->X, itersum->Y,
								RGB(0, 255, 0));
							sPoints1.insert(Point(itersum->X, itersum->Y));
						}
						*/
					/*newimage.SetPixel(itersum->X, itersum->Y,
						RGB(0, 0, 255));
					newimage.SetPixel(itersum->X-1, itersum->Y-1,
						RGB(0, 0, 255));
					newimage.SetPixel(itersum->X-1, itersum->Y,
						RGB(0, 0, 255));
					newimage.SetPixel(itersum->X, itersum->Y-1,
						RGB(0, 0, 255));
					newimage.SetPixel(itersum->X, itersum->Y+1,
						RGB(0, 0, 255));
					newimage.SetPixel(itersum->X+1, itersum->Y,
						RGB(0, 0, 255));
					newimage.SetPixel(itersum->X+1, itersum->Y+1,
						RGB(0, 0, 255));*/
				}
				itersum = itersum->getNext();
				if (!itersum && holen < Componentsum->HoleList.size()) {
					itersum = Componentsum->HoleList[holen++]->ListPoints->first();
				}
			}
			Componentsum = Componentsum->getNext();
		}
		std::cout << sPoints.size() << " " << sPoints1.size() << " no " << sPointsNo.size() << " e " << sPointsEdge.size() << " ";
		

		newimage.Save(_T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\WithPoints01.png"), Gdiplus::ImageFormatPNG);
		newimage.Destroy();

	/*	std::vector<Point> sPointsCommon;
		TConnected* Component1 = fig1->Components->first();
		while (Component1) {
			Point* iterfirst = Component1->Border->ListPoints->first();//Componentsum->HoleList[0]->ListPoints->first();
			while (iterfirst)
			{
				found = false;
				TConnected* Component2 = fig2->Components->first();            //if connected fig2
				while (Component2) {
					itersecond = Component2->Border->ListPoints->first();
					while (itersecond) {
						if (check(itersecond, iterfirst)) {
							sPointsCommon.push_back(Point(iterfirst->X, iterfirst->Y));
							found = true;
							break;
						}
						itersecond = itersecond->getNext();
					}
					if (found)
						break;
					Component2 = Component2->getNext();
				}
				iterfirst = iterfirst->getNext();
			}
			Component1 = Component1->getNext();
		}*/


		//  PaintBorders2(fig, imageF, sPoints);
		//  PaintInFile(imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\figwithmet.png"));

		Reconstruct* rec = new Reconstruct(source, firstH, secondH);

	/*	if (cPoint.size() > 0)
		{
			TConnected*  Component = rec->skeleton->Components->first();
			double mindist, distt;
			int zmin = 0;
			while (Component) {
				Element* el;
				if (Component->HoleList.size() > 0)
					el = Component->HoleList[0]->Elements[0];
				else
					el = Component->Border->Elements[0];
				while (el)
				{
					mindist = 1000000;
					Point* p;
					if (el->isVertex)
						p = ((Vertex*)el)->p;
					else
						p = ((Edge*)el)->org;
					for (int k = 0; k < cPoint.size(); ++k)
					{
						if ((distt = dist(p, &cPoint[k])) < mindist)
						{
							mindist = distt;
							zmin = k;
						}
					}
					p->X = cPoint[zmin].X;
					p->Y = cPoint[zmin].Y;

					if (!(el->isVertex))
					{
						p = ((Edge*)el)->dest;
						mindist = 1000000;
						for (int k = 0; k < cPoint.size(); ++k)
						{
							if ((distt = dist(p, &cPoint[k])) < mindist)
							{
								mindist = distt;
								zmin = k;
							}
						}
						p->X = cPoint[zmin].X;
						p->Y = cPoint[zmin].Y;
					}


					el = el->getNext();
				}

				Component = Component->getNext();
			}
		}
		*/
		rec->skeleton->MakeTriangDel();
		rec->skeleton->CutSkeleton(0);
		TConnected* Component = rec->skeleton->Components->first();
		rec->im = imcolors;
		while (Component) {
			rec->SetHeightforBorders(Component, sPoints/*No*//*Common*/, sPoints1/*sPointsNo*/, firstH, secondH);

		///	rec->SetHeightforBorders(Component, std::set<Point>()/*sPoints*//*No*//*Common*/, std::set<Point>()/*sPoints1*//*sPointsNo*/, firstH, secondH);


			TBone * Bone = Component->Bones->first();// Nodes->first();
			double mindist = 1000000;
			std::vector<double> mindists(ppps.size(), 10000);
			std::vector<TNode*> Nodes(ppps.size());
			TBone * minBone = Bone;
			TNode *minNode = Bone->org;
			while (Bone) {
				for (int i = 0; i < ppps.size(); ++i)             //1
				{
					double curdist = (pow(Bone->dest->X() - ppps[i].first, 2) + pow(Bone->dest->Y() - ppps[i].second, 2));
					double curdist2 = (pow(Bone->org->X() - ppps[i].first, 2) + pow(Bone->org->Y() - ppps[i].second, 2));
					if (curdist2 < curdist)
						curdist = curdist2;
					if ( curdist < mindists[i])
					{
						minBone = Bone;
						mindists[i] = curdist;
						//minNode = Bone->dest;
						Nodes[i] = Bone->dest;
					}
				}

				Bone = Bone->getNext();
			}

			for (int i = 0; i < Nodes.size(); ++i) {                 //1
				Nodes[i]->Disc->Rad = -1;
				//minNode->Disc->Rad = -1;
			}

			/*{
				TBone * bone = new TBone();
				TNode * newnode = new TNode();

				Point pp2 = get_perpendicular_pt_from_pt_to_line(Point(minBone->dest->X(), minBone->dest->Y()), Point(minBone->org->X(), minBone->org->Y()), Point(ppps[0].first, ppps[0].second));

				newnode->Disc = NewDisc(pp2.X, pp2.Y, -1); //NewDisc(ppps[0].first, ppps[0].second, -1); //without sites and other things now

				//TBone * boneepred = skeleton->Components->first()->Bones->first();
				int found = 0;
				for (int i = 0; i < 3; ++i) {
					if (!minBone->org->Sites[i])
						break;
					for (int j = 0; j < 3; ++j) {
						if (minBone->dest->Sites[j] == minBone->org->Sites[i]) {
							newnode->Sites[found++] = minBone->dest->Sites[j];
							//found = true;
						}
					}
				}
				//newnode->Sites[0] = minBone->dest->Sites[0]; //which similar
				//newnode->Sites[1] = minBone->dest->Sites[1]; // with org

				//bone->Com = skeleton->Components->first();
				bone->dest = minBone->dest;
				bone->org = newnode;
				minBone->dest = newnode;

				bone->dest->DetachBone(minBone);
				bone->dest->AddBone(bone);
				newnode->AddBone(bone);
				newnode->AddBone(minBone);

				bone->moveAsNextFor(minBone);
				newnode->moveAsNextFor(minBone->org);
			}*/
			
			Component = Component->getNext();
		}


		//rec->skeleton->MakeTriangDel();
		//rec->skeleton->CutSkeleton(0);
		rec->mainPart();
		//rec->vertPart(source, cvpoints);

		//rec->vertPart(source, sPointsCommon);
		//rec->vertPart(source, cont1[0]);
		//rec->vertPart(source, cont2[0]);

		//std::vector<std::vector<cv::Point>> cont3 = selectContours();
		//for (int i = 0; i < cont3.size(); ++i)
		//	rec->vertPart(source, cont3[i]);

		for (int i = 0; i < source1.GetWidth(); ++i)
			for (int j = 0; j < source1.GetHeight(); ++j) {
				//if (ISBlack(source2.GetPixel(i, j)))
				//	rec->imageF[i][j] = 1;
				//if (ISBlack(source1.GetPixel(i, j), true))
				//	rec->imageF[i][j] = 0;

				if (!srcim->getBit(i, j))
					rec->imageF[i][j] = 1;
				else
					rec->imageF[i][j] = 0.3;
			}
			PaintSkeletBones(rec->skeleton, rec->imageF, true);
			PaintInFile(rec->imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\fig03.png"));//after_cur_out
			//pr->selectPivot(0, 0);
		//source.Destroy();
		free(srcim1);
		free(srcim2);
		free(srcim);
		free(fig1);
		free(fig2);
		free(fig);
			//PaintInFile(rec->imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\after_cur_out.png"));
		return rec->cells;
	}

	int main(int argc, char *argv[]) {
		//static HBITMAP bmpSource = NULL;
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

		/*fs::recursive_directory_iterator it("C:\\Users\\Alexandra\\My\\Downloads\\Concha\\Concha\\NewSliceSeries\\Concha_020_x1y1x0y90z0x0y0z0copy\\bin_contours\\bins\\"), end;//("./")
		int si = 0;
		while (it != end) {
			if (it->path().extension() != ".jpg") {
				++it;
				continue;
			}
			std::string p = it->path().string();
			sources[si].Load(p.c_str());
			++si;
			++it;
		}*/
	    std::vector<std::vector<Cell>> cvec;
		std::vector<Point> cpoint;
		/*source0.Create(sources[0].first.GetWidth(), sources[0].first.GetHeight(), 32);
		for (int i = 0; i < sources[0].first.GetWidth(); ++i)
			for (int j = 0; j < sources[0].first.GetHeight(); ++j)
				source0.SetPixel(i, j, RGB(0, 0, 0));
		file3.open("C:/Users/Alexandra/My/Shape_reconstruction/data/Grid3.txt", ios::out | ios::app);*/
		for (int i = 0; i < sources.size() - 1; ++i)
		{
			//if success - better to divide
			cvec.push_back(main12(sources[i].first, sources[i + 1].first, sources[i + 1].second, sources[i].second, cpoint)); //tmp
			//cvec.push_back(main12(sources[i+1].first, sources[i].first, sources[i].second, sources[i+1].second, cpoint)); //tmp

			//if (i > 0)
			//	contours_matching(cvec[i-1], cvec[i], sources[i].second);
			//cpoint = calc_upper_border(cvec[i], sources[i+1].second);
		}

		//file3.close();

		/*ofstream file;
		file.open("C:/Users/Alexandra/My/Shape_reconstruction/data/Grid2.txt", ios::out | ios::app);

		std::pair<Point, Point> prevPoints = std::pair<Point, Point>(Point(-1,-1), Point(-1, -1));
		bool firstwas = false;
		for (int k = 0; k < cvec.size(); ++k) {
			for (int i = 0; i < cvec[k].size(); ++i) {
				int j = 0;
				for (auto bord = cvec[k][i].borders.begin(); bord != cvec[k][i].borders.end(); ++bord) {
					if (bord->first.X == bord->second.X && bord->first.Y == bord->second.Y) {
						++j;
						continue;
					}

					if (cvec[k][i].borders_color[j].first != 0) {
						++j;
						continue;
					}
					if (prevPoints.first.X == -1) {
						prevPoints = *bord;
						++j;
						continue;
					}

					if (!firstwas)
					{

						if (prevPoints.first.X == bord->first.X && prevPoints.first.Y == bord->first.Y
							|| prevPoints.first.X == bord->second.X && prevPoints.first.Y == bord->second.Y)
							prevPoints = std::pair<Point, Point>(prevPoints.second, prevPoints.first);
						firstwas = true;
						file << prevPoints.first.X << " " << prevPoints.first.Y << " " << 0 << std::endl;
					}

					if (prevPoints.second.X == bord->first.X && prevPoints.second.Y == bord->first.Y)
						prevPoints = *bord;
					else
						prevPoints = std::pair<Point, Point>(bord->second, bord->first);
					file << prevPoints.second.X << " " << prevPoints.second.Y << " " << 0 << std::endl;
					//file << "(" << bord->first.X << "," << bord->first.Y << "," << cvec[k][i].borders_color[j].first << ")";
					//file << "(" << bord->second.X << "," << bord->second.Y << "," << cvec[k][i].borders_color[j].second << ")";
					//file << "]";
					++j;
				}
				//file << "[";
				//file << "(" << cvec[k][i].skeletbone->dest->X() << "," << cvec[k][i].skeletbone->dest->Y() << "," << cvec[k][i].skeletbone->dest->f << ")";
				//file << "(" << cvec[k][i].skeletbone->org->X() << "," << cvec[k][i].skeletbone->org->Y() << "," << cvec[k][i].skeletbone->org->f << ")";
				//file << "]";
			}
			//file << std::endl;
		}
		file <<  prevPoints.second.X << " " << prevPoints.second.Y << " " << 0 << std::endl;

		prevPoints = std::pair<Point, Point>(Point(-1, -1), Point(-1, -1));
		firstwas = false;
		for (int k = 0; k < cvec.size() ; ++k) {
			for (int i = 0; i < cvec[k].size(); ++i) {
				int j = 0;
				for (auto bord = cvec[k][i].borders.begin(); bord != cvec[k][i].borders.end(); ++bord) {
					if (bord->first.X == bord->second.X && bord->first.Y == bord->second.Y) {
						++j;
						continue;
					}

					if (cvec[k][i].borders_color[j].first != 1) {
						++j;
						continue;
					}
					if (prevPoints.first.X == -1) {
						prevPoints = *bord;
						++j;
						continue;
					}

					if (!firstwas)
					{

						if (prevPoints.first.X == bord->first.X && prevPoints.first.Y == bord->first.Y
							|| prevPoints.first.X == bord->second.X && prevPoints.first.Y == bord->second.Y)
							prevPoints = std::pair<Point, Point>(prevPoints.second, prevPoints.first);
						firstwas = true;
						file << prevPoints.first.X << " " << prevPoints.first.Y << " " << 1 << std::endl;
					}

					if (prevPoints.second.X == bord->first.X && prevPoints.second.Y == bord->first.Y)
						prevPoints = *bord;
					else
						prevPoints = std::pair<Point, Point>(bord->second, bord->first);
					file <<  prevPoints.second.X << " " << prevPoints.second.Y << " " << 1 << std::endl;
					//file << "(" << bord->first.X << "," << bord->first.Y << "," << cvec[k][i].borders_color[j].first << ")";
					//file << "(" << bord->second.X << "," << bord->second.Y << "," << cvec[k][i].borders_color[j].second << ")";
					//file << "]";
					++j;
				}
				//file << "[";
				//file << "(" << cvec[k][i].skeletbone->dest->X() << "," << cvec[k][i].skeletbone->dest->Y() << "," << cvec[k][i].skeletbone->dest->f << ")";
				//file << "(" << cvec[k][i].skeletbone->org->X() << "," << cvec[k][i].skeletbone->org->Y() << "," << cvec[k][i].skeletbone->org->f << ")";
				//file << "]";
			}
			//file << std::endl;
		}
		file <<  prevPoints.second.X << " " << prevPoints.second.Y << " " << 1 << std::endl;
		*/
	/*	int k = cvec.size() - 1;
			for (int i = 0; i < cvec[k].size(); ++i) {
				int j = 0;
				for (auto bord = cvec[k][i].borders.begin(); bord != cvec[k][i].borders.end(); ++bord) {
					if (bord->first.X == bord->second.X && bord->first.Y == bord->second.Y) {
						++j;
						continue;
					}
					if (cvec[k][i].borders_color[j].first != 5) {
						++j;
						continue;
					}
					file << "[";
					file << "(" << bord->first.X << "," << bord->first.Y << "," << cvec[k][i].borders_color[j].first << ")";
					file << "(" << bord->second.X << "," << bord->second.Y << "," << cvec[k][i].borders_color[j].second << ")";
					file << "]";
					++j;
				}
				//file << "[";
				//file << "(" << cvec[k][i].skeletbone->dest->X() << "," << cvec[k][i].skeletbone->dest->Y() << "," << cvec[k][i].skeletbone->dest->f << ")";
				//file << "(" << cvec[k][i].skeletbone->org->X() << "," << cvec[k][i].skeletbone->org->Y() << "," << cvec[k][i].skeletbone->org->f << ")";
				//file << "]";
			}
			file << std::endl;
			for (int i = 0; i < cvec[k].size(); ++i) {
				int j = 0;
				for (auto bord = cvec[k][i].borders.begin(); bord != cvec[k][i].borders.end(); ++bord) {
					if (bord->first.X == bord->second.X && bord->first.Y == bord->second.Y) {
						++j;
						continue;
					}
					if (cvec[k][i].borders_color[j].first != 6) {
						++j;
						continue;
					}
					file << "[";
					file << "(" << bord->first.X << "," << bord->first.Y << "," << cvec[k][i].borders_color[j].first << ")";
					file << "(" << bord->second.X << "," << bord->second.Y << "," << cvec[k][i].borders_color[j].second << ")";
					file << "]";
					++j;
				}
				//file << "[";
				//file << "(" << cvec[k][i].skeletbone->dest->X() << "," << cvec[k][i].skeletbone->dest->Y() << "," << cvec[k][i].skeletbone->dest->f << ")";
				//file << "(" << cvec[k][i].skeletbone->org->X() << "," << cvec[k][i].skeletbone->org->Y() << "," << cvec[k][i].skeletbone->org->f << ")";
				//file << "]";
			}
			file << std::endl;*/
		main2(argc, argv, cvec);

		for (int i = 0; i < sources.size(); ++i)
		{
			sources[i].first.Destroy();
		}
		
		return 0;
	}
//}