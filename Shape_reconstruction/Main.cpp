
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
	//#define inverted true
//namespace old_main
//{
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
				source.SetPixel(i, j, /*RGB(255, 255, 255)-*/(source1.GetPixel(i, j) ^ source2.GetPixel(i, j)));

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
				}
				srcim1->setBit(i, j, ISBlack(source1.GetPixel(i, j)));
				srcim2->setBit(i, j, ISBlack(source2.GetPixel(i, j)));
				srcim->setBit(i, j, ISBlack(source.GetPixel(i, j), true));
				srcimboth->setBit(i, j, ISBlack(sourceboth.GetPixel(i, j), true));
				srcimedge->setBit(i, j, ISBlack(edge.GetPixel(i, j), true));
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
				if (/*srcim*/srcimedge->getBit(i, j))
					imageF[i][j] = 1;
			}
		}
		//PaintBorders(fig, imageF);
		PaintInFile(imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\figoutline.png"));

		std::vector<cv::Point> cvpoints; //map with 
		std::vector<std::vector<cv::Point>> cont1 = selectContours(std::string("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\fig1.png"));
		std::vector<std::vector<cv::Point>> cont2 = selectContours(std::string("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\fig2.png"));
		for (int i = 0; i < cont1.size(); ++i)
			for (int j = 0; j < cont1[i].size(); ++j)
				for (int i2 = 0; i2 < cont2.size(); ++i2)
					for (int j2 = 0; j2 < cont2[i2].size(); ++j2)
						if (cont1[i][j].x == cont2[i2][j2].x && cont1[i][j].y == cont2[i2][j2].y)
							cvpoints.push_back(cont1[i][j]);


		Point* itersum = fig->Components->first()->Border->ListPoints->first();
		bool found = false;
		Point* itersecond = fig1->Components->first()->Border->ListPoints->first();
		/*Point* iterfirst = fig2->Components->first()->Border->ListPoints->first();

		//dull intersect
		std::vector<Point> anglepoints;
		while (itersecond) {
			Point* iterfirst = fig2->Components->first()->Border->ListPoints->first();
			while (iterfirst)
			{
				Point res;
				bool b = intersect(*iterfirst, *iterfirst->getNextLooped(),
					*itersecond, *itersecond->getNextLooped(), res);
				if (b) {
					anglepoints.push_back(res);
					break;
				}
				iterfirst = iterfirst->getNext();
			}
			itersecond = itersecond->getNext();
		}*/

		std::vector<Point> borderPoints;
		for (int i = 0; i < edge.GetWidth(); ++i)
			for (int j = 0; j < edge.GetHeight(); ++j)
				if (srcimedge->getBit(i, j) == 1 &&
					srcimedge->getBit(i - 1, j - 1) + srcimedge->getBit(i - 1, j) + srcimedge->getBit(i, j - 1) +
					srcimedge->getBit(i + 1, j + 1) + srcimedge->getBit(i + 1, j) + srcimedge->getBit(i, j + 1) == 1)
				{
					borderPoints.push_back(Point(i, j));
				}

		std::set<Point> sPoints;
		std::set<Point> sPoints1;
		std::set<Point> sPointsNo;

		std::vector<Point> sPointsEdge;
		TConnected* Componentsum = fig->Components->first();

		while (Componentsum) {
			Point* itersum = Componentsum->Border->ListPoints->first();//Componentsum->HoleList[0]->ListPoints->first();
			while (itersum)
			{
				found = false;
				TConnected* Component2 = fig2->Components->first();            //if connected fig2
				while (Component2) {
					itersecond = Component2->Border->ListPoints->first();
					while (itersecond) {
						if (check(itersecond, itersum)) {
							//if (itersecond->X == itersum->X && itersecond->Y == itersum->Y) {
							sPoints.insert(Point(itersum->X, itersum->Y));
							found = true;
							//if (check (itersecond->getNextLooped(), itersum->getNextLooped() ) || check(itersecond->getNextLooped(), itersum->getNextLooped()) )//check(  itersum->getPrevLooped() ) )
							//sPoints.insert(Point(itersum->getNextLooped()->X, itersum->getNextLooped()->Y));
							//sPoints.insert(Point(itersum->getPrevLooped()->X, itersum->getPrevLooped()->Y));

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
					sPointsNo.insert(Point(itersum->X, itersum->Y));
					//sPointsEdge.push_back(Point(itersum->X, itersum->Y));

					/*for (int i = 0; i < borderPoints.size(); ++i)
						if (check_eps(&borderPoints[i], itersum) )
						{
							//found = true;
							sPointsEdge.push_back(*itersum);
							break;
						}
					found = false;
					TConnected* Component3 = figedge->Components->first();            //if connected fig2
					while (Component3) {
						itersecond = Component3->Border->ListPoints->first();
						while (itersecond) {
							if (check(itersecond, itersum)) {
								//if (itersecond->X == itersum->X && itersecond->Y == itersum->Y) {
								sPointsEdge.push_back(Point(itersum->X, itersum->Y));
								found = true;
								//if (check (itersecond->getNextLooped(), itersum->getNextLooped() ) || check(itersecond->getNextLooped(), itersum->getNextLooped()) )//check(  itersum->getPrevLooped() ) )
								//sPoints.insert(Point(itersum->getNextLooped()->X, itersum->getNextLooped()->Y));
								//sPoints.insert(Point(itersum->getPrevLooped()->X, itersum->getPrevLooped()->Y));

								break;
							}
							itersecond = itersecond->getNext();
						}
						if (found)
							break;
						Component3 = Component3->getNext();
					}*/

				}


				itersum = itersum->getNext();
			}
			Componentsum = Componentsum->getNext();
		}
		std::cout << sPoints.size() << " " << sPoints1.size() << " no " << sPointsNo.size() << " e " << sPointsEdge.size() << " ";





		std::vector<Point> sPointsCommon;
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
		}


		//  PaintBorders2(fig, imageF, sPoints);
		//  PaintInFile(imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\figwithmet.png"));

		Reconstruct* rec = new Reconstruct(source, firstH, secondH);
		rec->skeleton->MakeTriangDel();
		rec->skeleton->CutSkeleton(0);
		TConnected* Component = rec->skeleton->Components->first();
		while (Component) {
			rec->SetHeightforBorders(Component, sPoints/*No*//*Common*/, sPoints1/*sPointsNo*/, firstH, secondH);
			Component = Component->getNext();
		}
		//rec->skeleton->MakeTriangDel();
		//rec->skeleton->CutSkeleton(0);
		rec->mainPart();
		//rec->vertPart(source, cvpoints);

		//rec->vertPart(source, sPointsCommon);
		//rec->vertPart(source, cont1[0]);
		//rec->vertPart(source, cont2[0]);

		std::vector<std::vector<cv::Point>> cont3 = selectContours();
		for (int i = 0; i < cont3.size(); ++i)
			rec->vertPart(source, cont3[i]);

		/*for (int i = 0; i < source1.GetWidth(); ++i)
			for (int j = 0; j < source1.GetHeight(); ++j) {
				if (ISBlack(source2.GetPixel(i, j)))
					rec->imageF[i][j] = 1;
				if (ISBlack(source1.GetPixel(i, j), true))
					rec->imageF[i][j] = 0;
			}

			PaintInFile(rec->imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\after_cur_out.png"));*/
			//pr->selectPivot(0, 0);
		source.Destroy();
		free(srcim1);
		free(srcim2);
		free(srcim);
		free(fig1);
		free(fig2);
		free(fig);
		//	PaintInFile(rec->imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\after_cur_out.png"));
		return rec->cells;
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

		for (int i = 0; i < sources.size() - 1; ++i)
		{
			//if success - better to divide
			cvec.push_back(main12(sources[i].first, sources[i + 1].first, sources[i + 1].second, sources[i].second)); //tmp
		}

		main2(argc, argv, cvec);

		for (int i = 0; i < sources.size(); ++i)
		{
			sources[i].first.Destroy();
		}

		return 0;
	}
//}