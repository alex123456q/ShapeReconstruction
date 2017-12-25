#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#include "atlimage.h"
#include "Processing.h"
#include <fstream>
//#include "Reconstruct.cpp"
bool intersect (Point p1, Point q1, Point p2, Point q2, Point& res);
#include "paint_utils.h"
//#define inverted true

bool	ISBlack( COLORREF color, bool inverted = false) {
	bool ans = GetBValue(color) < 128
		|| GetRValue(color) < 128
		|| GetGValue(color) < 128;
	if (inverted)
		return !ans;
	return ans;
}


std::vector<Cell> main12(CImage& source1, CImage& source2, int firstH, int secondH)
{
	CImage source, sourceboth;

	int width = min(source1.GetWidth(), source2.GetWidth()), height = min(source1.GetHeight(), source2.GetHeight());
	
	source.Create(width, height, 32);
	//sourceboth.Create(width, height, 32);

	BitRaster* srcim1 = new BitRaster(width, height);
	BitRaster* srcim2 = new BitRaster(width, height);
	BitRaster* srcim = new BitRaster(width, height);

	for (int i = 0; i < width; ++i)
		for (int j = 0; j < height; ++j) {
			source.SetPixel(i, j, /*RGB(255, 255, 255)-*/(source1.GetPixel(i, j) ^ source2.GetPixel(i, j)));

			//sourceboth.SetPixel(i, j, (RGB(255, 255, 255) - source1.GetPixel(i, j)) & (RGB(255, 255, 255) - source2.GetPixel(i, j)));
			/*if (i > 0 && j > 0 && i < width - 1 && j < height - 1) {
				if (!ISBlack(source1.GetPixel(i, j)) &&
					(ISBlack(source1.GetPixel(i - 1, j)) ||
			     		ISBlack(source1.GetPixel(i, j - 1)) ||
						ISBlack(source1.GetPixel(i + 1, j)) ||
						ISBlack(source1.GetPixel(i, j + 1)) ||
						ISBlack(source1.GetPixel(i - 1, j - 1)) ||
						ISBlack(source1.GetPixel(i + 1, j + 1))
						)
					)
					source.SetPixel(i, j, RGB(255, 255, 255));
			}*/
			srcim1->setBit(i, j, ISBlack(source1.GetPixel(i, j)));
			srcim2->setBit(i, j, ISBlack(source2.GetPixel(i, j)));
			srcim->setBit(i, j, ISBlack(source.GetPixel(i, j), true));
		}

	TPolFigure* fig1 = new TPolFigure(srcim1, 10);// AreaIgnore  /*ïëîùàäü èãíîðèðóåìûõ êîíòóðîâ*/
	TPolFigure* fig2 = new TPolFigure(srcim2, 10);
	TPolFigure* fig = new TPolFigure(srcim, 10);
	std::vector<std::vector<double>> imageF;
	imageF.resize(width);
	for (int i = 0; i < width; ++i) {
		imageF[i].resize(height);
		std::fill(imageF[i].begin(), imageF[i].begin() + height, -1);
	}
	for (int i = 0; i < width; ++i) {
		for (int j = 0; j < height; ++j) {
			if (srcim->getBit(i, j))
				imageF[i][j] = 1;
		}
	}
	//PaintBorders(fig1, imageF);
	PaintInFile(imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\figoutline.png"));

	//    for (int i = 0; i < source1.GetWidth(); ++i){
	//        imageF[i].resize(source1.GetHeight());
	//        std::fill (imageF[i].begin(),imageF[i].begin()+source1.GetHeight(), -1);
	//    }
	//    PaintBorders(fig2, imageF);
	//   PaintInFile(imageF, _T("D:\\My documents\\Shape_reconstruction\\data\\fig2.png"));
	/* for (int i = 0; i < source1.GetWidth(); ++i){
	imageF[i].resize(source1.GetHeight());
	std::fill (imageF[i].begin(),imageF[i].begin()+source1.GetHeight(), -1);
	}
	PaintBorders(fig, imageF);
	PaintInFile(imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\figoutline.png"));*/

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

	std::set<Point> sPoints;
	TConnected* Componentsum = fig->Components->first();
	while (Componentsum) {
		Point* itersum = Componentsum->Border->ListPoints->first();
		while (itersum)
		{
			found = false;
			TConnected* Component2 = fig2->Components->first();
			while (Component2){
				itersecond = Component2->Border->ListPoints->first();
				while (itersecond) {
					if (itersecond->X == itersum->X && itersecond->Y == itersum->Y) {
						sPoints.insert(Point(itersum->X, itersum->Y));
						found = true;
						break;
					}
					itersecond = itersecond->getNext();
				}
				if (found)
					break;
				Component2 = Component2->getNext();
			}
			itersum = itersum->getNext();
		}
		Componentsum = Componentsum->getNext();
	}
	std::cout << sPoints.size();
	//  PaintBorders2(fig, imageF, sPoints);
	//  PaintInFile(imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\figwithmet.png"));

	Reconstruct*rec = new Reconstruct(source);
	TConnected* Component = rec->skeleton->Components->first();
	while (Component) {
		rec->SetHeightforBorders(Component, sPoints, firstH, secondH);
		Component = Component->getNext();
	}
	rec->mainPart();

	/*for (int i = 0; i < source1.GetWidth(); ++i)
		for (int j = 0; j < source1.GetHeight(); ++j) {
			if (ISBlack(source2.GetPixel(i, j)))
				rec->imageF[i][j] = 1;
			if (ISBlack(source1.GetPixel(i, j), true))
				rec->imageF[i][j] = 0;
		}*/

	//	PaintInFile(rec->imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\after_cur_out.png"));
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
		cvec.push_back( main12(sources[i].first, sources[i+1].first, sources[i+1].second, sources[i].second) ); //tmp
	}

	main2(argc, argv, cvec);

	for (int i = 0; i < sources.size(); ++i)
	{
		sources[i].first.Destroy();
	}

    return 0;
}
