#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#include "atlimage.h"
#include "Processing.h"

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


int main(int argc, char *argv[]) {
    //static HBITMAP bmpSource = NULL;
    CImage source1;
    CImage source2, source3;


/*	fs::recursive_directory_iterator it("C:\\Users\\Alexandra\\My\\Downloads\\Concha\\Concha\\NewSliceSeries\\Concha_020_x1y1x0y90z0x0y0z0copy\\bin_contours\\"), end;//("./")
	while (it != end) {
		if (it->path().extension() != ".jpg") {
			++it;
			continue;
		}
		std::string s  = it->path().string();
		source1.Load(s.c_str());
		++it;
		source2.Load(it->path().string().c_str());
		++it;
		source1.Destroy();
		source2.Destroy();
	}
	*/

    //CImage sm_img;
    //bmpSource = (HBITMAP)LoadImage(NULL, "D:\\Shape_reconstruction\\build-untitled1-QT4_8_6-Debug\\cur_out.png", IMAGE_BITMAP, 0, 0, LR_LOADFROMFILE);
    //source = (CImage)(bmpSource);
    //source1.Load("D:\\My documents\\Shape_reconstruction\\data\\input2.png");
    //source2.Load("C:\\Users\\Alexandra\\My\\Downloads\\Concha\\Concha\\NewSliceSeries\\Concha_020_x1y1x0y90z0x0y0z0copy\\bin_contours\\15_lu.jpg");
	source1.Load("C:\\Users\\Alexandra\\My\\Downloads\\Concha\\Concha\\NewSliceSeries\\Concha_020_x1y1x0y90z0x0y0z0copy\\bin_contours\\35_lu.jpg");
    source2.Load("C:\\Users\\Alexandra\\My\\Downloads\\Concha\\Concha\\NewSliceSeries\\Concha_020_x1y1x0y90z0x0y0z0copy\\bin_contours\\25_lu.jpg");
	//source3.Load("C:\\Users\\Alexandra\\My\\My documents\\Shape_reconstruction\\data\\difference.png");//dif input2

//	source2.Load("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\white.png");
//	source1.Load("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\white2.png");

    CImage source, sourceboth;
    source.Create(source1.GetWidth(), source1.GetHeight(), 32);
	sourceboth.Create(source1.GetWidth(), source1.GetHeight(), 32);

    BitRaster* srcim1 = new BitRaster(source1.GetWidth(), source1.GetHeight());
    BitRaster* srcim2 = new BitRaster(source1.GetWidth(), source1.GetHeight());
    BitRaster* srcim = new BitRaster(source1.GetWidth(), source1.GetHeight());

    for (int i = 0; i < source1.GetWidth(); ++i)
        for (int j = 0; j < source1.GetHeight(); ++j){
            source.SetPixel(i, j, /*RGB(255, 255, 255)-*/(source1.GetPixel(i, j) ^ source2.GetPixel(i, j)));
		
			sourceboth.SetPixel(i, j, (RGB(255, 255, 255)-source1.GetPixel(i, j)) & (RGB(255, 255, 255) - source2.GetPixel(i, j)));
			if (i > 0 && j > 0 && i < source1.GetWidth() - 1 && j < source1.GetHeight() - 1) {
				if (!ISBlack(source2.GetPixel(i, j)) &&
					(ISBlack(source2.GetPixel(i - 1, j)) ||
						ISBlack(source2.GetPixel(i, j - 1)) ||
						ISBlack(source2.GetPixel(i + 1, j)) ||
						ISBlack(source2.GetPixel(i, j + 1)) ||
						ISBlack(source2.GetPixel(i - 1, j - 1)) ||
						ISBlack(source2.GetPixel(i + 1, j + 1))
						)
					)
					source.SetPixel(i, j, RGB(255, 255, 255));
			}
            srcim1->setBit(i, j, ISBlack(source1.GetPixel(i, j)));
			srcim2->setBit(i, j, ISBlack(source2.GetPixel(i, j)));
            srcim->setBit(i, j, ISBlack(source.GetPixel(i, j), true));
   }

   TPolFigure* fig1 = new TPolFigure(srcim1, 10);// AreaIgnore  /*ïëîùàäü èãíîðèðóåìûõ êîíòóðîâ*/
   TPolFigure* fig2 = new TPolFigure(srcim2, 10);
   TPolFigure* fig = new TPolFigure(srcim, 10);
   std::vector<std::vector<double>> imageF;
   imageF.resize(source1.GetWidth());
			     for (int i = 0; i < source1.GetWidth(); ++i){
			         imageF[i].resize(source1.GetHeight());
			         std::fill (imageF[i].begin(),imageF[i].begin()+source1.GetHeight(), -1);
			     }
				 for (int i = 0; i < source1.GetWidth(); ++i) {
					 for (int j = 0; j < source1.GetHeight(); ++j) {
						 if (srcim->getBit(i, j))
							 imageF[i][j] = 1;
					 }
				 }
			    //PaintBorders(fig1, imageF);
			    PaintInFile( imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\figoutline.png"));
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
   std::set<Point> sPoints;
   while (itersum)
   {
//       if (!found)
      itersecond = fig1->Components->first()->Border->ListPoints->first();
//        else
//           itersecond = itersecond->getNext();
       while (itersecond){
           if (itersecond->X == itersum->X && itersecond->Y == itersum->Y){
               sPoints.insert(Point(itersum->X, itersum->Y));
               found = true;
               break;
           }
           itersecond = itersecond->getNext();
       }
       itersum = itersum->getNext();
   }
   std::cout << sPoints.size();
 //  PaintBorders2(fig, imageF, sPoints);
 //  PaintInFile(imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\figwithmet.png"));

    Reconstruct*rec = new Reconstruct(source);
	rec->SetHeightforBorders(sPoints);
    rec->mainPart();
	
 	for (int i = 0; i < source1.GetWidth(); ++i)
 		for (int j = 0; j < source1.GetHeight(); ++j) {
 			if (ISBlack(source2.GetPixel(i,j)))
 				rec->imageF[i][j] = 1;
 			if (ISBlack(source1.GetPixel(i, j), true))
 				rec->imageF[i][j] = 0;
 		}
	
//	PaintInFile(rec->imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\after_cur_out.png"));
    //pr->selectPivot(0, 0);
    source.Destroy();
	main2(argc, argv, rec->imageF, rec->cells);
//	PaintInFile(rec->imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\after_cur_out.png"));
    return 0;
}
