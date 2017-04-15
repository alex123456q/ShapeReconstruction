#include "atlimage.h"
#include "Processing.h"
bool intersect (Point p1, Point q1, Point p2, Point q2, Point& res);
int main(int argc, char *argv[]) {
    //static HBITMAP bmpSource = NULL;
    CImage source;
    //CImage sm_img;
    //bmpSource = (HBITMAP)LoadImage(NULL, "D:\\Shape_reconstruction\\build-untitled1-QT4_8_6-Debug\\cur_out.png", IMAGE_BITMAP, 0, 0, LR_LOADFROMFILE);
    //source = (CImage)(bmpSource);
    source.Load("D:\\My documents\\Shape_reconstruction\\data\\input.png");
    //HANDLE hBitMap = ::LoadImage(0, "c:\\mybmp.bmp",IMAGE_BITMAP, 0, 0, LR_LOADFROMFILE);
    //CBitmap bmp;
    //bmp.Attach((HBITMAP)hBitMap); 
    //sm_img.Create(500,500,source.GetBPP());
    //CImage cimage;
    //cimage.Create(500, 500, 32);
    //HDC imageHDC = cimage.GetDC();
    ///::BitBlt(imageHDC, 0, 0, 100, 100, bmpSource, 0, 0, SRCCOPY);
    //cimage.Save(L"c:\\test\\fileName.jpg", GUID_NULL);
    //cimage.ReleaseDC();
    //Point *res;
    //intersect(&Point(), &Point(), &Point(), &Point(), res);
    Processing* pr = new Processing(source);
    pr->selectPivot(0, 0);
    source.Destroy();
    return 0;
}