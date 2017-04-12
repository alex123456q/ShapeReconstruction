#ifndef MOVESKELETW_H
#define MOVESKELETW_H

#include <vector>
#include <stdio.h>
#include <iostream>
#include "SkeletonDemoGUI/SkeletonLib/BSTrans.h"
#include "atlimage.h"

class Processing
{
public:
    Processing(CImage im);
    ~Processing();
private:
    TPolFigure *skeleton;
    //HBITMAP image;
    BitRaster *srcimg;
    //int* imageF;
    vector<vector<double> > imageF;//(RR, std::vector<int>(CC));
    CImage image;
    bool point;
    int xPr, yPr;
    TNode* curNode;
    TNode* endPoint;//TNode*
    TNode* endPointst;
    TPoint* start;
    TPoint* finish;
    TPoint* tmpCorn;
    int tmp;
    std::vector<Point> circPoint;
    std::vector<TNode* > curRot;
    std::vector<TNode *> vertices;
public:
    void selectPivot(int px, int py);

    void renewSkelet();
    void Circles(double x, double y);
    void dfs(TNode* curNode);
    void changeSkelet(TNode* curNode, double x, double y);
};


#endif // MOVESKELETW_H
