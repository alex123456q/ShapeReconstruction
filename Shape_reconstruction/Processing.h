#ifndef MOVESKELETW_H
#define MOVESKELETW_H

#include <vector>
#include <stdio.h>
#include <iostream>
#include "SkeletonDemoGUI/SkeletonLib/BSTrans.h"
#include "atlimage.h"
#include <set>
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

struct point_compare {
    bool operator() (const std::pair<Point, Point>& lhs, const std::pair<Point, Point>& rhs) const{
        if (lhs.first.X < rhs.first.X)
            return true;
        if (lhs.first.X > rhs.first.X)
            return false;
        if (lhs.first.Y < rhs.first.Y)
            return true;
        if (lhs.first.Y > rhs.first.Y)
            return false;
        if (lhs.second.X < rhs.second.X)
            return true;
        if (lhs.second.X > rhs.second.X)
            return false;
        return (lhs.second.Y < rhs.second.Y);
    }
};
struct one_point_compare {
    bool operator() (const Point& lhs, const Point& rhs) const{
        if (lhs.X < rhs.X)
            return true;
        if (lhs.X > rhs.X)
            return false;
        return (lhs.Y < rhs.Y);
    }
};
struct Cell{
    int leftx, rightx, upy, downy;
    //vector<TBone> bones; 
    TNode* skeletnode;
    TBone* skeletbone;
    vector<Point> nodes;
    set/*vector*/<std::pair<Point, Point>, point_compare> borders;
    vector<Element*> cellels;

//     bool operator< (const Cell& another){
//         return leftx < another.leftx;
//     }
};
class Reconstruct{
    TPolFigure *skeleton;
    BitRaster *srcimg;
    vector<vector<double> > imageF;
    //vector< vector <Cell> > cells; 
    vector<Cell> cells; 
    CImage image;
public:
    Reconstruct(CImage im);
    ~Reconstruct();
    void makeSkelet();
    void makeCells();
    int mainPart();
    void SetInnerPointsofSkelet();
    void findClosestBorder(Cell& curcell, int i, int j, double&x, double&y, double& d1, double& f, std::vector<std::vector<double>>& imageF);
    void findClosestBone(Cell& curcell, int i, int j, double x, double y, double& x2, double& y2, double& d3, double& f, std::vector<std::vector<double>>& imageF);
};