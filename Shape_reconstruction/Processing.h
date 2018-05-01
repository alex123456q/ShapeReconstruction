#ifndef MOVESKELETW_H
#define MOVESKELETW_H

#include <vector>
#include <stdio.h>
#include <iostream>
#include "SkeletonDemoGUI/SkeletonLib/BSTrans.h"
#include "atlimage.h"
#include <set>
#include <map>
#include <unordered_set>
#include "opencv2/imgproc/imgproc.hpp"

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

enum Borders { Floor, Wall, Roof };

struct Cell{
    double leftx, rightx, upy, downy;
    //vector<TBone> bones; 
    TNode* skeletnode;
    TBone* skeletbone;
    vector<Point> nodes;
	std::vector<bool> paired;
	std::vector<double> color;
    std::/*unordered_*//*set*/vector<std::pair<Point, Point>/*, point_compare*/> borders;
	std::vector<std::pair<double, double> > borders_color;

	std::map<Borders, std::vector<std::pair<Point, Point>>> bords;
	std::map<Borders, std::vector<std::pair<double, double> >> bords_color;

	bool skeletboneex;
    vector<Element*> cellels;
	Cell()
	{
		bords[Floor] = std::vector<std::pair<Point, Point>>();
		bords[Wall] = std::vector<std::pair<Point, Point>>();
		bords[Roof] = std::vector<std::pair<Point, Point>>();

		bords_color[Floor] = std::vector<std::pair<double, double>>();
		bords_color[Wall] = std::vector<std::pair<double, double>>();
		bords_color[Roof] = std::vector<std::pair<double, double>>();
	}
//     bool operator< (const Cell& another){
//         return leftx < another.leftx;
//     }
};
class Reconstruct{
public:
    TPolFigure *skeleton;
	TPolFigure *skeletonVert;
public:
    vector<vector<double> > imageF;
	int firstH, secondH;
private:
    //vector< vector <Cell> > cells; 
//     vector<Cell> cells; 
    CImage image;
public:
	vector<vector<int> > im;
	std::vector<Point*> pointsvert;
	std::vector<int> partpointsvert;
	vector<Cell> cells;
    Reconstruct(CImage im, int col1, int col2);
	Reconstruct();
    ~Reconstruct();
    void makeSkelet(/*CImage im, TPolFigure *skeletonb*/);
    void makeCells(TConnected* Component);
    void addVerticalCells(TConnected* Component, int color1, int color2);

    int mainPart();
	int vertPart(CImage imVert, std::vector<cv::Point>& cvpoints);
	int vertPart(CImage imVert, std::vector<Point>& cvpoints);
    void SetInnerPointsofSkelet(TConnected* Component);
	void SetHeightforBorders(TConnected* Component, std::set<Point>& sPoints, std::set<Point>& sPointsEdge, int firstH, int secondH);
	void SetHeightforBorders2(TConnected* Component, std::set<Point>& sPoints2, int firstH_, int secondH_, double curcolor);
    void findClosestBorder(Cell& curcell, int i, int j, double&x, double&y, double& d1, double& f, std::vector<std::vector<double>>& imageF);
    void findClosestBone(Cell& curcell, int i, int j, double x, double y, double& x2, double& y2, double& d3, double& f, std::vector<std::vector<double>>& imageF);
};