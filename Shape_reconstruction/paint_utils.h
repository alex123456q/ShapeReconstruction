#pragma once;

#include <vector>
#include <stdio.h>
#include <iostream>
#include "SkeletonDemoGUI/SkeletonLib/BSTrans.h"
#include "atlimage.h"
#include "geom_utils.h"
#include <set>
struct Cell;
//#include <QtSvg>

static void PaintLine(std::vector<std::vector<double>>& imageF, Point* Node, Point*nextNode, double color){
    line m;
    m.a = Node->Y - nextNode->Y;
    m.b = nextNode->X - Node->X;
    m.c = -m.a*Node->X - m.b*Node->Y;
    double d1, d2;
    if (fabs(m.b) > fabs(m.a)){
        for (int x = min(Node->X, nextNode->X); x < max(Node->X, nextNode->X); ++x){
            int y = (-m.a*x-m.c)/m.b;
            if (m.b == 0)
                y = Node->Y;
            imageF[x][y] =  color ;
        }
    } else {
        for (int y = min(Node->Y, nextNode->Y); y < max(Node->Y, nextNode->Y); ++y){
            int x = (-m.b*y-m.c)/m.a;
            if (m.a == 0)
                x = Node->X;
            imageF[x][y] = color;
        }
    }
	//imageF[min(Node->X, nextNode->X)][min(Node->Y, nextNode->Y)] = color;
}

static void PaintSkeletBones(TPolFigure*skeleton,std::vector<std::vector<double>>& imageF, bool onlywhite){
    //TNode * Node = skeleton->Components->first()->Nodes->first();
	TConnected* Comp = skeleton->Components->first();
	while (Comp)
	{
		TBone* Bone = Comp->Bones->first();
		while (Bone) {
			//TBone* Bone = Node->Bones[i];
			TNode* Node = Bone->dest;
			TNode* nextNode = Bone->org;//GetNextNode(Node);
			line m;
			m.a = Node->Y() - nextNode->Y();
			m.b = nextNode->X() - Node->X();
			m.c = -m.a*Node->X() - m.b*Node->Y();
			double d1, d2;
			Point res;
			if (fabs(m.b) > fabs(m.a)) {
				for (int x = min(Node->X(), nextNode->X()); x < max(Node->X(), nextNode->X()); ++x) {
					int y;
					y = (-m.a*x - m.c) / m.b;
					if (m.b == 0)
						y = Node->Y();
					res.X = x;
					res.Y = y;
					if (Bone->Virt)
						FindParabolaPoint(Bone, x, y, res);
					d1 = DistPoint(&res, &Point(Node->Disc->X, Node->Disc->Y));//&bone->dest);
					d2 = DistPoint(&res, &Point(nextNode->Disc->X, nextNode->Disc->Y));//&bone->org);
					if (onlywhite)
						imageF[res.X][res.Y] = 1;// Node->f * (d2/(d1+d2)) + nextNode->f * (d1/(d1+d2)) ;    
					else
						imageF[res.X][res.Y] = Node->f * (d2 / (d1 + d2)) + nextNode->f * (d1 / (d1 + d2));
				}
			}
			else {
				for (int y = min(Node->Y(), nextNode->Y()); y < max(Node->Y(), nextNode->Y()); ++y) {
					int x = (-m.b*y - m.c) / m.a;
					if (m.a == 0)
						x = Node->X();
					res.X = x;
					res.Y = y;
					if (Bone->Virt)
						FindParabolaPoint(Bone, x, y, res);
					d1 = DistPoint(&res, &Point(Node->Disc->X, Node->Disc->Y));//&bone->dest);
					d2 = DistPoint(&res, &Point(nextNode->Disc->X, nextNode->Disc->Y));//&bone->org);
					if (onlywhite)
						imageF[res.X][res.Y] = 1;
					else
						imageF[res.X][res.Y] = Node->f * (d2 / (d1 + d2)) + nextNode->f * (d1 / (d1 + d2));
				}
			}
			//imageF[Node->X()][Node->Y()] = 0.5;
			//imageF[nextNode->X()][nextNode->Y()] = 0.5;
			Bone = Bone->getNext();
		}
		Comp = Comp->getNext();
	}
}

static void PaintBorders(TPolFigure*skeleton,std::vector<std::vector<double>>& imageF, double color = 1.0){
	TConnected* comp = skeleton->Components->first();
	while (comp) {
		Point* p = /*skeleton->Components->first()*/comp->Border->ListPoints->first();
		while (p) {
			Point* Node = p;
			Point* nextNode = p->getNextLooped();//GetNextNode(Node);
			PaintLine(imageF, Node, nextNode, color);
			p = p->getNext();
		}
		comp = comp->getNext();
	}
}

static void PaintBorders2(TPolFigure*skeleton,std::vector<std::vector<double>>& imageF, std::set<Point>& z){
    Point* p = skeleton->Components->first()->Border->ListPoints->first();
    while (p){    
        Point* Node = p;
        Point* nextNode = p->getNextLooped();//GetNextNode(Node);
        if (z.find(*Node)!=z.end() || z.find(*nextNode)!=z.end())
            PaintLine(imageF, Node, nextNode, 0.5);
        else
            PaintLine(imageF, Node, nextNode, 1);
        p = p->getNext();
    }
}

static void PaintInFile(/*CImage& image, */vector<vector<double> >& imageF, /*std::wstring&*/const char* filename){
    CImage newimage;
    newimage.Create(imageF.size(), imageF[0].size(), 32);
    for (int i = 0; i < imageF.size()/*image.GetWidth()*/; ++i)
        for (int j = 0; j < imageF[0].size()/*image.GetHeight()*/; ++j){
            if (imageF[i][j] < -0.5){
                imageF[i][j] = 0;
            }
            int val = imageF[i][j]*255;
            newimage.SetPixel(i, j,
                RGB(val, val, val));
        }
        newimage.Save(filename, Gdiplus::ImageFormatPNG );
        newimage.Destroy();
}

int main2(int argc, char** argv,/* std::vector<std::vector<double>> imageF,*/ std::vector<std::vector<Cell>>& cvec);