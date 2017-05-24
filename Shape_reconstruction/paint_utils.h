#pragma once;

#include <vector>
#include <stdio.h>
#include <iostream>
#include "SkeletonDemoGUI/SkeletonLib/BSTrans.h"
#include "atlimage.h"
#include "geom_utils.h"

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
            imageF[x][y] =  color;
        }
    }
}

static void PaintSkeletBones(TPolFigure*skeleton,std::vector<std::vector<double>>& imageF){
    //TNode * Node = skeleton->Components->first()->Nodes->first();
    TBone* Bone = skeleton->Components->first()->Bones->first();
    while (Bone){    
        //TBone* Bone = Node->Bones[i];
        TNode* Node = Bone->dest;
        TNode* nextNode = Bone->org;//GetNextNode(Node);
        line m;
        m.a = Node->Y() - nextNode->Y();
        m.b = nextNode->X() - Node->X();
        m.c = -m.a*Node->X() - m.b*Node->Y();
        double d1, d2;
        Point res;
        if (fabs(m.b) > fabs(m.a)){
            for (int x = min(Node->X(), nextNode->X()); x < max(Node->X(), nextNode->X()); ++x){
                int y;
                y = (-m.a*x-m.c)/m.b;
                if (m.b == 0)
                    y = Node->Y();
                res.X = x;
                res.Y = y;
                if (Bone->Virt)
                    FindParabolaPoint(Bone, x, y, res);
                d1 = DistPoint(&res, &Point(Node->Disc->X,Node->Disc->Y));//&bone->dest);
                d2 = DistPoint(&res, &Point(nextNode->Disc->X, nextNode->Disc->Y));//&bone->org);
                imageF[res.X][res.Y] =  Node->f * (d2/(d1+d2)) + nextNode->f * (d1/(d1+d2)) ;                
            }
        } else {
            for (int y = min(Node->Y(), nextNode->Y()); y < max(Node->Y(), nextNode->Y()); ++y){
                int x = (-m.b*y-m.c)/m.a;
                if (m.a == 0)
                    x = Node->X();
                res.X = x;
                res.Y = y;
                if (Bone->Virt)
                    FindParabolaPoint(Bone, x, y, res);
                d1 = DistPoint(&res, &Point(Node->Disc->X,Node->Disc->Y));//&bone->dest);
                d2 = DistPoint(&res, &Point(nextNode->Disc->X, nextNode->Disc->Y));//&bone->org);
                imageF[res.X][res.Y] =  Node->f * (d2/(d1+d2)) + nextNode->f * (d1/(d1+d2)) ;
            }
        }
        Bone = Bone->getNext();
    }
}

static void PaintBorders(TPolFigure*skeleton,std::vector<std::vector<double>>& imageF){
    Point* p = skeleton->Components->first()->Border->ListPoints->first();
    while (p){    
        Point* Node = p;
        Point* nextNode = p->getNextLooped();//GetNextNode(Node);
        PaintLine(imageF, Node, nextNode, 1);
        p = p->getNext();
    }
}

static void PaintInFile(CImage& image, vector<vector<double> >& imageF){
    CImage newimage;
    newimage.Create(image.GetWidth(), image.GetHeight(), 32);
    for (int i = 0; i < image.GetWidth(); ++i)
        for (int j = 0; j < image.GetHeight(); ++j){
            if (imageF[i][j] < -0.5){
                imageF[i][j] = 0;
            }
            int val = imageF[i][j]*255;
            newimage.SetPixel(i, j,
                RGB(val, val, val));
        }
        newimage.Save(_T("D:\\My documents\\Shape_reconstruction\\data\\after_cur_out.png"), Gdiplus::ImageFormatPNG );
        newimage.Destroy();
}