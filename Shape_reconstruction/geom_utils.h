#pragma once;

#include <vector>
#include <stdio.h>
#include <iostream>
#include "SkeletonDemoGUI/SkeletonLib/BSTrans.h"
#include "atlimage.h" 

static Point get_perpendicular_pt_from_pt_to_line(
    Point& pta,
    Point& ptb,
    Point& pt_from){
        Point pt_to;
//         double b1 = pt_from.X * (pta.X - ptb.X) + pt_from.Y * (pta.Y - ptb.Y);
//         double b2 = pta.X * ptb.Y - pta.Y * ptb.X;
//         pt_to.Y = (pta.X - ptb.X) * (pta.X - ptb.X) + (pta.Y - ptb.Y) * (pta.Y - ptb.Y);
//         double det_k = b1 * (pta.X - ptb.X) - b2 * (pta.Y - ptb.Y);
// 
//         pt_to.X = det_k/pt_to.Y;
//         det_k = (pta.X - ptb.X) * b2 + (pta.Y - ptb.Y) * b1;
//         pt_to.Y = det_k/pt_to.Y;



		double k = ((ptb.Y - pta.Y) * (pt_from.X - pta.X) - (ptb.X - pta.X) * (pt_from.Y - pta.Y)) / (pow((ptb.Y - pta.Y), 2) + pow((ptb.X - pta.X), 2));
		pt_to.X = pt_from.X - k * (ptb.Y - pta.Y);
		pt_to.Y = pt_from.Y + k * (ptb.X - pta.X);

        return pt_to;
}

static struct line {
    double a, b, c;
};

static double det (double a, double b, double c, double d) {
    return a * d - b * c;
}

static bool intersect (Point p1, Point q1, Point p2, Point q2, Point& res) {
    line m, n;
    m.a = p1.Y- q1.Y;
    m.b = q1.X - p1.X;
    m.c = -m.a*p1.X - m.b*p1.Y;
    n.a = p2.Y - q2.Y;
    n.b = q2.X - p2.X;
    n.c = -n.a*p2.X - n.b*p2.Y;
    double zn = -det (m.a, m.b, n.a, n.b);
    if (fabs (zn) < 0.000001)
        return false;
    res = Point( (det (m.c, m.b, n.c, n.b) / zn),  (det (m.a, m.c, n.a, n.c) / zn));
    return true;
} 
static double CalculateDistanceToBorder(TNode* Node){
    return Node->Disc->Rad;
}

static void FindParabolaPoint(TBone* bone, int x, int y, Point& res){
    int p = 0;
    bool fl = false;
    for (p = 0; p < 3/*bone->org->Kind()*/; ++p){
        if (!bone->org->Sites[p]->isVertex)
            continue;
        for (int z = 0; z < 3/*bone->dest->Kind()*/; ++z)
            if (bone->org->Sites[p] == bone->dest->Sites[z]){
                fl = true;
                break;
            }
            if (fl)
                break;
    }
    Vertex* n = (Vertex*)bone->org->Sites[p];
    fl = false;
    for (p = 0; p < 3/*bone->org->Kind()*/; ++p){
        if (bone->org->Sites[p]->isVertex)
            continue;
        for (int z = 0; z < 3/*bone->dest->Kind()*/; ++z)
            if (bone->org->Sites[p] == bone->dest->Sites[z]){
                fl = true;
                break;
            }
            if (fl)
                break;
    }
    double xnew=100000000;
    bool bbad = false;
    Point pp = get_perpendicular_pt_from_pt_to_line(
        Point ( ((Edge*)bone->org->Sites[p])->org->X, ((Edge*)bone->org->Sites[p])->org->Y ), 
        Point ( ((Edge*)bone->org->Sites[p])->dest->X, ((Edge*)bone->org->Sites[p])->dest->Y ),
        *n->p);
    Point pp2 = get_perpendicular_pt_from_pt_to_line(
        Point ( ((Edge*)bone->org->Sites[p])->org->X, ((Edge*)bone->org->Sites[p])->org->Y ), 
        Point ( ((Edge*)bone->org->Sites[p])->dest->X, ((Edge*)bone->org->Sites[p])->dest->Y ),
        Point(x,y));
    double a =  ( pow(pp.X - pp2.X, 2) + pow(pp.Y - pp2.Y, 2) ); 
    double b =  ( pow(pp.X - n->p->X, 2) + pow(pp.Y - n->p->Y, 2) ); 
    xnew = (a + b)/ (2*sqrt(b));//-d1;
    double d1 = DistPoint(&pp, n->p);
    res.X = pp2.X+(-pp.X+n->p->X)*(xnew/d1); 
    res.Y =  pp2.Y+(-pp.Y+n->p->Y)*(xnew/d1);
}