#pragma once;

#include <vector>
#include <set>
//#include <pair>
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


// namespace intersectintervals
// {
// 	const double EPS = 1E-9;
// 
// 	struct pt {
// 		double x, y;
// 	};
// 
// 	struct seg {
// 		pt p, q;
// 		int id;
// 
// 		double get_y(double x) const {
// 			if (abs(p.x - q.x) < EPS)  return p.y;
// 			return p.y + (q.y - p.y) * (x - p.x) / (q.x - p.x);
// 		}
// 	};
// 
// 
// 	inline bool intersect1d(double l1, double r1, double l2, double r2) {
// 		if (l1 > r1)  swap(l1, r1);
// 		if (l2 > r2)  swap(l2, r2);
// 		return max(l1, l2) <= min(r1, r2) + EPS;
// 	}
// 
// 	inline int vec(const pt & a, const pt & b, const pt & c) {
// 		double s = (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
// 		return abs(s) < EPS ? 0 : s > 0 ? +1 : -1;
// 	}
// 
// 	bool intersect(const seg & a, const seg & b) {
// 		return intersect1d(a.p.x, a.q.x, b.p.x, b.q.x)
// 			&& intersect1d(a.p.y, a.q.y, b.p.y, b.q.y)
// 			&& vec(a.p, a.q, b.p) * vec(a.p, a.q, b.q) <= 0
// 			&& vec(b.p, b.q, a.p) * vec(b.p, b.q, a.q) <= 0;
// 	}
// 
// 
// 	bool operator< (const seg & a, const seg & b) {
// 		double x = max(min(a.p.x, a.q.x), min(b.p.x, b.q.x));
// 		return a.get_y(x) < b.get_y(x) - EPS;
// 	}
// 
// 
// 	struct event {
// 		double x;
// 		int tp, id;
// 
// 		event() { }
// 		event(double x, int tp, int id)
// 			: x(x), tp(tp), id(id)
// 		{ }
// 
// 		bool operator< (const event & e) const {
// 			if (abs(x - e.x) > EPS)  return x < e.x;
// 			return tp > e.tp;
// 		}
// 	};
// 
// 	set<seg> s;
// 	vector < std::set<seg>::iterator > where;
// 
// 	inline std::set<seg>::iterator prev(std::set<seg>::iterator it) {
// 		return it == s.begin() ? s.end() : --it;
// 	}
// 
// 	inline std::set<seg>::iterator next(std::set<seg>::iterator it) {
// 		return ++it;
// 	}
// 
// 	std::vector<pair<int, int>> solve(const std::vector<seg> & a) {
// 		std::vector<pair<int, int>> outvector;
// 		int n = (int)a.size();
// 		vector<event> e;
// 		for (int i = 0; i < n; ++i) {
// 			e.push_back(event(min(a[i].p.x, a[i].q.x), +1, i));
// 			e.push_back(event(max(a[i].p.x, a[i].q.x), -1, i));
// 		}
// 		sort(e.begin(), e.end());
// 
// 		s.clear();
// 		where.resize(a.size());
// 		for (size_t i = 0; i < e.size(); ++i) {
// 			int id = e[i].id;
// 			if (e[i].tp == +1) {
// 				set<seg>::iterator
// 					nxt = s.lower_bound(a[id]),
// 					prv = prev(nxt);
// 				if (nxt != s.end() && intersect(*nxt, a[id]))
// 					outvector.push_back(make_pair(nxt->id, id)); //return 
// 				if (prv != s.end() && intersect(*prv, a[id]))
// 					outvector.push_back(make_pair(prv->id, id));//return make_pair(prv->id, id);
// 				where[id] = s.insert(nxt, a[id]);
// 			}
// 			else {
// 				set<seg>::iterator
// 					nxt = next(where[id]),
// 					prv = prev(where[id]);
// 				if (nxt != s.end() && prv != s.end() && intersect(*nxt, *prv))
// 					outvector.push_back(make_pair(prv->id, nxt->id));//return make_pair(prv->id, nxt->id);
// 				s.erase(where[id]);
// 			}
// 		}
// 
// 		return outvector;//make_pair(-1, -1);
// 	}
// }