#include <cmath>
#include "Processing.h"
#include "geom_utils.h"
#include "paint_utils.h"
// #include <boost/math/special_functions/round.hpp>

// struct Cell{
//     int leftx, rightx, upy, downy;
//     //vector<TBone> bones; 
//     TNode* skeletnode;
//     TBone* skeletbone;
//     vector<Point> nodes;
//     vector<std::pair<Point, Point>> borders;
//     vector<Element*> cellels;
// 
//     bool operator< (const Cell& another){
//         return leftx < another.x;
//     }
// };

// auto CompCells = [] (const Cell& e1, const Cell& e2)
// {
//     if (e1.rightx <= e2.leftx)
//         return true;
//     if (e1.leftx >= e2.rightx)
//         return false;
//     if (e1.downy >= e2.upy)
//         return true;
//     if (e1.upy <= e2.downy)
//         return false;
//     return (e1.leftx < e2.leftx);
// };

struct CellLessThan
{
    bool operator() (const Cell& e1, const Cell& e2) const
    {
//          if (e1.rightx < e2.leftx)
//              return true;
//          if (e1.leftx > e2.rightx)
//              return false;
// //          if (e1.downy > e2.upy)
// //              return true;
//          if (e1.upy < e2.downy)
//              return true;
        if (e1.leftx == e2.leftx)
            return (e1.downy < e2.downy);
        return (e1.leftx < e2.leftx);
    }
//     bool operator() (const Point& e1, const Cell& e2) const
//     {
//         //         if (e1.rightx <= e2.leftx)
//         //             return true;
//         //         if (e1.leftx >= e2.rightx)
//         //             return false;
//         //         if (e1.downy >= e2.upy)
//         //             return true;
//         //         if (e1.upy <= e2.downy)
//         //             return false;
//         if (e1.X < e2.leftx)
//             return true;
//         if (e1.X > e2.rightx)
//             return false;
// //         if (e1.Y > e2.upy)
// //             return false;
//         return (e1.Y < e2.downy);
// //         return (e1.leftx < e2.leftx);
//     }
//     bool operator() (const Cell& e1, const Point& e2) const
//     {
//         if (e2.X < e1.leftx)
//             return false;
//         if (e2.X > e1.rightx)
//             return true;
//         return (e2.Y > e1.upy);
//     }
} CellComp;


// class Reconstruct{
//     TPolFigure *skeleton;
//     BitRaster *srcimg;
//     vector<vector<double> > imageF;
//     //vector< vector <Cell> > cells; 
//     vector<Cell> cells; 
//     CImage image;
// public:
//     Reconstruct(CImage im);
//     ~Reconstruct();
//     void makeSkelet();
//     void makeCells();
//     int mainPart();
//     void SetInnerPointsofSkelet();
// };

// namespace std {
// 
//     template<class RandomIt, class T>
//     RandomIt binary_locate(RandomIt first, RandomIt last, const T& val, ) {
//         if(val == *first) return first;
//         auto d = std::distance(first, last);  
//         if(d==1) return first;
//         auto center = (first + (d/2));
//         if(val < *center) return binary_locate(first, center, val);
//         return binary_locate(center, last, val);
//     }  
// 
// }
// int pnpoly(int npol, float * xp, float * yp, float x, float y)
// {
//     int c = 0;
//     for (int i = 0, j = npol - 1; i < npol; j = i++) 
//     {
//         if ((
//             (yp[i]<yp[j]) && (yp[i]<=y) && (y<=yp[j]) &&
//             ((yp[j] - yp[i]) * (x - xp[i]) > (xp[j] - xp[i]) * (y - yp[i]))
//             ) || (
//             (yp[i]>yp[j]) && (yp[j]<=y) && (y<=yp[i]) &&
//             ((yp[j] - yp[i]) * (x - xp[i]) < (xp[j] - xp[i]) * (y - yp[i]))
//             ))
//             c = !c;
//     }
//     return c;
// }

bool pointINcell(Cell& cell, Point& p)
{
    bool f = false;
    std::set<Point, one_point_compare> points;
    for (auto iterat = cell.borders.begin(); iterat != cell.borders.end(); ++iterat){
        Point firstp = iterat->first;
        Point secondp = /*cell.borders[i]*/iterat->second;
        Point res;
        if ( intersect(firstp, secondp, p, Point(10, p.Y), res) ){
            if (abs(int(std::floor(res.X + 0.5)) - res.X) < 1e-3)
                res.X = int(std::floor(res.X + 0.5));
            if (abs(int(std::floor(res.Y + 0.5)) - res.Y) < 1e-3)
                res.Y = int(std::floor(res.Y + 0.5));
            if (res.X >= min(secondp.X, firstp.X) && res.X <= max(secondp.X, firstp.X) && 
                res.Y >= min(secondp.Y, firstp.Y) && res.Y <= max(secondp.Y, firstp.Y) &&
                (res.X > p.X) == (10 > p.X) ){
            if (points.find(res) != points.end())
                continue;
            points.insert(res);
//         if (firstp.X ==)
//         if ( Classify(&firstp, &secondp, &p) == LEFT_POS && ((firstp.Y<p.Y)&&(p.Y<=secondp.Y)) ||
//            Classify(&firstp, &secondp, &p) == RIGHT_POS && ((secondp.Y<p.Y)&&(p.Y<=firstp.Y)) )
       // if ( (p.Y <= max(firstp.Y, secondp.Y) && p.Y >= min(firstp.Y, secondp.Y) )
            f = !f;
            }
//         if (   (firstp.X - secondp.X)*(p.Y - secondp.Y) - (firstp.Y - secondp.Y)*(p.X - secondp.X) < 0 )
//             f  = false;
        }
    }
//     for (int i = 0; i < cell.borders.size(); ++i){
//         Point firstp = cell.borders[i].first;
//         Point secondp = cell.borders[i].second;
//         if
//         ((((firstp.Y<=p.Y) && (p.Y<secondp.Y)) || ((secondp.Y<=p.Y) && (p.Y<firstp.Y))) &&
//             (p.X > (secondp.X - firstp.X) * (p.Y - firstp.Y) / (secondp.Y - firstp.Y) + firstp.X))
//             c = !c;
// //         if ((
// //             (firstp.Y<secondp.Y) && (firstp.Y<=p.Y) && (p.Y<=secondp.Y) &&
// //             ((secondp.Y - firstp.Y) * (p.X - firstp.X) >= (secondp.X - firstp.X) * (p.Y - firstp.Y))
// //             ) || (
// //             (firstp.Y>secondp.Y) && (secondp.Y<=p.Y) && (p.Y<=firstp.Y) &&
// //             ((secondp.Y - firstp.Y) * (p.X - firstp.X) <= (secondp.X - firstp.X) * (p.Y - firstp.Y))
// //             ))
// //             c = !c;
//     }
    return f;
}

void FindCell(int x, int y, vector<Cell>& cells, Cell& outcell){
//     Cell value;
//     value.leftx = value.rightx = x;
//     value.upy = value.downy = y;
    int idx = 0;
//     std::pair<std::vector<Cell>::iterator,std::vector<Cell>::iterator> bounds;
//     bounds=std::equal_range (cells.begin(), cells.end(), value/*Point(x,y)*/, CellComp);
//     for (auto iterat = bounds.first; iterat != bounds.second; ++iterat){
//         if (iterat->rightx >= x && iterat->leftx <= x && 
//            iterat->upy >= y && iterat->downy <= y){
//              idx = iterat - cells.begin();
//              break;
//        }
//     }
    //long idx = (std::lower_bound(cells.begin(), cells.end(), Point(x, y), CellComp) - cells.begin());
    //if (idx > 0)
    //    --idx
    while (cells[idx].rightx < x)
        ++idx;
    for (int i = idx; i < cells.size(); ++i)
        if (pointINcell(cells[i], Point(x, y))){//
           // (cells[i].rightx >= x && cells[i].leftx <= x && 
          //  cells[i].upy >= y && cells[i].downy <= y){
            idx = i;
            break;
        }
    outcell = cells[idx];
    //std::binary_search(cells.begin(), cells.end(), value, CompCells());
}

void Reconstruct::makeCells(){
    //TNode * Node = skeleton->Components->first()->Nodes->first();
    TBone* Bone = skeleton->Components->first()->Bones->first();
    while (Bone){
        Cell newcell = Cell();
        newcell.skeletbone = Bone;
        TNode* orgnode = Bone->org;
        TNode* destnode = Bone->dest;
        for (int i = 0; i < 3; ++i){
            if (!orgnode->Sites[i])
                break;
            for (int j = 0; j < 3; ++j){
                if (destnode->Sites[j] == orgnode->Sites[i]){
                    if (destnode->Sites[j]->isVertex){
                        newcell.nodes.push_back(*((Vertex*)orgnode->Sites[i])->p);
                        Point pp1 = *((Vertex*)orgnode->Sites[i])->p;
                        newcell.borders.insert(std::pair<Point, Point>(pp1,Point(orgnode->X(), orgnode->Y())));
                        newcell.borders.insert(std::pair<Point, Point>(pp1,Point(destnode->X(), destnode->Y())));
                        newcell.cellels.push_back(orgnode->Sites[i]);
                        //PaintLine(imageF, &pp1, &Point(orgnode->X(), orgnode->Y()), 1);
                        //PaintLine(imageF, &pp1, &Point(destnode->X(), destnode->Y()), 1);
                    } else {
                        Edge* edge = (Edge*)destnode->Sites[j];
                        Point pp1 = get_perpendicular_pt_from_pt_to_line( *edge->dest, *edge->org, Point(destnode->X(), destnode->Y()));
                        newcell.nodes.push_back(pp1);
//                         if (pp1.X == destnode->X() && pp1.Y == destnode->Y())
//                             break;
                        Point pp2 = get_perpendicular_pt_from_pt_to_line( *edge->dest, *edge->org, Point(orgnode->X(), orgnode->Y()));
//                         if (pp2.X == orgnode->X() && pp2.Y == orgnode->Y())
//                             break;
                        newcell.nodes.push_back(pp2);
                        //if ( edge->WestDirect() )
                        newcell.borders.insert(std::pair<Point, Point>(pp1,pp2));
                        newcell.borders.insert(std::pair<Point, Point>(pp2,Point(orgnode->X(), orgnode->Y())));
                        newcell.borders.insert(std::pair<Point, Point>(pp1,Point(destnode->X(), destnode->Y())));

                        newcell.cellels.push_back(destnode->Sites[j]);
                        //PaintLine(imageF, &pp2, &Point(orgnode->X(), orgnode->Y()), 1);
                        //PaintLine(imageF, &pp1, &Point(destnode->X(), destnode->Y()), 1);
                    }
                    if (newcell.borders.size() > 5)
                        break;
                    break;
                }
            }
        }
        newcell.leftx = min(destnode->X(), orgnode->X());
        newcell.rightx = max(destnode->X(), orgnode->X());
        newcell.downy = min(destnode->Y(), orgnode->Y());
        newcell.upy = max(destnode->Y(), orgnode->Y());
        for (int i = 0; i < newcell.nodes.size(); ++i){
            if (newcell.nodes[i].X < newcell.leftx)
                newcell.leftx = newcell.nodes[i].X;
            if (newcell.nodes[i].X > newcell.rightx)
                newcell.rightx = newcell.nodes[i].X;
            if (newcell.nodes[i].Y < newcell.downy)
                newcell.downy = newcell.nodes[i].Y;
            if (newcell.nodes[i].Y > newcell.upy)
                newcell.upy = newcell.nodes[i].Y;
//            PaintLine(imageF, &newcell.nodes[i], &newcell.nodes[(i+1)%newcell.nodes.size()], 1);
        }

        cells.push_back(newcell); //here we want a vector of vectors and insert with conditions
        Bone = Bone->getNext();
    }
}

void Reconstruct::makeSkelet(){
    srcimg = new BitRaster(image.GetWidth(), image.GetHeight());
    bool inverted = false;
    for (int i = 0; i < image.GetHeight(); i++) {
        for(int j = 0; j < image.GetWidth(); j++) {
            bool isBlack = (GetBValue(image.GetPixel(j, i)) < 128
                || GetRValue(image.GetPixel(j, i)) < 128
                || GetGValue(image.GetPixel(j, i)) < 128);
            if (!inverted) {
                if (isBlack) {
                    srcimg->setBit(j, i, isBlack);
                }
            }
            else {
                if (!isBlack) {
                    srcimg->setBit(j, i, !isBlack);
                }
            }
        }
    }

    BondSkeletTrans(srcimg, 0, 100, skeleton);
    //skeleton->CutSkeleton(1);
    //skeleton->setFakeKind();
    //skeleton->fakeCutSkeleton(1);
    //this->update();
}

Reconstruct::Reconstruct(CImage im)
    : image(im)
{
    skeleton = NULL;
    //image = im;
    makeSkelet();
    imageF.resize(image.GetWidth());
    for (int i = 0; i < image.GetWidth(); ++i){
        imageF[i].resize(image.GetHeight());
        std::fill (imageF[i].begin(),imageF[i].begin()+image.GetHeight(), -1);
    }
    //imageF = std::vector<std::vector<doub>>(image.GetWidth(), std::vector<int>(image.GetHeight()));
}

Reconstruct::~Reconstruct()
{
}

void Reconstruct::SetInnerPointsofSkelet(){//(TPolFigure* skeleton, std::vector<std::vector<double>>& imageF){
    TNode * Node = skeleton->Components->first()->Nodes->first();
    while (Node){
        int sumint = 0;
        int sumnotint = 0;
        for (int i = 0; i < Node->Kind(); ++i){
            if ( Node->Sites[i]->Cont->Internal ){
                sumint++;
            } else {
                sumnotint++;
            }
        }
        int l = 0;
        int r = 1;
        if (sumint != 0){//(sumint != 0 && sumnotint != 0){//(sumint%2 == 1){
            Node->f = (l+r)/2.0;
            imageF[Node->Disc->X][Node->Disc->Y] = Node->f;
        }
        Node = Node->getNext();
    }
}

void SetOuterPointsofSkelet(TNode* Node, std::vector<std::vector<double>>& imageF){
    if (!Node)
        return;
    int count_bones = sizeof(Node->Bones)/sizeof(Node->Bones[0]);
    for (int i = 0; i < count_bones; ++i){
        TBone* Bone = Node->Bones[i];
        if (!Bone || Bone->GetNextNode(Node)->f > -1)
            continue;
        TNode* nextNode = Bone->GetNextNode(Node);
        double d1 = /*CalculateDistance*/DistPoint(&Point(nextNode->Disc->X,nextNode->Disc->Y),
            &Point(Node->Disc->X, Node->Disc->Y));
        double d2 = CalculateDistanceToBorder(nextNode);
        double val = 0;
        if (nextNode->Sites[0]->Cont->Internal)
            val = 1;
        nextNode->f = Node->f * (d2/(d2+d1)) + val/*nearest*/ *(d1/(d1+d2));
        imageF[int(nextNode->X())][int(nextNode->Y())]  = nextNode->f;
        SetOuterPointsofSkelet(nextNode, imageF);
    }
}

void Reconstruct::findClosestBorder(Cell& curcell, int i, int j, double&x, double&y, double& d1, double& f, std::vector<std::vector<double>>& imageF)
{
    d1 = 1000000;
    double curd;
    Point pp;
    for (int k = 0; k < curcell.cellels.size(); ++k){
        if (curcell.cellels[k]->isVertex){
            pp = *((Vertex*)curcell.cellels[k])->p;
            //curd = DistPoint(&Point(i, j), &((Vertex*)curcell.cellels[i]).p);
        } else {
            Edge* edge = (Edge*)curcell.cellels[k];
            pp = get_perpendicular_pt_from_pt_to_line(*edge->dest, *edge->org, Point(i,j));
//             curd = DistPoint(&Point(i, j), &pp);
        }
        curd = DistPoint(&Point(i, j), &pp);
        if (curd < d1){
            d1 = curd;
            f = curcell.cellels[k]->Cont->Internal? 1:0;
            x = pp.X;
            y = pp.Y;
        }
    }
    imageF[i][j] = f;
}
void Reconstruct::findClosestBone(Cell& curcell, int i, int j, double x, double y, double& x2, double& y2, double& d3, double& f, std::vector<std::vector<double>>& imageF)
{
    d3 = 10000000;
    Point res;
    TBone* bone = curcell.skeletbone;
    bool b = intersect(Point(x, y), Point(i, j),
        Point((bone->dest->Disc->X),(bone->dest->Disc->Y)) ,
        Point((bone->org->Disc->X), (bone->org->Disc->Y)) , res);
    if (bone->Virt){
        FindParabolaPoint(bone, i, j, res);
    }
    if (res.X < 0 || res.X > 500 || res.Y<0 || res.Y > 500)
        return;
    d3 = DistPoint(&Point(res.X, res.Y), &Point(i, j));
    x2 = res.X;
    y2 = res.Y;
    if (imageF[int(res.X)][int(res.Y)] > -0.5){
        f = imageF[int(res.X)][int(res.Y)];
    } else {
        double d1, d2;
        d1 = DistPoint(&Point(res.X, res.Y), &Point(bone->dest->Disc->X,bone->dest->Disc->Y));//&bone->dest);
        d2 = DistPoint(&Point(res.X, res.Y), &Point(bone->org->Disc->X, bone->org->Disc->Y));//&bone->org);
        f = bone->dest->f * (d2/(d1+d2)) + bone->org->f * (d1/(d1+d2));
    }
}
    

int Reconstruct::mainPart(){
    //TConnected* Com = skeleton->Components->first();
    SetInnerPointsofSkelet();   //SetInnerPoints(skeleton, imageF);
    //set outer points
    TNode * Node = skeleton->Components->first()->Nodes->first();
    while (Node){
        SetOuterPointsofSkelet(Node, imageF);
        Node = Node->getNext();
    }
    //////////////////////////////////////////////////////////////////////////

    PaintSkeletBones(skeleton, imageF);
    //PaintBorders(skeleton, imageF);
    //PaintInnerBorders(skeleton, imageF);

    makeCells();
    sort(cells.begin(), cells.end(), CellComp); //cellcompare()
    double d1, d2, h1, h2;
    double x, y;
    for (int i = 0; i < image.GetWidth(); i+=1)
        for (int j = 0; j < image.GetHeight(); ++j){
            if (!srcimg->getBit(i, j)){
                imageF[i][j] = 1;
            }
            if (imageF[i][j] > -0.5)
                continue;
            Cell curcell;
            FindCell(i, j, cells, curcell);
            findClosestBorder(curcell, i, j, x, y, d1, h1, imageF);
            if (i == x && y == j){
                imageF[i][j] = h1;
                continue;
            }
            double x2, y2;
            findClosestBone(curcell, i, j, x, y, x2, y2, d2, h2, imageF);
            //PaintLine(imageF, &Point(x, y), &Point(i, j), 1);
            //PaintLine(imageF, &Point(x2, y2), &Point(i, j), 1);

            //d1 = DistPoint(&Point(i, j), &Point(x, y));
            //d2 = DistPoint(&Point(i, j), &Point(x2, y2));
            imageF[i][j] = h1*(d2/(d1+d2)) + h2*(d1/(d1+d2));
        }
        //RePaintSkeletBones(skeleton, imageF);

    PaintInFile(image, imageF);
    return 0;
}
