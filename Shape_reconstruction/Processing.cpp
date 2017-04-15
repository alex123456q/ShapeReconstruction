#include "Processing.h"
#include<cmath>

void Processing::renewSkelet(){
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

Processing::Processing(CImage im)
    : image(im)
{
    skeleton = NULL;
    point = true;
    //image = im;
    tmpCorn = NULL;
    renewSkelet();
    imageF.resize(image.GetWidth());
    for (int i = 0; i < image.GetWidth(); ++i){
        imageF[i].resize(image.GetHeight());
        std::fill (imageF[i].begin(),imageF[i].begin()+image.GetHeight(), -1);
    }
    //imageF = std::vector<std::vector<doub>>(image.GetWidth(), std::vector<int>(image.GetHeight()));

    endPoint = NULL;
    curNode = NULL;
    xPr = 0;
    yPr = 0;
}

Processing::~Processing()
{
}

Point get_perpendicular_pt_from_pt_to_line(
       Point& pta,
       Point& ptb,
       Point& pt_from){
    Point pt_to;
    double b1 = pt_from.X * (pta.X - ptb.X) + pt_from.Y * (pta.Y - ptb.Y);
    double b2 = pta.X * ptb.Y - pta.Y * ptb.X;
    pt_to.Y = (pta.X - ptb.X) * (pta.X - ptb.X) + (pta.Y - ptb.Y) * (pta.Y - ptb.Y);
    double det_k = b1 * (pta.X - ptb.X) - b2 * (pta.Y - ptb.Y);

    pt_to.X = det_k/pt_to.Y;
    det_k = (pta.X - ptb.X) * b2 + (pta.Y - ptb.Y) * b1;
    pt_to.Y = det_k/pt_to.Y;
    return pt_to;
}

/*
 * dfs for all bones from chosen bone
 */
void Processing::dfs(TNode* curNode){
    if (curNode == NULL)
        return;
    for (int i = 0; i < curRot.size(); ++i)
        if (curNode == curRot[i])
            return;
    int i = 0;
    //curNode->Depth = -1;
    curRot.push_back(curNode);
    TBone* Bone;
    if (curNode->Kind() == TailNode){
        endPoint = curNode;
        endPointst = curNode->Bones[0]->dest;
        if (endPointst == endPoint)
            endPointst = curNode->Bones[0]->org;
    }
    while (i < curNode->Kind()){
        Bone = curNode->Bones[i];
        dfs(Bone->GetNextNode(curNode));
        ++i;
    }
}

/*
 * Choose the base bone for rotating and dfs
 */


void Processing::changeSkelet(TNode* curNode, double x, double y){
    TBone* Bone = curNode->Bones[0];
    double min = 1000000;
    int mini = 0;
    double tmp0;
    for (int i = 0; i < curNode->Kind(); ++i){
        Bone = curNode->Bones[i];

        Point pt_to = get_perpendicular_pt_from_pt_to_line(
                    Point(Bone->org->Disc->X,
                    Bone->org->Disc->Y),
                    Point(Bone->dest->Disc->X,
                    Bone->dest->Disc->Y),
                    Point(x, y));

        //if codirect
        tmp0 = pow(pt_to.X - x, 2) + pow(pt_to.Y - y, 2);
        /*if (!Codirect(&Point(Bone->org->Disc->X,
                              Bone->org->Disc->Y), &pt_to,
                              &Point(Bone->dest->Disc->X,
                              Bone->dest->Disc->Y),
                              )){*/
        if ( ((Bone->org->Disc->X - Bone->dest->Disc->X) *
           (curNode->Disc->X - x) +
           (Bone->org->Disc->Y - Bone->dest->Disc->Y) *
           (curNode->Disc->Y - y))*
           ((2*(curNode ==  Bone->org)-1)) < 0)
        {
            tmp0 = 1000000;
        }
        /*tmp0 = DistEdge(
                        &Point(Bone->org->Disc->X,
                        Bone->org->Disc->Y),
                        &Point(Bone->dest->Disc->X,
                        Bone->dest->Disc->Y),
                        &Point(x, y));*/
        std::cout<<tmp0<<" ";
        if (min > tmp0){
            min = tmp0;
            mini = i;
        }
    }
    std::cout << "min distance " << min << " for i=" << mini <<std::endl;
    curRot.push_back(curNode);
    tmp = mini;
    dfs(curNode->Bones[mini]->GetNextNode(curNode));
}


double det (double a, double b, double c, double d) {
    return a * d - b * c;
}

struct line {
    double a, b, c;
};

bool intersect (Point p1, Point q1, Point p2, Point q2, Point& res) {
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
    //if (res.x() < 0 || res.y() < 0)
    //    return false;
    return true;
}

void Processing::Circles(double x, double y){
    std::vector<TSite*> a;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            if (curRot[0]->Sites[i] == curRot[1]->Sites[j])
                a.push_back(curRot[0]->Sites[i]);

    for (int j=0; j < a.size(); ++j){
        TSite* i = a[j];

        //Element* i = k->Cont->Elements[i->NEl];
        if (i->IsVertex()){
            Vertex* n = (Vertex*)i;
            circPoint.push_back(Point(n->p->X, n->p->Y));
        }
        else {
            Edge* edge = (Edge*) i;
            Point newPoint = (get_perpendicular_pt_from_pt_to_line(
                                Point(edge->org->X,
                                edge->org->Y),
                                Point(edge->dest->X,
                                edge->dest->Y),
                                Point(curNode->Disc->X, curNode->Disc->Y))
                                );

            circPoint.push_back(newPoint);
        }
    }
}


double CalculateDistanceToBorder(TNode* Node){
    //size = sizeof(Node->Sites)/sizeof(Node->Sites[0]);
    return Node->Disc->Rad;
}

void SkeletonVertexIntrepolation(TNode* Node, std::vector<std::vector<double>>& imageF){
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
        SkeletonVertexIntrepolation(nextNode, imageF);
    }
}

void SetInnerPoints(TPolFigure* skeleton, std::vector<std::vector<double>>& imageF){
    TNode * Node = skeleton->Components->first()->Nodes->first();
    while (Node){
        int sumint = 0;
        for (int i = 0; i < 3; ++i){
            if (!Node->Sites[i])
                break;
            if ( Node->Sites[i]->Cont->Internal )
                sumint++;
        }
        int l = 0;
        int r = 1;
        if (sumint^1 == 1){
            Node->f = (l+r)/2.0;
            imageF[Node->Disc->X][Node->Disc->Y] = Node->f;
        }
        Node = Node->getNext();
    }
}
void PaintSkeletBones(TPolFigure*skeleton,std::vector<std::vector<double>>& imageF){
    //TNode * Node = skeleton->Components->first()->Nodes->first();
    TBone* Bone = skeleton->Components->first()->Bones->first();
    //int count_bones = sizeof(Node->Bones)/sizeof(Node->Bones[0]);
    //for (int i = 0; i < count_bones; ++i){
    while (Bone){    
        //TBone* Bone = Node->Bones[i];
        TNode* Node = Bone->dest;
        TNode* nextNode = Bone->org;//GetNextNode(Node);
        line m;
        m.a = Node->Y() - nextNode->Y();
        m.b = nextNode->X() - Node->X();
        m.c = -m.a*Node->X() - m.b*Node->Y();
        double d1, d2;
        if (fabs(m.b) > fabs(m.a)){
            for (int x = min(Node->X(), nextNode->X()); x < max(Node->X(), nextNode->X()); ++x){
                int y = (-m.a*x-m.c)/m.b;
                if (m.b == 0)
                    y = Node->Y();
                d1 = DistPoint(&Point(x, y), &Point(Node->Disc->X,Node->Disc->Y));//&bone->dest);
                d2 = DistPoint(&Point(x, y), &Point(nextNode->Disc->X, nextNode->Disc->Y));//&bone->org);
                imageF[x][y] =  Node->f * (d2/(d1+d2)) + nextNode->f * (d1/(d1+d2)) ;
            }
        } else {
            for (int y = min(Node->Y(), nextNode->Y()); y < max(Node->Y(), nextNode->Y()); ++y){
                int x = (-m.b*y-m.c)/m.a;
                if (m.a == 0)
                    x = Node->X();
                d1 = DistPoint(&Point(x, y), &Point(Node->Disc->X,Node->Disc->Y));//&bone->dest);
                d2 = DistPoint(&Point(x, y), &Point(nextNode->Disc->X, nextNode->Disc->Y));//&bone->org);
                imageF[x][y] =  Node->f * (d2/(d1+d2)) + nextNode->f * (d1/(d1+d2)) ;
            }
        }
        Bone = Bone->getNext();
    }
}
void findClosestBone(TPolFigure* skeleton, int i, int j, double x, double y, double& x2, double& y2, double& d3, double& f, std::vector<std::vector<double>>& imageF)
{
    TBone* bone = skeleton->Components->first()->Bones->first();
    double curd, d1, d2;
    d3 = 10000000;
    Point res;
    x2= y2 = 0;
    while (bone){
        bool b = intersect(Point(x, y), Point(i, j),
                  Point((bone->dest->Disc->X),(bone->dest->Disc->Y)) ,
                  Point((bone->org->Disc->X), (bone->org->Disc->Y)) , res);
        if (!b || res.Y < 0 || res.Y >= imageF[0].size() || res.X < 0 || res.X >= imageF.size()){
            bone = bone->getNext();
            continue;
        }
        if (Classify(&Point(bone->dest->Disc->X,bone->dest->Disc->Y), 
            &Point(bone->org->Disc->X, bone->org->Disc->Y), &res) != Between
            && Classify(&Point(bone->dest->Disc->X,bone->dest->Disc->Y), 
            &Point(bone->org->Disc->X, bone->org->Disc->Y), &res) != Origin
            && Classify(&Point(bone->dest->Disc->X,bone->dest->Disc->Y), 
            &Point(bone->org->Disc->X, bone->org->Disc->Y), &res) != Destination){
                //std::cout << "bone " << Classify(&Point(bone->dest->Disc->X,bone->dest->Disc->Y), 
               //     &Point(bone->org->Disc->X, bone->org->Disc->Y), &res);
            bone = bone->getNext();
            continue;
        }
        curd = DistPoint(&Point(res.X, res.Y), &Point(i, j));
        if (curd < d3){
            d3 = curd;
            x2 = res.X;
            y2 = res.Y;
            if (imageF[int(res.X)][int(res.Y)] > -0.5){
                f = imageF[int(res.X)][int(res.Y)];
                bone = bone->getNext();
                continue;
            }
            d1 = DistPoint(&Point(res.X, res.Y), &Point(bone->dest->Disc->X,bone->dest->Disc->Y));//&bone->dest);
            d2 = DistPoint(&Point(res.X, res.Y), &Point(bone->org->Disc->X, bone->org->Disc->Y));//&bone->org);
            f = bone->dest->f * (d2/(d1+d2)) + bone->org->f * (d1/(d1+d2));
            imageF[int(res.X)][int(res.Y)] = f;
        }
        bone = bone->getNext();
    }
}
void PaintLine(std::vector<std::vector<double>>& imageF, Point* Node, Point*nextNode, double color);
void PaintBorders(TPolFigure*skeleton,std::vector<std::vector<double>>& imageF){
    Point* p = skeleton->Components->first()->Border->ListPoints->first();
    while (p){    
        Point* Node = p;
        Point* nextNode = p->getNextLooped();//GetNextNode(Node);
        PaintLine(imageF, Node, nextNode, 1);
        p = p->getNext();
    }
}

void PaintLine(std::vector<std::vector<double>>& imageF, Point* Node, Point*nextNode, double color){
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


void PaintInnerBorders(TPolFigure*skeleton,std::vector<std::vector<double>>& imageF){
    Point* p = skeleton->Components->first()->HoleList[0]->ListPoints->first();
    while (p){    
        Point* Node = p;
        Point* nextNode = p->getNextLooped();//GetNextNode(Node);
        PaintLine(imageF, Node, nextNode, 1);
        p = p->getNext();
    }
}


void findClosestBorder(TPolFigure* skeleton, int i, int j, double&x, double&y, double& d1, double& f, std::vector<std::vector<double>>& imageF)
{
    //skeleton->Components->first()->Border->Elements[i]->isVertex;//listelements
    //skeleton->Components->first()->Border->Elements[i]->Cont->ListPoints;//external
    //skeleton->Components->first()->HoleList;//internal
    //curmin2, f2;//l r
    Point* p = skeleton->Components->first()->Border->ListPoints->first();
    double curd;
    d1 = 1000000;
    while (p){
        curd = DistPoint(p, &Point(i, j));
        if (curd < d1){
            d1 = curd;
            x = p->X;
            y = p->Y;
            f = 0;
            imageF[int(x)][int(y)] = 0;
        }
        Point pt_to = get_perpendicular_pt_from_pt_to_line(
                    *p,
                    *(p->getNextLooped()),
                    Point(i, j));
        if (Classify(p, p->getNextLooped(), &pt_to) != Between
            && Classify(p, p->getNextLooped(), &pt_to) != Origin
            && Classify(p, p->getNextLooped(), &pt_to) != Destination){
            p = p->getNext();
            std::cout << Classify(p, p->getNextLooped(), &pt_to);
            continue;
        }
        curd = DistPoint(&pt_to, &Point(i, j));
        if (curd < d1){
            d1 = curd;
            x = pt_to.X;
            y = pt_to.Y;
            f = 0;
            imageF[int(x)][int(y)] = 0;
        }
        p = p->getNext();
    }
    p = skeleton->Components->first()->HoleList[0]->ListPoints->first();
    while (p){
        curd = DistPoint(p, &Point(i, j));
        if (curd < d1){
            d1 = curd;
            x = p->X;
            y = p->Y;
            f = 1;
            imageF[int(x)][int(y)] = 1;
        }
        Point pt_to = get_perpendicular_pt_from_pt_to_line(
            *p,
            *(p->getNextLooped()),
            Point(i, j));
        if ((Classify(p, p->getNextLooped(), &pt_to) != Between
            && Classify(p, p->getNextLooped(), &pt_to) != Origin
            && Classify(p, p->getNextLooped(), &pt_to) != Destination)){
            p = p->getNext();
            continue;
        }

        //if (pt_to)
        curd = DistPoint(&pt_to, &Point(i, j));
        if (curd < d1){
            d1 = curd;
            x = pt_to.X;
            y = pt_to.Y;
            f = 1;
            imageF[int(x)][int(y)] = 1;
        }
        p = p->getNext();
    }
}

void Processing::selectPivot(int px, int py){

    //TConnected* Com = skeleton->Components->first();
    //TBone *minBone = NULL;

    SetInnerPoints(skeleton, imageF);

    TNode * Node = skeleton->Components->first()->Nodes->first();
    while (Node){
        SkeletonVertexIntrepolation(Node, imageF);
        Node = Node->getNext();
    }

    PaintSkeletBones(skeleton, imageF);
    //PaintBorders(skeleton, imageF);
    //PaintInnerBorders(skeleton, imageF);
    double d1, d2, h1, h2;
    double x, y;
    CImage newimage;
    for (int i = 0; i < image.GetWidth(); i+=1)
        for (int j = 0; j < image.GetHeight(); ++j){
            if (!srcimg->getBit(i, j)){
                imageF[i][j] = 1;
            }
            if (imageF[i][j] > -0.5)
                continue;
            findClosestBorder(skeleton, i, j, x, y, d1, h1, imageF);
            if (i == x && y == j){
                imageF[i][j] = h1;
                continue;
            }
            double x2, y2;
            findClosestBone(skeleton, i, j, x, y, x2, y2, d2, h2, imageF);
            //PaintLine(imageF, &Point(x, y), &Point(i, j), 1);
            //PaintLine(imageF, &Point(x2, y2), &Point(i, j), 1);

            //d1 = DistPoint(QPoint(i, j), QPoint(x, y));
            //d2 = DistPoint(QPoint(i, j), QPoint(x2, y2));
            imageF[i][j] = h1*(d2/(d1+d2)) + h2*(d1/(d1+d2));
        }
    newimage.Create(500, 500, 32);
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
