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
    n.c = -m.a*p2.X - m.b*p2.Y;
    double zn = -det (m.a, m.b, n.a, n.b);
    if (fabs (zn) < 0.000001)
        return false;
    res = Point( int(det (m.c, m.b, n.c, n.b) / zn),  int(det (m.a, m.c, n.a, n.c) / zn));
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

void findClosestBone(TPolFigure* skeleton, int i, int j, int x, int y, double& d3, double& f, std::vector<std::vector<double>>& imageF)
{
    TBone* bone = skeleton->Components->first()->Bones->first();
    double curd, d1, d2;
    d3 = 1000000;
    Point res;
    while (bone){
        bool b = intersect(Point(x, y), Point(i, j),
                  Point((bone->dest->Disc->X),(bone->dest->Disc->Y)) ,
                  Point((bone->org->Disc->X), (bone->org->Disc->Y)) , res);
        if (!b){
            bone = bone->getNext();
            continue;
        }
        if (Classify(&Point(bone->dest->Disc->X,bone->dest->Disc->Y), 
            &Point(bone->org->Disc->X, bone->org->Disc->Y), &res) != Between){
            bone = bone->getNext();
            continue;
        }
        if (res.Y < 0 || res.Y >= imageF[0].size() || res.X < 0 || res.X >= imageF.size()){
            bone = bone->getNext();
            continue;
        }
        curd = DistPoint(&Point(res.X, res.Y), &Point(i, j));
        if (curd < d3){
            d3 = curd;
            if (imageF[res.X][res.Y] > -1){
                f = imageF[res.X][res.Y];
                bone = bone->getNext();
                continue;
            }
            d1 = DistPoint(&Point(res.X, res.Y), &Point(bone->dest->Disc->X,bone->dest->Disc->Y));//&bone->dest);
            d2 = DistPoint(&Point(res.X, res.Y), &Point(bone->org->Disc->X, bone->org->Disc->Y));//&bone->org);
            f = bone->dest->f * (d2/(d1+d2)) + bone->org->f * (d1/(d1+d2));
            imageF[res.X][res.Y] = f;
        }
        bone = bone->getNext();
    }
}

void findClosestBorder(TPolFigure* skeleton, int i, int j, int&x, int&y, double& d1, double& f)
{
    //skeleton->Components->first()->Border->Elements[i]->isVertex;//listelements
    //skeleton->Components->first()->Border->Elements[i]->Cont->ListPoints;//external
    //skeleton->Components->first()->HoleList;//internal
    //curmin2, f2;//l r
    Point* p = skeleton->Components->first()->Border->ListPoints->first();
    double curd;
    d1 = 1000000;
    while (p){
        Point pt_to = get_perpendicular_pt_from_pt_to_line(
                    *p,
                    *(p->getNextLooped()),
                    Point(i, j));
        if (Classify(p, p->getNextLooped(), &pt_to) != Between){
            p = p->getNext();
            continue;
        }
        curd = DistPoint(&pt_to, &Point(i, j));
        if (curd < d1){
            d1 = curd;
            x = pt_to.X;
            y = pt_to.Y;
            f = 0;
        }
        p = p->getNext();
    }
    p = skeleton->Components->first()->HoleList[0]->ListPoints->first();
    while (p){
        Point pt_to = get_perpendicular_pt_from_pt_to_line(
            *p,
            *(p->getNextLooped()),
            Point(i, j));
        if (Classify(p, p->getNextLooped(), &pt_to) != Between){
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

    double d1, d2, h1, h2;
    int x, y;
    CImage newimage;
    newimage.Create(500, 500, 32);
    for (int i = 0; i < image.GetWidth(); ++i)
        for (int j = 0; j < image.GetHeight(); ++j)
            if (!srcimg->getBit(i, j))
                newimage.SetPixel(i, j,
                    RGB(255, 255, 255));
            else 
                newimage.SetPixel(i, j,
                RGB(0, 0, 0));
            
    newimage.Save(_T("D:\\My documents\\Shape_reconstruction\\data\\before_cur_out.png"), Gdiplus::ImageFormatPNG);
//for inner points
    for (int i = 0; i < image.GetWidth(); ++i)
        for (int j = 0; j < image.GetHeight(); ++j){
            if (!srcimg->getBit(i, j)){
                newimage.SetPixel(i, j,
                    RGB(255, 255, 255));
                continue;
            }
            if (imageF[i][j] >= 0){
                newimage.SetPixel(i, j, RGB(imageF[i][j], imageF[i][j], imageF[i][j]));
                continue;
            } else {
                //newimage.SetPixel(i, j, RGB(255, 255, 255));
                //continue;
            }
            findClosestBorder(skeleton, i, j, x, y, d1, h1);
            findClosestBone(skeleton, i, j, x, y, d2, h2, imageF);
            //d1 = DistPoint(QPoint(i, j), QPoint(x, y));
            //d2 = DistPoint(QPoint(i, j), QPoint(x2, y2));
            int val = 255*(h1*(d2/(d1+d2)) + h2*(d1/(d1+d2)));
            //tmp
            newimage.SetPixel(i, j,
                              RGB(val, val, val));
            //image.pixel(i, j) = /*im[i,j] =*/
        }
//    for (int i = 0; i < image.width(); ++i)
//        for (int j = image.height()/2; j < image.height(); ++j){
//            newimage.setPixel(QPoint(i, j), qRgb(0, 0, 0));

  //      }
    newimage.Save(_T("D:\\My documents\\Shape_reconstruction\\data\\after_cur_out.png"), Gdiplus::ImageFormatPNG );
    newimage.Destroy();

    //skeleton->Boundary->
    //skeleton->Components->first()->Border->Elements[i]->isVertex;//listelements
    //skeleton->Components->first()->Border->Elements[i]->Cont->ListPoints;//external
    //skeleton->Components->first()->HoleList;//internal
    //curmin2, f2;//l r
    //for bones
    //intersect(perp_on_border, current_point, bone.dist, bone.origin, res);
    //curmin3 = DistPoint(res,current_point);
    //f3 = current_point.f; // h1 * (d2/(d1+d2)) + h2 * (d1/(d1+d2)) hh - dist org
    //f =  h1 * (d2/(d1+d2)) + h2 * (d1/(d1+d2));
}


/*
void MoveSkeletW::paintEvent(QPaintEvent *)
{
    QPainter painter(this);

//    QImage mainImage(500, 300, 32);
//    mainImage.fill(qRgb(255, 255, 255));
//    painter.drawImage(0, 0, mainImage);
    path = QPainterPath();
    path.setFillRule(Qt::WindingFill);
    //painter.fillRect(0, 0, image.width(), image.height(), Qt::white);
    painter.drawImage(0, 0, image);
    //repaint( image.hasAlphaBuffer() );
    //QLabel myLabel;
    //myLabel.setPixmap(QPixmap::fromImage(myImage));
    //myLabel.show();
    image.save("cur_out.png");
    painter.setRenderHint(QPainter::Antialiasing);

    if (skeleton == NULL) {
        return;
    }
    //calloc()

    //memcpy(skeletonPaint, skeleton, sizeof(skeleton));
    //skeletonPaint->CutSkeleton(1);
    //skeleton->CutSkeleton(1);

    TPolFigure *skeletonPaint(skeleton);// = (TPolFigure *)malloc(sizeof(skeleton));//->copy();
    QPolygonF contour;
        painter.setPen(Qt::NoPen);//QPen(Qt::blue, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

        TContour* S = skeletonPaint->Boundary->first();

        while (S != NULL) {
            int cornersCount = S->ListPoints->cardinal();
            TPoint** points = new TPoint*[cornersCount];
            int i = 0;
            TPoint* Corn = S->ListPoints->first();
            while (Corn != NULL)
            {
                points[i++] = Corn;
                Corn = Corn->getNext();
            }
            for (int j = 0; j < cornersCount - 1; j++) {
                contour << QPointF(points[j]->X, points[j]->Y) << QPointF(points[j + 1]->X, points[j + 1]->Y);
         //c       painter.drawLine(points[j]->X, points[j]->Y, points[j + 1]->X, points[j + 1]->Y);
       //         painter.setPen(QPen(Qt::blue, 4.5, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
      //          painter.drawPoint(points[j]->X, points[j]->Y); ///
         //c       painter.setPen(QPen(Qt::blue, 1, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

            }
            contour << QPointF(points[cornersCount - 1]->X, points[cornersCount - 1]->Y) << QPointF(points[0]->X, points[0]->Y);

        //c    painter.drawLine(points[cornersCount - 1]->X, points[cornersCount - 1]->Y, points[0]->X, points[0]->Y);

            delete points;
            S = S->getNext();
        }
    path.addPolygon(contour);
  //c
    painter.setPen(QPen(Qt::blue, 0.5, Qt::SolidLine, Qt::RoundCap,
                        Qt::MiterJoin));

   // painter.setPen(QPen(Qt::green, 1, Qt::DotLine, Qt::RoundCap,
  //                      Qt::MiterJoin));
    //painter.setBrush(Qt::NoBrush);
    painter.setBrush(QBrush(Qt::green, Qt::SolidPattern)); //Qt::TransparentMode Qt::SolidPattern Qt::NoBrush  Dense1Pattern
    painter.drawPath(path);
    if (!curRot.empty()){
 //c
 //       path.addEllipse(QPointF(curRot[0]->Disc->X,
 //                               curRot[0]->Disc->Y),
 //                               curRot[0]->Disc->Rad, curRot[0]->Disc->Rad);
       painter.drawEllipse(QPointF(curRot[0]->Disc->X,
                           curRot[0]->Disc->Y),
               curRot[0]->Disc->Rad, curRot[0]->Disc->Rad);
        painter.drawEllipse(QPointF(curRot[1]->Disc->X,
                            curRot[1]->Disc->Y),
                curRot[1]->Disc->Rad, curRot[1]->Disc->Rad);
   }
   if (tmpCorn){
        painter.setPen(QPen(Qt::magenta, 10, Qt::SolidLine, Qt::RoundCap,
                            Qt::MiterJoin));
        painter.drawPoint(tmpCorn->X, tmpCorn->Y);
    }
   // рисуем скелеты
   TConnected* Com = skeletonPaint->Components->first();
   bool drawBones = 1;
   bool drawCircles = 0;
   while (Com != NULL) {
       painter.setPen(QPen(Qt::red, 1.5, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

       if (drawBones) {
           TBone* Bone = Com->Bones->first();
           while (Bone != NULL) {
               if (Bone->fake != 1){
                   painter.setPen(QPen(Qt::red, 1.5, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
               painter.drawLine(
                   Bone->org->Disc->X,
                   Bone->org->Disc->Y,
                   Bone->dest->Disc->X,
                   Bone->dest->Disc->Y
               );} else {
                   //painter.setPen(QPen(Qt::blue, 1.5, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
                  // painter.drawLine(
                  //     Bone->org->Disc->X,
                  //     Bone->org->Disc->Y,
                 //      Bone->dest->Disc->X,
                 //      Bone->dest->Disc->Y
                 //  );
               }

               Bone = Bone->getNext();
           }
       }

       painter.setPen(QPen(Qt::green, 1, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

       if (drawCircles) {
           TNode* Node = Com->Nodes->first();

           while (Node != NULL) {
               painter.drawEllipse(Node->X() + 1 - Node->r(),
                   Node->Y() + 1 - Node->r(),
                   2 * Node->r(),
                   2 * Node->r());

               Node = Node->getNext();
           }
       }

       Com = Com->getNext();
   }
   painter.setPen(QPen(Qt::red, 8, Qt::SolidLine, Qt::RoundCap,
                       Qt::MiterJoin));
   if (curNode){
//      painter.drawPoint(curNode->Disc->X, curNode->Disc->Y);
  //    painter.drawPoint(curNode->Bones[tmp]->GetNextNode(curNode)->Disc->X,
  //                      curNode->Bones[tmp]->GetNextNode(curNode)->Disc->Y);
   }
//   for (int i = 0; i < vertices.size(); ++i){
//       painter.drawPoint(vertices[i]->Disc->X, vertices[i]->Disc->Y);
       //painter.drawLine(vertices[i]->Disc->X, vertices[i]->Disc->Y, );
//   }
   //painter.setPen(QPen(Qt::black, 8, Qt::SolidLine, Qt::RoundCap,
    //                   Qt::MiterJoin));
   //for (int i = 0; i < circPoint.size(); ++i){
   //   painter.drawPoint(circPoint[i].X, circPoint[i].Y);
       //painter.drawLine(circPoint[i].X, circPoint[i].Y
   //}
   //painter.end();
   painter.setPen(QPen(Qt::magenta, 10, Qt::SolidLine, Qt::RoundCap,
                       Qt::MiterJoin));
   if (endPoint){
       painter.drawPoint(endPoint->Disc->X, endPoint->Disc->Y);
       painter.setPen(QPen(Qt::black, 8, Qt::SolidLine, Qt::RoundCap,
                           Qt::MiterJoin));
       painter.drawPoint(endPointst->Disc->X, endPointst->Disc->Y);
       painter.drawPoint(tmpCorn->X, tmpCorn->Y);
   }
}
*/
