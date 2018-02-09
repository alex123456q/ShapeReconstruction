#include <cmath>
#include "Processing.h"
#include "geom_utils.h"
#include "paint_utils.h"

struct CellLessThan
{
    bool operator() (const Cell& e1, const Cell& e2) const
    {
        if (e1.rightx == e2.rightx)
            return (e1.downy < e2.downy);
        return (e1.rightx < e2.rightx);  
    
    }
} CellComp;

bool pointINcell(Cell& cell, Point& p)
{
    bool f = false;
    if (cell.leftx > p.X || cell.upy < p.Y || cell.downy > p.Y)
        return false;
    std::set<Point, one_point_compare> points;
    for (auto iterat = cell.borders.begin(); iterat != cell.borders.end(); ++iterat){
        Point firstp = iterat->first;
        Point secondp = /*cell.borders[i]*/iterat->second;

		if (firstp.Y == secondp.Y)
			continue;
        Point res;
        if ( intersect(firstp, secondp, p, Point(0, p.Y), res) ){
            if (abs(int(std::floor(res.X + 0.5)) - res.X) < 1e-3)
                res.X = int(std::floor(res.X + 0.5));
            if (abs(int(std::floor(res.Y + 0.5)) - res.Y) < 1e-3)
                res.Y = int(std::floor(res.Y + 0.5));
            if ((res.X - min(secondp.X, firstp.X)) > -1e-5 && (max(secondp.X, firstp.X)-res.X) > -1e-5 &&
				(res.Y - min(secondp.Y, firstp.Y)) > -1e-5 && (max(secondp.Y, firstp.Y) - res.Y) > -1e-5 &&//res.Y >= min(secondp.Y, firstp.Y) && res.Y <= max(secondp.Y, firstp.Y) &&
                (res.X > p.X) == (0 > p.X) ){
            //if (points.find(res) != points.end())
            //    continue;
			if (abs(res.X - secondp.X) < 1e-5 && abs(res.Y - secondp.Y) < 1e-5 && secondp.Y > firstp.Y ||
				abs(res.X - firstp.X) < 1e-5 && abs(res.Y - firstp.Y) < 1e-5 && secondp.Y < firstp.Y){
				continue;
				}

            f = !f;
            }
        }
    }

    return f;
}

bool FindCell(int x, int y, vector<Cell>& cells, Cell& outcell){
    int idx = 0;

    while (cells[idx].rightx < x)
        ++idx;
	for (int i = idx; i < cells.size(); ++i)
        if (pointINcell(cells[i], Point(x, y))){
            idx = i;
            break;
		}
		else idx = -1;
	if (idx == -1)
		return false;
    outcell = cells[idx];
	return true;
}

void Reconstruct::addVerticalCells(TConnected* Component, int color1, int color2){
    TBone* Bone = Component->Bones->first();
    while (Bone){
        Cell newcell = Cell();
        newcell.skeletbone = Bone;
        TNode* orgnode = Bone->org;
        TNode* destnode = Bone->dest;

		destnode->f = orgnode->f = (color1 + color2) / 2;
		
		bool found = false;
		std::vector<int> bothsites;
        for (int i = 0; i < 3; ++i){
            if (!orgnode->Sites[i])
                break;
            for (int j = 0; j < 3; ++j){
				if (destnode->Sites[j] == orgnode->Sites[i]) {
					bothsites.push_back(j);
					found = true;
					break;
				}
// 				if (newcell.borders.size() > 5)
// 					break;
            }
        }
		int color = color1;
		for (int i = 0; i < bothsites.size(); ++i) {
			if (orgnode->Sites[bothsites[i]]->isVertex) {
				newcell.nodes.push_back(*((Vertex*)orgnode->Sites[bothsites[i]])->p);

				Point pp1 = *((Vertex*)orgnode->Sites[bothsites[i]])->p;
				newcell.borders.push_back(std::pair<Point, Point>(pp1, Point(orgnode->X(), orgnode->Y())));
				newcell.borders_color.push_back(std::pair<double, double>(color, orgnode->f));

				newcell.borders.push_back(std::pair<Point, Point>(pp1, Point(destnode->X(), destnode->Y())));
				newcell.borders_color.push_back(std::pair<double, double>(color, destnode->f));

				newcell.cellels.push_back(orgnode->Sites[i]);
			}
			else {
				Edge* edge = (Edge*)orgnode->Sites[bothsites[i]];
				Point pp1 = Point(destnode->X(), destnode->Y());//get_perpendicular_pt_from_pt_to_line( *edge->dest, *edge->org, Point(destnode->X(), destnode->Y()));		
				Point pp2 = Point(orgnode->X(), orgnode->Y());//get_perpendicular_pt_from_pt_to_line( *edge->dest, *edge->org, Point(orgnode->X(), orgnode->Y()));

				newcell.nodes.push_back(pp1);
				newcell.nodes.push_back(pp2);

				newcell.borders.push_back(std::pair<Point, Point>(pp1, pp2));
				newcell.borders.push_back(std::pair<Point, Point>(pp2, Point(orgnode->X(), orgnode->Y())));
				newcell.borders.push_back(std::pair<Point, Point>(pp1, Point(destnode->X(), destnode->Y())));


				newcell.borders_color.push_back(std::pair<double, double>(color, color));
				newcell.borders_color.push_back(std::pair<double, double>(color, orgnode->f));
				newcell.borders_color.push_back(std::pair<double, double>(color, destnode->f));

				newcell.cellels.push_back(orgnode->Sites[bothsites[i]]);
			}
			color = color2;
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
            //PaintLine(imageF, &newcell.nodes[i], &newcell.nodes[(i+1)%newcell.nodes.size()], 1);
        }

        cells.push_back(newcell); //here we want a vector of vectors and insert with conditions
        Bone = Bone->getNext();
		//++z;
		//if (z == 64)18
//		break;
    }
}

void Reconstruct::makeCells(TConnected* Component) {
	TBone* Bone = Component->Bones->first();
	while (Bone) {
		Cell newcell = Cell();
		newcell.skeletbone = Bone;
		TNode* orgnode = Bone->org;
		TNode* destnode = Bone->dest;
		bool found = false;
		for (int i = 0; i < 3; ++i) {
			if (!orgnode->Sites[i])
				break;
			if (orgnode->f == 0)
				std::cout << "0";
			for (int j = 0; j < 3; ++j) {
				if (destnode->Sites[j] == orgnode->Sites[i]) {
					found = true;
					if (destnode->Sites[j]->isVertex) {						
						Point pp1 = *((Vertex*)orgnode->Sites[i])->p;


						/*if (orgnode->Sites[i]->f == (firstH+secondH)/2.0)
						{
							std::cout << "ovk " << orgnode->f << std::endl;

							newcell.borders.push_back(std::pair<Point, Point>(pp1, Point(orgnode->X(), orgnode->Y())));
							newcell.borders_color.push_back(std::pair<double, double>(firstH, orgnode->f));

							newcell.borders.push_back(std::pair<Point, Point>(pp1, Point(destnode->X(), destnode->Y())));
							newcell.borders_color.push_back(std::pair<double, double>(firstH, destnode->f));

							newcell.borders.push_back(std::pair<Point, Point>(pp1, Point(orgnode->X(), orgnode->Y())));
							newcell.borders_color.push_back(std::pair<double, double>(secondH, orgnode->f));

							newcell.borders.push_back(std::pair<Point, Point>(pp1, Point(destnode->X(), destnode->Y())));
							newcell.borders_color.push_back(std::pair<double, double>(secondH, destnode->f));


							break;
						}*/

						newcell.nodes.push_back(pp1);
						newcell.borders./*insert*/push_back(std::pair<Point, Point>(pp1, Point(orgnode->X(), orgnode->Y())));
						newcell.borders_color.push_back(std::pair<double, double>(orgnode->Sites[i]->f, orgnode->f));
						newcell.borders.push_back/*insert*/(std::pair<Point, Point>(pp1, Point(destnode->X(), destnode->Y())));
						//newcell.borders_color.push_back(orgnode->Sites[i]->f);
						newcell.borders_color.push_back(std::pair<double, double>(orgnode->Sites[i]->f, destnode->f));

						// 						newcell.bordersNodes.insert(std::pair<TSite*, TNode*>(orgnode->Sites[i], orgnode));
						// 						newcell.bordersNodes.insert(std::pair<TSite*, TNode*>(orgnode->Sites[i], orgnode));

						newcell.cellels.push_back(orgnode->Sites[i]);
						//                         PaintLine(imageF, &pp1, &Point(orgnode->X(), orgnode->Y()), 1);
						//                         PaintLine(imageF, &pp1, &Point(destnode->X(), destnode->Y()), 1);
					}
					else {
						Edge* edge = (Edge*)destnode->Sites[j];
						
						Point pp1 = get_perpendicular_pt_from_pt_to_line(*edge->dest, *edge->org, Point(destnode->X(), destnode->Y()));
						//                         if (pp1.X == destnode->X() && pp1.Y == destnode->Y())
						//                             break;

						Point pp2 = get_perpendicular_pt_from_pt_to_line(*edge->dest, *edge->org, Point(orgnode->X(), orgnode->Y()));

						newcell.nodes.push_back(pp1);
						newcell.nodes.push_back(pp2);
						//if ( edge->WestDirect() )

						newcell.borders.push_back/*insert*/(std::pair<Point, Point>(pp1, pp2));
						newcell.borders.push_back/*insert*/(std::pair<Point, Point>(pp2, Point(orgnode->X(), orgnode->Y())));
						newcell.borders.push_back(std::pair<Point, Point>(pp1, Point(destnode->X(), destnode->Y())));



						/*if (edge->f == (firstH+secondH)/2)
						{
							//std::cout << "ok";
							
							if (edge->dest->X == destnode->X() && edge->dest->Y == destnode->Y() || edge->dest->X == orgnode->X() && edge->dest->Y == orgnode->Y())
							{
								std::cout << "yed";
								destnode->f = 0.5;
								orgnode->f = 0.5; //one of them
								edge->f = edge->getPrevLooped()->f;// orgnode->f;
							}

							if (edge->org->X == orgnode->X() && edge->org->Y == orgnode->Y() || edge->org->X == destnode->X() && edge->org->Y == destnode->Y())
							{
								std::cout << "yeo";
								orgnode->f = 0.5;
								destnode->f = 0.5;
								edge->f = edge->getNextLooped()->f;
							}
							newcell.borders_color.push_back(std::pair<double, double>(edge->f, edge->f));
							newcell.borders_color.push_back(std::pair<double, double>(edge->f, orgnode->f));
							newcell.borders_color.push_back(std::pair<double, double>(edge->f, destnode->f));
						}
						else*/
						{
							newcell.borders_color.push_back(std::pair<double, double>(edge->f, edge->f));
							newcell.borders_color.push_back(std::pair<double, double>(edge->f, orgnode->f));
							newcell.borders_color.push_back(std::pair<double, double>(edge->f, destnode->f));
						}
						


						//newcell.bordersNodes.insert(std::pair<TSite*, TNode*>((TSite*)destnode, orgnode));
						//newcell.bordersNodes.insert(std::pair<TSite*, TNode*>(destnode->Sites[j], orgnode));
						//newcell.bordersNodes.insert(std::pair<TSite*, TNode*>(destnode->Sites[j], destnode));

						newcell.cellels.push_back(destnode->Sites[j]);
						// 						PaintLine(imageF, &pp1, &pp2, 1);
						// 						PaintLine(imageF,  edge->dest, edge->org, 1);
						// 						PaintLine(imageF, &pp2, &Point(orgnode->X(), orgnode->Y()), 1);
						// 						PaintLine(imageF, &pp1, &Point(destnode->X(), destnode->Y()), 1);
					}
					break;
				}
				// 				if (newcell.borders.size() > 5)
				// 					break;
			}
			if (!found && orgnode->Sites[i]->f == 1e+3)
				std::cout << "n";
		}
		newcell.leftx = min(destnode->X(), orgnode->X());
		newcell.rightx = max(destnode->X(), orgnode->X());
		newcell.downy = min(destnode->Y(), orgnode->Y());
		newcell.upy = max(destnode->Y(), orgnode->Y());
		for (int i = 0; i < newcell.nodes.size(); ++i) {
			if (newcell.nodes[i].X < newcell.leftx)
				newcell.leftx = newcell.nodes[i].X;
			if (newcell.nodes[i].X > newcell.rightx)
				newcell.rightx = newcell.nodes[i].X;
			if (newcell.nodes[i].Y < newcell.downy)
				newcell.downy = newcell.nodes[i].Y;
			if (newcell.nodes[i].Y > newcell.upy)
				newcell.upy = newcell.nodes[i].Y;
			//PaintLine(imageF, &newcell.nodes[i], &newcell.nodes[(i+1)%newcell.nodes.size()], 1);
		}

		cells.push_back(newcell); //here we want a vector of vectors and insert with conditions
		Bone = Bone->getNext();
		//++z;
		//if (z == 64)18
		//		break;
	}
}

void Reconstruct::makeSkelet(/*CImage image, TPolFigure* skelet*/){
	BitRaster *srcimg = new BitRaster(image.GetWidth(), image.GetHeight());
    bool inverted = true/*false*//*true*/;
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

	
	skeleton = new TPolFigure(srcimg, 0);
    //skeleton->MakeTriangDel();
	//skeleton->CutSkeleton(0);

     //BondSkeletTrans(srcimg, 0, 0/*10*//*100*/, skeleton);
    
	 //skeleton->CutSkeleton(1);
    //skeleton->setFakeKind();
    //skeleton->fakeCutSkeleton(1);
    //this->update();
}

void Reconstruct::SetHeightforBorders(TConnected* Component, std::set<Point>& sPoints, std::set<Point>& sPointsEdge, int firstH_, int secondH_) {

	firstH = firstH_;
	secondH = secondH_;
	//Point * Node = skeleton->Components->first()->Border->ListPoints->first();//->Nodes->first();
	Element* el = Component->Border->Elements[0];
	std::vector<double> colors;
	while (el) {
		//el->f = firstH;
		if (el->isVertex) {
			//*((Vertex*)el)->p;
			//Point p = *((Vertex*)el)->p;
			if (sPointsEdge.find(*((Vertex*)el)->p) != sPointsEdge.end())
				el->f = secondH;//colors.push_back(secondH);//
			else if (sPoints.find(*((Vertex*)el)->p) != sPoints.end())
				el->f = firstH;//colors.push_back(firstH);//
			else {
				el->f = (secondH + firstH) / 2.0;//colors.push_back((secondH + firstH) / 2.0);//
				pointsvert.push_back( ((Vertex*)el)->p );
			}
			/*if (sPoints.find(p) != sPoints.end())
				el->f = firstH;                       //here
			else if (sPointsEdge.find(p) != sPointsEdge.end())
				el->f = secondH;
			//el->f = (firstH+secondH)/2.0;//secondH;//1e+3;
			else
				el->f = (firstH + secondH) / 2.0;//secondH;//1e+3;*/
				//el->f = secondH;
			//imageF[p.X][p.Y] = el->f;
		}
		el = el->getNext();
	}
	
	/*for (auto i = sPoints.begin(); i != sPoints.end(); ++i) {
		Point* p = new Point(*i);
		pointsvert.push_back(p);
	}*/

	el = Component->Border->Elements[1];
	int i = 0;
	while (el) {

		if (el->isVertex) {
			el = el->getNext();
			continue;
		}

		if (el->getNextLooped()->f*1.0 != (firstH + secondH) / 2.0 &&
			el->getPrevLooped()->f*1.0 != (firstH + secondH) / 2.0 &&
			el->getNextLooped()->f != el->getPrevLooped()->f)
		{
			pointsvert.push_back(((Vertex*)el->getNextLooped())->p);
			pointsvert.push_back(((Vertex*)el->getPrevLooped())->p);

			/*Point* p1 = ((Vertex*)el->getPrevLooped())->p;
			Point* p2 = ((Vertex*)el->getNextLooped())->p;
			Point* p3 = new Point((p1->X + p2->X) / 2, (p1->Y + p2->Y) / 2);
			Vertex* newVert = new Vertex(p3);
			Edge* newEdge1 = new Edge(p1, p3);
			Edge* newEdge2 = new Edge(p3, p2);
			newEdge1->moveAsNextFor(el);
			newVert->f = 0.5;
			newVert->moveAsNextFor(newEdge1);
			newEdge2->moveAsNextFor(newVert);
			el->removeFromCurrentList();
			el = newEdge1;*/
		}
	
		el->f = el->getPrevLooped()->f;//(el->getPrevLooped()->f + el->getNextLooped()->f) / 2.0; //min????????el->getPrevLooped()->f;
		
		if (el->getNextLooped()->f*1.0 == (firstH + secondH) / 2.0) {
			el->f = el->getPrevLooped()->f;
			//pointsvert.push_back( ((Vertex*)el->getNextLooped())->p );
		}
		if (el->getPrevLooped()->f*1.0 == (firstH + secondH) / 2.0) {
			el->f = el->getNextLooped()->f;
			//pointsvert.push_back(((Vertex*)el->getPrevLooped())->p);
		}

		//if (el->f == 0.5)
		//	pointsvert.push_back(new Point( ( ((Vertex*)el->getPrevLooped())->p->X + ((Vertex*)el->getNextLooped())->p->X) / 2 ,  ( ((Vertex*)el->getPrevLooped())->p->Y  + ((Vertex*)el->getNextLooped())->p->Y ) / 2)) ;

		el = el->getNext();
	}

	if (Component->HoleList.size() > 0)
		el = Component->HoleList[0]->Elements[0];//skeleton->Components->first()->Border->Elements[1];
	while (el) {

		el->f = secondH;
		el = el->getNext();
	}
}
Reconstruct::Reconstruct(CImage im, int col1, int col2)
    : image(im)//, firstH(col1), secondH(col2)
{
    skeleton = NULL;
    //image = im;
    makeSkelet(/*image,skeleton*/ );
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

void Reconstruct::SetInnerPointsofSkelet(TConnected* Component){//(TPolFigure* skeleton, std::vector<std::vector<double>>& imageF){
    TNode * Node = Component->Nodes->first();
    while (Node){
        int sumint = 0;
        int sumnotint = 0;
        for (int i = 0; i < 3/*Node->Kind()*/; ++i){
			if (!Node->Sites[i])
				break;
            if ( Node->Sites[i]->f == firstH/*Cont->Internal*/ ){
                sumint++;
            } else if (Node->Sites[i]->f == secondH){
                sumnotint++;
			}
			else
			{
				sumint++;
				sumnotint++;
			}
        }
        int l = firstH;
        int r = secondH;
        if (sumint != 0 && sumnotint != 0){//(sumint != 0 && sumnotint != 0){//(sumint%2 == 1){
            Node->f = (l+r)/2.0;
            imageF[Node->Disc->X][Node->Disc->Y] = Node->f;
        }
        Node = Node->getNext();
    }
}

void SetOuterPointsofSkelet(TNode* Node, std::vector<std::vector<double>>& imageF){
    if (!Node || Node->f < 0)
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
		double val = nextNode->Sites[0]->f;//0;
        //if (nextNode->Sites[0]->Cont->Internal)
        //    val = 1;
		if (val < 0)
			std::cout << "Kosyak";
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
            f = curcell.cellels[k]->f/*Cont->Internal? 1:0*/;
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
		//f = 0;
		//return;
		FindParabolaPoint(bone, i, j, res);
    }
    if (res.X < 0 || res.X > imageF.size() || res.Y<0 || res.Y > imageF[0].size())
        return;
    d3 = DistPoint(&Point(res.X, res.Y), &Point(i, j));
    x2 = res.X;
    y2 = res.Y;
//     if (imageF[int(res.X)][int(res.Y)] > -0.5){
//         f = imageF[int(res.X)][int(res.Y)];
//     } else {
        double d1, d2;
        d1 = DistPoint(&Point(res.X, res.Y), &Point(bone->dest->Disc->X,bone->dest->Disc->Y));//&bone->dest);
        d2 = DistPoint(&Point(res.X, res.Y), &Point(bone->org->Disc->X, bone->org->Disc->Y));//&bone->org);
		if (bone->dest->f < 0 || bone->org->f < 0)
			f = 0;
        f = bone->dest->f * (d2/(d1+d2)) + bone->org->f * (d1/(d1+d2));
//    }
}
    
int Reconstruct::vertPart(CImage imVert) 
{
	//std::sort(pointsvert.begin(), pointsvert.end());
	/*for (int i = 1; i < pointsvert.size(); i+=2)
	{
		Cell newcell = Cell();
		//if (pointsvert[i] == pointsvert[i-1])
		newcell.borders.push_back(std::pair<Point, Point>( *pointsvert[i], *pointsvert[i-1] ));
		//newcell.borders.push_back(std::pair<Point, Point>(pp2, Point(orgnode->X(), orgnode->Y())));
		//newcell.borders.push_back(std::pair<Point, Point>(pp1, Point(destnode->X(), destnode->Y())));
		newcell.borders_color.push_back(std::pair<double, double>((firstH+secondH)/2.0, (firstH + secondH) / 2.0) );
		cells.push_back(newcell);
	}*/
	for (int i = 1; i < pointsvert.size(); i += 2)
	{
		Cell newcell = Cell();
		newcell.borders.push_back(std::pair<Point, Point>(*pointsvert[i], *pointsvert[i - 1]));                        //if there are not one
		newcell.borders_color.push_back(std::pair<double, double>((firstH + secondH) / 2.0, (firstH + secondH) / 2.0));
		cells.push_back(newcell);
	}
	/*makeSkelet();
	TConnected* Com = skeletonVert->Components->first();
	while (Com)
	{
		TNode * Node = Com->Nodes->first();
		addVerticalCells(Com, firstH, secondH);
		Com = Com->getNext();
	}*/
	return 0;
}

int Reconstruct::mainPart(){
	//skeleton->MakeTriangDel();
	//skeleton->CutSkeleton(0);
    TConnected* Com = skeleton->Components->first();
	while (Com)
	{
		SetInnerPointsofSkelet(Com);   //SetInnerPoints(skeleton, imageF);
		//set outer points
		TNode * Node = Com->Nodes->first();
		while (Node) {
			SetOuterPointsofSkelet(Node, imageF);
			Node = Node->getNext();
		}
		//////////////////////////////////////////////////////////////////////////
		//PaintSkeletBones(skeleton, imageF, false);
		//PaintBorders(skeleton, imageF);//, std::set<Point>());
		//PaintInnerBorders(skeleton, imageF);

		makeCells(Com);
		Com = Com->getNext();
	}
	
   /* sort(cells.begin(), cells.end(), CellComp); //cellcompare()
     double d1, d2, h1, h2;
     double x, y;
 	std::cout << image.GetWidth() << " " << image.GetHeight();

     for (int i = 0; i < image.GetWidth(); ++i)
         for (int j = 0; j < image.GetHeight(); ++j){
             if (!srcimg->getBit(i, j)){
                 imageF[i][j] = 0;
 			}
 			//else imageF[i][j] = 1;
             if (imageF[i][j] > -0.5)
                 continue;
             Cell curcell;
 			if (!FindCell(i, j, cells, curcell)) {
 				imageF[i][j] = 0;
 				continue;
 			}
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
         }*/
//         //RePaintSkeletBones(skeleton, imageF);

//     PaintInFile(imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\after_cur_out.png"));
    return 0;
}
