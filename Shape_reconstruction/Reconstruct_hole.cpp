#include <cmath>
#include "Processing.h"
#include "geom_utils.h"
#include "paint_utils.h"
//namespace T {

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
		for (auto iterat = cell.borders.begin(); iterat != cell.borders.end(); ++iterat) {
			Point firstp = iterat->first;
			Point secondp = /*cell.borders[i]*/iterat->second;

			if (firstp.Y == secondp.Y)
				continue;
			Point res;
			if (intersect(firstp, secondp, p, Point(0, p.Y), res)) {
				if (abs(int(std::floor(res.X + 0.5)) - res.X) < 1e-3)
					res.X = int(std::floor(res.X + 0.5));
				if (abs(int(std::floor(res.Y + 0.5)) - res.Y) < 1e-3)
					res.Y = int(std::floor(res.Y + 0.5));
				if ((res.X - min(secondp.X, firstp.X)) > -1e-5 && (max(secondp.X, firstp.X) - res.X) > -1e-5 &&
					(res.Y - min(secondp.Y, firstp.Y)) > -1e-5 && (max(secondp.Y, firstp.Y) - res.Y) > -1e-5 &&//res.Y >= min(secondp.Y, firstp.Y) && res.Y <= max(secondp.Y, firstp.Y) &&
					(res.X > p.X) == (0 > p.X)) {
					//if (points.find(res) != points.end())
					//    continue;
					if (abs(res.X - secondp.X) < 1e-5 && abs(res.Y - secondp.Y) < 1e-5 && secondp.Y > firstp.Y ||
						abs(res.X - firstp.X) < 1e-5 && abs(res.Y - firstp.Y) < 1e-5 && secondp.Y < firstp.Y) {
						continue;
					}

					f = !f;
				}
			}
		}

		return f;
	}

	bool FindCell(int x, int y, vector<Cell>& cells, Cell& outcell) {
		int idx = 0;

		while (idx < cells.size() && cells[idx].rightx < x)
			++idx;
		for (int i = idx; i < cells.size(); ++i)
			if (pointINcell(cells[i], Point(x, y))) {
				idx = i;
				break;
			}
			else idx = -1;
		if (idx == -1 || idx == cells.size())
			return false;
		outcell = cells[idx];
		return true;
	}

	void Reconstruct::addVerticalCells(TConnected* Component, int color1, int color2) {
		TBone* Bone = Component->Bones->first();
		while (Bone) {
			Cell newcell = Cell();
			newcell.skeletbone = Bone;
			TNode* orgnode = Bone->org;
			TNode* destnode = Bone->dest;

			destnode->f = orgnode->f = (color1 + color2) / 2;

			bool found = false;
			std::vector<int> bothsites;
			for (int i = 0; i < 3; ++i) {
				if (!orgnode->Sites[i])
					break;
				for (int j = 0; j < 3; ++j) {
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

							Edge* edge = (Edge*)destnode->Sites[j];
							if (edge->f == (firstH + secondH) / 2.0)
							{
								continue;
							}
							//	continue;
							newcell.nodes.push_back(pp1);
							newcell.borders.push_back(std::pair<Point, Point>(pp1, Point(orgnode->X(), orgnode->Y())));
							newcell.borders_color.push_back(std::pair<double, double>(edge->f, orgnode->f));
							newcell.borders.push_back(std::pair<Point, Point>(pp1, Point(destnode->X(), destnode->Y())));
							newcell.borders_color.push_back(std::pair<double, double>(edge->f, destnode->f));

							newcell.bords[Borders::Wall].push_back(std::pair<Point, Point>(pp1, Point(orgnode->X(), orgnode->Y())));
							newcell.bords[Borders::Wall].push_back(std::pair<Point, Point>(pp1, Point(destnode->X(), destnode->Y())));
							newcell.bords_color[Borders::Wall].push_back(std::pair<double, double>(orgnode->Sites[i]->f, orgnode->f));
							newcell.bords_color[Borders::Wall].push_back(std::pair<double, double>(orgnode->Sites[i]->f, destnode->f));


							// 						newcell.bordersNodes.insert(std::pair<TSite*, TNode*>(orgnode->Sites[i], orgnode));
							// 						newcell.bordersNodes.insert(std::pair<TSite*, TNode*>(orgnode->Sites[i], orgnode));

							newcell.cellels.push_back(orgnode->Sites[i]);
							//                         PaintLine(imageF, &pp1, &Point(orgnode->X(), orgnode->Y()), 1);
							//                         PaintLine(imageF, &pp1, &Point(destnode->X(), destnode->Y()), 1);
						}
						else {
							Edge* edge = (Edge*)destnode->Sites[j];

							bool bachangecolor = false;

							//if (orgnode->Disc->Rad == -1 || destnode->Disc->Rad == -1)


							Point pp1 = get_perpendicular_pt_from_pt_to_line(*edge->dest, *edge->org, Point(destnode->X(), destnode->Y()));
							//                         if (pp1.X == destnode->X() && pp1.Y == destnode->Y())
							//                             break;
						
							if (destnode->Disc->Rad == -1 && edge->Cont->Internal || orgnode->Disc->Rad == -1 && !edge->Cont->Internal) {
								//continue;
								edge->f = edge->getNextLooped()->f; // org or dest it is point
								bachangecolor = true;
							}
							if (orgnode->Disc->Rad == -1 && edge->Cont->Internal || destnode->Disc->Rad == -1 && !edge->Cont->Internal) {
					//			continue;
								edge->f = edge->getNextLooped()->f; // org or dest
								bachangecolor = true;
							}
							

							Point pp2 = get_perpendicular_pt_from_pt_to_line(*edge->dest, *edge->org, Point(orgnode->X(), orgnode->Y()));

							newcell.nodes.push_back(pp1);
							newcell.nodes.push_back(pp2);
							//if ( edge->WestDirect() )
					//		if (edge->f != secondH && edge->f != 6)
					//			continue;
							if (edge->f == (firstH + secondH) / 2.0 && !(orgnode->Disc->Rad == -1) && destnode->Disc->Rad != -1)
							{

								/*int val = -1;
								int x = pp2.X;
								int y = pp2.Y;
								//set+vector
								std::map<int, int> mvals;
								//for (int i = x - 1; i <= x + 1; ++i)
								//	for (int j = y - 1; j <= y + 1; ++j)
								//	{
								int i = x, j = y;
								mvals[im[i][j]]++;
								mvals[im[i - 1][j]]++;
								mvals[im[i + 1][j]]++;
								mvals[im[i][j - 1]]++;
								mvals[im[i][j + 1]]++;

								if (mvals.size() != 4) {
									//max
									if ((*mvals.begin()).first == 0)
										val = (*mvals.rbegin()).first;
									else
										val = ((*mvals.begin()).first == firstH ? secondH : firstH);
								}
								else
									val = -1;
								edge->f = val;*/
								TBone* leftBone = Bone;
								TBone* rightBone = Bone;

								while (leftBone->org->Disc->Rad != -1 && rightBone->dest->Disc->Rad != -1)
								{
									leftBone = leftBone->getPrevLooped();
									rightBone = rightBone->getNextLooped();
								}
								if (rightBone->dest->Disc->Rad == -1 && edge->Cont->Internal || leftBone->org->Disc->Rad == -1 && !edge->Cont->Internal)
								{
									Element* leftedge = edge->getPrevLooped();                                  //not edge may be hole
									//edge->f = edge->getNextLooped()->f;
									while (leftedge->f == (firstH + secondH) / 2.0)
										leftedge = leftedge->getPrevLooped();
									edge->f = leftedge->f;

								}
								if (leftBone->org->Disc->Rad == -1 && edge->Cont->Internal || rightBone->dest->Disc->Rad == -1 && !edge->Cont->Internal)//(rightBone->dest->Disc->Rad == -1)
								{
									Element* rightedge = edge->getNextLooped();								 //not edge may be hole
									//edge->f = edge->getNextLooped()->f;
									while (rightedge->f == (firstH + secondH) / 2.0)
										rightedge = rightedge->getNextLooped();
									edge->f = rightedge->f;

								}
								bachangecolor = true;
							}
							//if (edge->f != (firstH+secondH)/2.0)
							{
								newcell.borders.push_back(std::pair<Point, Point>(pp2, Point(orgnode->X(), orgnode->Y())));
								newcell.borders.push_back(std::pair<Point, Point>(pp1, Point(destnode->X(), destnode->Y())));
								newcell.borders.push_back(std::pair<Point, Point>(pp2, pp1));
							}

							//if (edge->f != (firstH + secondH) / 2.0)
							{
								newcell.borders_color.push_back(std::pair<double, double>(edge->f, orgnode->f));
								newcell.borders_color.push_back(std::pair<double, double>(edge->f, destnode->f));
								newcell.borders_color.push_back(std::pair<double, double>(edge->f, edge->f));

							}
							/*if (edge->f == (firstH + secondH) / 2.0)
							{
								//pps
								
								newcell.borders.push_back(std::pair<Point, Point>(pp2, Point(orgnode->X(), orgnode->Y())));
								newcell.borders.push_back(std::pair<Point, Point>(pp1, Point(destnode->X(), destnode->Y())));
								newcell.borders.push_back(std::pair<Point, Point>(pp2, pp1));
								newcell.borders_color.push_back(std::pair<double, double>(firstH, orgnode->f));
								newcell.borders_color.push_back(std::pair<double, double>(firstH, destnode->f));
								newcell.borders_color.push_back(std::pair<double, double>(firstH, firstH));

								newcell.borders.push_back(std::pair<Point, Point>(pp2, Point(orgnode->X(), orgnode->Y())));
								newcell.borders.push_back(std::pair<Point, Point>(pp1, Point(destnode->X(), destnode->Y())));
								newcell.borders.push_back(std::pair<Point, Point>(pp2, pp1));
								newcell.borders_color.push_back(std::pair<double, double>(secondH, orgnode->f));
								newcell.borders_color.push_back(std::pair<double, double>(secondH, destnode->f));
								newcell.borders_color.push_back(std::pair<double, double>(secondH, secondH));

							}*/

							{
								//if (edge->f < secondH)//?????????????
								newcell.bords[Borders::Floor].push_back(std::pair<Point, Point>(pp1, pp2));
								newcell.bords_color[Borders::Floor].push_back(std::pair<double, double>(edge->f, edge->f));

								newcell.bords[Borders::Wall].push_back(std::pair<Point, Point>(pp2, Point(orgnode->X(), orgnode->Y())));
								newcell.bords[Borders::Wall].push_back(std::pair<Point, Point>(pp1, Point(destnode->X(), destnode->Y())));
								newcell.bords_color[Borders::Wall].push_back(std::pair<double, double>(edge->f, orgnode->f));
								newcell.bords_color[Borders::Wall].push_back(std::pair<double, double>(edge->f, destnode->f));

							}
							if (bachangecolor)
							{
								edge->f = (firstH + secondH) / 2.0;
							}

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

	void Reconstruct::makeSkelet(/*CImage image, TPolFigure* skelet*/) {
		BitRaster *srcimg = new BitRaster(image.GetWidth(), image.GetHeight());
		bool inverted = true/*false*//*true*/;
		for (int i = 0; i < image.GetHeight(); i++) {
			for (int j = 0; j < image.GetWidth(); j++) {
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


		skeleton = new TPolFigure(srcimg, -1.0);
		//skeleton->MakeTriangDel();
		//skeleton->CutSkeleton(0);

		 //BondSkeletTrans(srcimg, 0, 0/*10*//*100*/, skeleton);

		 //skeleton->CutSkeleton(1);
		//skeleton->setFakeKind();
		//skeleton->fakeCutSkeleton(1);
		//this->update();
	}

	/*void Reconstruct::SetHeightforBorders(TConnected* Component, std::set<Point>& sPoints, std::set<Point>& sPointsEdge, int firstH_, int secondH_) {

		static int comp = 0;
		firstH = firstH_;
		secondH = secondH_;
		//Point * Node = skeleton->Components->first()->Border->ListPoints->first();//->Nodes->first();
		Element* el = Component->Border->Elements[0];
		std::vector<double> colors;
		while (el) {
			el->f = firstH;
			//if (el->isVertex) {
			//	}

			//}
			el = el->getNext();
		}

		if (Component->HoleList.size() > 0)
			el = Component->HoleList[0]->Elements[0];//skeleton->Components->first()->Border->Elements[1];
		while (el) {

			el->f = secondH;
			el = el->getNext();
		}
		++comp;
	}*/

/*	void Reconstruct::SetHeightforBorders(TConnected* Component, std::set<Point>& sPoints, std::set<Point>& sPointsEdge, int firstH_, int secondH_) {

		static int comp = 0;
		firstH = firstH_;
		secondH = secondH_;
		//Point * Node = skeleton->Components->first()->Border->ListPoints->first();//->Nodes->first();
		Element* el = Component->Border->Elements[0];
		std::vector<double> colors;
		while (el) {
			el->f = firstH;
			el = el->getNext();
			}

		if (Component->HoleList.size() > 0)
			el = Component->HoleList[0]->Elements[0];//skeleton->Components->first()->Border->Elements[1];
		while (el) {

			el->f = secondH;
			el = el->getNext();
		}	
	}
	*/
	void Reconstruct::SetHeightforBorders(TConnected* Component, std::set<Point>& sPoints1, std::set<Point>& sPoints2, int firstH_, int secondH_) {

		//static int comp = 0;
		firstH = firstH_;
		secondH = secondH_;
		Element* el = Component->Border->Elements[0];//Points
		//std::vector<double> colors;
		int numgrey = 0;
		while (el) {
			if (el->isVertex) {



			/*	int x = ((Vertex*)el)->p->X;
				int y = ((Vertex*)el)->p->Y;
				//set+vector
				std::map<int, int> mvals;
				for (int i = x - 1; i <= x + 1; ++i)
					for (int j = y - 1; j <= y + 1; ++j)
					{
						if (mvals.find(this->im[i][j]) != mvals.end())
							mvals[this->im[i][j]]++;
						else
							mvals[this->im[i][j]] = 1;
					}
				if (mvals.size() == 2) {
					//max
					if ((*mvals.begin()).first == 0)
						el->f = (*mvals.rbegin()).first;
					else
						el->f = ( (*mvals.begin()).first == firstH ? secondH : firstH );
				}
				else
					el->f = 500;
			}*/

				if (sPoints2.find(*((Vertex*)el)->p) != sPoints2.end())
					if (sPoints1.find(*((Vertex*)el)->p) != sPoints1.end()) {
						el->f = (secondH + firstH) / 2.0;
						++numgrey;
					}
					else {
						/*int halfnumgrey = numgrey/2;
						Element *tmpel = el->getPrevLooped()->getPrevLooped();
						while (numgrey > 0)
						{
							if (halfnumgrey > 0)
								tmpel->f = secondH;
							else
								tmpel->f = firstH;
							halfnumgrey--;
							numgrey--;
							tmpel = tmpel->getPrevLooped()->getPrevLooped();
						}*/
						el->f = secondH;
					}
				else if (sPoints1.find(*((Vertex*)el)->p) != sPoints1.end())
				{
					/*int halfnumgrey = numgrey / 2;
					Element *tmpel = el->getPrevLooped()->getPrevLooped();
					while (numgrey > 0)
					{
						if (halfnumgrey > 0)
							tmpel->f = firstH;
						else
							tmpel->f = secondH;
						halfnumgrey--;
						numgrey--;
						tmpel = tmpel->getPrevLooped()->getPrevLooped();
					}*/
					el->f = firstH;
				}
				else {
					el->f = (secondH + firstH) / 2.0;
					numgrey++;
					//el->f = el->getPrevLooped()->getPrevLooped()->f;
					}
				}
			el = el->getNext();
		}

		el = Component->Border->Elements[1];
		//int i = 0;
		while (el) {
			if (el->isVertex) {
				el = el->getNext();
				continue;
			}
			el->f = el->getPrevLooped()->f;
			if (el->getNextLooped()->f*1.0 != (firstH + secondH) / 2.0 &&
				el->getPrevLooped()->f*1.0 != (firstH + secondH) / 2.0 &&
				el->getNextLooped()->f != el->getPrevLooped()->f)
			{
				el->f = (firstH + secondH) / 2.0;


				/*Point* p1 = ((Vertex*)el->getPrevLooped())->p;
				Point* p2 = ((Vertex*)el->getNextLooped())->p;
				Point* p3 = new Point((p1->X + p2->X) / 2, (p1->Y + p2->Y) / 2);
				Vertex* newVert = new Vertex(p3);
				Edge* newEdge1 = new Edge(p1, p3);
				Edge* newEdge2 = new Edge(p3, p2);
				newEdge1->f = ((Vertex*)el->getPrevLooped())->f;
				newEdge2->f = ((Vertex*)el->getNextLooped())->f;
				newEdge1->moveAsNextFor(el);
				newVert->f = 0.5;
				newVert->moveAsNextFor(newEdge1);
				newEdge2->moveAsNextFor(newVert);
				el->removeFromCurrentList();
				el = newEdge2; */
			/*	{	
					TBone * bone = new TBone();
					TNode * newnode = new TNode();
					newnode->Disc = NewDisc(5, 5, -1); //without sites and other things now
					TBone * boneepred = skeleton->Components->first()->Bones->first();
					//newnode->Bones[0] = boneepred;
					newnode->Sites[0] = boneepred->dest->Sites[0]; //which similar
					newnode->Sites[1] = boneepred->dest->Sites[1];
					//bone->Com = skeleton->Components->first();
					bone->dest = boneepred->dest;
					boneepred->dest = newnode;
					bone->org = newnode;
					bone->dest->DetachBone(boneepred);
					bone->dest->AddBone(bone);
					newnode->AddBone(bone);
					newnode->AddBone(boneepred);
					skeleton->Components->first()->Bones->moveAsNextFor(boneepred);
					newnode->moveAsNextFor(skeleton->Components->first()->Nodes->first());
				}*/
				
			} else 
			if (el->getNextLooped()->f*1.0 == (firstH + secondH) / 2.0) {
				el->f = el->getPrevLooped()->f;
			} else 

			if (el->getPrevLooped()->f*1.0 == (firstH + secondH) / 2.0) {
				el->f = el->getNextLooped()->f;
			}
			//el->f = el->getPrevLooped()->f;

//			el->f = (el->getPrevLooped()->f + el->getNextLooped()->f) / 2.0; //min????????el->getPrevLooped()->f; el->getPrevLooped()->f;//
			el = el->getNext();
		}

		
		int holen = 0;
		while (Component->HoleList.size() > holen) {
			el = Component->HoleList[holen]->Elements[0];//skeleton->Components->first()->Border->Elements[1];
			while (el) {
				if (el->isVertex) {
					if (sPoints2.find(*((Vertex*)el)->p) != sPoints2.end())
						if (sPoints1.find(*((Vertex*)el)->p) != sPoints1.end())
							el->f = (secondH + firstH) / 2.0;
						else
							el->f = secondH;
					else if (sPoints1.find(*((Vertex*)el)->p) != sPoints1.end())
						el->f = firstH;
					else {
						el->f = (secondH + firstH) / 2.0;
					}
				}
				/*if (el->isVertex) {
					int x = ((Vertex*)el)->p->X;
					int y = ((Vertex*)el)->p->Y;
					//set+vector
					std::map<int, int> mvals;
					for (int i = x - 1; i <= x + 1; ++i)
						for (int j = y - 1; j <= y + 1; ++j)
						{
							if (mvals.find(this->im[i][j]) != mvals.end())
								mvals[this->im[i][j]]++;
							else
								mvals[this->im[i][j]] = 1;
						}
					if (mvals.size() == 2) {
						//max
						if ((*mvals.begin()).first == 0)
							el->f = (*mvals.rbegin()).first;
						else
							el->f = ((*mvals.begin()).first == firstH ? secondH : firstH);
					}
					else
						el->f = 500;
				}*/
				el = el->getNext();
			}
			el = Component->HoleList[holen]->Elements[1];
			//int i = 0;
			while (el) {

				if (el->isVertex) {
					el = el->getNext();
					continue;
				}
				el->f = el->getPrevLooped()->f;

				if (el->getNextLooped()->f*1.0 != (firstH + secondH) / 2.0 &&
					el->getPrevLooped()->f*1.0 != (firstH + secondH) / 2.0 &&
					el->getNextLooped()->f != el->getPrevLooped()->f)
				{
					el->f = (firstH + secondH) / 2.0;
				}
				if (el->getNextLooped()->f*1.0 == (firstH + secondH) / 2.0) {
					el->f = el->getPrevLooped()->f;
				}

				if (el->getPrevLooped()->f*1.0 == (firstH + secondH) / 2.0) {
					el->f = el->getNextLooped()->f;
				}
				//el->f = (el->getPrevLooped()->f + el->getNextLooped()->f) / 2.0;
				//el->f = el->getPrevLooped()->f;

				el = el->getNext();
			}
			//++comp;
			++holen;
		}
	}

	Reconstruct::Reconstruct(CImage im, int col1, int col2)
		: image(im)//, firstH(col1), secondH(col2)
	{
		skeleton = NULL;
		//image = im;
		makeSkelet(/*image,skeleton*/);
		imageF.resize(image.GetWidth());
		for (int i = 0; i < image.GetWidth(); ++i) {
			imageF[i].resize(image.GetHeight());
			std::fill(imageF[i].begin(), imageF[i].begin() + image.GetHeight(), -1);
		}
		//imageF = std::vector<std::vector<doub>>(image.GetWidth(), std::vector<int>(image.GetHeight()));
	}

	Reconstruct::~Reconstruct()
	{
	}

	void Reconstruct::SetInnerPointsofSkelet(TConnected* Component) {//(TPolFigure* skeleton, std::vector<std::vector<double>>& imageF){
		TNode * Node = Component->Nodes->first();
		while (Node) {
			int sumint = 0;
			int sumnotint = 0;
			for (int i = 0; i < 3/*Node->Kind()*/; ++i) {
				if (!Node->Sites[i])
					break;
				if (Node->Sites[i]->f == firstH/*Cont->Internal*/ && !Node->Sites[i]->Cont->Internal ) {
					sumint++;
				}
				else if (Node->Sites[i]->f == secondH && !Node->Sites[i]->Cont->Internal) {
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
			if (sumint != 0 && sumnotint != 0) {//(sumint != 0 && sumnotint != 0){//(sumint%2 == 1){
				Node->f = (l + r) / 2.0;
				imageF[Node->Disc->X][Node->Disc->Y] = Node->f;
			}
			Node = Node->getNext();
		}
	}

	void SetOuterPointsofSkelet(TNode* Node, std::vector<std::vector<double>>& imageF) {
		if (!Node || Node->f < 0)
			return;
		int count_bones = sizeof(Node->Bones) / sizeof(Node->Bones[0]);
		for (int i = 0; i < count_bones; ++i) {
			TBone* Bone = Node->Bones[i];
			if (!Bone || Bone->GetNextNode(Node)->f > -1)
				continue;
			TNode* nextNode = Bone->GetNextNode(Node);
			double d1 = /*CalculateDistance*/DistPoint(&Point(nextNode->Disc->X, nextNode->Disc->Y),
				&Point(Node->Disc->X, Node->Disc->Y));
			double d2 = CalculateDistanceToBorder(nextNode);
			double val = nextNode->Sites[0]->f;//0;
			//if (nextNode->Sites[0]->Cont->Internal)
			//    val = 1;
			if (val < 0)
				std::cout << "Kosyak";
			nextNode->f = Node->f * (d2 / (d2 + d1)) + val/*nearest*/ *(d1 / (d1 + d2));
			imageF[int(nextNode->X())][int(nextNode->Y())] = nextNode->f;
			SetOuterPointsofSkelet(nextNode, imageF);
		}
	}

	void Reconstruct::findClosestBorder(Cell& curcell, int i, int j, double&x, double&y, double& d1, double& f, std::vector<std::vector<double>>& imageF)
	{
		d1 = 1000000;
		double curd;
		Point pp;
		for (int k = 0; k < curcell.cellels.size(); ++k) {
			if (curcell.cellels[k]->isVertex) {
				pp = *((Vertex*)curcell.cellels[k])->p;
				//curd = DistPoint(&Point(i, j), &((Vertex*)curcell.cellels[i]).p);
			}
			else {
				Edge* edge = (Edge*)curcell.cellels[k];
				pp = get_perpendicular_pt_from_pt_to_line(*edge->dest, *edge->org, Point(i, j));
				//             curd = DistPoint(&Point(i, j), &pp);
			}
			curd = DistPoint(&Point(i, j), &pp);
			if (curd < d1) {
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
			Point((bone->dest->Disc->X), (bone->dest->Disc->Y)),
			Point((bone->org->Disc->X), (bone->org->Disc->Y)), res);
		if (bone->Virt) {
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
		d1 = DistPoint(&Point(res.X, res.Y), &Point(bone->dest->Disc->X, bone->dest->Disc->Y));//&bone->dest);
		d2 = DistPoint(&Point(res.X, res.Y), &Point(bone->org->Disc->X, bone->org->Disc->Y));//&bone->org);
		if (bone->dest->f < 0 || bone->org->f < 0)
			f = 0;
		f = bone->dest->f * (d2 / (d1 + d2)) + bone->org->f * (d1 / (d1 + d2));
		//    }
	}
	/*#include "opencv2/highgui/highgui.hpp"
	#include "opencv2/imgproc/imgproc.hpp"
	std::vector<std::vector<cv::Point>> selectContours(std::string& filename = std::string("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\figoutline.png"))
	{
		std::vector<std::vector<cv::Point> > contours;
		std::vector<cv::Vec4i> hierarchy;
		cv::Mat src;
		src = cv::imread(filename, 1);
		cvtColor(src, src, CV_BGR2GRAY);
		cv::findContours(src, contours, hierarchy, CV_RETR_LIST, CV_CHAIN_APPROX_SIMPLE, cv::Point(0, 0));
		return contours;
	}*/
#include "opencv2/imgproc/imgproc.hpp"
	int Reconstruct::vertPart(CImage imVert, std::vector<Point>& cvpoints)
	{
		for (int i = 0; i < cvpoints.size(); ++i)
		{
			Cell newcell = Cell();
			Point p = cvpoints[i];
			Point p2 = cvpoints[(i + 1) % cvpoints.size()];
			if (pow(p.X - p2.X, 2) + pow(p.Y - p2.Y, 2) > 1000)
				continue;

			newcell.borders.push_back(std::pair<Point, Point>(p, p2));
			newcell.borders_color.push_back(std::pair<double, double>((firstH + secondH) / 2.0, (firstH + secondH) / 2.0));

			newcell.borders.push_back(std::pair<Point, Point>(p, p2));
			newcell.borders_color.push_back(std::pair<double, double>((firstH), (firstH)));

			newcell.borders.push_back(std::pair<Point, Point>(p, p2));
			newcell.borders_color.push_back(std::pair<double, double>((secondH), (secondH)));

			cells.push_back(newcell);
		}
		return 0;
	}
	int Reconstruct::vertPart(CImage imVert, std::vector<cv::Point>& cvpoints)
	{
		//std::sort(pointsvert.begin(), pointsvert.end());
		/*for (int i = 0; i < pointsvert.size(); ++i)
			for (int j = i + 1; j < pointsvert.size(); ++j)
			{
				if (partpointsvert[i] != partpointsvert[j])
				{
					Cell newcell = Cell();
					newcell.borders.push_back(std::pair<Point, Point>(*pointsvert[i], *pointsvert[j]));
					newcell.borders_color.push_back(std::pair<double, double>((firstH + secondH) / 2.0, (firstH + secondH) / 2.0));
					cells.push_back(newcell);
				}
			}*/
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

			//for (auto i = cvpoints.begin(); i != cvpoints.end(); ++i) 

		for (int i = 0; i < cvpoints.size(); ++i)
		{
			Cell newcell = Cell();
			cv::Point p = cvpoints[i];
			cv::Point p2 = cvpoints[(i + 1) % cvpoints.size()];

			newcell.borders.push_back(std::pair<Point, Point>(Point(p.x, p.y), Point(p2.x, p2.y)));
			newcell.borders_color.push_back(std::pair<double, double>((firstH + secondH) / 2.0, (firstH + secondH) / 2.0));

			newcell.borders.push_back(std::pair<Point, Point>(Point(p.x, p.y), Point(p2.x, p2.y)));
			newcell.borders_color.push_back(std::pair<double, double>((firstH), (firstH)));

			newcell.borders.push_back(std::pair<Point, Point>(Point(p.x, p.y), Point(p2.x, p2.y)));
			newcell.borders_color.push_back(std::pair<double, double>((secondH), (secondH)));

			cells.push_back(newcell);
		}

		/*std::vector<std::vector<cv::Point>> ps = selectContours();

		for (int i = 0; i < ps.size(); ++i)
		{
			for (int j = 0; j < ps[i].size()-1; ++j) {
				Cell newcell = Cell();
				cv::Point p = ps[i][j];//*ps[i].begin();
				cv::Point p2 = ps[i][j + 1];//*ps[i].rbegin();
				newcell.borders.push_back(std::pair<Point, Point>(Point(p.x, p.y), Point(p2.x, p2.y)));
				newcell.borders_color.push_back(std::pair<double, double>((firstH + secondH) / 2.0, (firstH + secondH) / 2.0));
				cells.push_back(newcell);
			}
		}*/
		/*for (int i = 1; i < pointsvert.size(); i += 2)
		{
			Cell newcell = Cell();
			newcell.borders.push_back(std::pair<Point, Point>(*pointsvert[i], *pointsvert[i - 1]));                        //if there are not one
			newcell.borders_color.push_back(std::pair<double, double>((firstH + secondH) / 2.0, (firstH + secondH) / 2.0));
			cells.push_back(newcell);
		}*/



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

	int Reconstruct::mainPart() {
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

	/*	sort(cells.begin(), cells.end(), CellComp); //cellcompare()
		  double d1, d2, h1, h2;
		  double x, y;
		 std::cout << image.GetWidth() << " " << image.GetHeight();

		  for (int i = 0; i < image.GetWidth(); ++i)
			  for (int j = 0; j < image.GetHeight(); ++j){
				  //if (!srcimg->getBit(i, j)){
				//	  imageF[i][j] = 0;
				 //}
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

			       PaintInFile(imageF, _T("C:\\Users\\Alexandra\\My\\Shape_reconstruction\\data\\after_cur_out_0.png"));
		return 0;
	}
//}