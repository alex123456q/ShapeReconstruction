#include <windows.h>
#include <GL/freeglut.h>   //Подключение библиотеки glut.h

#include "paint_utils.h"
#include "Processing.h"

//std::vector<std::vector<double>> imagetopaint;
std::vector<std::vector<Cell>> cellstopaint;
// float angle = 0.0f;
// float lx = 0, ly = 0;
// угол поворота камеры
float angle = -0.0133;// 0.0;
// координаты вектора направления движения камеры
float lx = 0.0133/*0*/ , lz = -1/*-1*/;
float scaleValue = 0.36f;//0.46f;//0.56f;//0.15f;
// XZ позиция камеры
float x = 40.0f, z = 205.0f;
//Ключи статуса камеры. Переменные инициализируются нулевыми значениями
//когда клавиши не нажаты
float deltaAngle = 0.0f;
float deltaMove = 0;
float deltaMovePerp = 0;
int xOrigin = -1;
// Коэффициент масштабирования по оси Y в соответствии с осями X и Z (НОВОЕ)
#define    HEIGHT_RATIO  1.f//1.5f
// void Initialize()
// {
// // 	//Выбрать фоновый (очищающий) цвет
// // 	glClearColor(1.0, 0.0, 1.0, 1.0);
// // 
// // 	//Установить проекцию
// // 	glMatrixMode(GL_PROJECTION);
// // 	glLoadIdentity();
// // 	glOrtho(0.0, 1.0, 0.0, 1.0, -1.0, 1.0);
// 	/* Initialize OpenGL Graphics */
// 		glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Set background color to black and opaque
// 		glClearDepth(1.0f);                   // Set background depth to farthest
// 		glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling
// 		glDepthFunc(GL_LEQUAL);    // Set the type of depth-test
// 		glShadeModel(GL_SMOOTH);   // Enable smooth shading
// 		glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);  // Nice perspective corrections
// }
// 
// void Draw()
// {
// // 	//Очищаем экран glClear(GL_COLOR_BUFFER_BIT);
// // 
// // 	//Отрисовка квадрата 
// // 	glColor3f(1.0, 1.0, 1.0); //Выбираем белый цвет
// // 	glBegin(GL_POLYGON);
// // 	glVertex3f(0.25, 0.25, 0.0); //Координаты квадрата
// // 	glVertex3f(0.75, 0.25, 0.0);
// // 	glVertex3f(0.75, 0.75, 0.0);
// // 	glVertex3f(0.25, 0.75, 0.0);
// // 
// // 	//another wqalls 
// // 	glEnd();
// // 	glFlush();
// 	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear color and depth buffers
// 	glMatrixMode(GL_MODELVIEW);     // To operate on model-view matrix
// 
// 									// Render a color-cube consisting of 6 quads with different colors
// 	glLoadIdentity();                 // Reset the model-view matrix
// 	glTranslatef(1.5f, 0.0f, -7.0f);  // Move right and into the screen
// 
// 	glBegin(GL_QUADS);                // Begin drawing the color cube with 6 quads
// 									  // Top face (y = 1.0f)
// 									  // Define vertices in counter-clockwise (CCW) order with normal pointing out
// 	glColor3f(0.0f, 1.0f, 0.0f);     // Green
// 	glVertex3f(1.0f, 1.0f, -1.0f);
// 	glVertex3f(-1.0f, 1.0f, -1.0f);
// 	glVertex3f(-1.0f, 1.0f, 1.0f);
// 	glVertex3f(1.0f, 1.0f, 1.0f);
// 
// 	// Bottom face (y = -1.0f)
// 	glColor3f(1.0f, 0.5f, 0.0f);     // Orange
// 	glVertex3f(1.0f, -1.0f, 1.0f);
// 	glVertex3f(-1.0f, -1.0f, 1.0f);
// 	glVertex3f(-1.0f, -1.0f, -1.0f);
// 	glVertex3f(1.0f, -1.0f, -1.0f);
// 
// 	// Front face  (z = 1.0f)
// 	glColor3f(1.0f, 0.0f, 0.0f);     // Red
// 	glVertex3f(1.0f, 1.0f, 1.0f);
// 	glVertex3f(-1.0f, 1.0f, 1.0f);
// 	glVertex3f(-1.0f, -1.0f, 1.0f);
// 	glVertex3f(1.0f, -1.0f, 1.0f);
// 
// 	// Back face (z = -1.0f)
// 	glColor3f(1.0f, 1.0f, 0.0f);     // Yellow
// 	glVertex3f(1.0f, -1.0f, -1.0f);
// 	glVertex3f(-1.0f, -1.0f, -1.0f);
// 	glVertex3f(-1.0f, 1.0f, -1.0f);
// 	glVertex3f(1.0f, 1.0f, -1.0f);
// 
// 	// Left face (x = -1.0f)
// 	glColor3f(0.0f, 0.0f, 1.0f);     // Blue
// 	glVertex3f(-1.0f, 1.0f, 1.0f);
// 	glVertex3f(-1.0f, 1.0f, -1.0f);
// 	glVertex3f(-1.0f, -1.0f, -1.0f);
// 	glVertex3f(-1.0f, -1.0f, 1.0f);
// 
// 	// Right face (x = 1.0f)
// 	glColor3f(1.0f, 0.0f, 1.0f);     // Magenta
// 	glVertex3f(1.0f, 1.0f, -1.0f);
// 	glVertex3f(1.0f, 1.0f, 1.0f);
// 	glVertex3f(1.0f, -1.0f, 1.0f);
// 	glVertex3f(1.0f, -1.0f, -1.0f);
// 	glEnd();  // End of drawing color-cube
// 
// 
// 	glutSwapBuffers();  // Swap the front and back frame buffers (double buffering)
// 
// }
// 
// /* Handler for window re-size event. Called back when the window first appears and
// whenever the window is re-sized with its new width and height */
// void reshape(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
// 											   // Compute aspect ratio of the new window
// 	if (height == 0) height = 1;                // To prevent divide by 0
// 	GLfloat aspect = (GLfloat)width / (GLfloat)height;
// 
// 	// Set the viewport to cover the new window
// 	glViewport(0, 0, width, height);
// 
// 	// Set the aspect ratio of the clipping volume to match the viewport
// 	glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
// 	glLoadIdentity();             // Reset
// 								  // Enable perspective projection with fovy, aspect, zNear and zFar
// 	gluPerspective(45.0f, aspect, 0.1f, 100.0f);
// }
// 
// // /* Main function: GLUT runs as a console application starting at main() */
// // int main(int argc, char** argv) {
// // 	glutInit(&argc, argv);            // Initialize GLUT
// // 	glutInitDisplayMode(GLUT_DOUBLE); // Enable double buffered mode
// // 	glutInitWindowSize(640, 480);   // Set the window's initial width & height
// // 	glutInitWindowPosition(50, 50); // Position the window's initial top-left corner
// // 	glutCreateWindow(title);          // Create window with the given title
// // 	glutDisplayFunc(display);       // Register callback handler for window re-paint event
// // 	glutReshapeFunc(reshape);       // Register callback handler for window re-size event
// // 	initGL();                       // Our own OpenGL initialization
// // 	glutMainLoop();                 // Enter the infinite event-processing loop
// // 	return 0;
// // }
// //Войти в главный цикл
// int main(int argc, char **argv)
// {
// 	glutInit(&argc, argv);
// 	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
// 	glutInitWindowSize(400, 400);		//Указываем размер окна
// 	glutInitWindowPosition(100, 100);	//Позиция окна
// 	glutCreateWindow("Polygon");		//Имя окна
// 	Initialize();						//Вызов функции Initialize
// 	glutDisplayFunc(Draw);				//Вызов функции отрисовки
// 	glutReshapeFunc(reshape);
// 	glutMainLoop();
// 
// 	return 0;
//}
//}
/* Global variables */
char title[] = "3D Shapes";

/* Initialize OpenGL Graphics */
void initGL() {
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Set background color to black and opaque
	glClearDepth(1.0f);                   // Set background depth to farthest
	glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling
	glDepthFunc(GL_LEQUAL);    // Set the type of depth-test
	glShadeModel(GL_SMOOTH);   // Enable smooth shading
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);  // Nice perspective corrections
}

void RenderingHeightMapLines(std::vector<std::vector<Cell>>& allcells/*, std::vector<std::vector<double>> imageF*/)
{
// 	glBegin(GL_QUADS);
// 	glColor3f(0.3f, 0.3f, 1.0f);
// 	glVertex3f(0, 0, 0);
// 	glVertex3f(0, 0, 582);
// 	glVertex3f(242, 0, 582);
// 	glVertex3f(242, 0, 0);
// 	glEnd();

// 	glBegin(GL_LINES); // Линии
// 	glColor3f(1.0f, 1.0f, 1.0f);
// 	glColor3f(1.0f, 1.0f, 1.0f);
	for (int k = 0; k < allcells.size(); ++k){
		for (int i = 0; i < allcells[k].size(); ++i) {
			glColor3f(1.0f, 1.0f, 1.0f);
			//auto cellborder = cells[i].bordersNodes.begin();
			//std::cout << " border size : " << cells[i].bordersNodes.size() <<  " " << cells[i].borders.size() << std::endl;
			int j = 0;
			for (auto bord = allcells[k][i].borders.begin(); bord != allcells[k][i].borders.end(); ++bord) {
				glBegin(GL_LINES); // Линии
				glColor3f(1.0f, 1.0f, 1.0f);

				if (!allcells[k][i].skeletbone) 
					glColor3f(1.0f, 0.0f, 0.0f);
				//if (fabs(imageF[int(bord->first.X)][int(bord->first.Y)]) <1e-3)
				//	glColor3f(1.0f, 0.7f, 0.7f);
				glVertex3f(bord->first.X, allcells[k][i].borders_color[j].first * 7 /*imageF[int(bord->first.X)][int(bord->first.Y)]*7*/, bord->first.Y);
				glColor3f(1.0f, 1.0f, 1.0f);
				if (!allcells[k][i].skeletbone)
					glColor3f(1.0f, 0.0f, 0.0f);
				//if (fabs(imageF[int(bord->second.X)][int(bord->second.Y)]) <1e-3)
				//	glColor3f(1.0f, 0.7f, 0.7f);
				glVertex3f(bord->second.X, allcells[k][i].borders_color[j].second * 7/*imageF[int(bord->second.X)][int(bord->second.Y)]*7*/, bord->second.Y);
				glEnd();;
				++j;
				//++cellborder;
			}
			glColor3f(0.7f, 1.0f, 0.7f);
			glBegin(GL_LINES); // Линии
			if (allcells[k][i].skeletbone) {
				glVertex3f(allcells[k][i].skeletbone->dest->X(), allcells[k][i].skeletbone->dest->f * 7, allcells[k][i].skeletbone->dest->Y());
				glVertex3f(allcells[k][i].skeletbone->org->X(), allcells[k][i].skeletbone->org->f * 7, allcells[k][i].skeletbone->org->Y());
			}
			glEnd();
		}
// 		for (int j = 0; j < cells[i].nodes.size()+1; ++j) {
// 			double x = cells[i].nodes[j%cells[i].nodes.size()].X;
// 			double y = cells[i].nodes[j%cells[i].nodes.size()].Y;
// 			glVertex3f(x, imageF[int(x)][int(y)]*7, y);  // Слева вверху
// // 			glVertex3f(1.0f, 1.0f, 0.0f);  // Справа вверху
// // 			glVertex3f(1.0f, -1.0f, 0.0f);  // Справа внизу
// // 			glVertex3f(-1.0f, -1.0f, 0.0f);  // Слева внизу
// 		}
		//glVertex3f(x, imageF[int(x)][int(y)], x);  // Слева вверху
	}
	//glEnd();
}

void RenderHeightMap(std::vector<std::vector<double>> pHeightMap) // Визуализация карты высоты с помощью квадратов
{
	int X = 0, Y = 0;    // Создаем пару переменных для перемещения по массиву
	int x, y, z;         // И еще три для удобства чтения

	if (pHeightMap.size() == 0) return;     // Данные корректны?
	//if (bRender)            // Что хотим визуализировать?
	glBegin(GL_QUADS); // Полигоны
	//else
	//	glBegin(GL_LINES); // Линии
	const int STEP_SIZE = 10;
	for (X = 0; X < pHeightMap.size()-STEP_SIZE; X += STEP_SIZE)
		for (Y = 0; Y < pHeightMap[0].size()-STEP_SIZE; Y += STEP_SIZE)
		{
			// Получаем (X, Y, Z) координаты нижней левой вершины
			x = X;
			y = pHeightMap[X][Y]*7;//Height(pHeightMap, X, Y);
			z = Y;

			// Устанавливаем цвет конкретной точки
			//SetVertexColor(pHeightMap, x, z);
			{
				float fColor = -0.15f + (pHeightMap[x][z]) /*/ 256.0f*/;
				//float fColor = 1.0;
				// Присвоить оттенок синего цвета для конкретной точки
				glColor4f(0.0f, 0.0f, fColor, 0.8);
			}

			glVertex3i(x, y, z);      // Визуализация ее

									  // Получаем (X, Y, Z) координаты верхней левой вершины
			x = X;
			y = pHeightMap[X][Y+STEP_SIZE]*7;//Height(pHeightMap, X, Y);
			z = Y + STEP_SIZE;

			// Устанавливаем цвет конкретной точки
			//SetVertexColor(pHeightMap, x, z);
			{
				float fColor = -0.15f + (pHeightMap[x][z]) /*/ 256.0f*/;
				//float fColor = 1.0;
				// Присвоить оттенок синего цвета для конкретной точки
				glColor4f(0.0f, 0.0f, fColor, 0.8);
			}

			glVertex3i(x, y, z);      // Визуализация ее

									  // Получаем (X, Y, Z) координаты верхней правой вершины
			x = X + STEP_SIZE;
			y = pHeightMap[X+STEP_SIZE][Y + STEP_SIZE]*7;//Height(pHeightMap, X, Y);
			z = Y + STEP_SIZE;

			// Устанавливаем цвет конкретной точки
			//SetVertexColor(pHeightMap, x, z);
			{
				float fColor = -0.15f + (pHeightMap[x][z]) /*/ 256.0f*/;
				//float fColor = 1.0;
				// Присвоить оттенок синего цвета для конкретной точки
				glColor4f(0.0f, 0.0f, fColor, 0.8);
			}

			glVertex3i(x, y, z);      // Визуализация ее

									  // Получаем (X, Y, Z) координаты нижней правой вершины
			x = X + STEP_SIZE;
			y = pHeightMap[X + STEP_SIZE][Y]*7;//Height(pHeightMap, X, Y);
			z = Y;

			// Устанавливаем цвет конкретной точки
			//SetVertexColor(pHeightMap, x, z);
			{
				float fColor = -0.15f + (pHeightMap[x][z]) /*/ 256.0f*/;
				//float fColor = 1.0;
				// Присвоить оттенок синего цвета для конкретной точки
				glColor4f(0.0f, 0.0f, fColor, 0.8);
			}

			glVertex3i(x, y, z);      // Визуализация ее
		}
	glEnd();
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);      // Сбрасываем цвет
}
/* Handler for window-repaint event. Called back when the window first appears and
whenever the window needs to be re-painted. */
void display() {
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear color and depth buffers
	glMatrixMode(GL_MODELVIEW);     // To operate on model-view matrix

									// Render a color-cube consisting of 6 quads with different colors
	glLoadIdentity();                 // Reset the model-view matrix
	glTranslatef(0.0f, 0.0f, 0.0f);  // Move right and into the screen


// 	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  // Очистка экрана и буфера глубины
// 	glLoadIdentity();          // Сброс просмотра
// 
// 							   //       Положение   Вид    Вектор вертикали
//  	gluLookAt(212, 60, 194, 186, 55, 171, 0, 1, 0);  // Определяет вид и положение камеры
	

     if (deltaMove || deltaMovePerp) {
		 x += deltaMove * lx * 0.1f - deltaMovePerp * lz * 0.1f;
		 z += deltaMove * lz * 0.1f + deltaMovePerp * lx * 0.1f;
	 }
	 //gluLookAt(50, 8, 225, 51, 2, 150, 0, 1, 0);  // Определяет вид и положение камеры


	gluLookAt(x, 5.0f, z,                    //0.8 1
		x + lx*75, 10.0f, z + lz*75,
		0.0f, 1.0f, 0.0f);

    glScalef(scaleValue, scaleValue * HEIGHT_RATIO, scaleValue);
	
	//RenderHeightMap(imagetopaint);        // Визализация карты высот
	RenderingHeightMapLines(cellstopaint/*, imagetopaint*/);

	glutSwapBuffers();  // Swap the front and back frame buffers (double buffering)
}
void processNormalKeys(unsigned char key, int xx, int yy) {

	if (key == 27)
		exit(0);
}

void pressKey(int key, int xx, int yy) {

	switch (key) {
	case GLUT_KEY_UP: deltaMove = 0.5f; break;
	case GLUT_KEY_DOWN: deltaMove = -0.5f; break;
	case GLUT_KEY_RIGHT: deltaMovePerp = 0.5f; break;
	case GLUT_KEY_LEFT: deltaMovePerp = -0.5f; break;
	}
}

void releaseKey(int key, int x, int y) {

	switch (key) {
	case GLUT_KEY_UP:
	case GLUT_KEY_DOWN: deltaMove = 0; break;
	case GLUT_KEY_RIGHT:
	case GLUT_KEY_LEFT: deltaMovePerp = 0; break;
	}
}

void mouseMove(int x, int y) {

	// this will only be true when the left button is down
	if (xOrigin >= 0) {

		// update deltaAngle
		deltaAngle = (x - xOrigin) * 0.001f;

		// update camera's direction
		lx = sin(angle + deltaAngle);
		lz = -cos(angle + deltaAngle);
	}
}

void mouseButton(int button, int state, int x, int y) {

	// only start motion if the left button is pressed
	if (button == GLUT_LEFT_BUTTON) {

		// when the button is released
		if (state == GLUT_UP) {
			angle += deltaAngle;
			xOrigin = -1;
		}
		else {// state = GLUT_DOWN
			xOrigin = x;
		}
	}
}
void processSpecialKeys(int key, int xx, int yy) {
	float fraction = 10.0f; //0.1
	switch (key) {
	case GLUT_KEY_LEFT:
		angle -= 0.01f;
		lx = sin(angle);
		lz = -cos(angle);
		break;
	case GLUT_KEY_RIGHT:
		angle += 0.01f;
		lx = sin(angle);
		lz = -cos(angle);
		break;
	case GLUT_KEY_UP:
		scaleValue += 0.001f;
		//x += lx * fraction;
		//z += lz * fraction;
		break;
	case GLUT_KEY_DOWN:
		scaleValue -= 0.001f;
		//x -= lx * fraction;
		//z -= lz * fraction;
		break;
	}
}
/* Handler for window re-size event. Called back when the window first appears and
whenever the window is re-sized with its new width and height */
void reshape(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
											   // Compute aspect ratio of the new window
	if (height == 0) height = 1;                // To prevent divide by 0
	GLfloat aspect = (GLfloat)width / (GLfloat)height;

	// Set the viewport to cover the new window
	glViewport(0, 0, width, height);

	// Set the aspect ratio of the clipping volume to match the viewport
	glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
	glLoadIdentity();             // Reset
								  // Enable perspective projection with fovy, aspect, zNear and zFar
	gluPerspective(45.0f, aspect, 0.1f, 100.0f);
}

/* Main function: GLUT runs as a console application starting at main() */
int main2(int argc, char** argv/*, std::vector<std::vector<double>> imageF*/, std::vector<std::vector<Cell>>& cvec) {

	//imagetopaint = imageF;
	cellstopaint = cvec;
// 	for (int i = 0; i < imagetopaint.size(); ++i)
// 		for (int j = 0; j < imagetopaint[0].size(); ++j)
// 			imagetopaint[i][j] *= 255.0;

	glutInit(&argc, argv);            // Initialize GLUT
	glutInitDisplayMode(GLUT_DOUBLE); // Enable double buffered mode
	glutInitWindowSize(640, 480);   // Set the window's initial width & height
	glutInitWindowPosition(50, 50); // Position the window's initial top-left corner
	glutCreateWindow(title);          // Create window with the given title
	
	glutDisplayFunc(display);       // Register callback handler for window re-paint event
	glutReshapeFunc(reshape);       // Register callback handler for window re-size event
	//glutSpecialFunc(processSpecialKeys);
	glutIdleFunc(display);

	glutIgnoreKeyRepeat(1);
	glutKeyboardFunc(processNormalKeys);
	glutSpecialFunc(pressKey);
	glutSpecialUpFunc(releaseKey);

	// регистрируем две новые функции
	glutMouseFunc(mouseButton);
	glutMotionFunc(mouseMove);

	initGL();                       // Our own OpenGL initialization
	glutMainLoop();                 // Enter the infinite event-processing loop
	return 0;
}