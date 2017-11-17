#include <windows.h>  // Заголовочный файл для Windows
#include <stdio.h>    // Заголовочный файл для стандартного ввода-вывода (НОВОЕ)
#include <gl\gl.h>    // Заголовочный файл для OpenGL32 библиотеки
#include <gl\glu.h>   // Заголовочный файл для GLu32 библиотеки
#include <gl\glaux.h> // Заголовочный файл для GLaux библиотеки

#pragma comment(lib, "opengl32.lib")     // Ссылка на OpenGL32.lib
#pragma comment(lib, "glu32.lib")        // Ссылка на Glu32.lib

#define    MAP_SIZE  1024     // Размер карты вершин (НОВОЕ)
#define    STEP_SIZE  16      // Ширина и высота каждого квадрата (НОВОЕ)
// Коэффициент масштабирования по оси Y в соответствии с осями X и Z (НОВОЕ)
#define    HEIGHT_RATIO  1.5f

HDC    hDC = NULL;        // Приватный контекст устройства GDI
HGLRC    hRC = NULL;      // Постоянный контекст рендеринга
HWND    hWnd = NULL;      // Указатель на наше окно
HINSTANCE  hInstance;   // Указывает на дескриптор текущего приложения

bool    keys[256];        // Массив для процедуры обработки клавиатуры
bool    active = TRUE;      // Флаг активности окна, по умолчанию=TRUE
bool    fullscreen = TRUE;  // Флаг полноэкранного режима, по умолчанию=TRUE
bool    bRender = TRUE;   // Флаг режима отображения полигонов,
						  // по умолчанию=TRUE (НОВОЕ)

BYTE g_HeightMap[MAP_SIZE*MAP_SIZE];    // Содержит карту вершин (НОВОЕ)

float scaleValue = 0.15f;               // Величина масштабирования поверхности (НОВОЕ)

LRESULT  CALLBACK WndProc(HWND, UINT, WPARAM, LPARAM); // Объявление для WndProc

GLvoid ReSizeGLScene(GLsizei width, GLsizei height)        // Изменить размер и инициализировать окно GL
{
	if (height == 0)              // Предотвращение деления на ноль
	{
		height = 1;
	}

	glViewport(0, 0, width, height);          // Сброс текущей области вывода
	glMatrixMode(GL_PROJECTION);            // Выбор матрицы проекций
	glLoadIdentity();              // Сброс матрицы проекции

								   // Вычисление соотношения геометрических размеров для окна
	gluPerspective(45.0f, (GLfloat)width / (GLfloat)height, 0.1f, 100.0f);

	glMatrixMode(GL_MODELVIEW);            // Выбор матрицы вида модели
	glLoadIdentity();              // Сброс матрицы вида модели
}

// Чтение и сохранение .RAW файла в pHeightMap
void LoadRawFile(LPSTR strName, int nSize, BYTE *pHeightMap)
{
// 	FILE *pFile = NULL;
// 
// 	// открытие файла в режиме бинарного чтения
// 	pFile = fopen(strName, "rb");
// 
// 	// Файл найден?
// 	if (pFile == NULL)
// 	{
// 		// Выводим сообщение об ошибке и выходим из процедуры
// 		MessageBox(NULL, "Can't Find The Height Map!", "Error", MB_OK);
// 		return;
// 	}
// 	// Загружаем .RAW файл в массив pHeightMap
// 	// Каждый раз читаем по одному байту, размер = ширина * высота
// 	fread(pHeightMap, 1, nSize, pFile);
// 
// 	// Проверяем на наличие ошибки
// 	int result = ferror(pFile);
// 
// 	// Если произошла ошибка
// 	if (result)
// 	{
// 		MessageBox(NULL, "Failed To Get Data!", "Error", MB_OK);
// 	}
// 
// 	// Закрываем файл
// 	fclose(pFile);
}

int InitGL(GLvoid)           // Инициализация OpenGL
{
	glShadeModel(GL_SMOOTH);   // Включить сглаживание
	glClearColor(0.0f, 0.0f, 0.0f, 0.5f);  // Очистка экрана черным цветом
	glClearDepth(1.0f);        // Установка буфера глубины
	glEnable(GL_DEPTH_TEST);   // Включить буфер глубины
	glDepthFunc(GL_LEQUAL);    // Тип теста глубины
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);  // Улучшенные вычисления перспективы

	// Читаем данные из файла и сохраняем их в массиве g_HeightMap array.
	// Также передаем размер файла (1024).

	LoadRawFile("Data/Terrain.raw", MAP_SIZE * MAP_SIZE, g_HeightMap);  // ( НОВОЕ )

	return TRUE;               // Инициализация прошла успешно
}

int Height(BYTE *pHeightMap, int X, int Y)      // Возвращает высоту из карты вершин (?)
{
	int x = X % MAP_SIZE;          // Проверка переменной х
	int y = Y % MAP_SIZE;          // Проверка переменной y

	if (!pHeightMap) return 0;      // Убедимся, что данные корректны
	return pHeightMap[x + (y * MAP_SIZE)];      // Возвращаем значение высоты
}
// Эта функция устанавливает значение цвета для конкретного номера, зависящего от номера высоты
void SetVertexColor(BYTE *pHeightMap, int x, int y)
{
	if (!pHeightMap) return;          // Данные корректны?

	float fColor = -0.15f + (Height(pHeightMap, x, y) / 256.0f);

	// Присвоить оттенок синего цвета для конкретной точки
	glColor3f(0.0f, 0.0f, fColor);
}
void RenderHeightMap(BYTE pHeightMap[]) // Визуализация карты высоты с помощью квадратов
{
	int X = 0, Y = 0;    // Создаем пару переменных для перемещения по массиву
	int x, y, z;         // И еще три для удобства чтения

	if (!pHeightMap) return;     // Данные корректны?
	if (bRender)            // Что хотим визуализировать?
		glBegin(GL_QUADS); // Полигоны
	else
		glBegin(GL_LINES); // Линии
	for (X = 0; X < MAP_SIZE; X += STEP_SIZE)
		for (Y = 0; Y < MAP_SIZE; Y += STEP_SIZE)
		{
			// Получаем (X, Y, Z) координаты нижней левой вершины
			x = X;
			y = Height(pHeightMap, X, Y);
			z = Y;

			// Устанавливаем цвет конкретной точки
			SetVertexColor(pHeightMap, x, z);

			glVertex3i(x, y, z);      // Визуализация ее

									  // Получаем (X, Y, Z) координаты верхней левой вершины
			x = X;
			y = Height(pHeightMap, X, Y + STEP_SIZE);
			z = Y + STEP_SIZE;

			// Устанавливаем цвет конкретной точки
			SetVertexColor(pHeightMap, x, z);

			glVertex3i(x, y, z);      // Визуализация ее

									  // Получаем (X, Y, Z) координаты верхней правой вершины
			x = X + STEP_SIZE;
			y = Height(pHeightMap, X + STEP_SIZE, Y + STEP_SIZE);
			z = Y + STEP_SIZE;

			// Устанавливаем цвет конкретной точки
			SetVertexColor(pHeightMap, x, z);

			glVertex3i(x, y, z);      // Визуализация ее

									  // Получаем (X, Y, Z) координаты нижней правой вершины
			x = X + STEP_SIZE;
			y = Height(pHeightMap, X + STEP_SIZE, Y);
			z = Y;

			// Устанавливаем цвет конкретной точки
			SetVertexColor(pHeightMap, x, z);

			glVertex3i(x, y, z);      // Визуализация ее
		}
	glEnd();
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);      // Сбрасываем цвет
}
int DrawGLScene(GLvoid)      // Здесь содержится код рисования
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  // Очистка экрана и буфера глубины
	glLoadIdentity();          // Сброс просмотра

							   //       Положение   Вид    Вектор вертикали
	gluLookAt(212, 60, 194, 186, 55, 171, 0, 1, 0);  // Определяет вид и положение камеры
	glScalef(scaleValue, scaleValue * HEIGHT_RATIO, scaleValue);
	RenderHeightMap(g_HeightMap);        // Визализация карты высот

	return TRUE;            // Идем дальше
}
GLvoid KillGLWindow(GLvoid)              // Корректное разрушение окна
{
	if (fullscreen)              // Мы в полноэкранном режиме?
	{
		ChangeDisplaySettings(NULL, 0);        // Если да, то переключаемся обратно в оконный режим
		ShowCursor(true);            // Показать курсор мышки
	}
	if (hRC)                // Существует ли Контекст Рендеринга?
	{
		if (!wglMakeCurrent(NULL, NULL))        // Возможно ли освободить RC и DC?
		{
			MessageBox(NULL, "Release Of DC And RC Failed.", "SHUTDOWN ERROR", MB_OK | MB_ICONINFORMATION);
		}
		if (!wglDeleteContext(hRC))        // Возможно ли удалить RC?
		{
			MessageBox(NULL, "Release Rendering Context Failed.", "SHUTDOWN ERROR", MB_OK | MB_ICONINFORMATION);
		}
		hRC = NULL;              // Установить RC в NULL
	}
	if (hDC && !ReleaseDC(hWnd, hDC))          // Возможно ли уничтожить DC?
	{
		MessageBox(NULL, "Release Device Context Failed.", "SHUTDOWN ERROR", MB_OK | MB_ICONINFORMATION);
		hDC = NULL;                // Установить DC в NULL
	}
	if (hWnd && !DestroyWindow(hWnd))            // Возможно ли уничтожить окно?
	{
		MessageBox(NULL, "Could Not Release hWnd.", "SHUTDOWN ERROR", MB_OK | MB_ICONINFORMATION);
		hWnd = NULL;                // Установить hWnd в NULL
	}
	if (!UnregisterClass("OpenGL", hInstance))        // Возможно ли разрегистрировать класс
	{
		MessageBox(NULL, "Could Not Unregister Class.", "SHUTDOWN ERROR", MB_OK | MB_ICONINFORMATION);
		hInstance = NULL;                // Установить hInstance в NULL
	}
}
BOOL CreateGLWindow(LPCWSTR title, int width, int height, int bits, bool fullscreenflag)
{
	GLuint    PixelFormat;              // Хранит результат после поиска
	WNDCLASS  wc;                // Структура класса окна
	DWORD    dwExStyle;              // Расширенный стиль окна
	DWORD    dwStyle;              // Обычный стиль окна
	RECT WindowRect;                // Grabs Rectangle Upper Left / Lower Right Values
	WindowRect.left = (long)0;              // Установить левую составляющую в 0
	WindowRect.right = (long)width;              // Установить правую составляющую в Width
	WindowRect.top = (long)0;                // Установить верхнюю составляющую в 0
	WindowRect.bottom = (long)height;              // Установить нижнюю составляющую в Height
	fullscreen = fullscreenflag;              // Устанавливаем значение глобальной переменной fullscreen
	hInstance = GetModuleHandle(NULL);        // Считаем дескриптор нашего приложения
	wc.style = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;      // Перерисуем при перемещении и создаём скрытый DC
	wc.lpfnWndProc = (WNDPROC)WndProc;          // Процедура обработки сообщений
	wc.cbClsExtra = 0;              // Нет дополнительной информации для окна
	wc.cbWndExtra = 0;              // Нет дополнительной информации для окна
	wc.hInstance = hInstance;            // Устанавливаем дескриптор
	wc.hIcon = LoadIcon(NULL, IDI_WINLOGO);        // Загружаем иконку по умолчанию
	wc.hCursor = LoadCursor(NULL, IDC_ARROW);        // Загружаем указатель мышки
	wc.hbrBackground = NULL;              // Фон не требуется для GL
	wc.lpszMenuName = NULL;              // Меню в окне не будет
	wc.lpszClassName = "OpenGL";            // Устанавливаем имя классу
	if (!RegisterClass(&wc))              // Пытаемся зарегистрировать класс окна
	{
		MessageBox(NULL, "Failed To Register The Window Class.", "ERROR", MB_OK | MB_ICONEXCLAMATION);
		return false;                // Выход и возвращение функцией значения false
	}
	if (fullscreen)                // Полноэкранный режим?
	{
		DEVMODE dmScreenSettings;            // Режим устройства
		memset(&dmScreenSettings, 0, sizeof(dmScreenSettings));    // Очистка для хранения установок
		dmScreenSettings.dmSize = sizeof(dmScreenSettings);      // Размер структуры Devmode
		dmScreenSettings.dmPelsWidth = width;        // Ширина экрана
		dmScreenSettings.dmPelsHeight = height;        // Высота экрана
		dmScreenSettings.dmBitsPerPel = bits;        // Глубина цвета
		dmScreenSettings.dmFields = DM_BITSPERPEL | DM_PELSWIDTH | DM_PELSHEIGHT;// Режим Пикселя
																				 // Пытаемся установить выбранный режим и получить результат.  Примечание: CDS_FULLSCREEN убирает панель управления.
		if (ChangeDisplaySettings(&dmScreenSettings, CDS_FULLSCREEN) != DISP_CHANGE_SUCCESSFUL)
		{
			// Если переключение в полноэкранный режим невозможно, будет предложено два варианта: оконный режим или выход.
			if (MessageBox(NULL, "The Requested Fullscreen Mode Is Not Supported By\nYour Video Card. Use Windowed Mode Instead?",
				"NeHe GL", MB_YESNO | MB_ICONEXCLAMATION) == IDYES)
			{
				fullscreen = false;          // Выбор оконного режима (fullscreen = false)
			}
			else
			{
				// Выскакивающее окно, сообщающее пользователю о закрытие окна.
				MessageBox(NULL, "Program Will Now Close.", "ERROR", MB_OK | MB_ICONSTOP);
				return false;            // Выход и возвращение функцией false
			}
		}
	}
	if (fullscreen)                  // Мы остались в полноэкранном режиме?
	{
		dwExStyle = WS_EX_APPWINDOW;          // Расширенный стиль окна
		dwStyle = WS_POPUP;            // Обычный стиль окна
		ShowCursor(false);              // Скрыть указатель мышки
	}
	else
	{
		dwExStyle = WS_EX_APPWINDOW | WS_EX_WINDOWEDGE;      // Расширенный стиль окна
		dwStyle = WS_OVERLAPPEDWINDOW;        // Обычный стиль окна
	}
	AdjustWindowRectEx(&WindowRect, dwStyle, false, dwExStyle);      // Подбирает окну подходящие размеры
	if (!(hWnd = CreateWindowEx(dwExStyle,          // Расширенный стиль для окна
		_T("OpenGL"),          // Имя класса
		title,            // Заголовок окна
		WS_CLIPSIBLINGS |        // Требуемый стиль для окна
		WS_CLIPCHILDREN |        // Требуемый стиль для окна
		dwStyle,          // Выбираемые стили для окна
		0, 0,            // Позиция окна
		WindowRect.right - WindowRect.left,    // Вычисление подходящей ширины
		WindowRect.bottom - WindowRect.top,    // Вычисление подходящей высоты
		NULL,            // Нет родительского
		NULL,            // Нет меню
		hInstance,          // Дескриптор приложения
		NULL)))          // Не передаём ничего до WM_CREATE (???)
	{
		KillGLWindow();                // Восстановить экран
		MessageBox(NULL, "Window Creation Error.", "ERROR", MB_OK | MB_ICONEXCLAMATION);
		return false;                // Вернуть false
	}
	static  PIXELFORMATDESCRIPTOR pfd =            // pfd сообщает Windows каким будет вывод на экран каждого пикселя
	{
		sizeof(PIXELFORMATDESCRIPTOR),            // Размер дескриптора данного формата пикселей
		1,                  // Номер версии
		PFD_DRAW_TO_WINDOW |              // Формат для Окна
		PFD_SUPPORT_OPENGL |              // Формат для OpenGL
		PFD_DOUBLEBUFFER,              // Формат для двойного буфера
		PFD_TYPE_RGBA,                // Требуется RGBA формат
		bits,                  // Выбирается бит глубины цвета
		0, 0, 0, 0, 0, 0,              // Игнорирование цветовых битов
		0,                  // Нет буфера прозрачности
		0,                  // Сдвиговый бит игнорируется
		0,                  // Нет буфера накопления
		0, 0, 0, 0,                // Биты накопления игнорируются
		32,                  // 32 битный Z-буфер (буфер глубины)
		0,                  // Нет буфера трафарета
		0,                  // Нет вспомогательных буферов
		PFD_MAIN_PLANE,                // Главный слой рисования
		0,                  // Зарезервировано
		0, 0, 0                  // Маски слоя игнорируются
	};
	if (!(hDC = GetDC(hWnd)))              // Можем ли мы получить Контекст Устройства?
	{
		KillGLWindow();                // Восстановить экран
		MessageBox(NULL, "Can't Create A GL Device Context.", "ERROR", MB_OK | MB_ICONEXCLAMATION);
		return false;                // Вернуть false
	}
	if (!(PixelFormat = ChoosePixelFormat(hDC, &pfd)))        // Найден ли подходящий формат пикселя?
	{
		KillGLWindow();                // Восстановить экран
		MessageBox(NULL, "Can't Find A Suitable PixelFormat.", "ERROR", MB_OK | MB_ICONEXCLAMATION);
		return false;                // Вернуть false
	}
	if (!SetPixelFormat(hDC, PixelFormat, &pfd))          // Возможно ли установить Формат Пикселя?
	{
		KillGLWindow();                // Восстановить экран
		MessageBox(NULL, "Can't Set The PixelFormat.", "ERROR", MB_OK | MB_ICONEXCLAMATION);
		return false;                // Вернуть false
	}
	if (!(hRC = wglCreateContext(hDC)))          // Возможно ли установить Контекст Рендеринга?
	{
		KillGLWindow();                // Восстановить экран
		MessageBox(NULL, "Can't Create A GL Rendering Context.", "ERROR", MB_OK | MB_ICONEXCLAMATION);
		return false;                // Вернуть false
	}
	if (!wglMakeCurrent(hDC, hRC))            // Попробовать активировать Контекст Рендеринга
	{
		KillGLWindow();                // Восстановить экран
		MessageBox(NULL, "Can't Activate The GL Rendering Context.", "ERROR", MB_OK | MB_ICONEXCLAMATION);
		return false;                // Вернуть false
	}
	ShowWindow(hWnd, SW_SHOW);              // Показать окно
	SetForegroundWindow(hWnd);              // Слегка повысим приоритет
	SetFocus(hWnd);                // Установить фокус клавиатуры на наше окно
	ReSizeGLScene(width, height);              // Настроим перспективу для нашего OpenGL экрана.
	if (!InitGL())                  // Инициализация только что созданного окна
	{
		KillGLWindow();                // Восстановить экран
		MessageBox(NULL, _T("Initialization Failed."), _T("ERROR"), MB_OK | MB_ICONEXCLAMATION);
		return false;                // Вернуть false
	}
	return true;                  // Всё в порядке!
}
LRESULT CALLBACK WndProc(
	HWND  hWnd,      // Указатель на окно
	UINT  uMsg,      // Сообщение для этого окна
	WPARAM  wParam,  // Параметры сообщения
	LPARAM  lParam)  // Параметры сообщения
{
	switch (uMsg)          // Проверим сообщения окна
	{
	case WM_ACTIVATE:    // Наблюдаем за сообщением об активизации окна
	{
		if (!HIWORD(wParam)) // Проверим состояние минимизации
		{
			active = TRUE;     // Программа активна
		}
		else
		{
			active = FALSE;    // Программа больше не активна
		}

		return 0;          // Вернуться к циклу сообщений
	}

	case WM_SYSCOMMAND:  // Перехватаем системную команду
	{
		switch (wParam)    // Проверка выбора системы
		{
		case SC_SCREENSAVE:   // Пытается включиться скринсэйвер?
		case SC_MONITORPOWER: // Монитор пытается переключиться в режим сохранения энергии?
			return 0;        // Не давать этому случиться
		}
		break;             // Выход
	}

	case WM_CLOSE:       // Мы получили сообщение о закрытии программы?
	{
		PostQuitMessage(0);// Послать сообщение о выходе
		return 0;          // Возврат обратно
	}

	case WM_LBUTTONDOWN: // Нажата левая клавиша мыши?
	{
		bRender = !bRender;// Изменить тип визуализации
		return 0;          // Вернуться
	}

	case WM_KEYDOWN:     // Клавиша была нажата?
	{
		keys[wParam] = TRUE; // Если так – пометим это TRUE
		return 0;          // Вернуться
	}

	case WM_KEYUP:       // Клавиша была отпущена?
	{
		keys[wParam] = FALSE; // Если так – пометим это FALSE
		return 0;          // Вернуться
	}

	case WM_SIZE:        // Изменились окна OpenGL
	{
		ReSizeGLScene(LOWORD(lParam), HIWORD(lParam));  // LoWord=ширина, HiWord=высота
		return 0;          // Вернуться
	}
	}
	// Пересылаем все прочие сообщения в DefWindowProc
	return DefWindowProc(hWnd, uMsg, wParam, lParam);
}
int WINAPI WinMain(
	HINSTANCE  hInstance,    // Экземпляр
	HINSTANCE  hPrevInstance,// Предыдущий экземпляр
	LPSTR    lpCmdLine,      // Параметры командной строки
	int    nCmdShow)         // Показать состояние окна
{
	MSG    msg;          // Структура сообщения окна
	BOOL  done = FALSE;    // Булевская переменная выхода из цикла

						   // Запросим пользователя, какой режим отображения он предпочитает
	if (MessageBox(NULL, "Would You Like To Run In Fullscreen Mode?",
		"Start FullScreen?", MB_YESNO | MB_ICONQUESTION) == IDNO)
	{
		fullscreen = FALSE;          // Оконный режим
	}

	// Создадим наше окно OpenGL
	if (!CreateGLWindow("NeHe & Ben Humphrey's Height Map Tutorial",
		640, 480, 16, fullscreen))
	{
		return 0;                  // Выходим если окно не было создано
	}

	while (!done)                 // Цикл, который продолжается пока done=FALSE
	{
		if (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))  // Есть ожидаемое сообщение?
		{
			if (msg.message == WM_QUIT)// Мы получили сообщение о выходе?
			{
				done = TRUE;             // Если так done=TRUE
			}
			else                     // Если нет, продолжаем работать с сообщениями окна
			{
				TranslateMessage(&msg);// Переводим сообщение
				DispatchMessage(&msg); // Отсылаем сообщение
			}
		}
		else                       // Если сообщений нет
		{
			// Рисуем сцену. Ожидаем нажатия кнопки ESC и сообщения о выходе от DrawGLScene()
			// Активно?  Было получено сообщение о выходе?
			if ((active && !DrawGLScene()) || keys[VK_ESCAPE])
			{
				done = TRUE;             // ESC или DrawGLScene просигналили "выход"
			}
			else if (active)         // Не время выходить, обновляем экран
			{
				SwapBuffers(hDC);      // Переключаем буферы (Двойная буферизация)
			}

			if (keys[VK_F1])         // Была нажата кнопка F1?
			{
				keys[VK_F1] = FALSE;     // Если так - установим значение FALSE
				KillGLWindow();        // Закроем текущее окно OpenGL 
				fullscreen = !fullscreen;// Переключим режим "Полный экран"/"Оконный"
										 // Заново создадим наше окно OpenGL
				if (!CreateGLWindow("NeHe & Ben Humphrey's Height Map Tutorial",
					640, 480, 16, fullscreen))
				{
					return 0;            // Выйти, если окно не было создано
				}
			}
			if (keys[VK_UP])           // Нажата клавиша ВВЕРХ?
				scaleValue += 0.001f;  // Увеличить переменную масштабирования

			if (keys[VK_DOWN])       // Нажата клавиша ВНИЗ?
				scaleValue -= 0.001f;  // Уменьшить переменную масштабирования
		}
	}

	// Shutdown
	KillGLWindow();              // Закроем окно
	return (msg.wParam);         // Выйдем из программы
} 