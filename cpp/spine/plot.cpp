/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#include "plot.h"

Plot * Plot::p = nullptr;

Plot::Plot(LPCWSTR name, unsigned int points)
{
	init(name, points);
}

Plot::Plot(LPCWSTR name, unsigned int points, double * xdata, double * ydata, BYTE r, BYTE g, BYTE b, double xmin, double xmax, double ymin, double ymax)
{
	init(name, points);
	add(r, g, b, xmin, xmax, ymin, ymax);
	PlotSeries * handle = get(0);
	for (unsigned int i = 0; i < points; i++)
		handle->add(i, ydata[i]);
	redraw();
}

Plot::Plot(LPCWSTR name, unsigned int points, double * ydata, BYTE r, BYTE g, BYTE b, double xmin, double xmax, double ymin, double ymax)
{
	init(name, points);
	add(r, g, b, xmin, xmax, ymin, ymax);
	PlotSeries * handle = get(0);
	for (unsigned int i = 0; i < points; i++)
		handle->add(i, ydata[i]);
	redraw();
}

void Plot::init(LPCWSTR name, unsigned int points)
{
	this->points = points;
	HINSTANCE hInstance = GetModuleHandle(NULL);

	WNDCLASS wc = { 0 };
	wc.lpfnWndProc = wndProc;
	wc.hInstance = hInstance;
	wc.hbrBackground = (HBRUSH)(COLOR_BACKGROUND);
	wc.lpszClassName = name;
	if (!RegisterClass(&wc))
		return;

	p = this;
	CreateWindow(wc.lpszClassName, name, WS_OVERLAPPEDWINDOW, 100, 100, 640, 480, 0, 0, hInstance, NULL);
	p = nullptr;

	if (handle) {
		ShowWindow(handle, SW_SHOWDEFAULT);
		UpdateWindow(handle);
	}
}

Plot::~Plot()
{
	for (vector<PlotSeries*>::iterator it = series.begin(); it != series.end(); ++it)
		delete *it;
}

int Plot::run(void)
{
	MSG msg = { 0 };
	while (GetMessage(&msg, NULL, WM_NULL, WM_NULL)) {
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}
	return (int)msg.wParam;
}

void Plot::runDraw()
{
	MSG msg = { 0 };
	while (PeekMessage(&msg, NULL, WM_NULL, WM_NULL, PM_REMOVE)) {
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}
}

void Plot::redraw(bool needs_clearing)
{
	if (needs_clearing)
		for (vector<PlotSeries*>::iterator it = series.begin(); it != series.end(); ++it)
			(*it)->clear();
	InvalidateRect(handle, NULL, TRUE);
	runDraw();
}

void Plot::add(BYTE r, BYTE g, BYTE b, double xmin, double xmax, double ymin, double ymax)
{
	PlotSeries * series = new PlotSeries(points, r, g, b, xmin, xmax, ymin, ymax);
	this->series.push_back(series);
}

PlotSeries * Plot::get(unsigned int n)
{
	return series[n];
}

void Plot::OnPaint(HDC hdc)
{
	HBRUSH brush = CreateSolidBrush(RGB(255, 255, 255));
	HBRUSH brushOld = (HBRUSH)SelectObject(hdc, brush);

	RECT clientRect;
	GetClientRect(handle, &clientRect);

	LONG plotxmin = clientRect.left + 10;
	LONG plotxmax = clientRect.right - 10;
	LONG plotymin = clientRect.top + 10;
	LONG plotymax = clientRect.bottom - 10;

	bool needs_clearing = false;
	for (vector<PlotSeries*>::iterator it = series.begin(); it != series.end(); ++it) {
		if ((*it)->needsClearing()) {
			needs_clearing = true;
			break;
		}
	}
	if (needs_clearing)
		Rectangle(hdc, plotxmin, plotymin, plotxmax, plotymax);
	for (vector<PlotSeries*>::iterator it = series.begin(); it != series.end(); ++it)
		(*it)->draw(hdc, plotxmin, plotxmax, plotymin, plotymax);

	SelectObject(hdc, brushOld);
	DeleteObject(brush);
}

LRESULT CALLBACK Plot::wndProc(_In_ HWND hWnd, _In_ UINT uMsg, _In_ WPARAM wParam, _In_ LPARAM lParam)
{
	Plot * p = (Plot*)(GetWindowLongPtr(hWnd, GWLP_USERDATA));

	if (p)
		return p->wndProcPlot(hWnd, uMsg, wParam, lParam);
	else if (Plot::p)
		return Plot::p->wndProcPlot(hWnd, uMsg, wParam, lParam);

	return DefWindowProc(hWnd, uMsg, wParam, lParam);
}

INT_PTR Plot::wndProcPlot(_In_ HWND hWnd, _In_ UINT uMsg, _In_ WPARAM wParam, _In_ LPARAM lParam)
{
	HDC          hdc;
	PAINTSTRUCT  ps;

	switch (uMsg)
	{
	case WM_CREATE:
		handle = hWnd;
		SetWindowLongPtr(hWnd, GWLP_USERDATA, (LONG_PTR)this);
		return 0;

	case WM_ERASEBKGND:
		return (LRESULT)1;

	case WM_PAINT:
		hdc = BeginPaint(hWnd, &ps);
		OnPaint(hdc);
		EndPaint(hWnd, &ps);
		return 0;

	case WM_CLOSE:
		PostQuitMessage(0);
		return 0;

	case WM_SIZE:
		redraw(true);
		return 0;

	default:
		return DefWindowProc(hWnd, uMsg, wParam, lParam);
	}
}

PlotSeries::PlotSeries(unsigned int points, BYTE r, BYTE g, BYTE b, double xmin, double xmax, double ymin, double ymax)
{
	this->r = r;
	this->g = g;
	this->b = b;
	x.reserve(points);
	y.reserve(points);
	this->points = 0;
	this->xmin = xmin;
	this->xmax = xmax;
	this->ymin = ymin;
	this->ymax = ymax;
	this->drawn_upto = 0;
	this->needs_clearing = true;
}

void PlotSeries::clear()
{
	needs_clearing = true;
}

bool PlotSeries::needsClearing()
{
	return needs_clearing;
}

void PlotSeries::add(double x, double y)
{
	if ((points == 0) && (xmin == 0) && (xmax == 0)) {
		xmin = x;
		xmax = x;
		needs_clearing = true;
	}
	if (x > xmax)
	{
		xmax = x;
		needs_clearing = true;
	}
	points++;
	this->x.push_back(x);
	this->y.push_back(y);
}

void PlotSeries::draw(HDC hdc, LONG plotxmin, LONG plotxmax, LONG plotymin, LONG plotymax)
{
	HPEN pen = CreatePen(PS_SOLID, 1, RGB(r, g, b));
	HPEN penOld = (HPEN)SelectObject(hdc, pen);

	if (needs_clearing)
		drawn_upto = 0;

	if (points)
	{
		LONG xcoor = (LONG) ((x[drawn_upto] - xmin) / (xmax - xmin) * (plotxmax - plotxmin) + plotxmin);
		LONG ycoor = (LONG) (plotymax - (y[drawn_upto] - ymin) / (ymax - ymin) * (plotymax - plotymin));
		MoveToEx(hdc, xcoor, ycoor, NULL);
		for (unsigned int i = drawn_upto + 1; i < points; i++) {
			xcoor = (LONG)((x[i] - xmin) / (xmax - xmin) * (plotxmax - plotxmin) + plotxmin);
			ycoor = (LONG)(plotymax - (y[i] - ymin) / (ymax - ymin) * (plotymax - plotymin));
			LineTo(hdc, xcoor, ycoor);
		}
		drawn_upto = points - 1;
	}
	needs_clearing = false;

	SelectObject(hdc, penOld);
	DeleteObject(pen);
}

void PlotSeries::xlim(double xmin, double xmax)
{
	this->xmin = xmin;
	this->xmax = xmax;
	this->needs_clearing = true;
}

void PlotSeries::ylim(double ymin, double ymax)
{
	this->ymin = ymin;
	this->ymax = ymax;
	this->needs_clearing = true;
}
