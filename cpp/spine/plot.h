/***************************************************************
* Author: J.P.G. van Dijk
***************************************************************/

#pragma once

#include <Windows.h>
#include <vector>

using namespace std;

class PlotSeries
{
private:
	unsigned int points;
	unsigned int drawn_upto;
	bool needs_clearing;
	double xmin, xmax, ymin, ymax;
	vector<double> x;
	vector<double> y;
	BYTE r, g, b;

public:
	PlotSeries(unsigned int points = 100, BYTE r = 0, BYTE g = 0, BYTE b = 0, double xmin = 0, double xmax = 0, double ymin = 0, double ymax = 1);
	void clear();
	bool needsClearing();
	void add(double x, double y);
	void draw(HDC hdc, LONG plotxmin, LONG plotxmax, LONG plotymin, LONG plotymax);
	void xlim(double xmin, double xmax);
	void ylim(double ymin, double ymax);
};

class Plot
{
private:
	unsigned int points;
	vector<PlotSeries*> series;

public:
	Plot(LPCWSTR name, unsigned int points = 100);
	Plot(LPCWSTR name, unsigned int points, double * xdata, double * ydata, BYTE r, BYTE g, BYTE b, double xmin, double xmax, double ymin, double ymax);
	Plot(LPCWSTR name, unsigned int points, double * ydata, BYTE r, BYTE g, BYTE b, double xmin, double xmax, double ymin, double ymax);
	~Plot();
	static int run(void);
	void runDraw();
	void redraw(bool needs_clearing = false);
	void add(BYTE r = 0, BYTE g = 0, BYTE b = 0, double xmin = 0, double xmax = 0, double ymin = 0, double ymax = 1);
	PlotSeries * get(unsigned int n);

private:
	void init(LPCWSTR name, unsigned int points = 100);
	void OnPaint(HDC hdc);

	// WIN32
public:
	static Plot * p;
	static LRESULT CALLBACK wndProc(_In_ HWND hWnd, _In_ UINT uMsg, _In_ WPARAM wParam, _In_ LPARAM lParam);
private:
	HWND handle = nullptr;
	INT_PTR wndProcPlot(_In_ HWND hWnd, _In_ UINT uMsg, _In_ WPARAM wParam, _In_ LPARAM lParam);
};
