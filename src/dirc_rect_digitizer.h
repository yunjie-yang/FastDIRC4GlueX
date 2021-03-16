#include "dirc_point.h"
#include "dirc_digitizer.h"
#include <vector>
#include <iostream>

//Root is infecting more files :(
#include <TRandom3.h>

#ifndef DIRC_RECT_DIGITIZER
#define DIRC_RECT_DIGITIZER
class DircRectDigitizer : public DircDigitizer
{
private:
	double minx,maxx,miny,maxy;
	double resx, resy;
	double t_unc;
	double t_bin_size;
	
	TRandom3 *dig_rand;
        int CHANNEL_PER_PMT = 64;
        int PMT_ROWS = 18;
        int PMT_COLUMNS = 6;
	
	double PIXEL_SIZE   = 6;
	double PIXEL_GAP    = 0.06;
	double PMT_GAP      = 4.58;
	double PMT_SIZE     = PIXEL_SIZE * 8 + PIXEL_GAP * 7;
	double PMT_SIZE_GAP = PMT_SIZE + PMT_GAP;

	int find_pixel(double x);
	double find_coordinate(int pixel_index);
	
public:
	DircRectDigitizer(\
		double iminx,\
		double imaxx,\
		double iresx,\
		double iminy,\
		double imaxy,\
		double iresy,\
		double it_unc,\
		double it_bin_size);

	int GetPmtID      ( int channel ) const;
	int GetPmtRow     ( int channel ) const;
	int GetPmtColumn  ( int channel ) const;
	int GetPixelID      ( int channel ) const;
	int GetPmtPixelRow     ( int channel ) const;
	int GetPmtPixelColumn  ( int channel ) const;
	int GetPixelRow     ( int channel ) const;
	int GetPixelColumn  ( int channel ) const;

	void digitize_point(dirc_point &pt);
	void digitize_points(std::vector<dirc_point> &points);

	void undigitize_point(dirc_point &pt);
	void undigitize_points(std::vector<dirc_point> &points);
	
	void get_random_point(dirc_point &pt);
};
#endif
