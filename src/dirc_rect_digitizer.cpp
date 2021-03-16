#include "dirc_point.h"
#include "dirc_rect_digitizer.h"
#include <vector>
DircRectDigitizer::DircRectDigitizer(\
		double iminx,\
		double imaxx,\
		double iresx,\
		double iminy,\
		double imaxy,\
		double iresy,\
		double it_unc,\
		double it_bin_size) : DircDigitizer()
{
//Digitizes to a rectagular grid
	minx = iminx;
	maxx = imaxx;
	resx = iresx;
	miny = iminy;
	maxy = imaxy;
	resy = iresy;
	t_unc = it_unc;
	t_bin_size = it_bin_size;
	
	dig_rand = new TRandom3();

}
void DircRectDigitizer::digitize_point(dirc_point &pt)
{
	/*
	int xdig = (x - minx)/resx;
	int ydig = (y - miny)/resy;
	
	double xout = resx*xdig + minx + resx/2;
	double yout = resy*ydig + miny + resy/2;
	
	if (x < minx)
	{
		xout = minx - resx/2;
	}
	else if (x > maxx)
	{
		xout = maxx + resx/2;
	}
	
	if (y < miny)
	{
		yout = miny - resy/2;
	}
	else if (y > maxy)
	{
		yout = maxy + resy/2;
	}
	
	pt.x = xout;
	pt.y = yout;
	pt.t += dig_rand->Gaus(0,t_unc);
	if (fabs(t_bin_size) > t_unc/100)
	{
		//Don't bother with binning if it's small
		//perhaps slow - a hard cut or a constant variable could be better here
		int tmp_t_bin = pt.t/(t_bin_size);
		pt.t = tmp_t_bin*t_bin_size + t_bin_size/2;
	}
	*/
	

	//Converting to as-built config
	double x = pt.x;
	double y = pt.y;
        int xhit = find_pixel(x);
        int yhit = find_pixel(y);

        //int yhit = (y - (-55.))/(PIXEL_SIZE + PIXEL_GAP);

	if (xhit == -1337 || yhit == -1337)
	{
		pt.ch = -999;
		pt.pixel_row = -999;
		pt.pixel_col = -999;
	}
	else
	{
		pt.pixel_row = xhit;
		pt.pixel_col = yhit;
		int pmt_id = (xhit/8) + (yhit/8) * PMT_ROWS;
		pt.ch = (pmt_id*CHANNEL_PER_PMT) + (yhit%8) * 8 + xhit%8;  
	}
	pt.x = find_coordinate(pt.pixel_row);
	pt.y = find_coordinate(pt.pixel_col);

}
int DircRectDigitizer::find_pixel(double x)
{
	int pmt_row    =  x / PMT_SIZE_GAP;
	double pixel_x = x - pmt_row*PMT_SIZE_GAP;
	int pmt_pixel_row = (pixel_x + PIXEL_GAP)/(PIXEL_SIZE+PIXEL_GAP);
	if (pixel_x > PMT_SIZE)
		return -1337;
	else
		return pmt_row*8 + pmt_pixel_row;
}

void DircRectDigitizer::undigitize_point(dirc_point &pt)
{
	pt.x = find_coordinate(pt.pixel_row);	
	pt.y = find_coordinate(pt.pixel_col);
}

double DircRectDigitizer::find_coordinate(int pixel_index)
{
	if (pixel_index < 0)
		return -999.;

	int pmt_i = pixel_index/8;
	int pmt_within_i = pixel_index%8;
	return pmt_i * PMT_SIZE_GAP + pmt_within_i * (PIXEL_SIZE + PIXEL_GAP) + 0.5 * PIXEL_SIZE;
}

void DircRectDigitizer::get_random_point(dirc_point &pt)
{

	pt.x = dig_rand->Rndm() * (PMT_SIZE_GAP * (PMT_ROWS - 1) + PMT_SIZE);
	pt.y = dig_rand->Rndm() * (PMT_SIZE_GAP * (PMT_COLUMNS - 1) + PMT_SIZE);
	pt.t = dig_rand->Rndm() * (350.);

	digitize_point(pt);
}

int
DircRectDigitizer::GetPmtID( int channel ) const 
{
	int MAX_BOX_CHANNEL = CHANNEL_PER_PMT * PMT_ROWS * PMT_COLUMNS;
	if(channel < MAX_BOX_CHANNEL) 
		return channel/CHANNEL_PER_PMT;
	else 
		return (channel-MAX_BOX_CHANNEL)/CHANNEL_PER_PMT;
}

int
DircRectDigitizer::GetPixelID( int channel ) const 
{
	return channel%CHANNEL_PER_PMT;
}

int DircRectDigitizer::GetPmtColumn( int channel ) const 
{
	int pmt = GetPmtID(channel);
	return pmt/PMT_ROWS; //  0 - 5
}

int DircRectDigitizer::GetPmtRow( int channel ) const 
{
	int pmt = GetPmtID(channel);
	return pmt%PMT_ROWS; //  0 - 17
}

int DircRectDigitizer::GetPmtPixelColumn( int channel ) const 
{
	int pix = GetPixelID(channel);
	return pix/8; //  0 - 7
}

int DircRectDigitizer::GetPmtPixelRow( int channel ) const 
{
	int pix = GetPixelID(channel);
	return pix%8; //  0 - 7
}

int DircRectDigitizer::GetPixelRow( int channel ) const 
{
	int pmt_row = GetPmtRow(channel);
	int pixel_row = GetPmtPixelRow(channel);
	return 8*pmt_row + pixel_row; //  0 - 143
}

int DircRectDigitizer::GetPixelColumn( int channel ) const 
{
	int pmt_column = GetPmtColumn(channel);
	int pixel_column = GetPmtPixelColumn(channel);
	return 8*pmt_column + pixel_column; // 0 - 47
}


void DircRectDigitizer::digitize_points(std::vector<dirc_point> &points)
{
	for (unsigned int i = 0; i < points.size(); i++)
	{
		//Wasting time Make it faster later TODO
		digitize_point(points[i]);
	}
}

void DircRectDigitizer::undigitize_points(std::vector<dirc_point> &points)
{
	for (unsigned int i = 0; i < points.size(); i++)
	{
		//Wasting time Make it faster later TODO
		undigitize_point(points[i]);
	}
}

	
