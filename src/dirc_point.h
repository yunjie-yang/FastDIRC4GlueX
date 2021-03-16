
#ifndef DIRC_POINT
#define DIRC_POINT

/*
struct dirc_point
{
	double x;
	double y;
	double t;
	int updown;
	int last_wall_x;
	int wedge_before_interface;
	double weight;
	double init_phi; //internal validation only
};
*/
class dirc_point
{
public:
        double x;
        double y;
        double t;
        int updown;
        int last_wall_x;
        int wedge_before_interface;
        double weight;
        double init_phi; //internal validation only
	int ch;
	int pixel_row;
	int pixel_col;
	bool cut;
	double min_dist[5];
	double frac_close[5];
	bool direct;
};

#endif
