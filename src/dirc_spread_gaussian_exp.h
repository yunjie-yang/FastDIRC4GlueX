#include "dirc_point.h"
#include "dirc_spread_radius.h"
#include <vector>
//#include <math.h>
#include <TRandom3.h>

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#ifndef DIRC_SPREAD_GAUSSIAN_EXP
#define DIRC_SPREAD_GAUSSIAN_EXP
class DircSpreadGaussianExp
{
private:
	double x_sig2inv,y_sig2inv,t_sig2inv;
	double spread_func_norm, spread_func_norm_inv;
	double lin_slope, r_trans, sigma2, sigma2inv,max_val;
	double min_probability;
	
	double sigma_cut2, prob_cut, prob_cut_per_point;
	double exp_norm_factor;
	double r2_cut;
	
	TFile* exp_hist_file;
	TH1D*  exp_hist;
	int bin_multiple;
	char exp_hist_filename[256];

	double norm_factor_xy;
public:
	DircSpreadGaussianExp(\
		double isigma, \
		double x_unc,\
		double y_unc,\
		double t_unc,\
		double isigma_cut = 5.,\
		double imin_prob  = 1e-6,\
		const char* iexp_hist_filename = "./cached_exp.root");
	void set_gaus_sigma(double isigma);

/*
	inline double radius_spread_function(double r2) __attribute__((always_inline))
	{
		//if (r2 < 5*sigma2)
		if (r2 < 25*sigma2)
		{
			return exp(-r2*sigma2inv);
		}
		else
		{
			return 0;
		}
	};
	inline double support_spread_function(dirc_point support, dirc_point test)__attribute__((always_inline))
	{
		double dx2,dy2,dt2;
		dx2 = support.x - test.x;
		dx2 *= dx2;
		dy2 = support.y - test.y;
		dy2 *= dy2;
		dt2 = support.t - test.t;
		dt2 *= dt2;
		return radius_spread_function(dx2*x_sig2inv+dy2*y_sig2inv+dt2*t_sig2inv);
	};
*/
	inline double support_spread_function(dirc_point support, dirc_point test, int &close_point)__attribute__((always_inline))
	{
		double dx2,dy2,dt2,r2;
		dx2 = support.x - test.x;
		dx2 *= dx2;
		dy2 = support.y - test.y;
		dy2 *= dy2;
		dt2 = support.t - test.t;
		dt2 *= dt2;
		r2 = dx2*x_sig2inv+dy2*y_sig2inv+dt2*t_sig2inv;
		if (r2 < sigma_cut2 * sigma2)
		{
			close_point++;
			return exp(-r2*sigma2inv);
		}
		else
		{
			return 0.;
		}
	};
	inline double support_spread_function_speedup(dirc_point support, dirc_point test, int &close_point)__attribute__((always_inline))
	{
		double dx2,dy2,dt2,r2;
		dx2 = support.x - test.x;
		dx2 *= dx2;
		dy2 = support.y - test.y;
		dy2 *= dy2;
		dt2 = support.t - test.t;
		dt2 *= dt2;
		r2 = dx2*x_sig2inv+dy2*y_sig2inv+dt2*t_sig2inv;

		//return 0.001;
		//return exp(-r2*sigma2inv);
		
		if (r2 < r2_cut)
		{
			close_point++;
			//return exp(-r2*sigma2inv);
			return exp_hist->GetBinContent(int(r2*sigma2inv*bin_multiple)+1);
			//return 0.1;
		}
		else
		{
			return 0.;
		}
		
	};
	double get_log_likelihood(std::vector<dirc_point> &support_points,std::vector<dirc_point> &inpoints);
	double get_log_likelihood_speedup(std::vector<dirc_point> &support_points,std::vector<dirc_point> &inpoints);

	// diagnostic
	inline double support_spread_function_diagnostic(dirc_point support, dirc_point test, int &close_point, double &prob_x, double &prob_y, double &prob_t)__attribute__((always_inline))
	{
		double dx2,dy2,dt2,r2_x,r2_y,r2_t,r2;
		dx2 = support.x - test.x;
		dx2 *= dx2;
		dy2 = support.y - test.y;
		dy2 *= dy2;
		dt2 = support.t - test.t;
		dt2 *= dt2;
		r2_x = dx2*x_sig2inv;
		r2_y = dy2*y_sig2inv;
		r2_t = dt2*t_sig2inv;
		//r2 = dx2*x_sig2inv+dy2*y_sig2inv+dt2*t_sig2inv;
		r2 = r2_x + r2_y + r2_t;
		if (r2 < sigma_cut2 * sigma2)
		{
			close_point++;
			prob_x = exp(-r2_x*sigma2inv);
			prob_y = exp(-r2_y*sigma2inv);
			prob_t = exp(-r2_t*sigma2inv);
/*
			printf("\n\n");
			printf("log(prob_x) = %12.04f\n",log(prob_x));
			printf("log(prob_y) = %12.04f\n",log(prob_y));
			printf("log(prob_t) = %12.04f\n",log(prob_t));
			printf("log(prob)   = %12.04f\n",log(exp(-r2*sigma2inv)));
*/

			return exp(-r2*sigma2inv);
		}
		else
		{
			return 0.;
		}
	};
	void calc_log_likelihoods_diagnostic(\
				std::vector<dirc_point> &inpoints,\
				std::vector<dirc_point> &support_points_pion,\
				std::vector<dirc_point> &support_points_kaon,\
				double &ll_pion, double &ll_kaon,\
				std::vector<double> &ll_points_pion,\
				std::vector<double> &ll_points_kaon,\
				std::vector<double> &ll_diff_points);
	double get_log_likelihood_cut(std::vector<dirc_point> &support_points,std::vector<dirc_point> &inpoints);

        void calc_PDF(\
                                std::vector<dirc_point> &support_points,\
                                TH2* hist_signalPDF);
        inline double calc_PDFval_xy(double &support_x, double &support_y, double &loc_x, double &loc_y)__attribute__((always_inline))
        {

		double dx2,dy2,r2;

                dx2 = support_x - loc_x;
                dx2 *= dx2;
                dy2 = support_y - loc_y;
                dy2 *= dy2;
                //dt2 = support.t - test.t;
                //dt2 *= dt2;
                r2 = dx2*x_sig2inv+dy2*y_sig2inv;

		return norm_factor_xy * exp(-r2*sigma2inv);
        };

        void calc_hit_distances(\
				std::vector<dirc_point> &inpoints,\
				TH2* hist_supports,\
				TH2* hist_signalPDF,\
				TH1* hist_hit_distance_xy,\
				TH1* hist_hit_distance_x,\
				TH1* hist_hit_distance_y,\
				bool verbose);

};
#endif
