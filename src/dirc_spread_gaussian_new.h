#include "dirc_point.h"
#include "dirc_spread_radius.h"
#include <vector>
//#include <math.h>
#include <TRandom3.h>

#include <TFile.h>
#include <TH1.h>

#ifndef DIRC_SPREAD_GAUSSIAN_NEW
#define DIRC_SPREAD_GAUSSIAN_NEW
class DircSpreadGaussianNew
{
private:
	double x_sig2inv,y_sig2inv,t_sig2inv;
	double spread_func_norm, spread_func_norm_inv;
	double lin_slope, r_trans, sigma2, sigma2inv,max_val;
	double min_probability;
	
	double sigma_cut2, prob_cut, prob_cut_per_point;
	double exp_norm_factor;
	double r2_cut;
	
public:
	DircSpreadGaussianNew(\
		double isigma, \
		double x_unc,\
		double y_unc,\
		double t_unc,\
		double isigma_cut = 5.,\
		double imin_prob  = 1e-6);
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




	inline double support_spread_function_distances(dirc_point support, dirc_point test,
				double &loc_min_dist_x, double &loc_min_dist_y, double &loc_min_dist_t,\
				 double &loc_min_dist_r, double &loc_min_dist_xy,\
				double &loc_close_supports_x,\
				double &loc_close_supports_y,\
				double &loc_close_supports_t,\
				double &loc_close_supports_r,\
				double &loc_close_supports_xy,\
				TH1F* hDelta_x,\
				TH1F* hDelta_y,\
				TH1F* hDelta_t,\
				TH1F* hDelta_r,\
				TH1F* hDelta_xy,\
				TH1F* hDelta_t_kinbin,\
				int    &Nsupports_t_slice)__attribute__((always_inline))
	{
		double dx,dy,dt,dr,dxy;
		double dx2,dy2,dt2,r2;

		//dx = std::abs(support.x - test.x);
		dx = test.x - support.x;
		dx2 = dx*dx;
		//dy = std::abs(support.y - test.y);
		dy = test.y - support.y;
		dy2 = dy*dy;
		//dt = std::abs(support.t - test.t);
		dt = test.t - support.t;
		dt2 = dt*dt;
		r2 = dx2*x_sig2inv+dy2*y_sig2inv+dt2*t_sig2inv;
		dr = std::sqrt(r2);	

		dxy = std::sqrt(dx2 + dy2);

		// perform analysis in a reasonably wide timing window
		if (std::abs(dt) < 20.)  
		{

			Nsupports_t_slice++;

			if ( dt2 *t_sig2inv < sigma_cut2 )
				loc_close_supports_t++;
			if ( dx2 * x_sig2inv < sigma_cut2 )
				loc_close_supports_x++;
			if ( dy2 * y_sig2inv < sigma_cut2 )
				loc_close_supports_y++;
			if ( dxy * dxy * x_sig2inv * y_sig2inv < 2.* sigma_cut2)
				loc_close_supports_xy++;
	

			//if (std::abs(dt) < loc_min_dist_t)
			if (std::abs(dt) < 20 && dxy < 20.)
			{
				hDelta_x->Fill(dx);
				hDelta_y->Fill(dy);
				hDelta_xy->Fill(dxy);
				hDelta_r->Fill(dr);
		
				hDelta_t->Fill(dt);

				hDelta_t_kinbin->Fill(dt);
			}

			loc_min_dist_x  = dx < loc_min_dist_x ? dx : loc_min_dist_x;
			loc_min_dist_y  = dy < loc_min_dist_y ? dy : loc_min_dist_y;
			loc_min_dist_t  = dt < loc_min_dist_t ? dt : loc_min_dist_t;

			if (dxy < loc_min_dist_xy)
			{
				loc_min_dist_xy = dxy;				
			}
			
			if (dr < loc_min_dist_r)
			{
				loc_min_dist_r = dr;
			}

		}
		
		if ( dt < 5. && dxy < 8.5)
			loc_close_supports_r++;

		if (r2 < sigma_cut2 * sigma2)
		{
			return exp(-r2*sigma2inv);
		}
		else
		{
			return 0.;
		}
	};

	double get_log_likelihood(std::vector<dirc_point> &support_points,std::vector<dirc_point> &inpoints);
	double get_log_likelihood_speedup(std::vector<dirc_point> &support_points,std::vector<dirc_point> &inpoints);

	double get_log_likelihood_distances(std::vector<dirc_point> &support_points,std::vector<dirc_point> &inpoints, TH1F* hDelta_x,\
			TH1F* hDelta_y, TH1F* hDelta_t, TH1F* hDelta_r, TH1F* hDelta_xy, TH1F* hDelta_t_kinbin_direct, TH1F* hDelta_t_kinbin_reflected);

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


};
#endif
