#include "dirc_point.h"
#include "dirc_spread_gaussian_new.h"
#include <vector>
#include <math.h>

DircSpreadGaussianNew::DircSpreadGaussianNew(\
	double isigma, \
	double x_unc,\
	double y_unc,\
	double t_unc,\
	double isigma_cut,\
	double imin_prob)
{

	sigma2 = isigma*isigma;
	sigma2inv = 1/sigma2;
	spread_func_norm = 1;
	spread_func_norm_inv=1/spread_func_norm;

	exp_norm_factor = 1/(pow(M_PI,3/2)*isigma*isigma*isigma*x_unc*y_unc*t_unc);
	
	
	x_sig2inv = 1/(x_unc*x_unc);
	y_sig2inv = 1/(y_unc*y_unc);
	t_sig2inv = 1/(t_unc*t_unc);
	
	sigma_cut2 = isigma_cut*isigma_cut;

	prob_cut_per_point = exp(-sigma_cut2*sigma2);

	min_probability = imin_prob;

	r2_cut = sigma_cut2 * sigma2;

}
/*
double DircSpreadGaussianNew::get_log_likelihood(std::vector<dirc_point> &support_points, std::vector<dirc_point> &inpoints)
{
	double tprob = 0;
	double rval = 0;
	int eval_count = 0;
	double log_mult = 1;
	double weight = 1;
	for (unsigned int i = 0; i < inpoints.size(); i++)
	{
		tprob = 0;
		for (unsigned int j = 0; j < support_points.size(); j++)
		{
			tprob += support_spread_function(support_points[j],inpoints[i]);
			eval_count++;
		}
		tprob /= support_points.size();
		tprob *= spread_func_norm_inv;
		
		rval += weight*log_mult*log(tprob+min_probability);
	}
	rval -= log(inpoints.size());
	
	return rval;
}
*/
double DircSpreadGaussianNew::get_log_likelihood(std::vector<dirc_point> &support_points, std::vector<dirc_point> &inpoints)
{
        double tprob = 0;
        double rval = 0;
        double log_mult = 1;
        double weight = 1;
        int close_point = 0; 
        int Npoints_passed = 0;
        for (unsigned int i = 0; i < inpoints.size(); i++)
        {
                tprob = 0;
                close_point = 0;
                for (unsigned int j = 0; j < support_points.size(); j++)
                {
                        tprob += support_spread_function(support_points[j],inpoints[i],close_point);
                }
                tprob /= support_points.size();
                tprob *= spread_func_norm_inv;

                rval += weight*log_mult*log(tprob+min_probability);

                if (close_point > 0)
                {
                        inpoints[i].cut = false;
                        Npoints_passed++;
                }
                else
                        inpoints[i].cut = true;
        }
        return rval;
}



double DircSpreadGaussianNew::get_log_likelihood_distances(std::vector<dirc_point> &support_points, std::vector<dirc_point> &inpoints,TH1F* hDelta_x, TH1F* hDelta_y, TH1F* hDelta_t, TH1F* hDelta_r, TH1F* hDelta_xy, TH1F* hDelta_t_kinbin_direct, TH1F* hDelta_t_kinbin_reflected)
{
	double tprob = 0;
	double rval = 0;
	double log_mult = 1;
	double weight = 1;
	int Npoints_passed = 0;
	//inpoints are data points

        double loc_min_dist_x  = 10000.;
        double loc_min_dist_y  = 10000.;
        double loc_min_dist_t  = 10000.;
        double loc_min_dist_r  = 10000.;
        double loc_min_dist_xy = 10000.;
     	double loc_close_supports_x = 0;
        double loc_close_supports_y = 0;
        double loc_close_supports_t = 0;
        double loc_close_supports_r = 0;
        double loc_close_supports_xy = 0;

	int Nsupports_t_slice = 0;

	for (unsigned int i = 0; i < inpoints.size(); i++)
	{
		tprob = 0;

		inpoints[i].min_dist[0] = 10000.;
		inpoints[i].min_dist[1] = 10000.;
		inpoints[i].min_dist[2] = 10000.;
		inpoints[i].min_dist[3] = 10000.;
		inpoints[i].min_dist[4] = 10000.;
		inpoints[i].frac_close[0] = -1;
		inpoints[i].frac_close[1] = -1;
		inpoints[i].frac_close[2] = -1;
		inpoints[i].frac_close[3] = -1;
		inpoints[i].frac_close[4] = -1;

	
		Nsupports_t_slice = 0;

        	loc_min_dist_x  = 10000.;
        	loc_min_dist_y  = 10000.;
        	loc_min_dist_t  = 10000.;
        	loc_min_dist_r  = 10000.;
        	loc_min_dist_xy = 10000.;

        	loc_close_supports_x  = 0;
        	loc_close_supports_y  = 0;
        	loc_close_supports_t  = 0;
        	loc_close_supports_r  = 0;
        	loc_close_supports_xy = 0;


		if (inpoints[i].direct == true)
		{
			for (unsigned int j = 0; j < support_points.size(); j++)
			{
				tprob += support_spread_function_distances(support_points[j],\
								inpoints[i],
								loc_min_dist_x,loc_min_dist_y,\
								loc_min_dist_t,loc_min_dist_r,loc_min_dist_xy,\
								loc_close_supports_x,loc_close_supports_y,\
								loc_close_supports_t,loc_close_supports_r,loc_close_supports_xy,\
								hDelta_x,hDelta_y,hDelta_t,hDelta_r, hDelta_xy, hDelta_t_kinbin_direct,\
								Nsupports_t_slice);
			}
		}
		else
		{
			for (unsigned int j = 0; j < support_points.size(); j++)
			{
				tprob += support_spread_function_distances(support_points[j],\
								inpoints[i],
								loc_min_dist_x,loc_min_dist_y,\
								loc_min_dist_t,loc_min_dist_r,loc_min_dist_xy,\
								loc_close_supports_x,loc_close_supports_y,\
								loc_close_supports_t,loc_close_supports_r,loc_close_supports_xy,\
								hDelta_x,hDelta_y,hDelta_t,hDelta_r, hDelta_xy, hDelta_t_kinbin_reflected,\
								Nsupports_t_slice);
			}
		}
		tprob /= support_points.size();
		tprob *= spread_func_norm_inv;
		
		rval += weight*log_mult*log(tprob+min_probability);

		if (loc_close_supports_r > 0)
		{
			inpoints[i].cut = false;
			Npoints_passed++;
		}
		else	
			inpoints[i].cut = true;

		inpoints[i].min_dist[0] = loc_min_dist_x;
		inpoints[i].min_dist[1] = loc_min_dist_y;
		inpoints[i].min_dist[2] = loc_min_dist_t;
		inpoints[i].min_dist[3] = loc_min_dist_r;
		inpoints[i].min_dist[4] = loc_min_dist_xy;

		inpoints[i].frac_close[0] = loc_close_supports_x  / Nsupports_t_slice;
		inpoints[i].frac_close[1] = loc_close_supports_y  / Nsupports_t_slice;
		inpoints[i].frac_close[2] = loc_close_supports_t  / Nsupports_t_slice;
		inpoints[i].frac_close[3] = loc_close_supports_r  / Nsupports_t_slice;
		inpoints[i].frac_close[4] = loc_close_supports_xy / Nsupports_t_slice;

	}
	return rval;
}

void DircSpreadGaussianNew::set_gaus_sigma(double isigma)
{
	sigma2 = isigma*isigma;
	sigma2inv = 1/sigma2;
	spread_func_norm = 1;
	spread_func_norm_inv=1/spread_func_norm;
}


void DircSpreadGaussianNew::calc_log_likelihoods_diagnostic(\
				std::vector<dirc_point> &inpoints,\
				std::vector<dirc_point> &support_points_pion,\
				std::vector<dirc_point> &support_points_kaon,\
				double &ll_pion, double &ll_kaon,\
				std::vector<double> &ll_points_pion,\
				std::vector<double> &ll_points_kaon,\
				std::vector<double> &ll_diff_points)
{
	double tprob_pion = 0;
	double tprob_kaon = 0;

	double rval_pion = 0;
	double rval_kaon = 0;

	int close_point_pion = 0;
	int close_point_kaon = 0;

	int Npoints_passed = 0;

	//double prob_mult_factor  = exp_norm_factor ;
	//double mean_support_points  = (double(support_points_pion.size())+double(support_points_kaon.size()))/2.;
	//double prob_mult_factor_pion = exp_norm_factor * mean_support_points/ double(support_points_pion.size());
	//double prob_mult_factor_kaon = exp_norm_factor * mean_support_points/ double(support_points_kaon.size());

	//double prob_mult_factor_pion = exp_norm_factor ;
	//double prob_mult_factor_kaon = exp_norm_factor ; 

	//double prob_cut_pion = prob_cut_per_point * int(support_points_pion.size());   
	//double prob_cut_kaon = prob_cut_per_point * int(support_points_kaon.size());   
	//double prob_cut_pion = exp(-8.);   
	//double prob_cut_kaon = exp(-8.);   

	double tmp_logprob_pion;
	double tmp_logprob_kaon;

	ll_points_pion.clear();
	ll_points_kaon.clear();
	ll_diff_points.clear();

	for (unsigned int i = 0; i < inpoints.size(); i++)
	{
		tprob_pion       = 0;
		tprob_kaon       = 0;
		close_point_pion = 0;
		close_point_kaon = 0;

		tmp_logprob_pion = 999;
		tmp_logprob_kaon = 999;

		for (unsigned int j = 0; j < support_points_pion.size(); j++)
		{
			tprob_pion += support_spread_function(support_points_pion[j],inpoints[i],close_point_pion);
		}
		for (unsigned int j = 0; j < support_points_kaon.size(); j++)
		{
			tprob_kaon += support_spread_function(support_points_kaon[j],inpoints[i],close_point_kaon);
		}



		// min prob idea
		tprob_pion /= double(support_points_pion.size());
		tprob_kaon /= double(support_points_kaon.size());

		tmp_logprob_pion = log(tprob_pion + min_probability);
		rval_pion += tmp_logprob_pion;

		tmp_logprob_kaon = log(tprob_kaon + min_probability);
		rval_kaon += tmp_logprob_kaon;

		ll_points_pion.push_back(tmp_logprob_pion);
		ll_points_kaon.push_back(tmp_logprob_kaon);
		ll_diff_points.push_back(tmp_logprob_pion - tmp_logprob_kaon);

/*
		// cut idea
		tprob_pion *= prob_mult_factor_pion;
		tprob_kaon *= prob_mult_factor_kaon;

		//if (tprob_pion > prob_cut_pion || tprob_kaon > prob_cut_kaon)
		if (tprob_pion > prob_cut_pion && tprob_kaon > prob_cut_kaon)
		{
			tmp_logprob_pion = log(tprob_pion);
			//tmp_logprob_pion = log(tprob_pion + prob_cut_pion);
			//tmp_logprob_pion = log(tprob_pion + min_probability);
			rval_pion  += tmp_logprob_pion;
			ll_points_pion.push_back(tmp_logprob_pion);


			tmp_logprob_kaon = log(tprob_kaon);
			//tmp_logprob_kaon = log(tprob_kaon + prob_cut_kaon);
			//tmp_logprob_kaon = log(tprob_kaon + min_probability);
			rval_kaon  += tmp_logprob_kaon;
			ll_points_kaon.push_back(tmp_logprob_kaon);

			ll_diff_points.push_back(tmp_logprob_pion - tmp_logprob_kaon);
		}
		//printf(" point #%2d, log(prob)_pion = %8.04f, log(prob)_kaon = %8.04f \n",i+1,tmp_logprob_pion,tmp_logprob_kaon);
*/



		if (close_point_pion > 0 && close_point_kaon > 0)
		{
			inpoints[i].cut = false;
			Npoints_passed++;
		}
		else	
			inpoints[i].cut = true;

	}



	ll_pion = rval_pion;
	ll_kaon = rval_kaon;

	//printf("ll_pion = %8.04f, ll_kaon = %8.04f \n",ll_pion,ll_kaon);
}

double DircSpreadGaussianNew::get_log_likelihood_cut(\
				std::vector<dirc_point> &support_points,\
				std::vector<dirc_point> &inpoints)
{
	double tprob = 0;

	double rval = 0;

	int close_point = 0;
	int Npoints_passed = 0;

	double prob_mult_factor  = exp_norm_factor ;
	double prob_cut = exp(-8.);   

	for (unsigned int i = 0; i < inpoints.size(); i++)
	{
		tprob       = 0;
		close_point = 0;

		for (unsigned int j = 0; j < support_points.size(); j++)
		{
			tprob += support_spread_function(support_points[j],inpoints[i],close_point);
		}

		tprob *= prob_mult_factor;

		if (tprob > prob_cut)
		{
			rval  += log(tprob);
		}

		if (close_point > 0 )
		{
			inpoints[i].cut = false;
			Npoints_passed++;
		}
		else	
			inpoints[i].cut = true;

	}

	return rval;
}
