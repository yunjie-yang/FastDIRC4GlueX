#include "dirc_point.h"
#include "dirc_spread_gaussian_exp.h"
#include <vector>
#include <math.h>

DircSpreadGaussianExp::DircSpreadGaussianExp(\
	double isigma, \
	double x_unc,\
	double y_unc,\
	double t_unc,\
	double isigma_cut,\
	double imin_prob,\
	const char* iexp_hist_filename)
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
	//min_probability = 1e-4;
	//min_probability = 1e-2;
	//min_probability = 1e-3;

	sprintf(exp_hist_filename,"/home/yunjiey/Documents/bayes_calib/FastDIRC_dev/cached_exp.root");
	exp_hist_file = new TFile(exp_hist_filename);
	exp_hist = (TH1D*) exp_hist_file -> Get("hist_exp");
	bin_multiple = exp_hist->GetNbinsX()/int(exp_hist->GetBinLowEdge(exp_hist->GetNbinsX()+1));

	r2_cut = sigma_cut2 * sigma2;

	norm_factor_xy = 1/(M_PI*sigma2*x_unc*y_unc);

}
/*
double DircSpreadGaussianExp::get_log_likelihood(std::vector<dirc_point> &support_points, std::vector<dirc_point> &inpoints)
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

double DircSpreadGaussianExp::get_log_likelihood(std::vector<dirc_point> &support_points, std::vector<dirc_point> &inpoints)
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
	if (Npoints_passed>0)
		rval -= log(Npoints_passed);
	else
		rval -= log(inpoints.size());
	
	return rval;
}

double DircSpreadGaussianExp::get_log_likelihood_speedup(std::vector<dirc_point> &support_points, std::vector<dirc_point> &inpoints)
{
	double tprob = 0;
	double rval = 0;
	double log_mult = 1;
	double weight = 1;
	int close_point = 0;
	int Npoints_passed = 0;
	//double rval_prod = 1.;
	//double rval_prod_log = 0.;
	for (unsigned int i = 0; i < inpoints.size(); i++)
	{
		tprob = 0;
		close_point = 0;
		for (unsigned int j = 0; j < support_points.size(); j++)
		{
			//tprob += support_spread_function(support_points[j],inpoints[i],close_point);
			tprob += support_spread_function_speedup(support_points[j],inpoints[i],close_point);
			//tprob += 0.001; 
		}
		tprob /= support_points.size();
		tprob *= spread_func_norm_inv;
		
		rval += weight*log_mult*log(tprob+min_probability);
		//rval_prod *= tprob+min_probability;

		if (close_point > 0)
		{
			inpoints[i].cut = false;
			Npoints_passed++;
		}
		else	
			inpoints[i].cut = true;

		//if (rval_prod < 1e-300)
		//{
		//	rval_prod_log += weight*log_mult*log(rval_prod);
		//	rval_prod = 1.;
		//}
	}
	/*
	if (Npoints_passed>0)
		rval -= log(Npoints_passed);
	else
		rval -= log(inpoints.size());
	*/
	//rval_prod_log += weight*log_mult*log(rval_prod);
	//printf("\n\n");
	//printf("rval          = %f\n",rval);	
	//printf("rval_prod_log = %f\n",rval_prod_log);	

	return rval;
	//return rval_prod_log;
}

void DircSpreadGaussianExp::set_gaus_sigma(double isigma)
{
	sigma2 = isigma*isigma;
	sigma2inv = 1/sigma2;
	spread_func_norm = 1;
	spread_func_norm_inv=1/spread_func_norm;
}


void DircSpreadGaussianExp::calc_log_likelihoods_diagnostic(\
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

double DircSpreadGaussianExp::get_log_likelihood_cut(\
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


void DircSpreadGaussianExp::calc_PDF(\
                                std::vector<dirc_point> &support_points,\
                                TH2* hist_signalPDF)
{
	int bin_x(-1),bin_y(-1);
	double loc_x(-999.), loc_y(-999.);
	double loc_i_x(-999.),loc_i_y(-999.);
	double loc_i_content(0.);
	for (unsigned int loc_i = 0 ; loc_i < support_points.size(); loc_i++ )
	{
		//if (support_points[loc_i].updown!=1)
		//	continue;

		loc_x = support_points[loc_i].x;
		loc_y = support_points[loc_i].y;
		bin_x = hist_signalPDF -> GetXaxis() -> FindBin(loc_x);
		bin_y = hist_signalPDF -> GetYaxis() -> FindBin(loc_y);
	
		for (int bin_i_x = bin_x - 5; bin_i_x <= bin_x + 5 ; bin_i_x++ )
		for (int bin_i_y = bin_y - 5; bin_i_y <= bin_y + 5 ; bin_i_y++ )
		{
			loc_i_x = hist_signalPDF->GetXaxis()->GetBinCenter(bin_i_x);
			loc_i_y = hist_signalPDF->GetYaxis()->GetBinCenter(bin_i_y);

			loc_i_content = hist_signalPDF->GetBinContent(bin_i_x,bin_i_y) + calc_PDFval_xy(loc_x,loc_y,loc_i_x,loc_i_y);
			hist_signalPDF->SetBinContent(bin_i_x,bin_i_y,loc_i_content);
		}

	}

	hist_signalPDF -> Scale(1./support_points.size());

}

void DircSpreadGaussianExp::calc_hit_distances(\
                                std::vector<dirc_point> &inpoints,\
                                TH2* hist_supports,\
                                TH2* hist_signalPDF,\
                                TH1* hist_hit_distance_xy,\
                                TH1* hist_hit_distance_x,\
                                TH1* hist_hit_distance_y,\
				bool verbose)
{

	double threshold = 5e-6;
	//double threshold = 1e-6;
	int layer_limit = 50;

	double loc_x(-999.),loc_y(-999.);
	int loc_x_bin(-1),loc_y_bin(-1);

	bool flag = true;
	int layer_counter = 0;

	int left_x_bin(0),right_x_bin(0),upper_y_bin(0),lower_y_bin(0);
	double loc_content = -999.;

	double loc_x_tmp(-1.),loc_y_tmp(-1);
	int loc_x_bin_tmp(-1),loc_y_bin_tmp(-1);
	double dist2_xy_tmp(-1.);

	int found_layer = -1;

	int loc_pixelrow = -1;
	int loc_pixelcol = -1;

	double min_dist2_xy = 1e10;
	double delta_x     = 999.;
	double delta_y     = 999.;

	int min_bin_x = 999;
	int min_bin_y = 999;

	for (unsigned int loc_i = 0 ; loc_i < inpoints.size() ; loc_i++)
	{
		//initialize

		flag = true;
		layer_counter = 0;
		found_layer = -1;
		min_dist2_xy = 1e10;
		delta_x     = 999.;
		delta_y     = 999.;
		min_bin_x = 999;
		min_bin_y = 999;



		loc_x = inpoints[loc_i].x;
		loc_y = inpoints[loc_i].y;
		loc_x_bin = hist_signalPDF->GetXaxis()->FindBin(loc_x); 
		loc_y_bin = hist_signalPDF->GetYaxis()->FindBin(loc_y); 
	
		if (verbose)
		{
			printf("\n\nhit #%d\n",loc_i);
			printf("(x,y)     = (%8.02f,%8.02f)\n",loc_x,loc_y);
			printf("(x,y) bin = (%8d,%8d)\n",loc_x_bin,loc_y_bin);
		}
		
		loc_pixelrow = hist_supports->GetXaxis()->FindBin(inpoints[loc_i].pixel_row);
		loc_pixelcol = hist_supports->GetYaxis()->FindBin(inpoints[loc_i].pixel_col);
	
		hist_supports  -> SetBinContent(loc_pixelrow,loc_pixelcol,100); 


		//searching from inside-out
		while (flag)
		{

			left_x_bin  = loc_x_bin - layer_counter;
			right_x_bin = loc_x_bin + layer_counter;

			upper_y_bin = loc_y_bin + layer_counter;
			lower_y_bin = loc_y_bin - layer_counter;

			if (verbose)
			{
				printf("\nlayer_counter = %d\n", layer_counter);
				printf("found_layer = %d, flag = %d\n",found_layer,flag);
			}
			//upper side
			if (verbose)
				printf("upper side...\n");
			loc_y_bin_tmp = upper_y_bin;
			for (loc_x_bin_tmp = left_x_bin ; loc_x_bin_tmp <= right_x_bin; loc_x_bin_tmp++)
			{
				loc_content = hist_signalPDF->GetBinContent(loc_x_bin_tmp,loc_y_bin_tmp);
				if (verbose)
					printf("loop: bin (x,y) = (%6d,%6d); content = %.2e\n",loc_x_bin_tmp,loc_y_bin_tmp,loc_content);
				if (loc_content > threshold)
				{
					found_layer = layer_counter;
					flag = false;
					loc_x_tmp = hist_signalPDF->GetXaxis()->GetBinCenter(loc_x_bin_tmp);
					loc_y_tmp = hist_signalPDF->GetYaxis()->GetBinCenter(loc_y_bin_tmp);
					dist2_xy_tmp = (loc_x_tmp - loc_x)*(loc_x_tmp - loc_x) + (loc_y_tmp - loc_y)*(loc_y_tmp-loc_y);
					if (dist2_xy_tmp < min_dist2_xy)
					{
						min_dist2_xy = dist2_xy_tmp;
						delta_x = loc_x_tmp - loc_x;
						delta_y = loc_y_tmp - loc_y;
						min_bin_x = loc_x_bin_tmp;
						min_bin_y = loc_y_bin_tmp;

					}
				}
			}
			//lower side
			if (verbose)
				printf("lower side...\n");
			loc_y_bin_tmp = lower_y_bin;
			for (int loc_x_bin_tmp = left_x_bin ; loc_x_bin_tmp <= right_x_bin; loc_x_bin_tmp++)
			{
				loc_content = hist_signalPDF->GetBinContent(loc_x_bin_tmp,loc_y_bin_tmp);
				if (verbose)
					printf("loop: bin (x,y) = (%6d,%6d); content = %.2e\n",loc_x_bin_tmp,loc_y_bin_tmp,loc_content);
				if (loc_content > threshold)
				{
					found_layer = layer_counter;
					flag = false;	
					loc_x_tmp = hist_signalPDF->GetXaxis()->GetBinCenter(loc_x_bin_tmp);
					loc_y_tmp = hist_signalPDF->GetYaxis()->GetBinCenter(loc_y_bin_tmp);
					dist2_xy_tmp = (loc_x_tmp - loc_x)*(loc_x_tmp - loc_x) + (loc_y_tmp - loc_y)*(loc_y_tmp-loc_y);
					if (dist2_xy_tmp < min_dist2_xy)
					{
						min_dist2_xy = dist2_xy_tmp;
						delta_x = loc_x_tmp - loc_x;
						delta_y = loc_y_tmp - loc_y;
						min_bin_x = loc_x_bin_tmp;
						min_bin_y = loc_y_bin_tmp;
					}
				}
			}
			//left side
			if (verbose)
				printf("left side...\n");
			loc_x_bin_tmp = left_x_bin;
			for (loc_y_bin_tmp = lower_y_bin + 1 ; loc_y_bin_tmp <= upper_y_bin -1 ; loc_y_bin_tmp++)
			{
				loc_content = hist_signalPDF->GetBinContent(loc_x_bin_tmp,loc_y_bin_tmp);
				if (verbose)
					printf("loop: bin (x,y) = (%6d,%6d); content = %.2e\n",loc_x_bin_tmp,loc_y_bin_tmp,loc_content);
				if (loc_content > threshold)
				{
					found_layer = layer_counter;
					flag = false;	
					loc_x_tmp = hist_signalPDF->GetXaxis()->GetBinCenter(loc_x_bin_tmp);
					loc_y_tmp = hist_signalPDF->GetYaxis()->GetBinCenter(loc_y_bin_tmp);
					dist2_xy_tmp = (loc_x_tmp - loc_x)*(loc_x_tmp - loc_x) + (loc_y_tmp - loc_y)*(loc_y_tmp-loc_y);
					if (dist2_xy_tmp < min_dist2_xy)
					{
						min_dist2_xy = dist2_xy_tmp;
						delta_x = loc_x_tmp - loc_x;
						delta_y = loc_y_tmp - loc_y;
						min_bin_x = loc_x_bin_tmp;
						min_bin_y = loc_y_bin_tmp;
					}
				}
			}
			//right side
			if (verbose)
				printf("right side...\n");
			loc_x_bin_tmp = right_x_bin;
			for (loc_y_bin_tmp = lower_y_bin + 1 ; loc_y_bin_tmp <= upper_y_bin -1 ; loc_y_bin_tmp++)
			{
				loc_content = hist_signalPDF->GetBinContent(loc_x_bin_tmp,loc_y_bin_tmp);
				if (verbose)
					printf("loop: bin (x,y) = (%6d,%6d); content = %.2e\n",loc_x_bin_tmp,loc_y_bin_tmp,loc_content);
				if (loc_content > threshold)
				{
					found_layer = layer_counter;
					flag = false;	
					loc_x_tmp = hist_signalPDF->GetXaxis()->GetBinCenter(loc_x_bin_tmp);
					loc_y_tmp = hist_signalPDF->GetYaxis()->GetBinCenter(loc_y_bin_tmp);
					dist2_xy_tmp = (loc_x_tmp - loc_x)*(loc_x_tmp - loc_x) + (loc_y_tmp - loc_y)*(loc_y_tmp-loc_y);
					if (dist2_xy_tmp < min_dist2_xy)
					{
						min_dist2_xy = dist2_xy_tmp;
						delta_x = loc_x_tmp - loc_x;
						delta_y = loc_y_tmp - loc_y;
						min_bin_x = loc_x_bin_tmp;
						min_bin_y = loc_y_bin_tmp;
					}
				}
			}

			if (layer_counter > layer_limit)
			{
				flag = false;
			}
			layer_counter++;

		}//end of search loop

		if (found_layer < 0)
			continue;
		
		if (verbose)
		//if (true)
		{
			printf("\n\nfound_layer = %d\n",found_layer);
			printf("(x,y)     = (%8.02f,%8.02f)\n",loc_x,loc_y);
			printf("(x,y) bin = (%8d,%8d)\n",loc_x_bin,loc_y_bin);
			printf("sqrt(min_dist2_xy) = %12.04f\n",sqrt(min_dist2_xy));
			printf("delta_x            = %12.04f\n",delta_x);
			printf("delta_y            = %12.04f\n",delta_y);
			printf("min_loc_x          = %12.04f\n",loc_x_tmp);
			printf("min_loc_y          = %12.04f\n",loc_y_tmp);
			printf("min_bin_x          = %3d\n",min_bin_x);
			printf("min_bin_y          = %3d\n",min_bin_y);

		}

		hist_hit_distance_xy->Fill(sqrt(min_dist2_xy));
		hist_hit_distance_x->Fill(delta_x);
		hist_hit_distance_y->Fill(delta_y);



	}//end hit loop




}

