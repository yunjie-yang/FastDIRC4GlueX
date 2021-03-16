#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <utility>

#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "src/dirc_point.h"
#include "src/dirc_threesegbox_sim.h"
#include "src/dirc_rect_digitizer.h"
#include "src/dirc_spread_gaussian.h"
#include "src/GlueXUserOptions.h"
#include "src/dirc_spread_gaussian_new.h"

#include <TFile.h>
#include <TTree.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TMinuit.h>


int main(int nargs, char* argv[])
{
        const char* config_str;
        char* geometry_infilename = new char[256];
        sprintf(geometry_infilename, "FastDIRC_geometry_input.csv");

        char* geometry_outfilename = new char[256];
        sprintf(geometry_outfilename,"dirc_model_geometry.csv");

        char* root_outfilename = new char[256];
        sprintf(root_outfilename,"output_hists.root");

        double kmass = .4937;
        double pimass = .1396;


        double particle_momentum = 4.;
        double particle_x = 0;
        double particle_y = 500;
        //double particle_theta = 4;
        //double particle_phi = 40;
        double particle_theta = 0;
        double particle_phi = 0;

	int num_runs = 1000;

        double ckov_unc = .003*57.3; //transport = 3mrad
        double pdf_unc_red_fac = 1;


        double mirror_r_difference = 400;

        int n_phi_phots = 150000;
        int n_z_phots = 4;


	int rseed = 1337;

        double pion_beta, kaon_beta/*, electron_beta:=1*/;
        pion_beta=kaon_beta=-1;


        double s_func_x = 6;
        double s_func_y = s_func_x;
        double s_func_t = 1.0;
        double sfunc_sig = 1;


	int n_sim_phots = 40;

        double resx = 6;
        double resy = 6;
        double t_unc = .27;
        double t_bin_size = 1;

        double digit_miny = -50;
        double digit_maxy = 300;

        double minx = -1500;
        double maxx = 1500;

        double ll_pion, ll_kaon;
        ll_pion = ll_kaon = 0.;

	//--------------- END of initialization ------------------//






	//-------------- READ IN CONFIG ------------------------//
        if (nargs==2) config_str = argv[1];
        else config_str = "config.in";
        printf("Running with config file: %s\n", config_str);

        std::map<int, int> opt_int;
        std::map<int, std::string> opt_str;
        std::map<int, double> opt_val;

        GlueXUserOptions user_opts;
        if (user_opts.ReadControl_in(config_str) == 0)
        {
                std::cerr << "Reading config.in failed" << std::endl;
                exit(-1);
        }

        if (user_opts.Find("OUTFILE", opt_str)) sprintf(root_outfilename,"%s",opt_str[1].c_str());
        if (user_opts.Find("GEOMETRY_INFILE", opt_str)) sprintf(geometry_infilename,"%s",opt_str[1].c_str());
        if (user_opts.Find("GEOMETRY_OUTFILE", opt_str)) sprintf(geometry_outfilename,"%s",opt_str[1].c_str());

        if (user_opts.Find("n", opt_val))  num_runs = int(opt_val[1]);
        if (user_opts.Find("n_phi_phots", opt_val))       n_phi_phots = int(opt_val[1]);

        if (user_opts.Find("particle_x", opt_val))            particle_x          = opt_val[1];
        if (user_opts.Find("particle_y", opt_val))            particle_y          = opt_val[1];
        if (user_opts.Find("particle_theta", opt_val))        particle_theta      = opt_val[1];
        if (user_opts.Find("particle_phi", opt_val))          particle_phi        = opt_val[1];
        if (user_opts.Find("particle_momentum", opt_val))     particle_momentum   = opt_val[1];


	
	//-------------- DEFINE DIGITIZER -------------//
        DircRectDigitizer digitizer(\
                        minx,\
                        maxx,\
                        resx,\
                        digit_miny,\
                        digit_maxy,\
                        resy,\
                        t_unc,\
                        t_bin_size);

	//-------------- DEFINE DIRC_MODEL -------------//
        double foc_mirror_size = 288;
        double main_mirror_angle_nominal = 74.0197;
        double mirror_r_nominal  = -1200 + mirror_r_difference;
        double pmt_angle_nominal = 47.87 ;

        DircThreeSegBoxSim *dirc_model = new DircThreeSegBoxSim(\
                        rseed,\
                        mirror_r_nominal,\
                        foc_mirror_size,\
                        main_mirror_angle_nominal,\
                        600,\
                        pmt_angle_nominal);


	//-------- DEFINE HISTOGRAMS ---------//
        TFile* tfile = new TFile(root_outfilename,"RECREATE");

        TH1F *ll_diff_pion = new TH1F("ll_diff_pion","Difference of log likelihood real = pion",200000,-200,200);
        TH1F *ll_diff_kaon = new TH1F("ll_diff_kaon","Difference of log likelihood real = kaon",200000,-200,200);

        TH2F *pion_dist_rowcol = new TH2F("pion_dist_rowcol","hit pattern - pion; Pixel Row ; Pixel Column",144,-0.5,143.5,48,-0.5,47.5);
	TH2F *kaon_dist_rowcol = new TH2F("kaon_dist_rowcol","hit pattern - kaon; Pixel Row ; Pixel Column",144,-0.5,143.5,48,-0.5,47.5);




	// ------- PREPARE FOR RUNNING  ----------//
	std::vector<dirc_point> support_points_pion;
	std::vector<dirc_point> support_points_kaon;

	std::vector<dirc_point> sim_points;

	DircSpreadGaussianNew* pdf_pion = new DircSpreadGaussianNew(\
			sfunc_sig,\
			s_func_x,\
			s_func_y,\
			s_func_t);
	DircSpreadGaussianNew* pdf_kaon = new DircSpreadGaussianNew(\
			sfunc_sig,\
			s_func_x,\
			s_func_y,\
			s_func_t);


	pion_beta = dirc_model->get_beta_from_p(particle_momentum,pimass);
	kaon_beta = dirc_model->get_beta_from_p(particle_momentum,kmass);

	printf("%12.04f%12.04f%12.04f%12.04f%12.04f\n",particle_momentum,particle_x,particle_y,particle_theta,particle_phi);


	for (int i = 0; i < num_runs; i++)
	{
                printf("\r                                                                      ");
                printf("\r Running particle %8d/%d  ",i+1,num_runs);
                fflush(stdout);


		//------ GENERATE PDF SUPPORTS ------//
		dirc_model->sim_reg_n_photons(\
				support_points_pion,\
				n_phi_phots,\
				n_z_phots,\
				-1,\
				1,\
				particle_x,\
				particle_y,\
				0.,\
				particle_theta,\
				particle_phi,\
				0,\
				ckov_unc/pdf_unc_red_fac,\
				pion_beta);

		dirc_model->sim_reg_n_photons(\
				support_points_kaon,\
				n_phi_phots,\
				n_z_phots,\
				-1,\
				1,\
				particle_x,\
				particle_y,\
				0.,\
				particle_theta,\
				particle_phi,\
				0,\
				ckov_unc/pdf_unc_red_fac,\
				kaon_beta);


		//------- SIMULATE POINTS ------//


		//pion hypothesis
		dirc_model->sim_rand_n_photons(\
				sim_points,\
				n_sim_phots,\
				-1,\
				1,\
				particle_x,\
				particle_y,\
				0.,\
				particle_theta,\
				particle_phi,\
				0.,\
				ckov_unc,\
				pion_beta);
		digitizer.digitize_points(sim_points);

                ll_pion = pdf_pion -> get_log_likelihood(support_points_pion,sim_points);
                ll_kaon = pdf_kaon -> get_log_likelihood(support_points_kaon,sim_points);
                ll_diff_pion->Fill(ll_pion-ll_kaon);

                for (unsigned int i = 0; i < sim_points.size(); i++)
                {
                        pion_dist_rowcol->Fill(sim_points[i].pixel_row,sim_points[i].pixel_col);
                }


                //kaon hypothesis
                dirc_model->sim_rand_n_photons(\
                                sim_points,\
                                n_sim_phots,\
                                -1,\
                                1,\
                                particle_x,\
                                particle_y,\
                                0.,\
                                particle_theta,\
                                particle_phi,\
                                0.,\
                                ckov_unc,\
                                kaon_beta);
                digitizer.digitize_points(sim_points);

                ll_pion = pdf_pion -> get_log_likelihood(support_points_pion,sim_points);
                ll_kaon = pdf_kaon -> get_log_likelihood(support_points_kaon,sim_points);
                ll_diff_kaon->Fill(ll_pion-ll_kaon);

                for (unsigned int i = 0; i < sim_points.size(); i++)
                {
                        kaon_dist_rowcol->Fill(sim_points[i].pixel_row,sim_points[i].pixel_col);
                }


	}//END OF num_runs loop

        printf("\n Run Completed. \n");


	//--------- WRITE HISTOGRAMS -----------//
        tfile->cd();

        ll_diff_pion->Write();
        ll_diff_kaon->Write();
	
	pion_dist_rowcol->Write();
	kaon_dist_rowcol->Write();

	tfile->Close();
	return 1;




}
