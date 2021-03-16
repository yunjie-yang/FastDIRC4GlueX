#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <utility>

#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "src/dirc_point.h"
#include "src/dirc_photon.h"
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


#include <TChain.h>
#include <TClonesArray.h>
#include "src/DrcHit.h"
#include "src/DrcEvent.h"

#define rad2deg 57.2958


int main(int nargs, char* argv[])
{
        const char* config_str;

        char* dirctree_filename = new char[256];
        sprintf(dirctree_filename,"dirc_tree.root");

        char* geometry_infilename = new char[256];
        sprintf(geometry_infilename, "FastDIRC_geometry_input.csv");

        char* geometry_outfilename = new char[256];
        sprintf(geometry_outfilename,"dirc_model_geometry.csv");

        char* root_outfilename = new char[256];
        sprintf(root_outfilename,"output_hists.root");

        double kmass = .493677;
        double pimass = .139571;


        double particle_momentum = 4.;
        double particle_x = 0;
        double particle_y = 0;
        double particle_theta = 4;
        double particle_phi = 40;
	double particle_bar = 1;

        double ckov_unc = .003*57.3; //transport = 3mrad
        double pdf_unc_red_fac = 1;


        double mirror_r_difference = 400;

        int n_phi_phots = 150000;
        int n_z_phots = 4;

	int n_sim_phots = 40;

	int rseed = 1337;

        double pion_beta, kaon_beta;
        pion_beta=kaon_beta=-1;


        double s_func_x = 6;
        double s_func_y = s_func_x;
        double s_func_t = 1;
        double sfunc_sig = 1;
        double sigma_cut = 5;


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

        int Ntracks = 1000;
        int Ntracks_pion = 1000;
        int Ntracks_kaon = 1000;


	int scale_support = 10;


	TRandom3* rand_gen = new TRandom3(1337);
	//bool smear_time = true;
	bool smear_time = false;

	double sigma_t_pion_smear = 0.;
	double sigma_t_kaon_smear = 0.;

	bool perform_selection = 1;

	//------------- OFFSETS --------------------//

        // PMT plane 
        double pmt_y_off_SouthLower     = 0;
        double pmt_z_off_SouthLower     = 0;
        double pmt_angle_off_SouthLower = 0;

        double pmt_y_off_NorthUpper     = 0;
        double pmt_z_off_NorthUpper     = 0;
        double pmt_angle_off_NorthUpper = 0;
        // 3-seg mirror
        std::vector<double> threeSeg_x_angle_offs_SouthLower = {0.,0.,0.};
        std::vector<double> threeSeg_y_angle_offs_SouthLower = {0.,0.,0.};
        std::vector<double> threeSeg_z_angle_offs_SouthLower = {0.,0.,0.};

        std::vector<double> threeSeg_y_pos_offs_SouthLower = {0.,0.,0.};
        std::vector<double> threeSeg_z_pos_offs_SouthLower = {0.,0.,0.};

        std::vector<double> threeSeg_x_angle_offs_NorthUpper = {0.,0.,0.};
        std::vector<double> threeSeg_y_angle_offs_NorthUpper = {0.,0.,0.};
        std::vector<double> threeSeg_z_angle_offs_NorthUpper = {0.,0.,0.};

        std::vector<double> threeSeg_y_pos_offs_NorthUpper = {0.,0.,0.};
        std::vector<double> threeSeg_z_pos_offs_NorthUpper = {0.,0.,0.};

        // bar box w.r.t. optical box
        double BB_OB_x_angle_off_SouthLower = 0;
        double BB_OB_y_angle_off_SouthLower = 0;
        double BB_OB_z_angle_off_SouthLower = 0;
        double BB_OB_x_pos_off_SouthLower = 0;
        double BB_OB_y_pos_off_SouthLower = 0;
        double BB_OB_z_pos_off_SouthLower = 0;

        double BB_OB_x_angle_off_NorthUpper = 0;
        double BB_OB_y_angle_off_NorthUpper = 0;
        double BB_OB_z_angle_off_NorthUpper = 0;
        double BB_OB_x_pos_off_NorthUpper = 0;
        double BB_OB_y_pos_off_NorthUpper = 0;
        double BB_OB_z_pos_off_NorthUpper = 0;

	// sidemirrors
	double sidemirror_xr_x_off_SouthLower = 0;
	double sidemirror_xl_x_off_SouthLower = 0;

	double sidemirror_xr_x_off_NorthUpper = 0;
	double sidemirror_xl_x_off_NorthUpper = 0;

	// large planar mirror
	double largePlanarMirrorY_y_off_SouthLower = 0;
	double largePlanarMirrorY_y_off_NorthUpper = 0;

	// large planar mirror
	double upperWedgeFarZ_z_off_SouthLower = 0;
	double upperWedgeFarZ_z_off_NorthUpper = 0;

	// support points offsets
	double support_x_off_SouthLower = 0.;
	double support_y_off_SouthLower = 0.;
	double support_t_off_SouthLower = 0.;

	double support_x_off_NorthUpper = 0.;
	double support_y_off_NorthUpper = 0.;
	double support_t_off_NorthUpper = 0.;


	// bar box w.r.t. tracking system

	double BB10_offsets_x_pos   = 0.;
	double BB10_offsets_y_pos   = 0.;
	double BB10_offsets_z_pos   = 0.;
	double BB10_offsets_x_angle = 0.;
	double BB10_offsets_y_angle = 0.;
	double BB10_offsets_z_angle = 0.;

	double BB11_offsets_x_pos   = 0.;
	double BB11_offsets_y_pos   = 0.;
	double BB11_offsets_z_pos   = 0.;
	double BB11_offsets_x_angle = 0.;
	double BB11_offsets_y_angle = 0.;
	double BB11_offsets_z_angle = 0.;

	double BB00_offsets_x_pos   = 0.;
	double BB00_offsets_y_pos   = 0.;
	double BB00_offsets_z_pos   = 0.;
	double BB00_offsets_x_angle = 0.;
	double BB00_offsets_y_angle = 0.;
	double BB00_offsets_z_angle = 0.;

	double BB01_offsets_x_pos   = 0.;
	double BB01_offsets_y_pos   = 0.;
	double BB01_offsets_z_pos   = 0.;
	double BB01_offsets_x_angle = 0.;
	double BB01_offsets_y_angle = 0.;
	double BB01_offsets_z_angle = 0.;


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
        if (user_opts.Find("DIRCTREE_INFILE", opt_str)) sprintf(dirctree_filename,"%s",opt_str[1].c_str());

        if (user_opts.Find("n_phi_phots", opt_val))       n_phi_phots = int(opt_val[1]);

        if (user_opts.Find("Ntracks", opt_val))            Ntracks = int(opt_val[1]);
        if (user_opts.Find("Ntracks_pion", opt_val))       Ntracks_pion = int(opt_val[1]);
        if (user_opts.Find("Ntracks_kaon", opt_val))       Ntracks_kaon = int(opt_val[1]);




	std::string tmp_str;

	// ------------- OFFSETS FOR SOUTHLOWER BOX --------------------//
        for (int seg_i = 0 ; seg_i < 3; seg_i++)
        {
                tmp_str = "threeSeg"+std::to_string(seg_i+1)+"_x_angle_off_SouthLower";
                if (user_opts.Find(tmp_str.c_str(), opt_val))  threeSeg_x_angle_offs_SouthLower[seg_i] = opt_val[1];
                tmp_str = "threeSeg"+std::to_string(seg_i+1)+"_y_angle_off_SouthLower";
                if (user_opts.Find(tmp_str.c_str(), opt_val))  threeSeg_y_angle_offs_SouthLower[seg_i] = opt_val[1];
                tmp_str = "threeSeg"+std::to_string(seg_i+1)+"_z_angle_off_SouthLower";
                if (user_opts.Find(tmp_str.c_str(), opt_val))  threeSeg_z_angle_offs_SouthLower[seg_i] = opt_val[1];
                tmp_str = "threeSeg"+std::to_string(seg_i+1)+"_y_pos_off_SouthLower";
                if (user_opts.Find(tmp_str.c_str(), opt_val))  threeSeg_y_pos_offs_SouthLower[seg_i] = opt_val[1];
                tmp_str = "threeSeg"+std::to_string(seg_i+1)+"_z_pos_off_SouthLower";
                if (user_opts.Find(tmp_str.c_str(), opt_val))  threeSeg_z_pos_offs_SouthLower[seg_i] = opt_val[1];
	}

        if (user_opts.Find("pmt_y_off_SouthLower", opt_val))          pmt_y_off_SouthLower = opt_val[1];
        if (user_opts.Find("pmt_z_off_SouthLower", opt_val))          pmt_z_off_SouthLower = opt_val[1];
        if (user_opts.Find("pmt_angle_off_SouthLower", opt_val))      pmt_angle_off_SouthLower = opt_val[1];

        if (user_opts.Find("BB_OB_x_pos_off_SouthLower", opt_val))       BB_OB_x_pos_off_SouthLower = opt_val[1];
        if (user_opts.Find("BB_OB_y_pos_off_SouthLower", opt_val))       BB_OB_y_pos_off_SouthLower = opt_val[1];
        if (user_opts.Find("BB_OB_z_pos_off_SouthLower", opt_val))       BB_OB_z_pos_off_SouthLower = opt_val[1];

        if (user_opts.Find("BB_OB_x_angle_off_SouthLower", opt_val))       BB_OB_x_angle_off_SouthLower = opt_val[1];
        if (user_opts.Find("BB_OB_y_angle_off_SouthLower", opt_val))       BB_OB_y_angle_off_SouthLower = opt_val[1];
        if (user_opts.Find("BB_OB_z_angle_off_SouthLower", opt_val))       BB_OB_z_angle_off_SouthLower = opt_val[1];

        if (user_opts.Find("sidemirror_xr_x_off_SouthLower", opt_val))     sidemirror_xr_x_off_SouthLower = opt_val[1];
        if (user_opts.Find("sidemirror_xl_x_off_SouthLower", opt_val))     sidemirror_xl_x_off_SouthLower = opt_val[1];

        if (user_opts.Find("largePlanarMirrorY_y_off_SouthLower", opt_val))    largePlanarMirrorY_y_off_SouthLower = opt_val[1];

        if (user_opts.Find("upperWedgeFarZ_z_off_SouthLower", opt_val))   upperWedgeFarZ_z_off_SouthLower = opt_val[1];

        if (user_opts.Find("support_x_off_SouthLower", opt_val))       support_x_off_SouthLower = opt_val[1];
        if (user_opts.Find("support_y_off_SouthLower", opt_val))       support_y_off_SouthLower = opt_val[1];
        if (user_opts.Find("support_t_off_SouthLower", opt_val))       support_t_off_SouthLower = opt_val[1];

        if (user_opts.Find("BB10_offsets_x_pos", opt_val))         BB10_offsets_x_pos   = opt_val[1];
        if (user_opts.Find("BB10_offsets_y_pos", opt_val))         BB10_offsets_y_pos   = opt_val[1];
        if (user_opts.Find("BB10_offsets_z_pos", opt_val))         BB10_offsets_z_pos   = opt_val[1];
        if (user_opts.Find("BB10_offsets_x_angle", opt_val))       BB10_offsets_x_angle = opt_val[1];
        if (user_opts.Find("BB10_offsets_y_angle", opt_val))       BB10_offsets_y_angle = opt_val[1];
        if (user_opts.Find("BB10_offsets_z_angle", opt_val))       BB10_offsets_z_angle = opt_val[1];

        if (user_opts.Find("BB11_offsets_x_pos", opt_val))         BB11_offsets_x_pos   = opt_val[1];
        if (user_opts.Find("BB11_offsets_y_pos", opt_val))         BB11_offsets_y_pos   = opt_val[1];
        if (user_opts.Find("BB11_offsets_z_pos", opt_val))         BB11_offsets_z_pos   = opt_val[1];
        if (user_opts.Find("BB11_offsets_x_angle", opt_val))       BB11_offsets_x_angle = opt_val[1];
        if (user_opts.Find("BB11_offsets_y_angle", opt_val))       BB11_offsets_y_angle = opt_val[1];
        if (user_opts.Find("BB11_offsets_z_angle", opt_val))       BB11_offsets_z_angle = opt_val[1];

	// ------------- OFFSETS FOR NORTHUPPER BOX --------------------//
	
        for (int seg_i = 0 ; seg_i < 3; seg_i++)
        {
                tmp_str = "threeSeg"+std::to_string(seg_i+1)+"_x_angle_off_NorthUpper";
                if (user_opts.Find(tmp_str.c_str(), opt_val))  threeSeg_x_angle_offs_NorthUpper[seg_i] = opt_val[1];
                tmp_str = "threeSeg"+std::to_string(seg_i+1)+"_y_angle_off_NorthUpper";
                if (user_opts.Find(tmp_str.c_str(), opt_val))  threeSeg_y_angle_offs_NorthUpper[seg_i] = opt_val[1];
                tmp_str = "threeSeg"+std::to_string(seg_i+1)+"_z_angle_off_NorthUpper";
                if (user_opts.Find(tmp_str.c_str(), opt_val))  threeSeg_z_angle_offs_NorthUpper[seg_i] = opt_val[1];
                tmp_str = "threeSeg"+std::to_string(seg_i+1)+"_y_pos_off_NorthUpper";
                if (user_opts.Find(tmp_str.c_str(), opt_val))  threeSeg_y_pos_offs_NorthUpper[seg_i] = opt_val[1];
                tmp_str = "threeSeg"+std::to_string(seg_i+1)+"_z_pos_off_NorthUpper";
                if (user_opts.Find(tmp_str.c_str(), opt_val))  threeSeg_z_pos_offs_NorthUpper[seg_i] = opt_val[1];
	}

        if (user_opts.Find("pmt_y_off_NorthUpper", opt_val))          pmt_y_off_NorthUpper = opt_val[1];
        if (user_opts.Find("pmt_z_off_NorthUpper", opt_val))          pmt_z_off_NorthUpper = opt_val[1];
        if (user_opts.Find("pmt_angle_off_NorthUpper", opt_val))      pmt_angle_off_NorthUpper = opt_val[1];

        if (user_opts.Find("BB_OB_x_pos_off_NorthUpper", opt_val))       BB_OB_x_pos_off_NorthUpper = opt_val[1];
        if (user_opts.Find("BB_OB_y_pos_off_NorthUpper", opt_val))       BB_OB_y_pos_off_NorthUpper = opt_val[1];
        if (user_opts.Find("BB_OB_z_pos_off_NorthUpper", opt_val))       BB_OB_z_pos_off_NorthUpper = opt_val[1];

        if (user_opts.Find("BB_OB_x_angle_off_NorthUpper", opt_val))       BB_OB_x_angle_off_NorthUpper = opt_val[1];
        if (user_opts.Find("BB_OB_y_angle_off_NorthUpper", opt_val))       BB_OB_y_angle_off_NorthUpper = opt_val[1];
        if (user_opts.Find("BB_OB_z_angle_off_NorthUpper", opt_val))       BB_OB_z_angle_off_NorthUpper = opt_val[1];

        if (user_opts.Find("sidemirror_xr_x_off_NorthUpper", opt_val))     sidemirror_xr_x_off_NorthUpper = opt_val[1];
        if (user_opts.Find("sidemirror_xl_x_off_NorthUpper", opt_val))     sidemirror_xl_x_off_NorthUpper = opt_val[1];

        if (user_opts.Find("largePlanarMirrorY_y_off_NorthUpper", opt_val))    largePlanarMirrorY_y_off_NorthUpper = opt_val[1];

        if (user_opts.Find("upperWedgeFarZ_z_off_NorthUpper", opt_val))   upperWedgeFarZ_z_off_NorthUpper = opt_val[1];

        if (user_opts.Find("support_x_off_NorthUpper", opt_val))       support_x_off_NorthUpper = opt_val[1];
        if (user_opts.Find("support_y_off_NorthUpper", opt_val))       support_y_off_NorthUpper = opt_val[1];
        if (user_opts.Find("support_t_off_NorthUpper", opt_val))       support_t_off_NorthUpper = opt_val[1];

        if (user_opts.Find("BB00_offsets_x_pos", opt_val))         BB00_offsets_x_pos   = opt_val[1];
        if (user_opts.Find("BB00_offsets_y_pos", opt_val))         BB00_offsets_y_pos   = opt_val[1];
        if (user_opts.Find("BB00_offsets_z_pos", opt_val))         BB00_offsets_z_pos   = opt_val[1];
        if (user_opts.Find("BB00_offsets_x_angle", opt_val))       BB00_offsets_x_angle = opt_val[1];
        if (user_opts.Find("BB00_offsets_y_angle", opt_val))       BB00_offsets_y_angle = opt_val[1];
        if (user_opts.Find("BB00_offsets_z_angle", opt_val))       BB00_offsets_z_angle = opt_val[1];

        if (user_opts.Find("BB01_offsets_x_pos", opt_val))         BB01_offsets_x_pos   = opt_val[1];
        if (user_opts.Find("BB01_offsets_y_pos", opt_val))         BB01_offsets_y_pos   = opt_val[1];
        if (user_opts.Find("BB01_offsets_z_pos", opt_val))         BB01_offsets_z_pos   = opt_val[1];
        if (user_opts.Find("BB01_offsets_x_angle", opt_val))       BB01_offsets_x_angle = opt_val[1];
        if (user_opts.Find("BB01_offsets_y_angle", opt_val))       BB01_offsets_y_angle = opt_val[1];
        if (user_opts.Find("BB01_offsets_z_angle", opt_val))       BB01_offsets_z_angle = opt_val[1];




        if (user_opts.Find("sigma_cut", opt_val))   sigma_cut = opt_val[1];

        if (user_opts.Find("sfunc_sig", opt_val))   sfunc_sig = opt_val[1];

        if (user_opts.Find("s_func_x", opt_val))   s_func_x = opt_val[1];
        if (user_opts.Find("s_func_y", opt_val))   s_func_y = opt_val[1];
        if (user_opts.Find("s_func_t", opt_val))   s_func_t = opt_val[1];

        if (user_opts.Find("sigma_t_pion_smear", opt_val))   sigma_t_pion_smear = opt_val[1];
        if (user_opts.Find("sigma_t_kaon_smear", opt_val))   sigma_t_kaon_smear = opt_val[1];

        if (user_opts.Find("PERFORM_SELECTION", opt_int))       perform_selection = bool(opt_int[1]);

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
                        pmt_angle_nominal,\
			geometry_infilename);


	//-------- DEFINE HISTOGRAMS ---------//
        TFile* tfile = new TFile(root_outfilename,"RECREATE");

	TDirectory* dir_integrated = tfile->mkdir("Integrated");
	dir_integrated->cd();

	// event/track selection diagnostics
	TH1F *hInvMass_rho = new TH1F("hInvMass_rho","Inv. Mass #rho",1500,0.,1.5);
	TH1F *hInvMass_phi = new TH1F("hInvMass_phi","Inv. Mass #phi",1500,0.5,2.0);

	TH1F *hMissingMassSquared_rho = new TH1F("hMissingMassSquared_rho","Missing Mass Squared",500,-0.1,0.1);
	TH1F *hMissingMassSquared_phi = new TH1F("hMissingMassSquared_phi","Missing Mass Squared",500,-0.1,0.1);

	TH1F *hChiSq_rho = new TH1F("hChiSq_rho","; #rho KinFit #chi^{2};",180,0.,45.);
	TH1F *hChiSq_phi = new TH1F("hChiSq_phi","; #phi KinFit #chi^{2};",180,0.,45.);

	TH1F *hTofTrackDist_pion = new TH1F("hTofTrackDist_pion","#pi; dist(extra. DIRC hit, matching TOF hit;",100,0.,10.);
	TH1F *hTofTrackDist_kaon = new TH1F("hTofTrackDist_kaon","K; dist(extra. DIRC hit, matching TOF hit;",100,0.,10.);

	TH1F *hTofTrackDeltaT_pion = new TH1F("hTofTrackDeltaT_pion","#pi; TOF #DeltaT (ns);",160,-2.,2.);	
	TH1F *hTofTrackDeltaT_kaon = new TH1F("hTofTrackDeltaT_kaon","K; TOF #DeltaT (ns);",160,-2.,2.);	

	TH1F *hDcHits_pion = new TH1F("hDcHits_pion","#pi; Num. DC Hits;",100,0,100);
	TH1F *hDcHits_kaon = new TH1F("hDcHits_kaon","K; Num. DC Hits;",100,0,100);

	TH2F *hTrackOccupancy_pion = new TH2F("hTrackOccupancy_pion","extrapolated point; x; bar",40,-100.,100.,48,-0.5,47.5);
	TH2F *hTrackOccupancy_kaon = new TH2F("hTrackOccupancy_kaon","extrapolated point; x; bar",40,-100.,100.,48,-0.5,47.5);


	// DIRC 
        TH1F *hit_dist_t_tree = new TH1F("hit_dist_t_tree","hit time ; hit time (ns);",800,0,800);
        TH1F *hit_dist_t_cut_tree = new TH1F("hit_dist_t_cut_tree","hit time ; hit time (ns);",800,0,800);


	// likelihood
        TH1F *ll_diff_pion_tree = new TH1F("ll_diff_pion_tree","Difference of log likelihood, pion selection",200000,-200,200);
        TH1F *ll_diff_kaon_tree = new TH1F("ll_diff_kaon_tree","Difference of log likelihood, kaon selection",200000,-200,200);
        TH1F *ll_diff_pion_tree_SouthLower = new TH1F("ll_diff_pion_tree_SouthLower",\
						"Difference of log likelihood, pion selection",200000,-200,200);
        TH1F *ll_diff_kaon_tree_SouthLower = new TH1F("ll_diff_kaon_tree_SouthLower",\
						"Difference of log likelihood, kaon selection",200000,-200,200);
        TH1F *ll_diff_pion_tree_NorthUpper = new TH1F("ll_diff_pion_tree_NorthUpper",\
						"Difference of log likelihood, pion selection",200000,-200,200);
        TH1F *ll_diff_kaon_tree_NorthUpper = new TH1F("ll_diff_kaon_tree_NorthUpper",\
						"Difference of log likelihood, kaon selection",200000,-200,200);


        TH1F *ll_diff_pion_sim = new TH1F("ll_diff_pion_sim","Difference of log likelihood, pion selection",200000,-200,200);
        TH1F *ll_diff_kaon_sim = new TH1F("ll_diff_kaon_sim","Difference of log likelihood, kaon selection",200000,-200,200);
        TH1F *ll_diff_pion_sim_SouthLower = new TH1F("ll_diff_pion_sim_SouthLower",\
						"Difference of log likelihood, pion selection",200000,-200,200);
        TH1F *ll_diff_kaon_sim_SouthLower = new TH1F("ll_diff_kaon_sim_SouthLower",\
						"Difference of log likelihood, kaon selection",200000,-200,200);
        TH1F *ll_diff_pion_sim_NorthUpper = new TH1F("ll_diff_pion_sim_NorthUpper",\
						"Difference of log likelihood, pion selection",200000,-200,200);
        TH1F *ll_diff_kaon_sim_NorthUpper = new TH1F("ll_diff_kaon_sim_NorthUpper",\
						"Difference of log likelihood, kaon selection",200000,-200,200);

	// photon yield 
        TH1F *Nph_tree          = new TH1F("Nph_tree","number of hits", 100,-.5,99.5);
        TH1F *Nph_cut_tree      = new TH1F("Nph_cut_tree","number of hits (selection)", 100,-.5,99.5);
        TH1F *Nph_cut_pion_tree = new TH1F("Nph_cut_pion_tree","number of photons (pion selection)", 100,-.5,99.5);
        TH1F *Nph_cut_kaon_tree = new TH1F("Nph_cut_kaon_tree","number of photons (kaon selection)", 100,-.5,99.5);

        TH1F *Nph_tree_SouthLower          = new TH1F("Nph_tree_SouthLower","number of hits", 100,-.5,99.5);
        TH1F *Nph_cut_tree_SouthLower      = new TH1F("Nph_cut_tree_SouthLower","number of hits (selection)", 100,-.5,99.5);
        TH1F *Nph_cut_pion_tree_SouthLower = new TH1F("Nph_cut_pion_tree_SouthLower","number of photons (pion selection)", 100,-.5,99.5);
        TH1F *Nph_cut_kaon_tree_SouthLower = new TH1F("Nph_cut_kaon_tree_SouthLower","number of photons (kaon selection)", 100,-.5,99.5);

        TH1F *Nph_tree_NorthUpper          = new TH1F("Nph_tree_NorthUpper","number of hits", 100,-.5,99.5);
        TH1F *Nph_cut_tree_NorthUpper      = new TH1F("Nph_cut_tree_NorthUpper","number of hits (selection)", 100,-.5,99.5);
        TH1F *Nph_cut_pion_tree_NorthUpper = new TH1F("Nph_cut_pion_tree_NorthUpper","number of photons (pion selection)", 100,-.5,99.5);
        TH1F *Nph_cut_kaon_tree_NorthUpper = new TH1F("Nph_cut_kaon_tree_NorthUpper","number of photons (kaon selection)", 100,-.5,99.5);


	// divide into bar boxes
	std::vector<int> BBNums = {10,11,1,0};

	TH1F *Nph_tree_BB[4];
	TH1F *Nph_cut_tree_BB[4];
	TH1F *Nph_cut_pion_tree_BB[4];
	TH1F *Nph_cut_kaon_tree_BB[4];

        TH1F *ll_diff_pion_tree_BB[4];
        TH1F *ll_diff_kaon_tree_BB[4];
        TH1F *ll_diff_pion_sim_BB[4];
        TH1F *ll_diff_kaon_sim_BB[4];


	for (int loc_i = 0; loc_i < 4; loc_i++ )
	{
		Nph_tree_BB[loc_i]          = new TH1F(Form("Nph_tree_BB%02d",BBNums[loc_i]),"number of hits", 200,-.5,199.5);
		Nph_cut_tree_BB[loc_i]      = new TH1F(Form("Nph_cut_tree_BB%02d",BBNums[loc_i]),"number of hits", 100,-.5,99.5);
		Nph_cut_pion_tree_BB[loc_i] = new TH1F(Form("Nph_cut_pion_tree_BB%02d",BBNums[loc_i]),"number of hits", 100,-.5,99.5);
		Nph_cut_kaon_tree_BB[loc_i] = new TH1F(Form("Nph_cut_kaon_tree_BB%02d",BBNums[loc_i]),"number of hits", 100,-.5,99.5);
        
		ll_diff_pion_tree_BB[loc_i] = new TH1F(Form("ll_diff_pion_tree_%02d",BBNums[loc_i]),\
                                                "Difference of log likelihood, pion selection",200000,-200,200);
		ll_diff_kaon_tree_BB[loc_i] = new TH1F(Form("ll_diff_kaon_tree_%02d",BBNums[loc_i]),\
                                                "Difference of log likelihood, kaon selection",200000,-200,200);
		ll_diff_pion_sim_BB[loc_i] = new TH1F(Form("ll_diff_pion_sim_%02d",BBNums[loc_i]),\
                                                "Difference of log likelihood, pion selection",200000,-200,200);
		ll_diff_kaon_sim_BB[loc_i] = new TH1F(Form("ll_diff_kaon_sim_%02d",BBNums[loc_i]),\
                                                "Difference of log likelihood, kaon selection",200000,-200,200);
	}


	TH2F *hDLLvsNph_pion_tree = new TH2F("hDLLvsNph_pion_tree",";Nph;DLL(#pi-K)",100,0,100,1000,-100,100);
	TH2F *hDLLvsNph_kaon_tree = new TH2F("hDLLvsNph_kaon_tree",";Nph;DLL(#pi-K)",100,0,100,1000,-100,100);

	TH2F *hDLLvsNphNear_pion_tree = new TH2F("hDLLvsNphNear_pion_tree",";Nph;DLL(#pi-K)",100,0,100,1000,-100,100);
	TH2F *hDLLvsNphNear_kaon_tree = new TH2F("hDLLvsNphNear_kaon_tree",";Nph;DLL(#pi-K)",100,0,100,1000,-100,100);

	TH2F *hDLLvsNph_pion_sim = new TH2F("hDLLvsNph_pion_sim",";Nph;DLL(#pi-K)",100,0,100,1000,-100,100);
	TH2F *hDLLvsNph_kaon_sim = new TH2F("hDLLvsNph_kaon_sim",";Nph;DLL(#pi-K)",100,0,100,1000,-100,100);

	TH1F *hTimePerTrack = new TH1F("hTimePerTrack",";Time;",4000,0,200.);


	// Min dist measure hists

	//std::vector<std::string> box_labels = {"SouthLower","NorthUpper"};
	std::vector<std::string> BB_labels = {"BB10","BB11","BB01","BB00"};
	std::vector<std::string> PID_labels = {"pion","kaon"};// PDF hypothesis
	std::vector<std::string> hit_types  = {"tree","sim"};
	std::vector<std::string> obs_labels = {"X","Y","T","R","XY"};
	std::vector<std::string> photon_types  = {"direct","reflected"};
	std::vector<std::vector<double>> obs_minmax = {{-20.,20.},{-20.,20.},{-20.,20.},{0.,20.},{0.,20.}};

	TH1F *hMap_MinDist[5][2][2][4];//obs, pion/kaon, tree/sim, BB
	TH1F *hMap_FracClose[5][2][2][4];
	TH1F *hMap_Delta[5][2][2][4];

	for (int loc_BB = 0 ; loc_BB < int(BB_labels.size()); loc_BB++)
	for (int loc_hittype = 0; loc_hittype < int(hit_types.size()); loc_hittype++)
	for (int loc_PID = 0 ; loc_PID < int(PID_labels.size()); loc_PID++)
	for (int loc_dim = 0 ; loc_dim < int(obs_labels.size()); loc_dim++)
	{
		hMap_MinDist[loc_dim][loc_PID][loc_hittype][loc_BB] = new TH1F(\
		Form("hMinDist%s_%s_%s_%s",obs_labels[loc_dim].c_str(),PID_labels[loc_PID].c_str(),hit_types[loc_hittype].c_str(),BB_labels[loc_BB].c_str()),Form("; min(#Delta %s)",obs_labels[loc_dim].c_str()),1000,0.,100.);

		hMap_FracClose[loc_dim][loc_PID][loc_hittype][loc_BB] = new TH1F(\
		Form("hFracClose%s_%s_%s_%s",obs_labels[loc_dim].c_str(),PID_labels[loc_PID].c_str(),hit_types[loc_hittype].c_str(),BB_labels[loc_BB].c_str()),Form("; frac. of support points that are \"close\" to a hit in %s",obs_labels[loc_dim].c_str()),1000,0.,1.);

		hMap_Delta[loc_dim][loc_PID][loc_hittype][loc_BB] = new TH1F(\
		Form("hDelta%s_%s_%s_%s",obs_labels[loc_dim].c_str(),PID_labels[loc_PID].c_str(),hit_types[loc_hittype].c_str(),BB_labels[loc_BB].c_str()),Form("; Delta %s (observed - support)",obs_labels[loc_dim].c_str()),1000,obs_minmax[loc_dim][0],obs_minmax[loc_dim][1]);
	}


        const int NumBars   = 48;
        const int NumXBins  = 40;

        const int NumPBins  = 1;
        //float PBinEdges[NumPBins][2] = {{2.9,3.1}};
        //float PBinEdges[NumPBins][2] = {{2.5,3.5}};
        //int PVals[NumPBins] = {3};

        //float PBinEdges[NumPBins][2] = {{2.9,3.1}};
        float PBinEdges[NumPBins][2] = {{4.,10.}};
        int PVals[NumPBins] = {4};


        double XBin_min, XBin_max;
        std::string locLabel_name,locLabel_title;

        TDirectory* PBin_dirs[NumPBins];
        TDirectory* KinBin_dirs[NumPBins][NumBars][NumXBins];

        TH2F* hMap_TrackOccupancy_pion[NumPBins];
        TH2F* hMap_TrackOccupancy_kaon[NumPBins];


        // (x,y) dist.
        TH2F* hMap_hit_dist_rowcol_pion_near_tree[NumPBins][NumBars][NumXBins];
        TH2F* hMap_hit_dist_rowcol_kaon_near_tree[NumPBins][NumBars][NumXBins];
        TH2F* hMap_hit_dist_rowcol_pion_away_tree[NumPBins][NumBars][NumXBins];
        TH2F* hMap_hit_dist_rowcol_kaon_away_tree[NumPBins][NumBars][NumXBins];

        TH2F* hMap_hit_dist_rowcol_pion_sim[NumPBins][NumBars][NumXBins];
        TH2F* hMap_hit_dist_rowcol_kaon_sim[NumPBins][NumBars][NumXBins];

        TH2F* hMap_support_dist_rowcol_pion_sim[NumPBins][NumBars][NumXBins];
        TH2F* hMap_support_dist_rowcol_kaon_sim[NumPBins][NumBars][NumXBins];


	// t dist.
        TH1F* hMap_hit_dist_t_pion_tree[NumPBins][NumBars][NumXBins];
        TH1F* hMap_hit_dist_t_kaon_tree[NumPBins][NumBars][NumXBins];

        TH1F* hMap_hit_dist_t_pion_sim[NumPBins][NumBars][NumXBins];
        TH1F* hMap_hit_dist_t_kaon_sim[NumPBins][NumBars][NumXBins];


        // Nph
        TH1F* hMap_Nph_pion_cut_tree[NumPBins][NumBars][NumXBins];
        TH1F* hMap_Nph_kaon_cut_tree[NumPBins][NumBars][NumXBins];

        TH1F* hMap_Nph_pion_cut_sim[NumPBins][NumBars][NumXBins];
        TH1F* hMap_Nph_kaon_cut_sim[NumPBins][NumBars][NumXBins];

	TH1F *hMap_Delta_Bar_XBin[NumPBins][NumBars][NumXBins][1][2][2][2];//DeltaT, pion/kaon, tree/sim, direct/reflected


        for (int locPBin = 0 ; locPBin < NumPBins ; locPBin++ )
        {
                PBin_dirs[locPBin] = tfile->mkdir(Form("P_%d",PVals[locPBin]));
                PBin_dirs[locPBin] -> cd();

                hMap_TrackOccupancy_pion[locPBin] = new TH2F(\
                                                Form("hTrackOccupancy_pion_%d",PVals[locPBin]),\
                                                "extrapolated point; x; bar ",\
                                                40,-100.,100.,48,-0.5,47.5);
                hMap_TrackOccupancy_kaon[locPBin] = new TH2F(\
                                                Form("hTrackOccupancy_kaon_%d",PVals[locPBin]),\
                                                "extrapolated point; x; bar ",\
                                                40,-100.,100.,48,-0.5,47.5);



                for (int locBar  = 0 ; locBar  < NumBars  ; locBar++  )
                for (int locXBin = 0 ; locXBin < NumXBins ; locXBin++ )
                {

                        locLabel_name  = Form("%d_%d_%d",PVals[locPBin],locBar,locXBin);
                        KinBin_dirs[locPBin][locBar][locXBin] = PBin_dirs[locPBin]->mkdir(locLabel_name.c_str());
                        KinBin_dirs[locPBin][locBar][locXBin] -> cd();

                        XBin_min = -100.0 + locXBin * 5.0;
                        XBin_max = XBin_min + 5.0;

                        locLabel_title = Form("Bar %d, x:(%3.01f, %3.01f), P:(%2.01f, %2.01f)",\
                                                locBar,XBin_min,XBin_max,PBinEdges[locPBin][0],PBinEdges[locPBin][1]);


			// (x,y) dist.
                        hMap_hit_dist_rowcol_pion_near_tree[locPBin][locBar][locXBin] = new TH2F(\
                                                        Form("hit_dist_rowcol_pion_near_tree_%s",locLabel_name.c_str()),\
                                                        Form("hit pattern %s; Pixel Row ; Pixel Column",locLabel_title.c_str()),\
                                                        144,-0.5,143.5,48,-0.5,47.5);
                        hMap_hit_dist_rowcol_kaon_near_tree[locPBin][locBar][locXBin] = new TH2F(\
                                                        Form("hit_dist_rowcol_kaon_near_tree_%s",locLabel_name.c_str()),\
                                                        Form("hit pattern %s; Pixel Row ; Pixel Column",locLabel_title.c_str()),\
                                                        144,-0.5,143.5,48,-0.5,47.5);
                        hMap_hit_dist_rowcol_pion_away_tree[locPBin][locBar][locXBin] = new TH2F(\
                                                        Form("hit_dist_rowcol_pion_away_tree_%s",locLabel_name.c_str()),\
                                                        Form("hit pattern %s; Pixel Row ; Pixel Column",locLabel_title.c_str()),\
                                                        144,-0.5,143.5,48,-0.5,47.5);
                        hMap_hit_dist_rowcol_kaon_away_tree[locPBin][locBar][locXBin] = new TH2F(\
                                                        Form("hit_dist_rowcol_kaon_away_tree_%s",locLabel_name.c_str()),\
                                                        Form("hit pattern %s; Pixel Row ; Pixel Column",locLabel_title.c_str()),\
                                                        144,-0.5,143.5,48,-0.5,47.5);

                        hMap_hit_dist_rowcol_pion_sim[locPBin][locBar][locXBin] = new TH2F(\
                                                        Form("hit_dist_rowcol_pion_sim_%s",locLabel_name.c_str()),\
                                                        Form("hit pattern %s; Pixel Row ; Pixel Column",locLabel_title.c_str()),\
                                                        144,-0.5,143.5,48,-0.5,47.5);
                        hMap_hit_dist_rowcol_kaon_sim[locPBin][locBar][locXBin] = new TH2F(\
                                                        Form("hit_dist_rowcol_kaon_sim_%s",locLabel_name.c_str()),\
                                                        Form("hit pattern %s; Pixel Row ; Pixel Column",locLabel_title.c_str()),\
                                                        144,-0.5,143.5,48,-0.5,47.5);

                        hMap_support_dist_rowcol_pion_sim[locPBin][locBar][locXBin] = new TH2F(\
                                                        Form("support_dist_rowcol_pion_sim_%s",locLabel_name.c_str()),\
                                                        Form("hit pattern %s; Pixel Row ; Pixel Column",locLabel_title.c_str()),\
                                                        144,-0.5,143.5,48,-0.5,47.5);
                        hMap_support_dist_rowcol_kaon_sim[locPBin][locBar][locXBin] = new TH2F(\
                                                        Form("support_dist_rowcol_kaon_sim_%s",locLabel_name.c_str()),\
                                                        Form("hit pattern %s; Pixel Row ; Pixel Column",locLabel_title.c_str()),\
                                                        144,-0.5,143.5,48,-0.5,47.5);


			// t dist.
                        hMap_hit_dist_t_pion_tree[locPBin][locBar][locXBin] = new TH1F(\
                                                                        Form("hit_dist_t_pion_tree_%s",locLabel_name.c_str()),\
                                                                        Form("hit time %s;t (ns);",locLabel_title.c_str()),\
                                                                        450,0,450);
                        hMap_hit_dist_t_kaon_tree[locPBin][locBar][locXBin] = new TH1F(\
                                                                        Form("hit_dist_t_kaon_tree_%s",locLabel_name.c_str()),\
                                                                        Form("hit time %s;t (ns);",locLabel_title.c_str()),\
                                                                        450,0,450);

                        hMap_hit_dist_t_pion_sim[locPBin][locBar][locXBin] = new TH1F(\
                                                                        Form("hit_dist_t_pion_sim_%s",locLabel_name.c_str()),\
                                                                        Form("hit time %s;t (ns);",locLabel_title.c_str()),\
                                                                        450,0,450);
                        hMap_hit_dist_t_kaon_sim[locPBin][locBar][locXBin] = new TH1F(\
                                                                        Form("hit_dist_t_kaon_sim_%s",locLabel_name.c_str()),\
                                                                        Form("hit time %s;t (ns);",locLabel_title.c_str()),\
                                                                        450,0,450);



			// Nph

                        hMap_Nph_pion_cut_tree[locPBin][locBar][locXBin] = new TH1F(\
                                                                        Form("Nph_pion_cut_tree_%s",locLabel_name.c_str()),\
                                                                        Form("%s, # photons (#pi data)",locLabel_title.c_str()),\
                                                                        100,-0.5,99.5);
                        hMap_Nph_kaon_cut_tree[locPBin][locBar][locXBin] = new TH1F(\
                                                                        Form("Nph_kaon_cut_tree_%s",locLabel_name.c_str()),\
                                                                        Form("%s, # photons (K data)",locLabel_title.c_str()),\
                                                                        100,-0.5,99.5);

                        hMap_Nph_pion_cut_sim[locPBin][locBar][locXBin] = new TH1F(\
                                                                        Form("Nph_pion_sim_cut_%s",locLabel_name.c_str()),\
                                                                        Form("%s, # photons (#pi sim)",locLabel_title.c_str()),\
                                                                        100,-0.5,99.5);
                        hMap_Nph_kaon_cut_sim[locPBin][locBar][locXBin] = new TH1F(\
                                                                        Form("Nph_kaon_sim_cut_%s",locLabel_name.c_str()),\
                                                                        Form("%s, # photons (K sim)",locLabel_title.c_str()),\
                                                                        100,-0.5,99.5);

			
			// DeltaT
			for (int loc_hittype = 0; loc_hittype < int(hit_types.size()); loc_hittype++)
			for (int loc_PID = 0 ; loc_PID < int(PID_labels.size()); loc_PID++)
			for (int loc_photontype = 0; loc_photontype < int(photon_types.size()); loc_photontype++)
			{
				hMap_Delta_Bar_XBin[locPBin][locBar][locXBin][0][loc_PID][loc_hittype][loc_photontype] = \ 
					new TH1F(Form("hDeltaT_%s_%s_%s",PID_labels[loc_PID].c_str(),hit_types[loc_hittype].c_str(),\
					photon_types[loc_photontype].c_str()),\
					Form("%s, #Delta t; Delta%s (observed - support)",locLabel_title.c_str(),obs_labels[2].c_str()),\
					500,obs_minmax[2][0],obs_minmax[2][1]);
			}


		}
	}

	//-------------- INITIALIZE TREE --------------------//

	bool locSouthLower = true;
	bool locSelectedRegion = false;


	int row_factor = 1;
	int multi      = 0;//pixelrow = 143 * multi + row_factor * data_pixel_row

	int locBin_PID = -1;
	int locBin_P   = -1;
	int locBin_X   = -1;
	int locBin_Bar = -1;
	int locBin_BB  = -1;
	int locBin_hittype = -1;

        bool locGoodParticle = false;

	
	double locInvMass = -1.;
	double locMissingMassSquared = -999.;

	double locChiSq = -1.;
	double locTofTrackDist = 999.;
	double locTofTrackDeltaT = 999.;
	int locDcHits = -1;


        TVector3 particle_p3_tree, particle_x3_tree;
        double particle_mom_tree   = 0.;
        double particle_x_tree     = 0.;
        double particle_t_tree     = 0.;
	int    particle_pid_tree   = -1;
	int    particle_bar_tree   = -1;


        int hit_ChannelId   = -1;
        int hit_pixelrow    = -1;
        int hit_pixelcol    = -1;
        double hit_t        = -1;

        std::vector<dirc_point> tree_points;

        int counter = 0;

	bool locPassMMS            = false;
	bool locPassInvMass        = false;
	bool locPassChiSq          = false;
	bool locPassTofTrackDist   = false;
	bool locPassTofTrackDeltaT = false;
	bool locPassDcHits         = false;
	

        TChain* event_chain = new TChain("dirc");
        event_chain -> Add(dirctree_filename);
        DrcEvent* evt = new DrcEvent();
        TClonesArray* particle_array = new TClonesArray("DrcEvent");
        event_chain->SetBranchAddress("DrcEvent", &particle_array);
        printf("\n\n Total Entries = %d\n",int(event_chain->GetEntries()));

	//------------- INITIALIZE FastDIRC -----------------//

	std::vector<dirc_point> support_points_pion;
	std::vector<dirc_point> support_points_kaon;

	DircSpreadGaussianNew* pdf_pion = new DircSpreadGaussianNew(\
			sfunc_sig,\
			s_func_x,\
			s_func_y,\
			s_func_t,\
			sigma_cut);

	DircSpreadGaussianNew* pdf_kaon = new DircSpreadGaussianNew(\
			sfunc_sig,\
			s_func_x,\
			s_func_y,\
			s_func_t,\
			sigma_cut);

	std::vector<dirc_point> sim_points_pion;
	std::vector<dirc_point> sim_points_kaon;

	int track_i      = 0;
	int track_pion_i = 0;
	int track_kaon_i = 0;


	int Nph_total_tree = -1;

	bool break_signal = false;

	clock_t timing_clock;
	double gen_time;

        for (int event_i = 0 ; event_i < event_chain->GetEntries() ; event_i++)
        {

		event_chain->GetEntry(event_i);
		if(event_i%100==0)
			printf("Event #%d Counter = %d\n",event_i,counter);

                if (Ntracks > 0 && track_i >= Ntracks)
                        break;

                if (Ntracks_pion > 0 && track_pion_i >= Ntracks_pion && \
                    Ntracks_kaon > 0 && track_kaon_i >= Ntracks_kaon)
                        break;

		if (break_signal)	
			break;
	
		for (int particle_i = 0; particle_i < particle_array->GetEntriesFast(); particle_i++)
		{
			locPassMMS            = false;
			locPassInvMass        = false;
			locPassChiSq          = false;
			locPassTofTrackDist   = false;
			locPassTofTrackDeltaT = false;
			locPassDcHits         = false;

			locGoodParticle = false;

			locSouthLower           = true;

			

			locSelectedRegion       = false;

	                evt      = (DrcEvent*) particle_array->At(particle_i);

               		locInvMass             = evt -> GetInvMass();
               		locMissingMassSquared  = evt -> GetMissMass();
			locChiSq               = evt -> GetChiSq();
			locTofTrackDist        = evt -> GetTofTrackDist();
			locTofTrackDeltaT      = evt -> GetTofTrackDeltaT();
			locDcHits              = evt -> GetDcHits();
	

	                particle_p3_tree       = evt->GetMomentum();
        	        particle_x3_tree       = evt->GetPosition();
                	particle_pid_tree      = evt->GetPdg();

                	particle_bar_tree      = int(evt->GetId());
                	particle_x_tree     = particle_x3_tree.X();

	                particle_mom_tree      = particle_p3_tree.Mag();


			if (abs(particle_pid_tree)==211)
				locBin_PID = 0;
			else if (abs(particle_pid_tree)==321)
				locBin_PID = 1;
			else
				continue;

			if (locBin_PID == 0)
			{
				hInvMass_rho            -> Fill(locInvMass);
				hMissingMassSquared_rho -> Fill(locMissingMassSquared);
				hChiSq_rho              -> Fill(locChiSq);

				hTofTrackDist_pion   -> Fill(locTofTrackDist);
				hTofTrackDeltaT_pion -> Fill(locTofTrackDeltaT);
				hDcHits_pion         -> Fill(locDcHits);

				if (particle_bar_tree <= 23)
					hTrackOccupancy_pion->Fill(particle_x_tree,23-particle_bar_tree);
				else
					hTrackOccupancy_pion->Fill(particle_x_tree,particle_bar_tree);
			}
			if (locBin_PID == 1)
			{
				hInvMass_phi            -> Fill(locInvMass);
				hMissingMassSquared_phi -> Fill(locMissingMassSquared);
				hChiSq_phi              -> Fill(locChiSq);

				hTofTrackDist_kaon   -> Fill(locTofTrackDist);
				hTofTrackDeltaT_kaon -> Fill(locTofTrackDeltaT);
				hDcHits_kaon         -> Fill(locDcHits);

				if (particle_bar_tree <= 23)
					hTrackOccupancy_kaon->Fill(particle_x_tree,23-particle_bar_tree);
				else
					hTrackOccupancy_kaon->Fill(particle_x_tree,particle_bar_tree);
					
			}


			if (locBin_PID == 0 && Ntracks_pion > 0 && track_pion_i >= Ntracks_pion)
				continue;
			if (locBin_PID == 1 && Ntracks_kaon > 0 && track_kaon_i >= Ntracks_kaon)
				continue;

			// ------- Applying Event/Track Quality Selection --------//
			if (locBin_PID == 0 && perform_selection)
			{
				locPassMMS = locMissingMassSquared > -0.004 && locMissingMassSquared < 0.003;
				if (!locPassMMS)
					continue;
			
				locPassInvMass = locInvMass > 0.66 && locInvMass < 0.84;
				if (!locPassInvMass)
					continue;

				locPassChiSq = locChiSq < 10.;
				if (!locPassChiSq)
					continue;

				locPassTofTrackDist = locTofTrackDist < 4.;
				if (!locPassTofTrackDist)
					continue;

				locPassTofTrackDeltaT = locTofTrackDeltaT > -0.5 && locTofTrackDeltaT < 0.5;
				if (!locPassTofTrackDeltaT)
					continue;

				locPassDcHits = locDcHits > 30;
				if (!locPassDcHits)
					continue;
			}
			if (locBin_PID == 1 && perform_selection)
			{
				locPassMMS = locMissingMassSquared > -0.003 && locMissingMassSquared < 0.002;
				if (!locPassMMS)
					continue;
			
				locPassInvMass = locInvMass > 1.012 && locInvMass < 1.028;
				if (!locPassInvMass)
					continue;

				locPassChiSq = locChiSq < 15.;
				if (!locPassChiSq)
					continue;

				locPassTofTrackDist = locTofTrackDist < 4.;
				if (!locPassTofTrackDist)
					continue;

				locPassTofTrackDeltaT = locTofTrackDeltaT > -0.5 && locTofTrackDeltaT < 0.5;
				if (!locPassTofTrackDeltaT)
					continue;

				locPassDcHits = locDcHits > 30;
				if (!locPassDcHits)
					continue;
			}

			// ------- Applying Track Kinematics Selection  --------//

			locBin_P = -1;
			for (int locPBin = 0 ; locPBin < NumPBins ; locPBin++ )
			{
				if (particle_mom_tree > PBinEdges[locPBin][0] && particle_mom_tree < PBinEdges[locPBin][1])
				{
					locBin_P = locPBin;
					break;
				}
			}
			locBin_P = 0; // temporary hack

			locBin_BB = -1;

			if (locBin_P == -1 && perform_selection)
				continue;


			locBin_Bar = particle_bar_tree;
			locBin_X   = hTrackOccupancy_pion -> GetXaxis() -> FindBin(particle_x_tree) - 1;

                        if (locBin_Bar >= 3 && locBin_Bar <= 8 && locBin_X >= 18 && locBin_X <= 21)
                                locSelectedRegion = true;
                        if (locBin_Bar >= 27 && locBin_Bar <= 32 && locBin_X >= 18 && locBin_X <= 21)
                                locSelectedRegion = true;


			if (particle_bar_tree <= 23)
				locSouthLower = true;
			else
				locSouthLower = false;


			if (0 <= particle_bar_tree && particle_bar_tree <= 11)
				locBin_BB = 0;
			if (12 <= particle_bar_tree && particle_bar_tree <= 23)
				locBin_BB = 1;
			if (24 <= particle_bar_tree && particle_bar_tree <= 35)
				locBin_BB = 2;
			if (36 <= particle_bar_tree && particle_bar_tree <= 47)
				locBin_BB = 3;


			locSelectedRegion = true;

	                locGoodParticle = dirc_model->convert_particle_kinematics_with_offsets(\
								particle_x,\
								particle_y,\
								particle_theta,\
								particle_phi,\
								particle_bar,\
								particle_x3_tree,\
								particle_p3_tree);

			particle_momentum = particle_mom_tree;

	                if (!locGoodParticle)
        	                continue;
			if (locSouthLower && particle_bar_tree != int(particle_bar))
				continue;
			if (!locSouthLower && particle_bar_tree != int(47-particle_bar))
				continue;

			// ------ Passed All Selection -------- //

			counter++;

			//------------ APPLY OFFSETS SOUTH BOX --------------//

			if (locSouthLower)
			{

				for (int seg_i = 0 ; seg_i < 3 ; seg_i++)
				{
					dirc_model->set_offsets_threeseg_mirror_angle(seg_i+1,\
										threeSeg_x_angle_offs_SouthLower[seg_i],\
										threeSeg_y_angle_offs_SouthLower[seg_i],\
										threeSeg_z_angle_offs_SouthLower[seg_i]);
				
					dirc_model->set_offsets_threeseg_mirror_position(seg_i+1,\
										threeSeg_y_pos_offs_SouthLower[seg_i],\
										threeSeg_z_pos_offs_SouthLower[seg_i]);

				}

				dirc_model->set_offsets_BB_OB_angle(\
									BB_OB_x_angle_off_SouthLower,\
									BB_OB_y_angle_off_SouthLower,\
									BB_OB_z_angle_off_SouthLower);
				dirc_model->set_offsets_BB_OB_position(\
									BB_OB_x_pos_off_SouthLower,\
									BB_OB_y_pos_off_SouthLower,\
									BB_OB_z_pos_off_SouthLower);

				dirc_model->set_offsets_pmt_plane_angle_y_z(\
									pmt_angle_off_SouthLower,\
									pmt_y_off_SouthLower,\
									pmt_z_off_SouthLower);

				dirc_model->set_offsets_sidemirror_x_offs(\
									sidemirror_xr_x_off_SouthLower,\
									sidemirror_xl_x_off_SouthLower);

				dirc_model->set_offsets_largePlanarMirrorY_y_off(\
									largePlanarMirrorY_y_off_SouthLower);

				dirc_model->set_offsets_upperWedgeFarZ_z_off(\
									upperWedgeFarZ_z_off_SouthLower);


				dirc_model->set_offsets_support_points_off(\
									support_x_off_SouthLower,\
									support_y_off_SouthLower,\
									support_t_off_SouthLower);

				dirc_model->set_offsets_BB_pos(\
								10,\
								BB10_offsets_x_pos,\
								BB10_offsets_y_pos,\
								BB10_offsets_z_pos);
				dirc_model->set_offsets_BB_pos(\
								11,\
								BB11_offsets_x_pos,\
								BB11_offsets_y_pos,\
								BB11_offsets_z_pos);
				dirc_model->set_offsets_BB_angle(\
								10,\
								BB10_offsets_x_angle,\
								BB10_offsets_y_angle,\
								BB10_offsets_z_angle);
				dirc_model->set_offsets_BB_angle(\
								11,\
								BB11_offsets_x_angle,\
								BB11_offsets_y_angle,\
								BB11_offsets_z_angle);

			}
			else
			{
				for (int seg_i = 0 ; seg_i < 3 ; seg_i++)
				{
					dirc_model->set_offsets_threeseg_mirror_angle(seg_i+1,\
										threeSeg_x_angle_offs_NorthUpper[seg_i],\
										threeSeg_y_angle_offs_NorthUpper[seg_i],\
										threeSeg_z_angle_offs_NorthUpper[seg_i]);
				
					dirc_model->set_offsets_threeseg_mirror_position(seg_i+1,\
										threeSeg_y_pos_offs_NorthUpper[seg_i],\
										threeSeg_z_pos_offs_NorthUpper[seg_i]);

				}

				dirc_model->set_offsets_BB_OB_angle(\
									BB_OB_x_angle_off_NorthUpper,\
									BB_OB_y_angle_off_NorthUpper,\
									BB_OB_z_angle_off_NorthUpper);
				dirc_model->set_offsets_BB_OB_position(\
									BB_OB_x_pos_off_NorthUpper,\
									BB_OB_y_pos_off_NorthUpper,\
									BB_OB_z_pos_off_NorthUpper);

				dirc_model->set_offsets_pmt_plane_angle_y_z(\
									pmt_angle_off_NorthUpper,\
									pmt_y_off_NorthUpper,\
									pmt_z_off_NorthUpper);

				dirc_model->set_offsets_sidemirror_x_offs(\
									sidemirror_xr_x_off_NorthUpper,\
									sidemirror_xl_x_off_NorthUpper);

				dirc_model->set_offsets_largePlanarMirrorY_y_off(\
									largePlanarMirrorY_y_off_NorthUpper);

				dirc_model->set_offsets_upperWedgeFarZ_z_off(\
									upperWedgeFarZ_z_off_NorthUpper);


				dirc_model->set_offsets_support_points_off(\
									support_x_off_NorthUpper,\
									support_y_off_NorthUpper,\
									support_t_off_NorthUpper);

				dirc_model->set_offsets_BB_pos(\
								0,\
								BB00_offsets_x_pos,\
								BB00_offsets_y_pos,\
								BB00_offsets_z_pos);
				dirc_model->set_offsets_BB_pos(\
								1,\
								BB01_offsets_x_pos,\
								BB01_offsets_y_pos,\
								BB01_offsets_z_pos);
				dirc_model->set_offsets_BB_angle(\
								0,\
								BB00_offsets_x_angle,\
								BB00_offsets_y_angle,\
								BB00_offsets_z_angle);
				dirc_model->set_offsets_BB_angle(\
								1,\
								BB01_offsets_x_angle,\
								BB01_offsets_y_angle,\
								BB01_offsets_z_angle);

			}


			// ----------- Get Hit info ------------ //

                	particle_t_tree        = evt->GetTime();


			if (locSouthLower)
			{
				row_factor = 1;
				multi = 0; 
			}
			else
			{
				row_factor = -1;
				multi = 1; 
			}

			for (auto hit : evt->GetHits())
			{
				dirc_point tree_point;
				hit_ChannelId   = hit.GetChannel();
				hit_t           = hit.GetLeadTime() - particle_t_tree;
				hit_pixelrow    = digitizer.GetPixelRow(hit_ChannelId);
				hit_pixelcol    = digitizer.GetPixelColumn(hit_ChannelId);

				hit_dist_t_tree      -> Fill(hit_t);

				tree_point.pixel_row = int(143 * multi + row_factor * hit_pixelrow);
				tree_point.pixel_col = hit_pixelcol;
				tree_point.t         = hit_t;

				if (hit_t < 38)
					tree_point.direct = true;
				else
					tree_point.direct = false;


				digitizer.undigitize_point(tree_point);

				tree_points.push_back(tree_point);

			}

			locBin_hittype = 0; // data

			
			Nph_total_tree = int(evt->GetHitSize());
	                Nph_tree -> Fill(Nph_total_tree);
	                Nph_tree_BB[locBin_BB] -> Fill(Nph_total_tree);
			if (locSouthLower)
	                	Nph_tree_SouthLower -> Fill(Nph_total_tree);
			else
	                	Nph_tree_NorthUpper -> Fill(Nph_total_tree);


			// -------- Generate PDFs -------------------//

			pion_beta = dirc_model->get_beta_from_p(particle_momentum,pimass);
			kaon_beta = dirc_model->get_beta_from_p(particle_momentum,kmass);


			timing_clock = clock();

			dirc_model->sim_reg_n_photons(\
					support_points_pion,\
					n_phi_phots,\
					n_z_phots,\
					-1,\
					particle_bar,\
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
					particle_bar,\
					particle_x,\
					particle_y,\
					0.,\
					particle_theta,\
					particle_phi,\
					0,\
					ckov_unc/pdf_unc_red_fac,\
					kaon_beta);

			timing_clock = clock() - timing_clock;
			gen_time = (float) timing_clock/(CLOCKS_PER_SEC)*1000.; // this is microsecond

			hTimePerTrack->Fill(gen_time);

			// ----------- Calulate DLL -------- //


			if (smear_time == true)
			{
				for (unsigned int loc_i = 0 ; loc_i < support_points_pion.size(); loc_i++ )
				{
					support_points_pion[loc_i].t += rand_gen->Gaus(0.,sigma_t_pion_smear);
				}
				for (unsigned int loc_i = 0 ; loc_i < support_points_kaon.size(); loc_i++ )
				{
					support_points_kaon[loc_i].t += rand_gen->Gaus(0.,sigma_t_kaon_smear);
				}
			}


	                ll_pion = pdf_pion->get_log_likelihood_distances(support_points_pion,tree_points,\
							hMap_Delta[0][0][locBin_hittype][locBin_BB],\
							hMap_Delta[1][0][locBin_hittype][locBin_BB],\
							hMap_Delta[2][0][locBin_hittype][locBin_BB],\
							hMap_Delta[3][0][locBin_hittype][locBin_BB],\
							hMap_Delta[4][0][locBin_hittype][locBin_BB],\
							hMap_Delta_Bar_XBin[locBin_P][locBin_Bar][locBin_X][0][0][locBin_hittype][0],\
							hMap_Delta_Bar_XBin[locBin_P][locBin_Bar][locBin_X][0][0][locBin_hittype][1]);

			int Ngoodhit_counter_pion = 0;

			for (unsigned int loc_i = 0; loc_i < tree_points.size(); loc_i++)
			{

				for (int loc_j = 0 ; loc_j < int(obs_labels.size()) ; loc_j++)
				{
					hMap_MinDist[loc_j][0][locBin_hittype][locBin_BB]->Fill(tree_points[loc_i].min_dist[loc_j]);
					hMap_FracClose[loc_j][0][locBin_hittype][locBin_BB]->Fill(tree_points[loc_i].frac_close[loc_j]);
				}


                                if (tree_points[loc_i].cut==false)
                                {
                                        Ngoodhit_counter_pion++;

                                        if (locBin_PID == 0)
                                        {
                                                hMap_hit_dist_t_pion_tree[locBin_P][locBin_Bar][locBin_X]->Fill(tree_points[loc_i].t);
                                                hMap_hit_dist_rowcol_pion_near_tree[locBin_P][locBin_Bar][locBin_X]->Fill(\
                                                                                                tree_points[loc_i].pixel_row,\
                                                                                                tree_points[loc_i].pixel_col);
                                        }
                                }
                                else
                                {
                                        if (locBin_PID == 0)
                                        {
                                                //hMap_hit_dist_t_pion_away_tree[locBin_P][locBin_Bar][locBin_X]->Fill(tree_points[i].t);
                                                hMap_hit_dist_rowcol_pion_away_tree[locBin_P][locBin_Bar][locBin_X]->Fill(\
                                                                                                tree_points[loc_i].pixel_row,\
                                                                                                tree_points[loc_i].pixel_col);
                                        }
                                }

			}

			if (locBin_PID == 0)
			{
                                Nph_cut_tree->Fill(Ngoodhit_counter_pion);
                                Nph_cut_pion_tree->Fill(Ngoodhit_counter_pion);

                                hMap_Nph_pion_cut_tree[locBin_P][locBin_Bar][locBin_X] -> Fill(Ngoodhit_counter_pion);
                                
                                Nph_cut_tree_BB[locBin_BB]->Fill(Ngoodhit_counter_pion);
				Nph_cut_pion_tree_BB[locBin_BB]->Fill(Ngoodhit_counter_pion);

				if (locSouthLower)
				{
                                	Nph_cut_tree_SouthLower->Fill(Ngoodhit_counter_pion);
                                	Nph_cut_pion_tree_SouthLower->Fill(Ngoodhit_counter_pion);
				}
				else
				{
                                	Nph_cut_tree_NorthUpper->Fill(Ngoodhit_counter_pion);
                                	Nph_cut_pion_tree_NorthUpper->Fill(Ngoodhit_counter_pion);
				}
			}

	                ll_kaon = pdf_kaon->get_log_likelihood_distances(support_points_kaon,tree_points,\
							hMap_Delta[0][1][locBin_hittype][locBin_BB],\
							hMap_Delta[1][1][locBin_hittype][locBin_BB],\
							hMap_Delta[2][1][locBin_hittype][locBin_BB],\
							hMap_Delta[3][1][locBin_hittype][locBin_BB],\
							hMap_Delta[4][1][locBin_hittype][locBin_BB],\
							hMap_Delta_Bar_XBin[locBin_P][locBin_Bar][locBin_X][0][1][locBin_hittype][0],\
							hMap_Delta_Bar_XBin[locBin_P][locBin_Bar][locBin_X][0][1][locBin_hittype][1]);
			int Ngoodhit_counter_kaon = 0;
			for (unsigned int loc_i = 0; loc_i < tree_points.size(); loc_i++)
			{
				for (int loc_j = 0 ; loc_j < int(obs_labels.size()) ; loc_j++)
				{
					hMap_MinDist[loc_j][1][locBin_hittype][locBin_BB]->Fill(tree_points[loc_i].min_dist[loc_j]);
					hMap_FracClose[loc_j][1][locBin_hittype][locBin_BB]->Fill(tree_points[loc_i].frac_close[loc_j]);
				
				}
				if (tree_points[loc_i].cut==false)
				{
                                        Ngoodhit_counter_kaon++;
                                        hit_dist_t_cut_tree->Fill(tree_points[loc_i].t);

                                        if (locBin_PID == 1)
                                        {
                                                hMap_hit_dist_t_kaon_tree[locBin_P][locBin_Bar][locBin_X]->Fill(tree_points[loc_i].t);
                                                hMap_hit_dist_rowcol_kaon_near_tree[locBin_P][locBin_Bar][locBin_X]->Fill(\
                                                                                                tree_points[loc_i].pixel_row,\
                                                                                                tree_points[loc_i].pixel_col);
                                        }
				}
                                else
                                {
                                        if (locBin_PID == 1)
                                        {
                                                //hMap_hit_dist_t_kaon_away_tree[locBin_P][locBin_Bar][locBin_X]->Fill(tree_points[i].t);
                                                hMap_hit_dist_rowcol_kaon_away_tree[locBin_P][locBin_Bar][locBin_X]->Fill(\
                                                                                                tree_points[loc_i].pixel_row,\
                                                                                                tree_points[loc_i].pixel_col);
                                        }
                                }

			}

                        if (locBin_PID == 1)
                        {
                                Nph_cut_tree->Fill(Ngoodhit_counter_kaon);
                                Nph_cut_kaon_tree->Fill(Ngoodhit_counter_kaon);

                                hMap_Nph_kaon_cut_tree[locBin_P][locBin_Bar][locBin_X] -> Fill(Ngoodhit_counter_kaon);

                                Nph_cut_tree_BB[locBin_BB]->Fill(Ngoodhit_counter_kaon);
				Nph_cut_kaon_tree_BB[locBin_BB]->Fill(Ngoodhit_counter_kaon);

				if (locSouthLower)
				{
                                	Nph_cut_tree_SouthLower->Fill(Ngoodhit_counter_kaon);
                                	Nph_cut_kaon_tree_SouthLower->Fill(Ngoodhit_counter_kaon);
				}
				else
				{
                                	Nph_cut_tree_NorthUpper->Fill(Ngoodhit_counter_kaon);
                                	Nph_cut_kaon_tree_NorthUpper->Fill(Ngoodhit_counter_kaon);
				}
                        }



			if (locBin_PID == 0)
			{
                		ll_diff_pion_tree->Fill(ll_pion-ll_kaon);

				hDLLvsNph_pion_tree -> Fill(Nph_total_tree,ll_pion-ll_kaon); 
				hDLLvsNphNear_pion_tree -> Fill(Ngoodhit_counter_pion,ll_pion-ll_kaon); 

				if (particle_bar_tree <= 23)
                                	hMap_TrackOccupancy_pion[locBin_P] -> Fill(particle_x_tree,23-particle_bar_tree);
				else
                                	hMap_TrackOccupancy_pion[locBin_P] -> Fill(particle_x_tree,particle_bar_tree);


                		ll_diff_pion_tree_BB[locBin_BB]->Fill(ll_pion-ll_kaon);

				if (locSouthLower)
				{
					if (locSelectedRegion)
                				ll_diff_pion_tree_SouthLower->Fill(ll_pion-ll_kaon);
				}
				else
				{
					if (locSelectedRegion)
                				ll_diff_pion_tree_NorthUpper->Fill(ll_pion-ll_kaon);
				}
			}
			if (locBin_PID == 1)
			{
                		ll_diff_kaon_tree->Fill(ll_pion-ll_kaon);

				hDLLvsNph_kaon_tree -> Fill(Nph_total_tree,ll_pion-ll_kaon); 
				hDLLvsNphNear_kaon_tree -> Fill(Ngoodhit_counter_kaon,ll_pion-ll_kaon); 

				if (particle_bar_tree <= 23)
                                	hMap_TrackOccupancy_kaon[locBin_P] -> Fill(particle_x_tree,23-particle_bar_tree);
				else
                                	hMap_TrackOccupancy_kaon[locBin_P] -> Fill(particle_x_tree,particle_bar_tree);

                		ll_diff_kaon_tree_BB[locBin_BB]->Fill(ll_pion-ll_kaon);

				if (locSouthLower)
				{
					if (locSelectedRegion)
                				ll_diff_kaon_tree_SouthLower->Fill(ll_pion-ll_kaon);
				}
				else
				{
					if (locSelectedRegion)
                				ll_diff_kaon_tree_NorthUpper->Fill(ll_pion-ll_kaon);
				}
			}


			// Studies with simulated hits

			locBin_hittype = 1;

                        if (locBin_PID == 0)
                        {
                                dirc_model->sim_rand_n_photons(\
                                                sim_points_pion,\
                                                n_sim_phots,\
                                                -1,\
                                                particle_bar,\
                                                particle_x,\
                                                particle_y,\
                                                0.,\
                                                particle_theta,\
                                                particle_phi,\
                                                0.,\
                                                ckov_unc,\
                                                pion_beta);

                                digitizer.digitize_points(sim_points_pion);


                                for (unsigned int loc_i = 0 ; loc_i < sim_points_pion.size(); loc_i++ )
                                {
                                        if (sim_points_pion[loc_i].pixel_row < 0 || sim_points_pion[loc_i].pixel_col < 0)
                                        {
                                                sim_points_pion.erase(sim_points_pion.begin() + loc_i);
                                        }
                                        else
                                        {
                                                hMap_hit_dist_t_pion_sim[locBin_P][locBin_Bar][locBin_X]->Fill(sim_points_pion[loc_i].t);
                                                hMap_hit_dist_rowcol_pion_sim[locBin_P][locBin_Bar][locBin_X]->Fill(\
                                                                                                sim_points_pion[loc_i].pixel_row,\
                                                                                                sim_points_pion[loc_i].pixel_col);
						if (sim_points_pion[loc_i].t < 38)
							sim_points_pion[loc_i].direct = true;
						else 
							sim_points_pion[loc_i].direct = false;


                                        }

                                }
	                
				Ngoodhit_counter_pion = 0;

				ll_pion = pdf_pion->get_log_likelihood_distances(support_points_pion,sim_points_pion,\
							hMap_Delta[0][0][locBin_hittype][locBin_BB],\
							hMap_Delta[1][0][locBin_hittype][locBin_BB],\
							hMap_Delta[2][0][locBin_hittype][locBin_BB],\
							hMap_Delta[3][0][locBin_hittype][locBin_BB],\
							hMap_Delta[4][0][locBin_hittype][locBin_BB],\
							hMap_Delta_Bar_XBin[locBin_P][locBin_Bar][locBin_X][0][0][locBin_hittype][0],\
							hMap_Delta_Bar_XBin[locBin_P][locBin_Bar][locBin_X][0][0][locBin_hittype][1]);

				for (int loc_i = 0 ; loc_i < int(sim_points_pion.size()); loc_i++)
				{
					for (int loc_j = 0 ; loc_j < int(obs_labels.size()) ; loc_j++)
					{
						hMap_MinDist[loc_j][0][locBin_hittype][locBin_BB]->Fill(sim_points_pion[loc_i].min_dist[loc_j]);
						hMap_FracClose[loc_j][0][locBin_hittype][locBin_BB]->Fill(sim_points_pion[loc_i].frac_close[loc_j]);
					}
					if (sim_points_pion[loc_i].cut == false)
						Ngoodhit_counter_pion++;

				}
                                hMap_Nph_pion_cut_sim[locBin_P][locBin_Bar][locBin_X] -> Fill(Ngoodhit_counter_pion);
				ll_kaon = pdf_kaon->get_log_likelihood_distances(support_points_pion,sim_points_pion,\
							hMap_Delta[0][1][locBin_hittype][locBin_BB],\
							hMap_Delta[1][1][locBin_hittype][locBin_BB],\
							hMap_Delta[2][1][locBin_hittype][locBin_BB],\
							hMap_Delta[3][1][locBin_hittype][locBin_BB],\
							hMap_Delta[4][1][locBin_hittype][locBin_BB],\
							hMap_Delta_Bar_XBin[locBin_P][locBin_Bar][locBin_X][0][1][locBin_hittype][0],\
							hMap_Delta_Bar_XBin[locBin_P][locBin_Bar][locBin_X][0][1][locBin_hittype][1]);

				for (int loc_i = 0 ; loc_i < int(sim_points_pion.size()); loc_i++)
				{
					for (int loc_j = 0 ; loc_j < int(obs_labels.size()) ; loc_j++)
					{
						hMap_MinDist[loc_j][1][locBin_hittype][locBin_BB]->Fill(sim_points_pion[loc_i].min_dist[loc_j]);
						hMap_FracClose[loc_j][1][locBin_hittype][locBin_BB]->Fill(sim_points_pion[loc_i].frac_close[loc_j]);
					}

				}

				hDLLvsNph_pion_sim -> Fill(int(sim_points_pion.size()),ll_pion-ll_kaon); 

                		ll_diff_pion_sim->Fill(ll_pion-ll_kaon);
				if (locSouthLower)
				{
					if (locSelectedRegion)
                				ll_diff_pion_sim_SouthLower->Fill(ll_pion-ll_kaon);
				}
				else
				{
					if (locSelectedRegion)
                				ll_diff_pion_sim_NorthUpper->Fill(ll_pion-ll_kaon);
				}
                		ll_diff_pion_sim_BB[locBin_BB]->Fill(ll_pion-ll_kaon);
                        }
 
                        if (locBin_PID == 1)
                        {
                                dirc_model->sim_rand_n_photons(\
                                                sim_points_kaon,\
                                                n_sim_phots,\
                                                -1,\
                                                particle_bar,\
                                                particle_x,\
                                                particle_y,\
                                                0.,\
                                                particle_theta,\
                                                particle_phi,\
                                                0.,\
                                                ckov_unc,\
                                                kaon_beta);
                                digitizer.digitize_points(sim_points_kaon);
                                for (unsigned int loc_i = 0 ; loc_i < sim_points_kaon.size(); loc_i++ )
                                {
                                        if (sim_points_kaon[loc_i].pixel_row < 0 || sim_points_kaon[loc_i].pixel_col < 0)
                                        {
                                                sim_points_kaon.erase(sim_points_kaon.begin() + loc_i);
                                        }
                                        else
                                        {
                                                hMap_hit_dist_t_kaon_sim[locBin_P][locBin_Bar][locBin_X]->Fill(sim_points_kaon[loc_i].t);
                                                hMap_hit_dist_rowcol_kaon_sim[locBin_P][locBin_Bar][locBin_X]->Fill(\
                                                                                                sim_points_kaon[loc_i].pixel_row,\
                                                                                                sim_points_kaon[loc_i].pixel_col);
						if (sim_points_kaon[loc_i].t < 38)
							sim_points_kaon[loc_i].direct = true;
						else 
							sim_points_kaon[loc_i].direct = false;
                                        }

                                }
	                
				ll_pion = pdf_pion->get_log_likelihood_distances(support_points_pion,sim_points_kaon,\
							hMap_Delta[0][0][locBin_hittype][locBin_BB],\
							hMap_Delta[1][0][locBin_hittype][locBin_BB],\
							hMap_Delta[2][0][locBin_hittype][locBin_BB],\
							hMap_Delta[3][0][locBin_hittype][locBin_BB],\
							hMap_Delta[4][0][locBin_hittype][locBin_BB],\
							hMap_Delta_Bar_XBin[locBin_P][locBin_Bar][locBin_X][0][0][locBin_hittype][0],\
							hMap_Delta_Bar_XBin[locBin_P][locBin_Bar][locBin_X][0][0][locBin_hittype][1]);

				for (int loc_i = 0 ; loc_i < int(sim_points_kaon.size()); loc_i++)
				{
					for (int loc_j = 0 ; loc_j < int(obs_labels.size()) ; loc_j++)
					{
						hMap_MinDist[loc_j][0][locBin_hittype][locBin_BB]->Fill(sim_points_kaon[loc_i].min_dist[loc_j]);
						hMap_FracClose[loc_j][0][locBin_hittype][locBin_BB]->Fill(sim_points_kaon[loc_i].frac_close[loc_j]);
					}

				}
				ll_kaon = pdf_kaon->get_log_likelihood_distances(support_points_pion,sim_points_kaon,\
							hMap_Delta[0][1][locBin_hittype][locBin_BB],\
							hMap_Delta[1][1][locBin_hittype][locBin_BB],\
							hMap_Delta[2][1][locBin_hittype][locBin_BB],\
							hMap_Delta[3][1][locBin_hittype][locBin_BB],\
							hMap_Delta[4][1][locBin_hittype][locBin_BB],\
							hMap_Delta_Bar_XBin[locBin_P][locBin_Bar][locBin_X][0][1][locBin_hittype][0],\
							hMap_Delta_Bar_XBin[locBin_P][locBin_Bar][locBin_X][0][1][locBin_hittype][1]);

				Ngoodhit_counter_kaon = 0;
				for (int loc_i = 0 ; loc_i < int(sim_points_kaon.size()); loc_i++)
				{
					for (int loc_j = 0 ; loc_j < int(obs_labels.size()) ; loc_j++)
					{
						hMap_MinDist[loc_j][1][locBin_hittype][locBin_BB]->Fill(sim_points_kaon[loc_i].min_dist[loc_j]);
						hMap_FracClose[loc_j][1][locBin_hittype][locBin_BB]->Fill(sim_points_kaon[loc_i].frac_close[loc_j]);
					}
					if (sim_points_kaon[loc_i].cut == false)
						Ngoodhit_counter_kaon++;

				}
                                hMap_Nph_kaon_cut_sim[locBin_P][locBin_Bar][locBin_X] -> Fill(Ngoodhit_counter_kaon);
				hDLLvsNph_kaon_sim -> Fill(int(sim_points_kaon.size()),ll_pion-ll_kaon); 
                		ll_diff_kaon_sim->Fill(ll_pion-ll_kaon);
				if (locSouthLower)
				{
					if (locSelectedRegion)
                				ll_diff_kaon_sim_SouthLower->Fill(ll_pion-ll_kaon);
				}
				else
				{
					if (locSelectedRegion)
                				ll_diff_kaon_sim_NorthUpper->Fill(ll_pion-ll_kaon);
				}
                		ll_diff_kaon_sim_BB[locBin_BB]->Fill(ll_pion-ll_kaon);
                        }




			//support points
                        for (unsigned int loc_i = 0 ; loc_i < support_points_pion.size(); loc_i++ )
                        {
                                if (loc_i % scale_support!=0)
                                        continue;
                                digitizer.digitize_point(support_points_pion[loc_i]);

                                if (support_points_pion[loc_i].pixel_row < 0 || support_points_pion[loc_i].pixel_col < 0)
                                {
                                        support_points_pion.erase(support_points_pion.begin() + loc_i);
                                }
                                else
                                {
                                        //hMap_support_dist_t_pion_sim[locBin_P][locBin_Bar][locBin_X]->Fill(support_points_pion[loc_i].t);
                                        hMap_support_dist_rowcol_pion_sim[locBin_P][locBin_Bar][locBin_X]->Fill(\
                                                             int(143*multi + row_factor*support_points_pion[loc_i].pixel_row),\
                                                                                        support_points_pion[loc_i].pixel_col);
                                }
                        }
                        for (unsigned int loc_i = 0 ; loc_i < support_points_kaon.size(); loc_i++ )
                        {
                                if (loc_i % scale_support!=0)
                                        continue;
                                digitizer.digitize_point(support_points_kaon[loc_i]);

                                if (support_points_kaon[loc_i].pixel_row < 0 || support_points_kaon[loc_i].pixel_col < 0)
                                {
                                        support_points_kaon.erase(support_points_kaon.begin() + loc_i);
                                }
                                else
                                {
                                        //hMap_support_dist_t_kaon_sim[locBin_P][locBin_Bar][locBin_X]->Fill(support_points_kaon[loc_i].t);
                                        hMap_support_dist_rowcol_kaon_sim[locBin_P][locBin_Bar][locBin_X]->Fill(\
                                                               int(143*multi+row_factor*support_points_kaon[loc_i].pixel_row),\
                                                                                        support_points_kaon[loc_i].pixel_col);
                                }

                        }



                        track_i++;
                        if (locBin_PID == 0)
                                track_pion_i++;
                        if (locBin_PID == 1)
                                track_kaon_i++;

			tree_points.clear();
			sim_points_pion.clear();
			sim_points_kaon.clear();
			support_points_pion.clear();
			support_points_kaon.clear();

		}//END of particle loop

	}//END of event loop


        tfile->cd();

	dir_integrated->Write();
        for (int locPBin = 0 ; locPBin < NumPBins; locPBin++)
        {
                PBin_dirs[locPBin] -> Write();
        }


	tfile->Close();

	return 1;


}
