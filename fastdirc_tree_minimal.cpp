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

  //int n_sim_phots = 40;
  int rseed = 1337;

  double pion_beta, kaon_beta;
  pion_beta=kaon_beta=-1;


  double s_func_x = 6;
  double s_func_y = s_func_x;
  double s_func_t = 1.0;
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

  //--------------- END of initialization ------------------//

  //-------------- READ IN CONFIG ------------------------//
  if (nargs==2) config_str = argv[1];
  else config_str = "config.in";
  printf("\nRunning with config file: %s\n", config_str);

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
  tfile->cd();

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

  TH1F *hMomentum = new TH1F("hMomentum","track momentum; p [GeV];", 100,0.,5.);

  // hit time
  TH1F *hit_dist_t_tree = new TH1F("hit_dist_t_tree","hit time ; hit time (ns);",800,0,800);
  TH1F *hit_dist_t_cut_pionhypo_tree = new TH1F("hit_dist_t_cut_pionhypo_tree","hit time ; hit time (ns);",800,0,800);
  TH1F *hit_dist_t_cut_kaonhypo_tree = new TH1F("hit_dist_t_cut_kaonhypo_tree","hit time ; hit time (ns);",800,0,800);

  // likelihood
  TH1F *ll_diff_pion_tree = new TH1F("ll_diff_pion_tree","Difference of log likelihood, pion selection",200000,-200,200);
  TH1F *ll_diff_kaon_tree = new TH1F("ll_diff_kaon_tree","Difference of log likelihood, kaon selection",200000,-200,200);

  // photon yield
  TH1F *Nph_tree          = new TH1F("Nph_tree","number of hits", 1001,-.5,1000.5);
  TH1F *Nph_cut_tree      = new TH1F("Nph_cut_tree","number of hits (selection)", 1001,-.5,1000.5);
  TH1F *Nph_cut_pion_tree = new TH1F("Nph_cut_pion_tree","number of photons (pion selection)", 1001,-.5,1000.5);
  TH1F *Nph_cut_kaon_tree = new TH1F("Nph_cut_kaon_tree","number of photons (kaon selection)", 1001,-.5,1000.5);


  //-------------- INITIALIZE TREE --------------------//

  bool locSouthLower = true;
  bool locSelectedRegion = false;

  int locBin_PID = -1;
  int locBin_X   = -1;
  int locBin_Bar = -1;

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
  printf("\nTotal Entries = %d\n\n",int(event_chain->GetEntries()));

  //------------- INITIALIZE FastDIRC -----------------//

  std::vector<dirc_point> support_points_pion;
  std::vector<dirc_point> support_points_kaon;

  DircSpreadGaussianNew* KDE_calc = new DircSpreadGaussianNew(\
    sfunc_sig,\
    s_func_x,\
    s_func_y,\
    s_func_t,\
    sigma_cut);

  int track_i      = 0;
  int track_pion_i = 0;
  int track_kaon_i = 0;

  int Nph_total_tree = -1;

  for (int event_i = 0 ; event_i < event_chain->GetEntries() ; event_i++){

    event_chain->GetEntry(event_i);
    if(event_i%500==0)
      printf("Event %d, # of pion tracks = %d, # of kaon tracks = %d\n",event_i,track_pion_i,track_kaon_i);

    if (Ntracks > 0 && track_i >= Ntracks)
            break;

    if (Ntracks_pion > 0 && track_pion_i >= Ntracks_pion && \
        Ntracks_kaon > 0 && track_kaon_i >= Ntracks_kaon)
            break;

    for (int particle_i = 0; particle_i < particle_array->GetEntriesFast(); particle_i++) {
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
      particle_x_tree        = particle_x3_tree.X();

      particle_mom_tree      = particle_p3_tree.Mag();

      if (abs(particle_pid_tree)==211)
      	locBin_PID = 0;
      else if (abs(particle_pid_tree)==321)
      	locBin_PID = 1;
      else
      	continue;

      if (locBin_PID == 0){
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
      if (locBin_PID == 1) {
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
      hMomentum -> Fill(particle_mom_tree);

      if (locBin_PID == 0 && Ntracks_pion > 0 && track_pion_i >= Ntracks_pion) continue;
      if (locBin_PID == 1 && Ntracks_kaon > 0 && track_kaon_i >= Ntracks_kaon) continue;

      // ------- Applying Selection  --------//
      locBin_Bar = particle_bar_tree;
      locBin_X   = hTrackOccupancy_pion -> GetXaxis() -> FindBin(particle_x_tree) - 1;

      if (locBin_Bar >= 3 && locBin_Bar <= 8 && locBin_X >= 18 && locBin_X <= 21) locSelectedRegion = true;
      if (locBin_Bar >= 27 && locBin_Bar <= 32 && locBin_X >= 18 && locBin_X <= 21) locSelectedRegion = true;

      if (particle_bar_tree <= 23) locSouthLower = true;
      else locSouthLower = false;

      locGoodParticle = dirc_model->convert_particle_kinematics_with_offsets(\
                                                      			particle_x,\
                                                      			particle_y,\
                                                      			particle_theta,\
                                                      			particle_phi,\
                                                      			particle_bar,\
                                                      			particle_x3_tree,\
                                                      			particle_p3_tree);
      particle_momentum = particle_mom_tree;

      if (!locGoodParticle) continue;
      if (locSouthLower && particle_bar_tree != int(particle_bar)) continue;
      if (!locSouthLower && particle_bar_tree != int(particle_bar+24)) continue;

      // ------ Passed All Selection -------- //
      counter++;

      // -------- Generate support points -------------------//

      pion_beta = dirc_model->get_beta_from_p(particle_momentum,pimass);
      kaon_beta = dirc_model->get_beta_from_p(particle_momentum,kmass);

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


      // ----------- Get Hit info ------------ //

      particle_t_tree        = evt->GetTime();
      for (auto hit : evt->GetHits()){
        dirc_point tree_point;
        hit_ChannelId   = hit.GetChannel();
        hit_t           = hit.GetLeadTime() - particle_t_tree;
        hit_pixelrow    = digitizer.GetPixelRow(hit_ChannelId);
        hit_pixelcol    = digitizer.GetPixelColumn(hit_ChannelId);

        hit_dist_t_tree      -> Fill(hit_t);

        tree_point.pixel_row = hit_pixelrow;
        tree_point.pixel_col = hit_pixelcol;
        tree_point.t         = hit_t;

        digitizer.undigitize_point(tree_point);

        tree_points.push_back(tree_point);
      }

      Nph_total_tree = int(evt->GetHitSize());
      Nph_tree -> Fill(Nph_total_tree);

      // ----------- Calulate DLL -------- //
      ll_pion = KDE_calc->get_log_likelihood(support_points_pion,tree_points); // ll under pion hypothesis
      int Ngoodhit_counter_pion = 0;
      for (unsigned int i = 0; i < tree_points.size(); i++){
        if (tree_points[i].cut==false){
          Ngoodhit_counter_pion++;
          hit_dist_t_cut_pionhypo_tree->Fill(tree_points[i].t);
        }
      }
      if (locBin_PID == 0){
        Nph_cut_tree->Fill(Ngoodhit_counter_pion);
        Nph_cut_pion_tree->Fill(Ngoodhit_counter_pion);
      }

      ll_kaon = KDE_calc->get_log_likelihood(support_points_kaon,tree_points); // LL under kaon hypothesis
      int Ngoodhit_counter_kaon = 0;
      for (unsigned int i = 0; i < tree_points.size(); i++){
        if (tree_points[i].cut==false){
          Ngoodhit_counter_kaon++;
          hit_dist_t_cut_kaonhypo_tree->Fill(tree_points[i].t);
        }
      }
      if (locBin_PID == 1) {
        Nph_cut_tree->Fill(Ngoodhit_counter_kaon);
        Nph_cut_kaon_tree->Fill(Ngoodhit_counter_kaon);
      }

      if (locBin_PID == 0) ll_diff_pion_tree->Fill(ll_pion-ll_kaon); // if the track's "true" label is pion
      if (locBin_PID == 1) ll_diff_kaon_tree->Fill(ll_pion-ll_kaon); // if the track's "true" label is kaon

      // ----------- Finishing -------- //
      track_i++;
      if (locBin_PID == 0) track_pion_i++;
      if (locBin_PID == 1) track_kaon_i++;

      tree_points.clear();
      support_points_pion.clear();
      support_points_kaon.clear();

    }//END of particle loop

  }//END of event loop

  tfile->Write();
  tfile->Close();

  return 1;
}
