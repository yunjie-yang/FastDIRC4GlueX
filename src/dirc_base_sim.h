#include <vector>
#include <utility>
#include <TRandom3.h>

#include <TVector3.h>
#include <TH1F.h>
#include <TH2F.h>

#include "dirc_point.h"
#include "dirc_photon.h"

#define rad2deg 57.2958

#ifndef DIRC_BASE_SIM
#define DIRC_BASE_SIM 
struct dirc_base_sim_tracking_step
{
	//position at start of step
	double x;
	double y;
	double z;

	//radians
	double sin_theta;
	double cos_theta;
	double sin_phi;
	double cos_phi;
};
class DircBaseSim
{
protected:
	double barLength;
	double barWidth;
	double barDepth;
	double windowThickness;

	double distDCBR11DCBR12;
        double B00A_x, B12A_x;
	double B00A_y,B11A_y,B12A_y,B23A_y;
	double B00A_z,B12A_z;

	double distDCBR35DCBR36;
        double B24A_x, B36A_x;
	double B24A_y,B35A_y,B36A_y,B47A_y;
	double B24A_z,B36A_z;

	double BB10_offsets_x_pos;
	double BB10_offsets_y_pos;
	double BB10_offsets_z_pos;
	double BB10_offsets_x_angle;
	double BB10_offsets_y_angle;
	double BB10_offsets_z_angle;

	double BB11_offsets_x_pos;
	double BB11_offsets_y_pos;
	double BB11_offsets_z_pos;
	double BB11_offsets_x_angle;
	double BB11_offsets_y_angle;
	double BB11_offsets_z_angle;

	double BB00_offsets_x_pos;
	double BB00_offsets_y_pos;
	double BB00_offsets_z_pos;
	double BB00_offsets_x_angle;
	double BB00_offsets_y_angle;
	double BB00_offsets_z_angle;

	double BB01_offsets_x_pos;
	double BB01_offsets_y_pos;
	double BB01_offsets_z_pos;
	double BB01_offsets_x_angle;
	double BB01_offsets_y_angle;
	double BB01_offsets_z_angle;

	TVector3 BB10_x0Vec, BB10_normVec, BB10_xaxisVec, BB10_yaxisVec;
	TVector3 BB11_x0Vec, BB11_normVec, BB11_xaxisVec, BB11_yaxisVec;
	TVector3 BB00_x0Vec, BB00_normVec, BB00_xaxisVec, BB00_yaxisVec;
	TVector3 BB01_x0Vec, BB01_normVec, BB01_xaxisVec, BB01_yaxisVec;


	double wedgeTop;
	double wedgeWidthOff;
	double wedgeDepthOff;
	double wedgeFarAngle;
	double wedgeCloseAngle;
	double wedgeWidth;
	double wedgeDepthHigh;
	double wedgeHeight;
	double upperWedgeDepthHigh;
	double upperWedgeTop;
	double upperWedgeHeight;
	double upperWedgeBottom;

	double wedgeClosePlaneNx;
	double wedgeClosePlaneNy;
	double wedgeClosePlaneNz;
	double wedgeClosePlaneD;

	double upperWedgeClosePlaneNx;
	double upperWedgeClosePlaneNy;
	double upperWedgeClosePlaneNz;
	double upperWedgeClosePlaneD;
	double lowerWedgeExtensionZ;

	double upperWedgeGap;
	
	bool upperWedgeNonUniform;
	double upperWedgeNonUniformSpread;
	
	double wedgeFarPlaneNx;
	double wedgeFarPlaneNy;
	double wedgeFarPlaneNz;
	double wedgeFarPlaneD;
	
	double upperWedgeFarPlaneNx;
	double upperWedgeFarPlaneNy;
	double upperWedgeFarPlaneNz;
	double upperWedgeFarPlaneD;
	
	double upperWedgeFarZ;
	double upperWedgeFarZ_input;
	double upperWedgeFarZ_z_off;
	
	double box_angle_off_cval;
	double box_angle_off_sval;
	double bar_box_xoff;
	double bar_box_yoff;
	double bar_box_zoff;

	double BB_OB_offsets_x_angle;
	double BB_OB_offsets_y_angle;
	double BB_OB_offsets_z_angle;
	
	double BB_OB_offsets_x_angle_cval;
	double BB_OB_offsets_x_angle_sval;
	double BB_OB_offsets_y_angle_cval;
	double BB_OB_offsets_y_angle_sval;
	double BB_OB_offsets_z_angle_cval;
	double BB_OB_offsets_z_angle_sval;

	double BB_OB_offsets_x_pos;
	double BB_OB_offsets_y_pos;
	double BB_OB_offsets_z_pos;

	double support_x_off;
	double support_y_off;
	double support_t_off;


	double quartzIndex;
	double quartzLiquidY;

	int wedge_bounces;
	int lastWallX;
	int wedgeBeforeInterface;

	double moliereP;
	bool useMoliere;

	double liquidIndex;
	double liquidAbsorbtion;
	bool use_liquid_n;
	bool use_quartz_n_for_liquid;
      
	std::vector<double> dist_traveled;
	bool kaleidoscope_plot;
	
	bool store_refraction;
	std::vector<double> refraction_before;
	std::vector<double> refraction_after;

	bool store_bounces;
	std::vector<int> x_bounces;
	std::vector<int> z_bounces;
	std::vector<int> x_direct_bounces;
	std::vector<int> z_direct_bounces;
	std::vector<int> x_indirect_bounces;
	std::vector<int> z_indirect_bounces;


	std::vector<double> generated_theta;
	std::vector<double> generated_phi;
	std::vector<double> generated_z;
	std::vector<double> generated_wavelength;

	
	double min_transmittance,max_transmittance,sep_transmittance;
	int num_transmittance;
	std::vector<double> quartz_transmittance;
	TRandom3 *rand_gen;

	bool midLineMode;
	int midLineWedgeWallFlip;

	bool upperWedgeAngleStore;
	std::vector<double> upper_wedge_incident;	
	
	void build_system();
	void spread_wedge_mirror();

	bool quartz_transmission_mc(double R, double lambda);
	bool absorbtion_mc(double dx, double dy);
	
	void rotate_2d(double &x, double &y, double cval, double sval);

	void build_BB_plane();


	void fill_rand_phi(\
		std::vector<dirc_point> &ovals,\
		int n_photons, \
		double ckov_theta /*= 47*/, \
        	double particle_bar /*=0*/, \
		double particle_x /*= 0*/, \
		double particle_y /*= 0*/, \
		double particle_t /*= 0*/, \
		double particle_theta /*= 0*/, \
		double particle_phi /*= 0*/,\
		double phi_theta_unc /*= .0015*57.3*/,\
		double ckov_theta_unc /* = .0055*57.3*/,\
		double beta /* = -1*/,\
		int save_kin /* =-1*/);

	void fill_reg_phi(\
		std::vector<dirc_point> &fill_points,\
		int n_photons_phi, \
		int n_photons_z,\
		double ckov_theta /*= 47*/, \
	        double particle_bar /*=0*/, \
		double particle_x /*= 0*/, \
		double particle_y /*= 0*/, \
		double particle_t /*= 0*/, \
		double particle_theta /*= 0*/, \
		double particle_phi /*= 0*/,\
		double phi_theta_unc, /*= 0*/
		double ckov_theta_unc /* = 0*/,\
		double beta /* = -1*/,\
		int save_kin /* =-1*/);

	void fill_reg_phi_exp(\
		std::vector<dirc_point> &fill_points,\
		int n_photons_phi, \
		int n_photons_z,\
		double ckov_theta /*= 47*/, \
	        double particle_bar /*=0*/, \
		double particle_x /*= 0*/, \
		double particle_y /*= 0*/, \
		double particle_t /*= 0*/, \
		double particle_theta /*= 0*/, \
		double particle_phi /*= 0*/,\
		double phi_theta_unc, /*= 0*/
		double ckov_theta_unc /* = 0*/,\
		double beta /* = -1*/,\
		TH1F* h_dx,\
		TH2F* h_2D_1,\
		TH2F* h_2D_2,\
		TH2F* h_2D_3,\
		int save_kin /* =-1*/);

	double generate_cos_moliere_angle(\
                double rad_length);

	void bar_box_interface(\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz);

	void fill_moliere_tracking_steps(\
		std::vector<dirc_base_sim_tracking_step> &rsteps,\
		double &travel_distance,\
		double step_length,\
		double start_theta,\
		double start_phi,\
		double start_x,\
		double start_y,\
		double start_z);
		

	
	double get_quartz_n(double lambda);
	double get_liquid_n(double lambda);
	bool optical_interface_z(\
		double n1,\
		double n2,\
		double &dx,\
		double &dy,\
		double &dz);

	bool interface_BB_OB(\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz,\
		double n1,\
		double n2);

	double warp_ray(\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz,\
		double cos_critical_angle);
	double warp_wedge(\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz);
	bool x_wedge_coerce_check(\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz,\
		double dt);
	//Utility function - should definitely be inlined
	void plane_reflect(\
		double Nx,\
		double Ny,\
		double Nz,\
		double D,\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz,\
		double &dt,\
		double offang = 0);
	//Utility function - should combine this with above for some speed
	double get_z_intercept(\
		double Nx,\
		double Ny,\
		double Nz,\
		double D,\
		double x,\
		double y,\
		double z,\
		double dx,\
		double dy,\
		double dz);
	//yet another inlinable utility function (make a second with no return?)
	double get_intercept_plane(\
		double Nx,\
		double Ny,\
		double Nz,\
		double D,\
		double &x,\
		double &y,\
		double &z,\
		double dx,\
		double dy,\
		double dz);


	//The compiler should be inlining this without our help
	double sgn(double val);
	

	//inherit and change this function for your implementation
	virtual void warp_readout_box(\
		dirc_point &out_val,\
		int particle_bar,\
		double &mm_index,\
		double &x,\
		double &y,\
		double &z,\
		double &dx,\
		double &dy,\
		double &dz) {}

	//additions to original FastDIRC
	char* geometry_infile   = new char[256];

public:
	//Also inherit and change this - allows for custom quantum efficiency.
	//See source/dirc_threeseg_box_sim.cpp for sample implementation
	virtual double get_cerenkov_angle_rand(double beta, double additional_spread, double &wavelength) {return -1;}
	

	void set_store_bounces(bool isb);
	void fill_bounces_vecs(\
		std::vector<int> &fxbounces,\
		std::vector<int> &fzbounces,\
		std::vector<int> &fxdirbounces,\
		std::vector<int> &fzdirbounces,\
		std::vector<int> &fxindirbounces,\
		std::vector<int> &fzindirbounces);
	void set_kaleidoscope_plot(bool ikp);
	std::vector<double> get_dist_traveled();
	void set_upper_wedge_angle_store(bool istore);
	std::vector<double> get_upper_wedge_incident();
	void set_liquid_index(double li);
	void set_bar_box_angle(double ang);
	void set_bar_box_offsets(double x, double y, double z);
	void set_wedge_mirror_rand(double ispread);
	double get_beta(double E, double m);
	double get_beta_from_p(double p, double m);
	void set_upper_wedge_angle_diff(double rads, double radsy_y = 0);	
	double get_bar_offset(int bar);
	int get_bar_from_x(double x);
	
	void set_use_quartz_n_for_liquid(bool iu);

	void set_moliere_p(double ip);
	void set_use_moliere(bool ium);

        bool convert_particle_kinematics(\
					double &particle_x,\
                                        double &particle_y,\
                                        double &particle_theta,\
                                        double &particle_phi,\
                                        double &particle_bar,\
                                        double particle_x_hall,\
                                        double particle_y_hall,\
                                        double particle_theta_hall,\
                                        double particle_phi_hall);
        bool convert_particle_kinematics_with_offsets(\
					double &particle_x,\
                                        double &particle_y,\
                                        double &particle_theta,\
                                        double &particle_phi,\
                                        double &particle_bar,\
                                        TVector3 particle_x3_hall,\
                                        TVector3 particle_p3_hall);

	void set_offsets_BB_OB_angle(double xang_off, double yang_off, double zang_off);
	void set_offsets_BB_OB_position(double xpos_off, double ypos_off, double zpos_off);
	void set_offsets_upperWedgeFarZ_z_off(double z_off);

	void set_offsets_support_points_off(double x_off, double y_off, double z_off);

	void set_offsets_BB_pos(int BB_number, double xpos_off, double ypos_off, double zpos_off);
	void set_offsets_BB_angle(int BB_number, double xang_off, double yang_off, double zang_off);


	//Random seed chosen arbitrarily
	//default parameters correspond to babar dirc bars
	//default upper wedge top is for gluex implementation.  Set to 0 to remove upper wedge
	DircBaseSim(\
		int rand_seed = 4357,\
		double ibarLength = 4900,\
		double ibarWidth = 35,\
		double ibarDepth = 17,\
		double iupperWedgeTop = 178.6,\
		const char* igeometry_file = "./geometry_files/FastDIRC_HDDS_Nominal.csv");

	std::vector<std::pair<double,double> > get_refraction_rand_phi(\
		std::vector<double> &before_interface,\
		std::vector<double> &after_interface,\
		std::vector<double> &pmt_incidence,\
		int n_photons, \
		double ckov_theta = 47, \
		double particle_x = 0, \
		double particle_y = 0, \
		double particle_theta = 0, \
		double particle_phi = 0,\
		double phi_theta_unc = .0015*57.3,\
		double ckov_theta_unc = .0055*57.3,\
		double beta = -1);
	void sim_rand_n_photons(\
		std::vector<dirc_point> &out_points,\
		int n_photons,\
		double ckov_theta = 47, \
        	double particle_bar= 0, \
		double particle_x = 0, \
		double particle_y = 0, \
		double particle_t = 0, \
		double particle_theta = 0, \
		double particle_phi = 0,\
		double phi_theta_unc = .08594,\
		double ckov_theta_unc = .3151,\
		double beta = -1,\
		int save_kin = -1);	
	void sim_reg_n_photons(\
		std::vector<dirc_point> &out_points,\
		int n_photons_phi,\
		int n_photons_z,\
		double ckov_theta = 47, \
        	double particle_bar= 0, \
		double particle_x = 0, \
		double particle_y = 0, \
		double particle_t = 0, \
		double particle_theta = 0, \
		double particle_phi = 0,\
		double phi_theta_unc = 0,\
		double ckov_theta_unc = 0,\
		double beta = -1,\
		int save_kin = -1);	
	void sim_reg_n_photons_exp(\
		std::vector<dirc_point> &out_points,\
		int n_photons_phi,\
		int n_photons_z,\
		double ckov_theta = 47, \
        	double particle_bar= 0, \
		double particle_x = 0, \
		double particle_y = 0, \
		double particle_t = 0, \
		double particle_theta = 0, \
		double particle_phi = 0,\
		double phi_theta_unc = 0,\
		double ckov_theta_unc = 0,\
		double beta = -1,\
		TH1F* h_dx = NULL,\
		TH2F* h_2D_1 = NULL,\
		TH2F* h_2D_2 = NULL,\
		TH2F* h_2D_3 = NULL,\
		int save_kin = -1);	
	bool track_single_photon(\
	        dirc_point &out_val,\
	        double emit_theta,\
	        double emit_phi,\
       		double particle_theta,\
	        double particle_phi,\
	        double particle_x,\
	        double particle_y,\
	        double particle_z,\
	        double particle_t,\
		int particle_bar);
	bool track_single_photon_beta(\
	        dirc_point &out_val,\
	        double particle_beta,\
	        double emit_phi,\
       		double particle_theta,\
	        double particle_phi,\
	        double particle_x,\
	        double particle_y,\
	        double particle_z,\
	        double particle_t,\
		int particle_bar);
	bool track_line_photon(\
	        dirc_point &out_val,\
	        double particle_beta,\
	        double emit_phi,\
       		double particle_theta,\
	        double particle_phi,\
	        double particle_x,\
	        double particle_y,\
	        double particle_z,\
	        double particle_t,\
		int particle_bar,\
		double z_at_top = 1);
	bool track_all_line_photons(\
                std::vector<dirc_point> &left_vals,\
                std::vector<dirc_point> &right_vals,\
                int points_per_side,\
                double emit_theta,\
                double particle_theta,\
                double particle_phi,\
                double particle_x,\
                double particle_y,\
                double particle_z,\
                double particle_t,\
                int particle_bar,\
                double z_at_top =1);
	void sim_lut_points(\
                std::vector<dirc_point> &ovals,\
                std::vector<double> &phis,\
                std::vector<double> &thetas,\
                int n_photons, \
                double particle_bar /*= 0*/);

	void test_from_wedge_top(\
                std::vector<dirc_point> &ovals,\
                int n_photons, \
                double particle_bar = 1, \
                double particle_x = 0, \
                double phot_theta = 0, \
                double phot_phi = 0,\
                double theta_unc = 0,\
                double phi_unc = 0,\
		double overall_theta = 0); 

	// ------ studies/additions by Y.Y. --------------//
	void sim_rand_n_photons_to_BarEnd(\
		std::vector<dirc_photon> &out_points,\
		int n_photons,\
		double ckov_theta = 47, \
        	double particle_bar= 0, \
		double particle_x = 0, \
		double particle_y = 0, \
		double particle_t = 0, \
		double particle_theta = 0, \
		double particle_phi = 0,\
		double phi_theta_unc = .08594,\
		double ckov_theta_unc = .3151,\
		double beta = -1,\
		int save_kin = -1);	

};
#endif
