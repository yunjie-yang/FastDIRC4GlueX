#include "plot_func.C"
void draw_DLLs_arg(string infilename, string outplot_dir, string plot_label_extra, double res, string hist_type, string run_type)
{
	gStyle->SetOptStat(0);

	string pion_hist,kaon_hist;

	if (run_type == "topdir")
	{
		pion_hist = "ll_diff_pion_"+hist_type;
		kaon_hist = "ll_diff_kaon_"+hist_type;
	}

	std::vector<std::vector<double>> PBins =  {
								{2.9, 3.1},
						  };
	vector<int> rebins = {2000};
	//vector<int> ranges = {85};
	//vector<int> ranges = {60};
	vector<int> ranges = {150};

	//TFile* infile = new TFile(Form("%s/%s.root",infile_dir.c_str(),infile_name.c_str()));
	TFile* infile = new TFile(infilename.c_str());
	string plot_title;
	for (size_t locPBin = 0; locPBin < PBins.size() ; locPBin++)
	{
		plot_title = Form("P: [%.1f, %.1f] GeV, %s",PBins[locPBin][0],PBins[locPBin][1],plot_label_extra.c_str());

		TH1F* hist_ll_diff_pion = (TH1F*) infile -> Get(pion_hist.c_str());
		TH1F* hist_ll_diff_kaon = (TH1F*) infile -> Get(kaon_hist.c_str());

		TCanvas* canv = new TCanvas("canv","canv",600,500);
		//draw_compare_DLL(hist_ll_diff_pion,hist_ll_diff_kaon,rebins[locPBin],ranges[locPBin],plot_title,res,true);
		draw_compare_DLL(hist_ll_diff_pion,hist_ll_diff_kaon,rebins[locPBin],ranges[locPBin],plot_title,res,false);

		canv->Print(Form("%s/dll_%s.pdf",outplot_dir.c_str(),plot_label_extra.c_str()));
		delete canv;

	}




}
