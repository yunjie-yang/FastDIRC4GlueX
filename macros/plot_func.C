void draw_compare_DLL(TH1* hist_dll_pion, TH1* hist_dll_kaon, int rebin, double range, string plot_title, double res = -1, bool doFit=false, string extra_label = "") 
{

        hist_dll_pion -> SetLineColor(kBlue);
        hist_dll_kaon -> SetLineColor(kRed);

        hist_dll_pion -> SetMarkerColor(kBlue);
        hist_dll_kaon -> SetMarkerColor(kRed);

        hist_dll_pion -> SetMarkerColor(kBlue);
        hist_dll_kaon -> SetMarkerColor(kRed);

        hist_dll_pion -> SetFillColorAlpha(kBlue,0.5);
        hist_dll_kaon -> SetFillColorAlpha(kRed,0.5);

        hist_dll_pion -> Rebin(rebin);
        hist_dll_kaon -> Rebin(rebin);

        hist_dll_pion -> Scale(1./hist_dll_pion->GetMaximum());
        hist_dll_kaon -> Scale(1./hist_dll_kaon->GetMaximum());

        //hist_dll_pion -> SetTitle(plot_title.c_str());
        hist_dll_pion -> SetTitle("");
        hist_dll_pion -> GetXaxis() -> SetTitle("#Deltalog(L)");
        hist_dll_pion -> GetXaxis() -> SetTitleSize(0.06);
        hist_dll_pion -> GetXaxis() -> SetTitleOffset(0.7);
        hist_dll_pion -> GetYaxis() -> SetTitleSize(0.06);
        hist_dll_pion -> GetYaxis() -> SetTitleOffset(0.45);
        hist_dll_pion -> GetYaxis() -> SetTitle("a.u.");
        hist_dll_pion -> GetXaxis() -> SetRangeUser(-range,range);


	TF1 *ff_pion, *ff_kaon;
	double sep=0,m1=0,m2=0,s1=0,s2=0;

	//if(hist_dll_pion->GetEntries()>100){
	if(doFit)
	{
		hist_dll_pion->Fit("gaus","S");
    		ff_pion = hist_dll_pion->GetFunction("gaus");
    		ff_pion->SetLineColor(1);
		m1=ff_pion->GetParameter(1);
		s1=ff_pion->GetParameter(2);
	}
	//if(hist_dll_kaon->GetEntries()>100){
	if(doFit)
	{
		hist_dll_kaon->Fit("gaus","S");
		ff_kaon = hist_dll_kaon->GetFunction("gaus");
		ff_kaon->SetLineColor(1);
		m2=ff_kaon->GetParameter(1);
		s2=ff_kaon->GetParameter(2);
	}

	if(doFit)
	{
		if(s1>0 && s2>0) sep = (fabs(m2-m1))/(0.5*(s1+s2));
		printf("m1  = %6.2f, s1 = %6.2f\n", m1, s1);
		printf("m2  = %6.2f, s2 = %6.2f\n", m2, s2);
		printf("sep = %6.2f\n", sep);
	}

        //TLegend* leg = new TLegend(0.725,0.675,0.9,0.9);
        TLegend* leg = new TLegend(0.85,0.675,0.975,0.9);
	leg->SetFillColor(0);
	leg->SetFillStyle(0);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);

	char buff[256];
	//sprintf(buff,"#pi (%d)",int(hist_dll_pion->GetEntries()));
	sprintf(buff," #pi");
        leg->AddEntry(hist_dll_pion,buff,"f");
	//sprintf(buff,"K (%d)",int(hist_dll_kaon->GetEntries()));
	sprintf(buff," K");
        leg->AddEntry(hist_dll_kaon,buff,"f");
	if (doFit)
		leg->AddEntry((TObject*)0, Form("%2.2f #sigma",sep), "");

        TLegend* leg2 = new TLegend(0.05,0.7,0.25,0.85);
	leg2->SetFillColor(0);
	leg2->SetFillStyle(0);
	leg2->SetBorderSize(0);
	leg2->SetFillStyle(0);
	leg2->SetTextAlign(13);
	leg2->SetTextSize(0.065);
	leg2->AddEntry((TObject*)0, Form("#sigma_{#theta_{C}} = %2.1f mrad",res), "");

	if (extra_label!="")
		hist_dll_pion->SetTitle(Form("%s",extra_label.c_str()));

        hist_dll_pion->Draw("HIST");
        hist_dll_kaon->Draw("HIST SAME");

	if (doFit)
	{
		ff_pion->Draw("SAME");
		ff_kaon->Draw("SAME");
	}
        leg->Draw();
	if (res > 0.)
        	leg2->Draw();


	gPad->SetTopMargin(0.08);
	gPad->SetLeftMargin(0.07);
	gPad->SetRightMargin(0.03);
        gPad->RedrawAxis();
}
