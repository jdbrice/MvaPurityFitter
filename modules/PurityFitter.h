#ifndef TOY_EFF_MAKER_H
#define TOY_EFF_MAKER_H

#include "HistoAnalyzer.h"
#include "XmlFunction.h"
#include "XmlHistogram.h"
#include "CutCollection.h"

#include "RooPlotLib.h"

#include "TRandom3.h"

#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TLatex.h"

#include "vendor/loguru.h"

TH1 * hSignalPDF = nullptr;
TH1 * hBackgroundPDF = nullptr;
TH1 * hKaonPDF = nullptr;
TH1 * hProtonPDF = nullptr;

Double_t fitfun( double *x, double *par ){

	double x0=x[0];
	double scale_sig = par[0];
	double scale_bg = par[1];
	double shift_sig = 0;//par[2];
	double shift_bg = 0;//par[3];

	int sig_bin = hSignalPDF->GetXaxis()->FindBin( x0 );
	int bg_bin = hBackgroundPDF->GetXaxis()->FindBin(x0);

	float sig_val = hSignalPDF->GetBinContent( sig_bin + shift_sig );
	float bg_val = hBackgroundPDF->GetBinContent( bg_bin + shift_bg );

	return sig_val * scale_sig + bg_val * scale_bg;
}

Double_t fitfun3( double *x, double *par ){

	double x0=x[0];
	double scale_sig = par[0];
	double scale_bg = par[1];
	double scale_kaon = par[2];
	

	int sig_bin  = hSignalPDF->GetXaxis()->FindBin( x0 );
	int bg_bin   = hBackgroundPDF->GetXaxis()->FindBin(x0);
	int kaon_bin = hKaonPDF->GetXaxis()->FindBin(x0);

	float sig_val = hSignalPDF->GetBinContent( sig_bin );
	float bg_val = hBackgroundPDF->GetBinContent( bg_bin );
	float kaon_val = hKaonPDF->GetBinContent( kaon_bin );

	return sig_val * scale_sig + bg_val * scale_bg + kaon_val * scale_kaon;
}

Double_t fitfun4( double *x, double *par ){

	double x0=x[0];
	double scale_sig = par[0];
	double scale_bg = par[1];
	double scale_kaon = par[2];
	double scale_proton = par[3];
	

	int sig_bin  = hSignalPDF->GetXaxis()->FindBin( x0 );
	int bg_bin   = hBackgroundPDF->GetXaxis()->FindBin(x0);
	int kaon_bin = hKaonPDF->GetXaxis()->FindBin(x0);
	int proton_bin = hProtonPDF->GetXaxis()->FindBin(x0);

	float sig_val = hSignalPDF->GetBinContent( sig_bin );
	float bg_val = hBackgroundPDF->GetBinContent( bg_bin );
	float kaon_val = hKaonPDF->GetBinContent( kaon_bin );
	float proton_val = hProtonPDF->GetBinContent( proton_bin );

	return sig_val * scale_sig + bg_val * scale_bg + kaon_val * scale_kaon + proton_val * scale_proton;
}




class PurityFitter : public HistoAnalyzer {
protected:

	HistoBins pid_bins;

	string var;
	string rpName;

	vector<string> var_names;
	map<string, HistoBins> bins;

	bool allow_proton;
	bool allow_kaon;

	bool export_img;

public:

	virtual void initialize(){
		HistoAnalyzer::initialize();

		
		var = config.getString( "var", "mlp" );
		pid_bins.load( config, var );

		var_names = config.getStringVector( "vars" );
		for ( string v : var_names ){
			bins[ v ].load( config, v );
		}

		bins[ "mlp" ].load( config, "mlp" );

		allow_kaon   = config.getBool( "fit:kaon", false );
		allow_proton = config.getBool( "fit:proton", false );

		export_img = config.getBool( "can:export", false );
	}


	TH1 * pt_projection( string name, TH2 * hin, float pt1, float pt2, HistoBins &rb ){

		TAxis *x = hin->GetXaxis();
		int ipt1 = x->FindBin( pt1 );
		int ipt2 = x->FindBin( pt2 );

		TH1 *hout =  hin->ProjectionY( name.c_str(), ipt1, ipt2 );
		if ( rb.nBins() > 2 )
			hout = hout->Rebin( rb.nBins(), (name + "_rb").c_str(), rb.bins.data() );
		return hout;
	}

	TH1 * subsample( TH1 * hin, float factor = 1.0 ){
		TH1 * hinsmooth = (TH1*)hin->Clone( TString::Format( "%s_smooth", hin->GetName() ) );
		float bw = hin->GetXaxis()->GetBinWidth(10);
		hinsmooth->Reset();
		TRandom3 r;
		r.SetSeed(0);
		size_t nSamples = 1000000;
		for ( int i = 0; i < nSamples; i ++ ){
			float v = hin->GetRandom();
			hinsmooth->Fill( v + r.Gaus( 0, bw * factor ) );
		}

		return hinsmooth;
	}

	void compare_other( TH1 * hsigFit, TH1 * hbgFit, TH1 * hkaonFit, TH1 * hprotonFit, float pt1, float pt2, TCanvas *can ){

		RooPlotLib rpl;

		bool use_kaon = false;
		if ( hkaonFit != nullptr )
			use_kaon = true;

		bool use_proton = false;
		if ( hprotonFit != nullptr )
			use_proton = true;

		can->Clear();
		if ( false == export_img )
			can->Divide( 2, 3, 0.001, 0.001 );
		float lw = 1;
		if ( export_img )
			lw = 2;
		int i = 1; 
		for ( string v : var_names ){
			int ix = i % 2;
			int iy = i / 3;
			can->cd(i);
			i++;
			TH2 * hdata_pos = get<TH2>( "pos_" + v + "_vs_pt", "data" );
			TH2 * hdata_neg = get<TH2>( "neg_" + v + "_vs_pt", "data" );

			TH2 * hmc_pos_sig = get<TH2>( "sig_pos_" + v, "mc" );
			TH2 * hmc_neg_sig = get<TH2>( "sig_neg_" + v, "mc" );

			TH2 * hmc_pos_bg = get<TH2>( "bg_pos_" + v, "mc" );
			TH2 * hmc_neg_bg = get<TH2>( "bg_neg_" + v, "mc" );

			TH2 * hmc_pos_kaon = get<TH2>( "kaon_pos_" + v, "kaon" );
			TH2 * hmc_neg_kaon = get<TH2>( "kaon_neg_" + v, "kaon" );

			TH2 * hmc_pos_proton = get<TH2>( "proton_pos_" + v, "proton" );
			TH2 * hmc_neg_proton = get<TH2>( "proton_neg_" + v, "proton" );

			TH2 * hdata_use = hdata_pos;
			TH2 * hsig_use = hmc_pos_sig;
			TH2 * hbg_use = hmc_pos_bg;
			TH2 * hkaon_use = hmc_pos_kaon;
			TH2 * hproton_use = hmc_pos_proton;
			if ( "neg" == config["charge"] ){
				hdata_use = hdata_neg;
				hsig_use  = hmc_neg_sig;
				hbg_use   = hmc_neg_bg;
				hkaon_use = hmc_neg_kaon;
				hproton_use = hmc_neg_proton;
			}

			string suffix="_" + v + "_"+dtes( pt1 ) +"_" + dtes(pt2);
			TH1 * hdata    = pt_projection( "data" + suffix, hdata_use, pt1, pt2, bins[ v ]  );
			TH1 * hsig     = pt_projection( "sig" + suffix, hsig_use, pt1, pt2, bins[ v ]  );
			TH1 * hbg      = pt_projection( "bg" + suffix, hbg_use, pt1, pt2, bins[ v ]  );
			TH1 * hkaon    = pt_projection( "kaon" + suffix, hkaon_use, pt1, pt2, bins[ v ]  );
			TH1 * hproton  = pt_projection( "proton" + suffix, hproton_use, pt1, pt2, bins[ v ]  );

			// hbg = subsample( hbg );
			// hkaon = subsample( hkaon );
			// hsig = subsample( hsig );

			hdata->Sumw2(true);
			hdata->Scale( 1.0 / hdata->Integral() );
			hsig->Scale( 1.0 / hsig->Integral() );
			hbg->Scale( 1.0 / hbg->Integral() );
			if ( use_kaon ) hkaon->Scale( 1.0 / hkaon->Integral() );
			if ( use_proton ) hproton->Scale( 1.0 / hproton->Integral() );

			rpl.style( hdata ).set( "xtitle", v + " " + config[ "units:" + v ] )
				.set( "title", TString::Format("%0.2f < pT < %0.2f (GeV/c)", pt1, pt2 ).Data() )
				.set( "ytitle", "dN/d" + v + " " + config[ "units" ] + "^{-1}" )
				.set( "lw", "1.5" )
				.set( "fc", "#999" )
				.set( "yto", 1.2 )
				.set( "logy", 1 )
				.set( "draw", "h" )
				// .set( "xr", -0.2, 1.1 )
				.set( "min", 1 )
				.set( "lc", "#222222" ).draw();

			hsig->Scale( hsigFit->Integral() / hsig->Integral() );
			hbg->Scale( hbgFit->Integral() / hbg->Integral() );
			if ( use_kaon ) hkaon->Scale( hkaonFit->Integral() / hkaon->Integral() );
			if ( use_proton ) hproton->Scale( hprotonFit->Integral() / hproton->Integral() );

			TH1 * hsum = (TH1*)hbg->Clone( ("hsum"+suffix).c_str() );
			hsum->Reset();

			hsum->Add( hsig );
			hsum->Add( hbg );
			if ( use_kaon )  hsum->Add( hkaon );
			if ( use_proton )  hsum->Add( hproton );

			rpl.style( hsig ).set( config, "style.sig" ).set( "lw", lw ).draw();
			rpl.style( hbg ).set( config, "style.bg" ).set( "lw", lw ).draw();
			if ( use_kaon ) rpl.style( hkaon ).set( config, "style.kaon" ).set( "lw", lw ).draw();
			if ( use_proton ) rpl.style( hproton ).set( config, "style.proton" ).set( "lw", lw ).draw();
			rpl.style( hsum ).set( config, "style.sum" ).set( "lw", lw ).draw();

			if ( export_img )
				can->Print( ("export/"+ config["mod"] +"/compare-" + v + "-" + dts( pt1 ) + "-" + dts(pt2) + ".pdf" ).c_str() );
		} // loop on vars
		can->Print( (rpName).c_str()  );
	}



	float calc_purity( TH1 * hsig, TH1 * hbg, TH1 * hkaonFit, TH1 * hprotonFit, float pt1, float pt2, float cut = 0.95 ){


		int ibsig = hsig->GetXaxis()->FindBin( cut );
		int ibbg = hbg->GetXaxis()->FindBin( cut );
		
		int ibsig2 = hsig->GetXaxis()->FindBin( 1.2 );
		int ibbg2 = hbg->GetXaxis()->FindBin( 1.2 );
		

		int ibkaonFit = -1;
		if ( nullptr != hkaonFit ) hkaonFit->GetXaxis()->FindBin( cut );
		int ibkaonFit2 = -1;
		if ( nullptr != hkaonFit ) hkaonFit->GetXaxis()->FindBin( 1.2 );

		float eff = hsig->Integral( ibsig, ibsig2 ) / hsig->Integral();
		float purity = 0;

		if ( nullptr != hkaonFit ) 
			purity = hsig->Integral( ibsig, ibsig2 ) / ( hbg->Integral( ibbg, ibbg2 ) + hkaonFit->Integral( ibkaonFit, ibkaonFit2 ) + hsig->Integral( ibsig, ibsig2 ) );
		else 
			purity = hsig->Integral( ibsig, ibsig2 ) / ( hbg->Integral( ibbg, ibbg2 ) + hsig->Integral( ibsig, ibsig2 ) );

		int ipt = book->get( "purity" )->GetXaxis()->FindBin( (pt1+pt2)/2.0 );
		book->get( "purity" )->SetBinContent( ipt, purity );
		book->get( "efficiency" )->SetBinContent( ipt, eff );
	}

	virtual void make(){
		LOG_SCOPE_FUNCTION( INFO );

		book->cd();
		book->makeAll( config, nodePath + ".histograms" );


		TH2 * hdata_pos = get<TH2>( "pos_" + var + "_vs_pt", "data" );
		TH2 * hdata_neg = get<TH2>( "neg_" + var + "_vs_pt", "data" );

		TH2 * hmc_pos_sig = get<TH2>( "sig_pos_" + var, "mc" );
		TH2 * hmc_neg_sig = get<TH2>( "sig_neg_" + var, "mc" );

		TH2 * hmc_pos_bg = get<TH2>( "bg_pos_" + var, "mc" );
		TH2 * hmc_neg_bg = get<TH2>( "bg_neg_" + var, "mc" );

		TH2 * hmc_pos_kaon = get<TH2>( "kaon_pos_" + var, "kaon" );
		TH2 * hmc_neg_kaon = get<TH2>( "kaon_neg_" + var, "kaon" );

		TH2 * hmc_pos_proton = get<TH2>( "proton_pos_" + var, "proton" );
		TH2 * hmc_neg_proton = get<TH2>( "proton_neg_" + var, "proton" );
		

		LOG_F( INFO, "hdata_pos=%p", hdata_pos );
		LOG_F( INFO, "hdata_neg=%p", hdata_neg );

		LOG_F( INFO, "hmc_pos_sig=%p", hmc_pos_sig );
		LOG_F( INFO, "hmc_neg_sig=%p", hmc_neg_sig );

		LOG_F( INFO, "hmc_pos_bg=%p", hmc_pos_bg );
		LOG_F( INFO, "hmc_neg_bg=%p", hmc_neg_bg );

		LOG_F( INFO, "hmc_pos_kaon=%p", hmc_pos_kaon );
		LOG_F( INFO, "hmc_neg_kaon=%p", hmc_neg_kaon );

		LOG_F( INFO, "hmc_pos_proton=%p", hmc_pos_proton );
		LOG_F( INFO, "hmc_neg_proton=%p", hmc_neg_proton );

		TCanvas *can = new TCanvas( "can", "can", config.get<int>( "can:w", 500 ), config.get<int>( "can:h", 500 ) );
		TCanvas *can2 = new TCanvas( "can2", "can2", config.get<int>( "can:w", 500 ), config.get<int>( "can:h2", 500 ) );
		can->SetTopMargin( 0.05 );
		can->SetRightMargin( 0.01 );
		can->SetBottomMargin( 0.11 );
		can->SetLeftMargin( 0.15 );

		can2->SetTopMargin( 0.05 );
		can2->SetRightMargin( 0.01 );
		can2->SetBottomMargin( 0.11 );
		can2->SetLeftMargin( 0.15 );


		rpName = config[nodePath + ".output.Report:url" ];
		can->Print( (rpName+"[").c_str() );
		
		gStyle->SetOptFit(111);
		gStyle->SetOptStat(0);

		RooPlotLib rpl;

		HistoBins pt_bins;
		pt_bins.load( config, "pt" );

		TH2 * hdata_use = hdata_pos;
		TH2 * hsig_use = hmc_pos_sig;
		TH2 * hbg_use = hmc_pos_bg;
		TH2 * hkaon_use = hmc_pos_kaon;
		TH2 * hproton_use = hmc_pos_proton;
		if ( "neg" == config["charge"] ){
			hdata_use = hdata_neg;
			hsig_use  = hmc_neg_sig;
			hbg_use   = hmc_neg_bg;
			hkaon_use = hmc_neg_kaon;
			hproton_use = hmc_neg_proton;
		}

		for ( int i = 0; i < pt_bins.nBins(); i++ ){
			float pt1 = pt_bins[i];
			float pt2 = pt_bins[(i+1)];
			LOG_F( INFO, "pt=(%f, %f)", pt1, pt2 );

			can->Clear();

			string suffix="_"+dtes( pt1 ) +"_" + dtes(pt2);
			TH1 * hdata    = pt_projection( "data" + suffix, hdata_use, pt1, pt2, bins[ var ] );
			TH1 * hsig     = pt_projection( "sig" + suffix, hsig_use, pt1, pt2, bins[ var ] );
			TH1 * hbg      = pt_projection( "bg" + suffix, hbg_use, pt1, pt2, bins[ var ] );
			TH1 * hkaon    = pt_projection( "kaon" + suffix, hkaon_use, pt1, pt2, bins[ var ] );
			TH1 * hproton  = pt_projection( "proton" + suffix, hproton_use, pt1, pt2, bins[ var ] );

			// TH1 * hbgsmooth = (TH1*)hbg->Clone( "hbgsmooth" );
			// float bw = hbg->GetXaxis()->GetBinWidth(10);
			// hbgsmooth->Reset();
			// TRandom3 r;
			// r.SetSeed(0);
			// size_t nSamples = 100000;
			// for ( int i = 0; i < nSamples; i ++ ){
			// 	float v = hbg->GetRandom();
			// 	hbgsmooth->Fill( v + r.Gaus( 0, bw ) );
			// }

			// hbg = hbgsmooth;

			hbg = subsample( hbg );
			// hproton = subsample( hproton, 1.5 );
			hkaon = subsample( hkaon );
			// hsig = subsample( hsig );

			
			LOG_F( INFO, "hdata   = %p", hdata );
			LOG_F( INFO, "hsig    = %p", hsig );
			LOG_F( INFO, "hbg     = %p", hbg );
			LOG_F( INFO, "hkaon   = %p", hkaon );
			LOG_F( INFO, "hproton = %p", hproton );

			hdata->Sumw2(true);
			float Idata = hdata->Integral();
			hdata->Scale( 1.0 / Idata );
			
			

			hsig->Scale( 1.0 / hsig->Integral() );
			hbg->Scale( 1.0 / hbg->Integral() );
			bool use_kaon = false;
			if ( pt1 < 3.0 && allow_kaon)
				use_kaon = true;

			bool use_proton = false;
			if ( pt1 < 2.0 && allow_proton)
				use_proton = true;
			
			if ( use_kaon )
				hkaon->Scale( 1.0 / hkaon->Integral() );
			if ( use_proton )
				hproton->Scale( 1.0 / hproton->Integral() );


			hSignalPDF     = (TH1*)hsig->Clone( "hSignalPDF" );
			hBackgroundPDF = (TH1*)hbg->Clone( "hBackgroundPDF" );
			if ( use_kaon )
				hKaonPDF = (TH1*)hkaon->Clone( "hKaonPDF" );
			else 
				hKaonPDF = nullptr;
			if ( use_proton )
				hProtonPDF = (TH1*)hproton->Clone( "hProtonPDF" );
			else 
				hProtonPDF = nullptr;

			
			TH1 * hsum = (TH1*)hbg->Clone( ("hsum"+suffix).c_str() );
			hsum->Reset();

			TF1 * ff = nullptr;
			if ( false == use_kaon && false == use_proton){
				ff = new TF1( "ff", fitfun, -1, 2, 2 );
				ff->SetParNames( "signal", "bg" );
				ff->SetParLimits( 0, 1e-4, 1e9 );
				ff->SetParLimits( 1, 1e-4, 1e9 );
				ff->SetParameters( 1.0, 100.0 );
			} else if ( true == use_kaon && false == use_proton ){
				ff = new TF1( "ff", fitfun3, -1, 2, 3 );
				ff->SetParNames( "signal", "bg", "kaon" );
				ff->SetParLimits( 0, 1e-4, 1e9 );
				ff->SetParLimits( 1, 1e-4, 1e9 );
				ff->SetParLimits( 2, 1e-4, 1e9 );
				ff->SetParameters( 1.0, 0.1, 0.01 );
				// ff->FixParameter( 1, 0.00001 );
			} else if ( false == use_kaon && true == use_proton ){
				ff = new TF1( "ff", fitfun3, -1, 2, 3 );
				ff->SetParNames( "signal", "bg", "proton" );
				// ff->SetParLimits( 0, 1e-4, 1e9 );
				// ff->SetParLimits( 1, 1e-4, 1e9 );
				// ff->SetParLimits( 2, 1e-4, 1e9 );
				ff->SetParameters( 1.0, 100.0, 50 );
			} else if ( use_kaon && use_proton ){
				ff = new TF1( "ff", fitfun4, -1, 2, 4 );
				ff->SetParNames( "signal", "bg", "kaon", "proton" );
				ff->SetParLimits( 0, 1e-4, 1e9 );
				ff->SetParLimits( 1, 1e-4, 1e9 );
				ff->SetParLimits( 2, 1e-4, 1e9 );
				ff->SetParLimits( 3, 1e-4, 1e9 );
				ff->SetParameters( 1.0, 100.0, 500.0, 1e-3 );
			}
			
			ff->SetNpx( 1000 );
			ff->SetLineColor(kBlack);
			ff->SetLineWidth(2);

			// hdata->Draw();
	
			rpl.style( hdata ).set( "xtitle", var + " " + config[ "units" ] )
				.set( "ytitle", "dN/d" + var + " " + config[ "units" ] + "^{-1}" )
				.set( "title", TString::Format("%s : %0.2f < pT < %0.2f (GeV/c)", config["charge"].c_str(), pt1, pt2 ).Data() )
				.set( "lw", "1.5" )
				.set( "fc", "#999" )
				.set( "yto", 1.3 )
				.set( "logy", 1 )
				.set( "draw", "h" )
				// .set( "xr", -0.2, 1.1 )
				// .set( "min", 1 )
				.set( "max", 1.0 )
				.set( "lc", "#222222" ).draw();

			string fitOpt = config[ "fit:opt" ];
			if ( config.exists( "fit:min" ) && config.exists( "fit:max" ) ){
				LOG_F( INFO, "Fitting (opt=%s) in range (%f, %f)", fitOpt.c_str(), config.get<float>( "fit:min" ), config.get<float>( "fit:max" ) );

				hdata->Fit( "ff", fitOpt.c_str(), "", config.get<float>( "fit:min" ), config.get<float>( "fit:max" ) );
				hdata->Fit( "ff", fitOpt.c_str(), "", config.get<float>( "fit:min" ), config.get<float>( "fit:max" ) );
				hdata->Fit( "ff", fitOpt.c_str(), "", config.get<float>( "fit:min" ), config.get<float>( "fit:max" ) );

			} else {
				hdata->Fit( "ff", fitOpt.c_str() );
				hdata->Fit( "ff", fitOpt.c_str() );
				hdata->Fit( "ff", fitOpt.c_str() );
			}

			

			hsig->Scale( ff->GetParameter( "signal" ) );
			hbg->Scale( ff->GetParameter( "bg" ) );
			if ( use_kaon )hkaon->Scale( ff->GetParameter( "kaon" ) );
			if ( use_proton ) hproton->Scale( ff->GetParameter( "proton" ) );

			hsum->Add( hsig );
			hsum->Add( hbg );
			if ( use_kaon )
				hsum->Add( hkaon );
			if ( use_proton )
				hsum->Add( hproton );

			rpl.style( hsig ).set( config, "style.sig" ).draw();
			rpl.style( hbg ).set( config, "style.bg" ).draw();
			if ( use_kaon ) rpl.style( hkaon ).set( config, "style.kaon" ).draw();
			if ( use_proton ) rpl.style( hproton ).set( config, "style.proton" ).draw();
			rpl.style( hsum ).set( config, "style.sum" ).draw();

			TLatex tl;
			tl.SetTextSize( 12.0 / 360.0 );
			tl.DrawLatexNDC( 0.5, 0.85, TString::Format("#chi^2 / ndf = %0.2f / %d = %0.2f", ff->GetChisquare(), ff->GetNDF(), ff->GetChisquare() / (float)ff->GetNDF() ) );
			tl.DrawLatexNDC( 0.5, 0.8, TString::Format("signal = %0.2f", ff->GetParameter( "signal" )) );
			tl.DrawLatexNDC( 0.5, 0.75, TString::Format("background (#pi) = %0.2f", ff->GetParameter( "bg" )) );
			if ( use_kaon ) tl.DrawLatexNDC( 0.5, 0.70, TString::Format("background (K) = %0.6f", ff->GetParameter( "kaon" )) );
			
			if ( use_proton && !use_kaon ) tl.DrawLatexNDC( 0.5, 0.70, TString::Format("background (P) = %0.6f", ff->GetParameter( "proton" )) );
			if ( use_proton && use_kaon ) tl.DrawLatexNDC( 0.5, 0.65, TString::Format("background (P) = %0.6f", ff->GetParameter( "proton" )) );
			

			int iptbin = book->get( "chi2ndf" )->GetXaxis()->FindBin( (pt1+pt2)/2.0 );
			book->get( "chi2ndf" )->SetBinContent( iptbin, ff->GetChisquare() / (float)ff->GetNDF() );

			TLegend *leg = new TLegend( 0.15, 0.8, 0.45, 0.95 );
			leg->SetNColumns(2);
			leg->AddEntry( hsig, "Signal" );
			leg->AddEntry( hbg, "#pi Background" );
			if ( use_kaon ) leg->AddEntry( hkaon, "K Background" );
			if ( use_proton ) leg->AddEntry( hproton, "p Background" );

			leg->Draw("same");

			can->Print( (rpName).c_str()  );
			if ( export_img ) can->Print( ("export/"+ config["mod"] + "/fit-" + dts( pt1 ) + "-" + dts(pt2) + ".pdf" ).c_str() );


			hsum->Scale( Idata );
			hdata->Scale( Idata );
			hsum->Divide( hdata );
			hsum->SetTitle( TString::Format("%0.2f < pT < %0.2f; %s; fit / data", pt1, pt2, var.c_str() ) );
			// gPad->SetLogy(0);
			
			if ( config.exists( "fit:min" ) && config.exists( "fit:max" ) ){
				hsum->Fit( "pol0", "R", "", config.get<float>( "fit:min" ), config.get<float>( "fit:max" ) );
			} else {
				hsum->Fit( "pol0" );
			}
			rpl.style( hsum )
				// .set( "xr", -0.2, 1.1 )
				.set( "logy", 0 )
				.set( "min", 0 )
				.set( "max", 2.0 )
				.set( "yto", 1.2 );
			hsum->Draw();

			can->Print( (rpName).c_str()  );
			if ( export_img ) {
				can2->cd();
				hsum->Draw();
				can2->Print( ("export/"+ config["mod"] +"/ratio-" + dts( pt1 ) + "-" + dts(pt2) + ".pdf" ).c_str() );
				can->cd();
			}


			if ( !use_kaon )
				hkaon = nullptr;
			if ( !use_proton )
				hproton = nullptr;
			
			compare_other( hsig, hbg, hkaon, hproton, pt1, pt2, can );
			calc_purity( hsig, hbg, hkaon, nullptr, pt1, pt2 );
			


			
			delete ff;
			delete hsum;
			delete hsig;
			delete hbg;
			delete hSignalPDF;
			delete hBackgroundPDF;

		}




		if ( 0 == config.getInt( "jobIndex" ) || -1 == config.getInt( "jobIndex" ) ){
			TNamed config_str( "config", config.toXml() );
			config_str.Write();
		}


		can->Print( (rpName+"]").c_str()  );

	} // make

};



#endif