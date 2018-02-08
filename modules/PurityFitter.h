#ifndef TOY_EFF_MAKER_H
#define TOY_EFF_MAKER_H

#include "HistoAnalyzer.h"
#include "XmlFunction.h"
#include "XmlHistogram.h"
#include "CutCollection.h"

#include "RooPlotLib.h"


#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TLatex.h"

#include "vendor/loguru.h"

TH1 * hSignalPDF = nullptr;
TH1 * hBackgroundPDF = nullptr;

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




class PurityFitter : public HistoAnalyzer {
protected:

	HistoBins pid_bins;

	string var;
	string rpName;

	vector<string> var_names;
	map<string, HistoBins> bins;

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

	void compare_other( TH1 * hsigFit, TH1 * hbgFit, float pt1, float pt2, TCanvas *can ){

		RooPlotLib rpl;

		can->Clear();
		can->Divide( 2, 3, 0.001, 0.001 );
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

			TH2 * hdata_use = hdata_pos;
			TH2 * hsig_use = hmc_pos_sig;
			TH2 * hbg_use = hmc_pos_bg;
			if ( "neg" == config["charge"] ){
				hdata_use = hdata_neg;
				hsig_use  = hmc_neg_sig;
				hbg_use   = hmc_neg_bg;
			}

			string suffix="_" + v + "_"+dtes( pt1 ) +"_" + dtes(pt2);
			TH1 * hdata    = pt_projection( "data" + suffix, hdata_use, pt1, pt2, bins[ v ]  );
			TH1 * hsig     = pt_projection( "sig" + suffix, hsig_use, pt1, pt2, bins[ v ]  );
			TH1 * hbg      = pt_projection( "bg" + suffix, hbg_use, pt1, pt2, bins[ v ]  );

			hdata->Sumw2(true);
			hdata->Scale( 1.0 / hdata->Integral() );
			hsig->Scale( 1.0 / hsig->Integral() );
			hbg->Scale( 1.0 / hbg->Integral() );

			rpl.style( hdata ).set( "xtitle", v + " " + config[ "units:" + v ] )
				.set( "title", TString::Format("%0.2f < pT < %0.2f", pt1, pt2 ).Data() )
				.set( "ytitle", "dN/d" + v + " " + config[ "units" ] + "^{-1}" )
				.set( "lw", "1.5" )
				.set( "fc", "#999" )
				.set( "yto", 1.5 )
				.set( "logy", 1 )
				.set( "draw", "h" )
				// .set( "xr", -0.2, 1.1 )
				.set( "min", 1 )
				.set( "lc", "#222222" ).draw();

			hsig->Scale( hsigFit->Integral() / hsig->Integral() );
			hbg->Scale( hbgFit->Integral() / hbg->Integral() );

			TH1 * hsum = (TH1*)hbg->Clone( ("hsum"+suffix).c_str() );
			hsum->Reset();

			hsum->Add( hsig );
			hsum->Add( hbg );

			rpl.style( hsig ).set( config, "style.sig" ).set( "lw", 1 ).draw();
			rpl.style( hbg ).set( config, "style.bg" ).set( "lw", 1 ).draw();
			rpl.style( hsum ).set( config, "style.sum" ).set( "lw", 1 ).draw();
		} // loop on vars
		can->Print( (rpName).c_str()  );
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

		

		LOG_F( INFO, "hdata_pos=%p", hdata_pos );
		LOG_F( INFO, "hdata_neg=%p", hdata_neg );

		LOG_F( INFO, "hmc_pos_sig=%p", hmc_pos_sig );
		LOG_F( INFO, "hmc_neg_sig=%p", hmc_neg_sig );

		LOG_F( INFO, "hmc_pos_bg=%p", hmc_pos_bg );
		LOG_F( INFO, "hmc_neg_bg=%p", hmc_neg_bg );

		TCanvas *can = new TCanvas( "can", "can", 500, 500 );
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
		if ( "neg" == config["charge"] ){
			hdata_use = hdata_neg;
			hsig_use  = hmc_neg_sig;
			hbg_use   = hmc_neg_bg;
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

			hdata->Sumw2(true);
			float Idata = hdata->Integral();
			hdata->Scale( 1.0 / Idata );
			
			

			hsig->Scale( 1.0 / hsig->Integral() );
			hbg->Scale( 1.0 / hbg->Integral() );

			hSignalPDF     = (TH1*)hsig->Clone( "hSignalPDF" );
			hBackgroundPDF = (TH1*)hbg->Clone( "hBackgroundPDF" );

			
			TH1 * hsum = (TH1*)hbg->Clone( ("hsum"+suffix).c_str() );
			hsum->Reset();

			TF1 * ff = new TF1( "ff", fitfun, -1, 1, 2 );
			ff->SetParNames( "signal", "bg" );
			ff->SetParLimits( 0, 1e-4, 1e9 );
			ff->SetParLimits( 1, 1e-4, 1e9 );
			ff->SetParameters( 1.0, 100.0 );
			ff->SetNpx( 1000 );
			ff->SetLineColor(kBlack);
			ff->SetLineWidth(2);

			// hdata->Draw();
	
			rpl.style( hdata ).set( "xtitle", var + " " + config[ "units" ] )
				.set( "ytitle", "dN/d" + var + " " + config[ "units" ] + "^{-1}" )
				.set( "title", TString::Format("%s : %0.2f < pT < %0.2f", config["charge"].c_str(), pt1, pt2 ).Data() )
				.set( "lw", "1.5" )
				.set( "fc", "#999" )
				.set( "yto", 1.3 )
				.set( "logy", 1 )
				.set( "draw", "h" )
				.set( "xr", -0.2, 1.1 )
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

			

			hsig->Scale( ff->GetParameter(0) );
			hbg->Scale( ff->GetParameter(1) );

			hsum->Add( hsig );
			hsum->Add( hbg );

			rpl.style( hsig ).set( config, "style.sig" ).draw();
			rpl.style( hbg ).set( config, "style.bg" ).draw();
			rpl.style( hsum ).set( config, "style.sum" ).draw();

			TLatex tl;
			tl.SetTextSize( 12.0 / 360.0 );
			tl.DrawLatexNDC( 0.4, 0.85, TString::Format("#chi^2 / ndf = %0.2f / %d = %0.2f", ff->GetChisquare(), ff->GetNDF(), ff->GetChisquare() / (float)ff->GetNDF() ) );
			tl.DrawLatexNDC( 0.4, 0.8, TString::Format("signal = %0.2f", ff->GetParameter(0)) );
			tl.DrawLatexNDC( 0.4, 0.75, TString::Format("background (#pi) = %0.2f", ff->GetParameter(1)) );
			

			can->Print( (rpName).c_str()  );


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
			rpl.style( hsum ).set( "logy", 0 ).set( "min", 0 ).set( "max", 2.0 ).set( "yto", 1.2 );
			hsum->Draw();

			can->Print( (rpName).c_str()  );

			compare_other( hsig, hbg, pt1, pt2, can );


			
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