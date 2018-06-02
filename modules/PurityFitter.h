#ifndef TOY_EFF_MAKER_H
#define TOY_EFF_MAKER_H

#include "HistoAnalyzer.h"
#include "XmlFunction.h"
#include "XmlHistogram.h"
#include "CutCollection.h"

#include "RooPlotLib.h"
#include "FitConfidence.h"

#include "TRandom3.h"

#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TLatex.h"

#include "vendor/loguru.h"
#include "TGraphAsymmErrors.h"

#include "TString.h"
// using TString;

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
	double scale_sig  = par[0];
	double scale_bg   = par[1];
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

	bool allow_proton, use_proton;
	bool allow_kaon, use_kaon;

	bool export_img;

	TCanvas *can = nullptr, *can2 = nullptr;

	float Idata = 0;

	map<string, vector<TH1 *>> nn_pt;

	TH2 * bgTemplate, *sigTemplate;

	size_t nSamples = 1000;
public:

	virtual void initialize(){
		HistoAnalyzer::initialize();

		
		var = config.getString( "var", "mlp" );
		pid_bins.load( config, var );

		var_names = config.getStringVector( "vars" );
	
		// bins[ "mlp" ].load( config, "mlp" );

		allow_kaon   = config.getBool( "fit:kaon", false );
		allow_proton = config.getBool( "fit:proton", false );

		export_img = config.getBool( "can:export", false );

		vector<string> paths = config.childrenOf( "bins" );
		for ( auto p : paths ){
			string tn = config.tagName(p);
			LOG_F(  INFO, "Loading Bins: %s", tn.c_str() );
			bins[ tn ].load( config, p );
			LOG_F( INFO, "%s", bins[tn].toString().c_str() );
		}


		book->cd();
		LOG_F( INFO, "var=%s", var.c_str() );
		build_pt_projections( "sig" );
		build_pt_projections( "bg" );
		build_pt_projections( "kaon" );
		build_pt_projections( "proton" );

		nSamples = config.getInt( "nSamples", 1000 );
	}

	void spliceIn( TH1 * h1, TH2 * h2, float x ){
		assert( h1 && h2 );

		TAxis * h2x = h2->GetXaxis();
		TAxis * h2y = h2->GetYaxis();
		int ix = h2x->FindBin( x );

		assert( h2y->GetNbins() == h1->GetXaxis()->GetNbins() );

		for ( int iy = 1; iy < h2y->GetNbins(); iy++ ){
			int ibin = h2->GetBin( ix, iy );
			h2->SetBinContent( ibin, h1->GetBinContent( iy ) );
		}

	}

	void build_pt_projections( string hname ){
		LOG_SCOPE_FUNCTION( INFO );
		LOG_F( INFO, "build_pt_projections( \"%s\", var=\"%s\" )", hname.c_str(), var.c_str() );
		TH2 * h2rawpos = get<TH2>( hname + "_pos_" + var, "mc" );
		TH2 * h2rbpos = HistoBins::rebin2D( hname + "posrb", h2rawpos, bins["ptTemplate"], bins[var] );

		TH2 * h2rawneg = get<TH2>( hname + "_neg_" + var, "mc" );
		TH2 * h2rbneg = HistoBins::rebin2D( hname + "negrb", h2rawpos, bins["ptTemplate"], bins[var] );

		vector<TH1 *> projections;
		LOG_F( INFO, "nBins(pos) = %d", h2rbpos->GetXaxis()->GetNbins() );
		for ( size_t i = 0; i <= h2rbpos->GetXaxis()->GetNbins(); i++ ){
			TH1 * h = h2rbpos->ProjectionY( TString::Format( "%s_pt_%lu", hname.c_str(), i ), i+1, i+1 );
			TH1 * hn = h2rbneg->ProjectionY( TString::Format( "%s_neg_pt_%lu", hname.c_str(), i ), i+1, i+1 );
			h->Add( hn );

			LOG_F( 9, "Integral pt[%lu] = %f", i, h->Integral() );
			projections.push_back( h );
		}
		nn_pt[ hname ] = projections;
	}

	TH1 * pt_projection( string name, TH2 * hin, double pt1, double pt2, HistoBins &rb ){
		LOG_SCOPE_FUNCTION(INFO);
		LOG_F( INFO, "args=(name=%s, hin=%p, pt1=%0.2f, pt2=%0.2f)", name.c_str(), hin, pt1, pt2 );
		TAxis *x = hin->GetXaxis();
		int ipt1 = x->FindBin( pt1 );
		int ipt2 = x->FindBin( pt2 );
		LOG_F( INFO, "(%0.2f -> %0.2f) = ( %d -> %d )", pt1, pt2, ipt1, ipt2 );

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
		LOG_SCOPE_FUNCTION( INFO );
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
			LOG_F( INFO, "compare to %s", v.c_str() );
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

			TH2 * hmc_pos_kaon = get<TH2>( "kaon_pos_" + v, "mc" );
			TH2 * hmc_neg_kaon = get<TH2>( "kaon_neg_" + v, "mc" );

			TH2 * hmc_pos_proton = get<TH2>( "proton_pos_" + v, "mc" );
			TH2 * hmc_neg_proton = get<TH2>( "proton_neg_" + v, "mc" );

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
			// hproton = subsample( hproton );
			// hsig = subsample( hsig );

			hdata->Sumw2(true);
			hdata->Scale( 1.0 / hdata->Integral() );
			hsig->Scale( 1.0 / hsig->Integral() );
			hbg->Scale( 1.0 / hbg->Integral() );
			if ( use_kaon ) hkaon->Scale( 1.0 / hkaon->Integral() );
			if ( use_proton ) hproton->Scale( 1.0 / hproton->Integral() );

			rpl.style( hdata ).set( "xtitle", v + " " + config[ "units:" + v ] )
				// .set( "title", TString::Format("%0.2f < pT < %0.2f (GeV/c)", pt1, pt2 ).Data() )
				.set( "ytitle", "dN/d" + v + " " + config[ "units:" + v ] + "^{-1}" )
				.set( "lw", "1.5" )
				.set( "fc", "#999" )
				.set( "yto", 1.2 )
				.set( "xtp", 16 )
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

	void fit_purity_pt( float pt1, float pt2, TH1 * hdata, TH1 * hsig, TH1 * hbg, TH1 * hkaon, TH1 * hproton, TH1 * hsum ){

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
			ff->SetParLimits( 0, 1e-6, 1e9 );
			ff->SetParLimits( 1, 1e-6, 1e9 );
			ff->SetParLimits( 2, 1e-6, 1e9 );
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
			ff->SetParLimits( 0, 1e-14, 1.0 );
			ff->SetParLimits( 1, 1e-14, 1.0 );
			ff->SetParLimits( 2, 1e-14, 1.0 );
			ff->SetParLimits( 3, 1e-14, 1.0 );
			ff->SetParameters( 0.1, 0.1, 0.1, 1e-3 );
		}
		
		ff->SetNpx( 1000 );
		ff->SetLineColor(kBlack);
		ff->SetLineWidth(2);


		// Execute fit
		string fitOpt = config[ "fit:opt" ] + "R";
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

		report_fit_pt( ff, pt1, pt2, hdata, hsig, hbg, hkaon, hproton, hsum );

		delete ff;
	}

	void report_fit_pt( TF1 * ff, float pt1, float pt2, TH1 * hdata, TH1 * hsig, TH1 * hbg, TH1 * hkaon, TH1 * hproton, TH1 * hsum){

		RooPlotLib rpl;

		rpl.style( hdata ).set( "xtitle", var + " " + config[ "units" ] )
			.set( "ytitle", "dN/d" + var + " " + config[ "units" ] + "^{-1}" )
			// .set( "title", TString::Format("%s : %0.2f < pT < %0.2f (GeV/c)", config["charge"].c_str(), pt1, pt2 ).Data() )
			.set( "lw", "1.5" )
			.set( "fc", "#999" )
			.set( "yto", 1.3 )
			.set( "logy", 1 )
			.set( "draw", "h" )
			// .set( "xr", -0.2, 1.1 )
			.set( "min", 1e-4 )
			.set( "max", 1.0 )
			.set( "lc", "#222222" )
			.set( config, "style.data" ).draw();

		hsig->Scale( ff->GetParameter( "signal" ) );
		hbg->Scale( ff->GetParameter( "bg" ) );
		if ( use_kaon ) hkaon->Scale( ff->GetParameter( "kaon" ) );
		if ( use_proton ) hproton->Scale( ff->GetParameter( "proton" ) );

		
		hsum->Add( hbg );
		if ( use_kaon )
			hsum->Add( hkaon );
		if ( use_proton )
			hsum->Add( hproton );

		spliceIn( hsum, bgTemplate, (pt1+pt2)/2.0 );

		hsum->Add( hsig );

		rpl.style( hsig ).set( config, "style.sig" ).draw();
		rpl.style( hbg ).set( config, "style.bg" ).draw();
		if ( use_kaon ) rpl.style( hkaon ).set( config, "style.kaon" ).draw();
		if ( use_proton ) rpl.style( hproton ).set( config, "style.proton" ).draw();
		rpl.style( hsum ).set( config, "style.sum" ).draw();


		TLatex tl;
		tl.SetTextFont( 42 );
		tl.SetTextSize( 14.0 / 360.0 );
		tl.DrawLatexNDC( config.get<float>("Text:x", 0.5), 0.90, TString::Format("%0.2f < p_{T} < %0.2f (GeV/c)", pt1, pt2 ) );
		tl.DrawLatexNDC( config.get<float>("Text:x", 0.5), 0.85, TString::Format("#chi^{2} / ndf = %0.2f / %d = %0.2f", ff->GetChisquare(), ff->GetNDF(), ff->GetChisquare() / (float)ff->GetNDF() ) );
		tl.DrawLatexNDC( config.get<float>("Text:x", 0.5), 0.8, TString::Format("Yield #mu = %0.3f #pm %0.3f", ff->GetParameter( "signal" ), ff->GetParError( ff->GetParNumber( "signal" )) ) );
		tl.DrawLatexNDC( config.get<float>("Text:x", 0.5), 0.75, TString::Format("Yield #pi = %0.3f #pm %0.3f", ff->GetParameter( "bg" ), ff->GetParError( ff->GetParNumber( "bg" ))));
		if ( use_kaon ) tl.DrawLatexNDC( config.get<float>("Text:x", 0.5), 0.70, TString::Format("Yield K = %0.3f #pm %0.3f", ff->GetParameter( "kaon" ), ff->GetParError( ff->GetParNumber( "kaon" ))) );
		if ( use_proton && !use_kaon ) tl.DrawLatexNDC( config.get<float>("Text:x", 0.5), 0.70, TString::Format("Yield P = %0.3f #pm %0.3f", ff->GetParameter( "proton" ), ff->GetParError( ff->GetParNumber( "proton" ))) );
		if ( use_proton && use_kaon ) tl.DrawLatexNDC( config.get<float>("Text:x", 0.5), 0.65, TString::Format("Yield P = %0.3f #pm %0.3f", ff->GetParameter( "proton" ), ff->GetParError( ff->GetParNumber( "proton" ))) );
		

		int iptbin = book->get( "chi2ndf" )->GetXaxis()->FindBin( (pt1+pt2)/2.0 );
		book->get( "chi2ndf" )->SetBinContent( iptbin, ff->GetChisquare() / (float)ff->GetNDF() );
		book->get( "frac_mu" )->SetBinContent( iptbin, ff->GetParameter( "signal" ) );
		book->get( "frac_mu" )->SetBinError( iptbin, ff->GetParError( ff->GetParNumber( "signal" )) );

		book->get( "frac_pi" )->SetBinContent( iptbin, ff->GetParameter( "bg" ) );
		book->get( "frac_pi" )->SetBinError( iptbin, ff->GetParError( ff->GetParNumber( "bg" )) );

		if ( use_kaon ){
			book->get( "frac_k" )->SetBinContent( iptbin, ff->GetParameter( "kaon" ) );
			book->get( "frac_k" )->SetBinError( iptbin, ff->GetParError( ff->GetParNumber( "kaon" )) );
		}
		if ( use_proton ){
			book->get( "frac_p" )->SetBinContent( iptbin, ff->GetParameter( "proton" ) );
			book->get( "frac_p" )->SetBinError( iptbin, ff->GetParError( ff->GetParNumber( "proton" )) );
		}
		


		vector<float> legPos = config.getFloatVector( "TLegend[0]:pos" );
		TLegend *leg0 = new TLegend( legPos[0], legPos[1], legPos[2], legPos[3] );
		leg0->SetBorderSize(0);
		leg0->SetNColumns( config.get<int>( "TLegend[0]:ncol", 2 ) );
		leg0->SetTextSize( config.get<float>( "TLegend[0]:point" ) / 360.0 );
		leg0->AddEntry( hdata, "Run15 p+p @ #sqrt{s}=200 GeV", "f" );
		

		legPos = config.getFloatVector( "TLegend[1]:pos" );
		TLegend *leg1 = new TLegend( legPos[0], legPos[1], legPos[2], legPos[3] );
		leg1->SetBorderSize(0);
		leg1->SetNColumns( config.get<int>( "TLegend[1]:ncol", 2 ) );
		leg1->SetTextSize( config.get<float>( "TLegend[1]:point" ) / 360.0 );
		leg1->AddEntry( hsig, "#mu", "l" );
		leg1->AddEntry( hbg, "#pi", "l" );
		if ( use_kaon ) leg1->AddEntry( hkaon, "K", "l" );
		if ( use_proton ) leg1->AddEntry( hproton, "p", "l" );
		

		leg0->Draw("same");
		leg1->Draw("same");

		can->Print( (rpName).c_str()  );
		if ( export_img ) can->Print( ("export/"+ config["mod"] + "/fit-" + dts( pt1 ) + "-" + dts(pt2) + ".pdf" ).c_str() );

		string suffix="_"+dtes( pt1 ) +"_" + dtes(pt2);
		TH1 * hratio = (TH1*)hsum->Clone( ("hratio_" + suffix).c_str() );
		TH1 * hdataratio = (TH1*)hdata->Clone( ("hdataratio_" + suffix).c_str() );


		hratio->Scale( Idata );
		hdataratio->Scale( Idata );
		hratio->Divide( hdataratio );
		hratio->SetTitle( TString::Format("%0.2f < pT < %0.2f; %s; fit / data", pt1, pt2, var.c_str() ) );
		// gPad->SetLogy(0);
		
		if ( config.exists( "fit:min" ) && config.exists( "fit:max" ) ){
			hratio->Fit( "pol0", "R", "", config.get<float>( "fit:min" ), config.get<float>( "fit:max" ) );
		} else {
			hratio->Fit( "pol0" );
		}
		rpl.style( hratio )
			// .set( "xr", -0.2, 1.1 )
			.set( "logy", 0 )
			.set( "min", 0 )
			.set( "max", 2.0 )
			.set( "yto", 1.2 )
			.set( config, "style.ratio" );
		hratio->Draw();

		can->Print( (rpName).c_str()  );
		if ( export_img ) {
			can2->cd();
			hratio->Draw();
			can2->Print( ("export/"+ config["mod"] +"/ratio-" + dts( pt1 ) + "-" + dts(pt2) + ".pdf" ).c_str() );
			can->cd();
		}


		if ( !use_kaon )
			hkaon = nullptr;
		if ( !use_proton )
			hproton = nullptr;

		// hsum->Scale( 1.0/Idata );
		// hdata->Scale( 1.0/Idata );
		
		compare_other( hsig, hbg, hkaon, hproton, pt1, pt2, can );
		calc_purity( hsig, hbg, hkaon, nullptr, pt1, pt2 );
	}

	void loop_on_pt(  ){
		LOG_SCOPE_FUNCTION( INFO );
		TH2 * hdata_pos = get<TH2>( "pos_" + var + "_vs_pt", "data" );
		TH2 * hdata_neg = get<TH2>( "neg_" + var + "_vs_pt", "data" );



		TH2 * hmc_pos_sig = get<TH2>( "sig_pos_" + var, "mc" );
		TH2 * hmc_neg_sig = get<TH2>( "sig_neg_" + var, "mc" );

		TH2 * hmc_pos_bg = get<TH2>( "bg_pos_" + var, "mc" );
		TH2 * hmc_neg_bg = get<TH2>( "bg_neg_" + var, "mc" );

		TH2 * hmc_pos_kaon = get<TH2>( "kaon_pos_" + var, "mc" );
		TH2 * hmc_neg_kaon = get<TH2>( "kaon_neg_" + var, "mc" );

		TH2 * hmc_pos_proton = get<TH2>( "proton_pos_" + var, "mc" );
		TH2 * hmc_neg_proton = get<TH2>( "proton_neg_" + var, "mc" );
		

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
		
		gStyle->SetOptFit(111);
		gStyle->SetOptStat(0);

		RooPlotLib rpl;

		HistoBins pt_bins = bins["pt"];
		

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

		assert( hdata_use && hsig_use && hbg_use );
		if ( use_kaon )
			assert( hkaon_use );
		if ( use_proton )
			assert( hproton_use );


		for ( int i = 0; i < pt_bins.nBins(); i++ ){
			double pt1 = pt_bins[i];
			double pt2 = pt_bins[(i+1)];
			LOG_F( INFO, "pt=(%f, %f)", pt1, pt2 );

			LOG_F( INFO, "%0.2f = %d, %0.2f = %d", pt1, hdata_use->GetXaxis()->FindBin( pt1 ), pt2, hdata_use->GetXaxis()->FindBin( pt2 ) );
			can->Clear();

			string suffix="_"+dtes( pt1 ) +"_" + dtes(pt2);
			TH1 * hdata    = pt_projection( "data"   + suffix, hdata_use, pt1, pt2, bins[ var ] );
			TH1 * hsig     = pt_projection( "sig"    + suffix, hsig_use, pt1, pt2, bins[ var ] );
			TH1 * hbg      = pt_projection( "bg"     + suffix, hbg_use, pt1, pt2, bins[ var ] );
			TH1 * hkaon    = pt_projection( "kaon"   + suffix, hkaon_use, pt1, pt2, bins[ var ] );
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

			// hbg = subsample( hbg );
			// hproton = subsample( hproton, 1.5 );
			// hkaon = subsample( hkaon );
			// hsig = subsample( hsig );

			
			LOG_F( INFO, "hdata   = %p", hdata );
			LOG_F( INFO, "hsig    = %p", hsig );
			LOG_F( INFO, "hbg     = %p", hbg );
			LOG_F( INFO, "hkaon   = %p", hkaon );
			LOG_F( INFO, "hproton = %p", hproton );

			hdata->Sumw2();
			Idata = hdata->Integral();
			hdata->Scale( 1.0 / Idata );

			// hsig->Sumw2();
			// hbg->Sumw2();

			hsig->Scale( 1.0 / hsig->Integral() );
			hbg->Scale( 1.0 / hbg->Integral() );
			use_kaon = false;
			if ( allow_kaon)
				use_kaon = true;

			use_proton = false;
			if ( allow_proton)
				use_proton = true;
			
			if ( use_kaon )
				hkaon->Scale( 1.0 / hkaon->Integral() );
			if ( use_proton )
				hproton->Scale( 1.0 / hproton->Integral() );


			hSignalPDF = nullptr;
			hBackgroundPDF = nullptr;
			hKaonPDF = nullptr;
			hProtonPDF = nullptr;

			// hSignalPDF     = (TH1*)hsig->Clone( "hSignalPDF" );
			// hBackgroundPDF = (TH1*)hbg->Clone( "hBackgroundPDF" );
			// if ( use_kaon )
			// 	hKaonPDF = (TH1*)hkaon->Clone( "hKaonPDF" );
			// else 
			// 	hKaonPDF = nullptr;
			// if ( use_proton )
			// 	hProtonPDF = (TH1*)hproton->Clone( "hProtonPDF" );
			// else 
			// 	hProtonPDF = nullptr;

			hSignalPDF     = hsig;
			hBackgroundPDF = hbg;
			if ( use_kaon )
				hKaonPDF   = hkaon;
			if ( use_proton )
				hProtonPDF = hproton;

			TH1 * hsum = (TH1*)hbg->Clone( ("hsum"+suffix).c_str() );
			hsum->Reset();

			fit_purity_pt( pt1, pt2, hdata, hsig, hbg, hkaon, hproton, hsum );
			
			// delete hsum;
			// delete hsig;
			// delete hbg;
			// delete hSignalPDF;
			// delete hBackgroundPDF;

		}
	}


	/* Generates template pid shapes from MC
	 * Samples the input pT distribution and builds a weighted template for the given kinematics
	 *
	 */
	void generate_templates( TH1 * hpt, TH1 * hmu, TH1 * hpi, TH1 * hk, TH1 * hp ){
		LOG_SCOPE_FUNCTION( INFO );
		
		assert( hpt != nullptr );
		assert( nullptr != hmu );
		assert( nullptr != hpi );
		assert( nullptr != hk );
		assert( nullptr != hp );

		if ( hpt->Integral() <= 0 ) {
			LOG_F( INFO, "pt is empty" );
			return;
		}

		assert( nn_pt.count( "sig" ) > 0 );
		assert( nn_pt.count( "bg" ) > 0 );
		assert( nn_pt.count( "kaon" ) > 0 );
		assert( nn_pt.count( "proton" ) > 0 );

		vector<TH1 * > sig_pt = nn_pt[ "sig" ];
		vector<TH1 * > pi_pt  = nn_pt[ "bg" ];
		vector<TH1 * > k_pt  = nn_pt[ "kaon" ];
		vector<TH1 * > p_pt  = nn_pt[ "proton" ];
		

		for ( size_t i = 0; i < nSamples; i++ ){
			
			float pt1 = hpt->GetRandom();
			if ( pt1 > 3 ) pt1 = 3;
			
			int ipt = bins[ "ptTemplate" ].findBin( pt1 );
			if ( ipt < 0 || ipt >= sig_pt.size() ) continue;

			float rMu = sig_pt[ipt]->GetRandom();
			float rPi = pi_pt[ipt]->GetRandom();
			float rK  =  k_pt[ipt]->GetRandom();
			
			float rP  =  -999;

			if ( p_pt.size() > ipt && hp != nullptr && pt1 < 3.0){
				rP = p_pt[ipt]->GetRandom();
			}

			// LOG_F( 9, "pt=%f, ipt=%d, rMu=%f, rPi=%f, rK=%f, rP=%f", pt1, ipt, rMu, rPi, rK, rP );

			if ( rMu <= 0 && "mlp" == var )
				rMu = -999;

			hmu->Fill( rMu );
			hpi->Fill( rPi );
			hk ->Fill( rK  );
			hp->Fill( rP );
		}
	} // generate_templates

	TF1 * fit_pid( TFitResultPtr * fr, TH1 * hvar, TH1 * hmu, TH1 * hpi, TH1 * hk = nullptr, TH1 * hp = nullptr ){

		TF1 * ff = nullptr;
		if ( nullptr == hk){
			ff = new TF1( "ff", fitfun, -1, 2, 2 );
			ff->SetParNames( "muon", "pion" );
			ff->SetParLimits( 0, 1e-4, 1e9 );
			ff->SetParLimits( 1, 1e-4, 1e9 );
			ff->SetParameters( 1.0, 100.0 );
		} else if ( nullptr == hp && nullptr != hk ) {
			ff = new TF1( "ff", fitfun3, -1, 2, 3 );
			ff->SetParNames( "muon", "pion", "kaon" );
			ff->SetParLimits( 0, 1e-16, 1.0 );
			ff->SetParLimits( 1, 1e-10, 1.0 );
			ff->SetParLimits( 2, 1e-10, 1.0 );
			ff->SetParameters( 1.0, 0.1, 0.01 );
			// ff->FixParameter( 1, 0.00001 );
		} else if ( nullptr == hk && nullptr != hp ){
			ff = new TF1( "ff", fitfun3, -1, 2, 3 );
			ff->SetParNames( "muon", "pion", "proton" );
			ff->SetParLimits( 0, 1e-10, 1.0 );
			ff->SetParLimits( 1, 1e-10, 1.0 );
			ff->SetParLimits( 2, 1e-10, 1.0 );
			ff->SetParameters( 0.1, 0.1, 0.50 );
		} else if ( nullptr != hk && nullptr != hp ){
			ff = new TF1( "ff", fitfun4, -1, 2, 4 );
			ff->SetParNames( "muon", "pion", "kaon", "proton" );
			ff->SetParLimits( 0, 1e-10, 1.0 );
			ff->SetParLimits( 1, 1e-10, 1.0 );
			ff->SetParLimits( 2, 1e-10, 1.0 );
			ff->SetParLimits( 3, 1e-10, 1.0 );
			ff->SetParameters( 1.0, 0.10, 0.50, 1e-3 );
		}
		
		ff->SetNpx( 1000 );
		ff->SetLineColor(kBlack);
		ff->SetLineWidth(2);


		// Execute fit
		string fitOpt = config[ "fit:opt" ] + "R S";
		if ( config.exists( "fit:min" ) && config.exists( "fit:max" ) ){
			LOG_F( INFO, "Fitting (opt=%s) in range (%f, %f)", fitOpt.c_str(), config.get<float>( "fit:min" ), config.get<float>( "fit:max" ) );

			(*fr) = hvar->Fit( ff, fitOpt.c_str(), "", config.get<float>( "fit:min" ), config.get<float>( "fit:max" ) );
		} else {
			(*fr) = hvar->Fit( ff, fitOpt.c_str() );
		}

		return ff;
	}

	void norm( TH1 * h ){
		assert( nullptr != h );
		h->Scale( 1.0 / h->Integral() );
	}

	void loop_on_mass( string scharge = "pos" ){
		LOG_SCOPE_FUNCTION( INFO );

		RooPlotLib rpl;

		TH2 * h2pt  = get<TH2>( scharge + "_pt_vs_mass", "data" );
		TH2 * h2var = get<TH2>( scharge + "_" + var + "_vs_mass", "data" );
		assert( nullptr != h2var );
		assert( nullptr != h2pt );

		TH2 * h2varrb = HistoBins::rebin2D( scharge + "_" + var + "_vs_mass_rb", h2var, bins["mass"], bins[ var ] );
		TH2 * h2ptrb  = HistoBins::rebin2D( scharge + "_pt_vs_mass_rb", h2pt, bins["mass"], bins["ptTemplate"] );
		assert( nullptr != h2varrb );
		assert( nullptr != h2ptrb );

		gStyle->SetOptFit(111);
		gStyle->SetOptStat(0);

		can->cd();
		for ( size_t iMass = 1; iMass <= h2varrb->GetXaxis()->GetNbins(); iMass ++  ){
			TH1 * hvar = (TH1*)h2varrb->ProjectionY( TString::Format( "%s_%s_vs_mass_%lu", scharge.c_str(), var.c_str(), iMass), iMass, iMass );
			TH1 * hpt  = (TH1*)h2ptrb ->ProjectionY( TString::Format( "%s_pt_vs_mass_%lu", scharge.c_str(), iMass), iMass, iMass );

			float m1 = h2varrb->GetXaxis()->GetBinLowEdge( iMass );
			float m2 = h2varrb->GetXaxis()->GetBinUpEdge( iMass );

			TH1 * hmu = new TH1F( TString::Format( "template_%s_mu_m%lu", scharge.c_str(), iMass ), "", bins[ var ].nBins(), bins[ var ].getBins().data() );
			TH1 * hpi = new TH1F( TString::Format( "template_%s_pi_m%lu", scharge.c_str(), iMass ), "", bins[ var ].nBins(), bins[ var ].getBins().data() );
			TH1 * hk  = new TH1F( TString::Format( "template_%s_k_m%lu", scharge.c_str(), iMass ), "", bins[ var ].nBins(), bins[ var ].getBins().data() );
			TH1 * hp  = new TH1F( TString::Format( "template_%s_p_m%lu", scharge.c_str(), iMass ), "", bins[ var ].nBins(), bins[ var ].getBins().data() );

			generate_templates( hpt, hmu, hpi, hk, hp );

			hvar->Sumw2();
			hmu->Sumw2();
			hpi->Sumw2();
			hk->Sumw2();
			hp->Sumw2();

			norm( hmu );
			norm( hpi );
			norm( hk );
			norm( hp );
			norm( hvar );

			hSignalPDF = hmu;
			hBackgroundPDF = hpi;
			hKaonPDF = hk;
			hProtonPDF = hp;

			TF1 * ff = nullptr;
			TFitResultPtr fr = nullptr;
			ff = fit_pid( &fr, hvar, hmu, hpi, hk, hp );

			int iptbin = book->get( "chi2ndf" )->GetXaxis()->FindBin( (m1+m2)/2.0 );
			book->get( "chi2ndf" )->SetBinContent( iptbin, ff->GetChisquare() / (float)ff->GetNDF() );
			book->get( "frac_mu" )->SetBinContent( iptbin, ff->GetParameter( "muon" ) );
			book->get( "frac_mu" )->SetBinError( iptbin, ff->GetParError( ff->GetParNumber( "muon" )) );

			book->get( "frac_pi" )->SetBinContent( iptbin, ff->GetParameter( "pion" ) );
			book->get( "frac_pi" )->SetBinError( iptbin, ff->GetParError( ff->GetParNumber( "pion" )) );

			book->get( "frac_k" )->SetBinContent( iptbin, ff->GetParameter( "kaon" ) );
			book->get( "frac_k" )->SetBinError( iptbin, ff->GetParError( ff->GetParNumber( "kaon" )) );

			book->get( "frac_p" )->SetBinContent( iptbin, ff->GetParameter( "proton" ) );
			book->get( "frac_p" )->SetBinError( iptbin, ff->GetParError( ff->GetParNumber( "proton" )) );
			

			rpl.style( hvar ).set( config, "style.data" ).draw();
			// gPad->SetLogy(1);

			hmu->Scale( ff->GetParameter( "muon" ) );
			hpi->Scale( ff->GetParameter( "pion" ) );
			hk->Scale( ff->GetParameter( "kaon" ) );
			hp->Scale( ff->GetParameter( "proton" ) );

			TH1 * hsum = (TH1*)hmu->Clone( TString::Format( "template_%s_sum_m%lu", scharge.c_str(), iMass ) );
			hsum->Add( hpi );
			hsum->Add( hk );
			hsum->Add( hp );

			rpl.style( hmu ).set(config, "style.sig" );
			rpl.style( hpi ).set(config, "style.bg" );
			rpl.style( hk ).set(config, "style.kaon" );
			rpl.style( hp ).set(config, "style.proton" );
			rpl.style( hsum ).set(config, "style.sum" );

			hmu->Draw( "same hist" );
			hpi->Draw( "same hist" );
			hk->Draw( "same hist" );
			hp->Draw( "same hist" );
			hsum->Draw( "same hist" );

			spliceIn( hsum, bgTemplate, (m1+m2)/2.0 );
			spliceIn( hmu,  sigTemplate, (m1+m2)/2.0 );

			TLatex tl;
			tl.SetTextSize( 12.0 / 360.0 );
			tl.DrawLatexNDC( 0.70, 0.90, TString::Format( "%0.2f < M < %0.2f (GeV/c^{2})", m1, m2 ) );
			tl.DrawLatexNDC( 0.18, 0.85, TString::Format("#chi^2 / ndf = %0.2f / %d = %0.2f", ff->GetChisquare(), ff->GetNDF(), ff->GetChisquare() / (float)ff->GetNDF() ) );
			tl.DrawLatexNDC( 0.18, 0.80, TString::Format("muon = %0.4f #pm %0.4f", ff->GetParameter( "muon" ), ff->GetParError( ff->GetParNumber( "muon" ) )) );
			tl.DrawLatexNDC( 0.58, 0.85, TString::Format("pion = %0.4f #pm %0.4f", ff->GetParameter( "pion" ), ff->GetParError( ff->GetParNumber( "pion" ) )) );
			tl.DrawLatexNDC( 0.58, 0.80, TString::Format("kaon = %0.4f #pm %0.4f", ff->GetParameter( "kaon" ), ff->GetParError( ff->GetParNumber( "kaon" ) )) );
			tl.DrawLatexNDC( 0.18, 0.75, TString::Format("proton = %0.4f #pm %0.4f", ff->GetParameter( "proton" ), ff->GetParError( ff->GetParNumber( "proton" ) )) );

			TLegend *leg = new TLegend( 0.18, 0.88, 0.75, 0.95 );
			leg->SetBorderSize( 0 );
			leg->SetNColumns(5);
			leg->AddEntry( hvar, "data" );
			leg->AddEntry( hmu, "#mu" );
			leg->AddEntry( hpi, "#pi" );
			leg->AddEntry( hk, "K" );
			leg->AddEntry( hsum, "Total Fit" );

			leg->Draw("same");

			can->Print( (rpName).c_str()  );
			if ( export_img )
				can->Print( ("export/"+ config["mod"] +"/fit-" + var + "-" + dts( m1 ) + "-" + dts(m2) + ".pdf" ).c_str() );

			delete ff;

		}
	}

	virtual void make(){
		LOG_SCOPE_FUNCTION( INFO );

		can = new TCanvas( "can", "can", config.get<int>( "can:w", 500 ), config.get<int>( "can:h", 500 ) );
		can2 = new TCanvas( "can2", "can2", config.get<int>( "can:w", 500 ), config.get<int>( "can:h2", 500 ) );

		can->SetTopMargin( 0.05 );
		can->SetRightMargin( 0.03 );
		can->SetBottomMargin( 0.11 );
		can->SetLeftMargin( 0.15 );

		can2->SetTopMargin( 0.05 );
		can2->SetRightMargin( 0.01 );
		can2->SetBottomMargin( 0.11 );
		can2->SetLeftMargin( 0.075 );


		rpName = config[nodePath + ".output.Report:url" ];
		can->Print( (rpName+"[").c_str() );

		book->cd();
		book->makeAll( config, nodePath + ".histograms" );
		bgTemplate  = (TH2*)book->get( "bg_template" );
		sigTemplate = (TH2*)book->get( "sig_template" );

		if ( "mass" == config.getString( "loop", "" ) ){
			loop_on_mass();
		} else {
			loop_on_pt();
		}

		if ( book->get("chi2ndf") ){
			RooPlotLib rpl;
			rpl.link( book );

			can->Clear();
			can->cd(0);
			can->SetLogy(0);
			rpl.style( "chi2ndf" ). set( config, "style.chi2ndf" ).draw();
			can->Print( (rpName).c_str()  );


			for ( string n : { "frac_mu", "frac_pi", "frac_k", "frac_p" } ){
				rpl.style( n ). set( config, "style." + n ) .set( config, "style.yield_frac" ).draw();
				can->Print( (rpName).c_str() );
			}

			rpl.style( "bg_template" ).set( config, "style.bg_template" ).draw();
			can->Print( (rpName).c_str() );

			rpl.style( "sig_template" ).set( config, "style.sig_template" ).draw();
			can->Print( (rpName).c_str() );
		}

		

		if ( 0 == config.getInt( "jobIndex" ) || -1 == config.getInt( "jobIndex" ) ){
			TNamed config_str( "config", config.toXml() );
			config_str.Write();
		}


		if ( nullptr != can )
			can->Print( (rpName+"]").c_str()  );
	} // make

};



#endif