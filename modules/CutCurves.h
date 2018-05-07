
#ifndef CUT_CURVES_H
#define CUT_CURVES_H


#include "HistoAnalyzer.h"
#include "vendor/loguru.h"

#include "RooPlotLib.h"


class CutCurves : public HistoAnalyzer {
protected:
public:

	CutCurves(){}
	~CutCurves(){}

	map<string, HistoBins> bins;


	virtual void initialize(){
		LOG_F( INFO, "" );
		book->cd();


		vector<string> paths = config.childrenOf( "bins" );
		for ( auto p : paths ){
			string tn = config.tagName(p);
			LOG_F(  INFO, "Loading Bins: %s", tn.c_str() );
			bins[ tn ].load( config, p );
			LOG_F( INFO, "%s", bins[tn].toString().c_str() );
		}

	}

	void make_purity_curve( string suffix, TH1 * hdata, TH1 * hmu, TH1 * hpion, TH1 * hkaon, TH1 * hprot, TH1 * hsum ){

		RooPlotLib rpl;

		TH1 * hpurity = (TH1*)hdata->Clone( ("purity" + suffix).c_str() );
		hpurity->Reset();

		TAxis * ax = hpurity->GetXaxis();
		int nBins = ax->GetNbins();
		for ( int i = 1; i <= nBins; i++ ){
			double Emu = 0, Epion = 0, Ekaon = 0, Eprot = 0;
			double Imu   = hmu->IntegralAndError( i, nBins, Emu );
			double Ipion = hpion->IntegralAndError( i, nBins, Epion );
			double Ikaon = 0;
			if ( nullptr != hkaon ) 
				hkaon->IntegralAndError( i, nBins, Ekaon );
			double Iprot = 0;
			if ( nullptr != hprot ) 
				hprot->IntegralAndError( i, nBins, Eprot );

			double Itotal = Imu + Ipion + Ikaon + Iprot;
			double p = Imu / Itotal;

			if ( Itotal <= 0 || p != p ) continue;

			hpurity->SetBinContent( i, p );
			hpurity->SetBinError( i, (Emu /Imu) * (p) );
			// LOG_F( INFO, "purity=%0.3f, Emu=%0.3f, Imu=%0.3f", p, Emu, Imu );
			rpl.style( hpurity ).set( config, "style.purity" );
		}
	} // make_purity_curve



	void make_eff_curve( string suffix, TH1 * h, string name ){
		RooPlotLib rpl;

		TH1 * heff = (TH1*)h->Clone( ("eff_" + name + suffix).c_str() );
		heff->Reset();

		TAxis * ax = heff->GetXaxis();
		int nBins = ax->GetNbins();
		for ( int i = 1; i <= nBins; i++ ){
			
			double Ecut   = 0, Etot = 0;
			double Icut   = h->IntegralAndError( i, nBins, Ecut );
			double Itotal = h->IntegralAndError( 1, nBins, Etot);
			double eff = Icut / Itotal;

			if ( Itotal <= 0 || eff != eff ) continue;

			heff->SetBinContent( i, eff );
			heff->SetBinError( i, (Ecut /Icut) * (eff) );
			
			rpl.style( heff ).set( config, "style.eff" );
		}
	} // make_eff_curve

	void make_roc_curve( string suffix, TH1 * hdata, TH1 * hmu, TH1 * hpion, TH1 * hkaon, TH1 * hprot, TH1 * hsum  ){

		TH1 * hroc = new TH1F( ("roc" + suffix).c_str(), ";#varepsilon_{signal}; 1 - #varepsilon_{bg}", 100, 0, 1.0 );

		TAxis * rax = hroc->GetXaxis();
		TAxis * nax = hmu->GetXaxis();
		for ( size_t i = 1; i <= rax->GetNbins(); i++ ){
			float targetEff = rax->GetBinCenter( i );
			double effBg = 1.0;

			for ( size_t j = 1; j <= nax->GetNbins(); j++ ){
				double Icut   = hmu->Integral( j, nax->GetNbins());
				double Itotal = hmu->Integral( 1, nax->GetNbins());
				double effMu = Icut / Itotal;
				if ( effMu < targetEff ) break;

				double IPicut   = hpion->Integral( j, nax->GetNbins());
				double IPitotal = hpion->Integral( 1, nax->GetNbins());

				double IKaoncut   = 0;
				double IKaontotal = 0;
				if ( nullptr != hkaon ){
					IKaoncut   = hkaon->Integral( j, nax->GetNbins());
					IKaontotal = hkaon->Integral( 1, nax->GetNbins());
				}
				
				double IProtoncut   = 0;
				double IProtontotal = 0;
				if ( nullptr != hprot ){
					IProtoncut   = hprot->Integral( j, nax->GetNbins());
					IProtontotal = hprot->Integral( 1, nax->GetNbins());
				}

				effBg = ( IPicut + IKaoncut + IProtoncut ) / ( IPitotal + IKaontotal + IProtontotal );

			} // loop to find bg efficiency at target signal efficiency

			// LOG_F( INFO, "(%0.3f, %0.3f)", targetEff, effBg );
			hroc->SetBinContent( i, 1.0 - effBg );

		}// loop rax bins
	} // make_roc_curve


	virtual void make(){
		LOG_F( INFO, "" );

		assert( bins.count("pt") && bins["pt"].nBins() >= 1 );

		for ( size_t i = 0; i < bins["pt"].nBins(); i++ ){

			double pt1 = bins["pt"][i];
			double pt2 = bins["pt"][(i+1)];

			string suffix= "_" + dtes( pt1 ) + "_" + dtes(pt2);

			TH1 * hdata = get<TH1>( "data" + suffix + "" );
			TH1 * hmu   = get<TH1>( "sig" + suffix + "_rb" );
			TH1 * hpion = get<TH1>( "bg" + suffix + "_rb" );
			TH1 * hkaon = get<TH1>( "kaon" + suffix + "_rb" );
			TH1 * hprot = get<TH1>( "proton" + suffix + "_rb" );
			TH1 * hsum  = get<TH1>( "hsum" + suffix );

			hmu->Scale( hdata->Integral() );
			hpion->Scale( hdata->Integral() );

			make_purity_curve( suffix, hdata, hmu, hpion, hkaon, hprot, hsum );
			make_eff_curve( suffix, hmu, "mu" );
			make_roc_curve( suffix, hdata, hmu, hpion, hkaon, hprot, hsum );

		}


	}
};

#endif