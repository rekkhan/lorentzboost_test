#include <stdio.h>
#include <math.h>
#include <iostream>

#include "TH1D.h"
#include "TCanvas.h"
#include "TPad.h"

#define higgs_mass 125
#define zboson_mass 91
#define pival 3.14159265




void VisualSetup ( TH1D* hist, int colour, int fillstyle, const char *xtitle )
{
	hist -> Sumw2( );
	hist -> SetFillStyle ( fillstyle );
	hist -> SetFillColor ( colour );
	hist -> SetLineWidth ( 2 );
	hist -> SetLineColor ( colour );

	hist -> GetXaxis( ) -> SetTitle ( xtitle );
	hist -> GetXaxis( ) -> SetTitleSize ( 0.05 );
	hist -> GetXaxis( ) -> SetTitleOffset ( 0.8 );

	hist -> GetYaxis( ) -> SetTitle ( "nEvent" );
	hist -> GetYaxis( ) -> SetTitleSize ( 0.05 );
	hist -> GetYaxis( ) -> SetTitleOffset ( 1.0 );
	hist -> GetYaxis( ) -> SetNdivisions ( 510 );
}


void PadSetup ( TPad *pad )
{
	pad -> SetLeftMargin ( 0.12 );
	pad -> SetRightMargin ( 0.05 );
	pad -> SetTopMargin (0.05 );
	pad -> SetBottomMargin (0.12 );
	pad -> SetGrid ( 0, 1 );
}





double func_decayangle (
		double *variables,
		double *parameters )
{
	double angle = variables[0];
	double exp_steepness = parameters[0];
	double lin_steepness = parameters[1];
	double vert_addition = parameters[2];

	double pdfvalue = 0;
	pdfvalue += exp ( exp_steepness * angle );
	pdfvalue += lin_steepness*angle + vert_addition;

	double result = ( (angle>=0) && (angle<pival) ) * pdfvalue;

	return result;
}





double func_initialeta (
		double *variables,
		double *parameters )
{
	double eta    = variables[0];
	double width     = parameters[0];
	double constrain = parameters[1] * width;

	double pdfvalue = 0;
	pdfvalue += ( abs(eta) < constrain ) * exp( - pow(eta/width, 2) / 2 );
	pdfvalue += ( abs(eta) < constrain ) * 0.02*width*pow( eta, 2 );

	double result = pdfvalue;

	return result;
}





int boostmultivector ( int method )
{
	// + NOTE on this code:
	//-------------------------------------------------------------------------------
	// [+] This code computes the cosine(Theta) by 2 different definition
	//  |-- The first definition compute the boosted angle between the Z boson
	//  |   and the Higgs boson.
	//  |   |-- This definition require both the Z boson and the Higgs boson
	//  |   |   being rotated before boosting,
	//  |   |-- This rotation rotate the momentum of the Higgs boson to the z-axis.
	//  |   `-- The boost vector is defined by the velocity of the Higgs boson
	//  |       after rotation
	//  |
	//  `-- The second definition consider the angle between the boosted Z boson
	//      and its new z-axis,
	//      `-- The boost vector is defined by the velocity of the Higgs boson
	//          after rotation
	// [+] Usage:
	//  |-- root -l -b -q boostmultivector.C\(<method>\)
	//  `-- method = 1 or 2
	//-------------------------------------------------------------------------------
	printf ( " [+] Program starts\n" );

	gErrorIgnoreLevel = kError, kBreak, kSysError, kFatal;
	gStyle -> SetOptStat( 0 );



	printf ( "  |-- Create histograms\n" );
	TH1D *hist_HiggsEta   = new TH1D ( "hist_HiggsEta", "", 30, -3.0, 3.0 );
	TH1D *hist_ZEta       = new TH1D ( "hist_ZEta",     "", 30, -3.0, 3.0 );
	TH1D *hist_PhoEta     = new TH1D ( "hist_PhoEta",   "", 30, -3.0, 3.0 );
	TH1D *hist_OrgCos     = new TH1D ( "hist_OrgCos",   "", 24, -1.2, 1.2 );
	TH1D *hist_BoostCos   = new TH1D ( "hist_BoostCos", "", 24, -1.2, 1.2 );
	TH1D *hist_OrgTheta   = new TH1D ( "hist_OrgTheta",   "", 36, -0.2, 3.4 );
	TH1D *hist_BoostTheta = new TH1D ( "hist_BoostTheta", "", 36, -0.2, 3.4 );

	VisualSetup ( hist_HiggsEta,   kOrange-3, 1001, "#eta(Higgs)" );
	VisualSetup ( hist_ZEta,       kAzure+2,  3345, "#eta(Z)" );
	VisualSetup ( hist_PhoEta,     kOrange,   3354, "#eta(#gamma)" );
	VisualSetup ( hist_OrgCos,     kBlack,    3345, "cos(#Theta)" );
	VisualSetup ( hist_BoostCos,   kTeal+2,   3354, "cos(#Theta)" );
	VisualSetup ( hist_OrgTheta,   kBlack,    3345, "#Theta" );
	VisualSetup ( hist_BoostTheta, kTeal+2,   3354, "#Theta" );



	// + Some function to determine the random angle
	//----------------------------------------------
	printf ( "  |-- Create functions\n" );

	// * The decay angle of Z boson in the Higgs' rest frame
	TF1 *pdf_decaytheta = new TF1 (
			"pdf_decaytheta",
			func_decayangle,
			0, pival, 3 );
	pdf_decaytheta -> SetParameters ( -0.5, -0.00, 0.05 );

	// * The polar angle of Higgs in the lab frame
	TF1 *pdf_higgseta = new TF1 (
			"pdf_higgseta",
			func_initialeta,
			-5, 5, 2 );
	pdf_higgseta -> SetParameters ( 1.8, 2.5 );



	TLorentzVector vec4_Higgs;
	TLorentzVector vec4_HiggsTmp;
	TLorentzVector vec4_System;
	TLorentzVector vec4_Zboson;
	TLorentzVector vec4_Photon;
	TLorentzVector vec4_Sys;
	TVector3       boostvector;

	float rotTheta;
	float rotPhi;


	TRandom3 *randomgen = new TRandom3( 0 );

	double higgs_momentum;
	double higgs_energy;

	double zboson_momentum;
	double zboson_energy;

	double photon_momentum;
	double photon_energy;

	double theta_higgsrestframe;
	double phi_higgsrestframe;

	double eta_HiggsInit;
	double theta_HiggsInit;
	double phi_HiggsInit;

	double cosine_theta;
	double boost_theta;



	printf ( "  |-- Start loop\n" );

	long nevent = 100000;

	for ( long eventloop=0; eventloop<nevent; eventloop++ )
	{
		// + Randomly generate particle
		//-----------------------------
		// * The Higgs boson
		higgs_momentum = abs ( randomgen->Gaus(30, 5) );
		higgs_energy   = sqrt (pow( higgs_mass, 2 ) + pow ( higgs_momentum, 2 ) );
		vec4_Higgs . SetPxPyPzE ( 0.0, 0.0, higgs_momentum, higgs_energy );

		// * the Z boson decay from Higgs at rest
		zboson_momentum = ( pow(higgs_mass, 2) - pow(zboson_mass, 2) ) / ( 2*higgs_mass );
		zboson_energy   = sqrt ( pow( zboson_mass, 2 ) + pow( zboson_momentum, 2 ) );
		vec4_Zboson . SetPxPyPzE ( 0.0, 0.0, zboson_momentum, zboson_energy );

		// * the Photon decay from Higgs at rest
		photon_momentum = -zboson_momentum;
		photon_energy   = abs ( photon_momentum );
		vec4_Photon . SetPxPyPzE ( 0.0, 0.0, photon_momentum, photon_energy );


		// + Rotate && boost the higgs rest frame
		//---------------------------------------
		// * Rotation by the decay angle of Z
		theta_higgsrestframe = pdf_decaytheta -> GetRandom ( 0, pival, randomgen );
		vec4_Zboson . RotateY ( theta_higgsrestframe );
		vec4_Photon . RotateY ( theta_higgsrestframe );

		// * Rotation randomly around z-axis
		phi_higgsrestframe = randomgen -> Uniform ( 0, 2*pival );
		vec4_Zboson . RotateZ ( phi_higgsrestframe );
		vec4_Photon . RotateZ ( phi_higgsrestframe );

		// * Lorentz boost to lab-frame
		vec4_Photon . Boost ( vec4_Higgs.BoostVector( ) );
		vec4_Zboson . Boost ( vec4_Higgs.BoostVector( ) );


		// * Rotate the whole system to mimic "real" Higgs decay
		eta_HiggsInit = pdf_higgseta -> GetRandom (-2.5, 2.5, randomgen);
		eta_HiggsInit += ( eta_HiggsInit >= 2.5 ) * abs ( 2.5 - abs(eta_HiggsInit) );
		eta_HiggsInit -= ( eta_HiggsInit <= -2.5 ) * abs ( 2.5 - abs(eta_HiggsInit) );
		theta_HiggsInit = 2 * atan ( exp( -eta_HiggsInit ) );
		phi_HiggsInit = randomgen -> Uniform ( 0, 2*pival );

		vec4_Higgs  . RotateY ( theta_HiggsInit );
		vec4_Zboson . RotateY ( theta_HiggsInit );
		vec4_Photon . RotateY ( theta_HiggsInit );

		vec4_Higgs  . RotateZ ( phi_HiggsInit );
		vec4_Zboson . RotateZ ( phi_HiggsInit );
		vec4_Photon . RotateZ ( phi_HiggsInit );

		double thetalab_Higgs = vec4_Higgs . Theta( );
		double philab_Higgs = vec4_Higgs . Phi( );



		// + Fill histogram
		//-----------------
		// * for full Z's eta, comment the line below
		//if ( vec4_Zboson.Eta( ) < 0 )   continue;
		hist_HiggsEta -> Fill ( vec4_Higgs.Eta() );
		hist_ZEta     -> Fill ( vec4_Zboson.Eta() );
		hist_PhoEta   -> Fill ( vec4_Photon.Eta() );
		hist_OrgCos   -> Fill ( cos( theta_higgsrestframe ) );
		hist_OrgTheta -> Fill ( ( theta_higgsrestframe ) );



		// + Boost procedure
		//------------------
		if ( method == 2 )
		{
			// * Rotate the ref frame of Z to make its new Z-axis match
			// the momentum of the Higgs boson
			// * The frame must be rotate around the Z axis by Higgs' phi first,
			// then around the Y axis by Higgs' theta
			vec4_Zboson . RotateZ ( -phi_HiggsInit );
			vec4_Zboson . RotateY ( -theta_HiggsInit );

			// * The Higgs boson must also be rotated by the same amount to get
			// the correct boost vector (pointing same direction as the new Z axis)
			vec4_Higgs . RotateZ ( -phi_HiggsInit );
			vec4_Higgs . RotateY ( -theta_HiggsInit );
		}

		vec4_Zboson . Boost ( -vec4_Higgs.BoostVector( ) );
		cosine_theta = cos ( vec4_Zboson.Theta( ) );
		boost_theta  = vec4_Zboson.Theta( );

		// * remove the if statement to calculate boost without rotation
		if ( method == 1 )
		{
			cosine_theta = cos ( vec4_Zboson.Angle( vec4_Higgs.Vect() ) );
			boost_theta  = vec4_Zboson.Angle( vec4_Higgs.Vect( ) );
		}

		hist_BoostCos   -> Fill ( cosine_theta );
		hist_BoostTheta -> Fill ( boost_theta );
	}



	TCanvas *canvas = new TCanvas ( "canvas", "", 1200, 800 );

	canvas -> cd( );
	TPad *pad1 = new TPad ( "pad1", "", 0.0, 0.90, 1.0, 1.0 );
	pad1 -> Draw( );
	pad1 -> cd( );
	std::string infotex = "";
	infotex += ( method==1 ) ? "No rotation before boost" : "";
	infotex += ( method==2 ) ? "Rotate the lab frame to match Higgs's momentum before boost" : "";
	TLatex *tex = new TLatex( );
	tex -> SetNDC( );
	tex -> SetTextSize  ( 0.3 );
	tex -> SetTextAlign ( 22 );
	tex -> DrawLatex ( 0.5, 0.75, infotex.data( ) );
	tex -> DrawLatex ( 0.5, 0.25, "#Theta: Polar angle of Z boson in the Higgs' rest frame" );


	canvas -> cd( );
	TPad *pad2 = new TPad ( "pad2", "", 0.0, 0.45, 0.5, 0.90 );
	PadSetup ( pad2 );
	pad2 -> Draw( );
	pad2 -> cd( );
	hist_HiggsEta -> Draw ( "hist" );


	canvas -> cd( );
	TPad *pad3 = new TPad ( "pad3", "", 0.5, 0.45, 1.0, 0.90 );
	PadSetup ( pad3 );
	pad3 -> Draw( );
	pad3 -> cd( );

	double maxEta = 1.3 * max (
			hist_ZEta->GetMaximum( ),
			hist_PhoEta->GetMaximum( ) );
	hist_ZEta -> SetMaximum ( maxEta );

	hist_ZEta   -> Draw ( "hist" );
	hist_PhoEta -> Draw ( "hist same" );

	TLegend *legendEta = new TLegend ( 0.13, 0.84, 0.94, 0.94 );
	legendEta -> SetNColumns ( 2 );
	legendEta -> SetTextSize ( 0.05 );
	legendEta -> AddEntry ( hist_ZEta, "Z boson" );
	legendEta -> AddEntry ( hist_PhoEta, "Photon" );
	legendEta -> Draw ( "same" );


	canvas -> cd( );
	TPad *pad4 = new TPad ( "pad4", "", 0.00, 0.0, 0.50, 0.45 );
	PadSetup ( pad4 );
	pad4 -> Draw( );
	pad4 -> cd( );

	double maxCos = 1.5 * max (
			hist_OrgCos->GetMaximum( ),
			hist_BoostCos->GetMaximum( ) );
	hist_OrgCos -> SetMaximum ( maxCos );
	hist_OrgCos  -> Draw ( "hist" );
	hist_BoostCos -> Draw ( "hist same" );

	TLegend *legendCos = new TLegend ( 0.13, 0.75, 0.94, 0.94 );
	legendCos -> SetTextSize ( 0.05 );
	legendCos -> AddEntry ( hist_OrgCos, "cos(decay angle of Z) (Higgs' rest frame)" );
	legendCos -> AddEntry ( hist_BoostCos, "cos(#Theta) (Higgs' rest frame)" );
	legendCos -> Draw ( "same" );


	canvas -> cd( );
	TPad *pad5 = new TPad ( "pad5", "", 0.50, 0.0, 1.00, 0.45 );
	PadSetup ( pad5 );
	pad5 -> Draw( );
	pad5 -> cd( );

	double maxTheta = 1.5 * max (
			hist_OrgTheta->GetMaximum( ),
			hist_BoostTheta->GetMaximum( ) );
	hist_OrgTheta -> SetMaximum ( maxTheta );
	hist_OrgTheta   -> Draw ( "hist" );
	hist_BoostTheta -> Draw ( "hist same" );

	TLegend *legendTheta = new TLegend ( 0.13, 0.75, 0.94, 0.94 );
	legendTheta -> SetTextSize ( 0.05 );
	legendTheta -> AddEntry ( hist_OrgTheta, "decay angle of Z (Higgs' rest frame)" );
	legendTheta -> AddEntry ( hist_BoostTheta, "#Theta (Higgs' rest frame)" );
	legendTheta -> Draw ( "same" );


	canvas -> SaveAs ( Form( "plot_method%d.png", method ) );
	printf ( "  |-- Save plot\n" );
	printf ( "  `-- Program exit\n" );


	return 0;
}
