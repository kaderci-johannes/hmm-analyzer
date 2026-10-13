#pragma once

#include "TString.h"
#include "TH1D.h"
#include "TDirectory.h"

namespace analysis
{
	namespace processing
	{
		struct DimuonSet
		{
			DimuonSet(TString const& postfix, const double scale = 1.)
			{
				_postfix = postfix;
			}

			void init()
			{
				gDirectory->mkdir(_postfix);
				gDirectory->cd(_postfix);
//				hNlead = new TH1D("hNlead", "hNlead", 11, -.5, 10.5);
				hDiJetMass = new TH1D("DiJetMass", "DiJetMass",
					1000, 0., 1000.);
				hDiJetMass->Sumw2(kTRUE);
				hDiJetdeta = new TH1D("DiJetdeta", "DiJetdeta",
					50, 0., 5.);
				hDiJetdeta->Sumw2(kTRUE);
				hDiJeteta = new TH1D("DiJeteta", "DiJeteta",
					100, -5., 5.);
				hDiJeteta->Sumw2(kTRUE);
				hDiJetPt = new TH1D("DiJetPt", "DiJetPt",
					400, 0., 400.);
				hDiJetPt->Sumw2(kTRUE);
				hDiMuonpt = new TH1D("DiMuonpt", "DiMuonpt",
					200, 0., 200.);
				hDiMuonpt->Sumw2(kTRUE);
				hDiMuonMass_full = new TH1D("DiMuonMass_full",
					"DiMuonMass_full", 200, 0, 200);
				hDiMuonMass_full->Sumw2(kTRUE);
				hDiMuonMass_Jpsi = new TH1D("DiMuonMass_Jpsi",
					"DiMuonMassJpsi", 10, 1, 6);
				hDiMuonMass_Jpsi->Sumw2(kTRUE);
				hDiMuonMass_Upsi = new TH1D("DiMuonMass_Upsi",
					"DiMuonMassUpsi", 8, 8, 12);
				hDiMuonMass_Upsi->Sumw2(kTRUE);
				hDiMuonMass = new TH1D("DiMuonMass",
					"DiMuonMass", 100, 50, 200);
				hDiMuonMass->Sumw2(kTRUE);
				hDiMuonMass_Higgs = new TH1D("DiMuonMass_Higgs",
					"DiMuonMass_Higgs", 50, 100, 150);
				hDiMuonMass_Higgs->Sumw2(kTRUE);
				hDiMuonMass_Higgs2 = new TH1D("DiMuonMass_Higgs2",
					"DiMuonMass_Higgs", 25, 100, 150);
				hDiMuonMass_Higgs2->Sumw2(kTRUE);
				hDiMuonMass_Higgs25 = new TH1D("DiMuonMass_Higgs25",
					"DiMuonMass_Higgs", 20, 100, 150);
				hDiMuonMass_Higgs25->Sumw2(kTRUE);
				hDiMuonMass_Higgs3 = new TH1D("DiMuonMass_Higgs3",
					"DiMuonMass_Higgs", 17, 100, 150);
				hDiMuonMass_Higgs3->Sumw2(kTRUE);
				hDiMuonMass_Higgs35 = new TH1D("DiMuonMass_Higgs35",
					"DiMuonMass_Higgs", 14, 100, 150);
				hDiMuonMass_Higgs35->Sumw2(kTRUE);
				hDiMuonMassFSR = new TH1D("DiMuonMassFSR_full",
					"DiMuonMassFSR", 200, 0, 200);
				hDiMuonMassFSR->Sumw2(kTRUE);
				hDiMuonMassFSR_Higgs = new TH1D("DiMuonMassFSR_Higgs",
					"DiMuonMassFSR_Higgs", 1000, 100, 150);
				hDiMuonMassFSR_Higgs->Sumw2(kTRUE);
				hDiMuonMassFSR_Higgs2 = new TH1D("DiMuonMassFSR_Higgs2",
					"DiMuonMassFSR_Higgs", 25, 100, 150);
				hDiMuonMassFSR_Higgs2->Sumw2(kTRUE);
				hDiMuonMassFSR_Higgs25 = new TH1D("DiMuonMassFSR_Higgs25",
					"DiMuonMassFSR_Higgs", 20, 100, 150);
				hDiMuonMassFSR_Higgs25->Sumw2(kTRUE);
				hDiMuonMassFSR_Higgs3 = new TH1D("DiMuonMassFSR_Higgs3",
					"DiMuonMassFSR_Higgs", 17, 100, 150);
				hDiMuonMassFSR_Higgs3->Sumw2(kTRUE);
				hDiMuonMassFSR_Higgs35 = new TH1D("DiMuonMassFSR_Higgs35",
					"DiMuonMassFSR_Higgs", 14, 100, 150);
				hDiMuonMassFSR_Higgs35->Sumw2(kTRUE);
				hDiMuoneta = new TH1D("DiMuoneta", "DiMuoneta",
					50, -2.5, 2.5);
				hDiMuoneta->Sumw2(kTRUE);
				hDiMuondphi = new TH1D("DiMuondphi",
					"DiMuondphi", 18, -3.6, 3.6);
				hDiMuondphi->Sumw2(kTRUE);
				hDiMuonpteta = new TH2D("DiMuonpteta", "DiMuonpteta",
					70, -3.5, 3.5, 100, 0, 200);
				hDiMuonpteta->Sumw2(kTRUE);
				hDiMuonptM = new TH2D("DiMuonptM", "Dimuon pt vs M",
					200, 0, 200, 200, 0, 200);
				hDiMuonptM->Sumw2(kTRUE);
				hDiMuonptM_Jpsi = new TH2D("DiMuonptM_Jpsi", "Dimuon pt vs M",
					10, 1, 6, 200, 0, 200);
				hDiMuonptM_Jpsi->Sumw2(kTRUE);
				hDiMuonptM_Upsi = new TH2D("DiMuonptM_Upsi", "Dimuon pt vs M",
					8, 8, 12, 200, 0, 200);
				hDiMuonptM_Upsi->Sumw2(kTRUE);
				hDiMuonptM_Z = new TH2D("DiMuonptM_Z", "Dimuon pt vs M",
					100, 50, 150, 200, 0, 200);
				hDiMuonptM_Z->Sumw2(kTRUE);
				hDiMuonptM_Higgs = new TH2D("DiMuonptM_Higgs", "Dimuon pt vs M",
					50, 100, 150, 200, 0, 200);
				hDiMuonptM_Higgs->Sumw2(kTRUE);
				hMuonpt_l = new TH1D("Muonpt_leading", "Muonpt_leading",
					100, 0, 100);
				hMuonpt_l->Sumw2(kTRUE);
				hMuonpt_s = new TH1D("Muonpt_subleading", "Muonpt_subleading",
									 100, 0, 100);
				hMuonpt_s->Sumw2(kTRUE);
				hMuoneta = new TH1D("Muoneta", "Muoneta", 
					50, -2.5, 2.5);
				hMuoneta->Sumw2(kTRUE);
				hMuonphi = new TH1D("Muonphi", "Muonphi", 
					36, -3.6, 3.6);
				hMuonphi->Sumw2(kTRUE);
				hMuonptpt = new TH2D("Muonptpt", "Muonptpt",
					200, 0., 200., 200, 0., 200.);
				hMuonptpt->Sumw2(kTRUE);
				hMuonpteta = new TH2D("Muonpteta", "Muonpteta",
					50, -2.5, 2.5, 100, 0, 200);
				hMuonpteta->Sumw2(kTRUE);
				hFsrPhotonpt = new TH1D("FsrPhotonpt", "FsrPhotonpt",
					75, 0, 75);
				hFsrPhotonpt->Sumw2(kTRUE);
				hFsrPhotoneta = new TH1D("FsrPhotoneta", "FsrPhotoneta",
					90, -4.5, 4.5);
				hFsrPhotoneta->Sumw2(kTRUE);
				hFsrPhotonphi = new TH1D("FsrPhotonphi", "FsrPhotonphi",
					36, -3.6, 3.6);
				hFsrPhotonphi->Sumw2(kTRUE);
				hFsrPhotonpteta = new TH2D("FsrPhotonpteta", "FsrPhotonpteta",
					50, -4.5, 4.5, 75, 0, 75);
				hFsrPhotonpteta->Sumw2(kTRUE);
				hZ = new TH1D("hZ", "hZ",
					75, -15, 15);
				hZ->Sumw2(kTRUE);
				hYvsX = new TH2D("hYvsX", "hYvsX",
					500, -.4, .4, 500, -.4, .4);
				hYvsX->Sumw2(kTRUE);
				gDirectory->cd("../");
			}

			void scale(double factor)
			{
				hDiJetMass->Scale(factor);
				hDiJeteta->Scale(factor);
				hDiJetdeta->Scale(factor);
				hDiJetPt->Scale(factor);

				hDiMuonMass_full->Scale(factor);
				hDiMuonMass_Jpsi->Scale(factor);
				hDiMuonMass_Upsi->Scale(factor);
				hDiMuonMass->Scale(factor);
				hDiMuonMass_Higgs->Scale(factor);
				hDiMuonMass_Higgs2->Scale(factor);
				hDiMuonMass_Higgs25->Scale(factor);
				hDiMuonMass_Higgs3->Scale(factor);
				hDiMuonMass_Higgs35->Scale(factor);

				hDiMuonMassFSR->Scale(factor);
				hDiMuonMassFSR_Higgs->Scale(factor);
				hDiMuonMassFSR_Higgs2->Scale(factor);
				hDiMuonMassFSR_Higgs25->Scale(factor);
				hDiMuonMassFSR_Higgs3->Scale(factor);
				hDiMuonMassFSR_Higgs35->Scale(factor);

				hDiMuonpt->Scale(factor);
				hDiMuoneta->Scale(factor);
				hDiMuondphi->Scale(factor);

				hDiMuonpteta->Scale(factor);
				hDiMuonptM->Scale(factor);
				hDiMuonptM_Jpsi->Scale(factor);
				hDiMuonptM_Upsi->Scale(factor);
				hDiMuonptM_Z->Scale(factor);
				hDiMuonptM_Higgs->Scale(factor);

				hMuonpt_l->Scale(factor);
				hMuonpt_s->Scale(factor);
				hMuoneta->Scale(factor);
				hMuonphi->Scale(factor);
				hMuonptpt->Scale(factor);
				hMuonpteta->Scale(factor);

				hFsrPhotonpt->Scale(factor);
				hFsrPhotoneta->Scale(factor);
				hFsrPhotonphi->Scale(factor);
				hFsrPhotonpteta->Scale(factor);

				hZ->Scale(factor);
				hYvsX->Scale(factor);
//				hNlead->Scale(factor);
			}

			~DimuonSet() {}

			TString _postfix;
			TH1D *hDiJetMass;
			TH1D *hDiJeteta;
			TH1D *hDiJetdeta;
			TH1D *hDiJetPt;

			TH1D *hDiMuonMass_full;
			TH1D *hDiMuonMass_Jpsi;
			TH1D *hDiMuonMass_Upsi;
			TH1D *hDiMuonMass;
			TH1D *hDiMuonMass_Higgs;
			TH1D *hDiMuonMass_Higgs2;
			TH1D *hDiMuonMass_Higgs25;
			TH1D *hDiMuonMass_Higgs3;
			TH1D *hDiMuonMass_Higgs35;

			TH1D *hDiMuonMassFSR;
			TH1D *hDiMuonMassFSR_Higgs;
			TH1D *hDiMuonMassFSR_Higgs2;
			TH1D *hDiMuonMassFSR_Higgs25;
			TH1D *hDiMuonMassFSR_Higgs3;
			TH1D *hDiMuonMassFSR_Higgs35;

			TH1D *hDiMuonpt;
			TH1D *hDiMuoneta;
			TH1D *hDiMuondphi;

			TH2D *hDiMuonpteta;
			TH2D *hDiMuonptM;
			TH2D *hDiMuonptM_Jpsi;
			TH2D *hDiMuonptM_Upsi;
			TH2D *hDiMuonptM_Z;
			TH2D *hDiMuonptM_Higgs;

			TH1D *hMuonpt_l;
			TH1D *hMuonpt_s;
			TH1D *hMuoneta;
			TH1D *hMuonphi;
			TH2D *hMuonptpt;
			TH2D *hMuonpteta;		//

			TH1D *hFsrPhotonpt;		//
			TH1D *hFsrPhotoneta;	//
			TH1D *hFsrPhotonphi;	//
			TH2D *hFsrPhotonpteta;	//

			TH1D *hZ;				//
			TH2D *hYvsX;			//
//			TH1D *hNlead;

		};
	}
}
