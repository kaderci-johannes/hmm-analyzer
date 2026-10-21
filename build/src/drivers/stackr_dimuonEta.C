void stackr_dimuonEta(const char *E_year="L = 7.98 fb^{-1}"){

    gStyle->SetOptStat(0);
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadRightMargin(0.05);

//  Histogram Ranges
    pair<double,double> hor = {-2.4, 2.4};  // Horizontal axis - min. and max.
    pair<double,double> rat = {0.9, 1.1};   // Ratio histogram vertical axis - min. and max.

//  double scaleFactor = 50 / 923.6;

    // --- Input files ---
    const char *dir = "./HTCondor/outputs/Sums/";

    TFile *fData = TFile::Open((TString(dir) + "data_2024_1_Muon.root").Data(),"r");
    TFile *fDY   = TFile::Open((TString(dir) + "mc_2024_1_DY.root").Data(),"r");
    TFile *fEWK  = TFile::Open((TString(dir) + "mc_2024_1_EWK.root").Data(),"r");
    TFile *fVV   = TFile::Open((TString(dir) + "mc_2024_1_VV.root").Data(),"r");
    TFile *fTT   = TFile::Open((TString(dir) + "mc_2024_1_TT.root").Data(),"r");
    TFile *fGGH  = TFile::Open((TString(dir) + "mc_2024_1_ggH.root").Data(),"r");
//  TFile *fVBF  = TFile::Open("","r");

    if(!fDY || !fEWK || !fTT || !fVV ||/*!fVBF ||*/ !fGGH || !fData){
        std::cerr << "ERROR: missing input files\n";
        return;
    }

    auto getHist = [&](TFile* f){
        TH1F* h = (TH1F*)f->Get("AllJets/DiMuoneta");
        return h ? (TH1F*)h->Clone() : nullptr;
    };

    TH1F *hDY   = getHist(fDY);
    TH1F *hEWK  = getHist(fEWK);
    TH1F *hTT   = getHist(fTT);
    TH1F *hVV   = getHist(fVV);
    TH1F *hGGH  = getHist(fGGH);
    TH1F *hData = getHist(fData);
//  TH1F *hVBF  = getHist(fVBF);

//  hTT->Scale(scaleFactor);

    // --- Bin-width normalization ---
    auto normBW = [](TH1F* h){
        for(int i=1;i<=h->GetNbinsX();++i){
            double w = h->GetBinWidth(i);
            h->SetBinContent(i, h->GetBinContent(i)/w);
            h->SetBinError(i,   h->GetBinError(i)/w);
        }
    };

    for(auto h : {hDY,hEWK,hTT,hVV,hGGH,hData/*,hVBF*/}) normBW(h);

    // --- Normalize MC to Data ---
    TH1F *hMC = (TH1F*)hDY->Clone("hMC");
    hMC->Add(hEWK); hMC->Add(hTT); hMC->Add(hVV);
    hMC->Add(hGGH);/* hMC->Add(hVBF);*/

    double mcScale = hData->Integral() / hMC->Integral();
    for(auto h : {hDY,hTT,hVV,hEWK,hGGH/*,hVBF*/}) h->Scale(mcScale);

    // ===================== DATA BLINDING (NUMERICAL) =====================
    const double blindLow  = 115.0;
    const double blindHigh = 135.0;

    TH1F *hDataBlind = (TH1F*)hData->Clone("hDataBlind");

    for(int i=1;i<=hDataBlind->GetNbinsX();++i){
        double x = hDataBlind->GetBinCenter(i);
        if(x > blindLow && x < blindHigh){
            hDataBlind->SetBinContent(i, 0.0);
            hDataBlind->SetBinError(i,   0.0);
        }
    }

    // ===================== Styles =====================
    hDY ->SetFillColor(kBlue-7);
    hEWK->SetFillColor(kTeal-7);
    hTT ->SetFillColor(kRed-7);
    hVV ->SetFillColor(kOrange-3);

    for(auto h : {hDY,hEWK,hTT,hVV}){
        h->SetFillStyle(1001);
        h->SetLineColor(kBlack);
        h->SetLineWidth(1);
    }

//  hVBF->SetLineColor(kGreen+2);
//  hVBF->SetLineWidth(3);
//  hVBF->SetFillStyle(0);

    hGGH->SetLineColor(kMagenta+2);
    hGGH->SetLineWidth(3);
    hGGH->SetFillStyle(0);

    hDataBlind->SetMarkerStyle(20);
    hDataBlind->SetMarkerSize(1.1);
    hDataBlind->SetLineWidth(2);

    // ===================== Canvas & pads =====================
    TCanvas *c = new TCanvas("c","Dimuon mass",900,800);

    TPad *pad1 = new TPad("pad1","",0,0.30,1,1);
    TPad *pad2 = new TPad("pad2","",0,0,1,0.30);

    pad1->SetBottomMargin(0.02);
    pad1->SetTopMargin(0.08);
    pad2->SetTopMargin(0.02);
    pad2->SetBottomMargin(0.35);

    pad1->Draw();
    pad2->Draw();

    // ===================== Upper pad =====================
    pad1->cd();
    pad1->SetLogy();

    THStack *hs = new THStack("hs","");

    hs->Add(hEWK);
    hs->Add(hVV);
    hs->Add(hTT);
    hs->Add(hDY);

    hs->Draw("HIST");
    hs->GetXaxis()->SetRangeUser(hor.first,hor.second);
    hs->GetYaxis()->SetTitle("Events");
    hs->SetMinimum(1e-2);

//  hVBF->Draw("HIST SAME");
    hGGH->Draw("HIST SAME");
    hDataBlind->Draw("E SAME");

    TLegend *leg = new TLegend(0.82,0.73,1.00,0.90);
    leg->AddEntry(hDataBlind,"Data","lep");
    leg->AddEntry(hDY,"DY","f");
    leg->AddEntry(hTT,"t#bar{t}","f");
    leg->AddEntry(hVV,"VV","f");
    leg->AddEntry(hEWK,"EWK","f");
    leg->AddEntry(hGGH,"ggH","l");
    //  leg->AddEntry(hVBF,"VBF","l");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();

    // --- CMS label ---
    char cms_pre[64];
    sprintf(cms_pre,"#bf{CMS} #it{Preliminary} #sqrt{s} = 13.6 TeV, %s",E_year);
    cout<<cms_pre<<endl;
    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.05);
    latex.DrawLatex(0.15,0.92,cms_pre);

    pad1->RedrawAxis();

    // ===================== Ratio pad =====================
    pad2->cd();

    TH1F *hRatio = (TH1F*)hData->Clone("hRatio");
    TH1F* hs_ = (TH1F*) ((TH1F*)hs->GetStack()->Last())->Clone("hsum");
    hRatio->Divide(hs_);

    for(int i=1;i<=hRatio->GetNbinsX();++i){
        double x = hRatio->GetBinCenter(i);
        if(x > blindLow && x < blindHigh){
            hRatio->SetBinContent(i, 0.0);
            hRatio->SetBinError(i,   0.0);
        }
    }

    hRatio->SetTitle("");
    hRatio->GetYaxis()->SetTitle("Data / MC");
    hRatio->GetYaxis()->SetRangeUser(rat.first,rat.second);
    hRatio->GetYaxis()->SetTitleSize(0.10);
    hRatio->GetYaxis()->SetLabelSize(0.10);
    hRatio->GetYaxis()->SetTitleOffset(0.45);

    hRatio->GetXaxis()->SetTitle("#eta_{#mu#mu}");
    hRatio->GetXaxis()->SetTitleSize(0.12);
    hRatio->GetXaxis()->SetLabelSize(0.10);
    hRatio->GetXaxis()->SetRangeUser(hor.first,hor.second);

    hRatio->Draw("E");

    // --- Vertical blinding lines ONLY in ratio pad ---
    TLine *lBlindLow  = new TLine(blindLow,  0.5, blindLow,  1.5);
    TLine *lBlindHigh = new TLine(blindHigh, 0.5, blindHigh, 1.5);

    lBlindLow ->SetLineStyle(2);
    lBlindHigh->SetLineStyle(2);
    lBlindLow ->SetLineColor(kGray+2);
    lBlindHigh->SetLineColor(kGray+2);
    lBlindLow ->SetLineWidth(2);
    lBlindHigh->SetLineWidth(2);

    lBlindLow ->Draw("SAME");
    lBlindHigh->Draw("SAME");

    TLine *l = new TLine(hor.first,1,hor.second,1);
    l->SetLineStyle(2);
    l->SetLineWidth(2);
    l->Draw();

    pad2->RedrawAxis();

    c->SaveAs("DimuonEta_log_withRatio_binWidth_ratioLines.pdf");
}
