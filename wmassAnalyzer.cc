//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Feb  9 09:23:08 2016 by ROOT version 5.34/32
// from TTree wmassAnalyzer/wmassAnalyzer
// found on file: /afs/cern.ch/work/m/mdunser/public/wmass/wmassAnalyzer_tree_6p2.root
//////////////////////////////////////////////////////////

#ifndef wmassAnalyzer_h
#define wmassAnalyzer_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TMath.h>
#include <TStopwatch.h>


#include <iostream>

using namespace std;

class wmassAnalyzer {
public :
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    // Declaration of leaf types
    // Int_t           run;
    // Int_t           lumi;
    // Int_t           evt;
    // Int_t           nvtx;
    Int_t           evtHasGoodVtx;
    // Int_t           Vtx_ndof;
    // Int_t           npu;
    // Double_t        scalePDF;
    Double_t        parton1_pdgId;
    Double_t        parton1_x;
    Double_t        parton2_pdgId;
    Double_t        parton2_x;
    Double_t        LHE_weight[466];
    // Double_t        pfmet;
    // Double_t        pfmet_phi;
    // Double_t        pfmet_sumEt;
    Double_t        tkmet;
    Double_t        tkmet_phi;
    Double_t        tkmet_sumEt;
    // Double_t        gentkmet;
    // Double_t        gentkmet_phi;
    // Double_t        gentkmet_sumEt;
    // Double_t        genpfmet;
    // Double_t        genpfmet_phi;
    // Double_t        genpfmet_sumEt;
    // Double_t        nopumet;
    // Double_t        nopumet_phi;
    // Double_t        nopumet_sumEt;
    // Double_t        pumet;
    // Double_t        pumet_phi;
    // Double_t        pumet_sumEt;
    // Double_t        pucmet;
    // Double_t        pucmet_phi;
    // Double_t        pucmet_sumEt;
    // Double_t        pfMetForRegression;
    // Double_t        pfMetForRegression_phi;
    // Double_t        pfMetForRegression_sumEt;
    Int_t           evtHasTrg;
    Int_t           njets;
    Int_t           evtWSel;
    Int_t           nMuons;
    Int_t           nTrgMuons;
    Int_t           noTrgMuonsLeadingPt;
    Double_t        W_pt;
    Double_t        W_phi;
    Double_t        W_mt;
    // Double_t        WGen_pt;
    // Double_t        WGen_rap;
    // Double_t        WGen_phi;
    // Double_t        WGen_mass;
    // Double_t        WGen_mt;
    // Double_t        NuGen_pt;
    // Double_t        NuGen_eta;
    // Double_t        NuGen_phi;
    // Double_t        NuGen_mass;
    Double_t        Mu_pt;
    Double_t        Mu_eta;
    Double_t        Mu_phi;
    Double_t        Mu_mass;
    Double_t        Mu_charge;
    // Double_t        MuRelIso;
    // Double_t        Mu_dxy;
    // Double_t        Mu_dz;
    Int_t           MuIsTightAndIso;
    // Int_t           MuIsTight;
    Double_t        pt_vis;
    Double_t        phi_vis;
    // Double_t        MuGen_pt;
    // Double_t        MuGen_eta;
    // Double_t        MuGen_phi;
    // Double_t        MuGen_mass;
    // Double_t        MuGen_charge;
    // Double_t        MuGenStatus1_pt;
    // Double_t        MuGenStatus1_eta;
    // Double_t        MuGenStatus1_phi;
    // Double_t        MuGenStatus1_mass;
    // Double_t        MuGenStatus1_charge;
    // Double_t        MuDRGenP;
    // Double_t        FSRWeight;
    // Double_t        MuCovMatrix[9];
    // Double_t        Jet_leading_pt;
    // Double_t        Jet_leading_eta;
    // Double_t        Jet_leading_phi;
    // Double_t        genWLept;

    // List of branches
    // TBranch        *b_run;   //!
    // TBranch        *b_lumi;   //!
    // TBranch        *b_evt;   //!
    // TBranch        *b_nvtx;   //!
    TBranch        *b_evtHasGoodVtx;   //!
    // TBranch        *b_Vtx_ndof;   //!
    // TBranch        *b_npu;   //!
    // TBranch        *b_scalePDF;   //!
    TBranch        *b_parton1_pdgId;   //!
    TBranch        *b_parton1_x;   //!
    TBranch        *b_parton2_pdgId;   //!
    TBranch        *b_parton2_x;   //!
    TBranch        *b_LHE_weight;   //!
    // TBranch        *b_pfmet;   //!
    // TBranch        *b_pfmet_phi;   //!
    // TBranch        *b_pfmet_sumEt;   //!
    TBranch        *b_tkmet;   //!
    TBranch        *b_tkmet_phi;   //!
    TBranch        *b_tkmet_sumEt;   //!
    // TBranch        *b_gentkmet;   //!
    // TBranch        *b_gentkmet_phi;   //!
    // TBranch        *b_gentkmet_sumEt;   //!
    // TBranch        *b_genpfmet;   //!
    // TBranch        *b_genpfmet_phi;   //!
    // TBranch        *b_genpfmet_sumEt;   //!
    // TBranch        *b_nopumet;   //!
    // TBranch        *b_nopumet_phi;   //!
    // TBranch        *b_nopumet_sumEt;   //!
    // TBranch        *b_pumet;   //!
    // TBranch        *b_pumet_phi;   //!
    // TBranch        *b_pumet_sumEt;   //!
    // TBranch        *b_pucmet;   //!
    // TBranch        *b_pucmet_phi;   //!
    // TBranch        *b_pucmet_sumEt;   //!
    // TBranch        *b_pfMetForRegression;   //!
    // TBranch        *b_pfMetForRegression_phi;   //!
    // TBranch        *b_pfMetForRegression_sumEt;   //!
    TBranch        *b_evtHasTrg;   //!
    TBranch        *b_njets;   //!
    TBranch        *b_evtWSel;   //!
    TBranch        *b_nMuons;   //!
    TBranch        *b_nTrgMuons;   //!
    TBranch        *b_noTrgMuonsLeadingPt;   //!
    TBranch        *b_W_pt;   //!
    TBranch        *b_W_phi;   //!
    TBranch        *b_W_mt;   //!
    // TBranch        *b_WGen_pt;   //!
    // TBranch        *b_WGen_rap;   //!
    // TBranch        *b_WGen_phi;   //!
    // TBranch        *b_WGen_mass;   //!
    // TBranch        *b_WGen_mt;   //!
    // TBranch        *b_NuGen_pt;   //!
    // TBranch        *b_NuGen_eta;   //!
    // TBranch        *b_NuGen_phi;   //!
    // TBranch        *b_NuGen_mass;   //!
    TBranch        *b_Mu_pt;   //!
    TBranch        *b_Mu_eta;   //!
    TBranch        *b_Mu_phi;   //!
    TBranch        *b_Mu_mass;   //!
    TBranch        *b_Mu_charge;   //!
    // TBranch        *b_MuRelIso;   //!
    // TBranch        *b_Mu_dxy;   //!
    // TBranch        *b_Mu_dz;   //!
    TBranch        *b_MuIsTightAndIso;   //!
    // TBranch        *b_MuIsTight;   //!
    TBranch        *b_pt_vis;   //!
    TBranch        *b_phi_vis;   //!
    // TBranch        *b_MuGen_pt;   //!
    // TBranch        *b_MuGen_eta;   //!
    // TBranch        *b_MuGen_phi;   //!
    // TBranch        *b_MuGen_mass;   //!
    // TBranch        *b_MuGen_charge;   //!
    // TBranch        *b_MuGenStatus1_pt;   //!
    // TBranch        *b_MuGenStatus1_eta;   //!
    // TBranch        *b_MuGenStatus1_phi;   //!
    // TBranch        *b_MuGenStatus1_mass;   //!
    // TBranch        *b_MuGenStatus1_charge;   //!
    // TBranch        *b_MuDRGenP;   //!
    // TBranch        *b_FSRWeight;   //!
    // TBranch        *b_MuCovMatrix;   //!
    // TBranch        *b_Jet_leading_pt;   //!
    // TBranch        *b_Jet_leading_eta;   //!
    // TBranch        *b_Jet_leading_phi;   //!
    // TBranch        *b_genWLept;   //!

    wmassAnalyzer(TTree *tree=0);
    virtual ~wmassAnalyzer();
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree *tree);
    virtual void     Loop();
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);


    // my functions hooray
    bool fIsData;
    Long64_t fNEntries;
    int fMaxEvents;
    virtual void RunJob(TString, bool);
    virtual void Begin(TFile *);
    virtual void End(TFile *);
    virtual bool IsGoodEvent();
    virtual void fillMTandETAHistograms();
    virtual void setMaxEvents(int me){fMaxEvents = me;};
    TH1F * nominal_mt       [7];
    TH1F * nominal_muEta    [7];
    TH1F * nominal_pos_mt   [7];
    TH1F * nominal_pos_muEta[7];
    TH1F * nominal_neg_mt   [7];
    TH1F * nominal_neg_muEta[7];
    TH1F * ct10_mt   [7][52]; // they have indices 309-361
    TH1F * ct10_muEta[7][52];
    TH1F * ct10_pos_mt   [7][52];
    TH1F * ct10_pos_muEta[7][52];
    TH1F * ct10_neg_mt   [7][52];
    TH1F * ct10_neg_muEta[7][52];

    

};

// #endif
// 
// #ifdef wmassAnalyzer_cxx
//
wmassAnalyzer::wmassAnalyzer(TTree *tree) : fChain(0) 
{
    fMaxEvents = -1;
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
    if (tree == 0) {
       TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/afs/cern.ch/work/m/mdunser/public/wmass/WTreeProducer_tree_6p2.root");
       if (!f || !f->IsOpen()) {
          f = new TFile("/afs/cern.ch/work/m/mdunser/public/wmass/WTreeProducer_tree_6p2.root");
       }
       f->GetObject("WTreeProducer",tree);

    }
    Init(tree);
}

wmassAnalyzer::~wmassAnalyzer()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

Int_t wmassAnalyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}
Long64_t wmassAnalyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
       fCurrent = fChain->GetTreeNumber();
       Notify();
    }
    return centry;
}

void wmassAnalyzer::Init(TTree *tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    // fChain->SetBranchAddress("run", &run, &b_run);
    // fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
    // fChain->SetBranchAddress("evt", &evt, &b_evt);
    // fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
    fChain->SetBranchAddress("evtHasGoodVtx", &evtHasGoodVtx, &b_evtHasGoodVtx);
    // fChain->SetBranchAddress("Vtx_ndof", &Vtx_ndof, &b_Vtx_ndof);
    // fChain->SetBranchAddress("npu", &npu, &b_npu);
    // fChain->SetBranchAddress("scalePDF", &scalePDF, &b_scalePDF);
    fChain->SetBranchAddress("parton1_pdgId", &parton1_pdgId, &b_parton1_pdgId);
    fChain->SetBranchAddress("parton1_x", &parton1_x, &b_parton1_x);
    fChain->SetBranchAddress("parton2_pdgId", &parton2_pdgId, &b_parton2_pdgId);
    fChain->SetBranchAddress("parton2_x", &parton2_x, &b_parton2_x);
    fChain->SetBranchAddress("LHE_weight", LHE_weight, &b_LHE_weight);
    // fChain->SetBranchAddress("pfmet", &pfmet, &b_pfmet);
    // fChain->SetBranchAddress("pfmet_phi", &pfmet_phi, &b_pfmet_phi);
    // fChain->SetBranchAddress("pfmet_sumEt", &pfmet_sumEt, &b_pfmet_sumEt);
    fChain->SetBranchAddress("tkmet", &tkmet, &b_tkmet);
    fChain->SetBranchAddress("tkmet_phi", &tkmet_phi, &b_tkmet_phi);
    fChain->SetBranchAddress("tkmet_sumEt", &tkmet_sumEt, &b_tkmet_sumEt);
    // fChain->SetBranchAddress("gentkmet", &gentkmet, &b_gentkmet);
    // fChain->SetBranchAddress("gentkmet_phi", &gentkmet_phi, &b_gentkmet_phi);
    // fChain->SetBranchAddress("gentkmet_sumEt", &gentkmet_sumEt, &b_gentkmet_sumEt);
    // fChain->SetBranchAddress("genpfmet", &genpfmet, &b_genpfmet);
    // fChain->SetBranchAddress("genpfmet_phi", &genpfmet_phi, &b_genpfmet_phi);
    // fChain->SetBranchAddress("genpfmet_sumEt", &genpfmet_sumEt, &b_genpfmet_sumEt);
    // fChain->SetBranchAddress("nopumet", &nopumet, &b_nopumet);
    // fChain->SetBranchAddress("nopumet_phi", &nopumet_phi, &b_nopumet_phi);
    // fChain->SetBranchAddress("nopumet_sumEt", &nopumet_sumEt, &b_nopumet_sumEt);
    // fChain->SetBranchAddress("pumet", &pumet, &b_pumet);
    // fChain->SetBranchAddress("pumet_phi", &pumet_phi, &b_pumet_phi);
    // fChain->SetBranchAddress("pumet_sumEt", &pumet_sumEt, &b_pumet_sumEt);
    // fChain->SetBranchAddress("pucmet", &pucmet, &b_pucmet);
    // fChain->SetBranchAddress("pucmet_phi", &pucmet_phi, &b_pucmet_phi);
    // fChain->SetBranchAddress("pucmet_sumEt", &pucmet_sumEt, &b_pucmet_sumEt);
    // fChain->SetBranchAddress("pfMetForRegression", &pfMetForRegression, &b_pfMetForRegression);
    // fChain->SetBranchAddress("pfMetForRegression_phi", &pfMetForRegression_phi, &b_pfMetForRegression_phi);
    // fChain->SetBranchAddress("pfMetForRegression_sumEt", &pfMetForRegression_sumEt, &b_pfMetForRegression_sumEt);
    fChain->SetBranchAddress("evtHasTrg", &evtHasTrg, &b_evtHasTrg);
    fChain->SetBranchAddress("njets", &njets, &b_njets);
    fChain->SetBranchAddress("evtWSel", &evtWSel, &b_evtWSel);
    fChain->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
    fChain->SetBranchAddress("nTrgMuons", &nTrgMuons, &b_nTrgMuons);
    fChain->SetBranchAddress("noTrgMuonsLeadingPt", &noTrgMuonsLeadingPt, &b_noTrgMuonsLeadingPt);
    fChain->SetBranchAddress("W_pt", &W_pt, &b_W_pt);
    fChain->SetBranchAddress("W_phi", &W_phi, &b_W_phi);
    fChain->SetBranchAddress("W_mt", &W_mt, &b_W_mt);
    // fChain->SetBranchAddress("WGen_pt", &WGen_pt, &b_WGen_pt);
    // fChain->SetBranchAddress("WGen_rap", &WGen_rap, &b_WGen_rap);
    // fChain->SetBranchAddress("WGen_phi", &WGen_phi, &b_WGen_phi);
    // fChain->SetBranchAddress("WGen_mass", &WGen_mass, &b_WGen_mass);
    // fChain->SetBranchAddress("WGen_mt", &WGen_mt, &b_WGen_mt);
    // fChain->SetBranchAddress("NuGen_pt", &NuGen_pt, &b_NuGen_pt);
    // fChain->SetBranchAddress("NuGen_eta", &NuGen_eta, &b_NuGen_eta);
    // fChain->SetBranchAddress("NuGen_phi", &NuGen_phi, &b_NuGen_phi);
    // fChain->SetBranchAddress("NuGen_mass", &NuGen_mass, &b_NuGen_mass);
    fChain->SetBranchAddress("Mu_pt", &Mu_pt, &b_Mu_pt);
    fChain->SetBranchAddress("Mu_eta", &Mu_eta, &b_Mu_eta);
    fChain->SetBranchAddress("Mu_phi", &Mu_phi, &b_Mu_phi);
    fChain->SetBranchAddress("Mu_mass", &Mu_mass, &b_Mu_mass);
    fChain->SetBranchAddress("Mu_charge", &Mu_charge, &b_Mu_charge);
    // fChain->SetBranchAddress("MuRelIso", &MuRelIso, &b_MuRelIso);
    // fChain->SetBranchAddress("Mu_dxy", &Mu_dxy, &b_Mu_dxy);
    // fChain->SetBranchAddress("Mu_dz", &Mu_dz, &b_Mu_dz);
    fChain->SetBranchAddress("MuIsTightAndIso", &MuIsTightAndIso, &b_MuIsTightAndIso);
    // fChain->SetBranchAddress("MuIsTight", &MuIsTight, &b_MuIsTight);
    fChain->SetBranchAddress("pt_vis", &pt_vis, &b_pt_vis);
    fChain->SetBranchAddress("phi_vis", &phi_vis, &b_phi_vis);
    // fChain->SetBranchAddress("MuGen_pt", &MuGen_pt, &b_MuGen_pt);
    // fChain->SetBranchAddress("MuGen_eta", &MuGen_eta, &b_MuGen_eta);
    // fChain->SetBranchAddress("MuGen_phi", &MuGen_phi, &b_MuGen_phi);
    // fChain->SetBranchAddress("MuGen_mass", &MuGen_mass, &b_MuGen_mass);
    // fChain->SetBranchAddress("MuGen_charge", &MuGen_charge, &b_MuGen_charge);
    // fChain->SetBranchAddress("MuGenStatus1_pt", &MuGenStatus1_pt, &b_MuGenStatus1_pt);
    // fChain->SetBranchAddress("MuGenStatus1_eta", &MuGenStatus1_eta, &b_MuGenStatus1_eta);
    // fChain->SetBranchAddress("MuGenStatus1_phi", &MuGenStatus1_phi, &b_MuGenStatus1_phi);
    // fChain->SetBranchAddress("MuGenStatus1_mass", &MuGenStatus1_mass, &b_MuGenStatus1_mass);
    // fChain->SetBranchAddress("MuGenStatus1_charge", &MuGenStatus1_charge, &b_MuGenStatus1_charge);
    // fChain->SetBranchAddress("MuDRGenP", &MuDRGenP, &b_MuDRGenP);
    // fChain->SetBranchAddress("FSRWeight", &FSRWeight, &b_FSRWeight);
    // fChain->SetBranchAddress("MuCovMatrix", MuCovMatrix, &b_MuCovMatrix);
    // fChain->SetBranchAddress("Jet_leading_pt", &Jet_leading_pt, &b_Jet_leading_pt);
    // fChain->SetBranchAddress("Jet_leading_eta", &Jet_leading_eta, &b_Jet_leading_eta);
    // fChain->SetBranchAddress("Jet_leading_phi", &Jet_leading_phi, &b_Jet_leading_phi);
    // fChain->SetBranchAddress("genWLept", &genWLept, &b_genWLept);
    Notify();
}

Bool_t wmassAnalyzer::Notify()
{
    return kTRUE;
}

void wmassAnalyzer::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}


void wmassAnalyzer::RunJob(TString filename, bool isData){
   fIsData = isData;

   TFile *file = TFile::Open(filename, "recreate");
   //do the analysis
   Begin(file);
   Loop();
   End(file);
}
  
void wmassAnalyzer::Begin(TFile *file){ // book the histograms and all
    file->cd();
    int nbins_mt = 250, nbins_muEta = 250;
    float mt_min = 30., mt_max = 200.;
    float muEta_min = -2.5, muEta_max = 2.5;

    // make histograms for mt and muEta
    for (int m = 0; m<7; m++){
        nominal_pos_mt   [m] = new TH1F(Form("nominal_pos_mt_m%d"   ,m), Form("nominal_pos_mt_m%d"   ,m), nbins_mt   , mt_min   , mt_max   );
        nominal_pos_muEta[m] = new TH1F(Form("nominal_pos_muEta_m%d",m), Form("nominal_pos_muEta_m%d",m), nbins_muEta, muEta_min, muEta_max);
        nominal_neg_mt   [m] = new TH1F(Form("nominal_neg_mt_m%d"   ,m), Form("nominal_neg_mt_m%d"   ,m), nbins_mt   , mt_min   , mt_max   );
        nominal_neg_muEta[m] = new TH1F(Form("nominal_neg_muEta_m%d",m), Form("nominal_neg_muEta_m%d",m), nbins_muEta, muEta_min, muEta_max);
        for (int p = 0; p<52; p++){
            ct10_pos_mt   [m][p] = new TH1F(Form("ct10_pos_mt_m%d_p%d",m,p)   , Form("ct10_pos_mt_m%d_p%d",m,p)   , nbins_mt   , mt_min   , mt_max   );
            ct10_pos_muEta[m][p] = new TH1F(Form("ct10_pos_muEta_m%d_p%d",m,p), Form("ct10_pos_muEta_m%d_p%d",m,p), nbins_muEta, muEta_min, muEta_max);
            ct10_neg_mt   [m][p] = new TH1F(Form("ct10_neg_mt_m%d_p%d",m,p)   , Form("ct10_neg_mt_m%d_p%d",m,p)   , nbins_mt   , mt_min   , mt_max   );
            ct10_neg_muEta[m][p] = new TH1F(Form("ct10_neg_muEta_m%d_p%d",m,p), Form("ct10_neg_muEta_m%d_p%d",m,p), nbins_muEta, muEta_min, muEta_max);
        }
    }
}

void wmassAnalyzer::End(TFile *file){
    file->cd();
    for (int m = 0; m<7; m++){
        nominal_mt       [m]= (TH1F*) nominal_pos_mt   [m]->Clone(Form("nominal_all_mt_m%d",m))   ; nominal_mt   [m]->Add(nominal_neg_mt   [m]); 
        nominal_muEta    [m]= (TH1F*) nominal_pos_muEta[m]->Clone(Form("nominal_all_muEta_m%d",m)); nominal_muEta[m]->Add(nominal_neg_muEta[m]);
        nominal_mt       [m]->Write();
        nominal_muEta    [m]->Write();
        nominal_pos_mt   [m]->Write(); 
        nominal_pos_muEta[m]->Write();
        nominal_neg_mt   [m]->Write(); 
        nominal_neg_muEta[m]->Write();
        for (int p = 0; p<52; p++){
            ct10_mt   [m][p] =(TH1F*) ct10_pos_mt   [m][p]->Clone(Form("ct10_all_mt_m%d_p%d",m,p))   ; ct10_mt   [m][p] ->Add(ct10_neg_mt   [m][p]);
            ct10_muEta[m][p] =(TH1F*) ct10_pos_muEta[m][p]->Clone(Form("ct10_all_muEta_m%d_p%d",m,p)); ct10_muEta[m][p] ->Add(ct10_neg_muEta[m][p]);
            ct10_mt       [m][p] -> Write();
            ct10_muEta    [m][p] -> Write();
            ct10_pos_mt   [m][p] -> Write();
            ct10_pos_muEta[m][p] -> Write();
            ct10_neg_mt   [m][p] -> Write();
            ct10_neg_muEta[m][p] -> Write();
        }
    }
    //file->Write();

    gROOT->SetBatch();
    TCanvas * c = new TCanvas("sanity_check_pdfs", "sanity_check_pdfs", 800, 600);
    c->cd();
    nominal_muEta[0] ->Draw("");
    for (int i = 2; i<10; i++){
        ct10_muEta[0][i] ->SetLineColor(i); ct10_muEta[0][i]->Draw("same");
    }
    c->Write();
    c->Clear();
    c->SetName("sanity_check_masses");
    for (int i = 0; i<7; i++){
        nominal_mt[i] ->SetLineColor(i+1); nominal_mt[i]->Draw("same");
    }
    c->Write();
    
    file->Close();
}

bool wmassAnalyzer::IsGoodEvent()
{
    // very basic event selection
    if(!evtHasGoodVtx)   return false;
    if(!evtHasTrg)       return false;
    if(tkmet  < 25.)     return false;
    if(Mu_pt  < 25.)     return false;
    if(W_mt   < 30.)     return false;
    if(!MuIsTightAndIso) return false;
    if(W_pt > 20.)       return false;
    return true;
}

void wmassAnalyzer::Loop()
{
    TStopwatch stopWatch;
    stopWatch.Start();

    //     This is the loop skeleton where:
    //    jentry is the global entry number in the chain
    //    ientry is the entry number in the current Tree
    //  Note that the argument to GetEntry must be:
    //    jentry for TChain::GetEntry
    //    ientry for TTree::GetEntry and TBranch::GetEntry

    if (fMaxEvents > 0) cout << " running on only " << fMaxEvents << " events " << endl;
    if (fChain == 0) return;

    fNEntries = fChain->GetEntriesFast();

    Long64_t nbytes = 0, nb = 0;
    Long64_t loopedEvents = 0;
    for (Long64_t jentry=0; jentry<fNEntries;jentry++) {
        if ( jentry > 0 && jentry%((int) pow(10, ((int)TMath::Log10(jentry))) )== 0) cout << ">>> Processing event # " << jentry << endl;
        if (fMaxEvents > 0 && jentry >= fMaxEvents) break;
        loopedEvents++;
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        if (!IsGoodEvent()) continue;

        // do relevant things here
        fillMTandETAHistograms();

    } // end loop on events
    double elapsed = stopWatch.RealTime();
    cout << "it took " << elapsed << " seconds to loop on " << loopedEvents << " events" << endl;
    cout << loopedEvents/elapsed << " ev/s" << endl;

} // end Loop() function

void wmassAnalyzer::fillMTandETAHistograms(){

    // loop on all ct10 variations
    for (int m = 0; m<7 ; m++){
    Mu_charge > 0 ? nominal_pos_mt   [m]->Fill(W_mt  , LHE_weight[202+m]) : nominal_neg_mt   [m]->Fill(W_mt  , LHE_weight[202+m]);
    Mu_charge > 0 ? nominal_pos_muEta[m]->Fill(Mu_eta, LHE_weight[202+m]) : nominal_neg_muEta[m]->Fill(Mu_eta, LHE_weight[202+m]);
    for (int p = 0; p<52; p++){
        if(Mu_charge > 0){
            ct10_pos_mt   [m][p] -> Fill(W_mt  , LHE_weight[309+p]*LHE_weight[202+m]) ;
            ct10_pos_muEta[m][p] -> Fill(Mu_eta, LHE_weight[309+p]*LHE_weight[202+m]) ;
        }
        else{
            ct10_neg_mt   [m][p] -> Fill(W_mt  , LHE_weight[309+p]*LHE_weight[202+m]) ;
            ct10_neg_muEta[m][p] -> Fill(Mu_eta, LHE_weight[309+p]*LHE_weight[202+m]) ;
        }
    
    } // end for loop over pdf variations
    } // end for loop over masses
} // end fillMTandETAHistograms

#endif // #ifdef wmassAnalyzer_cxx
