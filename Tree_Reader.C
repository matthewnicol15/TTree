#include <cstdlib>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TH1.h>
#include <TF1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <vector>
#include <TLorentzVector.h>


void Tree_Reader(){



  gROOT->ProcessLine(".L ./Loader.C+");
  TFile *f = new TFile("skim4_Tree1.root");
  TTree *t1 = (TTree*)f->Get("skim4_Tree");

  vector<TLorentzVector> *v_p4=0;


  TLorentzVector *readbeam=NULL;
  TLorentzVector *readtarget=NULL;


  vector<TLorentzVector> *v_vertex=0;

  vector<double> *v_beta=0;

  Double_t start_time;
  vector<double> *energy=0;
  vector<double> *charge=0;
  vector<double> *PID=0;
  vector<double> *chi2PID=0;

  Int_t readchargetracks;
  Int_t readprotonno;
  Int_t readpipno;
  Int_t readpimno;
  Int_t readelno;
  Int_t readkaonpno;



  t1->SetBranchAddress("p4",&v_p4);

  t1->SetBranchAddress("vertex",&v_vertex);

  t1->SetBranchAddress("beta",&v_beta);

  t1->SetBranchAddress("beam",&readbeam);
  t1->SetBranchAddress("target",&readtarget);

  t1->SetBranchAddress("start_time",&start_time);
  t1->SetBranchAddress("energy",&energy);
  t1->SetBranchAddress("charge",&charge);
  t1->SetBranchAddress("PID",&PID);
  t1->SetBranchAddress("chi2PID",&chi2PID);
  t1->SetBranchAddress("chargetracks",&readchargetracks);
  t1->SetBranchAddress("protonno",&readprotonno);
  t1->SetBranchAddress("pipno",&readpipno);
  t1->SetBranchAddress("pimno",&readpimno);
  t1->SetBranchAddress("elno",&readelno);
  t1->SetBranchAddress("kaonpno",&readkaonpno);

  TFile fileOutput1("Reduction_Factor.root","recreate");

  //Creating histograms
  auto* hmiss_1=new TH1F("miss_1","MM(e' K^{+});P [GeV];Counts",400,-1,4);
  auto* hmiss_2=new TH1F("miss_2","MM^2(e' K^{+} p);P [GeV];Counts",400,-2,2);
  auto* hmiss_3=new TH1F("miss_3","MM(e' K^{+} K^{+});P [GeV];Counts",300,0,3);
  auto* hmiss_4=new TH1F("miss_4","MM^2(e' K^{+} K^{+} p);P [GeV];Counts",400,-2,2);
  auto* hmiss_5=new TH1F("miss_5","MM(e' K^{+} K^{+} K^{+});P [GeV];Counts",400,-1,3);

  auto* hmiss_1c=new TH1F("miss_1c","MM(e' K^{+}) after cuts on delta beta and momentum of K^{+}; [GeV];Counts",400,-1,4);
  auto* hmiss_2c=new TH1F("miss_2c","MM^2(e' K^{+} p) after cuts on delta beta and momentum of K^{+}; [GeV];Counts",400,-2,2);
  auto* hmiss_3c=new TH1F("miss_3c","MM(e' K^{+} K^{+}) after cuts on delta beta and momentum of K^{+}; [GeV];Counts",300,0,3);
  auto* hmiss_4c=new TH1F("miss_4c","MM^2(e' K^{+} K^{+} p) after cuts on delta beta and momentum of K^{+}; [GeV];Counts",400,-2,2);
  auto* hmiss_5c=new TH1F("miss_5c","MM(e' K^{+} K^{+} K^{+}) after cuts on delta beta and momentum of K^{+}; [GeV];Counts",400,-1,3);

  auto* hmiss_1t=new TH1F("miss_1t","MM(e' K^{+}) after cuts on delta beta and momentum of K^{+} and MM^2(e' K^{+} p) cut); [GeV];Counts",400,-1,4);
  auto* hmiss_2t=new TH1F("miss_2t","MM^2(e' K^{+} p) after cuts on delta beta and momentum of K^{+} and MM(e' K^{+}) cut); [GeV];Counts",400,-2,2);
  auto* hmiss_3t=new TH1F("miss_3t","MM(e' K^{+} K^{+}) after cuts on delta beta and momentum of K^{+} and MM^2(e' K^{+} K^{+} p) cut); [GeV];Counts",300,0,3);
  auto* hmiss_4t=new TH1F("miss_4t","MM^2(e' K^{+} K^{+} p) after cuts on delta beta and momentum of K^{+} and MM(e' K^{+} K^{+}) cut); [GeV];Counts",400,-2,2);

  auto* hK=new TH2F("K","MM(e' K^{+}) against MM^2(e' K^{+} p));MM(e' K^{+}) [GeV];MM^2(e' K^{+} p)) [GeV]",400,-1,4,400,-2,2);
  auto* hKK=new TH2F("KK","MM(e' K^{+} K^{+}) against MM^2(e' K^{+} K^{+} p));MM(e' K^{+} K^{+}) [GeV];MM^2(e' K^{+} K^{+} p)) [GeV]",400,-1,4,400,-2,2);


  auto* hl0_1=new TH1F("l0_1","M(p #pi^{-});P [GeV];Counts",200,0,5);

  auto* hl0_1c=new TH1F("l0_1c","M(p #pi^{-}) after cuts on delta beta and momentum of K^{+};P [GeV];Counts",200,0,5);

  auto* hkaonpno=new TH1F("kaonpno","number of K^{+} in an event",200,0,6);

  vector<TLorentzVector> v_kp;

  TLorentzVector miss1;
  TLorentzVector miss2;
  TLorentzVector miss3;
  TLorentzVector miss4;
  TLorentzVector miss5;
  TLorentzVector miss6;


  TLorentzVector l0;

  Double_t beta_tof_pip;
  Double_t P_pip;
  Double_t beta_calc_pip;
  Double_t delta_beta_pip;

  Double_t beta_tof_pim;
  Double_t P_pim;
  Double_t beta_calc_pim;
  Double_t delta_beta_pim;

  Double_t beta_tof_pr;
  Double_t P_pr;
  Double_t beta_calc_pr;
  Double_t delta_beta_pr;

  Double_t P_el;

  vector<Double_t> v_beta_tof_kp;
  Double_t beta_tof_kp;
  vector<Double_t> v_P_kp;
  Double_t P_kp;
  vector<Double_t> v_beta_calc_kp;
  Double_t beta_calc_kp;
  vector<Double_t> v_delta_beta_kp;
  Double_t delta_beta_kp;

  TLorentzVector el;
  TLorentzVector pip;
  TLorentzVector pim;
  TLorentzVector pr;
  TLorentzVector kp;

  Long64_t nentries = t1->GetEntries();

  for(Long64_t i=0; i<nentries;i++){
    t1->GetEntry(i);
    if (i % 1000 == 0){
      fprintf (stderr, "%lld\r", i/1000);
      fflush (stderr);
        }
        v_kp.clear();
        v_beta_tof_kp.clear();
        v_P_kp.clear();
        v_beta_calc_kp.clear();
        v_delta_beta_kp.clear();


      Int_t Nparticles = v_p4->size();
      for(Int_t j=0; j<Nparticles; j++){

        if(PID->at(j)==11){
          el.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
        }
        else if(PID->at(j)==211){
          pip.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
          beta_tof_pip = v_beta->at(j);
        }
        else if(PID->at(j)==-211){
          pim.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
          beta_tof_pim = v_beta->at(j);
        }
        else if(PID->at(j)==2212){
          pr.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
          beta_tof_pr = v_beta->at(j);
        }
        else if(PID->at(j)==321){

          kp.SetXYZM(v_p4->at(j).Px(),v_p4->at(j).Py(),v_p4->at(j).Pz(),v_p4->at(j).M());
          P_kp = sqrt((pow(v_p4->at(j).Px(),2))+(pow(v_p4->at(j).Py(),2))+(pow(v_p4->at(j).Pz(),2)));
          beta_tof_kp = v_beta->at(j);
          beta_calc_kp = P_kp/(sqrt((pow(P_kp,2))+(pow(kp.M(),2))));
          delta_beta_kp = beta_calc_kp-beta_tof_kp;

          v_kp.push_back(kp);
          v_beta_tof_kp.push_back(beta_tof_kp);
          v_P_kp.push_back(P_kp);
          v_beta_calc_kp.push_back(beta_calc_kp);
          v_delta_beta_kp.push_back(delta_beta_kp);
        }
      }

      P_pip = sqrt((pow(pip.Px(),2))+(pow(pip.Py(),2))+(pow(pip.Pz(),2)));
      beta_calc_pip = P_pip/(sqrt((pow(P_pip,2))+(pow(pip.M(),2))));
      delta_beta_pip = beta_calc_pip-beta_tof_pip;

      P_pim = sqrt((pow(pim.Px(),2))+(pow(pim.Py(),2))+(pow(pim.Pz(),2)));
      beta_calc_pim = P_pim/(sqrt((pow(P_pim,2))+(pow(pim.M(),2))));
      delta_beta_pim = beta_calc_pim-beta_tof_pim;

      P_pr = sqrt((pow(pr.Px(),2))+(pow(pr.Py(),2))+(pow(pr.Pz(),2)));
      beta_calc_pr = P_pr/(sqrt((pow(P_pr,2))+(pow(pr.M(),2))));
      delta_beta_pr = beta_calc_pr-beta_tof_pr;

      P_el = sqrt((pow(el.Px(),2))+(pow(el.Py(),2))+(pow(el.Pz(),2)));
      for(int k=0;k<readkaonpno;k++){

    }

      //K^+ channel
      miss1 = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - el - kp;
      miss2 = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - el - kp - pr;
      l0 = pr + pim;

      //K^+ K^+ channel
      miss3 = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - el - v_kp[0] - v_kp[1];
      miss4 = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - el - v_kp[0] - v_kp[1] - pr;
      // cascade = pr + pim;
      //K^+ K^+ K^+ channel

      miss5 = (TLorentzVector)*readbeam + (TLorentzVector)*readtarget - el - v_kp[0] - v_kp[1] - v_kp[2];

      hkaonpno->Fill(readkaonpno);

      if(readkaonpno==1 && readelno==1 && readprotonno==1){
          hl0_1->Fill(l0.M());
          hmiss_1->Fill(miss1.M());
          hmiss_2->Fill(miss2.M2());
          if(abs(v_delta_beta_kp[0])<0.01 && v_P_kp[0]>0.2 && abs(delta_beta_pr)<0.01 && P_pr>0.2){
            hl0_1c->Fill(l0.M());
            hmiss_1c->Fill(miss1.M());
            hmiss_2c->Fill(miss2.M2());
            hK->Fill(miss1.M(),miss2.M2());
            if(abs(miss2.M2()<0.1))hmiss_1t->Fill(miss1.M());
            if(miss1.M()>1 && miss1.M()<1.3)hmiss_2t->Fill(miss2.M2());
          }
      }

      if(readkaonpno==2 && readelno==1 && readprotonno==1){
          hmiss_3->Fill(miss3.M());
          hmiss_4->Fill(miss4.M2());
          if(abs(v_delta_beta_kp[0])<0.01 && v_P_kp[0]>0.2  && abs(v_delta_beta_kp[1])<0.01 && v_P_kp[1]>0.2 && abs(delta_beta_pr)<0.01 && P_pr>0.2){
            hmiss_3c->Fill(miss3.M());
            hmiss_4c->Fill(miss4.M2());
            hKK->Fill(miss3.M(),miss4.M2());
            if(abs(miss4.M2()<0.1))hmiss_3t->Fill(miss3.M());
            if(miss3.M()>1.2 && miss3.M()<1.4)hmiss_4t->Fill(miss4.M2());
          }
      }

      if(readkaonpno==3 && readelno==1){
          hmiss_5->Fill(miss5.M());
          if(abs(v_delta_beta_kp[0])<0.01 && v_P_kp[0]>0.2  && abs(v_delta_beta_kp[1])<0.01 && v_P_kp[1]>0.2 && abs(v_delta_beta_kp[2])<0.01 && v_P_kp[2]>0.2){
            hmiss_5c->Fill(miss5.M());
          }
      }


      }

fileOutput1.Write();


}
