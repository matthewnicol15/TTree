#include <cstdlib>
#include <iostream>
#include <chrono>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TDatabasePDG.h>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include "clas12reader.h"
#include <vector>
using namespace clas12;


void Tree_Maker_New(){
  gROOT->ProcessLine(".L ./Loader.C+"); // Uses Loader.C file, make sure Loader.C is in this file path

  // Filepath to input files, use *.hipo to analyse all hipo files in a directory
  TString inputFile("/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v0/dst/train/skim4/*.hipo");

  // Creating a chain for the data from different files
  TChain fake("hipo");
  fake.Add(inputFile.Data());

  // Shortcut to a list of all the input file names
  auto files=fake.GetListOfFiles();

  // Create root file to save TTree in
  TFile f("/volatile/clas12/matthewn/skim4_Pass1_Tree_280520_01.root","recreate");
  // Creating TTree object
  TTree skim4_Tree_280520_01("skim4_Tree_280520_01","it's a tree!");

  // Access the PDG database to get information usind PID (e.g db->GetParticle(211)->Mass() for pi^{+} mass)
  auto db=TDatabasePDG::Instance();

  // Information to save to the tree
  // Any information specific to an event
  Int_t eventno; // Records the event number
  Int_t runno; // Records the run number
  Int_t triggerno; // Records the trigger number
  Double_t start_time; // Records the start time for the event

  // Any information specific to an individual particle
  TLorentzVector p4; // TLorentzVector for four vector of each particle (Px,Py,Pz,E)
  TLorentzVector vertex;   // Vertex position and time (Vx,Vy,Vz,Vt)
  Double_t beta; // beta from time of flight (TOF)
  Double_t status; // Gives information on which detectors and how many picked it up
  Double_t energy; // Energy measured by detector
  Double_t charge; // Charge measured by detector
  Double_t PID; // Records PID deterimined by TOF
  Double_t chi2PID; // Chi^2 for PID of TOF
  Double_t time; // Time recorded by by FTCAL,FTOF or CTOF
  Double_t path; // Path measured by FTCAL,FTOF or CTOF
  Double_t vertex_time; // Calculated vertex time from TOF information

  // Vectors of particle measurables for when you have more than one of a type of particle (e.g 2 pi^{-})
  vector<TLorentzVector> v_p4;
  vector<TLorentzVector> v_vertex;
  vector<double> v_beta;
  vector<double> v_status;
  vector<double> v_energy;
  vector<double> v_charge;
  vector<double> v_PID;
  vector<double> v_chi2PID;
  vector<double> v_time;
  vector<double> v_path;

  // Record the number of particles in an event
  Int_t elno; // electrons
  Int_t positronno; // positrons
  Int_t protonno; // protons
  Int_t antiprotonno; // antiprotons
  Int_t neutronno; // neutrons
  Int_t photonno; // photons
  Int_t pipno; // pi^{+}
  Int_t pimno; // pi^{-}
  Int_t pi0no; // pi^{0}
  Int_t kaonpno; // K^{+}
  Int_t kaonmno; // K^{-}

  // Setting TLorentzVectors for beam and target, for final analysis beam energy
  // will have to be changed depending on which run/runs you are analysing. Here
  // it is just set to 10.6 GeV
  TLorentzVector beam(0,0,10.6,10.6); // Set 4 vector four the beam, all momentum in z-direction
  TLorentzVector target(0,0,0,0.93827); // Set 4 vector for target, stationary so no momentum

  Double_t c=3*10^8; // Speed of light, used to calculate vertex time

  // Assign a branch to each measurable and name it
  skim4_Tree_280520_01.Branch("eventno",&eventno);
  skim4_Tree_280520_01.Branch("runno",&runno,"runno/I");
  skim4_Tree_280520_01.Branch("triggerno",&triggerno,"triggerno/I");
  skim4_Tree_280520_01.Branch("start_time",&start_time);
  skim4_Tree_280520_01.Branch("p4",&v_p4);
  skim4_Tree_280520_01.Branch("vertex",&v_vertex);
  skim4_Tree_280520_01.Branch("beta",&v_beta);
  skim4_Tree_280520_01.Branch("status",&v_status);
  skim4_Tree_280520_01.Branch("energy",&v_energy);
  skim4_Tree_280520_01.Branch("charge",&v_charge);
  skim4_Tree_280520_01.Branch("PID",&v_PID);
  skim4_Tree_280520_01.Branch("chi2PID",&v_chi2PID);
  skim4_Tree_280520_01.Branch("time",&v_time);
  skim4_Tree_280520_01.Branch("path",&v_path);
  skim4_Tree_280520_01.Branch("beam",&beam);
  skim4_Tree_280520_01.Branch("target",&target);
  skim4_Tree_280520_01.Branch("elno",&elno);
  skim4_Tree_280520_01.Branch("positronno",&positronno);
  skim4_Tree_280520_01.Branch("protonno",&protonno);
  skim4_Tree_280520_01.Branch("antiprotonno",&antiprotonno);
  skim4_Tree_280520_01.Branch("neutronno",&neutronno);
  skim4_Tree_280520_01.Branch("photonno",&photonno);
  skim4_Tree_280520_01.Branch("pipno",&pipno);
  skim4_Tree_280520_01.Branch("pimno",&pimno);
  skim4_Tree_280520_01.Branch("pi0no",&pi0no);
  skim4_Tree_280520_01.Branch("kaonpno",&kaonpno);
  skim4_Tree_280520_01.Branch("kaonmno",&kaonmno);

  // Going over all the input files listed above
  for(Int_t i=0;i<files->GetEntries();i++){

    // Create the CLAS12 event reader
    clas12reader c12(files->At(i)->GetTitle());

    // Prints out the file currently being analysed
    cout<<"file: "<<i<<endl;

    // This loop goes over the events within each file
    while(c12.next()==true){
      // Clear the vectors from the previous event
      v_p4.clear();
      v_vertex.clear();
      v_beta.clear();
      v_status.clear();
      v_energy.clear();
      v_charge.clear();
      v_PID.clear();
      v_chi2PID.clear();
      v_time.clear();
      v_path.clear();

      // Define how to access the information for each event
      eventno = c12.runconfig()->getEvent(); // Getting the event number
      runno = c12.runconfig()->getRun(); // Getting the run number
      triggerno = c12.runconfig()->getTrigger(); // Getting the trigger bit
      start_time = c12.event()->getStartTime(); // Getting start time for each event

      elno = 0;
      positronno = 0;
      protonno = 0;
      antiprotonno = 0;
      neutronno = 0;
      photonno = 0;
      pipno = 0;
      pimno = 0;
      pi0no = 0;
      kaonpno = 0;
      kaonmno = 0;

      // This loop goes over each particle within the current event
      for(auto& p : c12.getDetParticles()){

        // Define how to access the information for each particle
        p4.SetXYZM(p->par()->getPx(), p->par()->getPy(), p->par()->getPz(),0); // Sets the four vector for each particle, mass will be set once PID determined
        time = p->getTime(); // Gets the TOF
        path = p->getPath(); // Gets the path
        beta = p->par()->getBeta(); // Gets the beta from TOF
        vertex_time = time - path / (beta*c); // Calculated vertex time
        vertex.SetXYZT(p->par()->getVx(), p->par()->getVy(), p->par()->getVz(),vertex_time); // Sets the vertex vector for each particle
        status = p->par()->getStatus(); // Information on which detectors are involved for this particle
        energy =  p->getDetEnergy(); // Energy measured by detector
        charge = p->par()->getCharge(); // Charge of the particle measured
        PID = p->par()->getPid();  // PID determined by TOF
        chi2PID = p->par()->getChi2Pid(); // Chi^2 for PID from TOF

        // Save the particle information in the vectors by pushing it back
        v_vertex.push_back(vertex);
        v_beta.push_back(beta);
        v_status.push_back(status);
        v_energy.push_back(energy);
        v_charge.push_back(charge);
        v_PID.push_back(PID);
        v_chi2PID.push_back(chi2PID);
        v_time.push_back(time);
        v_path.push_back(path);


        // To assign particle mass we check the PID and use PDG database
        // Neutral particles first
        if(charge==0){
          // neutrons
          if(p->par()->getPid()==2112){
            p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(2112)->Mass());
            v_p4.push_back(p4); // Recording the 4 vector for each neutron
            neutronno++; // Increasing the count of neutrons
          }
          // pi^{0}
          else if(p->par()->getPid()==111){
            p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(111)->Mass());
            v_p4.push_back(p4); // Recording the 4 vector for each pi^{0}
            pi0no++; // Increasing the count of pi^{0}
          }
          // Photons
          else if(p->par()->getPid()==22){
            p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),0);
            v_p4.push_back(p4); // Recording the 4 vector for each photon
            photonno++; // Increasing the count of photons
          }
        }

        // Only records charged particles within a certain chi^2PID to reduce
        // cuts needed later. Can't use for neutrals as they are all set to 9999
        else if(abs(chi2PID)<5){
          // Electrons
          if(p->par()->getPid()==11){
            p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(11)->Mass());
            v_p4.push_back(p4); // Recording the 4 vector for each electron
            elno++; // Increasing the count of electrons
          }
          // Positrons
          else if(p->par()->getPid()==-11){
            p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(-11)->Mass());
            v_p4.push_back(p4); // Recording the 4 vector for each electron
            positronno++; // Increasing the count of positrons
          }
          // pi^{+}
          else if(p->par()->getPid()==211){
            p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(211)->Mass());
            v_p4.push_back(p4); // Recording the 4 vector for each pi^{+}
            pipno++; // Increasing the count of pi^{+}
          }
          // pi^{-}
          else if(p->par()->getPid()==-211){
            p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(-211)->Mass());
            v_p4.push_back(p4); // Recording the 4 vector for each pi^{-}
            pimno++; // Increasing the count of pi^{-}
          }
          // protons
          else if(p->par()->getPid()==2212){
            p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(2212)->Mass());
            v_p4.push_back(p4); // Recording the 4 vector for each proton
            protonno++; // Increasing the count of protons
          }
          // antiprotons
          else if(p->par()->getPid()==-2212){
            p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(-2212)->Mass());
            v_p4.push_back(p4); // Recording the 4 vector for each antiproton
            antiprotonno++; // Increasing the count of antiprotons
          }
          // K^{+}
          else if(p->par()->getPid()==321){
            p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(321)->Mass());
            v_p4.push_back(p4); // Recording the 4 vector for each K^{+}
            kaonpno++; // Increasing the count of K^{+}
          }
          // K^{-}
          else if(p->par()->getPid()==-321){
            p4.SetXYZM(p4.Px(),p4.Py(),p4.Pz(),db->GetParticle(-321)->Mass());
            v_p4.push_back(p4); // Recording the 4 vector for each K^{-}
            kaonmno++; // Increasing the count of K^{-}
          }
        }
      }

      // Here you can apply a basic skim for events you want to save in your tree
      if(kaonpno>0){
        //Fill the TTree
        skim4_Tree_280520_01.Fill();
      }
    }
  }
  skim4_Tree_280520_01.Write();
  f.Close();
}
