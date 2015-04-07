//STL
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

//Root
#include <TSystem.h>
#include <TString.h>
#include <TChain.h>
#include <TROOT.h> 
#include <TH1F.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1.h> 
#include <TProfile.h>
#include <TH2.h>
#include <TF2.h> 
#include <TMarker.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TLegend.h>
#include <TMath.h>
#include <TVirtualFitter.h>
#include <TBranch.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TNtuple.h>

//HESS
#include "utilities/StringTools.hh"
#include "background/BgMaps.hh"


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Global Variable
Plotters::SkyHist *onmap=0;
std::vector<int> fRunList;
std::vector<unsigned int> fAllowedTelList;
std::map<int, Stash::Coordinate> fObservationList; //!
std::map<int, float> fObservationZenithList; //!
double fEvtOffsetMax=2.7;
double fEvtEnergyMax=100;
double fEvtEnergyMin=0.1;

std::string PmBgRunInfoTreeName = "PmBgRunInfoTree";

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Example : 
//MergeOffRunsNewEnergieMaker("ngc253v14qc70.list","ngc253v14qc70_pmbgE","std_south_1b","ngc253",0.02,120,13)
//MergeOffRunsNewEnergieViewer("RadialLookup_ngc253","std_south_1b",0.02,120,13)
//MergeOffRunsNew("OFFlist_north_wholeTR","OFFlist_pmbgE/","hard_north_1b","RadialLookup_4TeVInf",4.,1000)
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int read_list(std::string listname);
void ReadList(std::vector<std::pair<int, std::string> >& List, 
	      std::string listname,
	      std::vector<double>& zen,
	      std::vector<double>& az);
void Residus(const char* filename,
		  TString path, const char* filename_RL,TString path_RL,TString Config,
		  const char* outname,
		  double EnergieMin=0.,
		  double EnergieMax=1000.);
///////////////////////////////////////////////////////////////////////////////////////////////////////////
int read_list(std::string listname)
{
  
  if (!listname.length()){
    std::cout << "Need list name !!" << std::endl;
    return 0;
  }

  std::ifstream runlist(listname.c_str());
  if (!runlist.is_open()){
    std::cout << "Problem opening list file " << listname << std::endl;
    return 0;
  }
  else std::cout << "Reading File list from: " << listname.c_str() << std::endl;
  
  std::string input;

  while (1){     

    std::getline(runlist, input, '\n');
    if((input.find("*") != std::string::npos) ||
       (input.find("h") != std::string::npos)){
      runlist >> input;
      continue;
    }
    
    if(runlist.eof()) break;
    std::vector<std::string> inlist = 
      split(input.c_str(),' ');
    
    if (inlist.size()==0) continue;
    
    int run = atoi(inlist[0].c_str());
    if(run == 0) continue;
    int allowedtels = 0;
    if(inlist.size()>1) allowedtels = atoi(inlist[1].c_str());
    
    if(allowedtels != 0) {
      fRunList.push_back(run);
      fAllowedTelList.push_back(allowedtels);
      std::cout << "Run " << run << " Tels: " << allowedtels << std::endl;
    } else std::cout << "Rejecting Run " << run 
		     << ", all telescopes failed run selection." << std::endl;
  }
  
  runlist.close();
  std::cout << "Number of files to process in the list: " << fRunList.size() << std::endl;
  return 1;
}

void ReadList(std::vector<std::pair<int, std::string> >& List, std::string listname,
	      std::vector<double>& zen,std::vector<double>& az)
{
  std::ifstream runlist(listname.c_str());
  if (!runlist.is_open()){
    std::cout << "Problem opening list file " << listname << std::endl;
    return;
  }
  else std::cout << "Reading File list from: " 
		 << listname.c_str() << std::endl;
  
  std::string input;

  while (1){     

    std::getline(runlist, input, '\n');
    if(input.find("*") != std::string::npos){
      runlist >> input;
      std::cout << input << std::endl;
      continue;
    }

    if(runlist.eof()) break;
    std::vector<std::string> inlist = 
      split(input.c_str(),' ');
    
    int run = atoi(inlist[0].c_str());
    if(run == 0) continue;
    std::string file = "/data/djannati/root/ph105Sq200/run";
    file += inlist[0];
    file += "_ph105Sq200.root";
    List.push_back(std::make_pair<int, std::string>(run, file));

    zen.push_back(atof(inlist[4].c_str()));
    az.push_back(atof(inlist[5].c_str()));

    std::cout <<file <<" "<< atof(inlist[4].c_str()) << " " << atof(inlist[5].c_str()) << std::endl;

  }

  runlist.close();
   
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////
void Residus(const char* filename,TString path,const char* filename_RL,TString path_RL,TString Config,const char* outname,double EnergieMin,double EnergieMax)
{
  /*
    filename:run list
    path:path where are the root file for the runs
    filename_RL: root file containing the histogram of the loockups
   path_RL: repertoire où se trouve les histos des radialoookups pour les differentes bande en energie, efficacite et zenith.
   Config: config used
   outname:path and filename for the results of the 2D residus maps.
  */

  gROOT->Clear();
  gROOT->Reset();

  // All Vector at zero: 
  fRunList.clear();
  fAllowedTelList.clear();
  fObservationList.clear();
  fObservationZenithList.clear();
  
  int ntotruns = 0;
  
  std::string str("../list/");
  str+= filename;
  read_list(str);
  
  
  //Define Zen Bands and histo names
  //Nband: donne le nombre de band en zenith
  Int_t Nzen=6;
  Int_t Nbands = Nzen-1;
  double zen[Nzen];
  zen[0] = 0.;
  zen[1] = 20.;
  zen[2] = 30.;
  zen[3] = 40.;
  zen[4] = 50.;
  zen[5] = 90.;
  
  Int_t Neff=4;
  Int_t Nbands_eff = Neff-1;
  double eff[Neff];
  eff[0] = 30.;
  eff[1] = 50.;
  eff[2] = 70.;
  eff[3] = 100.;
  
  Int_t Nbands_E=1;
  double E[Nbands_E+1];
  //Int_t Emax_bin=10;
  E[0]=fEvtEnergyMin;
  E[1]=fEvtEnergyMax;

  /*Int_t Nbands_E=50;
  double E[Nbands_E+1];
  Int_t Emax_bin=10;
  double bin_E=  (log10(Emax_bin)-log10(fEvtEnergyMin))/(Nbands_E-1);
  std::cout << bin_E<< endl;
  for(int ie=0; ie<Nbands_E;ie++){
      E[ie]=log10(fEvtEnergyMin)+ie*bin_E;
      std::cout << ie << " " << E[ie]<< endl;
  }
  E[Nbands_E]=log10(fEvtEnergyMax);
  std::cout<<Nbands_E << " " << E[Nbands_E] << endl;*/
  
  Int_t Nband_Eresidus=1;
  double E_residus[Nband_Eresidus+1];
  E_residus[0]=fEvtEnergyMin;
  E_residus[1]=fEvtEnergyMax;
  
  /*Int_t Nband_Eresidus=3;
  double E_residus[Nband_Eresidus+1];
  E_residus[0]=0.1;
  E_residus[1]=5;
  E_residus[2]=10;
  E_residus[3]=100;*/


  //Define the 2D histogram containing 
  std::vector< std::vector<std::vector<TH2F*> > > hist_nomevents;
  std::vector< std::vector<std::vector<TH2F*> > > hist_nomradiallookup;
  std::vector< std::vector<std::vector<TH2F*> > > hist_residus;
  hist_nomevents.resize(Nbands);
  hist_nomradiallookup.resize(Nbands);
  hist_residus.resize(Nbands);
  for (int i_zen=0;i_zen<Nbands;i_zen++){
    hist_nomevents[i_zen].resize(Nbands_eff);
    hist_nomradiallookup[i_zen].resize(Nbands_eff);
    hist_residus[i_zen].resize(Nbands_eff);
    for(int i_eff=0;i_eff<Nbands_eff;i_eff++) {
      hist_nomevents[i_zen][i_eff].resize(Nband_Eresidus);
      hist_nomradiallookup[i_zen][i_eff].resize(Nband_Eresidus);
      hist_residus[i_zen][i_eff].resize(Nband_Eresidus);
    }
  }
  std::vector<std::vector<std::string> > histnames;
  histnames.clear();
  
  
  std::string tabname_events;
  std::string tabname_RL;
  std::string tabname_plot;
  std::string zenmin, zenmax, effmin, effmax;
  std::string Emin, Emax;
  //besoin de ce string quand on ira chercher les histos
  std::string E_number;
  double nom_max=2.7;
  int bin_nom=700;
  for (int izen=0;izen<Nbands;izen++) 
    {
      for (int ieff=0;ieff<Nbands_eff;ieff++) 
	{
	  ostringstream str_convert1, str_convert2 ,str_convert3, str_convert4;
	  str_convert1 << zen[izen];
	  zenmin=str_convert1.str();
	  str_convert2 << zen[izen+1];
	  zenmax=str_convert2.str();
	  str_convert3 << eff[ieff];
	  effmin=str_convert3.str();
	  str_convert4 << eff[ieff+1];
	  effmax=str_convert4.str();
	  for (int i_E=0;i_E<Nband_Eresidus;i_E++) 
	    {
	      //ostringstream str_convert5, str_convert6;
	      //str_convert5 << E_residus[i_E];
	      //str_convert6 << E_residus[i_E+1];
	      ostringstream str_convert5;
	      str_convert5 << i_E;
	      Emin=str_convert5.str();
	      //Emax=str_convert6.str();
	      
	      tabname_events="Nominalposition_events_zen"+zenmin+"_"+zenmax+"deg_eff"+effmin+"_"+effmax+"_bin"+Emin;
	      tabname_RL="Nominalposition_RL_zen"+zenmin+"_"+zenmax+"deg_eff"+effmin+"_"+effmax+"_bin"+Emin;
	       tabname_plot="zen"+zenmin+"_"+zenmax+"deg_eff"+effmin+"_"+effmax+"_Emin_"+Emin+"Emax"+Emax+"TeV";
	     
	       hist_residus[izen][ieff][i_E] = new TH2F(tabname_events.c_str(),tabname_plot.c_str(),bin_nom, -nom_max, nom_max, bin_nom, -nom_max, nom_max);
	       hist_nomradiallookup[izen][ieff][i_E] = new TH2F(tabname_RL.c_str(),tabname_plot.c_str(),bin_nom, -nom_max, nom_max, bin_nom, -nom_max, nom_max);
	       
	      
	    }
	}
    }
  //Je pense que ca veut dire qu'il selectionne la methode du packman
  TString bgmakername("PMBgMaker");
  TString rootfilename_RL;
  rootfilename_RL=path_RL;
  rootfilename_RL+=filename_RL;  
  TFile *file_RL=new TFile(rootfilename_RL);
  TH1F* hist = NULL;
  std::string name_hist;
  std::string radname = "RadialLookup";
  double theta2;
  //Loop on runlist
  //fRunList: defini au debut du programme:vecteur d'entier
  for(unsigned int irun=0; irun< fRunList.size() ; irun++) 
    {
      TString rootfilename;
      rootfilename = "../root/";
      rootfilename += path;
      rootfilename+="/run";
      rootfilename+= fRunList[irun];
      std::cout << "\033[1;34;40m" << "---------------------------------Starting run :" << fRunList[irun] << "\033[0m" << std::endl;
      rootfilename+= "_"; 
      rootfilename+=Config; 
      if (fAllowedTelList.at(irun)!=30) 
	rootfilename+="_3Tel";
      rootfilename+=".root";
     
      //string.data: fname et rootfilename ont le meme content
      std::string fname=rootfilename.Data();
      std::cout << "Reading from file : " << rootfilename << std::endl;
      
      //En fait dans les rootfilename il y a les TTre que anne a cree
      TFile *filehandler=new TFile(rootfilename);
      //In case the file does not exist or is not a valid ROOT file, it is made a Zombie
      //donc check si il a pas de problem avec le file
      if(filehandler->IsZombie()) 
	{
	  std::cout << "\033[1;31;40m" << "Problem opening root file " << rootfilename << " Skipping it! " << "\033[0m" << std::endl;
	  continue;
	}
      
      gROOT->cd();	      
    
      
      
      //take each events, read leaves the and project it in a file
      // la lis les deux TTree present dans le TFile de nom RunInfoTree et EventsTree_BgMakerOff et les recupere avec get
      TTree *treeinfo = NULL;
      treeinfo = (TTree*)filehandler->Get("RunInfoTree");
      TTree *treeoff = NULL;
      treeoff = (TTree*)filehandler->Get("EventsTree_BgMakerOff");
      //treeinfo->Print();
      //treeoff->Print();
      if (!treeoff || !treeinfo) 
	{
	  std::cout << "Pb with the file !  " << std::endl;
	  continue;
	}
      
      TBranch *treeinfo_muoneff  = NULL;
      TBranch *treeinfo_alt  = NULL;
      TBranch *treeoff_evtenergy = NULL;
      TBranch *treeoff_evtoffset = NULL;
      TBranch *treeoff_evtXnominal = NULL;
      TBranch *treeoff_evtYnominal = NULL;
      //recupere les branch des deux TTree precedent:1) recupere l efficacite du run, 2)recupere l energie et l offset de chaque evenement
      treeinfo_muoneff  = treeinfo->GetBranch("fMuonEff");
      treeinfo_alt  = treeinfo->GetBranch("fPointALT");
      treeoff_evtenergy = treeoff->GetBranch("fEvtEnergy");
      treeoff_evtoffset = treeoff->GetBranch("fEvtOffset");
      treeoff_evtXnominal = treeoff->GetBranch("fEvtNomX");;
      treeoff_evtYnominal = treeoff->GetBranch("fEvtNomY");;
      
      Double_t fAlt;
      Double_t zenit;
      Double_t fMuonEff;
      Double_t fOffEvtEnergy;
      Double_t fOffEvtOffset;
      Double_t Xnom;
      Double_t Ynom;
      //voir si ces lignes permettent de lier les valeurs dans lesTTRee et les doubles crees, si on modifie l'un on omodifie l'autre, je comprend pas a quoi ca sert.
      // ca va lui servir juste en dessous dans la boucle sur les evenements, des qu'il va lire le TTree evenement par evenement avec GetEntry(iev), ca donnera la valeur lu de l'offset a fOffEvtOffset et la valeur lu de l'energy a fOffEvtEnergy
      treeinfo_muoneff->SetAddress(&fMuonEff);
      treeinfo_alt->SetAddress(&fAlt);
      treeoff_evtenergy->SetAddress(&fOffEvtEnergy);    
      treeoff_evtoffset->SetAddress(&fOffEvtOffset);
      treeoff_evtXnominal->SetAddress(&Xnom);
      treeoff_evtYnominal->SetAddress(&Ynom);
      // Je pense que premiere permet de recuperer efficacite du run et comme avec les lignes precedentes on l a relier a fMuonEff, fMuonEff=efficacité du run
      treeinfo_alt->GetEntry(0);
      zenit=90-fAlt;
      treeinfo_muoneff->GetEntry(0);
      fMuonEff=fMuonEff*100.;
      //SANITY CHECK
      if (zenit>=65.) {
	//std::cout << "\033[1;31;40m" << "The zenith angle is " << zenit << " deg, we skip it ! " << "\033[0m" << std::endl; 
	continue;
      }           
      //std::cout << "Zenith Angle : " << zenit  << " deg" << std::endl;
      int index = 0;
      for(int izen=0; izen<Nbands;izen++){
	if(zenit>zen[izen] && zenit<zen[izen+1]){
	    index=izen;
	    std::cout << izen << " " << zenit << endl;
	    break;
	  }
	  else if(zenit>zen[Nbands]){
	    //std::cout << " zenith superieur a zenith max" << std::endl;
	    break;
	  }
      }

      int index_eff = 0;
      for(int ieff=0; ieff<Nbands_eff;ieff++){
	if(fMuonEff>eff[ieff] && fMuonEff<eff[ieff+1]){
	    index_eff=ieff;
	    //std::cout << ieff << " " << fMuonEff << endl;
	    break;
	  }
	  else if(fMuonEff>eff[Nbands_eff]){
	    //std::cout << " efficacite superieur a efficacite max" << std::endl;
	    break;
	  }
      }
      
      // Fill histo.
      TH1F *fLookupEthresh = NULL;
      double offsetmax = 2.5;
      int nbinoff = 250;
      fLookupEthresh = new TH1F("lookupethresh","lookupethresh",nbinoff,0.,offsetmax);
    
      //J'ai l'impression que la ca lui permet de remplir un histo de 250 bins de Ethreshold a zero partout
       for (int il=1;il<=fLookupEthresh->GetNbinsX();++il) 
	{	
	  double myoff = fLookupEthresh->GetBinCenter(il);
	  double ethresh = 0.;
	  fLookupEthresh->SetBinContent(il,ethresh);
	  //std::cout << "il = " << il << " offset = " << myoff << " meanzenith = " << zenit << " fMuonEff = " << fMuonEff << " ethresh = " << ethresh << std::endl;
	}
      
      //Loop over the events 
      for (Int_t ievt=0;ievt<treeoff->GetEntries();++ievt) 
      {
	
	//Recupere pour chaque evenement energie et offset
	treeoff_evtenergy->GetEntry(ievt);    
	treeoff_evtoffset->GetEntry(ievt);
	treeoff_evtXnominal->GetEntry(ievt);
	treeoff_evtYnominal->GetEntry(ievt);
	// Ici fOffEvtOffset=a la valeur de l'offset qui vient d etre lu dans le TTree pour l evenement ievt.
	//Given a point fOffEvtOffset, approximates the value via linear interpolation based on the two nearest bin centers. Du coup vu que precedement pour chaque offset on a cree une valeur du ethreshol, pour chaque offset des evenements il peut remonter en interpolant entre les deux bins les plus proches, la valeur du ethreshold correspondant
	//seul truc je vois pas a quoi ca sert sachant qu'il a tout mis a zero avant dans l histo fLookupEthresh donc la valeur interpollé est aussi zero
	double ethresh = fLookupEthresh->Interpolate(fOffEvtOffset);
	//fEvtOffsetMax mis a 2.7 et defini au debut du programme
	Int_t index_E=0;
	Int_t index_E_residus=0;
	
	for(int ie=0; ie<Nbands_E;ie++){
	  if(fOffEvtEnergy>E[ie] && fOffEvtEnergy<E[ie+1]){
	    index_E=ie;
	    //std::cout << ie << " " << fOffEvtEnergy << endl;
	    break;
	       }
	}
	
	for(int ie_res=0; ie_res<Nband_Eresidus;ie_res++){
	  if(fOffEvtEnergy>E_residus[ie_res] && fOffEvtEnergy<E_residus[ie_res+1]){
	    index_E_residus=ie_res;	    
	    break;
	  }
	}
	if ( (fOffEvtEnergy >= EnergieMin) && (fOffEvtEnergy < EnergieMax) && (fOffEvtOffset < fEvtOffsetMax ) && (fOffEvtEnergy >= ethresh) ) 
	  {
	    //la on rempli les deux vecteur d'histo quon avait cree precedment en offset**2 et offset de 700 et 250 bins et de valeur max fEvtOffsetMax**2 et fEvtOffsetMax* respectivement. index donne la bande en zenith dans laquelle on est, pour chaque histo de bande en zenoth, on rempli l'histogramme avec l'offset de chaque evenement. Donc pour chaque bande en zenith, on a pour les 750 ou 250 bandes offset defini pour les histo hists et histsR le nombre d evenement correspondant.On peut donc avoir les radial lookup
	    hist_residus[index][index_eff][index_E_residus]->Fill(Xnom*180/TMath::Pi(),Ynom*180/TMath::Pi());
	    //hist_nomradiallookup[index][index_eff][index_E]->Fill(Xnom,Ynom);
	   
	    //std::cout << "Added the " << ievt << " Energy = " << fOffEvtEnergy << " Offset = " << fOffEvtOffset << " ethresh = " << ethresh << std::endl;
	  }
      }
      
      std::cout << "Run " << fRunList[irun] << " index : " << index <<  " Total Added Run = " << ++ntotruns << std::endl;
      
      delete fLookupEthresh;
      if (filehandler) delete filehandler;
      
      std::cout << "\033[1;32;40m" << "----------------> Run " << fRunList.at(irun) << " Sucessfully Added !  Total : " << ++ntotruns   << "\033[0m" << std::endl; 
      }
  double events_RL;
  int bin_number_RL, bin_number_2D;
  TCanvas *c1; 
  TCanvas *c2;
  string name_canvas;
  string name_hist_residus;
  TString foutname(outname);
  foutname+="_";
  foutname+=Config;
  foutname+=".root";
  //TCanvas *c1;
  //TCanvas *c2;
  TFile* out = new TFile(foutname,"recreate");
  out->cd();
  //for (int izen=0;izen<Nbands;izen++) 
  for (int izen=0;izen<1;izen++) 
    {
      //for (int ieff=1;ieff<Nbands_eff;ieff++) 
       for (int ieff=1;ieff<2;ieff++) 
	{
	  ostringstream str_convert1, str_convert2 ,str_convert3, str_convert4;
	  str_convert1 << zen[izen];
	  zenmin=str_convert1.str();
	  str_convert2 << zen[izen+1];
	  zenmax=str_convert2.str();
	  str_convert3 << eff[ieff];
	  effmin=str_convert3.str();
	  str_convert4 << eff[ieff+1];
	  effmax=str_convert4.str();
	  for (int i_Eresidus=0;i_Eresidus<Nband_Eresidus;i_Eresidus++) 
	    {
	      //std::cout << "I_Eresidus" << i_Eresidus << std::endl;
	      //std::cout << std::endl;
	      string binE;
	      ostringstream str_convert10;
	      str_convert10 << i_Eresidus;
	      binE=str_convert10.str();
	      ostringstream str_convert5,str_convert6;
	      str_convert5 << E_residus[i_Eresidus];
	      str_convert6 << E_residus[i_Eresidus+1];
	      Emin=str_convert5.str();
	      Emax=str_convert6.str();
	      name_canvas="zen"+zenmin+"_"+zenmax+"deg_eff"+effmin+"_"+effmax+"_Emin_"+Emin+"Emax"+Emax+"TeV.jpg";
	      name_hist_residus="Residus_zen"+zenmin+"_"+zenmax+"deg_eff"+effmin+"_"+effmax+"_bin"+binE;
	      for (int i_E=0;i_E<Nbands_E;i_E++) 
		{ 
		  if((E[i_E] >= E_residus[i_Eresidus]) & (E[i_E] <= E_residus[i_Eresidus+1])){
		    //juste pour verifier que pour la grosse bande les i_E s'enchaine bien de maniere continue
		    std::cout << "iE " << i_E << " min "<< E[i_E]<< " max "<< E[i_E+1]<< std::endl;
		    std::cout << "i_Eresidus " << i_Eresidus << " min "<< E_residus[i_Eresidus]<< " max "<< E_residus[i_Eresidus+1]<< std::endl;
		    //std ::cout << i_E << endl;
		    ostringstream str_convert7;
		    str_convert7 << i_E;
		    E_number=str_convert7.str();
		    name_hist=radname+"_zen"+zenmin+"_"+zenmax+"deg_eff"+effmin+"_"+effmax+"_BinE_"+E_number;
		    hist = (TH1F*)file_RL->Get(name_hist.c_str());
		    //c2=new TCanvas("pp");
		    //hist->Draw();
		    //c2->SaveAs("hist");
		    double Xaxis, Yaxis;
		    for (int iX=1; iX <= bin_nom; iX++){
		      for (int iY=1; iY <= bin_nom; iY++){
			//std::cout << "IXIY" << iX<<" " << iY<< std::endl;
			//std::cout << std::endl;
			Xaxis=hist_nomradiallookup[izen][ieff][i_Eresidus]->GetXaxis()->GetBinCenter(iX);
			Yaxis=hist_nomradiallookup[izen][ieff][i_Eresidus]->GetYaxis()->GetBinCenter(iY);
			theta2=(Xaxis*Xaxis+Yaxis*Yaxis);
			//std::cout << "theta2 " << theta2 << std::endl;
			bin_number_RL=hist->FindBin(theta2);
			//std::cout << "bin_number_RL " <<bin_number_RL <<std::endl;
			events_RL=hist->GetBinContent(bin_number_RL);
			//std::cout << "event_RL " << events_RL << std::endl;
			bin_number_2D=hist_nomradiallookup[izen][ieff][i_Eresidus]->FindBin(Xaxis,Yaxis);
			//std::cout << "bin_number_2D " << bin_number_2D << std::endl;
			//hist_nomradiallookup[izen][ieff][i_Eresidus]->AddBinContent(bin_number_2D,events_RL);
			std::cout << hist_nomradiallookup[izen][ieff][i_Eresidus]->GetBinContent(bin_number_2D) << std::endl;
			std::cout << "event_RL " << events_RL << std::endl;
			hist_nomradiallookup[izen][ieff][i_Eresidus]->SetBinContent(bin_number_2D,hist_nomradiallookup[izen][ieff][i_Eresidus]->GetBinContent(bin_number_2D)+events_RL);
			std::cout << hist_nomradiallookup[izen][ieff][i_Eresidus]->GetBinContent(bin_number_2D) << std::endl;
		      }
		    }
		    
		  }
		  
		}
	      //trouver koverwrite et comment changer nom suaver de histo pour pas remplacer
	      
	      hist_residus[izen][ieff][i_Eresidus]->Scale(1/ hist_residus[izen][ieff][i_Eresidus]->GetMaximum());
	      hist_nomevents[izen][ieff][i_Eresidus]=(TH2F*)hist_residus[izen][ieff][i_Eresidus]->Clone();
	      
	     
	      hist_nomradiallookup[izen][ieff][i_Eresidus]->Scale(1/hist_nomradiallookup[izen][ieff][i_Eresidus]->GetMaximum());
	      //std::cout << hist_residus[izen][ieff][i_Eresidus]->Integral() << std::endl;
	      //std::cout << hist_nomradiallookup[izen][ieff][i_Eresidus] << std::endl;
	      std::cout << std::endl;
	      hist_nomevents[izen][ieff][i_Eresidus]->Write();
	      hist_nomradiallookup[izen][ieff][i_Eresidus]->Write();
	      c1=new TCanvas(name_canvas.c_str());
	      c1->Divide(2,2);
	      c1->cd(1);
	      hist_nomevents[izen][ieff][i_Eresidus]->Draw("colz");
	      c1->cd(2);
	      hist_nomradiallookup[izen][ieff][i_Eresidus]->Draw("colz");
	      hist_residus[izen][ieff][i_Eresidus]->Add(hist_nomradiallookup[izen][ieff][i_Eresidus],-1);
	      
	      hist_residus[izen][ieff][i_Eresidus]->SetName(name_hist_residus.c_str());
	      hist_residus[izen][ieff][i_Eresidus]->Write();
	      //surement fire un TPAD
	      
	      c1->cd(3);
	      //hist_nomevents[izen][ieff][i_Eresidus]->Draw("colz");
	      hist_residus[izen][ieff][i_Eresidus]->Draw("colz");
	      c1->SaveAs(name_canvas.c_str());
	      //c1->SaveAs("test.jpg");
	    }
	}
    }
  out->Close();
  std::cout<< "c1"<< std::endl;
  delete c1;
  delete c2;
  std::cout<< "c1"<< std::endl;
  delete file_RL;
  std::cout<< "file"<< std::endl;
  //il ne faut pas delete le hist car comme formé a partir de file_RL, quand je delete file_RL ca delete hist
  //delete hist;
  std::cout<< "hist"<< std::endl;
  delete out;
  std::cout<< "out"<< std::endl;

  /*double events_RL;
  int bin_number_RL, bin_number_2D;
  TCanvas *c1; 
  string name_canvas;
  string name_hist_residus;
  TString foutname(outname);
  foutname+="_";
  foutname+=Config;
  foutname+=".root";
  //TCanvas *c1;
  //TCanvas *c2;
  TFile* out = new TFile(foutname,"recreate");
  out->cd();
  for (int izen=0;izen<Nbands;izen++) 
    {
      for (int ieff=0;ieff<Nbands_eff;ieff++) 
	{
	  ostringstream str_convert1, str_convert2 ,str_convert3, str_convert4;
	  str_convert1 << zen[izen];
	  zenmin=str_convert1.str();
	  str_convert2 << zen[izen+1];
	  zenmax=str_convert2.str();
	  str_convert3 << eff[ieff];
	  effmin=str_convert3.str();
	  str_convert4 << eff[ieff+1];
	  effmax=str_convert4.str();
	  for (int i_Eresidus=0;i_Eresidus<Nband_Eresidus;i_Eresidus++) 
	    {
	      //std::cout << "I_Eresidus" << i_Eresidus << std::endl;
	      //std::cout << std::endl;
	      ostringstream str_convert5,str_convert6;
	      str_convert5 << E_residus[i_Eresidus];
	      str_convert6 << E_residus[i_Eresidus+1];
	      Emin=str_convert5.str();
	      Emax=str_convert6.str();
	      name_canvas="zen"+zenmin+"_"+zenmax+"deg_eff"+effmin+"_"+effmax+"_Emin_"+Emin+"Emax"+Emax+"TeV.jpg";
	      name_hist_residus="Residus_zen"+zenmin+"_"+zenmax+"deg_eff"+effmin+"_"+effmax+"_Emin_"+Emin+"Emax"+Emax+"TeV.jpg";
	      for (int i_E=0;i_E<Nbands_E;i_E++) 
		{ 
		  if((E[i_E] > E_residus[i_Eresidus]) & (E[i_E] < E_residus[i_Eresidus+1])){
		    //juste pour verifier que pour la grosse bande les i_E s'enchaine bien de maniere continue
		     
		    //std ::cout << i_E << endl;
		    ostringstream str_convert7;
		    str_convert7 << i_E;
		    E_number=str_convert7.str();
		    name_hist=radname+"_zen"+zenmin+"_"+zenmax+"deg_eff"+effmin+"_"+effmax+"_BinE_"+E_number;
		    hist = (TH1F*)file_RL->Get(name_hist.c_str());
		    double Xaxis, Yaxis;
		    for (int iX=1; iX <= bin_nom; iX++){
		      for (int iY=1; iY <= bin_nom; iY++){
			//std::cout << "IXIY" << iX<<" " << iY<< std::endl;
			//std::cout << std::endl;
			Xaxis=hist_nomradiallookup[izen][ieff][i_Eresidus]->GetXaxis()->GetBinCenter(iX);
			Yaxis=hist_nomradiallookup[izen][ieff][i_Eresidus]->GetYaxis()->GetBinCenter(iY);
			theta2=(Xaxis*Xaxis+Yaxis*Yaxis);
			bin_number_RL=hist->FindBin(theta2);
			events_RL=hist->GetBinContent(bin_number_RL);
			bin_number_2D=hist_nomradiallookup[izen][ieff][i_Eresidus]->FindBin(Xaxis,Yaxis);
			hist_nomradiallookup[izen][ieff][i_Eresidus]->AddBinContent(bin_number_2D,events_RL);
		      }
		    }

		  }

		}
	      //trouver koverwrite et comment changer nom suaver de histo pour pas remplacer
	      hist_nomevents[izen][ieff][i_Eresidus]->Scale(1/hist_nomevents[izen][ieff][i_Eresidus]->Integral());
	      hist_nomradiallookup[izen][ieff][i_Eresidus]->Scale(1/hist_nomradiallookup[izen][ieff][i_Eresidus]->Integral());
	      hist_nomevents[izen][ieff][i_Eresidus]->Write("",TObject::kOverwrite);
	      hist_nomradiallookup[izen][ieff][i_Eresidus]->Write("",TObject::kOverwrite);
	      c1=new TCanvas(name_canvas.c_str());
	      c1->Divide(2,2);
	      c1->cd(1);
	      hist_nomevents[izen][ieff][i_Eresidus]->Draw("colz");
	      c1->cd(2);
	      hist_nomradiallookup[izen][ieff][i_Eresidus]->Draw("colz");
	      
	      hist_nomevents[izen][ieff][i_Eresidus]->Add(hist_nomradiallookup[izen][ieff][i_Eresidus],-1);
	      hist_residus[izen][ieff][i_Eresidus]=(TH2F*)hist_nomevents[izen][ieff][i_Eresidus]->Clone();
	      hist_residus[izen][ieff][i_Eresidus]->SetName(name_hist_residus.c_str());
	      hist_residus[izen][ieff][i_Eresidus]->Write("",TObject::kOverwrite);
	      //surement fire un TPAD
	      
	      c1->cd(3);
	      hist_residus[izen][ieff][i_Eresidus]->Draw("colz");
	      c1->SaveAs(name_canvas.c_str());
	     
	    }
	}
    }
  out->Close();
  std::cout<< "c1"<< std::endl;
  delete c1;
  std::cout<< "c1"<< std::endl;
  delete file_RL;
  std::cout<< "file"<< std::endl;
  //il ne faut pas delete le hist car comme formé a partir de file_RL, quand je delete file_RL ca delete hist
  //delete hist;
  std::cout<< "hist"<< std::endl;
  delete out;
  std::cout<< "out"<< std::endl;*/
}

