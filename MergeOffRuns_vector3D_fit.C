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
void MergeOffRuns(const char* filename,
		  TString path,
		  TString Config,
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
void MergeOffRuns(const char* filename,TString path,TString Config,const char* outname,double EnergieMin,double EnergieMax)
{
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
  /*Int_t Neff=2;
  Int_t Nbands_eff = Neff-1;
  double eff[Neff];
  eff[0] = 0.;
  eff[1] = 100.;*/

  Int_t Nbands_E=50;
  double E[Nbands_E+1];
  Int_t Emax_bin=10;
  double bin_E=  (log10(Emax_bin)-log10(fEvtEnergyMin))/(Nbands_E-1);
  std::cout << bin_E<< endl;
  for(int ie=0; ie<Nbands_E;ie++){
      E[ie]=log10(fEvtEnergyMin)+ie*bin_E;
      std::cout << ie << " " << E[ie]<< endl;
  }
  E[Nbands_E]=log10(fEvtEnergyMax);
  std::cout<<Nbands_E << " " << E[Nbands_E] << endl;
  TNtuple distrib_zenith("distrib_zenith", "distribution des zenith des runs", "zenith");
  TNtuple distrib_eff("distrib_efficacite", "distribution des efficacite des runs", "efficacite");
  TNtuple distrib_energy("distrib_energy", "distribution des energy des evenements", "energy");
  


  std::vector<std::vector<std::string> > histnames;
  histnames.clear();
  //histnames.push_back(std::string("Zen0_20deg"));
  //histnames.push_back(std::string("Zen20_30deg"));
  //histnames.push_back(std::string("Zen30_40deg"));
  //histnames.push_back(std::string("Zen40_45deg"));
  //histnames.push_back(std::string("Zen45_55deg"));
  //histnames.push_back(std::string("Zen55_90deg"));
  
  
  std::string radname = "RadialLookup";
  std::string radnameR = "RadialLookup_R"; 

  std::vector< std::vector<std::vector<TH1F*> > > hists;
  std::vector < std::vector<std::vector<TH1F*> > > histsR;
  std::vector< std::vector<std::vector<TH1F*> > > histssmooth;
  std::vector < std::vector<std::vector<TH1F*> > > histsfit;
  std::vector<std::vector<TH1F*> > Nevents_E_band_zen_eff;
  
  hists.resize(Nbands);
  histsR.resize(Nbands);
  histssmooth.resize(Nbands);
  histsfit.resize(Nbands);
  Nevents_E_band_zen_eff.resize(Nbands);
  //histssmooth.resize(Nbands);
  //histsfit.resize(Nbands);
  
  for (int i_zen=0;i_zen<Nbands;i_zen++){
    Nevents_E_band_zen_eff[i_zen].resize(Nbands_eff);
    hists[i_zen].resize(Nbands_eff);
    histsR[i_zen].resize(Nbands_eff);
    histssmooth[i_zen].resize(Nbands_eff);
    histsfit[i_zen].resize(Nbands_eff);
    for(int i_eff=0;i_eff<Nbands_eff;i_eff++) {
      hists[i_zen][i_eff].resize(Nbands_E);
      histsR[i_zen][i_eff].resize(Nbands_E);
      histssmooth[i_zen][i_eff].resize(Nbands_E);
      histsfit[i_zen][i_eff].resize(Nbands_E);
    }
  }
  
  
  TH1F *Nrun_zenith;
  TH1F *Nevent_zenith;
  Nrun_zenith=new TH1F("Nrun_bin_zenith","Nrun_bin_zenith", Nbands,0, Nbands); 
  Nevent_zenith=new TH1F("Nevent_bin_zenith","Nrun_event_zenith", Nbands,0, Nbands);
  TH1F *Nrun_eff;
  TH1F *Nevent_eff;
  Nrun_eff=new TH1F("Nrun_bin_efficacite","Nrun_bin_efficacite", Nbands_eff,0,Nbands_eff); 
  Nevent_eff=new TH1F("Nevent_bin_efficacite","Nrun_event_efficacite",Nbands_eff,0,Nbands_eff);
  TH1F *Nevent_E;
  Nevent_E=new TH1F("Nevent_bin_energy","Nrun_event_energy",Nbands_E,0,Nbands_E);
  // Initializing Tables : 
  //loop sur les N bands en zenith et remplie les histogrammes qui ont le nom de radia_lookup+_nom de chaque band en zenith.
  for (int izen=0;izen<Nbands;izen++) 
    {
      for (int ieff=0;ieff<Nbands_eff;ieff++) 
	{
	  std::string tabname;
	  std::string tabnameR;
	  std::string zenmin, zenmax, effmin, effmax;
	  ostringstream str_convert1, str_convert2 ,str_convert3, str_convert4;
	  str_convert1 << zen[izen];
	  zenmin=str_convert1.str();
	  str_convert2 << zen[izen+1];
	  zenmax=str_convert2.str();
	  str_convert3 << eff[ieff];
	  effmin=str_convert3.str();
	  str_convert4 << eff[ieff+1];
	  effmax=str_convert4.str();
	  std::string histoNeventE_2D_name="Nevent_bin_energy_zen"+zenmin+"_"+zenmax+"deg_eff"+effmin+"_"+effmax;
	  Nevents_E_band_zen_eff[izen][ieff]=new TH1F(histoNeventE_2D_name.c_str(),histoNeventE_2D_name.c_str(),Nbands_E,0,Nbands_E);
	  for (int i_E=0;i_E<Nbands_E;i_E++) 
	    {
	      std::string E_number;
	      ostringstream str_convert5;
	      str_convert5 << i_E;
	      E_number=str_convert5.str();
	      
	      tabname=radname+"_zen"+zenmin+"_"+zenmax+"deg_eff"+effmin+"_"+effmax+"_BinE_"+E_number;
	      tabnameR=radnameR+"_zen"+zenmin+"_"+zenmax+"deg_eff"+effmin+"_"+effmax+"_BinE_"+E_number;
	    //hists et histsR c'est des vecteurs d'histogram de root. 
	    //fEvtOffsetMax: defini au debut du program, donne l'offset max ici a 2.7
	    //Pour chaque band vecteur donc chaque band de zenith cree un histos de 700 bins de O a fEvtOffsetMax**2:histogramme des theta**2
	    //et un histogram de 250 bins de 0 a fEvtOffsetMax: histogram des theta
	      hists[izen][ieff][i_E] = new TH1F(tabname.c_str(),tabname.c_str(),700,0,TMath::Power(fEvtOffsetMax,2.));
	      histsR[izen][ieff][i_E] = new TH1F(tabnameR.c_str(),tabnameR.c_str(),250,0,fEvtOffsetMax);
	    }
	}
    }
  //Je pense que ca veut dire qu'il selectionne la methode du packman
  TString bgmakername("PMBgMaker");
  
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
      std::string zenithhistname(bgmakername);
      zenithhistname+="_Current_ZenithDistOn";
      TH1F* hz = NULL;
      // cree un histogramme hz avec un des histogrmme du filehander censé contenir les donnees d'un TTree. Il recupere l'histo du Bgmaker
      //voir avec vincent ce qu'il contient ce truc zenithhistname.c_str()?
      hz = (TH1F*)filehandler->Get(zenithhistname.c_str());
      //SANITY CHECK
      if (hz==NULL) {
	std::cout  << "\033[1;31;40m" << "Can't find " << zenithhistname << " histogram in the root file, skipping it ! " << "\033[0m" << std::endl;
	//il trouve pas l'histogramme du bgmaker dansc ce TTre donc la passe au run d'apres sur la boucle for avec le continue
	continue;
      }
      float zenit=0;
      zenit=hz->GetMean();
      delete hz;
      //SANITY CHECK
      if (zenit>=65.) {
	std::cout << "\033[1;31;40m" << "The zenith angle is " << zenit << " deg, we skip it ! " << "\033[0m" << std::endl; 
	continue;
      }           
      std::cout << "Zenith Angle : " << zenit  << " deg" << std::endl;
      int index = 0;
      for(int izen=0; izen<Nbands;izen++){
	if(zenit>zen[izen] && zenit<zen[izen+1]){
	    index=izen;
	    std::cout << izen << " " << zenit << endl;
	    break;
	  }
	  else if(zenit>zen[Nbands]){
	    std::cout << " zenith superieur a zenith max" << std::endl;
	    break;
	  }
      }
      
      
      //take each events, read leaves the and project it in a file
      // la lis les deux TTree present dans le TFile de nom RunInfoTree et EventsTree_BgMakerOff et les recupere avec get
      TTree *treeinfo = NULL;
      treeinfo = (TTree*)filehandler->Get("RunInfoTree");
      TTree *treeoff = NULL;
      treeoff = (TTree*)filehandler->Get("EventsTree_BgMakerOff");
      treeinfo->Print();
      treeoff->Print();
      if (!treeoff || !treeinfo) 
	{
	  std::cout << "Pb with the file !  " << std::endl;
	  continue;
	}
      
      TBranch *treeinfo_muoneff  = NULL;
      TBranch *treeoff_evtenergy = NULL;
      TBranch *treeoff_evtoffset = NULL;
      //recupere les branch des deux TTree precedent:1) recupere l efficacite du run, 2)recupere l energie et l offset de chaque evenement
      treeinfo_muoneff  = treeinfo->GetBranch("fMuonEff");
      treeoff_evtenergy = treeoff->GetBranch("fEvtEnergy");
      treeoff_evtoffset = treeoff->GetBranch("fEvtOffset");
      
      Double_t fMuonEff;
      Double_t fOffEvtEnergy;
      Double_t fOffEvtOffset;
      
      //voir si ces lignes permettent de lier les valeurs dans lesTTRee et les doubles crees, si on modifie l'un on omodifie l'autre, je comprend pas a quoi ca sert.
      // ca va lui servir juste en dessous dans la boucle sur les evenements, des qu'il va lire le TTree evenement par evenement avec GetEntry(iev), ca donnera la valeur lu de l'offset a fOffEvtOffset et la valeur lu de l'energy a fOffEvtEnergy
      treeinfo_muoneff->SetAddress(&fMuonEff);
      treeoff_evtenergy->SetAddress(&fOffEvtEnergy);    
      treeoff_evtoffset->SetAddress(&fOffEvtOffset);
     
      // Je pense que premiere permet de recuperer efficacite du run et comme avec les lignes precedentes on l a relier a fMuonEff, fMuonEff=efficacité du run
      treeinfo_muoneff->GetEntry(0);
      fMuonEff=fMuonEff*100.;
      
      int index_eff = 0;
      for(int ieff=0; ieff<Nbands_eff;ieff++){
	if(fMuonEff>eff[ieff] && fMuonEff<eff[ieff+1]){
	    index_eff=ieff;
	    std::cout << ieff << " " << fMuonEff << endl;
	    break;
	  }
	  else if(fMuonEff>eff[Nbands_eff]){
	    std::cout << " efficacite superieur a efficacite max" << std::endl;
	    break;
	  }
      }
      distrib_zenith.Fill(zenit);
      Nrun_zenith->Fill(index);
      distrib_eff.Fill(fMuonEff);
      Nrun_eff->Fill(index_eff);
      std::cout << zenit << " " << index << endl;
      std::cout << fMuonEff << " " << index_eff <<endl;
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
	  std::cout << "il = " << il << " offset = " << myoff << " meanzenith = " << zenit << " fMuonEff = " << fMuonEff << " ethresh = " << ethresh << std::endl;
	}
      
      //Loop over the events 
      for (Int_t ievt=0;ievt<treeoff->GetEntries();++ievt) 
      {
	Nevent_zenith->Fill(index);
	Nevent_eff->Fill(index_eff);
	//Recupere pour chaque evenement energie et offset
	treeoff_evtenergy->GetEntry(ievt);    
	treeoff_evtoffset->GetEntry(ievt);
	// Ici fOffEvtOffset=a la valeur de l'offset qui vient d etre lu dans le TTree pour l evenement ievt.
	//Given a point fOffEvtOffset, approximates the value via linear interpolation based on the two nearest bin centers. Du coup vu que precedement pour chaque offset on a cree une valeur du ethreshol, pour chaque offset des evenements il peut remonter en interpolant entre les deux bins les plus proches, la valeur du ethreshold correspondant
	//seul truc je vois pas a quoi ca sert sachant qu'il a tout mis a zero avant dans l histo fLookupEthresh donc la valeur interpollé est aussi zero
	double ethresh = fLookupEthresh->Interpolate(fOffEvtOffset);
	distrib_energy.Fill(fOffEvtEnergy);
	//fEvtOffsetMax mis a 2.7 et defini au debut du programme
	Int_t index_E=0;
	for(int ie=0; ie<Nbands_E;ie++){
	  if(log10(fOffEvtEnergy)>E[ie] && log10(fOffEvtEnergy)<E[ie+1]){
	    index_E=ie;
	    std::cout << ie << " " << log10(fOffEvtEnergy) << endl;
	    break;
	       }
	}
	Nevent_E->Fill(index_E);
	Nevents_E_band_zen_eff[index][index_eff]->Fill(index_E);
	if ( (fOffEvtEnergy >= EnergieMin) && (fOffEvtEnergy < EnergieMax) && (fOffEvtOffset < fEvtOffsetMax ) && (fOffEvtEnergy >= ethresh) ) 
	  {
	    //la on rempli les deux vecteur d'histo quon avait cree precedment en offset**2 et offset de 700 et 250 bins et de valeur max fEvtOffsetMax**2 et fEvtOffsetMax* respectivement. index donne la bande en zenith dans laquelle on est, pour chaque histo de bande en zenoth, on rempli l'histogramme avec l'offset de chaque evenement. Donc pour chaque bande en zenith, on a pour les 750 ou 250 bandes offset defini pour les histo hists et histsR le nombre d evenement correspondant.On peut donc avoir les radial lookup
	    hists[index][index_eff][index_E]->Fill(fOffEvtOffset*fOffEvtOffset);
	    histsR[index][index_eff][index_E]->Fill(fOffEvtOffset);
	    std::cout << "Added the " << ievt << " Energy = " << fOffEvtEnergy << " Offset = " << fOffEvtOffset << " ethresh = " << ethresh << std::endl;
	  }
      }
      
      std::cout << "Run " << fRunList[irun] << " index : " << index <<  " Total Added Run = " << ++ntotruns << std::endl;
      
      delete fLookupEthresh;
      if (filehandler) delete filehandler;
      
      std::cout << "\033[1;32;40m" << "----------------> Run " << fRunList.at(irun) << " Sucessfully Added !  Total : " << ++ntotruns   << "\033[0m" << std::endl; 
      }
  
  double statistic;
  double threshold_stat=80;
  for (int izen=0;izen<Nbands;izen++) 
    {
      for (int ieff=0;ieff<Nbands_eff;ieff++) 
	{
	  for (int i_E=0;i_E<Nbands_E;i_E++) 
	    {
	      //GetBincontent va de 1 a Nbre max bin, pas de 0 a Nbre max bin-1 c'est pour ca qu'on met i_E+1
	      statistic=Nevents_E_band_zen_eff[izen][ieff]->GetBinContent(i_E+1);
		if(statistic< threshold_stat){
		  std::cout << "WARNING: statictic inferieur a" << threshold_stat << ": il y a" << statistic  << "evenements in the zenithal band:"<< zen[izen]<< "-" << zen[izen+1] <<"degrees, efficacite band: " << eff[ieff] << "-" << eff[ieff+1] << " for the band in energy number "<< i_E<< ": E:"<< TMath::Power(10,E[i_E]) << "-" << TMath::Power(10,E[i_E+1]) << " TeV" << endl;
		}

	    }
	}
    }
  TString foutname(outname);
  foutname+="_";
  foutname+=Config;
  foutname+=".root";
  //TCanvas *c1;
  //TCanvas *c2;
  TFile* out = new TFile(foutname,"recreate");
  out->cd();
  std::string outnamebis(outname);
  for (int index=0;index<Nbands;index++) 
    {
      for (int index_eff=0;index_eff<Nbands_eff;index_eff++) 
	{
	  Nevents_E_band_zen_eff[index][index_eff]->Write("",TObject::kOverwrite);
	  for (int index_E=0;index_E<Nbands_E;index_E++) 
	    {
	      std::cout << "########################## BEGINNING TO  TREAT THE BAND #######################################" << std::endl;
	      std::cout << "zenithal band:"<< zen[index]<< "-" << zen[index+1] <<"degrees, efficacite band: " << eff[index_eff] << "-" << eff[index_eff+1] << ", band in energy number "<< index_E<< ": E:"<< TMath::Power(10,E[index_E]) << "-" << TMath::Power(10,E[index_E+1]) << " TeV" << endl;
	       
	      if(hists[index][index_eff][index_E]) {
		//hits est un vecteur d'histogram donc .at() est une fonction C++ qui retourne une reference de l'element a la position index. Contrairement a juste l'operateur [], at verifie que ca depapsee la limit de taille du vector.
		std::string histsmoothname = hists[index][index_eff][index_E]->GetName();
		histsmoothname+="_Smooth";
		
		TH1F *histtmprb      =(TH1F*)hists[index][index_eff][index_E]->Clone();
		TH1F *histtmpfit     = (TH1F*)hists[index][index_eff][index_E]->Clone();
		TH1F *histtmpsm      = (TH1F*)histtmpfit->Rebin(5,"histtmpsm");
		histssmooth[index][index_eff][index_E] = histtmpsm;
		
		// VIM : Test : 
		histssmooth[index][index_eff][index_E]->SetTitle(histsmoothname.c_str());
		histssmooth[index][index_eff][index_E]->SetName(histsmoothname.c_str());
		//comprend pas tres bien ce que ca fait la fonction Smooth?
		histssmooth[index][index_eff][index_E]->Smooth(3);
		
 
		//Rebinning and/or smoothing for the fit !!! BE CARREFULL !!!!!      
		
		int RebinTotal=2; 
		if (index<=Nbands-2) {
		  //getminimum donne la valeur min de l'histo
		  double minimum_bin = histtmprb->GetMinimum();
		  while (minimum_bin<20 && RebinTotal<=32) {
		    histtmprb->Rebin(2);
		    minimum_bin=histtmprb->GetMinimum();
		    RebinTotal=RebinTotal*2;
		  }
		  RebinTotal=0.5*RebinTotal;
		  //si c'est pas pasé dans le while rebintotal=1 apres avoir mumtiplié par 0.5 donc il le remet a 2 pour rebinner un peu
		  if (RebinTotal==1) {RebinTotal=2;}
		}
		if (index==Nbands-1) {RebinTotal=64;}
		std::cout<<"RebinTotal="<<RebinTotal<<"\n";
		histtmpfit->Rebin(RebinTotal);
		
		//Rebinning histo per 10 -> for the fit 
		TString pols[8];
		
		pols[0]="pol4";pols[1]="pol5";pols[2]="pol6";pols[3]="pol7";pols[4]="pol8";pols[5]="pol9";pols[6]="pol10";
		pols[7]="pol11";
		double DevMin=1;
		int polmin=1;
		double DevTest=0;
		bool flag=0,flagTailP=0;
		double DevMinTailP=1;
		int polminTailP=1;
		TVirtualFitter::Fitter(histtmpfit)->SetMaxIterations(20000);
		if (index<=Nbands-1) {
		  DevMin=1;
		  polmin=1;
		  DevTest=0;
		  for (int poltest=0;poltest<=7;poltest++) {
		    TF1 *fitfunc = new TF1("mypol",pols[poltest],0.,TMath::Power(fEvtOffsetMax,2.));
		    //dit que en x=0, le polynome=valeur de l'histo en zero. Donne une condition pour le fit
		    fitfunc->SetParameter(0,hists[index][index_eff][index_E]->GetBinContent(1));
		    histtmpfit->Fit("mypol","RINQ","",0.,TMath::Power(fEvtOffsetMax,2.));
		    int NDF=fitfunc->GetNDF();
		    double chi2=fitfunc->GetChisquare();
		    DevTest=TMath::Abs((chi2/NDF)-1);
		    //Dit que le nombre de point ou sera calculé les valeurs de la fonction c'est 700=nbre de barre de l'histogramme non rebinné. Il doit faire ca car ensuite il calcule le minimum de la fonction fitmin donc ce minimum est caculé entre 0 et 2.7**2 divisé en 700 intervalles
		    fitfunc->SetNpx(700);
		    std::cout<<"poltest  "<<pols[poltest]<<"  DevTest= "<<DevTest<<"NDF= "<<NDF<<" chi2 = "<< fitfunc->GetChisquare() <<"\t";
		    double fitmin=fitfunc->GetMinimum();
		    //double x_minimum=fitfunc->GetXaxis()->GetBinCenter(fitmin);
		    // cout<<"x_mim=  "<<x_minimum<<"\n"; 
		    
		    int MinPos=histtmpfit->GetMinimumBin();
		    //PErmet de savoir vu qu'on a rebinné l'histo sous certaine condition combien il reste de bin
		    int Nx=histtmpfit->GetNbinsX();
		    //Minposval lui donne la valeur en offser**2 du bon qui correspond au minimum de l histo
		    double MinPosVal=MinPos*2.7*2.7/Nx;
		    
		    if (DevTest<DevMin)
		      {
			if (fitmin>=0 && NDF>1)
			  { 		
			    //	  if (x_minimum > 6.25) 
			    {
			      DevMin=DevTest;
			      polmin=poltest;
			      flag=1;
			      std::cout<<"\n";
			    }
			  }
		      }
		  }
		  std::cout<<"polselect  "<<pols[polmin]<<"\n";
		  TF1 *fitfunc = new TF1("mypol",pols[polmin],0.,TMath::Power(fEvtOffsetMax,2.));
		  fitfunc->SetParameter(0,hists[index][index_eff][index_E]->GetBinContent(1));
		  histtmpfit->Fit("mypol","RINQ","",0.,TMath::Power(fEvtOffsetMax,2.));
		  //histtmpfit est l'histogramme rebinné sur lequel on fit, ensuite apres avoir determiné le meilleur polynome on remplie histsfit[index] de la fonction fitté sur 700 points.
		  std::cout << "FITTING RESULTS FOR THE zenithal band:"<< zen[index]<< "-" << zen[index+1] <<"degrees, efficacite band: " << eff[index_eff] << "-" << eff[index_eff+1] << ", band in energy number "<< index_E<< ": E:"<< TMath::Power(10,E[index_E]) << "-" << TMath::Power(10,E[index_E+1]) << " TeV" << endl;
		  std::cout << "probfit = " << fitfunc->GetProb() << " chi2 = " << fitfunc->GetChisquare() << " ndf = " << fitfunc->GetNDF() << std::endl;
		  std::cout << "before set 700 pix" << std::endl;
		  //dis que la fonction polynome choisie est calculé sur 700 points et du coup histsfit[index] est un histogramme de 700 intervalles correspondants a la fonction fittée
		  fitfunc->SetNpx(700);
		  histsfit[index][index_eff][index_E] = new TH1F(*((TH1F*)fitfunc->GetHistogram()));
		  
		  /*if (index==Nbands-1)
		    {
		      histsfit[index][index_eff][index_E]->Scale(0);
		      std::cout<<"----use previous fitted zenith part for this zenith part---------------------------"<<"\n";
		      
		      for (int ient=0; ient<=histsfit[index-1][index_eff][index_E]->GetNbinsX();ient++) {
			double yval=histsfit[index-1][index_eff][index_E]->GetBinContent(ient);
			double xval=histsfit[index-1][index_eff][index_E]->GetBinCenter(ient);
			histsfit[index][index_eff][index_E]->Fill(xval,yval);
		      }
		      
		      }*/
		  double scaling;
		  
		  if (index!=10) {
		    scaling=hists[index][index_eff][index_E] ->GetSum()/ histsfit[index][index_eff][index_E]->GetSum();
		  }
		  
		  histsfit[index][index_eff][index_E]->Scale(scaling);
		  double mxm=hists[index][index_eff][index_E]->GetMaximum();
		  double mxmft=histsfit[index][index_eff][index_E]->GetMaximum();
		  
		  std::cout<<"scaling  "<<scaling<<"\t"<<mxm<<"\t"<<mxmft<<"\n";
		  
		  std::cout << "After Getting Hist" << histsfit[index][index_eff][index_E]  << std::endl;
		  
		  std::string histsfitname = hists[index][index_eff][index_E]->GetName();
		  histsfitname+="_Fit";
		  
		  histsfit[index][index_eff][index_E]->SetName(histsfitname.c_str());
		  histsfit[index][index_eff][index_E]->SetTitle(histsfitname.c_str());
		  
		  
		  
		  std::cout << "writting hists and histR" << std::endl;
		  hists[index][index_eff][index_E]->Write("",TObject::kOverwrite);	
		  histsR[index][index_eff][index_E]->Write("",TObject::kOverwrite);
		  std::cout << "writting histssmooth" << std::endl; 
		  histssmooth[index][index_eff][index_E]->Write("",TObject::kOverwrite);
		  std::cout << "writting histsfit" << std::endl; 
		  histsfit[index][index_eff][index_E]->Write("",TObject::kOverwrite);	
		  
		  
		  //string histname;
		  //string histnameR;
		  //histname=hists[izen][ieff][i_E]->GetName();
		  //histnameR=histsR[izen][ieff][i_E]->GetName();
		  
		
		}		
	      }	      
	    }   
	}
    }
  //delete c1;
  //delete c2;
  distrib_zenith.Write();
  distrib_eff.Write();
  distrib_energy.Write();
  Nrun_zenith->Write();
  Nevent_zenith->Write();
  Nrun_eff->Write();
  Nevent_eff->Write();
  Nevent_E->Write();
  out->Close();
  delete out;
  
  
}

