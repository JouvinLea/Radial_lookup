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
double fEvtEnergyMin=0.2;

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
  Int_t Nzen=7;
  Int_t Nbands = Nzen-1;
  double zen[Nzen];
  zen[0] = 0.;
  zen[1] = 20.;
  zen[2] = 30.;
  zen[3] = 40.;
  zen[4] = 45.;
  zen[5] = 55.;
  zen[6] = 90.;
  
  Int_t Neff=4;
  Int_t Nbands_eff = Neff-1;
  double eff[Neff];
  eff[0] = 50.;
  eff[1] = 60.;
  eff[2] = 65.;
  eff[3] = 100.;

  Int_t Nbands_E=16;
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

  std::vector<std::vector<TH1F*> > hists;
  std::vector<std::vector<TH1F*> > histsR;
  //std::vector<std::vector<TH1F*> > histssmooth;
  //std::vector<std::vector<TH1F*> > histsfit;
  hists.resize(Nbands);
  histsR.resize(Nbands);
  //histssmooth.resize(Nbands);
  //histsfit.resize(Nbands);
  
  for (int i_zen=0;i_zen<Nbands;i_zen++){

    hists[i_zen].resize(Nbands_eff);
    histsR[i_zen].resize(Nbands_eff);
    //histssmooth[i_zen].resize(Nbands_eff);
    //histsfit[i_zen].resize(Nbands_eff);

  }
  
  
  TH1F *Nrun_zenith;
  TH1F *Nevent_zenith;
  Nrun_zenith=new TH1F("Nrun_bin_zenith","Nrun_bin_zenith",6,0,6); 
  Nevent_zenith=new TH1F("Nevent_bin_zenith","Nrun_event_zenith",6,0,6);
  TH1F *Nrun_eff;
  TH1F *Nevent_eff;
  Nrun_eff=new TH1F("Nrun_bin_efficacite","Nrun_bin_efficacite",3,0,3); 
  Nevent_eff=new TH1F("Nevent_bin_efficacite","Nrun_event_efficacite",3,0,3);
  // Initializing Tables : 
  //loop sur les N bands en zenith et remplie les histogrammes qui ont le nom de radia_lookup+_nom de chaque band en zenith.
  for (int izen=0;izen<Nbands;izen++) 
    {
      for (int ieff=0;ieff<Nbands_eff;ieff++) {
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
	
      tabname=radname+"_zen"+zenmin+"_"+zenmax+"deg_eff"+effmin+"_"+effmax;
      tabnameR=radnameR+"_zen"+zenmin+"_"+zenmax+"deg_eff"+effmin+"_"+effmax;
      //hists et histsR c'est des vecteurs d'histogram de root. 
      //fEvtOffsetMax: defini au debut du program, donne l'offset max ici a 2.7
      //Pour chaque band vecteur donc chaque band de zenith cree un histos de 700 bins de O a fEvtOffsetMax**2:histogramme des theta**2
      //et un histogram de 250 bins de 0 a fEvtOffsetMax: histogram des theta
      hists[izen][ieff] = new TH1F(tabname.c_str(),tabname.c_str(),700,0,TMath::Power(fEvtOffsetMax,2.));
      histsR[izen][ieff] = new TH1F(tabnameR.c_str(),tabnameR.c_str(),250,0,fEvtOffsetMax);
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
    
      if((zenit) >= zen[0] && (zenit) < zen[1]) index = 0;
      else if((zenit) >= zen[1] && (zenit) < zen[2]) index = 1;
      else if((zenit) >= zen[2] && (zenit) < zen[3]) index = 2;
      else if((zenit) >= zen[3] && (zenit) < zen[4]) index = 3;
      else if((zenit) >= zen[4] && (zenit) < zen[5]) index = 4;
      else if((zenit) >= zen[5]) index = 5;
      else 
	{
	  std::cout << " ERROR!" << std::endl;
	  continue;
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
    
      if((fMuonEff) >= eff[0] && (fMuonEff) < eff[1]) index_eff  = 0;
      else if((fMuonEff) >= eff[1] && (fMuonEff) < eff[2]) index_eff  = 1;
      else if((fMuonEff) >= eff[2]) index_eff = 2;
      else 
	{
	  std::cout << " ERROR!" << std::endl;
	  continue;
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
	       }
	}
	if ( (fOffEvtEnergy >= EnergieMin) && (fOffEvtEnergy < EnergieMax) && (fOffEvtOffset < fEvtOffsetMax ) && (fOffEvtEnergy >= ethresh) ) 
	  {
	    //la on rempli les deux vecteur d'histo quon avait cree precedment en offset**2 et offset de 700 et 250 bins et de valeur max fEvtOffsetMax**2 et fEvtOffsetMax* respectivement. index donne la bande en zenith dans laquelle on est, pour chaque histo de bande en zenoth, on rempli l'histogramme avec l'offset de chaque evenement. Donc pour chaque bande en zenith, on a pour les 750 ou 250 bandes offset defini pour les histo hists et histsR le nombre d evenement correspondant.On peut donc avoir les radial lookup
	    hists[index][index_eff]->Fill(fOffEvtOffset*fOffEvtOffset);
	    histsR[index][index_eff]->Fill(fOffEvtOffset);
	    std::cout << "Added the " << ievt << " Energy = " << fOffEvtEnergy << " Offset = " << fOffEvtOffset << " ethresh = " << ethresh << std::endl;
	  }
      }
      
      std::cout << "Run " << fRunList[irun] << " index : " << index <<  " Total Added Run = " << ++ntotruns << std::endl;
      
      delete fLookupEthresh;
      if (filehandler) delete filehandler;
      
      std::cout << "\033[1;32;40m" << "----------------> Run " << fRunList.at(irun) << " Sucessfully Added !  Total : " << ++ntotruns   << "\033[0m" << std::endl; 
    }
  
  
  /*TString foutname(outname);
  foutname+="_";
  foutname+=Config;
  foutname+=".root";
  TCanvas *c1;
  TCanvas *c2;
  TFile* out = new TFile(foutname,"recreate");
  out->cd();
  std::string outnamebis(outname);
  for (int izen=0;izen<Nbands;izen++) 
    {
      for (int ieff=0;ieff<Nbands_eff;ieff++) {
	string histname;
	string histnameR;
	histname=hists[izen][ieff]->GetName();
	histnameR=histsR[izen][ieff]->GetName();
	hists[izen][ieff]->Write("",TObject::kOverwrite);	
	histsR[izen][ieff]->Write("",TObject::kOverwrite);
	string name_canvas1;
	string name_canvas2;
	name_canvas1=outnamebis+"_"+histname+".jpg";
	name_canvas2=outnamebis+"_"+histnameR+".jpg";
	c1=new TCanvas("theta**2");
	hists[izen][ieff]->Draw();
	c1->SaveAs(name_canvas1.c_str());
	c2=new TCanvas("theta");
	histsR[izen][ieff]->Draw();
	c2->SaveAs(name_canvas2.c_str());
    }
  }
  delete c1;
  delete c2;
  distrib_zenith.Write();
  distrib_eff.Write();
  distrib_energy.Write();
  Nrun_zenith->Write();
  Nevent_zenith->Write();
  Nrun_eff->Write();
  Nevent_eff->Write();
  out->Close();
  delete out;*/
 
  
}

