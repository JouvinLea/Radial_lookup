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
  Int_t Nbands = 6;
  double zen1 = 0.;
  double zen2 = 20.;
  double zen3 = 30.;
  double zen4 = 40.;
  double zen5 = 45.;
  double zen6 = 55.;
  
  std::vector<std::string> histnames;
  histnames.clear();
  histnames.push_back(std::string("Zen0_20deg"));
  histnames.push_back(std::string("Zen20_30deg"));
  histnames.push_back(std::string("Zen30_40deg"));
  histnames.push_back(std::string("Zen40_45deg"));
  histnames.push_back(std::string("Zen45_55deg"));
  histnames.push_back(std::string("Zen55_90deg"));
  
  std::string radname = "RadialLookup_";
  std::string radnameR = "RadialLookup_R"; 

  std::vector<TH1F*> hists;
  hists.resize(Nbands);
  
  std::vector<TH1F*> histsR;
  histsR.resize(Nbands);
  
  std::vector<TH1F*> histssmooth;
  histssmooth.resize(Nbands);
  
  std::vector<TH1F*> histsfit;
  histsfit.resize(Nbands);
  
  TH1F *Nrun;
  TH1F *Nevent;
  Nrun=new TH1F("Nrun_bin_zenith","Nrun_bin_zenith",6,0,6); 
  Nevent=new TH1F("Nevent_bin_zenith","Nrun_event_zenith",6,0,6);
  // Initializing Tables : 
  //loop sur les N bands en zenith et remplie les histogrammes qui ont le nom de radia_lookup+_nom de chaque band en zenith.
  for (int itab=0;itab<Nbands;itab++) 
    {
      std::string tabname;
      std::string tabnameR;
      
      tabname=radname+histnames.at(itab);
      tabnameR=radnameR+histnames.at(itab);
      //hists et histsR c'est des vecteurs d'histogram de root. 
      //fEvtOffsetMax: defini au debut du program, donne l'offset max ici a 2.7
      //Pour chaque band vecteur donc chaque band de zenith cree un histos de 700 bins de O a fEvtOffsetMax**2:histogramme des theta**2
      //et un histogram de 250 bins de 0 a fEvtOffsetMax: histogram des theta
      hists[itab] = new TH1F(tabname.c_str(),tabname.c_str(),700,0,TMath::Power(fEvtOffsetMax,2.));
      histsR[itab] = new TH1F(tabnameR.c_str(),tabnameR.c_str(),250,0,fEvtOffsetMax);
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
    
      if((zenit) >= zen1 && (zenit) < zen2) index = 0;
      else if((zenit) >= zen2 && (zenit) < zen3) index = 1;
      else if((zenit) >= zen3 && (zenit) < zen4) index = 2;
      else if((zenit) >= zen4 && (zenit) < zen5) index = 3;
      else if((zenit) >= zen5 && (zenit) < zen6) index = 4;
      else if((zenit) >= zen6) index = 5;
      else 
	{
	  std::cout << " ERROR!" << std::endl;
	  continue;
	}

      Nrun->Fill(index);
      //take each events, read leaves the and project it in a file
      // la lis les deux TTree present dans le TFile de nom RunInfoTree et EventsTree_BgMakerOff et les recupere avec get
      TTree *treeinfo = NULL;
      treeinfo = (TTree*)filehandler->Get("RunInfoTree");
      TTree *treeoff = NULL;
      treeoff = (TTree*)filehandler->Get("EventsTree_BgMakerOff");
    
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
	Nevent->Fill(index);
	//Recupere pour chaque evenement energie et offset
	treeoff_evtenergy->GetEntry(ievt);    
	treeoff_evtoffset->GetEntry(ievt);
	// Ici fOffEvtOffset=a la valeur de l'offset qui vient d etre lu dans le TTree pour l evenement ievt.
	//Given a point fOffEvtOffset, approximates the value via linear interpolation based on the two nearest bin centers. Du coup vu que precedement pour chaque offset on a cree une valeur du ethreshol, pour chaque offset des evenements il peut remonter en interpolant entre les deux bins les plus proches, la valeur du ethreshold correspondant
	//seul truc je vois pas a quoi ca sert sachant qu'il a tout mis a zero avant dans l histo fLookupEthresh donc la valeur interpollé est aussi zero
	double ethresh = fLookupEthresh->Interpolate(fOffEvtOffset);
	//fEvtOffsetMax mis a 2.7 et defini au debut du programme
	if ( (fOffEvtEnergy >= EnergieMin) && (fOffEvtEnergy < EnergieMax) && (fOffEvtOffset < fEvtOffsetMax ) && (fOffEvtEnergy >= ethresh) ) 
	  {
	    //la on rempli les deux vecteur d'histo quon avait cree precedment en offset**2 et offset de 700 et 250 bins et de valeur max fEvtOffsetMax**2 et fEvtOffsetMax* respectivement. index donne la bande en zenith dans laquelle on est, pour chaque histo de bande en zenoth, on rempli l'histogramme avec l'offset de chaque evenement. Donc pour chaque bande en zenith, on a pour les 750 ou 250 bandes offset defini pour les histo hists et histsR le nombre d evenement correspondant.On peut donc avoir les radial lookup
	    hists[index]->Fill(fOffEvtOffset*fOffEvtOffset);
	    histsR[index]->Fill(fOffEvtOffset);
	    std::cout << "Added the " << ievt << " Energy = " << fOffEvtEnergy << " Offset = " << fOffEvtOffset << " ethresh = " << ethresh << std::endl;
	  }
      }
      
      std::cout << "Run " << fRunList[irun] << " index : " << index <<  " Total Added Run = " << ++ntotruns << std::endl;
      
      delete fLookupEthresh;
      if (filehandler) delete filehandler;
      
      std::cout << "\033[1;32;40m" << "----------------> Run " << fRunList.at(irun) << " Sucessfully Added !  Total : " << ++ntotruns   << "\033[0m" << std::endl; 
    }
  
  
  TString foutname(outname);
  foutname+="_";
  foutname+=Config;
  foutname+=".root";
  
  TFile* out = new TFile(foutname,"recreate");
  out->cd();
  
  for(unsigned int index = 0; index < hists.size(); ++index){
    std::cout << "########################## BEGINNING TO  TREAT THE BAND " << histnames.at(index) << " #######################################" << std::endl;
    if(hists[index]) {
      //hits est un vecteur d'histogram donc .at() est une fonction C++ qui retourne une reference de l'element a la position index. Contrairement a juste l'operateur [], at verifie que ca depapsee la limit de taille du vector.
      std::string histsmoothname = (hists.at(index))->GetName();
      histsmoothname+="_Smooth";
      
      TH1F *histtmprb      =(TH1F*)hists[index]->Clone();
      TH1F *histtmpfit     = (TH1F*)hists[index]->Clone();
      TH1F *histtmpsm      = (TH1F*)histtmpfit->Rebin(5,"histtmpsm");
      histssmooth.at(index) = histtmpsm;
      
      // VIM : Test : 
      (histssmooth.at(index))->SetTitle(histsmoothname.c_str());
      (histssmooth.at(index))->SetName(histsmoothname.c_str());
      //comprend pas tres bien ce que ca fait la fonction Smooth?
      (histssmooth.at(index))->Smooth(3);
      
 
      //Rebinning and/or smoothing for the fit !!! BE CARREFULL !!!!!      
      int RebinTotal=2; 
      if (index<=4) {
	
	double minimum_bin = histtmprb->GetMinimum();
	while (minimum_bin<20 && RebinTotal<=32) {
	  histtmprb->Rebin(2);
	  minimum_bin=histtmprb->GetMinimum();
	  RebinTotal=RebinTotal*2;
	}
	RebinTotal=0.5*RebinTotal;
	if (RebinTotal==1) {RebinTotal=2;}
      }
      if (index==5) {RebinTotal=64;}
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
      if (index<=5) {
	DevMin=1;
	polmin=1;
	DevTest=0;
	for (int poltest=0;poltest<=7;poltest++) {
	  TF1 *fitfunc = new TF1("mypol",pols[poltest],0.,TMath::Power(fEvtOffsetMax,2.));
	  fitfunc->SetParameter(0,hists[index]->GetBinContent(1));
	  histtmpfit->Fit("mypol","RINQ","",0.,TMath::Power(fEvtOffsetMax,2.));
	  int NDF=fitfunc->GetNDF();
	  double chi2=fitfunc->GetChisquare();
	  DevTest=TMath::Abs((chi2/NDF)-1);
	  fitfunc->SetNpx(700);
	  std::cout<<"poltest  "<<pols[poltest]<<"  DevTest= "<<DevTest<<"NDF= "<<NDF<<" chi2 = "<< fitfunc->GetChisquare() <<"\t";
	  double fitmin=fitfunc->GetMinimum();
	  //double x_minimum=fitfunc->GetXaxis()->GetBinCenter(fitmin);
	  // cout<<"x_mim=  "<<x_minimum<<"\n"; 
	  int MinPos=histtmpfit->GetMinimumBin();

	  int Nx=histtmpfit->GetNbinsX();
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
	fitfunc->SetParameter(0,hists[index]->GetBinContent(1));
	histtmpfit->Fit("mypol","RINQ","",0.,TMath::Power(fEvtOffsetMax,2.));
	std::cout << "FITTING RESULTS FOR THE " << histnames.at(index) << " BANDS :" << std::endl; 
	std::cout << "probfit = " << fitfunc->GetProb() << " chi2 = " << fitfunc->GetChisquare() << " ndf = " << fitfunc->GetNDF() << std::endl;
 	std::cout << "before set 700 pix" << std::endl;
	fitfunc->SetNpx(700);
	histsfit[index] = new TH1F(*((TH1F*)fitfunc->GetHistogram()));

	if (index==5)
	  {
	    histsfit[index]->Scale(0);
	    std::cout<<"----use previous fitted zenith part for this zenith part---------------------------"<<"\n";
	    
	    for (int ient=0; ient<=histsfit[index-1]->GetNbinsX();ient++) {
	      double yval=histsfit[index-1]->GetBinContent(ient);
	      double xval=histsfit[index-1]->GetBinCenter(ient);
	      histsfit[index]->Fill(xval,yval);
	   }
	    
	  }
	double scaling;
	
	if (index!=10) {
	  scaling=hists[index]->GetSum()/ histsfit[index]->GetSum();
	}

	histsfit[index]->Scale(scaling);
	double mxm=hists[index]->GetMaximum();
	double mxmft=histsfit[index]->GetMaximum();
	
	std::cout<<"scaling  "<<scaling<<"\t"<<mxm<<"\t"<<mxmft<<"\n";
	
	std::cout << "After Getting Hist" << histsfit[index] << std::endl;
	
	std::string histsfitname = (hists.at(index))->GetName();
	histsfitname+="_Fit";
	
	histsfit[index]->SetName(histsfitname.c_str());
	histsfit[index]->SetTitle(histsfitname.c_str());
	

	
	std::cout << "writting hists" << std::endl; 
	(hists.at(index))->Write("",TObject::kOverwrite);
	std::cout << "writting histssmooth" << std::endl; 
	(histssmooth.at(index))->Write("",TObject::kOverwrite);
	std::cout << "writting histsfit" << std::endl; 
	(histsfit.at(index))->Write("",TObject::kOverwrite);	
	(histsR.at(index))->Write("",TObject::kOverwrite);
      }
     
    }
  }
  Nrun->Write();
  Nevent->Write();
  out->Close();
  delete out;
 
  
}

