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
#include <TNtuple>

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
  eff[0] = 0.;
  eff[1] = 30.;
  eff[2] = 60.;
  eff[3] = 100.;
  
  
  TNtuple distrib_zenith("distrib_zenith", "distribution des zenith des runs", "zenith");
  TNtuple distrib_eff("dsitrib_efficacite", "distribution des efficacite des runs", "efficacite");
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
  
  
  TH1F *Nrun;
  TH1F *Nevent;
  Nrun=new TH1F("Nrun_bin_zenith","Nrun_bin_zenith",6,0,6); 
  Nevent=new TH1F("Nevent_bin_zenith","Nrun_event_zenith",6,0,6);
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
      effmin=str_convert3.str();histsnameR+".jpg";
      hists[index][index_eff]->Write("",TObject::kOverwrite);	
      histsR[index][index_eff]->Write("",TObject::kOverwrite);
      c1=new TCanvas("theta**2");
      hists[index][index_eff]->Draw();
      c1->SaveAs(name_canvas1.c_str());
      c2=new TCanvas("theta");
      histsR[index][index_eff]->Draw();
      c2->SaveAs(name_canvas2.c_str());
    }
  }
  delete c1;
  delete c2;
  /*for(unsigned int index = 0; index < Nbands; ++index){
    for(unsigned int index_eff = 0; index < Nbands_eff; ++index_eff){
      std::cout << "########################## BEGINNING TO  TREAT THE BAND: Zenith: " << zen[index] << "-"<< zen[index+1] << ", Efficacite:" << eff[index_eff] << "-" << eff[index_eff+1] << " #######################################" << std::endl;
      if(hists[index][index_eff]) {
	//hits est un vecteur d'histogram donc .at() est une fonction C++ qui retourne une reference de l'element a la position index. Contrairement a juste l'operateur [], at verifie que ca depapsee la limit de taille du vector.
	std::string histsmoothname = (hists[index][index_eff])->GetName();
	histsmoothname+="_Smooth";
	
	TH1F *histtmprb      =(TH1F*)hists[index][index_eff]->Clone();
	TH1F *histtmpfit     = (TH1F*)hists[index][index_eff]->Clone();
	TH1F *histtmpsm      = (TH1F*)histtmpfit->Rebin(5,"histtmpsm");
	histssmooth[index][index_eff]= histtmpsm;
	
      // VIM : Test : 
	histssmooth[index][index_eff]->SetTitle(histsmoothname.c_str());
	histssmooth[index][index_eff]->SetName(histsmoothname.c_str());
	//comprend pas tres bien ce que ca fait la fonction Smooth?
	histssmooth[index][index_eff]->Smooth(3);
      
 
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
	    fitfunc->SetParameter(0,hists[index][index_eff]->GetBinContent(1));
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
	  fitfunc->SetParameter(0,hists[index][index_eff]->GetBinContent(1));
	  histtmpfit->Fit("mypol","RINQ","",0.,TMath::Power(fEvtOffsetMax,2.));
	  std::cout << "FITTING RESULTS FOR THE " << histnames[index][index_eff]<< " BANDS :" << std::endl; 
	  std::cout << "probfit = " << fitfunc->GetProb() << " chi2 = " << fitfunc->GetChisquare() << " ndf = " << fitfunc->GetNDF() << std::endl;
	  std::cout << "before set 700 pix" << std::endl;
	  fitfunc->SetNpx(700);
	  histsfit[index][index_eff] = new TH1F(*((TH1F*)fitfunc->GetHistogram()));
	  
	  //Comprend pas tres bien ce qu'il fait la pour index=5
	  if (index==5)
	    {
	      histsfit[index][index_eff]->Scale(0);
	      std::cout<<"----use previous fitted zenith part for this zenith part---------------------------"<<"\n";
	      
	      for (int ient=0; ient<=histsfit[index-1][index_eff-1]->GetNbinsX();ient++) {
		double yval=histsfit[index-1][index_eff-1]->GetBinContent(ient);
		double xval=histsfit[index-1][index_eff-1]->GetBinCenter(ient);
		histsfit[index][index_eff]->Fill(xval,yval);
	      }
	      
	    }
	  double scaling;
	  
	  if (index!=10) {
	    scaling=hists[index][index_eff]->GetSum()/ histsfit[index][index_eff]->GetSum();
	  }
	  
	  histsfit[index][index_eff]->Scale(scaling);
	  double mxm=hists[index][index_eff]->GetMaximum();
	  double mxmft=histsfit[index][index_eff]->GetMaximum();
	  
	  std::cout<<"scaling  "<<scaling<<"\t"<<mxm<<"\t"<<mxmft<<"\n";
	  
	  std::cout << "After Getting Hist" << histsfit[index][index_eff] << std::endl;
	  
	  std::string histsfitname = hists.[index][index_eff]->GetName();
	  histsfitname+="_Fit";
	  
	  histsfit[index][index_eff]->SetName(histsfitname.c_str());
	  histsfit[index][index_eff]->SetTitle(histsfitname.c_str());
	  
	  
	  
	  std::cout << "writting hists" << std::endl; 
	  hists[index][index_eff]->Write("",TObject::kOverwrite);
	  std::cout << "writting histssmooth" << std::endl; 
	  histssmooth[index][index_eff]->Write("",TObject::kOverwrite);
	  std::cout << "writting histsfit" << std::endl; 
	  histsfit[index][index_eff]->Write("",TObject::kOverwrite);	
	  histsR[index][index_eff]->Write("",TObject::kOverwrite);
	}
	
      }
    }
    }*/
  distrib_zenith->Write();
  distrib_eff->Write();
  distrib_eff->Write();
  Nrun->Write();
  Nevent->Write();
  out->Close();
  delete out;
 
  
}

