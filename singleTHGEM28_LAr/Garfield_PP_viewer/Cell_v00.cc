#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
//#include <curses.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#if defined(__WIN32__)
#define NOMINMAX
#include <windows.h>
#include <direct.h>
#else
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#endif

#include <TCanvas.h>
#include <TROOT.h>
#include <TApplication.h>
#include <TH1F.h>
#include <TH2D.h>

#include "MediumMagboltz.hh"
#include "SolidBox.hh"
#include "GeometrySimple.hh"
#include "ComponentElmer.hh"
#include "ViewGeometry.hh"
#include "ViewField.hh"
#include "ViewDrift.hh"
#include "ViewFEMesh_My.hh"
#include "ComponentConstant.hh"
#include "Sensor.hh"
#include "TrackHeed.hh"
#include "Plotting.hh"
#include "AvalancheMC_My.hh"
#include "LAr_Mat.hh"
#include "Mat_Vd_const.hh"
#include "Random.hh"
#include "GarfieldConstants.hh"

using namespace Garfield;

//Garfield sizes are in cm
#define _VISUAL 0
#define _X_SIZE 0.0225
#define _Y_SIZE 0.03897114317
#define MESH_ std::string("../v00.01_THGEM1/")
#define RESULT_ std::string("../Elmer_v00.01/case_v01.result")

void ensure_file(std::string fname); //makes sure file can be created later on
void ensure_folder(std::string folder);

void ensure_file(std::string fname)
{
	std::string folder = fname;
	while ((folder.back() != '\\') &&(folder.back()!='/') &&!folder.empty())
		folder.pop_back();
	if (!folder.empty())
		folder.pop_back();
	ensure_folder(folder);
}

void ensure_folder(std::string folder)
{
#if defined(__WIN32__)
	if (!folder.empty()) {
		DWORD ftyp = GetFileAttributesA(folder.c_str());
		if (!(ftyp & FILE_ATTRIBUTE_DIRECTORY) || ftyp == INVALID_FILE_ATTRIBUTES) {
			int code = system(("mkdir \"" + folder + "\"").c_str());
			if (code)
				std::cout << "mkdir error: " << GetLastError() << std::endl;
		}
	}
#else
	struct stat st;
	stat(folder.c_str(), &st);
	if (!S_ISDIR(st.st_mode)) {
		int code = system(("mkdir \"" + folder + "\"").c_str());
		if (code)
			std::cout << "mkdir error: " << code << std::endl;
	}
#endif //_WIN32__
}

void plot_field(ComponentFieldMap* fm, std::string filename, double x0 ,double y0,double z0 ,double x1,double y1 ,double z1, int Num,std::string name="", double L_fine=0, int Num_fine=0)
{
  ensure_file(filename);
  double dphi=0;
  double L = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0));
  double s=0,x=x0,y=y0,z=z0;
  double ex,ey,ez,ee;
  Medium* m;
  int status;
  int LN=(Num-Num_fine)/2.0;
  int RN=Num-LN;
  double ds=0;
  if ((L_fine<=0)||(Num_fine>=Num)||(L_fine>=s)||(RN<LN)||(Num_fine<=0)||(LN>Num)||(LN<0))
  {
    Num_fine=0;
    L_fine=0;
    LN=Num;
    RN=-Num;
  }
  std::ofstream flll(filename,std::ios_base::out);
	//flll<<"//z\tx\ty\tL\tEabs\tEx\tEy\tEz"
  for (int hh=0;hh<=Num;hh++) //hh<=Num in order to cover [x0;x1], not [x0;x1)
  {
     x=(x1-x0)*s/L+x0;
     y=(y1-y0)*s/L+y0;
     z=(z1-z0)*s/L+z0;
     fm->ElectricField(x,y,z,ex,ey,ez,m,status);
     ee=std::sqrt(ex*ex+ey*ey+ez*ez);
     if (ee>0) //there was a bug: because of -0.0 comparison in elmer component
     //resulting in appearing of 0 field in arbitrary points - fixed by me, yet better to leave this condition
     {
        if((x0==x1)&&(y0==y1))
           flll<<z<<'\t';
        else if ((z0==z1)&&(y0==y1))
           flll<<x<<'\t';
        else if ((x0==x1)&&(z0==z1))
           flll<<y<<'\t';
        else
           flll<<s<<'\t';
        flll<<ee<<'\t'<<ex<<'\t'<<ey<<'\t'<<ez<<std::endl;
     }
     if ((hh<LN)||(hh>RN))
     {
        ds=(L-L_fine)/(Num-LN-RN-1);
        dphi+=(ds/L)*(ex*(x1-x0)+ey*(y1-y0)+ez*(z1-z0));
        s+=ds;
     }
     else
     {
        ds=L_fine/(RN-LN+1);
        dphi+=ex*ds*(x1-x0)/L+ey*ds*(y1-y0)/L+ez*ds*(z1-z0)/L;
        s+=ds;
     }
  }
  if (name.size())
      std::cout<<"Int[E(r)*dr] on "<<name<<" = "<<dphi<<std::endl;
  flll.close();
}

int main(int argc, char * argv[]) {

  TApplication app("app", &argc, argv);
  plottingEngine.SetDefaultStyle();
  ComponentElmer* fm = new ComponentElmer();
  bool is_loaded = fm->Initialise(MESH_+"mesh.header", MESH_+"mesh.elements",MESH_+ "mesh.nodes",MESH_+ "diels.dat",RESULT_,"mm");
  if (is_loaded)
  	std::cout<<"load sucess"<<std::endl;
  else
  {
	std::cout<<"load failure, quiting"<<std::endl;
        return 0;
  }
  //fm->EnableXHexPeriodicity();
  fm->EnableMirrorPeriodicityX();
  fm->EnableMirrorPeriodicityY();
  for (int ff=0;ff<fm->GetNumberOfMaterials();ff++)
     fm->DriftMedium(ff);
    // Make a medium
  Mat_Vd_const* LAr = new Mat_Vd_const(5e4);// LAr_Mat("LAr_drift_data.txt", "LAr_diffusionL_data.txt","LAr_diffusionT_data.txt");
  //LAr->SetComposition("ar", 100.);
  LAr->SetTemperature(89);
  LAr->EnableDrift();
	LAr->SetPressure(760.);

  MediumMagboltz* met = new MediumMagboltz();
  met->SetComposition("co2", 100.);
  met->SetTemperature(293.15);
  met->SetPressure(760.);
  met->DisableDrift();
  MediumMagboltz* diel = new MediumMagboltz();
  diel->SetComposition("co2", 20.,"h2",80.);
  diel->SetTemperature(293.15);
  diel->SetPressure(760.);
  diel->DisableDrift();
  std::cout<<"Number of materials: "<<fm->GetNumberOfMedia()<<std::endl;

  fm->SetMedium(0,LAr);
  fm->SetMedium(1,diel);  //FR4
  fm->SetMedium(2,met);
  fm->SetMedium(3,met);
   // Make a sensor
  Sensor* sensor = new Sensor();
  sensor->AddComponent(fm);
  sensor->SetArea(-5.5*_X_SIZE,-5.5*_Y_SIZE,-0.32, 5.5*_X_SIZE, 5.5*_Y_SIZE, 0.32);
  TCanvas* c3 = new TCanvas("drift lines xz", "drift lines xz");

  plot_field(fm,"v00.01.01/center_axis_field.txt",_X_SIZE,_Y_SIZE,-0.3,_X_SIZE,_Y_SIZE,0.3, 3000, "axis z");
  plot_field(fm,"v00.01.01/x_p0.05_axis_field.txt",_X_SIZE+0.005, _Y_SIZE, -0.3,_X_SIZE+0.005,_Y_SIZE, 0.3, 2000);
	plot_field(fm,"v00.01.01/y_p0.05_axis_field.txt",_X_SIZE, _Y_SIZE+0.005, -0.3,_X_SIZE,_Y_SIZE+0.005, 0.3, 2000);

  //TESTING MATERIAL
  /*
  std::ofstream file("A1", std::ios_base::trunc);
  file<<"E[V]\tVel[cm/sec]\tDL[cm^2/sec]\tDT[cm^2/sec]\n";
  for (int h = 0;h<200;h++)
  {
     double e = 100+h*(7000-100.0)/(199);
     double vx,vy,vz;
     LAr->ElectronVelocity(e,0,0,0,0,0,vx,vy,vz);
     file<<e<<"\t"<<-vx*1e9<<"\t";
     LAr->ElectronDiffusion(e,0,0,0,0,0,vy,vz);
     file<<-vx*1e9*vy*vy/2<<"\t"<<-vx*1e9*vz*vz/2<<std::endl;
  }
  file.close();
  */
  //FINISHED MATERIAL TESTING

  fm->NotDriftMedium(1);
  fm->NotDriftMedium(2);
  fm->NotDriftMedium(3);

  ViewDrift* viewDr = new ViewDrift();
  viewDr->SetArea(-5.5*_X_SIZE, -5.5*_Y_SIZE, -0.32, 5.5*_X_SIZE, 5.5*_Y_SIZE, 0.32);
  viewDr->SetCanvas(c3);
  //viewDr->EnableDebugging();

  AvalancheMC_My *aval = new AvalancheMC_My();
  aval->EnableDiffusion();
  //aval->DisableDiffusion();
  aval->DisableDebugging();
  aval->SetSensor(sensor);
  aval->SetTimeSteps(0.01);

  ViewFEMesh_My* vMsh = new ViewFEMesh_My();
  vMsh->SetCanvas(c3);
  vMsh->SetComponent(fm);
  vMsh->SetArea(-5.5*_X_SIZE, -5.5*_Y_SIZE, -0.15, 5.5*_X_SIZE, 5.5*_Y_SIZE, 0.15);
  vMsh->SetPlane(0,-1,0,0,_Y_SIZE,0);
  vMsh->SetFillMesh(true);
  vMsh->SetColor(1,kGray);
  vMsh->SetColor(2,kYellow);
  vMsh->SetColor(3,kYellow);
  vMsh->EnableAxes();
  vMsh->EnableDebugging();
  vMsh->SetViewDrift(viewDr);
	vMsh->SetDrawViewRegion(false);
  //aval->EnablePlotting(viewDr);

  int num_pass=0;
  int num_norm=0;
  int num_reached_anode=0;
  int tmp=0, stat=0;
  double x1,y1,z1,t1,x2,y2,z2,t2;
  int NUM=10000;
  NUM = 0;
  int interval = NUM/100;
  int counter=1;

  TCanvas* c4 = new TCanvas("electron distr before_","electron distr before_");
  TH2D *hist = new TH2D("electron distribution below", "electron distribution below",100, -_X_SIZE, _X_SIZE,100, -_Y_SIZE, _Y_SIZE);

  TH2D *hist2 = new TH2D("electron distr after", "electron distribution after",100,-_X_SIZE, _X_SIZE, 100, -_Y_SIZE, _Y_SIZE);

//#endif
  for (int gg=0;gg<NUM;gg++)
  {
     //std::cout<<"#"<<gg<<std::endl;
     double x=-_X_SIZE+RndmUniform()*(_X_SIZE*2),y=-_Y_SIZE+_Y_SIZE*2*RndmUniform();
     hist->Fill(x,y);
     aval->DriftElectron(x, y, -0.129, 0);
     tmp=aval->GetNumberOfElectronEndpoints();
     if (tmp>0) {
        aval->GetElectronEndpoint(tmp-1, x1, y1, z1, t1, x2, y2, z2, t2, stat);
        if (StatusLeftDriftMedium==stat)
           num_norm++;
        if (StatusLeftDriftArea==stat) {
           num_norm++;
           num_pass++;
           hist2->Fill(x2,y2);
        }
        if (std::abs(z2-0.15)<0.01)
           ++num_reached_anode;
     } else {
	 	std::cout<<"e#"<<gg<<": no endpoint"<<std::endl;
	 }
     if (interval*counter<gg) {
         std::cout<<"done "<<(int)(100.0*gg/NUM)<<"%\n";
         ++counter;
     }
  }
  std::cout<<"TRANSPARENCY: "<<(double)num_pass/num_norm<<std::endl;
  std::cout<<"\t"<<num_pass<<"/"<<num_norm<<std::endl;
  std::cout<<"\tNum reached anone: "<<num_reached_anode<<std::endl;
  hist->SetMarkerSize(15);
  hist->Draw();
  c4->Update();
  TCanvas* c5 = new TCanvas("drift lines after","drift lines after");
  hist2->SetMarkerSize(15);
  hist2->Draw();
  c5->Update();

  AvalancheMC_My *aval2 = new AvalancheMC_My();
  aval2->EnableDiffusion();
  //aval2->DisableDiffusion();
  aval2->DisableDebugging();
  aval2->SetSensor(sensor);
  aval2->SetTimeSteps(0.01);
  aval2->EnablePlotting(viewDr);
	double ex,ey,ez,ee_above, ee_below;
  Medium* m;
  int status;
	fm->ElectricField(_X_SIZE,_Y_SIZE, 0.149,ex,ey,ez,m,status);
	ee_above=std::sqrt(ex*ex+ey*ey+ez*ez);
	fm->ElectricField(_X_SIZE,_Y_SIZE, -0.149,ex,ey,ez,m,status);
	ee_below=std::sqrt(ex*ex+ey*ey+ez*ez);
	double E_ratio = ee_above/ee_below;
	std::cout<<"Drift field: "<<ee_below<<"V/cm"<<std::endl;
	std::cout<<"Induction field: "<<ee_above<<"V/cm"<<std::endl;
	std::cout<<"Field ratio: "<<E_ratio<<std::endl;
  for (int gg=55;gg<55;gg++) {
     aval2->DriftElectron((-5.2*_X_SIZE+gg*(_X_SIZE*5.2*2)/54),_Y_SIZE, -0.149, 0);
  }
	for (int gg=55;gg<(int) (55.0 / E_ratio);gg++) {
     aval2->DriftElectron((-5.2*_X_SIZE+gg*(_X_SIZE*5.2*2)/54 * E_ratio),_Y_SIZE, 0.149, 0);
  }
  std::cout<<"vMsh->Plot(): ";
  std::cout<<vMsh->Plot()<<std::endl;
  std::cout<<"app.Run"<<std::endl;
  app.Run(kTRUE);
  return 0;
}
