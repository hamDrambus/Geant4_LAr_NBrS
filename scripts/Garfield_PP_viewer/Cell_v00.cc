#include <iostream>
#include <sstream>
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
#include "ViewFEMesh.hh"
#include "ComponentConstant.hh"
#include "Sensor.hh"
#include "TrackHeed.hh"
#include "Plotting.hh"
#include "AvalancheMC_My.hh"
#include "Mat_Vd_const.hh"
#include "Random.hh"
#include "GarfieldConstants.hh"

using namespace Garfield;

// Garfield sizes are in cm
//Thin GEM (standard)
//#define _X_SIZE 0.0035
//#define _Y_SIZE 0.00606217783
//#define _Z_CATHODE -0.103
//#define _Z_ANODE 0.103
#define _TIME_STEP 0.2
//28% CERN THGEM
#define _X_SIZE 0.0225
#define _Y_SIZE 0.038971143170
#define _Z_CATHODE -0.503
#define _Z_ANODE 0.503

//qewr
//2D and 3D can be viewed with the same code due to selected geometries!
#define MESH_ std::string("../../tests/19_field_THGEM1_220113/v00.01_THGEM1/")
#define RESULT_ std::string("../../tests/19_field_THGEM1_220113/Elmer_v00.01/case_2650v.result")
#define RESULT_FOLDER std::string("../../tests/19_field_THGEM1_220113/Elmer_v00.01/")

void ensure_file(std::string fname); //makes sure file can be created later on
void ensure_folder(std::string folder);
std::string dbl_to_str (double val, int precision);

std::string dbl_to_str (double val, int precision)
{
	std::stringstream ss;
	ss<<std::fixed<<std::setprecision(precision)<<val;
	return ss.str();
}

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

void plot_field(ComponentFieldMap* fm, std::string filename, double x0, double y0, double z0, double x1, double y1, double z1, std::size_t Num, std::string name="", double L_fine=0, std::size_t Num_fine=0)
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
  if ((L_fine<=0)||(Num_fine>=Num)||(L_fine>=s)||(RN<LN)||(Num_fine<=0)||(LN>Num)||(LN<0)) {
    Num_fine=0;
    L_fine=0;
    LN=Num;
    RN=-Num;
  }
  std::ofstream flll(filename,std::ios_base::out);
	//flll<<"//z\tx\ty\tL\tEabs\tEx\tEy\tEz"
  for (int hh=0;hh<=Num;hh++) { //hh<=Num in order to cover [x0;x1], not [x0;x1)
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
     if ((hh<LN)||(hh>RN)) {
        ds=(L-L_fine)/(Num-LN-RN-1);
        dphi+=(ds/L)*(ex*(x1-x0)+ey*(y1-y0)+ez*(z1-z0));
        s+=ds;
     } else {
        ds=L_fine/(RN-LN+1);
        dphi+=ex*ds*(x1-x0)/L+ey*ds*(y1-y0)/L+ez*ds*(z1-z0)/L;
        s+=ds;
     }
  }
  if (name.size())
		std::cout<<"Int[E(r)*dr] on "<<name<<" = "<<dphi<<std::endl;
  flll.close();
}

class PointE3D {
public:
	PointE3D() : E(0), Ex(0), Ey(0), Ez(0), x(DBL_MAX), y(DBL_MAX), z(DBL_MAX)
	{}
	double E;
	double Ex;
	double Ey;
	double Ez;
	double x;
	double y;
	double z;
};

PointE3D find_max_E_along_line(ComponentFieldMap* fm, double x0, double y0, double z0, double x1, double y1, double z1, std::size_t Num_pts)
{
	PointE3D out;
  double L = sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0));
  double x=x0,y=y0,z=z0;
  double ex,ey,ez,ee;
  Medium* m;
  int status;
  for (int hh=0; hh<= Num_pts; ++hh) { //hh<=Num in order to cover [x0;x1], not [x0;x1)
		double s = ((double)hh)/Num_pts;
		x=(x1-x0)*s+x0;
		y=(y1-y0)*s+y0;
		z=(z1-z0)*s+z0;
		fm->ElectricField(x,y,z,ex,ey,ez,m,status);
		ee=std::sqrt(ex*ex+ey*ey+ez*ez);
		if (ee > out.E) {
			out.x = x;
			out.y = y;
			out.z = z;
			out.E = ee;
			out.Ex = ex;
			out.Ey = ey;
			out.Ez = ez;
		}
	}
  return out;
}

double E_at_hole_center(ComponentFieldMap* fm)
{
  double ex,ey,ez,ee;
  Medium* m;
  int status;
	fm->ElectricField(_X_SIZE* 0.999999, _Y_SIZE* 0.999999, 0, ex, ey, ez, m, status);
	ee=std::sqrt(ex*ex+ey*ey+ez*ez);
  return ee;
}

ComponentElmer* LoadFieldMap(std::string mesh_folder, std::string elmer_result) {
	ComponentElmer* fm = new ComponentElmer();
  bool is_loaded = fm->Initialise(mesh_folder+"mesh.header", mesh_folder+"mesh.elements", mesh_folder+ "mesh.nodes", mesh_folder+ "diels.dat", elmer_result, "mm");
	if (is_loaded)
  	std::cout<<"Loaded successfully"<<std::endl;
  else {
		std::cout<<"Load failure, quiting"<<std::endl;
		return nullptr;
  }
  fm->EnableMirrorPeriodicityX();
  fm->EnableMirrorPeriodicityY();
  for (int i=0; i < fm->GetNumberOfMaterials(); i++)
		fm->NotDriftMedium(i);
	fm->DriftMedium(0);
  Mat_Vd_const* LAr = new Mat_Vd_const(5e4);
  //LAr->SetComposition("ar", 100.);
  LAr->SetTemperature(89);
  LAr->EnableDrift();
	LAr->SetPressure(760.);

  MediumMagboltz* met = new MediumMagboltz();
  met->SetComposition("co2", 100.);
  met->SetTemperature(293.15);
  met->SetPressure(760.);
  met->EnableDrift(false);
  MediumMagboltz* diel = new MediumMagboltz();
  diel->SetComposition("co2", 20.,"h2",80.);
  diel->SetTemperature(293.15);
  diel->SetPressure(760.);
  diel->EnableDrift(false);

  fm->SetMedium(0,LAr);
  fm->SetMedium(1,diel);  //FR4
  fm->SetMedium(2,met);
  fm->SetMedium(3,met);
	return fm;
}

void DriftLinesXZ(ComponentElmer* fm, std::size_t line_n, bool use_above = true) {
	if (0 != line_n)
		--line_n;
	else
		return;
	Sensor* sensor = new Sensor();
  sensor->AddComponent(fm);
  sensor->SetArea(-5.5*_X_SIZE,-5.5*_Y_SIZE, _Z_CATHODE, 5.5*_X_SIZE, 5.5*_Y_SIZE, _Z_ANODE);
  TCanvas* c3 = new TCanvas("drift lines xz", "drift lines xz");

	ViewDrift* viewDr = new ViewDrift();
	viewDr->SetArea(-5.5*_X_SIZE, -5.5*_Y_SIZE, _Z_CATHODE*0.2, 5.5*_X_SIZE, 5.5*_Y_SIZE, _Z_ANODE*0.2);
	viewDr->SetCanvas(c3);

	ViewFEMesh* vMsh = new ViewFEMesh();
  vMsh->SetCanvas(c3);
  vMsh->SetComponent(fm);
  vMsh->SetArea(-5.5*_X_SIZE, -5.5*_Y_SIZE, _Z_CATHODE*0.2, 5.5*_X_SIZE, 5.5*_Y_SIZE, _Z_ANODE*0.2);
  vMsh->SetPlane(0,-1,0,0,_Y_SIZE,0);
  vMsh->SetFillMesh(true);
  vMsh->SetColor(1,kGray);
  vMsh->SetColor(2,kYellow);
  vMsh->SetColor(3,kYellow);
  vMsh->EnableAxes();
  vMsh->EnableDebugging(false);
  vMsh->SetViewDrift(viewDr);
	vMsh->SetDrawViewRegion(false);

	AvalancheMC_My *aval = new AvalancheMC_My();
  //aval->EnableDiffusion();
  aval->DisableDiffusion();
	aval->EnableDebugging(false);
  aval->SetSensor(sensor);
  aval->SetTimeSteps(_TIME_STEP);
  aval->EnablePlotting(viewDr);
	double ex,ey,ez_above, ez_below,ee_above, ee_below;
  Medium* m;
  int status;
	fm->ElectricField(_X_SIZE* 0.999999,_Y_SIZE* 0.999999, _Z_CATHODE*0.2,ex,ey,ez_below,m,status);
	ee_below=std::sqrt(ex*ex+ey*ey+ez_below*ez_below);
	fm->ElectricField(_X_SIZE* 0.999999,_Y_SIZE* 0.999999, _Z_ANODE*0.2,ex,ey,ez_above,m,status);
	ee_above=std::sqrt(ex*ex+ey*ey+ez_above*ez_above);
	double E_ratio = ee_above/ee_below;
	std::cout<<"Drift field: "<<ee_below<<"V/cm"<<std::endl;
	std::cout<<"Induction field: "<<ee_above<<"V/cm"<<std::endl;
	std::cout<<"Field ratio: "<<E_ratio<<std::endl;
  for (std::size_t gg=0; gg <= line_n; ++gg) {
		if (ez_below < 0)
    	aval->DriftElectron((-5.0*_X_SIZE+gg*(_X_SIZE*5.0*2)/line_n),_Y_SIZE, _Z_CATHODE*0.2, 0);
		else
			aval->DriftPositron((-5.0*_X_SIZE+gg*(_X_SIZE*5.0*2)/line_n),_Y_SIZE, _Z_CATHODE*0.2, 0);
  }
	std::size_t N = (std::size_t) (line_n * E_ratio);
	for (std::size_t gg=0; gg <= N && use_above; ++gg) {
		if (ez_above > 0)
    	aval->DriftElectron((-5.0*_X_SIZE+gg*(_X_SIZE*5.0*2)/N),_Y_SIZE, _Z_ANODE*0.2, 0);
		else
			aval->DriftPositron((-5.0*_X_SIZE+gg*(_X_SIZE*5.0*2)/N),_Y_SIZE, _Z_ANODE*0.2, 0);
  }
  std::cout<<"DriftLinesX: Calling vMsh->Plot()."<<std::endl;
  std::cout<<vMsh->Plot()<<std::endl;
	std::cout<<"DriftLinesX: Finished vMsh->Plot()."<<std::endl;
}

double THGEM_transparency(ComponentElmer* fm, std::size_t N, bool do_draw = false)
{
	if (0 == N)
		return 0;
	Sensor* sensor = new Sensor();
  sensor->AddComponent(fm);
  sensor->SetArea(-5.5*_X_SIZE, -5.5*_Y_SIZE, _Z_CATHODE, 5.5*_X_SIZE, 5.5*_Y_SIZE, _Z_ANODE);

	AvalancheMC_My *aval = new AvalancheMC_My();
  //aval->EnableDiffusion();
  aval->DisableDiffusion();
  aval->EnableDebugging(false);
  aval->SetSensor(sensor);
  aval->SetTimeSteps(_TIME_STEP);

	TH2D *hist = nullptr, *hist2 = nullptr;
	if (do_draw) {
		hist = new TH2D("electron distribution below", "electron distribution below", 100, -_X_SIZE, _X_SIZE, 100, -_Y_SIZE, _Y_SIZE);
		hist2 = new TH2D("electron distr after", "electron distribution after", 100,-_X_SIZE, _X_SIZE, 100, -_Y_SIZE, _Y_SIZE);
	}

	std::size_t num_pass=0;
  std::size_t num_norm=0;
  std::size_t num_reached_anode=0;
  int tmp=0, stat=0;
  double x1,y1,z1,t1,x2,y2,z2,t2;
  std::size_t interval = N/100;
  std::size_t counter=1;

	std::cout<<"Calculating THGEM transparency..."<<std::endl;
	for (std::size_t gg=0; gg < N; ++gg) {
		//std::cout<<"#"<<gg<<std::endl;
		double x=-_X_SIZE+RndmUniform()*(_X_SIZE*2), y=-_Y_SIZE+_Y_SIZE*2*RndmUniform();
		hist->Fill(x,y);
		aval->DriftElectron(x, y, _Z_CATHODE*0.3, 0);
		tmp=aval->GetNumberOfElectronEndpoints();
		if (tmp > 0) {
			aval->GetElectronEndpoint(tmp-1, x1, y1, z1, t1, x2, y2, z2, t2, stat);
			if (StatusLeftDriftMedium==stat)
				++num_norm;
			if (StatusLeftDriftArea==stat) {
				++num_norm;
				++num_pass;
				hist2->Fill(x2,y2);
			}
			if (std::abs(z2-_Z_ANODE) < 0.0001)
				++num_reached_anode;
		} else {
			std::cout<<"e#"<<gg<<": no endpoint"<<std::endl;
		}
		if (interval*counter < gg) {
			std::cout<<"done "<<(int)(100.0*gg/N)<<"%\n";
			++counter;
		}
  }
	double result = (double)num_pass/num_norm;
  std::cout<<"TRANSPARENCY: "<<result<<std::endl;
  std::cout<<"\t"<<num_pass<<"/"<<num_norm<<std::endl;
  std::cout<<"\tNumber of e which reached anode: "<<num_reached_anode<<std::endl;
	if (do_draw) {
		TCanvas* c4 = new TCanvas("electron distr before_", "electron distr before_");
	  hist->SetMarkerSize(0.2);
	  hist->Draw();
	  c4->Update();
	  TCanvas* c5 = new TCanvas("electron distr after", "electron distr after");
	  hist2->SetMarkerSize(0.2);
	  hist2->Draw();
	  c5->Update();
	}
	return result;
}

int main(int argc, char * argv[]) {
	// Last modified on 2022.11.22
	// Plots drift lines and calculates THGEM transparency (no diffusion)
  TApplication app("app", &argc, argv);
	gStyle->SetCanvasDefH(800);
	gStyle->SetCanvasDefW(800);
  plottingEngine.SetDefaultStyle();

	if (true) {
	  ComponentElmer* fm = LoadFieldMap(MESH_, RESULT_);
		if (nullptr == fm)
			return -1;
		//adsf
	  plot_field(fm,"../../tests/19_field_THGEM1_220113/center_axis_field_2650V_16.47kV.txt",_X_SIZE * 0.999999, _Y_SIZE * 0.999999, _Z_CATHODE,_X_SIZE * 0.999999,_Y_SIZE * 0.999999, _Z_ANODE, 3000, "axis z");
		PointE3D pt_max = find_max_E_along_line(fm, _X_SIZE * 0.999999, _Y_SIZE * 0.999999, _Z_CATHODE, _X_SIZE * 0.999999, _Y_SIZE * 0.999999, _Z_ANODE, 3000);
		std::cout<<"Max field is "<<pt_max.E/1000 <<" kV/cm at z = "<< pt_max.z * 10 << " mm"<<std::endl;
		double E_center = E_at_hole_center(fm);
		std::cout<<"E at hole center = "<<E_center/1000 <<" kV/cm"<<std::endl;
		//THGEM_transparency(fm, 1000, true);
	  AvalancheMC_My *aval = new AvalancheMC_My();
		DriftLinesXZ(fm, 30, true);
	}
	if (false) {
		std::vector<double> Vs = {0, 178, 298, 373, 447, 521, 595, 746, 930, 1130, 1257, 1506, 1757, 2009, 2260, 2350};
		std::string plot_file = "../../tests/19_field_THGEM1_220113/hole_field_vs_Vt_1.0atm.txt";
		ensure_file(plot_file);
		std::ofstream str(plot_file, std::ios_base::trunc);
		str<<"Vdivider[V]\tEmax[kV/cm]\tE at hole center[kV/cm]"<<std::endl;
		for (std::size_t i = 0; i!=Vs.size(); ++i) {
			std::string fieldmap = RESULT_FOLDER + "case_"+dbl_to_str(Vs[i], 0) +"v.result";
			ComponentElmer* fm = LoadFieldMap(MESH_, fieldmap);
			if (nullptr == fm) {
				std::cout<<"Elmer file: \""<<fieldmap<<"\""<<std::endl;
				continue;
			}
			PointE3D pt_max = find_max_E_along_line(fm, _X_SIZE * 0.999999, _Y_SIZE * 0.999999, _Z_CATHODE*0.2, _X_SIZE * 0.999999, _Y_SIZE * 0.999999, _Z_ANODE*0.2, 2000);
			double E_center = E_at_hole_center(fm);
			str<<Vs[i]<<"\t"<<pt_max.E/1000.0<<"\t"<<E_center/1000.0<<std::endl;
			delete fm;
		}
		str.close();
	}

  std::cout<<"app.Run"<<std::endl;
  app.Run(kTRUE);
  return 0;
}
