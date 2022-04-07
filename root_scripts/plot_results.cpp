// Run .L init.cpp before executing this script

void plot_results(void)
{
  gStyle->SetCanvasDefH(800);
	gStyle->SetCanvasDefW(800);
  gErrorIgnoreLevel = 1001;
  gStyle->SetStatY(0.9);
	gStyle->SetStatX(0.9);
  gStyle->SetGridStyle(3);
	gStyle->SetGridColor(14);
	gStyle->SetGridWidth(1);
	gStyle->SetOptStat("e");

  std::string plot_name = "SiPM matrix photon hits x-y";
  bool linear_x = true, linear_y = true, linear_z = true;
  std::string X_axis_title = "X [mm]";
  std::string Y_axis_title = "Y [mm]";
  std::pair<double, double> X_axis_range(-5.0, 5.0);
  std::pair<double, double> Y_axis_range(-5.0, 5.0);
  std::pair<double, double> X_axis_zoom(0, 0);
  std::pair<double, double> Y_axis_zoom(0, 0);
  int X_n_bins = 1000, Y_n_bins = 1000;

  std::string input_file = "../tests/test_00_SiPM_THGEM1_shading/generated.dat";
  //std::string input_file = "../tests/test_00_SiPM_THGEM1_shading/recorded.dat";
  PlotInfo plot_info;
  plot_info.selection = SelectSiPMs;
  plot_info.histogram = NULL;
  plot_info.plot_par_x = PlotParameter::gammaXpos;
  plot_info.plot_par_y = PlotParameter::gammaYpos;

  TCanvas *c_00 = new TCanvas ((std::string("0") + plot_name).c_str(), plot_name.c_str());
	c_00->SetGrid(); c_00->SetTicks(); c_00->ToggleEventStatus(); c_00->ToggleToolBar();
  if (!linear_x)
    c_00->SetLogx();
  if (!linear_y)
    c_00->SetLogy();
  if (!linear_z)
    c_00->SetLogz();
	TH2F* frame_00 = new TH2F("frame_00", plot_name.c_str(), 500, X_axis_range.first, X_axis_range.second, 500, Y_axis_range.first, Y_axis_range.second);
	frame_00->GetXaxis()->SetTitle(X_axis_title.c_str());
  frame_00->GetYaxis()->SetTitle(Y_axis_title.c_str());
	if (X_axis_zoom.first!=X_axis_zoom.second)
		frame_00->GetXaxis()->SetRangeUser(X_axis_zoom.first, X_axis_zoom.second);
  if (Y_axis_zoom.first!=Y_axis_zoom.second)
		frame_00->GetYaxis()->SetRangeUser(Y_axis_zoom.first, Y_axis_zoom.second);
	frame_00->Draw();
  frame_00->SetStats(false);

  TH2D* hist_00 = new TH2D("hist_00", plot_name.c_str(), X_n_bins, X_axis_range.first, X_axis_range.second,
                        Y_n_bins, Y_axis_range.first, Y_axis_range.second);
  plot_info.histogram = hist_00;

  std::ifstream str;
  str.open(input_file);
  if (!str.is_open()) {
    std::cerr<<"Could not open file \""<<input_file<<"\"!"<<std::endl;
  }
  FileToHist(plot_info, str);
  hist_00->Draw("SCATsames");
  c_00->Update();
}
