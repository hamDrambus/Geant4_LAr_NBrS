// Run .L init.cpp before executing this script

void plot_Npe_spectrum(void)
{
  gStyle->SetCanvasDefH(800);
	gStyle->SetCanvasDefW(1000);
  gErrorIgnoreLevel = 1001;
  gStyle->SetStatY(0.9);
	gStyle->SetStatX(0.9);
  gStyle->SetStatH(0.3);
  gStyle->SetGridStyle(3);
	gStyle->SetGridColor(14);
	gStyle->SetGridWidth(1);
	gStyle->SetOptStat("e");

  double NBrS_yield_factor = 10;
  bool isPMT = 1;
  //std::string input_file = "../results/v8/6180V/generated.dat";
  std::string input_file = "../results/v9.0_transfer_with_diffusion/1765V/recorded.dat";

  std::string plot_name = isPMT ? "4PMTs photon hits spectrum" : "SiPM matrix photon hits spectrum";
  //std::string plot_name = "Generated photon spectrum";
  std::string QE_filename = isPMT ? "../data/quantum_efficiency/PMT_R6041_506MOD.dat" : "../data/quantum_efficiency/SiPM_s13360-6050pe_46V.dat";
  bool linear_x = true, linear_y = true, linear_z = true;
  std::string X_axis_title = "#lambda [nm]";
  std::string Y_axis_title = "Counts [PE]";
  std::pair<double, double> X_axis_range(0, 1000.0);
  std::pair<double, double> Y_axis_range(0.0, 1e3);
  std::pair<double, double> Y2_axis_range(0.0, 100.0);
  std::pair<double, double> X_axis_zoom(0, 0);
  std::pair<double, double> Y_axis_zoom(0, 0);
  int X_n_bins = 100, Y_n_bins = 1000;

  PlotInfo plot_info;
  plot_info.selection = isPMT ? SelectPMTs : SelectSiPMs;
  plot_info.histogram = NULL;
  plot_info.plot_par_x = PlotParameter::gammaSpectrum;
  plot_info.plot_par_y = PlotParameter::None;

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
  double title_offset_00 = frame_00->GetYaxis()->GetTitleOffset() + 1.3;
  frame_00->GetYaxis()->SetTitleOffset(title_offset_00);
	if (X_axis_zoom.first!=X_axis_zoom.second)
		frame_00->GetXaxis()->SetRangeUser(X_axis_zoom.first, X_axis_zoom.second);
  if (Y_axis_zoom.first!=Y_axis_zoom.second)
		frame_00->GetYaxis()->SetRangeUser(Y_axis_zoom.first, Y_axis_zoom.second);
	frame_00->Draw();
  frame_00->SetStats(false);

  TH1D* hist_00 = new TH1D("hist_00", plot_name.c_str(), X_n_bins, X_axis_range.first, X_axis_range.second);
  plot_info.histogram = hist_00;

  std::ifstream str;
  str.open(input_file);
  if (!str.is_open()) {
    std::cerr<<"Could not open file \""<<input_file<<"\"!"<<std::endl;
  }
  FileToHist(plot_info, str);
  hist_00->Draw("sames");
  c_00->Update();

  ExperimentalXY* SiPM_QE = LoadQE(QE_filename, true);
  // scale graph to the pad coordinates
  float rightmax = Y2_axis_range.second - Y2_axis_range.first;
  float scale = gPad->GetUymax()/rightmax;
  TGraph *graph_00 = GraphFromData(SiPM_QE, 1.0, 100*scale); //100 is for percents
  if (nullptr != graph_00) {
    graph_00->SetLineColor(kRed);
    graph_00->SetLineWidth(2);
    //graph_00->Scale(100 * scale, "y"); // Not present in ROOT v6.14
    graph_00->Draw("same");
    // draw an axis on the right side
    auto axis = new TGaxis(gPad->GetUxmax(),gPad->GetUymin(),
    gPad->GetUxmax(), gPad->GetUymax(), Y2_axis_range.first, Y2_axis_range.second, 510, "+L");
    double title_offset = axis->GetTitleOffset() + 0.1;
    axis->ImportAxisAttributes(frame_00->GetYaxis());
    axis->SetTitle("SiPM PDE [%]");
    axis->SetTitleOffset(title_offset);
    axis->Draw();
  }
  std::cout<<"Input PE number = "<<CalculateNpe(hist_00)<<std::endl;
  double N_real = CalculateNpe(hist_00, SiPM_QE);
  std::cout<<"Resulting PE number = "<<N_real<<std::endl;
  double QE = ((double)N_real)/CalculateNpe(hist_00);
  std::cout<<"<QE> = "<<QE<<std::endl;
  std::cout<<"*************************************"<<std::endl;
  std::cout<<"Npe[-1] = "<<GetNpeCh(plot_info, -1)<<std::endl;
  std::cout<<"NBrS real yield was multiplied by: "<<NBrS_yield_factor<<std::endl;
  std::cout<<"Number of electrons = "<<plot_info.N_electrons<<std::endl;
  if (isPMT) {
    std::cout<<"PMT recorded Npe raw = "<<GetNpePMTraw(plot_info)<<std::endl;
    std::cout<<"Average Npe per 1 PMT (grid is accounted for) = "<<GetNpePMTavg(plot_info)<<std::endl;
    std::cout<<"Average real Npe per 1 PMT per 1 e = "<<GetNpePMTavg(plot_info) * QE / plot_info.N_electrons / NBrS_yield_factor <<std::endl;
  } else {
    std::cout<<"SiPMs recorded Npe raw = "<<GetNpeSiPMs(plot_info)<<std::endl;
    std::cout<<"Npe for 23 SiPMs = "<<GetNpeSiPMs23(plot_info)<<std::endl;
    std::cout<<"Average real Npe for 23 SiPMs per 1 e = "<<GetNpeSiPMs23(plot_info) * QE / plot_info.N_electrons / NBrS_yield_factor <<std::endl;
  }

  TPaveStats *ps = (TPaveStats*)c_00->GetPrimitive("stats");
  AddStatValue(ps, "Detected", N_real);
  // the following line is needed to avoid that the automatic redrawing of stats
  hist_00->SetStats(false);

  c_00->Update();
}
