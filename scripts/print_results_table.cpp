// Run .L init.cpp
// before executing this script

// Prints all results obtained in simulation vs voltage (or drift step):
// 1. Total generated N electrons
// 2. Total generated N photons (divided by NBrS factor) per single electron
// 3. Total N photons that reached SiPM-matrix in geant4 (divided by NBrS factor, some channels my be turned off)
// 4. Total N photons that reached 1 PMT in geant4 (average, divided by NBrS factor)
// 5. Average light collection efficiency for SiPM-matrix (some channels my be turned off) (= [3]/[2])
// 6. Average light collection efficiency (LCE) for single PMT (= [4]/[2])
// 7. <QE> for SiPM-matrix using spectrum of photons that reached PMTs
// 8. <QE> for PMTs using spectrum of photons that reached PMTs
// 9. <QE> for SiPM-matrix using spectrum of generated photons
// 10. <QE> for PMTs using spectrum of generated photons
// 11. N photoelectrons per electron for SiPM-matrix (some channels my be turned off)
// 12. N photoelectrons per electron for 1 PMT (= Total generated N photons [2] * LCE * <QE>)

void print_results_table(void)
{
  std::string this_str = "print_Npe_vs_V";

  //std::vector<double> Vs = {6180, 5993, 5728, 5297, 4856, 4413, 3972, 3531, 3090, 2648, 2206, 1765};
  //std::vector<double> Vs = {200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, \
  1600, 1765, 1900, 2100, 2206, 2400, 2648, 2800, 3090, 3150, 3300, 3531, 3700, 3972, \
  4200, 4413, 4600, 4856, 5050, 5297, 5500, 5728, 5993, 6180, 6400, 6600, 6800, 7000, \
  7300, 7600, 7900, 8100, 8400, 8700, 9000};
  std::vector<double> Vs = {80, 90, 100, 110, 120, 130, 139, 160, 185, 200, 231, 250, 277, 300, 323, \
    350, 369, 390, 416, 450, 462, 480, 508, 525, 554, 580, 600, 620, 646, 670, 693, \
    620, 739, 760, 785, 800, 831, 850, 877, 900, 923, 950, 1000, 1050, 1100};
  //std::vector<double> Vs = {11.0, 12.0, 13.0, 14.0, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5};
  //Drift steps:
  //std::vector<double> Vs = {0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 1.0, 2.0, 4.0, 7.0, 10.0};
  std::sort(Vs.begin(), Vs.end());
  std::string folder = "../results/v15_GEM1/elastic_XS_with_diff/";
  double NBrS_yield_factor = 100;
  std::vector<int> exclude_channles = {44};
  std::vector<std::string> inputs_gen, inputs_rec;
  for (std::size_t i = 0, i_end_ = Vs.size(); i!=i_end_; ++i) {
    inputs_gen.push_back(folder+dbl_to_str(Vs[i], 0) +"V/generated.dat");
    inputs_rec.push_back(folder+dbl_to_str(Vs[i], 0) +"V/recorded.dat");
  }
  std::string output_fname = folder + "Results_vs_V.txt";
  ExperimentalXY* PMTs_QE_data = LoadQE("../data/quantum_efficiency/PMT_R6041_506MOD.dat", true);
  ExperimentalXY* SiPMs_QE_data = LoadQE("../data/quantum_efficiency/SiPM_s13360-6050pe_46V.dat", true);
  PlotInfo plot_info_generated, plot_info_PMTs, plot_info_SiPMs;
  plot_info_generated.selection = SelectAll;
  plot_info_PMTs.selection = SelectPMTs;
  plot_info_SiPMs.selection = SelectSiPMs;
  plot_info_generated.plot_par_x = plot_info_PMTs.plot_par_x = plot_info_SiPMs.plot_par_x = PlotParameter::gammaSpectrum;
  plot_info_generated.plot_par_y = plot_info_PMTs.plot_par_y = plot_info_SiPMs.plot_par_y = PlotParameter::None;
  int X_n_bins = 100;
  std::pair<double, double> X_axis_range(0, 1000.0);
  TH1D* hist_gen = new TH1D("hist_gen", "Generated photon spectrum", X_n_bins, X_axis_range.first, X_axis_range.second);
  plot_info_generated.histogram = hist_gen;
  TH1D* hist_PMTs = new TH1D("hist_PMTs", "PMTs photon hit spectrum", X_n_bins, X_axis_range.first, X_axis_range.second);
  plot_info_PMTs.histogram = hist_PMTs;
  TH1D* hist_SiPMs = new TH1D("hist_SiPMs", "SiPMs photon hit spectrum", X_n_bins, X_axis_range.first, X_axis_range.second);
  plot_info_SiPMs.histogram = hist_SiPMs;

  std::ofstream str;
  str.open(output_fname, std::ios_base::trunc);
  if (!str.is_open()) {
    std::cerr<<this_str<<":Error: Could not open output file \""<<output_fname<<"\"!"<<std::endl;
    return;
  }
  std::string str_SiPMs = "24SiPMs";
  str<<"//Vgem1[V]\tN_electrons\tN_photons_per_e\tNph_"<<str_SiPMs<<"\tNph_1PMT\tLCE_"<<str_SiPMs<<"\tLCE_1PMT\tQE_"<<str_SiPMs<<"_recorded"
      "\tQE_1PMT_recorded\tQE_"<<str_SiPMs<<"_generated\tQE_1PMT_generated\tNpe_per_e_"<<str_SiPMs<<"\tNpe_per_e_1PMT"<<std::endl;
  for (std::size_t i = 0, i_end_ = Vs.size(); i!=i_end_; ++i) {
    plot_info_generated.histogram->Reset("ICESM");
    plot_info_SiPMs.histogram->Reset("ICESM");
    plot_info_PMTs.histogram->Reset("ICESM");
    std::ifstream inp;
    inp.open(inputs_gen[i]);
    if (!inp.is_open()) {
      std::cerr<<"Could not open file \""<<inputs_gen[i]<<"\"! Skipping."<<std::endl;
      continue;
    }
    FileToHist(plot_info_generated, inp);
    inp.close();
    inp.open(inputs_rec[i]);
    if (!inp.is_open()) {
      std::cerr<<"Could not open file \""<<inputs_rec[i]<<"\"! Skipping."<<std::endl;
      continue;
    }
    FileToHist(plot_info_SiPMs, inp);
    inp.clear();
    inp.seekg(0, ios::beg);
    FileToHist(plot_info_PMTs, inp);
    inp.close();

    double Nphotons_per_e = GetNpeCh(plot_info_generated, -1) / NBrS_yield_factor / plot_info_generated.N_electrons;
    double SiPMs_ph_rec = GetNpeSiPMsSome(plot_info_SiPMs, exclude_channles) / NBrS_yield_factor; // Total number of photoelectrons recorded in geant4 adjusted for NBrS yield factor.
    double PMTs_ph_rec = GetNpePMTavg(plot_info_PMTs) / NBrS_yield_factor;
    double LCE_SiPMs = SiPMs_ph_rec / Nphotons_per_e / plot_info_generated.N_electrons; // Light collection efficiency
    double LCE_PMTs = PMTs_ph_rec / Nphotons_per_e / plot_info_generated.N_electrons;
    double SiPMs_QE_rec = CalculateNpe(hist_SiPMs, SiPMs_QE_data)/CalculateNpe(hist_SiPMs); //QE = quantum efficiency
    double PMTs_QE_rec = CalculateNpe(hist_PMTs, PMTs_QE_data)/CalculateNpe(hist_PMTs);
    double SiPMs_QE_gen = CalculateNpe(hist_gen, SiPMs_QE_data)/CalculateNpe(hist_gen);
    double PMTs_QE_gen = CalculateNpe(hist_gen, PMTs_QE_data)/CalculateNpe(hist_gen);
    double SiPMs_Npe_per_e = SiPMs_QE_rec * LCE_SiPMs * Nphotons_per_e;
    double PMTs_Npe_per_e = PMTs_QE_rec * LCE_PMTs * Nphotons_per_e;
    // Diference between SiPMs_QE_gen and SiPMs_QE_rec reflects how different are
    // spectrum of generated photons and spectrum of photons that actually reached SiPM matrix.
    str<<dbl_to_str(Vs[i], 2)<<"\t"
       <<plot_info_generated.N_electrons<<"\t"
       <<Nphotons_per_e<<"\t"
       <<SiPMs_ph_rec<<"\t"
       <<PMTs_ph_rec<<"\t"
       <<LCE_SiPMs<<"\t"
       <<LCE_PMTs<<"\t"
       <<SiPMs_QE_rec<<"\t"
       <<PMTs_QE_rec<<"\t"
       <<SiPMs_QE_gen<<"\t"
       <<PMTs_QE_gen<<"\t"
       <<SiPMs_Npe_per_e<<"\t"
       <<PMTs_Npe_per_e<<std::endl;
  }
  str.close();
}
