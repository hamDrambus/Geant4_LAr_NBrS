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

struct definitions {
  std::vector<double> Vs;
  std::string folder;
  double NBrS_yield_factor;
  std::vector<int> exclude_channles;
  std::string out_file_header;
  std::vector<std::string> inputs_generated;
  std::vector<std::string> inputs_recorded;
  std::string output_fname;
};

definitions use_v10_old_setup(void) {
  definitions defs;
  defs.folder = "../results/v10_old_setup/exact_XS_with_diffusion_T11_Tr_Boyle15/";
  defs.NBrS_yield_factor = 20;
  // parameters below should be rarely changed
  defs.Vs = {200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, \
  1600, 1762, 1765, 1900, 2100, 2206, 2400, 2644, 2648, 2800, 3090, 3150, 3300, 3525, 3531, 3700, 3972, \
  4200, 4413, 4600, 4847, 4856, 5050, 5288, 5297, 5500, 5728, 5738, 5993, 6169, 6180, 6400, 6600, 6800, \
  7000, 7400, 7800, 8200, 8600, 9200};
  std::sort(defs.Vs.begin(), defs.Vs.end());
  defs.exclude_channles = {43, 44};
  std::string str_SiPMs = "23SiPMs";
  defs.out_file_header = "//V1[V]\tN_electrons\tN_photons_per_e\tNph_"+str_SiPMs+"\tNph_1PMT\tLCE_"+str_SiPMs+"\tLCE_1PMT\tQE_"+str_SiPMs+
        "_recorded\tQE_1PMT_recorded\tQE_"+str_SiPMs+"_generated\tQE_1PMT_generated\tNpe_per_e_"+str_SiPMs+"\tNpe_per_e_1PMT";
  for (std::size_t i = 0, i_end_ = defs.Vs.size(); i!=i_end_; ++i) {
    defs.inputs_generated.push_back(defs.folder+dbl_to_str(defs.Vs[i], 0) +"V/generated.dat");
    defs.inputs_recorded.push_back(defs.folder+dbl_to_str(defs.Vs[i], 0) +"V/recorded.dat");
  }
  defs.output_fname = defs.folder + "Results_vs_V1.txt";
  return defs;
}

definitions use_v11_setup_y2022(void) {
  definitions defs;
  defs.folder = "../results/v11_setup_y2022/transfer_XS_with_diff/";
  defs.NBrS_yield_factor = 1000;
  // parameters below should be rarely changed
  defs.Vs = {6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 21.0, 22.0};
  std::sort(defs.Vs.begin(), defs.Vs.end());
  defs.exclude_channles = {};
  std::string str_SiPMs = "25SiPMs";
  defs.out_file_header = "//V0[kV]\tN_electrons\tN_photons_per_e\tNph_"+str_SiPMs+"\tNph_1PMT\tLCE_"+str_SiPMs+"\tLCE_1PMT\tQE_"+str_SiPMs+
        "_recorded\tQE_1PMT_recorded\tQE_"+str_SiPMs+"_generated\tQE_1PMT_generated\tNpe_per_e_"+str_SiPMs+"\tNpe_per_e_1PMT";
  for (std::size_t i = 0, i_end_ = defs.Vs.size(); i!=i_end_; ++i) {
    defs.inputs_generated.push_back(defs.folder+dbl_to_str(defs.Vs[i], 1) +"kV/generated.dat");
    defs.inputs_recorded.push_back(defs.folder+dbl_to_str(defs.Vs[i], 1) +"kV/recorded.dat");
  }
  defs.output_fname = defs.folder + "Results_vs_V0.txt";
  return defs;
}

definitions use_v15_GEM1(void) {
  definitions defs;
  defs.folder = "../results/v15_GEM1/transfer_XS_with_diff/";
  defs.NBrS_yield_factor = 500;
  // parameters below should be rarely changed
  defs.Vs = {80, 90, 100, 110, 120, 130, 139, 160, 185, 200, 231, 250, 277, 300, 323, \
    350, 369, 390, 416, 450, 462, 480, 508, 525, 554, 580, 600, 620, 646, 670, 693, \
    739, 760, 785, 800, 831, 850, 877, 900, 923, 950, 1000, 1050, 1100};
  std::sort(defs.Vs.begin(), defs.Vs.end());
  defs.exclude_channles = {44};
  std::string str_SiPMs = "24SiPMs";
  defs.out_file_header = "//V1[V]\tN_electrons\tN_photons_per_e\tNph_"+str_SiPMs+"\tNph_1PMT\tLCE_"+str_SiPMs+"\tLCE_1PMT\tQE_"+str_SiPMs+
        "_recorded\tQE_1PMT_recorded\tQE_"+str_SiPMs+"_generated\tQE_1PMT_generated\tNpe_per_e_"+str_SiPMs+"\tNpe_per_e_1PMT";
  for (std::size_t i = 0, i_end_ = defs.Vs.size(); i!=i_end_; ++i) {
    defs.inputs_generated.push_back(defs.folder+dbl_to_str(defs.Vs[i], 0) +"V/generated.dat");
    defs.inputs_recorded.push_back(defs.folder+dbl_to_str(defs.Vs[i], 0) +"V/recorded.dat");
  }
  defs.output_fname = defs.folder + "Results_vs_V1.txt";
  return defs;
}


void print_results_table(void)
{
  std::string this_str = "print_Npe_vs_V";
  definitions defs = use_v10_old_setup();
  //definitions defs = use_v11_setup_y2022();
  //definitions defs = use_v15_GEM1();

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
  str.open(defs.output_fname, std::ios_base::trunc);
  if (!str.is_open()) {
    std::cerr<<this_str<<":Error: Could not open output file \""<<defs.output_fname<<"\"!"<<std::endl;
    return;
  }
  str<<defs.out_file_header<<std::endl;
  for (std::size_t i = 0, i_end_ = defs.Vs.size(); i!=i_end_; ++i) {
    plot_info_generated.histogram->Reset("ICESM");
    plot_info_SiPMs.histogram->Reset("ICESM");
    plot_info_PMTs.histogram->Reset("ICESM");
    std::ifstream inp;
    inp.open(defs.inputs_generated[i]);
    if (!inp.is_open()) {
      std::cerr<<"Could not open file \""<<defs.inputs_generated[i]<<"\"! Skipping."<<std::endl;
      continue;
    }
    FileToHist(plot_info_generated, inp);
    inp.close();
    inp.open(defs.inputs_recorded[i]);
    if (!inp.is_open()) {
      std::cerr<<"Could not open file \""<<defs.inputs_recorded[i]<<"\"! Skipping."<<std::endl;
      continue;
    }
    FileToHist(plot_info_SiPMs, inp);
    FileToHist(plot_info_PMTs, inp);
    inp.close();

    double Nphotons_per_e = GetNpeCh(plot_info_generated, -1) / defs.NBrS_yield_factor / plot_info_generated.N_electrons;
    double SiPMs_ph_rec = GetNpeSiPMsSome(plot_info_SiPMs, defs.exclude_channles) / defs.NBrS_yield_factor; // Total number of photoelectrons recorded in geant4 adjusted for NBrS yield factor.
    double PMTs_ph_rec = GetNpePMTavg(plot_info_PMTs) / defs.NBrS_yield_factor;
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
    str<<dbl_to_str(defs.Vs[i], 2)<<"\t"
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
