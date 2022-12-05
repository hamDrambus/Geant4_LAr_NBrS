// Run .L init.cpp and .L plot_Npe_spectrum.cpp
// before executing this script

void print_Npe_vs_V (void)
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
  std::string folder = "../results/v15_GEM1/transfer_XS_with_diff/";
  double NBrS_yield_factor = 100;
  std::vector<std::string> inputs;
  for (std::size_t i = 0, i_end_ = Vs.size(); i!=i_end_; ++i) {
    inputs.push_back(folder+dbl_to_str(Vs[i], 0) +"V/generated.dat");
  }
  std::string output_fname = folder + "Npe_vs_V.txt";

  std::ofstream str;
  str.open(output_fname, std::ios_base::trunc);
  if (!str.is_open()) {
    std::cerr<<this_str<<":Error: Could not open output file \""<<output_fname<<"\"!"<<std::endl;
    return;
  }
  str<<"//Vgem1[V]	Npe_24SiPM	Npe_1PMT"<<std::endl;
  str<<"//Npe is from NBrS and per one electron"<<std::endl;
  QE_result PMT_QE, SiPM_QE;
  for (std::size_t i = 0, i_end_ = Vs.size(); i!=i_end_; ++i) {
    PMT_QE = plot_Npe_spectrum(inputs[i], true, NBrS_yield_factor, {}, false);
    SiPM_QE = plot_Npe_spectrum(inputs[i], false, NBrS_yield_factor, {44}, false);
    str<<dbl_to_str(Vs[i], 2)<<"\t"
       <<SiPM_QE.Npe_per_e<<"\t"<<PMT_QE.Npe_per_e<<std::endl;
  }
  str.close();
}
