// Run .L init.cpp and .L plot_Npe_spectrum.cpp
// before executing this script

void print_Npe_vs_V (void)
{
  std::string this_str = "print_Npe_vs_V";

  //std::vector<double> Vs = {6180, 5993, 5728, 5297, 4856, 4413, 3972, 3531, 3090, 2648, 2206, 1765};
  std::vector<double> Vs = {200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, \
  1600, 1765, 1900, 2100, 2206, 2400, 2648, 2800, 3090, 3150, 3300, 3531, 3700, 3972, \
  4200, 4413, 4600, 4856, 5050, 5297, 5500, 5728, 5993, 6180, 6400, 6600, 6800, 7000, \
  7300, 7600, 7900, 8100, 8400, 8700, 9000};
  //std::vector<double> Vs = {11.0, 12.0, 13.0, 14.0, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5};
  std::sort(Vs.begin(), Vs.end());
  std::string folder = "../results/v10_old_setup/transfer_XS_more_fields/";
  double NBrS_yield_factor = 10;
  bool is_23_SiPMS = true;
  std::vector<std::string> inputs;
  for (std::size_t i = 0, i_end_ = Vs.size(); i!=i_end_; ++i) {
    inputs.push_back(folder + dbl_to_str(Vs[i], 0) +"V/recorded.dat");
  }
  std::string output_fname = folder + "Npe_vs_V.txt";

  std::ofstream str;
  str.open(output_fname, std::ios_base::trunc);
  if (!str.is_open()) {
    std::cerr<<this_str<<":Error: Could not open output file \""<<output_fname<<"\"!"<<std::endl;
    return;
  }
  if (is_23_SiPMS)
    str<<"//Voltage[V]	Npe_23SiPM	Npe_1PMT"<<std::endl;
  else
    str<<"//Voltage[V]	Npe_25SiPM	Npe_1PMT"<<std::endl;
  str<<"//Npe is from NBrS and per one electron"<<std::endl;
  for (std::size_t i = 0, i_end_ = Vs.size(); i!=i_end_; ++i) {
    str<<dbl_to_str(Vs[i], 1)<<"\t"
       <<plot_Npe_spectrum(inputs[i], false, NBrS_yield_factor, !is_23_SiPMS)<<"\t"
       <<plot_Npe_spectrum(inputs[i], true, NBrS_yield_factor, !is_23_SiPMS)<<std::endl;
  }
  str.close();
}
