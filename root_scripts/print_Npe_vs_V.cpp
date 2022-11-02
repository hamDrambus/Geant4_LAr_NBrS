// Run .L init.cpp and .L plot_Npe_spectrum.cpp
// before executing this script

void print_Npe_vs_V (void)
{
  std::string this_str = "print_Npe_vs_V";

  std::vector<double> Vs = {6180, 5993, 5728, 5297, 4856, 4413, 3972, 3531, 3090, 2648, 2206, 1765};
  //std::vector<double> Vs = {11.0, 12.0, 13.0, 14.0, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5};
  std::sort(Vs.begin(), Vs.end());
  std::string folder = "../results/v10_old_setup/transfer_XS_with_diffusion/";
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
