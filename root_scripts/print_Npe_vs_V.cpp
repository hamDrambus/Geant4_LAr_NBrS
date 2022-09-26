// Run .L init.cpp and .L plot_Npe_spectrum.cpp
// before executing this script

void print_Npe_vs_V (void)
{
  std::string this_str = "print_Npe_vs_V";

  std::vector<int> Vs = {6180, 5993, 5728, 5297, 4856, 4413, 3972, 3531, 3090, 2648, 2206, 1765};
  std::sort(Vs.begin(), Vs.end());
  std::string folder = "../results/v10.1_transfer/";
  double NBrS_yield_factor = 10;
  std::vector<std::string> inputs;
  for (std::size_t i = 0, i_end_ = Vs.size(); i!=i_end_; ++i) {
    inputs.push_back(folder + std::to_string(Vs[i]) +"V/recorded.dat");
  }
  std::string output_fname = folder + "Npe_vs_V.txt";

  std::ofstream str;
  str.open(output_fname, std::ios_base::trunc);
  if (!str.is_open()) {
    std::cerr<<this_str<<":Error: Could not open output file \""<<output_fname<<"\"!"<<std::endl;
    return;
  }
  str<<"//Voltage[V]	Npe_23SiPM	Npe_1PMT"<<std::endl
     <<"//Npe is from NBrS and per one electron"<<std::endl;
  for (std::size_t i = 0, i_end_ = Vs.size(); i!=i_end_; ++i) {
    str<<Vs[i]<<"\t"
       <<plot_Npe_spectrum(inputs[i], false, NBrS_yield_factor)<<"\t"
       <<plot_Npe_spectrum(inputs[i], true, NBrS_yield_factor)<<std::endl;
  }
  str.close();
}
