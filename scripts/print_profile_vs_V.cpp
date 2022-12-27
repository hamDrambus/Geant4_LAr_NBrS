// Run .L init.cpp and .L plot_Npe_profile.cpp
// before executing this script

void print_profile_vs_V (void)
{
  std::string this_str = "print_profile_vs_V";
  std::vector<double> Vs;
  std::string simulation_folder;
  std::string experiment_results;
  std::string output_fname;
  double NBrS_yield_factor;
  auto fname_builder = +[](std::string folder, double voltage) -> std::string { return folder + dbl_to_str(voltage, 0); };
  auto get_charge = +[](double voltage) -> double { return 1.0; };

  int check = 0;
  if (false) {
    Vs = {1765, 2206, 2648, 3090, 3531, 3972, 4413, 4856, 5297, 5738, 5993, 6180};
    simulation_folder = "../results/v10_old_setup/transfer_XS_with_diffusion/";
    experiment_results = "../../../Post_processing/220203/results_v5/SiPM_Npes.txt";
    output_fname = "../../../Post_processing/220203/results_v5/Source_XY.txt";
    NBrS_yield_factor = 10;
    get_charge = +[](double voltage) -> double { return 1.27e6; };
    fname_builder = +[](std::string folder, double voltage) -> std::string { return folder+dbl_to_str(voltage, 0) +"V/recorded.dat"; };
    ++check;
  }
  if (false) {
    Vs = {1762, 2644, 3525, 3966, 4847, 5288, 5728, 6169};
    simulation_folder = "../results/v10_old_setup/transfer_XS_with_diffusion/";
    experiment_results = "../../../Post_processing/220804/results_v5/SiPM_Npes_Q1.00.txt";
    output_fname = "../../../Post_processing/220804/results_v5/Source_XY_Q1.00.txt";
    NBrS_yield_factor = 10;
    get_charge = +[](double voltage) -> double { return 1.27e6; };
    fname_builder = +[](std::string folder, double voltage) -> std::string { return folder+dbl_to_str(voltage, 0) +"V/recorded.dat"; };
    ++check;
  }
  if (false) {
    Vs = {1762, 2644, 3525, 3966, 4847, 5288, 5728, 6169};
    simulation_folder = "../results/v10_old_setup/transfer_XS_with_diffusion/";
    experiment_results = "../../../Post_processing/220804/results_v5/SiPM_Npes_Q0.20.txt";
    output_fname = "../../../Post_processing/220804/results_v5/Source_XY_Q0.20.txt";
    NBrS_yield_factor = 10;
    get_charge = +[](double voltage) -> double { return 1.575e5; };
    fname_builder = +[](std::string folder, double voltage) -> std::string { return folder+dbl_to_str(voltage, 0) +"V/recorded.dat"; };
    ++check;
  }
  if (false) {
    Vs = {2644, 3966, 5288, 6169};
    simulation_folder = "../results/v10_old_setup/transfer_XS_with_diffusion/";
    experiment_results = "../../../Post_processing/220804/results_v5/SiPM_Npes_Q0.04.txt";
    output_fname = "../../../Post_processing/220804/results_v5/Source_XY_Q0.04.txt";
    NBrS_yield_factor = 10;
    get_charge = +[](double voltage) -> double { return 5.39e4; };
    fname_builder = +[](std::string folder, double voltage) -> std::string { return folder+dbl_to_str(voltage, 0) +"V/recorded.dat"; };
    ++check;
  }
  if (true) {
    Vs = {11.0, 12.0, 13.0, 14.0, 15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5};
    simulation_folder = "../results/v11_setup_y2022/transfer_XS_with_diff/";
    experiment_results = "../../../Post_processing/221020/results_v1/SiPM_Npes.txt";
    output_fname = "../../../Post_processing/221020/results_v1/Source_XY.txt";
    NBrS_yield_factor = 100;
    get_charge = +[](double voltage) -> double {
      double Edrift = voltage * 3.0 * 80/800/4.8; // kV/cm
      double Qnominal = 5.5e6/23.6/(1 + 23.1/Edrift); //5.5 MeV alpha-particle. W=23.6+-0.3 eV, Krecombination = 23.1+-1.1 kV/cm
      return Qnominal * 0.93; // Average between approximate adjustment for 150 us electron lifetime and nominal charge.
    };
    fname_builder = +[](std::string folder, double voltage) -> std::string { return folder+dbl_to_str(voltage, 1) +"kV/recorded.dat"; };
    ++check;
  }

  if (check != 1) {
    std::cerr<<"Invalid parameter configureation! Aborting."<<std::endl;
    return;
  }
  std::sort(Vs.begin(), Vs.end());
  std::vector<std::string> inputs;
  for (std::size_t i = 0, i_end_ = Vs.size(); i!=i_end_; ++i) {
    inputs.push_back(fname_builder(simulation_folder, Vs[i]));
  }
  std::ofstream str;
  str.open(output_fname, std::ios_base::trunc);
  if (!str.is_open()) {
    std::cerr<<this_str<<":Error: Could not open output file \""<<output_fname<<"\"! Aborting."<<std::endl;
    return;
  }
  str<<"//Voltage[kV]\tSource x [mm]\tSource y [mm]"<<std::endl;
  str<<"//Source coordinates are from SiPM weighted average and are adjusted for systematic shift."<<std::endl;
  Npe_profile_result result;
  for (std::size_t i = 0, i_end_ = Vs.size(); i!=i_end_; ++i) {
    result = plot_Npe_profile(inputs[i], experiment_results, get_charge(Vs[i]), NBrS_yield_factor, false);
    str<<dbl_to_str(Vs[i], 2)<<"\t"
       <<result.source_xy_exp.first<<"\t"<<result.source_xy_exp.second<<std::endl;
  }
  str.close();
}
