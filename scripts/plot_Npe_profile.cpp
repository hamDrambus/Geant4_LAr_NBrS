// Run .L init.cpp before executing this script
// Created on 2022.12.12

// plot_Npe_profile("../results/v10_old_setup/transfer_XS_with_diffusion/6169V/recorded.dat", "../../../Post_processing/220804/results_v5/SiPM_Npes_Q1.00.txt", 1.27e6, 10);
// plot_Npe_profile("../results/v10_old_setup/transfer_XS_with_diffusion/6169V/recorded.dat", "../../../Post_processing/220804/results_v5/SiPM_Npes_Q0.20.txt", 1.58e5, 10);
// plot_Npe_profile("../results/v10_old_setup/transfer_XS_with_diffusion/6169V/recorded.dat", "../../../Post_processing/220804/results_v5/SiPM_Npes_Q0.04.txt", 5.4e4, 10);
// plot_Npe_profile("../results/v10_old_setup/MT_default/6169V/recorded.dat", "../../../Post_processing/220804/results_v5/SiPM_Npes_Q0.20.txt", 1.58e5, 10);

double SiPM_Npe_stat_error = 1e-3;

struct Npe_profile_result {
  // This block has the same size.
  std::vector<std::vector<int>> equidistant_ch_groups;
  std::vector<int> n_ch_exp; //[i] is != equidistant_ch_groups[i].size() because some channels are turned off in the experiment
  std::vector<double> Radii;
  std::vector<double> Npes_per_e_sim;
  std::vector<double> Npes_per_e_sim_indiv_QE;
  std::vector<double> Npes_per_e_exp;

  std::map<int, double> Npe_per_e_per_ch_exp;
  std::pair<double, double> source_xy_exp; // true-ish coordinate of source obtained from experiment.
};

std::map<int, double> ReadExperimentData(std::string filename, std::string voltage) {
  std::map<int, double> result;
  std::ifstream str;
  str.open(filename);
  if (!str.is_open()) {
    std::cerr<<"Error:ReadExperimentData: could not open file \""<<filename<<"\"!"<<std::endl;
    return result;
  }
  double Voltage;
  try {
    Voltage = stod(voltage);
  } catch (std::exception &e) {
    std::cerr << "Error:ReadExperimentData: could not get voltage from string \"" << voltage << "\"." << std::endl;
    std::cerr << e.what() << std::endl;
    return result;
  }
  std::vector<int> chs;
  std::vector<double> Npes;
  std::size_t line_n = 0;
  std::string line, word;
  bool found = false;
  while (!str.eof() && str.is_open()) {
		std::getline(str, line);
		++line_n;
    try {
  		if (line_n == 1) { // Header line is "V 32 33 .... 59" where numbers are channel numbers
        word = strtoken(line, "\t "); // ignore first word
        while (!line.empty()) {
          word = strtoken(line, "\t ");
          if (!word.empty())
            chs.push_back(std::stoi(word));
        }
        continue;
      }
      word = strtoken(line, "\t ");
      double V = stod(word);
      if ((V - Voltage) < 1e-10*max(V, Voltage)) { // extract only necessary voltage from table in the file
        found = true;
        while (!line.empty()) {
          word = strtoken(line, "\t ");
          if (!word.empty())
            Npes.push_back(std::max(std::stod(word), 0.0));
        }
        break;
      }
		} catch (std::invalid_argument &e) {
    } catch (std::out_of_range &e) {
			std::cerr << "Warning:ReadExperimentData: Invalid data on line " << line_n << std::endl;
			std::cerr << e.what() << std::endl;
			return result;
		} catch (std::exception &e) {
			std::cerr << "Warning:ReadExperimentData: Unforeseen exception on line " << line_n << std::endl;
			std::cerr << e.what() << std::endl;
			return result;
		}
	}
  if (!found) {
    std::cerr<<"Warning:ReadExperimentData: voltage "<<voltage<<" V was not found\n";
    std::cerr<<"\tin the table \""<<filename<<"\"."<<std::endl;
    return result;
  }
  if (chs.size()!=Npes.size()) {
    std::cerr<<"Warning:ReadExperimentData: number of channels and number of table entries are inconsistent!\n";
    std::cerr<<"\tCheck input file \""<<filename<<"\""<<std::endl;
  }
  for (std::size_t i = 0, i_end_ = std::min(chs.size(), Npes.size()); i!=i_end_; ++i)
    result[chs[i]] = Npes[i];
  return result;
}

double get_SiPM_QE(std::ifstream &sim_data, ExperimentalXY* QE, int channel = -1) {
  PlotInfo plot_info;
  plot_info.selection = channel < 0 ? SelectSiPMs : SelectChannel;
  TH1D* hist_00 = new TH1D("hist_QE", "hist_QE", 1000, 0, 1000.0);
  plot_info.histogram = hist_00;
  plot_info.plot_par_x = PlotParameter::gammaSpectrum;
  plot_info.plot_par_y = PlotParameter::None;
  plot_info.selection_info = &channel;
  FileToHist(plot_info, sim_data);
  double Npe_geant4 = CalculateNpe(hist_00);
  double N_real = CalculateNpe(hist_00, QE);

  if (nullptr != plot_info.histogram)
    plot_info.histogram->Delete();
  return Npe_geant4 < 0.1 ? 0 : N_real/Npe_geant4;
}

// Return S2 (source) center based on weighted average.
std::pair<double, double> plot_experiment_XY(std::map<int, double> &Npe_exp, double charge, bool do_draw = true) {
  std::pair<double, double> xy_avg(0, 0);
  double weight = 0;
  TH2D* hist_01 = new TH2D("SiPM experimental S2 profile", "SiPM experimental S2 profile", 2*SiPM_n_rows+1,
						-SiPM_pitch * (SiPM_n_rows + 1) / 2 + SiPM_size / 2.0, SiPM_pitch * (SiPM_n_rows + 1) / 2 - SiPM_size / 2.0,
						2*SiPM_n_rows+1,
						-SiPM_pitch * (SiPM_n_rows + 1) / 2 + SiPM_size / 2.0, SiPM_pitch * (SiPM_n_rows + 1) / 2 - SiPM_size / 2.0);
  for (auto const& it : Npe_exp) {
    std::pair<double, double> XY = SiPM_channel_XY(it.first);
    if (XY.first == DBL_MAX || XY.second == DBL_MAX)
      continue;
    hist_01->Fill(XY.first, XY.second,  it.second == 0 ? 0.1*SiPM_Npe_stat_error/charge : it.second/charge);
    xy_avg.first += XY.first * it.second;
    xy_avg.second += XY.second * it.second;
    weight += it.second;
  }
  // Fix missing channels to avoid systematic x,y shift. Only cases of missing ch 43 and 44 are considered.
  bool is43 = (Npe_exp.find(43) != Npe_exp.end());
  bool is44 = (Npe_exp.find(44) != Npe_exp.end());
  if (!is43 && (Npe_exp.size() == 24 || (Npe_exp.size() == 23 && !is44))) { // mising only ch 43 or only channels 43 and 44
    double Npe43 = (Npe_exp[40] * Npe_exp[42] / Npe_exp[55] +
        Npe_exp[40] * Npe_exp[58] / Npe_exp[41]) / 2.0;
    std::pair<double, double> XY = SiPM_channel_XY(43);
    xy_avg.first += XY.first * Npe43;
    xy_avg.second += XY.second * Npe43;
    weight += Npe43;
  }
  if (!is44 && (Npe_exp.size() == 24 || (Npe_exp.size() == 23 && !is43))) { // mising only ch 44 or only channels 44 and 43
    double Npe44 = (Npe_exp[57] * Npe_exp[59] / Npe_exp[56] +
        0.5*(Npe_exp[57] + Npe_exp[59])) / 2.0;
    std::pair<double, double> XY = SiPM_channel_XY(44);
    xy_avg.first += XY.first * Npe44;
    xy_avg.second += XY.second * Npe44;
    weight += Npe44;
  }
  xy_avg.first /= weight;
  xy_avg.second /= weight;
  if (do_draw) {
    TCanvas *c_00 = new TCanvas ((std::string("01_") + "SiPM experimental S2 profile").c_str(), "SiPM experimental S2 profile");
    c_00->SetGrid(); c_00->SetTicks(); c_00->ToggleEventStatus(); c_00->ToggleToolBar();
    hist_01->SetStats(false);
    hist_01->GetXaxis()->SetTitle("X (mm)");
    hist_01->GetYaxis()->SetTitle("Y (mm)");
    double title_offset_01 = hist_01->GetYaxis()->GetTitleOffset() + 1.3;
    hist_01->GetYaxis()->SetTitleOffset(title_offset_01);
    //hist_01->Draw("lego");
    hist_01->Draw("colz");
    c_00->Update();
  }
  if (!do_draw)
    hist_01->Delete();
  return xy_avg;
}

std::pair<double, double> plot_simulation_XY(std::map<int, long int> Npe_exp, double charge, bool do_draw = true) {
  std::pair<double, double> xy_avg(0, 0);
  double weight = 0;
  TH2D* hist_01 = new TH2D("SiPM simulation S2 profile", "SiPM simulation S2 profile", 2*SiPM_n_rows+1,
						-SiPM_pitch * (SiPM_n_rows + 1) / 2 + SiPM_size / 2.0, SiPM_pitch * (SiPM_n_rows + 1) / 2 - SiPM_size / 2.0,
						2*SiPM_n_rows+1,
						-SiPM_pitch * (SiPM_n_rows + 1) / 2 + SiPM_size / 2.0, SiPM_pitch * (SiPM_n_rows + 1) / 2 - SiPM_size / 2.0);
  for (auto const& it : Npe_exp) {
    std::pair<double, double> XY = SiPM_channel_XY(it.first);
    if (XY.first == DBL_MAX || XY.second == DBL_MAX)
      continue;
    hist_01->Fill(XY.first, XY.second, it.second/charge);
    xy_avg.first += XY.first * it.second;
    xy_avg.second += XY.second * it.second;
    weight += it.second;
  }
  // Fix missing channels to avoid systematic x,y shift. Only cases of missing ch 43 and 44 are considered.
  bool is43 = (Npe_exp.find(43) != Npe_exp.end());
  bool is44 = (Npe_exp.find(44) != Npe_exp.end());
  if (!is43 && (Npe_exp.size() == 24 || (Npe_exp.size() == 23 && !is44))) { // mising only ch 43 or only channels 43 and 44
    double Npe43 = ((double)Npe_exp[40] * Npe_exp[42] / Npe_exp[55] +
        (double)Npe_exp[40] * Npe_exp[58] / Npe_exp[41]) / 2.0;
    std::pair<double, double> XY = SiPM_channel_XY(43);
    xy_avg.first += XY.first * Npe43;
    xy_avg.second += XY.second * Npe43;
    weight += Npe43;
  }
  if (!is44 && (Npe_exp.size() == 24 || (Npe_exp.size() == 23 && !is43))) { // mising only ch 44 or only channels 44 and 43
    double Npe44 = ((double)Npe_exp[57] * Npe_exp[59] / Npe_exp[56] +
        0.5*(Npe_exp[57] + Npe_exp[59])) / 2.0;
    std::pair<double, double> XY = SiPM_channel_XY(44);
    xy_avg.first += XY.first * Npe44;
    xy_avg.second += XY.second * Npe44;
    weight += Npe44;
  }
  xy_avg.first = (weight < SiPM_Npe_stat_error ? 0 : xy_avg.first/weight);
  xy_avg.second = (weight < SiPM_Npe_stat_error ? 0 : xy_avg.second/weight);
  if (do_draw) {
    TCanvas *c_00 = new TCanvas ((std::string("01_") + "SiPM simulation S2 profile").c_str(), "SiPM simulation S2 profile");
    c_00->SetGrid(); c_00->SetTicks(); c_00->ToggleEventStatus(); c_00->ToggleToolBar();
    hist_01->SetStats(false);
    hist_01->GetXaxis()->SetTitle("X (mm)");
    hist_01->GetYaxis()->SetTitle("Y (mm)");
    double title_offset_01 = hist_01->GetYaxis()->GetTitleOffset() + 1.3;
    hist_01->GetYaxis()->SetTitleOffset(title_offset_01);
    //hist_01->Draw("lego");
    hist_01->Draw("colz");
    c_00->Update();
  }
  if (!do_draw)
    hist_01->Delete();
  return xy_avg;
}

void plot_simulation_projection(std::map<int, long int> Npe_exp, double charge, bool project_on_X = true) {
  std::string name = std::string("SiPM simulation ") + (project_on_X ? "X" : "Y") + " projection";
  TH1D* hist_01 = new TH1D(name.c_str(), name.c_str(), 2*SiPM_n_rows+1,
						-SiPM_pitch * (SiPM_n_rows + 1) / 2 + SiPM_size / 2.0, SiPM_pitch * (SiPM_n_rows + 1) / 2 - SiPM_size / 2.0);
  for (auto const& it : Npe_exp) {
    std::pair<double, double> XY = SiPM_channel_XY(it.first);
    if (XY.first == DBL_MAX || XY.second == DBL_MAX)
      continue;
    hist_01->Fill(project_on_X ? XY.first : XY.second, it.second/charge);
  }
  TCanvas *c_00 = new TCanvas ((std::string("01_") + name).c_str(), name.c_str());
  c_00->SetGrid(); c_00->SetTicks(); c_00->ToggleEventStatus(); c_00->ToggleToolBar();
  hist_01->SetLineWidth(2);
  hist_01->GetXaxis()->SetTitle(project_on_X ? "X (mm)" : "Y (mm)");
  hist_01->GetYaxis()->SetTitle("Average S2 (PE/e)");
  double title_offset_01 = hist_01->GetYaxis()->GetTitleOffset() + 1.3;
  hist_01->GetYaxis()->SetTitleOffset(title_offset_01);
  hist_01->Draw("HIST");
  c_00->Update();
  return;
}

void plot_experiment_projection(std::map<int, double> Npe_exp, double charge, bool project_on_X = true) {
  std::string name = std::string("SiPM experiment ") + (project_on_X ? "X" : "Y") + " projection";
  TH1D* hist_01 = new TH1D(name.c_str(), name.c_str(), 2*SiPM_n_rows+1,
						-SiPM_pitch * (SiPM_n_rows + 1) / 2 + SiPM_size / 2.0, SiPM_pitch * (SiPM_n_rows + 1) / 2 - SiPM_size / 2.0);
  for (auto const& it : Npe_exp) {
    std::pair<double, double> XY = SiPM_channel_XY(it.first);
    if (XY.first == DBL_MAX || XY.second == DBL_MAX)
      continue;
    hist_01->Fill(project_on_X ? XY.first : XY.second, it.second/charge);
  }
  TCanvas *c_00 = new TCanvas ((std::string("01_") + name).c_str(), name.c_str());
  c_00->SetGrid(); c_00->SetTicks(); c_00->ToggleEventStatus(); c_00->ToggleToolBar();
  hist_01->SetLineWidth(2);
  hist_01->GetXaxis()->SetTitle(project_on_X ? "X (mm)" : "Y (mm)");
  hist_01->GetYaxis()->SetTitle("Average S2 (PE/e)");
  double title_offset_01 = hist_01->GetYaxis()->GetTitleOffset() + 1.3;
  hist_01->GetYaxis()->SetTitleOffset(title_offset_01);
  hist_01->Draw("HIST");
  c_00->Update();
  return;
}

// Using Oleynikov's position reconstruction results.
std::pair<double, double> weighted_xy_to_real(std::pair<double, double> xy_in) {
  std::pair<double, double> xy_out;
  xy_out.first = xy_in.first * 1.86737 + std::pow(xy_in.first, 2) * 0.00231 + std::pow(xy_in.first, 3) * 0.00291;
  xy_out.second = xy_in.second * 1.86737 + std::pow(xy_in.second, 2) * 0.00231 + std::pow(xy_in.second, 3) * 0.00291;
  return xy_out;
}

Npe_profile_result plot_Npe_profile(std::string sim_input_file = "../results/v10.0_elastic/2206V/recorded.dat", std::string exp_input_file = "../../../Post_processing/220804/results_v5/SiPM_Npes_Q0.20.txt",
  double experiment_charge = 1.27e6, double NBrS_yield_factor = 10, bool do_draw = true)
{
  Npe_profile_result res;

  gStyle->SetCanvasDefH(800);
	gStyle->SetCanvasDefW(1000);
  gErrorIgnoreLevel = 1001;
  gStyle->SetStatY(0.9);
	gStyle->SetStatX(0.9);
  gStyle->SetStatH(0.3);
  gStyle->SetGridStyle(3);
	gStyle->SetGridColor(14);
	gStyle->SetGridWidth(1);
	gStyle->SetOptStat("mr");

  TString input_str(sim_input_file);
  TObjArray *subStrL = TPRegexp("/[0-9]+.{0,1}[0-9]+k{0,1}[vV]/").MatchS(input_str);
  TObjString *match = (TObjString *)subStrL->Last();
  if (match == nullptr) {
    std::cerr<<"Could not determine voltage! Aborting."<<std::endl;
    return res;
  }
  std::string Voltage(match->GetString().Data());

  Voltage.erase(Voltage.end()-2, Voltage.end());
  Voltage.erase(Voltage.begin(), Voltage.begin()+1);
  if (Voltage.back() == 'k')
    Voltage.erase(Voltage.end()-1, Voltage.end());

  std::string plot_name = "SiPM matrix S2 radius profile";
  std::string X_axis_title = "R (mm)";
  std::string Y_axis_title = "S2 amplitude (PE/e)";
  std::pair<double, double> X_axis_range(0, 40.0);
  std::pair<double, double> Y_axis_range(0.0, 0.0);
  std::pair<double, double> X_axis_zoom(0, 0);
  std::pair<double, double> Y_axis_zoom(0, 0);
  int X_n_bins = 60;
  ExperimentalXY* SiPM_QE = LoadQE("../data/quantum_efficiency/SiPM_s13360-6050pe_46V.dat", true);

  PlotInfo plot_info;
  plot_info.selection = SelectSiPMs;
  plot_info.histogram = NULL;
  plot_info.plot_par_x = PlotParameter::None; //PlotParameter::gammaNProfile;
  plot_info.plot_par_y = PlotParameter::None;
  //TH1D* hist_00 = new TH1D("hist_00", plot_name.c_str(), X_n_bins, X_axis_range.first, X_axis_range.second);
  //plot_info.histogram = hist_00;

  std::ifstream str;
  str.open(sim_input_file);
  if (!str.is_open()) {
    std::cerr<<"Could not open file \""<<sim_input_file<<"\"!"<<std::endl;
  }
  FileToHist(plot_info, str);
  res.Npe_per_e_per_ch_exp = ReadExperimentData(exp_input_file, Voltage); // Naturally, channels in experiment must be consistent with definitions in init.cpp
  double QE_avg = get_SiPM_QE(str, SiPM_QE);

  std::vector<int> used_channels;
  for (std::size_t i = 0; i < SiPM_n_rows*SiPM_n_rows; ++i) {
    int channel = SiPM_index_to_channel(i);
    if (channel < 0)
      continue;
    std::vector<int> eq_chs = SiPM_equidistant_channels(channel);
    bool intersection = false;
    for (std::size_t i1 = 0; i1 != used_channels.size() && !intersection; ++i1) {
      for (std::size_t i2 = 0; i2 != eq_chs.size() && !intersection; ++i2)
        if (used_channels[i1] == eq_chs[i2])
          intersection = true;
    }
    if (!intersection) { // Creating unique equidistant channel groups.
      for (std::size_t i2 = 0; i2 != eq_chs.size() && !intersection; ++i2)
        used_channels.push_back(eq_chs[i2]);
      res.equidistant_ch_groups.push_back(eq_chs);
      res.Radii.push_back(SiPM_channel_R(channel));
      double Npe_avg_sim = 0;
      double Npe_avg_sim_with_QE = 0;
      double Npe_avg_exp = 0;
      std::size_t N_exp = 0;
      for (std::size_t ch_i = 0; ch_i != eq_chs.size(); ++ch_i) {
        channel = eq_chs[ch_i];
        int Npe_ch = GetNpeCh(plot_info, channel);
        Npe_avg_sim += Npe_ch;
        if (do_draw)
          Npe_avg_sim_with_QE += Npe_ch * get_SiPM_QE(str, SiPM_QE, channel);
        auto e = res.Npe_per_e_per_ch_exp.find(channel);
        if (e != res.Npe_per_e_per_ch_exp.end()) {
          Npe_avg_exp += res.Npe_per_e_per_ch_exp[channel];
          N_exp += 1;
        }
      }
      Npe_avg_sim /= eq_chs.size();
      if (do_draw)
        Npe_avg_sim_with_QE /= eq_chs.size();
      else
        Npe_avg_sim_with_QE = Npe_avg_sim;
      Npe_avg_exp = (N_exp == 0 ? -1 : Npe_avg_exp/N_exp);
      res.Npes_per_e_sim.push_back(Npe_avg_sim * QE_avg / plot_info.N_electrons / NBrS_yield_factor);
      res.Npes_per_e_sim_indiv_QE.push_back(Npe_avg_sim_with_QE / plot_info.N_electrons / NBrS_yield_factor);
      res.Npes_per_e_exp.push_back(Npe_avg_exp / experiment_charge);
      res.n_ch_exp.push_back(N_exp);
      Y_axis_range.second = std::max(Npe_avg_sim * QE_avg / plot_info.N_electrons / NBrS_yield_factor * 1.3, Y_axis_range.second);
      Y_axis_range.second = std::max(Npe_avg_exp / experiment_charge * 1.3, Y_axis_range.second);
    }
  }
  //
  //double normalization = 1;
  //for (std::size_t i = 0, i_end_ = Radii.size(); i!=i_end_; ++i) {
  //  if (Radii[i] < SiPM_pitch*1e-3) {
  //    normalization = Npes_sim[i]/Npes_exp[i];
  //    break;
  //  }
  //}
  if (do_draw) {
    std::cout<<"Total Npe per e:\n";
    double total_sim = 0, total_sim_indiv_QE = 0, total_exp = 0;
    for (std::size_t i = 0, i_end_ = res.Radii.size(); i!=i_end_; ++i) {
      total_sim += res.Npes_per_e_sim[i] * res.n_ch_exp[i];
      total_sim_indiv_QE += res.Npes_per_e_sim_indiv_QE[i] * res.n_ch_exp[i];
      total_exp += res.Npes_per_e_exp[i] * res.n_ch_exp[i];
    }
    std::cout<<"\tTheory, average PDE: "<<total_sim<<"\n";
    std::cout<<"\tTheory, individual PDE: "<<total_sim_indiv_QE<<"\n";
    std::cout<<"\tExperiment: "<<total_exp<<std::endl;

    std::pair<double, double> S2_center = plot_experiment_XY(res.Npe_per_e_per_ch_exp, experiment_charge, true);
    std::cout<<"Experimental source center, not adjusted (x, y) = ("<<S2_center.first <<", "<<S2_center.second <<")"<<std::endl;
    res.source_xy_exp = weighted_xy_to_real(S2_center);
    std::cout<<"Experimental source center, adjusted (x, y) = ("<<res.source_xy_exp.first <<", "<<res.source_xy_exp.second <<")"<<std::endl;

    S2_center = plot_simulation_XY(plot_info.Npes, plot_info.N_electrons * NBrS_yield_factor / QE_avg, true);
    std::cout<<"Simulation source center, not adjusted (x, y) = ("<<S2_center.first <<", "<<S2_center.second <<")"<<std::endl;
    S2_center = weighted_xy_to_real(S2_center);
    std::cout<<"Simulation source center, adjusted (x, y) = ("<<S2_center.first <<", "<<S2_center.second <<")"<<std::endl;

    plot_simulation_projection(plot_info.Npes, plot_info.N_electrons * NBrS_yield_factor / QE_avg, true);
    plot_experiment_projection(res.Npe_per_e_per_ch_exp, experiment_charge, true);

    TCanvas *c_00 = new TCanvas ((std::string("00_") + plot_name).c_str(), plot_name.c_str());
  	c_00->SetGrid(); c_00->SetTicks(); c_00->ToggleEventStatus(); c_00->ToggleToolBar();
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
    if (plot_info.histogram)
      plot_info.histogram->Draw("sames");
    c_00->Update();

    TLegend *legend = new TLegend(0.60, 0.75, 0.9, 0.9);
  	legend->SetHeader((std::string("V_{T} = ") + Voltage + "V").c_str());
  	legend->SetMargin(0.25);
    TGraph *graph_00 = GraphFromData(res.Radii, res.Npes_per_e_sim);
    if (nullptr != graph_00) {
      graph_00->SetMarkerColor(kRed);
      graph_00->SetMarkerSize(2);
      graph_00->SetMarkerStyle(8);
      graph_00->Draw("Psame");
      legend->AddEntry(graph_00, (std::string("Theory, average PDE")).c_str(), "p");
    }
    TGraph *graph_02 = GraphFromData(res.Radii, res.Npes_per_e_sim_indiv_QE);
    if (nullptr != graph_02) {
      graph_02->SetMarkerColor(kBlue);
      graph_02->SetMarkerSize(2);
      graph_02->SetMarkerStyle(22);
      graph_02->Draw("Psame");
      legend->AddEntry(graph_02, (std::string("Theory, individual PDE")).c_str(), "p");
    }
    TGraph *graph_01 = GraphFromData(res.Radii, res.Npes_per_e_exp);
    if (nullptr != graph_01) {
      graph_01->SetMarkerColor(kBlack);
      graph_01->SetMarkerSize(2);
      graph_01->SetMarkerStyle(25);
      graph_01->Draw("Psame");
      legend->AddEntry(graph_01, (std::string("Experiment")).c_str(), "p");
    }
    // Place to modify displayed stats.
    // the following line is needed to avoid that the automatic redrawing of stats
    if (plot_info.histogram)
      plot_info.histogram->SetStats(false);
    if (graph_01 || graph_00)
      legend->Draw("same");
    c_00->Update();
  } else {
    std::pair<double, double> S2_center = plot_experiment_XY(res.Npe_per_e_per_ch_exp, experiment_charge, false);
    res.source_xy_exp = weighted_xy_to_real(S2_center);
    if (plot_info.histogram)
      plot_info.histogram->Delete();
  }
  if (SiPM_QE)
    delete SiPM_QE;
  return res;
}
