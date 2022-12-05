enum PlotParameter
{
  None,
  gammaXpos,
  gammaYpos,
  gammaZpos,
  gammaEnergy,
  gammaSpectrum,
  gammaTime,
  gammaAngle,
  gammaNProfile
};

const int SiPM_n_rows = 5;
const double SiPM_size = 6; // mm
const double SiPM_pitch = 10; // mm
const double SiPM_matrix_center_x = 0, SiPM_matrix_center_y = 0; // mm

std::string dbl_to_str (double val, int precision)
{
	std::stringstream ss;
	ss<<std::fixed<<std::setprecision(precision)<<val;
	return ss.str();
}

double ev_to_nm(double energy) {
  return 1239.84198 / energy;
}

int SiPM_index_to_channel (int index) {
  std::map<int, int> index_to_ch = {
    {0, 42},
    {1, 43},
    {2, 58},
    {3, 59},
    {4, 44},
    {5, 55},
    {6, 40},
    {7, 41},
    {8, 56},
    {9, 57},
    {10, 52},
    {11, 53},
    {12, 38},
    {13, 39},
    {14, 54},
    {15, 35},
    {16, 50},
    {17, 51},
    {18, 36},
    {19, 37},
    {20, 32},
    {21, 33},
    {22, 48},
    {23, 49},
    {24, 34} };
  if (index >=0 && index < 25) {
    return index_to_ch[index];
  }
  return -1;
}

// prunes list of channels to uniques and only those that exist
std::vector<int> SiPM_valid_channels_only (std::vector<int> channels) {
  std::sort(channels.begin(), channels.end());
  auto last = std::unique(channels.begin(), channels.end());
  channels.erase(last, channels.end());
  std::vector<int> valid_channels;
  for (std::size_t i = 0; i < SiPM_n_rows*SiPM_n_rows; ++i) {
    int channel = SiPM_index_to_channel(i);
    if (channel < 0)
      continue;
    for (int ch : channels) {
      if (ch == channel) {
        valid_channels.push_back(channel);
        break;
      }
    }
  }
  return valid_channels;
}

// Returns distatance from (0, 0) to SiPM #channel
double SiPM_channel_R(int channel) {
  int index = -1;
  for (std::size_t i = 0; i < SiPM_n_rows*SiPM_n_rows; ++i) {
    if (channel == SiPM_index_to_channel(i)) {
      index = i;
      break;
    }
  }
  if (index < 0)
    return -1.0;
  // Same as in source/geant4/detector/DetectorParameterisation.cpp (DetectorParameterisation::ComputeTransformation)
  double X_position = (-SiPM_n_rows / 2.0 + index % SiPM_n_rows + 0.5) * SiPM_pitch;
	double Y_position = (-SiPM_n_rows / 2.0 + (index / SiPM_n_rows) % SiPM_n_rows + 0.5) * SiPM_pitch;
  return sqrt(X_position*X_position + Y_position*Y_position);
}

std::vector<int> SiPM_equidistant_channels(int channel) {
  std::vector<int> out;
  double R = SiPM_channel_R(channel);
  if (R < 0)
    return out;
  for (std::size_t i = 0; i < SiPM_n_rows*SiPM_n_rows; ++i) {
    int ch = SiPM_index_to_channel(i);
    double r = SiPM_channel_R(ch);
    if (fabs(r-R) < 0.0001 * SiPM_pitch)
      out.push_back(ch);
  }
  return out;
}

class ElectronInfo {
public:
  int index;
  double pos_x; //position
  double pos_y;
  double pos_z;
  long seed;
  unsigned long photon_number;
  ElectronInfo() : index(-1), pos_x(0), pos_y(0), pos_z(0),
    seed(-1), photon_number(0) {}
};

class PhotonInfo {
public:
  double energy;
  double pos_x; //position
  double pos_y;
  double pos_z;
  double time;
  double mom_x; //momentum
  double mom_y;
  double mom_z;
  PhotonInfo() : energy(0), pos_x(0), pos_y(0), pos_z(0),
    time(0), mom_x(0), mom_y(0), mom_z(0) {}
};

class ExperimentalXY {
public:
  ROOT::Math::Interpolator interpolator;
  std::vector<double> Xs;
  std::vector<double> Ys;
  ExperimentalXY(ROOT::Math::Interpolation::Type type = ROOT::Math::Interpolation::Type::kLINEAR) : interpolator(0, type) {}
  double Eval(double x) {
    if (Xs.empty() || Ys.size()!=Ys.size())
      return 0.0;
    if (x < Xs[0] || x > Xs[Xs.size()-1])
      return 0.0;
    double out = interpolator.Eval(x);
    if (isnan(out) || isinf(out)) {
      std::cerr<<"ExperimentalXY::Eval("<<x<<") = "<<out<<std::endl;
      out = 0.0;
    }
    return out;
  }
};

std::string strtoken(std::string &in, std::string break_symbs)
{
  std::string out;
  while (!in.empty())
  {
    char a = in.front();
    in.erase(in.begin());
    bool break_ = false;
    for (auto h = break_symbs.begin(); h != break_symbs.end(); ++h)
      if (a == *h) {
        break_ = true;
        break;
      }
    if ((break_) && (out.empty()))
      continue;
    if (break_)
      return out;
    out.push_back(a);
  }
  return out;
}

typedef bool (*Selector)(ElectronInfo& electron, PhotonInfo& photon); // Add to plot if true

bool SelectPMTs (ElectronInfo& electron, PhotonInfo& photon)
{
  if (photon.pos_x > 75 || photon.pos_x < -75 || photon.pos_y > 75 || photon.pos_y < -75) {
    return true;
  }
  return false;
}

bool SelectSiPMs (ElectronInfo& electron, PhotonInfo& photon)
{
  if (photon.pos_x > 75 || photon.pos_x < -75 || photon.pos_y > 75 || photon.pos_y < -75) {
    return false;
  }
  return true;
}

bool SelectAll (ElectronInfo& electron, PhotonInfo& photon)
{
  return true;
}

int PhotonToChannel (ElectronInfo& electron, PhotonInfo& photon) {
  const double max_XY = 0.5 * (SiPM_n_rows * SiPM_pitch);
  const double PMMA_plate_center_z = 78.15;

  if (SelectPMTs(electron, photon)) { // channels as in geant4
    if (photon.pos_x > 75 && photon.pos_y < 75 && photon.pos_y > -75)
      return 1;
    if (photon.pos_x < -75 && photon.pos_y < 75 && photon.pos_y > -75)
      return 0;
    if (photon.pos_y < -75 && photon.pos_x < 75 && photon.pos_x > -75)
      return 2;
    if (photon.pos_y > 75 && photon.pos_x < 75 && photon.pos_x > -75)
      return 3;
  }
  if (SelectSiPMs(electron, photon)) { // channels as in experiment
    if (!((photon.pos_y - SiPM_matrix_center_y) > max_XY || (photon.pos_y - SiPM_matrix_center_y) < -max_XY
        || (photon.pos_x - SiPM_matrix_center_x)> max_XY || (photon.pos_x - SiPM_matrix_center_x) < -max_XY
        || photon.pos_z < PMMA_plate_center_z)) {
      int index_x = (int)((photon.pos_x - SiPM_matrix_center_x + max_XY) / SiPM_pitch);
      int index_y = (int)((photon.pos_y - SiPM_matrix_center_x + max_XY) / SiPM_pitch);
      int index = index_y * SiPM_n_rows + index_x;
      return SiPM_index_to_channel(index);
    }
  }
  return -1;
}

struct PlotInfo {
  TH1* histogram;
  PlotParameter plot_par_x;
  PlotParameter plot_par_y;
  Selector selection;
  std::map<int, long int> Npes; //per channel. ch 0-3 - PMTs, ch 32-63 - SiPMs
  unsigned int N_electrons;
};

double PickValue(PlotParameter par, const ElectronInfo& electron, const PhotonInfo& photon)
{
  switch(par) {
  case PlotParameter::gammaXpos: {
    return photon.pos_x;
  }
  case PlotParameter::gammaYpos: {
    return photon.pos_y;
  }
  case PlotParameter::gammaZpos: {
    return photon.pos_z;
  }
  case PlotParameter::gammaEnergy: {
    return photon.energy;
  }
  case PlotParameter::gammaSpectrum: {
    return ev_to_nm(photon.energy);
  }
  case PlotParameter::gammaTime: {
    return photon.time;
  }
  case PlotParameter::gammaAngle: {
    double mom = photon.mom_x*photon.mom_x + photon.mom_y*photon.mom_y + photon.mom_z*photon.mom_z;
    mom = std::sqrt(mom);
    if (photon.pos_x > 75 || photon.pos_x < -75) { // hit PMT, X axis normal
      double mom_yz =photon.mom_y*photon.mom_y + photon.mom_z*photon.mom_z;
      mom_yz = std::sqrt(mom_yz);
      return std::asin(mom_yz/mom);
    }
    if (photon.pos_y > 75 || photon.pos_y < -75) { // hit PMT, Y axis normal
      double mom_xz =photon.mom_x*photon.mom_x + photon.mom_z*photon.mom_z;
      mom_xz = std::sqrt(mom_xz);
      return std::asin(mom_xz/mom);
    }
    // hit SiPM, Z axis normal
    double mom_xy =photon.mom_x*photon.mom_x + photon.mom_y*photon.mom_y;
    mom_xy = std::sqrt(mom_xy);
    return std::asin(mom_xy/mom);
  }
  case PlotParameter::gammaNProfile: { //returns distance from hit SiPM to SiPM-matrix center
    if (photon.pos_x > 75 || photon.pos_x < -75 || photon.pos_y > 75 || photon.pos_y < -75) {
      return -DBL_MAX; // Hit PMT
    }
    int x_index = (int) (photon.pos_x)/(SiPM_pitch/2.0);
    int y_index = (int) (photon.pos_y)/(SiPM_pitch/2.0);
    double x_pos = SiPM_pitch*x_index;
    double y_pos = SiPM_pitch*y_index;
    return std::sqrt(x_pos*x_pos + y_pos*y_pos);
  }
  default:
    return -DBL_MAX;
  }
}

void FillHist(PlotInfo& plot_info, ElectronInfo& electron, PhotonInfo& photon)
{
  int ch = PhotonToChannel(electron, photon);
  auto e = plot_info.Npes.find(ch);
  if (e == plot_info.Npes.end()) {
    plot_info.Npes[ch] = 0;
  }
  ++plot_info.Npes[ch];

  if (plot_info.selection==NULL || plot_info.selection(electron, photon)) {
    if (PlotParameter::None == plot_info.plot_par_x)
      plot_info.plot_par_x = plot_info.plot_par_y;
    if (PlotParameter::None == plot_info.plot_par_x)
      return;
    if (PlotParameter::None == plot_info.plot_par_y) { // 1D histogram
      TH1D* hist = (TH1D*) plot_info.histogram;
      hist->Fill(PickValue(plot_info.plot_par_x, electron, photon), 1.0);
      return;
    }
    // 2D histogram
    TH2D* hist = (TH2D*) plot_info.histogram;
    hist->Fill(PickValue(plot_info.plot_par_x, electron, photon),
          PickValue(plot_info.plot_par_y, electron, photon));
  }
}

bool FileToHist(PlotInfo& plot_info, std::ifstream &str)
{
  if (!str.is_open())
    return false;
  plot_info.Npes.clear();
  plot_info.N_electrons = 0;
  std::string line, word;
  ElectronInfo electron;
  PhotonInfo photon;
	int line_n = 0;
  unsigned long photon_N = 0;
	while (!str.eof() && str.is_open()) {
		std::getline(str, line);
		++line_n;
		if (line.size() >= 2) //Ignore simple c style comment
			if ((line[0] == '/') && (line[1] == '/'))
				continue;
		try {
      if (0 >= photon_N) { //read electron info
        ElectronInfo electron_temp;
        electron = electron_temp;
        photon_N = 0;

  			word = strtoken(line, "\t ");
  			electron_temp.index = std::stoi(word);
  			word = strtoken(line, "\t ");
        electron_temp.photon_number = std::stoul(word);
        word = strtoken(line, "\t ");
  			electron_temp.pos_x = std::stod(word);
        word = strtoken(line, "\t ");
  			electron_temp.pos_y = std::stod(word);
        word = strtoken(line, "\t ");
  			electron_temp.pos_z = std::stod(word);
        word = strtoken(line, "\t ");
        word.erase(std::remove(word.begin(), word.end(), '"'), word.end());
        electron_temp.seed = std::stol(word);

        electron = electron_temp;
        photon_N = electron.photon_number;
        ++plot_info.N_electrons;
      } else { //read each photon for single electron
        photon_N--;
        PhotonInfo photon_temp;
        photon = photon_temp;

        word = strtoken(line, "\t ");
  			photon_temp.energy = std::stod(word);
        word = strtoken(line, "\t ");
  			photon_temp.pos_x = std::stod(word);
        word = strtoken(line, "\t ");
  			photon_temp.pos_y = std::stod(word);
        word = strtoken(line, "\t ");
  			photon_temp.pos_z = std::stod(word);
        word = strtoken(line, "\t ");
  			photon_temp.time = std::stod(word);
        word = strtoken(line, "\t ");
  			photon_temp.mom_x = std::stod(word);
        word = strtoken(line, "\t ");
  			photon_temp.mom_y = std::stod(word);
        word = strtoken(line, "\t ");
  			photon_temp.mom_z = std::stod(word);

        photon = photon_temp;
        FillHist(plot_info, electron, photon);
      }
		} catch (std::invalid_argument &e) {
    } catch (std::out_of_range &e) {
			std::cerr << "FileToHist: Invalid data on line " << line_n << std::endl;
			std::cerr << e.what() << std::endl;
			return false;
		} catch (std::exception &e) {
			std::cerr << "FileToHist: Unforeseen exception on line " << line_n << std::endl;
			std::cerr << e.what() << std::endl;
			return false;
		}
	}
  return true;
}

//column number starts from 0
bool LoadColumnData(std::string file, std::vector<double> &Vs, std::size_t y_column)
{
	Vs.clear();
	std::ifstream str;
	str.open(file);
	if (!str.is_open()) {
		std::cerr << "Error: Failed to open file \"" << file << "\"!" << std::endl;
		return false;
	}
	std::string line, word;
	int line_n = 0;
	while (!str.eof()) {
		std::getline(str, line);
		++line_n;
		if (line.size() >= 2) //Ignore simple c style comment
			if ((line[0] == '/') && (line[1] == '/'))
				continue;
		std::size_t column = 0;
		word = strtoken(line, "\t ");
		while (column < y_column && !word.empty()) {
			word = strtoken(line, "\t ");
			++column;
		}
		if (word.empty())
			continue;
		double val = std::stod(word);
		Vs.push_back(val);
	}
	return true;
}

ExperimentalXY* LoadQE (std::string file, bool in_eVs) {
  ExperimentalXY* out = new ExperimentalXY(ROOT::Math::Interpolation::Type::kLINEAR);
  if (!LoadColumnData(file, out->Xs, 0) || !LoadColumnData(file, out->Ys, 1)) {
    delete out;
    return nullptr;
  }
  if (in_eVs) {
    for (std::size_t i = 0, i_end_=out->Xs.size(); i!=i_end_; ++i) {
      out->Xs[i] = ev_to_nm(out->Xs[i]);
    }
    for (std::size_t i = 1, i_end_=out->Xs.size(); i!=i_end_ && i_end_!=0; ++i) {
      if (out->Xs[i] <= out->Xs[i-1]) {
        std::cerr<<"LoadQE:Non increasing Xs: i="<<i<<" Xs[i]="<<out->Xs[i]<<" Xs[i-1]="<<out->Xs[i-1]<<std::endl;
      }
    }
  }
  out->interpolator.SetData(out->Xs, out->Ys);
  return out;
}

TGraph* GraphFromData(ExperimentalXY* data, double scaleX = 1.0, double scaleY = 1.0) {
  TGraph* out = nullptr;
  Int_t size = std::min(data->Xs.size(), data->Ys.size());
	Double_t *xs = NULL, *ys = NULL;
	if (size > 0) {
		xs = new Double_t[size];
		ys = new Double_t[size];
		for (Int_t i = 0; i != size; ++i) {
			xs[i] = data->Xs[i] * scaleX;
			ys[i] = data->Ys[i] * scaleY;
		}
		out = new TGraph(size, xs, ys);
		delete [] xs;
		delete [] ys;
	}
  return out;
}

double CalculateNpe(TH1D* spectrum_counts, ExperimentalXY* QEdata = nullptr) {
  double out = 0.0;
  if (spectrum_counts == nullptr) {
    std::cout<<"CalculateNpe: Invalid input"<<std::endl;
    return out;
  }
  if (QEdata == nullptr || QEdata->Xs.empty()) {
    for (int bin = 1, bin_end = spectrum_counts->GetNbinsX()+1; bin!=bin_end; ++bin) {
      double N_counts = spectrum_counts->GetBinContent(bin);
      out += N_counts;
    }
    return out;
  }
  for (int bin = 1, bin_end = spectrum_counts->GetNbinsX()+1; bin!=bin_end; ++bin) {
    double lambda = spectrum_counts->GetBinCenter(bin);
    double N_counts = spectrum_counts->GetBinContent(bin);
    double QE = QEdata->Eval(lambda);
    out += QE * N_counts;
  }
  return out;
}

void AddStatValue(TPaveStats *ps, std::string name, double value) {
  if (nullptr == ps)
    return;
  ps->SetName("mystats");
  TList *listOfLines = ps->GetListOfLines();
  // Add a new line in the stat box.
  // Note that "=" is a control character
  //std::string format = gStyle->GetStatFormat();
  std::string format = ps->GetStatFormat();
  TString line = TString::Format((name + "=%" + format).c_str(), value);
  TLatex *myt = new TLatex(0,0,line.Data());
  myt->SetTextColor(gStyle->GetStatColor());
  myt->SetTextFont(gStyle->GetStatFont());
  myt->SetTextSize(gStyle->GetStatFontSize());
  listOfLines->Add(myt);
  return;
}

int GetNpeCh(PlotInfo& plot_info, int ch) {
  auto e = plot_info.Npes.find(ch);
  if (e == plot_info.Npes.end()) {
    return 0;
  }
  return plot_info.Npes[ch];
}

int GetNpePMTraw (PlotInfo& plot_info) {
  std::size_t out = GetNpeCh(plot_info, 0) + GetNpeCh(plot_info, 2) + GetNpeCh(plot_info, 3) + GetNpeCh(plot_info, 1);
  return out;
}

int GetNpePMTavg (PlotInfo& plot_info) {
  std::size_t PMT3 = GetNpeCh(plot_info, 0) + GetNpeCh(plot_info, 2) + GetNpeCh(plot_info, 3);
  double grid_fraction = (double) GetNpeCh(plot_info, 1) * 3.0 / PMT3;
  if (PMT3==0)
    grid_fraction = 1.0;
  std::size_t out = std::round((PMT3 * grid_fraction + GetNpeCh(plot_info, 1)) / 4.0);
  return out;
}

int GetNpeSiPMs (PlotInfo& plot_info) {
  std::size_t out = 0;
  for (int i = 0; i != SiPM_n_rows*SiPM_n_rows; ++i) {
    out += GetNpeCh(plot_info, SiPM_index_to_channel(i));
  }
  return out;
}

int GetNpeSiPMsSome (PlotInfo& plot_info, std::vector<int> exclude_ch) {
  std::size_t out = GetNpeSiPMs(plot_info);
  exclude_ch = SiPM_valid_channels_only(exclude_ch);
  double rem = 0;
  for (int channel : exclude_ch) {
    std::vector<int> chs = SiPM_equidistant_channels(channel);
    double average = 0;
    for (int ch : chs)
      average += GetNpeCh(plot_info, ch);
    rem += average/chs.size();
  }
  out -= std::round(rem);
  return out;
}

int GetNpeSiPMs23 (PlotInfo& plot_info) {
  std::vector<int> exclude_chs = {44, 43};
  return GetNpeSiPMsSome(plot_info, exclude_chs);
}
