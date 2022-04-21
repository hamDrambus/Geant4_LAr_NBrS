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

double ev_to_nm(double energy) {
  return 1239.84198 / energy;
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

struct PlotInfo {
  TH1* histogram;
  PlotParameter plot_par_x;
  PlotParameter plot_par_y;
  Selector selection;
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
    double size_SiPM = 6.0; // [mm]
	  double SiPM_spacing = 10; // [mm]
    int x_index = (int) (photon.pos_x)/(SiPM_spacing/2.0);
    int y_index = (int) (photon.pos_y)/(SiPM_spacing/2.0);
    double x_pos = SiPM_spacing*x_index;
    double y_pos = SiPM_spacing*y_index;
    return std::sqrt(x_pos*x_pos + y_pos*y_pos);
  }
  default:
    return -DBL_MAX;
  }
}

void FillHist(PlotInfo& plot_info, ElectronInfo& electron, PhotonInfo& photon)
{
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
