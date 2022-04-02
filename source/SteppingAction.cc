#include "SteppingAction.hh"

SteppingAction::SteppingAction(DetectorConstruction* myDC)
	: myDetector(myDC), is_photon_inside_single_THGEM_hole(false)
{
	delta_z = 0.01;
}

void SteppingAction::UserSteppingAction(const G4Step* theStep)
{
	G4Track* theTrack = theStep->GetTrack();
	G4int trackID = theTrack->GetTrackID();
	G4ParticleDefinition* particleType = theTrack->GetDefinition();

	G4StepPoint* thePrePoint = theStep->GetPreStepPoint();
	G4VPhysicalVolume* thePrePV = thePrePoint->GetPhysicalVolume();

	G4StepPoint* thePostPoint = theStep->GetPostStepPoint();
	G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();

	G4ThreeVector pos;
	std::string postVolName;
	std::string prevVolName;
	G4ThreeVector MomentumDirection;
	if (thePostPV != NULL) {
		pos = theStep->GetPostStepPoint()->GetPosition();
		postVolName = theStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();
		MomentumDirection = theStep->GetTrack()->GetMomentumDirection();
	}
	if (thePrePV != nullptr)
	  prevVolName = thePrePV->GetName();

	//std::cout<<"************"<<std::endl;
	//std::cout<<"PrevVolume: \""<<prevVolName<<"\""<<std::endl;
	//std::cout<<"PostVolume: \""<<postVolName<<"\""<<std::endl;
	return G4UserSteppingAction::UserSteppingAction(theStep);
}

void SteppingAction::ReturnFromSingleTHGEMHole(const G4Step* theStep)
{
	G4double z_size;
	if (Volume_init_THGEM == gPars::det_dims.THGEM1_cell_container_name)
		z_size = gPars::det_dims.THGEM1_width_total;

	//final position inside SingleTHGEMHole
	const G4ThreeVector& pos_f_h = theStep->GetPostStepPoint()->GetPosition();
	const G4ThreeVector& mom_f_h = theStep->GetTrack()->GetMomentumDirection();
	bool is_bottom_to_top = (pos_f_h.z() == gPars::det_dims.xyz_position_SingleTHGEMHole + z_size / 2.0) && mom_f_h.z() > 0;
	bool is_top_to_bottom = (pos_f_h.z() == gPars::det_dims.xyz_position_SingleTHGEMHole - z_size / 2.0) && mom_f_h.z() < 0;
	//if photon go out from SingleTHGEMHole
	if (is_photon_inside_single_THGEM_hole && (is_bottom_to_top || is_top_to_bottom)) {
		is_photon_inside_single_THGEM_hole = false;

		//find the x,y shift in the SingleTHGEMHole
		double dx = pos_f_h.x() - Position_init_SingleTHGEMHole.x();
		double dy = pos_f_h.y() - Position_init_SingleTHGEMHole.y();

		const G4ThreeVector& MomentumDirection_f_THGEM = theStep->GetTrack()->GetMomentumDirection();
		double zm = MomentumDirection_f_THGEM.z();
		double zm_sign = zm > 0 ? 1 : -1;

		bool is_changed_direction = (MomentumDirection_init_THGEM.z() > 0) == !(MomentumDirection_f_THGEM.z()>0); /* a==!b means XOR*/
		double dz = is_changed_direction ? 0 : z_size*zm_sign;

		const G4ThreeVector& pos_f_THGEM =
			G4ThreeVector(Position_init_THGEM.x() + dx, Position_init_THGEM.y() + dy,
				Position_init_THGEM.z() + dz - zm_sign*delta_z);

		/*cout << "zm = " << zm << endl;
		cout << "zm_sign = " << zm_sign << endl;
		cout << "is_changed_direction = " << is_changed_direction << endl;
		cout << "dz = " << dz << endl;*/

		theStep->GetTrack()->SetPosition(pos_f_THGEM);
	}
}

void SteppingAction::PassThroughGEM(const G4Step* theStep, G4double z_pos, G4double z_size)
{
	//change MomentumDirection and position for optical photon incidents on THGEM
	const G4ThreeVector& pos_i = theStep->GetPostStepPoint()->GetPosition();
	if (pos_i.z() > (z_pos - delta_z) && pos_i.z() < (z_pos + delta_z)) {
		//get
		const G4ThreeVector& MomentumDirection_i = theStep->GetTrack()->GetMomentumDirection();
		const G4ThreeVector& Position_i = theStep->GetTrack()->GetPosition();

		//THGEM CERN 28%
		//const double step = 0.9 * mm;
		//const double radius = 0.25 * mm;

		////THGEM ELECROCONNECT 75%
		//const double step = 1.1 * mm;
		//const double radius = 0.5 * mm;

		double step = gPars::det_dims.THGEM1_hole_pitch;

		const double step_y = step * sqrt(3) / 2;
		double y_center;
		double x_center;

		double x = Position_i.x();
		double y = Position_i.y();
		double x_abs = fabs(Position_i.x());
		double y_abs = fabs(Position_i.y());

		// the nearest y_center
		if (y_abs < step_y / 2) {
			y_center = 0;
		} else {
			y_center = ((int)((y_abs - step_y / 2) / step_y) + 1) * step_y * y / y_abs;
		}
		//the nearest x_center
		if (y_abs < step_y / 2) {
			int n_x = (int)(x / step);
			int sign = x > 0 ? 1 : -1;
			x_center = (n_x + sign * 0.5) * step;
		} else {
			if (((int)((y_abs - (step_y / 2)) / (step_y)) + 1) % 2 == 1 /*if y row is odd*/) {
				if (x_abs < step / 2.0) {
					x_center = 0;
				} else {
					int n_x = (x_abs - step / 2) / step;
					int sign = x > 0 ? 1 : -1;
					x_center = (n_x + 1) * step * sign;
				}
			} else {
				int n_x = (int)(x / step);
				int sign = x > 0 ? 1 : -1;
				x_center = (n_x + sign * 0.5) * step;
			}
		}

		//kill or allow photon to pass
		double distance = sqrt(pow(x - x_center, 2.0) + pow(y - y_center, 2.0));
		if (distance < gPars::det_dims.THGEM1_hole_radius) {
			//set
			//const G4ThreeVector& MomentumDirection_f = G4ThreeVector(0, 1, 0);
			double xm = MomentumDirection_i.x();
			double ym = MomentumDirection_i.y();
			double zm = MomentumDirection_i.z();

			/*const double dx = xm * 0.5*mm / sqrt(xm * xm + zm * zm);//probably, incorrect
			const double dy = ym * 0.5*mm / sqrt(ym * ym + zm * zm);*/

			const double sign = zm > 0 ? 1 : -1;
			const double dx = z_size * xm / (zm*sign);//it should be correct
			const double dy = z_size * ym / (zm*sign);

			const G4ThreeVector& Position_f = G4ThreeVector(x + dx, y + dy, (z_pos + z_size*sign));

			distance = sqrt(pow(Position_f.x() - x_center, 2.0) + pow(Position_f.y() - y_center, 2.0));

			//cout << "center: " << x_center << "\t" << y_center << "\t" << endl;
			//cout << "dxy: " << dx << "\t" << dy << endl;
			//cout << "Position_f: " << Position_f.x() << "\t" << Position_f.y() << endl;
			//cout << "distance: " << distance << endl;
			if (distance < gPars::det_dims.THGEM1_hole_radius) {
				theStep->GetTrack()->SetPosition(Position_f);
				//theStep->GetTrack()->SetMomentumDirection(MomentumDirection_i); // I see strange result and don't know why
				//cout << "InHole !" << "x_center = " << x_center << " \t y_center = " << y_center << endl;
				//cout << "InHole !" << endl;
				//cout << "InHole !" << "Position_f.x() = " << Position_f.x() << " \t Position_f.y()  = " << Position_f.y() <<
					//" \t Position_f.z()  = " << Position_f.z() << endl;
			} else {
				theStep->GetTrack()->SetTrackStatus(fStopAndKill);
				//system("pause");
			}
		} else { //distance >= gPars::det_dims.radius_THGEM_hole
			//cout << "OutHole !" << "x_center = " << x_center << " \t y_center = " << y_center << endl;
			//cout << "distance: " << distance << endl;
			//#ifndef Cu_REFLECTION
			//theStep->GetTrack()->SetTrackStatus(fStopAndKill);
			//#endif //Cu_REFLECTION
		}
	}
}

/*
void SteppingAction::PassThroughGEM(const G4Step* theStep, G4double z_pos, G4double z_size)
{
	//change MomentumDirection and position for optical photon incidents on THGEM

	const G4ThreeVector& pos_i = theStep->GetPostStepPoint()->GetPosition();
	if (pos_i.z() > (z_pos - delta_z) && pos_i.z() < (z_pos + delta_z)) {
		//get
		const G4ThreeVector& MomentumDirection_i = theStep->GetTrack()->GetMomentumDirection();
		const G4ThreeVector& Position_i = theStep->GetTrack()->GetPosition();
		Position_init_THGEM = Position_i;

		double step = gPars::det_dims.step_THGEM_hole;
		const double step_y = step * sqrt(3) / 2;
		double y_center;
		double x_center;

		double x = Position_i.x();
		double y = Position_i.y();
		double x_abs = fabs(Position_i.x());
		double y_abs = fabs(Position_i.y());

		// the nearest y_center
		if (y_abs < step_y / 2) {
			y_center = 0;
		} else {
			y_center = ((int)((y_abs - step_y / 2) / step_y) + 1) * step_y * y / y_abs;
		}
		//the nearest x_center
		if (y_abs < step_y / 2) {
			int n_x = (int)(x / step);
			int sign = x > 0 ? 1 : -1;
			x_center = (n_x + sign * 0.5) * step;
		} else {
			if (((int)((y_abs - (step_y / 2)) / (step_y)) + 1) % 2 == 1) {  //if y row is odd
				if (x_abs < step / 2.0) {
					x_center = 0;
				} else {
					int n_x = (x_abs - step / 2) / step;
					int sign = x > 0 ? 1 : -1;
					x_center = (n_x + 1) * step * sign;
				}
			} else {
				int n_x = (int)(x / step);
				int sign = x > 0 ? 1 : -1;
				x_center = (n_x + sign * 0.5) * step;
			}
		}

		//kill or allow photon to pass
		double distance = sqrt(pow(x - x_center, 2.0) + pow(y - y_center, 2.0));
		//cout << "x - x_center = " << x - x_center << " \t y - y_center = " << y - y_center << endl;
		if (distance < gPars::det_dims.radius_THGEM_hole) {
			//set
			//const G4ThreeVector& MomentumDirection_f = G4ThreeVector(0, 1, 0);
			double xm = MomentumDirection_i.x();
			double ym = MomentumDirection_i.y();
			double zm = MomentumDirection_i.z();

			//set new initial position inside SingleTHGEMHole
			double zm_sign = zm > 0 ? 1 : -1;
			const G4ThreeVector& Position_new_i = G4ThreeVector(gPars::det_dims.xyz_position_SingleTHGEMHole + x - x_center,
				gPars::det_dims.xyz_position_SingleTHGEMHole + y - y_center, gPars::det_dims.xyz_position_SingleTHGEMHole - zm_sign*z_size/2.0);
			Position_init_SingleTHGEMHole = Position_new_i;
			if (distance < gPars::det_dims.radius_THGEM_hole) {
				MomentumDirection_init_THGEM = MomentumDirection_i;
				Volume_init_THGEM = theStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();
				//inside SingleTHGEMHole, start propagation
				is_photon_inside_single_THGEM_hole = true;
				theStep->GetTrack()->SetPosition(Position_new_i);
				//theStep->GetTrack()->SetMomentumDirection(MomentumDirection_i); // I see strange result and don't know why
				//cout << "InHole !" << "x_center = " << x_center << " \t y_center = " << y_center << endl;
				//cout << "x - x_center = " << x - x_center << "\t y - y_center = " << y - y_center << endl;
			} else {
				theStep->GetTrack()->SetTrackStatus(fStopAndKill);
				//system("pause");
			}
		} else {
			//cout << "OutHole !" << "x_center = " << x_center << " \t y_center = " << y_center << endl;
			//theStep->GetTrack()->SetTrackStatus(fStopAndKill);
		}
	}
}
*/
