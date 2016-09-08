#include "HosTracks.hpp"

HosData::HosData() {}

HosTrack::HosTrack()
		: mdot(0)
{

}

HosTrack::HosTrack(double md, const std::vector<double>& m, const std::vector<double>& r, const std::vector<double>& rph, const std::vector<double>& l, const std::vector<double>& lt, const std::vector<double>& t)
		: mdot(md)
{
	data.mass = m;
	data.rstar = r;
	data.rphot = rph;
	data.lstar = l;
	data.ltot = lt;
	data.time = t;
}

HosTracks::HosTracks(std::string filename)
{
	std::vector<double> mdot, m, r, rph, l, lt, t, dum;
	std::vector<int> id;
	std::string delimiter = "\t";
	readcols(filename, delimiter, mdot, id, m, r, rph, l, lt, t, dum);
	int prev = 0;
	int ntracks = mdot.size();
	for (int i = 0; i < ntracks; i++) {
		if (id[i] == 0 && i != 0) {
			std::vector<double> mpart(m.begin() + prev, m.begin()+i);
			std::vector<double> rpart(r.begin() + prev, r.begin()+i);
			std::vector<double> rphpart(rph.begin() + prev, rph.begin()+i);
			std::vector<double> lpart(l.begin() + prev, l.begin()+i);
			std::vector<double> ltpart(lt.begin() + prev, lt.begin()+i);
			std::vector<double> tpart(t.begin() + prev, t.begin()+i);
			HosTrack a_track(mdot[prev], mpart, rpart, rphpart, lpart, ltpart, tpart);
			tracks.push_back(a_track);
			prev = i;
		}
	}
	std::vector<double> mpart(m.begin() + prev, m.end());
	std::vector<double> rpart(r.begin() + prev, r.end());
	std::vector<double> rphpart(rph.begin() + prev, rph.end());
	std::vector<double> lpart(l.begin() + prev, l.end());
	std::vector<double> ltpart(lt.begin() + prev, lt.end());
	std::vector<double> tpart(t.begin() + prev, t.end());
	HosTrack a_track(mdot[prev], mpart, rpart, rphpart, lpart, ltpart, tpart);
	tracks.push_back(a_track);
}
