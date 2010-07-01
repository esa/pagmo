#include<fstream>
#include"boost/algorithm/string.hpp"

#include"planet_mpcorb.h"
#include"exceptions.h"

static const int mpcorb_format[8][2] =
{
	{92,11},	// a (AU)
	{70,9},		// e
	{59,9},		// i (deg)
	{48,9},		// Omega (deg)
	{37,9},		// omega (deg)
	{26,9},		// M (deg)
	{20,5},		// Epoch (shitty format)
	{166,28}	// Asteroid readable name
};

namespace kep_toolbox{

planet_mpcorb::planet_mpcorb(const std::string& name)
{
	array6D elem;
	std::string lower_case_name = name;
	boost::algorithm::to_lower(lower_case_name);

	std::ifstream mpcorbfile("MPCORB.DAT");
	if (!mpcorbfile) {
		throw_value_error("Could not find MPCORB.DAT is it in the current directory?");
	}
	else {
		std::string tmp;
		std::string line;

		//skipping the first lines of the file
		do {
			std::getline(mpcorbfile,line);
		} while (!boost::algorithm::find_first(line,"-----------------"));

		//reading through the whole file
		while (!mpcorbfile.eof()) {
			//read line
			std::getline(mpcorbfile,line);
			//to lowercase
			boost::algorithm::to_lower(line);
			//if the asteroid name is there
			if(boost::algorithm::find_first(line,lower_case_name)) {
				//read keplerian elements from MPCORB.DAT
				for (int i = 0; i < 6; ++i) {
					tmp.clear();
					tmp.append(&line[mpcorb_format[i][0]],mpcorb_format[i][1]);
					boost::algorithm::trim(tmp);
					elem[i] = boost::lexical_cast<double>(tmp);
				}
				// Converting orbital elements to the dictatorial PaGMO units.
				elem[0] *= ASTRO_AU;
				for (int i = 2; i < 6; ++i) {
					elem[i] *= ASTRO_DEG2RAD;
				}
				// Deal with MPCORB data format
				tmp.clear();
				tmp.append(&line[mpcorb_format[6][0]],mpcorb_format[6][1]);
				boost::algorithm::trim(tmp);
				boost::gregorian::greg_year anno = packed_date2number(tmp[0]) * 100 + boost::lexical_cast<int>(std::string(&tmp[1],&tmp[3]));
				boost::gregorian::greg_month mese = packed_date2number(tmp[3]);
				boost::gregorian::greg_day giorno = packed_date2number(tmp[4]);
				epoch epoch(anno,mese,giorno);
				// Record asteroid name.
				tmp.clear();
				tmp.append(&line[mpcorb_format[7][0]],mpcorb_format[7][1]);
				boost::algorithm::trim(tmp);
				build_planet(epoch,elem,ASTRO_MU_SUN,100,100,100,tmp);
				mpcorbfile.close();
				break;
			}
		}
		if (mpcorbfile.eof()) {
			mpcorbfile.close();
			throw_value_error("Could not find minor planet in MPCORB.DAT");
		}
	}
}

planet_mpcorb::planet_mpcorb(int row)
{
	array6D elem;

	std::ifstream mpcorbfile("MPCORB.DAT");
	if (!mpcorbfile) {
		throw_value_error("Could not find MPCORB.DAT is it in the current directory?");
	}
	else {
		std::string tmp;
		std::string line;

		//skipping the first lines of the file
		do {
			std::getline(mpcorbfile,line);
		} while (!boost::algorithm::find_first(line,"-----------------"));

		//skipping to row
		int i=0;
		do {
			//read line
			std::getline(mpcorbfile,line);
			++i;
		} while (!mpcorbfile.eof() && !((i-1)==row));

		if (mpcorbfile.eof()) {
			mpcorbfile.close();
			throw_value_error("Could not find row in MPCORB.DAT");
		}
		//to lowercase
		boost::algorithm::to_lower(line);

		//read keplerian elements from MPCORB.DAT
		for (int i = 0; i < 6; ++i) {
			tmp.clear();
			tmp.append(&line[mpcorb_format[i][0]],mpcorb_format[i][1]);
			boost::algorithm::trim(tmp);
			try {
				elem[i] = boost::lexical_cast<double>(tmp);
			} catch (boost::bad_lexical_cast) {
				std::cout << "Could not construct planet_mpcorb("<< row << "). Empty line?\n";
				throw boost::bad_lexical_cast();
			}
		}
		// Converting orbital elements to the dictatorial PaGMO units.
		elem[0] *= ASTRO_AU;
		for (int i = 2; i < 6; ++i) {
			elem[i] *= ASTRO_DEG2RAD;
		}
		// Deal with MPCORB data format
		tmp.clear();
		tmp.append(&line[mpcorb_format[6][0]],mpcorb_format[6][1]);
		boost::algorithm::trim(tmp);
		boost::gregorian::greg_year anno = packed_date2number(tmp[0]) * 100 + boost::lexical_cast<int>(std::string(&tmp[1],&tmp[3]));
		boost::gregorian::greg_month mese = packed_date2number(tmp[3]);
		boost::gregorian::greg_day giorno = packed_date2number(tmp[4]);
		epoch epoch(anno,mese,giorno);
		// Record asteroid name.
		tmp.clear();
		tmp.append(&line[mpcorb_format[7][0]],mpcorb_format[7][1]);
		boost::algorithm::trim(tmp);
		build_planet(epoch,elem,ASTRO_MU_SUN,100,100,100,tmp);
		mpcorbfile.close();
	}

}

// Convert mpcorb packed dates conventio into number. (lower case assumed)
// TODO: check locale ASCII.
int planet_mpcorb::packed_date2number(char c)
{
	return static_cast<int>(c) - (boost::algorithm::is_alpha()(c) ? 87 : 48);
}

} //namespace
