#include <iostream> //std::cout
#include <fstream> //std::ofstream
#include <iomanip> //std::setprecision
#include <vector> //std::vector
#include <random>
#include <functional>

#include "selene/include/selene.h"

#include "Galaxy.hpp"
#include "Logger.hpp"
#include "ParseLua.hpp"
#include "FileManagement.hpp"

void showUsage();
void parseParameters(const std::string& filename, GalaxyParameters& p);
std::string parseOutputDirectory(const std::string& paramfilename);

int main(int argc, char** argv) {
	bool silent = false;
	std::string paramFile = "config/usr/galsim-config.lua";

	// Parse parameters
	if (argc > 3) {
		showUsage();
		return 0;
	}
	if (argc > 1) {
		for (int iarg = 1; iarg < argc; ++iarg) {
			std::string arg = argv[iarg];
			std::string prefix1("--config=");
			std::string prefix2("-s");

			if (!arg.compare(0, prefix1.size(), prefix1))
				paramFile = arg.substr(prefix1.size()).c_str();
			else if (!arg.compare(0, prefix2.size(), prefix2))
				silent = true;
			else {
				showUsage();
				return 0;
			}
		}
	}

	std::unique_ptr<LogPolicyInterface> consoleLogPolicy 
		= std::unique_ptr<ConsoleLogPolicy>(new ConsoleLogPolicy());
	consoleLogPolicy->setLogLevel(silent ? SeverityType::ERROR : SeverityType::DEBUG);
	Logger::Instance().registerLogPolicy("console", std::move(consoleLogPolicy));

	GalaxyParameters gal_params;

	try {
		gal_params.outputDirectory = parseOutputDirectory(paramFile);

		FileManagement::makeDirectoryPath(gal_params.outputDirectory);
		FileManagement::makeDirectoryPath(gal_params.outputDirectory + "/log");

		std::unique_ptr<LogPolicyInterface> fileLogPolicy 
			= std::unique_ptr<FileLogPolicy>(
				  new FileLogPolicy(gal_params.outputDirectory + "/log/galsim.log")
			  );
		fileLogPolicy->setLogLevel(SeverityType::NOTICE);
		Logger::Instance().registerLogPolicy("file", std::move(fileLogPolicy));
	}
	catch (std::exception& e) {
		Logger::Instance().print<SeverityType::FATAL_ERROR>(e.what());
		return 0;
	}

	try {
		parseParameters(paramFile, gal_params);
		Galaxy galaxy(gal_params);
	}
	catch (std::exception& e) {
		Logger::Instance().print<SeverityType::FATAL_ERROR>(e.what());
		std::cout << "See " << gal_params.outputDirectory << "/log/galsim.log for more details." << std::endl;
		return 0;
	}

	return 0;
}

std::string parseOutputDirectory(const std::string& paramfilename) {
	std::string fullfilename = paramfilename;
	std::string outputDir = "";
	// Create new Lua state and load the lua libraries
	sel::State luaState{true};

	if (!luaState.Load(fullfilename)) {
		throw std::runtime_error("ParseParameters: could not open lua file: " + fullfilename);
	}
	else {
		parseLuaVariable(luaState["Parameters"]["output_directory"], outputDir);
	}

	return outputDir;
}

void parseParameters(const std::string& filename, GalaxyParameters& p) {
	// Create new Lua state and load the lua libraries
	sel::State luaState{true};

	if (!luaState.Load(filename)) {
		throw std::runtime_error("ParseParameters: could not open lua file: " + filename);
	}
	else {
		std::ostringstream num[4];
		num[0] << 1;
		num[1] << 2;
		num[2] << 3;
		num[3] << 4;

		parseLuaVariable(luaState["Parameters"]["densityfile"], p.densityfile);
		parseLuaVariable(luaState["Parameters"]["starpopfile"], p.starpopfile);

		parseLuaVariable(luaState["Parameters"]["output_directory"], p.outputDirectory);
		const size_t len = p.outputDirectory.size();
		if ( p.outputDirectory[len-1] != '\\' && p.outputDirectory[len-1] != '/' )
			p.outputDirectory.push_back('/');

		parseLuaVariable(luaState["Parameters"]["sfr"], p.sfr);
		parseLuaVariable(luaState["Parameters"]["total_time"], p.total_time);

		parseLuaVariable(luaState["Parameters"]["radius"], p.radius);
		parseLuaVariable(luaState["Parameters"]["height"], p.height);
		parseLuaVariable(luaState["Parameters"]["resolution"], p.resolution);
		parseLuaVariable(luaState["Parameters"]["solar_position_x"], p.solar_position[0]);
		parseLuaVariable(luaState["Parameters"]["solar_position_y"], p.solar_position[1]);
		parseLuaVariable(luaState["Parameters"]["solar_position_z"], p.solar_position[2]);
		parseLuaVariable(luaState["Parameters"]["A_inner"], p.A_inner);
		parseLuaVariable(luaState["Parameters"]["A_a"], p.Aa);
		parseLuaVariable(luaState["Parameters"]["A_1"], p.A1);
		parseLuaVariable(luaState["Parameters"]["A_2"], p.A2);
		parseLuaVariable(luaState["Parameters"]["H_1"], p.H1);
		parseLuaVariable(luaState["Parameters"]["H_2"], p.H2);
		parseLuaVariable(luaState["Parameters"]["n_1"], p.n1);
		parseLuaVariable(luaState["Parameters"]["n_2"], p.n2);
		parseLuaVariable(luaState["Parameters"]["h_1"], p.h1);
		parseLuaVariable(luaState["Parameters"]["arm_min_radius"], p.rmin);
		for (int i = 0; i < 4; ++i) {
			parseLuaVariable(luaState["Parameters"][("wa_w"+num[i].str()).c_str()], p.wa_wj[i]);
			parseLuaVariable(luaState["Parameters"][("na_f"+num[i].str()).c_str()], p.na_fj[i]);
			parseLuaVariable(luaState["Parameters"][("ha_h"+num[i].str()).c_str()], p.ha_hj[i]);
			parseLuaVariable(luaState["Parameters"][("arm_joins_"+num[i].str()).c_str()], p.sp_joins[i]);
			parseLuaVariable(luaState["Parameters"][("logarm_rmin_"+num[i].str()).c_str()], p.sp_rmin[i]);
			parseLuaVariable(luaState["Parameters"][("logarm_thmin_"+num[i].str()).c_str()], p.sp_thmin[i]);
			parseLuaVariable(luaState["Parameters"][("logarm_asp_"+num[i].str()).c_str()], p.sp_asp[i]);
			parseLuaVariable(luaState["Parameters"][("logarm_rmax_"+num[i].str()).c_str()], p.sp_rmax[i]);
		}
	}
}

void showUsage() {
	std::cout << "galsim [-s] [--config=<filename>]" << std::endl;
}
