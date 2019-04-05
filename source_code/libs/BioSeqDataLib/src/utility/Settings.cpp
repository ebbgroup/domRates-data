#include "Settings.hpp"

namespace BioSeqDataLib
{

	void
	Settings::readSettings()
	{
		const char *domainWorldPath = std::getenv("DOMAINWORLD_DATA");
		fs::path settingsPath;
		if (domainWorldPath == nullptr)
		{
			settingsPath = std::getenv("HOME");
			settingsPath /= ".domainWorld";
		}
		else
			settingsPath = domainWorldPath;


		this->mapping["domainworld_path"] = settingsPath;
		this->mapping["dsm"] = settingsPath / "dsm";
		this->mapping["rads_db"] = settingsPath / "rads_db";

		settingsPath /= "general.conf";
		AP::Input settingsFile(settingsPath);
		std::string line;
		while(getline(settingsFile, line))
		{
			if (line.empty() || (line[0] == '#'))
				continue;
			auto tokens = BioSeqDataLib::split(line, "=");
			std::transform(tokens[0].begin(), tokens[0].end(), tokens[0].begin(), ::tolower);
			tokens[1].erase (std::remove(tokens[1].begin(), tokens[1].end(), '"'), tokens[1].end());
			this->mapping.emplace(std::make_pair(std::move(tokens[0]), std::move(tokens[1])));
		}
	}
}
