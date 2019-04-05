/*
 * Settings.hpp
 *
 *  Created on: October 26, 2016
 *      Author: Carsten Kemena
 *	 Copyright: 2016
 *
 *  This file is part of BioSeqDataLib.
 *
 *  BioSeqDataLib is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  BioSeqDataLib is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with BioSeqDataLib.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file Settings.hpp
 * \brief This file contains the settings class for Domain World.
 */

#ifndef SETTINGS_HPP_
#define SETTINGS_HPP_

#include <algorithm>
#include <string>
#include <map>

#include <boost/filesystem.hpp>

#include "../external/Input.hpp"
#include "../utility/stringHelpers.hpp"

namespace fs=boost::filesystem;
namespace AP=AlgorithmPack;

namespace BioSeqDataLib
{

	/**
	 * \brief The settings class is a small class to allow simply read the settings for the configuration file.
	 *
	 */
	class Settings
	{
	private:
		std::map<std::string, fs::path> mapping;

	public:

		Settings()
		{}

		Settings(const Settings&) = default;
		Settings(Settings &&) = default;
		virtual ~Settings() = default;

		/**
		 * \brief Return the settings for a specific path
		 * \param name The name of the configuration that is needed.
		 */
		const fs::path
		&operator[](const std::string &name) const
		{
			return mapping.at(name);
		}

		/*
		fs::path
		&operator[](const std::string &name)
		{
			return mapping[name];
		}*/

		/**
		 * \brief Reads the settings from the system.
		 */
		void
		readSettings();


	};
}

#endif /* SETTINGS_HPP_ */
