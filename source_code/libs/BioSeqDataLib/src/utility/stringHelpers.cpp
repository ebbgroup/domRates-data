#include "stringHelpers.hpp"

namespace BioSeqDataLib
{


std::vector<std::string>
split(const std::string &s, const std::string &delimiters, bool keepEmpty)
{
	size_t start=0;
	size_t end=s.find_first_of(delimiters);
	std::vector<std::string> output;
	size_t len=s.size();
	while (start != len)
	{
		if (keepEmpty || (end-start >0))
			output.emplace_back(s.substr(start, end-start));

		if (end == std::string::npos)
			break;
		start=end+1;
		end = s.find_first_of(delimiters, start);
	}

	return output;
}

void
split(const std::string &s, const std::string &delimiters, std::vector<std::string> &tokens, bool keepEmpty)
{
	size_t start=0;
	size_t end=s.find_first_of(delimiters);
	size_t len=s.size();
	while (start != len)
	{
		if (keepEmpty || (end-start >0))
			tokens.emplace_back(s.substr(start, end-start));

		if (end == std::string::npos)
			break;

		start=end+1;
		end = s.find_first_of(delimiters, start);
	}
}


void
trimRight(std::string &str)
{
	size_t i= str.size();
	while ( (--i != 0) && (std::isspace(str[i])) )
		;
	str.resize(i+1);
}

void
removeSpaces(std::string &str, int startPos)
{
	size_t pos=startPos;
	size_t length = str.size();
	for (size_t i=startPos; i<length; ++i)
	{
		if (!std::isspace(str[i]))
			str[pos++] = str[i];
	}
	str.resize(pos);
}

void
strip(std::string &line, size_t startPos)
{
	line.erase(std::remove_if(line.begin()+startPos, line.end(), [] (std::string::value_type ch){ return !isalpha(ch); }), line.end());
}


}
