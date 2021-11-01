//freeman_arg_parse.h
//by Jared Freeman

//utility functions for extracting formatted arguments from passed-in arg vector


#ifndef FREEMAN_ARG_PARSE_H
#define FREEMAN_ARG_PARSE_H

#include <iostream>
#include <string>
#include <vector>

std::vector<std::string> ExtractArgs(std::string match_criterion, std::vector<std::string> arg_vec)
{
	std::vector<std::string> r_vec;

	/*
	for (std::vector<std::string>::iterator i = arg_vec.begin() ; i != arg_vec.end(); ++i)
	{
		std::cout << *i << "\n";
	}
	*/

	for(int j = 0; j < arg_vec.size(); j++)
	{
		if(match_criterion == arg_vec[j] && j < arg_vec.size() - 1)
		{
			//std::cout << "match: " << match_criterion << "\n";
			r_vec.push_back(arg_vec[j + 1]);
		}
	}

	return r_vec;
}

#endif
