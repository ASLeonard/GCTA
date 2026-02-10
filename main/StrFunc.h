/*
 * Interface to the string-operation functions
 *
 * 2010 by Jian Yang <jian.yang.qt@gmail.com>
 *
 * This file is distributed under the GNU General Public
 * License, Version 3.  Please see the file LICENSE for more
 * details
 */

#ifndef _STRFUNC_H
#define _STRFUNC_H

#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <iostream>

namespace StrFunc
{
	bool str_within_quto(const std::string &str, std::string &str_buf);
	int split_string(const std::string &str, std::vector<std::string> &vec_str, std::string separator=" ,\t;\n");
	std::string first_string(const std::string &str, const char separator);
	std::string last_string(const std::string &str, const char separator);
	void to_upper(std::string &str);
	void to_lower(std::string &str);
	std::string get_sub_str(const std::string & rst, int pos);
	bool StrEqual(const std::string &StrA, const std::string &StrB, bool NoCaseSens=true);
	bool StrVecEqual(const std::vector<std::string> &VsBufA, const std::vector<std::string> &VsBufB, int Pos);

	// find a string in a string vector ignoring upper or lower case
	std::vector<std::string>::iterator find(std::vector<std::string> &target_vs, const std::string &target_str);

	// find a char in a string ignoring upper or lower case
	std::string::iterator find(std::string &target_str, const char target_ch);

	// go to the postion of a give string in a stream ignoring upper or lower case
	bool goto_str(std::istream &in_file, const std::string &str);

	// rewind a stream
	void rewind_if(std::istream &in_file);

	// match two vectors
	void match(const std::vector<std::string> &VecA, const std::vector<std::string> &VecB, std::vector<int> &VecC);
        // compare two string ignore the case
        bool i_compare(std::string const& a, std::string const& b);
}

#endif
