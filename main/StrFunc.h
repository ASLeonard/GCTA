/*
 * Interface to the string-operation functions
 *
 * 2010 by Jian Yang <jian.yang.qt@gmail.com>
 * Modernized 2026 with C++20/23/26 features
 *
 * This file is distributed under the GNU General Public
 * License, Version 3.  Please see the file LICENSE for more
 * details
 */

#ifndef _STRFUNC_H
#define _STRFUNC_H

#include <string>
#include <string_view>
#include <vector>
#include <algorithm>
#include <ranges>
#include <map>
#include <iostream>
#include <cctype>

namespace StrFunc
{
	[[nodiscard]] bool str_within_quto(std::string_view str, std::string &str_buf);
	[[nodiscard]] int split_string(std::string_view str, std::vector<std::string> &vec_str, std::string_view separator=" ,\t;\n");
	[[nodiscard]] std::string first_string(std::string_view str, char separator);
	[[nodiscard]] std::string last_string(std::string_view str, char separator);
	void to_upper(std::string &str) noexcept;
	void to_lower(std::string &str) noexcept;
	[[nodiscard]] std::string get_sub_str(std::string_view rst, int pos);
	[[nodiscard]] bool StrEqual(std::string_view StrA, std::string_view StrB, bool NoCaseSens=true) noexcept;
	[[nodiscard]] bool StrVecEqual(const std::vector<std::string> &VsBufA, const std::vector<std::string> &VsBufB, int Pos);

	// find a string in a string vector ignoring upper or lower case
	[[nodiscard]] std::vector<std::string>::iterator find(std::vector<std::string> &target_vs, std::string_view target_str);

	// find a char in a string ignoring upper or lower case
	[[nodiscard]] std::string::iterator find(std::string &target_str, char target_ch);

	// go to the position of a given string in a stream ignoring upper or lower case
	[[nodiscard]] bool goto_str(std::istream &in_file, std::string_view str);

	// rewind a stream
	void rewind_if(std::istream &in_file) noexcept;

	// match two vectors
	void match(const std::vector<std::string> &VecA, const std::vector<std::string> &VecB, std::vector<int> &VecC);
	
	// compare two strings ignoring case (modern C++20 constexpr implementation)
	[[nodiscard]] constexpr bool i_compare(std::string_view a, std::string_view b) noexcept
	{
		if (a.length() != b.length()) {
			return false;
		}
		
		return std::ranges::equal(a, b, [](char c1, char c2) {
			return std::toupper(static_cast<unsigned char>(c1)) == 
			       std::toupper(static_cast<unsigned char>(c2));
		});
	}
}

#endif
