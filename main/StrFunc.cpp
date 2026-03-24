/*
 * Implementations of the string-operation functions
 *
 * 2010 by Jian Yang <jian.yang.qt@gmail.com>
 * Modernized 2026 with C++20/23/26 features
 *
 * This file is distributed under the GNU General Public
 * License, Version 3.  Please see the file LICENSE for more
 * details
 */

#include "StrFunc.h"
#include <cctype>
#include <ranges>
#include "Logger.h"

int StrFunc::split_string(std::string_view str, std::vector<std::string> &vec_str, std::string_view separator)
{
	if(str.empty()) return 0;
	vec_str.clear();

	// Build set of characters to keep (not separators)
	constexpr std::string_view all_symbols = 
		"`1234567890-=~!@#$%^&*()_+qwertyuiop[]\\asdfghjkl;'zxcvbnm,./"
		"QWERTYUIOP{}|ASDFGHJKL:\"ZXCVBNM<>? \t\n";
	
	std::string symbol_pool{all_symbols};
	
	// Remove separator characters from the pool
	for (const char sep : separator) {
		if (auto pos = symbol_pool.find(sep); pos != std::string::npos) {
			symbol_pool.erase(pos, 1);
		}
	}

	// Extract tokens
	std::string str_buf;
	str_buf.reserve(str.length());  // Optimize allocation
	bool in_token = false;

	for (const char ch : str) {
		if (symbol_pool.find(ch) != std::string::npos) {
			in_token = true;
			str_buf += ch;
		} else {
			if (in_token) {
				vec_str.push_back(std::move(str_buf));
				str_buf.clear();
				in_token = false;
			}
		}
	}
	
	if (in_token) {
		vec_str.push_back(std::move(str_buf));
	}

	return static_cast<int>(vec_str.size());
}

std::string StrFunc::first_string(std::string_view str, char separator)
{
	if (auto pos = str.find(separator); pos != std::string_view::npos) {
		return std::string{str.substr(0, pos)};
	}
	return {};
}

std::string StrFunc::last_string(std::string_view str, char separator)
{
	if (auto pos = str.find_last_of(separator); pos != std::string_view::npos) {
		return std::string{str.substr(pos + 1)};
	}
	return {};
}

void StrFunc::to_upper(std::string &str) noexcept
{
	std::ranges::transform(str, str.begin(), 
		[](unsigned char c) { return std::toupper(c); });
}

void StrFunc::to_lower(std::string &str) noexcept
{
	std::ranges::transform(str, str.begin(), 
		[](unsigned char c) { return std::tolower(c); });
}

std::string StrFunc::get_sub_str(std::string_view rst, int pos)
{
	std::vector<std::string> vs_buf;
	if (split_string(rst, vs_buf) > pos) {
		return vs_buf[pos];
	}
	return {};
}

bool StrFunc::StrEqual(std::string_view StrA, std::string_view StrB, bool NoCaseSens) noexcept
{
	if (!NoCaseSens) {
		return StrA == StrB;
	}
	return i_compare(StrA, StrB);
}

bool StrFunc::StrVecEqual(const std::vector<std::string> &VsBufA, const std::vector<std::string> &VsBufB, int Pos)
{
	if (VsBufA.size() != VsBufB.size()) {
		return false;
	}
	if (Pos >= static_cast<int>(VsBufA.size())) {
		LOGGER.e(0, "invalid Pos. StrFunc::StrVecEqual");
	}
	return std::ranges::equal(
		VsBufA | std::views::drop(Pos),
		VsBufB | std::views::drop(Pos)
	);
}

bool StrFunc::str_within_quto(std::string_view str, std::string &str_buf)
{
	auto begin = str.find_first_of('"');
	auto end = str.find_last_of('"');
	
	if (begin == std::string_view::npos || end == std::string_view::npos || begin == end) {
		return false;
	}

	str_buf = str.substr(begin + 1, end - begin - 1);
	return true;
}

std::vector<std::string>::iterator StrFunc::find(std::vector<std::string> &target_vs, std::string_view target_str)
{
	auto it = std::ranges::find_if(target_vs, [target_str](const std::string& s) {
		return i_compare(s, target_str);
	});
	return it;
}

std::string::iterator StrFunc::find(std::string &target_str, char target_ch)
{
	auto it = std::ranges::find_if(target_str, [target_ch](char c) {
		return std::toupper(static_cast<unsigned char>(c)) == 
		       std::toupper(static_cast<unsigned char>(target_ch));
	});
	return it;
}

bool StrFunc::goto_str(std::istream &in_file, std::string_view str)
{
	std::string str_buf;
	std::string query_str{str};
	std::vector<std::string> vs_buf;
	to_upper(query_str);
	
	while (in_file >> str_buf) {
		if (split_string(str_buf, vs_buf) > 0) {
			str_buf = vs_buf[0];
		} else {
			continue;
		}
		
		to_upper(str_buf);
		if (str_buf == "#") {
			getline(in_file, str_buf);
			continue;
		}
		if (str_buf == query_str) {
			return true;
		}
	}

	return false;
}

void StrFunc::rewind_if(std::istream &in_file) noexcept
{
	in_file.clear(std::ios::goodbit);
	in_file.seekg(std::ios::beg);
}

void StrFunc::match(const std::vector<std::string> &VecA, const std::vector<std::string> &VecB, std::vector<int> &VecC)
{
	// Build lookup map from VecB
	std::map<std::string, int> id_map;
	for (int i = 0; i < static_cast<int>(VecB.size()); ++i) {
		id_map.emplace(VecB[i], i);
	}
	
	// Match elements from VecA
	VecC.clear();
	VecC.reserve(VecA.size());
	for (const auto& item : VecA) {
		if (auto it = id_map.find(item); it != id_map.end()) {
			VecC.push_back(it->second);
		} else {
			VecC.push_back(-9);
		}
	}
}
