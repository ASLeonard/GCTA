#pragma once

#include <map>
#include <string>
#include <vector>

class PCAStream {
public:
    static std::map<std::string, std::string> options;
    static std::map<std::string, double> options_d;
    static std::vector<std::string> processFunctions;

    static int registerOption(std::map<std::string, std::vector<std::string>>& options_in);
    static void processMain();
};
