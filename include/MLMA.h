#pragma once
#include <map>
#include <string>
#include <vector>

class MLMA {
public:
    static int  registerOption(std::map<std::string, std::vector<std::string>>& options_in);
    static void processMain();

private:
    static std::map<std::string, std::string> options;
    static std::map<std::string, double>      options_d;
    static std::map<std::string, std::vector<double>> options_vd;
    static std::vector<std::string>           processFunctions;
};
