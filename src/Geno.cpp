/*
   GCTA: a tool for Genome-wide Complex Trait Analysis

   New implementation: read and process genotype of plink format in block way.

   Depends on the class of marker and phenotype

   Developed by Zhili Zheng<zhilizheng@outlook.com>

   This file is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   A copy of the GNU General Public License is attached along with this program.
   If not, see <http://www.gnu.org/licenses/>.
*/

#define NOMINMAX
#include "Geno.h"
#include "constants.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <thread>
#include <chrono>
#include <ctime>
#include <iostream>
#include <iterator>
#include <cmath>
#include <sstream>
#include <iomanip>
#include "utils.hpp"
#include "omp.h"
#include "GenoBackendFactory.h"
#include <cstring>
#include <boost/algorithm/string.hpp>
#include "OptionIO.h"
#include "zlib.h"
#include "zstd.h"
#include <cstring>
#include "cpu.h"
#include <Eigen/Eigen>
#include <algorithm>
#include "third_party/Pgenlib/PgenReader.h"
#include <numeric>

#ifdef _WIN64
  #include <intrin.h>
  uint32_t __inline CTZ64U(uint64_t value){
      unsigned long tz = 0;
      _BitScanForward64(&tz, value);
      return tz;
  }
  
  uint32_t __inline CLZ64U(uint64_t value){
      unsigned long lz = 0;
      _BitScanReverse64(&lz, value);
      return 63 - lz;
  }
#else
  //#define CTZU __builtin_ctz
  //#define CLZU __builtin_clz
  #if defined(__linux__) && GCTA_CPU_x86
  __attribute__((target("default")))
  #endif
  uint32_t CTZ64U(uint64_t value){
      return __builtin_ctzll(value);
  }
  #if defined(__linux__) && GCTA_CPU_x86
  __attribute__((target("popcnt")))
  uint32_t CTZ64U(uint64_t value){
      return __builtin_ctzll(value);
  }
  #endif
 
#endif

#if defined(__linux__) && GCTA_CPU_x86
__attribute__((target("default")))
#endif
uint64_t fill_inter_zero(uint64_t x) {
   uint64_t t;
   t = (x ^ (x >> 16)) & 0x00000000FFFF0000;
   x = x ^ t ^ (t << 16);
   t = (x ^ (x >> 8)) & 0x0000FF000000FF00;
   x = x ^ t ^ (t << 8);
   t = (x ^ (x >> 4)) & 0x00F000F000F000F0;
   x = x ^ t ^ (t << 4);
   t = (x ^ (x >> 2)) & 0x0C0C0C0C0C0C0C0C;
   x = x ^ t ^ (t << 2);
   t = (x ^ (x >> 1)) & 0x2222222222222222;
   x = x ^ t ^ (t << 1);
   return x;
}
#if defined(__linux__) && GCTA_CPU_x86
#include <x86intrin.h>
__attribute__((target("bmi2")))
uint64_t fill_inter_zero(uint64_t x) {
    return _pdep_u64(x, 0x5555555555555555U);
}
#endif


typedef uint32_t halfword_t;
const uintptr_t k1LU = (uintptr_t)1;

using std::thread;
using std::to_string;

map<string, string> Geno::options;
map<string, double> Geno::options_d;
vector<string> Geno::processFunctions;

Geno::Geno(Pheno* pheno, Marker* marker) {
    geno_files.clear();
    int num_geno = 0;
    if(options.find("geno_file") != options.end()){
        genoFormat = "BED";
        geno_files.push_back(options["geno_file"]);
        num_geno++;
        hasInfo = false;
    }

    if(options.find("m_file") != options.end()){
        genoFormat = "BED";
        boost::split(geno_files, options["m_file"], boost::is_any_of("\t "));
        //std::transform(geno_files.begin(), geno_files.end(), geno_files.begin(), [](string r){return r + ".bed";});
        num_geno++;
        hasInfo = false;
    }

    if(options.find("pgen_file") != options.end()){
        genoFormat = "PGEN";
        geno_files.push_back(options["pgen_file"]);
        boost::split(geno_files, options["pgen_file"], boost::is_any_of("\t "));
        num_geno++;
        hasInfo = false;
    }

    if(options.find("mpgen_file") != options.end()){
        genoFormat = "PGEN";
        boost::split(geno_files, options["mpgen_file"], boost::is_any_of("\t "));
        num_geno++;
        hasInfo = false;
    }

    if(options.find("bgen_file") != options.end()){
        genoFormat = "BGEN";
        geno_files.push_back(options["bgen_file"]);
        num_geno++;
        hasInfo = true;
    }

    if(options.find("mbgen_file") != options.end()){
        genoFormat = "BGEN";
        boost::split(geno_files, options["mbgen_file"], boost::is_any_of("\t "));
        num_geno++;
        hasInfo = true;
    }

    if(num_geno == 0){
        LOGGER.e(0, "no genotype file is specified");
    }

    //open and check genotype files

    this->pheno = pheno;
    this->marker = marker;

    //this->sampleKeepIndex = pheno->get_index_keep();

    //register format handlers



    // register getGenoDouble dispatch (fix PGEN bug: was aliased to BED)
    getGenoDoubleFuncs["BED"]  = &Geno::getGenoDouble_bed;
    getGenoDoubleFuncs["PGEN"] = &Geno::getGenoDouble_pgen;
    getGenoDoubleFuncs["BGEN"] = &Geno::getGenoDouble_bgen;
    // preGenoDouble*/endGenoDouble*/readGeno* maps removed:
    // those are now handled by BedBackend / BgenBackend / PgenBackend.

    string alleleFileName = "";
    if(options.find("update_freq_file") != options.end()){
        alleleFileName = options["update_freq_file"];
    }

    setMAF(options_d["min_maf"]);
    setMaxMAF(options_d["max_maf"]);
    setFilterInfo(options_d["info_score"]);
    setFilterMiss(1.0 - options_d["geno_rate"]);

    string filterprompt = "Threshold to filter variants:";
    bool outFilterPrompt = false;
    if(options_d["min_maf"] != 0.0){
        filterprompt += " MAF > " + to_string(options_d["min_maf"]);
        outFilterPrompt = true;
    }
    if(options_d["max_maf"] != 0.5){
        filterprompt += string(outFilterPrompt ? "," : "") + " MAF < " + to_string(options_d["max_maf"]);
        outFilterPrompt = true;
    }
    if(options_d["info_score"] != 0.0){
        filterprompt += string(outFilterPrompt ? "," : "") + " imputation INFO score > " + to_string(options_d["info_score"]);
        outFilterPrompt = true;
    }

    if(options_d["geno_rate"] != 1.0){
        filterprompt += string(outFilterPrompt ? "," : "") + " missingness rate < " + to_string(options_d["geno_rate"]);
        outFilterPrompt = true;
    }
    if(outFilterPrompt){
        LOGGER << filterprompt << "." << std::endl;
    }
    if(options_d["dos_dc"] == 1.0){
        iGRMdc = 1;
        iDC = 1;
        LOGGER << "Switch to the full dosage compensation mode." << std::endl;
    }else if(options_d["dos_dc"] == 0.0){
        iGRMdc = 0;
        iDC = 0;
        LOGGER << "Switch to the no dosage compensation mode." << std::endl;
    }else{
        // default equal variance mode for grm
        iGRMdc = -1;
        // compensation mode for x
        iDC = 1;
    }

    init_AF(alleleFileName);

    //olds
    //init_AsyncBuffer();

    //num_keep_sample = 0;
    //init_keep();
    //olds;
}

Geno::~Geno(){
}

uint32_t Geno::getTotalMarker(){
    return total_markers;
}

void Geno::setSexMode(){
    std::map<string, vector<string>> t_option;
    t_option["--chrx"] = {};
    t_option["--filter-sex"] = {}; 
    Pheno::registerOption(t_option);
    Marker::registerOption(t_option);
    Geno::registerOption(t_option);
}


//
//true:  filtered; false: not necessary to filter
bool Geno::filterMAF(){
    if((options_d["min_maf"] != 0.0) || (options_d["max_maf"] != 0.5)){
        LOGGER.i(0, "Computing allele frequencies...");

        int N = static_cast<int>(marker->count_extract());
        AFA1.assign(N, 0.0);
        countMarkers.assign(N, 0);
        vector<uint32_t> extractIdx(N);
        std::iota(extractIdx.begin(), extractIdx.end(), 0);

        // Temporarily open all filters so getGenoDouble computes AF for every
        // marker regardless of MAF/missingness thresholds.
        double sv_min = min_maf, sv_max = max_maf, sv_miss = dFilterMiss;
        min_maf = 0.0; max_maf = 1.0; dFilterMiss = 0.0;

        loopDouble(extractIdx, Constants::NUM_MARKER_READ, false, false, false, false,
            {[this](uintptr_t *buf, const vector<uint32_t> &exIdx) {
                int n = static_cast<int>(exIdx.size());
                #pragma omp parallel for schedule(static)
                for(int i = 0; i < n; ++i){
                    GenoBufItem item;
                    item.extractedMarkerIndex = exIdx[i];
                    getGenoDouble(buf, i, &item);
                    if(item.valid){
                        AFA1[exIdx[i]]         = item.af;
                        countMarkers[exIdx[i]] = item.nValidAllele;
                    }
                }
            }}, false);

        min_maf = sv_min; max_maf = sv_max; dFilterMiss = sv_miss;

        // Apply MAF filter
        double min_maf_eps = options_d["min_maf"] * (1.0 - Constants::SMALL_EPSILON);
        double max_maf_eps = options_d["max_maf"] * (1.0 + Constants::SMALL_EPSILON);
        LOGGER.d(0, "min_maf: " + to_string(min_maf_eps) + " max_maf: " + to_string(max_maf_eps));
        vector<uint32_t> extract_index;
        double cur_AF;

        for(int index = 0; index != static_cast<int>(AFA1.size()); index++){
            cur_AF = AFA1[index];
            if(cur_AF > 0.5) cur_AF = 1.0 - cur_AF;
            if((cur_AF > min_maf_eps) && (cur_AF < max_maf_eps)){
                extract_index.push_back(index);
                LOGGER.d(0, to_string(index) + ": " + to_string(cur_AF));
            }
        }

        vector<double> AFA1o = AFA1;
        vector<uint32_t> countMarkerso = countMarkers;

        AFA1.resize(extract_index.size());
        countMarkers.resize(extract_index.size());

        #pragma omp parallel for
        for(uint32_t index = 0; index < extract_index.size(); index++){
            uint32_t cur_index = extract_index[index];
            AFA1[index]        = AFA1o[cur_index];
            countMarkers[index] = countMarkerso[cur_index];
        }

        marker->keep_extracted_index(extract_index);

        LOGGER.i(0, to_string(extract_index.size()) + " SNPs remain from --maf or --max-maf,  ");
        return true;
    }else{
        return false;
    }
}

void Geno::init_AF(string alleleFileName) {
    AFA1.clear();
    //countA1A2.clear();
    //countA1A1.clear();
    //countA2A2.clear();
    countMarkers.clear();
    //RDev.clear();
    if(!alleleFileName.empty()){
        LOGGER.i(0, "Reading frequencies from [" + alleleFileName + "]...");
        vector<int> field_return = {2};
        vector<string> fields;
        vector<bool> a_rev;
        marker->matchSNPListFile(alleleFileName, 3, field_return, fields, a_rev, false);

        AFA1.resize(a_rev.size());
        vector<uint32_t> extract_index;
        bool filterByMaf = false;
        for(int i = 0; i < a_rev.size(); i++){
            double af;
            try{
                af = stod(fields[i]);
            }catch(std::out_of_range &){
                LOGGER.e(0, "the third column should be numeric");
            }
            if(af < 0 || af > 1.0){
                LOGGER.e(0, "frequency values should range from 0 to 1");
            }
            double maf = std::min(af, 1.0 - af);
            if(maf > min_maf && maf < max_maf){
                extract_index.push_back(i);
            }else{
                filterByMaf = true;
            }

            if(a_rev[i]){
                AFA1[i] = 1.0 - af;
            }else{
                AFA1[i] = af;
            }
        }
        LOGGER.i(0, "Frequencies of " + to_string(AFA1.size()) + " SNPs are updated.");

        marker->keep_extracted_index(extract_index);
        vector<double> AFA1o = AFA1;
        //vector<uint32_t> countMarkerso = countMarkers;

        AFA1.resize(extract_index.size());
        countMarkers.resize(extract_index.size());

        #pragma omp parallel for
        for(uint32_t index = 0; index < extract_index.size(); index++){
            uint32_t cur_index = extract_index[index];
            AFA1[index] = AFA1o[cur_index];
           // countMarkers[index] = countMarkerso[cur_index];
        }
        if(filterByMaf)LOGGER << "  " << extract_index.size() << " SNPs remain after MAF filtering." << std::endl;
        bHasPreAF = true;
    }
}

void Geno::out_freq(string filename){
    string name_frq = filename + ".frq";
    LOGGER.i(0, "Saving allele frequencies...");
    std::ofstream o_freq(name_frq.c_str());
    if (!o_freq) { LOGGER.e(0, "cannot open the file [" + name_frq + "] to write"); }
    vector<string> out_contents;
    out_contents.reserve(AFA1.size() + 1);
    out_contents.push_back("CHR\tSNP\tPOS\tA1\tA2\tAF\tNCHROBS");
    for(int i = 0; i != AFA1.size(); i++){
        out_contents.push_back(marker->get_marker(marker->getRawIndex(i)) + "\t" + to_string(AFA1[i])
                               + "\t" + to_string(countMarkers[i]));
    }
    std::copy(out_contents.begin(), out_contents.end(), std::ostream_iterator<string>(o_freq, "\n"));
    o_freq.close();
    LOGGER.i(0, "Allele frequencies of " + to_string(AFA1.size()) + " SNPs have been saved in the file [" + name_frq + "]");
}

bool Geno::getGenoHasInfo(){
    return hasInfo;
}

void Geno::setGRMMode(bool grm, bool dominace){
    this->bGRM = grm;
    this->bGRMDom = dominace;
}

void Geno::setGenoItemSize(uint32_t &genoSize, uint32_t &missSize){
    genoSize = keepSampleCT;
    missSize = missPtrSize;
}

void Geno::getGenoDouble(uintptr_t *buf, int bufIndex, GenoBufItem* gbuf){
    (this->*getGenoDoubleFuncs[genoFormat])(buf, bufIndex, gbuf);
}

void Geno::setMaleWeight(double &weight, bool &needWeight){
    weight = 1.0;
    if(bGRM){
        weight = sqrt(0.5);
        if(iGRMdc == 1){
            weight *= sqrt(2.0);
        }else if(iGRMdc == 0){
            weight *= sqrt(0.5);
        }
    }else{
        if(iDC == 0){
            weight = 0.5;
        }
    }
    if(std::abs(weight - 1.0) > 1e-6){
        needWeight = true;
    }else{
        needWeight = false;
    }
}

void Geno::getGenoDouble_pgen(uintptr_t *buf, int idx, GenoBufItem* gbuf){
    SNPInfo snpinfo;
    uintptr_t *cur_buf = buf + idx * pgenGenoBuf1PtrSize;
    uintptr_t *geno_buf = cur_buf;
    uintptr_t *dosage_present = cur_buf + pgenGenoPtrSize;
    uint16_t *dosage_main = reinterpret_cast<uint16_t*>(cur_buf + pgenGenoPtrSize + pgenDosagePresentPtrSize);
    uint32_t dosage_ct = cur_buf[pgenGenoBuf1PtrSize - 1];
    uint8_t isSexXY = isMarkersSexXYs[curBufferIndex];

    const vector<uint32_t> *curMaleIndex = NULL;
    if(isSexXY == 1){
        curMaleIndex = &keepMaleExtractIndex;
    }

    string err;
    if(!PgenReader::CountHardDosage(cur_buf, dosage_main, dosage_present, curMaleIndex, keepSampleCT, dosage_ct, &snpinfo, err)){
        LOGGER.e(0, err);
    }

    double af = snpinfo.af;
    double std = snpinfo.std;
    if(bHasPreAF){
        af = AFA1[gbuf->extractedMarkerIndex];
        std = 2.0 * af * (1.0 - af);
    }
    double maf = std::min(af, 1.0 - af);
    if(maf >= min_maf && maf <= max_maf && snpinfo.nMissRate >= dFilterMiss){
        gbuf->valid = true;
        gbuf->af = af;
        gbuf->nValidN = snpinfo.N;
        gbuf->nValidAllele = snpinfo.AlCount;
        gbuf->mean = 2.0 * af;
        gbuf->sd = std;

        if(bMakeGeno){
            if(isSexXY == 1){
                double weight;
                bool needWeight;
                setMaleWeight(weight, needWeight);
                if(needWeight){
                    for(int i = 0 ; i < keepMaleSampleCT; i++){
                        gbuf->geno[keepMaleExtractIndex[i]] *= weight;
                    }
                }
            }
        }
        if(bMakeMiss){
            gbuf->missing.resize(missPtrSize, 0);
        }
    }
}

void Geno::getGenoDouble_bed(uintptr_t *buf, int idx, GenoBufItem* gbuf){
    SNPInfo snpinfo;
    uintptr_t *cur_buf = buf + idx * bedRawGenoBuf1PtrSize;
    uint8_t isSexXY = isMarkersSexXYs[curBufferIndex];
    bool hasNoHET = true;
    if(isSexXY != 1){
        PgenReader::CountHardFreqMissExt(cur_buf, keepMaskInterPtr, rawSampleCT, keepSampleCT, &snpinfo, f_std);
    }else{
        string errmsg;
        hasNoHET = PgenReader::CountHardFreqMissExtX(cur_buf, keepMaskInterPtr, maleMaskInterPtr, rawSampleCT, keepSampleCT, keepMaleSampleCT, &snpinfo, errmsg, iDC==1, f_std);
    }
    uint32_t curExtractIndex = gbuf->extractedMarkerIndex;
    bool isEffRev = marker->isEffecRev(curExtractIndex);
    double af = isEffRev ? (1.0 - snpinfo.af) : snpinfo.af;
    if(bHasPreAF){
        af = AFA1[curExtractIndex];
        snpinfo.mean = 2 * af;
        snpinfo.std = 2 * af * (1.0 - af); 
    }
    double maf = std::min(af, 1.0 - af);
    if(maf >= min_maf && maf <= max_maf){
        if(snpinfo.nMissRate >= dFilterMiss){
            gbuf->valid = true;
            gbuf->af = af;
            gbuf->nValidN = snpinfo.N;
            gbuf->nValidAllele = snpinfo.AlCount;

            if(bGRM){
                gbuf->mean = 2.0 * af;
            }else{
                gbuf->mean = snpinfo.mean;
            }
            if(f_std){
                gbuf->sd = snpinfo.std;
            }else{
                gbuf->sd = gbuf->mean * (1.0 - af);
            }
            if(bMakeGeno){
                double mu = gbuf->mean;
                double sd = gbuf->sd;
                if(sd < 1.0e-50){
                    gbuf->valid = false;
                    return;
                }

                double center_value = 0.0;
                double rdev = 1.0;
                double a0, a1, a2, na;

                if(!bGRMDom){
                    if(bGenoCenter){
                        center_value = mu;
                    }
                    if(bGenoStd){
                        rdev = sqrt(1.0 / sd);
                    }
                    double aa0 = 0.0, aa2 = 2.0;
                    if(isEffRev){
                        double temp = aa0;
                        aa0 = aa2;
                        aa2 = aa0;
                    }

                    a0 = (aa0 - center_value) * rdev;
                    a1 = (1.0 - center_value) * rdev;
                    a2 = (aa2 - center_value) * rdev;
                    na = (mu - center_value) * rdev;
               }else{
                   double psq = 0.5 * mu * mu;
                   if(bGenoCenter)center_value = psq;
                   if(bGenoStd){
                       rdev = 1.0 / sd;
                   }
                   double aa0 = 0.0, aa2 = 2.0 * mu - 2.0;
                   if(isEffRev){
                       double temp = aa0;
                       aa0 = aa2;
                       aa2 = temp;
                   }
                   a0 = (aa0 - center_value) * rdev;
                   a1 = (mu - center_value) * rdev;
                   a2 = (aa2 - center_value) * rdev;
                   na = (psq - center_value)*rdev;
                }

                const double lookup[32] __attribute__ ((aligned (16))) = GET_TABLE16(a0, a1, a2, na);
                gbuf->geno.resize(keepSampleCT);
                uintptr_t * pmiss = NULL;
                if(bMakeMiss){
                    gbuf->missing.resize(missPtrSize); 
                    pmiss = gbuf->missing.data();
                }
                PgenReader::ExtractDoubleExt(cur_buf, keepMaskPtr, rawSampleCT, keepSampleCT, lookup, gbuf->geno.data(), pmiss); 
                if(isSexXY == 1){
                    double weight;
                    bool needWeight;
                    setMaleWeight(weight, needWeight);
                    if(needWeight){
                        if(bGRM){
                            for(int i = 0 ; i < keepMaleSampleCT; i++){
                                gbuf->geno[keepMaleExtractIndex[i]] *= weight;
                            }
                        }else{
                            double correctWeight = (weight - 1) * rdev * center_value;
                            for(int i = 0 ; i < keepMaleSampleCT; i++){
                                uint32_t curIndex = keepMaleExtractIndex[i];
                                gbuf->geno[curIndex] *= weight;
                                gbuf->geno[curIndex] += correctWeight;
                            }
                        }
                    }
                }
            }
            return;
        }
    }
    gbuf->valid = false;
}

void calDosage_bgen(uint32_t prob1, uint32_t prob2, uint64_t &dosage, uint32_t &prob1d){
    prob1d = prob1 * 2;
    dosage = prob1d + prob2;
}

void calDosagePhase_bgen(uint32_t prob1, uint32_t prob2, uint64_t &dosage, uint32_t &prob1d){
    dosage = prob1 + prob2;
    prob1d = 2 * prob1 * prob2;
}

void Geno::getGenoDouble_bgen(uintptr_t *buf, int idx, GenoBufItem* gbuf){
    SNPInfo snpinfo;
    uintptr_t *cur_buf = buf + idx * bgenRawGenoBuf1PtrSize;
    uint8_t *curbuf = (uint8_t*)cur_buf;
    int fileIndex = fileIndexBuf[curBufferIndex];

    int compressFormat = compressFormats[fileIndex];

    uint16_t L16;
    uint32_t L32;
    memcpy(&L16, curbuf, sizeof(L16));
    curbuf += sizeof(L16) + L16;
    memcpy(&L16, curbuf, sizeof(L16));
    curbuf += sizeof(L16) + L16;
    memcpy(&L16, curbuf, sizeof(L16));
    curbuf += sizeof(L16) + L16;
    curbuf += sizeof(uint32_t);
    memcpy(&L16, curbuf, sizeof(L16));
    curbuf += sizeof(L16);
    for(int i = 0; i < L16; i++){
        memcpy(&L32, curbuf, sizeof(L32));
        curbuf += sizeof(L32) + L32;
    }

    uint32_t len_comp, len_decomp;
    memcpy(&len_comp, curbuf, sizeof(len_comp));
    curbuf += sizeof(len_comp);

    if(compressFormat == 0){
        len_decomp = len_comp;
    }else{
        len_comp -= 4;
        memcpy(&len_decomp, curbuf, sizeof(len_decomp));
        curbuf += sizeof(len_decomp);
    }

    string error_promp = to_string(gbuf->extractedMarkerIndex) + "th SNP of [" + geno_files[fileIndex] + "]."; 
    uint8_t *dec_data;
    if(compressFormat != 0){
        dec_data = new uint8_t[len_decomp + 8];
        uint32_t curCompSize = len_comp;
        if(compressFormat == 1){
            uint32_t Ldecomp = len_decomp;
            int z_result = uncompress((Bytef*)dec_data, (uLongf*)&Ldecomp, (Bytef*)curbuf, curCompSize);
            if(z_result != Z_OK || len_decomp != Ldecomp){
                LOGGER.e(0, "decompressing genotype data error in " + error_promp); 
            }
        }else if(compressFormat == 2){
            uint64_t const rSize = ZSTD_getFrameContentSize((void*)curbuf, curCompSize);
            switch(rSize){
                case ZSTD_CONTENTSIZE_ERROR:
                    LOGGER.e(0, "not compressed by zstd in " + error_promp);
                    break;
                case ZSTD_CONTENTSIZE_UNKNOWN:
                    LOGGER.e(0, "original size unknown in " + error_promp);
                    break;
            }
            if(rSize != len_decomp){
                LOGGER.e(0, "size stated in the compressed file is different from " + error_promp);
            }
            size_t const dSize = ZSTD_decompress((void *)dec_data, len_decomp, (void*)curbuf, curCompSize); 
            if(ZSTD_isError(dSize)){
                LOGGER.e(0, "decompressing genotype error: " + string(ZSTD_getErrorName(dSize)) + " in " + error_promp);
            }
        }else{
            LOGGER.e(0, "unknown compress format in " + error_promp);
        }
    }else{
        dec_data = curbuf;
    }

    uint32_t n_sample = *(uint32_t *)dec_data;
    if(n_sample != rawCountSamples[fileIndex]){
        LOGGER.e(0, "inconsistent number of individuals in " + error_promp);
    }
    uint16_t num_alleles = *(uint16_t *)(dec_data + 4);
    if(num_alleles != 2){
        LOGGER.e(0, "multi-allelic SNPs detected in " + error_promp);
    }

    uint8_t min_ploidy = *(uint8_t *)(dec_data + 6);
    uint8_t max_ploidy = *(uint8_t *)(dec_data + 7);
    uint8_t * sample_ploidy = (uint8_t *)(dec_data + 8);
    if(min_ploidy != 2){
        LOGGER.e(0, "multiploidy detected in " + error_promp);
    }

    uint8_t *geno_prob = sample_ploidy + n_sample;
    uint8_t is_phased = *(geno_prob);
    uint8_t bits_prob = *(geno_prob+1);
    uint8_t* X_prob = geno_prob + 2;
    uint32_t len_prob = len_decomp - n_sample - 10;
    void (*calFunc)(uint32_t, uint32_t, uint64_t&, uint32_t &prob1d);
    if(is_phased){
        calFunc = &calDosagePhase_bgen;
    }else{
        calFunc = &calDosage_bgen;
    }

    uint8_t double_bits_prob = bits_prob * 2;
    vector<uint32_t> miss_index;

    uint8_t isSexXY = isMarkersSexXYs[curBufferIndex];

    uint64_t mask = (1U << bits_prob) - 1;
    uint64_t dosage_sum = 0, fij_sum = 0, dosage2_sum = 0;
    uint32_t validN = 0;
    uint32_t validAllele = 0;
    vector<uint32_t> dosages(keepSampleCT);
    uint32_t max_dos = mask * 2 + 1;
    bool has_miss = false;

    uint32_t curSampleCT = keepSampleCT;
    vector<uint32_t> *curSampleIndexPtr = &sampleKeepIndex;

    for(int j = 0; j < curSampleCT; j++){
        uint32_t sindex = (*curSampleIndexPtr)[j];
        uint8_t item_ploidy = sample_ploidy[sindex];
        if(item_ploidy > 128){
            miss_index.push_back(sindex);
            has_miss = true;
            dosages[j] = max_dos;
        }else if(item_ploidy == 2){
            uint32_t start_bits = sindex * double_bits_prob;
            uint64_t geno_temp;
            memcpy(&geno_temp, &(X_prob[start_bits/CHAR_BIT]), sizeof(geno_temp));
            geno_temp = geno_temp >> (start_bits % CHAR_BIT);
            uint32_t prob1 = geno_temp & mask;
            uint32_t prob2 = (geno_temp >> bits_prob) & mask;
            uint32_t prob1d;
            uint64_t dosage;
            calFunc(prob1, prob2, dosage, prob1d);
            dosages[j] = dosage;
            dosage_sum += dosage;
            dosage2_sum += dosage * dosage;
            fij_sum += prob1d;
            validN++;
            validAllele += 2;
        }else{
            LOGGER.e(0, "multiploidy detected in " + error_promp);
        }
    }

    double dosage_sum_half = dosage_sum;
    double dosage2_sum_half = dosage2_sum;
    if(isSexXY == 1){
        for(int j = 0; j < keepMaleSampleCT; j++){
            uint32_t sindex = keepMaleIndex[j];
            uint8_t item_ploidy = sample_ploidy[sindex];
            if(item_ploidy == 2){
                uint32_t start_bits = sindex * double_bits_prob;
                uint64_t geno_temp;
                memcpy(&geno_temp, &(X_prob[start_bits/CHAR_BIT]), sizeof(geno_temp));
                geno_temp = geno_temp >> (start_bits % CHAR_BIT);
                uint32_t prob1 = geno_temp & mask;
                uint32_t prob2 = (geno_temp >> bits_prob) & mask;
                uint32_t prob1d;
                uint64_t dosage;
                calFunc(prob1, prob2, dosage, prob1d);
                uint32_t prob1_true = prob1d / 2;
                dosage_sum_half -= prob1_true;
                dosage2_sum_half -= (dosage * dosage - (uint64_t)prob1_true * prob1_true);
                validAllele--;
            }
        }
    }

    if(compressFormat != 0){
        delete[] dec_data;
    }

    double maskd = (double)mask;
    double af = (double)dosage_sum_half / maskd / validAllele;
    double mean;
    bool bEffRev = this->marker->isEffecRev(gbuf->extractedMarkerIndex);
    if(bEffRev){
        af = 1.0 - af;
    }
    double std = 2.0 * af * (1.0 - af);
    double info = 0.0;
    double mask2 = mask * mask;
    if(std < 1e-50){
        info = 1.0;
    }else{
        double dos2_fij_sum = (double)(dosage_sum + fij_sum)/ maskd - (double)dosage2_sum / mask2; 
        info = 1.0 - dos2_fij_sum / (std * validN);
    }
    if(is_phased){
        info = 1;
    }

    if(bHasPreAF){
        af = AFA1[gbuf->extractedMarkerIndex];
        mean = 2.0 * af;
        std = 2.0 * af * (1.0 - af);
    }else{
        if(iDC == 1 && (!bGRM)){
            double dos_double = (double)dosage_sum / maskd;
            mean = dos_double / validN;
            std = ((double)dosage2_sum / mask2 - dos_double * mean)/(validN - 1);
        }else{
            double dos_double = (double)dosage_sum_half / maskd;
            mean = dos_double / validN;
            std = ((double)dosage2_sum_half / mask2 - dos_double * mean)/(validN - 1);
        }
    }

    double maf = std::min(af, 1.0 - af);
    if(maf >= min_maf && maf <= max_maf){
        double nMissRate = 1.0*validN / curSampleCT;
        if(nMissRate >= dFilterMiss && info >= dFilterInfo){
            gbuf->valid = true;
            gbuf->af = af;
            gbuf->nValidN = validN;
            gbuf->nValidAllele = validAllele;
            gbuf->info = is_phased ? (std::numeric_limits<double>::quiet_NaN()) : info;
            gbuf->mean = mean;
            gbuf->sd = std;
            if(bMakeGeno){
                double mu = gbuf->mean;
                if(std < 1.0e-50){
                    gbuf->valid = false;
                    return;
                }

                double* dos_lookup = new double[max_dos + 2];
                double center_value = 0.0;
                double rdev = 1.0;
                if(!bGRMDom){
                    if(bGenoCenter){
                        center_value = mu;
                    }
                    if(bGenoStd){
                        rdev = sqrt(1.0 / std);
                    }
                    for(uint32_t i = 0; i < max_dos; i++){
                        double tdos = (double)i / mask;
                        tdos = bEffRev ? (2.0 - tdos) : tdos;
                        dos_lookup[i] = (tdos - center_value) *rdev;
                    }
                    dos_lookup[max_dos] = ( mu - center_value) * rdev;
                }else{
                    uint32_t cut05 = ceil(0.5 * mask);
                    uint32_t cut15 = ceil(1.5 * mask);
                    double psq = 0.5 * mu * mu;
                    double dos00 = bEffRev ? (2.0 * mu - 2.0) : 0.0;
                    double dos01 = mu;
                    double dos10 = bEffRev ? 0.0 : (2.0 * mu - 2.0);
                    double dosna = psq;

                    if(bGenoCenter)center_value = psq;
                    if(bGenoStd){
                        rdev = 1.0 / std;
                    }
                    double a0 = (dos00 - center_value) * rdev;
                    double a1 = (dos01 - center_value) * rdev;
                    double a2 = (dos10 - center_value) * rdev;
                    double na = (dosna - center_value) * rdev;
                    
                    for(uint32_t i = 0; i < cut05; i++){
                        dos_lookup[i] = a0;
                    }
                    for(uint32_t i = cut05; i < cut15; i++){
                        dos_lookup[i] = a1;
                    }
                    for(uint32_t i = cut15; i < max_dos; i++){
                        dos_lookup[i] = a2;
                    }
                    dos_lookup[max_dos] = (dosna - center_value) * rdev;
                }

                gbuf->geno.resize(curSampleCT);
                for(int j = 0; j < curSampleCT; j++){
                    gbuf->geno[j] = dos_lookup[dosages[j]];
                }
                delete[] dos_lookup;
                if(isSexXY == 1){
                    double weight;
                    bool needWeight;
                    setMaleWeight(weight, needWeight);
                    if(needWeight){
                        if(bGRM){
                            for(int i = 0 ; i < keepMaleSampleCT; i++){
                                gbuf->geno[keepMaleExtractIndex[i]] *= weight;
                            }
                        }else{
                            double correctWeight = (weight - 1) * rdev * center_value;
                            for(int i = 0 ; i < keepMaleSampleCT; i++){
                                uint32_t curIndex = keepMaleExtractIndex[i];
                                gbuf->geno[curIndex] *= weight;
                                gbuf->geno[curIndex] += correctWeight;
                            }
                        }
                    }
                }
            }
            if(bMakeMiss){
                gbuf->missing.resize(missPtrSize, 0); 
                const int ptrsize = sizeof(uintptr_t) * CHAR_BIT;
                for(int j = 0; j < (int)miss_index.size(); j++){
                    int cur_index = miss_index[j];
                    gbuf->missing[cur_index / ptrsize] |= (1UL << (cur_index % ptrsize));
                }
            }
            return;
        }
    }
    gbuf->valid = false;
}

void Geno::loopDouble(const vector<uint32_t> &extractIndex, int numMarkerBuf,
                      bool bMakeGeno, bool bGenoCenter, bool bGenoStd, bool bMakeMiss,
                      vector<function<void(uintptr_t *buf, const vector<uint32_t> &exIndex)>> callbacks,
                      bool showLog)
{
    // ── Common pre-logic (was in preGenoDouble) ──────────────────────────
    sampleKeepIndex    = pheno->get_index_keep();
    keepSampleCT       = sampleKeepIndex.size();
    rawSampleCT        = pheno->count_raw();
    numMarkerBlock     = numMarkerBuf;
    keepSexIndex       = pheno->getSexValidRawIndex();
    keepMaleIndex      = pheno->getMaleRawIndex();
    keepMaleExtractIndex = pheno->getMaleExtractIndex();
    keepSexSampleCT    = keepSexIndex.size();
    keepMaleSampleCT   = keepMaleIndex.size();
    this->bMakeGeno    = bMakeGeno;
    this->bGenoCenter  = bGenoCenter;
    this->bGenoStd     = bGenoStd;
    this->bMakeMiss    = bMakeMiss;

    // ── Create format-specific backend (runs format pre-logic) ───────────
    std::unique_ptr<GenoBackend> backend;
    if (genoFormat == "BED") {
        backend = makeBedBackend(*this);
    } else if (genoFormat == "PGEN") {
        backend = makePgenBackend(*this);
    } else if (genoFormat == "BGEN") {
        backend = makeBgenBackend(*this);
    } else {
        LOGGER.e(0, "loopDouble: unknown genotype format: " + genoFormat);
    }

    // ── Build baseIndexLookup (depends on rawCountSNPs set by backend) ───
    baseIndexLookup.clear();
    baseIndexLookup.push_back(0);
    int32_t sumIndex = 0;
    for (int i = 0; i < (int)geno_files.size() - 1; ++i) {
        sumIndex += rawCountSNPs[i];
        baseIndexLookup.push_back(sumIndex);
    }
    missPtrSize = PgenReader::GetSubsetMaskSize(keepSampleCT);

    // ── Single-slot block metadata (index 0 is always current block) ─────
    isMarkersSexXYs.assign(1, 0);
    fileIndexBuf.assign(1, 0);
    curBufferIndex = 0;

    // ── Thread pools ──────────────────────────────────────────────────────
    //    io_pool  : 1 thread — file handles and PgenReader are not thread-safe.
    //    cpu_pool : 1 thread — Option A from the coroutine migration plan.
    //      Callbacks use #pragma omp parallel for internally, so OMP expands
    //      to omp_get_max_threads() workers inside each block.  The full
    //      machine is therefore used; cpu_pool only acts as the scheduler
    //      that launches each OMP region.
    //
    //    *** Do NOT raise cpu_pool above 1 without first removing or
    //        nesting-guarding every OMP pragma in GRM.cpp and FastFAM.cpp.
    //        Concurrent OMP regions from multiple cpu_pool threads would
    //        oversubscribe the machine and degrade performance. ***
    exec::static_thread_pool io_pool{1};
    exec::static_thread_pool cpu_pool{1};   // OMP handles intra-block parallelism

    LOGGER.ts("LOOP_GENO_TOT");
    LOGGER.ts("LOOP_GENO_PRE");
    uint32_t nFinishedMarker = 0;
    const uint32_t nTMarker  = static_cast<uint32_t>(extractIndex.size());
    int pre_block = 0;

    // ── Adapter: GenoBlock → legacy (uintptr_t*, vector<uint32_t>) ───────
    auto blockCallback = [&](GenoBlock &block) {
        isMarkersSexXYs[0] = block.isSexXY;
        fileIndexBuf[0]    = block.fileIndex;
        curBufferIndex     = 0;

        for (auto &cb : callbacks) {
            cb(block.buf.data(), block.extractIndex);
        }

        nFinishedMarker += block.numMarkers;
        if (showLog) {
            int cur_block = nFinishedMarker >> 14;
            if (cur_block > pre_block) {
                pre_block = cur_block;
                float time_p = LOGGER.tp("LOOP_GENO_PRE");
                if (time_p > 300) {
                    LOGGER.ts("LOOP_GENO_PRE");
                    float elapse_time      = LOGGER.tp("LOOP_GENO_TOT");
                    float finished_percent = (float)nFinishedMarker / nTMarker;
                    float remain_time      = (1.0f / finished_percent - 1) * elapse_time / 60;
                    std::ostringstream ss;
                    ss << std::fixed << std::setprecision(1)
                       << finished_percent * 100
                       << "% Estimated time remaining " << remain_time << " min";
                    LOGGER.i(1, ss.str());
                }
            }
        }
    };

    // ── Run coroutine pipeline (blocks until all markers processed) ───────
    stdexec::sync_wait(backend->stream(
        io_pool.get_scheduler(),
        cpu_pool.get_scheduler(),
        extractIndex,
        std::move(blockCallback)));

    if (showLog) {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(1)
           << "100% finished in " << LOGGER.tp("LOOP_GENO_TOT") << " sec";
        LOGGER.i(1, ss.str());
        LOGGER << nFinishedMarker << " SNPs have been processed." << std::endl;
    }
    // backend destructor runs end-logic (masks freed, files closed).
}

//function for BGEN
inline void bgen12ExtractVal(uint64_t val, uint8_t bit_prob, uint64_t mask, uint64_t &v1, uint64_t &v2){
    v1 = val & mask;
    v2 = (val >> bit_prob) & mask;
}

// call genotype
void dosageFunc(uint32_t prob1, uint32_t prob2, uint64_t mask, uint32_t cutVal, uint32_t A1U, uint32_t A1L, double &gval, bool &miss){
    gval = (double)(prob1 * 2 + prob2) / mask;
}


void dosageCallFunc(uint32_t prob1, uint32_t prob2, uint64_t mask, uint32_t cutVal, uint32_t A1U, uint32_t A1L, double &gval, bool &miss){
    uint32_t dos = prob1 * 2 + prob2;
    if(dos > A1U){
        gval = 2.0;
    }else if(dos < A1L){
        gval = 0.0;
    }else{
        gval = 1.0;
    }
}

void hardCallFunc(uint32_t prob1, uint32_t prob2, uint64_t mask, uint32_t cutVal, uint32_t A1U, uint32_t A1L, double &gval, bool &miss){
    miss = false;
    if(prob1 >= cutVal){
        gval = 2.0;
    }else if(prob2 >= cutVal){
        gval = 1.0;
    }else if(mask - prob1 - prob2 >= cutVal){
        gval = 0.0;
    }else{
        gval = 0.0;
        miss = true;
    }
};
        /* obsoleted functions hard calling
        uint8_t double_bits_prob = bits_prob * 2;
        uint32_t A1U = floor(mask * 1.5);
        uint32_t A1L = ceil(mask * 0.5);
        uint32_t cutVal = ceil(mask * options_d["hard_call_thresh"]);

        void (*callFunc)(uint32_t prob1, uint32_t prob2, uint64_t mask, uint32_t cutVal, uint32_t A1U, uint32_t A1L, double &gval, bool &miss);
        if(options.find("dosage_call") != options.end()){
            callFunc = dosageCallFunc;
        }else if(options.find("dosage") != options.end()){
            callFunc = dosageFunc;
        }else{
            callFunc = hardCallFunc;
        }

        auto callFuncB = [callFunc, mask, cutVal, A1U, A1L](uint32_t prob1, uint32_t prob2, double &gval, bool&miss){
            callFunc(prob1, prob2, mask, cutVal, A1U, A1L, gval, miss);
        };
    
        //vector<uint32_t> prob1s(n_res), prob2s(n_res);
        vector<uint32_t> miss_index;
        double sum_geno = 0;
        double *base_geno = (gbuf->geno).data() + i * n_res;
        uint64_t *base_miss;
        uint32_t n_valid = 0;
        if(gbuf->saveMiss) base_miss = (gbuf->miss).data() + i * n_bmiss; 

        uint64_t dosage_sum = 0, fij_dosage2_sum = 0;
        bool *lookup_miss = new bool[mask + 1];


        for(int j = 0; j < n_res; j++){
            uint32_t sindex = sampleKeepIndex[j];
            uint8_t item_ploidy = sample_ploidy[sindex];
            if(item_ploidy > 128){
                miss_index.push_back(sindex);
            }else if(item_ploidy == 2){
                uint32_t start_bits = sindex * double_bits_prob;
                uint64_t geno_temp;
                memcpy(&geno_temp, &(X_prob[start_bits/CHAR_BIT]), sizeof(uint64_t));
                geno_temp = geno_temp >> (start_bits % CHAR_BIT);
                uint32_t prob1 = geno_temp & mask;
                uint32_t prob2 = (geno_temp >> bits_prob) & mask;
                uint32_t prob1d = prob1 * 2;
                uint32_t dosage = prob1d + prob2;

                

                dosage_sum += dosage;

                uint32_t fij = dosage + prob1d;
                fij_dosage2_sum += (fij - dosage * dosage);
            }else{
               LOGGER.e(0, "multi-allelic SNPs detected in " + error_promp);
            }
        }

        for(int j = 0; j < n_res; j++){
            uint32_t sindex = sampleKeepIndex[j];
            uint8_t item_ploidy = sample_ploidy[sindex];
            if(item_ploidy > 128){
                miss_index.push_back(sindex);
            }else if(item_ploidy == 2){
                uint32_t start_bits = sindex * double_bits_prob;
                uint64_t geno_temp;
                memcpy(&geno_temp, &(X_prob[start_bits/CHAR_BIT]), sizeof(uint64_t));
                geno_temp = geno_temp >> (start_bits % CHAR_BIT);
                uint32_t prob1 = geno_temp & mask;
                uint32_t prob2 = (geno_temp >> bits_prob) & mask;

                //total_prob += ((prob1 << 1) + prob2);
                //prob1s[i] = prob1;
                //prob2s[i] = prob2;
                double gval;
                bool miss;
                callFuncB(prob1, prob2, gval, miss);
                base_geno[j] = gval;
                if(miss){
                    miss_index.push_back(j);
                    if(gbuf->saveMiss) base_miss[j / 64] |= (1UL << (j % 64)); 
                }else{
                    n_valid++;
                }
           }else{
               LOGGER.e(0, "multi-allelic SNPs detected in " + error_promp);
           }
        }
        delete[] dec_data;
        Eigen::VectorXd Vgeno = Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,1> > (base_geno, n_res);
        double af = Vgeno.sum() / (2*n_valid);
        (gbuf->validN)[i] = n_valid;
        (gbuf->af)[i] = af;
        */

union Geno_prob{
    char byte[4];
    uint32_t value = 0;
};


void Geno::bgen2bed(const vector<uint32_t> &raw_marker_index){
    LOGGER << "Old bgen 2 bed" << std::endl;
    LOGGER.ts("LOOP_BGEN_BED");
    LOGGER.ts("LOOP_BGEN_TOT");
    vector<uint32_t>& index_keep = pheno->get_index_keep();

    /*
    std::ofstream out((options["out"] + "_sample_index.txt").c_str());
    for(auto & item : index_keep){
        out << item << std::endl;
    }
    out.close();
    */
    
    auto buf_size = (num_raw_sample + 31) / 32;
    size_t buf_size_byte = buf_size * 8;

    int num_marker = 1;
    int num_markers = raw_marker_index.size();
    LOGGER << "samples: " << num_raw_sample << ", keep_sample: " << index_keep.size() << std::endl;
    LOGGER << "Markers: " << num_markers << std::endl;

    FILE * h_bgen = fopen(options["bgen_file"].c_str(), "rb");
    /*
    std::ofstream infos("exmaple.txt");
    infos << "index\traw_index\tpos\tLen_comp\tLen_decomp" << std::endl;
    */
    #pragma omp parallel for schedule(static) ordered
    for(uint32_t index = 0; index < num_markers; index++){
        //LOGGER.i(0, to_string(index) + "NUM_thread: " + to_string(omp_get_max_threads()));
        auto raw_index = raw_marker_index[index];
        uint64_t *buf = new uint64_t[buf_size]();
        uint64_t byte_pos, byte_size;
        this->marker->getStartPosSize(raw_index, byte_pos, byte_size);
        uint32_t len_comp, len_decomp;
        char * snp_data;

        #pragma omp ordered
        {
            fseek(h_bgen, byte_pos, SEEK_SET);
            len_comp = read1Byte<uint32_t>(h_bgen) - 4;
            len_decomp = read1Byte<uint32_t>(h_bgen);
            snp_data = new char[len_comp];
            readBytes(h_bgen, len_comp, snp_data);
            //infos << index << "\t" << raw_index << "\t" << byte_pos << "\t" << len_comp << "\t" << len_decomp << std::endl;
        }
        uLongf dec_size = len_decomp;

        char * dec_data =  new char[len_decomp];
        int z_result = uncompress((Bytef*)dec_data, &dec_size, (Bytef*)snp_data, len_comp);
        delete[] snp_data;
        if(z_result == Z_MEM_ERROR || z_result == Z_BUF_ERROR || dec_size != len_decomp){
            LOGGER.e(0, "decompressing genotype data error in " + to_string(raw_index) + "th SNP."); 
        }

        uint32_t n_sample = *(uint32_t *)dec_data;
        if(n_sample != num_raw_sample){
            LOGGER.e(0, "inconsistent number of samples in " + to_string(raw_index) + "th SNP." );
        }
        uint16_t num_alleles = *(uint16_t *)(dec_data + 4);
        if(num_alleles != 2){
            LOGGER.e(0, "multi-allelic SNPs detected likely because the bgen file is malformed.");
        }

        uint8_t min_ploidy = *(uint8_t *)(dec_data + 6);//2
        uint8_t max_ploidy = *(uint8_t *)(dec_data + 7); //2
        uint8_t * sample_ploidy = (uint8_t *)(dec_data + 8);

        uint8_t *geno_prob = sample_ploidy + n_sample;
        uint8_t is_phased = *(geno_prob);
        uint8_t bits_prob = *(geno_prob+1);
        uint8_t* X_prob = geno_prob + 2;
        uint32_t len_prob = len_decomp - n_sample - 10;
        if(is_phased){
            LOGGER.e(0, "GCTA does not support phased data currently.");
        }

        int byte_per_prob = bits_prob / 8;
        int double_byte_per_prob = byte_per_prob * 2;
        if(bits_prob % 8 != 0){
            LOGGER.e(0, "GCTA does not support probability bits other than in byte units.");
        }

        if(len_prob != double_byte_per_prob * n_sample){
            LOGGER.e(0, "malformed data in " + to_string(raw_index) + "th SNP.");
        }
        /*
        infos << index << "_2\t" << raw_index << "\t" << byte_pos << "\t" << len_comp << "\t" << len_prob << std::endl;
        FILE *obgen = fopen((to_string(index) + ".bin").c_str(), "wb");
        fwrite(snp_data, sizeof(char), len_comp, obgen);
        fwrite(X_prob, sizeof(char), len_prob, obgen);
        fclose(obgen);
        */

        uint32_t base_value = (1 << bits_prob) - 1;

        uint8_t *buf_ptr = (uint8_t *)buf;
        if(options.find("dosage_call") == options.end()){
            uint32_t cut_value = ceil(base_value * options_d["hard_call_thresh"]);
            for(uint32_t i = 0; i < num_keep_sample; i++){
                uint32_t item_byte = i >> 2;
                uint32_t move_byte = (i & 3) << 1;

                uint32_t sindex = index_keep[i];
                uint8_t item_ploidy = sample_ploidy[sindex];

                uint8_t geno_value;
                if(item_ploidy > 128){
                    geno_value = 1;
                }else if(item_ploidy == 2){
                    auto base = sindex * double_byte_per_prob;
                    auto base1 = base + byte_per_prob;
                    Geno_prob prob_item;
                    Geno_prob prob_item1;
                    /*
                       memcpy(prob_item.byte, X_prob + base,  byte_per_prob); 
                       memcpy(prob_item1.byte, X_prob + base1, byte_per_prob); 
                       */
                    for(int i = 0 ; i != byte_per_prob; i++){
                        prob_item.byte[i] = X_prob[base + i];
                        prob_item1.byte[i] = X_prob[base1 + i];
                    }

                    uint32_t t1 = prob_item.value;
                    uint32_t t2 = prob_item1.value;
                    uint32_t t3 = base_value - t1 - t2;
                    if(t1 >= cut_value){
                        geno_value = 0;
                    }else if(t2 >= cut_value){
                        geno_value = 2;
                    }else if(t3 >= cut_value){
                        geno_value = 3;
                    }else{
                        geno_value = 1;
                    }
                }else{
                    LOGGER.e(0, "multi-allelic SNPs detected in the " + to_string(raw_index) + "th SNP.");
                }
                buf_ptr[item_byte] += geno_value << move_byte;
            }
        }else{
            uint32_t A1U = floor(base_value * 1.5);
            uint32_t A1L = ceil(base_value * 0.5);


            for(uint32_t i = 0; i < num_keep_sample; i++){
                uint32_t item_byte = i >> 2;
                uint32_t move_byte = (i & 3) << 1;

                uint32_t sindex = index_keep[i];
                uint8_t item_ploidy = sample_ploidy[sindex];

                uint8_t geno_value;
                if(item_ploidy > 128){
                    //missing
                    geno_value = 1;
                }else if(item_ploidy == 2){
                    auto base = sindex * double_byte_per_prob;
                    auto base1 = base + byte_per_prob;
                    Geno_prob prob_item;
                    Geno_prob prob_item1;
                    /*
                       memcpy(prob_item.byte, X_prob + base,  byte_per_prob); 
                       memcpy(prob_item1.byte, X_prob + base1, byte_per_prob); 
                       */
                    for(int i = 0 ; i != byte_per_prob; i++){
                        prob_item.byte[i] = X_prob[base + i];
                        prob_item1.byte[i] = X_prob[base1 + i];
                    }

                    uint32_t t1 = prob_item.value;
                    uint32_t t2 = prob_item1.value;
                    uint32_t dosageA = 2 * t1 + t2;
                    if(dosageA > A1U){
                        geno_value = 0;
                    }else if(dosageA < A1L){
                        geno_value = 3;
                    }else{
                        geno_value = 2;
                    }
                }else{
                    LOGGER.e(0, " multi-allelic SNPs detected in the " + to_string(raw_index) + "th SNP.");
                }
                buf_ptr[item_byte] += geno_value << move_byte;
            }
 
        }
        //LOGGER.i(0, "MIDDLE: " + to_string(index) + "NUM_thread: " + to_string(omp_get_max_threads()));

        // exactly one 'ordered' directive must appear in the loop body of an enclosing directive
        // #pragma omp ordered
        save_bed(buf, num_marker);
        delete[] buf;
        delete[] dec_data;
        //#pragma omp ordered
        //LOGGER.i(0, "Finished " + to_string(index) + "NUM_thread: " + to_string(omp_get_max_threads()));
        if(index % 10000 == 0){
            float time_p = LOGGER.tp("LOOP_BGEN_BED");
            if(time_p > 300){
                LOGGER.ts("LOOP_BGEN_BED");
                float elapse_time = LOGGER.tp("LOOP_BGEN_TOT");
                float finished_percent = (float) index / num_markers;
                float remain_time = (1.0 / finished_percent - 1) * elapse_time / 60;

                std::ostringstream ss;
                ss << std::fixed << std::setprecision(1) << finished_percent * 100 << "% Estimated time remaining " << remain_time << " min"; 
                
                LOGGER.i(1, ss.str());
            }
        }

    }
    //infos.close();
    closeOut();
    fclose(h_bgen);
}

// extracted and revised from plink2.0
// GPL v3, license detailed on github
// https://github.com/chrchang/plink-ng

void Geno::save_bed(uint64_t *buf, int num_marker){
    static string err_string = "can't write to [" + options["out"] + ".bed].";
    static bool inited = false;
    if(!inited){
        hOut = fopen((options["out"] + ".bed").c_str(), "wb");
        if(hOut == NULL){
            LOGGER.e(0, err_string);
        }
        uint8_t header[3] = {0x6c, 0x1b, 0x01};
        if(3 != fwrite(header, sizeof(uint8_t), 3, hOut)){
            LOGGER.e(0, err_string);
        }
        inited = true;
    }
    uint64_t base_buffer = 0;
    for(int i = 0; i < num_marker; i++){
        uint8_t *buffer = (uint8_t *)(buf + base_buffer);
        if(fwrite(buffer, sizeof(uint8_t), num_byte_keep_geno1, hOut) != num_byte_keep_geno1){
            LOGGER.e(0, err_string);
        }
        base_buffer += num_item_1geno;
    }
}

void Geno::closeOut(){
    fclose(hOut);
}

void copy_quaterarr_nonempty_subset(uint64_t* raw_quaterarr[], const uint64_t* subset_mask, uint32_t raw_quaterarr_entry_ct, uint32_t subset_entry_ct, uint64_t* output_quaterarr[], const int num_marker) {
    // in plink 2.0, we probably want (0-based) bit raw_quaterarr_entry_ct of
    // subset_mask to be always allocated and unset.  This removes a few special
    // cases re: iterating past the end of arrays.
    static const uint32_t kBitsPerWordD2 = 32;

    uint64_t cur_output_word[num_marker];
    memset(cur_output_word, 0, num_marker * 8);

    uint64_t* output_quaterarr_iter[num_marker];
    uint64_t* output_quaterarr_last[num_marker];
    for(int i = 0; i != num_marker; i++){
        output_quaterarr_iter[i] = output_quaterarr[i];
        output_quaterarr_last[i] = &(output_quaterarr[i][subset_entry_ct / kBitsPerWordD2]);
    }
    const uint32_t word_write_halfshift_end = subset_entry_ct % kBitsPerWordD2;
    uint32_t word_write_halfshift = 0;
    // if <= 2/3-filled, use sparse copy algorithm
    // (tried copy_bitarr_subset() approach, that actually worsened things)
    if (subset_entry_ct * (3 * k1LU) <= raw_quaterarr_entry_ct * (2 * k1LU)) {
        uint32_t subset_mask_widx = 0;
        while (1) {
            const uint64_t cur_include_word = subset_mask[subset_mask_widx];
            if (cur_include_word) {
                uint32_t wordhalf_idx = 0;
                uint32_t cur_include_halfword = (halfword_t)cur_include_word;
                while (1) {
                    if (cur_include_halfword) {
                        uint64_t raw_quaterarr_word[num_marker];
                        uint32_t temp_index = subset_mask_widx * 2 + wordhalf_idx;
                        for(int i = 0; i != num_marker; i++){
                            raw_quaterarr_word[i] = raw_quaterarr[i][temp_index];
                        }
                        do {
                            uint32_t rqa_idx_lowbits = CTZ64U(cur_include_halfword);
                            uint32_t lshift = word_write_halfshift * 2; 
                            uint32_t rshift = rqa_idx_lowbits * 2;
                            for(int i = 0; i != num_marker; i++){
                                cur_output_word[i] |= ((raw_quaterarr_word[i] >> rshift) & 3) << lshift;
                            }
                            if (++word_write_halfshift == kBitsPerWordD2) {
                                for(int i = 0; i != num_marker; i++){
                                    *(output_quaterarr_iter[i])++ = cur_output_word[i];
                                }
                                word_write_halfshift = 0;
                                //cur_output_word = 0;
                                memset(cur_output_word, 0, num_marker * 8);
                            }
                            cur_include_halfword &= cur_include_halfword - 1;
                        } while (cur_include_halfword);
                    }
                    if (wordhalf_idx) {
                        break;
                    }
                    ++wordhalf_idx;
                    cur_include_halfword = cur_include_word >> kBitsPerWordD2;
                }
                if (output_quaterarr_iter[0] == output_quaterarr_last[0]) {
                    if (word_write_halfshift == word_write_halfshift_end) {
                        if (word_write_halfshift_end) {
                            for(int i = 0; i != num_marker; i++){
                                *(output_quaterarr_last[i]) = cur_output_word[i];
                            }
                        }
                        return;
                    }
                }
            }
            ++subset_mask_widx;
        }
    }

    const uint64_t* raw_quaterarr_iter[num_marker];
    for(int i = 0; i != num_marker; i++){
        raw_quaterarr_iter[i] = raw_quaterarr[i];
    }
    //const uint64_t* raw_quaterarr_iter = raw_quaterarr;
    while (1) {
        const uint64_t cur_include_word = *subset_mask++;
        uint32_t wordhalf_idx = 0;
        uint64_t cur_include_halfword = (halfword_t)cur_include_word;
        while (1) {
           // uintptr_t raw_quaterarr_word = *raw_quaterarr_iter++;
            uint64_t raw_quaterarr_word[num_marker];
            for(int i = 0; i != num_marker; i++){
                raw_quaterarr_word[i] = *(raw_quaterarr_iter[i]++);
            }
            while (cur_include_halfword) {
                uint32_t rqa_idx_lowbits = CTZ64U(cur_include_halfword); // tailing zero
                uint64_t halfword_invshifted = (~cur_include_halfword) >> rqa_idx_lowbits;
                uint64_t raw_quaterarr_curblock_unmasked[num_marker];
                int m_bit = rqa_idx_lowbits * 2;
                for(int i = 0; i != num_marker; i++){
                    raw_quaterarr_curblock_unmasked[i] = raw_quaterarr_word[i] >> m_bit; 
                }
                //uintptr_t raw_quaterarr_curblock_unmasked = raw_quaterarr_word >> (rqa_idx_lowbits * 2); // remove mask bit tailing zero, not to keep
                uint32_t rqa_block_len = CTZ64U(halfword_invshifted);  // find another keep
                uint32_t block_len_limit = kBitsPerWordD2 - word_write_halfshift;
                m_bit = 2 * word_write_halfshift;
                for(int i = 0; i != num_marker; i++){
                    cur_output_word[i] |= raw_quaterarr_curblock_unmasked[i] << m_bit; // avoid overwrite current saved bits
                }
                if (rqa_block_len < block_len_limit) { //2  16
                    word_write_halfshift += rqa_block_len; // 0 2
                    m_bit = 2 * word_write_halfshift;
                    uint64_t temp_mask = (k1LU << m_bit) - k1LU;
                    for(int i = 0; i != num_marker; i++){
                        cur_output_word[i] &= temp_mask; // mask high end, and keep low needed bits
                    }
                } else {
                    // no need to mask, extra bits vanish off the high end
                    for(int i = 0; i != num_marker; i++){
                        *(output_quaterarr_iter[i]++) = cur_output_word[i];
                    }
                    word_write_halfshift = rqa_block_len - block_len_limit;
                    if (word_write_halfshift) {
                        uint64_t t_mask = ((k1LU << (2 * word_write_halfshift)) - k1LU), mi_bit = 2 * block_len_limit;
                        for(int i = 0; i != num_marker; i++){
                            cur_output_word[i] = (raw_quaterarr_curblock_unmasked[i] >> mi_bit) & t_mask;
                        }
                    } else {
                        // avoid potential right-shift-[word length]
                        //cur_output_word = 0;
                        memset(cur_output_word, 0, num_marker * 8);
                    }
                }
                cur_include_halfword &= (~(k1LU << (rqa_block_len + rqa_idx_lowbits))) + k1LU;
            }
            if (wordhalf_idx) {
                break;
            }
            ++wordhalf_idx;
            cur_include_halfword = cur_include_word >> kBitsPerWordD2;
        }
        if (output_quaterarr_iter[0] == output_quaterarr_last[0]) {
            if (word_write_halfshift == word_write_halfshift_end) {
                if (word_write_halfshift_end) {
                    for(int i = 0; i != num_marker; i++){
                        *(output_quaterarr_last[i]) = cur_output_word[i];
                    }
                }
                return;
            }
        }
    }
}


// === Restored functions ===

void Geno::setMAF(double val){
    if(val < 0){
        LOGGER.e(0, "MAF can't be negative: " + to_string(val));
    }
    this->min_maf = val * (1.0 - Constants::SMALL_EPSILON);
    deterFilterMAF();
}

double Geno::getMAF(){
    return this->min_maf;
}

double Geno::getFilterInfo(){
    return this->dFilterInfo;
}

double Geno::getFilterMiss(){
    return this->dFilterMiss;
}

void Geno::deterFilterMAF(){
    if(max_maf < min_maf){
        LOGGER.e(0, "the value specified for --max-maf can't be smaller than that for --min-maf");
    }

    if(std::abs(min_maf) <= 1e-10 && std::abs(max_maf - 0.5) <= 1e-10){
        this->bFilterMAF = false;
    }else{
        this->bFilterMAF = true;
    }
}

void Geno::setMaxMAF(double val){
    if(val <= 0 || val > 0.5){
        LOGGER.e(0, "the value specified for --max-maf can't be negative or larger than 0.5");
    }
    this->max_maf = val * (1.0 + Constants::SMALL_EPSILON);
    deterFilterMAF();
}

void Geno::setFilterInfo(double val){
    this->dFilterInfo = val;
}

void Geno::setFilterMiss(double val){
    this->dFilterMiss = val;
}


void Geno::addOneFileOption(string key_store, string append_string, string key_name,
                                     map<string, vector<string>> options_in) {
    if(options_in.find(key_name) != options_in.end()){
        if(options_in[key_name].size() == 1){
            options[key_store] = options_in[key_name][0] + append_string;
        }else if(options_in[key_name].size() > 1){
            options[key_store] = options_in[key_name][0] + append_string;
            LOGGER.w(0, "Geno: multiple " + key_name + ", use the first one only" );
        }else{
            LOGGER.e(0, "no " + key_name + " parameter found");
        }
        std::ifstream f(options[key_store].c_str());
        if(!f.good()){
            LOGGER.e(0, key_name + " " + options[key_store] + " not found");
        }
        f.close();
    }
}

int Geno::registerOption(map<string, vector<string>>& options_in) {
    int return_value = 0;
    addOneFileOption("geno_file", ".bed", "--bfile", options_in);
    addOneFileOption("bgen_file", "", "--bgen", options_in);
    addOneFileOption("pgen_file", ".pgen", "--pfile", options_in);
    addOneFileOption("pgen_file", ".pgen", "--bpfile", options_in);

    addMFileListsOption("m_file", ".bed", "--mbfile", options_in, options);
    addMFileListsOption("mbgen_file", ".bgen", "--mbgen", options_in, options);
    addMFileListsOption("mpgen_file", ".pgen", "--mbpfile", options_in, options);
    addMFileListsOption("mpgen_file", ".pgen", "--mpfile", options_in, options);

    options_d["min_maf"] = 0.0;
    options_d["max_maf"] = 0.5;
    if(options_in.find("--maf") != options_in.end()){
        auto option = options_in["--maf"];
        if(option.size() == 1){
            try{
                options_d["min_maf"] = std::stod(option[0]);
            }catch(std::invalid_argument&){
                LOGGER.e(0, "invalid value for --maf");
            }
            if(options_d["min_maf"]<0.0 || options_d["max_maf"]>0.5){
                LOGGER.e(0, "value specified for--maf can't be smaller than 0 or larger than 0.5");
            }else if(options_d["min_maf"] == 0.0){
                options_in["--nofilter"] = {};
            }

        }else{
            LOGGER.e(0, "GCTA does not support multiple values for --maf currently");
        }
        options_in.erase("--maf");
    }

     if(options_in.find("--max-maf") != options_in.end()){
        auto option = options_in["--max-maf"];
        if(option.size() == 1){
            try{
                options_d["max_maf"] = std::stod(option[0]);
           }catch(std::invalid_argument&){
                LOGGER.e(0, "invalid value for --maf");
           }
           if(options_d["max_maf"] < 0.0 || options_d["max_maf"] > 0.5){
               LOGGER.e(0, "the value specified for --max-maf can't be smaller than 0 or larger than 0.5");
           }
        }else{
            LOGGER.e(0, " GCTA does not support multiple values for --maf currently ");
        }
        options_in.erase("--max-maf");
    }

    if(options_d["min_maf"] > options_d["max_maf"]){
        LOGGER.e(0, "value specified for --max-maf can't be smaller than that for --min-maf");
    }


    addOneValOption<double>("geno_rate", "--geno", options_in, options_d, 1.0, 0.0, 1.0);
    if(options_in.find("--geno") != options_in.end()){
        if(options_d["geno_rate"] == 0.0){
            options_in["--nofilter"] = {};
        }
    }

    addOneValOption<double>("info_score", "--info", options_in, options_d, 0.0, 0.0, 1.0);
    addOneValOption<double>("dos_dc", "--dc", options_in, options_d, -1.0, -1.0, 1.0);



    options_d["hard_call_thresh"] = 0.9;
    string flag = "--hard-call-thresh";
    if(options_in.find(flag) != options_in.end()){
        auto option = options_in[flag];
        if(options.size() == 1){
            try{
                options_d["hard_call_thresh"] = std::stod(option[0]);
            }catch(std::invalid_argument&){
                LOGGER.e(0, "invalid value in " + flag);
            }
        }else{
            LOGGER.e(0, " GCTA does not support multiple values for " + flag + " currently.");
        }
        options_in.erase(flag);
    }

    flag = "--dosage-call";
    if(options_in.find(flag) != options_in.end()){
        options["dosage_call"] = "true";
        options_in.erase(flag);
    }

    if(options_in.find("--freq") != options_in.end()){
        processFunctions.push_back("freq");
        if(options_in["--freq"].size() != 0){
            LOGGER.w(0, "--freq should not be followed by any parameter, if you want to calculate the allele frequencies in founders only, "
                    "please specify by --founders option");
        }
        options_in.erase("--freq");

        options["out"] = options_in["--out"][0];

        return_value++;
    }

    if(options_in.find("--freqx") != options_in.end()){
        processFunctions.push_back("freqx");
        if(options_in["--freqx"].size() != 0){
            LOGGER.w(0, "--freq should not be followed by any other parameter, if you want to calculate the allele frequencies in founders only, "
                    "please specify by --founders option");
        }
        options_in.erase("--freqx");

        options["out"] = options_in["--out"][0];

        return_value++;
    }

    if(options_in.find("--make-bed") != options_in.end()){
        if(options.find("bgen_file") == options.end()){
            processFunctions.push_back("make_bed");
        }else{
            processFunctions.push_back("make_bed_bgen");
        } 

        options_in.erase("--make-bed");
        options["out"] = options_in["--out"][0];

        return_value++;
    }

    if(options_in.find("--recodet") != options_in.end()){
        processFunctions.push_back("recodet");
        options["recode_method"] = "nomiss";
        if(options_in["--recodet"].size() > 0){
            string cop = options_in["--recodet"][0];
            vector<string> ops = {"nomiss", "raw", "std"};
            if(std::find(ops.begin(), ops.end(), cop) != ops.end()){
                options["recode_method"] = cop;
            }else{
                LOGGER.e(0, "can't recognize recode method: " + options_in["--recodet"][0]);
            }
        }
        options_in.erase("--recodet");
        options["out"] = options_in["--out"][0];
        return_value++;
    }


    addOneFileOption("update_freq_file", "", "--update-freq", options_in);

    if(options_in.find("--filter-sex") != options_in.end()){
        options["sex"] = "yes";
    }

    if(options_in.find("--sum-geno-x") != options_in.end()){
        processFunctions.push_back("sum_geno_x");
        options["sex"] = "yes";
        std::map<string, vector<string>> t_option;
        t_option["--chrx"] = {};
        t_option["--filter-sex"] = {};
        Pheno::registerOption(t_option);
        Marker::registerOption(t_option);
        options["out"] = options_in["--out"][0];
        return_value++;
    }



    return return_value;
}

void Geno::processFreq(){
    string name_out = options["out"] + ".frq";
    int buf_size = 23068672;
    osBuf.resize(buf_size);
    osOut.rdbuf()->pubsetbuf(&osBuf[0], buf_size);
 
    osOut.open(name_out.c_str());
    if (!osOut) { LOGGER.e(0, "cannot open the file [" + name_out + "] to write."); }
    osOut << "CHR\tSNP\tPOS\tA1\tA2\tAF\tNCHROBS";
    if(hasInfo){
        osOut << "\tINFO";
    }
    osOut << "\n";

    LOGGER << "Computing allele frequencies and saving them to [" << name_out << "]..." << std::endl;

    int nMarker = 128;
    vector<uint32_t> extractIndex(marker->count_extract());
    std::iota(extractIndex.begin(), extractIndex.end(), 0);
    
    vector<function<void (uintptr_t *, const vector<uint32_t> &)>> callBacks;
    callBacks.push_back(bind(&Geno::freq_func, this, _1, _2));

    numMarkerOutput = 0;
    loopDouble(extractIndex, nMarker, false, false, false, false, callBacks);

    osOut.flush();
    osOut.close();
    LOGGER << "Saved " << numMarkerOutput << " SNPs." << std::endl;
}

void Geno::processRecodet(){
    string name_out = options["out"] + ".xmat";
    int buf_size = 23068672;
    osBuf.resize(buf_size);
    osOut.rdbuf()->pubsetbuf(&osBuf[0], buf_size);
 
    osOut.open(name_out.c_str());
    if (!osOut) { LOGGER.e(0, "cannot open the file [" + name_out + "] to write."); }
    osOut << "CHR\tSNP\tPOS\tA1\tA2\tAF\tNCHROBS";
    if(hasInfo){
        osOut << "\tINFO";
    }

    LOGGER << "Recoding genotypes and saving them to [" << name_out << "]..." << std::endl;
    uint32_t n_sample = pheno->count_keep();
    vector<string> phenoID = pheno->get_id(0, n_sample - 1, "|");

    for(auto & phenItem : phenoID){
        osOut << "\t" << phenItem;
    }
    osOut << "\n";

    int nMarker = 128;
    bool center, std, saveMiss;
    if(options["recode_method"] == "std"){
        center = true;
        std = true;
        saveMiss = false;
    }else if(options["recode_method"] == "nomiss"){
        center = false;
        std = false;
        saveMiss = false;
    }else if(options["recode_method"] == "raw"){
        center = false;
        std = false;
        saveMiss = true;
    }
    bRecodeSaveMiss = saveMiss;

    vector<uint32_t> extractIndex(marker->count_extract());
    std::iota(extractIndex.begin(), extractIndex.end(), 0);
    
    vector<function<void (uintptr_t *, const vector<uint32_t> &)>> callBacks;
    callBacks.push_back(bind(&Geno::recode_func, this, _1, _2));

    numMarkerOutput = 0;
    loopDouble(extractIndex, nMarker, true, center, std, saveMiss, callBacks);

    osOut.flush();
    osOut.close();
    LOGGER << "Saved " << numMarkerOutput << " SNPs." << std::endl;

}

void Geno::freq_func(uintptr_t* genobuf, const vector<uint32_t> &markerIndex){
    int num_marker = markerIndex.size();
    vector<uint8_t> isValids(num_marker);
    vector<double> af(num_marker);
    vector<uint32_t> nValidAllele(num_marker);
    vector<double> info(num_marker);
    /*
    vector<string> outs(num_marker);
    int bufsize = keepSampleCT * 10;
    for(int i = 0; i < num_marker; i++){
        outs[i].reserve(bufsize);
    }
    */
    #pragma omp parallel for schedule(dynamic)
    for(int i = 0; i < num_marker; i++){
        uint32_t cur_marker = markerIndex[i];
        GenoBufItem item;
        item.extractedMarkerIndex = cur_marker;

        getGenoDouble(genobuf, i, &item);
        isValids[i] = item.valid;
        if(item.valid){
            af[i] = item.af;
            nValidAllele[i] = item.nValidAllele;
            info[i] = item.info;
        }
    }
    //output
    for(int i = 0; i != num_marker; i++){
        if(isValids[i]){
            numMarkerOutput++;
            osOut << marker->getMarkerStrExtract(markerIndex[i]) << "\t"  << af[i]
                <<"\t" << nValidAllele[i];
            if(hasInfo)osOut << "\t" << info[i];
            osOut << "\n";
        }
    }

}

void Geno::recode_func(uintptr_t* genobuf, const vector<uint32_t> &markerIndex){
    int num_marker = markerIndex.size();
    //vector<uint8_t> isValids(num_marker);
    /*
    vector<string> outs(num_marker);
    int bufsize = keepSampleCT * 10;
    for(int i = 0; i < num_marker; i++){
        outs[i].reserve(bufsize);
    }
    */
    #pragma omp parallel for ordered schedule(static,1)
    for(int i = 0; i < num_marker; i++){
        uint32_t cur_marker = markerIndex[i];
        GenoBufItem item;
        item.extractedMarkerIndex = cur_marker;

        getGenoDouble(genobuf, i, &item);

        #pragma omp ordered
        {
            if(item.valid) {
                numMarkerOutput++;
                osOut << marker->getMarkerStrExtract(cur_marker) << "\t"
                    << item.af << "\t" << item.nValidAllele;
                if(hasInfo) osOut << "\t" << item.info;
                for(int j = 0; j < keepSampleCT; j++){
                    osOut << "\t";
                    if(bRecodeSaveMiss && (item.missing[j/64] & (1UL << (j %64)))){
                        osOut << "NA";
                    }else{
                        osOut << item.geno[j];
                    }
                }
                osOut << "\n";
            }
        }
    }
}



void Geno::processMain() {
    for(auto &process_function : processFunctions){
        if(process_function == "freq"){
            Pheno pheno;
            Marker marker;
            Geno geno(&pheno, &marker);
            geno.processFreq();
        }

        if(process_function == "recodet"){
            Pheno pheno;
            Marker marker;
            Geno geno(&pheno, &marker);
            geno.processRecodet();
        }

        if(process_function == "make_bed"){
            LOGGER.e(0, "make_bed has not been migrated to the coroutine pipeline yet.");
        }

        if(process_function == "make_bed_bgen"){
            Pheno pheno;
            Marker marker;
            Geno geno(&pheno, &marker);
            string filename = options["out"];
            pheno.save_pheno(filename + ".fam");
            marker.save_marker(filename + ".bim");
            LOGGER.i(0, "Converting bgen to PLINK binary PED format [" + filename + ".bed]...");
            geno.bgen2bed(marker.get_extract_index());
            LOGGER.i(0, "Genotype has been saved.");
        }

        if(process_function == "sum_geno_x"){
            LOGGER.e(0, "sum_geno_x has not been migrated to the coroutine pipeline yet.");
        }
    }
}
