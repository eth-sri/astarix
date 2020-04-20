#include "utils.h"

namespace astarix {

void process_mem_usage(double& vm_usage, double& resident_set) {
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entrees in stat that we don't care about
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}

double b2gb(size_t bytes) {
    return bytes / 1024.0 / 1024.0 / 1024.0;
}

const char* edgeType2str(EdgeType type) {
    switch(type) {
        case ORIG:  return "ORIG";
        case INS:   return "ins";
        case DEL:   return "del";
        case SUBST: return "subst";
        case JUMP:  return "jump";
        case EdgeType_after_type: assert(false);
    }
    return "?";
}

char compl_nucl(char c) {
    if (c == 'A') return 'T';
    if (c == 'C') return 'G';
    if (c == 'G') return 'C';
    if (c == 'T') return 'A';
    std::cerr << "Bad nucleotide '" << c << "'" << std::endl;
    assert(false);
    return '!';
}

int nucl2num(char c) {
    if (c == 'A') return 0;
    if (c == 'C') return 1;
    if (c == 'G') return 2;
    if (c == 'T') return 3;
    assert(false);
    return '!';
}

bool is_nucl(char c) {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N';
}

std::string to_lower(std::string s) {
    for (char &c: s)
        if (c >= 'A' && c <= 'Z')
            c -= 'A'-'a';
    return s;
}

bool is_letter(char c) {
    return is_nucl(c) || c == EPS;
}

bool is_extended_nucl(char nucl) {
    return extended_nucls.find(nucl) != std::string::npos;
}

bool are_all_nucls(const std::string &s) {
    for (char c: s) {
        if (!is_nucl(toupper(c))) {
            std::cerr << "Not a nucleotide: [" << c << "]" << std::endl;
            return false;
        }
    }
    return true;
}

void write(std::ostream& os, int from, const edge_t &e) {
    os << from << " " << e.to << " " << e.label << " " << edgeType2str(e.type); 
}

std::string extended2orignucls(char nucl) {
    if (nucl == 'R') return "AG";
    if (nucl == 'Y') return "CT";
    if (nucl == 'K') return "GT";
    if (nucl == 'M') return "AC";
    if (nucl == 'S') return "CG";
    if (nucl == 'W') return "AT";
    if (nucl == 'B') return "CGT";
    if (nucl == 'D') return "AGT";
    if (nucl == 'H') return "ACT";
    if (nucl == 'V') return "ACG";
    if (nucl == 'N') return "ACGT";
    assert(is_nucl(nucl));
    return "!";
}

double sample() {
    return 1.0 * rand() / INT_MAX;
}

bool hasEnding(const std::string &fullString, const std::string &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

}
