#include "../../include/skrf_cpp/io/general.h"
#include "../../include/skrf_cpp/touchstone.h"
#include <filesystem>
#include <stdexcept>

namespace skrf_cpp {
namespace io {

void write_network(const std::string &path, const Network &net, bool overwrite) {
    std::string out = path;
    // if no extension, append .s2p
    if(out.find('.') == std::string::npos) out += ".s2p";
    if(!overwrite) {
        if(std::filesystem::exists(out)) return;
    }
    Touchstone::save(out, net);
}

Network read_network(const std::string &path) {
    // only support touchstone files for now
    return Touchstone::load(path);
}

std::vector<std::pair<std::string, Network>> read_all_networks(const std::string &dir) {
    std::vector<std::pair<std::string, Network>> out;
    for(auto &p : std::filesystem::directory_iterator(dir)) {
        if(!p.is_regular_file()) continue;
        auto name = p.path().filename().string();
        auto ext = p.path().extension().string();
        // match .s1p .s2p .s3p .sNp (case-insensitive)
        if(ext.size() >= 3 && (ext[0]=='.' || ext[0]=='.')) {
            std::string lext = ext;
            std::transform(lext.begin(), lext.end(), lext.begin(), [](unsigned char c){ return std::tolower(c); });
            if(lext.size()>=3 && lext[1]=='s') { /* skip */ }
        }
        // simple check: look for pattern like .s?p
        if(name.size() >= 3) {
            std::string lower = name;
            std::transform(lower.begin(), lower.end(), lower.begin(), [](unsigned char c){ return std::tolower(c); });
            if(lower.find(".s1p") != std::string::npos || lower.find(".s2p") != std::string::npos || lower.find(".s3p") != std::string::npos || lower.find(".s4p") != std::string::npos) {
                try {
                    Network net = Touchstone::load(p.path().string());
                    std::string key = p.path().stem().string();
                    out.emplace_back(key, std::move(net));
                } catch(...) {
                    // skip
                }
            }
        }
    }
    return out;
}

} // namespace io
} // namespace skrf_cpp
