#include "../../include/skrf_cpp/io/metas.h"
#include "../../include/skrf_cpp/mathFunctions.h"
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace skrf_cpp {
namespace io {

void ns_to_sdatcv(const NetworkSet &ns, const std::string &fname, bool polar) {
    if(ns.size() == 0) throw std::runtime_error("empty NetworkSet");
    const Network &ntwk = ns.ntwk_set.front();
    int nports = ntwk.n_ports;
    // prepare freq column
    std::ofstream ofs(fname);
    if(!ofs) throw std::runtime_error("cannot open metas output: " + fname);
    // minimal header
    ofs << "SDATCV\nPorts\n";
    for(int k=0;k<nports;++k) ofs << (k+1) << (k+1==nports?"\n":"\t");
    // z0
    ofs << "Zr\n";
    // write frequencies and S (flattened mean)
    // compute mean network
    Network mean = ns.mean_s();
    size_t npts = mean.freqs.size();
    for(size_t i=0;i<npts;++i) {
        ofs << mean.freqs[i].hz;
        const auto &flat = mean.sparams[i];
        for(const auto &c : flat) ofs << '\t' << std::real(c) << '\t' << std::imag(c);
        ofs << '\n';
    }
    ofs.close();
}

} // namespace io
} // namespace skrf_cpp
