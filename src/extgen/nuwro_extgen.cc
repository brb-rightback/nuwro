// [uboone patch] NuWro request mode for MicroBooNE ExtGen (ubsim/ExtGen).
//
// ExtGen lets GENIE choose the neutrino, its energy, the target nucleus and
// the current (CC/NC) from the real flux and geometry, and asks an external
// generator for the final state of each such interaction ("request"). This
// program is NuWro's side of that: it serves requests one by one, each at its
// exact energy, with NuWro's own physics (makeevent + finishevent), and
// writes the generator-neutral ExtGen files. Stock NuWro is not changed: this
// is a separate executable (bin/nuwro_extgen); nuwro.h only grants it access
// (friend class ExtGenDriver).
//
// One process serves one group: a target nucleus, a neutrino flavour and a
// current, all fixed by the card and -p overrides (nucleus_p, nucleus_n,
// beam_particle, dyn_* of that current). Only the energy changes from request
// to request.
//
//   nuwro_extgen -m serve -r requests.txt -s served_events.txt -t trial_stats.txt
//                -i card.txt -o PREFIX [-p "param=value"] ...
//   nuwro_extgen -m scan -n TRIALS -r energies.txt -t scan.txt
//                -i card.txt -o PREFIX [-p "param=value"] ...
//
// -i, -o and -p are NuWro's own options (PREFIX.par lists the parameters).
// random_seed must be > 1 (set it with -p); Pythia6 is seeded from it too.
//
// Serve mode. Request file lines: "req_id E_GeV M_1 ... M_K", one Max per
// active channel, in the order of the header line
//   # extgen-nuwro-requests v1 channels=<dyn>,<dyn>,...
// which must equal this process's active channels (dyn codes of NuWro:
// 0/1 QEL, 2/3 RES, 4/5 DIS, 6/7 COH, 8/9 MEC, 10 HYP; even CC, odd NC).
// Max values in 1e-38 cm^2 per nucleon. Each trial draws channel k with
// probability M_k / sum(M), makes a NuWro event of channel k at the
// request's energy (weight w = dsigma estimate) and accepts it with
// probability w / M_k; a rejected trial starts again with a new channel.
// Accepted events are therefore distributed as the cross section summed over
// channels, exactly, as long as w <= M_k always. A weight above its Max
// stops the program (exit code 3), so that the caller retries with larger
// Max values; a channel with M_k = 0 is never generated.
// The trial weight u = w * sum(M) / M_k has mean sigma_total (the sum of the
// channels), which the in-job check (V4) compares with the spline table.
//
// Scan mode. Request lines "req_id E_GeV" (more columns are ignored). For
// each request and each active channel, TRIALS trials (no FSI): mean weight
// (sigma_k), sum of squares, largest weight. Feeds the spline and Max tables.
//
// Output (ExtGen formats, ubsim/ExtGen ExtGenIO.h):
//   served_events.txt  EVENT req_id E_GeV dyn simb_mode simb_inttype nparts
//                      P status pdg px py pz E   (GeV, target rest frame, nu along +z)
//                      status 11: struck nucleon (QEL, RES, DIS, HYP: NuWro's
//                      in[1]; MEC 2p2h: the nucleon pair in[1], in[2] as one cluster,
//                      pdg 2000000200 + number of protons, as GENIE; none for
//                      COH and for 3-nucleon MEC); status 1: every particle
//                      leaving the nucleus (post), the primary lepton first.
//                      gen_proc = NuWro dyn, + 100 for 3-nucleon MEC (3p3h).
//   trial_stats.txt    req_id tgt_pdg nu_pdg ccnc E_GeV n_trials sum_u sum_u2 max_u
//                      n_dropped   (u per nucleus, 1e-38 cm^2)
//   scan file          req_id E_GeV dyn n_trials sum_w sum_w2 max_w n_nan
//                      (w per nucleon, 1e-38 cm^2; n_trials valid trials,
//                      n_nan discarded NaN trials on top)
// NaN weights (from numerical failures inside NuWro's models; seen in the
// hybrid RES model for antineutrinos above ~1.5 GeV, ~0.3% of those trials)
// are discarded as stock NuWro's chooser does (Dyn::add refuses them): such
// a trial does not count in any statistic, so cross sections are means over
// valid trials, as NuWro reports them; the number is reported. An accepted
// event without its primary lepton after FSI is dropped (n_dropped) and the
// request continues.

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "nuwro.h"
#include "event1.h"
#include "generatormt.h"
#include "particle.h"
#include "beam.h"

NuWro nuwro;   // the stock main.cc defines it; the library code refers to it

extern "C" {
// Fortran stdout to /dev/null, as stock NuWro::main does first: Pythia6 prints
// an error block for each fragmentation retry in DIS (PYSTRF "caught in
// infinite loop", ~2 per DIS event at 3 GeV), which stock NuWro never shows
void shhpythiaitokay_(void);
// Pythia6 (hadronization in DIS) has its own RANMAR generator, which NuWro
// never seeds: every process would replay the same sequence. RANMAR takes
// seeds 0..900000000; MRPY(2) = 0 makes PYR re-initialize on its next call.
extern struct { int mrpy[6]; double rrpy[100]; } pydatr_;
}

namespace {

const char* kProg = "nuwro_extgen";
const double kUnit = 1e-38;   // cm^2: ExtGen's cross-section unit
const long kMaxTrials = 100000000L;   // per request: a closed channel set never ends

void die(const std::string& msg, int code = 2)
{
  std::cerr << kProg << ": ERROR: " << msg << std::endl;
  std::exit(code);
}

// a fixed-energy neutrino along +z (the requests' frame)
class FixedBeam : public beam
{
 public:
  FixedBeam(int pdg, double eMeV) : pdg_(pdg), e_(eMeV) {}
  particle shoot(bool) override
  {
    particle p(pdg_, 0.);
    p.set_momentum(vec(0, 0, 1));
    p.set_energy(e_);
    return p;
  }
 private:
  int pdg_;
  double e_;
};

struct Request {
  long id;
  double E;               // GeV, as read (written back unchanged)
  std::vector<double> M;  // serve: Max per active channel, 1e-38 cm^2 per nucleon
};

struct Channel {
  int dyn;
  std::string label;
};

int SimbMode(int dyn)
{
  switch (dyn) {
    case 0: case 1: case 10: return 0;   // kQE (hyperon production: QE, as GENIE)
    case 2: case 3: return 1;            // kRes (NuWro's single-pion region)
    case 4: case 5: return 2;            // kDIS
    case 6: case 7: return 3;            // kCoh
    case 8: case 9: return 10;           // kMEC
  }
  return -1;
}

// GENIE2ART interaction types (kNuanceOffset + NUANCE code); resonances and
// MEC 1000, as GENIE2ART writes outside its single-pion channels
int SimbIntType(int dyn)
{
  switch (dyn) {
    case 0: case 10: return 1001;   // kCCQE
    case 1: return 1002;            // kNCQE
    case 4: return 1091;            // kCCDIS
    case 5: return 1092;            // kNCDIS
    case 6: return 1097;            // kCCCOH
    case 7: return 1096;            // kNCCOH
  }
  return 1000;
}

struct ServedParticle {
  int status, pdg;
  double px, py, pz, E;   // GeV
};

}  // namespace

class ExtGenDriver
{
 public:
  explicit ExtGenDriver(NuWro& n) : nw(n) {}

  void Setup()
  {
    params& p = nw.p;
    if (p.random_seed <= 1) die("random_seed must be > 1 (set it with -p random_seed=N)");
    pydatr_.mrpy[0] = p.random_seed % 900000000;
    pydatr_.mrpy[1] = 0;
    if (p.target_type != 0) die("target_type must be 0 (one nucleus per process)");
    if (p.beam_test_only) die("beam_test_only must be 0");
    int pdg = p.beam_particle;
    if (!(std::abs(pdg) == 12 || std::abs(pdg) == 14 || std::abs(pdg) == 16))
      die("beam_particle must be a neutrino");
    if (p.dyn_lep) die("dyn_lep must be 0 (neutrino-electron scattering has no nucleus target)");
    if (p.dyn_qel_el || p.dyn_res_el) die("electron channels must be off");
    const bool cc = p.dyn_qel_cc || p.dyn_res_cc || p.dyn_dis_cc || p.dyn_coh_cc ||
                    p.dyn_mec_cc || p.dyn_hyp_cc;
    const bool nc = p.dyn_qel_nc || p.dyn_res_nc || p.dyn_dis_nc || p.dyn_coh_nc || p.dyn_mec_nc;
    if (cc == nc) die("enable the channels of exactly one current (dyn_*_cc or dyn_*_nc)");
    isCC = cc;
    const struct { bool on; int dyn; const char* label; } all[] = {
      {p.dyn_qel_cc, 0, "QELcc"}, {p.dyn_qel_nc, 1, "QELnc"}, {p.dyn_res_cc, 2, "REScc"},
      {p.dyn_res_nc, 3, "RESnc"}, {p.dyn_dis_cc, 4, "DIScc"}, {p.dyn_dis_nc, 5, "DISnc"},
      {p.dyn_coh_cc, 6, "COHcc"}, {p.dyn_coh_nc, 7, "COHnc"}, {p.dyn_mec_cc, 8, "MECcc"},
      {p.dyn_mec_nc, 9, "MECnc"}, {p.dyn_hyp_cc, 10, "HYPcc"}};
    for (auto const& c : all)
      if (c.on) channels.push_back({c.dyn, c.label});
    Z = p.nucleus_p;
    A = p.nucleus_p + p.nucleus_n;
    if (A < 1) die("no target nucleus (nucleus_p + nucleus_n = 0)");
    // sfevent (QEL with sf_method > 0) has spectral functions for C12, O16,
    // Ar40 and Fe56 only, and returns weight 0 for any other nucleus, silently
    const int N = p.nucleus_n;
    const bool sfNucleus = (Z == 6 && N == 6) || (Z == 8 && N == 8) || (Z == 18 && N == 22) ||
                           (Z == 26 && N == 30);
    if ((p.dyn_qel_cc || p.dyn_qel_nc) && p.sf_method > 0 && !sfNucleus)
      die("sf_method " + std::to_string((int)p.sf_method) + " has no spectral function for Z=" +
          std::to_string(Z) + " N=" + std::to_string(N) +
          " (QEL would have weight 0): use sf_method=0 for this nucleus");
    tgtPdg = 1000000000 + Z * 10000 + A * 10;
    nuPdg = pdg;
    const int sign = pdg > 0 ? 1 : -1;
    leptonPdg = isCC ? pdg - sign : pdg;
  }

  std::string ChannelList() const
  {
    std::string s;
    for (size_t k = 0; k < channels.size(); ++k)
      s += (k ? "," : "") + std::to_string(channels[k].dyn);
    return s;
  }

  std::vector<Request> ReadRequests(const std::string& path, bool serve)
  {
    std::ifstream in(path.c_str());
    if (!in) die("cannot read " + path);
    std::vector<Request> reqs;
    std::string line;
    bool header = false;
    long lineNo = 0;
    while (std::getline(in, line)) {
      ++lineNo;
      const std::string where = path + ":" + std::to_string(lineNo);
      if (line.compare(0, 30, "# extgen-nuwro-requests v1 cha") == 0) {
        const size_t pos = line.find("channels=");
        std::string list = line.substr(pos + 9);
        list = list.substr(0, list.find_first_of(" \t"));
        if (list != ChannelList())
          die(where + ": channels " + list + ", this process has " + ChannelList());
        header = true;
        continue;
      }
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      Request r;
      if (!(ss >> r.id >> r.E)) die(where + ": expected req_id E_GeV");
      if (!(r.E > 0)) die(where + ": energy must be positive");
      if (serve) {
        double m;
        while (ss >> m) {
          if (!(m >= 0)) die(where + ": Max must be >= 0");
          r.M.push_back(m);
        }
        if (r.M.size() != channels.size())
          die(where + ": " + std::to_string(r.M.size()) + " Max values for " +
              std::to_string(channels.size()) + " channels");
        double sum = 0;
        for (double m : r.M) sum += m;
        if (!(sum > 0)) die(where + ": all Max values are 0");
      }
      reqs.push_back(r);
    }
    if (serve && !header) die(path + ": no '# extgen-nuwro-requests v1 channels=...' line");
    if (reqs.empty()) die(path + ": no requests");
    return reqs;
  }

  // a fresh event of channel dyn at energy E (GeV) through NuWro's makeevent;
  // returns the weight in 1e-38 cm^2 per nucleon (NaN as -1)
  double MakeEvent(event*& e, int dyn)
  {
    e = new event();
    e->dyn = dyn;
    nw.makeevent(e, nw.p);
    const double w = e->weight / kUnit;
    return std::isnan(w) ? -1. : w;
  }

  void SetEnergy(double eGeV)
  {
    delete nw._beam;
    nw._beam = new FixedBeam(nuPdg, eGeV * 1000.);   // NuWro works in MeV
  }

  void Scan(const std::vector<Request>& reqs, long trials, const std::string& outPath)
  {
    FILE* out = std::fopen(outPath.c_str(), "w");
    if (!out) die("cannot write " + outPath);
    std::fprintf(out, "# extgen-nuwro-scan v1 version=%s tgt_pdg=%d nu_pdg=%d ccnc=%s "
                      "channels=%s trials=%ld seed=%d\n",
                 VERSION, tgtPdg, nuPdg, isCC ? "CC" : "NC", ChannelList().c_str(), trials,
                 (int)nw.p.random_seed);
    std::fprintf(out, "# req_id E_GeV dyn n_trials sum_w sum_w2 max_w n_nan"
                      "   (w: 1e-38 cm^2 per nucleon)\n");
    for (auto const& r : reqs) {
      SetEnergy(r.E);
      for (auto const& c : channels) {
        double sw = 0, sw2 = 0, mw = 0;
        long nnan = 0;
        for (long i = 0; i < trials; ) {
          event* e = nullptr;
          const double w = MakeEvent(e, c.dyn);
          delete e;
          if (w < 0) {   // NaN: discarded, as stock NuWro does
            if (++nnan > trials) die("more NaN weights than trials in channel " +
                                     std::to_string(c.dyn), 4);
            continue;
          }
          ++i;
          sw += w; sw2 += w * w;
          if (w > mw) mw = w;
        }
        std::fprintf(out, "%ld %.17g %d %ld %.17g %.17g %.17g %ld\n", r.id, r.E, c.dyn, trials,
                     sw, sw2, mw, nnan);
      }
      std::fflush(out);
      std::cout << kProg << ": scanned request " << r.id << " (E " << r.E << " GeV)" << std::endl;
    }
    if (std::fclose(out) != 0) die("error writing " + outPath);
  }

  void Serve(const std::vector<Request>& reqs, const std::string& evPath,
             const std::string& statPath, const std::string& reqPath)
  {
    FILE* ev = std::fopen(evPath.c_str(), "w");
    if (!ev) die("cannot write " + evPath);
    FILE* ts = std::fopen(statPath.c_str(), "w");
    if (!ts) die("cannot write " + statPath);
    std::fprintf(ev, "# extgen-served v1 generator=NuWro version=%s requests=%s tgt_pdg=%d "
                     "nu_pdg=%d ccnc=%s channels=%s seed=%d\n",
                 VERSION, reqPath.c_str(), tgtPdg, nuPdg, isCC ? "CC" : "NC",
                 ChannelList().c_str(), (int)nw.p.random_seed);
    std::fprintf(ev, "# EVENT req_id E_GeV gen_proc simb_mode simb_inttype nparts\n");
    std::fprintf(ev, "# P status pdg px py pz E   (GeV, neutrino along +z)\n");
    std::fprintf(ev, "# gen_proc = NuWro dyn (+ 100: 3-nucleon MEC)\n");
    std::fprintf(ts, "# extgen-trial-stats v1 generator=NuWro\n");
    std::fprintf(ts, "# req_id tgt_pdg nu_pdg ccnc E_GeV n_trials sum_w sum_w2 max_w n_dropped"
                     "   (w: 1e-38 cm^2 per nucleus)\n");
    long nNan = 0;
    for (auto const& r : reqs) {
      SetEnergy(r.E);
      double msum = 0;
      for (double m : r.M) msum += m;
      long n = 0, nDropped = 0;
      double su = 0, su2 = 0, mu = 0;
      bool served = false;
      while (!served) {
        if (++n > kMaxTrials)
          die("request " + std::to_string(r.id) + ": no event in " + std::to_string(kMaxTrials) +
              " trials", 4);
        // channel k with probability M_k / sum(M); channels with M_k = 0 never
        // (the last one with M_k > 0 if rounding puts x at the very end)
        const double x = frandom() * msum;
        size_t k = r.M.size();
        double acc = 0;
        for (size_t j = 0; j < r.M.size(); ++j) {
          if (r.M[j] == 0) continue;
          k = j;
          acc += r.M[j];
          if (x < acc) break;
        }
        event* e = nullptr;
        const double w = MakeEvent(e, channels[k].dyn);
        if (w < 0) {   // NaN: discarded, as stock NuWro does (not a trial)
          delete e;
          --n;
          if (++nNan > kMaxTrials) die("too many NaN weights", 4);
          continue;
        }
        if (w > r.M[k]) {
          char buf[300];
          std::snprintf(buf, sizeof buf,
                        "OVERWEIGHT request %ld E %.17g GeV channel %d: weight %.6g > Max %.6g "
                        "(1e-38 cm^2 per nucleon)",
                        r.id, r.E, channels[k].dyn, w, r.M[k]);
          delete e;
          die(buf, 3);
        }
        const double u = w * msum / r.M[k] * A;   // per nucleus
        su += u; su2 += u * u;
        if (u > mu) mu = u;
        if (frandom() * r.M[k] < w) {
          nw.finishevent(e, nw.p);
          if (WriteEvent(ev, r, *e)) served = true;
          else ++nDropped;
        }
        delete e;
      }
      std::fprintf(ts, "%ld %d %d %s %.17g %ld %.17g %.17g %.17g %ld\n", r.id, tgtPdg, nuPdg,
                   isCC ? "CC" : "NC", r.E, n, su, su2, mu, nDropped);
    }
    if (std::fclose(ev) != 0) die("error writing " + evPath);
    if (std::fclose(ts) != 0) die("error writing " + statPath);
    std::cout << kProg << ": served " << reqs.size() << " requests";
    if (nNan) std::cout << " (" << nNan << " NaN weights discarded)";
    std::cout << std::endl;
  }

 private:
  static ServedParticle Convert(int status, int pdg, particle& p)
  {
    return {status, pdg, p.x / 1000., p.y / 1000., p.z / 1000., p.t / 1000.};
  }

  // false if the event has no primary lepton among the particles leaving the nucleus
  bool WriteEvent(FILE* ev, const Request& r, event& e)
  {
    std::vector<ServedParticle> parts;
    const int dyn = e.dyn;
    if (dyn == 8 || dyn == 9) {   // MEC: in = nu + the initial nucleons (makeevent adds none)
      if (e.in.size() == 3) {
        particle& a = e.in[1];
        particle& b = e.in[2];
        const int np = (a.pdg == 2212) + (b.pdg == 2212);
        parts.push_back({11, 2000000200 + np, (a.x + b.x) / 1000., (a.y + b.y) / 1000.,
                         (a.z + b.z) / 1000., (a.t + b.t) / 1000.});
      }
    } else if (dyn != 6 && dyn != 7 && e.in.size() > 1) {
      parts.push_back(Convert(11, e.in[1].pdg, e.in[1]));
    }
    int iLepton = -1;
    for (size_t i = 0; i < e.post.size(); ++i)
      if (e.post[i].pdg == leptonPdg) { iLepton = i; break; }
    if (iLepton < 0) return false;
    parts.push_back(Convert(1, leptonPdg, e.post[iLepton]));
    for (size_t i = 0; i < e.post.size(); ++i)
      if ((int)i != iLepton) parts.push_back(Convert(1, e.post[i].pdg, e.post[i]));
    // gen_proc: NuWro's dyn, + 100 for 3-nucleon MEC (Valencia 3p3h)
    const int genProc = ((dyn == 8 || dyn == 9) && e.in.size() == 4) ? dyn + 100 : dyn;
    std::fprintf(ev, "EVENT %ld %.17g %d %d %d %zu\n", r.id, r.E, genProc, SimbMode(dyn),
                 SimbIntType(dyn), parts.size());
    for (auto const& q : parts)
      std::fprintf(ev, "P %d %d %.17g %.17g %.17g %.17g\n", q.status, q.pdg, q.px, q.py, q.pz, q.E);
    return true;
  }

  NuWro& nw;
  std::vector<Channel> channels;
  bool isCC = true;
  int Z = 0, A = 0, tgtPdg = 0, nuPdg = 0, leptonPdg = 0;
};

int main(int argc, char** argv)
{
  std::string mode, reqPath, evPath, statPath;
  long trials = 0;
  std::vector<char*> nwArgs{argv[0]};   // passed on to NuWro::init (-i, -o, -p)
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    if (i + 1 >= argc) die("option " + a + " needs a value");
    if (a == "-m") mode = argv[++i];
    else if (a == "-r") reqPath = argv[++i];
    else if (a == "-s") evPath = argv[++i];
    else if (a == "-t") statPath = argv[++i];
    else if (a == "-n") trials = std::atol(argv[++i]);
    else if (a == "-i" || a == "-o" || a == "-p") {
      nwArgs.push_back(argv[i]);
      nwArgs.push_back(argv[++i]);
    }
    else die("unknown option " + a + " (see src/extgen/nuwro_extgen.cc)");
  }
  if (mode != "serve" && mode != "scan") die("-m serve or -m scan");
  if (reqPath.empty() || statPath.empty()) die("give -r and -t");
  if (mode == "serve" && evPath.empty()) die("serve mode needs -s");
  if (mode == "scan" && trials < 1) die("scan mode needs -n TRIALS");

  try {
    shhpythiaitokay_();
    nuwro.init(nwArgs.size(), nwArgs.data());
    ExtGenDriver d(nuwro);
    d.Setup();
    const auto reqs = d.ReadRequests(reqPath, mode == "serve");
    if (mode == "scan") d.Scan(reqs, trials, statPath);
    else d.Serve(reqs, evPath, statPath, reqPath);
  }
  catch (std::string const& s) { die("NuWro: " + s, 4); }
  catch (char const* s) { die(std::string("NuWro: ") + s, 4); }
  catch (std::exception const& x) { die(x.what(), 4); }
  return 0;
}
