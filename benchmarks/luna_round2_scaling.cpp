// Standalone, sparse length-scaling benchmark for the Luna round2 decoders.
//
// Compile this file in one of two modes:
//   -DLUNA_ROUND2_BEAM  and link BeamAlign.cpp
//   -DLUNA_ROUND2_NUSS  and link nussinov.cpp
//
// The Nussinov path deliberately calls only the sparse p/q overloads.  It
// never constructs a dense L-by-L probability matrix.  The benchmark is an
// empirical fixed-positive-beam experiment; it is not an asymptotic proof.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <sys/resource.h>

#if defined(LUNA_ROUND2_BEAM) && defined(LUNA_ROUND2_NUSS)
#error "select exactly one Luna round2 scaling implementation"
#elif defined(LUNA_ROUND2_BEAM)
#include "BeamAlign.h"
#include "nussinov.h"
#elif defined(LUNA_ROUND2_NUSS)
#include "nussinov.h"
#else
#error "define LUNA_ROUND2_BEAM or LUNA_ROUND2_NUSS"
#endif

namespace {

using Clock = std::chrono::steady_clock;

struct Hash64
{
  std::uint64_t value = 14695981039346656037ull;

  void bytes(const void* data, std::size_t size)
  {
    const auto* input = static_cast<const std::uint8_t*>(data);
    for (std::size_t i = 0; i < size; ++i) {
      value ^= input[i];
      value *= 1099511628211ull;
    }
  }

  template <typename T>
  void scalar(const T& item)
  {
    bytes(&item, sizeof(item));
  }

  void text(const std::string& item)
  {
    scalar<std::uint64_t>(item.size());
    bytes(item.data(), item.size());
  }
};

std::string
hex_digest(const Hash64& hash)
{
  std::ostringstream stream;
  stream << std::hex << std::setw(16) << std::setfill('0') << hash.value;
  return stream.str();
}

std::uint32_t
float_bits(float value)
{
  std::uint32_t bits = 0;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

std::uint64_t
double_bits(double value)
{
  std::uint64_t bits = 0;
  static_assert(sizeof(bits) == sizeof(value));
  std::memcpy(&bits, &value, sizeof(bits));
  return bits;
}

struct Usage
{
  double wall_seconds = 0.0;
  double user_seconds = 0.0;
  double system_seconds = 0.0;
  long max_rss_kb = 0;
};

double
time_value(const timeval& value)
{
  return static_cast<double>(value.tv_sec) +
         static_cast<double>(value.tv_usec) / 1000000.0;
}

Usage
usage_between(const Clock::time_point& start, const rusage& before,
              const rusage& after)
{
  Usage usage;
  usage.wall_seconds =
      std::chrono::duration<double>(Clock::now() - start).count();
  usage.user_seconds = time_value(after.ru_utime) - time_value(before.ru_utime);
  usage.system_seconds =
      time_value(after.ru_stime) - time_value(before.ru_stime);
  usage.max_rss_kb = after.ru_maxrss;
  return usage;
}

struct Options
{
  std::string variant;
  std::string pattern = "concentrated";
  std::vector<unsigned> lengths = {250, 500, 1000, 2000, 4000};
  unsigned beam = 100;
  unsigned repetitions = 3;
  unsigned warmup = 1;
  std::uint64_t seed = 42;
  std::string output;
};

std::vector<unsigned>
parse_lengths(const std::string& text)
{
  std::vector<unsigned> lengths;
  std::size_t begin = 0;
  while (begin < text.size()) {
    const std::size_t end = text.find(',', begin);
    const std::string token = text.substr(
        begin, end == std::string::npos ? std::string::npos : end - begin);
    if (token.empty())
      throw std::invalid_argument("empty length in --lengths");
    const unsigned length = static_cast<unsigned>(std::stoul(token));
    if (length < 3)
      throw std::invalid_argument("length must be at least 3");
    lengths.push_back(length);
    if (end == std::string::npos)
      break;
    begin = end + 1;
  }
  if (lengths.empty())
    throw std::invalid_argument("--lengths is empty");
  return lengths;
}

void
usage(const char* program)
{
  std::cerr
      << "usage: " << program
      << " --variant beam|nuss-certified|nuss-uncertified"
      << " [--lengths 250,500,1000,2000,4000] [--beam 100]"
      << " [--repetitions 3] [--warmup 1] [--seed 42]"
      << " [--pattern concentrated] [--output FILE]\n";
}

Options
parse_options(int argc, char** argv)
{
  Options options;
  for (int i = 1; i < argc; ++i) {
    const std::string argument = argv[i];
    const auto value = [&](const char* name) -> std::string {
      if (i + 1 >= argc)
        throw std::invalid_argument(std::string(name) + " requires a value");
      return argv[++i];
    };
    if (argument == "--variant")
      options.variant = value("--variant");
    else if (argument == "--lengths")
      options.lengths = parse_lengths(value("--lengths"));
    else if (argument == "--beam")
      options.beam = static_cast<unsigned>(std::stoul(value("--beam")));
    else if (argument == "--repetitions")
      options.repetitions =
          static_cast<unsigned>(std::stoul(value("--repetitions")));
    else if (argument == "--warmup")
      options.warmup = static_cast<unsigned>(std::stoul(value("--warmup")));
    else if (argument == "--seed")
      options.seed = std::stoull(value("--seed"));
    else if (argument == "--pattern")
      options.pattern = value("--pattern");
    else if (argument == "--output")
      options.output = value("--output");
    else if (argument == "--help" || argument == "-h") {
      usage(argv[0]);
      std::exit(0);
    } else {
      throw std::invalid_argument("unknown argument: " + argument);
    }
  }
  if (options.variant.empty())
    throw std::invalid_argument("--variant is required");
  if (options.beam == 0)
    throw std::invalid_argument("beam must be positive; beam=0 is excluded");
  if (options.repetitions == 0)
    throw std::invalid_argument("repetitions must be positive");
  if (options.pattern != "concentrated")
    throw std::invalid_argument("only --pattern concentrated is supported");
#if defined(LUNA_ROUND2_BEAM)
  if (options.variant != "beam")
    throw std::invalid_argument("Beam executable accepts only --variant beam");
#else
  if (options.variant != "nuss-certified" &&
      options.variant != "nuss-uncertified")
    throw std::invalid_argument(
        "Nussinov executable accepts only certified or uncertified variants");
#endif
  return options;
}

struct SparseEntry
{
  unsigned left = 0;
  unsigned right = 0;
  float value = 0.0f;
};

struct NussInput
{
  SparseFloatMatrix probability;
  SparseFloatMatrix multiplier;
  std::vector<std::vector<uint>> pairs_by_right;
  std::vector<SparseEntry> probability_entries;
  std::vector<SparseEntry> multiplier_entries;
  std::size_t raw_support = 0;
  std::size_t unique_support = 0;
  std::string input_digest;
};

std::uint64_t
pair_key(unsigned left, unsigned right)
{
  return (static_cast<std::uint64_t>(left) << 32) | right;
}

void
add_entry(SparseFloatMatrix& matrix, std::vector<SparseEntry>& entries,
          std::unordered_set<std::uint64_t>& seen, unsigned left,
          unsigned right, float value)
{
  const std::uint64_t key = pair_key(left, right);
  if (!seen.insert(key).second)
    return;
  matrix.set(left, right, value);
  entries.push_back({left, right, value});
}

NussInput
make_nuss_input(unsigned length, std::uint64_t seed, bool preassembled_support)
{
  NussInput input;
  input.probability.assign(length, length);
  input.multiplier.assign(length, length);
  std::unordered_set<std::uint64_t> probability_seen;
  std::unordered_set<std::uint64_t> multiplier_seen;
  probability_seen.reserve(static_cast<std::size_t>(length) * 3);
  multiplier_seen.reserve(static_cast<std::size_t>(length) * 2);

  // The sparse matrix wrapper accepts span >= 3 (left + 2 < right), while
  // the preassembled cached-support API accepts span >= 2.  Keep these paths
  // explicit so their support semantics are not accidentally conflated.
  const unsigned backbone_span = preassembled_support ? 2 : 3;
  for (unsigned right = backbone_span; right < length; ++right)
    add_entry(input.probability, input.probability_entries,
              probability_seen, right - backbone_span, right, 0.86f);

  // All additional probability edges share the final right endpoint.  This
  // is intentionally concentrated but still O(L) total support.
  const unsigned last_right_max_left = preassembled_support
      ? length - 3 : (length >= 4 ? length - 4 : 0);
  for (unsigned left = 0; left <= last_right_max_left; ++left) {
    const float fraction = static_cast<float>(left) /
                           static_cast<float>(std::max(1u, last_right_max_left));
    const float star_value = 0.78f - 0.20f * fraction;
    add_entry(input.probability, input.probability_entries,
              probability_seen, left, length - 1, star_value);
  }

  // Duplicate a deterministic subset in q so the decoder's support-union
  // path sees raw duplicates as well as distinct q-only support.
  for (unsigned right = backbone_span; right < length; right += 2)
    add_entry(input.multiplier, input.multiplier_entries,
              multiplier_seen, right - backbone_span, right, 0.04f);

  // q-only edges have a negative multiplier and therefore remain profitable
  // at p=0.  They occupy one other concentrated endpoint and are O(L).
  const unsigned q_right = length - 2;
  const unsigned q_limit = std::max(1u, length / 4);
  for (unsigned left = 0; left < q_limit && left + 2 < q_right; ++left) {
    const float fraction = static_cast<float>(left) /
                           static_cast<float>(std::max(1u, q_limit - 1));
    add_entry(input.multiplier, input.multiplier_entries,
              multiplier_seen, left, q_right, -0.35f + 0.02f * fraction);
  }

  std::map<std::pair<unsigned, unsigned>, std::pair<float, float>> support;
  for (const auto& entry : input.probability_entries)
    support[{entry.left, entry.right}].first = entry.value;
  for (const auto& entry : input.multiplier_entries)
    support[{entry.left, entry.right}].second = entry.value;
  input.raw_support = input.probability_entries.size() +
                      input.multiplier_entries.size();
  input.unique_support = support.size();
  input.pairs_by_right.assign(length, {});
  for (const auto& [pair, values] : support)
    input.pairs_by_right[pair.second].push_back(pair.first);

  Hash64 digest;
  digest.text("luna-round2-sparse-nuss-input-v1");
  digest.scalar(length);
  digest.scalar(seed);
  for (const auto& [pair, values] : support) {
    digest.scalar(pair.first);
    digest.scalar(pair.second);
    digest.scalar(float_bits(values.first));
    digest.scalar(float_bits(values.second));
  }
  input.input_digest = hex_digest(digest);
  return input;
}

struct Digests
{
  std::string score;
  std::string structure;
  std::string bounds;
};

struct Row
{
  std::string variant;
  std::string pattern;
  std::string support_mode;
  unsigned length = 0;
  unsigned beam = 0;
  unsigned warmup = 0;
  unsigned repetitions = 0;
  std::uint64_t seed = 0;
  std::size_t p_nnz = 0;
  std::size_t q_nnz = 0;
  std::size_t raw_support = 0;
  std::size_t unique_support = 0;
  std::string params_digest;
  std::string input_digest;
  std::string score;
  Digests digests;
  std::size_t structure_pairs = 0;
  float upper_bound = std::numeric_limits<float>::quiet_NaN();
  float pruned_upper_bound = std::numeric_limits<float>::quiet_NaN();
  float additive_upper_bound = std::numeric_limits<float>::quiet_NaN();
  std::size_t pruned_states = 0;
  double wall_seconds = 0.0;
  double user_seconds = 0.0;
  double system_seconds = 0.0;
  long max_rss_kb = 0;
  bool stable = true;
};

std::string
params_digest(const Options& options, unsigned length)
{
  Hash64 digest;
  digest.text("luna-round2-scaling-v1");
  digest.text(options.variant);
  digest.text(options.pattern);
  digest.scalar(length);
  digest.scalar(options.beam);
  digest.scalar(options.warmup);
  digest.scalar(options.repetitions);
  digest.scalar(options.seed);
  return hex_digest(digest);
}

Digests
beam_digests(double score, const std::vector<unsigned>& mapping)
{
  Hash64 score_hash;
  score_hash.text("score");
  score_hash.scalar(double_bits(score));
  Hash64 structure_hash;
  structure_hash.text("mapping");
  for (const unsigned value : mapping)
    structure_hash.scalar(value);
  Hash64 bounds_hash;
  bounds_hash.text("not-applicable");
  return {hex_digest(score_hash), hex_digest(structure_hash),
          hex_digest(bounds_hash)};
}

template <typename T>
bool
same_bits(T lhs, T rhs)
{
  return std::memcmp(&lhs, &rhs, sizeof(T)) == 0;
}

bool
same_nuss_result(const LinearNussinovResult& lhs,
                 const LinearNussinovResult& rhs)
{
  return same_bits(lhs.score, rhs.score) &&
         same_bits(lhs.upper_bound, rhs.upper_bound) &&
         same_bits(lhs.pruned_upper_bound, rhs.pruned_upper_bound) &&
         same_bits(lhs.additive_upper_bound, rhs.additive_upper_bound) &&
         lhs.pruned_states == rhs.pruned_states;
}

Digests
nuss_digests(const LinearNussinovResult& result, const VU& structure)
{
  Hash64 score_hash;
  score_hash.text("score");
  score_hash.scalar(float_bits(result.score));
  Hash64 structure_hash;
  structure_hash.text("structure");
  for (const uint value : structure)
    structure_hash.scalar(value);
  Hash64 bounds_hash;
  bounds_hash.text("bounds");
  bounds_hash.scalar(float_bits(result.upper_bound));
  bounds_hash.scalar(float_bits(result.pruned_upper_bound));
  bounds_hash.scalar(float_bits(result.additive_upper_bound));
  bounds_hash.scalar(result.pruned_states);
  return {hex_digest(score_hash), hex_digest(structure_hash),
          hex_digest(bounds_hash)};
}

std::size_t
count_pairs(const VU& structure)
{
  std::size_t pairs = 0;
  for (unsigned left = 0; left < structure.size(); ++left)
    if (structure[left] != static_cast<unsigned>(-1) &&
        structure[left] > left)
      ++pairs;
  return pairs;
}

void
write_csv_header(std::ostream& output)
{
  output
      << "schema_version,variant,pattern,support_mode,length,beam,warmup,repetitions,seed,"
      << "p_nnz,q_nnz,raw_support,unique_support,params_digest,input_digest,"
      << "score,score_digest,structure_pairs,structure_digest,upper_bound,"
      << "pruned_upper_bound,additive_upper_bound,pruned_states,bounds_digest,"
      << "wall_seconds,user_seconds,system_seconds,cpu_seconds,max_rss_kb,"
      << "output_stable,status\n";
}

void
write_optional_float(std::ostream& output, float value)
{
  if (std::isnan(value))
    return;
  output << std::setprecision(9) << value;
}

void
write_row(std::ostream& output, const Row& row)
{
  output << "1," << row.variant << ',' << row.pattern << ','
         << row.support_mode << ',' << row.length << ',' << row.beam << ','
         << row.warmup << ',' << row.repetitions
         << ',' << row.seed << ',' << row.p_nnz << ',' << row.q_nnz << ','
         << row.raw_support << ',' << row.unique_support << ','
         << row.params_digest << ',' << row.input_digest << ',' << std::setprecision(12)
         << row.score << ',' << row.digests.score << ',' << row.structure_pairs
         << ',' << row.digests.structure << ',';
  write_optional_float(output, row.upper_bound);
  output << ',';
  write_optional_float(output, row.pruned_upper_bound);
  output << ',';
  write_optional_float(output, row.additive_upper_bound);
  output << ',' << row.pruned_states << ',' << row.digests.bounds << ','
         << std::setprecision(12) << row.wall_seconds << ','
         << row.user_seconds << ',' << row.system_seconds << ','
         << row.user_seconds + row.system_seconds << ',' << row.max_rss_kb
         << ',' << (row.stable ? "true" : "false") << ",ok\n";
}

#if defined(LUNA_ROUND2_BEAM)

Row
run_beam(const Options& options, unsigned length)
{
  std::unordered_map<std::uint64_t, double> match_scores;
  match_scores.reserve(static_cast<std::size_t>(length) * 2);
  const std::uint64_t stride = static_cast<std::uint64_t>(length) + 1;
  for (unsigned i = 0; i < length; ++i)
    match_scores.emplace(static_cast<std::uint64_t>(i) * stride + i, 1.0);
  for (unsigned i = 0; i + 2 < length; ++i)
    match_scores.emplace(static_cast<std::uint64_t>(i) * stride + length - 1,
                         0.25);

  BeamAlign decoder(static_cast<int>(options.beam));
  std::vector<unsigned> mapping;
  const auto run_once = [&]() {
    return decoder.max_alignment(
        length, length, mapping,
        [&](int i, int k) {
          const auto it = match_scores.find(
              static_cast<std::uint64_t>(i) * stride + k);
          return it == match_scores.end() ? 0.0 : it->second;
        });
  };
  for (unsigned i = 0; i < options.warmup; ++i)
    run_once();

  double score = 0.0;
  Usage total;
  Digests first_digests;
  bool have_first = false;
  bool stable = true;
  for (unsigned i = 0; i < options.repetitions; ++i) {
    rusage before{};
    rusage after{};
    getrusage(RUSAGE_SELF, &before);
    const auto start = Clock::now();
    score = run_once();
    const Usage usage = usage_between(start, before,
                                      (getrusage(RUSAGE_SELF, &after), after));
    total.wall_seconds += usage.wall_seconds;
    total.user_seconds += usage.user_seconds;
    total.system_seconds += usage.system_seconds;
    total.max_rss_kb = std::max(total.max_rss_kb, usage.max_rss_kb);
    const Digests digests = beam_digests(score, mapping);
    if (!have_first) {
      first_digests = digests;
      have_first = true;
    } else if (digests.score != first_digests.score ||
               digests.structure != first_digests.structure ||
               digests.bounds != first_digests.bounds) {
      stable = false;
    }
  }

  Row row;
  row.variant = options.variant;
  row.pattern = options.pattern;
  row.support_mode = "beam-callback";
  row.length = length;
  row.beam = options.beam;
  row.warmup = options.warmup;
  row.repetitions = options.repetitions;
  row.seed = options.seed;
  row.raw_support = match_scores.size();
  row.unique_support = match_scores.size();
  row.params_digest = params_digest(options, length);
  Hash64 input_hash;
  input_hash.text("luna-round2-beam-input-v1");
  input_hash.scalar(length);
  input_hash.scalar(options.seed);
  std::vector<std::pair<std::uint64_t, double>> ordered_scores(
      match_scores.begin(), match_scores.end());
  std::sort(ordered_scores.begin(), ordered_scores.end());
  for (const auto& [key, value] : ordered_scores) {
    input_hash.scalar(key);
    input_hash.scalar(double_bits(value));
  }
  row.input_digest = hex_digest(input_hash);
  row.score = std::to_string(score);
  row.digests = beam_digests(score, mapping);
  row.structure_pairs = 0;
  row.wall_seconds = total.wall_seconds / options.repetitions;
  row.user_seconds = total.user_seconds / options.repetitions;
  row.system_seconds = total.system_seconds / options.repetitions;
  row.max_rss_kb = total.max_rss_kb;
  row.stable = stable;
  return row;
}

#else

Row
run_nuss(const Options& options, unsigned length)
{
  const bool preassembled_support = options.variant == "nuss-certified";
  NussInput input = make_nuss_input(length, options.seed,
                                    preassembled_support);
  LinearNussinov decoder(0.2f, options.beam);
  VU structure;
  LinearNussinovResult result;
  bool stable = true;
  Digests first_digests;
  bool have_first = false;

  const auto run_once = [&]() {
    if (options.variant == "nuss-certified")
      return decoder.decode_certified(1.0f, input.probability,
                                      input.multiplier, input.pairs_by_right,
                                      structure);
    const float score = decoder.decode(1.0f, input.probability,
                                       input.multiplier, structure);
    LinearNussinovResult uncertified;
    uncertified.score = score;
    uncertified.upper_bound = score;
    uncertified.pruned_upper_bound = score;
    uncertified.additive_upper_bound = score;
    return uncertified;
  };
  for (unsigned i = 0; i < options.warmup; ++i)
    result = run_once();

  Usage total;
  for (unsigned i = 0; i < options.repetitions; ++i) {
    rusage before{};
    rusage after{};
    getrusage(RUSAGE_SELF, &before);
    const auto start = Clock::now();
    result = run_once();
    const Usage usage = usage_between(start, before,
                                      (getrusage(RUSAGE_SELF, &after), after));
    total.wall_seconds += usage.wall_seconds;
    total.user_seconds += usage.user_seconds;
    total.system_seconds += usage.system_seconds;
    total.max_rss_kb = std::max(total.max_rss_kb, usage.max_rss_kb);
    const Digests digests = nuss_digests(result, structure);
    if (!have_first) {
      first_digests = digests;
      have_first = true;
    } else if (digests.score != first_digests.score ||
               digests.structure != first_digests.structure ||
               digests.bounds != first_digests.bounds) {
      stable = false;
    }
  }

  Row row;
  row.variant = options.variant;
  row.pattern = options.pattern;
  row.support_mode = preassembled_support
      ? "preassembled-cached-support" : "sparse-matrix-wrapper";
  row.length = length;
  row.beam = options.beam;
  row.warmup = options.warmup;
  row.repetitions = options.repetitions;
  row.seed = options.seed;
  row.p_nnz = input.probability_entries.size();
  row.q_nnz = input.multiplier_entries.size();
  row.raw_support = input.raw_support;
  row.unique_support = input.unique_support;
  row.params_digest = params_digest(options, length);
  row.input_digest = input.input_digest;
  row.score = std::to_string(result.score);
  row.digests = nuss_digests(result, structure);
  row.structure_pairs = count_pairs(structure);
  row.upper_bound = result.upper_bound;
  row.pruned_upper_bound = result.pruned_upper_bound;
  row.additive_upper_bound = result.additive_upper_bound;
  row.pruned_states = result.pruned_states;
  row.wall_seconds = total.wall_seconds / options.repetitions;
  row.user_seconds = total.user_seconds / options.repetitions;
  row.system_seconds = total.system_seconds / options.repetitions;
  row.max_rss_kb = total.max_rss_kb;
  row.stable = stable;
  return row;
}

#endif

} // namespace

#if defined(LUNA_ROUND2_NUSS)
// nussinov.cpp uses these decoder-base symbols but the standalone benchmark
// intentionally does not link the full folding stack.
const uint Fold::Decoder::n_support_brackets = 4 + 26;
const char* Fold::Decoder::left_brackets = "([{<ABCDEFGHIJKLMNOPQRSTUVWXYZ";
const char* Fold::Decoder::right_brackets = ")]}>abcdefghijklmnopqrstuvwxyz";
#endif

int
main(int argc, char** argv)
{
  try {
    const Options options = parse_options(argc, argv);
    std::ofstream file;
    std::ostream* output = &std::cout;
    if (!options.output.empty()) {
      file.open(options.output);
      if (!file)
        throw std::runtime_error("cannot open --output: " + options.output);
      output = &file;
    }
    write_csv_header(*output);
    for (const unsigned length : options.lengths) {
#if defined(LUNA_ROUND2_BEAM)
      write_row(*output, run_beam(options, length));
#else
      write_row(*output, run_nuss(options, length));
#endif
    }
    return 0;
  } catch (const std::exception& error) {
    std::cerr << "ERROR: " << error.what() << '\n';
    usage(argv[0]);
    return 2;
  }
}
