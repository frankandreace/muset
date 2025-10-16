#ifndef AGGREGATOR_H
#define AGGREGATOR_H

#include <vector>
#include <stdint.h>
#include <iostream>

#include <fmt/format.h>

#include <kmat_tools/utils.h>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

// This interface is used to abstract the possible aggregating statistics for per-sample k-mer count in unitigs
class Aggregator {
public:

    struct AbundanceFractions {
    std::vector<double> abundances;
    std::vector<double> fractions;
    };

    virtual ~Aggregator() = default;
    // Process the counts for a single k-mer belonging to a unitig
    virtual void process_kmer(size_t bucket_id, const std::vector<uint32_t>& kmer_counts) = 0;
    // Calculate the final abundance for a given unitig and sample
    virtual std::pair<double, double> get_abundance_fraction(size_t bucket_id, size_t sample_id) const = 0;
    virtual AbundanceFractions get_all_abundance_fractions(size_t bucket_id) const = 0;
    virtual AbundanceFractions get_connected_component_abundance_fraction(const std::vector<size_t>& bucket_ids) const = 0;
};

// Reimplementation of Riccardo's mean computation logic
class MeanAggregator: public Aggregator {
  private:
  std::vector<std::vector<std::pair<uint32_t,uint32_t>>> m_samples_count;
  std::vector<size_t> m_kmer_counts;
  size_t m_num_samples;
  double m_min_fraction;
  // the pair contains <nb_present, abundance_sum>

  public:

  MeanAggregator(size_t num_samples, size_t num_buckets, double min_fraction):
  m_num_samples(num_samples), m_min_fraction(min_fraction) {
    m_samples_count.resize(num_buckets);
    for(size_t i {0}; i < num_buckets; i++){
      m_samples_count[i].resize(m_num_samples, {0,0});
    }
    m_kmer_counts.resize(num_buckets);
  }

  ~MeanAggregator() override = default;

  void process_kmer(size_t bucket_id, const std::vector<uint32_t>& kmer_counts) override {
    if (bucket_id >= m_samples_count.size()) {
      spdlog::debug(fmt::format("ERROR: GOT UTG ID {} BUT MAX POSSIBLE IS {}", bucket_id, m_samples_count.size()));
      return;
    }
    if (kmer_counts.size() != m_num_samples) {
      spdlog::debug(fmt::format("ERROR: GOT THESE NUMBER OF SAMPLES {}. EXPECTED {}", kmer_counts.size(), m_num_samples));
        return;
    }

    auto& samples = m_samples_count[bucket_id];
    size_t idx{0};
    for(auto& [nb_present,abundance_sum]: samples) {
        uint32_t num = kmer_counts[idx];
        // std::cout << "[INSIDE-" << bucket_id << "]: I see " << num << " and report " << uint32_t{num > 0} << ". old sum: " << nb_present << "\t";
        nb_present = kmat::add_sat(nb_present, uint32_t{num > 0});
        // std::cout << "New sum:" << nb_present << std::endl;
        abundance_sum = kmat::add_sat(abundance_sum, num);
        idx++;
    }
    m_kmer_counts[bucket_id]++;
  }

  std::pair<double,double> get_abundance_fraction(size_t bucket_id, size_t sample_id) const override{
    const size_t unitig_num_kmers = m_kmer_counts[bucket_id];
    if (bucket_id >= m_samples_count.size() || sample_id >= m_num_samples || unitig_num_kmers == 0) {
            return std::pair(0.0,0.0);
        }

    const auto& [nb_present, abundance_sum] = m_samples_count[bucket_id][sample_id];
    const double inv_kmers = 1.0 / unitig_num_kmers;
    // std::cout << "[INSIDE] GOT " << nb_present << " present out of " << unitig_num_kmers << std::endl;
    // return std::pair(static_cast<double>(abundance_sum) / unitig_num_kmers,static_cast<double>(nb_present)/ unitig_num_kmers);
    double fraction = static_cast<double>(nb_present) * inv_kmers;
    double abundance = fraction >= m_min_fraction ? static_cast<double>(abundance_sum) * inv_kmers : 0.0;
    return std::pair(abundance, fraction);
  }

  AbundanceFractions get_all_abundance_fractions(size_t bucket_id) const override {
    const size_t unitig_num_kmers = m_kmer_counts[bucket_id];
    AbundanceFractions result;
    result.abundances.resize(m_num_samples);
    result.fractions.resize(m_num_samples);
    if (bucket_id >= m_samples_count.size() || unitig_num_kmers == 0) {
      return result;
    }

    const auto& samples = m_samples_count[bucket_id];
    const double inv_kmers = 1.0 / unitig_num_kmers;

    for(size_t i = 0; i < m_num_samples; i++) {
      const auto& [nb_present, abundance_sum] = samples[i];
      result.fractions[i] = nb_present * inv_kmers;
      result.abundances[i] = result.fractions[i] >= m_min_fraction ? abundance_sum * inv_kmers : 0.0;
    }

    return result;
  }

  virtual AbundanceFractions get_connected_component_abundance_fraction(const std::vector<size_t>& bucket_ids) const override {
    AbundanceFractions result;
    size_t unitig_num_kmers {0};

    //checking for errors and computing total number k-mers
    for (const auto& bucket_id : bucket_ids){
      unitig_num_kmers += m_kmer_counts[bucket_id];
      if (bucket_id >= m_samples_count.size()) {
      return result;
      }
    }
    if (unitig_num_kmers == 0 || bucket_ids.size() == 0) return result;
    result.abundances.resize(m_num_samples);
    result.fractions.resize(m_num_samples);

    //iteration over all buckets
    for (const auto& bucket_id: bucket_ids){
      // summing fraction and abundance
      const auto& samples = m_samples_count[bucket_id];
      for(size_t i = 0; i < m_num_samples; i++) {
        const auto& [nb_present, abundance_sum] = samples[i];
        result.fractions[i] += nb_present;
        result.abundances[i] += abundance_sum;
      }
    }
    // dividing by total number of k-mer
    const double inv_kmers = 1.0 / unitig_num_kmers;
    for(size_t i = 0; i < m_num_samples; i++) {
        result.fractions[i] *= inv_kmers;
        result.abundances[i] = result.fractions[i] >= m_min_fraction ?  result.abundances[i] * inv_kmers : 0.0;
    }
    return result;
  }

};


// Median computation logic
class MedianAggregator: public Aggregator {
  private:
  std::vector<std::vector<std::map<uint32_t, uint32_t>>> m_samples_count;
  size_t m_num_samples;
  double m_min_fraction;
  std::vector<size_t> m_kmer_counts;
  // trying to reduce memory footprint, instead of having a vector of counts, I store the possible counts in a dictionary
  // I hope that most of the counts per sample will be the same so that I can save space.

  public:

  MedianAggregator(size_t num_samples, size_t num_buckets, double min_fraction):
  m_num_samples(num_samples), m_min_fraction(min_fraction) {
    m_samples_count.resize(num_buckets);
    for(int i {0}; i < num_buckets; i++){
      m_samples_count[i].resize(m_num_samples);
    }
    m_kmer_counts.resize(num_buckets);
  }

  ~MedianAggregator() override = default;

  void process_kmer(size_t bucket_id, const std::vector<uint32_t>& kmer_counts) override {
    if (bucket_id >= m_samples_count.size()) {
      spdlog::debug(fmt::format("ERROR: GOT UTG ID {} BUT MAX POSSIBLE IS {}", bucket_id, m_samples_count.size()));
      return;
    }
    if (kmer_counts.size() != m_num_samples) {
      spdlog::debug(fmt::format("ERROR: GOT THESE NUMBER OF SAMPLES {}. EXPECTED {}", kmer_counts.size(), m_num_samples));
        return;
    }

    auto& samples = m_samples_count[bucket_id];
    for(size_t idx {0}; idx < samples.size(); idx++) {
        samples[idx][kmer_counts[idx]]++;
    }
    m_kmer_counts[bucket_id]++;
  }

  std::pair<double, double> get_abundance_fraction(size_t bucket_id, size_t sample_id) const override{
    const size_t unitig_num_kmers = m_kmer_counts[bucket_id];
    const auto& sample_map = m_samples_count[bucket_id][sample_id];

    // 1 - computing fraction by counting number of occurrences.
    size_t nb_present {0};
    for (const auto& [key, value] : sample_map) {
        if (key > 0) {
            nb_present += value;
        }
    }

    // 2 - computing abundance by finding the median value from the dictionoary of occurencies.
    // Since I am using a std::map, keys are already sorted in ascending order. I just need to traverse the map
    // and stop at the middle value(s). In case odd, take the value, in case even, take
    double abundance {0.0};

    if (nb_present != 0) {
      const size_t idx1 = (nb_present - 1) / 2;
      const size_t idx2 = nb_present / 2;

      size_t running = 0;
      uint32_t m1 = 0;
      uint32_t m2 = 0;
      bool have_m1 = false;

      // Traverse sorted map to find median value(s)
      for (const auto& [count, freq] : sample_map) {
          if (count == 0) continue;  // Skip zeros for median calculation

          size_t next_running = running + freq;

          if (!have_m1 && next_running > idx1) {
              m1 = count;
              have_m1 = true;
          }

          if (next_running > idx2) {
              m2 = count;
              break;
          }

          running = next_running;
      }

      abundance = (idx1 == idx2) ? static_cast<double>(m2)
                                  : (static_cast<double>(m1) + m2) * 0.5;
    }

    // returning
    const double fraction = static_cast<double>(nb_present) / static_cast<double>(unitig_num_kmers);
    if (fraction < m_min_fraction) abundance = 0.0;
    return {abundance, fraction};
  }

  AbundanceFractions get_all_abundance_fractions(size_t bucket_id) const override {
    const size_t unitig_num_kmers = m_kmer_counts[bucket_id];
    AbundanceFractions result;
    result.abundances.resize(m_num_samples);
    result.fractions.resize(m_num_samples);

    if (bucket_id >= m_samples_count.size() || unitig_num_kmers == 0) {
        return result;
    }

    const auto& all_samples = m_samples_count[bucket_id];
    const double inv_kmers = 1.0 / unitig_num_kmers;

    for(size_t sample_id = 0; sample_id < m_num_samples; sample_id++) {
        const auto& sample_map = all_samples[sample_id];

        // 1 - Count present k-mers (non-zero)
        size_t nb_present = 0;
        for (const auto& [key, value] : sample_map) {
            if (key > 0) {
                nb_present += value;
            }
        }

        // 2 - Compute median
        double abundance = 0.0;

        if (nb_present != 0) {
            const size_t idx1 = (nb_present - 1) / 2;
            const size_t idx2 = nb_present / 2;

            size_t running = 0;
            uint32_t m1 = 0;
            uint32_t m2 = 0;
            bool have_m1 = false;

            // Traverse sorted map to find median value(s)
            for (const auto& [count, freq] : sample_map) {
                if (count == 0) continue;  // Skip zeros for median calculation

                size_t next_running = running + freq;

                if (!have_m1 && next_running > idx1) {
                    m1 = count;
                    have_m1 = true;
                }

                if (next_running > idx2) {
                    m2 = count;
                    break;
                }

                running = next_running;
            }

            abundance = (idx1 == idx2) ? static_cast<double>(m2)
                                       : (static_cast<double>(m1) + m2) * 0.5;
        }

        // 3 - Store results
        result.fractions[sample_id] = nb_present * inv_kmers;
        result.abundances[sample_id] = (result.fractions[sample_id] >= m_min_fraction) ? abundance : 0.0;
    }

    return result;
  }

  virtual AbundanceFractions get_connected_component_abundance_fraction(const std::vector<size_t>& bucket_ids) const override {
    size_t unitig_num_kmers {0};
    AbundanceFractions result;

    //checking for errors and computing total number k-mers
    for (const auto& bucket_id : bucket_ids){
      unitig_num_kmers += m_kmer_counts[bucket_id];
      if (bucket_id >= m_samples_count.size()) {
      return result;
      }
    }
    if (unitig_num_kmers == 0 || bucket_ids.size() == 0) return result;

    result.abundances.resize(m_num_samples);
    result.fractions.resize(m_num_samples);

    // Prepare merged maps for all samples
    std::vector<std::map<uint32_t, uint32_t>> merged_maps(m_num_samples);

    // CACHE-FRIENDLY: Iterate buckets outer, samples inner
    for (const auto& bucket_id : bucket_ids) {
        const auto& all_samples = m_samples_count[bucket_id];  // Sequential access coming

        for (size_t sample_id = 0; sample_id < m_num_samples; sample_id++) {
            const auto& sample_map = all_samples[sample_id];  // Contiguous!

            for (const auto& [count, freq] : sample_map) {
                merged_maps[sample_id][count] += freq;
            }
        }
    }

    // Compute statistics from merged maps
    const double inv_kmers = 1.0 / unitig_num_kmers;

    for (size_t sample_id = 0; sample_id < m_num_samples; sample_id++) {
        const auto& merged_map = merged_maps[sample_id];

        // Count present k-mers
        size_t nb_present = 0;
        for (const auto& [key, value] : merged_map) {
            if (key > 0) {
                nb_present += value;
            }
        }

        // Compute median
        double abundance = 0.0;

        if (nb_present != 0) {
            const size_t idx1 = (nb_present - 1) / 2;
            const size_t idx2 = nb_present / 2;

            size_t running = 0;
            uint32_t m1 = 0;
            uint32_t m2 = 0;
            bool have_m1 = false;

            for (const auto& [count, freq] : merged_map) {
                if (count == 0) continue;

                size_t next_running = running + freq;

                if (!have_m1 && next_running > idx1) {
                    m1 = count;
                    have_m1 = true;
                }

                if (next_running > idx2) {
                    m2 = count;
                    break;
                }

                running = next_running;
            }

            abundance = (idx1 == idx2) ? static_cast<double>(m2)
                                       : (static_cast<double>(m1) + m2) * 0.5;
        }

        result.fractions[sample_id] = nb_present * inv_kmers;
        result.abundances[sample_id] = (result.fractions[sample_id] >= m_min_fraction) ? abundance : 0.0;
    }

  return result;// Prepare merged maps for all samples
  }

};

#endif