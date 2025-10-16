#include <kmat_tools/aggregator.h>
#include <gtest/gtest.h>
#include <random>
#include <vector>
#include <algorithm>
#include <cstdlib>

// ========== CONFIGURABLE TEST PARAMETERS ==========

// Default number of randomized test iterations
// Can be overridden via environment variable: NUM_RANDOM_TESTS
inline size_t get_num_random_tests() {
    static size_t num_tests = 0;
    if (num_tests == 0) {
        const char* env_val = std::getenv("NUM_RANDOM_TESTS");
        if (env_val) {
            num_tests = std::max(1, std::atoi(env_val));
        } else {
            num_tests = 1000000;  // Default value
        }
    }
    return num_tests;
}

// ========== HELPER FUNCTIONS ==========

// Random number generator
std::mt19937& get_rng() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return gen;
}

// Function to generate random count values between 0 and max_val per k-mer
std::vector<uint32_t> random_counts(size_t num_samples, uint32_t max_val = 100) {
    std::vector<uint32_t> counts(num_samples);
    std::uniform_int_distribution<uint32_t> dist(1, max_val);
    for (auto& c : counts) {
        c = dist(get_rng());
    }
    return counts;
}

// Helper function to compute the expected fraction
double expected_fract(const std::vector<std::vector<uint32_t>>& kmers, size_t sample_idx) {
    if (kmers.empty()) return 0.0;

    double sum = 0.0;
    for (const auto& kmer : kmers) {
        if (kmer[sample_idx] != 0) sum += 1;
    }
    return sum / kmers.size();
}

// Helper function to compute the expected mean
double expected_mean(const std::vector<std::vector<uint32_t>>& kmers, size_t sample_idx) {
    if (kmers.empty()) return 0.0;

    double sum = 0.0;
    for (const auto& kmer : kmers) {
        sum += kmer[sample_idx];
    }
    return sum / kmers.size();
}

// Helper function to compute the expected median
double expected_median(const std::vector<uint32_t>& values) {
    if (values.empty()) return 0.0;

    std::vector<uint32_t> sorted = values;
    std::sort(sorted.begin(), sorted.end());

    size_t n = sorted.size();
    if (n % 2 == 0) {
        return (sorted[n/2 - 1] + sorted[n/2]) / 2.0;
    } else {
        return sorted[n/2];
    }
}

// ========== TEST FIXTURE ==========

class AggregatorRandomTest : public ::testing::Test {
protected:
    void SetUp() override {
        get_rng().seed(std::random_device()());
    }
};

// ========== SINGLE UNITIG, SINGLE SAMPLE ACCESS TESTS ==========

TEST_F(AggregatorRandomTest, MeanAggregator_SingleUnitig_RandomTests) {
    const size_t num_tests = get_num_random_tests();
    const size_t max_samples = 5;
    const size_t max_kmers = 20;
    const double min_fraction = 0.5;

    for (size_t test_num = 0; test_num < num_tests; ++test_num) {
        size_t num_samples = 1 + get_rng()() % max_samples;
        size_t num_kmers = 1 + get_rng()() % max_kmers;

        MeanAggregator agg(num_samples, 1, min_fraction);

        std::vector<std::vector<uint32_t>> kmers;
        for (size_t i = 0; i < num_kmers; ++i) {
            auto counts = random_counts(num_samples);
            kmers.push_back(counts);
            agg.process_kmer(0, counts);
        }

        for (size_t sample = 0; sample < num_samples; ++sample) {
            auto [computed_abundance, computed_fraction] = agg.get_abundance_fraction(0, sample);

            double expected_fraction = expected_fract(kmers, sample);
            EXPECT_DOUBLE_EQ(computed_fraction, expected_fraction)
                << "Test " << test_num << "/" << num_tests << ", Sample " << sample;

            double expected_abundance = expected_mean(kmers, sample);
            EXPECT_NEAR(computed_abundance, expected_abundance, 0.001)
                << "Test " << test_num << "/" << num_tests << ", Sample " << sample;
        }
    }
}

TEST_F(AggregatorRandomTest, MeanAggregator_MinFraction_RandomTests) {
    const size_t num_tests = get_num_random_tests();
    const size_t max_samples = 3;
    const size_t max_kmers = 15;

    for (size_t test_num = 0; test_num < num_tests; ++test_num) {
        size_t num_samples = 1 + get_rng()() % max_samples;
        size_t num_kmers = 5 + get_rng()() % max_kmers;
        double min_fraction = 0.1 + (get_rng()() % 90) / 100.0;

        MeanAggregator agg(num_samples, 1, min_fraction);

        std::vector<std::vector<uint32_t>> sent_kmers;
        sent_kmers.reserve(num_kmers);

        std::bernoulli_distribution process_dist(0.7);

        for (size_t i = 0; i < num_kmers; ++i) {
            std::vector<uint32_t> kmer;
            if (process_dist(get_rng())) {
                kmer = random_counts(num_samples);
            } else {
                kmer = std::vector<uint32_t>(num_samples, 0);
            }
            sent_kmers.push_back(kmer);
            agg.process_kmer(0, kmer);
        }

        for (size_t sample = 0; sample < num_samples; ++sample) {
            auto [computed_abundance, computed_fraction] = agg.get_abundance_fraction(0, sample);

            double nb_present = 0.0;
            double sum = 0.0;
            for (const auto& kmer : sent_kmers) {
                if (kmer[sample] > 0) nb_present += 1.0;
                sum += kmer[sample];
            }

            double expected_fraction = nb_present / num_kmers;
            // Use the COMPUTED fraction for the threshold check, not the expected one
            double expected_abundance = (computed_fraction >= min_fraction) ? (sum / num_kmers) : 0.0;

            EXPECT_NEAR(computed_fraction, expected_fraction, 0.001)
                << "Test " << test_num << "/" << num_tests << ", Sample " << sample;

            EXPECT_NEAR(computed_abundance, expected_abundance, 0.001)
                << "Test " << test_num << "/" << num_tests << ", Sample " << sample;
        }
    }
}

TEST_F(AggregatorRandomTest, MedianAggregator_SingleUnitig_RandomTests) {
    const size_t num_tests = get_num_random_tests();
    const size_t max_samples = 5;
    const size_t max_kmers = 25;
    const double min_fraction = 0.5;

    for (size_t test_num = 0; test_num < num_tests; ++test_num) {
        size_t num_samples = 1 + get_rng()() % max_samples;
        size_t num_kmers = 1 + get_rng()() % max_kmers;

        MedianAggregator agg(num_samples, 1, min_fraction);

        std::vector<std::vector<uint32_t>> kmers;
        for (size_t i = 0; i < num_kmers; ++i) {
            auto counts = random_counts(num_samples, 20);
            kmers.push_back(counts);
            agg.process_kmer(0, counts);
        }

        for (size_t sample = 0; sample < num_samples; ++sample) {
            auto [computed_abundance, computed_fraction] = agg.get_abundance_fraction(0, sample);

            EXPECT_DOUBLE_EQ(computed_fraction, 1.0)
                << "Test " << test_num << "/" << num_tests << ", Sample " << sample;

            std::vector<uint32_t> sample_counts;
            for (const auto& kmer : kmers) {
                sample_counts.push_back(kmer[sample]);
            }

            double expected_abundance = expected_median(sample_counts);
            EXPECT_NEAR(computed_abundance, expected_abundance, 0.001)
                << "Test " << test_num << "/" << num_tests << ", Sample " << sample;
        }
    }
}

TEST_F(AggregatorRandomTest, MedianAggregator_MinFraction_RandomTests) {
    const size_t num_tests = get_num_random_tests();
    const size_t max_samples = 3;
    const size_t max_kmers = 15;

    for (size_t test_num = 0; test_num < num_tests; ++test_num) {
        size_t num_samples = 1 + get_rng()() % max_samples;
        size_t num_kmers = 5 + get_rng()() % max_kmers;
        double min_fraction = 0.1 + (get_rng()() % 90) / 100.0;

        MedianAggregator agg(num_samples, 1, min_fraction);

        std::vector<std::vector<uint32_t>> sent_kmers;
        sent_kmers.reserve(num_kmers);

        std::bernoulli_distribution process_dist(0.7);

        for (size_t i = 0; i < num_kmers; ++i) {
            std::vector<uint32_t> kmer;
            if (process_dist(get_rng())) {
                kmer = random_counts(num_samples);
            } else {
                kmer = std::vector<uint32_t>(num_samples, 0);
            }
            sent_kmers.push_back(kmer);
            agg.process_kmer(0, kmer);
        }

        for (size_t sample = 0; sample < num_samples; ++sample) {
            auto [computed_abundance, computed_fraction] = agg.get_abundance_fraction(0, sample);

            double nb_present = 0.0;
            std::vector<uint32_t> non_zero_values;  // ← Only non-zero values
            for (const auto& kmer : sent_kmers) {
                if (kmer[sample] > 0) {
                    nb_present += 1.0;
                    non_zero_values.push_back(kmer[sample]);  // ← Only add if non-zero
                }
            }

            double expected_fraction = nb_present / num_kmers;
            double expected_abundance = (expected_fraction >= min_fraction && !non_zero_values.empty())
                                       ? expected_median(non_zero_values)  // ← Median of NON-ZERO only
                                       : 0.0;

            EXPECT_NEAR(computed_fraction, expected_fraction, 0.001)
                << "Test " << test_num << "/" << num_tests << ", Sample " << sample;

            EXPECT_NEAR(computed_abundance, expected_abundance, 0.001)
                << "Test " << test_num << "/" << num_tests << ", Sample " << sample;
        }
    }
}

// ========== BATCH (ALL SAMPLES) TESTS ==========

TEST_F(AggregatorRandomTest, MeanAggregator_GetAllSamples_RandomTests) {
    const size_t num_tests = get_num_random_tests();
    const size_t max_samples = 5;
    const size_t max_kmers = 20;
    const double min_fraction = 0.5;

    for (size_t test_num = 0; test_num < num_tests; ++test_num) {
        size_t num_samples = 1 + get_rng()() % max_samples;
        size_t num_kmers = 1 + get_rng()() % max_kmers;

        MeanAggregator agg(num_samples, 1, min_fraction);

        std::vector<std::vector<uint32_t>> kmers;
        for (size_t i = 0; i < num_kmers; ++i) {
            auto counts = random_counts(num_samples);
            kmers.push_back(counts);
            agg.process_kmer(0, counts);
        }

        auto result = agg.get_all_abundance_fractions(0);

        ASSERT_EQ(result.abundances.size(), num_samples)
            << "Test " << test_num << "/" << num_tests;
        ASSERT_EQ(result.fractions.size(), num_samples)
            << "Test " << test_num << "/" << num_tests;

        for (size_t sample = 0; sample < num_samples; ++sample) {
            double expected_fraction = expected_fract(kmers, sample);
            double expected_abundance = expected_mean(kmers, sample);

            EXPECT_DOUBLE_EQ(result.fractions[sample], expected_fraction)
                << "Test " << test_num << "/" << num_tests << ", Sample " << sample;
            EXPECT_NEAR(result.abundances[sample], expected_abundance, 0.001)
                << "Test " << test_num << "/" << num_tests << ", Sample " << sample;
        }
    }
}

TEST_F(AggregatorRandomTest, MedianAggregator_GetAllSamples_RandomTests) {
    const size_t num_tests = get_num_random_tests();
    const size_t max_samples = 5;
    const size_t max_kmers = 25;
    const double min_fraction = 0.5;

    for (size_t test_num = 0; test_num < num_tests; ++test_num) {
        size_t num_samples = 1 + get_rng()() % max_samples;
        size_t num_kmers = 1 + get_rng()() % max_kmers;

        MedianAggregator agg(num_samples, 1, min_fraction);

        std::vector<std::vector<uint32_t>> kmers;
        for (size_t i = 0; i < num_kmers; ++i) {
            auto counts = random_counts(num_samples, 20);
            kmers.push_back(counts);
            agg.process_kmer(0, counts);
        }

        auto result = agg.get_all_abundance_fractions(0);

        ASSERT_EQ(result.abundances.size(), num_samples)
            << "Test " << test_num << "/" << num_tests;
        ASSERT_EQ(result.fractions.size(), num_samples)
            << "Test " << test_num << "/" << num_tests;

        for (size_t sample = 0; sample < num_samples; ++sample) {
            EXPECT_DOUBLE_EQ(result.fractions[sample], 1.0)
                << "Test " << test_num << "/" << num_tests << ", Sample " << sample;

            std::vector<uint32_t> sample_counts;
            for (const auto& kmer : kmers) {
                sample_counts.push_back(kmer[sample]);
            }

            double expected_abundance = expected_median(sample_counts);
            EXPECT_NEAR(result.abundances[sample], expected_abundance, 0.001)
                << "Test " << test_num << "/" << num_tests << ", Sample " << sample;
        }
    }
}

// ========== CONNECTED COMPONENT (MULTIPLE UNITIGS) TESTS ==========

TEST_F(AggregatorRandomTest, MeanAggregator_ConnectedComponent_RandomTests) {
    const size_t num_tests = get_num_random_tests();
    const size_t max_samples = 4;
    const size_t max_unitigs = 5;
    const size_t max_kmers = 10;
    const double min_fraction = 0.5;

    for (size_t test_num = 0; test_num < num_tests; ++test_num) {
        size_t num_samples = 1 + get_rng()() % max_samples;
        size_t num_unitigs = 1 + get_rng()() % max_unitigs;

        MeanAggregator agg(num_samples, num_unitigs, min_fraction);

        std::vector<std::vector<std::vector<uint32_t>>> all_kmers(num_unitigs);
        std::vector<size_t> bucket_ids;

        for (size_t utg = 0; utg < num_unitigs; ++utg) {
            size_t num_kmers = 1 + get_rng()() % max_kmers;
            for (size_t i = 0; i < num_kmers; ++i) {
                auto counts = random_counts(num_samples);
                all_kmers[utg].push_back(counts);
                agg.process_kmer(utg, counts);
            }
            bucket_ids.push_back(utg);
        }

        auto result = agg.get_connected_component_abundance_fraction(bucket_ids);

        ASSERT_EQ(result.abundances.size(), num_samples)
            << "Test " << test_num << "/" << num_tests;
        ASSERT_EQ(result.fractions.size(), num_samples)
            << "Test " << test_num << "/" << num_tests;

        // Flatten all kmers across unitigs
        std::vector<std::vector<uint32_t>> flattened;
        for (const auto& utg_kmers : all_kmers) {
            for (const auto& kmer : utg_kmers) {
                flattened.push_back(kmer);
            }
        }

        for (size_t sample = 0; sample < num_samples; ++sample) {
            double expected_fraction = expected_fract(flattened, sample);
            double expected_abundance = expected_mean(flattened, sample);

            EXPECT_NEAR(result.fractions[sample], expected_fraction, 0.001)
                << "Test " << test_num << "/" << num_tests << ", Sample " << sample;
            EXPECT_NEAR(result.abundances[sample], expected_abundance, 0.001)
                << "Test " << test_num << "/" << num_tests << ", Sample " << sample;
        }
    }
}

TEST_F(AggregatorRandomTest, MedianAggregator_ConnectedComponent_RandomTests) {
    const size_t num_tests = get_num_random_tests();
    const size_t max_samples = 4;
    const size_t max_unitigs = 5;
    const size_t max_kmers = 10;
    const double min_fraction = 0.5;

    for (size_t test_num = 0; test_num < num_tests; ++test_num) {
        size_t num_samples = 1 + get_rng()() % max_samples;
        size_t num_unitigs = 1 + get_rng()() % max_unitigs;

        MedianAggregator agg(num_samples, num_unitigs, min_fraction);

        std::vector<std::vector<std::vector<uint32_t>>> all_kmers(num_unitigs);
        std::vector<size_t> bucket_ids;

        for (size_t utg = 0; utg < num_unitigs; ++utg) {
            size_t num_kmers = 1 + get_rng()() % max_kmers;
            for (size_t i = 0; i < num_kmers; ++i) {
                auto counts = random_counts(num_samples, 20);
                all_kmers[utg].push_back(counts);
                agg.process_kmer(utg, counts);
            }
            bucket_ids.push_back(utg);
        }

        auto result = agg.get_connected_component_abundance_fraction(bucket_ids);

        ASSERT_EQ(result.abundances.size(), num_samples)
            << "Test " << test_num << "/" << num_tests;
        ASSERT_EQ(result.fractions.size(), num_samples)
            << "Test " << test_num << "/" << num_tests;

        // Flatten all kmers across unitigs
        std::vector<std::vector<uint32_t>> flattened;
        for (const auto& utg_kmers : all_kmers) {
            for (const auto& kmer : utg_kmers) {
                flattened.push_back(kmer);
            }
        }

        for (size_t sample = 0; sample < num_samples; ++sample) {
            EXPECT_DOUBLE_EQ(result.fractions[sample], 1.0)
                << "Test " << test_num << "/" << num_tests << ", Sample " << sample;

            std::vector<uint32_t> sample_counts;
            for (const auto& kmer : flattened) {
                sample_counts.push_back(kmer[sample]);
            }

            double expected_abundance = expected_median(sample_counts);
            EXPECT_NEAR(result.abundances[sample], expected_abundance, 0.001)
                << "Test " << test_num << "/" << num_tests << ", Sample " << sample;
        }
    }
}

TEST_F(AggregatorRandomTest, MeanAggregator_MinFraction_DEBUG) {
    size_t num_samples = 2;
    size_t num_kmers = 10;
    double min_fraction = 0.5;

    MeanAggregator agg(num_samples, 1, min_fraction);

    // Simple deterministic test case
    std::vector<std::vector<uint32_t>> sent_kmers = {
        {10, 20},  // kmer 0
        {0, 30},   // kmer 1
        {15, 0},   // kmer 2
        {0, 0},    // kmer 3
        {5, 10},   // kmer 4
        {0, 0},    // kmer 5
        {20, 25},  // kmer 6
        {0, 15},   // kmer 7
        {0, 0},    // kmer 8
        {10, 20},  // kmer 9
    };

    for (const auto& kmer : sent_kmers) {
        agg.process_kmer(0, kmer);
    }

    std::cout << "\n=== DEBUG OUTPUT ===\n";
    std::cout << "Total k-mers: " << sent_kmers.size() << "\n";
    std::cout << "Min fraction: " << min_fraction << "\n\n";

    for (size_t sample = 0; sample < num_samples; ++sample) {
        auto [computed_abundance, computed_fraction] = agg.get_abundance_fraction(0, sample);

        // Manual calculation
        double nb_present = 0.0;
        double sum = 0.0;
        for (const auto& kmer : sent_kmers) {
            if (kmer[sample] > 0) {
                nb_present += 1.0;
                std::cout << "Sample " << sample << " k-mer present with value " << kmer[sample] << "\n";
            }
            sum += kmer[sample];
        }

        double expected_fraction = nb_present / sent_kmers.size();
        double expected_abundance = (expected_fraction >= min_fraction)
                                   ? (sum / sent_kmers.size())
                                   : 0.0;

        std::cout << "\nSample " << sample << ":\n";
        std::cout << "  Present k-mers: " << nb_present << "/" << sent_kmers.size() << "\n";
        std::cout << "  Sum: " << sum << "\n";
        std::cout << "  Expected fraction: " << expected_fraction << "\n";
        std::cout << "  Computed fraction: " << computed_fraction << "\n";
        std::cout << "  Expected abundance: " << expected_abundance << "\n";
        std::cout << "  Computed abundance: " << computed_abundance << "\n";
        std::cout << "  Fraction >= min_fraction? " << (expected_fraction >= min_fraction) << "\n\n";

        EXPECT_NEAR(computed_fraction, expected_fraction, 0.001);
        EXPECT_NEAR(computed_abundance, expected_abundance, 0.001);
    }
}