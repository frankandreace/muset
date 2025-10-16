#include "kmat_tools/graph_reader.h"
#include <gtest/gtest.h>

#include <fstream>
#include <random>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <filesystem>
#include <cstring>
#include <unistd.h>
#include <sys/stat.h>
#include <queue>
#include <set>

namespace fs = std::filesystem;

int MAX_NUM_RANDOM_ITERATIONS {1000};

std::string createTempDirectory(const std::string& name) {
    fs::path temp_root = fs::temp_directory_path();
    fs::path temp_dir = temp_root / name;
    std::cerr << "Creating temp directory: " << temp_dir << std::endl;
    if (fs::exists(temp_dir)) {std::cerr << "Already exists" << std::endl; return temp_dir.string();}
    std::error_code ec;
    if (!fs::exists(temp_dir) && !fs::create_directory(temp_dir, ec)) {
        std::cerr << "Error creating temp directory: " << ec.message() << std::endl;
        throw std::runtime_error("Failed to create temporary directory");
    }

    return temp_dir.string();
}

std::string createTempFastaFile(const std::vector<std::pair<uint64_t, std::vector<std::pair<uint64_t, bool>>>>& graph) {
    std::string temp_dir_name = "muset_unit_tests";
    std::string dir_path = createTempDirectory(temp_dir_name);
    char template_name[] = "/tmp/muset_unit_tests/test_graph_XXXXXX";

    // Create temporary file using mkstemp
    int fd = mkstemp(template_name);
    if (fd == -1) {
        throw std::runtime_error("Failed to create temporary file");
    }

    // Use file descriptor to open a FILE* for writing
    FILE* file = fdopen(fd, "w");
    if (!file) {
        close(fd);
        unlink(template_name);  // Clean up
        throw std::runtime_error("Failed to open temporary file for writing");
    }

    // Write to temporary file
    for (const auto& [node, links] : graph) {
        fprintf(file, ">%lu LN:i:100", node);  // Dummy length
        for (const auto& [linked_node, orientation] : links) {
            fprintf(file, " L:%c:%lu:%c", orientation ? '+' : '-', linked_node, orientation ? '+' : '-');
        }
        fprintf(file, "\n");
        fprintf(file, "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"); // Dummy sequence
    }

    // Close the file descriptor to ensure all data is written
    fclose(file);

    std::cerr << "Temporary file created at: " << template_name << std::endl;
    std::cout << "Absolute path: " << std::filesystem::absolute(template_name) << std::endl;

    return std::string(template_name);
}


// BFS-based verification of connected components
std::vector<std::vector<uint64_t>> computeConnectedComponentsBFS(
    const std::vector<std::pair<uint64_t, std::vector<std::pair<uint64_t, bool>>>>& graph) {

    // First build adjacency list
    std::vector<std::vector<uint64_t>> adj(graph.size());
    for (uint64_t i = 0; i < graph.size(); ++i) {
        for (const auto& [j, _] : graph[i].second) {
            adj[i].push_back(j);
            // Since it's an undirected graph, add the reverse edge if not already present
            if (std::find(adj[j].begin(), adj[j].end(), i) == adj[j].end()) {
                adj[j].push_back(i);
            }
        }
    }

    std::vector<std::vector<uint64_t>> components;
    std::vector<bool> visited(graph.size(), false);

    for (uint64_t i = 0; i < graph.size(); ++i) {
        if (!visited[i]) {
            std::vector<uint64_t> component;
            std::queue<uint64_t> q;
            q.push(i);
            visited[i] = true;

            while (!q.empty()) {
                uint64_t current = q.front();
                q.pop();
                component.push_back(current);

                for (uint64_t neighbor : adj[current]) {
                    if (!visited[neighbor]) {
                        visited[neighbor] = true;
                        q.push(neighbor);
                    }
                }
            }

            std::sort(component.begin(), component.end());
            components.push_back(component);
        }
    }

    return components;
}

// Helper function to generate random graphs
std::vector<std::pair<uint64_t, std::vector<std::pair<uint64_t, bool>>>> generateRandomGraph(
    size_t num_nodes,
    double edge_probability,
    uint64_t seed = 42,
    bool directed = false) {

    std::mt19937 gen(seed);
    std::uniform_real_distribution<> dis(0.0, 1.0);
    std::uniform_int_distribution<int> orient_dis(0, 1);  // Changed from bool to int

    std::vector<std::pair<uint64_t, std::vector<std::pair<uint64_t, bool>>>> graph(num_nodes);

    for (uint64_t i = 0; i < num_nodes; ++i) {
        graph[i].first = i;
        for (uint64_t j = i + 1; j < num_nodes; ++j) {
            if (dis(gen) < edge_probability) {
                // Add edge in both directions with random orientation
                graph[i].second.emplace_back(j, static_cast<bool>(orient_dis(gen)));
                if (!directed) {
                    graph[j].second.emplace_back(i, static_cast<bool>(orient_dis(gen)));
                }
            }
        }
    }

    return graph;
}


class GraphHandlerTest : public ::testing::Test {
protected:
    std::vector<std::string> temp_files;

    void SetUp() override {
        spdlog::set_level(spdlog::level::err);
        spdlog::set_pattern("[%^%l%$] %v");
    }

    // void TearDown() override {
    //     for (const auto& file : temp_files) {
    //         if (fs::exists(file)) {
    //             fs::remove(file);
    //         }
    //     }
    //     temp_files.clear();
    // }

    std::string registerTempFile(const std::string& filename) {
        temp_files.push_back(filename);
        return filename;
    }

    // Helper to compare components (order-independent)
    void assertComponentsEqual(
        const std::vector<std::vector<uint64_t>>& actual,
        const std::vector<std::vector<uint64_t>>& expected) {

        ASSERT_EQ(actual.size(), expected.size());

        // Make copies to sort
        auto sorted_actual = actual;
        auto sorted_expected = expected;

        for (auto& comp : sorted_actual) {
            std::sort(comp.begin(), comp.end());
        }
        for (auto& comp : sorted_expected) {
            std::sort(comp.begin(), comp.end());
        }

        std::sort(sorted_actual.begin(), sorted_actual.end());
        std::sort(sorted_expected.begin(), sorted_expected.end());

        for (size_t i = 0; i < sorted_actual.size(); ++i) {
            EXPECT_EQ(sorted_actual[i], sorted_expected[i])
                << "Component mismatch at index " << i
                << "\nActual:   " << vectorToString(sorted_actual[i])
                << "\nExpected: " << vectorToString(sorted_expected[i]);
        }
    }

    // Helper to convert vector to string for error messages
    std::string vectorToString(const std::vector<uint64_t>& vec) {
        std::ostringstream oss;
        oss << "[";
        for (size_t i = 0; i < vec.size(); ++i) {
            if (i > 0) oss << ", ";
            oss << vec[i];
        }
        oss << "]";
        return oss.str();
    }

    // Helper to verify all nodes are present exactly once
    void verifyAllNodesPresent(
        const std::vector<std::vector<uint64_t>>& components,
        size_t expected_node_count) {

        std::set<uint64_t> nodes;
        for (const auto& comp : components) {
            for (uint64_t node : comp) {
                EXPECT_LT(node, expected_node_count) << "Node ID out of range: " << node;
                nodes.insert(node);
            }
        }

        EXPECT_EQ(nodes.size(), expected_node_count)
            << "Not all nodes are present in components. Expected "
            << expected_node_count << " nodes, found " << nodes.size();

        for (uint64_t i = 0; i < expected_node_count; ++i) {
            EXPECT_TRUE(nodes.count(i))
                << "Node " << i << " is missing from components";
        }
    }

    // Helper to verify component sizes match expected distribution
    void verifyComponentSizes(
        const std::vector<std::vector<uint64_t>>& components,
        const std::vector<size_t>& expected_sizes) {

        std::vector<size_t> actual_sizes;
        for (const auto& comp : components) {
            actual_sizes.push_back(comp.size());
        }

        std::sort(actual_sizes.begin(), actual_sizes.end());
        std::vector<size_t> sorted_expected = expected_sizes;
        std::sort(sorted_expected.begin(), sorted_expected.end());

        EXPECT_EQ(actual_sizes, sorted_expected)
            << "Component size distribution mismatch.\n"
            << "Actual sizes:   " << vectorToString(actual_sizes) << "\n"
            << "Expected sizes: " << vectorToString(sorted_expected);
    }
};

TEST_F(GraphHandlerTest, EmptyGraph) {
    auto graph = std::vector<std::pair<uint64_t, std::vector<std::pair<uint64_t, bool>>>>();
    std::string filename = registerTempFile(createTempFastaFile(graph));

    GraphHandler handler(filename);
    handler.read_graph_into_connected_components();

    auto components = handler.get_components();
    EXPECT_TRUE(components.empty());
}

TEST_F(GraphHandlerTest, SingleNode) {
    auto graph = std::vector<std::pair<uint64_t, std::vector<std::pair<uint64_t, bool>>>>{
        {0, {}}
    };
    std::string filename = registerTempFile(createTempFastaFile(graph));


    GraphHandler handler(filename);
    handler.read_graph_into_connected_components();

    auto components = handler.get_components();
    auto expected = computeConnectedComponentsBFS(graph);

    assertComponentsEqual(components, expected);
    verifyAllNodesPresent(components, 1);
}

TEST_F(GraphHandlerTest, TwoConnectedNodes) {
    auto graph = std::vector<std::pair<uint64_t, std::vector<std::pair<uint64_t, bool>>>>{
        {0, {{1, true}}},
        {1, {{0, false}}}
    };
    std::string filename = registerTempFile(createTempFastaFile(graph));

    GraphHandler handler(filename);
    handler.read_graph_into_connected_components();

    auto components = handler.get_components();
    auto expected = computeConnectedComponentsBFS(graph);

    assertComponentsEqual(components, expected);
    verifyAllNodesPresent(components, 2);
}

TEST_F(GraphHandlerTest, TwoDisconnectedNodes) {
    auto graph = std::vector<std::pair<uint64_t, std::vector<std::pair<uint64_t, bool>>>>{
        {0, {}},
        {1, {}}
    };
    std::string filename = registerTempFile(createTempFastaFile(graph));

    GraphHandler handler(filename);
    handler.read_graph_into_connected_components();

    auto components = handler.get_components();
    auto expected = computeConnectedComponentsBFS(graph);

    assertComponentsEqual(components, expected);
    verifyAllNodesPresent(components, 2);
    verifyComponentSizes(components, {1, 1});
}

TEST_F(GraphHandlerTest, ThreeNodesLine) {
    auto graph = std::vector<std::pair<uint64_t, std::vector<std::pair<uint64_t, bool>>>>{
        {0, {{1, true}}},
        {1, {{0, false}, {2, true}}},
        {2, {{1, false}}}
    };
    std::string filename = registerTempFile(createTempFastaFile(graph));

    GraphHandler handler(filename);
    handler.read_graph_into_connected_components();

    auto components = handler.get_components();
    auto expected = computeConnectedComponentsBFS(graph);

    assertComponentsEqual(components, expected);
    verifyAllNodesPresent(components, 3);
    verifyComponentSizes(components, {3});
}

TEST_F(GraphHandlerTest, ThreeNodesStar) {
    auto graph = std::vector<std::pair<uint64_t, std::vector<std::pair<uint64_t, bool>>>>{
        {0, {{1, true}, {2, true}}},
        {1, {{0, false}}},
        {2, {{0, false}}}
    };
    std::string filename = registerTempFile(createTempFastaFile(graph));

    GraphHandler handler(filename);
    handler.read_graph_into_connected_components();

    auto components = handler.get_components();
    auto expected = computeConnectedComponentsBFS(graph);

    assertComponentsEqual(components, expected);
    verifyAllNodesPresent(components, 3);
    verifyComponentSizes(components, {3});
}

TEST_F(GraphHandlerTest, ComplexGraphWithMultipleComponents) {
    // Graph with three components:
    // Component 1: 0-1-2 (line)
    // Component 2: 3-4 (pair)
    // Component 3: 5 (isolated)
    auto graph = std::vector<std::pair<uint64_t, std::vector<std::pair<uint64_t, bool>>>>{
        {0, {{1, true}}},
        {1, {{0, false}, {2, true}}},
        {2, {{1, false}}},
        {3, {{4, true}}},
        {4, {{3, false}}},
        {5, {}}
    };
    std::string filename = registerTempFile(createTempFastaFile(graph));

    GraphHandler handler(filename);
    handler.read_graph_into_connected_components();

    auto components = handler.get_components();
    auto expected = computeConnectedComponentsBFS(graph);

    assertComponentsEqual(components, expected);
    verifyAllNodesPresent(components, 6);
    verifyComponentSizes(components, {1, 2, 3});
}

TEST_F(GraphHandlerTest, RandomSmallGraph) {
    for (int seed = 0; seed < 10; ++seed) {
        auto graph = generateRandomGraph(20, 0.2, seed);
        std::string filename = registerTempFile(createTempFastaFile(graph));

        GraphHandler handler(filename);
        handler.read_graph_into_connected_components();

        auto actual_components = handler.get_components();
        auto expected_components = computeConnectedComponentsBFS(graph);

        assertComponentsEqual(actual_components, expected_components);
        verifyAllNodesPresent(actual_components, graph.size());
    }
}

TEST_F(GraphHandlerTest, RandomMediumGraph) {
    const size_t num_nodes = 100;
    const double edge_prob = 0.05;
    const uint64_t seed_start = 12345;

    for (int seed {seed_start}; seed < seed_start + MAX_NUM_RANDOM_ITERATIONS; seed++){

        auto graph = generateRandomGraph(num_nodes, edge_prob, seed);
        std::string filename = registerTempFile(createTempFastaFile(graph));

        GraphHandler handler(filename);
        handler.read_graph_into_connected_components();

        auto actual_components = handler.get_components();
        auto expected_components = computeConnectedComponentsBFS(graph);

        assertComponentsEqual(actual_components, expected_components);
        verifyAllNodesPresent(actual_components, num_nodes);
    }
}

TEST_F(GraphHandlerTest, PerformanceLargeGraph) {
    const size_t num_nodes = 1000;
    const double edge_prob = 0.01;
    const uint64_t seed = 67890;

    auto graph = generateRandomGraph(num_nodes, edge_prob, seed);
    std::string filename = registerTempFile(createTempFastaFile(graph));

    GraphHandler handler(filename);

    // Time the connected components computation
    auto start = std::chrono::high_resolution_clock::now();
    handler.read_graph_into_connected_components();
    auto components = handler.get_components();
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    // Verify correctness with BFS
    auto expected_components = computeConnectedComponentsBFS(graph);
    assertComponentsEqual(components, expected_components);
    verifyAllNodesPresent(components, num_nodes);

    // Performance expectation - should complete in reasonable time
    EXPECT_LT(duration, 2000) << "Performance test took too long: " << duration << "ms";

    std::cout << "\nPerformance test completed in " << duration << "ms for "
              << num_nodes << " nodes with edge probability " << edge_prob << "\n";
}

TEST_F(GraphHandlerTest, RealExampleGraph) {
    std::string content = R"(
>0 LN:i:75 L:+:1:- L:+:5:-
TCATCGCCTCCGGGATTTGCGGATCTCGTCCACGACTTTGGACGTTTCGCTGCGGAGCATGGTCATATACGTGAT
>1 LN:i:83 L:+:0:- L:+:6:-
TTCGGGCGGCGGTAGATCGCCTTGATGCCCATCTTCTTCATCAGCGTGGCGACGTGTCGCCGCCCGGTCTCCAGACCTTCTCC
>2 LN:i:68 L:+:6:-
TTAAAGGTTTTTGATAATGGCTTTCGAAGTGCTGAAAAGTAATCCGTAATCAAAACCGTATACATGAG
>3 LN:i:72 L:+:8:-
GAGACCGGGCGGCTGCACGTCGCCACGCTGATGAAGAAGATGGGCATCGAGGCGATCTACCGTCGCCCGAAC
>4 LN:i:86
GTGTTCGGGCGGCGGTAGATCGCCTCGATGCCCATCTTCTTCATCAGCGTGGCGACGTGTCGCCGCCCGGTCTCCAGACCTTCTCC
>5 LN:i:69 L:+:0:- L:+:7:-
GCCGGTCTTATCACGTATATGACCATGCTCCGCAGCGAAACGTCCAAAGTCGTGGACGAGATCCGCAAA
>6 LN:i:62 L:+:2:- L:+:8:-
ATCTCATGTATACGGTTTTGATTACGGATTACTTTTCAGCACTTCGAAAGCCATTATCAAAA
>7 LN:i:64 L:+:7:- L:-:5:+
CTCGTCCACGACTTTGGACGTTTCGCTGCGGAGCATGGTCATATACGTGATAAGACCGGCTACA
>8 LN:i:62 L:+:3:- L:+:6:+
GTGTTCGGGCGACGGTAGATCGCCTCGATGCCCATCTTCTTCATCAGCGTGGCGACGTGCAG
>9 LN:i:70
AGCTTTTCTTCCTGTGCCTTGATTTTGCTGTCTAGTGTTTCAGTAACGCTCATTAATCCAATACCATCAT
>10 LN:i:85
TTGCTCCACTGTAATTCTGCAGACTACGGGCTTTATTCATCGCGCAGCTGAACGCATTCGCACGATCGCCTGCACTTCCGTTATG
)";

    std::string filename = "real_example.fa";
    std::ofstream out(filename);
    out << content;
    out.close();
    registerTempFile(filename);

    GraphHandler handler(filename);
    handler.read_graph_into_connected_components();

    auto components = handler.get_components();

    // Build graph structure for BFS verification
    std::vector<std::pair<uint64_t, std::vector<std::pair<uint64_t, bool>>>> graph(11);
    for (uint64_t i = 0; i < 11; ++i) {
        graph[i].first = i;
    }

    // Add edges based on the FASTA content
    graph[0].second = {{1, true}, {5, true}};
    graph[1].second = {{0, false}, {6, true}};
    graph[2].second = {{6, true}};
    graph[3].second = {{8, true}};
    // 4 has no links
    graph[5].second = {{0, false}, {7, true}};
    graph[6].second = {{2, false}, {8, true}};
    graph[7].second = {{5, false}};
    graph[8].second = {{3, false}, {6, false}};
    // 9 and 10 have no links

    auto expected_components = computeConnectedComponentsBFS(graph);
    assertComponentsEqual(components, expected_components);
    verifyAllNodesPresent(components, 11);
}

TEST_F(GraphHandlerTest, TestFileHandling) {
    // Test with a non-existent file
    EXPECT_NO_THROW({
        GraphHandler handler("nonexistent_file.fa");
        handler.read_graph_into_connected_components();
        auto components = handler.get_components();
        EXPECT_TRUE(components.empty());
    });

    // Test with an empty file
    std::string empty_filename = "empty_test.fa";
    std::ofstream empty_file(empty_filename);
    empty_file.close();
    registerTempFile(empty_filename);

    EXPECT_NO_THROW({
        GraphHandler handler(empty_filename);
        handler.read_graph_into_connected_components();
        auto components = handler.get_components();
        EXPECT_TRUE(components.empty());
    });
}

TEST_F(GraphHandlerTest, TestMalformedInput) {
    // Test with malformed node IDs
    std::string malformed_filename = "malformed_test.fa";
    std::ofstream malformed_file(malformed_filename);
    malformed_file << ">not_a_number LN:i:10\n";
    malformed_file << "ACGT\n";
    malformed_file.close();
    registerTempFile(malformed_filename);

    EXPECT_NO_THROW({
        GraphHandler handler(malformed_filename);
        handler.read_graph_into_connected_components();
    });

    // Test with malformed link
    std::string malformed_link_filename = "malformed_link_test.fa";
    std::ofstream malformed_link_file(malformed_link_filename);
    malformed_link_file << ">0 LN:i:10 L:+:not_a_number:1\n";
    malformed_link_file << "ACGT\n";
    malformed_link_file.close();
    registerTempFile(malformed_link_filename);

    EXPECT_NO_THROW({
        GraphHandler handler(malformed_link_filename);
        handler.read_graph_into_connected_components();
    });
}

TEST_F(GraphHandlerTest, StressTest) {
    const size_t num_nodes = 500000;
    const double edge_prob = 0.002; // Sparse graph to keep it manageable
    const uint64_t seed = 424242;

    auto graph = generateRandomGraph(num_nodes, edge_prob, seed);
    std::string filename = registerTempFile(createTempFastaFile(graph));

    GraphHandler handler(filename);

    // Time the connected components computation
    auto start = std::chrono::high_resolution_clock::now();
    handler.read_graph_into_connected_components();
    auto components = handler.get_components();
    auto end = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    // Verify correctness with BFS on a subset (for performance)
    if (num_nodes <= 1000) {
        auto expected_components = computeConnectedComponentsBFS(graph);
        assertComponentsEqual(components, expected_components);
    }

    verifyAllNodesPresent(components, num_nodes);

    std::cout << "\nStress test completed in " << duration << "ms for "
              << num_nodes << " nodes with edge probability " << edge_prob << "\n";
}