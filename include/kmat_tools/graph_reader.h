#ifndef GRAPH_READER_H
#define GRAPH_READER_H

#include <vector>
#include <stdint.h>
#include <sstream>
#include <string>
#include <iostream>

#include <fmt/format.h>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <kseq++/seqio.hpp>

std::vector<std::string> split(const std::string &s, char delimiter)
{
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);

  while (std::getline(tokenStream, token, delimiter))
  {
    tokens.push_back(token);
  }

  return tokens;
}

class GraphHandler
{
private:
  uint64_t m_unitigs_count;
  std::vector<uint64_t> m_group;
  std::vector<uint32_t> m_size;
  std::string m_filename;

  void union_set(uint64_t node_id1, uint64_t node_id2)
  {
    bool node_id1_in_group {is_in_group(node_id1)};
    bool node_id2_in_group {is_in_group(node_id2)};
    // std::cerr << fmt::format(" {}:{} ; {}:{}", node_id1, node_id1_in_group, node_id2, node_id2_in_group) << std::endl;
    // if both are in a group
    if (node_id1_in_group && node_id2_in_group){
      // find the set leader of both
      uint64_t set_id1 = find_set(node_id1);
      uint64_t set_id2 = find_set(node_id2);
      // std::cerr << fmt::format("[2] n:{} -> s:{}; n:{} -> s:{}", node_id1, set_id1, node_id2, set_id2) << std::endl;

      // if the same do nothing, if different add the smallest to the larger
      if (set_id1 != set_id2){
        // using ternary operator to assign smaller and larger
        uint64_t larger_set = (m_size[set_id1] >= m_size[set_id2]) ? set_id1 : set_id2 ;
        uint64_t smaller_set = (larger_set == set_id1) ? set_id2 : set_id1;

        // assigning smaller to larger
        m_group[smaller_set] = larger_set;
        m_size[larger_set] += m_size[smaller_set];
        // std::cerr << fmt::format("[2a] n:{} -> s:{}", smaller_set, larger_set) << std::endl;

      }
    }

    // else if both not in a group
    else if (!(node_id1_in_group || node_id2_in_group)) {
      // assigning node_id2 to the group of node_id1
      m_group[node_id2] = node_id1;
      m_size[node_id1] = 2;
      // std::cerr << fmt::format("[0] n:{} -> s:{}", node_id2, node_id1) << std::endl;
    }

    // else if just one in a group
    else {
      // understanding which one to assign and which one is in group (ternary operator)
      uint64_t to_add = (node_id1_in_group == true) ? node_id2 : node_id1;
      uint64_t in_group = (node_id1_in_group == true) ? node_id1 : node_id2;

      // assigning it to the one already in group.
      uint64_t set_id = find_set(in_group);
      m_group[to_add] = set_id;
      m_size[set_id] += 1;
      // std::cerr << fmt::format("[1] n:{} -> s:{}", to_add, set_id) << std::endl;
    }
  }

  uint64_t find_set(uint64_t node_id)
  {
    // find the root
    uint64_t current_node = node_id;
    while(m_group[current_node] != current_node){
      current_node = m_group[current_node];
    }
    uint64_t root = current_node;

    // reassign as group 'leader' the root to make the tree flat and not deep
    current_node = node_id;
    while(current_node != root){
      current_node = m_group[current_node];
      m_group[current_node] = root;
    }

    // return the root as the set leader
    return root;
  }

  bool is_in_group(uint64_t node_id)
  {
    if (m_group[node_id] == node_id && m_size[node_id] == 0) {return false;}
    else {return true;}
  }

public:
  GraphHandler(std::string filename)
  {
    // Opening file
    spdlog::info(fmt::format("[Graph-Handler] Opening {}.", filename.c_str()));
    klibpp::SeqStreamIn seq_stream(filename.c_str());
    if (!seq_stream)
    {
      spdlog::error(fmt::format("Error. Cannot open {} to read unitig link information.", filename));
    }
    m_filename = filename;
    // Count number of records
    klibpp::KSeq record;
    m_unitigs_count = 0;
    while (seq_stream >> record)
    {
      m_unitigs_count++;
    }
    spdlog::info(fmt::format("[Graph-Handler] Counted {} unitigs.", m_unitigs_count));
    // Resize vectors
    spdlog::info("[Graph-Handler] Resizing vectors.");
    m_group.resize(m_unitigs_count);
    m_size.resize(m_unitigs_count, 0);

    // Fill group vector
    spdlog::info("[Graph-Handler] Filling group vector.");
    for (uint64_t utg{0}; utg < m_unitigs_count; utg++)
    {
      m_group[utg] = utg;
    }
    spdlog::info("[Graph-Handler] DONE.");
  }

  ~GraphHandler() = default;

  void read_graph_into_connected_components()
  {
    // Open file
    spdlog::info(fmt::format("[Graph-Handler] Re-Opening {}.", m_filename.c_str()));
    klibpp::SeqStreamIn seq_stream(m_filename.c_str());
    if (!seq_stream)
    {
      spdlog::error(fmt::format("Error. Cannot open {} to read untig link information.", m_filename));
    }

    // Read sequences and computes connected components using union find

    std::string utg_id;
    klibpp::KSeq record;
    char fasta_comment_delimiter{' '};
    char link_delimiter{':'};
    spdlog::info(fmt::format("[Graph-Handler] Opening {} for unitig link extraction.", m_filename.c_str()));
    while (seq_stream >> record)
    {
      uint64_t utg_id, linked_node;
      bool is_used {false};
      spdlog::debug(fmt::format("Reading utg {} with metadata {} and sequence {}", record.name, record.comment, record.seq));
      try{
      utg_id = std::stoull(record.name); // record.seq  record.seq.size()
      }
      catch (const std::invalid_argument& ia){
        spdlog::error(fmt::format("Error. Cannot convert {} to integer in unitig name parsing.", record.name));
        continue;
      }
      std::vector<std::string> unitig_records = split(record.comment, fasta_comment_delimiter);
      for (const auto &info_field : unitig_records)
      {
        if (info_field.compare(0, 2, "L:") == 0)
        { // IT IS A LINK RECORD
          std::vector<std::string> link_records = split(info_field, link_delimiter);
          try{
            linked_node = std::stoull(link_records[2]);
          }
          catch (const std::invalid_argument& ia){
            spdlog::error(fmt::format("Error. Cannot convert {} to integer in link parsing.", link_records[1]));
            continue;
          }
          if(utg_id != linked_node){
            union_set(utg_id, linked_node);
            // spdlog::error(fmt::format("UNION OPERATION BETWEEN {} and {}", utg_id, linked_node));
            is_used = true;
          }
        }
      }
      if (!is_used){m_size[utg_id] = 1; spdlog::debug(fmt::format("NODE {} had no links. It is still its own component.", utg_id));}
    }
  spdlog::info("[Graph-Handler] Done. Now getting components.");
  }

  std::vector<std::vector<uint64_t>> get_components(){
    // In order to avoid many reallocations, I will do this in 2 passes.
    // Pass 1 is to count the # of connected components and their size;
    // Then I allocate the vector of vectors accordingly to then avoid any reallocations
    // Finally I scan again and put the nodes into their cc
    spdlog::info("[Graph-Handler] Dumping components.");
    std::unordered_map<uint64_t, uint64_t> rootCount;
    for(uint64_t utg {0}; utg < m_unitigs_count; utg++){
      uint64_t set_repr {find_set(utg)};
      // std::cerr << fmt::format("UTG {} -> {}", utg, set_repr) << std::endl;
      rootCount[set_repr]++;
    }

    std::vector<std::vector<uint64_t>> components;
    components.reserve(rootCount.size());

    std::unordered_map<uint64_t, std::vector<uint64_t>> rootToNodes;
    for (uint64_t utg {0}; utg < m_unitigs_count; utg++) {
        uint64_t root = find_set(utg);
        if (rootToNodes[root].empty()) {
            rootToNodes[root].reserve(rootCount[root]); // Pre-allocating
        }
        rootToNodes[root].push_back(utg);
    }

    for (auto& [root, nodes] : rootToNodes) {
        components.push_back(move(nodes));
    }

    return components;
  }

  std::vector<uint64_t> get_component_of(uint64_t node);
};

#endif