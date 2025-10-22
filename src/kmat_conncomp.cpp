#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>

#include <kseq++/seqio.hpp>
#include "../external/sshash/dictionary.hpp"

#include <kmat_tools/cmd/conncomp.h>
#include <kmat_tools/matrix.h>
#include <kmat_tools/task.h>
#include <kmat_tools/utils.h>

#include <kmat_tools/aggregator.h>
#include <kmat_tools/matrix_writer.h>
#include <kmat_tools/graph_reader.h>

namespace fs = std::filesystem;


namespace kmat {

int main_connected_component(conncomp_opt_t opt)
{
    // input validation
    fs::path unitig_path = opt->inputs[0];
    if(!fs::is_regular_file(unitig_path)) {
        throw std::runtime_error(fmt::format("unitig file \"{}\" does not exist", unitig_path.c_str()));
    }

    fs::path matrix_path = opt->inputs[1];
    if(!fs::is_regular_file(matrix_path)) {
        throw std::runtime_error(fmt::format("k-mer matrix file \"{}\" does not exist", matrix_path.c_str()));
    }

    if(opt->mini_size >= opt->kmer_size) {
        throw std::runtime_error("minimizer size must be smaller than k-mer size");
    }

    spdlog::info("building k-mer dictionary.");
    spdlog::info(fmt::format("Storing info for {}.sshash.log", opt->prefix));


    sshash::dictionary kmer_dict;
    {
        std::string sshash_logfile = fmt::format("{}.sshash.log", opt->prefix);
        std::ofstream ofs(sshash_logfile, std::ios::out);
        std::streambuf *coutbuf = std::cout.rdbuf();
        if (ofs.good()) {
          std::cout.rdbuf(ofs.rdbuf());
        }

        sshash::build_configuration build_config;
        build_config.k = opt->kmer_size;
        build_config.m = opt->mini_size;
        build_config.c = 5.0;
        build_config.pthash_threads = opt->nb_threads;
        build_config.canonical_parsing = true;
        build_config.verbose = false;
        kmer_dict.build(unitig_path, build_config);

        std::cout.rdbuf(coutbuf);
    }
    spdlog::info("Finished building k-mer dictionary");
    size_t max_utgs {kmer_dict.num_contigs()};
    spdlog::debug(fmt::format("k-mer processed: {}", kmer_dict.size()));
    spdlog::debug(fmt::format("unitigs processed: {}", max_utgs));

    ////////////////////////////////////////////////

    // COMPUTING CONNECTED COMPONENTS
    spdlog::info(fmt::format("Handling {} graph.", unitig_path.c_str()));
    GraphHandler handler(unitig_path);
    spdlog::info("Computing connected components.");
    handler.read_graph_into_connected_components();
    std::vector<std::vector<uint64_t>> connected_components = handler.get_components();
    spdlog::info(fmt::format("Computed {} connected components.", connected_components.size()));
    // REORGANIZING UNITIGS IDS FOR CONNECTED COMPONENTS METRICS
    // TODO: I WILL DO IT WHEN ALL THE OTHER LOGICS ARE IMPLEMENTED

    // COMPUTING STATS

    spdlog::info("aggregating k-mer counts");

    TextMatrixReader mat(matrix_path);

    std::string kmer;
    std::vector<uint32_t> kmer_counts;
    bool has_kmer = mat.read_kmer_counts(kmer,kmer_counts);

    std::size_t nb_samples {has_kmer ? kmer_counts.size() : 0};
    spdlog::debug(fmt::format("samples: {}", nb_samples));

    size_t number_unitigs {kmer_dict.num_contigs()};

    std::unique_ptr<Aggregator> aggregator;
    if (opt->abundance_metric == "mean") {
        spdlog::info(fmt::format("Computing mean ({}) for components.", opt->abundance_metric));
        aggregator = std::make_unique<MeanAggregator>(nb_samples, number_unitigs, opt->min_frac );
    }
    else { // median
        spdlog::info(fmt::format("Computing median ({}) for components.", opt->abundance_metric));
        aggregator = std::make_unique<MedianAggregator>(nb_samples, number_unitigs, opt->min_frac );
    }

    while(has_kmer) {
        auto res = kmer_dict.lookup_advanced(kmer.c_str());
        if (res.kmer_id == sshash::constants::invalid_uint64) {
            has_kmer = mat.read_kmer_counts(kmer,kmer_counts);
            continue;
        }

        std::size_t utg_id = res.contig_id;
        aggregator->process_kmer(utg_id, kmer_counts);

        has_kmer = mat.read_kmer_counts(kmer,kmer_counts);
    }

    ////////////////////////////////////////////////

    // WRITING UNITIG MATRIX TO DISK
    spdlog::info("writing unitig matrix");
    std:string utg_prefix = opt->prefix + "unitigs";
    std::unique_ptr<MatrixWriter> writer;
    if (opt->output_format == "txt") {
        // NORMAL TEXT FILE
        writer = std::make_unique<TextMatrixWriter>(utg_prefix, opt->write_frac_matrix );
    } else if (opt->output_format == "tsv") {
        // TSV COMPRESSED FOR DOWNSTREAM IN MEMORY
        writer = std::make_unique<CompressedTSVMatrixWriter>(utg_prefix, opt->write_frac_matrix, nb_samples);
    } else {
        // IF NOT RECOGNIZED DEFAULT TO TXT FOR BACKWARD COMPATIBILITY AND AVOID DISRUPTION
        spdlog::info(fmt::format("OUTPUT FORMAT {} NOT RECOGNIZED. DEFAULT TO TXT.", opt->output_format));
        writer = std::make_unique<TextMatrixWriter>(utg_prefix, opt->write_frac_matrix );
    }

    std::vector<double> utg_abundances(nb_samples);
    std::vector<double> utg_fractions(nb_samples);
    Aggregator::AbundanceFractions result;
    for(uint64_t utg_id=0; utg_id < max_utgs; utg_id++) {
        // FILL SAMPLES VECTOR WITH ABUNDANCE FRACTION FOR THE UTG
        result = aggregator->get_all_abundance_fractions(utg_id);
        utg_abundances = std::move(result.abundances);
        if (opt->write_frac_matrix) {
            utg_fractions = std::move(result.fractions);
        }
        // DUMP IT TO DISK
        std::string utg_identifier {">" + std::to_string(utg_id)};
        writer->write_row(utg_identifier, utg_abundances, utg_fractions);
        // REPEAT
    }

    // WRITING CONNECTED COMPONENTS MATRIX'
    spdlog::info("writing connected_components matrix");
    std::string cc_prefix = opt->prefix + "connected_components";
    if (opt->output_format == "txt") {
        // NORMAL TEXT FILE
        writer = std::make_unique<TextMatrixWriter>(cc_prefix, opt->write_frac_matrix );
    } else if (opt->output_format == "tsv") {
        // TSV COMPRESSED FOR DOWNSTREAM IN MEMORY
        writer = std::make_unique<CompressedTSVMatrixWriter>(cc_prefix, opt->write_frac_matrix, nb_samples);
    } else {
        // IF NOT RECOGNIZED DEFAULT TO TXT FOR BACKWARD COMPATIBILITY AND AVOID DISRUPTION
        spdlog::info(fmt::format("OUTPUT FORMAT {} NOT RECOGNIZED. DEFAULT TO TXT.", opt->output_format));
        writer = std::make_unique<TextMatrixWriter>(cc_prefix, opt->write_frac_matrix );
    }


    for(uint64_t cc_id=0; cc_id < connected_components.size(); cc_id++) {
        // FILL SAMPLES VECTOR WITH ABUNDANCE FRACTION FOR THE UTG
        result = aggregator->get_connected_component_abundance_fraction(connected_components[cc_id]);
        // DUMP IT TO DISK
        std::string utg_identifier {">CC_" + std::to_string(cc_id)};
        writer->write_row(utg_identifier, utg_abundances, utg_fractions);
        // REPEAT
    }

    // WRITING THE CONNECTED COMPONENT LIST IN GZIPPED TSV
    spdlog::info("writing connected_components list in tsv");

    std::string cc_nodes_prefix = opt->prefix + "connected_components_to_nodes.tsv.gz";
    CompressedTSVComponentWriter cc_list_writer(cc_nodes_prefix);
    for(uint64_t cc_id=0; cc_id < connected_components.size(); cc_id++) {
        std::string comp_id = fmt::format("CC_{}",cc_id);
        cc_list_writer.write_row(comp_id, connected_components[cc_id]);
    }

    return 0;
}


kmat_opt_t conncomp_cli(std::shared_ptr<bc::Parser<1>> cli, conncomp_opt_t opt)
{
    bc::cmd_t conncomp = cli->add_command("conncomp", "Create a connected component matrix, a unitig matrix and connected component list.");

    conncomp->add_group("main options", "");

    conncomp->add_param("-k/--kmer-size", fmt::format("k-mer size [8,{}].", KL[MUSET_KMER_N-1]-1))
        ->meta("INT")
        ->def("31")
        ->checker(bc::check::f::range(8, KL[MUSET_KMER_N-1]-1))
        ->setter(opt->kmer_size);

    conncomp->add_param("-p/--prefix", "output files prefix.")
        ->meta("FILE")
        ->def("out")
        ->setter(opt->prefix);

    conncomp->add_param("-f/--min-frac", "set unitig average abundance to 0 if its k-mer fraction is below this threshold [0,1].")
        ->meta("FLOAT")
        ->def("0.0")
        ->checker(bc::check::f::range(0.0, 1.0))
        ->setter(opt->min_frac);

    conncomp->add_param("--out-frac", "output an additional matrix containing k-mer fractions.")
        ->as_flag()
        ->setter(opt->write_frac_matrix);

    conncomp->add_group("other options", "");

    conncomp->add_param("-m/--minimizer-size", "minimizer size")
        ->meta("INT")
        ->def("15")
        ->setter(opt->mini_size);

    conncomp->add_param("-t/--threads", "number of threads.")
        ->meta("INT")
        ->def("4")
        ->checker(bc::check::is_number)
        ->setter(opt->nb_threads);

    conncomp->add_param("--abundance-metric", "metric to use for abundance: mean or median.")
    ->meta("STRING")
    ->def("mean") // Default to mean for backward compatibility
    ->checker(bc::check::f::in("mean|median"))
    ->setter(opt->abundance_metric);

    conncomp->add_param("--output-format", "Output format can be either 'txt' or 'tsv' (tsv is gzip compressed).")
    ->meta("STRING")
    ->def("txt") // Default to txt for backward compatibility
    ->checker(bc::check::f::in("txt|tsv"))
    ->setter(opt->output_format);

    conncomp->add_param("-h/--help", "show this message and exit.")
         ->as_flag()
         ->action(bc::Action::ShowHelp);

    conncomp->add_param("-v/--version", "show version and exit.")
         ->as_flag()
         ->action(bc::Action::ShowVersion);

    conncomp->set_positionals(2, "<unitigs.fasta> <kmer_matrix>", "a unitig fasta file and a text-based k-mer matrix");

    return opt;
}

};