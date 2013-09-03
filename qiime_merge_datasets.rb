#!/usr/bin/env ruby

require 'pp'
require 'optparse'
require_relative 'lib/qiime'

ARGV << "-h" if ARGV.empty?

cmd_init = File.basename($0) + " " + ARGV.join(" ")

options = {}

OptionParser.new do |opts|
  opts.banner = "Usage: #{File.basename(__FILE__)} [options]"

  opts.on("-h", "--help", "Display this screen" ) do
    $stderr.puts opts
    exit
  end

  opts.on("-m", "--mapping_files <files>", Array, "Mapping file to process") do |o|
    options[:mapping_files] = o
  end

  opts.on("-f", "--fasta_files <files>", Array, "inflated.fna FASTA files") do |o|
    options[:fasta_files] = o
  end

  opts.on("-o", "--dir_out <dir>", String, "Output directory") do |o|
    options[:dir_out] = o
  end

  opts.on("--force", "Force overwrite log file and output directory") do |o|
    options[:force] = o
  end

  opts.on("-M", "--catagory <string>", String, "Mapping catagory (mapping file column label)") do |o| 
    options[:catagory] = o
  end

  options[:cpus] = Qiime::DEFAULT_CPUS
  opts.on("-c", "--cpus <int>", Integer, "Number of CPUs to use (#{Qiime::DEFAULT_CPUS})") do |o|
    options[:cpus] = o
  end

  opts.on("-e", "--email <string>", String, "Send email alert") do |o| 
    options[:email] = o
  end
end.parse!

raise OptionParser::MissingArgument, "--mapping_files" if options[:mapping_files].nil?
raise OptionParser::MissingArgument, "--fasta_files"   if options[:fasta_files].nil?
raise OptionParser::MissingArgument, "--dir_out"       if options[:dir_out].nil?
raise OptionParser::InvalidOption, "Only 1 mapping file specified" if options[:mapping_files].size == 1
raise OptionParser::InvalidOption, "Only 1 fasta   file specified" if options[:fasta_files].size   == 1
raise OptionParser::InvalidOption, "Number of mapping and fasta files don't match" if options[:mapping_files].size != options[:fasta_files].size

q = Qiime::Pipeline.new(options)
q.log_delete                   if options[:force]
q.log_init(cmd_init)
q.dir_delete                   if options[:force]
q.dir_create
q.print_qiime_config
q.merge_id_maps
#q.check_id_map
q.merge_fasta_files
q.pick_de_novo_otus
q.print_biom_table_summary
q.make_otu_heatmap_html
q.make_otu_network
q.wf_taxa_summary
q.alpha_diversity
q.beta_diversity_through_plots
q.jackknifed_beta_diversity
q.make_bootstrapped_tree
q.make_3d_plots
q.send_mail("Finished merging datasets") if options[:email]

puts "All done."

END { q.send_mail("Interrupted") if options[:email] }
