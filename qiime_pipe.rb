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

  opts.on("-s", "--file_sff <file>", String, "SFF file to process") do |o|
    options[:file_sff] = o
  end

  opts.on("-i", "--dir_illumina <dir>", String, "Illumina directory to process") do |o|
    options[:dir_illumina] = o
  end

  opts.on("-m", "--file_map <file>", String, "Mapping file to process") do |o|
    options[:file_map] = o
  end

  opts.on("-r", "--remote_map <key>", String, "Remote mapping file to process") do |o|
    options[:remote_map] = o
  end

  opts.on("-o", "--dir_out <dir>", String, "Output directory") do |o|
    options[:dir_out] = o
  end

  opts.on("-f", "--force", "Force overwrite log file and output directory") do |o|
    options[:force] = o
  end

  opts.on("-d", "--denoise", "Denoise data") do |o|
    options[:denoise] = o
  end

  opts.on("-c", "--chimera", "Chimera filter data") do |o|
    options[:chimera] = o
  end

  opts.on("-M", "--catagory <string>", String, "Mapping catagory (mapping file column label") do |o| 
    options[:catagory] = o
  end

  options[:chimera_db] = Qiime::DEFAULT_CHIMERA_DB
  opts.on("-D", "--chimera_db <file>", String, "Chimere database (#{Qiime::DEFAULT_CHIMERA_DB})") do |o|
    options[:chimera_db] = o || Qiime::DEFAULT_CHIMERA_DB
  end

  options[:barcode_size] = Qiime::DEFAULT_BARCODE_SIZE
  opts.on("-b", "--barcode_size <int>", Integer, "Size of barcodes used (#{Qiime::DEFAULT_BARCODE_SIZE})") do |o|
    options[:barcode_size] = o
  end

  options[:cpus] = Qiime::DEFAULT_CPUS
  opts.on("-C", "--cpus <int>", Integer, "Number of CPUs to use (#{Qiime::DEFAULT_CPUS})") do |o|
    options[:cpus] = o
  end

  opts.on("-e", "--email <string>", String, "Send email alert") do |o| 
    options[:email] = o
  end
end.parse!

unless options[:file_map] or options[:remote_map]
  raise OptionParser::MissingArgument, "--file_map or --remote_map"
end

if options[:file_map] and not File.file? options[:file_map]
  raise OptionParser::InvalidOption, "no such file: #{options[:file_map]}"
end

unless options[:file_sff] or options[:dir_illumina]
  raise OptionParser::MissingArgument, "--file_sff or --dir_illumina must be specified"
end

if options[:sff_file] and not File.file?(options[:file_sff])
  raise OptionParser::InvalidOption, "no such file: #{options[:file_sff]}"
end

if options[:dir_illumina] and not File.directory?(options[:dir_illumina])
  raise OptionParser::InvalidOption, "no such directory: #{options[:dir_illumina]}"
end

raise OptionParser::MissingArgument, "--dir_out"  if options[:dir_out].nil?

q = Qiime::Pipeline.new(options)
q.log_delete                   if options[:force]
q.log_init(cmd_init)
q.dir_delete                   if options[:force]
q.dir_create
q.print_qiime_config
q.load_remote_mapping_file     if options[:remote_map]
q.check_id_map
q.process_sff                  if options[:file_sff]
q.process_illumina             if options[:dir_illumina]
q.split_libraries              if options[:file_sff]
q.denoiser_preprocessor        if options[:denoise]
q.denoiser                     if options[:denoise]
q.inflate_denoiser_output      if options[:denoise]
q.identify_chimeric_seq        if options[:chimera]
q.filter_fasta                 if options[:chimera]
#q.chimera_check                if options[:chimera]
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

if options[:email]
  if options[:file_sff]
    project = File.basename(options[:file_sff])
  else
    project = options[:dir_illumina]
  end

  q.send_mail("Finished: #{project}")
end

puts "All done."

END { q.send_mail("Interrupted") if options[:email] }
