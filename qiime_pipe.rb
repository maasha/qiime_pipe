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

  opts.on("-m", "--file_map <file>", String, "Mapping file to process") do |o|
    options[:file_map] = o
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

raise OptionParser::MissingArgument, "--file_sff" if options[:file_sff].nil?
raise OptionParser::MissingArgument, "--file_map" if options[:file_map].nil?
raise OptionParser::MissingArgument, "--dir_out"  if options[:dir_out].nil?
raise OptionParser::InvalidOption,   "no such file: #{options[:file_sff]}" unless File.file?(options[:file_sff])
raise OptionParser::InvalidOption,   "no such file: #{options[:file_map]}" unless File.file?(options[:file_map])

q = Qiime::Pipeline.new(options)
q.log_delete                   if options[:force]
q.log_init(cmd_init)
q.dir_delete                   if options[:force]
q.dir_create
q.print_qiime_config
q.check_id_map
q.process_sff
q.split_libraries
q.denoise_wrapper              if options[:denoise]
q.inflate_denoiser_output      if options[:denoise]
q.chimera_check                if options[:chimera]
q.pick_otus_through_otu_table
q.per_library_stats
q.make_otu_heatmap_html
q.make_otu_network
q.wf_taxa_summary
q.alpha_diversity
q.beta_diversity_through_plots
q.jackknifed_beta_diversity
q.make_bootstrapped_tree
q.make_3d_plots
q.send_mail("Finished: " + File.basename(options[:file_sff])) if options[:email]

puts "All done."

END { q.send_mail("Interrupted") if options[:email] }
