#!/usr/bin/env ruby

require 'pp'
require 'optparse'
require_relative 'lib/qiime'

ARGV << "-h" if ARGV.empty?

cmd_init = File.basename($0) + " " + ARGV.join(" ")

def option_parser(args)
  options = {}

  OptionParser.new do |opts|
    opts.banner = <<USAGE
Description: #{File.basename(__FILE__)} is used to process a 16S RNA amplicon dataset using QIIME.
  the dataset can either be a 454 dataset where the starting point is a SFF file or an Illumina
  dataset where the starting point is a directory with demultiplexed FASTQ files.

  Mapping files can be supplied either as text files or as remote mapping files using Google Docs:
  http://qiime.org/tutorials/remote_mapping_files.html

  Note that for processing Illumina data the SampleID column must match the FASTQ file prefix. Since
  underscores '_' are not allowed in mapping files you must use dots '.' instead, however the #{File.basename(__FILE__)}
  script will translate these when matching file prefixes.

  The progress of the analysis will be recording the the <--dir_out>.log file that will also contain
  the command and options used with #{File.basename(__FILE__)}.

Usage: #{File.basename(__FILE__)} [options] <--illumina_dirs|--file_sff> <--remote_map|--file_map> <--dir_out>

Examples:
  Processing a 454 dataset:    
    
    #{File.basename(__FILE__)} -s GXS0P3T01.sff -m map_file.txt -o Result
  
  Processing an Illumina dataset:

    #{File.basename(__FILE__)} -i Fastq/ -r 0AnzomiBiZW0ddDVrdENlNG5lTWpBTm5kNjRGbjVpQmc -o Result
    
  Restarting a stopped job:

    #{File.basename(__FILE__)} --restart Result.log

USAGE

    opts.on("-h", "--help", "Display this screen" ) do
      $stderr.puts opts
      exit
    end

    opts.on("--nofigures", "Skip figure generation") do |o|
      options[:nofigures] = o
    end

    opts.on("--notree", "Skip alignment and tree building") do |o|
      options[:notree] = o
    end

    opts.on("--restart <log file>", "Restart job") do |o|
      options[:restart] = o
    end

    opts.on("-s", "--file_sff <file>", String, "SFF file to process") do |o|
      options[:file_sff] = o
    end

    opts.on("-i", "--illumina_dirs <dir>[,<dir>[,<dir]] ...", Array, "Illumina directories to process") do |o|
      options[:illumina_dirs] = o
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

    opts.on("--parameter_file <file>", String, "QIIME parameters file") do |o|
      options[:parameter_file] = o
    end

    opts.on("--trim_primers", "Trim primers from reads prior to assembly") do |o|
      options[:trim_primers] = o
    end

    opts.on("-M", "--category <string>", String, "Mapping category (mapping file column label") do |o| 
      options[:catagory] = o
    end

    opts.on("-D", "--chimera_db <file>", String, "Chimera database (#{Qiime::DEFAULT_CHIMERA_DB})") do |o|
      options[:chimera_db] = o || Qiime::DEFAULT_CHIMERA_DB #TODO: removing this default doesnt work for some reason
    end

    opts.on("-b", "--barcode_size <int>", Integer, "Size of barcodes used (#{Qiime::DEFAULT_BARCODE_SIZE})") do |o|
      options[:barcode_size] = o
    end

    opts.on("-C", "--cpus <int>", Integer, "Number of CPUs to use (#{Qiime::DEFAULT_CPUS})") do |o|
      options[:cpus] = o
    end

    opts.on("-R", "--r_starter <git url>", String, "Git URL for cloning the R starter pack (#{Qiime::DEFAULT_R_STARTER})") do |o|
      options[:r_starter] = o || Qiime::DEFAULT_R_STARTER
    end

    opts.on("-e", "--email <string>", String, "Send email alert") do |o| 
      options[:email] = o
    end
  end.parse!(args)

  options
end

options = option_parser(ARGV)

if options[:restart]
  raise OptionParser::InvalidOption, "no such file: #{options[:restart]}" unless File.file? options[:restart]

  File.open(options[:restart]) do |ios|
    line = ios.gets

    if line[0] == '#'
      options = option_parser(line[1 .. -1].split " ")
    else
      raise "Failed to parse log file: #{options[:restart]}"
    end
  end
end

unless options[:file_map] or options[:remote_map]
  raise OptionParser::MissingArgument, "--file_map or --remote_map"
end

if options[:file_map] and not File.file? options[:file_map]
  raise OptionParser::InvalidOption, "no such file: #{options[:file_map]}"
end

unless options[:file_sff] or options[:illumina_dirs]
  raise OptionParser::MissingArgument, "--file_sff or --illumina_dirs must be specified"
end

if options[:sff_file] and not File.file?(options[:file_sff])
  raise OptionParser::InvalidOption, "no such file: #{options[:file_sff]}"
end

if options[:file_parameters] and not File.file?(options[:file_parameters])
  raise OptionParser::InvalidOption, "no such file: #{options[:file_parameters]}"
end

if options[:illumina_dirs]
  options[:illumina_dirs].each do |dir|
    raise OptionParser::InvalidOption, "no such directory: #{dir}" unless File.directory? dir
  end
end

#This option cannot be overwritten by parameter files, as in process_illumina.rb
options[:cpus] ||= Qiime::DEFAULT_CPUS 

raise OptionParser::MissingArgument, "--dir_out"  if options[:dir_out].nil?
if options[:parameter_file]
  raise OptionParser::InvalidArgument, "No no such parameter_file: #{options[:parameter_file]}" unless File.file? options[:parameter_file]

  File.open options[:parameter_file] do |ios|
    ios.each do |line|
      line.chomp!

      script, leftover = line.split ':'
      option, value    = leftover.split /\s+/

      if script == File.basename(__FILE__).sub(/.rb$/, '')
        options[option.to_sym] ||= value
      end
    end
  end
end
options[:chimera_db] ||= Qiime::DEFAULT_CHIMERA_DB
options[:barcode_size] ||= Qiime::DEFAULT_BARCODE_SIZE
options[:r_starter] ||= Qiime::DEFAULT_R_STARTER


q = Qiime::Pipeline.new(options)
q.log_delete                   if options[:force]
q.log_init(cmd_init)
q.dir_delete                   if options[:force]
q.dir_create
q.print_qiime_config
q.load_remote_mapping_file     if options[:remote_map]
q.check_id_map
q.process_sff                  if options[:file_sff]
q.process_illumina             if options[:illumina_dirs]
q.split_libraries              if options[:file_sff]
q.denoiser_preprocessor        if options[:denoise]
q.denoiser                     if options[:denoise]
q.inflate_denoiser_output      if options[:denoise]
q.identify_chimeric_seq        if options[:chimera]
q.filter_fasta                 if options[:chimera]
#q.chimera_check                if options[:chimera]
q.pick_de_novo_otus
if !options[:nofigures]
  q.print_biom_table_summary
  q.make_otu_heatmap_html
  q.make_otu_network
  q.wf_taxa_summary
  q.alpha_diversity
  q.beta_diversity_through_plots
  q.jackknifed_beta_diversity
  q.make_bootstrapped_tree
  q.make_3d_plots
end

q.initialize_R_starter
if options[:email]
  if options[:file_sff]
    project = File.basename(options[:file_sff])
  else
    project = options[:illumina_dirs].inspect
  end

  q.send_mail("Finished: #{project}")
end

puts "All done."

END { q.send_mail("Interrupted") if options[:email] }
