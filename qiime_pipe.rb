#!/usr/bin/env ruby

require 'pp'
require 'optparse'

ARGV << "-h" if ARGV.empty?

default_chimera_db  = "/home/maasha/Install/QIIME1.6/data/Gold/gold.fa"
default_barcode_size = 10
default_cpus         = 1

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

  options[:chimera_db] = default_chimera_db
  opts.on("-D", "--chimera_db <file>", String, "Chimere database (#{default_chimera_db})") do |o|
    options[:chimera_db] = o || default_chimera_db
  end

  options[:barcode_size] = default_barcode_size
  opts.on("-b", "--barcode_size <int>", Integer, "Size of barcodes used (#{default_barcode_size})") do |o|
    options[:barcode_size] = o
  end

  options[:cpus] = default_cpus
  opts.on("-C", "--cpus <int>", Integer, "Number of CPUs to use (#{default_cpus})") do |o|
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
# raise OptionParser::InvalidOption,   "directory exists: #{options[:dir_out]} - use --force to overwrite" if File.directory?(options[:dir_out]) and not options[:force]

class Qiime
  def initialize(options)
    @options     = options
    @file_log    = @options[:dir_out] + ".log"
    @min_samples = 0
  end

  def log_delete
    File.delete(@file_log) if File.file? @file_log
  end

  def dir_delete
    require 'fileutils'

    FileUtils.rm_rf @options[:dir_out] if File.directory? @options[:dir_out]
  end

  def dir_create
    run "mkdir #{@options[:dir_out]}"
  end

  def print_qiime_config
    log = "#{@options[:dir_out]}/print_qiime_config.log"
    run "print_qiime_config.py -t > #{log} 2>&1"
  end

  def check_id_map
    dir_out = "#{@options[:dir_out]}/mapping_output"
    run "check_id_map.py -m #{@options[:file_map]} -o #{dir_out} > /dev/null"
  end

  def process_sff
    base         = File.basename(@options[:file_sff], ".sff")
    file_sff_txt = "#{@options[:dir_out]}/#{base}.sff.txt"
    file_fasta   = "#{@options[:dir_out]}/#{base}.fna"
    file_qual    = "#{@options[:dir_out]}/#{base}.qual"

    if @options[:denoise]
      run "sffinfo #{@options[:file_sff]} > #{file_sff_txt}"
      run "sffinfo -s #{@options[:file_sff]} > #{file_fasta}"
      run "sffinfo -q #{@options[:file_sff]} > #{file_qual}"
    else
      run "process_sff.py -i #{@options[:file_sff]} -o #{@options[:dir_out]}"
    end
  end

  def split_libraries
    dir_out    = "#{@options[:dir_out]}/split_library_output"
    base       = File.basename(@options[:file_sff], ".sff")
    file_fasta = "#{@options[:dir_out]}/#{base}.fna"
    file_qual  = "#{@options[:dir_out]}/#{base}.qual"
    run "split_libraries.py -b #{@options[:barcode_size]} -m #{@options[:file_map]} -f #{file_fasta} -q #{file_qual} -o #{dir_out}"
  end

  def denoise_wrapper
    dir_out      = "#{@options[:dir_out]}/denoised"
    base         = File.basename(@options[:file_sff], ".sff")
    file_sff_txt = "#{@options[:dir_out]}/#{base}.sff.txt"
    file_fasta   = "#{@options[:dir_out]}/split_library_output/seqs.fna"
    run "denoise_wrapper.py --titanium -i #{file_sff_txt} -f #{file_fasta} -o #{dir_out} -m #{@options[:file_map]} -n #{@options[:cpus]} --force_overwrite"
  end

  def inflate_denoiser_output
    file_centroids  = "#{@options[:dir_out]}/denoised/centroids.fasta"
    file_singletons = "#{@options[:dir_out]}/denoised/singletons.fasta"
    file_fasta      = "#{@options[:dir_out]}/split_library_output/seqs.fna"
    file_map        = "#{@options[:dir_out]}/denoised/denoiser_mapping.txt"
    file_out        = "#{@options[:dir_out]}/denoised/inflated.fna"
    run "inflate_denoiser_output.py -c #{file_centroids} -s #{file_singletons} -f #{file_fasta} -d #{file_map} -o #{file_out}"
  end

  def chimera_check
    if @options[:denoise]
      file_fasta = "#{@options[:dir_out]}/denoised/inflated.fna"
    else
      file_fasta = "#{@options[:dir_out]}/split_library_output/seqs.fna"
    end

    Dir.mkdir("#{@options[:dir_out]}/chimera") unless File.directory? "#{@options[:dir_out]}/chimera"

    file_ref         = @options[:chimera_db]
    file_chimeras    = "#{@options[:dir_out]}/chimera/chimeras.fasta"
    file_nonchimeras = "#{@options[:dir_out]}/chimera/nonchimeras.fasta"

    run "usearch -quiet -uchime #{file_fasta} -db #{file_ref} -chimeras #{file_chimeras} -nonchimeras #{file_nonchimeras}"
  end

  def pick_otus_through_otu_table
    if @options[:chimera]
      file_fasta = "#{@options[:dir_out]}/chimera/nonchimeras.fasta"
    else
      if @options[:denoise]
        file_fasta = "#{@options[:dir_out]}/denoised/inflated.fna"
      else
        file_fasta = "#{@options[:dir_out]}/split_library_output/seqs.fna"
      end
    end

    dir_out = "#{@options[:dir_out]}/otus"
    run "pick_otus_through_otu_table.py -i #{file_fasta} -o #{dir_out} -a -f"
  end

  def per_library_stats
    file_biom  = "#{@options[:dir_out]}/otus/otu_table.biom"
    file_stats = "#{file_biom}.stats"
    run "per_library_stats.py -i #{file_biom} > #{file_stats}"

    File.open(file_stats, 'r') do |ios|
      ios.each_line do |line|
        if line.match(/Min: (\d+)/)
          @min_samples = $1.to_i

          break
        end
      end
    end

    if @min_samples == 0
      error = "Fail: Failed to parse min samples."
      self.send_mail(error) if @options[:email]
      raise error
    end
  end

  def make_otu_heatmap_html
    file_biom = "#{@options[:dir_out]}/otus/otu_table.biom"
    dir_out   = "#{@options[:dir_out]}/otus/OTU_Heatmap"
    run "make_otu_heatmap_html.py -i #{file_biom} -o #{dir_out}"
  end

  def make_otu_network
    file_biom = "#{@options[:dir_out]}/otus/otu_table.biom"
    dir_out   = "#{@options[:dir_out]}/otus/OTU_Network"
    run "make_otu_network.py -m #{@options[:file_map]} -i #{file_biom} -o #{dir_out}"
  end

  # TODO consider adding  -c Description
  def wf_taxa_summary
    file_biom = "#{@options[:dir_out]}/otus/otu_table.biom"
    dir_out   = "#{@options[:dir_out]}/wf_taxa_summary"
    run "summarize_taxa_through_plots.py -i #{file_biom} -o #{dir_out} -m #{@options[:file_map]} -f"
  end

  def alpha_diversity
    file_biom   = "#{@options[:dir_out]}/otus/otu_table.biom"
    file_tree   = "#{@options[:dir_out]}/otus/rep_set.tre"
    file_params = "#{@options[:dir_out]}/alpha_params.txt"
    dir_out     = "#{@options[:dir_out]}/wf_arare"

    File.open(file_params, 'w') do |ios|
      ios.puts "alpha_diversity:metrics shannon,PD_whole_tree,chao1,observed_species"
    end

    run "alpha_rarefaction.py -i #{file_biom} -m #{@options[:file_map]} -o #{dir_out} -p #{file_params} -t #{file_tree} -a -f"
  end

  def beta_diversity_through_plots
    file_biom = "#{@options[:dir_out]}/otus/otu_table.biom"
    file_tree = "#{@options[:dir_out]}/otus/rep_set.tre"
    dir_out   = "#{@options[:dir_out]}/wf_bdiv_even"
    run "beta_diversity_through_plots.py -i #{file_biom} -m #{@options[:file_map]} -o #{dir_out} -t #{file_tree} -e #{@min_samples} -a -f"
  end

  def jackknifed_beta_diversity
    file_biom = "#{@options[:dir_out]}/otus/otu_table.biom"
    file_tree = "#{@options[:dir_out]}/otus/rep_set.tre"
    dir_out   = "#{@options[:dir_out]}/wf_jack"
    samples   = (@min_samples * 0.75).to_i
    run "jackknifed_beta_diversity.py -i #{file_biom} -t #{file_tree} -m #{@options[:file_map]} -o #{dir_out} -e #{samples} -a -f"
  end

  def make_bootstrapped_tree
    file_tree = "#{@options[:dir_out]}/wf_jack/unweighted_unifrac/upgma_cmp/master_tree.tre"
    file_sup  = "#{@options[:dir_out]}/wf_jack/unweighted_unifrac/upgma_cmp/jackknife_support.txt"
    file_out  = "#{@options[:dir_out]}/wf_jack/unweighted_unifrac/upgma_cmp/jackknife_named_nodes.pdf"
    run "make_bootstrapped_tree.py -m #{file_tree} -s #{file_sup} -o #{file_out}"
  end

  def make_3d_plots
    file_weight = "#{@options[:dir_out]}/wf_bdiv_even/unweighted_unifrac_pc.txt"
    file_table  = "#{@options[:dir_out]}/wf_taxa_summary/otu_table_L3.txt"
    dir_out     = "#{@options[:dir_out]}/3d_biplot"
    run "make_3d_plots.py -i #{file_weight} -m #{@options[:file_map]} -t #{file_table} --n_taxa_keep 5 -o #{dir_out}"
  end

  def send_mail(subject)
    cmd = %{cat #{@file_log} | mail -s "#{subject}" #{@options[:email]}}

    system(cmd)

    raise "Command: #{cmd} failed." unless $?.success?
  end

  private

  def run(cmd)
    print "Running: #{cmd} ... "

    if run? cmd
      puts "SKIPPING"
    else
      log "INIT", cmd

      system(cmd)

      if $?.success?
        puts "OK"
        log "OK", cmd
      else
        puts "FAIL"
        log "FAIL", cmd

        self.send_mail("Fail: " + File.basename(@options[:file_sff])) if @options[:email]

        raise "FAIL"
      end
    end
  end

  def run?(cmd)
    if cmd =~ /^sffinfo/   # sffinfo is used multiple times so we use the whole string to check
      cmd_str = cmd
    else
      cmd_str = cmd.split(" ").first
    end

    if File.readable? @file_log
      File.open(@file_log, "r") do |ios|
        ios.each_line do |line|
          return true if line.match(cmd_str) and line =~ /OK$/
        end
      end
    end

    false
  end

  def log(status, cmd)
    time = Time.now

    File.open(@file_log, "a") do |ios|
      ios.puts [time, cmd, status].join("\t")
    end
  end
end

q = Qiime.new(options)
q.log_delete                   if options[:force]
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

