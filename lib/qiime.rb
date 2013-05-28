module Qiime
  DEFAULT_CHIMERA_DB   = "/home/maasha/Install/QIIME1.6/data/Gold/gold.fa"
  DEFAULT_BARCODE_SIZE = 10
  DEFAULT_CPUS         = 1

  class QiimeError < StandardError; end

  class MapFile
    def initialize
      @mapping_table = []
      @count         = 0
    end

    # Method to parse a QIIME mapping file. The first line starting with '#' is
    # kept as the header while the remaing are skipped.
    def parse_mapping_file(file)
      got_header = false

      File.open(file, 'r') do |ios|
        ios.each do |line|
          if line[0] == '#'
            next if got_header

            raise QiimeError "Mapping file must start with '#SampleID'" unless line =~ /^#SampleID/

            got_header = true
          end

          @mapping_table << line.chomp.split("\t")
        end
      end
    end

    # Method to merge a QIIME mapping file into existing mapping table to
    # @mapping_table. 
    def merge_mapping_file(file)
      got_header = false

      File.open(file, 'r') do |ios|
        ios.each do |line|
          if line[0] == '#'
            next if got_header

            raise QiimeError "Mapping file must start with '#SampleID'" unless line =~ /^#SampleID/

            cols1 = @mapping_table.first.size
            cols2 = line.chomp.split("\t").size
            
            raise QiimeError "Column count differ between mapping files: #{cols1} != #{cols2}" if cols1 != cols2

            got_header = true
          else
            @mapping_table << line.chomp.split("\t")
          end
        end
      end
    end

    # Method to convert a QIIME mapping table to a string.
    def to_s
      table = ""

      @mapping_table.each do |row|
        table << row.join("\t") + $/
      end

      table.chomp
    end
  end

  class Pipeline
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
      run "pick_otus_through_otu_table.py -i #{file_fasta} -o #{dir_out} -a -O #{@options[:cpus]} -f"
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

      @options[:email] = nil
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
end
