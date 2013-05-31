module Qiime
  DEFAULT_CHIMERA_DB   = "/home/maasha/Install/QIIME1.6/data/Gold/gold.fa"
  DEFAULT_BARCODE_SIZE = 10
  DEFAULT_CPUS         = 1

  class QiimeError < StandardError; end

  class MapFile
    def initialize
      @mapping_header = []
      @mapping_table  = []
      @file_count     = 0
    end

    # Method to parse and add a QIIME mapping file to the mapping table.
    def <<(file)
      got_header = false

      File.open(file, 'r') do |ios|
        ios.each do |line|
          if line[0] == '#'
            if line =~ /^#SampleID/
              got_header = true

              check_header(line)

              @mapping_header = line.chomp.split("\t") if @mapping_header.empty?
            end
          else
            if got_header
              @mapping_table << (@file_count.to_s + line.chomp).split("\t")
            else
              raise QiimeError, "Mapping file must start with '#SampleID'"
            end
          end
        end
      end

      @file_count += 1
    end

    # Method to convert a QIIME mapping table to a string.
    def to_s
      table = @mapping_header.join("\t") + $/

      @mapping_table.each do |row|
        table << row.join("\t") + $/
      end

      table.chomp
    end

    private

    def check_header(line)
      unless @mapping_header.empty?
        fields = line.chomp.split("\t")

        unless fields.size == @mapping_header.size
          raise QiimeError, "Mapping file column count mismatch"
        end

        line.chomp.split("\t").each_with_index do |field, i|
          unless field == @mapping_header[i]
            raise QiimeError, "Mapping file fields name mismatch: #{field} != #{@mapping_header[i]}"
          end
        end
      end
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

    def log_init(cmd_init)
      unless File.file? @file_log
        File.open(@file_log, "w") do |ios|
          ios.puts "#" + cmd_init
        end
      end
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

    def merge_id_maps
      mapping_files       = @options[:mapping_files].join(",")
      output_file         = "#{@options[:dir_out]}/merged.map"
      @options[:file_map] = output_file
      run "qiime_merge_mapping_files.rb -m #{mapping_files} -o #{output_file}"
    end

    def merge_fasta_files
      fasta_files      = @options[:fasta_files].join(",")
      output_file      = "#{@options[:dir_out]}/merged.fasta"
      @options[:merge] = output_file

      run "qiime_merge_fasta_files.rb -f #{fasta_files} -o #{output_file}"
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

      cmd = "denoise_wrapper.py --titanium -i #{file_sff_txt} -f #{file_fasta} -o #{dir_out} -m #{@options[:file_map]} -n #{@options[:cpus]}"

      interrupted? cmd ? self.denoiser : run(cmd)
    end

    def denoiser
      dir_out      = "#{@options[:dir_out]}/denoised"
      dir_resume   = dir_out + "_resumed"
      base         = File.basename(@options[:file_sff], ".sff")
      file_sff_txt = "#{@options[:dir_out]}/#{base}.sff.txt"
      file_fasta   = "#{@options[:dir_out]}/split_library_output/seqs.fna"

      if checkpoint = get_checkpoint(dir_out)
        run "denoiser.py --titanium -i #{file_sff_txt} -f #{file_fasta} -o #{dir_resume} -p #{dir_out} --checkpoint #{checkpoint} -n #{@options[:cpus]}"
        File.rename(dir_out, dir_out + "_bak")
        File.rename(dir_resume, dir_out)
      else
        raise QiimeError, "No checkpoint found"
      end
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
      if @options[:merge]   # merging datasets
        file_fasta = "#{@options[:dir_out]}/merged.fasta"
      elsif @options[:chimera]
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
        raise QiimeError, error
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

    def wf_taxa_summary
      file_biom = "#{@options[:dir_out]}/otus/otu_table.biom"
      dir_out   = "#{@options[:dir_out]}/wf_taxa_summary"

      if @options[:catagory]
        run "summarize_taxa_through_plots.py -i #{file_biom} -o #{dir_out} -m #{@options[:file_map]} -c #{@options[:catagory]} -f"
      else
        run "summarize_taxa_through_plots.py -i #{file_biom} -o #{dir_out} -m #{@options[:file_map]} -f"
      end
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
      dir_out     = "#{@options[:dir_out]}/3d_biplot"

      if @options[:catagory]
        file_table = "#{@options[:dir_out]}/wf_taxa_summary/#{@options[:catagory]}_otu_table_L3.txt"
      else
        file_table = "#{@options[:dir_out]}/wf_taxa_summary/otu_table_L3.txt"
      end

      run "make_3d_plots.py -i #{file_weight} -m #{@options[:file_map]} -t #{file_table} --n_taxa_keep 5 -o #{dir_out}"
    end

    def send_mail(subject)
      cmd = %{cat #{@file_log} | mail -s "#{subject}" #{@options[:email]}}

      system(cmd)

      raise QiimeError, "Command: #{cmd} failed." unless $?.success?

      @options[:email] = nil
    end

    private

    def run(cmd)
      print "Running: #{cmd} ... "

      if ok? cmd
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

          raise QiimeError, "FAIL"
        end
      end
    end

    def ok?(cmd)
      run_status(cmd) == :OK
    end

    def interrupted?(cmd)
      status = run_status(cmd)
      status == :FAIL or status == :INIT
    end

    def run_status(cmd)
      if cmd =~ /^sffinfo/   # sffinfo is used multiple times so we use the whole string to check
        cmd_str = cmd
      else
        cmd_str = cmd.split(" ").first
      end

      lines = []

      if File.readable? @file_log
        File.open(@file_log, "r") do |ios|
          ios.each_line do |line|
            next if line[0] == '#'
            lines << line.chomp
          end
        end
      end

      lines.reverse.each do |line|
        if line.match(cmd_str)
          return line.split("\t").last.to_sym
        end
      end
    end

    def log(status, cmd)
      time = Time.now

      File.open(@file_log, "a") do |ios|
        ios.puts [time, cmd, status].join("\t")
      end
    end

    def get_checkpoint(dir)
      checkpoints = Dir.glob "#{dir}/checkpoints/*.pickle"

      checkpoints.sort! do |a, b|
        a_num = a.match('\d+')
        b_num = b.match('\d+')

        a_num.to_s.to_i <=> b_num.to_s.to_i
      end
      
      checkpoints.last
    end
  end
end
