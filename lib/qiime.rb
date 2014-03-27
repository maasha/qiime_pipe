module Qiime
  require 'fileutils'

  # DEFAULT_CHIMERA_DB   = "/home/maasha/install/QIIME1.7/data/Gold/gold.fa"
  # DEFAULT_CHIMERA_DB   = "/home/maasha/install/QIIME1.7/data/gg_otus_4feb2011/rep_set/gg_97_otus_4feb2011.fasta"
  DEFAULT_CHIMERA_DB   = "/home/maasha/install/QIIME1.7/data/gg_otus_4feb2011/rep_set/v4_slice_97rep.fasta"
  DEFAULT_BARCODE_SIZE = 10
  DEFAULT_CPUS         = 1

  class QiimeError < StandardError; end

  class ParameterFile
    def initialize(options)
      @options = options
    end

    def expand_relevant_parameters(command)
      ret = ""
      File.open @options[:parameter_file] do |ios|
        ios.each do |line|
          line.chomp!

          script, leftover = line.split ':'
          option, value    = leftover.split /\s+/

          if script == command.chomp(File.extname(command))
            ret = ret + " --" + option + " " + value + " "  
          end
        end
      end
      return ret
      exit
    end
  end

  class MapFile
    def initialize
      @mapping_header = []
      @mapping_table  = []
      @file_count     = 0
    end

    # Method to parse a QIIME mapping file to the mapping table.
    def parse(file)
      got_header = false

      File.open(file, 'r') do |ios|
        ios.each do |line|
          line.gsub!(/\r/, "\n")

          if line[0] == '#'
            if line =~ /^#SampleID/
              got_header = true

              check_header(line)

              @mapping_header = line[1 .. line.length].chomp.split("\t") if @mapping_header.empty?
            end
          else
            if got_header
              @mapping_table << line.chomp.split("\t")
            else
              raise QiimeError, "Mapping file must start with '#SampleID'"
            end
          end
        end
      end

      self
    end

    # Method to parse and add a QIIME mapping file to the mapping table.
    def <<(file)
      got_header = false

      File.open(file, 'r') do |ios|
        ios.each do |line|
          line.gsub!(/\r/, "\n")

          if line[0] == '#'
            if line =~ /^#SampleID/
              got_header = true

              check_header(line)

              @mapping_header = line[1 .. line.length].chomp.split("\t") if @mapping_header.empty?
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
      table = "#" + @mapping_header.join("\t") + $/

      @mapping_table.each do |row|
        table << row.join("\t") + $/
      end

      table.chomp
    end

    # Method to return a mapping table as an array.
    def to_a
      @mapping_header + @mapping_table
    end

    # Method to return the data for a specified column as an array.
    def column(name)
      col  = nil
      data = []

      @mapping_header.each_with_index do |header, i|
        if header.to_sym == name
          col = i
          break
        end
      end

      @mapping_table.each do |row|
        data << row[col]
      end

      data.empty? ? nil : data
    end

    private

    def check_header(line)
      unless @mapping_header.empty?
        line = line[1 .. -1]
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

      if @options[:parameter_file]
        @params = ParameterFile.new(@options)
      end

      if @options[:file_sff]
        @options[:dataset_name] = File.basename(@options[:file_sff])
      elsif @options[:illumina_dirs]
        @options[:dataset_name] = @options[:illumina_dirs].map { |dir| File.basename dir }.join(';')
      else
        @options[:dataset_name] = "unknown dataset"
      end
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

      case ENV['SHELL']
      when /bash/ then run "print_qiime_config.py -t > #{log} 2>&1"
      when /tcsh/ then run "print_qiime_config.py -t >& #{log}"
      else raise "Unknown shell in environment"
      end
    end

    def load_remote_mapping_file
      output_file         = "#{@options[:dir_out]}/mapping_file.txt"
      @options[:file_map] = output_file
      run "load_remote_mapping_file.py -k #{@options[:remote_map]} -o #{output_file}"
    end

    def merge_id_maps
      mapping_files       = @options[:mapping_files].join(",")
      output_file         = "#{@options[:dir_out]}/mapping_file_merged.txt"
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
      run "check_id_map.py -m #{@options[:file_map]} -o #{dir_out} -v"

      log_file = File.join(dir_out, File.basename(@options[:file_map], ".txt") + ".log")

      File.open(log_file) do |ios|
        ios.each do |line|
          if line.chomp == "No errors or warnings found in mapping file."
            break
          else
            error = "Fail: #{@options[:dataset_name]} Errors and warnings found in mapping file."
            self.send_mail(error) if @options[:email]
            raise QiimeError, error
          end
        end
      end
    end

    def process_illumina
      dirs_in = @options[:illumina_dirs].join(',')
      dir_out = "#{@options[:dir_out]}/split_library_output"
      if @options[:trim_primers]
        run "process_illumina.rb --trim_primers -i #{dirs_in} -m #{@options[:file_map]} -o #{dir_out} -C #{@options[:cpus]} -f"
      else
        run "process_illumina.rb -i #{dirs_in} -m #{@options[:file_map]} -o #{dir_out} -C #{@options[:cpus]} -f"
      end
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

    def denoiser_preprocessor
      dir_out      = "#{@options[:dir_out]}/denoiser_preprocessor"
      base         = File.basename(@options[:file_sff], ".sff")
      file_sff_txt = "#{@options[:dir_out]}/#{base}.sff.txt"
      file_fasta   = "#{@options[:dir_out]}/split_library_output/seqs.fna"

      map = MapFile.new
      map << @options[:file_map]

      primer_seq = map.column(:LinkerPrimerSequence).first

      run "denoiser_preprocess.py -i #{file_sff_txt} -f #{file_fasta} -o #{dir_out} -p #{primer_seq}"
    end

    def denoiser
      dir_in       = "#{@options[:dir_out]}/denoiser_preprocessor"
      dir_out      = "#{@options[:dir_out]}/denoised"
      dir_resume   = dir_out + "_resumed"
      base         = File.basename(@options[:file_sff], ".sff")
      file_sff_txt = "#{@options[:dir_out]}/#{base}.sff.txt"
      file_fasta   = "#{@options[:dir_out]}/split_library_output/seqs.fna"

      if interrupted? "denoiser.py"
        if checkpoint = get_checkpoint(dir_resume)
          File.rename(dir_out, dir_out + "." + Time.now.to_f.to_s)

          run "denoiser.py --titanium -i #{file_sff_txt} -f #{file_fasta} -o #{dir_out} -p #{dir_in} --checkpoint #{checkpoint} -c -n #{@options[:cpus]}"
        elsif checkpoint = get_checkpoint(dir_out)
          File.rename(dir_resume, dir_resume + "." + Time.now.to_f.to_s) if File.directory? dir_resume
          run "denoiser.py --titanium -i #{file_sff_txt} -f #{file_fasta} -o #{dir_resume} -p #{dir_in} --checkpoint #{checkpoint} -c -n #{@options[:cpus]}"

          File.rename(dir_out, dir_out + "." + Time.now.to_f.to_s)
          File.rename(dir_resume, dir_out)
        else
          run "denoiser.py --titanium -i #{file_sff_txt} -f #{file_fasta} -o #{dir_out} -p #{dir_in} -c -n #{@options[:cpus]} --force"
        end
      else
        run "denoiser.py --titanium -i #{file_sff_txt} -f #{file_fasta} -o #{dir_out} -p #{dir_in} -c -n #{@options[:cpus]}"
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

    def identify_chimeric_seq
      if @options[:denoise]
        file_fasta = "#{@options[:dir_out]}/denoised/inflated.fna"
      else
        file_fasta = "#{@options[:dir_out]}/split_library_output/seqs.fna"
      end

      Dir.mkdir("#{@options[:dir_out]}/chimera") unless File.directory? "#{@options[:dir_out]}/chimera"

      file_ref      = @options[:chimera_db]
      dir_chimeras  = "#{@options[:dir_out]}/chimera/"

      run "identify_chimeric_seqs.py -m usearch61 -i #{file_fasta} -o #{dir_chimeras} -r #{file_ref} --suppress_usearch61_denovo"  # --suppress_usearch61_ref
    end

    def filter_fasta
      if @options[:denoise]
        file_fasta = "#{@options[:dir_out]}/denoised/inflated.fna"
      else
        file_fasta = "#{@options[:dir_out]}/split_library_output/seqs.fna"
      end

      file_chimeras    = "#{@options[:dir_out]}/chimera/chimeras.txt"
      file_nonchimeras = "#{@options[:dir_out]}/chimera/nonchimeras.fasta"

      run "filter_fasta.py -f #{file_fasta} -o #{file_nonchimeras} -s #{file_chimeras} -n"
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

    def pick_de_novo_otus
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
      picking_method = "uclust"
      alignment_method = "pynast"
      classification_method = "rdp"
      filename = File.basename(file_fasta)
      filename = filename.chomp(File.extname(filename))

      run "pick_otus.py -i #{file_fasta} -m #{picking_method} -o #{dir_out}/#{picking_method}_picked_otus/"
      run "mkdir #{dir_out}/rep_set/"
      run "pick_rep_set.py -i #{dir_out}/#{picking_method}_picked_otus/#{filename}_otus.txt -f #{file_fasta} -o #{dir_out}/rep_set/rep_set.fasta" 
        
        if !@options[:notree]
          run "parallel_align_seqs_pynast.py -i #{dir_out}/rep_set/rep_set.fasta  -o #{dir_out}/#{alignment_method}_aligned/ -O #{@options[:cpus]}" #FIXME: This is not awesome because it limits alignment options to pynast, but align_seqs.py does not support parallelism, and we usually do pynast or no alignments anyway, so the clumsy logic needed to support this is postponed 
          run "filter_alignment.py -i #{dir_out}/#{alignment_method}_aligned/rep_set_aligned.fasta   -o #{dir_out}/#{alignment_method}_aligned/ "
          run "make_phylogeny.py -i #{dir_out}/#{alignment_method}_aligned/rep_set_aligned_pfiltered.fasta -o #{dir_out}/rep_set.tre"
        end
      run "parallel_assign_taxonomy_rdp.py -i #{dir_out}/rep_set/rep_set.fasta  -o #{dir_out}/rdp_assigned_taxonomy/ -O #{@options[:cpus]}" #FIXME: This has the same problem as above but is much more problematic
      run "make_otu_table.py -i #{dir_out}/#{picking_method}_picked_otus/#{filename}_otus.txt -t #{dir_out}/#{classification_method}_assigned_taxonomy/rep_set_tax_assignments.txt -o #{dir_out}/otu_table.biom"
    end

    def print_biom_table_summary
      file_biom  = "#{@options[:dir_out]}/otus/otu_table.biom"
      file_stats = "#{file_biom}.stats"
      run "print_biom_table_summary.py -i #{file_biom} > #{file_stats}"

      File.open(file_stats, 'r') do |ios|
        ios.each_line do |line|
          if line.match(/Min: (\d+)/)
            @min_samples = $1.to_i

            break
          end
        end
      end

      if @min_samples == 0
        error = "Fail: #{@options[:dataset_name]} Failed to parse min samples."
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

      run "alpha_rarefaction.py -i #{file_biom} -m #{@options[:file_map]} -o #{dir_out} -p #{file_params} -t #{file_tree} -a -O #{@options[:cpus]} -f"
    end

    def beta_diversity_through_plots
      file_biom = "#{@options[:dir_out]}/otus/otu_table.biom"
      file_tree = "#{@options[:dir_out]}/otus/rep_set.tre"
      dir_out   = "#{@options[:dir_out]}/wf_bdiv_even"
      run "beta_diversity_through_plots.py -i #{file_biom} -m #{@options[:file_map]} -o #{dir_out} -t #{file_tree} -e #{@min_samples} -a -O #{@options[:cpus]} -f"
    end

    def jackknifed_beta_diversity
      file_biom = "#{@options[:dir_out]}/otus/otu_table.biom"
      file_tree = "#{@options[:dir_out]}/otus/rep_set.tre"
      dir_out   = "#{@options[:dir_out]}/wf_jack"
      samples   = (@min_samples * 0.75).to_i
      run "jackknifed_beta_diversity.py -i #{file_biom} -t #{file_tree} -m #{@options[:file_map]} -o #{dir_out} -e #{samples} -a -O #{@options[:cpus]} -f"
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
      if @params
        #Expands parameters to the commandline 
        #TODO: This needs translation rules that can pass parameters to the parallel version of the programs and vice versa
        program = cmd.split(" ").first
        cmd += @params.expand_relevant_parameters(program)
      end

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

          self.send_mail("Fail: #{@options[:dataset_name]}") if @options[:email]

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
      elsif cmd =~ /^mkdir/  # mkdir may be used multiuple times so same as above
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

      nil
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
        a_num = a.match('(\d+)\.pickle$').to_a.last.to_i
        b_num = b.match('(\d+)\.pickle$').to_a.last.to_i

        a_num <=> b_num
      end
      
      checkpoints.last
    end
  end
end
