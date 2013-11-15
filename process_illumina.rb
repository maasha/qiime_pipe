#!/usr/bin/env ruby

require 'pp'
require 'optparse'
require 'parallel'
require 'maasha/fastq'
require 'maasha/seq/assemble'
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

  opts.on("-i", "--input_dir <dir>", String, "Input directoy with FASTQ files") do |o|
    options[:input_dir] = o
  end

  opts.on("-m", "--map_file <file>", String, "Mapping file to process") do |o|
    options[:map_file] = o
  end

  opts.on("-o", "--output_dir <dir>", String, "Output directoy with FASTQ files (default / )") do |o|
    options[:output_dir] = o
  end

  opts.on("-f", "--force", "Force overwrite output directory") do |o|
    options[:force] = o
  end

  opts.on("-C", "--cpus <int>", Integer, "Number of CPUs to use (default 1)") do |o|
    options[:cpus] = o
  end

  opts.on("--trim_qual <int>", Integer, "Minimum quality (default 20)") do |o|
    options[:trim_qual] = o
  end

  opts.on("--trim_len <int>", Integer, "Minimum stretch length (default 3)") do |o|
    options[:trim_len] = o
  end

  opts.on("--trim_primers", "Trim primers from reads prior to assembly") do |o|
    options[:trim_forward] = o
  end

  opts.on("--min_len <int>", Integer, "Minimum sequence length (default 40)") do |o|
    options[:min_len] = o
  end

  opts.on("--mismatches_max <int>", Integer, "Maximum number of mismatches in percent (default 5)") do |o|
    options[:mismatches_max] = o
  end

  opts.on("--overlap_min <int>", Integer, "Minimum assembly overlap (default 15)") do |o|
    options[:overlap_min] = o
  end
end.parse!

raise OptionParser::MissingArgument, "No input directory specified."  unless options[:input_dir]
raise OptionParser::MissingArgument, "No mapping file specified."     unless options[:map_file]
raise OptionParser::MissingArgument, "No output directory specified." unless options[:output_dir]
raise OptionParser::InvalidArgument, "No no such directory: #{options[:input_dir]}" unless File.directory? options[:input_dir]

if File.directory? options[:output_dir]
  if options[:force]
    FileUtils.rm_rf options[:output_dir]
    FileUtils.mkdir options[:output_dir]
  else
    raise OptionParser::InvalidArgument, "Output directory exists. Use --force to overwrite"
  end
else
  FileUtils.mkdir options[:output_dir]
end

options[:log_dir] = File.join(options[:output_dir], "log")
FileUtils.mkdir options[:log_dir]

options[:seq_dir] = File.join(options[:output_dir], "seq")
FileUtils.mkdir options[:seq_dir]

options[:overlap_min]    ||= 15
options[:mismatches_max] ||= 5
options[:trim_qual]      ||= 20
options[:trim_len]       ||= 3
options[:cpus]           ||= 1
options[:min_len]        ||= 40

name_hash = {}
m = Qiime::MapFile.new
m.parse(options[:map_file]).column(:SampleID).each { |n| name_hash[n] = true }

forward_primer = m.column(:LinkerPrimerSequence).first
reverse_primer = m.column(:ReversePrimerSequence).first

fastq_files = Dir.glob("#{options[:input_dir]}/*").select { |f| name_hash[File.basename(f).split('_').first] }

def sample_names(fastq_files)
  samples = Hash.new { |h, k| h[k] = {} }
  file1   = nil
  file2   = nil
  prefix1 = nil
  prefix2 = nil

  fastq_files.each do |file|
    base = File.basename file

    if base.match(/(.+)_(R[1-2])_/)
      prefix = $1
      pair   = $2

      case pair
      when "R1" then samples[prefix.to_sym][:file1] = file
      when "R2" then samples[prefix.to_sym][:file2] = file
      else raise "Bad pair value: #{pair}"
      end
    else
      raise "Failed to find base prefix: #{base}"
    end
  end

  samples
end

samples = sample_names(fastq_files)

samples.map { |sample| raise "Bad sample #{sample}" if sample.size != 2 }

Parallel.each(samples, in_processes: options[:cpus]) do |sample, files|
  stats = Hash.new(0)
  count = 0

  in1 = Fastq.open(files[:file1])
  in2 = Fastq.open(files[:file2])
  out = Fastq.open(File.join(options[:seq_dir], "#{sample}.fna"), 'w')

  while entry1 = in1.get_entry and entry2 = in2.get_entry
    stats[:reads_total] += 2
    stats[:bases_total] += entry1.length + entry2.length

    if options[:trim_primers]
      if entry1.patmatch_trim_left!(forward_primer, max_mismatches: 2, max_insertions: 1, max_deletions: 1)
        stats[:e1_fprimer_found] += 1

        trim1 = entry1.quality_trim(options[:trim_qual], options[:trim_len])
        stats[:e1_bases_ok]   += trim1.length
        stats[:e1_bases_trim] += entry1.length - trim1.length

        if trim1.length >= options[:min_len]
          stats[:e1_length_ok] += 1

          if entry2.patmatch_trim_left!(reverse_primer, max_mismatches: 2, max_insertions: 1, max_deletions: 1)
            stats[:e2_rprimer_found] += 1

            trim2 = entry2.quality_trim(options[:trim_qual], options[:trim_len])
            stats[:e2_bases_ok]   += trim2.length
            stats[:e2_bases_trim] += entry2.length - trim2.length

            if trim2.length >= options[:min_len]
              stats[:e2_length_ok] += 1
              trim2.type = :dna
              trim2.reverse!.complement!

              if assembly = Assemble.pair(trim1, trim2, options)
                stats[:reads_assembled_ok] += 1
                stats[:bases_assembled] += assembly.length

                assembly.seq_name = sample.to_s.sub(/_S\d+_L\d{3}/, "") + "_#{count} " + assembly.seq_name
                out.puts assembly.to_fasta

                count += 1
              else
                stats[:reads_assembled_fail] += 1
              end
            else
              stats[:e2_length_bad] += 1
            end
          else
            stats[:e2_rprimer_miss] += 1
          end
        else
          stats[:e1_length_bad] += 1
        end
      else
        stats[:e1_fprimer_miss] += 1
      end
    else
      trim1 = entry1.quality_trim(options[:trim_qual], options[:trim_len])
      stats[:e1_bases_ok]   += trim1.length
      stats[:e1_bases_trim] += entry1.length - trim1.length

      if trim1.length >= options[:min_len]
        stats[:e1_length_ok] += 1

        trim2 = entry2.quality_trim(options[:trim_qual], options[:trim_len])
        stats[:e2_bases_ok]   += trim2.length
        stats[:e2_bases_trim] += entry2.length - trim2.length

        if trim2.length >= options[:min_len]
          stats[:e2_length_ok] += 1
          trim2.type = :dna
          trim2.reverse!.complement!

          if assembly = Assemble.pair(trim1, trim2, options)
            stats[:reads_assembled_ok] += 1
            stats[:bases_assembled] += assembly.length

            assembly.seq_name = sample.to_s.sub(/_S\d+_L\d{3}/, "") + "_#{count} " + assembly.seq_name
            out.puts assembly.to_fasta

            count += 1
          else
            stats[:reads_assembled_fail] += 1
          end
        else
          stats[:e2_length_bad] += 1
        end
      else
        stats[:e1_length_bad] += 1
      end
    end
  end

  in1.close
  in2.close
  out.close

  File.open(File.join(options[:log_dir], "#{sample}.log"), 'w') do |ios|
    ios.puts [
      sample,
      stats[:e1_fprimer_found],
      stats[:e1_fprimer_miss],
      stats[:e1_bases_ok],
      stats[:e1_bases_trim],
      stats[:e1_length_ok],
      stats[:e1_length_bad],
      stats[:e2_rprimer_found],
      stats[:e2_rprimer_miss],
      stats[:e2_bases_ok],
      stats[:e2_bases_trim],
      stats[:e2_length_ok],
      stats[:e2_length_bad],
      stats[:reads_assembled_ok],
      stats[:reads_assembled_fail],
      stats[:bases_assembled]
    ].join("\t")
  end
end

system("cat #{options[:seq_dir]}/* > #{options[:output_dir]}/seqs.fna")

stats = []

log_files = Dir.glob("#{options[:log_dir]}/*")

log_files.each do |file|
  File.open(file) do |ios|
    stats << ios.gets
  end
end

File.open(File.join(options[:output_dir], "log.txt"), 'w') do |ios|
  ios.puts "#" + %w{sample
    e1_fprimer_found
    e1_fprimer_miss
    e1_bases_ok
    e1_bases_trim
    e1_length_ok
    e1_length_bad
    e2_rprimer_found
    e2_rprimer_miss
    e2_bases_ok
    e2_bases_trim
    e2_length_ok
    e2_length_bad
    reads_assembled_ok
    reads_assembled_fail
    bases_assembled
  }.join("\t")
  stats.each { |s| ios.puts s }
end

FileUtils.rm_rf options[:seq_dir]
FileUtils.rm_rf options[:log_dir]



__END__

