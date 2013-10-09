#!/usr/bin/env ruby

require 'pp'
require 'optparse'
require 'parallel'
require 'maasha/fastq'
require 'maasha/seq/assemble'

ARGV << "-h" if ARGV.empty?

cmd_init = File.basename($0) + " " + ARGV.join(" ")

options = {}

OptionParser.new do |opts|
  opts.banner = "Usage: #{File.basename(__FILE__)} [options]"

  opts.on("-h", "--help", "Display this screen" ) do
    $stderr.puts opts
    exit
  end

  opts.on("--trim_qual <int>", Integer, "Minimum quality (default 20)") do |o|
    options[:trim_qual] = o
  end

  opts.on("--trim_len <int>", Integer, "Minimum stretch length (default 3)") do |o|
    options[:trim_len] = o
  end

  opts.on("--mismatches_max <int>", Integer, "Maximum number of mismatches in percent (default 5)") do |o|
    options[:mismatches_max] = o
  end

  opts.on("--overlap_min <int>", Integer, "Minimum asembly overlap (default 15)") do |o|
    options[:overlap_min] = o
  end

  opts.on("-C", "--cpus <int>", Integer, "Number of CPUs to use (default 1)") do |o|
    options[:cpus] = o
  end
end.parse!

options[:overlap_min]    ||= 15
options[:mismatches_max] ||= 5
options[:trim_qual]      ||= 20
options[:trim_len]       ||= 3
options[:cpus]           ||= 1

fastq_files = ARGV

def sample_names(fastq_files)
  samples = {}
  file1   = nil
  file2   = nil
  prefix1 = nil
  prefix2 = nil

  fastq_files.each do |file|
    case file
    when /(.+)_R1_/ then
      file1   = file
      prefix1 = File.basename $1 
    when /(.+)_R2_/ then
      file2   = file
      prefix2 = File.basename $1
    else
      puts "unmatched file: #{file}"
    end

    if prefix1 == prefix2
      samples[prefix1.to_sym] = {file1: file1, file2: file2}
    end
  end

  samples
end

log = STDERR

samples = sample_names(fastq_files)

log.puts "#" + %w{sample reads_total bases_total reads_assembled_ok reads_assembled_fail bases_assembled bases_ok bases_trimmed }.join("\t")

Parallel.each(samples, in_processes: options[:cpus]) do |sample, files|
  stats = Hash.new(0)
  count = 0

  in1 = Fastq.open(files[:file1])
  in2 = Fastq.open(files[:file2])

  while entry1 = in1.get_entry and entry2 = in2.get_entry
    stats[:reads_total] += 2
    stats[:bases_total] += entry1.length + entry2.length

    entry2.type = :dna

    if assembly = Assemble.pair(entry1, entry2.reverse.complement, options)
      stats[:reads_assembled_ok] += 1
      stats[:bases_assembled] += assembly.length

      trim = assembly.quality_trim(options[:trim_qual], options[:trim_len])

      stats[:bases_ok]      += trim.length
      stats[:bases_trimmed] += assembly.length - trim.length

      trim.seq_name = sample.to_s.sub(/_L\d{3}/, "") + "_#{count} " + trim.seq_name
      puts trim.to_fastq

      count += 1
    else
      stats[:reads_assembled_fail] += 1
    end
  end

  in1.close
  in2.close

  log.puts [
    sample,
    stats[:reads_total],
    stats[:bases_total],
    stats[:reads_assembled_ok],
    stats[:reads_assembled_fail],
    stats[:bases_assembled],
    stats[:bases_ok],
    stats[:bases_trimmed]
  ].join("\t")
end
