#!/usr/bin/env ruby

require 'pp'
require 'optparse'
require 'google_hash'
require 'maasha/fasta'
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

  opts.on("-i", "--input <file>", String, "input file in FASTA format to depreplicate") do |o|
    options[:input] = o
  end

  opts.on("-o", "--output <file>", String, "output file in FASTA format after depreplication") do |o|
    options[:output] = o
  end

  opts.on("-f", "--force", "force overwrite output file") do
    options[:force] = true
  end

  opts.on("-v", "--verbose", "verbose output") do
    options[:verbose] = true
  end
end.parse!

raise OptionParser::MissingArgument, "No input file specified."  unless options[:input]
raise OptionParser::MissingArgument, "No output file specified." unless options[:output]
raise OptionParser::InvalidArgument, "No input file: #{options[:input]}" unless File.file? options[:input]

if File.file? options[:output] and not options[:force]
  raise OptionParser::InvalidArgument, "Output file exists: #{options[:output]} - use --force to overwrite"
end

hash = GoogleHashDenseLongToInt.new

Fasta.open(options[:output], 'w') do |output|
  Fasta.open(options[:input]) do |input|
    input.each_with_index do |entry, i|
      $stderr.puts "Read #{i} entries" if (i % 10000) == 0 and options[:verbose]
      key = entry.seq.upcase.hash
      unless hash[key]
        output.puts entry.to_fasta
        hash[key] = 1
      end
    end
  end
end
