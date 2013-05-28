#!/usr/bin/env ruby

require 'pp'
require 'optparse'
require_relative 'lib/qiime'

ARGV << "-h" if ARGV.empty?

options = {}

OptionParser.new do |opts|
  opts.banner = "Usage: #{File.basename(__FILE__)} [options]"

  opts.on("-h", "--help", "Display this screen" ) do
    $stderr.puts opts
    exit
  end

  opts.on("-f", "--fasta_files <files>", Array, "Fasta files to merge") do |o|
    options[:fasta_files] = o
  end

  opts.on("-o", "--output_file <file>", String, "Output file") do |o|
    options[:output_file] = o
  end
end.parse!

raise OptionParser::MissingArgument, "--fasta_files" if options[:fasta_files].nil?
raise OptionParser::MissingArgument, "--output_file" if options[:output_file].nil?
raise OptionParser::InvalidOption, "Only 1 fasta file specified" if options[:fasta_files].size == 1

file_count = 0

File.open(options[:output_file], 'w') do |o|
  options[:fasta_files].each do |file|
    raise OptionParser::InvalidOption, "no such file: #{file}" unless File.file?(file)

    File.open(file, 'r') do |i|
      i.each do |line|
        line = ">" + file_count.to_s + line[1 .. -1] if line[0] == '>'
        o.puts line
      end
    end

    file_count += 1
  end
end

