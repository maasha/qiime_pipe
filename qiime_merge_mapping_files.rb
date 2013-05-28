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

  opts.on("-m", "--mapping_files <files>", Array, "Mapping files to merge") do |o|
    options[:mapping_files] = o
  end

  opts.on("-o", "--output_file <file>", String, "Output file") do |o|
    options[:output_file] = o
  end
end.parse!

raise OptionParser::MissingArgument, "--mapping_files" if options[:mapping_files].nil?
raise OptionParser::MissingArgument, "--output_file"   if options[:output_file].nil?
raise OptionParser::InvalidOption, "Only 1 mapping file specified" if options[:mapping_files].size == 1

m = Qiime::MapFile.new

options[:mapping_files].each do |file|
  raise OptionParser::InvalidOption, "no such file: #{file}" unless File.file?(file)

  m << file
end

File.open options[:output_file], 'w' do |ios|
  ios.puts m
end
