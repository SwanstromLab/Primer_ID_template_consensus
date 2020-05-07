# Calculate all postition error rate compared to reference from sample consensus,
# list type of substitution

# require RUBY package 'viral_seq'

# example:
# ruby mut_table.rb [target_directory]

# [target_directory] contains (end-joined) TCS files of the same sequencing region,
# one file per library. Make sure all TCS files are from the same region.
# example directory structure:
# target_directory
# ├───lib1.fasta
# ├───lib2.fasta
# ...



require "viral_seq"
require 'colorize'

mut = []
bases = %w[A T C G]
bases.permutation(2).each do |c|
  m = c[0] + "\s->\s" + c[1]
  mut << m
end

indir = ARGV[0]

if ARGV[1]
  pos = ARGV[1].to_i
else
  pos = 1
end

p = pos

files = []
Dir.chdir(indir){files = Dir.glob("*")}

outdir = indir + "_error"
Dir.mkdir(outdir) unless File.directory?(outdir)
out_file2 = indir + "_type.csv"
out2 = File.open(out_file2,"w")
out2.puts "lib,TCS,Length," + mut.join(",") + ",Total Errors"

merged_seq = ViralSeq::SeqHash.new

files.each do |f|
  infile = indir + "/" + f
  sequences = ViralSeq::SeqHash.fa(infile)
  merged_seq += sequences
end

ref_seq = merged_seq.consensus
l = ref_seq.size
ref_file = indir + "_ref"

out_ref = File.open(ref_file, "w")

out_ref.puts "reference\t" + ref_seq + "\n"
%w{A C G T}.each do |nt|
  out_ref.puts nt + "\t" + ref_seq.count(nt).to_s
end

out_ref.close

files.sort_by{|f| f.split("_")[0].to_i}.each do |f|
  infile = indir + "/" + f
  sh = ViralSeq::SeqHash.fa(infile)
  sequences = sh.dna_hash
  outfile = outdir + "/" + f.split("_")[0] + '.csv'
  out = File.open(outfile, 'w')
  out.puts "ref_poisition,positions,consensus,TCS#,A,C,G,T"

  sequences.each do |name,seq|
    count_e = 0
    (0..(l-1)).each do |n|
      count_e += 1 unless ref_seq[n] == seq[n]
    end
    if count_e > l * 0.02
      puts "This sequences has more than #{(l * 0.01).round.to_s.black} errors per sequence: " + name.magenta
      aligned = ViralSeq::Muscle.align(ref_seq, seq)
      if (aligned[0] + aligned[1]).include?("-")
        puts "It is removed because of indels.".red
        sequences.delete(name)
        next
      else
        puts "No indel detected. It is NOT removed.".blue
      end
    end
  end

  sh.dna_hash = sequences

  sh.error_table.each do |error|
    out.puts p.to_s + ',' + error.join(',')
    p += 1
  end

  p = pos

  pattern = []
  (0..(l-1)).each do |n|
    sequences.each do |_name,seq|
      if ref_seq[n] != seq[n]
        pattern_seq = ref_seq[n] + "\s->\s" + seq[n]
        pattern << pattern_seq
      end
    end
  end

  out.close

  pattern_c = pattern.count_freq

  out2.print f + "," + sequences.size.to_s + "," + l.to_s + ","

  mut.each do |m|
    out2.print pattern_c[m].to_s + ","
  end
  out2.print pattern_c.values.sum.to_s
  out2.print "\n"
end

out2.close
