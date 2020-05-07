# end join PID sequences using consensus model
# require RUBY package 'viral_seq'

# example:
# ruby end_join.rb [target_directory]


# [target_directory] contains library subdirectories, each has one r1 and one r2 file in it.
# example:
# target_directory
# ├───lib1
# │     lib1_r1.txt
# │     lib1_r2.txt
# ├───lib2
# │     lib1_r2.txt
# │     lib1_r2.txt
# ...


require 'viral_seq'

indir = ARGV[0]

libs = []

Dir.chdir(indir) {libs = Dir.glob('*')}

libs.each do |lib|
  puts "processing #{lib}..."
  regions = []
  Dir.chdir(File.join(indir,lib)) {regions = Dir.glob('*')}
  regions.each do |r|
    puts "\tprocessing #{r}..."
    seq_dir = File.join(indir,lib,r,'consensus')
    sh_pair = ViralSeq::SeqHashPair.fa(seq_dir)
    joined_seq = sh_pair.join2(model: :con)
    joined_file = File.join(seq_dir, "combined.txt")
    joined_seq.write_nt_fa(joined_file)
  end
end
