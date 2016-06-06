=begin
Version 1.20-05JUN2016
Create Primer ID template consensus sequences from raw MiSeq FASTq file
Input = directory of raw sequences of two ends (R1 and R2 fasta files, unzipped)
Require parameters:
  list of Primer Sequence of cDNA primer and 1st round PCR forward Primer, including a tag for the pair name
  ignore the first nucleotide of Primer ID: Yes/No
=end
ver = "1.20-05JUN2016"
#############Patch Note#############
=begin
1.Now allow multiplexed Primer ID sequencing system. Input primers in pairs for all sets.
2.Add option to ignore the 1st nucleotide of the Primer ID. 
=end



#############cDNA gene-specific region and forward primer region needs to be defined#######

#mutilple cDNA primers
primers = {}

#change set_name, forward primer sequence and cDNA primer sequence. both forward primer sequence and cDNA primer sequence should include the entire sequence, not just biological sequence
#primers["set_name"] = ["forward primer sequence", "cDNA primer sequence"]
#example of primer 
primers["V1V3"] = ["GCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAGNNNNTTATGGGATCAAAGCCTAAAGCCATGTGTA","GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNNNNNCAGTCCATTTTGCTCTACTAATGTTACAATGTGC"]

#ignore the first nucleotide of the PID, default value true
$ignore_first_nt = true

#input file is the directory containing sequences from both ends of one library
indir = ARGV[0]

#####################General Methods

#convert array to hash in a memory-saving way
def array_to_hash(array)
  count = 0
  hash = Hash.new
  (array.length / 2).times do
    hash[array[count]] = array[count+1]
    count += 2
  end
  return hash
end

#count frequencies of elements in a array.
def count(array)
  hash = Hash.new(0)
  array.each do |element|
    hash[element] +=1
  end
  return hash
end

#calculate consensus cutoff, based on 8-nucleotide long
def calculate_cut_off(m)
  n = 0
  if m <= 10
    n = 2
  elsif m <= 8500
    n = -1.24*10**-21*m**6 + 3.53*10**-17*m**5 - 3.90*10**-13*m**4 + 2.12*10**-9*m**3 - 6.06*10**-6*m**2 + 1.80*10**-2*m + 3.15
  else
    n = 0.0079 * m + 9.4869
  end
  n = n.round
  n = 2 if n < 3
  return n
end

#obtain a consensus sequences
def consensus_without_alignment(seq_array,gap_treatment = 1)
  length = seq_array[0].size
  consensus_bases = []
  (0..(length-1)).each do |n|
    bases = []
    seq_array.each do |seq|
      bases << seq[n]
    end
    if gap_treatment == 1
      consensus_bases << creat_consensus_base_non_gap(bases)
    else
      consensus_bases << creat_consensus_base_gap(bases)
    end
  end
  consensus_seq = consensus_bases.join('')
end

#create a consensus base call at a position. 
def creat_consensus_base_non_gap(base_array_input)
  base_array = Array.new(base_array_input)
  consensus_base = '-'
  number_of_bases = base_array.size
  h = Hash.new(0)
  if base_array.size >0
    base_array.each do |base|
      h[base] += 1
    end
    max_number = h.values.max
    max_list = []
    h.each do |k,v|
      if v == max_number
        max_list << k
      end
    end
    maxi_list_size = max_list.size
    if maxi_list_size == 1
      consensus_base = max_list.shift
    elsif maxi_list_size >= 3
      consensus_base = "N"
    elsif maxi_list_size == 2
      if max_list.include?("A") and max_list.include?("T")
        consensus_base = "W"
      elsif max_list.include?("A") and max_list.include?("C")
        consensus_base = "M"
      elsif max_list.include?("A") and max_list.include?("G")
        consensus_base = "R"
      elsif max_list.include?("T") and max_list.include?("C")
        consensus_base = "Y"
      elsif max_list.include?("G") and max_list.include?("C")
        consensus_base = "S"
      elsif max_list.include?("T") and max_list.include?("G")
        consensus_base = "K"
      elsif max_list.include?('-')
        max_list.delete('-')
        consensus_base = max_list.shift
      end
    end
  end
  return consensus_base.chr
end

module Enumerable
  def median     
    len = self.length
    sorted = self.sort
    median = len % 2 == 1 ? sorted[len/2] : (sorted[len/2 - 1] + sorted[len/2]).to_f / 2
  end 
  
  def sum
     self.inject(0){|accum, i| accum + i }
  end
  
  def mean
    self.sum/self.length.to_f
  end
  
  def sample_variance
    m = self.mean
    sum = self.inject(0){|accum, i| accum + (i-m)**2 }
    sum/(self.length - 1).to_f
  end
  
  def stdev
    return Math.sqrt(self.sample_variance)
  end

end

#####################End of General Methods

#obtain files for two ends for the input directory
indir = ARGV[0]
libname = File.basename(indir)

files = []
Dir.chdir(indir) do
  files = Dir.glob("*")
end
r1_f = ""
r2_f = ""
files.each do |f|
  if f =~ /R1_001\.fastq/
    r1_f = indir + "/" + f
  elsif f =~ /R2_001\.fastq/
    r2_f = indir + "/" + f
  end
end

t = Time.now
outdir = indir + "/consensus_out" + "_" + t.year.to_s + "_" + t.month.to_s + "_" + t.day.to_s + "_" + t.hour.to_s + "_" + t.min.to_s
Dir.mkdir(outdir) unless File.directory?(outdir)



primers.each do |setname,primer_pair|
  puts "Processing " + setname
  n_all_seq = 0
  n_filter_r1 = 0
  n_filter_r2 = 0
  n_paired = 0
  forward_primer = primer_pair[0]
  reverse_primer = primer_pair[1]
  forward_primer.match(/(N+)(\w+)$/)
  forward_n = $1.size
  $forward_bio_primer = $2
  forward_bio_primer_size = $forward_bio_primer.size
  forward_starting_number = forward_n + forward_bio_primer_size
  reverse_primer.match(/(N+)(\w+)$/)
  reverse_n = $1.size
  $reverse_bio_primer = $2
  reverse_bio_primer_size = $reverse_bio_primer.size
  $ignore_first_nt ? id_l = reverse_n - 1 : id_l = reverse_n
  reverse_starting_number = reverse_n + reverse_bio_primer_size
  
  def filter_r2(input_file,id_l=8)
    ref = $reverse_bio_primer
    l = ref.size
    count = 0
    sequence_a = []
    sequence_h = {}
    
    File.open(input_file,'r') do |file|
      file.readlines.collect do |line|
        count +=1
        count_m = count % 4
        if count_m == 1
          line.tr!('@','>')
          sequence_a << line.chomp
        elsif count_m == 2
          sequence_a << line.chomp
        end
      end
    end
    
    sequence_h = array_to_hash(sequence_a)
    sequence_passed = {}
    sequence_h.each do |name,seq|
      next if seq[1..-2] =~ /N/
      next if seq =~ /A{11}/
      next if seq =~ /T{11}/
     
      primer = seq[($ignore_first_nt ? id_l + 1 : id_l),l]
      if primer == ref
        sequence_passed[name] = seq
      end
    end
    return sequence_passed
  end
  
  def filter_r1(input_file)
    ref = $forward_bio_primer
    l = ref.size
    count = 0
    sequence_a = []
    sequence_h = {}
    
    File.open(input_file,'r') do |file|
      file.readlines.collect do |line|
        count +=1
        count_m = count % 4
        if count_m == 1
          line.tr!('@','>')
          sequence_a << line.chomp
        elsif count_m == 2
          sequence_a << line.chomp
        end
      end
    end
    
    sequence_h = array_to_hash(sequence_a)
    n = sequence_h.size
    sequence_passed = {}
    sequence_h.each do |name,seq|
      next if seq[1..-2] =~ /N/
      next if seq =~ /A{11}/
      next if seq =~ /T{11}/
      
      primer = seq[4,l]
      if primer == ref
        sequence_passed[name] = seq
      end
    end
    return [sequence_passed,n]
  end
  
  puts "Filtering R1...."
  r1_temp = filter_r1(r1_f)
  
  filtered_r1_h = r1_temp[0]
  n_all_seq = r1_temp[1]
  print "The number of raw sequences is #{n_all_seq.to_s}\n"
  
  n_filter_r1 = filtered_r1_h.size
  puts "Filtering R2...."
  filtered_r2_h = filter_r2(r2_f,id_l)
  n_filter_r2 = filtered_r2_h.size
  
  print "R1: #{n_filter_r1}\n"
  print "R2: #{n_filter_r2}\n"
  
  puts "Pairing...."
  sequence_rtag1 = {}
  sequence_rtag2 = {}
  
  filtered_r1_h.each do |k,v|
    k =~ /\s/
    k2 = $`
    sequence_rtag1[k2]= v
  end
  
  filtered_r2_h.each do |k,v|
    k =~ /\s/
    k2 = $`
    sequence_rtag2[k2]= v
  end
  
  keys = sequence_rtag1.keys & sequence_rtag2.keys
  
  paired_r1 = {}
  paired_r2 = {}
  
  keys.each do |k|
    paired_r1[k] = sequence_rtag1[k]
    paired_r2[k] = sequence_rtag2[k]
  end
  
  n_paired = keys.size
  puts "Paired raw sequences are : #{n_paired.to_s}"
  
  #create a temp file. Temp file contains sequence names, primer ids, and sequences from two ends
  puts "Create Temp File...."
  temp_out = indir + "/temp_seq"
  Dir.mkdir(temp_out) unless File.directory?(temp_out)
  temp_file = temp_out + "/temp_file_" + setname
  temp_file_out = File.open(temp_file,'w')
  
  #building hashes for Primer ID, and two end sequences
  id = {}
  bio_id = {}
  bio_non_id = {}
  
  paired_r2.each do |k,r2_seq|
    r1 = paired_r1[k]
    id[k] = r2_seq[1,id_l]
    bio_id[k] = r2_seq[reverse_starting_number..-2]
    bio_non_id[k] = r1[forward_starting_number..-2]
    temp_file_out.print k+ "\n" + id[k] + "\n" + bio_id[k] + "\n"+bio_non_id[k] + "\n"
  end
  temp_file_out.close
  
  #hashes of Primer ID list and Primer ID distribution 
  primer_id_list = {}
  primer_id_dis = {}  
  
  puts "Calculate consensus cutoff...."
  #count primer ID
  primer_id_list = id.values
  primer_id_count = count(primer_id_list)
  #Primer ID distribution
  primer_id_dis = count(primer_id_count.values)
  primer_id_in_use = {}
  
 #calculate distinct_to_raw
  distinct_to_raw = (primer_id_count.size/primer_id_list.size.to_f).round(3)
  #define consensus cutoff
  #in case very little raw sequences, i.e. less than 5 unique PIDs. ignore this set and move to the next set.
  if primer_id_dis.keys.size < 5
    File.unlink(temp_file)
    next
  end
  max_id = primer_id_dis.keys.sort[-5..-1].mean 
  n = calculate_cut_off(max_id)
  puts "Consensus cutoff is #{n}"
  puts "Creating consensus..."
  
  #Pick primer ID over threshold n
  primer_id_count_over_n = []
  primer_id_count.each do |primer_id,count|
    primer_id_count_over_n << primer_id if count > n
  end
  nn = primer_id_count_over_n.size
  puts "Number of consensus to process: #{nn}"

  #output part 1
  out_dir_set = outdir + "/" + setname
  Dir.mkdir(out_dir_set)
  out_dir_consensus = out_dir_set + "/consensus" 
  Dir.mkdir(out_dir_consensus)
  
  outfile_id = out_dir_consensus + "/r2.txt"
  outfile_non_id = out_dir_consensus + "/r1.txt"
  
  f1 = File.open(outfile_id,'w')
  f2 = File.open(outfile_non_id,'w')
  
  outdir_primer_id = out_dir_set + "/primer_id"
  Dir.mkdir(outdir_primer_id)
  
  outfile_primer_id_count = outdir_primer_id + "/primer_id_count"
  outfile_primer_id_dis = outdir_primer_id + "/primer_id_dis"
  outfile_primer_id_in_use = outdir_primer_id + "/primer_id_in_use"
  
  f3 = File.open(outfile_primer_id_count,'w')
  f4 = File.open(outfile_primer_id_dis,'w')
  f5 = File.open(outfile_primer_id_in_use,'w')
  
  f3.print "Primer ID List and Counts\n\n"
  f3.print "Primer ID\tCounts\n"
  
  primer_id_count.each do |k,v|
    f3.print k + "\t" + v.to_s + "\n"
  end
  f3.close
  
  f4.print "Primer ID Frequence\n\n"
  f4.print "Frequence\tCounts\n"
  primer_id_dis.keys.sort.each do |c|
    w = primer_id_dis[c]
    f4.print c.to_s + "\t" + w.to_s + "\n"
  end
  f4.close
  #output part 2

  #List of sequences with same primer ID over n.Create consensus
  consensus = {}
  m = 0
  primer_id_count_over_n.each do |primer_id|
    m += 1
    puts "Now processing number #{m}" if m%100 == 0
    seq_with_same_primer_id = []
    id.each do |seq_name,primer_id_test|
      seq_with_same_primer_id << seq_name if primer_id_test == primer_id
    end
    list_id_part = []
    list_non_id_part = []
    seq_with_same_primer_id.each do |seq_name|
      id_part = bio_id[seq_name]
      non_id_part = bio_non_id[seq_name]
      list_id_part << id_part
      list_non_id_part << non_id_part
    end
    #consensus name including the Primer ID and number of raw sequences of that Primer ID, library name and setname.
    consensus_name = ">" + primer_id + "_" + seq_with_same_primer_id.size.to_s + "_" + libname + "_" + setname
    consensus_id_part = consensus_without_alignment(list_id_part)
    consensus_non_id_part = consensus_without_alignment(list_non_id_part)
    #consensus name including the Primer ID and number of raw sequences of that Primer ID
    next if consensus_id_part =~ /[^ATCG]/
    next if consensus_non_id_part =~ /[^ATCG]/
    #get reverse complement sequence of the R2 region
    consensus_id_part.reverse!.tr!('ATCG','TAGC')
    primer_id_in_use[primer_id] = seq_with_same_primer_id.size
    consensus[consensus_name] = [consensus_id_part,consensus_non_id_part]
  end
  
  n_con = consensus.size
  puts "Number of consensus sequences:\t" + n_con.to_s
  #output part 2
  consensus.each do |seq_name,seq|
    f1.print seq_name + "_r2\n" + seq[0] + "\n"
    f2.print seq_name + "_r1\n" + seq[1] + "\n"
  end
  
  f1.close
  f2.close
  
  
  f5.print "Primer ID used to create consensus\n\n"
  f5.print "Primer ID\tCounts\n"
  primer_id_in_use.each do |k,v|
    f5.print k + "\t" + v.to_s + "\n"
  end
  f5.close
  
  #output log file
  log = out_dir_set + "/log.txt"
  
  log_f = File.open(log,'w')
  
  log_f.print "Primer ID pair-end consensus creator Version #{ver}\n\n"
  
  log_f.print "Primer ID pair-end consensus creator\n\n"
  
  log_f.print "Runtime: #{t}\n\n"
  
  log_f.print "Primer set name:\n#{setname}\n\n"
  
  log_f.puts "Forward primer sequence:\t" + forward_primer
  log_f.puts "Reverse primer sequence:\t" + reverse_primer
  
  log_f.print "\nNumber of Raw Sequences for each end is: #{n_all_seq}\n\n"
  
  log_f.print "Number of R1 passed filtered is: #{n_filter_r1}\n\n"
  
  log_f.print "Number of R2 passed filtered is: #{n_filter_r2}\n\n"
  
  log_f.print "Number of sequences paired is: #{n_paired}\n\n"
  
  log_f.print "The consensus threshold is #{n}.\n\n"
  
  log_f.print "Length of Primer ID is #{id_l.to_s}.\n\n"
  
  log_f.print "The number of consensus sequences process (including ambiguities) is #{nn}\n\n"
  
  log_f.print "The number of consensus sequences is #{n_con}\n\n"
  
  log_f.print "Distinct Primer ID to raw is #{distinct_to_raw}\n\n"

  log_f.print "Resampling Parameter is #{(n_con/nn.to_f).round(3)}\n\n"
  
  if distinct_to_raw > 0.1
    log_f.print "WARNING: NOT ENOUGH RAW SEQUENCES, SAMPLING DEPTH MAY NOT BE REVEALED!!!"
    print "\t\t\t****************************\nWARNING: NOT ENOUGH RAW SEQUENCES, SAMPLING DEPTH MAY NOT BE REVEALED!!!\n\t\t\t****************************\n"
  end
  
  log_f.close
  print `rm -rf #{temp_out}`
end

print `rm -rf #{r1_f}`
print `rm -rf #{r2_f}`