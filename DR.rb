=begin
DR version 1.00-13JUN2017
based on TCS Pipeline Version 1.33-19FEB2017 for HIV-1 Multiplexing Drug Resistance Tesing
Regions include
Protease, RT, INT, V1/V3

Create Primer ID template consensus sequences from raw MiSeq FASTq file
Input = directory of raw sequences of two ends (R1 and R2 fasta files, unzipped)
Require parameters:
  list of Primer Sequence of cDNA primer and 1st round PCR forward Primer, including a tag for the pair name
  ignore the first nucleotide of Primer ID: Yes/No
=end
ver = "1.00-13JUN2017 based on TCS 1.33-19FEB2017"
#############Patch Note#############
=begin
  1. sequence locator filter for regions.
  2. PR r1 and r2 regions match need alignment. 
=end


#ignore the first nucleotide of the PID, default value true, remove the # before use

#$ignore_first_nt = true unless defined? $ignore_first_nt

primers = {}

primers["RT"] = ["GCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAGNNNNGGCCATTGACAGAAGAAAAAATAAAAGC","GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNNNNNNNCAGTCACTATAGGCTGTACTGTCCATTTATC"]
primers["V1V3"] = ["GCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAGNNNNTTATGGGATCAAAGCCTAAAGCCATGTGTA","GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNNNNNNNCAGTCCATTTTGCTCTACTAATGTTACAATGTGC"]
primers["INT"] = ["GCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAGNNNNAAAAGGAGAAGCCATGCATG","GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNNNNNNNATCGAATACTGCCATTTGTACTGC"]
primers["PR"] = ["GCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAGNNNNTCAGAGCAGACCAGAGCCAACAGCCCCA", "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNNNNNCAGTTTAACTTTTGGGCCATCCATTCC"]

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

#calculate consensus cutoff
#error rate = 0.02
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

=begin
#error rate = 0.01
def calculate_cut_off(m)
  n = 0
  if m <= 10
    n = 2
  else
    n = 1.09*10**-26*m**6 + 7.82*10**-22*m**5 - 1.93*10**-16*m**4 + 1.01*10**-11*m**3 - 2.31*10**-7*m**2 + 0.00645*m + 2.872
  end
  n = n.round
  n = 2 if n < 3
  return n
end
=end

=begin
#error rate = 0.005
def calculate_cut_off(m)
  n = 0
  if m <= 10
    n = 2
  else
    n = -9.59*10**-27*m**6 + 3.27*10**-21*m**5 - 3.05*10**-16*m**4 + 1.2*10**-11*m**3 - 2.19*10**-7*m**2 + 0.004044*m + 2.273
  end
  n = n.round
  n = 2 if n < 3
  return n
end
=end

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


#primer with ambiguities to match
def primer_match (primer = "")
  match = ""
  primer.each_char.each do |base|
    base_array = to_list(base)
    if base_array.size == 1
      match += base_array[0]
    else
      pattern = "[" + base_array.join("|") + "]"
      match += pattern
    end
  end
  return match
end

def to_list(base = "")
  list = []
  case base
  when /[A|T|C|G]/
    list << base
  when "W"
    list = ['A','T']
  when "S"
    list = ['C','G']
  when "M"
    list = ['A','C']
  when 'K'
    list = ['G','C']
  when 'R'
    list = ['A','G']
  when 'Y'
    list = ['C','T']
  when 'B'
    list = ['C','G','T']
  when 'D'
    list = ['A','G','T']
  when 'H'
    list = ['A','C','T']
  when 'V'
    list = ['A','C','G']
  when 'N'
    list = ['A','T','C','G']
  end
  return list
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


#compare PID with sequences which have identical sequences. 
#PIDs differ by 1 base will be recognized. 
#if PID1 is x time greater than PID2, PID2 will be disgarded

def filter_similar_pid(sequence_hash = {}, cutoff = 10)
  seq = sequence_hash
  uni_seq = seq.values.uniq
  uni_seq_pid = {}
  uni_seq.each do |k|
    seq.each do |name,s|
      name = name[1..-1]
      if k == s
        if uni_seq_pid[k]
          uni_seq_pid[k] << [name.split("_")[0],name.split("_")[1]]
        else
          uni_seq_pid[k] = []
          uni_seq_pid[k] << [name.split("_")[0],name.split("_")[1]]
        end
      end
    end 
  end
  
  dup_pid = []
  uni_seq_pid.values.each do |v|
    next if v.size == 1
    pid_hash = Hash[v]
    list = pid_hash.keys
    list2 = Array.new(list)
    pairs = []
  
    list.each do |k|
      list2.delete(k)
      list2.each do |k1|
        pairs << [k,k1]
      end
    end
    
    pairs.each do |p|
      pid1 = p[0]
      pid2 = p[1]
      if two_pid_x_base_different(pid1,pid2,1)
        n1 = pid_hash[pid1].to_i
        n2 = pid_hash[pid2].to_i
        if n1 >= cutoff * n2
          dup_pid << pid2
          #puts pid1 + "\t" + n1.to_s + "\t" + pid2 + "\t" + n2.to_s
        elsif n2 >= cutoff * n1
          dup_pid << pid1
          #puts pid2 + "\t" + n2.to_s + "\t" + pid1 + "\t" + n1.to_s
        end
      end
    end
  end

  new_seq = {}
  seq.each do |name,s|
    pid = name.split("_")[0][1..-1]
    unless dup_pid.include?(pid)
      new_seq[name] = s
    end
  end
  return new_seq
end

#compare two primer ID sequences. If they differ in x base, return boolean value "TURE", else, return boolean value "FALSE"
def two_pid_x_base_different(pid1="",pid2="", x=0)
  l = pid1.size
  m = l - x
  n = 0
  if pid1.size != pid2.size
    return false
  else
    (0..(pid1.size - 1)).each do |k|
      if pid1[k] == pid2[k]
        n += 1
      end
    end
    if n >= m
      return true
    else
      return false
    end
  end
end

#sequence locator:
#align with HXB2, return location, percentage of similarity
def sequence_locator(seq="",temp_dir=File.dirname($0))

  muscledir = "/nas02/apps/muscle-3.8.31/bin/muscle"

  hxb2_ref = "TGGAAGGGCTAATTCACTCCCAACGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGGGCCAGGGATCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGAGAAGTTAGAAGAAGCCAACAAAGGAGAGAACACCAGCTTGTTACACCCTGTGAGCCTGCATGGAATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACATGGCCCGAGAGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAGGGACTTTCCGCTGGGGACTTTCCAGGGAGGCGTGGCCTGGGCGGGACTGGGGAGTGGCGAGCCCTCAGATCCTGCATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCAGTGGCGCCCGAACAGGGACCTGAAAGCGAAAGGGAAACCAGAGGAGCTCTCTCGACGCAGGACTCGGCTTGCTGAAGCGCGCACGGCAAGAGGCGAGGGGCGGCGACTGGTGAGTACGCCAAAAATTTTGACTAGCGGAGGCTAGAAGGAGAGAGATGGGTGCGAGAGCGTCAGTATTAAGCGGGGGAGAATTAGATCGATGGGAAAAAATTCGGTTAAGGCCAGGGGGAAAGAAAAAATATAAATTAAAACATATAGTATGGGCAAGCAGGGAGCTAGAACGATTCGCAGTTAATCCTGGCCTGTTAGAAACATCAGAAGGCTGTAGACAAATACTGGGACAGCTACAACCATCCCTTCAGACAGGATCAGAAGAACTTAGATCATTATATAATACAGTAGCAACCCTCTATTGTGTGCATCAAAGGATAGAGATAAAAGACACCAAGGAAGCTTTAGACAAGATAGAGGAAGAGCAAAACAAAAGTAAGAAAAAAGCACAGCAAGCAGCAGCTGACACAGGACACAGCAATCAGGTCAGCCAAAATTACCCTATAGTGCAGAACATCCAGGGGCAAATGGTACATCAGGCCATATCACCTAGAACTTTAAATGCATGGGTAAAAGTAGTAGAAGAGAAGGCTTTCAGCCCAGAAGTGATACCCATGTTTTCAGCATTATCAGAAGGAGCCACCCCACAAGATTTAAACACCATGCTAAACACAGTGGGGGGACATCAAGCAGCCATGCAAATGTTAAAAGAGACCATCAATGAGGAAGCTGCAGAATGGGATAGAGTGCATCCAGTGCATGCAGGGCCTATTGCACCAGGCCAGATGAGAGAACCAAGGGGAAGTGACATAGCAGGAACTACTAGTACCCTTCAGGAACAAATAGGATGGATGACAAATAATCCACCTATCCCAGTAGGAGAAATTTATAAAAGATGGATAATCCTGGGATTAAATAAAATAGTAAGAATGTATAGCCCTACCAGCATTCTGGACATAAGACAAGGACCAAAGGAACCCTTTAGAGACTATGTAGACCGGTTCTATAAAACTCTAAGAGCCGAGCAAGCTTCACAGGAGGTAAAAAATTGGATGACAGAAACCTTGTTGGTCCAAAATGCGAACCCAGATTGTAAGACTATTTTAAAAGCATTGGGACCAGCGGCTACACTAGAAGAAATGATGACAGCATGTCAGGGAGTAGGAGGACCCGGCCATAAGGCAAGAGTTTTGGCTGAAGCAATGAGCCAAGTAACAAATTCAGCTACCATAATGATGCAGAGAGGCAATTTTAGGAACCAAAGAAAGATTGTTAAGTGTTTCAATTGTGGCAAAGAAGGGCACACAGCCAGAAATTGCAGGGCCCCTAGGAAAAAGGGCTGTTGGAAATGTGGAAAGGAAGGACACCAAATGAAAGATTGTACTGAGAGACAGGCTAATTTTTTAGGGAAGATCTGGCCTTCCTACAAGGGAAGGCCAGGGAATTTTCTTCAGAGCAGACCAGAGCCAACAGCCCCACCAGAAGAGAGCTTCAGGTCTGGGGTAGAGACAACAACTCCCCCTCAGAAGCAGGAGCCGATAGACAAGGAACTGTATCCTTTAACTTCCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTTCCCATTAGCCCTATTGAGACTGTACCAGTAAAATTAAAGCCAGGAATGGATGGCCCAAAAGTTAAACAATGGCCATTGACAGAAGAAAAAATAAAAGCATTAGTAGAAATTTGTACAGAGATGGAAAAGGAAGGGAAAATTTCAAAAATTGGGCCTGAAAATCCATACAATACTCCAGTATTTGCCATAAAGAAAAAAGACAGTACTAAATGGAGAAAATTAGTAGATTTCAGAGAACTTAATAAGAGAACTCAAGACTTCTGGGAAGTTCAATTAGGAATACCACATCCCGCAGGGTTAAAAAAGAAAAAATCAGTAACAGTACTGGATGTGGGTGATGCATATTTTTCAGTTCCCTTAGATGAAGACTTCAGGAAGTATACTGCATTTACCATACCTAGTATAAACAATGAGACACCAGGGATTAGATATCAGTACAATGTGCTTCCACAGGGATGGAAAGGATCACCAGCAATATTCCAAAGTAGCATGACAAAAATCTTAGAGCCTTTTAGAAAACAAAATCCAGACATAGTTATCTATCAATACATGGATGATTTGTATGTAGGATCTGACTTAGAAATAGGGCAGCATAGAACAAAAATAGAGGAGCTGAGACAACATCTGTTGAGGTGGGGACTTACCACACCAGACAAAAAACATCAGAAAGAACCTCCATTCCTTTGGATGGGTTATGAACTCCATCCTGATAAATGGACAGTACAGCCTATAGTGCTGCCAGAAAAAGACAGCTGGACTGTCAATGACATACAGAAGTTAGTGGGGAAATTGAATTGGGCAAGTCAGATTTACCCAGGGATTAAAGTAAGGCAATTATGTAAACTCCTTAGAGGAACCAAAGCACTAACAGAAGTAATACCACTAACAGAAGAAGCAGAGCTAGAACTGGCAGAAAACAGAGAGATTCTAAAAGAACCAGTACATGGAGTGTATTATGACCCATCAAAAGACTTAATAGCAGAAATACAGAAGCAGGGGCAAGGCCAATGGACATATCAAATTTATCAAGAGCCATTTAAAAATCTGAAAACAGGAAAATATGCAAGAATGAGGGGTGCCCACACTAATGATGTAAAACAATTAACAGAGGCAGTGCAAAAAATAACCACAGAAAGCATAGTAATATGGGGAAAGACTCCTAAATTTAAACTGCCCATACAAAAGGAAACATGGGAAACATGGTGGACAGAGTATTGGCAAGCCACCTGGATTCCTGAGTGGGAGTTTGTTAATACCCCTCCCTTAGTGAAATTATGGTACCAGTTAGAGAAAGAACCCATAGTAGGAGCAGAAACCTTCTATGTAGATGGGGCAGCTAACAGGGAGACTAAATTAGGAAAAGCAGGATATGTTACTAATAGAGGAAGACAAAAAGTTGTCACCCTAACTGACACAACAAATCAGAAGACTGAGTTACAAGCAATTTATCTAGCTTTGCAGGATTCGGGATTAGAAGTAAACATAGTAACAGACTCACAATATGCATTAGGAATCATTCAAGCACAACCAGATCAAAGTGAATCAGAGTTAGTCAATCAAATAATAGAGCAGTTAATAAAAAAGGAAAAGGTCTATCTGGCATGGGTACCAGCACACAAAGGAATTGGAGGAAATGAACAAGTAGATAAATTAGTCAGTGCTGGAATCAGGAAAGTACTATTTTTAGATGGAATAGATAAGGCCCAAGATGAACATGAGAAATATCACAGTAATTGGAGAGCAATGGCTAGTGATTTTAACCTGCCACCTGTAGTAGCAAAAGAAATAGTAGCCAGCTGTGATAAATGTCAGCTAAAAGGAGAAGCCATGCATGGACAAGTAGACTGTAGTCCAGGAATATGGCAACTAGATTGTACACATTTAGAAGGAAAAGTTATCCTGGTAGCAGTTCATGTAGCCAGTGGATATATAGAAGCAGAAGTTATTCCAGCAGAAACAGGGCAGGAAACAGCATATTTTCTTTTAAAATTAGCAGGAAGATGGCCAGTAAAAACAATACATACTGACAATGGCAGCAATTTCACCGGTGCTACGGTTAGGGCCGCCTGTTGGTGGGCGGGAATCAAGCAGGAATTTGGAATTCCCTACAATCCCCAAAGTCAAGGAGTAGTAGAATCTATGAATAAAGAATTAAAGAAAATTATAGGACAGGTAAGAGATCAGGCTGAACATCTTAAGACAGCAGTACAAATGGCAGTATTCATCCACAATTTTAAAAGAAAAGGGGGGATTGGGGGGTACAGTGCAGGGGAAAGAATAGTAGACATAATAGCAACAGACATACAAACTAAAGAATTACAAAAACAAATTACAAAAATTCAAAATTTTCGGGTTTATTACAGGGACAGCAGAAATCCACTTTGGAAAGGACCAGCAAAGCTCCTCTGGAAAGGTGAAGGGGCAGTAGTAATACAAGATAATAGTGACATAAAAGTAGTGCCAAGAAGAAAAGCAAAGATCATTAGGGATTATGGAAAACAGATGGCAGGTGATGATTGTGTGGCAAGTAGACAGGATGAGGATTAGAACATGGAAAAGTTTAGTAAAACACCATATGTATGTTTCAGGGAAAGCTAGGGGATGGTTTTATAGACATCACTATGAAAGCCCTCATCCAAGAATAAGTTCAGAAGTACACATCCCACTAGGGGATGCTAGATTGGTAATAACAACATATTGGGGTCTGCATACAGGAGAAAGAGACTGGCATTTGGGTCAGGGAGTCTCCATAGAATGGAGGAAAAAGAGATATAGCACACAAGTAGACCCTGAACTAGCAGACCAACTAATTCATCTGTATTACTTTGACTGTTTTTCAGACTCTGCTATAAGAAAGGCCTTATTAGGACACATAGTTAGCCCTAGGTGTGAATATCAAGCAGGACATAACAAGGTAGGATCTCTACAATACTTGGCACTAGCAGCATTAATAACACCAAAAAAGATAAAGCCACCTTTGCCTAGTGTTACGAAACTGACAGAGGATAGATGGAACAAGCCCCAGAAGACCAAGGGCCACAGAGGGAGCCACACAATGAATGGACACTAGAGCTTTTAGAGGAGCTTAAGAATGAAGCTGTTAGACATTTTCCTAGGATTTGGCTCCATGGCTTAGGGCAACATATCTATGAAACTTATGGGGATACTTGGGCAGGAGTGGAAGCCATAATAAGAATTCTGCAACAACTGCTGTTTATCCATTTTCAGAATTGGGTGTCGACATAGCAGAATAGGCGTTACTCGACAGAGGAGAGCAAGAAATGGAGCCAGTAGATCCTAGACTAGAGCCCTGGAAGCATCCAGGAAGTCAGCCTAAAACTGCTTGTACCAATTGCTATTGTAAAAAGTGTTGCTTTCATTGCCAAGTTTGTTTCATAACAAAAGCCTTAGGCATCTCCTATGGCAGGAAGAAGCGGAGACAGCGACGAAGAGCTCATCAGAACAGTCAGACTCATCAAGCTTCTCTATCAAAGCAGTAAGTAGTACATGTAACGCAACCTATACCAATAGTAGCAATAGTAGCATTAGTAGTAGCAATAATAATAGCAATAGTTGTGTGGTCCATAGTAATCATAGAATATAGGAAAATATTAAGACAAAGAAAAATAGACAGGTTAATTGATAGACTAATAGAAAGAGCAGAAGACAGTGGCAATGAGAGTGAAGGAGAAATATCAGCACTTGTGGAGATGGGGGTGGAGATGGGGCACCATGCTCCTTGGGATGTTGATGATCTGTAGTGCTACAGAAAAATTGTGGGTCACAGTCTATTATGGGGTACCTGTGTGGAAGGAAGCAACCACCACTCTATTTTGTGCATCAGATGCTAAAGCATATGATACAGAGGTACATAATGTTTGGGCCACACATGCCTGTGTACCCACAGACCCCAACCCACAAGAAGTAGTATTGGTAAATGTGACAGAAAATTTTAACATGTGGAAAAATGACATGGTAGAACAGATGCATGAGGATATAATCAGTTTATGGGATCAAAGCCTAAAGCCATGTGTAAAATTAACCCCACTCTGTGTTAGTTTAAAGTGCACTGATTTGAAGAATGATACTAATACCAATAGTAGTAGCGGGAGAATGATAATGGAGAAAGGAGAGATAAAAAACTGCTCTTTCAATATCAGCACAAGCATAAGAGGTAAGGTGCAGAAAGAATATGCATTTTTTTATAAACTTGATATAATACCAATAGATAATGATACTACCAGCTATAAGTTGACAAGTTGTAACACCTCAGTCATTACACAGGCCTGTCCAAAGGTATCCTTTGAGCCAATTCCCATACATTATTGTGCCCCGGCTGGTTTTGCGATTCTAAAATGTAATAATAAGACGTTCAATGGAACAGGACCATGTACAAATGTCAGCACAGTACAATGTACACATGGAATTAGGCCAGTAGTATCAACTCAACTGCTGTTAAATGGCAGTCTAGCAGAAGAAGAGGTAGTAATTAGATCTGTCAATTTCACGGACAATGCTAAAACCATAATAGTACAGCTGAACACATCTGTAGAAATTAATTGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGAGAGCATTTGTTACAATAGGAAAAATAGGAAATATGAGACAAGCACATTGTAACATTAGTAGAGCAAAATGGAATAACACTTTAAAACAGATAGCTAGCAAATTAAGAGAACAATTTGGAAATAATAAAACAATAATCTTTAAGCAATCCTCAGGAGGGGACCCAGAAATTGTAACGCACAGTTTTAATTGTGGAGGGGAATTTTTCTACTGTAATTCAACACAACTGTTTAATAGTACTTGGTTTAATAGTACTTGGAGTACTGAAGGGTCAAATAACACTGAAGGAAGTGACACAATCACCCTCCCATGCAGAATAAAACAAATTATAAACATGTGGCAGAAAGTAGGAAAAGCAATGTATGCCCCTCCCATCAGTGGACAAATTAGATGTTCATCAAATATTACAGGGCTGCTATTAACAAGAGATGGTGGTAATAGCAACAATGAGTCCGAGATCTTCAGACCTGGAGGAGGAGATATGAGGGACAATTGGAGAAGTGAATTATATAAATATAAAGTAGTAAAAATTGAACCATTAGGAGTAGCACCCACCAAGGCAAAGAGAAGAGTGGTGCAGAGAGAAAAAAGAGCAGTGGGAATAGGAGCTTTGTTCCTTGGGTTCTTGGGAGCAGCAGGAAGCACTATGGGCGCAGCCTCAATGACGCTGACGGTACAGGCCAGACAATTATTGTCTGGTATAGTGCAGCAGCAGAACAATTTGCTGAGGGCTATTGAGGCGCAACAGCATCTGTTGCAACTCACAGTCTGGGGCATCAAGCAGCTCCAGGCAAGAATCCTGGCTGTGGAAAGATACCTAAAGGATCAACAGCTCCTGGGGATTTGGGGTTGCTCTGGAAAACTCATTTGCACCACTGCTGTGCCTTGGAATGCTAGTTGGAGTAATAAATCTCTGGAACAGATTTGGAATCACACGACCTGGATGGAGTGGGACAGAGAAATTAACAATTACACAAGCTTAATACACTCCTTAATTGAAGAATCGCAAAACCAGCAAGAAAAGAATGAACAAGAATTATTGGAATTAGATAAATGGGCAAGTTTGTGGAATTGGTTTAACATAACAAATTGGCTGTGGTATATAAAATTATTCATAATGATAGTAGGAGGCTTGGTAGGTTTAAGAATAGTTTTTGCTGTACTTTCTATAGTGAATAGAGTTAGGCAGGGATATTCACCATTATCGTTTCAGACCCACCTCCCAACCCCGAGGGGACCCGACAGGCCCGAAGGAATAGAAGAAGAAGGTGGAGAGAGAGACAGAGACAGATCCATTCGATTAGTGAACGGATCCTTGGCACTTATCTGGGACGATCTGCGGAGCCTGTGCCTCTTCAGCTACCACCGCTTGAGAGACTTACTCTTGATTGTAACGAGGATTGTGGAACTTCTGGGACGCAGGGGGTGGGAAGCCCTCAAATATTGGTGGAATCTCCTACAGTATTGGAGTCAGGAACTAAAGAATAGTGCTGTTAGCTTGCTCAATGCCACAGCCATAGCAGTAGCTGAGGGGACAGATAGGGTTATAGAAGTAGTACAAGGAGCTTGTAGAGCTATTCGCCACATACCTAGAAGAATAAGACAGGGCTTGGAAAGGATTTTGCTATAAGATGGGTGGCAAGTGGTCAAAAAGTAGTGTGATTGGATGGCCTACTGTAAGGGAAAGAATGAGACGAGCTGAGCCAGCAGCAGATAGGGTGGGAGCAGCATCTCGAGACCTGGAAAAACATGGAGCAATCACAAGTAGCAATACAGCAGCTACCAATGCTGCTTGTGCCTGGCTAGAAGCACAAGAGGAGGAGGAGGTGGGTTTTCCAGTCACACCTCAGGTACCTTTAAGACCAATGACTTACAAGGCAGCTGTAGATCTTAGCCACTTTTTAAAAGAAAAGGGGGGACTGGAAGGGCTAATTCACTCCCAAAGAAGACAAGATATCCTTGATCTGTGGATCTACCACACACAAGGCTACTTCCCTGATTAGCAGAACTACACACCAGGGCCAGGGGTCAGATATCCACTGACCTTTGGATGGTGCTACAAGCTAGTACCAGTTGAGCCAGATAAGATAGAAGAGGCCAATAAAGGAGAGAACACCAGCTTGTTACACCCTGTGAGCCTGCATGGGATGGATGACCCGGAGAGAGAAGTGTTAGAGTGGAGGTTTGACAGCCGCCTAGCATTTCATCACGTGGCCCGAGAGCTGCATCCGGAGTACTTCAAGAACTGCTGACATCGAGCTTGCTACAAGGGACTTTCCGCTGGGGACTTTCCAGGGAGGCGTGGCCTGGGCGGGACTGGGGAGTGGCGAGCCCTCAGATCCTGCATATAAGCAGCTGCTTTTTGCCTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAGCCTGGGAGCTCTCTGGCTAACTAGGGAACCCACTGCTTAAGCCTCAATAAAGCTTGCCTTGAGTGCTTCAAGTAGTGTGTGCCCGTCTGTTGTGTGACTCTGGTAACTAGAGATCCCTCAGACCCTTTTAGTCAGTGTGGAAAATCTCTAGCA"
  hxb2_l = hxb2_ref.size
  temp_file = temp_dir + "/temp"
  temp_aln = temp_dir + "/temp_aln"
  
  l1 = 0
  l2 = 0
  name = ">test"
  temp_in = File.open(temp_file,"w")
  temp_in.puts ">ref"
  temp_in.puts hxb2_ref
  temp_in.puts name
  temp_in.puts seq
  temp_in.close
  
  print `#{muscledir} -in #{temp_file} -out #{temp_aln} -quiet`
  aln_seq = fasta_to_hash(temp_aln)
  aln_test = aln_seq[name]
  aln_test =~ /^(\-*)(\w.*\w)(\-*)$/
  gap_begin = $1.size
  gap_end = $3.size
  aln_test2 = $2
  ref = aln_seq[">ref"]
  ref = ref[gap_begin..(-gap_end-1)]
  ref_size = ref.size
  if ref_size > 1.3*(seq.size)
    l1 = l1 + gap_begin
    l2 = l2 + gap_end
    max_seq = aln_test2.scan(/[ACGT]+/).max_by(&:length)
    aln_test2 =~ /#{max_seq}/
    before_aln_seq = $`
    before_aln = $`.size
    post_aln_seq = $'
    post_aln = $'.size
    before_aln_seq_size = before_aln_seq.scan(/[ACGT]+/).join("").size
    b1 = (1.3 * before_aln_seq_size).to_i
    post_aln_seq_size = post_aln_seq.scan(/[ACGT]+/).join("").size
    b2 = (1.3 * post_aln_seq_size).to_i
    if (before_aln > seq.size) and (post_aln <= seq.size)
      ref = ref[(before_aln - b1)..(ref_size - post_aln - 1)]
      l1 = l1 + (before_aln - b1)
    elsif (post_aln > seq.size) and (before_aln <= seq.size)
      ref = ref[before_aln..(ref_size - post_aln - 1 + b2)]
      l2 = l2 + post_aln - b2
    elsif (post_aln > seq.size) and (before_aln > seq.size)
      ref = ref[(before_aln - b1)..(ref_size - post_aln - 1 + b2)]
      l1 = l1 + (before_aln - b1)
      l2 = l2 + (post_aln - b2)
    end
    temp_in = File.open(temp_file,"w")
    temp_in.puts ">ref"
    temp_in.puts ref
    temp_in.puts name
    temp_in.puts seq
    temp_in.close
    print `#{muscledir} -in #{temp_file} -out #{temp_aln} -quiet`
    aln_seq = fasta_to_hash(temp_aln)
    aln_test = aln_seq[name]
    aln_test =~ /^(\-*)(\w.*\w)(\-*)$/
    gap_begin = $1.size
    gap_end = $3.size
    aln_test2 = $2
    ref = aln_seq[">ref"]
    ref = ref[gap_begin..(-gap_end-1)]
    ref_size = ref.size
  end
  aln_seq = fasta_to_hash(temp_aln)
  aln_test = aln_seq[name]
  aln_test =~ /^(\-*)(\w.*\w)(\-*)$/
  gap_begin = $1.size
  gap_end = $3.size
  aln_test = $2
  aln_test =~ /^(\w+)(\-*)\w/
  s1 = $1.size
  g1 = $2.size
  aln_test =~ /\w(\-*)(\w+)$/
  s2 = $2.size
  g2 = $1.size
  ref = aln_seq[">ref"]
  ref = ref[gap_begin..(-gap_end-1)]

  l1 = l1 + gap_begin
  l2 = l2 + gap_end
  repeat = 0
  
  if g1 == g2 and (s1 + g1 + s2) == ref.size
    if s1 > s2 and g2 > 2*s2
      ref = ref[0..(-g2-1)]
      repeat = 1
      l2 = l2 + g2
    elsif s1 < s2 and g1 > 2*s1
      ref = ref[g1..-1]
      repeat = 1
      l1 = l1 + g1
    end
  else
    if g1 > 2*s1
      ref = ref[g1..-1]
      repeat = 1
      l1 = l1 + g1
    end
    if g2 > 2*s2
      ref = ref[0..(-g2 - 1)]
      repeat = 1
      l2 = l2 + g2
    end 
  end

  while repeat == 1
    temp_in = File.open(temp_file,"w")
    temp_in.puts ">ref"
    temp_in.puts ref
    temp_in.puts name
    temp_in.puts seq
    temp_in.close
    print `#{muscledir} -in #{temp_file} -out #{temp_aln} -quiet`
    aln_seq = fasta_to_hash(temp_aln)
    aln_test = aln_seq[name]
    aln_test =~ /^(\-*)(\w.*\w)(\-*)$/
    gap_begin = $1.size
    gap_end = $3.size
    aln_test = $2
    aln_test =~ /^(\w+)(\-*)\w/
    s1 = $1.size
    g1 = $2.size
    aln_test =~ /\w(\-*)(\w+)$/
    s2 = $2.size
    g2 = $1.size
    ref = aln_seq[">ref"]
    ref = ref[gap_begin..(-gap_end-1)]
    l1 = l1 + gap_begin
    l2 = l2 + gap_end
    repeat = 0
    if g1 > 2*s1
      ref = ref[g1..-1]
      repeat = 1
      l1 = l1 + g1
    end
    if g2 > 2*s2
      ref = ref[0..(-g2 - 1)]
      repeat = 1
      l2 = l2 + g2
    end
  end
  aln_seq = fasta_to_hash(temp_aln)
  aln_test = aln_seq[name]
  aln_test =~ /^(\-*)(\w.*\w)(\-*)$/
  gap_begin = $1.size
  gap_end = $3.size
  aln_test2 = $2
  ref = aln_seq[">ref"]
  ref = ref[gap_begin..(-gap_end-1)]
  ref_size = ref.size
  sim_count = 0
  (0..(ref_size-1)).each do |n|
    ref_base = ref[n]
    test_base = aln_test2[n]
    sim_count += 1 if ref_base == test_base
  end
  similarity = (sim_count/ref_size.to_f*100).round(1)
  print `rm -f #{temp_file}`
  print `rm -f #{temp_aln}`
  return [(l1+1),(hxb2_l-l2),similarity]
end

def fasta_to_hash(infile)
  f=File.open(infile,"r")
  return_hash = {}
  name = ""
  while line = f.gets do
    if line =~ /^\>/
      name = line.chomp
      return_hash[name] = ""
    else
      return_hash[name] += line.chomp
    end
  end
  f.close
  return return_hash
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
  if f =~ /r1/i
    r1_f = indir + "/" + f
  elsif f =~ /r2/i
    r2_f = indir + "/" + f
  end
end

t = Time.now
#outdir = indir + "/consensus_out" + "_" + t.year.to_s + "_" + t.month.to_s + "_" + t.day.to_s + "_" + t.hour.to_s + "_" + t.min.to_s
outdir = indir + "/" + File.basename(indir)
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
    ref = primer_match($reverse_bio_primer)
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
      if primer =~ /#{ref}/
        sequence_passed[name] = seq
      end
    end
    return sequence_passed
  end
  
  def filter_r1(input_file)
    ref = primer_match($forward_bio_primer)
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
      if primer =~ /#{ref}/
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
    id[k] = r2_seq[($ignore_first_nt ? 1 : 0),id_l]
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
  Dir.mkdir(out_dir_set) unless File.directory?(out_dir_set)
  out_dir_consensus = out_dir_set + "/consensus" 
  Dir.mkdir(out_dir_consensus) unless File.directory?(out_dir_consensus)
  
  outfile_id = out_dir_consensus + "/r2.txt"
  outfile_non_id = out_dir_consensus + "/r1.txt"
  outfile_combined = out_dir_consensus + "/combined.txt"
  
  f1 = File.open(outfile_id,'w')
  f2 = File.open(outfile_non_id,'w')
  f6 = File.open(outfile_combined,'w')
  
  outdir_primer_id = out_dir_set + "/primer_id"
  Dir.mkdir(outdir_primer_id) unless File.directory?(outdir_primer_id)
  
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
  
  consensus_filtered = {}
  r1_consensus = {}
  r2_consensus = {}
  consensus.each do |seq_name,seq|
    r1_consensus[seq_name] = seq[1]
    r2_consensus[seq_name] = seq[0]
  end
  consensus_number_temp = consensus.size
  
  max_pid_comb = 4**id_l
  
  if consensus_number_temp < 0.003*max_pid_comb
    puts "Applying PID post consensus filter..."
    r1_consensus_filtered = filter_similar_pid(r1_consensus,10)
    r2_consensus_filtered = filter_similar_pid(r2_consensus,10)
    common_pid = r1_consensus_filtered.keys & r2_consensus_filtered.keys
    common_pid.each do |pid|
      consensus_filtered[pid] = [r2_consensus_filtered[pid],r1_consensus_filtered[pid]]
    end
  else
    consensus_filtered = consensus
  end
  
  n_con = consensus_filtered.size
  puts "Number of consensus sequences:\t" + n_con.to_s
  #output part 3
  
  r1 = {}
  r2 = {}
  consensus_filtered.each do |seq_name,seq|
    r1[seq_name] = seq[1]
    r2[seq_name] = seq[0]
    f1.print seq_name + "_r2\n" + seq[0] + "\n"
    f2.print seq_name + "_r1\n" + seq[1] + "\n"
  end
  
  r1_seq_uni = r1.values.uniq
  r2_seq_uni = r2.values.uniq
  r1_uni_pass = []
  r2_uni_pass = []
  r1_seq_uni.each do |k|
    loc = sequence_locator(k,temp_out)
    case setname
    when "RT"
      next if loc[3]
      if (loc[0] == 2648) and (loc[1] == 2914)
        r1_uni_pass << k
      end
    when "PR"
      if (loc[0] < 2253) and (loc[1] > 2253)
        r1_uni_pass << k
      end 
    when "INT"
      next if loc[3]
      if (loc[0] == 4384) and (loc[1] == 4658)
        r1_uni_pass << k
      end 
    when "V1V3"
      if loc[0] == 6585
        r1_uni_pass << k
      end
    end
  end
  r2_seq_uni.each do |k|
    loc = sequence_locator(k,temp_out)
    case setname
    when "RT"
      next if loc[3]
      if (loc[0] == 3001) and (loc[1] == 3257)
        r2_uni_pass << k
      end
    when "PR"
      next if loc[3]
      if (loc[0] == 2329) and (loc[1] == 2591)
        r2_uni_pass << k
      end      
    when "INT"
      next if loc[3]
      if (loc[0] == 4488) and (loc[1] == 4751)
        r2_uni_pass << k
      end      
    when "V1V3"
      if loc[1] == 7208
        r2_uni_pass << k
      end
    end
    
  end
  r1_pass = {}
  r1_uni_pass.each do |seq|
    r1.each do |k,v|
      if v == seq
        r1_pass[k] = v
      end
    end
  end
  r2_pass = {}
  r2_uni_pass.each do |seq|
    r2.each do |k,v|
      if v == seq
        r2_pass[k] = v
      end
    end
  end
  r1_pass_size = r1_pass.size
  r2_pass_size = r2_pass.size
  
  combined_k = r1_pass.keys & r2_pass.keys
  combined_seq = {}
  
  combined_k.each do |k|
    case setname
    when /V1V3|RT/
      s = r1_pass[k] + r2_pass[k]
      combined_seq[k] = s
    when "PR"
      r1s = r1_pass[k]
      r2s = r2_pass[k]
      loc_r1 = sequence_locator(r1s,temp_out)
      x = loc_r1[1]
      s = r1s[-(x-2252)..-1] + r2s[(x-2328),(2549-x)]
      next unless s
      combined_seq[k] = s
    when "INT"
      s = r1_pass[k] + r2_pass[k][171..-1]
      next unless s
      combined_seq[k] = s
    end
    f6.puts k + "_combined"
    f6.puts s
  end
  f6.close
  combined_size = combined_seq.size
  print "Number of combined sequences: #{combined_size}\n"
  
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
  
  log_f.print "TCS-DR pipeline #{ver}\n\n"
  
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
  
  log_f.print "R1 Location filter: #{r1_pass_size}\n\n"

  log_f.print "R2 Location filter: #{r2_pass_size}\n\n"

  log_f.print "The number of combined sequences is #{combined_size}\n\n"
  
  if distinct_to_raw > 0.1
    log_f.print "WARNING: NOT ENOUGH RAW SEQUENCES, SAMPLING DEPTH MAY NOT BE REVEALED!!!"
    print "\t\t\t****************************\nWARNING: NOT ENOUGH RAW SEQUENCES, SAMPLING DEPTH MAY NOT BE REVEALED!!!\n\t\t\t****************************\n"
  end
  
  log_f.close
  print `rm -rf #{temp_out}`
end

outdir_tar = outdir + ".tar.gz"

if File.exists?(outdir_tar)
  File.unlink(outdir_tar)
end
Dir.chdir(indir) {print `tar -czf #{File.basename(outdir_tar)} #{File.basename(outdir)}`}

print `rm -rf #{outdir}`
print `rm -rf #{r1_f}`
print `rm -rf #{r2_f}`