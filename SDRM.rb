load "sequence.rb"
#SDRM analysis for PR, RT and IN regions.
#input argv is the directory combined_TCS_per_lib dir from log_multi.rb script. 

indir = ARGV[0]
libs = Dir[indir + "/*"]
outdir = indir + "_SDRM"
Dir.mkdir(outdir) unless File.directory?(outdir)

libs.each do |lib|
  next unless File.directory?(lib)
  lib_name = File.basename(lib)
  out_lib_dir = outdir + "/" + lib_name
  Dir.mkdir(out_lib_dir) unless File.directory?(out_lib_dir)
  sub_seq_files = Dir[lib + "/*"]
  seq_summary_file = out_lib_dir + "/" + lib_name + "_summary.csv"
  seq_summary_out = File.open(seq_summary_file, "w")
  seq_summary_out.puts "Region,TCS,TCS with A3G/F hypermutation,TCS with stop codon,TCS w/o hypermutation and stop codon, Poisson cutoff for minority mutation (>=)"
  
  point_mutation_file = out_lib_dir + "/" + lib_name + "_substitution.csv"
  point_mutation_out = File.open(point_mutation_file, "w")
  point_mutation_out.puts "region,TCS,AA position,wild type,mutation,number,percentage,95% CI low, 95% CI high, notes"
  
  linkage_file = out_lib_dir + "/" + lib_name + "_linkage.csv"
  linkage_out = File.open(linkage_file, "w")
  linkage_out.puts "region,TCS,mutation linkage,number,percentage,95% CI low, 95% CI high, notes"
  
  aa_report_file = out_lib_dir + "/" + lib_name + "_aa.csv"
  aa_report_out = File.open(aa_report_file, "w")
  aa_report_out.puts "region,ref.aa.positions,TCS.number," + $amino_acid_list.join(",")
  
  filtered_seq_dir = out_lib_dir + "/" + lib_name + "_filtered_seq"
  Dir.mkdir(filtered_seq_dir) unless File.directory?(filtered_seq_dir)
  
  point_mutation_list = []
  linkage_list = []
  aa_report_list = []
  
  sub_seq_files.each do |sub_seq|
    seq_basename = File.basename(sub_seq)
    seqs = fasta_to_hash(sub_seq)
    if seq_basename =~ /V1V3/i
      seq_summary_out.puts "V1V3,#{seqs.size.to_s},NA,NA,NA,NA,"
      print `cp #{sub_seq} #{filtered_seq_dir}`
    elsif seq_basename =~ /PR/i
      hypermut_seq = a3g_hypermut_seq_hash(seqs)[0]
      stop_codon_seq = stop_codon_seq_hash(seqs, 0)
      filtered_seq = seqs.difference(hypermut_seq).difference(stop_codon_seq)
      p_cutoff = poisson_minority_cutoff(filtered_seq.values, 0.0001, 20)
      seq_summary_out.puts "PR,#{seqs.size.to_s},#{hypermut_seq.size.to_s},#{stop_codon_seq.size.to_s},#{filtered_seq.size.to_s},#{p_cutoff.to_s}"
      filtered_out = File.open((filtered_seq_dir + "/" + seq_basename), "w")
      filtered_seq.each {|k,v| filtered_out.puts k; filtered_out.puts v}
      sdrm = sdrm_pr_bulk(filtered_seq, p_cutoff)
      point_mutation_list += sdrm[0]
      linkage_list += sdrm[1]
      aa_report_list += sdrm[2]

    elsif seq_basename =~/IN/i
      hypermut_seq = a3g_hypermut_seq_hash(seqs)[0]
      stop_codon_seq = stop_codon_seq_hash(seqs, 2)
      filtered_seq = seqs.difference(hypermut_seq).difference(stop_codon_seq)
      p_cutoff = poisson_minority_cutoff(filtered_seq.values, 0.0001, 20)
      seq_summary_out.puts "IN,#{seqs.size.to_s},#{hypermut_seq.size.to_s},#{stop_codon_seq.size.to_s},#{filtered_seq.size.to_s},#{p_cutoff.to_s}"
      filtered_out = File.open((filtered_seq_dir + "/" + seq_basename), "w")
      filtered_seq.each {|k,v| filtered_out.puts k; filtered_out.puts v}
      sdrm = sdrm_in_bulk(filtered_seq, p_cutoff)
      point_mutation_list += sdrm[0]
      linkage_list += sdrm[1]
      aa_report_list += sdrm[2]
      
    elsif seq_basename =~/RT/i
      rt_seq1 = {}
      rt_seq2 = {}
      seqs.each do |k,v|
        rt_seq1[k] = v[0,267]
        rt_seq2[k] = v[267..-1]
      end
      hypermut_seq_rt1 = a3g_hypermut_seq_hash(rt_seq1)[0]
      hypermut_seq_rt2 = a3g_hypermut_seq_hash(rt_seq2)[0]
      stop_codon_seq_rt1 = stop_codon_seq_hash(rt_seq1, 1)
      stop_codon_seq_rt2 = stop_codon_seq_hash(rt_seq2, 2)
      hypermut_seq_keys = (hypermut_seq_rt1.keys | hypermut_seq_rt2.keys)
      stop_codon_seq_keys = (stop_codon_seq_rt1.keys | stop_codon_seq_rt2.keys)
      reject_keys = (hypermut_seq_keys | stop_codon_seq_keys)
      filtered_seq = seqs.reject {|k,v| reject_keys.include?(k) }
      p_cutoff = poisson_minority_cutoff(filtered_seq.values, 0.0001, 20)
      seq_summary_out.puts "RT,#{seqs.size.to_s},#{hypermut_seq_keys.size.to_s},#{stop_codon_seq_keys.size.to_s},#{filtered_seq.size.to_s},#{p_cutoff.to_s}"
      filtered_out = File.open((filtered_seq_dir + "/" + seq_basename), "w")
      filtered_seq.each {|k,v| filtered_out.puts k; filtered_out.puts v}
      sdrm = sdrm_rt_bulk(filtered_seq, p_cutoff)
      point_mutation_list += sdrm[0]
      linkage_list += sdrm[1]
      aa_report_list += sdrm[2]
      
    end
  end
 
  point_mutation_list.each do |record|
    point_mutation_out.puts record.join(",")
  end
  linkage_list.each do |record|
    linkage_out.puts record.join(",")
  end
  aa_report_list.each do |record|
    aa_report_out.puts record.join(",")
  end
end
