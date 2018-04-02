#seperate consensus sequences for regions;

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

indir = ARGV[0].chomp
if indir[-1] == "/"
  indir = indir[0..-2]
end

rename = true if ARGV[1] == "-rename"

libs = []

Dir.chdir(indir) {libs = Dir.glob("*")}
outdir = indir + "_consensus"
Dir.mkdir(outdir) unless File.directory?(outdir)
outdir1 = outdir + "/consensus1"
outdir2 = outdir + "/consensus2"
outdir3 = outdir + "/consensus3"
Dir.mkdir(outdir1) unless File.directory?(outdir1)
Dir.mkdir(outdir2) unless File.directory?(outdir2)
Dir.mkdir(outdir3) unless File.directory?(outdir3)

combined_exist = false

libs.each do |lib|
  dir = []
  Dir.chdir((indir + "/" + lib)) {dir = Dir.glob("*")}
  if dir[0] == lib
    dir2 = indir + "/" + lib + "/" + lib
  else
    dir2 = indir + "/" + lib
  end
  bars = []
  Dir.chdir(dir2) {bars = Dir.glob("*")}
  bars.each do |bar|
    r1_file = ""
    r2_file = ""
    combined_file = ""
    dir3 = dir2 + "/" + bar + "/consensus"
    
    files = []
    Dir.chdir(dir3) {files = Dir.glob("*")}
    files.each do |f|
      if (f =~ /r2/i)
        r2_file = f
      elsif (f =~ /r1/i)
        r1_file = f
      elsif f =~ /combined/
        combined_file = f
        combined_exist = true
      end
    end
    r1_file = dir3 + "/" + r1_file
    r2_file = dir3 + "/" + r2_file
    
    r1 = fasta_to_hash(r1_file)
    r2 = fasta_to_hash(r2_file)
    r1_new = {}
    r2_new = {}
    if rename
      r1.each do |k,v|
        k = k + "_" + lib + "_" + bar + "_r1"
        r1_new[k] = v
      end
      r2.each do |k,v|
        k = k + "_" + lib + "_" + bar + "_r2"
        r2_new[k] = v
      end
    else
      r1_new = r1
      r2_new = r2
    end
    
    c1_d1 = outdir1 + "/" + bar
    Dir.mkdir(c1_d1) unless File.directory?(c1_d1)
    c1_d2_r1 = c1_d1 + "/R1"
    c1_d2_r2 = c1_d1 + "/R2"
    Dir.mkdir(c1_d2_r1) unless File.directory?(c1_d2_r1)
    Dir.mkdir(c1_d2_r2) unless File.directory?(c1_d2_r2)
    p1 = c1_d2_r1 + "/" + lib + "_" + bar + "_r1"
    p2 = c1_d2_r2 + "/" + lib + "_" + bar + "_r2"
    out1_r1 = File.open(p1,"w")
    out1_r2 = File.open(p2,"w")
    r1_new.each do |k,v|
      out1_r1.puts k
      out1_r1.puts v
    end
    r2_new.each do |k,v|
      out1_r2.puts k
      out1_r2.puts v
    end
    out1_r1.close
    out1_r2.close
    
    c2_d1 = outdir2 + "/" + lib
    Dir.mkdir(c2_d1) unless File.directory?(c2_d1)
    c2_d2 = c2_d1 + "/" + bar
    Dir.mkdir(c2_d2) unless File.directory?(c2_d2)
    p1 = c2_d2 + "/" + lib + "_" + bar + "_r1"
    p2 = c2_d2 + "/" + lib + "_" + bar + "_r2"
    out1_r1 = File.open(p1,"w")
    out1_r2 = File.open(p2,"w")
    r1_new.each do |k,v|
      out1_r1.puts k
      out1_r1.puts v
    end
    r2_new.each do |k,v|
      out1_r2.puts k
      out1_r2.puts v
    end
    out1_r1.close
    out1_r2.close
    
    c3_d1 = outdir3 + "/" + bar
    Dir.mkdir(c3_d1) unless File.directory?(c3_d1)
    c3_d2 = c3_d1 + "/" + lib
    Dir.mkdir(c3_d2) unless File.directory?(c3_d2)
    p1 = c3_d2 + "/" + lib + "_" + bar + "_r1"
    p2 = c3_d2 + "/" + lib + "_" + bar + "_r2"
    out1_r1 = File.open(p1,"w")
    out1_r2 = File.open(p2,"w")
    r1_new.each do |k,v|
      out1_r1.puts k
      out1_r1.puts v
    end
    r2_new.each do |k,v|
      out1_r2.puts k
      out1_r2.puts v
    end
    out1_r1.close
    out1_r2.close
    
    next if combined_file == ""
    combined_exist = true
  
    outdir4 = outdir + "/combined_consensus"
    Dir.mkdir(outdir4) unless File.directory?(outdir4)
    outdir_bar = outdir4 + "/" + bar
    Dir.mkdir(outdir_bar) unless File.directory?(outdir_bar)
    
    combined_file = dir3 + "/" + combined_file
    path_combined = outdir_bar + "/" + lib + "_" + bar + "_combined"
    print `cp #{combined_file} #{path_combined}` 
  end
end

indir = ARGV[0]
pools = []

outdir = indir + "_log"

Dir.mkdir(outdir) unless File.directory?(outdir)

log_file2 = outdir + "/log.csv"

log_out2 = File.open(log_file2,"w")

if combined_exist
  log_out2.puts "Poolname,Barcode,Raw Sequences per barcode,R1 Raw,R2 Raw,Paired Raw,Cutoff,PID Length,Consensus1,Consensus2,Distinct to Raw,Resampling index,R1_loc,R2_loc,Combined"
else
  log_out2.puts "Poolname,Barcode,Raw Sequences per barcode,R1 Raw,R2 Raw,Paired Raw,Cutoff,PID Length,Consensus1,Consensus2,Distinct to Raw,Resampling index"
end

a1 = []

Dir.chdir(indir){pools = Dir.glob("*")}

pools.each do |poolname|
  pool_dir = indir + "/" + poolname
  files1 = []
  Dir.chdir(pool_dir){files1 = Dir.glob("*")}
  next if files1.size == 0
  if poolname == files1[0]
    consensus_dir = pool_dir + "/" + poolname
  else
    consensus_dir = pool_dir
  end
 
  barcodes = []
  Dir.chdir(consensus_dir){barcodes = Dir.glob("*")}
  barcodes.each do |b1|
    a2 = []
    a2 << poolname
    a2 << b1
    log_f = consensus_dir + "/" + b1 + "/log.txt"
    lines = File.readlines(log_f).delete_if{|k|k == "\n"}
    (7..14).each do |line|
      lines[line].chomp[-15..-1] =~ /\d+/
      a2 << $&
    end
    (15..16).each do |line|
      a2 << lines[line].chomp.split("\s")[-1]
    end
    if combined_exist
      (17..19).each do |line|
         a2 << lines[line].chomp.split("\s")[-1]
      end
    end
    
    a1 << a2
  end
end

a1.each do |v|
  v.each do |n|
    log_out2.print n.to_s + ","
  end
  log_out2.print "\n"
end
log_out2.close
