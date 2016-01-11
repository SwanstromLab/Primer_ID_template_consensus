#input sequence file
infile = ARGV[0]

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

sequences = fasta_to_hash(infile)
puts sequences

outfile = File.dirname(infile) + "/" + File.basename(infile,".*") +  "_converted.txt"
out = File.open(outfile,"w")

sequences.each do |name,seq|
  out.puts name
  out.puts seq.scan(/../).collect {|c|c.to_i.chr}.join
end

out.close