type Seq
	name::String
	seq::String
end

function read_fasta_file(filename::String)
	a = open(filename,"r")
	tname = ""
	tseq = ""
  seqs = Seq[]
	first = true
	for i in eachline(a)
    i = strip(i)
    if length(i) == 0
      continue
    end
		if i[1]=='>'
			if first == false
        ts = Seq(tname,tseq)
				push!(seqs,ts)
			else
				first = false
			end
			tname = i[2:length(i)]
			tseq = ""
			continue
		else
			tseq = uppercase(string(tseq,i[1:length(i)]))
		end
	end
  ts = Seq(tname,tseq)
	push!(seqs,ts)
  return seqs
end

#read_fasta_file("test.fasta")
