include("fasta.jl")
include("node.jl")
include("treelike.jl")

function run_like(ARGS)
  if length(ARGS) != 2
    println("arguments: treefile alnfile")
    exit()
  end
  #read the tree
  a = open(ARGS[1],"r")
  s = ""
  for i in eachline(a)
    s = i
    break
  end
  tree = read_tree_string(s)

  #read the seqs
  seqs = read_fasta_file(ARGS[2])

  match_tips_and_seqs(tree,seqs)
  #basefreq = [0.1,0.4,0.2,0.3]
  #rmatrix = [0 .3 .4 .3 ; .3 0 .3 .4 ; .4 .3 0 .3 ; 0.3 .4 .3 0]
  basefreq = [0.25,0.25,0.25,0.25]
  rmatrix = [0. 1. 1. 1. ; 1. 0. 1. 1. ; 1. 1. 0. 1. ; 1. 1. 1. 0.]
  @time x = calc_nuc_tree_likelihood(tree,rmatrix,basefreq,length(seqs[1].seq),true)
  println("loglike: $(x)")
  #Profile.print(format=:flat)
end

run_like(ARGS)
