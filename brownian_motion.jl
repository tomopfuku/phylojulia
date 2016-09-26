include("node.jl")

function read_tree(infl)
    treefl = open(infl)
    nwk = readall(treefl)
    return read_tree_string(nwk)
end
    
function read_traits(infl)
    traits = readdlm(infl,'\t')
    traits = traits[1:2,2:10]
    dic = Dict()
    for i in [1:9]
        dic[traits[1,i]]=traits[2,i]
    end
    return traits 
end

function brlen(node)
    println(node.label)
    #println(node.length)
    #println(node.children)
end

#get shared history between i and j (distance of mrca from root)
function mrca_height(tree,tax1::Node,tax2::Node)
  t1_path = Set()
  curnode = tax1
  while curnode.parent != None
    push!(t1_path,curnode.parent)
    curnode = curnode.parent
  end
  curnode = tax2
  while in(curnode,t1_path) == false
    curnode = curnode.parent
  end
  return curnode.height
end

#get distance between i and j
function pairwise_dist(tree,tax1::Node,tax2::Node)
  t1_path = Set()
  curnode = tax1
  while curnode.parent != None
    push!(t1_path,curnode.parent)
    curnode = curnode.parent
  end
  curnode = tax2
  while in(curnode,t1_path) == false
    curnode = curnode.parent
  end
  return tax1.height-curnode.height
end

function ou_vcv(alpha,sigma,tree = "/home/carolioness/documents/yeast_expression/expression_evol_yeast/data/RAxML_bipartitions.fa100.rr.treepl")
  tree = read_tree(tree)
  node_height(tree)
  mat = Array(Float64,9,9)
  icount = 0
  jcount = 0
  for i in return_preorder_array(tree)
    if i.istip == true
      icount += 1
    end
    for j in return_preorder_array(tree)
      if i.istip == true & j.istip == true
        if jcount < 9
          jcount += 1
        else
          jcount = 1
        end
        dist = pairwise_dist(tree,i,j)
        mat[icount,jcount] = (sigma/(2*alpha))*exp(-alpha*dist) 
      end
    end
  end
  return mat
end

function bm_vcv(sigma,tree = "/home/carolioness/documents/yeast_expression/expression_evol_yeast/data/RAxML_bipartitions.fa100.rr.treepl")
  tree = read_tree(tree)
  node_height(tree)
  mat = Array(Float64,9,9)
  icount = 0
  jcount = 0
  for i in return_preorder_array(tree) 
    if i.istip == true
      icount += 1
    end
    for j in return_preorder_array(tree)
      if i.istip == true & j.istip == true
        if jcount < 9
          jcount += 1
        else
          jcount = 1
        end
        shared_brlen = mrca_height(tree,i,j)
        mat[icount,jcount] = shared_brlen
      end
    end
  end
  return mat*sigma 
end

function fit_bm()
  tips = read_traits("/home/carolioness/documents/yeast_expression/expression_evol_yeast/data/expr.tsv")
  tree = read_tree("/home/carolioness/documents/yeast_expression/expression_evol_yeast/data/RAxML_bipartitions.fa100.rr.treepl")
  node_height(tree)
  taxa = tips[1,1:end]
end










