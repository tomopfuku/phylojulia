include("node.jl")

function calc_nuc_q_matrix(rmatrix::Array{Float64,2},basefreq::Array{Float64,1})
  t = ones(4,4)
  for i in 1:4
    for j in 1:4
      if j == i
        t[i,j] = basefreq[i] * rmatrix[i,j]
      else
        t[i,j] = basefreq[i] * basefreq[j] * rmatrix[i,j]
      end
    end
  end
  tscale = sum(t)-sum(diag(t))
  t = t/tscale
  @inbounds @simd for i in 1:4
    t[i,i] = 0 - sum(t[i,:]) - t[i,i]
    t[i,:] = t[i,:] / basefreq[i]
  end
  return t
end

function match_tips_and_seqs(tree::Node,seqs::Array{Seq,1})
  y(x::Node,seqs::Array{Seq,1}) = if length(x.children) == 0
    test = false
    for j in seqs
      if j.name == x.label
        x.data["seq"] = j.seq
        test = true
        break
      end
    end
    if test == false
      println("can't find $(x.label) in seqs")
      exit()
    end
  end
  postorder(tree,x -> y(x,seqs))
end

function calc_p_matrix(q::Array{Float64,2},t::Float64)
  return expm(q*t)
end

function calc_nuc_tree_likelihood(tree::Node,rmatrix::Array{Float64,2},basefreq::Array{Float64,1},sites::Int,verbose::Bool)
  position = {'A'=>[1],'C'=>[2],'G'=>[3],'T'=>[4],'-'=>[1,2,3,4],'N'=>[1,2,3,4],'Y'=>[2,4],'R'=>[1,3],'W'=>[1,4],'M'=>[1,2],'B'=>[3,2,4],'V'=>[3,1,2],'S'=>[3,2],'K'=>[3,4],'H'=>[1,4,2],'D'=>[1,3,4]}
  q = calc_nuc_q_matrix(rmatrix,basefreq)
  #trying the deconstructed type !SLOWER!
  #evals,evecs = eig(q)
  #ievecs = inv(evecs)
  ps = Dict{Float64,Array{Float64,2}}()
  poa = return_postorder_array(tree)
  loglike = 0
  for s in 1:sites
    probs = []
    for x in poa #postorder traversal nodes in a list
       if length(x.children) == 0
        x.probs = zeros(4)
        positions = position[x.data["seq"][s]]
        for j in positions
          x.probs[j] = 1.
        end
      else
        x.probs = ones(4)
        for j in 1:4
          for i in x.children
            if in(i.length,keys(ps))
              p = ps[i.length]
            else
              p = calc_p_matrix(q,i.length)
              #using precalculated evecs and evals !SLOWER!
              #p = evecs*diagm(exp(evals*i.length))*ievecs
              ps[i.length] = p
            end
            templike = 0
            #templike = sum([(p[j,k] * i.probs[k]) for k in 1:4]) !SLOWER!
            @inbounds @simd for k in 1:4
              templike += (p[j,k] * i.probs[k])
            end
            x.probs[j] *= templike
          end
        end
        if x == tree
          probs = x.probs
        end
      end
    end
    for k in 1:4
      probs[k] *= basefreq[k]
    end
    loglike += (log(sum(probs)))
  end
  return loglike
end

#basefreq = [0.1,0.4,0.2,0.3]
#rmatrix = [0 .3 .4 .3 ; .3 0 .3 .4 ; .4 .3 0 .3 ; 0.3 .4 .3 0]
#calc_nuc_q_matrix(rmatrix,basefreq)
