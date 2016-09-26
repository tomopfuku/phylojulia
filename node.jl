type Node
  label::String
  length::Float64
  parent
  children::Array{}
  data::Dict
  probs::Array{Float64,1}
  height::Float64
  istip::Bool
end

function return_preorder_array(n::Node)
  out = Node[]
  st = Node[]
  push!(st,n)
  while length(st) > 0
    c = pop!(st)
    push!(out,c)
    for i in c.children
      push!(st,i)
    end
  end
  return out
end

function return_postorder_array(n::Node)
  out = Node[]
  st = Node[]
  push!(st,n)
  while length(st) > 0
    c = pop!(st)
    push!(out,c)
    for i in c.children
      push!(st,i)
    end
  end
  return reverse(out)
end

function postorder (n::Node,f)
  for i in n.children
    postorder(i,f)
  end
  return f(n)
end

function node_height (n::Node)
  if n.parent == None
    n.height = 0.0
  end
  for i in n.children
    i.height = i.parent.height + i.length
    node_height(i)
  end
end

function add_child(p::Node,c::Node)
    push!(p.children,c)
    c.parent = p
end

function read_tree_string(instr::String)
    index = 1
    nextchar = instr[index]
    start = true
    keepgoing = true
    curnode = None
    root = None
    while keepgoing == true
        if nextchar == '\('
            if start == true
                root = Node("",0.0,None,Node[],Dict{}(),[],0.0,false)
                curnode = root
                start = false
            else
                newnode = Node("",0.0,None,Node[],Dict{}(),[],0.0,false)
                add_child(curnode,newnode)
                curnode = newnode
            end
        elseif nextchar == ','
            curnode = curnode.parent
        elseif nextchar == '\)'
            curnode = curnode.parent
            index += 1
            nextchar = instr[index]
            name = ""
            while true
                if nextchar == ',' || nextchar == '\)' || nextchar == ':' || nextchar == ';' || nextchar == '['
                    break
                end
                name = "$name$nextchar"
                index += 1
                nextchar = instr[index]
            end
            curnode.label = name
            index -= 1
        elseif nextchar == ';'
            keepgoing = false
            break
        elseif nextchar == ':'
            index += 1
            nextchar = instr[index]
            brlen = ""
            while true
                if nextchar == ',' || nextchar == '\)' || nextchar == ':' || nextchar == ';' || nextchar == '['
                    break
                end
                brlen = "$brlen$nextchar"
                index += 1
                nextchar = instr[index]
            end
            curnode.length = float(brlen)
            index -= 1
        elseif nextchar == ' '
            index += 1
            nextchar = instr[index]
        else
            newnode = Node("",0.0,None,Node[],Dict{}(),[],0.0,true) #TODO fix this
            add_child(curnode,newnode)
            curnode = newnode
            #istip
            name = ""
            while true
                if nextchar == ',' || nextchar == '\)' || nextchar == ':' || nextchar == ';' || nextchar == '['
                    break
                end
                name = "$name$nextchar"
                index += 1
                nextchar = instr[index]
            end
            curnode.label = name
            index -= 1
        end
        if index < length(instr)
            index += 1
        else
            keepgoing = false
        end
        nextchar = instr[index]
    end
    return root
end

function get_newick_repr(x::Node,bl=false)
    ret = ""
    for i in 1:length(x.children)
        if i == 1
            ret = string(ret,"(")
        end
        ret = string(ret,get_newick_repr(x.children[i],bl))
        if i == length(x.children)
            ret = string(ret,")")
        else
            ret = string(ret,",")
        end
    end
    if length(x.label) > 0
        ret = string(ret,x.label)
    end
    if bl==true
        ret = string(ret,":",string(x.length))
    end
    return ret
end

#s = "(a:3,(b:0.1,c:1.3)int_|_and_33.5:5)root;"
#tree = read_tree_string(s)
#println(get_newick_repr(tree,true))
#y(x) = if length(x.children) > 0
#  print(x.label," ")
#end
#postorder(tree,x -> y(x))
#for i in return_postorder_array(tree)
#  println(i.label)
#end
