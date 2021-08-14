function hi()
  print("hi")
end


# The MIT License (MIT)
#
# Copyright (c) January, 2014 michael otte
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.


# example of a node that can be used in the heap, where T is the type of data
# other nodes can also be used as long as they have the fields: heapIndex::Int
# inHeap::Bool and these fields are initialized the same way. Further, some
# way of calculating the key from the fields in the nodes must be provided
# (see the following functions)
mutable struct HeapNode{T}
  data::T
  heapIndex::Int
  inHeap::Bool

  # constructor functions
  function HeapNode{T}(D::T) where {T}
    new{T}(D,-1, false)
  end
end

############################# some helper functions ######################
# note that these can be changed via the constructors to allow different #
# ways of doing things, e.g., in case the same node is supposed to be in #
# multiple heaps at the same time.                                       #
##########################################################################

# returns the key value of the node
keyDefault(node::T) where {T} = node.data

# default less than function
lessThanDefault(a::T, b::T) where {T} = (a.data < b.data)

# default greater than function
greaterThanDefault(a::T, b::T) where {T} = (a.data > b.data)

# default heap marker function (marks when a node is in the heap)
markDefault(node::T) where {T} = (node.inHeap = true)

# default heap unmarker function (un marks when a node is removed)
unmarkDefault(node::T) where {T} = (node.inHeap = false)

# default heap check marker function (checks if the node is marked)
markedDefault(node::T) where {T} = node.inHeap

# sets the heap index to the value
setIndexDefault(node::T, val::Int) where {T} = (node.heapIndex = val)

# set the heap index to the unused value
unsetIndexDefault(node::T) where {T} = (node.heapIndex = -1)

# returns the heap index
getIndexDefault(node::T) where {T} = node.heapIndex


# a binary heap data structure that stores nodes of type T
mutable struct BinaryHeap{T, TK}
  heapNode::Array{T}         # stores the things that are in the heap
  indexOfLast::Int           # the index of the last node in the heap array
  parentOfLast::Int          # stores the index of the parent of the last node

  # functions used for interacting with marks and indicies
  key::Function              # returns the key value of a node
  lessThan::Function         # the less than function for type key type TK
  greaterThan::Function      # the greater than function for type key type TK
  mark::Function             # marks a node as being in the heap
  unmark::Function           # unmarks a node (as not being in the heap)
  marked::Function           # returns true if the node is in the heap
  setIndex::Function         # function used to set the index
  unsetIndex::Function       # function used to unset the index
  getIndex::Function         # function that returns the index

  # constructor functions
  function BinaryHeap{T, TK}() where {T, TK}
    maxSize::Int = 64;
    H = Array{T}(undef, maxSize, 1)
    il = 0;
    ipl = -1;
    new{T, TK}(H,il,ipl, keyDefault, lessThanDefault, greaterThanDefault,
    markDefault, unmarkDefault, markedDefault, setIndexDefault,
    unsetIndexDefault, getIndexDefault)
  end

  function BinaryHeap{T, TK}(maxSize::Int) where {T, TK}
    H = Array{T}(undef, maxSize, 1)
    il = 0;
    ipl = -1;
    new{T, TK}(H,il,ipl, keyDefault, lessThanDefault, greaterThanDefault,
    markDefault, unmarkDefault, markedDefault, setIndexDefault,
    unsetIndexDefault, getIndexDefault)
  end

  function BinaryHeap{T, TK}(key::Function, lessThan::Function,
    greaterThan::Function, mark::Function, unmark::Function,
    marked::Function, setIndex::Function, unsetIndex::Function,
    getIndex::Function) where {T, TK}
    maxSize::Int = 64;
    H = Array{T}(undef, maxSize, 1)
    il = 0;
    ipl = -1;
    new(H,il,ipl, key, lessThan, greaterThan, mark, unmark, marked, setIndex, unsetIndex, getIndex)
  end

  function BinaryHeap{T, TK}(maxSize::Int, key::Function, lessThan::Function,
    greaterThan::Function, mark::Function, unmark::Function, marked::Function,
    setIndex::Function, unsetIndex::Function, getIndex::Function) where {T, TK}
    H = Array{T}(undef, maxSize, 1)
    il = 0;
    ipl = -1;
    new(H,il,ipl, key, lessThan, greaterThan, mark, unmark, marked, setIndex, unsetIndex, getIndex)
  end
end

# compares a node n with its parent, and switches them if the parent's
# cost is more than the node's cost. Repeats if a switch happens.
function bubbleUp(H::THeap, n::Int) where {THeap}
  if(n == 1)
    return
  end

  parent::Int = div(n,2);
  while n != 1 && H.greaterThan(H.heapNode[parent,1], H.heapNode[n,1])
    # swap graph node pointers
    tempNode = H.heapNode[parent,1];
    H.heapNode[parent,1] = H.heapNode[n,1];
    H.heapNode[n,1] = tempNode;

    # update graph node heap index values
    H.setIndex(H.heapNode[parent,1], parent);
    H.setIndex(H.heapNode[n,1], n);

    # get new node and parent indicies
    n = parent;
    parent = div(n,2);
  end
end


# compares a node n with its children, and switches them if a child's cost
# is less than the node's cost. Repeats if a switch happens.
function bubbleDown(H::THeap, n::Int) where {THeap}
  child::Int = 0

  if 2*n == H.indexOfLast
    child = 2*n
  elseif 2*n+1 > H.indexOfLast
    return
  elseif H.lessThan(H.heapNode[2*n,1], H.heapNode[2*n+1,1])
    child = 2*n
  else
    child = 2*n+1
  end

  while n <= H.parentOfLast && H.lessThan(H.heapNode[child,1], H.heapNode[n,1])
    # swap graph node pointers
    tempNode = H.heapNode[child,1]
    H.heapNode[child,1] = H.heapNode[n,1]
    H.heapNode[n,1] = tempNode

    # update graph node heap index values
    H.setIndex(H.heapNode[child,1], child)
    H.setIndex(H.heapNode[n,1], n)

    # get new node and child indicies
    n = child;

    if 2*n == H.indexOfLast
      child = 2*n
    elseif 2*n+1 > H.indexOfLast
      return
    elseif H.lessThan(H.heapNode[2*n], H.heapNode[2*n+1])
       child = 2*n
    else
       child = 2*n+1
    end
  end
end

# add thisNode to the heap
function addToHeap(H::THeap, thisNode::T) where {THeap, T}

  if H.indexOfLast == size(H.heapNode,1)
    # out of space, so double size
    maxSize = size(H.heapNode,1)
    heapNodeTemp = copy(H.heapNode)
    H.heapNode = Array{T}(undef, maxSize*2, 1)
    H.heapNode[1:maxSize,1] = heapNodeTemp
  end

  if !H.marked(thisNode)
    H.indexOfLast += 1
    H.parentOfLast = div(H.indexOfLast,2)
    H.heapNode[H.indexOfLast,1] = thisNode
    H.setIndex(thisNode, H.indexOfLast)
    bubbleUp(H, H.indexOfLast)
    H.mark(thisNode)
  else
    println("problems")
a = b
  end
end

# returns the node that is on the top of the heap
function topHeap(H::THeap) where {THeap}
  if H.indexOfLast < 1
    return false
  end
  return H.heapNode[1,1]
end

# removes the top valued node from the heap and returns it
function popHeap(H::THeap) where {THeap}
  if H.indexOfLast < 1
    return false
  end
  oldTopNode = H.heapNode[1,1]
  H.heapNode[1,1] = H.heapNode[H.indexOfLast,1]
  H.setIndex(H.heapNode[1,1], 1)
  H.indexOfLast -= 1
  H.parentOfLast = div((H.indexOfLast),2)
  bubbleDown(H, 1)
  H.unmark(oldTopNode)
  H.unsetIndex(oldTopNode)
  return oldTopNode
end

# removes the node from the heap, assumes that it is in the heap
function removeFromHeap(H::THeap, thisNode::T) where {THeap, T}
  n::Int = H.getIndex(thisNode)

  movedNode = H.heapNode[H.indexOfLast,1]
  H.heapNode[n,1] = movedNode
  H.setIndex(movedNode, n)
  H.indexOfLast -= 1
  H.parentOfLast = div((H.indexOfLast),2)
  bubbleUp(H, n)
  bubbleDown(H, H.getIndex(movedNode))
  H.unmark(thisNode)
  H.unsetIndex(thisNode)
end

# updates a node that is already in the heap
function updateHeap(H::THeap, thisNode::T) where {THeap, T}
  if !H.marked(thisNode)
    error("trying to update a node that is not in the heap\n")
  end

  bubbleUp(H, H.getIndex(thisNode))
  bubbleDown(H, H.getIndex(thisNode))
end

# prints the heap values on the command line
function printHeap(H::THeap) where {THeap}
  i::Int = 1
  p::Int = 1

  while i <= H.indexOfLast
    print("$(H.key(H.heapNode[i,1])) ->  ");

    if 2*i <= H.indexOfLast
      print("$(H.key(H.heapNode[2*i,1]))");
    else
      print("NULL");
    end

    print("    ");

    if 2*i+1 <= H.indexOfLast
      println("$(H.key(H.heapNode[2*i+1,1]))");
    else
      println("NULL");
    end

    i+=1
  end
  print("\n\n")

end


# pops each item out and displays it
function printPopAllHeap(H::THeap) where {THeap}
  while H.indexOfLast >= 1
    node = popHeap(H)
    println("$(H.key(node))")
  end
end

# returns true if heap is good, false if bad, also prints a command line message
function checkHeap(H::THeap) where {THeap}
  i::Int = 2
  if H.indexOfLast < 1
    print("Heap is empty\n")
    return true
  elseif H.getIndex(H.heapNode[1,1]) != 1
    print("There is a problem with the heap (root)\n")
    return false
  end
  while i <= H.indexOfLast
    if(H.lessThan(H.heapNode[i,1], H.heapNode[div(i,2),1]))
      print("There is a problem with the heap order\n")
      return false;
    elseif(H.getIndex(H.heapNode[i,1]) != i)
      print("There is a problem with the heap node data  $(H.getIndex(H.heapNode[i,1])) != $(i) \n")
      return false;
    end
    i+=1
  end
  println("The heap is OK")
  return true;
end

# removes all items from the heap (quicker than popping everything)
# returns an array contining the heap items (unsorted)
function cleanHeap(H::THeap) where {THeap}
  retNodes = H.heapNode[1:H.indexOfLast,:]

  for i = 1:H.indexOfLast
    H.unmark(H.heapNode[i,1])
    H.unsetIndex(H.heapNode[i,1])
  end

  H.indexOfLast = 0
  H.parentOfLast = -1

  return retNodes
end

###############################################################################
# Use these functions instead if we want the heap to return the biggest thing #
###############################################################################

# compares a node n with its parent, and switches them if the parent's
# cost is less than the node's cost. Repeats if a switch happens.
function bubbleUpB(H::THeap, n::Int) where {THeap}
  if(n == 1)
    return
  end

  parent::Int = div(n,2);
  while n != 1 && H.lessThan(H.heapNode[parent,1], H.heapNode[n,1])
    # swap graph node pointers
    tempNode = H.heapNode[parent,1];
    H.heapNode[parent,1] = H.heapNode[n,1];
    H.heapNode[n,1] = tempNode;

    # update graph node heap index values
    H.setIndex(H.heapNode[parent,1], parent);
    H.setIndex(H.heapNode[n,1], n);

    # get new node and parent indicies
    n = parent;
    parent = div(n,2);
  end
end


# compares a node n with its children, and switches them if a child's cost
# is more than the node's cost. Repeats if a switch happens.
function bubbleDownB(H::THeap, n::Int) where {THeap}
   child::Int

  if 2*n == H.indexOfLast
    child = 2*n
  elseif 2*n+1 > H.indexOfLast
    return
  elseif H.greaterThan(H.heapNode[2*n,1], H.heapNode[2*n+1,1])
    child = 2*n
  else
    child = 2*n+1
  end

  while n <= H.parentOfLast && H.greaterThan(H.heapNode[child,1], H.heapNode[n,1])
    # swap graph node pointers
    tempNode = H.heapNode[child,1]
    H.heapNode[child,1] = H.heapNode[n,1]
    H.heapNode[n,1] = tempNode

    # update graph node heap index values
    H.setIndex(H.heapNode[child,1], child)
    H.setIndex(H.heapNode[n,1], n)

    # get new node and child indicies
    n = child;

    if 2*n == H.indexOfLast
      child = 2*n
    elseif 2*n+1 > H.indexOfLast
      return
    elseif H.greaterThan(H.heapNode[2*n], H.heapNode[2*n+1])
       child = 2*n
    else
       child = 2*n+1
    end
  end
end

# add thisNode to the heap
function addToHeapB(H::THeap, thisNode::T) where {THeap, T}

  if H.indexOfLast == size(H.heapNode,1)
    # out of space, so double size
    maxSize = size(H.heapNode,1)
    heapNodeTemp = copy(H.heapNode)
    H.heapNode = Array{T}(undef, maxSize*2, 1)
    H.heapNode[1:maxSize,1] = heapNodeTemp
  end

  if !H.marked(thisNode)
    H.indexOfLast += 1
    H.parentOfLast = div(H.indexOfLast,2)
    H.heapNode[H.indexOfLast,1] = thisNode
    H.setIndex(thisNode, H.indexOfLast)
    bubbleUpB(H, H.indexOfLast)
    H.mark(thisNode)
  end
end

# returns the node that is on the top of the heap
topHeapB(H::THeap) where {THeap} = topHeap(H::THeap)


# removes the top valued node from the heap and returns it
function popHeapB(H::THeap) where {THeap}
  if H.indexOfLast < 1
    return false
  end
  oldTopNode = H.heapNode[1,1]
  H.heapNode[1,1] = H.heapNode[H.indexOfLast,1]
  H.setIndex(H.heapNode[1,1], 1)
  H.indexOfLast -= 1
  H.parentOfLast = div((H.indexOfLast),2)
  bubbleDownB(H, 1)
  H.unmark(oldTopNode)
  H.unsetIndex(oldTopNode)

  return oldTopNode
end


# removes the node from the heap, assumes that it is in the heap
function removeFromHeapB(H::THeap, thisNode::T) where {THeap, T}
  n::Int = H.getIndex(thisNode)

  movedNode = H.heapNode[H.indexOfLast,1]
  H.heapNode[n,1] = movedNode
  H.setIndex(movedNode, n)
  H.indexOfLast -= 1
  H.parentOfLast = div((H.indexOfLast),2)
  bubbleUpB(H, n)
  bubbleDownB(H, H.getIndex(movedNode))
  H.unmark(thisNode)
  H.unsetIndex(thisNode)
end

# updates a node that is already in the heap
function updateHeapB(H::THeap, thisNode::T) where {THeap, T}
  if !H.marked(thisNode)
    error("trying to update a node that is not in the heap\n")
  end

  bubbleUpB(H, H.getIndex(thisNode))
  bubbleDownB(H, H.getIndex(thisNode))
end

# prints the heap values on the command line
printHeapB(H::THeap) where {THeap} = printHeap(H::THeap)

# pops each item out and displays it
printPopAllHeapB(H::THeap) where {THeap} = printPopAllHeap(H::THeap)


# returns 1 if heap is good, 0 if bad, also prints a command line message
function checkHeapB(H::THeap) where {THeap}
  i::Int = 2
  if H.indexOfLast < 1
    print("Heap is empty\n")
    return true
  elseif H.getIndex(H.heapNode[1,1]) != 1
    print("There is a problem with the heap (root)\n")
    return false
  end
  while i <= H.indexOfLast
    if(H.greaterThan(H.heapNode[i,1], H.heapNode[div(i,2),1]) || H.getIndex(H.heapNode[i,1]) != i)
      print("There is a problem with the heap \n")
      return false;
    end
    i+=1
  end
  println("The heap is OK")
  return true;
end

# removes all items from the heap (quicker than popping everything)
# returns an array contining the heap items (unsorted) and another
# array containing their key values
cleanHeapB(H::THeap) where {THeap} = cleanHeap(H::THeap)

#function testCase()
#  H = BinaryHeap{HeapNode{Float64}, Float64}(64)
#
#  for i = 1:100
#
#    if rand() > .2
#      #println("add")
#      node = HeapNode{Float64}(rand())
#      addToHeap(H, node)
#    elseif rand() > .5
#      #println("pop")
#      popHeap(H)
#    elseif H.indexOfLast > 1
#      #println("remove")
#      randN = randi(H.indexOfLast)
#      removeFromHeap(H, H.heapNode[randN])
#    end
#
#    if rand() > .5 && H.indexOfLast > 1
#      #println("update")
#      randN = randi(H.indexOfLast)
#      H.heapNode[randN].data = rand()
#      updateHeap(H, H.heapNode[randN])
#    end
#
#    checkHeap(H)
#  end
#  printHeap(H)
#  printPopAllHeap(H)
#end
#
#
#function testCaseB()
#  H = BinaryHeap{HeapNode{Float64}, Float64}(64)
#
#  for i = 1:100
#
#    if rand() > .2
#      #println("add")
#      node = HeapNode{Float64}(rand())
#      addToHeapB(H, node)
#    elseif rand() > .5
#      #println("pop")
#      popHeapB(H)
#    elseif H.indexOfLast > 1
#      #println("remove")
#      randN = randi(H.indexOfLast)
#      removeFromHeapB(H, H.heapNode[randN])
#    end
#
#    if rand() > .5 && H.indexOfLast > 1
#      #println("update")
#      randN = randi(H.indexOfLast)
#      H.heapNode[randN].data = rand()
#      updateHeapB(H, H.heapNode[randN])
#    end
#
#    checkHeapB(H)
#  end
#  printHeapB(H)
#  printPopAllHeapB(H)
#end
