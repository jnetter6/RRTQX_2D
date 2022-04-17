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



# saves the data to the file
#function saveData(data, fileName)
#  fptr = open(fileName, "w")
#  for i = 1:size(data,1)
#    writecsv(fptr, data[i,:])
#  end
#  close(fptr)
#end

include("ghostPoint.jl")
# require("jlist.jl")

# note that this implimentaion of the KD tree should not be used for cartesian spaces in practice (since I've written a faster implimentation that saves dist squared instad of distance for that particular case)

# also note that, by default, the space is assumed to be R^d
# however, it is also possible to allow the space to be some cartesian product of
# R^n and S^m subspaces, more codeing is required if other subspaces are desired,
# e.g. quarternions. Make sure to use an approperiate distance function.
# The particular topology of the space is defined in KDTree{T},
# see notes. The way that wrapping subspaces are handles, currently, is to run
# a check with the original point, and then to run an additional run 3^(w) - 1
# checks, where w is the number of wrapping dimensions, for each possibility of
# wrapping the querery point around some/all of the wrappable dimensions.

# in practice, many of these can be skipped if the distance from the wrapped point
# to the closest point in the unwrapped space > the distance to the nearest point
# found so far. FINALLY NOTE that the kd-tree distance function used to help
# find nearest neighbors should actually treat S dimesnions as R dimensions, since
# ghost points are used for wrapping!! so you'll need a swperate kd-tree distance
# function than the one used in the rest of your code where S dimension are
# wrapped inside the distance function.



# example of a node that can be used in the KDTree, where T is the type of data
# used to measure distance along each dimension
# other nodes can also be used as long as they have these fields that
# are initialized as follows by a default constructor and the parent and
# children types are the same as the node's type itself
mutable struct KDTreeNode{T}

  # data used for KD Tree
  kdInTree::Bool           # set to true if this node is in the kd-tree
  kdParentExist::Bool      # set to true if parent in the tree is used
  kdChildLExist::Bool      # set to true if left child in the tree is used
  kdChildRExist::Bool      # set to true if right child in the tree is used

  # data used for heap in KNN-search
  heapIndex::Int           # named such to allow the use of default heap functions
  inHeap::Bool             # ditto
  data::Float64            # ditto, this will hold the distance

  # more data used for KD Tree
  position::Array{T}       # a dX1 array where d is the dimesnions of the space
  kdSplit::Int             # the dimension used for splitting at this node
  kdParent::KDTreeNode{T}  # parent in the tree
  kdChildL::KDTreeNode{T}  # left child in the tree
  kdChildR::KDTreeNode{T}  # right child in the tree

  # constructors
  KDTreeNode{T}() where {T} = new{T}(false, false, false, false, -1, false)

end

# a KD-Tree data structure that stores nodes of type T
mutable struct KDTree{T}
  d::Int                      # the number of dimensions in the space
  distanceFunction::Function  # the distance function to use: f(poseA, poseB)
  treeSize::Int               # the number of nodes in the KD-Tree

  numWraps::Int               # the total number of dimensions that wrap
  wraps::Array{Int}           # a vector of length d containing a list of
                              # all the dimensions that wrapAround
  wrapPoints::Array{Float64}  # space is assumed to start at 0 and end at
                              # wrapPoints[i] along dimension wraps[i]

  root::T                     # the root node

  # constructors
  KDTree{T}(d::Int, f::Function, wraps::Array{Int}, wrapPoints::Array{Float64}) where {T} = new{T}(d,
  f,0,length(wraps), wraps, wrapPoints)
  KDTree{T}(d::Int, f::Function) where {T} = new{T}(d,f,0,0)
  KDTree{T}() where {T} = new{T}(0, none, 0, 0)
end

function KDTreeInit(K::KDTree, d::Int, f::Function)
  K.d = d
  K.f = f
  K.treeSize = 0
end

# inserts a new node into the tree
function kdInsert(tree::TKD, node::T) where {TKD, T}
  if node.kdInTree
    return
  end
  node.kdInTree = true

  if tree.treeSize == 0
    tree.root = node
    tree.root.kdSplit = 1
    tree.treeSize = 1
    return
  end

  # figure out where to put this node
  parent::T = tree.root
  while true
    if node.position[parent.kdSplit] < parent.position[parent.kdSplit]
      # traverse tree to the left
      if !parent.kdChildLExist
        # the node gets inserted as the left child of the parent
        parent.kdChildL = node
        parent.kdChildLExist = true
        break
      end

      parent = parent.kdChildL
      continue
    else
      # traverse tree to the right
      if !parent.kdChildRExist
        # the node gets inserted as the right child of the parent
        parent.kdChildR = node
        parent.kdChildRExist = true
        break
      end

      parent = parent.kdChildR
      continue
    end
  end

  node.kdParent = parent
  node.kdParentExist = true
  if parent.kdSplit == tree.d
    node.kdSplit = 1
  else
    node.kdSplit = parent.kdSplit + 1
  end
  tree.treeSize += 1
end

# prints the sub-tree on the command line using dfs
function kdPrintSubTree(node::T) where {T}
   print("[$(node.kdSplit)]:$(node.position[node.kdSplit]) ->  ")

   if node.kdChildLExist
     print("$(node.kdChildL.position[node.kdSplit])")
   else
     print("NULL")
   end

   print("   |   ")

   if node.kdChildRExist
     println("$(node.kdChildR.position[node.kdSplit])")
   else
     println("NULL")
   end

   if node.kdChildLExist
     kdPrintSubTree(node.kdChildL)
   end
   if node.kdChildRExist
     kdPrintSubTree(node.kdChildR)
   end

end

# prints the tree on the command line
function kdPrintTree(tree::TKD) where {TKD}
   if tree.treeSize == 0
     println("tree is empty")
     return
   end

   kdPrintSubTree(tree.root)
end



############################### Nearest ##################################


# naive implimentation used for error checking DO NOT USE
function kdFindNearestInSubtreeNaive(distanceFunction::Function, root::T,
  queryPoint::Array{Float64}) where {T}

  # calculate this dist
  bestDist = distanceFunction(queryPoint,root.position)
  bestNode = root

  # find left subtree dist
  if root.kdChildLExist
    (Lnode, Ldist) = kdFindNearestInSubtreeNaive(distanceFunction, root.kdChildL, queryPoint)
    if Ldist < bestDist
      bestDist = Ldist
      bestNode = Lnode
    end
  end

  # find right subtree dist
  if root.kdChildRExist
    (Rnode, Rdist) = kdFindNearestInSubtreeNaive(distanceFunction, root.kdChildR, queryPoint)
    if Rdist < bestDist
      bestDist = Rdist
      bestNode = Rnode
    end
  end

  return (bestNode, bestDist)
end


# naive implimentation used for error checking DO NOT USE
function kdFindNearestNaive(tree::TKD, queryPoint::Array{Float64}) where {TKD}
  return kdFindNearestInSubtreeNaive(tree.distanceFunction, tree.root, queryPoint)
end



# returns the nearest node to queryPoint in the subtree starting at root
# and also its distance, it takes also takes a suggestion for a
# possible closest node (and uses that if it is best)
function kdFindNearestInSubtree(distanceFunction::Function, root::T,
  queryPoint::Array{Float64}, suggestedClosestNode::T,
  suggestedClosestDist::Float64) where {T}

  # walk down the tree as if the node would be inserted
  parent::T = root
  currentClosestNode::T = suggestedClosestNode
  currentClosestDist::Float64 = suggestedClosestDist
  while true
    if queryPoint[parent.kdSplit] < parent.position[parent.kdSplit]
      # traverse tree to the left
      if !parent.kdChildLExist
        # the queryPoint would be inserted as the left child of the parent
        break
      end
      parent = parent.kdChildL
      continue
    else
      # traverse tree to the right
      if !parent.kdChildRExist
        # the queryPoint would be inserted as the right child of the parent
        break
      end
      parent = parent.kdChildR
      continue
    end
  end

  newDist::Float64 = distanceFunction(queryPoint,parent.position)
  if newDist < currentClosestDist
    currentClosestNode = parent
    currentClosestDist = newDist
  end

  # now walk back up the tree (will break out when done)
  while true
    # now check if there could possibly be any closer nodes on the other
    # side of the parent, if not then check grandparent etc.

    parentHyperPlaneDist = (queryPoint[parent.kdSplit] - parent.position[parent.kdSplit])

    if parentHyperPlaneDist > currentClosestDist
      # then there could not be any closer nodes on the other side of the parent
      # (and the parent itself is also too far away

      if parent == root
        # the parent is the root and we are done
        return (currentClosestNode, currentClosestDist)
      end

      parent = parent.kdParent
      continue
    end

    # if we are here, then there could be a closer node on the other side of the
    # parent (including the parent itself)

    # first check the parent itself (if it is not already the closest node)
    if currentClosestNode != parent
      newDist = distanceFunction(queryPoint,parent.position)
      if newDist < currentClosestDist
        currentClosestNode = parent
        currentClosestDist = newDist
      end
    end

    # now check on the other side of the parent
    if queryPoint[parent.kdSplit] < parent.position[parent.kdSplit] && parent.kdChildRExist
      # queryPoint is on the left side of the parent, so we need to look
      # at the right side of it (if it exists)

      # find right subtree dist

      (Rnode, Rdist) = kdFindNearestInSubtree(distanceFunction, parent.kdChildR, queryPoint, currentClosestNode, currentClosestDist)

      if Rdist < currentClosestDist
        currentClosestDist = Rdist
        currentClosestNode = Rnode
      end

    elseif parent.position[parent.kdSplit] <= queryPoint[parent.kdSplit] && parent.kdChildLExist
      # queryPoint is on the right side of the parent, so we need to look
      # at the left side of it (if it exists)

      # find left subtree dist

      (Lnode, Ldist) = kdFindNearestInSubtree(distanceFunction, parent.kdChildL, queryPoint, currentClosestNode, currentClosestDist)
      if Ldist < currentClosestDist
        currentClosestDist = Ldist
        currentClosestNode = Lnode
      end
    end

    if parent == root
        # the parent is the root and we are done
      return (currentClosestNode, currentClosestDist)
    end

    parent = parent.kdParent
  end
end

function kdFindNearestInSubtreeTerrain(distanceFunction::Function, root::T,
  queryPoint::Array{Float64}, suggestedClosestNode::T,
  suggestedClosestDist::Float64) where {T}

  # walk down the tree as if the node would be inserted
  parent::T = root
  currentClosestNode::T = suggestedClosestNode
  currentClosestDist::Float64 = suggestedClosestDist
  while true
    if queryPoint[parent.kdSplit] < parent.position[parent.kdSplit]
      # traverse tree to the left
      if !parent.kdChildLExist
        # the queryPoint would be inserted as the left child of the parent
        break
      end
      parent = parent.kdChildL
      continue
    else
      # traverse tree to the right
      if !parent.kdChildRExist
        # the queryPoint would be inserted as the right child of the parent
        break
      end
      parent = parent.kdChildR
      continue
    end
  end

  newDist::Float64 = distanceFunction(queryPoint,parent.position)
  midNode = [((queryPoint[1]+parent.position[1])/2), ((queryPoint[2]+parent.position[2])/2)]
  if ((Wdist(midNode, [0.0, 0.0]) < 5.0) || (Wdist(midNode, [-2.0, -2.0]) < 5.0) || (Wdist(midNode, [-4.0, -4.0]) < 5.0))
    newDist = newDist *8.0
  end
  if newDist < currentClosestDist
    currentClosestNode = parent
    currentClosestDist = newDist
  end

  # now walk back up the tree (will break out when done)
  while true
    # now check if there could possibly be any closer nodes on the other
    # side of the parent, if not then check grandparent etc.

    parentHyperPlaneDist = (queryPoint[parent.kdSplit] - parent.position[parent.kdSplit])

    if parentHyperPlaneDist > currentClosestDist
      # then there could not be any closer nodes on the other side of the parent
      # (and the parent itself is also too far away

      if parent == root
        # the parent is the root and we are done
        return (currentClosestNode, currentClosestDist)
      end

      parent = parent.kdParent
      continue
    end

    # if we are here, then there could be a closer node on the other side of the
    # parent (including the parent itself)

    # first check the parent itself (if it is not already the closest node)
    if currentClosestNode != parent
      newDist = distanceFunction(queryPoint,parent.position)
      midNode = [((queryPoint[1]+parent.position[1])/2), ((queryPoint[2]+parent.position[2])/2)]
      if ((Wdist(midNode, [0.0, 0.0]) < 5.0) || (Wdist(midNode, [-2.0, -2.0]) < 5.0) || (Wdist(midNode, [-4.0, -4.0]) < 5.0))
        newDist = newDist *8.0
      end
      if newDist < currentClosestDist
        currentClosestNode = parent
        currentClosestDist = newDist
      end
    end

    # now check on the other side of the parent
    if queryPoint[parent.kdSplit] < parent.position[parent.kdSplit] && parent.kdChildRExist
      # queryPoint is on the left side of the parent, so we need to look
      # at the right side of it (if it exists)

      # find right subtree dist

      (Rnode, Rdist) = kdFindNearestInSubtreeTerrain(distanceFunction, parent.kdChildR, queryPoint, currentClosestNode, currentClosestDist)

      if Rdist < currentClosestDist
        currentClosestDist = Rdist
        currentClosestNode = Rnode
      end

    elseif parent.position[parent.kdSplit] <= queryPoint[parent.kdSplit] && parent.kdChildLExist
      # queryPoint is on the right side of the parent, so we need to look
      # at the left side of it (if it exists)

      # find left subtree dist

      (Lnode, Ldist) = kdFindNearestInSubtreeTerrain(distanceFunction, parent.kdChildL, queryPoint, currentClosestNode, currentClosestDist)
      if Ldist < currentClosestDist
        currentClosestDist = Ldist
        currentClosestNode = Lnode
      end
    end

    if parent == root
        # the parent is the root and we are done
      return (currentClosestNode, currentClosestDist)
    end

    parent = parent.kdParent
  end
end

# returns the nearest node to queryPoint and also its distance
function kdFindNearest(tree::TKD, queryPoint::Array{Float64}) where {TKD}
  # initial search (only search if the space does not wrap around)
  distToRoot::Float64 = tree.distanceFunction(queryPoint, tree.root.position)
  (Lnode, Ldist) = kdFindNearestInSubtree(tree.distanceFunction, tree.root, queryPoint, tree.root, distToRoot)

  if tree.numWraps > 0
    # if dimensions wrap around, we need to search vs. identities (ghosts)

    pointIterator = ghostPointIterator{TKD}(tree, queryPoint)
    while true
      thisGhostPoint = getNextGhostPoint(pointIterator, Ldist)
      if thisGhostPoint == nothing
        break
      end

      # now see if any points in the space are closer to this ghost
      distGhostToRoot = tree.distanceFunction(thisGhostPoint, tree.root.position)
      (thisLnode, thisLdist) = kdFindNearestInSubtree(tree.distanceFunction, tree.root, thisGhostPoint, tree.root, distGhostToRoot)

      if thisLdist < Ldist
        # found closer point
        Ldist = thisLdist
        Lnode = thisLnode
      end
    end
  end

  return (Lnode, Ldist)
end

function kdFindNearestTerrain(tree::TKD, queryPoint::Array{Float64}) where {TKD}
  # initial search (only search if the space does not wrap around)
  distToRoot::Float64 = tree.distanceFunction(queryPoint, tree.root.position)
  (Lnode, Ldist) = kdFindNearestInSubtreeTerrain(tree.distanceFunction, tree.root, queryPoint, tree.root, distToRoot)

  if tree.numWraps > 0
    # if dimensions wrap around, we need to search vs. identities (ghosts)

    pointIterator = ghostPointIterator{TKD}(tree, queryPoint)
    while true
      thisGhostPoint = getNextGhostPoint(pointIterator, Ldist)
      if thisGhostPoint == nothing
        break
      end

      # now see if any points in the space are closer to this ghost
      distGhostToRoot = tree.distanceFunction(thisGhostPoint, tree.root.position)
      (thisLnode, thisLdist) = kdFindNearestInSubtreeTerrain(tree.distanceFunction, tree.root, thisGhostPoint, tree.root, distGhostToRoot)

      if thisLdist < Ldist
        # found closer point
        Ldist = thisLdist
        Lnode = thisLnode
      end
    end
  end

  return (Lnode, Ldist)
end



# returns the nearest node to queryPoint in the subtree starting at root
# and also its distance, it takes also takes a suggestion for a
# possible closest node (and uses that if it is best)
function kdFindNearestInSubtreeWithGuess(distanceFunction::Function,
  root::T, queryPoint::Array{Float64}, suggestedClosestNode::T,
  suggestedClosestDist::Float64) where {T}

  # walk down the tree as if the node would be inserted
  parent::T = suggestedClosestNode.kdParent
  currentClosestNode::T = suggestedClosestNode
  currentClosestDist::Float64 = suggestedClosestDist
  while true
    if queryPoint[parent.kdSplit] < parent.position[parent.kdSplit]
      # traverse tree to the left
      if !parent.kdChildLExist
        # the queryPoint would be inserted as the left child of the parent
        break
      end
      parent = parent.kdChildL
      continue
    else
      # traverse tree to the right
      if !parent.kdChildRExist
        # the queryPoint would be inserted as the right child of the parent
        break
      end
      parent = parent.kdChildR
      continue
    end
  end

  newDist::Float64 = distanceFunction(queryPoint,parent.position)
  if newDist < currentClosestDist
    currentClosestNode = parent
    currentClosestDist = newDist
  end

  # now walk back up the tree (will break out when done)
  while true
    # now check if there could possibly be any closer nodes on the other
    # side of the parent, if not then check grandparent etc.

    parentHyperPlaneDist = (queryPoint[parent.kdSplit] - parent.position[parent.kdSplit])

    if parentHyperPlaneDist > currentClosestDist
      # then there could not be any closer nodes on the other side of the parent
      # (and the parent itself is also too far away

      if parent == root
        # the parent is the root and we are done
        return (currentClosestNode, currentClosestDist)
      end

      parent = parent.kdParent
      continue
    end

    # if we are here, then there could be a closer node on the other side of the
    # parent (including the parent itself)

    # first check the parent itself (if it is not already the closest node)
    if currentClosestNode != parent
      newDist = distanceFunction(queryPoint,parent.position)
      if newDist < currentClosestDist
        currentClosestNode = parent
        currentClosestDist = newDist
      end
    end

    # now check on the other side of the parent
    if queryPoint[parent.kdSplit] < parent.position[parent.kdSplit] && parent.kdChildRExist
      # queryPoint is on the left side of the parent, so we need to look
      # at the right side of it (if it exists)

      # find right subtree dist

      (Rnode, Rdist) = kdFindNearestInSubtree(distanceFunction, parent.kdChildR, queryPoint, currentClosestNode, currentClosestDist)

      if Rdist < currentClosestDist
        currentClosestDist = Rdist
        currentClosestNode = Rnode
      end

    elseif parent.position[parent.kdSplit] <= queryPoint[parent.kdSplit] && parent.kdChildLExist
      # queryPoint is on the right side of the parent, so we need to look
      # at the left side of it (if it exists)

      # find left subtree dist

      (Lnode, Ldist) = kdFindNearestInSubtree(distanceFunction, parent.kdChildL, queryPoint, currentClosestNode, currentClosestDist)
      if Ldist < currentClosestDist
        currentClosestDist = Ldist
        currentClosestNode = Lnode
      end
    end

    if parent == root
        # the parent is the root and we are done
        # need to do one last check vs the root
        thisDist = distanceFunction(queryPoint,parent.position)
        if thisDist < currentClosestDist
          return (parent, thisDist)
        end

      return (currentClosestNode, currentClosestDist)
    end

    parent = parent.kdParent
  end
end


# returns the nearest node to queryPoint and also its distance squared,
# instead of starting at the root, it starts at guess
function kdFindNearestWithGuesstree(tree::TKD, queryPoint::Array{Float64}, guess::T) where {TKD, T}
  distToGuess::Float64 = tree.distanceFunction(queryPoint, guess.position)
  if guess == tree.root
    return  kdFindNearestInSubtree(tree.distanceFunction, tree.root, queryPoint, tree.root, distToGuess)
  end

  (Lnode, Ldist) = kdFindNearestInSubtreeWithGuess(tree.distanceFunction, tree.root, queryPoint, guess, distToGuess)

  if tree.numWraps > 0
    # if dimensions wrap around, we need to search vs. identities (ghosts)

    pointIterator = ghostPointIterator{TKD}(tree, queryPoint)
    while true
      thisGhostPoint = getNextGhostPoint(pointIterator, Ldist)
      if thisGhostPoint == nothing
        break
      end

      # now see if any points in the space are closer to this ghost
      distGhostToGuess = tree.distanceFunction(thisGhostPoint, guess.position)
      (thisLnode, thisLdist) = kdFindNearestInSubtreeWithGuess(tree.distanceFunction, tree.root, thisGhostPoint, guess, distGhostToGuess)

      if thisLdist < Ldist
        # found closer point
        Ldist = thisLdist
        Lnode = thisLnode
      end
    end
  end

  return (Lnode, Ldist)
end

############################### K Nearest #############################



# naive implimentation used for error checking DO NOT USE
function kdFindKNearestInSubtreeNaive(distanceFunction::Function, root::T,
  queryPoint::Array{Float64}, nearestHeap::BinaryHeap) where {T}

  # basically, just add all nodes in the tree to the heap
  if !root.inHeap
    root.data = distanceFunction(queryPoint, root.position)
    addToHeapB(nearestHeap, root)
  end

  # go down the left subtree
  if root.kdChildLExist
    kdFindKNearestInSubtreeNaive(distanceFunction, root.kdChildL, queryPoint, nearestHeap)
  end

  # go down the right subtree
  if root.kdChildRExist
    kdFindKNearestInSubtreeNaive(distanceFunction, root.kdChildR, queryPoint, nearestHeap)
  end
end


# naive implimentation used for error checking DO NOT USE
function kdFindKNearestNaive(tree::TKD, k::Int, queryPoint::Array{Float64}) where {TKD}
  H = BinaryHeap{typeof(tree.root), Float64}(64)

  # add all nodes in tree to a heap
  kdFindKNearestInSubtreeNaive(tree.distanceFunction, tree.root, queryPoint, H)

  # now remove far away nodes until we have the k that we were looking for
  while H.indexOfLast > k
    popHeapB(H)
  end
  return cleanHeapB(H)
end



# adds the node to the heap IF there is space in the current heap without growing
# past k, othewise the curent top is removed first, returns the current top
function addToKNNHeap(H::THeap, thisNode::T, key::Float64, k::Int) where {THeap, T}
  if thisNode.inHeap
    return topHeapB(H)
  elseif H.indexOfLast < k
    # just insert
    thisNode.data = key
    addToHeapB(H, thisNode)
  elseif H.heapNode[1,1].data > key
    popHeapB(H)
    thisNode.data = key
    addToHeapB(H, thisNode)
  end
  return topHeapB(H)
end


# finds the K nearest nodes to queryPoint in the subtree starting at root
# note that this data is stored in the nearestHeap, and the heap may also contain
# nodes before this function is called explicitly returns the node of the nearest
# set that is FURTHEREST from the querery along with its distance
# assumes that there is at least one node in the heap to begin with
# IF this is the first call to this function (e.g., from kdFindKNearest) then
# it is also assumed that a dummy node has been added that has an Inf key value
# this makes things easier with checking that all K slots are used during the
# recursion
function kdFindKNearestInSubtree(distanceFunction::Function, root::T, k::Int,
  queryPoint::Array{Float64}, nearestHeap::BinaryHeap) where {T}

  # walk down the tree as if the node would be inserted
  parent::T = root
  currentWorstClosestNode = topHeapB(nearestHeap)
  currentWorstClosestDist = currentWorstClosestNode.data
  while true
    if queryPoint[parent.kdSplit] < parent.position[parent.kdSplit]
      # traverse tree to the left
      if !parent.kdChildLExist
        # the queryPoint would be inserted as the left child of the parent
        break
      end
      parent = parent.kdChildL
      continue
    else
      # traverse tree to the right
      if !parent.kdChildRExist
        # the queryPoint would be inserted as the right child of the parent
        break
      end
      parent = parent.kdChildR
      continue
    end
  end

  newDist::Float64 = distanceFunction(queryPoint,parent.position)
  if newDist < currentWorstClosestDist
    currentWorstClosestNode = addToKNNHeap(nearestHeap, parent, newDist, k)
    currentWorstClosestDist = currentWorstClosestNode.data
  end

  # now walk back up the tree (will break out when done)
  while true

    # now check if there could possibly be any closer nodes on the other
    # side of the parent, if not then check grandparent etc.

    parentHyperPlaneDist = (queryPoint[parent.kdSplit] - parent.position[parent.kdSplit])

    if parentHyperPlaneDist >  currentWorstClosestDist
      # then there could not be any closer nodes on the other side of the parent
      # (and the parent itself is also too far away

      if parent == root
        # the parent is the root and we are done
        return currentWorstClosestNode
      end

      parent = parent.kdParent
      continue
    end

    # if we are here, then there could be a closer node on the other side of the
    # parent (including the parent itself)

    # first check the parent itself (if it is not already one of the closest nodes)
    if !parent.inHeap
      newDist = distanceFunction(queryPoint,parent.position)
      if newDist < currentWorstClosestDist
        currentWorstClosestNode = addToKNNHeap(nearestHeap, parent, newDist, k)
        currentWorstClosestDist = currentWorstClosestNode.data
      end
    end

    # now check on the other side of the parent
    if queryPoint[parent.kdSplit] < parent.position[parent.kdSplit] && parent.kdChildRExist
      # queryPoint is on the left side of the parent, so we need to look
      # at the right side of it (if it exists)
      currentWorstClosestNode = kdFindKNearestInSubtree(distanceFunction, parent.kdChildR, k, queryPoint, nearestHeap)
      currentWorstClosestDist = currentWorstClosestNode.data
    elseif parent.position[parent.kdSplit] <= queryPoint[parent.kdSplit] && parent.kdChildLExist
      # queryPoint is on the right side of the parent, so we need to look
      # at the left side of it (if it exists)
      currentWorstClosestNode = kdFindKNearestInSubtree(distanceFunction, parent.kdChildL, k, queryPoint, nearestHeap)
      currentWorstClosestDist = currentWorstClosestNode.data
    end

    if parent == root
      # the parent is the root and we are done

      return currentWorstClosestNode
    end

    parent = parent.kdParent
  end
end

# returns the K nearest nodes to queryPoint and also their distance
# (note that they are not sorted, but they are in (reverse) heap order)
function kdFindKNearest(tree::TKD, k::Int, queryPoint::Array{Float64}) where {TKD}
  H = BinaryHeap{typeof(tree.root), Float64}(k)

  # insert root node in heap
  tree.root.data = tree.distanceFunction(queryPoint, tree.root.position)
  addToHeapB(H, tree.root)

  # insert a dummy node in the heap with inf key
  dummyNode = (typeof(tree.root))()
  dummyNode.data = Inf
  addToHeapB(H, dummyNode)

  # find k nearest neighbors
  kdFindKNearestInSubtree(tree.distanceFunction, tree.root, k, queryPoint, H)

  if tree.numWraps > 0
    error("knn search has not been implimented for wrapped space")
  end


  # remove the dummy node if still there (guarenteed to be on top, due to Inf key)
  topNode = topHeapB(H)
  if topNode == dummyNode
    popHeapB(H)
  end

  return cleanHeapB(H)
end



############################# Within Range #############################



# naive implimentation used for error checking DO NOT USE
function kdFindWithinRangeInSubtreeNaive(distanceFunction::Function,
  root::T, queryPoint::Array{Float64}, nearestHeap::BinaryHeap, range::Float64) where {T}

  # basically, just add all nodes in the tree to the heap if they are within range
  if !root.inHeap && distanceFunction(queryPoint, root.position) <= range
    root.data = distanceFunction(queryPoint, root.position)
    addToHeapB(nearestHeap, root)
  end

  # go down the left subtree
  if root.kdChildLExist
    kdFindWithinRangeInSubtreeNaive(distanceFunction, root.kdChildL, queryPoint, nearestHeap, range)
  end

  # go down the right subtree
  if root.kdChildRExist
    kdFindWithinRangeInSubtreeNaive(distanceFunction, root.kdChildR, queryPoint, nearestHeap, range)
  end
end


# naive implimentation used for error checking DO NOT USE
function kdFindWithinRangeNaive(tree::TKD, range::Float64, queryPoint::Array{Float64}) where {TKD}
  H = BinaryHeap{typeof(tree.root), Float64}(64)

  # add all nodes in tree to a heap
  kdFindWithinRangeInSubtreeNaive(tree.distanceFunction, tree.root, queryPoint, H, range)

  return cleanHeapB(H)
end


# adds the node to the list if it is not already there
function addToRangeList(S::Tlist, thisNode::T, key::Float64) where {Tlist, T}
  if thisNode.inHeap
    return
  end
  thisNode.inHeap = true      # note that inHeap is a misnomer since this is a list
  JlistPush(S, thisNode, key)
end

# pops the range list
function popFromRangeList(S::Tlist) where {Tlist}
  (thisNode, key) = JlistPopKey(S)
  thisNode.inHeap = false     # note that inHeap is a misnomer since this is a list

  return (thisNode, key)
end

# empty the range list
function emptyRangeList(S::Tlist) where {Tlist}
  while S.length > 0
    (thisNode, key) = JlistPopKey(S)
    thisNode.inHeap = false  # note that inHeap is a misnomer since this is a list
  end
end

# empty the range list and print its key values
function emptyAndPrintRangeList(S::Tlist) where {Tlist}
  while S.length > 0
    node = popFromRangeList(S)
    println(node.data)
  end
end

# finds all nodes within range of the queryPoint in the subtree starting at root
# and also their distance squared, note that this data is stored in the nodeList
# the nodeList may also contain nodes before this function is called
function kdFindWithinRangeInSubtree(distanceFunction::Function,
  root::T, range::Float64, queryPoint::Array{Float64}, nodeList::ST) where {T,ST}

  # walk down the tree as if the node would be inserted
  parent::T = root
  while true
    if queryPoint[parent.kdSplit] < parent.position[parent.kdSplit]
      # traverse tree to the left
      if !parent.kdChildLExist
        # the queryPoint would be inserted as the left child of the parent
        break
      end
      parent = parent.kdChildL
      continue
    else
      # traverse tree to the right
      if !parent.kdChildRExist
        # the queryPoint would be inserted as the right child of the parent
        break
      end
      parent = parent.kdChildR
      continue
    end
  end

  #println(queryPoint)

#println(parent.position)

  newDist::Float64 = distanceFunction(queryPoint,parent.position)
  if newDist < range
    addToRangeList(nodeList, parent, newDist)
  end

  # now walk back up the tree (will break out when done)
  while true

    # now check if there could possibly be any nodes on the other
    # side of the parent within range, if not then check grandparent etc.

    parentHyperPlaneDist = queryPoint[parent.kdSplit] - parent.position[parent.kdSplit]

    if parentHyperPlaneDist >  range
      # then there could not be any closer nodes within range on the other
      # side of the parent (and the parent itself is also too far away

      if parent == root
        # the parent is the root and we are done
        return
      end

      parent = parent.kdParent
      continue
    end

    # if we are here, then there could be a closer node on the other side of the
    # parent (including the parent itself) that is within range

    # first check the parent itself (if it is not already one of the closest nodes)
    if !parent.inHeap    # note that inHeap is a misnomer since this is a list
      newDist = distanceFunction(queryPoint,parent.position)
      if newDist < range
        addToRangeList(nodeList, parent, newDist)
      end
    end

    # now check on the other side of the parent
    if queryPoint[parent.kdSplit] < parent.position[parent.kdSplit] && parent.kdChildRExist
      # queryPoint is on the left side of the parent, so we need to look
      # at the right side of it (if it exists)
      kdFindWithinRangeInSubtree(distanceFunction, parent.kdChildR, range, queryPoint, nodeList)
    elseif parent.position[parent.kdSplit] <= queryPoint[parent.kdSplit] && parent.kdChildLExist
      # queryPoint is on the right side of the parent, so we need to look
      # at the left side of it (if it exists)
      kdFindWithinRangeInSubtree(distanceFunction, parent.kdChildL, range, queryPoint, nodeList)
    end

    if parent == root
      # the parent is the root and we are done
      return
    end

    parent = parent.kdParent
  end
end


# returns all nodes within range of queryPoint and also their distance
# they are returned in a list with elements of type "(node, Float64)"
function kdFindWithinRange(tree::TKD, range::Float64, queryPoint::Array{Float64}) where {TKD}

  # make a list of type "(node, Float64)"
  L = JList{typeof(tree.root)}()

  # insert root node in list if it is within range
  distToRoot::Float64 = tree.distanceFunction(queryPoint, tree.root.position)
  if distToRoot <= range
    addToRangeList(L, tree.root, distToRoot)
  end

  # find nodes within range
  kdFindWithinRangeInSubtree(tree.distanceFunction, tree.root, range, queryPoint, L)

  if tree.numWraps > 0
    # if dimensions wrap around, we need to search vs. identities (ghosts)

    pointIterator = ghostPointIterator{TKD}(tree, queryPoint)
    while true
      thisGhostPoint = getNextGhostPoint(pointIterator, range)
      if thisGhostPoint == nothing
        break
      end

      # now see if any points in the space are closer to this ghost
      kdFindWithinRangeInSubtree(tree.distanceFunction, tree.root, range, thisGhostPoint, L)
    end
  end

  return L
end

# returns all nodes within range of queryPoint and also their distance
# they are returned in a list with elements of type "(node, Float64)"
# the list to be used in passed in, so that additional points can be added
# to it (e.g., if we want to have one list containing the points that
# are close to a couple of different points X1 .. Xn, then call this
# for X2 ... Xn after first calling kdFindWithinRange for X1
function kdFindMoreWithinRange(tree::TKD, range::Float64,
  queryPoint::Array{Float64}, L::TL) where {TKD, TL}

  # insert root node in list if it is within range
  distToRoot::Float64 = tree.distanceFunction(queryPoint, tree.root.position)
  if distToRoot <= range
    addToRangeList(L, tree.root, distToRoot)
  end

  # find nodes within range
  kdFindWithinRangeInSubtree(tree.distanceFunction, tree.root, range, queryPoint, L)

  if tree.numWraps > 0
    # if dimensions wrap around, we need to search vs. identities (ghosts)

    pointIterator = ghostPointIterator{TKD}(tree, queryPoint)
    while true
      thisGhostPoint = getNextGhostPoint(pointIterator, range)
      if thisGhostPoint == nothing
        break
      end

      # now see if any points in the space are closer to this ghost
      kdFindWithinRangeInSubtree(tree.distanceFunction, tree.root, range, thisGhostPoint, L)
    end
  end

  return L
end

# inserts a new point into the tree (used only for debugging)
function kdInsert(tree::TKD, A::Array{Float64}) where {TKD}
  N = (typeof(tree.root))()
  N.position = A
  kdInsert(tree::TKD, N)
end

# examples of a few distance functions:

# returns distance between two points from the space of a cartesian product of

# 2 two dimensional point robots
oneTwoDRobotDist(poseA::Array{T}, poseB::Array{T}) where {T} = (
sqrt((poseA[1]-poseB[1])^2 + (poseA[2]-poseB[2])^2))

# 2 two dimensional point robots
twoTwoDRobotDist(poseA::Array{T}, poseB::Array{T}) where {T} = (
sqrt((poseA[1]-poseB[1])^2 + (poseA[2]-poseB[2])^2) +
sqrt((poseA[3]-poseB[3])^2 + (poseA[4]-poseB[4])^2))

# returns distance between two points from the space of a cartesian product of
# 3 two dimensional point robots
threeTwoDRobotDist(poseA::Array{T}, poseB::Array{T}) where {T} = (
sqrt((poseA[1]-poseB[1])^2 + (poseA[2]-poseB[2])^2) +
sqrt((poseA[3]-poseB[3])^2 + (poseA[4]-poseB[4])^2) +
sqrt((poseA[5]-poseB[5])^2 + (poseA[6]-poseB[6])^2))

# returns distance between two points from the space of a cartesian product of
# 4 two dimensional point robots
fourTwoDRobotDist(poseA::Array{T}, poseB::Array{T}) where {T} = (
sqrt((poseA[1]-poseB[1])^2 + (poseA[2]-poseB[2])^2) +
sqrt((poseA[3]-poseB[3])^2 + (poseA[4]-poseB[4])^2) +
sqrt((poseA[5]-poseB[5])^2 + (poseA[6]-poseB[6])^2) +
sqrt((poseA[7]-poseB[7])^2 + (poseA[8]-poseB[8])^2))

# returns distance between two points from the space of a cartesian product of
# 5 two dimensional point robots
fiveTwoDRobotDist(poseA::Array{T}, poseB::Array{T}) where {T} = (
sqrt((poseA[1]-poseB[1])^2 + (poseA[2]-poseB[2])^2) +
sqrt((poseA[3]-poseB[3])^2 + (poseA[4]-poseB[4])^2) +
sqrt((poseA[5]-poseB[5])^2 + (poseA[6]-poseB[6])^2) +
sqrt((poseA[7]-poseB[7])^2 + (poseA[8]-poseB[8])^2) +
sqrt((poseA[9]-poseB[9])^2 + (poseA[10]-poseB[10])^2))

## distance in a 3Dimension space that wraps in all 3 dimensions
# Note if using this distance, then the dist function actually
# used for kd-tree search should be the 3 dimensional euclidian metric
# since wrapping is dealt with using ghost points
# (i.e., the next function down)
# note that 1.0 in the following gunction means that the length of each
# dimension is 1.0
S3Dist(poseA::Array{T}, poseB::Array{T}) where {T} = (sqrt(
min(abs(poseA[1]-poseB[1]), min(poseA[1], poseB[1])+1.0-max(poseA[1], poseB[1]))^2 +
min(abs(poseA[2]-poseB[2]), min(poseA[2], poseB[2])+1.0-max(poseA[2], poseB[2]))^2 +
min(abs(poseA[3]-poseB[3]), min(poseA[3], poseB[3])+1.0-max(poseA[3], poseB[3]))^2 ))


S3KDSearchDist(poseA::Array{T}, poseB::Array{T}) where {T} = (sqrt(
(poseA[1]-poseB[1])^2 +
(poseA[2]-poseB[2])^2 +
(poseA[3]-poseB[3])^2 ))


## distance in a R2 X S
# Note if using this distance, then the dist function actually
# used for kd-tree search should be the 3 dimensional euclidian metric
# since wrapping is dealt with using ghost points
# (i.e., the next function down)
# note that 1.0 in the following gunction means that the length of
# wrap dimension is 1.0
R2SDist(poseA::Array{T}, poseB::Array{T}) where {T} = (sqrt(
(poseA[1]-poseB[1])^2 +
(poseA[2]-poseB[2])^2 +
min(abs(poseA[3]-poseB[3]), min(poseA[3], poseB[3])+1.0-max(poseA[3], poseB[3]))^2 ))


R2SKDSearchDist(poseA::Array{T}, poseB::Array{T}) where {T} = (sqrt(
(poseA[1]-poseB[1])^2 +
(poseA[2]-poseB[2])^2 +
(poseA[3]-poseB[3])^2 ))


#function testCase()
#  T=KDTree{KDTreeNode{Float64}}(6, threeTwoDRobotDist)
#
#tic()
#  for j = 1:100
#    for i = 1:10000
#      kdInsert(T, rand(6,1))
#    end
#
#    queryPoint = rand(6,1)
#    (nearestPoint , nearestDist) = kdFindNearest(T, queryPoint)
#    (nearestPointNaive , nearestDistNaive) = kdFindNearestNaive(T, queryPoint)
#
#    if nearestDist != nearestDistNaive
#      println(queryPoint)
#      println(nearestPoint.position)
#      println(nearestPointNaive.position)
#
#      println(nearestDist)
#      println(nearestDistNaive)
#      println("problems")
#      return
#    end
#
#   nodesNaive = kdFindKNearestNaive(T, 10, queryPoint)
#    nodes = kdFindKNearest(T, 10, queryPoint)
#
#    if nodesNaive[1,1] != nodes[1,1]
#      println(nodesNaive[1,1].data)
#      println(nodes[1,1].data)
#      println("problems")
#      return
#    end
#
#    nodesNaive = kdFindWithinRangeNaive(T, 0.1, queryPoint)
#    L = kdFindWithinRange(T, 0.1, queryPoint)
#
#    if size(nodesNaive,1) != L.length
#      println(nodesNaive[1,1].data)
#      emptyAndPrintRangeList(L)
#      println("problems: $(size(nodesNaive,1))  $(L.length)")
#      return
#    end
#
#    emptyRangeList(L)
#
#  end
#toc()
#end

#function testGhost()
#  total = 1
#  for t = 1:total
#    wraps = [3]
#    wrapPoints = [1.0]
#    T = KDTree{KDTreeNode{Float64}}(3, R2SKDSearchDist, wraps, wrapPoints)
#    TNaive = KDTree{KDTreeNode{Float64}}(3, R2SDist, wraps, wrapPoints)
#
#    fptr = open("temp/points.txt", "w")
#    for i = 1:500
#      p = rand(1,3)
#      kdInsert(T, p)
#      kdInsert(TNaive, p)
#      writecsv(fptr, p)
#    end
#
#    queryPoint = rand(1,3)
#    queryPoint = [0.5 0.5 0.0 ]
#    writecsv(fptr, queryPoint)
#    close(fptr)
#
# #   (nearestPoint , nearestDist) = kdFindNearest(T, queryPoint)
# #   (nearestPointNaive , nearestDistNaive) = kdFindNearestNaive(TNaive, queryPoint)
#
# #  if nearestPoint.position != nearestPointNaive.position
# #     println(" ----------- problems -------------    ")
# #    println(queryPoint)
# #     println(nearestPoint.position)
# #     println(nearestPointNaive.position)
# #     println(nearestDist)
# #     println(nearestDistNaive)
# #     println("----")
#
# #     println("$(S3Dist(nearestPoint.position, queryPoint))")
# #     println("$(S3Dist(nearestPointNaive.position, queryPoint))")
# #   end
#
#    nodesNaive = kdFindWithinRangeNaive(TNaive, 0.2, queryPoint)
#    L = kdFindWithinRange(T, 0.2, queryPoint)
#
#    if size(nodesNaive,1) != L.length
#      println(nodesNaive[1,1].data)
#      emptyAndPrintRangeList(L)
#      println("problems: $(size(nodesNaive,1))  $(L.length)")
#      return
#    end
#
#
#    fptr = open("temp/closestPoints.txt", "w")
#    thisNode = L.front
#    while thisNode != thisNode.child
#      writecsv(fptr, thisNode.data.position)
#      thisNode = thisNode.child
#    end
#    close(fptr)
#
#    emptyRangeList(L)
#  end
#
#end
