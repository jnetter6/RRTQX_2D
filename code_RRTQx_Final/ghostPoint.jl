function a()
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

# this contains an iterator that efficiently gives the next ghost Point needed
# for searching in a kd tree that has some dimensions the wrap around
# note it starts at the first ghost and does not return the original point

include("heap.jl")

mutable struct ghostPointIterator{TKD}
  kdTree::TKD                # the kd tree hold that this is being used with

  queryPoint::Array{Float64} # the actual point that all the ghosts are "identical" with

  wrapDimFlags::Array{Int}   # wrapDimFlags indicates the current permutation
                             # of all ghost points that we are currently looking at

  ghostTreeDepth::Int        # a pointer to the current depth of the "gost tree"
                             # note that this "tree" is only a theoretical construct that
                             # determines the order in which the ghosts are returned
                             # it should not be confused with the kdTree

  currentGhost::Array{Float64} # the current ghost that we are returning

  closestUnwrappedPoint::Array{Float64} # the closest point in the normal space to the
                                        # currentGhost, dist between this and ghost can
                                        # be used as uristic to skip unhelpfull ghosts


  ghostPointIterator{TKD}(kdTree::TKD, queryPoint::Array{Float64}) where {TKD} =
  new{TKD}(kdTree, queryPoint, zeros(Int, kdTree.numWraps), kdTree.numWraps,
  copy(queryPoint), copy(queryPoint))

end

# this returns the next ghost point, note it starts at the first -ghost- and does
# not return the original point
function getNextGhostPoint(G::TG, bestDist::Float64) where {TG}

  while true # will return out when done
    # go up tree until we find a wrapDimFlags[ghostTreeDepth] == 0 (this
    # indicates that we need to try permutations where ghost is wrapped
    # around tree.wraps[ghostTreeDepth]

    while G.ghostTreeDepth > 0 && G.wrapDimFlags[G.ghostTreeDepth] != 0
      G.ghostTreeDepth -= 1
    end

    if G.ghostTreeDepth == 0
      # we are finished, no more ghosts
      return nothing
    end

    # otherwise we are at at depth where wrapDimFlags[ghostTreeDepth] == 0
    G.wrapDimFlags[G.ghostTreeDepth] = 1

    # calculate this (wrapped) dimension of the ghost
    dimVal = G.queryPoint[G.kdTree.wraps[G.ghostTreeDepth]]
    dimClosest = 0.0
    if G.queryPoint[G.kdTree.wraps[G.ghostTreeDepth]] < G.kdTree.wrapPoints[G.ghostTreeDepth]/2.0
      # wrap to the right
      dimVal += G.kdTree.wrapPoints[G.ghostTreeDepth]
      dimClosest += G.kdTree.wrapPoints[G.ghostTreeDepth]
    else
      # wrap to the left
      dimVal -= G.kdTree.wrapPoints[G.ghostTreeDepth]
    end
    G.currentGhost[G.kdTree.wraps[G.ghostTreeDepth]] = dimVal
    G.closestUnwrappedPoint[G.kdTree.wraps[G.ghostTreeDepth]] = dimClosest

    # finally, move back down the tree to the left-most possible leaf,
    # marking the path with 0s, and populating approperiate dimension of
    # ghost point with values
    while G.ghostTreeDepth < G.kdTree.numWraps
      G.ghostTreeDepth += 1
      G.wrapDimFlags[G.ghostTreeDepth] = 0
      G.currentGhost[G.kdTree.wraps[G.ghostTreeDepth]] = G.queryPoint[G.kdTree.wraps[G.ghostTreeDepth]]
      G.closestUnwrappedPoint[G.kdTree.wraps[G.ghostTreeDepth]] = G.currentGhost[G.kdTree.wraps[G.ghostTreeDepth]]
    end

    # check if closest point in unwrapped space is further than best distance
    if G.kdTree.distanceFunction(G.closestUnwrappedPoint, G.currentGhost) > bestDist
      continue
    end

    return copy(G.currentGhost)
  end

end
