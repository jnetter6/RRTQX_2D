function makeMyTextEditorDisplayNice2() # thanks!
  return
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



### an edge for a Dubins state space
mutable struct DubinsEdge{T}
  startNode::T
  endNode::T
  dist::Float64  # the distance between startNode and endNode (i.e., how far the
                 # robot must travel through the -configuration- space to get from
                 # startNode to endNode. This is the distance that is used to
                 # calculate RRTTreeCost and RRTLMC

  distOriginal::Float64 # saves the original value of dist, so we don't need to
                        # recalculate if this edge is removed and then added again


  listItemInStartNode::JListNode{DubinsEdge{T}} # pointer to this edges location
                                                 # in startNode

  listItemInEndNode::JListNode{DubinsEdge{T}}   # pointer to this edges location
                                                 # in endNode

  Wdist::Float64 # this contains the distance that the robot must travel through
                 # the -workspace- along the edge (so far only used for time based
                 # c-spaces)


  dubinsType::String     # one of the following types depending on edge:
                               # lsl, rsr, lsr, rsl, lrl, rlr


  trajectory::Array{Float64,2} # stores a descritized version of the dubins path
                               # [x y] or [x y time] depending of it time is being used

  velocity::Float64 # the velocity that this robot travels along this edge
                    # only used if time is part of the state space

  DubinsEdge{T}() where {T} = new{T}()
end
