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


### this is an a simple edge between two nodes of type T
### ALL OTHER EDGE TYPES MUST HAVE THESE FIELDS, INCLUDING DEFAULT CONSTRUCTOR.
### FURTHERMORE, TYPEALIASING SHOULD BE USED TO MAKE Edge{T} <: OF WHATEVER
### EDGE IS ACTUALLY BEING USED. SEE ADDITIONAL NOTES BELOW REGARDING newEdge()
### this particular (example) type should be used when trajectories between nodes
### follow straigh-line paths between those nodes in the C-space
mutable struct SimpleEdge{T}
  startNode::T
  endNode::T
  dist::Float64  # the distance between startNode and endNode (i.e., how far the
                 # robot must travel through the -configuration- space to get from
                 # startNode to endNode. This is the distance that is used to
                 # calculate RRTTreeCost and RRTLMC

  distOriginal::Float64 # saves the original value of dist, so we don't need to
                        # recalculate if this edge is removed and then added again

  listItemInStartNode::JListNode{SimpleEdge{T}} # pointer to this edges location
                                                # in startNode

  listItemInEndNode::JListNode{SimpleEdge{T}}   # pointer to this edges location
                                                # in endNode

  Wdist::Float64 # this contains the distance that the robot must travel through
                 # the -workspace- along the edge (so far only used for time based
                 # c-spaces)

  SimpleEdge{T}() where {T} = new{T}()
end
