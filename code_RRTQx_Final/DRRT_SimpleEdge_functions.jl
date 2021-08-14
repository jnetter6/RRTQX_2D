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



######################### !!!!! critical functions !!!!! ##########################
###  !  !  !  !  !  !  !  !  !  !  !  !  !  !  !  !  !  !  !  !  !  !  !  !  !  ###
## ! functions that must be modified vs. the particular C-Space and Workspace  ! ##
## ! that are being used !!!                                                   ! ##
###  !  !  !  !  !  !  !  !  !  !  !  !  !  !  !  !  !  !  !  !  !  !  !  !  !  ###
###################################################################################


## returns the distance between two points, assuming they are in the C-Space, this
## should obey the triangle inequality, it is used for
#
# use the following version of dist for Euclidian space
dist(x::Array,y::Array) = euclidianDist(x,y)


## This returns the distance that is used internally in the kd-tree (see notes
## there!!! the kd-implimentation handles wrapping dimensions (e.g. a circle or
## taurus) by running multipe searches, 1 for each wrapped identity of the
## querery point. So this is the distance function that is the -non-wrapping-
## version of what the actual distance function between points in the KD space is
#
# use the following version of KDdist for Euclidian space (with and without time)
KDdist(x::Array,y::Array) = euclidianDist(x,y)


## returns the workspace distance between two points, this should obey the triangle
## inequality. e.g., in the current version of this code, it is assumed that the
## workspace is the first two dimensions of the C-space, it is used for calculating
## the "real world" distance that the robot covers (e.g., in a particualr ammount
## of time when determining if a move is possible given robot max speed).
#
# use the following version of Wdist for Euclidian space and Dubbin's space
Wdist(x::Array,y::Array) = euclidianDist(x[1:2],y[1:2])



# moves newPoint toward closestPoint such that each robot is no
# further than delta, points repesent the cartesian product of R robots
#
# use the following version of saturate for Euclidian space
function saturate(newPoint::Array{Float64}, closestPoint::Array{Float64}, delta)
  thisDist = dist(newPoint, closestPoint)
  if thisDist > delta
    newPoint = closestPoint  + (newPoint - closestPoint)*delta/thisDist
  end
end


################################### edge functions ################################
### Functions that interact with edges and not much else.                       ###
###################################################################################

### allocates a new edge
function newEdge(startNode::T, endNode::T) where {T}
  E = Edge{T}()
  E.startNode = startNode
  E.endNode = endNode
  return E
end



# returns true if the dynamics of the robot and space allow a robot to follow the edge
#
## SimpleEdge version
function validMove(S::TS, edge::SimpleEdge{T}) where {T, TS}
  if S.spaceHasTime
    # note that planning happens in reverse time, i.e., time = 0 is at the root
    # of the search tree, and thus the time of startNode must be greater than
    # the time of the end node.
    return (edge.Wdist <= (edge.startNode.position[3] - edge.endNode.position[3])*S.robotVelocity)
  end

  # if space does not have time then we assume that the move is always valid
  return true
end



### returns the pose of a robot that is located dist along the edge
### NOTE THAT 'dist' and 'far' are with respect to whatever type
### of distance is stored in edge.dist
#
## simple edge version
function poseAtDistAlongEdge(edge::SimpleEdge, distAlongEdge::Float64)
  if edge.dist == 0.0
    return copy(edge.endNode.position)
  end

  ratioAlongEdge = distAlongEdge/edge.dist
  ret = edge.startNode.position + ratioAlongEdge*(edge.endNode.position - edge.startNode.position)
  return ret
end



### returns the pose of a robot that is located time along the edge
#
## simple edge version
function poseAtTimeAlongEdge(edge::SimpleEdge, timeAlongEdge::Float64)
  ratioAlongEdge = timeAlongEdge/(edge.startNode.position[3] - edge.endNode.position[3])
  ret = edge.startNode.position + ratioAlongEdge*(edge.endNode.position - edge.startNode.position)
  return ret
end


### the next function saves the last part of the trajectory into the file filePtr,
### starting distFromFront from the front. NOTE THAT 'distFromFront' is with
### respect to whatever type of distance is stored in edge.dist
#
## Simple Edge
function saveEndOfTrajectory(filePtr::TFP, edge::SimpleEdge, distFromFront::Float64) where {TFP}
  ratioAlongEdge = distFromFront/edge.dist
  ret = edge.startNode.position + ratioAlongEdge*(edge.endNode.position - edge.startNode.position)

  writedlm(filePtr, ret, ',')
  writedlm(filePtr, edge.endNode.position, ',')
end


### the next function saves the last part of the trajectory into the file filePtr,
# starting timeFromFront from the front
#
## simple edge version
function saveEndOfTrajectoryTime(filePtr::TFP, edge::SimpleEdge,
  timeFromFront::Float64) where {TFP}
  ratioAlongEdge = timeFromFront/(edge.startNode.position[3] - edge.endNode.position[3])
  ret = edge.startNode.position + ratioAlongEdge*(edge.endNode.position - edge.startNode.position)

  writedlm(filePtr, ret, ',')
  writedlm(filePtr, edge.endNode.position, ',')
end




### NOTE: A version of calculateTrajectory() must be created for whatever edge
### type is being used. FURTHERMORE, the trajectory can be stored in any usefull
### form, but there are a number of functions (e.g., for collision checking) that
### must be modified to use whatever form the trajectory is stored in. For example,
### For a holonomic robot in euclidian space this is simply a straight line and
### the trajectory is implicitly defined by the edge between the nodes. For more
### complex systems the trajectory is not a straight line and is stored as a local
### path (called a trajectory) within the edge. AT THE VERY LEAST this function
### should populate the dist field of edge
#
## basic version, trajectory is implicitly defined as straight path along line
## segment between edge.startNode and edge.endNode, Euclidian space is assumed
function calculateTrajectory(S::TS, edge::SimpleEdge) where {TS}
  edge.dist = dist(edge.startNode, edge.endNode)
  edge.distOriginal = edge.dist
  edge.Wdist = Wdist(edge.startNode, edge.endNode)
end


### this calculates a trajectory of what the robot is supposed to do when it is
### hovering "in place"
#
## simple edge version
calculateHoverTrajectory(S::TS, edge::SimpleEdge) where {TS} = calculateTrajectory(S, edge)

### this saves the trajectory that is stored in the edge to the file
#
# simple edge version
function saveEdgeTrajectory(filePtr::TFP, edge::SimpleEdge) where {TFP}
  writedlm(filePtr, edge.startNode.position, ',')
  writedlm(filePtr, edge.endNode.position, ',')
end



########################## collision checking functions ###########################
### Note, these are only collision checking functions that depend on edge type, ###
### more general collision checking functions appear in DRRT.jl                 ###
###################################################################################


### checks if the edge is in collision with a particular obstacle
### (one version for each edge type)
#
## SimpleEdge version
explicitEdgeCheck(S::CSpace{T}, edge::SimpleEdge,
obstacle::OT) where {T, OT} = explicitEdgeCheck2D(obstacle,
edge.startNode.position, edge.endNode.position, S.robotRadius)
