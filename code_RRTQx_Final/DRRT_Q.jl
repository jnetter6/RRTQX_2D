using DelimitedFiles

function makeMyTextEditorDisplayNice() # thanks!
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


### This version of the code uses an explicit notion of edges, which allows
### edges to represent arbitrary trajectories. It also uses directed edges,
### which allows non symmetric connection between nodes


################################# simple functions ################################
### random simple functions.                                                    ###
###################################################################################


# saves the data to the file
function saveData(data, fileName)
  fptr = open(fileName, "w")

  for i = 1:size(data,1)
    writedlm(fptr, hcat(data[i,:]'), ',')
  end
  close(fptr)
end

# prints the length of the goal's rrtLMC
function printPathLengths(goalNode::T) where {T}
  print("\n goal.rrtLMC: $(goalNode.rrtLMC) goal.rrtTreeCost: $(goalNode.rrtLMC) \n")
  println()
end

# design is limited in scope to help with the readObstaclesFromfile() function
function str2array(S::String)

  # count number of elements to put in the return structure
  # ignoring special cases of missing numbers or extra commas
  # and assumes no comma after the last number
  p = 1
  L = 0
  Lmax = length(S)

  while p <= Lmax
    if p == Lmax || S[p] == ',' || S[p] == '\n' || S[p] == '\0'
      L += 1
      if S[p] != ','
        break
      end
    end
    p += 1
  end

  p = 1
  q = 1
  i = 1

  retVal = Array{Float64}(undef, 1, L)
  while i <= L
    while q < Lmax && S[q] != ',' && S[q] != '\n' && S[q] != '\0'
      q += 1
    end
    if S[q] == ',' || S[q] == '\n' || S[q] == '\0'
      retVal[i] = parse(Float64, S[p:q-1])
      p = q = q+1
    else
      retVal[i] = parse(Float64, S[p:q])
    end

    i += 1

  end
  return retVal
end


################################### node functions ################################
### Functions that interact with nodes and not much else. This includes a       ###
### handful of functions that traverse the search tree in one way or another.   ###
###################################################################################

# returns distance between two nodes (e.g., in the C-space), this should
# obey the triangle inequality
function dist(x::RRTNode,y::RRTNode)
  dist(x.position, y.position)
end

# returns distance between two nodes (e.g., in the workspace), this should
# obey the triangle inequality. Note that it is assumed that the first two
# dimensions of the C-space are position in the workspace.
function Wdist(x::RRTNode,y::RRTNode)
  Wdist(x.position, y.position)
end


# returns true if nodeA certifies nodeB (only used in Alg A)
function certifierOf(nodeA::T, nodeB::T) where {T}

  certifyingNode = nodeB
  if !nodeB.hasCertificate
    certifyingNode = nodeB.certifyingNode
  end

  if certifyingNode == nodeA
    return true
  end


  if nodeA.hasCertificate && Wdist(nodeA.position, nodeB.position) < nodeA.certificateValue
    return true
  end
  return false
end


# prints each node in the subtree of the KD-tree followed by its
# parent in the RRT-tree, assumes that fptr is already opened
function saveRRTSubTree(node::T, root::T, fptr::IOStream) where {T}
   if node.rrtParentUsed
     writedlm(fptr, [node.position node.rrtTreeCost], ',')
     writedlm(fptr, [node.rrtParentEdge.endNode.position  node.rrtParentEdge.endNode.rrtTreeCost], ',')
   end

   if node.kdChildLExist
     saveRRTSubTree(node.kdChildL, root, fptr::IOStream)
   end
   if node.kdChildRExist
     saveRRTSubTree(node.kdChildR, root, fptr::IOStream)
   end
end


# saves the RRT tree to a file (in edge form)
# it uses the KD-Tree structure to traverse the tree
function saveRRTTree(tree::TKD, fileName) where {TKD}
    fptr = open(fileName, "w")
    if tree.treeSize > 0
      saveRRTSubTree(tree.root, tree.root, fptr)
    end
    close(fptr)
end


# prints each RRT-graph edge from the sub-tree of the
# KD-tree, assumes that fptr is already opened
function saveRRTSubGraph(node::T, fptr::IOStream) where {T}

   listItem = node.rrtNeighborsOut.front
   for i = 1:node.rrtNeighborsOut.length
     writedlm(fptr, node.position, ',')
     writedlm(fptr, listItem.data.position, ',')

     listItem = listItem.child   # iterate thorugh list
   end
   if node.kdChildLExist
     saveRRTSubGraph(node.kdChildL, fptr::IOStream)
   end
   if node.kdChildRExist
     saveRRTSubGraph(node.kdChildR, fptr::IOStream)
   end
end


# saves all RRT edges to a file (in edge form)
# it uses the KD-Tree structure to traverse the tree
function saveRRTGraph(tree::TKD, fileName) where {TKD}
    fptr = open(fileName, "w")
    if tree.treeSize > 0
      saveRRTSubGraph(tree.root, fptr)
    end
    close(fptr)
end


# prints each node in the subtree of the KD-tree,
# assumes that fptr is already opened
function saveRRTSubNodes(node::T, fptr::IOStream) where {T}
   writedlm(fptr, [node.position node.rrtTreeCost node.rrtLMC], ',')

   if node.kdChildLExist
     saveRRTSubNodes(node.kdChildL, fptr::IOStream)
   end
   if node.kdChildRExist
     saveRRTSubNodes(node.kdChildR, fptr::IOStream)
   end
end


# saves the RRT nodes to a file
# it uses the KD-Tree structure to traverse the tree
function saveRRTNodes(tree::TKD, fileName) where {TKD}
    fptr = open(fileName, "w")
    if tree.treeSize > 0
      saveRRTSubNodes(tree.root, fptr)
    else
      println("tree size = 0")
    end
    close(fptr)
end


# prints each node that is in collision in the subtree of the KD-tree,
# assumes that fptr is already opened
function saveRRTSubNodesCollision(node::T, fptr::IOStream) where {T}
  writedlm(fptr, [node.position min(node.rrtTreeCost, node.rrtLMC)], ',')

  if node.kdChildLExist
    saveRRTSubNodesCollision(node.kdChildL, fptr::IOStream)
  end
  if node.kdChildRExist
    saveRRTSubNodesCollision(node.kdChildR, fptr::IOStream)
  end
end


# saves the RRT nodes to a file, uses the KD-Tree structure to traverse the tree
function saveRRTNodesCollision(tree::TKD, fileName) where {TKD}
  fptr = open(fileName, "w")
  if tree.treeSize > 0
    saveRRTSubNodesCollision(tree.root, fptr)
  else
    println("tree size = 0")
  end
  close(fptr)
end


# saves the path (sequence of trajectories from edges) between robot and the first
# target node and then between the first target node and the root
function saveRRTPath(S::TS, node::T, root::T, robot::RobotData{T}, fileName) where {T, TS}
  fptr = open(fileName, "w")
  thisNode = node

  if robot.robotEdgeForPlottingUsed && robot.robotEdgeForPlotting.startNode != thisNode
    # save trajectory from first edge (connecting robot to graph)
    # (unless the entire thing will be handled as part of nextNode)
    if S.spaceHasTime
      saveEndOfTrajectoryTime(fptr, robot.robotEdgeForPlotting, robot.timeAlongRobotEdgeForPlotting)
    else
      saveEndOfTrajectory(fptr, robot.robotEdgeForPlotting, robot.distAlongRobotEdgeForPlotting)
    end
  end

  i = 0 # added to detect for inf paths
  while thisNode != root && thisNode.rrtParentUsed && i < 1000
    saveEdgeTrajectory(fptr, thisNode.rrtParentEdge)
    thisNode = thisNode.rrtParentEdge.endNode
    i += 1 # added to detect for inf paths
  end

  if S.spaceHasTime
    writedlm(fptr, hcat(thisNode.position[1,1:3]'), ',')
  else
    writedlm(fptr, hcat(thisNode.position[1,1:2]'), ',')
  end
  close(fptr)
end


# prints each node in the subtree of the KD-tree that has a certificate
# followed by its certificate
function saveCertificates(node::T, fptr::IOStream) where {T}
   if node.hasCertificate
       writedlm(fptr, [node.position], ',')
   end

   if node.kdChildLExist
     saveCertificates(node.kdChildL, fptr::IOStream)
   end
   if node.kdChildRExist
     saveCertificates(node.kdChildR, fptr::IOStream)
   end
end


# saves certificates to a file (along with nodes they are from)
function saveCertificates(tree::TKD, fileName) where {TKD}
    fptr = open(fileName, "w")
    if tree.treeSize > 0
      saveCertificates(tree.root, fptr)
    end
    close(fptr)
end


# prints each node in the subtree of the KD-tree that has a certificate
# followed by its certificate
function saveCertificateEdges(node::T, fptr::IOStream) where {T}
   if(!node.hasCertificate)
     writedlm(fptr, [node.position node.certifyingNode.position], ',')
   end

   if node.kdChildLExist
     saveCertificateEdges(node.kdChildL, fptr::IOStream)
   end
   if node.kdChildRExist
     saveCertificateEdges(node.kdChildR, fptr::IOStream)
   end
end


# saves certificate pointer edges to a file
function saveCertificateEdges(tree::TKD, fileName) where {TKD}
    fptr = open(fileName, "w")
    if tree.treeSize > 0
      saveCertificateEdges(tree.root, fptr)
    end
    close(fptr)
end

# extracts the cost of the graph path from the node to the start by adding up
# edge lengths between node and the root node
function extractPathLength(node::T, root::T) where {T}
  #println("extracting path length")

  pathLength = 0.0
  thisNode = node
  while thisNode != root
    if !thisNode.rrtParentUsed
      pathLength = Inf
      break
    end

    pathLength += thisNode.rrtParentEdge.dist
    thisNode = thisNode.rrtParentEdge.endNode
  end

  #println("done extracting path length")
  return pathLength
end



################################# obstacle functions ##############################
### This does -not- include collision checking, which appears lower in its own  ###
### section. Functions involving obstacles that require CSpace access also      ###
### appear lower down in the CSpace section or in the RRTx section.             ###
###################################################################################


# decreases the life of the obstacle
decreaseLife(O::Obstacle) = (O.lifeSpan -= 1.0)

# for obstacle type 7 (time obstacles with paths unknown to the robot)
# this is used to calculate the current path (used for collision checking)
# from unknownPath (which describes how the obstace moves vs time)
# based on the current time. Note that times in the future are closer to S.start[3]
function changeObstacleDirection(S::TS, O::Obstacle, currentTime::Float64) where {TS}

  endTime = S.start[3]
  pathTimeStep = 3.0 # edges in path are no longer than this in the time dimension

  # find the path segment of unknownPath that the obstacle is now moving along
  # note that future is represented by times closer to S.start[3] and unknwonPath
  # is stored from the future (low) to the past (high)
  # (yes, this is nonintuitive, perhaps I will change the convention someday)
  while (O.unknownPath[O.nextDirectionChangeInd,3] > currentTime &&
         O.nextDirectionChangeInd > 1)
    O.nextDirectionChangeInd -= 1
  end

  # calculate the start and end points of the segment to place in path
  # and remember nextDirectionChangeTime
  highPoint = zeros(Float64, 3)
  lowPoint = zeros(Float64, 3)
  if (O.unknownPath[O.nextDirectionChangeInd, 3] <= currentTime &&
      O.nextDirectionChangeInd == size(O.unknownPath,1))
  #  # obstacle has not started moving yet, so we assume that it is at the start
  #  # of the "earliest" (greatest time) edge of unknownPath
  #
  #  highPoint =  [O.unknownPath[end,1:2] currentTime]
  #  lowPoint = [O.unknownPath[end,1:2] (currentTime - pathTimeStep)]
  #  O.nextDirectionChangeTime = O.unknownPath[end,3]

    # obstacle has not started moving yet so
    # assume that we know the first movement segment of the robot

    indWeCareAbout = size(O.unknownPath, 1)

    highPoint = O.unknownPath[indWeCareAbout,:]
    lowPoint = O.unknownPath[indWeCareAbout-1,:]
    O.nextDirectionChangeTime = O.unknownPath[indWeCareAbout-1,3]


  elseif (O.unknownPath[O.nextDirectionChangeInd, 3] > currentTime &&
      O.nextDirectionChangeInd <= 1)
    # time has progressed further than this obstacle has a path for
    # so we assume that it remains at the end of its path forever

    highPoint = vcat(O.unknownPath[1,1:2], currentTime)
    lowPoint =  vcat(O.unknownPath[1,1:2], (currentTime - pathTimeStep))
    O.nextDirectionChangeTime = -Inf
  else

    highPoint = O.unknownPath[O.nextDirectionChangeInd + 1,:]
    lowPoint = O.unknownPath[O.nextDirectionChangeInd,:]
    O.nextDirectionChangeTime = O.unknownPath[O.nextDirectionChangeInd,3]
  end


  # calculate path, which is a line parallel to the edge {highPoint,lowPoint}
  # and that goes from currentTime to endTime in the time dimension.
  # The required time must be calculated based on the obstac's speed.
  # (Saving as a series of segments increases collision checking speed because
  # we can ignore those that have no time overlap)

  mx = (highPoint[1] - lowPoint[1])/(highPoint[3] - lowPoint[3])
  my = (highPoint[2] - lowPoint[2])/(highPoint[3] - lowPoint[3])

  ts = collect(currentTime:-pathTimeStep:endTime)
  L = length(ts)

  O.path = zeros(Float64, length(ts),3)

  for t = 1:length(ts)
    O.path[t,1] = lowPoint[1] + (ts[L-t+1] - lowPoint[3])*mx
    O.path[t,2] = lowPoint[2] + (ts[L-t+1] - lowPoint[3])*my
    O.path[t,3] = lowPoint[3] + ts[L-t+1] - lowPoint[3]
  end

end



# saves obstacle data to a file
function saveObstacleLocations(obstacles::List{Obstacle}, fileName)

  fptr = open(fileName, "w")

  listNode = obstacles.front
  while listNode != listNode.child
    ob = listNode.data

    if ob.kind == 6 || ob.kind == 7
      # for time obstacles we save data about path [x, y, time, radius]
      # then we put NaNs rows between obstacls

      for i = 1:size(ob.path,1)
        writedlm(fptr, [(reshape(ob.path[i,:], 1, length(ob.path[i,:])) + [ob.position 0.0]) ob.radius], ',')
      end
      writedlm(fptr, [NaN NaN NaN NaN], ',')


      listNode = listNode.child
      continue
    end

    if ob.kind != 3 && ob.kind != 4  && ob.kind != 5
      println("warning cannot save non-polygon obstacle to file (not implimented) ")
      listNode = listNode.child
      continue
    end

    if ob.obstacleUnused
      listNode = listNode.child
      continue
    end

    for i = 1:size(ob.polygon,1)
      writedlm(fptr, reshape(ob.polygon[i,:], 1, length(ob.polygon[i,:])), ',')
    end
    writedlm(fptr, reshape(ob.polygon[1,:], 1, length(ob.polygon[1, :])), ',')
    writedlm(fptr, [NaN NaN], ',')
    listNode = listNode.child
  end

  close(fptr)

end


################################## CSpace functions ###############################
### functions that interact CSpace, this includes sampling functions.           ###
###################################################################################


# returns a random point in S
randPointDefault(S::CSpace{Float64}) = ( S.lowerBounds + ( rand(1, S.d) .* S.width ))

# returns a random node from S
randNodeDefault(S::CSpace{Float64}) = RRTNode{Float64}(S.lowerBounds + ( rand(1, S.d) .* S.width ))

# returns a random node from S, or the goal with probability pGoal
function randNodeOrGoal(S::CSpace)
  if rand() > S.pGoal
    return randNodeDefault(S)
  end
  return S.goalNode
end

# returns a random node from S, but when itsSamplePoint == 0 it returns
# itsSamplePoint instead
function randNodeIts(S::CSpace{T}) where {T}
  if S.itsUntilSample == 0
    S.itsUntilSample -= 1
    return RRTNode{T}(S.itsSamplePoint)
  end
  S.itsUntilSample -= 1
  return randNodeOrGoal(S)
end

# returns a random node from S, but when waitTime has passed it returns
# timeSamplePoint instead
function randNodeTime(S::CSpace{T}) where {T}
  if S.waitTime != Inf && S.timeElapsed >= S.waitTime
    S.waitTime = Inf
    return RRTNode{T}(S.timeSamplePoint)
  end
  return randNodeOrGoal(S)
end


# returns a random node from S, but when waitTime has passed it returns
# timeSamplePoint instead, also sets the first obstacle to unused
function randNodeTimeWithObstacleRemove(S::CSpace{T}) where {T}
  if S.waitTime != Inf && S.elapsedTime >= S.waitTime
    S.waitTime = Inf
    S.obstaclesToRemove.obstacleUnused = true
    S.obstaclesToRemove.startTime = -Inf
    S.obstaclesToRemove.lifeSpan = 0.0
    S.obstaclesToRemove.senseableObstacle = false
    S.obstaclesToRemove.obstacleUnusedAfterSense = true

    #println(S.timeSamplePoint)
    println("xxxxxxxxxxxxxx sample xxxxxxxxxxxxxxxx")
    a = copy(S.timeSamplePoint)
    return RRTNode{T}(a)
  end

  return randNodeOrGoal(S)
end


# returns a random node from S, but when waitTime has passed it returns
# timeSamplePoint instead, also sets the first obstacle to unused
function randNodeItsWithObstacleRemove(S::CSpace{T}) where {T}

  S.itsUntilSample -= 1
  if S.itsUntilSample == 0
    S.waitTime = Inf
    S.obstaclesToRemove.obstacleUnused = true
    S.obstaclesToRemove.startTime = -Inf
    S.obstaclesToRemove.lifeSpan = 0.0
    S.obstaclesToRemove.senseableObstacle = false
    S.obstaclesToRemove.obstacleUnusedAfterSense = true

    #println(S.timeSamplePoint)
    println("xxxxxxxxxxxxxx sample xxxxxxxxxxxxxxxx")
    a = copy(S.timeSamplePoint)
    return RRTNode{T}(a)
  end

  return randNodeOrGoal(S)
end


# this returns a random node unless there are points in the sample
# stack, in which case it returns the first one of those
function randNodeOrFromStack(S::CSpace{T}) where {T}
  if S.sampleStack.length > 0
    return RRTNode{T}(JlistPop(S.sampleStack))
  else
    return randNodeOrGoal(S)
  end
end

# this returns a random node where the time dimension is drawn uniformly at
# random from {the min time the robot could reach the point in an obstacle-less
# environment travelinga at max speed} and {current move time}---unless there are
# points in the sample stack, in which case it returns the first one of those.
function randNodeInTimeOrFromStack(S::CSpace{T}) where {T}
  if S.sampleStack.length > 0
    return RRTNode{T}(JlistPop(S.sampleStack))
  else
    newNode = randNodeOrGoal(S)
    if newNode == S.goalNode
      return newNode
    end

    minTimeToReachNode = S.start[3] + sqrt((newNode.position[1] - S.root.position[1])^2 + (newNode.position[2] - S.root.position[2])^2)/S.robotVelocity

    # if point is too soon vs robot's available speed
    # or it is in the "past" and the robot is moving
    if (newNode.position[3] < minTimeToReachNode ||
        (newNode.position[3] > S.moveGoal.position[3] && S.moveGoal != S.goalNode))

      # resample time in ok range
      newNode.position[3] = minTimeToReachNode + rand()*(S.moveGoal.position[3]-minTimeToReachNode)
    end

    return newNode
  end
end


# returns a random point from within the obstacle
function randomSampleObs(S::CSpace, KD::TKD, ob::Obstacle) where {TKD}
  # calculation of the number of samples to use could be made more accurate
  # by also multiplying by the ratio of free vs total random samples we
  # have observed up to this point over the total space)

  # sample from a hypercube that contains the entire obstacle, and
  # reject if the the point is not within the obstacle


  # calculate volume of the bounding hypercube
  obHypervolumeBound::Float64 = 0.0
  if !S.spaceHasTime && !S.spaceHasTheta
    # euclidian space, no time dimension

    obHypervolumeBound = (2.0*ob.radius)^S.d

    if S.hypervolume == 0.0 # need to calculate hypervolume of space
      S.hypervolume = prod(S.width)
    end
  elseif !S.spaceHasTime && S.spaceHasTheta
    # dubins car, no time dimension

    obHypervolumeBound = (2.0*ob.radius)^2 # * S.width[4]

    if S.hypervolume == 0.0 # need to calculate hypervolume of space
      S.hypervolume = prod(S.width[1:2])   # * S.width[4]
    end
  else
    error("not coded yet")
  end

  numObsSamples = KD.treeSize * obHypervolumeBound/S.hypervolume + 1.0

  for smp = 1.0:1.0:numObsSamples
    newPoint = rand(1, S.d)
    newPoint[1:2] = ob.position[1:2] .- ob.radius + newPoint[1:2] * ob.radius*2.0

    if quickCheck2D(ob, newPoint)
      JlistPush(S.sampleStack, newPoint)
    end
  end
end


# adds the obstacle to the C space
addObsToCSpace(C::CSpace, Ob::Obstacle) = listPush(C.obstacles, Ob)


# returns a random obstacle of type 1 with radius, makes sure that it is not
# on the start or goal
function randObstacle1(S::CSpace, radius::Float64)
  obs = Obstacle(1, randPointDefault(S), radius)
  while inCollision(obs, S.start) || inCollision(obs, S.goal)
    obs = Obstacle(1, randPointDefault(S), radius)
  end
  return obs
end

# returns a random obstacle of type 2 with span, makes sure that it is not
# on the start or goal
function randObstacle2(S::CSpace, span::Array{Float64})
  obs = Obstacle(2, randPointDefault(S), span)
  (Ba, junk) = inCollision(obs, S.start)
  (Bb, junk) = inCollision(obs, S.goal)
  while Ba || Bb
    obs = Obstacle(2, randPointDefault(S), span)
  end
  return obs
end



# reads a polygon from a file
# stores each obstacle obsMult times to simulate more complex environments
# whenver this is used for non-experimental pourposes obsMult should be 1
function readObstaclesFromfile(S::CSpace, filename, obsMult)
  a = open(filename, "r")

  # get the number of polygons
  P = parse(Int, readline(a))

  for p = 1:P
    # get the numer of point in this polygon
    N = parse(Int, readline(a))

    polygon = Array{Float64}(undef, N, 2)
    for n = 1:N
      polygon[n,:] = str2array(readline(a))
    end

    for i = 1:obsMult
      addObsToCSpace(S, Obstacle(3, polygon))
    end
  end

  close(a)
end


# reads a polygon from a file that is dynamic (start time and lifespan)
# stores each obstacle obsMult times to simulate more complex environments
# whenver this is used for non-experimental pourposes obsMult should be 1
function readDynamicObstaclesFromfile(S::CSpace, filename, obsMult)
  a = open(filename, "r")

  # get the number of polygons
  P = parse(Int, readline(a))

  for p = 1:P
    # get the numer of point in this polygon
    N = parse(Int, readline(a))

    polygon = Array{Float64}(undef, N, 2)
    for n = 1:N
      polygon[n,:] = str2array(readline(a))
    end

    startTimeAndLifeSpan = str2array(readline(a))

    for i = 1:obsMult
      ob = Obstacle(3, polygon)
      ob.startTime = startTimeAndLifeSpan[1]
      ob.lifeSpan = startTimeAndLifeSpan[2]
      ob.obstacleUnused = true

      addObsToCSpace(S, ob)
    end
  end

  close(a)
end


# reads a polygon from a file that is discoverable (always obs, or vanishes
# or appears once the robot gets within sensor range).
# stores each obstacle obsMult times to simulate more complex environments
# whenver this is used for non-experimental pourposes obsMult should be 1
function readDiscoverablecObstaclesFromfile(S::CSpace, filename, obsMult)
  a = open(filename, "r")

  # get the number of polygons
  P = parse(Int, readline(a))

  for p = 1:P
    # get the numer of point in this polygon
    N = parse(Int, readline(a))

    polygon = Array{Float64}(undef, N, 2)
    for n = 1:N
      polygon[n,:] = str2array(readline(a))
    end

    obsBehavourType = parse(Int, readline(a))

    for i = 1:obsMult
      ob = Obstacle(3, polygon)

      if obsBehavourType == 0 # normal obstacle
        ob.senseableObstacle = false
        ob.obstacleUnusedAfterSense = false
        ob.obstacleUnused = false
        println("normal")
      elseif obsBehavourType == -1 # vanishes when robot is within range
        ob.senseableObstacle = true
        ob.obstacleUnusedAfterSense = true
        ob.obstacleUnused = false
        println("vanishing")
      elseif obsBehavourType == 1  # appears when robot is within range
        ob.senseableObstacle = true
        ob.obstacleUnusedAfterSense = false
        ob.obstacleUnused = true
        println("appearing")
      else
        error("unknown behavoiur type")
      end

      addObsToCSpace(S, ob)
    end


  end

  close(a)
end



# reads a polygon from a file, which represents a prismatic obstacle vs
# all dimensions higher than the first two, where span in higher dimensions
# is stored in the fields PrismSpanMin::Array{Float64} and
# PrismSpanMax::Array{Float64}
# this was done because I am lazy, I needed a high dimensional experiment to
# evaluate convergence in a high dimensional bottleneck and did't feel like
# coding up a full-blown D-dimensioal obstacle checking thing
function readPrismaticObstaclesFromfile(S::CSpace, filename, totalDims)
  a = open(filename, "r")

  # get the number of polygons
  P = parse(Int, readline(a))

  for p = 1:P
    # get the numer of point in this polygon
    N = parse(Int, readline(a))

    polygon = Array{Float64}(undef, N, 2)
    PrismSpanMin = Array{Float64}(undef, 1, totalDims)
    PrismSpanMax = Array{Float64}(undef, 1, totalDims)
    for n = 1:N
      data = str2array(readline(a))
      polygon[n,:] = data[1:2]

      PrismSpanMin[3:totalDims] = data[3:2:length(data)]
      PrismSpanMax[3:totalDims] = data[4:2:length(data)]
    end
    addObsToCSpace(S, Obstacle(5, polygon, PrismSpanMin, PrismSpanMax))
  end

  close(a)
end


# reads a directional polygon from a file
# stores each obstacle obsMult times to simulate more complex environments
# whenver this is used for non-experimental pourposes obsMult should be 1
function readDirectionalObstaclesFromfile(S::CSpace, filename, obsMult)
  a = open(filename, "r")

  # get the number of polygons
  P = parse(Int, readline(a))

  for p = 1:P
    # get the direction of this obstacle
    D = readline(a)

    # get the numer of point in this polygon
    N = parse(Int, readline(a))

    polygon = Array{Float64}(undef, N, 2)
    for n = 1:N
      polygon[n,:] = str2array(readline(a))
    end

    for i = 1:obsMult
      if D[1] != 'X'
        addObsToCSpace(S, Obstacle(4, polygon, D[1]))
      else
        addObsToCSpace(S, Obstacle(3, polygon))
      end
    end
  end

  close(a)
end



# reads time obstacles (i.e., a polygon that moves vs. time) from the file
# movement is defiend by a 2d path
function readTimeObstaclesFromfile(S::CSpace, filename, obsMult)
  a = open(filename, "r")

  # get the number of polygons
  P = parse(Int, readline(a))

  for p = 1:P
    # get the numer of point in this polygon
    N = parse(Int, readline(a))

    polygon = Array{Float64}(undef, N, 2)
    for n = 1:N
      polygon[n,:] = str2array(readline(a))
    end

    obsSpeed = parse(Float64, readline(a))
    M = parse(Int, readline(a))
    movepath = Array{Float64}(undef, M, 3)
    for m = 1:M
      movepath[m,:] = str2array(readline(a))
    end


    for i = 1:obsMult
      ob = Obstacle(6, polygon) ######
      ob.senseableObstacle = false
      ob.obstacleUnusedAfterSense = false
      ob.obstacleUnused = false
      ob.velocity = obsSpeed
      ob.path = movepath
      ob.originalPolygon = polygon
      ob.nextDirectionChangeTime = -Inf


      addObsToCSpace(S, ob)
    end
  end

  close(a)
end


# reads time obstacles (i.e., a polygon that moves vs. time) from the file
# movement is defiend by a 2d path but that path is unknown to the robot,
# which must assume that the obstacle will continue moving along its current
# vector at any given time
function readDynamicTimeObstaclesFromfile(S::CSpace, filename, obsMult)
  a = open(filename, "r")

  # get the number of polygons
  P = parse(Int, readline(a))

  for p = 1:P
    # get the numer of point in this polygon
    N = parse(Int, readline(a))

    polygon = Array{Float64}(undef, N, 2)
    for n = 1:N
      polygon[n,:] = str2array(readline(a))
    end

    obsSpeed = parse(Float64, readline(a))
    M = parse(Int, readline(a))
    movepath = Array{Float64}(undef, M, 3)
    for m = 1:M
      movepath[m,:] = str2array(readline(a))
    end


    for i = 1:obsMult
      ob = Obstacle(7, polygon) ######
      ob.senseableObstacle = false
      ob.obstacleUnusedAfterSense = false
      ob.obstacleUnused = false
      ob.velocity = obsSpeed
      ob.unknownPath = movepath
      ob.nextDirectionChangeInd = size(movepath, 1)
      ob.nextDirectionChangeTime = movepath[ob.nextDirectionChangeInd,3]
      ob.originalPolygon = polygon
      ob.lastDirectionChangeTime = Inf
      changeObstacleDirection(S, ob, S.goal[3])

      addObsToCSpace(S, ob)
    end
  end

  close(a)
end


# this function is used to move obstacles around or remove them
function updateObstacles(C::CSpace)
  obstacleListNode = C.obstacles.front
  for i = 1:C.obstacles.length
    thisObstacle = obstacleListNode.data
    decreaseLife(thisObstacle)
    obstacleListNode = obstacleListNode.child #iterate
  end
end

# makes sure that thisNode can, in fact, reach all neighbors
# (used for error checking)
function checkNeighborsForEdgeProblems(S::TS, thisNode::T) where {T, TS}

  if thisNode.rrtParentUsed
    if explicitEdgeCheck(S, thisNode, thisNode.rrtParentEdge.endNode)
      error("problems $(thisNode.inPriorityQueue) $(thisNode.inOSQueue)")
    end
  end
  listItem = thisNode.rrtNeighborsOut.front
  while listItem != listItem.child
    neighborNode = listItem.data

    if neighborNode.rrtParentUsed && explicitEdgeCheck(S, neighborNode, neighborNode.rrtParentEdge.endNode)
      error("problems $(thisNode.inPriorityQueue) $(thisNode.inOSQueue)")
    end

    listItem = listItem.child
  end
end



##############################  geometric functions ###############################
### Slightly more more complex geometric functions.                             ###
###################################################################################

# return true if the point is in the polygon (open set of it, at least)
# point is a 1X2 point, polygon is pX2 where every row represents a point
# and rows i and i+1 are assumed to define an edge of the polygon (also
# top and bottom rows are connected) Note that the polygon DOES NOT have to be
# convex, but it should probably be simple (or unexpected behaviour may occour)
function pointInPolygon(point::Array{Float64}, polygon::Array{Float64})
  # this uses the MacMartin version of the crossings test

  P = size(polygon, 1)
  if P < 2
    return false
  end

  numCrossings = 0

  # start with the last point vs the first point
  startPoint = polygon[P, 1:2]
  for i = 1:P
    endPoint = polygon[i, 1:2]

    # check if this edge crosses the y value of point
    if (startPoint[2] > point[2] && endPoint[2] < point[2]) ||
       (startPoint[2] < point[2] && endPoint[2] > point[2])

       # it does, so now we need to see if the ray from point -> [inf, 0]
       # intersects the line.
       if startPoint[1] > point[1] && endPoint[1] > point[1]
         # definatly yes if both x coordiantes are right of the point
         numCrossings += 1
       elseif startPoint[1] < point[1] && endPoint[1] < point[1]
         # definatly no if both x coordiantes are left of the point
       else
         # otherwise we need to do the "expensive" calculation


         T = 2*max(startPoint[1], endPoint[1])

         x = (-((startPoint[1]*endPoint[2]-startPoint[2]*endPoint[1])*(point[1]-T))+((startPoint[1]-endPoint[1])*(point[1]*point[2]-point[2]*T)))/((startPoint[2]-endPoint[2])*(point[1]-T));

         if x > point[1]
           numCrossings += 1
         end
       end
    end
    startPoint = copy(endPoint)
  end

  # now check number of crossings (odd means point inside polygon)
  if mod(numCrossings,2) == 0 # even number of corssings
    return false
  end
  return true
end

# returns the min distance squared between the point and the segment
# [startPoint, endPoint] assumes a 2d space
function distanceSqrdPointToSegment(point::Array{Float64}, startPoint::Array{Float64}, endPoint::Array{Float64})
  vx = point[1]-startPoint[1];
  vy = point[2]-startPoint[2];
  ux = endPoint[1]-startPoint[1];
  uy = endPoint[2]-startPoint[2] ;
  determinate = vx*ux + vy*uy;

  if determinate <= 0
    #thisClosestPoint = startPoint;
    return vx*vx + vy*vy;
  else
    len = ux*ux + uy*uy;
    if determinate >= len
      #thisClosestPoint = endPoint;
      return (endPoint[1]-point[1])^2 + (endPoint[2]-point[2])^2;
    else
      #ex = ux / sqrt(len);
      #ey = uy / sqrt(len);
      #f = ex * vx + ey * vy;
      #thisClosestPoint = [startPoint[1] + f * ex,  startPoint[2] + f * ey];
      return (ux*vy-uy*vx)^2 / len;
    end
  end
end

# this returns the distance of closest point on the bouandary of polygon to point
# assumes a 2D space
function distToPolygonSqrd(point::Array{Float64}, polygon::Array{Float64})
  minDistanceSqrd = convert(Float64,Inf)

  # start with the last point vs the first point
  P = size(polygon, 1)
  startPoint = polygon[P, 1:2]

  for i = 1:P
    endPoint = polygon[i, 1:2]

    thisDistanceSqrd = distanceSqrdPointToSegment(point, startPoint, endPoint)

    if thisDistanceSqrd < minDistanceSqrd
      minDistanceSqrd = thisDistanceSqrd
    end
    startPoint = copy(endPoint)
  end

  return minDistanceSqrd
end


# this returns the distance of closest point on the bouandary of polygon to point
# assumes a 2D space
# note: hacky, see notes on high d bottleneck experiment elseshere in code
function distToPolygonPrismSqrd(point::Array{Float64}, polygon::Array{Float64}, thisObstacle::Obstacle)
  minDistanceSqrd = convert(Float64,Inf)

  # start with the last point vs the first point
  P = size(polygon, 1)
  startPoint = polygon[P, 1:2]

  for i = 1:P
    endPoint = polygon[i, 1:2]

    # first check that are actually in prism with respect to other dims
    # this is hacky (brittle)
    # see notes on stuff coded for high d bottleneck experiment elsewhere in code
    if (! segInPrism(thisObstacle, point, point))
      # hacky, only works for what I coded it for (binary detection)
      continue
    end
    thisDistanceSqrd = distanceSqrdPointToSegment(point, startPoint, endPoint)

    if thisDistanceSqrd < minDistanceSqrd
      minDistanceSqrd = thisDistanceSqrd
    end
    startPoint = copy(endPoint)
  end

  return minDistanceSqrd
end



# all input args represent points, this returns the minimum distance between line
# segments [PA PB] and [QA QB] and assumes 2D space
function segmentDistSqrd(PA::Array{Float64}, PB::Array{Float64}, QA::Array{Float64}, QB::Array{Float64})
  # check of the points are definatly not in collision by seeing if both
  # points of Q are on the same side of the line containing P and vise verse

  possibleIntersect = true

  # first check if P is close to vertical
  if abs(PB[1] - PA[1]) < .000001
    if (QA[1] >= PA[1] && QB[1] >= PA[1]) || (QA[1] <= PA[1] && QB[1] <= PA[1])
      # Q is on one side of P
      possibleIntersect = false
    end
  else
    # P is not close to vertical
    m = (PB[2] - PA[2])/(PB[1] - PA[1])

    # equation for points on P: y = m*(x - PA[1]) + PA[2]
    diffA = (m*(QA[1] - PA[1]) + PA[2]) - QA[2]
    diffB = (m*(QB[1] - PA[1]) + PA[2]) - QB[2]
    if (diffA > 0.0 && diffB > 0.0) || (diffA < 0.0 && diffB < 0.0)
      # Q is either fully above or below the line containing P
      possibleIntersect = false
    end
  end

  if possibleIntersect
    # first check if Q is close to vertical
    if abs(QB[1] - QA[1]) < .000001
      if (PA[1] >= QA[1] && PB[1] >= QA[1]) || (PA[1] <= QA[1] && PB[1] <= QA[1])
        # P is on one side of Q
        possibleIntersect = false
      end
    else
      # Q is not close to vertical
      m = (QB[2] - QA[2])/(QB[1] - QA[1])

      # equation for points on Q: y = m*(x - QA[1]) + QA[2]
      diffA = (m*(PA[1] - QA[1]) + QA[2]) - PA[2]
      diffB = (m*(PB[1] - QA[1]) + QA[2]) - PB[2]
      if (diffA > 0.0 && diffB > 0.0) || (diffA < 0.0 && diffB < 0.0)
        # P is either fully above or below the line containing Q
        possibleIntersect = false
      end
    end
  end

  if possibleIntersect
    # then there is an intersection for sure
    return 0.0
  end

  # when the lines do not intersect,
  # in 2D the min distance must be between one segment's end point
  # and the other segment (assuming lines are not parallel
  return min(distanceSqrdPointToSegment(PA, QA, QB),
                  distanceSqrdPointToSegment(PB, QA, QB),
                  distanceSqrdPointToSegment(QA, PA, PB),
                  distanceSqrdPointToSegment(QB, PA, PB))
end


# returns the index of the first time coord (3rd dimension) smaller than the time
# NOTE THIS SHOULD BE REPLACE WITH BINARY SEARCH WHEN I GET A CHANCE
function findIndexBeforeTime(path::Array{Float64,2}, timeToFind::Float64)

  if size(path, 1) < 1
    return -1
  end

  i = 0
  while i+1 <= size(path, 1) && path[i+1,3] < timeToFind
    i += 1
  end
  return i
end

# finds the transform of polygon to the approperiate position at
# the time of the point, based on the former's path through time
# assumes obstacle type 6
function findTransformObsToTimeOfPoint(thisObstacle::Obstacle, point::Array{Float64})
  # start by finding the indicies of the path edge containing time
  indBefore = findIndexBeforeTime(thisObstacle.path, point[1,3])

  # if before or after times we know for obstace, assume that it is
  # stationary at the relivant path end before/after this time
  if indBefore < 1
    dx = thisObstacle.path[1,1]
    dy = thisObstacle.path[1,2]
    return (dx, dy)
  elseif indBefore == size(thisObstacle.path,1)
    dx = thisObstacle.path[end,1]
    dy = thisObstacle.path[end,2]
    return (dx, dy)
  end

  indAfter = indBefore + 1
  perportionAlongEdge = ((point[3] - thisObstacle.path[indBefore,3])/
                         (thisObstacle.path[indAfter,3] - thisObstacle.path[indBefore,3]))

  dx = thisObstacle.path[indBefore,1] + perportionAlongEdge*(thisObstacle.path[indAfter,1] - thisObstacle.path[indBefore,1])
  dy = thisObstacle.path[indBefore,2] + perportionAlongEdge*(thisObstacle.path[indAfter,2] - thisObstacle.path[indBefore,2])

  return (dx, dy)
end



########################## collision checking functions ###########################
### Collision checking, etc. This includes certificate stuff that is currently  ###
### unused, but could be added later with little difficulty.                    ###
###################################################################################

# checks if the 2D point is inside the obstacle
# (does not account for robot radius)
function quickCheck2D(thisObstacle::Obstacle, point::Array{Float64})
  pointInsideObstacle = false

  if thisObstacle.obstacleUnused || thisObstacle.lifeSpan <= 0
    return false
  end

  # do a quick check based on a lower bound of the distance to the obstacle
  if (1 <= thisObstacle.kind <= 5) && Wdist(thisObstacle.position, point) > thisObstacle.radius
    return false
  end

  if thisObstacle.kind == 1
    # lower bound is actual distance in this case
    return true
  elseif thisObstacle.kind == 2 || thisObstacle.kind == 4 # added latter recently, might not supposed to be here
    error("need to impliment this")

  elseif thisObstacle.kind == 3

    # if that doesn't collide then need to check vs the polygon itself
    if pointInPolygon(point, thisObstacle.polygon)
      return true
    end
  elseif thisObstacle.kind == 5
    # if that doesn't collide then need to check vs the polygon itself
    if pointInPolygon(point[1:2], thisObstacle.polygon)
      if segInPrism(thisObstacle, point, point)
        return true
      end
    end
  elseif thisObstacle.kind == 6 || thisObstacle.kind == 7

    # first need to transform position and polygon to the approperiate position at
    # the time of the point, based on the former's path through time
    (dx, dy) = findTransformObsToTimeOfPoint(thisObstacle, point)

    # do a quick check based on a lower bound of the distance to the obstacle
    if Wdist(thisObstacle.position + [dx dy] , point) > thisObstacle.radius
      return false
    end

    # transform polygon and then do a normal check vs it
    thisObstacle.polygon[:,1] = thisObstacle.originalPolygon[:,1] .+ dx
    thisObstacle.polygon[:,2] = thisObstacle.originalPolygon[:,2] .+ dy
    if pointInPolygon(point, thisObstacle.polygon)
      return true
    end

  end

  return false
end

# checks if the point is inside any obstacles
# (does not account for robot radii)
function quickCheck(C::CSpace{Float64}, point::Array{Float64})

  R = div(C.d,2)
  retCert = ones(R)*convert(Float64,Inf)

  obstacleListNode = C.obstacles.front
  for i = 1:C.obstacles.length

    if quickCheck2D(obstacleListNode.data, point)
      return true
    end

    obstacleListNode = obstacleListNode.child #iterate
  end


  return false
end


# checks if the node is in collision with an obstacle
quickCheck(C::CSpace{T}, node::RRTNode{T}) where {T} = quickCheck(C, node.position)



# checks if the 2D point is in collision with the obstacle
# if -not- it retruns the distance to the closest point on that obstacle
# minDist is the closest point found so far (init to inf before calling)
# robot radius is the radius of this robot
function explicitPointCheck2D(thisObstacle::Obstacle, point::Array{Float64}, minDist::Float64, robotRadius::Float64)

  thisDist::Float64 = convert(Float64,Inf)

  if thisObstacle.obstacleUnused || thisObstacle.lifeSpan <= 0
    return (false, minDist)
  end

  if (1 <= thisObstacle.kind <= 5)
    # do a quick check to see if any points on the obstacle might be closer
    # to point than minDist based on the ball around the obstacle

    # calculate distance from robot boundary to obstacle center
    thisDist = Wdist(thisObstacle.position, point) - robotRadius
    if thisDist - thisObstacle.radius > minDist
      return (false, minDist)
    end
  end

  if thisObstacle.kind == 1
    # ball is the actual obstacle, so we have a new minimum
    thisDist = thisDist - thisObstacle.radius
    if thisDist < 0.0
      return (true, 0.0)
    end

  elseif thisObstacle.kind == 2
    error("need to impliment this")

  elseif thisObstacle.kind == 3

    if pointInPolygon(point, thisObstacle.polygon)
       return (true, 0.0)
    end
    thisDist = sqrt(distToPolygonSqrd(point, thisObstacle.polygon)) - robotRadius

    if thisDist < 0.0
      return (true, 0.0)
    end

  elseif thisObstacle.kind == 5

    if pointInPolygon(point[1:2], thisObstacle.polygon)
      if segInPrism(thisObstacle, point, point)
        return (true, 0.0)
      end
    end
    thisDist = sqrt(distToPolygonPrismSqrd(point, thisObstacle.polygon, thisObstacle)) - robotRadius

    if thisDist < 0.0
      return (true, 0.0)
    end
  elseif thisObstacle.kind == 6 || thisObstacle.kind == 7

    # first need to transform position and polygon to the approperiate position at
    # the time of the point, based on the former's path through time
    (dx, dy) = findTransformObsToTimeOfPoint(thisObstacle, point)

    # do a quick check to see if any points on the obstacle might be closer
    # to point than minDist based on the ball around the obstacle

    # calculate distance from robot boundary to obstacle center
    thisDist = Wdist(thisObstacle.position + [dx dy], point) - robotRadius
    if thisDist - thisObstacle.radius > minDist
      return (false, minDist)
    end

    # transform polygon and then do the rest of a normal check
    thisObstacle.polygon[:,1] = thisObstacle.originalPolygon[:,1] .+ dx
    thisObstacle.polygon[:,2] = thisObstacle.originalPolygon[:,2] .+ dy

    if pointInPolygon(point, thisObstacle.polygon)
       return (true, 0.0)
    end
    thisDist = sqrt(distToPolygonSqrd(point, thisObstacle.polygon)) - robotRadius
    if thisDist < 0.0
      return (true, 0.0)
    end

  else
    error("need to impliment this")
  end

  return (false, min(minDist, thisDist))
end


# checks if the point is in collision with any obstacles
# (DOES account for robot radii)
# returns a bool and also, IF the point is NOT in collision, the values that
# should be used for the certificate (dist of closest obstacle per robot)
function explicitPointCheck(C::CSpace{T}, point::Array{Float64}) where {T}

  # if we are ignoring obstacles
  if C.inWarmupTime
    return (false, Inf)
  end

  # first do a quick check to see if the point can be determined
  # in collision with minimal work. NOTE a quick check is NOT an implicit check
  if quickCheck(C, point)
    return (true, 0.0)
  end


  # the point in not inside any of the obstacles, but it may still be in
  # collision due to the robot radius
  R = div(C.d,2)

  retCert = convert(Float64,Inf)

  obstacleListNode = C.obstacles.front
  for i = 1:C.obstacles.length

    (thisRetVal, thisCert) = explicitPointCheck2D(obstacleListNode.data, point, retCert, C.robotRadius)

    if thisRetVal
      return (true, 0.0)
    end
    if thisCert < retCert
      retCert = thisCert
    end
    obstacleListNode = obstacleListNode.child #iterate
  end


  return (false, retCert)
end


# checks if the node is in collision with an obstacle
# (DOES account for robot radii)
explicitNodeCheck(C::CSpace{T}, node::RRTNode{T}) where {T} = explicitPointCheck(C, node.position)


# returns true if some part of the segment [startPoint endPoint]
# is inside the x = X[3:D] projection (highest d-2 dimensions)
# this is coded sloppy, however, it is not expected to be called much,
# It was created only for use with the high d bottlenect experiment
function segInPrism(thisObstacle::Obstacle, startPoint::Array{Float64}, endPoint::Array{Float64})

  # find thinest dimesnion
  delta = Inf
  for d = 3:length(thisObstacle.PrismSpanMax)
    if delta > thisObstacle.PrismSpanMax[d] - thisObstacle.PrismSpanMin[d]
      delta = thisObstacle.PrismSpanMax[d] - thisObstacle.PrismSpanMin[d]
    end
  end
  delta /= 20.0
  steps = dist(startPoint, endPoint)/delta

  if startPoint == endPoint
    delta = 1
    steps = 1
  end

  for s = 0.0:1.0:steps
    testPose = startPoint + (endPoint - startPoint)*s/steps

    allIn = true
    for d = 3:length(thisObstacle.PrismSpanMax)

      if (testPose[d] < thisObstacle.PrismSpanMin[d] ||
          testPose[d] > thisObstacle.PrismSpanMax[d])

          allIn = false
          break
      end
    end

    if allIn
      return true
    end
  end
  return false
end


# checks if the edge between the points is in collision with the obstacle
# (the point is closer than robotRadius to the edge)
function explicitEdgeCheck2D(thisObstacle::Obstacle, startPoint::Array{Float64}, endPoint::Array{Float64}, robotRadius::Float64)

  if thisObstacle.obstacleUnused || thisObstacle.lifeSpan <= 0
    return false
  end


  # do a quick check to see if any points on the obstacle might be closer
  # to the edge than robot radius
  if 1 <= thisObstacle.kind <= 5

    # calculate distance squared from center of the obstacle to the edge
    # (projected into first two dims)
    distSqrd = distanceSqrdPointToSegment(thisObstacle.position, startPoint[1:2], endPoint[1:2])
    if distSqrd > (robotRadius + thisObstacle.radius)^2
      return false
    end
  end

  if thisObstacle.kind == 1
    # ball is the actual obstacle, so in collision
    return true
  elseif thisObstacle.kind == 2
    error("need to impliment this")

  elseif thisObstacle.kind == 3 || thisObstacle.kind == 5
    # need to check vs all edges in the polygon

    P = size(thisObstacle.polygon, 1)
    if P < 2
      return false
    end

    # start with the last point vs the first point
    A = thisObstacle.polygon[P, 1:2]
    for i = 1:P
      B = thisObstacle.polygon[i, 1:2]

      if segmentDistSqrd(startPoint[1:2], endPoint[1:2], A, B) < robotRadius^2
        # there is a colision (with the 2d projection of the obstacle)

        if thisObstacle.kind == 5
          # if this is a high-D prismatic obstacle,
          # then still need to make sure that the segment is
          # in the obstacle with respect to the other dims,
          # also calculate the thinest dimension of the ob

          if segInPrism(thisObstacle, startPoint, endPoint)
            return true
          end
        else # thisObstacle.kind == 3
          return true
        end
      end
      A = copy(B)
    end
  elseif thisObstacle.kind == 6 || thisObstacle.kind == 7
    # must check all edges of obstacle path that overlap with the robot edge in time

    # make life easier by always checking past to future
    if startPoint[3] < endPoint[3]
      earlyPoint = startPoint
      latePoint = endPoint
    else
      latePoint = startPoint
      earlyPoint = endPoint
    end


    firstObsInd = max(findIndexBeforeTime(thisObstacle.path, earlyPoint[3]), 1)
    lastObsInd = min(1+findIndexBeforeTime(thisObstacle.path, latePoint[3]), size(thisObstacle.path,1))

    if lastObsInd <= firstObsInd
      # obstacle does not overlap in time
      return false
    end



    for i_start = firstObsInd:(lastObsInd - 1)
      i_end = i_start + 1

      x_1 = earlyPoint[1]; # robot start x
      y_1 = earlyPoint[2]; # robot start y
      T_1 = earlyPoint[3]; # robot start time



      # calculate the time of minimum approach of the centers of obstacle and robot
      x_2 = thisObstacle.path[i_start,1] + thisObstacle.position[1]; # obs start x
      y_2 = thisObstacle.path[i_start,2] + thisObstacle.position[2]; # obs start y
      T_2 = thisObstacle.path[i_start,3];                            # obs start time




      # calculate intermediate quantities (parametric slopes)
      m_x1 = (latePoint[1] - x_1)/(latePoint[3]-T_1);
      m_y1 = (latePoint[2] - y_1)/(latePoint[3]-T_1);
      m_x2 = (thisObstacle.path[i_end,1] + thisObstacle.position[1] - x_2)/(thisObstacle.path[i_end,3]-T_2);
      m_y2 = (thisObstacle.path[i_end,2] + thisObstacle.position[2] - y_2)/(thisObstacle.path[i_end,3]-T_2);


      # solve for time of closest pass of centers:
      T_c = ((m_x1^2 * T_1 + m_x2 * (m_x2 * T_2 + x_1 - x_2) -
              m_x1 * (m_x2 * (T_1 + T_2) + x_1 - x_2) +
              (m_y1 - m_y2) * (m_y1 * T_1 - m_y2 * T_2 - y_1 + y_2)) /
             ((m_x1 - m_x2)^2 + (m_y1 - m_y2)^2));

      # now bound T_c by the allowable times of the robot and the obstacle
      if T_c < max(T_1, T_2)
        T_c = max(T_1, T_2);
      elseif T_c > min(latePoint[3], thisObstacle.path[i_end,3])
        T_c = min(latePoint[3], thisObstacle.path[i_end,3]);
      end

      # finally see if the distance between the robot and the obstacle at T_c
      # is close enough to cause a conflict
      r_x = m_x1*(T_c - T_1) + x_1;  # robot x at T_c
      r_y = m_y1*(T_c - T_1) + y_1;  # robot y at T_c
      o_x = m_x2*(T_c - T_2) + x_2;  # obstacle x at T_c
      o_y = m_y2*(T_c - T_2) + y_2;  # obstacle y at T_c

      if (r_x - o_x)^2 + (r_y - o_y)^2 < (thisObstacle.radius + robotRadius)^2
        # then there is a collision
        return true
      end
    end
  end
  return false
end




# this checks if the edge is in collision with any obstacles in the C-space
# recall that Edge is aliased to the particular type of edge being used
function explicitEdgeCheck(C::CSpace{T}, edge::Edge) where {T}

  # if we are ignoring obstacles
  if C.inWarmupTime
    return false
  end


  obstacleListNode = C.obstacles.front
  for i = 1:C.obstacles.length

    if explicitEdgeCheck(C, edge, obstacleListNode.data)
      return true
    end

    obstacleListNode = obstacleListNode.child #iterate
  end
  return false
end


# implicitly checks if this the edge between the point and the node
# is safe based on the (nearest neighbor) node, it returns flags as follows:
# 0 = safe
# 2 = unknown
# as well as the node that certifies node::T
function implicitEdgeCheck(point::Array{Float64}, node::T) where {T}
  certifyingNode = node
  if !node.hasCertificate
    certifyingNode = node.certifyingNode
  end

  if Wdist(point, certifyingNode.position) > certifyingNode.certificateValue
    # need explicit check for sure
    return (2, certifyingNode)
  end


  # implicitly safe
  return (0, certifyingNode)
end

# implicitly checks if the edge from nodeA to nodeB is safe based on
# nodeB's certificate, return values are described above
implicitEdgeCheck(nodeA::T, nodeB::T) where {T} = implicitEdgeCheck(nodeA.position, nodeB)


# implicitly checks if nodeA is safe based on nodeB's certificate, return values are described above, note that this assumes triangle inequality holds
implicitPointCheck(pointA::P, nodeB::T) where {P,T} = implicitEdgeCheck(pointA, nodeB)

# implicitly checks if pointA is safe based on nodeB's certificate, return values are described above, note that this assumes triangle inequality holds
implicitNodeCheck(nodeA::T, nodeB::T) where {T} = implicitEdgeCheck(nodeA.position, nodeB)

# implicitly checks if this the edge between the point and the node
# is safe based on the (nearest neighbor) node (its certificate is passed
# seperatly, it returns flags as follows:
# 0 = safe
# 2 = unknown
# as well as the node that certifies node::T
function implicitEdgeCheckCert(point::Array{Float64}, node::T, cert::C) where {T,C}
  certifyingNode = node

  # check per robot
  if Wdist(point, certifyingNode.position) > cert
    # need explicit check for sure
    return (2, certifyingNode)
  end

  # implicitly safe
  return (0, certifyingNode)
end

# implicitly checks if the edge from nodeA to nodeB is safe based on
# nodeB's certificate (passed seperatly), return values are described above
implicitEdgeCheckCert(nodeA::T, nodeB::T, nodeBCert::C) where {T,C} = implicitEdgeCheckCert(nodeA.position,
  nodeB, nodeBCert)

############################### RRT Functions #####################################
### Functions used for RRT. Some of these are also used in RRT*, RRT#, and RRTx.###
###################################################################################

# takes care of inserting a new node
function extend(S::TCS, KD::TKD, Q::rrtQueue, newNode::T,
  closestNode::T, delta::Float64, hyberBallRad::Float64, moveGoal::T) where {T, TCS, TKD}


  # first calculate the shortest trajectory (and its distance) that gets from
  # newNode to closestNode while obeying the constraints of the state space and the
  # dynamics of the robot
  thisEdge = newEdge(newNode, closestNode)
  calculateTrajectory(S, thisEdge)

  # figure out if we can link to the nearest node
  if !validMove(S, thisEdge) || explicitEdgeCheck(S, thisEdge)
    # we cannot link to nearest neighbor
    return
  end

  # otherwise we can link to the nearest neighbor
  newNode.rrtParentEdge = thisEdge
  newNode.rrtTreeCost = closestNode.rrtLMC + newNode.rrtParentEdge.dist
  newNode.rrtLMC = newNode.rrtTreeCost # only for compatability with visualization
  newNode.rrtParentUsed = true


  # insert the new node into the KDTree
  kdInsert(KD, newNode)
end


############################### RRT* Functions ####################################
### Functions used for RRT*. Some of these are also used in RRT#, and RRTx.     ###
###################################################################################

# this looks through the list to try to find the best parent for newNode
# if the list is empty then it uses closestNode instead. If a parent is found
# then newNode is linked to its parent. If saveAllEdges is true then it saves
# a copy of each edge in each node (this is used to avoud work duplication in
# RRT# and RRTx)
function findBestParent(S::TCS, newNode::T, nodeList::TL, closestNode::T,
  saveAllEdges::Bool) where {T, TCS, TL}

  # if the list is empty
  if nodeList.length == 0
    if S.goalNode != newNode
      JlistPush(nodeList, closestNode)
    end
  end

  # update LMC value based on nodes in the list
  newNode.rrtLMC = Inf
  newNode.rrtTreeCost = Inf
  newNode.rrtParentUsed = false

  # find best parent (and even if one exists)
  listItem = nodeList.front
  for i = 1:nodeList.length
    nearNode = listItem.data

    # first calculate the shortest trajectory (and its distance) that gets from
    # newNode to nearNode while obeying the constraints of the state space and the
    # dynamics of the robot
    thisEdge = newEdge(newNode, nearNode)
    calculateTrajectory(S, thisEdge)

    if saveAllEdges
      nearNode.tempEdge = thisEdge
    end

    # check for validity vs edge collisions vs obstacles, and vs the time-dynamics
    # of the robot and space

    if (explicitEdgeCheck(S, thisEdge) || !validMove(S, thisEdge))

      if saveAllEdges
        nearNode.tempEdge.dist = Inf
      end

      listItem = listItem.child         # iterate thorugh list
      continue
    end

    if newNode.rrtLMC > nearNode.rrtLMC + thisEdge.dist
      # found a potentially better parrent
      newNode.rrtLMC = nearNode.rrtLMC + thisEdge.dist
      newNode.rrtParentEdge = thisEdge
      newNode.rrtParentUsed = true
    end

    listItem = listItem.child   # iterate thorugh list
  end
end


# takes care of inserting a new node
function extend(S::TCS, KD::TKD, Q::rrtStarQueue, newNode::T,
  closestNode::T, delta::Float64, hyberBallRad::Float64, moveGoal::T) where {T, TCS, TKD}


  # find all nodes within the (shrinking) hyperball of (saturated) newNode
  nodeList = kdFindWithinRange(KD, hyberBallRad, newNode.position)

  # try to find and link to best parent
  findBestParent(S, newNode, nodeList, closestNode, false)
  if !newNode.rrtParentUsed
    emptyRangeList(nodeList)     # clean up
    return
  end

  newNode.rrtTreeCost = newNode.rrtLMC

  # insert the new node into the KDTree
  kdInsert(KD, newNode)

  # if this is inserted in an unhelpful part of the C-space then
  # we don't waste time rewiring (assumes triangle inequality, added
  # by MO, not technically part of RRT* but can only improve it)
  if newNode.rrtLMC > moveGoal.rrtLMC
    emptyRangeList(nodeList)     # clean up
    return
  end


  # now rewire neighbors that should use newNode as their parent
  listItem = nodeList.front
  for i = 1:nodeList.length
    nearNode = listItem.data

    # watch out for cycles
    if newNode.rrtParentEdge.endNode == nearNode
      listItem = listItem.child   # iterate thorugh list
      continue
    end

    # calculate the shortest trajectory (and its distance) that gets from
    # nearNode to newNode while obeying the constraints of the state space
    # and the dynamics of the robot
    thisEdge = newEdge(nearNode, newNode)
    calculateTrajectory(S, thisEdge)

    # rewire neighbors that would do better to use this node as their parent
    # unless they are in collision or impossible due to dynamics of robot/space
    if (nearNode.rrtLMC > newNode.rrtLMC + thisEdge.dist &&
        validMove(S, thisEdge) && !explicitEdgeCheck(S, thisEdge))

      # make this node the parent of the neighbor node
      nearNode.rrtParentEdge = thisEdge
      nearNode.rrtParentUsed = true

      # recalculate tree cost of neighbor
      nearNode.rrtTreeCost = nearNode.rrtLMC = newNode.rrtLMC + thisEdge.dist
    end
    listItem = listItem.child   # iterate thorugh list
  end
  emptyRangeList(nodeList)     # clean up
end


############################### RRT# Functions ####################################
### Functions used for RRT#. Some of these are also used in RRTx.               ###
### This includes the priority heap related key functions, etc.                 ###
###################################################################################

# returns the key value of the node
function keyQ(node::T) where {T}
  g_min = min(node.rrtTreeCost, node.rrtLMC)
  #return (g_min + node.rrtH, g_min) ################!!!!!!!!!!!!!!!!!!!!!
  return (g_min + 0.0, g_min)
end

# less than function for key values
function lessQ(a::T, b::T) where {T}
  (a_key_first, a_key_second) = keyQ(a)
  (b_key_first, b_key_second) = keyQ(b)
  if (a_key_first < b_key_first) || (a_key_first == b_key_first && a_key_second < b_key_second) || (a_key_first == b_key_first && a_key_second == b_key_second && a.isMoveGoal)
    return true
  end
  return false
end

# greater than function for key values
function greaterQ(a::T, b::T) where {T}
  (a_key_first, a_key_second) = keyQ(a)
  (b_key_first, b_key_second) = keyQ(b)
  if (a_key_first > b_key_first) || (a_key_first == b_key_first && a_key_second > b_key_second) || (a_key_first == b_key_first && a_key_second == b_key_second && b.isMoveGoal)
    return true
  end
  return false
end

# priority queue marker function (marks when a node is in the queue)
markQ(node::T) where {T} = (node.inPriorityQueue = true)

# priority queue unmarker function (un marks when a node is removed)
unmarkQ(node::T) where {T} = (node.inPriorityQueue = false)

# priority queue check marker function (checks if the node is marked)
markedQ(node::T) where {T} = node.inPriorityQueue

# sets the priority queue index to the value
setIndexQ(node::T, val::Int) where {T} = (node.priorityQueueIndex = val)

# set the priority queue index to the unused value
unsetIndexQ(node::T) where {T} = (node.priorityQueueIndex = -1)

# returns the priority queue index
getIndexQ(node::T) where {T} = node.priorityQueueIndex

# checks all nodes in the heap to see if there are edge problems
function checkHeapForEdgeProblems(Q::TQ) where {TQ}
  for i = 1:Q.Q.indexOfLast
    node = Q.Q.heapNode[i]
    checkNeighborsForEdgeProblems(Q.S, node)
  end
end

# resets the neighbor iterator
function resetNeighborIterator(It::RRTNodeNeighborIterator{T, TE}) where {T, TE}
  It.listFlag = 0
end

# rrt# based version (note, this is technically not used in the classical
# version of RRT#, but it can be used to allow propogate successors to work with
# RRT# so that obstacles can be added
# this returns the JlistNode containing the next neighbor of the node
# for which this iterator was created
function nextOutNeighbor(It::RRTNodeNeighborIterator{T, TE}, Q::rrtSharpQueue) where {T, TE}

  if It.listFlag == 0
    It.listItem = It.thisNode.rrtNeighborsOut.front
    It.listFlag = 1
  else
    It.listItem = It.listItem.child
  end
  if It.listItem == It.listItem.child
    # done with all neighbors
    return nothing
  end

  return It.listItem
end

# rrt# based version (note, this is technically not used in the classical
# version of RRT#, but it can be used to allow propogate successors to work with
# RRT# so that obstacles can be added
# this returns the JlistNode containing the next neighbor of the node
# for which this iterator was created
function nextInNeighbor(It::RRTNodeNeighborIterator{T, TE}, Q::rrtSharpQueue) where {T, TE}


  if It.listFlag == 0
    It.listItem = It.thisNode.rrtNeighborsIn.front
    It.listFlag = 1
  else
    It.listItem = It.listItem.child
  end
  if It.listItem == It.listItem.child
    # done with all neighbors
    return nothing
  end

  return It.listItem
end


# links an edge -from- node -to- newNeighbor, edge should already
# be populated correctly
function makeNeighborOf(newNeighbor::T, node::T, edge::Edge) where {T}

  JlistPush(node.rrtNeighborsOut, edge)
  edge.listItemInStartNode = node.rrtNeighborsOut.front

  JlistPush(newNeighbor.rrtNeighborsIn, edge)
  edge.listItemInEndNode = newNeighbor.rrtNeighborsIn.front

  return 0
end

# links an "initial" "out" edge -from- node -to- newNeighbor, edge should
# already be populated correctly. This is actually only used for RRTx
# but is included here because of its similarity to the function above
function makeInitialOutNeighborOf(newNeighbor::T, node::T, edge::Edge) where {T}
  JlistPush(node.InitialNeighborListOut, edge)
end

# links an "initial" "in" edge -"from"- node -"to"- newNeighbor (i.e.,
# the edge is only stored on the recieving node and not on the sending node.
# edge should already be populated correctly. This is actually only used for
# RRTx but is included here because of its similarity to the functions above
function makeInitialInNeighborOf(newNeighbor::T, node::T, edge::Edge) where {T}
  JlistPush(node.InitialNeighborListIn, edge)
end


# recalculates LMC based on neighbors that this node can reach
# Note the first argument is unused but necessary for the multiple dispatch
# that us used to differentiate between RRT* and RRT#
function recalculateLMC(Q::rrtSharpQueue, thisNode::T, root::T) where {T}

  if thisNode == root
    return
  end

  oldrrtLMC = thisNode.rrtLMC

  listItem = thisNode.rrtNeighborsOut.front
  for i = 1:thisNode.rrtNeighborsOut.length
    neighborEdge = listItem.data
    neighborNode = neighborEdge.endNode
    neighborDist = neighborEdge.dist

    if thisNode.rrtLMC > neighborNode.rrtLMC + neighborDist  && validMove(Q.S, neighborEdge)
      # found a potentially better parrent
      thisNode.rrtLMC = neighborNode.rrtLMC + neighborDist
      thisNode.rrtParentEdge = listItem.data
      thisNode.rrtParentUsed = true
    end

    listItem = listItem.child   # iterate through list
  end
end


# updates the priority queue (adds node if necessary, does not if not)
function updateQueue(Q::rrtSharpQueue, newNode::TN, root::TN, hyberBallRad) where {TN}

  recalculateLMC(Q, newNode, root)  # internally ignors root

  if markedQ(newNode)
    updateHeap(Q.Q, newNode)
    removeFromHeap(Q.Q, newNode)
  end
  if newNode.rrtTreeCost != newNode.rrtLMC
      addToHeap(Q.Q, newNode)
  end
end

# takes care of inserting a new node
#function extend{T, TCS, TKD}(S::TCS, KD::TKD, Q::rrtSharpQueue, newNode::T, closestNode::T, delta::Float64, hyberBallRad::Float64, moveGoal::T)

function extend(S::TCS, KD::TKD, Q::rrtSharpQueue, newNode::RRTNode{Float64},
  closestNode::RRTNode{Float64}, delta::Float64,
  hyberBallRad::Float64, moveGoal::RRTNode{Float64}) where {TCS, TKD}


  # find all nodes within the (shrinking) hyperball of (saturated) newNode
  nodeList = kdFindWithinRange(KD, hyberBallRad, newNode.position)

  # try to find and link to best parent, this also saves the edges from newNode
  # to the neighbors in the field "tempEdge" of the neighbors. This saves time in
  # the case that trajectory calculation is complicated.
  findBestParent(S, newNode, nodeList, closestNode, true)

  # if no parent was found then ignore this node
  if !newNode.rrtParentUsed
    emptyRangeList(nodeList)     # clean up
    return
  end

  # insert the new node into the KDTree
  kdInsert(KD, newNode)

  # first pass, make edges between newNode and all of its valid neighbors
  # Note that the edges have been stored in "tempEdge" field of the neighbors
  listItem = nodeList.front
  while listItem != listItem.child
    nearNode = listItem.data

    if nearNode.tempEdge.dist == Inf  # obstacle, edge invalid
      listItem = listItem.child      # iterate through list
      continue
    end

    # make nearNode a neighbor (can be reached from) newNode
    makeNeighborOf(nearNode, newNode, nearNode.tempEdge)

    listItem = listItem.child      # iterate through list
  end


  # second pass, make edges (if possible) between all valid nodes in D-ball
  # and newNode, also rewire neighbors that should use newNode as their parent
  listItem = nodeList.front
  while listItem != listItem.child
    nearNode = listItem.data

    # in the general case the trajectories along edges are not simply
    # the reverse of each other, therefore we need to calculate and
    # check the trajectory along the edge from nearNode to newNode.
    thisEdge = newEdge(nearNode, newNode)
    calculateTrajectory(S, thisEdge)

    if validMove(S, thisEdge) && !explicitEdgeCheck(S, thisEdge)
      makeNeighborOf(newNode, nearNode, thisEdge)
    else
      # edge cannot be created
      listItem = listItem.child   # iterate thorugh list
      continue
    end

    # rewire neighbors that would do better to use this node as their parent
    # unless they are not in the relivant portion of the space vs. moveGoal
    if (nearNode.rrtLMC > newNode.rrtLMC + thisEdge.dist &&
        newNode.rrtParentEdge.endNode != nearNode &&
        newNode.rrtLMC + thisEdge.dist < moveGoal.rrtLMC)

      # make this node the parent of the neighbor node
      nearNode.rrtParentEdge = thisEdge
      nearNode.rrtParentUsed = true

      # recalculate tree cost of neighbor
      nearNode.rrtLMC = newNode.rrtLMC + thisEdge.dist

      # insert the neighbor into priority queue if it is not consistant
      if nearNode.rrtLMC != nearNode.rrtTreeCost && markedQ(nearNode)
        updateHeap(Q.Q, nearNode)
      elseif nearNode.rrtTreeCost != nearNode.rrtLMC
        addToHeap(Q.Q, nearNode)
      elseif newNode.rrtTreeCost == newNode.rrtLMC  && markedQ(nearNode)
        updateHeap(Q.Q, nearNode)
        removeFromHeap(Q.Q, nearNode)
      end

    end
    listItem = listItem.child   # iterate thorugh list
  end
  emptyRangeList(nodeList)     # clean up

  # insert the node into the piority queue
  addToHeap(Q.Q, newNode)


end


# propogates cost information through the graph
function reduceInconsistency(Q::rrtSharpQueue, goalNode::T,
  robotRad::Float64, root::T, hyberBallRad::Float64) where {T}


  while Q.Q.indexOfLast > 0 && (lessQ(topHeap(Q.Q), goalNode) || goalNode.rrtLMC == Inf || goalNode.rrtTreeCost == Inf || markedQ(goalNode))



#Q.numReduces += 1;

    thisNode = popHeap(Q.Q)
    thisNode.rrtTreeCost = thisNode.rrtLMC

    listItem = thisNode.rrtNeighborsIn.front
    for i = 1:thisNode.rrtNeighborsIn.length
      neighborNode = listItem.data.startNode

      updateQueue(Q, neighborNode, root, hyberBallRad)

      listItem = listItem.child   # iterate through list
    end
  end
end


############################### RRTx Functions ####################################
### Functions used for RRTx.                                                    ###
###################################################################################


# successor stack marker function (marks when a node is in the successor stack OS)
markOS(node::T) where {T} = (node.inOSQueue = true)

# successor stack unmarker function (un marks when a node is removed from OS)
unmarkOS(node::T) where {T} = (node.inOSQueue = false)

# successor stack queue check marker function (checks if the node is marked OS)
markedOS(node::T) where {T} = node.inOSQueue

# makes sure the node is in the priority queue
function verifyInQueue(Q::QT, node::TN) where {QT, TN}
  if markedQ(node)
    updateHeap(Q.Q, node)
  else
    addToHeap(Q.Q, node)
  end
end

# makes sure the node is in the OS Queue
# removes it from the normal Q if necessary
function verifyInOSQueue(Q::QT, node::TN) where {QT, TN}

  if markedQ(node)
    updateHeap(Q.Q, node)
    removeFromHeap(Q.Q, node)
  end

  if !markedOS(node)
    markOS(node)
    JlistPush(Q.OS, node)
  end
end

# removes members of the current neighbor list of thisNode that are too far away
function cullCurrentNeighbors(thisNode::T, hyberBallRad::Float64) where {T}

  # remove outgoing edges from thisNode that are now too long
  listItem = thisNode.rrtNeighborsOut.front
  while listItem != listItem.child
    nextItem = listItem.child # since we may remove listItem from the list
    if listItem.data.dist > hyberBallRad
      neighborEdge = listItem.data
      neighorNode = neighborEdge.endNode
      JlistRemove(thisNode.rrtNeighborsOut, neighborEdge.listItemInStartNode)
      JlistRemove(neighorNode.rrtNeighborsIn, neighborEdge.listItemInEndNode)
    end
    listItem = nextItem
  end
end


# RRTx based version
# this returns the JlistNode containing the next outgoing neighbor edge of the
# node for which this iterator was created
function nextOutNeighbor(It::RRTNodeNeighborIterator{T, TE}, Q::rrtXQueue) where {T, TE}

  if It.listFlag == 0
    It.listItem = It.thisNode.InitialNeighborListOut.front
    It.listFlag = 1
  else
    It.listItem = It.listItem.child
  end
  while It.listItem == It.listItem.child
    # go to the next place that neighbors are stored
    if It.listFlag == 1
      It.listItem = It.thisNode.rrtNeighborsOut.front
    else
      # done with all neighbors
      return nothing
    end
    It.listFlag += 1
  end

  return It.listItem
end


# RRTx based version
# this returns the JlistNode containing the next outgoing neighbor edge of the
# node for which this iterator was created
function nextInNeighbor(It::RRTNodeNeighborIterator{T, TE}, Q::rrtXQueue) where {T, TE}

  if It.listFlag == 0
    It.listItem = It.thisNode.InitialNeighborListIn.front
    It.listFlag = 1
  else
    It.listItem = It.listItem.child
  end
  while It.listItem == It.listItem.child
    # go to the next place that neighbors are stored
    if It.listFlag == 1
      It.listItem = It.thisNode.rrtNeighborsIn.front
    else
      # done with all neighbors
      return nothing
    end
    It.listFlag += 1
  end

  return It.listItem
end



# makes newParent the parent of node via the edge
function makeParentOf(newParent::T, node::T, edge::Edge, root::T) where {T}

  ## error checking
  #if newParent.rrtParentUsed && newParent.rrtParentEdge.endNode == node
  #  error("trying to make node's parent use node as its parent")
  #end
  #if newParent == node
  #  error("trying to make node its own parent")
  #end

  # remove the node from its old parent's successor list
  if node.rrtParentUsed
    JlistRemove(node.rrtParentEdge.endNode.SuccessorList, node.successorListItemInParent)
  end

  # make newParent the parent of node
  node.rrtParentEdge = edge
  node.rrtParentUsed = true

  # place a (non-trajectory) reverse edge into newParent's sucessor
  # list and save a pointer to its postion in that list. This edge
  # is used to help keep track of successors and not for movement.
  backEdge = newEdge(newParent, node)
  backEdge.dist = Inf
  JlistPush(newParent.SuccessorList, backEdge, Inf)
  node.successorListItemInParent = newParent.SuccessorList.front

end


# recalculates LMC based on neighbors
function recalculateLMCMineVTwo(Q::TQ, thisNode::T, root::T,
  hyberBallRad::Float64) where {TQ, T}

  if thisNode == root
    return
  end

  newParentFound = false
  rrtParent = nothing
  rrtParentEdge = nothing

  # remove outdated nodes from current neighbors list
  cullCurrentNeighbors(thisNode, hyberBallRad)

  # get an iterator for this nodes neighbors
  thisNodeOutNeighbors = RRTNodeNeighborIterator{T, Edge{T}}(thisNode)

  # set the iterator to the first neighbor
  listItem = nextOutNeighbor(thisNodeOutNeighbors, Q)

  while listItem != nothing
    neighborEdge = listItem.data
    neighborNode = neighborEdge.endNode
    neighborDist = neighborEdge.dist
    nextItem = listItem.child

    if markedOS(neighborNode) #|| !neighborNode.rrtParentUsed
      #neighborNode already in OS queue (orphaned) or unwired

      listItem = nextOutNeighbor(thisNodeOutNeighbors, Q)
      continue
    end


    if (thisNode.rrtLMC > neighborNode.rrtLMC + neighborDist &&
        (!neighborNode.rrtParentUsed ||
         neighborNode.rrtParentEdge.endNode != thisNode) &&
        validMove(Q.S, neighborEdge)) # not sure why this last part of this check is necessary!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! if it is, we need modify the paper algorithm writeup to include it

      # found a better parrent
      thisNode.rrtLMC = neighborNode.rrtLMC + neighborDist
      rrtParent = neighborNode
      rrtParentEdge = listItem.data
      newParentFound = true
    end

    listItem = nextOutNeighbor(thisNodeOutNeighbors, Q)
  end

  if newParentFound # this node found a viable parent
    makeParentOf(rrtParent, thisNode, rrtParentEdge, root)
  end
end


# takes care of inserting a new node
function extend(S::TCS, KD::TKD, Q::rrtXQueue,
  newNode::T, closestNode::T, delta::Float64, hyberBallRad::Float64,
  moveGoal::T) where {T, TCS, TKD}

  # find all nodes within the (shrinking) hyperball of (saturated) newNode
  nodeList = kdFindWithinRange(KD, hyberBallRad, newNode.position)

  # try to find and link to best parent, this also saves the edges from newNode
  # to the neighbors in the field "tempEdge" of the neighbors. This saves time in
  # the case that trajectory calculation is complicated.
  findBestParent(S, newNode, nodeList, closestNode, true)

  # if no parent was found then ignore this node
  if !newNode.rrtParentUsed
    emptyRangeList(nodeList)     # clean up
    return
  end

  # for RRTx, need to add the new node to its parent's successor list
  # place a (non-trajectory) reverse edge into newParent's sucessor
  # list and save a pointer to its postion in that list. This edge
  # is used to help keep track of successors and not for movement.
  parentNode = newNode.rrtParentEdge.endNode
  backEdge = newEdge(parentNode, newNode)
  backEdge.dist = Inf
  JlistPush(parentNode.SuccessorList, backEdge, Inf)
  newNode.successorListItemInParent = parentNode.SuccessorList.front

  # insert the new node into the KDTree
  kdInsert(KD, newNode)

  # second pass, if there was a parent, then link with neighbors  and rewire
  # neighbors that would do better to use newNode as their parent.
  # Note that the edges -from- newNode -to- its neighbors have been stored
  # in "tempEdge" field of the neighbors
  listItem = nodeList.front
  for i = 1:nodeList.length
    nearNode = listItem.data

    # if edge from newNode to nearNode was valid
    if listItem.key != Inf
      # add to initial out neighbor list of newNode
      # (allows information propogation from newNode to nearNode always)
      makeInitialOutNeighborOf(nearNode, newNode, nearNode.tempEdge)

      # add to current neighbor list of newNode (allows information propogation
      # from newNode to nearNode and vise verse, but only while they are in D-Ball)
      makeNeighborOf(nearNode, newNode, nearNode.tempEdge)
    end

    # in the general case the trajectories along edges are not simply
    # the reverse of each other, therefore we need to calculate and
    # check the trajectory along the edge from nearNode to newNode.
    thisEdge = newEdge(nearNode, newNode)
    calculateTrajectory(S, thisEdge)

    if validMove(S, thisEdge) && !explicitEdgeCheck(S, thisEdge)
      # add to initial in neighbor list of newNode
      # (allows information propogation from newNode to nearNode always)
      makeInitialInNeighborOf(newNode, nearNode, thisEdge)


      # add to current neighbor list of newNode (allows information propogation
      # from newNode to nearNode and vise verse, but only while they are in D-Ball)
      makeNeighborOf(newNode, nearNode, thisEdge)
    else
      # edge cannot be created
      listItem = listItem.child   # iterate thorugh list
      continue
    end

    # rewire neighbors that would do better to use this node as their parent
    # unless they are not in the relivant portion of the space vs. moveGoal
    if (nearNode.rrtLMC > newNode.rrtLMC + thisEdge.dist &&
        newNode.rrtParentEdge.endNode != nearNode &&
        newNode.rrtLMC + thisEdge.dist < moveGoal.rrtLMC)

      # make this node the parent of the neighbor node
      makeParentOf(newNode, nearNode, thisEdge, KD.root)

      # recalculate tree cost of neighbor
      oldLmc = nearNode.rrtLMC
      nearNode.rrtLMC = newNode.rrtLMC + thisEdge.dist

      # insert neighbor into priority queue if cost reduction is great enough
      if oldLmc - nearNode.rrtLMC > Q.changeThresh && nearNode != KD.root
        verifyInQueue(Q, nearNode)
      end
    end

    listItem = listItem.child   # iterate thorugh list
  end
  emptyRangeList(nodeList)     # clean up

  # insert the node into the piority queue
  addToHeap(Q.Q, newNode)
end

function extendTerrain(S::TCS, KD::TKD, Q::rrtXQueue,
  newNode::T, closestNode::T, delta::Float64, hyberBallRad::Float64,
  moveGoal::T) where {T, TCS, TKD}

  # find all nodes within the (shrinking) hyperball of (saturated) newNode
  nodeList = kdFindWithinRange(KD, hyberBallRad, newNode.position)

  # try to find and link to best parent, this also saves the edges from newNode
  # to the neighbors in the field "tempEdge" of the neighbors. This saves time in
  # the case that trajectory calculation is complicated.
  findBestParent(S, newNode, nodeList, closestNode, true)

  # if no parent was found then ignore this node
  if !newNode.rrtParentUsed
    emptyRangeList(nodeList)     # clean up
    return
  end

  # for RRTx, need to add the new node to its parent's successor list
  # place a (non-trajectory) reverse edge into newParent's sucessor
  # list and save a pointer to its postion in that list. This edge
  # is used to help keep track of successors and not for movement.
  parentNode = newNode.rrtParentEdge.endNode
  backEdge = newEdge(parentNode, newNode)
  backEdge.dist = Inf
  JlistPush(parentNode.SuccessorList, backEdge, Inf)
  newNode.successorListItemInParent = parentNode.SuccessorList.front

  # insert the new node into the KDTree
  kdInsert(KD, newNode)

  # second pass, if there was a parent, then link with neighbors  and rewire
  # neighbors that would do better to use newNode as their parent.
  # Note that the edges -from- newNode -to- its neighbors have been stored
  # in "tempEdge" field of the neighbors
  listItem = nodeList.front
  for i = 1:nodeList.length
    nearNode = listItem.data

    # if edge from newNode to nearNode was valid
    if listItem.key != Inf
      # add to initial out neighbor list of newNode
      # (allows information propogation from newNode to nearNode always)
      makeInitialOutNeighborOf(nearNode, newNode, nearNode.tempEdge)

      # add to current neighbor list of newNode (allows information propogation
      # from newNode to nearNode and vise verse, but only while they are in D-Ball)
      makeNeighborOf(nearNode, newNode, nearNode.tempEdge)
    end

    # in the general case the trajectories along edges are not simply
    # the reverse of each other, therefore we need to calculate and
    # check the trajectory along the edge from nearNode to newNode.
    thisEdge = newEdge(nearNode, newNode)
    calculateTrajectoryTerrain(S, thisEdge)

    if validMove(S, thisEdge) && !explicitEdgeCheck(S, thisEdge)
      # add to initial in neighbor list of newNode
      # (allows information propogation from newNode to nearNode always)
      makeInitialInNeighborOf(newNode, nearNode, thisEdge)


      # add to current neighbor list of newNode (allows information propogation
      # from newNode to nearNode and vise verse, but only while they are in D-Ball)
      makeNeighborOf(newNode, nearNode, thisEdge)
    else
      # edge cannot be created
      listItem = listItem.child   # iterate thorugh list
      continue
    end

    # rewire neighbors that would do better to use this node as their parent
    # unless they are not in the relivant portion of the space vs. moveGoal
    if (nearNode.rrtLMC > newNode.rrtLMC + thisEdge.dist &&
        newNode.rrtParentEdge.endNode != nearNode &&
        newNode.rrtLMC + thisEdge.dist < moveGoal.rrtLMC)

      # make this node the parent of the neighbor node
      makeParentOf(newNode, nearNode, thisEdge, KD.root)

      # recalculate tree cost of neighbor
      oldLmc = nearNode.rrtLMC
      nearNode.rrtLMC = newNode.rrtLMC + thisEdge.dist

      # insert neighbor into priority queue if cost reduction is great enough
      if oldLmc - nearNode.rrtLMC > Q.changeThresh && nearNode != KD.root
        verifyInQueue(Q, nearNode)
      end
    end

    listItem = listItem.child   # iterate thorugh list
  end
  emptyRangeList(nodeList)     # clean up

  # insert the node into the piority queue
  addToHeap(Q.Q, newNode)
end


## this is the (non-initial) rewire function used by RRTx that is responsible
## for propogating changes through the graph
function rewire(Q::TQ, thisNode::T, root::T, hyberBallRad::Float64,
  changeThresh::Float64) where {T, TQ}

  # only explicitly propogate changes if they are large enough
  deltaCost = thisNode.rrtTreeCost - thisNode.rrtLMC
  if deltaCost <= changeThresh
    # note that using simply "<" (and not "<=") causes problems
    # actaully, I think that the above note may now be outdated
    #thisNode.rrtTreeCost = thisNode.rrtLMC !!!! Now happens after return
    return
  end

  # remove outdated nodes from current neighbors list
  cullCurrentNeighbors(thisNode, hyberBallRad)

  # get an iterator for this node's neighbors and iterate through list
  thisNodeInNeighbors = RRTNodeNeighborIterator{T, Edge{T}}(thisNode)
  listItem = nextInNeighbor(thisNodeInNeighbors, Q)
  while listItem != nothing
    neighborEdge = listItem.data
    neighborNode = neighborEdge.startNode

    # ignor this node's parent and also nodes that [cannot reach thisNode due to
    # dynamics of robot or space] -- not sure if we need [...] part since neighbors
    # are not initially created that cannot reach this node
    if ((thisNode.rrtParentUsed && thisNode.rrtParentEdge.endNode == neighborNode) ||
        !validMove(Q.S, neighborEdge))

      listItem = nextInNeighbor(thisNodeInNeighbors, Q)
      continue
    end

    neighborEdge = listItem.data

    deltaCostNeighbor = neighborNode.rrtLMC - (thisNode.rrtLMC + neighborEdge.dist)
    if deltaCostNeighbor > 0
      # neighborNode should use thisNode as its parent (note it may already be)

      neighborNode.rrtLMC = thisNode.rrtLMC + neighborEdge.dist

      if(!neighborNode.rrtParentUsed || neighborNode.rrtParentEdge.endNode != thisNode)
        makeParentOf(thisNode, neighborNode, neighborEdge, root)
      end

      # if the reduction is great enough, then propogate through the neighbor
      if neighborNode.rrtTreeCost - neighborNode.rrtLMC > changeThresh
        verifyInQueue(Q, neighborNode)
      end
    end

    listItem = nextInNeighbor(thisNodeInNeighbors, Q)
  end
end


# propogates cost information through the graph
function reduceInconsistency(Q::rrtXQueue, goalNode::T,
  robotRad::Float64, root::T, hyberBallRad::Float64) where {T}

  while Q.Q.indexOfLast > 0 && (lessQ(topHeap(Q.Q), goalNode) || goalNode.rrtLMC == Inf || goalNode.rrtTreeCost == Inf || markedQ(goalNode))

    thisNode = popHeap(Q.Q)

    # update neighbors of thisNode if it has changed more than change thresh
    if thisNode.rrtTreeCost - thisNode.rrtLMC > Q.changeThresh
      recalculateLMCMineVTwo(Q, thisNode, root, hyberBallRad)
      rewire(Q, thisNode, root, hyberBallRad, Q.changeThresh)
    end
    thisNode.rrtTreeCost = thisNode.rrtLMC
  end
end


# propogate orphan status to all nodes in the basin(s) of attraction of
# the nodes in Q.OS (that have higher cost). This also takes the robot
# to remember if node the robot was moving at is one of the nodes that
# has become an orphan
function propogateDescendants(Q::T, R::RT) where {T, RT}

  if Q.OS.length <= 0
    return
  end

  # first pass, accumulate all such nodes in a single list, and mark them as
  # belonging in that list, we'll just use the OS stack that we've been using
  # adding nodes to the front while moving from back to front
  OS_list_item = Q.OS.back
  while OS_list_item != OS_list_item.parent
    thisNode = OS_list_item.data

    # add all of this' node's successors to OS Stack
    SuccessorList_item = thisNode.SuccessorList.front
    while SuccessorList_item != SuccessorList_item.child
      sucessorNode = SuccessorList_item.data.endNode

      verifyInOSQueue(Q, sucessorNode) # pushes to front of OS

      SuccessorList_item = SuccessorList_item.child
    end

    OS_list_item = OS_list_item.parent
  end


  # second pass, put all -out neighbors- of the nodes in OS (not including nodes
  # in OS) into Q and tell them to force rewire. Note going back to front makes Q
  # adjustements slightly faster, since nodes near the front tend to have higher
  # costs
  OS_list_item = Q.OS.back
  while OS_list_item != OS_list_item.parent
    thisNode = OS_list_item.data

    # get an iterator for this node's neighbors
    thisNodeOutNeighbors = RRTNodeNeighborIterator{RRTNode{Float64}, Edge{RRTNode{Float64}}}(thisNode)

    # now iterate through list (add all neighbors to the Q, except those in OS)
    listItem = nextOutNeighbor(thisNodeOutNeighbors, Q)
    while listItem != nothing
      neighborNode = listItem.data.endNode

      if markedOS(neighborNode) #|| !neighborNode.rrtParentUsed
        #neighborNode already in OS queue (orphaned) or unwired
        listItem = nextOutNeighbor(thisNodeOutNeighbors, Q)
        continue
      end

      # otherwise, make sure that neighborNode is in normal queue
      neighborNode.rrtTreeCost = Inf # node will be inserted with LMC key and then
                                     # guarenteed to propogate cost forward since
                                     # usefull nodes have rrtLMC < rrtTreeCost
      verifyInQueue(Q, neighborNode)

      listItem = nextOutNeighbor(thisNodeOutNeighbors, Q)
    end

    # add parent to the Q, unless it is in OS
    if thisNode.rrtParentUsed && !markedOS(thisNode.rrtParentEdge.endNode)
      # thisNode.rrtParent.rrtTreeCost = Inf
      thisNode.rrtParentEdge.endNode.rrtTreeCost = Inf
      verifyInQueue(Q, thisNode.rrtParentEdge.endNode)
    end

    OS_list_item = OS_list_item.parent
  end

  # third pass, remove all nodes from OS, unmark them, and remove
  # their connections to their parents, and if one was the robot's target
  # then take approperiate measures
  while Q.OS.length > 0

    thisNode = JlistPop(Q.OS)
    unmarkOS(thisNode)

    if thisNode == R.nextMoveTarget
      R.currentMoveInvalid = true
    end

    if thisNode.rrtParentUsed
      # remove thisNode from its parent's successor list
      JlistRemove(thisNode.rrtParentEdge.endNode.SuccessorList, thisNode.successorListItemInParent)

      # thisNode now has no parent
      thisNode.rrtParentEdge = newEdge(thisNode, thisNode)
      thisNode.rrtParentEdge.dist = Inf
      thisNode.rrtParentUsed = false
    end

    thisNode.rrtTreeCost = Inf
    thisNode.rrtLMC = Inf
  end
end




# if C-Space has a time dimension, add a sequence of descendents
# to the root, where each {(great)^n}-grandchild is at the same
# position as root, but at a sequence of times from 0 to the last ("earliest")
# time that the robot could arrive at that position (assuming no obstacles)
# but at a tree distance defined to be 0 from the root.
# this helps the robot to reach the goal location as quickly as possible instead
# of burning time
function addOtherTimesToRoot(S::TS, KD::TKD, goal::RRTNode{Float64},
  root::RRTNode{Float64}, searchType::String) where {TS, TKD}

  insertStep = 2.0

  lastTimeToInsert = goal.position[3] - Wdist(root,goal)/S.robotVelocity
  firstTimeToInsert = S.start[3] + insertStep
  previousNode = root
  safeToGoal = true
  for timeToInsert = firstTimeToInsert:insertStep:lastTimeToInsert
    newPose = copy(root.position)
    newPose[3] = timeToInsert

    newNode = RRTNode{Float64}(newPose)

    # edge from newNode to previousNode
    thisEdge = newEdge(newNode, previousNode)
    calculateHoverTrajectory(S, thisEdge)

    if searchType == "RRT*"

      # make this node the parent of the neighbor node
      newNode.rrtParentEdge = thisEdge
      newNode.rrtParentUsed = true

    elseif searchType == "RRT#"

      # make this node the parent of the neighbor node
      newNode.rrtParentEdge = thisEdge
      newNode.rrtParentUsed = true

      makeNeighborOf(newNode, previousNode, thisEdge)

    elseif searchType == "RRTx"

      makeParentOf(previousNode, newNode, thisEdge, root)
      makeInitialOutNeighborOf(previousNode, newNode, thisEdge)

      # initial neighbor list edge
      makeInitialInNeighborOf(newNode, previousNode, thisEdge)
    end


    # make sure that this edge is safe
    if explicitEdgeCheck(S, thisEdge)
      # not safe
      thisEdge.dist = Inf
      safeToGoal = false
      newNode.rrtLMC = Inf;
      newNode.rrtTreeCost = Inf;
    elseif safeToGoal
      # if the edge has safe path all the way to the "real" goal then
      # we make the cost of reaching "real" goal 0 from newNode
      thisEdge.dist = 0.0
      newNode.rrtLMC = 0.0;
      newNode.rrtTreeCost = 0.0;
    else
      thisEdge.dist = Inf
      newNode.rrtLMC = Inf;
      newNode.rrtTreeCost = Inf;
    end


    kdInsert(KD, newNode)


    previousNode = newNode
  end
end

# attempts to find a new move target for the robot, places it into Robot data
# (used when the old target has become invalid)
function findNewTarget(S::TS, KD::TKD, R::RobotData, hyberBallRad::Float64) where {TS, TKD}
  # start by looking at a hyperball of possible targets with a radius
  # determined by max of {hyberBallRad} and {the previous edge
  # legth (the edge that is now invalud)}. If this fails, then we try larger
  # and larger balls until we have looked at all nodes, if still cannot
  # find a target then the robot is out of luck, and we do nothing

  R.robotEdgeUsed = false
  R.distAlongRobotEdge = 0.0
  R.timeAlongRobotEdge = 0.0
  R.robotEdgeForPlottingUsed = false
  R.distAlongRobotEdgeForPlotting = 0.0
  R.timeAlongRobotEdgeForPlotting = 0.0

  println("move target has become invalid")
  searchBallRad = max(hyberBallRad, dist(R.robotPose, R.nextMoveTarget.position))
  maxSearchBallRad = dist(S.lowerBounds, S.upperBounds)
  searchBallRad = min(searchBallRad, maxSearchBallRad)
  L = kdFindWithinRange(KD, searchBallRad, R.robotPose)
  dummyRobotNode = RRTNode{Float64}(R.robotPose) # a temp node at robot pose
  edgeToBestNeighbor = Edge{RRTNode{Float64}}

  while true # will break out when done
    println("searching for new target within radius $(searchBallRad)")

    bestDistToNeighbor = Inf
    bestDistToGoal = Inf
    bestNeighbor = nothing

    ptr = L.front
    while ptr != ptr.child
      neighborNode = ptr.data

      thisEdge = newEdge(dummyRobotNode, neighborNode)
      calculateTrajectory(S, thisEdge)

      if validMove(S, thisEdge) && !explicitEdgeCheck(S, thisEdge)
        # a safe point was found, see if it is the best so far

        distToGoal = neighborNode.rrtLMC + thisEdge.dist
        if distToGoal < bestDistToGoal && validMove(S, thisEdge)

          # found a new and better neighbor
          bestDistToGoal = distToGoal
          bestDistToNeighbor = thisEdge.dist
          bestNeighbor = neighborNode
          edgeToBestNeighbor = thisEdge
        end
      end

      ptr = ptr.child
    end
    # done trying to find a target within ball of searchBallRad

    # if a valid neighbor was found, then use it
    if bestDistToGoal != Inf
      R.nextMoveTarget = bestNeighbor
      R.distanceFromNextRobotPoseToNextMoveTarget = bestDistToNeighbor
      R.currentMoveInvalid = false
      println("found a valid move target")


      R.robotEdge = edgeToBestNeighbor
      R.robotEdgeForPlotting = edgeToBestNeighbor
      R.robotEdgeUsed = true
      R.robotEdgeForPlottingUsed = true

      if S.spaceHasTime
        R.timeAlongRobotEdge = 0.0            # note this is updated before robot moves
        R.timeAlongRobotEdgeForPlotting = 0.0
      else
        R.distAlongRobotEdge = 0.0            # note this is updated before robot moves
        R.distAlongRobotEdgeForPlotting = 0.0
      end


      # set moveGoal to be nextMoveTarget NOTE, may want to actually insert a new
      # node at the robot's position and use that instead, since these
      # "edges" created between robot pose and R.nextMoveTarget may be lengthy
      S.moveGoal.isMoveGoal = false
      S.moveGoal = R.nextMoveTarget
      S.moveGoal.isMoveGoal = true

      break
    end

    searchBallRad *= 2
    if searchBallRad > maxSearchBallRad
      println("unable to find a valid move target")
      error("unable to find a valid move target")
      break
    end
    kdFindMoreWithinRange(KD, searchBallRad, R.robotPose, L)
  end
  emptyRangeList(L)    # cleanup
end


# move robot the distance that it would move in slice_time time
# if time is not a dimension of the C-Space, then a constant velocity is assumed.
# This also updates the moveGoal in the event that the robot has lost connectivity
# with the graph due to dynamic obstacles breaking the first edge of its path
function moveRobot(S::TS, Q::TQ, KD::TKD, slice_time::Float64, root::RRTNode{Float64},
  hyberBallRad::Float64, R::RobotData) where {TS, TQ, TKD}

  # start by updating the location of the robot based on how it moved
  # since the last update (as well as the total path that it has followed)
  if R.moving
    R.robotPose = R.nextRobotPose
    R.robotMovePath[R.numRobotMovePoints+1:R.numRobotMovePoints+R.numLocalMovePoints,:] = R.robotLocalPath[1:R.numLocalMovePoints,:]
    R.numRobotMovePoints += R.numLocalMovePoints

    if !S.spaceHasTime
      println("new robot pose: [$(R.robotPose[1]) $(R.robotPose[2])]")
    else
      println("new robot pose: [$(R.robotPose[1]) $(R.robotPose[2]) $(R.robotPose[3])]")
    end

    # save stuff for plotting
    if S.spaceHasTime
      R.timeAlongRobotEdgeForPlotting = R.timeAlongRobotEdge
    else
      R.distAlongRobotEdgeForPlotting = R.distAlongRobotEdge
    end
    R.robotEdgeForPlotting = R.robotEdge
    R.robotEdgeForPlottingUsed = true


  else
    # movement has just started, so remember that the robot is now moving
    R.moving = true

    if !S.moveGoal.rrtParentUsed
      # no parent has been found for the node at the robot's position
      R.currentMoveInvalid = true
    else
      R.robotEdge = S.moveGoal.rrtParentEdge
      R.robotEdgeForPlotting = R.robotEdge
      R.robotEdgeUsed = true
      R.robotEdgeForPlottingUsed = true

      if S.spaceHasTime
        R.timeAlongRobotEdge = 0.0
        R.timeAlongRobotEdgeForPlotting = 0.0
      else
        R.distAlongRobotEdge = 0.0
        R.distAlongRobotEdgeForPlotting = 0.0
      end

    end
  end


  # if the robot's current move target has been invalidated due to
  # dynamic obstacles then we need to attempt to find a new (safe) move target
  # NOTE we handle newly invalid moveTarget after moving the robot (since
  # the robot has already moved this time slice)
  if R.currentMoveInvalid
    findNewTarget(S, KD, R, hyberBallRad)
  else
    # recall that moveGoal is the node who's key is used to determing the
    # level set of cost propogation (this should theoretically be further
    # than the robot from the root of the tree, which will happen here
    # assuming that robot moves at least one edge each slice time. even if
    # that does not happen, things will still be ok in practice as long as
    # robot is "close" to moveGoal
    S.moveGoal.isMoveGoal = false
    S.moveGoal = R.nextMoveTarget
    S.moveGoal.isMoveGoal = true
  end


  # finally, we calculate the point to which the robot will move to in slice_time
  # and remember it for the next time this function is called. we also remember
  # all the nodes that it will visit along the way in the local path
  # and the part of the edge trajectory that takes the robot to the first local
  # point (the latter two things are used for visuialization)

  if !S.spaceHasTime
    # not using time dimension, so assume speed is equal to robotVelocity

    nextNode = R.nextMoveTarget

    # calculate distance from robot to the end of the current edge it is following
    nextDist = R.robotEdge.dist - R.distAlongRobotEdge

    distRemaining = S.robotVelocity*slice_time   # dist left for robot to move

    # save first local path point
    R.numLocalMovePoints = 1
    R.robotLocalPath[R.numLocalMovePoints,:] = copy(R.robotPose)

    # starting at current location (and looking ahead to nextNode), follow parent
    # pointers back for the approperiate distance (or root or dead end)
    while (nextDist <= distRemaining && nextNode != root &&
           nextNode.rrtParentUsed && nextNode != nextNode.rrtParentEdge.endNode)

      # can go all the way to nextNode and still have some distance left to spare

      # remember the robot will move through this point
      R.numLocalMovePoints += 1
      R.robotLocalPath[R.numLocalMovePoints,:] = copy(nextNode.position)

      # recalculate remaining distance
      distRemaining -= nextDist

      # reset distance along edge
      R.distAlongRobotEdge = 0.0

      # update trajectory that the robot will be in the middle of
      R.robotEdge = nextNode.rrtParentEdge
      R.robotEdgeUsed = true

      # calculate the dist of that trajectory
      nextDist = R.robotEdge.dist

      # update the next node (at the end of that trajectory)
      nextNode = R.robotEdge.endNode

    end


    # either: 1) nextDist > distRemaining
    # (or)    2) the path we were following now ends at nextNode

    # calculate next pose of the robot
    if nextDist > distRemaining
      R.distAlongRobotEdge += distRemaining
      R.nextRobotPose = poseAtDistAlongEdge(R.robotEdge, R.distAlongRobotEdge)
    else
      # the next node is the end of this tree and we reach it
      R.nextRobotPose = nextNode.position
      R.distAlongRobotEdge = R.robotEdge.dist
    end

    R.nextMoveTarget = R.robotEdge.endNode

    # remember last point in local path
    R.numLocalMovePoints += 1
    R.robotLocalPath[R.numLocalMovePoints,:] = R.nextRobotPose

  else # S.spaceHasTime
    # space has time, so path is parameterized by time as well
    nextNode = R.nextMoveTarget

    # save first local path point
    R.numLocalMovePoints = 1
    R.robotLocalPath[R.numLocalMovePoints,:] = copy(R.robotPose)

    targetTime = R.robotPose[3] - slice_time
    while (targetTime < R.robotEdge.endNode.position[3] &&
           nextNode != root && nextNode.rrtParentUsed &&
           nextNode != nextNode.rrtParentEdge.endNode)

      # can go all the way to nextNode and still have some time left to spare

      # remember the robot will move through this point
      R.numLocalMovePoints += 1
      R.robotLocalPath[R.numLocalMovePoints,:] = copy(nextNode.position)

      # update trajectory that the robot will be in the middle of
      R.robotEdge = nextNode.rrtParentEdge
      R.robotEdgeUsed = true

      # update the next node (at the end of that trajectory)
      nextNode = nextNode.rrtParentEdge.endNode
    end


    # either: 1) targetTime >= nextNode.position[3]
    # (or)    2) the path we were following now ends at nextNode

    # calculate next pose of the robot
    if targetTime >= nextNode.position[3]
      R.timeAlongRobotEdge = R.robotEdge.startNode.position[3] - targetTime
      R.nextRobotPose = poseAtTimeAlongEdge(R.robotEdge, R.timeAlongRobotEdge)

    else
      # the next node is the end of this tree and we reach it
      R.nextRobotPose = nextNode.position
      R.timeAlongRobotEdge = R.robotEdge.startNode.position[3] - R.robotEdge.endNode.position[3]
    end

    R.nextMoveTarget = R.robotEdge.endNode

    # remember last point in local path
    R.numLocalMovePoints += 1
    R.robotLocalPath[R.numLocalMovePoints,:] = R.nextRobotPose
  end
end



# this returns a -rangeList- (see the KD tree code) containing all points
# that are in conflict with the obstacle, NOTE that the rangeList
# must be DESTROYED PROPERLY using emptyRangeList(L) to avoid problems
function findPointsInConflictWithObstacle(S::TS, KD::TKD, ob::Obstacle, root::T) where{TS, TKD, T}

  # find points that are in conflict with the obstacle, start by over-estimate,
  # finding all within the following bounding hyper-sphere(s)
  L::JList{T} = JList{typeof(KD.root)}()
  if 1 <= ob.kind <= 5
    # 2D obstacle
    if !S.spaceHasTime && !S.spaceHasTheta
      # euclidian space without time
      searchRange = S.robotRadius + S.delta + ob.radius
      L = kdFindWithinRange(KD, searchRange, ob.position)
    elseif !S.spaceHasTime && S.spaceHasTheta
      # Dubin's robot without time [x,y, 0.0, theta]
      searchRange = S.robotRadius + S.delta + ob.radius + pi
      obsCenterDubins = [ob.position[1:2]' 0.0 pi]
      L = kdFindWithinRange(KD, searchRange, obsCenterDubins)
    else
      error( "this type of obstacle not coded for this type of space")
    end
  elseif 6 <= ob.kind <= 7
    # 2D obstacle with time, find points within range of each point along the
    # time path, accumulating all points that are in any of the bounding
    # hyper-spheres

    baseSearchRange = S.robotRadius + S.delta + ob.radius

    ## for debugging only:
    #AllqueryPoses = zeros(Float64, size(ob.path,1)-1, 3)

    for i = 1:size(ob.path,1)
      # make query pose the middle of the edge, and add 1/2 edge length
      # through the C-space to the baseSearchRange (this is an overestimate)

      if size(ob.path,1) == 1
        j = 1
      else
        j = i+1
      end

      val = reshape(ob.path[i,:] + ob.path[j,:], 1, length(ob.path[i,:]))
      queryPose = [ob.position 0.0] + val/2.0

      if S.spaceHasTheta
        # Dubins Car
        queryPose = [queryPose pi]
      end

      ## for debugging only:
      #AllqueryPoses[i-1,:] = queryPose

      searchRange = baseSearchRange + euclidianDist(ob.path[i,:], ob.path[j,:])/2.0

      if S.spaceHasTheta
        searchRange += pi
      end


      if i == 1
        L = kdFindWithinRange(KD, searchRange, queryPose)
      else
        kdFindMoreWithinRange(KD, searchRange, queryPose, L)
      end

      if j == size(ob.path,1)
        break
      end
    end

    ## for debugging only:
    #saveData(AllqueryPoses, "temp/Qnodes_$(fileCounter).txt")
  else
    error("this case not coded yet")
  end
  return L
end


# This adds the obstacle (checks for edge conflicts with the obstacle and then
# puts the affected nodes into the approperiate heaps)
function addNewObstacle(S::TS, KD::TKD, Q::TQ,
  ob::Obstacle, root::T, fileCounter::Int, R::TR) where {TS, TKD, TQ, T, TR}
  ob.obstacleUnused = false

  # find all point that are in conflict with the obstacle
  L = findPointsInConflictWithObstacle(S, KD, ob, root)

  ## save data for debugging:
  #points = zeros(Float64, L.length, 3)
  #ptct = 0


  #for all nodes that might be in conflict
  while L.length > 0
    (thisNode, key) = popFromRangeList(L)
    # check all of their edges

    #### see if this node's neighbors can be reached

    # get an iterator for this node's (out) neighbors
    thisNodeOutNeighbors = RRTNodeNeighborIterator{T, Edge{T}}(thisNode)

    # now iterate through list
    listItem = nextOutNeighbor(thisNodeOutNeighbors, Q)
    while listItem != nothing
      neighborEdge = listItem.data
      nextItem = nextOutNeighbor(thisNodeOutNeighbors, Q) # we may alter neighborlists

      if explicitEdgeCheck(S, neighborEdge, ob)
        listItem.data.dist = Inf # mark edge to neighbor at inf cost
      end
      listItem = nextItem
    end

    #### done with seeing if this node's neighbors can be reached

    #### now see if this node's parent can be reached
    if thisNode.rrtParentUsed && explicitEdgeCheck(S, thisNode.rrtParentEdge, ob)

      # for debugging:
      #points[ptct + 1,:] = thisNode.position

      # remove thisNode from its parent's successor list
      JlistRemove(thisNode.rrtParentEdge.endNode.SuccessorList, thisNode.successorListItemInParent)

      # this node now has no parent
      thisNode.rrtParentEdge.endNode = thisNode
      thisNode.rrtParentEdge.dist = Inf
      thisNode.rrtParentUsed = false

      verifyInOSQueue(Q, thisNode)
    end
    #### done with seeing if this node's parent can be reached

    ## for debugging only:
    #if points[ptct+1,:] != [0.0 0.0 0.0]
    #  ptct +=1
    #end
  end

  ## for debugging only:
  #saveData(points[1:ptct,:], "temp/Cnodes_$(fileCounter).txt")

  # cleanup
  emptyRangeList(L)

  # now check the robot's current move to its target
  if R.robotEdgeUsed && explicitEdgeCheck(S, R.robotEdge, ob)
    R.currentMoveInvalid = true
  end
end


# This removes the obstacle (checks for edge conflicts with the obstacle and then
# puts the affected nodes into the approperiate heaps)
function removeObstacle(S::TS, KD::TKD, Q::TQ,
  ob::Obstacle, root::T, hyberBallRad::Float64, timeElapsed::Float64,
  moveGoal::T) where{T, TS, TKD, TQ}

  # find all point that were in conflict with the obstacle
  L = findPointsInConflictWithObstacle(S, KD, ob, root)

  #for all nodes that might be in conflict
  while L.length > 0
    (thisNode, key) = popFromRangeList(L)
    # check all of their edges

    #### see if this node's (out) neighbors were blocked by the obstacle

    # get an iterator for this nodes neighbors
    thisNodeOutNeighbors = RRTNodeNeighborIterator{T, Edge{T}}(thisNode)
    neighborsWereBlocked = false

    # now iterate through list
    listItem = nextOutNeighbor(thisNodeOutNeighbors, Q)
    while listItem != nothing
      neighborEdge = listItem.data
      neighborNode = listItem.data.endNode
      nextItem = nextOutNeighbor(thisNodeOutNeighbors, Q) # we may alter neighborlists

      if neighborEdge.dist == Inf && explicitEdgeCheck(S, neighborEdge, ob)
        # this edge used to be in collison with at least one obstacle,
        # and in particular, the obstacle in question, however it still may be in
        # conflict with another obstacle, so need to check

        list_item = S.obstacles.front
        conflictsWithOtherObs = false
        while list_item != list_item.child
          obOther = list_item.data
          if obOther != ob && !obOther.obstacleUnused && obOther.startTime <= timeElapsed <= (obOther.startTime + obOther.lifeSpan)
            if explicitEdgeCheck(S, neighborEdge, obOther)
              conflictsWithOtherObs = true
              break
            end
          end
          list_item = list_item.child
        end


        if !conflictsWithOtherObs
          # reset edge length to actual cost
          listItem.data.dist = listItem.data.distOriginal

          neighborsWereBlocked = true # remember that this node had neighbors that
                                      # were blocked
        end
      end
      listItem = nextItem
    end

    #### done with seeing if this node's neighbors were blocked by the obstacle
    if neighborsWereBlocked
      recalculateLMCMineVTwo(Q, thisNode, root, hyberBallRad)
      if thisNode.rrtTreeCost != thisNode.rrtLMC && lessQ(thisNode, moveGoal)
        verifyInQueue(Q, thisNode)
      end
    end
  end
  emptyRangeList(L)

  ob.obstacleUnused = true
end





############################# main (despite the name) ###########################
# the following now calls RRT, RRT*, RRT#, and RRTx. Behaviour is determined by #
# paramiters passed in, note that RRT# vs (RRT* and RRT) helper functions above #
# are called correctly using julia's multiple disbatch, where the type of       #
# queue being used is different for each algorithm                              #
#################################################################################


# S is the CSpace, the algorithm runs until either N nodes have been sampled
# or TimeOut seconds pass, delta is the saturation distance, ballConstant is
# the ball constant
function RRTX(S::TS, total_planning_time::Float64, slice_time::Float64,
  delta::Float64, ballConstant::Float64, changeThresh::Float64,
  searchType::String, MoveRobotFlag::Bool, saveVideoData::Bool,
  statsArgs...) where {TS}


  T = RRTNode{Float64}

  # NOTE THIS IS HARD CODED HERE (SHOULD PROBABLY MAKE INPUT ARGUMENT)
  robotSensorRange = 20.0 # used for "sensing" obstacles

  if length(statsArgs) >= 2
    dataFileName = statsArgs[2]
  else
    dataFileName = "data.txt"
  end

  startTime = time_ns()
  save_elapsed_time = 0.0 # will hold save time to correct for time spent
                          # writing to files (only used for visualization video)

  ### do initialization stuff:

  # init a KDtree that will be used
  # MAKE SURE THIS USES APPROPERIATE DISTANCE FUNCTION !!!!!!!!!
  if S.spaceHasTheta
    # using dubins car, 4th dimension wraps around at 0 == 2pi
    KD = KDTree{RRTNode{Float64}}(S.d, KDdist, [4], [2*pi])
  else
    KD = KDTree{RRTNode{Float64}}(S.d, KDdist)
  end

  # init queue. Note that for algorithms that do not require a queue this
  # is essentially an empty variable that is used to help with polymorphism,
  # which makes life easy, the code fast, and allows much of the same code
  # to be used for all algorithms. Note about RRTx: I decided to forgoe RRT*+
  # and focus only on #+ sice the latter will likely be more useful than the
  # former in practice
  if searchType == "RRT"
    Q = rrtQueue{Bool}()
  elseif searchType == "RRT*"
    Q = rrtStarQueue{Bool}()
  elseif searchType == "RRT#"
    Q = rrtSharpQueue{RRTNode{Float64}, typeof((Float64, Float64))}()
    Q.Q = BinaryHeap{RRTNode{Float64}, typeof((Float64, Float64))}(keyQ, lessQ, greaterQ, markQ, unmarkQ, markedQ, setIndexQ, unsetIndexQ, getIndexQ)
    Q.S = S

#Q.numReduces = 0 # for testing, take out when done

  elseif searchType == "RRTx"
    Q = rrtXQueue{RRTNode{Float64}, typeof((Float64, Float64))}()
    Q.Q = BinaryHeap{RRTNode{Float64}, typeof((Float64, Float64))}(keyQ, lessQ, greaterQ, markQ, unmarkQ, markedQ, setIndexQ, unsetIndexQ, getIndexQ)
    Q.OS = JList{RRTNode{Float64}}()
    Q.S = S
    Q.changeThresh = changeThresh
  else
    error("unknown search type: $(searchType)")
  end

  S.sampleStack = JList{Array{Float64,2}}() # stores a stack of points that we
                                            # desire to insert in the future in
                                            # (used when an obstacle is removed
  S.delta = delta


  robotRads = S.robotRadius




  # define root node in the search tree graph
  root = RRTNode{Float64}(S.start)

  # explicit check root
  (explicitlyUnSafe, unused) = explicitNodeCheck(S, root)
  if explicitlyUnSafe
    error("root is not safe")
  end

  root.rrtTreeCost = 0.0
  root.rrtLMC = 0.0

  # insert the root into the KDtree
  kdInsert(KD, root)


  # define a goal node
  goal = RRTNode{Float64}(S.goal)
  goal.rrtTreeCost = Inf
  goal.rrtLMC = Inf
  S.goalNode = goal
  S.root = root

  S.moveGoal = goal # this will store a node at least as far from the root as the robot
                    # during movement it key is used to limit propogation beyond the
                    # region we care about

  S.moveGoal.isMoveGoal = true

  # paramiters that have to do with the robot path following simulation
  R = RobotData{RRTNode{Float64}}(copy(S.goal), goal, 20000)

  vCounter = 0; # helps with visuilizing data
  S.fileCtr = vCounter

  sliceCounter = 0; # helps with saving accurate time data


  if S.spaceHasTime
    # add other "times" to root of tree
    addOtherTimesToRoot(S, KD, goal, root, searchType)
  end
  ### end of initialization stuff


  # if saving stats about run, then allocate memory to store data
  if (length(statsArgs) >= 1 && statsArgs[1])

    savingStats = true
    estimated_number_of_data_points = 4*Int(ceil(total_planning_time/slice_time))

    checkPtr = 1      # array location to save stats

    itOfCheck = Array{Int64}(undef, estimated_number_of_data_points)
    itOfCheck[1] = 0

    elapsedTime = Array{Float64}(undef, estimated_number_of_data_points)
    elapsedTime[1] = 0.0

    nodesInGraph = Array{Int64}(undef, estimated_number_of_data_points)
    nodesInGraph[1] = 1

    costOfGoal = Array{Float64}(undef, estimated_number_of_data_points)
    costOfGoal[1] = Inf

    #numReduces = Array{Int64}(undef, estimated_number_of_data_points)
    #numReduces[1] = 0
  else
    savingStats = false
  end

  # while planning time left, plan. (will break out when done)
  robot_slice_start = time_ns()
  S.startTimeNs = robot_slice_start
  S.elapsedTime = 0.0

  oldrrtLMC = Inf

timeForGC = 0

  # while v_bot != v_goal
  while (true)
    hyberBallRad = min(delta, ballConstant*((log(1+KD.treeSize)/(KD.treeSize))^(1/S.d)))
    itOfCheck[checkPtr] += 1
    now_time = time_ns()

    # calculate the end time of the first slice
    slice_end_time = (1+sliceCounter)*slice_time

    # see if warmup time has ended
    warmUpTimeJustEnded = false
    if S.inWarmupTime && S.warmupTime < S.elapsedTime
      warmUpTimeJustEnded = true
      S.inWarmupTime = false
    end

    ### add/remove newly "detected" obstacles ###

    ### beginning of remove obstacle

    # remove obstacles at the required time
    S.elapsedTime = (time_ns() - S.startTimeNs)/1000000000 - save_elapsed_time
    list_item = S.obstacles.front
    removedObstacle = false
    while list_item != list_item.child
      ob = list_item.data

      if !ob.senseableObstacle && !ob.obstacleUnused && (ob.startTime + ob.lifeSpan <= S.elapsedTime)
        # time to remove obstacle
        removeObstacle(S, KD, Q, ob, root, hyberBallRad, S.elapsedTime, S.moveGoal)
        removedObstacle = true
      elseif ob.senseableObstacle && ob.obstacleUnusedAfterSense && Wdist(R.robotPose, ob.position) < robotSensorRange + ob.radius
        # place to remove obstacle

        # because the space that used to be in this obstacle was never sampled
        # there will be a hole in the graph where it used to be. The following
        # attempts to mitigate this problem by requiring that the next few samples
        # come from the space that used to be inside the obstacle
        randomSampleObs(S, KD, ob) # stores samples in the sample stack
        removeObstacle(S, KD, Q, ob, root, hyberBallRad, S.elapsedTime, S.moveGoal)
        ob.senseableObstacle = false
        ob.startTime = Inf
        removedObstacle = true
      elseif S.spaceHasTime && ob.nextDirectionChangeTime > R.robotPose[3] && ob.lastDirectionChangeTime != R.robotPose[3]
        # a moving obstacle with unknown path is changing direction, so remove
        # its old anticipated trajectory

        removeObstacle(S, KD, Q, ob, root, hyberBallRad, S.elapsedTime, S.moveGoal)
        ob.obstacleUnused = false # obstacle is still used
        removedObstacle = true
      end

      list_item = list_item.child
    end
    if removedObstacle
      println("Removed obstacle")

      reduceInconsistency(Q, S.moveGoal, robotRads, root, hyberBallRad)
    end
    ### end of remove obstacle

    ### beginning of add obstacle
    # add obstacles at the required time
    list_item = S.obstacles.front
    addedObstacle = false
    while list_item != list_item.child
      ob = list_item.data

      if !ob.senseableObstacle && ob.obstacleUnused && (ob.startTime <= S.elapsedTime <= ob.startTime + ob.lifeSpan)
        # time to add
        addNewObstacle(S, KD, Q, ob, root, vCounter, R)
        addedObstacle = true
      elseif ob.senseableObstacle && !ob.obstacleUnusedAfterSense && Wdist(R.robotPose, ob.position) < robotSensorRange + ob.radius
        # place to add obstacle
        addNewObstacle(S, KD, Q, ob, root, vCounter, R)
        ob.senseableObstacle = false
        addedObstacle = true
      elseif S.spaceHasTime && ob.nextDirectionChangeTime > R.robotPose[3] && ob.lastDirectionChangeTime != R.robotPose[3]
        # time that a moving obstacle with unknown path changes direction
        ob.obstacleUnused = false
        changeObstacleDirection(S, ob, R.robotPose[3])
        addNewObstacle(S, KD, Q, ob, root, vCounter, R)
        ob.lastDirectionChangeTime = copy(R.robotPose[3])
        #println("$(ob.nextDirectionChangeTime)  $(S.moveGoal.position[3]) ")
        addedObstacle = true
      elseif warmUpTimeJustEnded && !ob.obstacleUnused
        # warm up time is over, so we need to treat all active obstacles
        # as if they have just been added
        addNewObstacle(S, KD, Q, ob, root, vCounter, R)
        addedObstacle = true
      end

      list_item = list_item.child
    end
    if addedObstacle

      # propogate inf cost to all nodes beyond the obstacle and in its
      # basin of attraction
      propogateDescendants(Q, R)

      if !markedOS(S.moveGoal) # I'm pretty sure this is always true, since OS is emopty here -- M.O.
        verifyInQueue(Q, S.moveGoal)
      end

      println("added obstacle")

      reduceInconsistency(Q, S.moveGoal, robotRads, root, hyberBallRad)
    end
    ### end of add obstacle


# uncomment to explicitly save data after add/remove:
#    if addedObstacle || removedObstacle
#     ## visualize graph #############
#      if saveVideoData
#        before_save_time = time_ns()
#        saveRRTTree(KD, "temp/edges_$(vCounter).txt")
#        saveRRTNodes(KD, "temp/nodes_$(vCounter).txt")
#        #saveRRTNodesCollision(KD, "temp/Cnodes_$(vCounter).txt")
#        saveRRTPath(S, S.moveGoal , root, R, "temp/path_$(vCounter).txt")
#        saveObstacleLocations(S.obstacles, "temp/obstacles_$(vCounter).txt")
#        saveData(R.robotMovePath[1:R.numRobotMovePoints,:], "temp/robotMovePath_$(vCounter).txt")
#        vCounter += 1
#        save_elapsed_time += (time_ns()-before_save_time)/1000000000
#      end
#      ## end of visualize graph ######
#    end

    ### done with add/remove newly "detected" obstacles ###


    # if this robot has used all of its allotted planning time of this slice
    # then make robot move + visualize + save data
    S.elapsedTime = (time_ns() - S.startTimeNs)/1000000000 - save_elapsed_time
    if S.elapsedTime >= slice_end_time

      # use current waypoints to solve tpbvp's
        # find current robot path (straight line)
        #


      # find max deviation
      # augment obstacles using max deviation
      # add obstacles
      # solve tpbvp's using new waypoints
      # save robot moving path (Q) and make robots move
      # after each slice, visualize
      # -> (robot moving for slice_time - visualize - run RRTQX - visualize) ->

      # calculate the end time of the next slice
      slice_end_time = (1+sliceCounter)*slice_time



      robot_slice_start = now_time

      sliceCounter += 1

      truncElapsedTime = floor(S.elapsedTime * 1000)/1000

      println("slice $(sliceCounter) --- $(truncElapsedTime)  --------  $(S.moveGoal.rrtTreeCost) $(S.moveGoal.rrtLMC)----")

      # if saving stats
      if length(statsArgs) >= 1 && statsArgs[1]
        # record data
        elapsedTime[checkPtr] = S.elapsedTime
      end

      ## move robot if the robot is allowed to move, otherwise planning is finished
      # so break out of the control loop
      if elapsedTime[checkPtr] > total_planning_time + slice_time
        if MoveRobotFlag
          moveRobot(S, Q, KD, slice_time, root, hyberBallRad, R)
        else
          println("done (not moving robot)")
          break
        end
      end

      if searchType == "RRT#" || searchType == "RRTx"
        reduceInconsistency(Q, S.moveGoal, robotRads, root, hyberBallRad)
        if(S.moveGoal.rrtLMC != oldrrtLMC)
          #printPathLengths(moveGoal)
          oldrrtLMC = S.moveGoal.rrtLMC
        end
      end

      ## visualize graph #############
      ## also run & visualize Q ##
      if saveVideoData

        before_save_time = time_ns()

        saveRRTTree(KD, "temp/edges_$(vCounter).txt")
        saveRRTNodes(KD, "temp/nodes_$(vCounter).txt")
        #saveRRTNodesCollision(KD, "temp/Cnodes_$(vCounter).txt")
        saveRRTPath(S, S.moveGoal, root, R, "temp/path_$(vCounter).txt")
        saveObstacleLocations(S.obstacles, "temp/obstacles_$(vCounter).txt")
        saveData(R.robotMovePath[1:R.numRobotMovePoints,:], "temp/robotMovePath_$(vCounter).txt")

        vCounter += 1
        S.fileCtr = vCounter

        save_elapsed_time += (time_ns()-before_save_time)/1000000000
      end
      ## end of visualize graph ######

      # check if the robot has reached its movement goal
      if R.robotPose == root.position
        break
      end

      # if saving stats
      if length(statsArgs) >= 1 && statsArgs[1]
        # update statistics about run, assuming that we are saving them

        if checkPtr < length(costOfGoal)
          checkPtr += 1
          itOfCheck[checkPtr] = itOfCheck[checkPtr-1] + 1

          nodesInGraph[checkPtr] = KD.treeSize
          costOfGoal[checkPtr] = min(goal.rrtTreeCost, goal.rrtLMC)
          #costOfGoal[checkPtr] = extractPathLength(goal , root)
          #numReduces[checkPtr] = Q.numReduces
        else
          #println("WARNING: out of space to save stats")
        end
      end
    end # if this robot has used all of its allotted planning time of this slice

    #### END of obstacle and robot pose update
    #### START of normal graph search stuff


    # pick a random node
    newNode = S.randNode(S)

    if newNode.kdInTree # happens when we explicitly sample the goal every so often
      # nodes will be updated automatically as information gets propogated to it
      continue
    end


    # find closest old node to the new node
    (closestNode, closestDist) = kdFindNearest(KD, newNode.position)

    # saturate
    #if closestDist > delta && newNode != S.goalNode
    #  newNode.position = closestNode.position  + (newNode.position - closestNode.position)*delta/closestDist
    #end

    if closestDist > delta && newNode != S.goalNode
      saturate(newNode.position, closestNode.position, delta)
    end



    # check for collisions vs static obstacles
    (explicitlyUnSafe, retCert) = explicitNodeCheck(S, newNode)

    if explicitlyUnSafe
      continue
    end

    #!!! Need to look into this
    GC.enable(false)

    # extend
    extend(S, KD, Q, newNode, closestNode, delta, hyberBallRad, S.moveGoal)



    # make graph consistant (RRT# and RRTx)
    if searchType == "RRT#" || searchType == "RRTx"
      reduceInconsistency(Q, S.moveGoal, robotRads, root, hyberBallRad)
      if(S.moveGoal.rrtLMC != oldrrtLMC)
        #printPathLengths(S.moveGoal)
        oldrrtLMC = S.moveGoal.rrtLMC
      end
    end

    GC.enable(true)
  end
  ## end of while(true)

  elapsedTime[checkPtr] = (time_ns()-startTime)/1000000000

  if (length(statsArgs) >= 1 && statsArgs[1])
    if (!goal.rrtParentUsed)
      print("goal has no parent ")
    end

    stats = hcat(elapsedTime, itOfCheck, nodesInGraph, costOfGoal)

    saveData(stats[1:checkPtr,:], dataFileName)

    #reduceData = [numReduces', nodesInGraph']'
    #saveData(reduceData[1:checkPtr,:], "temp/reduceStats.txt")
  end

  moveLength = sum(sqrt, sum((R.robotMovePath[1:R.numRobotMovePoints-1, :] - R.robotMovePath[2:R.numRobotMovePoints, :]).^2, dims=2))

  println("distance traveled by robot: $(moveLength[1])")
  return true
end



################################################################################
############################ end of file #######################################
################################################################################

# saves obstacle data to a file
function saveOriginalObstacleLocations_Q(obstacles::List{Obstacle}, fileName)

  fptr = open(fileName, "w")

  listNode = obstacles.front
  while listNode != listNode.child
    ob = listNode.data

    if ob.kind == 6 || ob.kind == 7
      # for time obstacles we save data about path [x, y, time, radius]
      # then we put NaNs rows between obstacls

      for i = 1:size(ob.path,1)
        writedlm(fptr, [(reshape(ob.path[i,:], 1, length(ob.path[i,:])) + [ob.position 0.0]) ob.radius], ',')
      end
      writedlm(fptr, [NaN NaN NaN NaN], ',')


      listNode = listNode.child
      continue
    end

    if ob.kind != 3 && ob.kind != 4  && ob.kind != 5
      println("warning cannot save non-polygon obstacle to file (not implimented) ")
      listNode = listNode.child
      continue
    end

    if ob.obstacleUnused
      listNode = listNode.child
      continue
    end

    for i = 1:size(ob.polygonWithoutAug,1)
      writedlm(fptr, reshape(ob.polygonWithoutAug[i,:], 1, length(ob.polygonWithoutAug[i,:])), ',')
    end
    writedlm(fptr, reshape(ob.polygonWithoutAug[1,:], 1, length(ob.polygonWithoutAug[1, :])), ',')
    writedlm(fptr, [NaN NaN], ',')
    listNode = listNode.child
  end

  close(fptr)

end

# saves the RRT nodes to a file
# it uses the KD-Tree structure to traverse the tree
function saveKds_Q(S::CSpace, fileName)
  fptr = open(fileName, "w")
  writedlm(fptr, S.kino_dist)
  writedlm(fptr, S.augDist)
  close(fptr)
end

function saveRRTPath_Q(S::TS, node::T, root::T, robot::RobotData{T}, fileName) where {T, TS}
  fptr = open(fileName, "w")
  thisNode = node
  i = 0 # added to detect for inf paths
  while thisNode != root && thisNode.rrtParentUsed && i < 1000
    saveEdgeTrajectory(fptr, thisNode.rrtParentEdge)
    thisNode = thisNode.rrtParentEdge.endNode
    i += 1 # added to detect for inf paths
  end
  writedlm(fptr, hcat(thisNode.position[1,1:2]'), ',')
  close(fptr)
end
