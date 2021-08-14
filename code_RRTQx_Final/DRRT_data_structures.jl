function makeMyTextEditorDisplayNice0() # thanks!
  return
end

################################# data structures ################################
### Note: other things like heap, kdtree, etc. are loaded from files.          ###
##################################################################################


### !!!!!                                                              !!!!!
### TYPEALIASING SHOULD BE USED TO MAKE Edge{T} BE WHATEVER EDGE IS USED.
### This has been moved to be external to this file, but still is required
### !!!!!                                                              !!!!!
# typealias Edge{T} DubinsEdge{T}
# typealias Edge{T} SimpleEdge{T}


### the node that is used to build the search graph, it will also be inserted
### into the KD tree to help find nearest neighbors, and various queues and
### lists depending on the algorithm that is being used.
### see above for info regarding edges
mutable struct RRTNode{T}

  # data used for KD Tree
  kdInTree::Bool           # set to true if this node is in the kd-tree
  kdParentExist::Bool      # set to true if parent in the tree is used
  kdChildLExist::Bool      # set to true if left child in the tree is used
  kdChildRExist::Bool      # set to true if right child in the tree is used

  # data used for heap in KNN-search
  heapIndex::Int
  inHeap::Bool

  # data used for RRT
  rrtParentUsed::Bool                        # flag for if this node has a parent

  # data used for RRT#
  rrtNeighborsOut::JList{Edge{RRTNode{T}}}   # edges in the graph that can be reached
                                             # from this node. In + this holds such
                                             # nodes in the current D-ball

  rrtNeighborsIn::JList{Edge{RRTNode{T}}}    # edges in the graph that reach this
                                             # node. In + this holds such nodes in the
                                             # current D-ball


  priorityQueueIndex::Int                    # index in the queue
  inPriorityQueue::Bool                      # flag for in the queue

  # data used for RRTx (my idea)
  SuccessorList::JList{Edge{RRTNode{T}}}     # edges to nodes that use this node as
                                             #their parent

  InitialNeighborListOut::JList{Edge{RRTNode{T}}}  # edges to nodes in the original
                                                   # ball that can be reached
                                                   # from this node.

  InitialNeighborListIn::JList{Edge{RRTNode{T}}}   # edges to nodes in the original
                                                   # ball that can reach this node.


  inOSQueue::Bool                            # flag for in the OS queue
  isMoveGoal::Bool                           # true if this is move goal (robot pose)

  # position (used for everything)
  position::Array{T}          # a dX1 array where d is the dimesnions of the space

  # more data used for KD Tree
  kdSplit::Int                # the dimension used for splitting at this node
  kdParent::RRTNode{T}        # parent in the tree
  kdChildL::RRTNode{T}        # left child in the tree
  kdChildR::RRTNode{T}        # right child in the tree

  # more data used for RRT
  rrtParentEdge::Edge{RRTNode{T}} # edge to the node that is this node's parent

  # # QEdge
  # qParentEdge::QEdge{RRTNode{T}} # Q edge to the node that is this node's parent
  #                                # the startNode and endNode are very near to this
  #                                # node and this node's parent (except for the robot
  #                                # initial position, that would be the same)

  # more data used for RRT*
  rrtTreeCost::Float64        # the cost to get to the root through the tree

  # more data used for RRT#
  rrtLMC::Float64             # locally minimum cost (1-step look ahead)
  rrtH::Float64               # the heuristic estimate of the distance to the goal !!!!!
  tempEdge::Edge{RRTNode{T}}  # this is a temporary storage location to avoid
                              # calculating the same trajectory multiple times

  # more data used for RRTx
  successorListItemInParent::JListNode{Edge{RRTNode{T}}} # pointer to the list node in
                                                         # the parent's successor list
                                                         # that holds parent's edge to
                                                         # this node

  ## data for certificate based collision detection (currently unused, but could be
  ## added in future without much difficulty
  #hasCertificate::Bool       # true if this node was explicitly checked
  #certificateValue::Float64  # value of certificate (if hasCertificate) 1 per robot
  #certifyingNode::RRTNode{T}   # the node that certified (if !hasCertificate)

  # constructors
  function RRTNode{T}() where {T}
    rrtNeighborsOut = JList{Edge{RRTNode{T}}}()
    rrtNeighborsIn = JList{Edge{RRTNode{T}}}()

    SuccessorList = JList{Edge{RRTNode{T}}}()
    InitialNeighborListOut = JList{Edge{RRTNode{T}}}()
    InitialNeighborListIn = JList{Edge{RRTNode{T}}}()

    new{T}(false, false, false, false, -1, false, false, rrtNeighborsOut,
    rrtNeighborsIn, -1, false, SuccessorList, InitialNeighborListOut,
    InitialNeighborListIn, false, false)
  end

  function RRTNode{T}(A::Array{Float64}) where {T}
    rrtNeighborsOut = JList{Edge{RRTNode{T}}}()
    rrtNeighborsIn = JList{Edge{RRTNode{T}}}()

    SuccessorList = JList{Edge{RRTNode{T}}}()
    InitialNeighborListOut = JList{Edge{RRTNode{T}}}()
    InitialNeighborListIn = JList{Edge{RRTNode{T}}}()

    new{T}(false, false, false, false, -1, false, false, rrtNeighborsOut,
    rrtNeighborsIn, -1, false, SuccessorList, InitialNeighborListOut,
    InitialNeighborListIn, false, false, A)
  end
end


### This holds an obstacle, there are many types of obstacles each with their
### own type of behaviour
mutable struct Obstacle
  kind::Int    # 1 = ball
               # 2 = axis aligned hyperrectangle,
               # 3 = polygon
               # 4 = polygon with safe direction
               # 5 = high dimensional prismatic polygon
               # 6 = polygon that moves in time along a predefined path
               # 7 = similar to 6 but the robot does not "know" obstacle path a priori

  startTime::Float64  # obstacle appears this long after the start of the run
                      # 0 by default

  lifeSpan::Float64   # the lifespan of the obstacle (defaults to Inf)
  obstacleUnused::Bool # if true, then this obstacle will not be checked


  senseableObstacle::Bool # true if this obstacle can be sensed by the robot
                          # i.e., may change state after the robot gets close enough
                          # is set to false after sensing by robot happens, the
                          # distance at which sensing occours is set in the part of
                          # the code pertaining to the robot (and not here)
                          # default is false

  obstacleUnusedAfterSense::Bool # obstacleUnused is set to this value after the
                                 # robot senses this obstacle, default is true


  position::Array{Float64} # initial position of obstacle

  # data for D-dimensional ball obstacles (kind = 1) and as bound on obstacle
  # (all kinds)
  radius::Float64
  radiusWithoutAug::Float64 # save radius before aug -------- rrtqx
  # data for an axis aligned D-dimensional hyperrectangle obstacle (kind = 2)
  span::Array{Float64} # distance away from the center that this obstacle spans

  # data used for a polygon (kind = 3,6)
  polygon::Array{Float64} # each row holds a point of the polygon

  polygonWithoutAug::Array{Float64} # save polygon before aug -------- rrtqx

  # direction is used for directional polygon (kind = 4)
  direction::Char

  # used only for obstacle kind 5
  PrismSpanMin::Array{Float64}
  PrismSpanMax::Array{Float64}

  # stuff for time obstacles (kind = 6,7)
  velocity::Float64       # speed that the obstacle moves at
  path::Array{Float64, 2} # were p[i,:] is the i-th point of an offset path this
                          # obstacle follows, note that p[1,:] should be [0, 0]
                          # this path the robot knows about

  # the following two field is used to keep track of original values, values
  # used for collision checking are calculated as needed vs. the relivant time
  # projection by adding the current offset path value to the following
  originalPolygon::Array{Float64}


  # stuff for time obstacles with "unknown" path (kind = 7)
  unknownPath::Array{Float64, 2} # this is the path the obstacle will actually follow
                                 # but that the robot does not know about. path,
                                 # which is what is used for collision checking,
                                 # is re-calculated from unknownPath each time the
                                 # obstacle changes direction

  nextDirectionChangeTime::Float64 # this is the next time that a direction change
                                   # of this obs will occour (unknown to robot).
                                   # note that time count's down to 0, so this is
                                   # generally less than the current time associated
                                   # with robot's state as it moves

  nextDirectionChangeInd::Int64     # unknownPath index of nextDirectionChangeTime
  lastDirectionChangeTime::Float64  # unknownPath time of last change

  # constructors:
  # ball:
  function Obstacle(kind::Int, position::Array{Float64}, radius::Float64)
    new(kind, 0.0, Inf, false, false, true, position, radius)
  end

  # hyperrectangle:
  function Obstacle(kind::Int, position::Array{Float64}, span::Array{Float64})
    # calculate position of center and distance of center to furtherest point
    # this is used to do quick obstacle checks

    radius = sqrt(sum(span.^2))
    new(kind, 0.0, Inf, false, false, true, position, radius, span)
  end



  # polygon:
  function Obstacle(kind::Int, polygon::Array{Float64})
    # calculate position of center and distance of center to furtherest point
    # this is used to do quick obstacle checks

    polygonWithoutAug = copy(polygon) # -------- rrtqx

    position = [(maximum(polygon[:,1]) + minimum(polygon[:,1]))/2.0 (maximum(polygon[:,2]) + minimum(polygon[:,2]))/2.0]

    radius = sqrt(maximum(sum((broadcast(-, polygon, position)).^2, dims=2)))
    radiusWithoutAug = copy(radius) # -------- rrtqx
    span = [-1.0 -1.0] # just a dummy (could combine span and polygon)
    new(kind, 0.0, Inf, false, false, true, position, radius, radiusWithoutAug, span, polygon, polygonWithoutAug) # ------ rrtqx
  end


  # polygon with direction:
  # (N, S, E, W), e.g., if N, then movement is allowed as long as dy is positive
  function Obstacle(kind::Int, polygon::Array{Float64}, direction::Char)
    # calculate position of center and distance of center to furtherest point
    # this is used to do quick obstacle checks

    position = [(maximum(polygon[:,1]) + minimum(polygon[:,1]))/2.0 (maximum(polygon[:,2]) + minimum(polygon[:,2]))/2.0]
    radius = sqrt(max(sum((broadcast(-, polygon, position)).^2,2)))
    span = [-1.0 -1.0] # just a dummy (could combine span and polygon)
    new(kind, 0.0, Inf, false, false, true, position, radius, span, polygon, direction)
  end


  # High-D prismatic polygon:
  function Obstacle(kind::Int, polygon::Array{Float64}, PrismSpanMin::Array{Float64}, PrismSpanMax::Array{Float64})
    error("probably not used any more")
    Ob = Obstacle(kind, polygon)
    Ob.PrismSpanMin = copy(PrismSpanMin)
    Ob.PrismSpanMax = copy(PrismSpanMax)
    return Ob
  end
end




### the configuration space. This holds data about the start, goal, obstacles, etc.
### and is used for sampling new points.
mutable struct CSpace{T}
  d::Int                          # dimensions
  obstacles::List{Obstacle}       # a list of obstacles
  # augObs::List{Obstacle}          # a list of augmented obstacles for visualization ----- rrtqx
  obsDelta::Float64               # the granularity of obstacle checks on edges
  lowerBounds::Array{Float64}     # 1XD array containing the lower bounds
  upperBounds::Array{Float64}     # 1XD array containing the upper bounds
  width::Array{Float64}           # 1XD array containing upper_bounds-lower_bounds
  start::Array{Float64}           # 1XD array containing the start location
  goal::Array{Float64}            # 1XD array containing the goal location

  # flags that indicate what type of search space we are using (these are mosly here
  # to reduce the ammoung of duplicate code for similar spaces, althouth they should
  # probably one day be replaced with a different approach that takes advantage of
  # Julia's multiple dispatch and polymophism)

  spaceHasTime::Bool              # if true then the 3rd dimension of the space is
                                  # time
  spaceHasTheta::Bool             # if true then the 4th dimension of the space is
                                  # theta, in particular a dubins system is used

  # stuff for sampling functions
  pGoal::Float64                  # the probability that the goal is sampled
  randNode::Function              # the sampling function to use (takes a Cspace)
  goalNode::RRTNode{T}            # the goal node
  root::RRTNode{T}                # the root node
  moveGoal::RRTNode{T}            # the current movegoal (robot position) node
  itsUntilSample::Int             # a count down to sample a particular point
  itsSamplePoint::Array{Float64}  # sample this when itsUntilSample == 0

  timeSamplePoint::Array{Float64} # sample this when waitTime has passed
  waitTime::Float64               # time to wait in seconds
  startTimeNs::UInt64             # time this started
  elapsedTime::Float64            # elapsed time since start (where time spent saving
                                  # experimental data has been removed)
  obstaclesToRemove::Obstacle     # an obstacle to remove

  robotRadius::Float64            # robot radius
  robotVelocity::Float64          # robot velocity (used for dubins w/o time)

  dubinsMinVelocity::Float64     # min velocity of dubins car (for dubins + time)
  dubinsMaxVelocity::Float64     # max velocity of dubins car (for dubins + time)


  sampleStack::JList{Array{Float64,2}} # points to sample in the future

  hypervolume::Float64            # hypervolume of space

  delta::Float64                  # RRT parameter delta

  minTurningRadius::Float64       # min truning radius, e.g., for dubins car

  fileCtr::Int64                  # file counter, used for debugging only

  warmupTime::Float64             # the ammount of warm up time allowed (obstacles are
                                  # are ignored for warm up time)
  inWarmupTime::Bool              # true if we are in the warmup time

  numCoveredLocal::Int64          # for Q visualization
	numLocal::Int64
  numEsqTrigLocal::Int64
  kino_dist::Float64
  localKd::Float64
  augDist::Float64
  NormEsqvec::Array{Float64}
  TrigCondvec::Array{Float64}
  lastVelocity::Array{Float64}
  # constructors

  function CSpace{T}(D::Int, ObsDelta, L, U, S, G) where {T}
    O = List{Obstacle}()
    # aO = List{Obstacle}() ------ rrtqx
    CS = new{T}(D, O, ObsDelta, L, U, U-L, S, G, false, false)
    # CS = new{T}(D, O, ObsDelta, L, U, U-L, S, G, false, false)
    CS.hypervolume = 0.0 # flag indicating that this needs to be calculated

    CS.inWarmupTime = false
    CS.warmupTime = 0.0 # default value for time for build graph with no obstacles

    return CS
  end
end

## queue data structure used for RRT, basically empty, used to make coding easier
mutable struct rrtQueue{T}
  rrtQueue{T}() where {T} = new{T}()
end


## queue data structure used for RRT*, basically empty, used to make coding easier
mutable struct rrtStarQueue{T}
  rrtStarQueue{T}() where {T} = new{T}()
end

## queue data structure used for RRT#
mutable struct rrtSharpQueue{A,B}
  Q::BinaryHeap{A,B}  # normal queue (sorted based on cost from goal)
  S::CSpace


  #numReduces::Int64 # for testing, take out when done

  rrtSharpQueue{A, B}() where {A, B} = new{A, B}()
end

## queue data structure used for RRTx
mutable struct rrtXQueue{A,B}
  Q::BinaryHeap{A,B}  # normal queue (sorted based on cost from goal)
  OS::JList{A}        # obstacle successor stack

  S::CSpace
  changeThresh::Float64  # the threshold of local changes that we care about

  rrtXQueue{A, B}() where {A, B} = new{A, B}()
end


# this is used to make iteration through a particualr node's neighbor edges easier
# given that each node stores all of its neighbor edges in three different places
mutable struct RRTNodeNeighborIterator{T, TE}
  thisNode::T                # the node who's neighbors we are iterating through

  listFlag::Int64            # flag with the following values:
                             #   0: uninitialized
                             #   1: successors
                             #   2: original neighbors
                             #   3: current neighbors

  listItem::JListNode{TE}    # a pointer to the position in the current neighbor list
                             # we are iterating through

  # constructors
  RRTNodeNeighborIterator{T, TE}(node::T) where {T, TE} = new{T, TE}(node::T, 0)
end



# This holds the stuff associated with the robot that is necessary for movement.
# Although some of the fields are used primarily for simulation of robot movement,
# currentMoveInvalid is important for the algorithm in general.
mutable struct RobotData{T}
  robotPose::Array{Float64,2}       # this is where the robot is (i.e., where it was
                                    # at the end of the last control loop)

  nextRobotPose::Array{Float64,2}   # this is where the robot will be at the end
                                    # of the current control loop

  nextMoveTarget::T                 # this is the node at the root-end of the edge
                                    # contains nextRobotPose

  distanceFromNextRobotPoseToNextMoveTarget::Float64
                                    # this holds the distance from nextRobotPose to
                                    # nextMoveTarget along the trajectory the robot
                                    # will be following at that time

  moving::Bool                      # set to true when the robot starts moving

  currentMoveInvalid::Bool          # this gets set to true if nextMoveTarget
                                    # has become invalid due to dynamic obstacles

  robotMovePath::Array{Float64,2}   # this holds the path the robot has followed
                                    # from the start of movement up through robotPose

  numRobotMovePoints::Int64       # the number of points in robotMovePath

  robotLocalPath::Array{Float64,2}  # this holds the path between robotPose and
                                    # nextRobotPose (not including the former)

  numLocalMovePoints::Int64       # the number of points in robotLocalPath

  robotEdge::Edge{T}                # this is the edge that contians the
                                    # trajectory that the robot is currently
                                    # following
  robotEdgeUsed::Bool               # true if the robotEdge is populated

  # note that currently only one of the two following paramiters is used at a time,
  # which one is used depends on if time is explicitly part of the state space

  distAlongRobotEdge::Float64       # the current distance that the robot "will be"
                                    # along robotEdge (i.e., next time slice)

  timeAlongRobotEdge::Float64       # the current time that the robot "will be"
                                    # along robotEdge (i.e., next time slice)

  # the following things are only used to help with visualization
  robotEdgeForPlottingUsed::Bool
  robotEdgeForPlotting::Edge{T}
  distAlongRobotEdgeForPlotting::Float64
  timeAlongRobotEdgeForPlotting::Float64

  function RobotData{T}(rP::Array{Float64,2}, nMT::T, maxPathNodes::Int) where {T}

    R = new{T}()

    R.robotPose = rP
    R.nextRobotPose = rP
    R.nextMoveTarget = nMT
    R.distanceFromNextRobotPoseToNextMoveTarget = 0.0

    R.moving = false
    R.currentMoveInvalid = false

    R.robotMovePath = Array{Float64}(undef, maxPathNodes, length(rP))
    R.robotMovePath[1,:] = rP
    R.numRobotMovePoints = 1

    R.robotLocalPath = Array{Float64}(undef, maxPathNodes, length(rP))
    R.numLocalMovePoints

    R.robotEdgeUsed = false
    R.distAlongRobotEdge = 0.0
    R.timeAlongRobotEdge = 0.0

    R.robotEdgeForPlottingUsed = false
    R.distAlongRobotEdgeForPlotting = 0.0
    R.timeAlongRobotEdgeForPlotting = 0.0

    return R
  end
end
