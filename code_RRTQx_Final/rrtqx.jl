function rrtqx(S::TS, total_planning_time::Float64, slice_time::Float64,
  delta::Float64, ballConstant::Float64, changeThresh::Float64,
  searchType::String, MoveRobotFlag::Bool, saveVideoData::Bool,
  statsArgs...) where {TS}


  T = RRTNode{Float64}

  # NOTE THIS IS HARD CODED HERE (SHOULD PROBABLY MAKE INPUT ARGUMENT)
  robotSensorRange = 20 # 20 # used for "sensing" obstacles

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

  KD = KDTree{RRTNode{Float64}}(S.d, KDdist)

  # init queue. Note that for algorithms that do not require a queue this
  # is essentially an empty variable that is used to help with polymorphism,
  # which makes life easy, the code fast, and allows much of the same code
  # to be used for all algorithms. Note about RRTx: I decided to forgoe RRT*+
  # and focus only on #+ sice the latter will likely be more useful than the
  # former in practice

  if searchType == "RRTx"
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

  vCounter = 0 # helps with visuilizing data
  S.fileCtr = vCounter

  sliceCounter = 0 # helps with saving accurate time data

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

  # firstTimeQLearning = true # to check the first time to get Q_Path. default to be true

  # qPathHasCollision = false # check collision. default to be false

  S.augDist = 0.0             # initialize augDist
  S.kino_dist = 0.0

  localPoseAndKd = Array{Tuple{Array{Float64,2},Float64}}(undef,1000)
  localNormEsq = Array{Float64}(undef,3000)
  localTrigCond = Array{Float64}(undef,3000)
  S.numCoveredLocal = 0
  S.numLocal = 0
  # S.numErrTrigCoveredLocal = 0
  # S.numErrTrigLocal = 0
  ### initialize learning process


  ### end of initializing learning process

  # NormEsqvec = Array{Float64}(undef, 0)
  # TrigCondvec = Array{Float64}(undef, 0)
  S.NormEsqvec = zeros(0,)
  S.TrigCondvec = zeros(0,)

  S.lastVelocity = [0.; 0.]


  # environmentChangeFinished = false
  # while v_bot != v_goal

  while true
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

	### deal with larger kino_dist

	if S.kino_dist > S.augDist

    # not count in time for Minkowski sum
    before_save_time = time_ns()

    println("--------------------------------------------------------------------- New KD: $(S.kino_dist)")

	  S.augDist = S.kino_dist
	  
	  obstacleAugmentation(S, S.augDist)

	  # add latest augmented obs
	  list_item = S.obstacles.front
	  while list_item != list_item.child
		ob = list_item.data
		if !ob.obstacleUnused
		  addNewObstacle(S, KD, Q, ob, root, vCounter, R)
	    end
		list_item = list_item.child
	  end
	  # println("-------------------------------------------------------------------------------R.currentMoveInvalid = $(R.currentMoveInvalid)")
	  propogateDescendants(Q, R)
	  # println("-------------------------------------------------------------------------------R.currentMoveInvalid = $(R.currentMoveInvalid)")

	  if R.currentMoveInvalid && explicitPointCheck(S, R.robotPose)[1] # R in aug obs && R not in ori obs
		R.currentMoveInvalid = false
		println("-------------------------------------------------------------------------------R.currentMoveInvalid = $(R.currentMoveInvalid) ----------- robot in aug obs")
	  end

	  if !markedOS(S.moveGoal) # I'm pretty sure this is always true, since OS is emopty here -- M.O.
		verifyInQueue(Q, S.moveGoal)
	  end
	  reduceInconsistency(Q, S.moveGoal, robotRads, root, hyberBallRad)

	  save_elapsed_time += (time_ns()-before_save_time)/1000000000

	  # println("obstacle enlarged")
	end

	### end of deal with larger kino_dist


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

	# if S.elapsedTime >= 13.0 && !environmentChangeFinished
	# 	ob = S.obstacles.front.data
	# 	randomSampleObs(S, KD, ob)
	# 	removeObstacle(S, KD, Q, ob, root, hyberBallRad, S.elapsedTime, S.moveGoal)
	# 	removedObstacle = true
	# 	environmentChangeFinished = true
	# end

	if removedObstacle
      println("----------------------------------------------------------------------------- Removed obstacle")
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
	  # println("-------------------------------------------------------------------------------R.currentMoveInvalid = $(R.currentMoveInvalid)")
	  propogateDescendants(Q, R)
	  # println("-------------------------------------------------------------------------------R.currentMoveInvalid = $(R.currentMoveInvalid)")
      if !markedOS(S.moveGoal) # I'm pretty sure this is always true, since OS is emopty here -- M.O.
        verifyInQueue(Q, S.moveGoal)
      end
      println("--------------------------------------------------------------------------------- Added obstacle")
      reduceInconsistency(Q, S.moveGoal, robotRads, root, hyberBallRad)
    end
    ### end of add obstacle

# 	if (S.kino_dist > augDist) || removedObstacle || addedObstacle
#
# 	  obstacleAugmented = false
#
# 	  if S.kino_dist > augDist
# 		println("kd increased")
# 		augDist = S.kino_dist
# 	  end
#
# 	  # pop augmented obs if needed
# 	  listEmpty(S.augObs)
# 	  # now S.obstacles only contains original obs, S.augObs is empty
# ##################################################################################################
# 	  # push augmented obs into S.obstacles and S.augObs
# 	  obstacleAugmentation(S, augDist)
#
# 	  list_item = S.augObs.front
# 	  while list_item != list_item.child
# 		  ob = list_item.data
# 		  addNewObstacle(S, KD, Q, ob, root, vCounter, R)
# 		  list_item = list_item.child
# 	  end
# 	  propogateDescendants(Q, R)
# 	  if !markedOS(S.moveGoal) # I'm pretty sure this is always true, since OS is emopty here -- M.O.
# 		  verifyInQueue(Q, S.moveGoal)
# 	  end
# 	  reduceInconsistency(Q, S.moveGoal, robotRads, root, hyberBallRad)
# 	  obstacleAugmented = true
# 	  println("added augmented obstacles")
#
#     end
    ### done with add/remove newly "detected" obstacles ###




    # if this robot has used all of its allotted planning time of this slice
    S.elapsedTime = (time_ns() - S.startTimeNs)/1000000000 - save_elapsed_time
    if S.elapsedTime >= slice_end_time

      # calculate the end time of the next slice
      slice_end_time = (1+sliceCounter)*slice_time

      robot_slice_start = now_time

      sliceCounter += 1

      truncElapsedTime = floor(S.elapsedTime * 1000)/1000

      println("slice $(sliceCounter) --- $(truncElapsedTime) -------- $(S.moveGoal.rrtTreeCost) $(S.moveGoal.rrtLMC) ----")

      # if saving stats
      if length(statsArgs) >= 1 && statsArgs[1]
        # record data
        elapsedTime[checkPtr] = S.elapsedTime
      end

      ## move robot if the robot is allowed to move, otherwise planning is finished
      # so break out of the control loop
      if elapsedTime[checkPtr] > total_planning_time + slice_time
        if MoveRobotFlag
          moveRobot_Q(S, Q, KD, slice_time, root, hyberBallRad, R, localPoseAndKd, localNormEsq, localTrigCond,save_elapsed_time)#, NormEsqvec, TrigCondvec) # 2 steps, update S.kino_dist
        else
          println("done (not moving robot)")
          break
        end
      end

      if searchType == "RRT#" || searchType == "RRTx"
        reduceInconsistency(Q, S.moveGoal, robotRads, root, hyberBallRad)
        if S.moveGoal.rrtLMC != oldrrtLMC
          #printPathLengths(moveGoal)
          oldrrtLMC = S.moveGoal.rrtLMC
        end
	  end

      ## visualize graph #############
      if saveVideoData
        before_save_time = time_ns()

		# for visualization, S.obstacles only contains original obs
		# if obstacleAugmented
		# 	saveObstacleLocations(S.augObs, "temp/augObs_$(vCounter).txt")
		# 	for i = 1:S.augObs.length
		# 		listPop(S.obstacles)
		# 	end
		# end

        saveRRTTree(KD, "temp/edges_$(vCounter).txt")
        saveRRTNodes(KD, "temp/nodes_$(vCounter).txt")
        #saveRRTNodesCollision(KD, "temp/Cnodes_$(vCounter).txt")
        saveRRTPath_Q(S, S.moveGoal, root, R, "temp/path_$(vCounter).txt")
        saveObstacleLocations(S.obstacles, "temp/obstacles_$(vCounter).txt")
		saveOriginalObstacleLocations_Q(S.obstacles, "temp/originalObs_$(vCounter).txt")
        saveData(R.robotMovePath[1:R.numRobotMovePoints,:], "temp/robotMovePath_$(vCounter).txt")
		saveKds_Q(S, "temp/kd_$(vCounter).txt")

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
    end

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
      print("goal has no parent")
    end

    stats = hcat(elapsedTime, itOfCheck, nodesInGraph, costOfGoal)

    saveData(stats[1:checkPtr,:], dataFileName)

    #reduceData = [numReduces', nodesInGraph']'
    #saveData(reduceData[1:checkPtr,:], "temp/reduceStats.txt")
  end

  moveLength = sum(sqrt, sum((R.robotMovePath[1:R.numRobotMovePoints-1, :] - R.robotMovePath[2:R.numRobotMovePoints, :]).^2, dims=2))

  println("distance traveled by robot: $(moveLength[1])")
  println("KD_max: $(S.augDist)")
  return (S.NormEsqvec, S.TrigCondvec)
  # saveData(tr, "temp/Trig.txt")
  # saveData(er, "temp/Esq.txt")
end


function multirrtqx(S::Array{TS}, N::Int64, total_planning_time::Float64, slice_time::Float64,
  delta::Float64, ballConstant::Float64, changeThresh::Float64,
  searchType::String, MoveRobotFlag::Bool, saveVideoData::Bool,
  statsArgs...) where {TS}

  T = RRTNode{Float64}

  # NOTE THIS IS HARD CODED HERE (SHOULD PROBABLY MAKE INPUT ARGUMENT)
  robotSensorRange = 20 # 20 # used for "sensing" obstacles

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
  KD = []

  allDone = false

  for i = 1:N
    push!(KD, KDTree{RRTNode{Float64}}(S[i].d, KDdist))
  end
  # init queue. Note that for algorithms that do not require a queue this
  # is essentially an empty variable that is used to help with polymorphism,
  # which makes life easy, the code fast, and allows much of the same code
  # to be used for all algorithms. Note about RRTx: I decided to forgoe RRT*+
  # and focus only on #+ sice the latter will likely be more useful than the
  # former in practice

  if searchType == "RRTx"
    Q = []
    for i = 1:N
      push!(Q, rrtXQueue{RRTNode{Float64}, typeof((Float64, Float64))}())
      Q[i].Q = BinaryHeap{RRTNode{Float64}, typeof((Float64, Float64))}(keyQ, lessQ, greaterQ, markQ, unmarkQ, markedQ, setIndexQ, unsetIndexQ, getIndexQ)
      Q[i].OS = JList{RRTNode{Float64}}()
      Q[i].S = S[i]
      Q[i].changeThresh = changeThresh
    end
  else
    error("unknown search type: $(searchType)")
  end
  for i = 1:N
    S[i].sampleStack = JList{Array{Float64,2}}() # stores a stack of points that we
                                              # desire to insert in the future in
                                              # (used when an obstacle is removed
    S[i].delta = delta
  end

  robotRads = S[1].robotRadius


  # define root node in the search tree graph
  root = []
  for i = 1:N
    push!(root, RRTNode{Float64}(S[i].start))
    
    # explicit check root
    (explicitlyUnSafe, unused) = explicitNodeCheck(S[i], root[i])
    if explicitlyUnSafe
      error("root is not safe")
    end
  

    root[i].rrtTreeCost = 0.0
    root[i].rrtLMC = 0.0

    # insert the root into the KDtree
    kdInsert(KD[i], root[i])
  end

  goal = []

  for i = 1:N
    # define a goal node
    push!(goal, RRTNode{Float64}(S[i].goal))
    goal[i].rrtTreeCost = Inf
    goal[i].rrtLMC = Inf
    S[i].goalNode = goal[i]
    S[i].root = root[i]

    S[i].moveGoal = goal[i] # this will store a node at least as far from the root as the robot
                    # during movement it key is used to limit propogation beyond the
                    # region we care about

    S[i].moveGoal.isMoveGoal = true
  end

  # paramiters that have to do with the robot path following simulation
  R = []
  for i = 1:N
    push!(R, RobotData{RRTNode{Float64}}(copy(S[i].goal), goal[i], 20000))
  end
  
  vCounter = []
  for i = 1:N
    push!(vCounter, 0) # helps with visuilizing data
    S[i].fileCtr = vCounter[i]
  end

  sliceCounter = 0 # helps with saving accurate time data

  ### end of initialization stuff


  # if saving stats about run, then allocate memory to store data
  if (length(statsArgs) >= 1 && statsArgs[1])

    savingStats = true
    estimated_number_of_data_points = 4*Int(ceil(total_planning_time/slice_time))
    checkPtr = []
    for i = 1:N
      push!(checkPtr, 1)
    end      # array location to save stats
    itOfCheck = []
    elapsedTime = []
    nodesInGraph = []
    costOfGoal = []
    for i = 1:N
      temp = Array{Int64}(undef, estimated_number_of_data_points)
      temp[1] = 0
      push!(itOfCheck, temp)

      temp = Array{Float64}(undef, estimated_number_of_data_points)
      temp[1] = 0.0
      push!(elapsedTime, temp)

      temp = Array{Int64}(undef, estimated_number_of_data_points)
      temp[1] = 1
      push!(nodesInGraph, temp)

      temp = Array{Float64}(undef, estimated_number_of_data_points)
      temp[1] = Inf
      push!(costOfGoal, temp)
    end

    #numReduces = Array{Int64}(undef, estimated_number_of_data_points)
    #numReduces[1] = 0
  else
    savingStats = false
  end

  # while planning time left, plan. (will break out when done)
  oldrrtLMC = []
  timeForGC = []
  robot_slice_start = time_ns()
  for i = 1:N
    S[i].startTimeNs = robot_slice_start
    S[i].elapsedTime = 0.0
    push!(oldrrtLMC, Inf)
    push!(timeForGC, 0)
  end

  # firstTimeQLearning = true # to check the first time to get Q_Path. default to be true

  # qPathHasCollision = false # check collision. default to be false
  for i = 1:N
    S[i].augDist = 0.0             # initialize augDist
    S[i].kino_dist = 0.0
  end

  localPoseAndKd = []
  localNormEsq = []
  localTrigCond = []
  for i = 1:N
    push!(localPoseAndKd, Array{Tuple{Array{Float64,2},Float64}}(undef,1000))
    push!(localNormEsq, Array{Float64}(undef,3000))
    push!(localTrigCond, Array{Float64}(undef,3000))
    S[i].numCoveredLocal = 0
    S[i].numLocal = 0
    S[i].numEsqTrigLocal = 0
  end
  # S.numErrTrigCoveredLocal = 0
  # S.numErrTrigLocal = 0
  ### initialize learning process


  ### end of initializing learning process

  # NormEsqvec = Array{Float64}(undef, 0)
  # TrigCondvec = Array{Float64}(undef, 0)
  lastVel = []

  for i = 1:N
    S[i].NormEsqvec = zeros(0,)
    S[i].TrigCondvec = zeros(0,)

    S[i].lastVelocity = [0.; 0.]
    push!(lastVel, [0.; 0.])
  end

  # environmentChangeFinished = false
  # while v_bot != v_goal
  hyberBallRad = []
  for i = 1:N
    push!(hyberBallRad, delta)
  end

  currPos = []
  prevPos = []
  for i = 1:N
    push!(currPos, R[i].robotPose)
    push!(prevPos, Array{Array{Float64,2}}(undef, 5))
    for j = 1:5
      prevPos[i][j] = currPos[i]
    end
  end

  
  currObsPos = []
  nextObsPos = []
  nextObsPos2 = []
  nextPos = []
  nextPos2 = []
  for i = 1:N
    push!(currObsPos, Array{Float64}(undef, 4, 2))
    push!(nextObsPos, Array{Float64}(undef, 4, 2))
    push!(nextObsPos2, Array{Float64}(undef, 4, 2))
    push!(nextPos, R[i].robotPose)
    push!(nextPos2, R[i].robotPose)
  end

  BVPEnds = []
  maxKDs = []
  tempKDs = [0.0]
  whichBVP = []
  level = []
  for i = 1:N
    tempEnds = [R[i].robotPose]
    push!(BVPEnds, tempEnds)
    push!(maxKDs, tempKDs)
    push!(whichBVP, 1)
    push!(level, 0)
  end

  BVPJustChanged = []
  BVPCounter = []
  NextBVPCheck = []
  for i = 1:N
    push!(BVPCounter, 0)
    push!(BVPJustChanged, false)
    push!(NextBVPCheck, true)
  end

  while true
    for i = 1:N
      if (NextBVPCheck[i] == false)
        if (BVPCounter[i] > 5)
          BVPCounter[i] = 0
          NextBVPCheck[i] = true
        else
          BVPCounter[i] += 1
        end
      end
      if (BVPJustChanged[i] == true)
        push!(BVPEnds[i], currPos[i])
        BVPJustChanged[i] = false
        whichBVP[i] += 1
        push!(maxKDs[i], 0.0)
        kdAndLV = sim_TNNLS_B_CT_Local_Max(BVPEnds[i][size(BVPEnds[i])[1] - 1][:], BVPEnds[i][size(BVPEnds[i])[1]][:], lastVel[i])
        lastVel[i] = kdAndLV[2]
        if ((maxKDs[i][size(maxKDs[i])[1] - 1]) > (kdAndLV[1] + .05))
          level[i] = 1
        end
      end
    end



    if (allDone == true)
      break
    end
    for i = 1:N
      hyberBallRad[i] = min(delta, ballConstant*((log(1+KD[i].treeSize)/(KD[i].treeSize))^(1/S[i].d)))
      itOfCheck[i][checkPtr[i]] += 1
    end
    now_time = time_ns()

    # calculate the end time of the first slice
    slice_end_time = (1+sliceCounter)*slice_time

    # see if warmup time has ended
    warmUpTimeJustEnded = false
    for i = 1:N
      if S[i].inWarmupTime && S[i].warmupTime < S[i].elapsedTime
        warmUpTimeJustEnded = true
        S[i].inWarmupTime = false
      end
    end

	### deal with larger kino_dist
  for i = 1:N
    if S[i].kino_dist > maxKDs[i][whichBVP[i]]
      maxKDs[i][whichBVP[i]] = S[i].kino_dist
    end

	  if S[i].kino_dist > S[i].augDist

      # not count in time for Minkowski sum
      before_save_time = time_ns()

      println("--------------------------------------------------------------------- New KD: $(S[i].kino_dist)")

	    S[i].augDist = S[i].kino_dist
	  
  	  obstacleAugmentation(S[i], S[i].augDist)

	    # add latest augmented obs
  	  list_item = S[i].obstacles.front
	    while list_item != list_item.child
		  ob = list_item.data
		  if !ob.obstacleUnused
		    addNewObstacle(S[i], KD[i], Q[i], ob, root[i], vCounter[i], R[i])
	      end
		  list_item = list_item.child
	    end
	  # println("-------------------------------------------------------------------------------R.currentMoveInvalid = $(R.currentMoveInvalid)")
	    propogateDescendants(Q[i], R[i])
	  # println("-------------------------------------------------------------------------------R.currentMoveInvalid = $(R.currentMoveInvalid)")

	    if R[i].currentMoveInvalid && explicitPointCheck(S[i], R[i].robotPose)[1] # R in aug obs && R not in ori obs
		  R[i].currentMoveInvalid = false
		  println("-------------------------------------------------------------------------------R.currentMoveInvalid = $(R[i].currentMoveInvalid) ----------- robot in aug obs")
	    end

  	  if !markedOS(S[i].moveGoal) # I'm pretty sure this is always true, since OS is emopty here -- M.O.
		  verifyInQueue(Q[i], S[i].moveGoal)
	    end
	    reduceInconsistency(Q[i], S[i].moveGoal, robotRads, root[i], hyberBallRad[i])

	    save_elapsed_time += (time_ns()-before_save_time)/1000000000

	  # println("obstacle enlarged")
	  end
  end

	### end of deal with larger kino_dist


    ### add/remove newly "detected" obstacles ###
    ### beginning of remove obstacle
    if vCounter[1] > 30
    for i = 1:N
      currObsPos[i][1,:] = [(currPos[i][1] - 1.0), (currPos[i][2] - 1.0)]
      currObsPos[i][2,:] = [(currPos[i][1] + 1.0), (currPos[i][2] - 1.0)]
      currObsPos[i][3,:] = [(currPos[i][1] + 1.0), (currPos[i][2] + 1.0)]
      currObsPos[i][4,:] = [(currPos[i][1] - 1.0), (currPos[i][2] + 1.0)]
      
      nextPos[i] = [(currPos[i][1] + 2*(currPos[i][1]-prevPos[i][5][1])), (currPos[i][2] + 2*(currPos[i][2]-prevPos[i][5][2]))]
      nextPos2[i] = [(currPos[i][1] + 4*(currPos[i][1]-prevPos[i][5][1])), (currPos[i][2] + 4*(currPos[i][2]-prevPos[i][5][2]))]

      nextObsPos[i][1,:] = [(nextPos[i][1] - 1.4), (nextPos[i][2] - 1.4)]
      nextObsPos[i][2,:] = [(nextPos[i][1] + 1.4), (nextPos[i][2] - 1.4)]
      nextObsPos[i][3,:] = [(nextPos[i][1] + 1.4), (nextPos[i][2] + 1.4)]
      nextObsPos[i][4,:] = [(nextPos[i][1] - 1.4), (nextPos[i][2] + 1.4)]

      nextObsPos2[i][1,:] = [(nextPos2[i][1] - 1.4), (nextPos2[i][2] - 1.4)]
      nextObsPos2[i][2,:] = [(nextPos2[i][1] + 1.4), (nextPos2[i][2] - 1.4)]
      nextObsPos2[i][3,:] = [(nextPos2[i][1] + 1.4), (nextPos2[i][2] + 1.4)]
      nextObsPos2[i][4,:] = [(nextPos2[i][1] - 1.4), (nextPos2[i][2] + 1.4)]

      currObs = Obstacle(3, currObsPos[i])
      nextObs = Obstacle(3, nextObsPos[i])
      nextObs2 = Obstacle(3, nextObsPos2[i])
    
      currObs.startTime = S[i].elapsedTime
      currObs.lifeSpan = slice_time*3
      currObs.obstacleUnused = true

      nextObs.startTime = S[i].elapsedTime
      nextObs.lifeSpan = slice_time*3
      nextObs.obstacleUnused = true

      nextObs2.startTime = S[i].elapsedTime
      nextObs2.lifeSpan = slice_time*3
      nextObs2.obstacleUnused = true

      #if (i != 2)
      #  if (Wdist(R[2].robotPose, currPos[i]) < robotSensorRange)
      #    addObsToCSpace(S[2], currObs)
      #    if (Wdist(R[2].robotPose, nextPos[i]) > 2.0)
      #      addObsToCSpace(S[2], nextObs)
      #    end
      #    if (Wdist(R[2].robotPose, nextPos2[i]) > 2.0)
      #      addObsToCSpace(S[2], nextObs2)
      #    end
      #  end
      #end
      if (i != 1)
        if (Wdist(R[1].robotPose, currPos[i]) < robotSensorRange)
          #if (level[i] == 0)
            addObsToCSpace(S[1], currObs)
          #end
          if (Wdist(R[1].robotPose, nextPos[i]) > 2.0)
            #if (level[i] == 0)
              addObsToCSpace(S[1], nextObs)
            #end
          end
          if (Wdist(R[1].robotPose, nextPos2[i]) > 2.0)
            #if (level[i] == 0)
              addObsToCSpace(S[1], nextObs2)
            #end
          end
        end
      end
    end
    end
    # remove obstacles at the required time
    for i = 1:N
    S[i].elapsedTime = (time_ns() - S[i].startTimeNs)/1000000000 - save_elapsed_time
    list_item = S[i].obstacles.front
    removedObstacle = false
    while list_item != list_item.child
      ob = list_item.data

      if !ob.senseableObstacle && !ob.obstacleUnused && (ob.startTime + ob.lifeSpan <= S[i].elapsedTime)
        # time to remove obstacle
        removeObstacle(S[i], KD[i], Q[i], ob, root[i], hyberBallRad[i], S[i].elapsedTime, S[i].moveGoal)
        removedObstacle = true
      elseif ob.senseableObstacle && ob.obstacleUnusedAfterSense && Wdist(R[i].robotPose, ob.position) < robotSensorRange + ob.radius
        # place to remove obstacle

        # because the space that used to be in this obstacle was never sampled
        # there will be a hole in the graph where it used to be. The following
        # attempts to mitigate this problem by requiring that the next few samples
        # come from the space that used to be inside the obstacle
        randomSampleObs(S[i], KD[i], ob) # stores samples in the sample stack
        removeObstacle(S[i], KD[i], Q[i], ob, root[i], hyberBallRad[i], S[i].elapsedTime, S[i].moveGoal)
        ob.senseableObstacle = false
        ob.startTime = Inf
        removedObstacle = true
      elseif S[i].spaceHasTime && ob.nextDirectionChangeTime > R[i].robotPose[3] && ob.lastDirectionChangeTime != R[i].robotPose[3]
        # a moving obstacle with unknown path is changing direction, so remove
        # its old anticipated trajectory

        removeObstacle(S[i], KD[i], Q[i], ob, root[i], hyberBallRad[i], S[i].elapsedTime, S[i].moveGoal)
        ob.obstacleUnused = false # obstacle is still used
        removedObstacle = true
      end

      list_item = list_item.child
    end
  

	  # if S.elapsedTime >= 13.0 && !environmentChangeFinished
	  # 	ob = S.obstacles.front.data
	  # 	randomSampleObs(S, KD, ob)
	  # 	removeObstacle(S, KD, Q, ob, root, hyberBallRad, S.elapsedTime, S.moveGoal)
	  # 	removedObstacle = true
	  # 	environmentChangeFinished = true
	  # end

	  if removedObstacle
      #println("----------------------------------------------------------------------------- Removed obstacle")
      reduceInconsistency(Q[i], S[i].moveGoal, robotRads, root[i], hyberBallRad[i])
    end
  end
    ### end of remove obstacle

    ### beginning of add obstacle
    # add obstacles at the required time
  for i = 1:N
    list_item = S[i].obstacles.front
    addedObstacle = false
    while list_item != list_item.child
      ob = list_item.data

      if !ob.senseableObstacle && ob.obstacleUnused && (ob.startTime <= S[i].elapsedTime <= ob.startTime + ob.lifeSpan)
        # time to add
        addNewObstacle(S[i], KD[i], Q[i], ob, root[i], vCounter[i], R[i])
        addedObstacle = true
      elseif ob.senseableObstacle && !ob.obstacleUnusedAfterSense && Wdist(R[i].robotPose, ob.position) < robotSensorRange + ob.radius
        # place to add obstacle
        addNewObstacle(S[i], KD[i], Q[i], ob, root[i], vCounter[i], R[i])
        ob.senseableObstacle = false
        addedObstacle = true
      elseif S[i].spaceHasTime && ob.nextDirectionChangeTime > R[i].robotPose[3] && ob.lastDirectionChangeTime != R[i].robotPose[3]
        # time that a moving obstacle with unknown path changes direction
        ob.obstacleUnused = false
        changeObstacleDirection(S[i], ob, R[i].robotPose[3])
        addNewObstacle(S[i], KD[i], Q[i], ob, root[i], vCounter[i], R[i])
        ob.lastDirectionChangeTime = copy(R[i].robotPose[3])
        #println("$(ob.nextDirectionChangeTime)  $(S.moveGoal.position[3]) ")
        addedObstacle = true
      elseif warmUpTimeJustEnded && !ob.obstacleUnused
        # warm up time is over, so we need to treat all active obstacles
        # as if they have just been added
        addNewObstacle(S[i], KD[i], Q[i], ob, root[i], vCounter[i], R[i])
        addedObstacle = true
      end

      list_item = list_item.child
    end
    if addedObstacle
      # propogate inf cost to all nodes beyond the obstacle and in its
      # basin of attraction
	  # println("-------------------------------------------------------------------------------R.currentMoveInvalid = $(R.currentMoveInvalid)")
	  propogateDescendants(Q[i], R[i])
	  # println("-------------------------------------------------------------------------------R.currentMoveInvalid = $(R.currentMoveInvalid)")
      if !markedOS(S[i].moveGoal) # I'm pretty sure this is always true, since OS is emopty here -- M.O.
        verifyInQueue(Q[i], S[i].moveGoal)
      end
      #println("--------------------------------------------------------------------------------- Added obstacle")
      reduceInconsistency(Q[i], S[i].moveGoal, robotRads, root[i], hyberBallRad[i])
    end
  end
    ### end of add obstacle

# 	if (S.kino_dist > augDist) || removedObstacle || addedObstacle
#
# 	  obstacleAugmented = false
#
# 	  if S.kino_dist > augDist
# 		println("kd increased")
# 		augDist = S.kino_dist
# 	  end
#
# 	  # pop augmented obs if needed
# 	  listEmpty(S.augObs)
# 	  # now S.obstacles only contains original obs, S.augObs is empty
# ##################################################################################################
# 	  # push augmented obs into S.obstacles and S.augObs
# 	  obstacleAugmentation(S, augDist)
#
# 	  list_item = S.augObs.front
# 	  while list_item != list_item.child
# 		  ob = list_item.data
# 		  addNewObstacle(S, KD, Q, ob, root, vCounter[i], R)
# 		  list_item = list_item.child
# 	  end
# 	  propogateDescendants(Q, R)
# 	  if !markedOS(S.moveGoal) # I'm pretty sure this is always true, since OS is emopty here -- M.O.
# 		  verifyInQueue(Q, S.moveGoal)
# 	  end
# 	  reduceInconsistency(Q, S.moveGoal, robotRads, root, hyberBallRad)
# 	  obstacleAugmented = true
# 	  println("added augmented obstacles")
#
#     end
    ### done with add/remove newly "detected" obstacles ###

    allDone = true
    for i = 1:N
      if R[i].robotPose == root[i].position
        break
      else
        allDone = false
      end
    end

    for i = 1:N
    # if this robot has used all of its allotted planning time of this slice
    S[i].elapsedTime = (time_ns() - S[i].startTimeNs)/1000000000 - save_elapsed_time
    if S[i].elapsedTime >= slice_end_time

      # calculate the end time of the next slice
      slice_end_time = (1+sliceCounter)*slice_time

      robot_slice_start = now_time
      if i == 1
        sliceCounter += 1
      end

      truncElapsedTime = floor(S[i].elapsedTime * 1000)/1000
      if i == 1
        println("slice $(sliceCounter) --- $(truncElapsedTime) -------- $(S[i].moveGoal.rrtTreeCost) $(S[i].moveGoal.rrtLMC) ----")
      end

      for j = 5:-1:2
        prevPos[i][j] = prevPos[i][j-1]
      end
      currPos[i] = R[i].robotPose
      prevPos[i][1] = currPos[i]

      if (vCounter[i] > 35)
        pastVec = [(prevPos[i][4][1] - prevPos[i][5][1]), (prevPos[i][4][2] - prevPos[i][5][2])]
        currVec = [(prevPos[i][1][1] - prevPos[i][2][1]), (prevPos[i][1][2] - prevPos[i][2][2])]
        angle = acos(((pastVec[1]*currVec[1]) + (pastVec[2]*currVec[2]))/(sqrt(pastVec[1]^2 + pastVec[2]^2)*sqrt(currVec[1]^2 + currVec[2]^2)))
        angle = angle*(360/(2*pi))
        if (abs(angle) > 18)
          BVPJustChanged[i] = true
          NextBVPCheck[i] = false
        end
      end

      # if saving stats
      if length(statsArgs) >= 1 && statsArgs[1]
        # record data
        elapsedTime[i][checkPtr[i]] = S[i].elapsedTime
      end

      ## move robot if the robot is allowed to move, otherwise planning is finished
      # so break out of the control loop
      if elapsedTime[i][checkPtr[i]] > total_planning_time + slice_time
        if MoveRobotFlag
          moveRobot_Q(S[i], Q[i], KD[i], slice_time, root[i], hyberBallRad[i], R[i], localPoseAndKd[i], localNormEsq[i], localTrigCond[i],save_elapsed_time)#, NormEsqvec, TrigCondvec) # 2 steps, update S.kino_dist
        else
          println("done (not moving robot)")
          break
        end
      end

      if searchType == "RRT#" || searchType == "RRTx"
        reduceInconsistency(Q[i], S[i].moveGoal, robotRads, root[i], hyberBallRad[i])
        if (S[i].moveGoal.rrtLMC != oldrrtLMC[i])
          oldrrtLMC[i] = (S[i].moveGoal.rrtLMC)
        end
	  end

      ## visualize graph #############
      if saveVideoData
        before_save_time = time_ns()

		# for visualization, S.obstacles only contains original obs
		# if obstacleAugmented
		# 	saveObstacleLocations(S.augObs, "temp/augObs_$(vCounter[i]).txt")
		# 	for i = 1:S.augObs.length
		# 		listPop(S.obstacles)
		# 	end
		# end

        saveRRTTree(KD[i], "temp/edges_$(i)_$(vCounter[i]).txt")
        saveRRTNodes(KD[i], "temp/nodes_$(i)_$(vCounter[i]).txt")
        #saveRRTNodesCollision(KD, "temp/Cnodes_$(vCounter[i]).txt")
        saveRRTPath_Q(S[i], S[i].moveGoal, root[i], R[i], "temp/path_$(i)_$(vCounter[i]).txt")
        saveObstacleLocations(S[i].obstacles, "temp/obstacles_$(i)_$(vCounter[i]).txt")
		    saveOriginalObstacleLocations_Q(S[i].obstacles, "temp/originalObs_$(i)_$(vCounter[i]).txt")
        saveData(R[i].robotMovePath[1:R[i].numRobotMovePoints,:], "temp/robotMovePath_$(i)_$(vCounter[i]).txt")
		    saveKds_Q(S[i], "temp/kd_$(i)_$(vCounter[i]).txt")

        
        vCounter[i] += 1
        S[i].fileCtr = vCounter[i]
        vCounter[i] = vCounter[1]


        save_elapsed_time += (time_ns()-before_save_time)/1000000000
      end
      ## end of visualize graph ######

      # check if the robot has reached its movement goal

      # if saving stats
      if length(statsArgs) >= 1 && statsArgs[1]
        # update statistics about run, assuming that we are saving them

        if checkPtr[i] < length(costOfGoal[i])
          checkPtr[i] += 1
          itOfCheck[i][checkPtr[i]] = itOfCheck[i][(checkPtr[i]-1)] + 1

          nodesInGraph[i][checkPtr[i]] = KD[i].treeSize
          costOfGoal[i][checkPtr[i]] = min(goal[i].rrtTreeCost, goal[i].rrtLMC)
          #costOfGoal[checkPtr[i]] = extractPathLength(goal , root)
          #numReduces[checkPtr[i]] = Q.numReduces
        else
          #println("WARNING: out of space to save stats")
        end
      end
    end
  end

    #### END of obstacle and robot pose update
    #### START of normal graph search stuff

    for i = 1:N
    # pick a random node
    newNode = S[i].randNode(S[i])

    if newNode.kdInTree # happens when we explicitly sample the goal every so often
      # nodes will be updated automatically as information gets propogated to it
      continue
    end


    # find closest old node to the new node
    #(closestNode, closestDist) = kdFindNearest(KD[i], newNode.position)
    (closestNode, closestDist) = kdFindNearestTerrain(KD[i], newNode.position)


    # saturate
    #if closestDist > delta && newNode != S.goalNode
    #  newNode.position = closestNode.position  + (newNode.position - closestNode.position)*delta/closestDist
    #end

    if closestDist > delta && newNode != S[i].goalNode
      #saturate(newNode.position, closestNode.position, delta)
      saturateTerrain(newNode.position, closestNode.position, delta)
    end



    # check for collisions vs static obstacles
    (explicitlyUnSafe, retCert) = explicitNodeCheck(S[i], newNode)

    if explicitlyUnSafe
      continue
    end

    #!!! Need to look into this
    GC.enable(false)

    # extend
    midNode = [((newNode.position[1]+closestNode.position[1])/2.0), ((newNode.position[2]+closestNode.position[2])/2.0)]
    if ((Wdist(midNode, [0.0, 0.0]) < 5.0) || (Wdist(midNode, [-2.0, -2.0]) < 5.0) || (Wdist(midNode, [-4.0, -4.0]) < 5.0))
      extendTerrain(S[i], KD[i], Q[i], newNode, closestNode, delta, hyberBallRad[i], S[i].moveGoal)
    else
      extend(S[i], KD[i], Q[i], newNode, closestNode, delta, hyberBallRad[i], S[i].moveGoal)
    end



    # make graph consistant (RRT# and RRTx)
    if searchType == "RRT#" || searchType == "RRTx"
      reduceInconsistency(Q[i], S[i].moveGoal, robotRads, root[i], hyberBallRad[i])
      if(S[i].moveGoal.rrtLMC != oldrrtLMC[i])
        #printPathLengths(S.moveGoal)
        oldrrtLMC[i] = S[i].moveGoal.rrtLMC
      end
    end

    GC.enable(true)
    end
  end

  ## end of while(true)
  for i = 1:N
  elapsedTime[i][checkPtr[i]] = (time_ns()-startTime)/1000000000

  if (length(statsArgs) >= 1 && statsArgs[1])
    if (!goal[i].rrtParentUsed)
      print("goal has no parent")
    end

    stats = hcat(elapsedTime[i], itOfCheck[i], nodesInGraph[i], costOfGoal[i])

    #saveData(stats[1:checkPtr[i],:], dataFileName[i])

    #reduceData = [numReduces', nodesInGraph']'
    #saveData(reduceData[1:checkPtr,:], "temp/reduceStats.txt")
  end
  moveLength = 0
  moveLength = sum(sqrt, sum((R[i].robotMovePath[1:R[i].numRobotMovePoints-1, :] - R[i].robotMovePath[2:R[i].numRobotMovePoints, :]).^2, dims=2))

  println("distance traveled by robot: $(moveLength[1])")
  println("KD_max: $(S[i].augDist)")
  end
  println(BVPEnds[1])
  println(maxKDs[1])
  println(level)
  return (S[1].NormEsqvec, S[1].TrigCondvec)
  # saveData(tr, "temp/Trig.txt")
  # saveData(er, "temp/Esq.txt")
end
