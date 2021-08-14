function makeDisplayNice()
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



# this contains the code that was used to run experiments and simulations for the paper

include("heap.jl")
include("list.jl")
include("jlist.jl")
include("kdTree_general.jl")
include("DRRT_distance_functions.jl")
include("DRRT_SimpleEdge.jl")
const Edge{T} = SimpleEdge{T}
include("DRRT_data_structures.jl")
include("DRRT_SimpleEdge_functions.jl")
include("DRRT.jl")

#--------------------------------------------------------------------------------#
#------------- after this line is code used for most recent version -------------#
#--------------------------------------------------------------------------------#

# this runs the specified algorithm in a static 2D environemnt and
# saves data in ./temp to allow debugging
# NOTE: plot resulting simulation with make "make_video.m"
function Static_2D_Debug()

  algorithmName = "RRTx"  # this is the algorithm to run, valid options are:
                         # "RRT", "RRT*", "RRT#", and "RRTx"

  changeThresh = 1.0     # ONLY USED FOR RRTx

  expName = "Debug"      # the name of this experiment (for naming output files)

  total_time = 8.0       # total planning time (move after this, and keep planning)
  slice_time = 1.0/10.0       # for saving data

  envRad = 50.0          # environment spans -envRad to envRad in each dimension
  robotRad = 0.5
  start = 5*[0.0 -8.0]   # robot goes to here (start location of search tree)
  goal = [-40.0 40.0]    # robot comes from here (goal location of search tree)

  obstacleFile = "environments/rand_Static.txt"
 # obstacleFile = "environments/empty.txt"
  success(`mkdir experiments/$(expName)`)
  # success(`cmd /c mkdir experiments/$(expName)`)

  MoveRobot = true
  saveVideoData = true

  d = 2                  # number of dimensions
  timeOut = Inf          # a backup timeout in seconds
  saveTree = true        # if true then saves the tree in out.txt

  lowerBounds = -envRad*ones(1,d)
  upperBounds = envRad*ones(1,d)

  C = CSpace{Float64}(d, -1.0, lowerBounds, upperBounds, start, goal)

  # init robot radii
  C.robotRadius =  robotRad

  # init robot velocity
  C.robotVelocity = 2.0

  # load obstacles
  readDiscoverablecObstaclesFromfile(C, obstacleFile, 1)

  # set up sampling function
  C.randNode = randNodeOrFromStack # use this function to return random nodes
  C.pGoal = .01

  # space paramiters
  C.spaceHasTime = false
  C.spaceHasTheta = false

  dataFile = "experiments/$(expName)/debug_data.txt"

  RRTX(C, total_time, slice_time, 5.0, 100.0, changeThresh, algorithmName, MoveRobot, saveVideoData, saveTree, dataFile);

end



# this runs the specified algorithm in a Dynamic 2D environemnt and
# saves data in ./temp to allow debugging
# NOTE: plot resulting simulation with make "make_video.m"
function Dynamic_2D_Debug()

  algorithmName = "RRTx" # this is the algorithm to run, valid options are:
                         # "RRT", "RRT*", "RRT#", and "RRTx"

  changeThresh = 1.0     # ONLY USED FOR RRTx

  expName = "Debug"      # the name of this experiment (for naming output files)

  total_time = 8.0       # total planning time (move after this, and keep planning)
  slice_time = 1.0/10.0  # how often we move robot (also for saving data =1/FPS)
  obstacleComplexity = 1 # insert each obstacle this many times in a row
                         # to simulate more complex environments

  # random obstacle test:
  envRad = 50.0          # environment spans -envRad to envRad in each dimension
  robotRad = 0.53
  start = [0.0 -40.0]   # robot goes to here (start location of search tree)
  goal = [-40.0 40.0]    # robot comes from here (goal location of search tree)

  obstacleFile = "environments/rand_DiscForest_.txt"
  # success(`mkdir experiments/$(expName)`)
  success(`cmd /c mkdir experiments/$(expName)`)

  MoveRobot = true
  saveVideoData = true

  d = 2                  # number of dimensions
  timeOut = Inf          # a backup timeout in seconds
  saveTree = true        # if true then saves the tree in out.txt

  lowerBounds = -envRad*ones(1,d)
  upperBounds = envRad*ones(1,d)

  C = CSpace{Float64}(d, -1.0, lowerBounds, upperBounds, start, goal)

  # init robot radii
  C.robotRadius =  robotRad

  # init robot velocity
  C.robotVelocity = 20.0
  # space paramiters
  C.spaceHasTime = false
  C.spaceHasTheta = false
  # load obstacles
  readDiscoverablecObstaclesFromfile(C, obstacleFile, 1)

  # set up sampling function
  C.randNode = randNodeOrFromStack # use this function to return random nodes
  C.pGoal = .01

  # space paramiters
  C.spaceHasTime = false
  C.spaceHasTheta = false

  dataFile = "experiments/$(expName)/debug_data.txt"

  RRTX(C, total_time, slice_time, 5.0, 100.0, changeThresh, algorithmName, MoveRobot, saveVideoData, saveTree, dataFile);

end


# this runs the specified algorithm in a static 2D + time environemnt and
# saves data in ./temp to allow debugging
# NOTE: plot resulting simulation with make "make_time_video_static.m"
function Static_2D_time_Debug()

  algorithmName = "RRTx" # "RRT", "RRT*", "RRT#", and "RRTx"
  changeThresh = 1.0     # ONLY USED FOR RRTx

  expName = "Debug"      # the name of this experiment (for naming output files)
  total_time = 3.0       # total planning time (move after this, and keep planning)
  slice_time = 1.0/10.0  # how often we move robot (also for saving data =1/FPS)

  envRad = 50.0          # environment spans -envRad to envRad in each dimension
  maxTime = 35.0         # this is the span along the time dimension that we care about
  robotRad = 5.0

  # note: time represents "time to arrive at the place the robot want to get to"
  start =[0.0 -40.0 0.0]        # robot goes to here (start location of search tree)
  goal = [0.0 40.0 maxTime]   # robot comes from here (goal location of search tree)

  obstacleFile = "environments/rand_StaticTime.txt"
  # success(`mkdir experiments/$(expName)`)
  success(`cmd /c mkdir experiments/$(expName)`)

  MoveRobot = true
  saveVideoData = true

  d = 3                  # number of dimensions [x, y, time]
  timeOut = Inf          # a backup timeout in seconds
  saveTree = true        # if true then saves the tree in out.txt

  lowerBounds = [-envRad*ones(1,2) 0.0]
  upperBounds = [envRad*ones(1,2) maxTime]

  C = CSpace{Float64}(d, -1.0, lowerBounds, upperBounds, start, goal)

  # init robot radii
  C.robotRadius =  robotRad

  # init robot velocity
  C.robotVelocity = 20.0

  # load obstacles
  readTimeObstaclesFromfile(C, obstacleFile, 1)

  # set up sampling function
  C.randNode = randNodeInTimeOrFromStack # use this function to return random nodes
  C.pGoal = .01

  # space paramiters
  C.spaceHasTime = true
  C.spaceHasTheta = false

  dataFile = "experiments/$(expName)/debug_data.txt"

  RRTX(C, total_time, slice_time, 10.0, 100.0, changeThresh, algorithmName, MoveRobot, saveVideoData, saveTree, dataFile);

end



# this runs the specified algorithm in a dynamic 2D + time environemnt and
# saves data in ./temp to allow debugging
# NOTE: plot resulting simulation with make "make_time_video_dynamic.m"
function Dynamic_2D_time_Debug()

 algorithmName = "RRTx" # this is the algorithm to run, valid options are:
                        # "RRT", "RRT*", "RRT#", and "RRTx"
  changeThresh = 1.0    # ONLY USED FOR RRTx

  expName = "Debug"     # the name of this experiment (for naming output files)

  total_time = 5.0      # total planning time (move after this, and keep planning)
  robotVelocity = 20.0  # m/s
  slice_time = 1.0/10.0 # how often we move robot (also for saving data =1/FPS)

  # random obstacle test:
  envRad = 50.0         # environment spans -envRad to envRad in each dimension
  maxTime = 35.0        # this is the span along the time dimension that we care about
  minTime = 0.0 #-15.0       # this is how far beyond time goal we allow planning

  robotRad = 2.0

  # note: time represents "time to arrive at the place the robot want to get to"
  start =[0.0 -40.0 minTime]  # robot goes to here (start location of search tree)
  goal = [0.0 40.0 maxTime]   # robot comes from here (goal location of search tree)

  obstacleFile = "environments/rand_DynamicTime.txt"
  # success(`mkdir experiments/$(expName)`)
  success(`cmd /c mkdir experiments/$(expName)`)

  MoveRobot = true
  saveVideoData = true

  d = 3                  # number of dimensions [x, y, time]
  timeOut = Inf          # a backup timeout in seconds
  saveTree = true        # if true then saves the tree in out.txt

  lowerBounds = [-envRad*ones(1,2) minTime]
  upperBounds = [envRad*ones(1,2) maxTime]

  C = CSpace{Float64}(d, -1.0, lowerBounds, upperBounds, start, goal)

  # init robot radii
  C.robotRadius =  robotRad

  # init robot velocity
  C.robotVelocity = 80.0

  # load obstacles
  readDynamicTimeObstaclesFromfile(C, obstacleFile, 1)

  # set up sampling function
  C.randNode = randNodeInTimeOrFromStack # use this function to return random nodes
  C.pGoal = .01

  # space paramiters
  C.spaceHasTime = true
  C.spaceHasTheta = false

  dataFile = "experiments/$(expName)/debug_data.txt"

  RRTX(C, total_time, slice_time, 10.0, 100.0, changeThresh, "RRTx" , MoveRobot, saveVideoData, saveTree, dataFile);

end



# ---------- end of "simulations", start of repeated trial "experiments" -------------

# this runs the specified algorithm in a Dynamic 2D environemnt and
# saves data in ./temp to allow debugging
function BottleNeckRepeats()
  firstTrial = 1
  totalTrials = 50  # run trials from firstTrial to totalTrials

  changeThresh = 0.5     # ONLY USED FOR RRTx
  expName = "repeatedBNeck"  # the name of this experiment (for naming output files)

  total_time = 60.0       # total planning time (move after this, and keep planning)
  slice_time = 1.0/10.0  # how often we move robot (also for saving data =1/FPS)
  obstacleComplexity = 1 # insert each obstacle this many times in a row
                         # to simulate more complex environments

  # bottleneck test:
  envRad = 50.0          # environment spans -envRad to envRad in each dimension
  robotRad = 0.5
  start = [0.0 -40.0]   # robot goes to here (start location of search tree)
  goal = [0.0 40.0]    # robot comes from here (goal location of search tree)

  obstacleFile = "environments/BNeck.txt"
  success(`mkdir experiments/$(expName)`)

  MoveRobot = false
  saveVideoData = false
  d = 2                  # number of dimensions
  timeOut = Inf          # a backup timeout in seconds
  saveTree = true        # if true then saves the tree in out.txt

  lowerBounds = -envRad*ones(1,d)
  upperBounds = envRad*ones(1,d)

  for trial = firstTrial:totalTrials
    for algType = 1:3

      if algType == 1
        algorithmName = "RRTx"
        algstr = "plus"
      elseif algType == 2
        algorithmName = "RRT#"
        algstr = "sharp"
      elseif algType == 3
        algorithmName = "RRT*"
        algstr = "star"
      end

      C = CSpace{Float64}(d, -1.0, lowerBounds, upperBounds, start, goal)

      # init robot radii
      C.robotRadius =  robotRad

      # init robot velocity (not really applicable for this experiment)
      C.robotVelocity = 20.0
      # space paramiters
      C.spaceHasTime = false
      C.spaceHasTheta = false
      # load obstacles
      readObstaclesFromfile(C, obstacleFile, 1)


      # set up sampling function (after 10 seconds removes bottleneck and sample it)
      C.randNode = randNodeTimeWithObstacleRemove # use this function to return random nodes
      C.pGoal = .01
      C.waitTime = 10.0
      C.obstaclesToRemove = C.obstacles.front.data
      C.timeSamplePoint = [0.0 0.0]

      # space paramiters
      C.spaceHasTime = false
      C.spaceHasTheta = false

      dataFile = "experiments/$(expName)/$(algstr)_$(trial).txt"

      println("#### running trial $(trial) of $(algorithmName)")

      RRTX(C, total_time, slice_time, 5.0, 100.0, changeThresh, algorithmName, MoveRobot, saveVideoData, saveTree, dataFile);

    end
  end
end

# this runs the specified algorithm in a Dynamic 2D environemnt and
# saves data in ./temp to allow debugging
function Static_2D_Repeats()
  firstTrial = 1
  totalTrials = 50  # run trials from firstTrial to totalTrials

  changeThresh = 0.5     # ONLY USED FOR RRTx
  expName = "repeated_Static_2D"  # the name of this experiment (for naming output files)

  total_time = 30.0       # total planning time (move after this, and keep planning)
  slice_time = 1.0/10.0  # how often we move robot (also for saving data =1/FPS)
  obstacleComplexity = 1 # insert each obstacle this many times in a row
                         # to simulate more complex environments

  # bottleneck test:
  envRad = 50.0          # environment spans -envRad to envRad in each dimension
  robotRad = 0.5
  start = [35.0 -40.0]   # robot goes to here (start location of search tree)
  goal = [-40.0 40.0]    # robot comes from here (goal location of search tree)

  obstacleFile = "environments/rand_Static.txt"
  success(`mkdir experiments/$(expName)`)

  MoveRobot = false
  saveVideoData = false
  d = 2                  # number of dimensions
  timeOut = Inf          # a backup timeout in seconds
  saveTree = true        # if true then saves the tree in out.txt

  lowerBounds = -envRad*ones(1,d)
  upperBounds = envRad*ones(1,d)

  for trial = firstTrial:totalTrials
    for algType = 1:3

      if algType == 1
        algorithmName = "RRTx"
        algstr = "plus"
      elseif algType == 2
        algorithmName = "RRT#"
        algstr = "sharp"
      elseif algType == 3
        algorithmName = "RRT*"
        algstr = "star"
      end

      C = CSpace{Float64}(d, -1.0, lowerBounds, upperBounds, start, goal)

      # init robot radii
      C.robotRadius =  robotRad

      # init robot velocity (not really applicable for this experiment)
      C.robotVelocity = 20.0
      # space paramiters
      C.spaceHasTime = false
      C.spaceHasTheta = false
      # load obstacles
      readDiscoverablecObstaclesFromfile(C, obstacleFile, 1)

      # set up sampling function
      C.randNode = randNodeOrFromStack # use this function to return random nodes
      C.pGoal = .01

      # space paramiters
      C.spaceHasTime = false
      C.spaceHasTheta = false

      dataFile = "experiments/$(expName)/$(algstr)_$(trial).txt"

      println("#### running trial $(trial) of $(algorithmName)")

      RRTX(C, total_time, slice_time, 5.0, 100.0, changeThresh, algorithmName, MoveRobot, saveVideoData, saveTree, dataFile);

    end
  end
end




# this runs the specified algorithm in a Dynamic 2D environemnt and
# saves data in ./temp to allow debugging
function BottleNeckRepeatsIts()
  firstTrial = 1
  totalTrials = 50  # run trials from firstTrial to totalTrials

  changeThresh = 0.5     # ONLY USED FOR RRTx
  #expName = "repeatedBNeckIts"  # the name of this experiment (for naming output files)
  expName = "repeatedBNeckIts2"  # the name of this experiment (for naming output files)


  total_time = 60.0       # total planning time (move after this, and keep planning)
  slice_time = 1.0/10.0  # how often we move robot (also for saving data =1/FPS)
  obstacleComplexity = 1 # insert each obstacle this many times in a row
                         # to simulate more complex environments

  # bottleneck test:
  envRad = 50.0          # environment spans -envRad to envRad in each dimension
  robotRad = 0.5
  start = [0.0 -40.0]   # robot goes to here (start location of search tree)
  goal = [0.0 40.0]    # robot comes from here (goal location of search tree)

  obstacleFile = "environments/BNeck.txt"
  success(`mkdir experiments/$(expName)`)

  MoveRobot = false
  saveVideoData = false
  d = 2                  # number of dimensions
  timeOut = Inf          # a backup timeout in seconds
  saveTree = true        # if true then saves the tree in out.txt

  lowerBounds = -envRad*ones(1,d)
  upperBounds = envRad*ones(1,d)

  for trial = firstTrial:totalTrials
    for algType = 1:3

      if algType == 1
        algorithmName = "RRTx"
        algstr = "plus"
      elseif algType == 2
        algorithmName = "RRT#"
        algstr = "sharp"
      elseif algType == 3
        algorithmName = "RRT*"
        algstr = "star"
      end

      C = CSpace{Float64}(d, -1.0, lowerBounds, upperBounds, start, goal)

      # init robot radii
      C.robotRadius =  robotRad

      # init robot velocity (not really applicable for this experiment)
      C.robotVelocity = 20.0
      # space paramiters
      C.spaceHasTime = false
      C.spaceHasTheta = false
      # load obstacles
      readObstaclesFromfile(C, obstacleFile, 1)


      # set up sampling function (after 10 seconds removes bottleneck and sample it)
      C.randNode = randNodeItsWithObstacleRemove # use this function to return random nodes
      C.pGoal = .01
      #C.itsUntilSample = 5000
      C.itsUntilSample = 6000
      C.obstaclesToRemove = C.obstacles.front.data
      C.timeSamplePoint = [0.0 0.0]

      # space paramiters
      C.spaceHasTime = false
      C.spaceHasTheta = false

      dataFile = "experiments/$(expName)/$(algstr)_$(trial).txt"

      println("#### running trial $(trial) of $(algorithmName)")

      RRTX(C, total_time, slice_time, 5.0, 100.0, changeThresh, algorithmName, MoveRobot, saveVideoData, saveTree, dataFile);

    end
  end
end




# this runs the specified algorithm in a Static 4D environemnt and
# saves data in ./temp to allow debugging
function Static_Time_Repeats()

  firstTrial = 6
  totalTrials = 50  # run trials from firstTrial to totalTrials

  changeThresh =  0.0 # 0.1 # 0.5   # ONLY USED FOR RRTx

  expName = "TimeRepeats4"  # the name of this experiment (for naming output files)


  total_time = 10.0*60       # total planning time (move after this, and keep planning)
  slice_time = 1.0  # how often we move robot (also for saving data =1/FPS)
  obstacleComplexity = 1 # insert each obstacle this many times in a row
                         # to simulate more complex environments

  envRad = 50.0          # environment spans -envRad to envRad in each dimension
  maxTime = 35.0         # this is the span along the time dimension that we care about
  minTime = -30.0
  robotRad = 2.0

  # note: time represents "time to arrive at the place the robot want to get to"
  start =[40.0 40.0 minTime pi/3]        # robot goes to here (start location of search tree)
  goal = [-40.0 -40.0 maxTime pi/3]   # robot comes from here (goal location of search tree)

  #obstacleFile = "environments/rand_StaticTime_7.txt"
  obstacleFile = "environments/rand_StaticTime_8.txt"

  success(`mkdir experiments/$(expName)`)

  MoveRobot = false
  saveVideoData = false

  d = 4                  # number of dimensions [x, y, 0.0, theta]
  timeOut = Inf          # a backup timeout in seconds
  saveTree = true        # if true then saves the tree in out.txt

  lowerBounds = [-envRad -envRad minTime 0.0]
  upperBounds = [envRad envRad maxTime 2*pi]

  for trial = firstTrial:totalTrials
    for algType =  1:3

      if algType == 1
        algorithmName = "RRTx"
        algstr = "plus"
      elseif algType == 2
        algorithmName = "RRT#"
        algstr = "sharp"
      elseif algType == 3
        algorithmName = "RRT*"
        algstr = "star"
      end

      C = CSpace{Float64}(d, -1.0, lowerBounds, upperBounds, start, goal)

      # init robot radii
      C.robotRadius =  robotRad

      # init robot velocity
      C.robotVelocity = 20.0
      C.dubbinsMinVelocity = 5.0
      C.dubbinsMaxVelocity = 20.0

      C.minTurningRadius = 2.0

      # load obstacles
      readTimeObstaclesFromfile(C, obstacleFile, 1)

      # set up sampling function
      C.randNode = randNodeOrFromStack # use this function to return random nodes
      C.pGoal = .01

      # space paramiters
      C.spaceHasTime = true
      C.spaceHasTheta = false

      dataFile = "experiments/$(expName)/$(algstr)_$(trial).txt"

      println("#### running trial $(trial) of $(algorithmName)")

      RRTX(C, total_time, slice_time, 10.0, 100.0, changeThresh, algorithmName, MoveRobot, saveVideoData, saveTree, dataFile);

    end
  end
end



# ----------------------------- more videos ---------------------------------


# this runs the specified algorithm in a Dynamic 2D environemnt and
# saves data in ./temp to allow debugging
function Dynamic_2D_Forest()

  algorithmName = "RRTx" # this is the algorithm to run, valid options are:
                         # "RRT", "RRT*", "RRT#", and "RRTx"

  changeThresh = 1.0     # ONLY USED FOR RRTx

  expName = "Debug"      # the name of this experiment (for naming output files)

  total_time = 10.0       # total planning time (move after this, and keep planning)
  slice_time = 1.0/10.0  # how often we move robot (also for saving data =1/FPS)
  obstacleComplexity = 1 # insert each obstacle this many times in a row
                         # to simulate more complex environments

  # random obstacle test:
  envRad = 50.0          # environment spans -envRad to envRad in each dimension
  robotRad = 1.5
  start = [40.0 40.0]   # robot goes to here (start location of search tree)
  goal = [-40.0 -40.0]    # robot comes from here (goal location of search tree)

  obstacleFile = "environments/rand_DiscForest_.txt"
  success(`mkdir experiments/$(expName)`)

  MoveRobot = true
  saveVideoData = true

  d = 2                  # number of dimensions
  timeOut = Inf          # a backup timeout in seconds
  saveTree = true        # if true then saves the tree in out.txt

  lowerBounds = -envRad*ones(1,d)
  upperBounds = envRad*ones(1,d)

  C = CSpace{Float64}(d, -1.0, lowerBounds, upperBounds, start, goal)

  # init robot radii
  C.robotRadius =  robotRad

  # init robot velocity
  C.robotVelocity = 20.0
  # space paramiters
  C.spaceHasTime = false
  C.spaceHasTheta = false
  # load obstacles
  readDiscoverablecObstaclesFromfile(C, obstacleFile, 1)

  # set up sampling function
  C.randNode = randNodeOrFromStack # use this function to return random nodes
  C.pGoal = .01

  # space paramiters
  C.spaceHasTime = false
  C.spaceHasTheta = false

  dataFile = "experiments/$(expName)/debug_data.txt"

  RRTX(C, total_time, slice_time, 5.0, 100.0, changeThresh, algorithmName, MoveRobot, saveVideoData, saveTree, dataFile);

end


# this runs the specified algorithm in a Dynamic 2D environemnt and
# saves data in ./temp to allow debugging
function Dynamic_2D_Fort()

  algorithmName = "RRTx" # this is the algorithm to run, valid options are:
                         # "RRT", "RRT*", "RRT#", and "RRTx"

  changeThresh = 1.0     # ONLY USED FOR RRTx

  expName = "Debug"      # the name of this experiment (for naming output files)

  total_time = 5.0       # total planning time (move after this, and keep planning)
  slice_time = 1.0/10.0  # how often we move robot (also for saving data =1/FPS)
  obstacleComplexity = 1 # insert each obstacle this many times in a row
                         # to simulate more complex environments

  # random obstacle test:
  envRad = 50.0          # environment spans -envRad to envRad in each dimension
  robotRad = 1.5
  start = [0.0 -45.0]   # robot goes to here (start location of search tree)
  goal = [0.0 0.0]    # robot comes from here (goal location of search tree)

  obstacleFile = "environments/rand_Disc_3.txt"
  success(`mkdir experiments/$(expName)`)

  MoveRobot = true
  saveVideoData = true

  d = 2                  # number of dimensions
  timeOut = Inf          # a backup timeout in seconds
  saveTree = true        # if true then saves the tree in out.txt

  lowerBounds = -envRad*ones(1,d)
  upperBounds = envRad*ones(1,d)

  C = CSpace{Float64}(d, -1.0, lowerBounds, upperBounds, start, goal)

  # init robot radii
  C.robotRadius =  robotRad

  # init robot velocity
  C.robotVelocity = 20.0
  # space paramiters
  C.spaceHasTime = false
  C.spaceHasTheta = false
  # load obstacles
  readDiscoverablecObstaclesFromfile(C, obstacleFile, 1)

  # set up sampling function
  C.randNode = randNodeOrFromStack # use this function to return random nodes
  C.pGoal = .01

  # space paramiters
  C.spaceHasTime = false
  C.spaceHasTheta = false

  dataFile = "experiments/$(expName)/debug_data.txt"

  RRTX(C, total_time, slice_time, 5.0, 100.0, changeThresh, algorithmName, MoveRobot, saveVideoData, saveTree, dataFile);

end


# this runs the specified algorithm in a static 2D + time environemnt and
# saves data in ./temp to allow debugging
function Static_2D_time_grid()

  algorithmName = "RRTx" # this is the algorithm to run, valid options are:
                         # "RRT", "RRT*", "RRT#", and "RRTx"
  changeThresh = 1.0     # ONLY USED FOR RRTx

  expName = "Debug"      # the name of this experiment (for naming output files)
  total_time = 3.0       # total planning time (move after this, and keep planning)
  slice_time = 1.0/10.0  # how often we move robot (also for saving data =1/FPS)

  envRad = 50.0          # environment spans -envRad to envRad in each dimension
  maxTime = 35.0         # this is the span along the time dimension that we care about
  minTime = 20
  robotRad = 2.0

  # note: time represents "time to arrive at the place the robot want to get to"
  start =[40.0 40.0 minTime]        # robot goes to here (start location of search tree)
  goal = [-40.0 -40.0 maxTime]   # robot comes from here (goal location of search tree)

  obstacleFile = "environments/rand_StaticTime_6.txt"
  success(`mkdir experiments/$(expName)`)

  MoveRobot = true
  saveVideoData = true

  d = 3                  # number of dimensions [x, y, time]
  timeOut = Inf          # a backup timeout in seconds
  saveTree = true        # if true then saves the tree in out.txt

  lowerBounds = [-envRad*ones(1,2) 0.0]
  upperBounds = [envRad*ones(1,2) maxTime]

  C = CSpace{Float64}(d, -1.0, lowerBounds, upperBounds, start, goal)

  # init robot radii
  C.robotRadius =  robotRad

  # init robot velocity
  C.robotVelocity = 20.0

  # load obstacles
  readTimeObstaclesFromfile(C, obstacleFile, 1)

  # set up sampling function
  C.randNode = randNodeInTimeOrFromStack # use this function to return random nodes
  C.pGoal = .01

  # space paramiters
  C.spaceHasTime = true
  C.spaceHasTheta = false

  dataFile = "experiments/$(expName)/debug_data.txt"

  RRTX(C, total_time, slice_time, 10.0, 100.0, changeThresh, algorithmName, MoveRobot, saveVideoData, saveTree, dataFile);

end


# this runs the specified algorithm in a dynamic 2D + time environemnt and
# saves data in ./temp to allow debugging
function Dynamic_2D_time_Busy()

 algorithmName = "RRTx" # this is the algorithm to run, valid options are:
                        # "RRT", "RRT*", "RRT#", and "RRTx"
  changeThresh = 1.0    # ONLY USED FOR RRTx

  expName = "Debug"     # the name of this experiment (for naming output files)

  total_time = 5.0      # total planning time (move after this, and keep planning)
  robotVelocity = 20.0  # m/s
  slice_time = 1.0/10.0 # how often we move robot (also for saving data =1/FPS)

  # random obstacle test:
  envRad = 50.0         # environment spans -envRad to envRad in each dimension
  maxTime = 35.0        # this is the span along the time dimension that we care about
  minTime = 0.0 #-15.0       # this is how far beyond time goal we allow planning

  robotRad = 2.0

  # note: time represents "time to arrive at the place the robot want to get to"
  start =[-40.0 -40.0 minTime]  # robot goes to here (start location of search tree)
  goal = [40.0 40.0 maxTime]   # robot comes from here (goal location of search tree)

  obstacleFile = "environments/rand_DynamicTime_3.txt"
  success(`mkdir experiments/$(expName)`)

  MoveRobot = true
  saveVideoData = true

  d = 3                  # number of dimensions [x, y, time]
  timeOut = Inf          # a backup timeout in seconds
  saveTree = true        # if true then saves the tree in out.txt

  lowerBounds = [-envRad*ones(1,2) minTime]
  upperBounds = [envRad*ones(1,2) maxTime]

  C = CSpace{Float64}(d, -1.0, lowerBounds, upperBounds, start, goal)

  # init robot radii
  C.robotRadius =  robotRad

  # init robot velocity
  C.robotVelocity = 80.0

  # load obstacles
  readDynamicTimeObstaclesFromfile(C, obstacleFile, 1)

  # set up sampling function
  C.randNode = randNodeInTimeOrFromStack # use this function to return random nodes
  C.pGoal = .01

  # space paramiters
  C.spaceHasTime = true
  C.spaceHasTheta = false

  dataFile = "experiments/$(expName)/debug_data.txt"

  RRTX(C, total_time, slice_time, 10.0, 100.0, changeThresh, "RRTx" , MoveRobot, saveVideoData, saveTree, dataFile);

end
