# this contains experiments of RRT-QX
include("heap.jl")
include("list.jl")
include("jlist.jl")
include("kdTree_general.jl")
include("DRRT_distance_functions.jl")
include("DRRT_SimpleEdge.jl")
const Edge{T} = SimpleEdge{T}
include("DRRT_data_structures.jl")
include("DRRT_SimpleEdge_functions.jl")
include("DRRT_Q.jl")
include("rrtqx.jl")
# include("movingrrtqx.jl")
include("Q-learning.jl")
# include("relaxedPE.jl")
include("obstacleAugmentation.jl")
# include("ACC2020.jl")
# include("kd.jl")



# this runs the specified algorithm in a 2D environemnt with unpredictably appearing/disapeearing obstacles
# saves data in ./temp to allow debugging
# NOTE: plot resulting simulation with "make_video_Q.m"
function sim_RRTQX()

  algorithmName = "RRTx"  # this is the algorithm to run, valid options are:
                         # "RRT", "RRT*", "RRT#", and "RRTx"

  changeThresh = 1.0     # ONLY USED FOR RRTx

  expName = "Debug"      # the name of this experiment (for naming output files)

  total_time = 3.0       # total planning time (move after this, and keep planning)
  slice_time = 1.0/10.0       # for saving data

  envRad = 18.0 # 50.0          # environment spans -envRad to envRad in each dimension
  robotRad = 0.5
  # robotSensorRange = 0.4
  start3 = [-3.0 14.0] # 5*[0.0 -8.0]   # robot goes to here (start location of search tree) [15.0 -10.0]
  goal3 = [-7.0 -7.0] # [-40.0 40.0]    # robot comes from here (goal location of search tree) [-18.0 18.0]
  start = [14.0 -5.0] # 5*[0.0 -8.0]   # robot goes to here (start location of search tree) [15.0 -10.0]
  goal = [-12.0 12.0] # [-40.0 40.0]    # robot comes from here (goal location of search tree) [-18.0 18.0]
  start4 = [14.0 -12.0] # 5*[0.0 -8.0]   # robot goes to here (start location of search tree) [15.0 -10.0]
  goal4 = [3.0 -6.0] # [-40.0 40.0]    # robot comes from here (goal location of search tree) [-18.0 18.0]
  start2 = [-12.0 17.0] # [-40.0 40.0]    # robot comes from here (goal location of search tree) [-18.0 18.0]
  goal2 = [12.0 -12.0] # 5*[0.0 -8.0]   # robot goes to here (start location of search tree) [15.0 -10.0]


  # start = [0.0 -8.]
  # goal = [-7.0 7.]
  obstacleFile = "environments/CDC2020Test.txt" # Dynamic_5: ACC 20; Static_5: small envir;
  # obstacleFile = "environments/empty.txt"
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
  C2 = CSpace{Float64}(d, -1.0, lowerBounds, upperBounds, start2, goal2)
  C3 = CSpace{Float64}(d, -1.0, lowerBounds, upperBounds, start3, goal3)
  C4 = CSpace{Float64}(d, -1.0, lowerBounds, upperBounds, start4, goal4)

  # init robot radii
  C.robotRadius = robotRad
  C2.robotRadius = robotRad
  C3.robotRadius = robotRad
  C4.robotRadius = robotRad

  # init robot velocity - useless in RRT-QX
  C.robotVelocity = 2.0
  C2.robotVelocity = 2.0
  C3.robotVelocity = 2.0
  C4.robotVelocity = 2.0

  # load obstacles
  readDiscoverablecObstaclesFromfile(C, obstacleFile, 1)
  readDiscoverablecObstaclesFromfile(C2, obstacleFile, 1)
  readDiscoverablecObstaclesFromfile(C3, obstacleFile, 1)
  readDiscoverablecObstaclesFromfile(C4, obstacleFile, 1)

  # set up sampling function
  C.randNode = randNodeOrFromStack # use this function to return random nodes
  C.pGoal = .01
  C2.randNode = randNodeOrFromStack # use this function to return random nodes
  C2.pGoal = .01
  C3.randNode = randNodeOrFromStack # use this function to return random nodes
  C3.pGoal = .01
  C4.randNode = randNodeOrFromStack # use this function to return random nodes
  C4.pGoal = .01

  # space paramiters
  C.spaceHasTime = false
  C.spaceHasTheta = false
  C2.spaceHasTime = false
  C2.spaceHasTheta = false
  C3.spaceHasTime = false
  C3.spaceHasTheta = false
  C4.spaceHasTime = false
  C4.spaceHasTheta = false

  dataFile = "experiments/$(expName)/debug_data.txt"

  T = []

  push!(T, C)
  push!(T, C2)
  push!(T, C3)
  push!(T, C4)
  N = (size(T)[1])

  multirrtqx(T, N, total_time, slice_time, 10.0, 100.0, changeThresh, algorithmName, MoveRobot, saveVideoData, saveTree, dataFile)

end