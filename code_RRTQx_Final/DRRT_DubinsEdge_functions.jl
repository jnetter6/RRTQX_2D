function makeMyTextEditorDisplayNice20() # thanks!
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
# use the following version of dist for Dubins space [X Y Theta Time]
dist(x::Array,y::Array) = R3SDist(x,y)


## This returns the distance that is used internally in the kd-tree (see notes
## there!!! the kd-implimentation handles wrapping dimensions (e.g. a circle or
## taurus) by running multipe searches, 1 for each wrapped identity of the
## querery point. So this is the distance function that is the -non-wrapping-
## version of what the actual distance function between points in the KD space is
#
# use the following version of KDdist for Dubins space (with and without time, since
# time === 0 in the latter case, but still exists as a dimension in points)
KDdist(x::Array,y::Array) = euclidianDist(x,y)


## returns the workspace distance between two points, this should obey the triangle
## inequality. e.g., in the current version of this code, it is assumed that the
## workspace is the first two dimensions of the C-space, it is used for calculating
## the "real world" distance that the robot covers (e.g., in a particualr ammount
## of time when determining if a move is possible given robot max speed).
#
# use the following version of Wdist for Euclidian space and Dubin's space
Wdist(x::Array,y::Array) = euclidianDist(x[1:2],y[1:2])



# moves newPoint toward closestPoint such that each robot is no
# further than delta, points repesent the cartesian product of R robots
#
# use the following version of saturate for Dubins space
function saturate(newPoint::Array{Float64}, closestPoint::Array{Float64}, delta)

  thisDist = dist(newPoint, closestPoint)
  if thisDist > delta

    # first scale non-theta dimensions
    newPoint[1:3] = closestPoint[1:3]  + (newPoint[1:3] - closestPoint[1:3])*delta/thisDist

    # saturate theta in the shorter of the two directions that it can go
    if abs(newPoint[4] - closestPoint[4]) < pi
      # saturate in the normal way
      newPoint[4] = closestPoint[4]  + (newPoint[4] - closestPoint[4])*delta/thisDist
    else
      # saturate in the opposite way

      # start at wrapped identity of new point theta
      newPoint[4] = (newPoint[4] < pi ? newPoint[4] + 2*pi : newPoint[4] - 2*pi)

      # now saturate
      newPoint[4] = closestPoint[4]  + (newPoint[4] - closestPoint[4])*delta/thisDist

      # finally, wrap back to the identity that is on [0 2*pi]
      newPoint[4] = max(min(newPoint[4], 2*pi), 0.0)
    end
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
## DubinsEdgeVersion
function validMove(S::TS, edge::DubinsEdge{T}) where {T, TS}
  if S.spaceHasTime
    # note that planning happens in reverse time, i.e., time = 0 is at the root
    # of the search tree, and thus the time of startNode must be greater than
    # the time of the end node.
    return ((edge.startNode.position[3] > edge.endNode.position[3]) && (S.dubinsMinVelocity <= edge.velocity <= S.dubinsMaxVelocity))
  end

  # if space does not have time then we assume that the move is always valid
  return true
end



### returns the pose of a robot that is located dist along the edge
### NOTE THAT 'dist' and 'far' are with respect to whatever type
### of distance is stored in edge.dist
#
## Dubins edge version, NOTE THAT THIS COULD BE MADE MORE EFFICIENT IN A
## NUMBER OF WAYS.
function poseAtDistAlongEdge(edge::DubinsEdge, distAlongEdge::Float64)
  distRemaining = copy(distAlongEdge)

  if size(edge.trajectory,1) < 2 || edge.dist <= distAlongEdge
    return copy(edge.endNode.position)
  end

  # find the piece of the trajectory that contains the point at the
  # desired distance
  i = 2
  thisDist = Inf
  timeInPath = (size(edge.trajectory,2) >= 3)
  while i <= size(edge.trajectory,1)
    thisDist = ( timeInPath ?
        dubinsDistAlongTimePath(edge.trajectory[i-1,:], edge.trajectory[i,:]) :
        dubinsDistAlongPath(edge.trajectory[i-1,:], edge.trajectory[i,:])       )

    if distRemaining - thisDist <= 0
      break
    end

    distRemaining -= thisDist
    i += 1
  end

  if distRemaining > thisDist
    # in case of rare subtraction based percision errors
    distRemaining = thisDist
  end

  # now calculate pose along that piece
  ratio = distRemaining/thisDist
  ret = edge.trajectory[i-1,:] + ratio*(edge.trajectory[i,:] - edge.trajectory[i-1,:])
  retTimeRatio = distAlongEdge/edge.dist
  retTime = edge.startNode.position[3] + retTimeRatio*(edge.endNode.position[3] - edge.startNode.position[3])
  retTheta = atan( edge.trajectory[i,2] - edge.trajectory[i-1,2] ,  edge.trajectory[i,1] - edge.trajectory[i-1,1])

  return [ret[1] ret[2] retTime retTheta]
end



### returns the pose of a robot that is located time along the edge
#
## Dubins edge version, NOTE THAT THIS COULD BE MADE MORE EFFICIENT IN A
## NUMBER OF WAYS
function poseAtTimeAlongEdge(edge::DubinsEdge, timeAlongEdge::Float64)
  if size(edge.trajectory,1) < 2 || (edge.startNode.position[3] - edge.endNode.position[3]) <= timeAlongEdge

    return copy(edge.endNode.position)
  end

  # find the piece of the trajectory that contains the time at the
  # desired distance
  i = 2
  while edge.trajectory[i,3] > edge.startNode.position[3] - timeAlongEdge
   i += 1
  end

  # now calculate pose along that piece
  ratio = (edge.trajectory[i-1,3] - (edge.startNode.position[3] - timeAlongEdge))/(edge.trajectory[i-1,3] - edge.trajectory[i,3])

  ret = edge.trajectory[i-1,:] + ratio*(edge.trajectory[i,:] - edge.trajectory[i-1,:])
  retTime = edge.startNode.position[3] - timeAlongEdge
  retTheta = atan( edge.trajectory[i,2] - edge.trajectory[i-1,2] ,  edge.trajectory[i,1] - edge.trajectory[i-1,1])

  return [ret[1] ret[2] retTime retTheta]
end


### the next function saves the last part of the trajectory into the file filePtr,
### starting distFromFront from the front. NOTE THAT 'distFromFront' is with
### respect to whatever type of distance is stored in edge.dist
#
## Dubins edge. NOTE edge.dist is the distance through the subspace
## defined by either [X Y Time] or [X Y] depending of if time is being
## used or not, respectively. Note that this is not, in general, the same
## distance that is retured by either dist() or Wdist() which is why
## I have hard-coded sqrt(sum((x-y).^2)) instead
function saveEndOfTrajectory(filePtr::TFP, edge::DubinsEdge,
  distFromFront::Float64) where {TFP}
  distRemaining = copy(distFromFront)

  if size(edge.trajectory,1) < 2 || edge.dist <= distFromFront
    writedlm(filePtr, reshape(edge.endNode.position[1,1:2], 1, 2), ',')
    return
  end

  # find the piece of the trajectory that contains the point at the
  # desired distance
  i = 2
  thisDist = Inf
  timeInPath = (size(edge.trajectory,2) >= 3)
  while i <= size(edge.trajectory,1)
    thisDist = ( timeInPath ?
        dubinsDistAlongTimePath(edge.trajectory[i-1,:], edge.trajectory[i,:]) :
        dubinsDistAlongPath(edge.trajectory[i-1,:], edge.trajectory[i,:])       )

    if distRemaining - thisDist <= 0
      break
    end

    distRemaining -= thisDist
    i += 1
  end

  if distRemaining > thisDist
    # in case of rare subtraction based percision errors
    distRemaining = thisDist
  end

  # now calculate pose along that piece
  ratio = distRemaining/thisDist
  ret = edge.trajectory[i-1,:] + ratio*(edge.trajectory[i,:] - edge.trajectory[i-1,:])

  # retTimeRatio = distFromFront/edge.dist
  # retTime = edge.startNode.position[3] + retTimeRatio*(edge.endNode.position[3] - edge.startNode.position[3])
  # retTheta = atan( edge.trajectory[i,2] - edge.trajectory[i-1,2] ,  edge.trajectory[i,1] - edge.trajectory[i-1,1])
  # writecsv(filePtr, [ret[1] ret[2] retTime retTheta])
  writedlm(filePtr, ret, ',')

  # finally, save the rest of the points of the trajectory into the file
  #totalDistAt_j = distFromFront - distRemaining + thisDist
  for j = i:size(edge.trajectory,1)

    positionAt_j = edge.trajectory[j,:]
    #timeRatio_j = positionAt_j/edge.dist
    #timeAt_j = edge.startNode.position[3] + timeRatio_j*(edge.endNode.position[3] - edge.startNode.position[3])
    #thetaAt_j = atan( edge.trajectory[j,2] - edge.trajectory[j-1,2] ,  edge.trajectory[j,1] - edge.trajectory[j-1,1])
    #writecsv(filePtr, [positionAt_j[1] positionAt_j[2] timeAt_j thetaAt_j])

    writedlm(filePtr, positionAt_j, ',')

    # account for dist that the next point in the trajectory will add
    #if j < size(edge.trajectory,1)
    #  totalDistAt_j += dist(edge.trajectory[j,:], edge.trajectory[j+1,:])
    #end
  end
end



### the next function saves the last part of the trajectory into the file filePtr,
# starting timeFromFront from the front
#
## Dubins edge version, NOTE THAT THIS COULD BE MADE MORE EFFICIENT IN A
## NUMBER OF WAYS
function saveEndOfTrajectoryTime(filePtr::TFP, edge::DubinsEdge,
  timeFromFront::Float64) where {TFP}
  if size(edge.trajectory,1) < 2 || (edge.startNode.position[3] - edge.endNode.position[3]) <= timeFromFront

    writedlm(filePtr, reshape(edge.endNode.position[1,1:3], 1, 3), ',')
    return
  end

  # find the piece of the trajectory that contains the time at the
  # desired distance
  i = 2
  while edge.trajectory[i,3] > edge.startNode.position[3] - timeFromFront
   i += 1
  end

  # now calculate pose along that piece
  ratio = (edge.trajectory[i-1,3] - (edge.startNode.position[3] - timeFromFront))/(edge.trajectory[i-1,3] - edge.trajectory[i,3])

  ret = edge.trajectory[i-1,:] + ratio*(edge.trajectory[i,:] - edge.trajectory[i-1,:])
  retTime = edge.startNode.position[3] - timeFromFront
  # retTheta = atan( edge.trajectory[i,2] - edge.trajectory[i-1,2] ,  edge.trajectory[i,1] - edge.trajectory[i-1,1])

  writedlm(filePtr, [ret[1] ret[2] retTime], ',')

  # finally, save the rest of the points of the trajectory into the file
  #totalDistAt_j = distFromFront - distRemaining + thisDist
  for j = i:size(edge.trajectory,1)
    writedlm(filePtr, reshape(edge.trajectory[j,1:3], 1, 3), ',')
  end
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
## dubins version, figures out which one of the 6 possibilities is the shortest
## (ignoring obstacles) subject to the robot's (constant) velocity and minimum
## turning radius
function calculateTrajectory(S::TS, edge::DubinsEdge) where {TS}

  r_min = S.minTurningRadius

  initial_location = edge.startNode.position[1:2];
  initial_theta = edge.startNode.position[4];

  goal_location = edge.endNode.position[1:2];
  goal_theta = edge.endNode.position[4];

  #println("initial_location = [$(initial_location[1]), $(initial_location[2])]")
  #println("initial_theta = $(initial_theta)")
  #println("goal_location = [$(goal_location[1]), $(goal_location[2])]")
  #println("goal_theta = $(goal_theta)")


  # calculate the center of the right-turn and left turn circles

  # ... right-turn initial_location circle
  irc_center = initial_location + r_min * [cos(initial_theta-pi/2.0), sin(initial_theta-pi/2.0)];

  # ... left-turn initial_location circle
  ilc_center = initial_location + r_min * [cos(initial_theta+pi/2.0), sin(initial_theta+pi/2.0)];

  # ... right-turn goal_location circle
  grc_center = goal_location + r_min * [cos(goal_theta-pi/2.0), sin(goal_theta-pi/2.0)];

  # ... the left-turn initial_location circle
  glc_center = goal_location + r_min * [cos(goal_theta+pi/2.0), sin(goal_theta+pi/2.0)];


  # remember the best distance and the associated type of trajectory
  bestDist = Inf
  bestTrajType = "xxx"


  # calculate tangent points for right-straight-left path
  # r-s-l requires "inner" tangent points
  D = sqrt(sum((glc_center - irc_center).^2.0));
  v = (glc_center - irc_center)./D;
  R = -2.0*r_min/D;
  if abs(R) > 1.0
    rsl_length = Inf
  else
    sq = sqrt(1.0-R^2);
    a = r_min*(R*v[1] + v[2]*sq);
    b = r_min*(R*v[2] - v[1]*sq);
    rsl_tangent_x = [irc_center[1] - a, glc_center[1] + a];
    rsl_tangent_y = [irc_center[2] - b, glc_center[2] + b];

    firstDist = rightTurnDist(initial_location, [rsl_tangent_x[1], rsl_tangent_y[1]], irc_center, r_min);
    secondDist = sqrt(sum(([rsl_tangent_x[2], rsl_tangent_y[2]] - [rsl_tangent_x[1], rsl_tangent_y[1]]).^2));
    thirdDist = leftTurnDist([rsl_tangent_x[2], rsl_tangent_y[2]], goal_location, glc_center, r_min);
    rsl_length = firstDist+secondDist+thirdDist

    if bestDist > rsl_length
      bestDist = rsl_length
      bestTrajType = "rsl"
    end
  end


  # calculate tangent points for right-straight-right path
  # r-s-r requires "outer" tangent points, and the tangent vector is parallel
  # to the vector from irc to grc
  D = sqrt(sum((grc_center - irc_center).^2));
  v = (grc_center - irc_center)./D;
  rsr_tangent_x = -r_min*v[2] .+ [irc_center[1], grc_center[1]];
  rsr_tangent_y = r_min*v[1] .+ [irc_center[2], grc_center[2]];

  firstDist = rightTurnDist(initial_location, [rsr_tangent_x[1], rsr_tangent_y[1]], irc_center, r_min);
  secondDist = sqrt(sum(([rsr_tangent_x[2], rsr_tangent_y[2]] - [rsr_tangent_x[1], rsr_tangent_y[1]]).^2));
  thirdDist = rightTurnDist([rsr_tangent_x[2], rsr_tangent_y[2]], goal_location, grc_center, r_min);
  rsr_length = firstDist+secondDist+thirdDist

  if bestDist > rsr_length
    bestDist = rsr_length
    bestTrajType = "rsr"
  end


  # calculate (if advantagous) the right-left-right path
  rlr_rl_tangent = [NaN NaN]
  rlr_lr_tangent = [NaN NaN]
  if D < 4.0*r_min
    # start by finding the center of the "left" circle
    theta = -acos(D/(4*r_min)) + atan(v[2], v[1]);
    rlr_l_circle_center = irc_center + 2*r_min*[cos(theta), sin(theta)];
    rlr_rl_tangent = (rlr_l_circle_center + irc_center)./2.0;
    rlr_lr_tangent = (rlr_l_circle_center + grc_center)./2.0;

    firstDist = rightTurnDist(initial_location, rlr_rl_tangent, irc_center, r_min);
    secondDist = leftTurnDist(rlr_rl_tangent, rlr_lr_tangent, rlr_l_circle_center, r_min);
    thirdDist = rightTurnDist(rlr_lr_tangent, goal_location, grc_center, r_min);
    rlr_length = firstDist+secondDist+thirdDist

    if bestDist > rlr_length
      bestDist = rlr_length
      bestTrajType = "rlr"
    end
  else
    rlr_length = Inf
  end


  # calculate tangent points for left-straight-right path
  # l-s-r requires "inner" tangent points
  D = sqrt(sum((grc_center - ilc_center).^2));
  v = (grc_center - ilc_center)./D;
  R = 2.0*r_min/D;
  if abs(R) > 1
    lsr_length = Inf
  else
    sq = sqrt(1-R^2);
    a = R*v[1] + v[2]*sq;
    b = R*v[2] - v[1]*sq;
    lsr_tangent_x = [ilc_center[1] + a*r_min, grc_center[1] - a*r_min];
    lsr_tangent_y = [ilc_center[2] + b*r_min, grc_center[2] - b*r_min];

    firstDist = leftTurnDist(initial_location, [lsr_tangent_x[1] lsr_tangent_y[1]], ilc_center, r_min);
    secondDist = sqrt(sum(([lsr_tangent_x[2], lsr_tangent_y[2]] - [lsr_tangent_x[1], lsr_tangent_y[1]]).^2));
    thirdDist = rightTurnDist([lsr_tangent_x[2], lsr_tangent_y[2]], goal_location, grc_center, r_min);
    lsr_length = firstDist+secondDist+thirdDist

    if bestDist > lsr_length
      bestDist = lsr_length
      bestTrajType = "lsr"
    end

  end


  # calculate tangent points for left-straight-left path
  # l-s-l requires "outer" tangent points, and the tangent vector is parallel
  # to the vector from irc to grc
  D = sqrt(sum((glc_center - ilc_center).^2));
  v = (glc_center - ilc_center)./D;
  lsl_tangent_x = r_min*v[2] .+ [ilc_center[1], glc_center[1]];
  lsl_tangent_y = -r_min*v[1] .+ [ilc_center[2], glc_center[2]];

  firstDist = leftTurnDist(initial_location, [lsl_tangent_x[1] lsl_tangent_y[1]], ilc_center, r_min);
  secondDist = sqrt(sum(([lsl_tangent_x[2], lsl_tangent_y[2]] - [lsl_tangent_x[1], lsl_tangent_y[1]]).^2));
  thirdDist = leftTurnDist([lsl_tangent_x[2], lsl_tangent_y[2]], goal_location, glc_center, r_min);
  lsl_length = firstDist+secondDist+thirdDist

  if bestDist > lsl_length
    bestDist = lsl_length
    bestTrajType = "lsl"
  end


  # calculate (if advantagous) the left-right-left path
  lrl_rl_tangent = [NaN NaN]
  lrl_lr_tangent = [NaN NaN]
  if D < 4.0*r_min
    # start by finding the center of the "right" circle
    theta = acos(D/(4*r_min)) + atan(v[2], v[1]);
    lrl_r_circle_center = ilc_center + 2.0*r_min*[cos(theta), sin(theta)];
    lrl_lr_tangent = (lrl_r_circle_center + ilc_center)./2.0;
    lrl_rl_tangent = (lrl_r_circle_center + glc_center)./2.0;

    firstDist = leftTurnDist(initial_location, lrl_lr_tangent, ilc_center, r_min);
    secondDist = rightTurnDist(lrl_lr_tangent, lrl_rl_tangent, lrl_r_circle_center, r_min);
    thirdDist = leftTurnDist(lrl_rl_tangent, goal_location, glc_center, r_min);
    lrl_length = firstDist+secondDist+thirdDist

    if bestDist > lrl_length
      bestDist = lrl_length
      bestTrajType = "lrl"
    end
  else
    lrl_length = Inf
  end


  # now save the best path in the trajectory field

  delta_phi = .1; # this is the angle granularity (in radians) used for
                  # discritizing the arcs of the paths (straight lines
                  # are saved as a single segment)

  # calculate first part of path
  if bestTrajType[1] == 'r'
    # first part of path is a right hand turn

    if bestTrajType == "rsl"
      p = [rsl_tangent_x[1], rsl_tangent_y[1]];
    elseif bestTrajType == "rsr"
      p = [rsr_tangent_x[1], rsr_tangent_y[1]];
    else #  bestTrajType == "rlr"
      p = rlr_rl_tangent;
    end

    phi_start = atan(initial_location[2]-irc_center[2], initial_location[1]-irc_center[1]);
    phi_end = atan(p[2]-irc_center[2], p[1]-irc_center[1]);

    if phi_end > phi_start
        phi_end = phi_end - 2.0*pi;
    end

    phis = (phi_end == phi_start) ? phi_start : collect(phi_start:-delta_phi:phi_end);

    first_path_x = irc_center[1] .+ r_min*cos.(phis);
    first_path_y = irc_center[2] .+ r_min*sin.(phis);

  elseif bestTrajType[1] == 'l'
    # first part of path is a left hand turn
    if bestTrajType == "lsl"
      p = [lsl_tangent_x[1], lsl_tangent_y[1]];
    elseif bestTrajType == "lsr"
      p = [lsr_tangent_x[1], lsr_tangent_y[1]];
    else #  bestTrajType == "rlr"
      p = lrl_lr_tangent;
    end

    phi_start = atan(initial_location[2]-ilc_center[2], initial_location[1]-ilc_center[1]);

    phi_end = atan(p[2]-ilc_center[2], p[1]-ilc_center[1]);

    if phi_end < phi_start
        phi_end = phi_end + 2.0*pi;
    end

    phis = (phi_end == phi_start) ? phi_start : collect(phi_start:delta_phi:phi_end);
    first_path_x = ilc_center[1] .+ r_min*cos.(phis);
    first_path_y = ilc_center[2] .+ r_min*sin.(phis);
  end


  # calculate second part of path
  if bestTrajType[2] == 's'
    # second part of path is a straight line

    if bestTrajType == "lsr"
      p1 = [lsr_tangent_x[1], lsr_tangent_y[1]];
      p2 = [lsr_tangent_x[2], lsr_tangent_y[2]];
    elseif bestTrajType == "lsl"
      p1 = [lsl_tangent_x[1], lsl_tangent_y[1]];
      p2 = [lsl_tangent_x[2], lsl_tangent_y[2]];
    elseif bestTrajType == "rsr"
      p1 = [rsr_tangent_x[1], rsr_tangent_y[1]];
      p2 = [rsr_tangent_x[2], rsr_tangent_y[2]];
    else # bestTrajType == "rsl"
      p1 = [rsl_tangent_x[1], rsl_tangent_y[1]];
      p2 = [rsl_tangent_x[2], rsl_tangent_y[2]];
    end

    second_path_x = [p1[1], p2[1]];
    second_path_y = [p1[2], p2[2]];

  elseif bestTrajType[2] == 'r'
    # second part of path is a right turn

    phi_start = atan(lrl_lr_tangent[2]-lrl_r_circle_center[2], lrl_lr_tangent[1]-lrl_r_circle_center[1]);
    phi_end = atan(lrl_rl_tangent[2]-lrl_r_circle_center[2], lrl_rl_tangent[1]-lrl_r_circle_center[1]);

    if phi_end > phi_start
      phi_end = phi_end - 2.0*pi;
    end

    phis = (phi_end == phi_start) ? phi_start : collect(phi_start:-delta_phi:phi_end);

    second_path_x = lrl_r_circle_center[1] .+ r_min*cos.(phis);
    second_path_y = lrl_r_circle_center[2] .+ r_min*sin.(phis);

  elseif bestTrajType[2] == 'l'
    # second part of path is a left turn

    phi_start = atan(rlr_rl_tangent[2]-rlr_l_circle_center[2], rlr_rl_tangent[1]-rlr_l_circle_center[1]);
    phi_end = atan(rlr_lr_tangent[2]-rlr_l_circle_center[2], rlr_lr_tangent[1]-rlr_l_circle_center[1]);

    if phi_end < phi_start
      phi_end = phi_end + 2.0*pi;
    end

    phis = (phi_end == phi_start) ? phi_start : collect(phi_start:delta_phi:phi_end);

    second_path_x = rlr_l_circle_center[1] .+ r_min*cos.(phis);
    second_path_y = rlr_l_circle_center[2] .+ r_min*sin.(phis);
  end

  # calculate third part of path
  if bestTrajType[3] == 'r'
    # third part of path is a right hand turn

    if bestTrajType == "rsr"
      p = [rsr_tangent_x[2], rsr_tangent_y[2]];
    elseif bestTrajType == "lsr"
      p = [lsr_tangent_x[2], lsr_tangent_y[2]];
    else #  bestTrajType == "rlr"
      p = rlr_lr_tangent;
    end

    phi_start = atan(p[2]-grc_center[2], p[1]-grc_center[1]);
    phi_end = atan(goal_location[2]-grc_center[2], goal_location[1]-grc_center[1]);

    if phi_end > phi_start
        phi_end = phi_end - 2.0*pi;
    end

    phis = (phi_end == phi_start) ? phi_start : collect(phi_start:-delta_phi:phi_end);

    third_path_x = grc_center[1] .+ r_min*cos.(phis);
    third_path_y = grc_center[2] .+ r_min*sin.(phis);
  elseif bestTrajType[3] == 'l'
    # third part of path is a left hand turn

    if bestTrajType == "lsl"
      p = [lsl_tangent_x[2] lsl_tangent_y[2]];
    elseif bestTrajType == "rsl"
      p = [rsl_tangent_x[2] rsl_tangent_y[2]];
    else #  bestTrajType == "lrl"
      p = lrl_rl_tangent;
    end

    phi_start = atan(p[2]-glc_center[2], p[1]-glc_center[1]);
    phi_end = atan(goal_location[2]-glc_center[2], goal_location[1]-glc_center[1]);

    if phi_end < phi_start
        phi_end = phi_end + 2.0*pi;
    end

    phis = (phi_end == phi_start) ? phi_start : collect(phi_start:delta_phi:phi_end);

    third_path_x = glc_center[1] .+ r_min*cos.(phis);
    third_path_y = glc_center[2] .+ r_min*sin.(phis);
  end


  edge.dubinsType = bestTrajType
  edge.Wdist = bestDist          # distance that the robot moves in workspace

  if edge.Wdist == Inf
    edge.dist = Inf
  elseif S.spaceHasTime
    # calcualte C-space edge length
    # NOTE THIS HARDCODED batch-version of dubinsDistAlongPath only works because
    # we assume constant speed along the edge
    edge.dist = sqrt(bestDist^2 + (edge.startNode.position[3] - edge.endNode.position[3])^2)

    # we need to calculate the time paramiterization for the robot along the path.
    # NOTE: we make the simplifying assumption that the robot travels at constant
    # velocity along the entire path. FURTHERMORE, the total time that it requires
    # to traverse the trajectory is determined by the difference in time positions
    # of the end-points. ALSO NOTE: This function is not responsible for determining
    # if it is possible for the robot to actually achieve this speed---which is done
    # in the function validMove(). FINALLY, I have also assumed that the trajectory
    # is sampled well enough in "curvy" parts that the distance between points
    # can be approximated by the straight-line distance between these points
    # i.e., the sum of the segment lengths is close enough to edge.Wdist that
    # problems will not occour

    edge.velocity = edge.Wdist/(edge.startNode.position[3] - edge.endNode.position[3])

    # build trajectory with 0 for times
    edge.trajectory = [[first_path_x  first_path_y zeros(size(first_path_x))  ];
                       [second_path_x second_path_y zeros(size(second_path_x))];
                       [third_path_x third_path_y zeros(size(third_path_x))   ] ]

    # now calculate times
    edge.trajectory[1,3] = edge.startNode.position[3]
    cumulativeDist = 0.0
    for i = 2:(size(edge.trajectory,1)-1)
      cumulativeDist += Wdist(edge.trajectory[i-1,1:2], edge.trajectory[i,1:2])
      edge.trajectory[i,3] = edge.startNode.position[3] - cumulativeDist/edge.velocity
    end
    edge.trajectory[end,1:3] = edge.endNode.position[1:3] # make end point exact

  else
    edge.dist = bestDist
    edge.trajectory = [[first_path_x  first_path_y ];
                       [second_path_x second_path_y];
                       [third_path_x third_path_y]   ]

  end

  edge.distOriginal = edge.dist

  # for debugging only:
  #saveData(edge.trajectory, "temp/DubinsTrajectory_$(S.fileCtr).txt")
end


### this calculates a trajectory of what the robot is supposed to do when it is
### hovering "in place"
#
## Dubins edge version
function calculateHoverTrajectory(S::TS, edge::DubinsEdge) where {TS}
  edge.dubinsType = "xxx"
  edge.Wdist = 0.0
  edge.dist = 0.0

  if S.spaceHasTime
    edge.velocity = S.dubinsMinVelocity

    edge.trajectory = vcat(edge.startNode.position[1,1:3]', edge.endNode.position[1,1:3]')
  else
    edge.trajectory = vcat(edge.startNode.position[1,1:2]', edge.endNode.position[1,1:2]')
  end
end


### this saves the trajectory that is stored in the edge to the file
#
# Dubins edge version
function saveEdgeTrajectory(filePtr::TFP, edge::DubinsEdge) where {TFP}
  writedlm(filePtr, edge.trajectory, ',')
end



########################## collision checking functions ###########################
### Note, these are only collision checking functions that depend on edge type, ###
### more general collision checking functions appear in DRRT.jl                 ###
###################################################################################


### checks if the edge is in collision with a particular obstacle
### (one version for each edge type)
#
## DubinsEdgeVersion
function explicitEdgeCheck(S::CSpace{T}, edge::DubinsEdge, obstacle::OT) where {T, OT}
  # We know that all points in the dubin's trajectory are within
  # r*S.minTurningRadius of the 2D segment from startNode to endNode,
  # thus, the edge is safe if a robot with this additional radius can traverse
  # that segment


  if !explicitEdgeCheck2D(obstacle, edge.startNode.position, edge.endNode.position,
    S.robotRadius + 2*S.minTurningRadius)
    return false
  end

  # if the segment was in conflict then we need to do the full check of the
  # trajectory segments. In the future this could be improved by making a
  # function that can check arcs of the dubin's path all at once instead
  # of just the line segments currently stored in the trajectory.

  for i = 2:size(edge.trajectory,1)
    if explicitEdgeCheck2D(obstacle, edge.trajectory[i-1,:], edge.trajectory[i,:], S.robotRadius)
      return true
    end
  end

  return false
end
