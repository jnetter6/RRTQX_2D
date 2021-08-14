function makeMyTextEditorDisplayNice4() # thanks!
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


########################## basic Distance Functions ###############################
### A number of basic distance functions. These are not called direcly unless   ###
### as part of code for a space in which the specific type of distance is known.###
### Instead, they are usually aliased to have other names depending on the      ###
### space and robot being used.                                                 ###
###################################################################################

# returns the Euclidian distance
euclidianDist(x::Array,y::Array) = sqrt(sum((x-y).^2))

# returns the "straight-line" distance in a space containing [X Y Time Theta]
# where theta exists on [0 2*pi] and wrapps around, i.e., 0 === 2*pi
R3SDist(x::Array,y::Array) = (sqrt( sum((x[1:3]-y[1:3]).^2) + min(abs(x[4]-y[4]), min(x[4], y[4])+2.0*pi-max(x[4], y[4]))^2 ))

# returns the Distance between two points in the projection of a dubins space
# that contains [X Y] depending on if time is being used or not
# this is useful for calculating the distance between two points that are close
# to each other on the same dubin's path, e.g. when change in theta is small
# and so the Theta component can be ignored
dubinsDistAlongPath(x::Array,y::Array)  = sqrt(sum((x[1:2]-y[1:2]).^2))

# returns the Distance between two points in the projection of a dubins space
# that contains [X Y T] depending on if time is being used or not
# this is useful for calculating the distance between two points that are close
# to each other on the same dubin's path, e.g. when change in theta is small
# and so the Theta component can be ignored
dubinsDistAlongTimePath(x::Array,y::Array)  = sqrt(sum((x[1:3]-y[1:3]).^2))


# helps with dubins car
# returns the distance that the car travels to get from
# pointA to pointB, where both are asumed to be on the circle of radoius r
# that is centered at circleCenter
function rightTurnDist(pointA, pointB, circleCenter, r)
    theta = atan(pointA[2]-circleCenter[2], pointA[1]-circleCenter[1]) - atan(pointB[2]-circleCenter[2], pointB[1]-circleCenter[1]);
    if theta < 0
       theta = theta + 2*pi;
    end
    return theta*r;
end

# helps with dubins car
# returns the distance that the car travels to get from
# pointA to pointB, where both are asumed to be on the circle of radoius r
# that is centered at circleCenter
function leftTurnDist(pointA, pointB, circleCenter, r)
    theta = atan(pointB[2]-circleCenter[2], pointB[1]-circleCenter[1]) - atan(pointA[2]-circleCenter[2], pointA[1]-circleCenter[1]);
    if theta < 0
       theta = theta + 2*pi;
    end
    return theta*r;
end
