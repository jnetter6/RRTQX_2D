function make_nice()

end


# saves the data to the file
function saveData(data, fileName)
  fptr = open(fileName, "w")
  for i = 1:size(data,1)
    writedlm(fptr, data[i,:], ',')
  end
  close(fptr)
end

# This solves the 1D two point bouandary value problem, given the start and
# end positions and velocites at those positions, this returns the sequence
# of accelerations necessary to make the system achieve that, as well as a
# resulting positions and velocites. a is the bound on acceleration magnitude
# so acceleration can be -a, a, or 0. in this simplified version we assume
# that there is one  acceleration at the begining and another at the end,
# with a possible coast in the middle. returns nans if cannot find a solution
# otherwise returns [t_1 t_2 a_1 a_2 v_coast]
function steering_1D(x_start, x_end, v_start, v_end, t_start, t_end, a)

  # we assume that the starting velocity is positive, if not then just need to
  # flip signs on all x and v and then solve, and then flip signs back

  flip_flag = false
  if v_start < 0
    # flip direction, solve and then flip back
    flip_flag = true

    x_start = -x_start
    x_end = -x_end
    v_start = -v_start
    v_end = -v_end
  end

  t_1 = t_2 = a_1 = a_2 = v_coast = 0.0 # all these are reset later


  delta_x = x_end - x_start
  delta_t = t_end - t_start

  # solve everything for between 0 and delta_t, then add in


  # first make sure that we have enough time to do the required velocity
  # change
  t_hat = abs(v_end - v_start)/a

  if t_hat > delta_t
    # warning('not enough time to change velocity')
    return fill(NaN,1,5)
  end

  v_max = (delta_t*a + v_start + v_end)/2
  v_min =  (v_start + v_end - delta_t*a)/2

  tau_1 = (v_start^2)/(a*2)   # distance to v = 0 from start
  tau_2 = (v_end^2)/(a*2)     # magnitude of distance from v = 0 to goal

  t_a = v_start/a             # time to zero v from start

  if v_end >= 0
    # overall case A: both start and end velocity is positive

    t_b = delta_t - v_end/a     # -time to zero v from goal


    # Figure out which one of the sub-cases we have,
    # first few are for non-overlapping tau tringles
    # second few are for overlapping tau triangles
    if t_a <= t_b


      # first calculte the max and min positions that can be reached at t_end
      # to make sure that a solution exists.

      delta_x_max = (delta_t + t_a + t_end-t_b)*v_max/2 - tau_1 - tau_2

      if x_start + delta_x_max < x_end
        # warning('cannot be solved: end is too high to get to in time')
        return fill(NaN,1,5)
      end

      delta_x_min = (delta_t - (delta_t-t_b) - t_a)*v_min/2 + tau_1 + tau_2

      if x_start + delta_x_min > x_end
        # warning('cannot be solved: end is too low to get to in time')
        return fill(NaN,1,5)
      end

      if delta_x < tau_1 + tau_2
        # sub-case 1, we need negative velocity at some point

        z = sqrt((delta_x - delta_x_min)*a)
        t_delta_min_v = (t_a+t_b)/2

        v_coast = v_min + z        # velocity to coast at
        t_1 = t_delta_min_v - z/a  # time at start of coast
        t_2 = t_delta_min_v + z/a  # time at end of coast

      elseif delta_x < delta_t*min(v_start, v_end) + ((v_start - v_end)^2)/(2*a)
        # sub-case 2, reduce velocity, wait, increase velocity

        z = sqrt((delta_x - delta_x_min)*a)
        t_delta_min_v = (t_a+t_b)/2

        v_coast = v_min + z        # velocity to coast at
        t_1 = t_delta_min_v - z/a  # time at start of coast
        t_2 = t_delta_min_v + z/a  # time at end of coast

      elseif delta_x < delta_t*max(v_start, v_end) - ((v_start - v_end)^2)/(2*a)
        # sub-case 3, reduce velocity, wait, reduce velocity OR
        #             increase velocity, wait, increase velocity

        if tau_1 < tau_2

          v_coast = (delta_x + tau_1 - tau_2)/(t_b + t_a)
          t_1 = (v_coast - v_start)/a
          t_2 = t_b + v_coast/a

        else

          v_coast = (delta_x - tau_1 + tau_2)/(2*delta_t - t_a - t_b)
          t_1 = (v_start - v_coast)/a
          t_2 = delta_t - (v_coast - v_end)/a
        end
      else
        # subcase 4, increase velocity, wait, decrease velocity

        z = sqrt((delta_x_max - delta_x)*a)
        t_delta_max_v = (v_max-v_start)/a

         v_coast = v_max - z         # velocity to coast at
         t_1 = t_delta_max_v - z/a   # time at start of coast
         t_2 = t_delta_max_v + z/a   # time at end of coast

      end
    else # t_b < t_a

      t_v_max = (v_max - v_start)/a
      delta_x_max = v_max*delta_t - (v_max-v_start)^2/(2*a) - (v_max-v_end)^2/(2*a)
      if delta_x_max < delta_x
        # warning('cannot solve: end is to high to get to in time')
        return fill(NaN,1,5)
      end

      delta_x_min = -(delta_t - (delta_t-t_b) - t_a)*v_min/2 + tau_1 + tau_2
      tau_3 = v_min^2/a

      if delta_x_min > delta_x
        # warning('cannot solve: end is to low to get to in time')
        return fill(NaN,1,5)
      end

      if delta_x < tau_1 + tau_2 - tau_3
        # sub-case 1, we need negative velocity at some point
        # warning('cannot solve: need negative velocity and cannot attain it')

        return fill(NaN,1,5)
      elseif (tau_1 < tau_2 && delta_x <= tau_1 + tau_2 + t_b*v_start) ||
             (tau_2 <= tau_1 && delta_x <= tau_1 + tau_2 + (delta_t-t_a)*v_end)
        # sub-case 2 reduce velocity, wait, increase velocity
        # but final velocity is more

        z= sqrt((delta_x-tau_1-tau_2+tau_3)*a)

        v_coast = (t_a - t_b)/2*a + z
        t_1 = (t_a+t_b)/2 - z/a
        t_2 = (t_a+t_b)/2 + z/a

      elseif delta_x < delta_t*max(v_end,v_start) - ((v_end-v_start)^2)/(2*a)
        # sub-case 3, reduce velocity, wait, reduce velocity OR
        #             increase velocity, wait, increase velocity

        if tau_1 < tau_2

          x_i = t_b + v_start/a
          z = (delta_x - tau_1 - tau_2 - v_start*t_b)/x_i

          v_coast = v_start + z
          t_1 = z/a
          t_2 = z/a + x_i

        else
          x_i = delta_t - (t_a - v_end/a)
          z = (delta_x - tau_1 - tau_2 - (delta_t-t_a)*v_end)/x_i

          v_coast = v_end + z
          t_1 = delta_t - z/a - x_i
          t_2 = delta_t - z/a
        end
      else
        # subcase 4, increase velocity, wait, decrease velocity

        z = sqrt((delta_x_max-delta_x)*a)

        v_coast = v_max - z
        t_1 = t_v_max - z/a
        t_2 = t_v_max + z/a


      end
    end

  else
    # overall case B: start velocity is positive, and end velocity is negative

    t_b = delta_t + v_end/a     # -time to zero v from goal

    # first calculte the max and min positions that can be reached at t_end
    # to make sure that a solution exists.

    delta_x_max = (t_a + t_b)*v_max/2 - tau_1 - tau_2

    if x_start + delta_x_max < x_end
      # warning('cannot be solved: end is too high to get to in time')
      return fill(NaN,1,5)
    end

    #delta_x_min = ((v_min + v_end)*(v_end - v_min) - v_min^2 + v_start^2)/(2*a)
    delta_x_min = (2*delta_t - t_b -t_a)*v_min/2 + tau_1 + tau_2

    if x_start + delta_x_min > x_end
      # warning('cannot be solved: end is too low to get to in time')
      return fill(NaN,1,5)
    end

    # the problem can be solved, so figure out which one of the 4 basic
    # sub-cases we need to solve

    if delta_x < tau_1 - tau_2 + v_end*(t_b-t_a)
      # sub-case 1: decelerate to negative velocity wait, then accelerate
      # (to still negative velocity)

      z = sqrt((delta_x - delta_x_min)*a)
      t_delta_min_v = (v_start-v_min)/a

      v_coast = v_min + z        # velocity to coast at
      t_1 = t_delta_min_v - z/a  # time at start of coast
      t_2 = t_delta_min_v + z/a  # time at end of coast

    elseif delta_x < tau_1 - tau_2 + v_start*(delta_t+(v_end-v_start)/a)
      # sub-case 2 and 3: decelerate to negative/positive velocity, wait, then
      # decelerate more

      v_coast = (delta_x - tau_1 + tau_2)/(t_b - t_a) # velocity to coast at
      t_1 = t_a - v_coast/a       # time at start of coast
      t_2 = t_b - v_coast/a       # time at end of coast

    else
      # sub-case 4: accelerate, wait, decelerate

      z = sqrt((delta_x_max - delta_x)*a)
      t_delta_max_v = (v_max-v_start)/a

      v_coast = v_max - z         # velocity to coast at
      t_1 = t_delta_max_v - z/a   # time at start of coast
      t_2 = t_delta_max_v + z/a   # time at end of coast
    end
  end

  if flip_flag == true
    # flip back

    flip_flag = true

    x_start = -x_start
    x_end = -x_end
    v_start = -v_start
    v_end = -v_end

    v_coast = -v_coast
  end


  # note:
  # a_1 is acceleration from t_start to t_1
  # a_2 is acceleration from t_2 to t_end
  # (acceleration from t_1 to t_2 is 0)
  # v_coast is the velocity we coast at from t_1 to t_2

  a_1 = (v_coast - v_start)/t_1
  a_2 = (v_end - v_coast)/(delta_t - t_2)

  # add back in start time
  t_1 += t_start
  t_2 += t_start

  return [t_1 t_2 a_1 a_2 v_coast]
end


# solves the ND steering function for the two point bouandary value problem
# all inputs are vectors of length D. a contaisn max acceleration magnitudes per
# dimension. x is position, v is velocity. This returns a vector of times and
# three matricies, one each for position, velocity, and acceleration at those times
# maticies are (t X D). Note that acceleration(i,:) gives the acceleration from time
# i to time i+1, while velocity and position give the exact velocities and positions
# a each time. This also returns a flag indicating (true) if the solve is successful
function steering_ND(x_start, x_end, v_start, v_end, t_start, t_end, a)
  # find dimensionality of space
  D = length(x_start)

  # first solve the 1D problem in each dimension (a is between t, while
  # v, x are at t), each D gets 4 times including start/goal that things
  # change
  raw_t = fill(NaN,2,D)  # note that first and last time are the same for all
                         # D and so left out of here
  raw_a = fill(NaN,3,D)

  for d = 1:D
    rets1D = steering_1D(x_start[d], x_end[d], v_start[d], v_end[d], t_start, t_end, a[d])

    # Note: rets1D = [t_1 t_2 a_1 a_2 v_coast]

    println("$(d)  ret   $(rets1D)")

    if isnan(rets1D[1])
      # this dimension is impossible to solve
      # return (false, Array(Float64,0,0), Array(Float64,0,0), Array(Float64,0,0), Array(Float64,0,0))
      return (false, Array{Float64,2}(undef,0,0), Array{Float64,2}(undef,0,0), Array{Float64,2}(undef,0,0), Array{Float64,2}(undef,0,0))
    end

    raw_t[:,d] = rets1D[1:2]
    raw_a[:,d] = [rets1D[3] 0.0 rets1D[4]]
  end

  # interleave different times (note that if this is ever used for
  # reallyhigh D systems, then this should be rewritten to use sorting
  # first). This next part makes times consistant per row, which will
  # rewuire adding more points

  # we will have D*2 + 2 timesteps max
  Time = fill(NaN,D*2+2,1)
  A = fill(NaN,D*2+1,D)  # A between, while V X at a particular Time
  V = fill(NaN,D*2+2,D)
  X = fill(NaN,D*2+2,D)
  # Time = Array{Float64,2}(undef,D*2+2,1)
  # A = Array{Float64,2}(undef,D*2+1,D)  # A between, while V X at a particular Time
  # V = Array{Float64,2}(undef,D*2+2,D)
  # X = Array{Float64,2}(undef,D*2+2,D)

  # do times first and accelerations first
  t_i = fill(1, 1, D)  # index into time along each dimension

  # Time[1,:] = t_start
  # Time[end,:] = t_end
  Time[1] = t_start
  Time[end] = t_end
  A[1,:] = raw_a[1,:]

  println("raw_t: $(raw_t)")

  for t = 2:D*2+1
    # find the next most recent time
    next_time = Inf
    cooresponding_D = -1

    for j = 1:D
      if t_i[j] < 3 && raw_t[t_i[j],j] < next_time
        next_time = raw_t[t_i[j],j]
        cooresponding_D = j
      end
    end

    # now put this next time into times
    if cooresponding_D == -1
      # at t_end
      Time[t] = t_end
    else
      Time[t] = raw_t[t_i[cooresponding_D],cooresponding_D]
      t_i[cooresponding_D] += 1
    end

    # take care of accelerations up to this time
    for j = 1:D
      A[t,j] = raw_a[t_i[j],j]
    end
  end

  # now calculate V and X based on A
  V[1,:] = v_start[:]
  X[1,:] = x_start[:]
  for t = 2:D*2+2
    Dt = (Time[t,:] .- Time[t-1,:])
    V[t,:] = V[t-1,:] .+ (Dt .* A[t-1,:])
    X[t,:] = X[t-1,:] .+ (Dt .* V[t-1,:]) .+ (Dt.^2 .* A[t-1,:] / 2.0)
  end

  return (true, Time, X, V, A)
end


# given T, A, V, X this calculates intermediate values at times at a resolution
# of dt
function fineGrain(Time, A, V, X, dt)

  dims = size(A, 2)
  Dt = Time[end] - Time[1]
  num_t = Int(div(Dt,dt)) + length(Time) + 1

  Time_fine = fill(NaN, num_t, 1)
  A_fine = fill(NaN, num_t, dims)
  V_fine = fill(NaN, num_t, dims)
  X_fine = fill(NaN, num_t, dims)

  i = 0;
  for t = 1:length(Time)-1
    i += 1
    Time_fine[i] = Time[t]
    A_fine[i,:] = A[t,:]
    V_fine[i,:] = V[t,:]
    X_fine[i,:] = X[t,:]

    while Time_fine[i] + dt <= Time[t+1]
      i += 1

      Time_fine[i] = Time_fine[i-1] + dt
      A_fine[i,:] = A[t,:]
      V_fine[i,:] = V_fine[i-1,:] .+ (dt * A_fine[i,:])
      X_fine[i,:] = X_fine[i-1,:] .+ (dt * V_fine[i-1,:]) .+ ((dt^2)/2.0 .* A_fine[i-1,:])
    end
  end

  Time_fine[end] = Time[end]
  A_fine[end,:] = A[end,:]
  V_fine[end,:] = V[end,:]
  X_fine[end,:] = X[end,:]

  return (Time_fine, A_fine, V_fine, X_fine)

end


# function for testing stuff in this file test
function test()

#  x_start = rand()
#  x_end  = rand()
#  v_start  = rand()
#  v_end  = rand()
#  t_start  = 0.0
#  t_end  = 10.0
#  a = rand()

#  ret = steering_1D(x_start, x_end, v_start, v_end, t_start, t_end, a)
#  println(ret)

  D = 3

  x_start = rand(1,D)
  x_end  = rand(1,D)
  v_start  = rand(1,D)
  v_end  = rand(1,D)
  t_start  = 0.0
  t_end  = 10.0
  a = rand(1,D)

  (foundSolution, Time, X, V, A) = steering_ND(x_start, x_end, v_start, v_end, t_start, t_end, a)

  if !foundSolution
    println("did not find solution")
    return
  end

  println("t:")
  println(Time)

  println("a:")
  println(A)

  println("v:")
  println(X)

  println("x:")
  println(V)


  saveData(Time, "temp/times_raw.txt")
  saveData(A, "temp/a_raw.txt")
  saveData(V, "temp/v_raw.txt")
  saveData(X, "temp/x_raw.txt")

  (Time_fine, A_fine, V_fine, X_fine) = fineGrain(Time, A, V, X, .01)

  saveData(Time_fine, "temp/times_fine.txt")
  saveData(A_fine, "temp/a_fine.txt")
  saveData(V_fine, "temp/v_fine.txt")
  saveData(X_fine, "temp/x_fine.txt")

end
