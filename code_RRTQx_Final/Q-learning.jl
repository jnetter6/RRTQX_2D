### this contains functions used for applying Q-learning to TPBVP between nodes
using LinearAlgebra
using ControlSystems
using DifferentialEquations
using Interpolations
using Plots

function sim_TNNLS_A_infinite_LQR()
	n, m = 3, 1
	A = [-1.01887 0.90506 -0.00215; 0.82225 -1.07741 -0.17555; 0 0 -1] # n x n
	B = reshape([0.0; 0.0; 1.0], (3, 1)) # n x m
	M = 1.0*Matrix(I, n, n) # n x n
	R = 0.1*Matrix(I, m, m) # m x m
	L = lqr(A, B, M, R)
	println(L)
end

# infinite simulation A in TNNLS with F-16 model
function sim_TNNLS_A_infinite(NN::Int)#,amplitude, percent) # NN is the # of iterations
	# start = time()
    n, m = 3, 1
    A = [-1.01887 0.90506 -0.00215; 0.82225 -1.07741 -0.17555; 0 0 -1] # n x n
    B = reshape([0.0; 0.0; 1.0], (3, 1)) # n x m
    M = 1.0*Matrix(I, n, n) # n x n
    R = 0.1*Matrix(I, m, m) # m x m

	# check controllability & observability
    Co = ctrb(A, B)
    unco = size(A, 1) - rank(Co)
    # println("unco = $unco")
    Ob = obsv(sqrt(M), A)
    unob = size(sqrt(M), 1) - rank(Ob)
    # println("unob = $unob")

	if unco + unob > 0
		println("system uncontrollable and/or unobservable")
		return
	end

    # ODE parameters
    Tf, T, N = 45, 0.05, 900 # finite horizon

    αc, αa = 90, 2.5

	amplitude, percent = 0.05, 50

    Wc0 = [10.0*ones(Int((n + m)*(n + m + 1)/2 - m*m), 1); reshape(R, (m*m, 1))][:] # (n+m)(n+m+1)/2 x 1; 10 x 1
	Wa0 = 0.5*ones(n,) # n x m; 3 x 1
    u0 = 5.0*ones(m,) # m x 1; 1 x 1
    uf = 0.001

    x0 = [1.0, 5.0, 1.0]
    xf = [2.0, 7.0, 3.0]

	Quu = 0
	Qxu = 0

    p0 = x0'*M*x0

    t_save = [0,]
    u_save = [u0; u0]
    # x_save = zeros(17, 1) # 3+10+3+1
	x_save = [[x0; Wc0; Wa0; p0],]
    # println("x_save = $x_save")

    uDelay  = interp2PWC(u_save[:], -1, 1) # return an interpolation function
	x1Delay = interp2PWC(getindex.(x_save, 1), -1, 1)
    x2Delay = interp2PWC(getindex.(x_save, 2), -1, 1)
    x3Delay = interp2PWC(getindex.(x_save, 3), -1, 1)

	uvec = 0
	xdist = norm(x0 - xf)
	error = 0
    # solve ODEs
	for i = 1:NN
		start = time()
		# global xf, uf, Tf, T, A, B, M, R, Pt, per, amp, αa, αc, Quu, Qxu, Wcfinal, uDelay, x1Delay, x2Delay, x3Delay, x_save, t_save, uvec
		# uvec = 1
		t = ((i-1)*T, i*T)
		# t = (0, 1)
		p = [xf, uf, Tf, T, A, B, M, R, percent, amplitude, αa, αc, Quu, Qxu, uDelay, x1Delay, x2Delay, x3Delay, x_save, t_save, uvec, i,NN,start]
		# if i == NN println(1) end
        # sol = solve(ODEProblem(babyETC_A, x_save[end], t, p), DP5())
		sol = solve(ODEProblem(babyCT_A_infinite, x_save[end], t, p), DP5())
		# if i == NN println(2) end
		t_save = [t_save; sol.t[2:end]] # vcat(t_save, sol.t) save time and kick out the duplicated initial value
		x_save = [x_save; sol.u[2:end]] # vcat(x_save, sol.u) save state
		# if i == NN println(3) end
		# println("sol.u[$i] = $(sol.u)")
		# println("sol.t[$i] = $(sol.t)")
		uvec = p[end-3]
		# uvec = p[end]
		# println("uvec[$i] = $(p[end])")
		# if i == NN println(4) end
        uDelay  = interp2PWC(uvec, -1, i*T+.01) # interpolate control input
		x1Delay = interp2PWC(getindex.(x_save, 1), -1, i*T+.01)
		x2Delay = interp2PWC(getindex.(x_save, 2), -1, i*T+.01)
		x3Delay = interp2PWC(getindex.(x_save, 3), -1, i*T+.01)
		# println(i)
		# if euclidean(x_save[end][1:3], xf) < error*xdist break end
		if norm(x_save[end][1:3] - xf) < error*xdist break end
	end

	x_save_1 = getindex.(x_save, 1)
	# x_save_2 = getindex.(x_save, 2)
	# x_save_3 = getindex.(x_save, 3)

	# println("t_save = $t_save")
	# println("length_t_save = $(length(t_save))")
	# println("length_x_save = $(length(x_save))")
	# println("x1 = $(x_save_1)")
	# println("x2 = $(x_save_2[end])")
	# println("x3 = $(x_save_3[end])")
	println("x1 = $(x_save[end][1])")
	println("x2 = $(x_save[end][2])")
	println("x3 = $(x_save[end][3])")
	plot(t_save, getindex.(x_save, 1), w = 2)
	plot!(t_save, getindex.(x_save, 2))#, lw = 2)
	plot!(t_save, getindex.(x_save, 3))#, lw = 2)
	# println("Wa = $(x_save[end][14:16])")
	# println(x_save[end])
	# elapsed = time() - start
end

# ode fcn for infinite simulation A in TNNLS
function babyCT_A_infinite(dotx, x, p, t)

	# global xf, uf, Tf, T, A, B, M, R, Pt, per, amp, αa, αc, Quu, Qxu, uDelay, x1Delay, x2Delay, x3Delay, x_save, t_save, uvec
	xf, uf, Tf, T, A, B, M, R, percent, amplitude, αa, αc, Quu, Qxu, uDelay, x1Delay, x2Delay, x3Delay, x_save, t_save, uu,i,NN,start = p # uu useless here
	# println(t[1][1])
    n, m = 3, 1
	Wc = x[(n+1) : Int(n+(n+m)*(n+m+1)/2)] # 4-13
    Wa = x[Int(n+(n+m)*(n+m+1)/2+1) : Int(n+(n+m)*(n+m+1)/2+n)] # 14-16
    # P = x[end] # 17
	P = x[end]

	# Update control
    ud = zeros(m,)
	ud[1] = Wa'*(x[1:n] - xf)

	# Delays
    uddelay = zeros(m,)
	uddelay[1] = uDelay(t - T)
    xdelay = zeros(n,)
	xdelay[1] = x1Delay(t - T)
    xdelay[2] = x2Delay(t - T)
    xdelay[3] = x3Delay(t - T)
	# Kronecker products
	U = vcat(x[1:n] - xf, ud) # augmented state
    UkU = vcat(U[1]^2, U[1]*U[2], U[1]*U[3], U[1]*ud,
               U[2]^2, U[2]*U[3], U[2]*ud,
               U[3]^2, U[3]*ud,
               ud[1]^2)
    UkUdelay = vcat(xdelay[1]^2, xdelay[1]*xdelay[2], xdelay[1]*xdelay[3], xdelay[1]*uddelay,
                    xdelay[2]^2, xdelay[2]*xdelay[3], xdelay[2]*uddelay,
                    xdelay[3]^2, xdelay[3]*uddelay,
                    uddelay[1]^2)
    Quu = Wc[end] # m x m
    Quu_inv = inv(Quu)
    Qxu = reshape(Wc[Int((n+m)*(n+m+1)/2-m^2-m*n+1):Int((n+m)*(n+m+1)/2-m^2)], (n, m)) # n x m
    Qux = Qxu' # m x n

    σ = UkU - UkUdelay

	# Integral reinforcement dynamics
	dP = 0.5*((x[1:n] - xf)'*M*(x[1:n] - xf) - xdelay'*M*xdelay + ud'*R*ud - uddelay'*R*uddelay)

	# approximation errors
    ec = P + Wc'*UkU - Wc'*UkUdelay
    ea = Wa'*U[1:n] + [Quu_inv*Qux[:,1]; Quu_inv*Qux[:,2]; Quu_inv*Qux[:,3]]'*U[1:n]

	# Critic dynamics
	dWc = -αc*(σ./(σ'*σ+1)^2)*ec

	# Actor dynamics
    dWa = -αa*U[1:n]*ea'

	# Persistence excitation
    unew = zeros(m,)
	# if t <= (p[10]/100)*p.Tf
	if t <= (percent/100)*Tf
	    # unew[1] = ud[1] + amplitude*exp(-.0005*t)*(sin(0.1*t)^2 + cos(0.7*t)^2) + sin(1.5*t)^2*cos(0.1*t) + sin(pi*t) + cos(0.1*t)
	    # unew[1] = ud[1] + p[11]*exp(-.0005*t)*(sin(0.1*t)^2+cos(0.7*t)^2)+sin(1.5*t)^2*cos(0.1*t)+sin(pi*t)+cos(0.1*t)
	    unew[1] = (ud[1]+amplitude*exp(-0.009*t)*2*(sin(t)^2*cos(t)+sin(2*t)^2*cos(0.1*t)+sin(-1.2*t)^2*cos(0.5*t)+sin(t)^5+sin(1.12*t)^2+cos(2.4*t)*sin(2.4*t)^3))
	    # unew[2] = (ud[2]+amplitude*exp(-0.009*t)*2*(sin(t)^2*cos(t)+sin(2*t)^2*cos(0.1*t)+sin(-1.2*t)^2*cos(0.5*t)+sin(t)^5+sin(1.12*t)^2+cos(2.4*t)*sin(2.4*t)^3))
	else
	    unew = ud
	end
	dx = A*U[1:n] + B*unew

	uu = unew
	# p[end-3] = unew
	# Augmented state
	# for j in 1:n dotx[j] = dx[j] end
	# for j in (n+1):convert(Int64, (n+(n+m)*(n+m+1)/2)) dotx[j] = dWc[j-n] end
	# for j in convert(Int64, (n+(n+m)*(n+m+1)/2+1)):convert(Int64, (n+(n+m)*(n+m+1)/2+n)) dotx[j] = dWa[j-convert(Int64, (n+(n+m)*(n+m+1)/2))] end
	# dotx[end] = dP
	if (time() - start) > 2
		# return dotx .= vcat(dx, dWc, dWa, dP)
		return
	end
	dotx .= vcat(dx, dWc, dWa, dP)
	# println("dotx = $dotx")

end

# simulation A in TNNLS with F-16 model
function sim_TNNLS_A(NN::Int)#,amplitude, percent) # NN is the # of iterations
	# start = time()
    n, m = 3, 1
    A = [-1.01887 0.90506 -0.00215; 0.82225 -1.07741 -0.17555; 0 0 -1] # n x n
    B = reshape([0.0; 0.0; 1.0], (3, 1)) # n x m
    M = 1.0*Matrix(I, n, n) # n x n
    R = 0.1*Matrix(I, m, m) # m x m
	Pt = 0.5*Matrix(I, n, n)

	# check controllability & observability
    Co = ctrb(A, B)
    unco = size(A, 1) - rank(Co)
    # println("unco = $unco")
    Ob = obsv(sqrt(M), A)
    unob = size(sqrt(M), 1) - rank(Ob)
    # println("unob = $unob")

	if unco + unob > 0
		println("system uncontrollable and/or unobservable")
		return
	end

    # ODE parameters
    Tf, T, N = 45, 0.05, 900 # finite horizon

    αc, αa = 90, 2.5

	amplitude, percent = 0.05, 50

    Wc0 = [10.0*ones(Int((n + m)*(n + m + 1)/2 - m*m), 1); reshape(R, (m*m, 1))][:] # (n+m)(n+m+1)/2 x 1; 10 x 1
	Wa0 = 0.5*ones(n,) # n x m; 3 x 1
    u0 = 5.0*ones(m,) # m x 1; 1 x 1

    x0 = [1.0, 5.0, 1.0]
    xf = [2.0, 7.0, 3.0]
	#
	# Quu = 0
	# Qxu = 0

    p0 = x0'*M*x0

    t_save = [0,]
    u_save = [u0]
    # x_save = zeros(17, 1) # 3+10+3+1
	x_save = [[x0; Wc0; Wa0; p0],]
    # println("x_save = $x_save")

    uDelay  = interp2PWC(u_save, -1, 1) # return an interpolation function
	x1Delay = interp2PWC(getindex.(x_save, 1), -1, 1)
    x2Delay = interp2PWC(getindex.(x_save, 2), -1, 1)
    x3Delay = interp2PWC(getindex.(x_save, 3), -1, 1)

	# uvec = 0
	xdist = norm(x0 - xf)
	error = 0
    # solve ODEs
	for i = 1:NN
		start = time()
		# global xf, uf, Tf, T, A, B, M, R, Pt, per, amp, αa, αc, Quu, Qxu, Wcfinal, uDelay, x1Delay, x2Delay, x3Delay, x_save, t_save, uvec
		# uvec = 1
		t = ((i-1)*T, i*T)
		# t = (0, 1)
		p = [xf, Tf, T, A, B, M, R, Pt, percent, amplitude, αa, αc, uDelay, x1Delay, x2Delay, x3Delay, x_save, t_save, u_save, NN,start]
		# if i == NN println(1) end
        # sol = solve(ODEProblem(babyETC_A, x_save[end], t, p), DP5())
		sol = solve(ODEProblem(babyCT_A, x_save[end], t, p), DP5())
		# if i == NN println(2) end
		t_save = [t_save; sol.t[2:end]] # vcat(t_save, sol.t) save time and kick out the duplicated initial value
		x_save = [x_save; sol.u[2:end]] # vcat(x_save, sol.u) save state
		# if i == NN println(3) end
		# println("sol.u[$i] = $(sol.u)")
		# println("sol.t[$i] = $(sol.t)")
		u_save = p[end-3]
		# uvec = p[end]
		# println("uvec[$i] = $(p[end])")
		# if i == NN println(4) end
        uDelay  = interp2PWC(u_save, -1, i*T+.01) # interpolate control input
		x1Delay = interp2PWC(getindex.(x_save, 1), -1, i*T+.01)
		x2Delay = interp2PWC(getindex.(x_save, 2), -1, i*T+.01)
		x3Delay = interp2PWC(getindex.(x_save, 3), -1, i*T+.01)
		# println(i)
		# if euclidean(x_save[end][1:3], xf) < error*xdist break end
		if norm(x_save[end][1:3] - xf) < error*xdist break end
	end

	x_save_1 = getindex.(x_save, 1)
	# x_save_2 = getindex.(x_save, 2)
	# x_save_3 = getindex.(x_save, 3)

	# println("t_save = $t_save")
	# println("u_save = $u_save")
	# println("length_t_save = $(length(t_save))")
	# println("length_x_save = $(length(x_save))")
	# println("x1 = $(x_save_1)")
	# println("x2 = $(x_save_2[end])")
	# println("x3 = $(x_save_3[end])")
	println("x1 = $(x_save[end][1])")
	println("x2 = $(x_save[end][2])")
	println("x3 = $(x_save[end][3])")
	# println("$(x_save[2][1:3]-x_save[1][1:3])")
	plot(t_save, getindex.(x_save, 1), width = 2)
	plot!(t_save, getindex.(x_save, 2), width = 2)
	plot!(t_save, getindex.(x_save, 3), width = 2)
	# plot(t_save, u_save, width = 2)
	# println("Wa = $(x_save[end][14:16])")
	# println(x_save[end])
	# elapsed = time() - start
end

# ode fcn for simulation A in TNNLS
function babyCT_A(dotx, x, p, t)

	# global xf, uf, Tf, T, A, B, M, R, Pt, per, amp, αa, αc, Quu, Qxu, uDelay, x1Delay, x2Delay, x3Delay, x_save, t_save, uvec
	xf, Tf, T, A, B, M, R, Pt, percent, amplitude, αa, αc, uDelay, x1Delay, x2Delay, x3Delay, x_save, t_save, u_save,NN,start = p # uu useless here
	# println(t[1][1])
    n, m = 3, 1
	Wc = x[(n+1) : Int(n+(n+m)*(n+m+1)/2)] # 4-13
    Wa = x[Int(n+(n+m)*(n+m+1)/2+1) : Int(n+(n+m)*(n+m+1)/2+n)] # 14-16
    P = x[end] # 17

	# Update control
    ud = zeros(m,)
	ud[1] = Wa'*(x[1:n] - xf)
	# g(i,NN,2)
	# Delays
    uddelay = zeros(m,)
	uddelay[1] = uDelay(t - T)
    xdelay = zeros(n,)
	xdelay[1] = x1Delay(t - T)
    xdelay[2] = x2Delay(t - T)
    xdelay[3] = x3Delay(t - T)
	# g(i,NN,3)
	# Kronecker products
	U = vcat(x[1:n] - xf, ud) # augmented state
    UkU = vcat(U[1]^2, U[1]*U[2], U[1]*U[3], U[1]*ud,
               U[2]^2, U[2]*U[3], U[2]*ud,
               U[3]^2, U[3]*ud,
               ud[1]^2)
    UkUdelay = vcat(xdelay[1]^2, xdelay[1]*xdelay[2], xdelay[1]*xdelay[3], xdelay[1]*uddelay,
                    xdelay[2]^2, xdelay[2]*xdelay[3], xdelay[2]*uddelay,
                    xdelay[3]^2, xdelay[3]*uddelay,
                    uddelay[1]^2)
	# g(i,NN,4)
    Quu = Wc[end] # m x m
    Quu_inv = inv(Quu)
    Qxu = reshape(Wc[Int((n+m)*(n+m+1)/2-m^2-m*n+1):Int((n+m)*(n+m+1)/2-m^2)], (n, m)) # n x m
    Qux = Qxu' # m x n

    σ = UkU - UkUdelay
	σ_f	= UkU

	# Integral reinforcement dynamics
	dP = 0.5*(U[1:n]'*M*U[1:n] - xdelay'*M*xdelay + ud'*R*ud - uddelay'*R*uddelay)

	# approximation errors
    ec = P + Wc'*σ
	ecfinal = 0.5*U[1:n]'*Pt*U[1:n] - Wc'*σ_f
    ea = Wa'*U[1:n] + [Quu_inv*Qux[:,1]; Quu_inv*Qux[:,2]; Quu_inv*Qux[:,3]]'*U[1:n]

	# Critic dynamics
	dWc = -αc*((σ./(σ'*σ+1)^2)*ec + (σ_f./(σ_f'*σ_f+1)^2)*ecfinal)
	# g(i,NN,5)
	# Actor dynamics
    dWa = -αa*U[1:n]*ea'

	# Persistence excitation
    unew = zeros(m,)
	# if t <= (p[10]/100)*p.Tf
	if t <= (percent/100)*Tf
	    # unew[1] = ud[1] + amplitude*exp(-.0005*t)*(sin(0.1*t)^2 + cos(0.7*t)^2) + sin(1.5*t)^2*cos(0.1*t) + sin(pi*t) + cos(0.1*t)
	    # unew[1] = ud[1] + p[11]*exp(-.0005*t)*(sin(0.1*t)^2+cos(0.7*t)^2)+sin(1.5*t)^2*cos(0.1*t)+sin(pi*t)+cos(0.1*t)
	    unew[1] = (ud[1]+amplitude*exp(-0.009*t)*2*(sin(t)^2*cos(t)+sin(2*t)^2*cos(0.1*t)+sin(-1.2*t)^2*cos(0.5*t)+sin(t)^5+sin(1.12*t)^2+cos(2.4*t)*sin(2.4*t)^3))
	    # unew[2] = (ud[2]+amplitude*exp(-0.009*t)*2*(sin(t)^2*cos(t)+sin(2*t)^2*cos(0.1*t)+sin(-1.2*t)^2*cos(0.5*t)+sin(t)^5+sin(1.12*t)^2+cos(2.4*t)*sin(2.4*t)^3))
	else
	    unew = ud
	end
	# g(i,NN,6)
	u_save = [u_save; unew]
	# p[end-3] = unew
	# Augmented state
	# for j in 1:n dotx[j] = dx[j] end
	# for j in (n+1):convert(Int64, (n+(n+m)*(n+m+1)/2)) dotx[j] = dWc[j-n] end
	# for j in convert(Int64, (n+(n+m)*(n+m+1)/2+1)):convert(Int64, (n+(n+m)*(n+m+1)/2+n)) dotx[j] = dWa[j-convert(Int64, (n+(n+m)*(n+m+1)/2))] end
	# dotx[end] = dP
	dx = A*U[1:n] + B*unew
	if (time() - start) > 2
		# return dotx .= vcat(dx, dWc, dWa, dP)
		return
	end
	dotx .= vcat(dx, dWc, dWa, dP)
	# println("dotx = $dotx")

end

# simulation B in TNNLS, also for book chapter
function sim_TNNLS_B_CT_Local(x1, x2, S)#, goalVelocity) # x1 -> x2

	# System Dynamics for feedback
    n, m = 4, 2
	m_net = 10 # 10 - WAFR
	m_fuel = 30 # 30 - WAFR
	alpha = 0.05 # 0.05 - WAFR
	m_total = m_net + m_fuel # 40 - WAFR
	kx = 20 # 10 - CDC, 20 - WAFR
	ky = 20 # 10 - CDC, 20 - WAFR
	cx = 45 # 10 - CDC, 45 - WAFR
	cy = 45 # 10 - CDC, 45 - WAFR

    A = [0 0 1 0; 0 0 0 1; -kx/m_total 0 -cx/m_total 0; 0 -ky/m_total 0 -cx/m_total] # n x n
    B = [0 0; 0 0; 1/m_total 0; 0 1/m_total] # n x m
    M = 10.0*Matrix(I, n, n) # n x n
    R = 2.0*Matrix(I, m, m) # m x m
	Pt = 0.5*Matrix(I, n, n)

	# check controllability & observability
    Co = ctrb(A, B)
    unco = size(A, 1) - rank(Co)
    # println("unco = $unco")
    Ob = obsv(sqrt(M), A)
    unob = size(sqrt(M), 1) - rank(Ob)
    # println("unob = $unob")

	if unco + unob > 0  error("system uncontrollable and/or unobservable")  end

    # ODE parametersYo
    Tf, T, N = 10, 0.05, 200 # finite horizon
    αc, αa = 90, 1.2
	amplitude, percent = 0.1, 50

	# start and goal
	xf = [x2; 0.0;0.0]
	x0 = [x1; S.lastVelocity]#0.0;0.0]

	# goal orientation
	theta = atan((x2-x1)[2]/(x2-x1)[1])

	# initialization
	Wc0 = [10.0*ones(4); 5cos(theta+1); 5sin(theta+1); 10.0*ones(3); 5cos(theta+1); 5sin(theta+1); 10.0*ones(2); 5cos(theta+1); 5sin(theta+1); 10.0; 5cos(theta+1); 5sin(theta+1); 1; 0; 1];
	Wa10 = 5*cos(theta+1)*ones(n,) # n x m; 4 x 2
	Wa20 = 5*sin(theta+1)*ones(n,)
    # Wc0 = [10.0*ones(Int((n+m)*(n+m+1)/2 - m*(m+1)/2), 1); 1; 0; 1][:] # (n+m)(n+m+1)/2 x 1; 10 x 1
	# Wa10 = 0.5*ones(n,) # n x m; 4 x 2
	# Wa20 = 0.5*ones(n,)

	u0 = [Wa10'*(x0 - xf); Wa20'*(x0 - xf)] # 5.0*ones(m,) # m x 1; 2 x 1

    t_save = [0,]
    u_save = [u0,] # u at every T
	x_save = [[x0; Wc0; Wa10; Wa20; 0],]
	uvec = [u0,] # all u
    # println("u_save = $u_save")

    u1Delay = interp2PWC(getindex.(u_save, 1), -1, 1) # return an interpolation function
	u2Delay = interp2PWC(getindex.(u_save, 2), -1, 1)
	x1Delay = interp2PWC(getindex.(x_save, 1), -1, 1)
    x2Delay = interp2PWC(getindex.(x_save, 2), -1, 1)
    x3Delay = interp2PWC(getindex.(x_save, 3), -1, 1)
	x4Delay = interp2PWC(getindex.(x_save, 4), -1, 1)
	pDelay = interp2PWC(getindex.(x_save, 34), -1, 1)

	# uvec1, uvec2 = 0, 0
	# trigCondvec = zeros(0)
	# normESqvec = zeros(0)
	xdist = norm(x0[1:2] - xf[1:2])
	error = 0.15
	# localKd = 0.0
	maxIter = 10000
	poseAndKd = Array{Tuple{Array{Float64,2},Float64}}(undef,0)

    # solve ODEs
	for i = 1:maxIter
		# global xf, uf, Tf, T, A, B, M, R, Pt, per, amp, αa, αc, Quu, Qxu, Wcfinal, uDelay, x1Delay, x2Delay, x3Delay, x_save, t_save, uvec
		t = ((i-1)*T, i*T)
		p = [xf, T, Tf, A, B, M, R, Pt, percent, amplitude, αa, αc, u1Delay, u2Delay, x1Delay, x2Delay, x3Delay, x4Delay, pDelay, uvec] #, normESqvec, trigCondvec]#, x_save, t_save, uvec1, uvec2]#, i,maxIter,start]
		# if i == NN println(1) end
        sol = solve(ODEProblem(babyCT_B, x_save[end], t, p), DP5())
		# if i == NN println(2) end
		t_save = [t_save; sol.t[2:end]] # vcat(t_save, sol.t) save time
		x_save = [x_save; sol.u[2:end]] # vcat(x_save, sol.u) save state
		u_save = [u_save; [uvec[end]]]
		# if i == NN println(3) end
		# println("sol.u[$i] = $(sol.u)")
		# println("sol.t[$i] = $(sol.t)")
		# uvec1 = p[end-1]
		# uvec2 = p[end]
		# println("uvec[$i] = $(p[end])")
		# if i == NN println(4) end
        u1Delay = interp2PWC(getindex.(u_save, 1), -1, i*T+.01) # interpolate control input
		u2Delay = interp2PWC(getindex.(u_save, 2), -1, i*T+.01)
		x1Delay = interp2PWC(getindex.(x_save, 1), -1, i*T+.01)
		x2Delay = interp2PWC(getindex.(x_save, 2), -1, i*T+.01)
		x3Delay = interp2PWC(getindex.(x_save, 3), -1, i*T+.01)
		x4Delay = interp2PWC(getindex.(x_save, 4), -1, i*T+.01)
		pDelay = interp2PWC(getindex.(x_save, 34), -1, i*T+.01)
		# println(i)

		localKd = distPoint2Line(x_save[end][1:2], x1, x2)
		poseAndKd = [poseAndKd; ([x_save[end][1] x_save[end][2]], copy(localKd))]
		# if distPoint2Line(x_save[end][1:2], x1, x2) > kino_dist
		# 	kino_dist = distPoint2Line(x_save[end][1:2], x1, x2)
		# end

		# iter += 1

		# if euclidean(x_save[end][1:3], xf) < error*xdist break end
		if norm(x_save[end][1:2] - xf[1:2]) < error*xdist break end
	end
	S.lastVelocity = x_save[end][3:4]
	(poseAndKd, [0; 0], [0; 0])
	# (poseAndKd, normESqvec, trigCondvec)#,t_save)#eventTimes)
	# plot(getindex.(getindex.(poseAndKd, 1), 1), getindex.(getindex.(poseAndKd, 1), 2))
	# x_save_1 = getindex.(x_save, 1)
	# x_save_2 = getindex.(x_save, 2)
	# x_save_3 = getindex.(x_save, 3)
	# println("t_save = $t_save")
	# println("x1 = $(x_save_1[end])")
	# println("x2 = $(x_save_2[end])")
	# println("x3 = $(x_save_3[end])")
	# println("x1 = $(x_save[end][1])")
	# println("x2 = $(x_save[end][2])")
	# println("x3 = $(x_save[end][3])")
	# println("x4 = $(x_save[end][4])")
	# (path = [getindex.(x_save, 1) getindex.(x_save, 2)], kino_dist = kino_dist, lastVelocity = x_save[end][3:4])#, time = t_save)
	# [getindex.(x_save, 1) getindex.(x_save, 2) getindex.(x_save, 3) getindex.(x_save, 4)]
	# plot(getindex.(x_save, 1), getindex.(x_save, 2), width = 2, aspectratio = 1)
	# plot(t_save, getindex.(x_save, 1), width = 2)
	# plot!(t_save, getindex.(x_save, 2), width = 2)
	# plot!(t_save, getindex.(x_save, 3), width = 2)
	# plot!(t_save, getindex.(x_save, 4), width = 2)
	# elapsed = time() - start
end

function sim_TNNLS_B_CT_Local_Max(x1, x2, lastVel)#, goalVelocity) # x1 -> x2

	# System Dynamics for feedback
    n, m = 4, 2
	m_net = 10 # 10 - WAFR
	m_fuel = 30 # 30 - WAFR
	alpha = 0.05 # 0.05 - WAFR
	m_total = m_net + m_fuel # 40 - WAFR
	kx = 20 # 10 - CDC, 20 - WAFR
	ky = 20 # 10 - CDC, 20 - WAFR
	cx = 45 # 10 - CDC, 45 - WAFR
	cy = 45 # 10 - CDC, 45 - WAFR

    A = [0 0 1 0; 0 0 0 1; -kx/m_total 0 -cx/m_total 0; 0 -ky/m_total 0 -cx/m_total] # n x n
    B = [0 0; 0 0; 1/m_total 0; 0 1/m_total] # n x m
    M = 10.0*Matrix(I, n, n) # n x n
    R = 2.0*Matrix(I, m, m) # m x m
	Pt = 0.5*Matrix(I, n, n)

	# check controllability & observability
    Co = ctrb(A, B)
    unco = size(A, 1) - rank(Co)
    # println("unco = $unco")
    Ob = obsv(sqrt(M), A)
    unob = size(sqrt(M), 1) - rank(Ob)
    # println("unob = $unob")

	if unco + unob > 0  error("system uncontrollable and/or unobservable")  end

    # ODE parameters
    Tf, T, N = 10, 0.05, 200 # finite horizon
    αc, αa = 90, 1.2
	amplitude, percent = 0.1, 50

	# start and goal
	xf = [x2; 0.0;0.0]
	x0 = [x1; lastVel]#0.0;0.0]

	# goal orientation
	theta = atan((x2-x1)[2]/(x2-x1)[1])

	# initialization
	Wc0 = [10.0*ones(4); 5cos(theta+1); 5sin(theta+1); 10.0*ones(3); 5cos(theta+1); 5sin(theta+1); 10.0*ones(2); 5cos(theta+1); 5sin(theta+1); 10.0; 5cos(theta+1); 5sin(theta+1); 1; 0; 1];
	Wa10 = 5*cos(theta+1)*ones(n,) # n x m; 4 x 2
	Wa20 = 5*sin(theta+1)*ones(n,)
    # Wc0 = [10.0*ones(Int((n+m)*(n+m+1)/2 - m*(m+1)/2), 1); 1; 0; 1][:] # (n+m)(n+m+1)/2 x 1; 10 x 1
	# Wa10 = 0.5*ones(n,) # n x m; 4 x 2
	# Wa20 = 0.5*ones(n,)

	u0 = [Wa10'*(x0 - xf); Wa20'*(x0 - xf)] # 5.0*ones(m,) # m x 1; 2 x 1

    t_save = [0,]
    u_save = [u0,] # u at every T
	x_save = [[x0; Wc0; Wa10; Wa20; 0],]
	uvec = [u0,] # all u
    # println("u_save = $u_save")

    u1Delay = interp2PWC(getindex.(u_save, 1), -1, 1) # return an interpolation function
	u2Delay = interp2PWC(getindex.(u_save, 2), -1, 1)
	x1Delay = interp2PWC(getindex.(x_save, 1), -1, 1)
    x2Delay = interp2PWC(getindex.(x_save, 2), -1, 1)
    x3Delay = interp2PWC(getindex.(x_save, 3), -1, 1)
	x4Delay = interp2PWC(getindex.(x_save, 4), -1, 1)
	pDelay = interp2PWC(getindex.(x_save, 34), -1, 1)

	# uvec1, uvec2 = 0, 0
	# trigCondvec = zeros(0)
	# normESqvec = zeros(0)
	xdist = norm(x0[1:2] - xf[1:2])
	error = 0.15
	# localKd = 0.0
	maxIter = 10000
	poseAndKd = Array{Tuple{Array{Float64,2},Float64}}(undef,0)

    # solve ODEs
	for i = 1:maxIter
		# global xf, uf, Tf, T, A, B, M, R, Pt, per, amp, αa, αc, Quu, Qxu, Wcfinal, uDelay, x1Delay, x2Delay, x3Delay, x_save, t_save, uvec
		t = ((i-1)*T, i*T)
		p = [xf, T, Tf, A, B, M, R, Pt, percent, amplitude, αa, αc, u1Delay, u2Delay, x1Delay, x2Delay, x3Delay, x4Delay, pDelay, uvec] #, normESqvec, trigCondvec]#, x_save, t_save, uvec1, uvec2]#, i,maxIter,start]
		# if i == NN println(1) end
        sol = solve(ODEProblem(babyCT_B, x_save[end], t, p), DP5())
		# if i == NN println(2) end
		t_save = [t_save; sol.t[2:end]] # vcat(t_save, sol.t) save time
		x_save = [x_save; sol.u[2:end]] # vcat(x_save, sol.u) save state
		u_save = [u_save; [uvec[end]]]
		# if i == NN println(3) end
		# println("sol.u[$i] = $(sol.u)")
		# println("sol.t[$i] = $(sol.t)")
		# uvec1 = p[end-1]
		# uvec2 = p[end]
		# println("uvec[$i] = $(p[end])")
		# if i == NN println(4) end
        u1Delay = interp2PWC(getindex.(u_save, 1), -1, i*T+.01) # interpolate control input
		u2Delay = interp2PWC(getindex.(u_save, 2), -1, i*T+.01)
		x1Delay = interp2PWC(getindex.(x_save, 1), -1, i*T+.01)
		x2Delay = interp2PWC(getindex.(x_save, 2), -1, i*T+.01)
		x3Delay = interp2PWC(getindex.(x_save, 3), -1, i*T+.01)
		x4Delay = interp2PWC(getindex.(x_save, 4), -1, i*T+.01)
		pDelay = interp2PWC(getindex.(x_save, 34), -1, i*T+.01)
		# println(i)

		localKd = distPoint2Line(x_save[end][1:2], x1, x2)
		poseAndKd = [poseAndKd; ([x_save[end][1] x_save[end][2]], copy(localKd))]
		# if distPoint2Line(x_save[end][1:2], x1, x2) > kino_dist
		# 	kino_dist = distPoint2Line(x_save[end][1:2], x1, x2)
		# end

		# iter += 1

		# if euclidean(x_save[end][1:3], xf) < error*xdist break end
		if norm(x_save[end][1:2] - xf[1:2]) < error*xdist break end
	end
	lastVel = x_save[end][3:4]
	maxKd = maximum(getindex.(poseAndKd,2))
	(maxKd, lastVel)
	# (poseAndKd, normESqvec, trigCondvec)#,t_save)#eventTimes)
	# plot(getindex.(getindex.(poseAndKd, 1), 1), getindex.(getindex.(poseAndKd, 1), 2))
	# x_save_1 = getindex.(x_save, 1)
	# x_save_2 = getindex.(x_save, 2)
	# x_save_3 = getindex.(x_save, 3)
	# println("t_save = $t_save")
	# println("x1 = $(x_save_1[end])")
	# println("x2 = $(x_save_2[end])")
	# println("x3 = $(x_save_3[end])")
	# println("x1 = $(x_save[end][1])")
	# println("x2 = $(x_save[end][2])")
	# println("x3 = $(x_save[end][3])")
	# println("x4 = $(x_save[end][4])")
	# (path = [getindex.(x_save, 1) getindex.(x_save, 2)], kino_dist = kino_dist, lastVelocity = x_save[end][3:4])#, time = t_save)
	# [getindex.(x_save, 1) getindex.(x_save, 2) getindex.(x_save, 3) getindex.(x_save, 4)]
	# plot(getindex.(x_save, 1), getindex.(x_save, 2), width = 2, aspectratio = 1)
	# plot(t_save, getindex.(x_save, 1), width = 2)
	# plot!(t_save, getindex.(x_save, 2), width = 2)
	# plot!(t_save, getindex.(x_save, 3), width = 2)
	# plot!(t_save, getindex.(x_save, 4), width = 2)
	# elapsed = time() - start
end

# ode fcn for simulation B in TNNLS
function babyCT_B(dotx, x, p, t)

	# global xf, uf, Tf, T, A, B, M, R, Pt, per, amp, αa, αc, Quu, Qxu, Wcfinal, uDelay, x1Delay, x2Delay, x3Delay, x_save, t_save, uvec
	xf, T, Tf, A, B, M, R, Pt, percent, amplitude, αa, αc, u1Delay, u2Delay, x1Delay, x2Delay, x3Delay, x4Delay, pDelay, uvec = p#,i,maxIter,start = p # uu, vv useless here
	# println(t[1][1])
    n, m = 4, 2
	Wc = x[(n+1) : Int(n+(n+m)*(n+m+1)/2)] # 4-13
    Wa1 = x[26:29]#x[convert(Int64, (n+(n+m)*(n+m+1)/2+1)) : convert(Int64, (n+(n+m)*(n+m+1)/2+n))] # 14-16
	Wa2 = x[30:33]
	# P = x[end] # 17
	P = x[end]

	# Update control
    ud = zeros(m,)
	ud[1] = Wa1'*(x[1:n] - xf)
	ud[2] = Wa2'*(x[1:n] - xf)

	# Delays
    uddelay = zeros(m,)
	uddelay[1] = u1Delay(t - T)
	uddelay[2] = u2Delay(t - T)
    xdelay = zeros(n,)
	xdelay[1] = x1Delay(t - T)
    xdelay[2] = x2Delay(t - T)
    xdelay[3] = x3Delay(t - T)
	xdelay[4] = x4Delay(t - T)
	pdelay = pDelay(t - T)

	# Kronecker products
	U = vcat(x[1:n] - xf, ud[1], ud[2]) # augmented state
    UkU = vcat(U[1]^2, U[1]*U[2], U[1]*U[3], U[1]*U[4], U[1]*ud[1], U[1]*ud[2],
               U[2]^2, U[2]*U[3], U[2]*U[4], U[2]*ud[1], U[2]*ud[2],
               U[3]^2, U[3]*U[4], U[3]*ud[1], U[3]*ud[2],
			   U[4]^2, U[4]*ud[1], U[4]*ud[2],
               ud[1]^2, ud[1]*ud[2],
			   ud[2]^2)
    UkUdelay = vcat(xdelay[1]^2, xdelay[1]*xdelay[2], xdelay[1]*xdelay[3], xdelay[1]*xdelay[4], xdelay[1]*uddelay[1], xdelay[1]*uddelay[2],
                    xdelay[2]^2, xdelay[2]*xdelay[3], xdelay[2]*xdelay[4], xdelay[2]*uddelay[1], xdelay[2]*uddelay[2],
                    xdelay[3]^2, xdelay[3]*xdelay[4], xdelay[3]*uddelay[1], xdelay[3]*uddelay[2],
					xdelay[4]^2, xdelay[4]*uddelay[1], xdelay[4]*uddelay[2],
                    uddelay[1]^2, uddelay[1]*uddelay[2],
					uddelay[2]^2)
	# g(i,maxIter,4)
    Quu = [Wc[end-2] Wc[end-1]; Wc[end-1] Wc[end]] # m x m
    Quu_inv = inv(Quu)
    Qux = [Wc[5] Wc[10] Wc[14] Wc[17]; Wc[6] Wc[11] Wc[15] Wc[18]] # n x m
    Qxu = Qux' # m x n

	# Integral reinforcement dynamics
	dP = 0.5*(U[1:n]'*M*U[1:n] + ud'*R*ud)

	# approximation errors
	σ = UkU - UkUdelay
	σ_f	= UkU
    ec = P - pdelay + Wc'*σ
	ecfinal = 0.5*U[1:n]'*Pt*U[1:n] - Wc'*σ_f
    ea1 = Wa1'*U[1:n] + [Quu_inv[1,:]'*Qux[:,1]; Quu_inv[1,:]'*Qux[:,2]; Quu_inv[1,:]'*Qux[:,3]; Quu_inv[1,:]'*Qux[:,4]]'*U[1:n]
	ea2 = Wa2'*U[1:n] + [Quu_inv[2,:]'*Qux[:,1]; Quu_inv[2,:]'*Qux[:,2]; Quu_inv[2,:]'*Qux[:,3]; Quu_inv[2,:]'*Qux[:,4]]'*U[1:n]

	# Critic dynamics
	dWc = -αc*((σ./(σ'*σ+1).^2)*ec + (σ_f./(σ_f'*σ_f+1).^2)*ecfinal)
	# g(i,maxIter,5)
	# Actor dynamics
    dWa1 = -αa*U[1:n]*ea1'
	dWa2 = -αa*U[1:n]*ea2'

	# Persistence excitation
    unew = zeros(m,)
	# if t <= (p[10]/100)*p.Tf
	if t <= (percent/100)*Tf
	    # unew[1] = ud[1] + amplitude*exp(-.0005*t)*(sin(0.1*t)^2 + cos(0.7*t)^2) + sin(1.5*t)^2*cos(0.1*t) + sin(pi*t) + cos(0.1*t)
	    # unew[1] = ud[1] + p[11]*exp(-.0005*t)*(sin(0.1*t)^2+cos(0.7*t)^2)+sin(1.5*t)^2*cos(0.1*t)+sin(pi*t)+cos(0.1*t)
	    unew[1] = (ud[1]+amplitude*exp(-0.009*t)*2*(sin(t)^2*cos(t)+sin(2*t)^2*cos(0.1*t)+sin(-1.2*t)^2*cos(0.5*t)+sin(t)^5+sin(1.12*t)^2+cos(2.4*t)*sin(2.4*t)^3))
	    unew[2] = (ud[2]+amplitude*exp(-0.009*t)*2*(sin(t)^2*cos(t)+sin(2*t)^2*cos(0.1*t)+sin(-1.2*t)^2*cos(0.5*t)+sin(t)^5+sin(1.12*t)^2+cos(2.4*t)*sin(2.4*t)^3))
	else
	    unew .= ud
	end
	dx = A*U[1:n] + B*unew
	# g(i,maxIter,6)
	# p[end] = unew
	# p[end-1:end] .= unew
	uvec = [uvec; [unew]]
	# normESqvec = [normESqvec; 0]
	# trigCondvec = [trigCondvec; 0]
	# if (time() - start) > 2
	# 	# return dotx .= vcat(dx, dWc, dWa, dP)
	# 	return
	# end
	dotx .= vcat(dx, dWc, dWa1, dWa2, dP)

end

# ode fcn for simulation B in TNNLS with ETC
function babyETC_B(dotx, x, p, t)

	global eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa10, Wa20, sizePlant, β, L, L1, Wavec, eventTimes

	# global xf, uf, Tf, T, A, B, M, R, Pt, per, amp, αa, αc, Quu, Qxu, Wcfinal, uDelay, x1Delay, x2Delay, x3Delay, x_save, t_save, uvec
	xf, uf, Tf, T, A, B, M, R, Pt, percent, amplitude, αa, αc, Quu, Qxu, Wcfinal, u1Delay, u2Delay, x1Delay, x2Delay, x3Delay, x4Delay, x_save, t_save, uu, vv=p#,eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa10, Wa20, sizePlant, β, L, L1, Wavec, eventTimes=p#,i,maxIter,start = p # uu, vv useless here
	# println(t[1][1])
    n, m = 4, 2
	Wc = x[(n+1) : Int(n+(n+m)*(n+m+1)/2)] # 4-13
    # Wa1 = x[26:29]#x[convert(Int64, (n+(n+m)*(n+m+1)/2+1)) : convert(Int64, (n+(n+m)*(n+m+1)/2+n))] # 14-16
	# Wa2 = x[30:33]
	# P = x[end] # 17
	P = x[end]

	# Augmented final state
	UkUfinal = vcat(xf[1]^2, xf[1]*xf[2], xf[1]*xf[3], xf[1]*xf[4], xf[1]*uf[1], xf[1]*uf[2],
					xf[2]^2, xf[2]*xf[3], xf[2]*xf[4], xf[2]*uf[1], xf[2]*uf[2],
					xf[3]^2, xf[3]*xf[4], xf[3]*uf[1], xf[3]*uf[2],
					xf[4]^2, xf[4]*uf[1], xf[4]*uf[2],
					uf[1]^2, uf[1]*uf[2],
					uf[2]^2)
    Pfinal = 0.5*xf'*Pt*xf

	# Update control
    ud = zeros(m,)
	# ud[1] = Wa10'*(x[1:n] - xf)
	# ud[2] = Wa20'*(x[1:n] - xf)
	ud[1] = Wa10'*Xhat
	ud[2] = Wa20'*Xhat
	# Delays
    uddelay = zeros(m,)
	uddelay[1] = u1Delay(t - T)
	uddelay[2] = u2Delay(t - T)
    xdelay = zeros(n,)
	xdelay[1] = x1Delay(t - T)
    xdelay[2] = x2Delay(t - T)
    xdelay[3] = x3Delay(t - T)
	xdelay[4] = x4Delay(t - T)
	# Kronecker products
	U = vcat(x[1:n] - xf, ud[1], ud[2]) # augmented state
    UkU = vcat(U[1]^2, U[1]*U[2], U[1]*U[3], U[1]*U[4], U[1]*ud[1], U[1]*ud[2],
               U[2]^2, U[2]*U[3], U[2]*U[4], U[2]*ud[1], U[2]*ud[2],
               U[3]^2, U[3]*U[4], U[3]*ud[1], U[3]*ud[2],
			   U[4]^2, U[4]*ud[1], U[4]*ud[2],
               ud[1]^2, ud[1]*ud[2],
			   ud[2]^2)
    UkUdelay = vcat(xdelay[1]^2, xdelay[1]*xdelay[2], xdelay[1]*xdelay[3], xdelay[1]*xdelay[4], xdelay[1]*uddelay[1], xdelay[1]*uddelay[2],
                    xdelay[2]^2, xdelay[2]*xdelay[3], xdelay[2]*xdelay[4], xdelay[2]*uddelay[1], xdelay[2]*uddelay[2],
                    xdelay[3]^2, xdelay[3]*xdelay[4], xdelay[3]*uddelay[1], xdelay[3]*uddelay[2],
					xdelay[4]^2, xdelay[4]*uddelay[1], xdelay[4]*uddelay[2],
                    uddelay[1]^2, uddelay[1]*uddelay[2],
					uddelay[2]^2)
    Quu = [Wc[end-3] Wc[end-2]; Wc[end-1] Wc[end]] # m x m
    Quu_inv = inv(Quu)
    Qxu = reshape(Wc[Int((n+m)*(n+m+1)/2-m^2-m*n+1):Int((n+m)*(n+m+1)/2-m^2)], (n, m)) # n x m
    Qux = Qxu' # m x n

    σ = UkU - UkUdelay

	# Triggering condition
	normXSq = abs(U[1])^2 + abs(U[2])^2 + U[3]^2 + U[4]^2
	trigCond = ((1 - β^2)*eigMinM/(4*(L^2+L1^2)*eigMaxR))*normXSq + (eigMinR/(4*(L^2+L1^2)*eigMaxR))*(abs(ud[1])^2+abs(ud[2])^2)

	# Integral reinforcement dynamics
	dP = 0.5*((x[1:n] - xf)'*M*(x[1:n] - xf) - xdelay'*M*xdelay + ud'*R*ud - uddelay'*R*uddelay)

	# approximation errors
    ec = P + Wc'*UkU - Wc'*UkUdelay
	ecfinal = Pfinal - Wcfinal'*UkUfinal
    # ea1 = Wa10'*U[1:n] + [Quu_inv[1,:]'*Qux[:,1]; Quu_inv[1,:]'*Qux[:,2]; Quu_inv[1,:]'*Qux[:,3]; Quu_inv[1,:]'*Qux[:,4]]'*U[1:n]
	# ea2 = Wa20'*U[1:n] + [Quu_inv[2,:]'*Qux[:,1]; Quu_inv[2,:]'*Qux[:,2]; Quu_inv[2,:]'*Qux[:,3]; Quu_inv[2,:]'*Qux[:,4]]'*U[1:n]

	# gap
	e = Xhat - U[1:n]
	normESq = norm(e)^2

	# Critic dynamics
	dWc = -αc*((σ./(σ'*σ+1).^2)*ec + (UkUfinal./(UkUfinal'*UkUfinal+1).^2)*ecfinal)
	# # Actor dynamics
    # dWa1 = -αa*U[1:n]*ea1'
	# dWa2 = -αa*U[1:n]*ea2'

	# Triggering updates
	if normESq >= 0.95*trigCond # flows
	    Xhat = U[1:n]
		ea1 = Wa10'*Xhat + [Quu_inv[1,:]'*Qux[:,1]; Quu_inv[1,:]'*Qux[:,2]; Quu_inv[1,:]'*Qux[:,3]; Quu_inv[1,:]'*Qux[:,4]]'*Xhat
		ea2 = Wa20'*Xhat + [Quu_inv[2,:]'*Qux[:,1]; Quu_inv[2,:]'*Qux[:,2]; Quu_inv[2,:]'*Qux[:,3]; Quu_inv[2,:]'*Qux[:,4]]'*Xhat
		Wa10 = Wa10 - αa*Xhat*ea1'
		Wa20 = Wa20 - αa*Xhat*ea2'
		# Wavec = [Wavec; Wa10'; Wa20']
		eventTimes = [eventTimes; t]
	end

	# Persistence excitation
    unew = zeros(m,)
	# if t <= (p[10]/100)*p.Tf
	if t <= (percent/100)*Tf
	    # unew[1] = ud[1] + amplitude*exp(-.0005*t)*(sin(0.1*t)^2 + cos(0.7*t)^2) + sin(1.5*t)^2*cos(0.1*t) + sin(pi*t) + cos(0.1*t)
	    # unew[1] = ud[1] + p[11]*exp(-.0005*t)*(sin(0.1*t)^2+cos(0.7*t)^2)+sin(1.5*t)^2*cos(0.1*t)+sin(pi*t)+cos(0.1*t)
	    unew[1] = (ud[1]+amplitude*exp(-0.009*t)*2*(sin(t)^2*cos(t)+sin(2*t)^2*cos(0.1*t)+sin(-1.2*t)^2*cos(0.5*t)+sin(t)^5+sin(1.12*t)^2+cos(2.4*t)*sin(2.4*t)^3))
	    unew[2] = (ud[2]+amplitude*exp(-0.009*t)*2*(sin(t)^2*cos(t)+sin(2*t)^2*cos(0.1*t)+sin(-1.2*t)^2*cos(0.5*t)+sin(t)^5+sin(1.12*t)^2+cos(2.4*t)*sin(2.4*t)^3))
	else
	    unew .= ud
	end
	dx = A*U[1:n] + B*unew
	# p[end] = unew
	# p[end-1:end] .= unew
	uu, vv = unew[1], unew[2]
	normESqvec = [normESqvec; normESq]
	# p[31] = [p[31]; normESq]
	trigCondvec = [trigCondvec; trigCond]
	# p[32] = [p[32]; trigCond]
	# println("$trigCondvec")

	dotx .= vcat(dx, dWc, dP)
end

# simulation B in TNNLS with no obstacle
function sim_TNNLS_B()
	A1 = [95
	95.1624763874382
	96.3500590637060
	95
	95.7491235875551
	95.8949548141154
	95.2662347545682
	94.5336038239492
	96.5842499393074
	94.6270087246292
	92.5781708327261
	87.8636823940887
	78.6153725115791
	75.8455106021162
	74.8657498993810
	73.2009390282228
	66.7554039660350
	61.6369123742990
	62.3376790194971
	66.1866185062487
	73.2496103174612
	74.9672816485974
	71.2384314193777
	62.9965487463500
	54.1592137828044
	45.0869966971835
	36.1550899086504
	30.3190139376311
	28.7757550839645
	31
	32.2432019856655
	37.2050843021910
	37.2629322061965
	31.2160663936524
	27.4654314433054
	26.5220950790295
	21.5694145274602
	19.8894748894053
	10.1877025150662
	6.76107477784817
	6.01791486042274
	4.67204153198610
	5.03990445537371
	4.76011337360909
	5.24515712193008
	5.29209680765824
	5.43568921537421
	5.39223227027637
	5]

	A2 = [5
	5.37253305758653
	14.5876961244558
	24
	32.7951116418459
	41.8142124861278
	51.1061226157655
	60.4983947935726
	68.6640489958543
	77.3180400565964
	85.4962489657631
	94.2604590177179
	91.3637628118781
	83.6540728490662
	74.1684988744362
	65.6149538046676
	59.4410996417406
	51.6573069198732
	41.9201699744068
	33.9985644913744
	28.6812739087495
	19.8223476444452
	10.8783283619548
	10.1174440439029
	8.37316040198836
	7.76038235196808
	10.0356183221714
	17.3971443366884
	26.4901805407076
	34
	42.0485057834142
	46.6395091707108
	56.3911233740648
	62.2766902498206
	69.4960648449329
	78.8283994953343
	86.9052055346931
	94.1996002583170
	93.4441679039921
	86.0653216398322
	76.5276437799853
	68.4130596316205
	58.5079038864695
	48.6025072470124
	38.6308531383844
	28.8655587027563
	18.9411561962995
	8.96116135138184
	5]
	lastVelocity = [0.; 0.]
	waypoint = zeros(0,2)
	numEdge = size(A1, 1) - 1
	for j = 1:numEdge
		x0 = [A1[j], A2[j]]
		xf = [A1[j+1], A2[j+1]]
		# if j == numEdge
		# 	goalVelocity = [A1[j+1] - A1[j]; A2[j+1] - A2[j]]
		# else
		# 	goalVelocity = [A1[j+2] - A1[j]; A2[j+2] - A2[j]]/2
		# end
		# goalVelocity = [A1[j+1] - A1[j]; A2[j+1] - A2[j]]
		localPath = sim_TNNLS_B_Single(x0, xf, lastVelocity)#, goalVelocity) sim_TNNLS_B_ETC_robust_relaxedPE_Local
		A1[j+1] = localPath[1][end, 1]
		A2[j+1] = localPath[1][end, 2]
		lastVelocity = localPath[3]
		waypoint = [waypoint; localPath[1]]
	end
	plot(waypoint[:, 1], waypoint[:, 2])
end

# interp2PWC is a function that approximates the input vector data to a
# piecewise continuous function given the initial and final time
function interp2PWC(y, xi, xf)
    row = size(y, 1)
	if row == 1
	    xdata = range(-1.0, stop = xf, length = row + 1) # linspace in MATLAB; collect also works
		itp = interpolate([y[1];y[end]], BSpline(Cubic(Interpolations.Line(OnGrid()))))
	else
		xdata = range(xi, stop = xf, length = row)
		itp = interpolate(y, BSpline(Cubic(Interpolations.Line(OnGrid()))))
		# xdata = range(xi, stop = xf, length = 2)
		# itp = interpolate([y[1];y[end]], BSpline(Cubic(Line(OnGrid()))))
	end
	Interpolations.scale(itp, xdata)
end

# for debugging
function g(i, NN, j)
	if i == NN
		println(j)
	end
end

# get kinodynamic distance
function distPoint2Line(p, v1, v2)
	a = [v1 - v2; 0]
	b = [p - v2; 0]
	d = norm(cross(a,b))/norm(a)
end

function readPath(filename)	readdlm(filename, ',', Float64) end

function testChangeTarget()
	x1 = [0.;0.]
	x2 = [2.;0.]
	x3 = [2.;2.]
	path = nothing
	path = sim_TNNLS_B_Single(x1, x2, 50)[1]
	xnew = path[end,:][:]
	path = [path; sim_TNNLS_B_Single(xnew, x3, 10000)[1]]
	plot(path[:,1], path[:,2], aspectratio = 1)
end

# change localTemp for different control scheme
function moveRobot_Q(S, Q, KD, slice_time, root, hyberBallRad, R, localPoseAndKd, localNormEsq, localTrigCond,save_elapsed_time)#, NormEsqvec, TrigCondvec)
	# global NormEsqvec, TrigCondvec
	if R.moving
	  R.robotPose = R.nextRobotPose
	  R.robotMovePath[R.numRobotMovePoints+1:R.numRobotMovePoints+R.numLocalMovePoints,:] = R.robotLocalPath[1:R.numLocalMovePoints,:]
	  R.numRobotMovePoints += R.numLocalMovePoints

	  # NormEsqVec[R.numEsqTrigPoints+1:R.numEsqTrigPoints+R.numLocalEsqTrigPoints] = R.localEsq[:]
	  # TrigCondVec[R.numEsqTrigPoints+1:R.numEsqTrigPoints+R.numLocalEsqTrigPoints] = R.localTrig[:]
	  # R.numEsqTrigPoints += R.numLocalEsqTrigPoints
	  # S.numEsqTrigCovered = S.numEsqTrigCoveredNext

	  println("new robot pose: [$(R.robotPose[1]) $(R.robotPose[2])]")# --- new kino dist: $(getindex(localPoseAndKd[S.numCoveredLocal],2))")
	else
	  # movement has just started, so remember that the robot is now moving
	  R.moving = true
	  # localPoseAndKd = [localPoseAndKd; ([1. 1.], 0.0)] # the goal node
	  if !S.moveGoal.rrtParentUsed
	    # no parent has been found for the node at the robot's position
	    R.currentMoveInvalid = true
	  end
	end

	if R.currentMoveInvalid
	  # save covered Esq & Trig
  	  S.NormEsqvec = [S.NormEsqvec; localNormEsq[1:Int(ceil(S.numEsqTrigLocal*S.numCoveredLocal/S.numLocal))]]
	  S.TrigCondvec = [S.TrigCondvec; localTrigCond[1:Int(ceil(S.numEsqTrigLocal*S.numCoveredLocal/S.numLocal))]]
	  # println("localNormEsq = $(localNormEsq[1:Int(ceil(S.numEsqTrigLocal*S.numCoveredLocal/S.numLocal))])")
	  # println("localTrigCond = $(localTrigCond[1:Int(ceil(S.numEsqTrigLocal*S.numCoveredLocal/S.numLocal))])")
	  # now start new TPBVP
	  before_save_time = time_ns()
	  findNewTarget(S, KD, R, hyberBallRad)
	  save_elapsed_time += (time_ns()-before_save_time)/1000000000
	  x0 = R.robotPose[:]
	  xf = S.moveGoal.position[:]
	  before_save_time = time_ns()
	  localTemp = sim_TNNLS_B_CT_Local(x0, xf, S) # sim_TNNLS_B_ETC_relaxedPE_Local#sim_TNNLS_B_Local(x0, xf), sim_TNNLS_B_ETC_robust_relaxedPE_Local
	  save_elapsed_time += (time_ns()-before_save_time)/1000000000
	  localPoseAndKd[1:length(localTemp[1])] = localTemp[1]
	  localNormEsq[1:length(localTemp[2])] = localTemp[2]
	  localTrigCond[1:length(localTemp[3])] = localTemp[3]

	  S.numCoveredLocal = 0
	  S.numLocal = length(localTemp[1])
	  # S.numEsqTrigCovered = 0
	  # S.numEsqTrigCoveredNext = 0
	  S.numEsqTrigLocal = length(localTemp[2])
	else
	  S.moveGoal.isMoveGoal = false
	  S.moveGoal = R.nextMoveTarget
	  S.moveGoal.isMoveGoal = true
	end

	if S.numLocal - S.numCoveredLocal == 0
	  if R.nextMoveTarget.position != S.root.position#R.nextMoveTarget.rrtParentEdge.endNode
	    S.NormEsqvec = [S.NormEsqvec; localNormEsq[1:S.numEsqTrigLocal]]
	    S.TrigCondvec = [S.TrigCondvec; localTrigCond[1:S.numEsqTrigLocal]]
		# println("localNormEsq = $(localNormEsq[1:S.numEsqTrigLocal])")
		# println("localTrigCond = $(localTrigCond[1:S.numEsqTrigLocal])")

		x0 = R.robotPose[:]
		R.nextMoveTarget = R.nextMoveTarget.rrtParentEdge.endNode
		S.moveGoal = R.nextMoveTarget
		xf = S.moveGoal.position[:]
		before_save_time = time_ns()
		localTemp = sim_TNNLS_B_CT_Local(x0, xf, S) #sim_TNNLS_B_ETC_relaxedPE_Local(x0, xf, S)#sim_TNNLS_B_ETC_robust_relaxedPE_Local
		save_elapsed_time += (time_ns()-before_save_time)/1000000000
		localPoseAndKd[1:length(localTemp[1])] = localTemp[1]
		localNormEsq[1:length(localTemp[2])] = localTemp[2]
		localTrigCond[1:length(localTemp[3])] = localTemp[3]
		S.numLocal = length(localTemp[1])#length(localPoseAndKd)
		S.numCoveredLocal = 2
		R.numLocalMovePoints = 2
		R.robotLocalPath[1:R.numLocalMovePoints,:] = [localPoseAndKd[1][1]; localPoseAndKd[2][1]]
		R.nextRobotPose = localPoseAndKd[2][1]

		# S.numEsqTrigCovered = 0
		S.numEsqTrigLocal = length(localTemp[2])
		# S.numEsqTrigCoveredNext = floor(S.numEsqTrigLocal*S.numCoveredLocal/S.numLocal)
		#
	    # R.numLocalEsqTrigPoints = S.numEsqTrigCoveredNext - S.numEsqTrigCovered
  	    # R.localEsq = localNormEsq[S.numEsqTrigCovered+1:S.numEsqTrigCoveredNext]
  	    # R.localTrig = localTrigCond[S.numEsqTrigCovered+1:S.numEsqTrigCoveredNext]

	  else # end
		R.numLocalMovePoints = 1
	    R.robotLocalPath[1:R.numLocalMovePoints,:] = [R.nextMoveTarget.position[1] R.nextMoveTarget.position[2]]
		R.nextRobotPose = R.nextMoveTarget.position

		# R.numLocalEsqTrigPoints = 1
		# R.localEsq = localNormEsq[end]
		# R.localTrig = localTrigCond[end]
		# S.numEsqTrigCoveredNext = S.numEsqTrigLocal
	  end
	elseif S.numLocal - S.numCoveredLocal == 1
	  if R.nextMoveTarget.position != S.root.position#R.nextMoveTarget.rrtParentEdge.endNode
		R.robotLocalPath[1,:] = localPoseAndKd[S.numLocal][1]

		S.NormEsqvec = [S.NormEsqvec; localNormEsq[1:S.numEsqTrigLocal]]
  	    S.TrigCondvec = [S.TrigCondvec; localTrigCond[1:S.numEsqTrigLocal]]
		# println("localNormEsq = $(localNormEsq[1:S.numEsqTrigLocal])")
		# println("localTrigCond = $(localTrigCond[1:S.numEsqTrigLocal])")

		# generate new trajectory
		x0 = localPoseAndKd[S.numLocal][1][:]
		R.nextMoveTarget = R.nextMoveTarget.rrtParentEdge.endNode
		S.moveGoal = R.nextMoveTarget
		xf = S.moveGoal.position[:]
		before_save_time = time_ns()
		localTemp = sim_TNNLS_B_CT_Local(x0, xf, S) #sim_TNNLS_B_ETC_relaxedPE_Local(x0, xf, S)#sim_TNNLS_B_ETC_robust_relaxedPE_Local
		save_elapsed_time += (time_ns()-before_save_time)/1000000000
		localPoseAndKd[1:length(localTemp[1])] = localTemp[1]
		localNormEsq[1:length(localTemp[2])] = localTemp[2]
		localTrigCond[1:length(localTemp[3])] = localTemp[3]
		S.numLocal = length(localTemp[1])#length(localPoseAndKd)
		S.numCoveredLocal = 1
		R.numLocalMovePoints = 2
		R.robotLocalPath[2,:] = localPoseAndKd[S.numCoveredLocal][1]
		R.nextRobotPose = localPoseAndKd[S.numCoveredLocal][1]
		# S.numEsqTrigCovered = floor(S.numEsqTrigLocal*S.numCoveredLocal/S.numLocal)
		S.numEsqTrigLocal = length(localTemp[2])
		# S.numEsqTrigCoveredNext = floor(S.numEsqTrigLocal*S.numCoveredLocal/S.numLocal)
		#
		# R.numLocalEsqTrigPoints = S.numEsqTrigCoveredNext - S.numEsqTrigCovered
		# R.localEsq = localNormEsq[S.numEsqTrigCovered+1:S.numEsqTrigCoveredNext]
		# R.localTrig = localTrigCond[S.numEsqTrigCovered+1:S.numEsqTrigCoveredNext]
	  else # end
		R.numLocalMovePoints = 2
	  	R.robotLocalPath[1,:] = localPoseAndKd[S.numLocal][1]
		R.robotLocalPath[2,:] = [R.nextMoveTarget.position[1] R.nextMoveTarget.position[2]]
	    R.nextRobotPose = R.nextMoveTarget.position
	  end
	elseif S.numLocal - S.numCoveredLocal >= 2
	  # println("length of localPoseAndKd = $(S.numLocal)")
	  R.numLocalMovePoints = 2
	  # R.robotLocalPath[1:R.numLocalMovePoints,:] = [localPoseAndKd[S.numCoveredLocal+1][1]; localPoseAndKd[S.numCoveredLocal+2][1]]
	  R.robotLocalPath[1,:] = localPoseAndKd[S.numCoveredLocal+1][1]
	  R.robotLocalPath[2,:] = localPoseAndKd[S.numCoveredLocal+2][1]
	  R.nextRobotPose = localPoseAndKd[S.numCoveredLocal+2][1]
	  S.numCoveredLocal += 2

	  # S.numEsqTrigCoveredNext = floor(S.numEsqTrigLocal*S.numCoveredLocal/S.numLocal)
	  # R.numLocalEsqTrigPoints = S.numEsqTrigCoveredNext - S.numEsqTrigCovered
	  # R.localEsq = localNormEsq[S.numEsqTrigCovered+1:S.numEsqTrigCoveredNext]
	  # R.localTrig = localTrigCond[S.numEsqTrigCovered+1:S.numEsqTrigCoveredNext]
	else
	  error("error in moveRobot")
	end

	S.kino_dist = localPoseAndKd[S.numCoveredLocal][2]

end

function sim_TNNLS_B_ETC_Local(x1, x2, S) # x1 -> x2, ACC 20

	global eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa10, Wa20, sizePlant, β, L, L1, Wavec, eventTimes

	# System Dynamics for feedback
    n, m = 4, 2
	m_net = 10 # 10 - WAFR
	m_fuel = 30 # 30 - WAFR
	alpha = 0.05 # 0.05 - WAFR
	m_total = m_net + m_fuel # 40 - WAFR
	kx = 20 # 10 - CDC, 20 - WAFR
	ky = 20 # 10 - CDC, 20 - WAFR
	cx = 45 # 10 - CDC, 45 - WAFR
	cy = 45 # 10 - CDC, 45 - WAFR

    A = [0 0 1 0; 0 0 0 1; -kx/m_total 0 -cx/m_total 0; 0 -ky/m_total 0 -cx/m_total] # n x n
    B = [0 0; 0 0; 1/m_total 0; 0 1/m_total] # n x m
    M = 15.0*Matrix(I, n, n) # n x n # 1, .1, 10
    R = .55*Matrix(I, m, m) # m x m # 1, 10, .1
	Pt = 0.5*Matrix(I, n, n)

	# check controllability & observability
    Co = ctrb(A, B)
    unco = size(A, 1) - rank(Co)
    Ob = obsv(sqrt(M), A)
    unob = size(sqrt(M), 1) - rank(Ob)

	if unco + unob > 0	error("system uncontrollable and/or unobservable")	end

    # ODE parameters
    Tf, T, N = 10, 0.05, 200 # finite horizon
    αc = 90 # αa = 1.2
	amplitude, percent = 0.1, 50

    Wc0 = [10.0*ones(Int((n + m)*(n + m + 1)/2 - m*m), 1); reshape(R, (m*m, 1))][:] # (n+m)(n+m+1)/2 x 1; 10 x 1
	Wa10 = 0.5*ones(n,) # n x m; 4 x 2
	Wa20 = 0.5*ones(n,)
    u0 = zeros(m,) # m x 1; 2 x 1
    uf = 0.001*ones(m,)
    Wcfinal = 12.0*ones(Int((n + m)*(n + m + 1)/2),)

	# Trigger parameters
	sizePlant, lixo = size(A)
	eigMinM = minimum(eigvals(M))
	eigMaxR = maximum(eigvals(R))
	eigMinR = minimum(eigvals(R))
	L, β = 30, .6
	L1 = .9*(β*sqrt(eigMinM/eigMaxR))
	αa_UB = (8*eigMinR-4)/(eigMinR+2)
	αa = .25*αa_UB
	finalTime = zeros(0,)

	# x0 = [x1;0.0;0.0]
	x0 = [x1; S.lastVelocity]
	xf = [x2;0.0;0.0]
	Quu, Qxu = 0, 0

    p0 = x0'*M*x0

	Xhat = x0 - xf

    t_save = [0,]
    u_save = u0
	x_save = [[x0; Wc0; p0],]

    u1Delay = interp2PWC(u_save[1], -1, 1) # return an interpolation function
	u2Delay = interp2PWC(u_save[2], -1, 1)
	x1Delay = interp2PWC(getindex.(x_save, 1), -1, 1)
    x2Delay = interp2PWC(getindex.(x_save, 2), -1, 1)
    x3Delay = interp2PWC(getindex.(x_save, 3), -1, 1)
	x4Delay = interp2PWC(getindex.(x_save, 4), -1, 1)

	uvec1, uvec2 = 0, 0
	xdist = norm(x0[1:2] - xf[1:2])
	error = 0.1
	localKd = 0.0
	maxIter = 10000
	poseAndKd = Array{Tuple{Array{Float64,2},Float64}}(undef,0)

	trigCondvec = zeros(0,)
	normESqvec = zeros(0,)
	Wavec = zeros(4,)
	eventTimes = zeros(0,)

    # solve ODEs
	for i = 1:maxIter
		t = ((i-1)*T, i*T)
		p = [xf, uf, Tf, T, A, B, M, R, Pt, percent, amplitude, αa, αc, Quu, Qxu, Wcfinal, u1Delay, u2Delay, x1Delay, x2Delay, x3Delay, x4Delay, x_save, t_save, uvec1, uvec2]#,eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa10, Wa20, sizePlant, β, L, L1, Wavec, eventTimes]#, i,maxIter,start]
        sol = solve(ODEProblem(babyETC_B, x_save[end], t, p), DP5(), reltol = 1e-4, abstol = 1e-4, dtmax = .05)
		t_save = [t_save; sol.t[2:end]] # vcat(t_save, sol.t) save time
		x_save = [x_save; sol.u[2:end]] # vcat(x_save, sol.u) save state
		# uvec1 = p[end-1]
		# uvec2 = p[end]
		uvec1 = p[25]
		uvec2 = p[26]
		# trigCondvec = p[32]
		# normESqvec = p[31]
		# println("$trigCondvec")
        u1Delay = interp2PWC(uvec1, -1, i*T+.01) # interpolate control input
		u2Delay = interp2PWC(uvec2, -1, i*T+.01)
		x1Delay = interp2PWC(getindex.(x_save, 1), -1, i*T+.01)
		x2Delay = interp2PWC(getindex.(x_save, 2), -1, i*T+.01)
		x3Delay = interp2PWC(getindex.(x_save, 3), -1, i*T+.01)
		x4Delay = interp2PWC(getindex.(x_save, 4), -1, i*T+.01)
		trigCondDelay = interp2PWC(trigCondvec, -1, i*T+.01)

		# S.kino_dist is always the max kd till now
		localKd = distPoint2Line(x_save[end][1:2], x1, x2)
		poseAndKd = [poseAndKd; ([x_save[end][1] x_save[end][2]], copy(localKd))]
		if norm(x_save[end][1:2] - xf[1:2]) < error*xdist break end
	end
	S.lastVelocity = x_save[end][3:4]
	(poseAndKd,normESqvec,trigCondvec)#,t_save)#eventTimes)
end

function sim_TNNLS_B_Local(x1, x2) # x1 -> x2
	# System Dynamics for feedback
    n, m = 4, 2
	m_net = 10 # 10 - WAFR
	m_fuel = 30 # 30 - WAFR
	alpha = 0.03 # 0.05 - WAFR
	m_total = m_net + m_fuel # 40 - WAFR
	kx = 20 # 10 - CDC, 20 - WAFR
	ky = 20 # 10 - CDC, 20 - WAFR
	cx = 45 # 10 - CDC, 45 - WAFR
	cy = 45 # 10 - CDC, 45 - WAFR

    A = [0 0 1 0; 0 0 0 1; -kx/m_total 0 -cx/m_total 0; 0 -ky/m_total 0 -cx/m_total] # n x n
    B = [0 0; 0 0; 1/m_total 0; 0 1/m_total] # n x m
    M = 1.0*Matrix(I, n, n) # n x n # 1, .1, 10
    R = 1.0*Matrix(I, m, m) # m x m # 1, 10, .1
	Pt = 0.5*Matrix(I, n, n)

	# check controllability & observability
    Co = ctrb(A, B)
    unco = size(A, 1) - rank(Co)
    Ob = obsv(sqrt(M), A)
    unob = size(sqrt(M), 1) - rank(Ob)

	if unco + unob > 0	error("system uncontrollable and/or unobservable")	end

    # ODE parameters
    Tf, T, N = 45, 0.05, 900 # finite horizon
    αc, αa = 50, 1.2
	amplitude, percent = 0.1, 50

    Wc0 = [10.0*ones(Int((n + m)*(n + m + 1)/2 - m*m), 1); reshape(R, (m*m, 1))][:] # (n+m)(n+m+1)/2 x 1; 10 x 1
	Wa10 = 0.5*ones(n,) # n x m; 4 x 2
	Wa20 = 0.5*ones(n,)
    u0 = 5.0*ones(m,) # m x 1; 2 x 1
    uf = 5.0*ones(m,)
    Wcfinal = 12.0*ones(Int((n + m)*(n + m + 1)/2),)

	x0 = [x1;0.0;0.0]
	xf = [x2;0.0;0.0]
	Quu = 0
	Qxu = 0

    p0 = x0'*M*x0

    t_save = [0,]
    u_save = u0
	x_save = [[x0; Wc0; Wa10; Wa20; p0],]

    u1Delay = interp2PWC(u_save[1], -1, 1) # return an interpolation function
	u2Delay = interp2PWC(u_save[2], -1, 1)
	x1Delay = interp2PWC(getindex.(x_save, 1), -1, 1)
    x2Delay = interp2PWC(getindex.(x_save, 2), -1, 1)
    x3Delay = interp2PWC(getindex.(x_save, 3), -1, 1)
	x4Delay = interp2PWC(getindex.(x_save, 4), -1, 1)

	uvec1, uvec2 = 0, 0
	xdist = norm(x0[1:2] - xf[1:2])
	error = 0.1
	# globalKd = 0.0
	localKd = 0.0
	maxIter = 10000
	poseAndKd = Array{Tuple{Array{Float64,2},Float64}}(undef,0)

    # solve ODEs
	for i = 1:maxIter
		t = ((i-1)*T, i*T)
		p = [xf, uf, Tf, T, A, B, M, R, Pt, percent, amplitude, αa, αc, Quu, Qxu, Wcfinal, u1Delay, u2Delay, x1Delay, x2Delay, x3Delay, x4Delay, x_save, t_save, uvec1, uvec2]#, i,maxIter,start]
        sol = solve(ODEProblem(babyCT_B, x_save[end], t, p), DP5())
		# t_save = [t_save; sol.t[2:end]] # vcat(t_save, sol.t) save time
		x_save = [x_save; sol.u[2:end]] # vcat(x_save, sol.u) save state
		uvec1 = p[end-1]
		uvec2 = p[end]
        u1Delay = interp2PWC(uvec1, -1, i*T+.01) # interpolate control input
		u2Delay = interp2PWC(uvec2, -1, i*T+.01)
		x1Delay = interp2PWC(getindex.(x_save, 1), -1, i*T+.01)
		x2Delay = interp2PWC(getindex.(x_save, 2), -1, i*T+.01)
		x3Delay = interp2PWC(getindex.(x_save, 3), -1, i*T+.01)
		x4Delay = interp2PWC(getindex.(x_save, 4), -1, i*T+.01)
		# S.kino_dist is always the max kd till now
		localKd = distPoint2Line(x_save[end][1:2], x1, x2)
		# if localKd > globalKd
		# 	globalKd = distPoint2Line(x_save[end][1:2], x1, x2)
		# end
		# kino_dist = distPoint2Line(x_save[end][1:2], x1, x2)
		poseAndKd = [poseAndKd; ([x_save[end][1] x_save[end][2]], copy(localKd))]
		# println(sqrt(x_save[end][3]^2 + x_save[end][4]^2))
		if norm(x_save[end][1:2] - xf[1:2]) < error*xdist break end
	end
	# println(length(poseAndKd))
	poseAndKd
end
