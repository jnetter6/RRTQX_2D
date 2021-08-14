# ode fcn for simulation A in journal20 with ETC, H-infinity, and relaxed PE
function babyETC_robust_relaxedPE_journal20_A(dotx, x, p, t)

	global eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa10, Wa20, sizePlant, β, L, L1, Ld, Wavec, Wcvec, Wdvec, numDataRecord, d_real, UkU, Wc, e_buff, Ω, mu, nu, nudelay, xi, basis_on, alpha_basis, epsilon_basis, time_basis, uvec, dvec

	# global xf, uf, Tf, T, A, B, M, R, Pt, per, amp, αa, αc, Quu, Qxu, Wcfinal, uDelay, x1Delay, x2Delay, x3Delay, x_save, t_save, uvec
	xf, T, A, B, F, M, R, γ, Pt, αa, αc, αd, u1Delay, u2Delay, x1Delay, x2Delay, x3Delay, x4Delay, dDelay, pDelay = p#,eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa10, Wa20, sizePlant, β, L, L1, Wavec, eventTimes=p#,i,maxIter,start = p # uu, vv useless here percent, amplitude,
    n, m, q = 4, 2, 1
	Wc = x[(n+1) : Int(n+(n+m+q)*(n+m+q+1)/2)] # 5-32
	Wd = x[Int(n+(n+m+q)*(n+m+q+1)/2+1) : Int(n+(n+m+q)*(n+m+q+1)/2+n)] # 33-36
	P = x[end] # 37
	d_real = 0#.3*sin(t)

	# Basis
	if basis_on == 1
		r = time_basis - t
		rdelay = time_basis - (t - T)
		for i = 1:Int((n+m+q)*(n+m+q+1)/2)
			for j = 1:Int((n+m+q)*(n+m+q+1)/2)
				if i+j < 23
					nu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-1) # alpha_basis*tanh(epsilon_basis*r)^(i+j-1); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-1); % cos(2*pi*(r)/10)^(i+j-1); %
					nudelay[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-1)
				elseif i+j == 23
					nu[i,j] = alpha_basis*tanh(epsilon_basis*r) # alpha_basis*tanh(epsilon_basis*r); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2); % cos(2*pi*(r)/10); %
					nudelay[i,j] = alpha_basis*tanh(epsilon_basis*r)
				else #if i+j > 23
					nu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-(((n+m+q)*(n+m+q+1)/2)+1)) # alpha_basis*tanh(epsilon_basis*r)^(i+j-22); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-22) ; % cos(2*pi*(r)/10)^(i+j-22); %
					nudelay[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-(((n+m+q)*(n+m+q+1)/2)+1))
				end
			end
		end
		for i = 1:n
			for j = 1:n
				if i+j < 6
					mu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-1) # alpha_basis*tanh(epsilon_basis*r)^(i+j-1); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-1); %cos(2*pi*(r)/10)^(i+j-1); %
					xi[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-1)
				elseif i+j == 6
					mu[i,j] = alpha_basis*tanh(epsilon_basis*r) # alpha_basis*tanh(epsilon_basis*r); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2); % cos(2*pi*(r)/100); %
					xi[i,j] = alpha_basis*tanh(epsilon_basis*r)
				else #if i+j > 6
					mu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-(n+1)) # alpha_basis*tanh(epsilon_basis*r)^(i+j-5); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-5) ; % cos(2*pi*(r)/10)^(i+j-5); %
					xi[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-(n+1))
				end
			end
		end
	end

	# Update control
    ud = zeros(m,)
	# ud[1] = Wa10'*(x[1:n] - xf)
	# ud[2] = Wa20'*(x[1:n] - xf)
	ud[1] = Wa10'*mu*Xhat
	ud[2] = Wa20'*mu*Xhat

	# Update disturbance
	d = Wd'*xi*(x[1:4] - xf)

	# Delays
    uddelay = zeros(m,)
	uddelay[1] = u1Delay(t - T)
	uddelay[2] = u2Delay(t - T)
    xdelay = zeros(n,)
	xdelay[1] = x1Delay(t - T)
    xdelay[2] = x2Delay(t - T)
    xdelay[3] = x3Delay(t - T)
	xdelay[4] = x4Delay(t - T)
	ddelay = dDelay(t - T)
	pdelay = pDelay(t - T)

	# Kronecker products
	U = vcat(x[1:n] - xf, ud[1], ud[2], d) # augmented state
    UkU = vcat(U[1]^2, U[1]*U[2], U[1]*U[3], U[1]*U[4], U[1]*ud[1], U[1]*ud[2], U[1]*d,
               U[2]^2, U[2]*U[3], U[2]*U[4], U[2]*ud[1], U[2]*ud[2], U[2]*d,
               U[3]^2, U[3]*U[4], U[3]*ud[1], U[3]*ud[2], U[3]*d,
			   U[4]^2, U[4]*ud[1], U[4]*ud[2], U[4]*d,
               ud[1]^2, ud[1]*ud[2], ud[1]*d,
			   ud[2]^2, ud[2]*d,
			   d^2)
    UkUdelay = vcat(xdelay[1]^2, xdelay[1]*xdelay[2], xdelay[1]*xdelay[3], xdelay[1]*xdelay[4], xdelay[1]*uddelay[1], xdelay[1]*uddelay[2], xdelay[1]*ddelay,
                    xdelay[2]^2, xdelay[2]*xdelay[3], xdelay[2]*xdelay[4], xdelay[2]*uddelay[1], xdelay[2]*uddelay[2], xdelay[2]*ddelay,
                    xdelay[3]^2, xdelay[3]*xdelay[4], xdelay[3]*uddelay[1], xdelay[3]*uddelay[2], xdelay[3]*ddelay,
					xdelay[4]^2, xdelay[4]*uddelay[1], xdelay[4]*uddelay[2], xdelay[4]*ddelay,
                    uddelay[1]^2, uddelay[1]*uddelay[2], uddelay[1]*ddelay,
					uddelay[2]^2, uddelay[2]*ddelay,
					ddelay^2)
    Quu = [Wc[23] Wc[24]; Wc[24] Wc[26]] # m x m 24-27
    Quu_inv = inv(Quu)
    Qux = [Wc[5] Wc[11] Wc[16] Wc[20]; Wc[6] Wc[12] Wc[17] Wc[21]] #reshape(Wc[Int((n+m+q)*(n+m+q+1)/2-m^2-q^2-m*n+1):Int((n+m+q)*(n+m+q+1)/2-m^2-q^2)], (n, m)) # n x m 16-23
    Qxu = Qux' # m x n
    Qdx = [Wc[7] Wc[13] Wc[18] Wc[22]]#reshape(Wc[Int((n+m+q)*(n+m+q+1)/2-m^2-q^2-m*n-n*q+1):Int((n+m+q)*(n+m+q+1)/2-m^2-q^2-m*n)], (n, q)) # n x q 12-15
    Qxd = Qdx'
    # Qdd = Wc[end] # q x q
	Qdd_inv = 1.0/Wc[end]

	# Triggering condition
	normXSq = abs(U[1])^2 + abs(U[2])^2 + U[3]^2 + U[4]^2
	trigCond = (((1 - β^2)*eigMinM - γ^2*Ld^2)/(2*L^2*eigMaxR))*normXSq + ((1-β^2)*eigMinR/(2*L^2*eigMaxR))*(abs(ud[1])^2+abs(ud[2])^2)

	# Integral reinforcement dynamics
	# dP = 0.5*((x[1:n] - xf)'*M*(x[1:n] - xf) - (xdelay - xf)'*M*(xdelay - xf) + ud'*R*ud - uddelay'*R*uddelay - γ^2*d^2 + γ^2*ddelay^2)
	dP = 0.5*(U[1:n]'*M*U[1:n] + ud'*R*ud - γ^2*d_real^2) # d_real

	# approximation errors
	σ = nu*UkU - nudelay*UkUdelay 
	σ_f	= nu*UkU
    ec = P - pdelay + Wc'*σ #P + Wc'*σ
	ecfinal = 0.5*U[1:n]'*Pt*U[1:n] - Wc'*σ_f
	ed = Wd'*xi*U[1:n] + (Qdd_inv*Qdx*U[1:n])[1] # q x 1

	# gap
	e = Xhat - U[1:n]
	normESq = norm(e)^2

	# Critic dynamics
	relaxedPETerm = zeros(Int((n+m+q)*(n+m+q+1)/2),)
	for rrr in 1:numDataRecord
		relaxedPETerm = relaxedPETerm + e_buff[rrr]*Ω[:,rrr]./((1 + Ω[:,rrr]'*Ω[:,rrr]).^2)
	end
	dWc = -αc*((σ./(σ'*σ+1).^2)*ec + (σ_f./(σ_f'*σ_f+1).^2)*ecfinal + relaxedPETerm)

	# Triggering updates
	if normESq >= 0.95*trigCond # flows
	    Xhat = U[1:n]
		ea1 = Wa10'*mu*Xhat + [Quu_inv[1,:]'*Qux[:,1]; Quu_inv[1,:]'*Qux[:,2]; Quu_inv[1,:]'*Qux[:,3]; Quu_inv[1,:]'*Qux[:,4]]'*Xhat
		ea2 = Wa20'*mu*Xhat + [Quu_inv[2,:]'*Qux[:,1]; Quu_inv[2,:]'*Qux[:,2]; Quu_inv[2,:]'*Qux[:,3]; Quu_inv[2,:]'*Qux[:,4]]'*Xhat
		Wa10 = Wa10 - αa*mu*Xhat*ea1'/(1 + (mu*Xhat)'*(mu*Xhat))
		Wa20 = Wa20 - αa*mu*Xhat*ea2'/(1 + (mu*Xhat)'*(mu*Xhat))
		Wavec = [Wavec [Wa10; Wa20]]
		eventTimes = [eventTimes; t]
	end

	# Disturbance dynamics
	dWd = -αd*xi*U[1:n]*ed'/(1 + (xi*U[1:n])'*(xi*U[1:n]))

	# System dynamics
	dx = A*U[1:n] + B*ud + F*d_real

	# save vecs
	# p[end] = unew
	# p[end-1:end] .= unew
	uvec = [uvec; [ud]]
	normESqvec = [normESqvec; normESq]
	# p[31] = [p[31]; normESq]
	trigCondvec = [trigCondvec; trigCond]
	# p[32] = [p[32]; trigCond]
	# println("$trigCondvec")
	dvec = [dvec; d]
	# dd = d

	dotx .= vcat(dx, dWc, dWd, dP)
end

function sim_journal20_A_ETC_robust_relaxedPE_Local(x0, xf)#, S) # lastVelocity # x1 -> x2, Journal 20

	global eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa10, Wa20, sizePlant, β, L, L1, Ld, Wavec, Wcvec, Wdvec, numDataRecord, d_real, UkU, Wc, e_buff, Ω, mu, nu, nudelay, xi, basis_on, alpha_basis, epsilon_basis, time_basis, uvec, dvec

	# System Dynamics for feedback
    n, m, q = 4, 2, 1
	m_net = 10 # 10 - WAFR
	m_fuel = 30 # 30 - WAFR
	alpha = 0.05 # 0.05 - WAFR
	m_total = m_net + m_fuel # 40 - WAFR
	kx = 20 # 10 - CDC, 20 - WAFR
	ky = 20 # 10 - CDC, 20 - WAFR
	cx = 45 # 10 - CDC, 45 - WAFR
	cy = 45 # 10 - CDC, 45 - WAFR

    A = [0 0 1 0; 0 0 0 1; -kx/m_total 0 -cx/m_total 0; 0 -ky/m_total 0 -cx/m_total]# [0.883537   1.10463   1.68383  -3.77859; -1.37899    1.51254  -4.4073    2.13532; 1.78866   -4.57991   3.03146  -0.864618; 3.80966    3.99251   1.75758  -2.51588] # n x n
    B = [0 0; 0 0; 1/m_total 0; 0 1/m_total] # n x m [0. 0.; 0. 0.; 1. 0.; 0. 1.]
    F = [0.; 0.; 1.; 1.] # n x q
    M = 1.0*Matrix(I, n, n) # n x n # 1, .1, 10
    R = 0.1*Matrix(I, m, m) # m x m # 1, 10, .1
	γ = .1
	Pt = 0.5*Matrix(I, n, n)

	# check controllability & observability
    Co = ctrb(A, B)
    unco = size(A, 1) - rank(Co)
    Ob = obsv(sqrt(M), A)
    unob = size(sqrt(M), 1) - rank(Ob)

	if unco + unob > 0	error("system uncontrollable and/or unobservable")	end

    # ODE parameters
    Tf, T, N = 10, 0.05, 200 # finite horizon
	αc, αd = 90, 3 # αa = 1.2
	
	# goal orientation
	# theta = atan((x2-x1)[2]/(x2-x1)[1])

	# Wc0 = 10*rand(28)[:]
	Wc0 = [10.0*ones(Int((n+m+q)*(n+m+q+1)/2 - m*(m+1)/2 - m*q - q^2), 1); .525; 0; 0; .525; 0; -γ^2][:] # (n+m+q)(n+m+q+1)/2 x 1; 28 x 1
	# Wc0 = 20*[0.328555; 0.756487; 0.273798; 0.309953; 0.438633; 0.307851; 0.335527; 0.0851001; 0.608546; 0.282202; 0.399089; 0.281325; 0.528826; 0.234641; 0.584131; 0.313394; 0.897723; 0.675292; 0.73439; 0.445563; 0.592254; 0.596938; 0.634189; 0.332015; 0.714059; 0.511094; 0.0298947; 0.0722483]
	# Wc0 = [10.0*ones(11); 6.0*ones(11); .525; 0; 0; .525; 0; -γ^2][:] # (n+m+q)(n+m+q+1)/2 x 1; 28 x 1
	# Wc0 = [10.0*ones(4); -5; 5; 10; 10.0*ones(3); -5; 5; 10; 10.0*ones(2); -5; 5; 10; 10; -5; 5; 10; .525; 0; 10; .525; 10; -1];
	# Wc0 = [10.0*ones(4); 5cos(theta); 5sin(theta); 0.5; 10.0*ones(3); 5cos(theta); 5sin(theta); 0.5; 10.0*ones(2); 5cos(theta); 5sin(theta); 0.5; 10.0; 5cos(theta); 5sin(theta); 0.5; .525; 0; 0; .525; 0; -γ^2];
	
	# Wa10 = 5*cos(theta)*ones(n,) # n x m; 4 x 2
	# Wa20 = 5*sin(theta)*ones(n,)
	Wa10 = .5*ones(n,) # n x m; 4 x 2
	Wa20 = .5*ones(n,)
	Wd0 = 0.5*ones(n,)
    u0 = zeros(m,) # m x 1; 2 x 1
    d0 = 0
    # d_real = 1

	# Trigger parameters
	sizePlant, lixo = size(A)
	eigMinM = minimum(eigvals(M))
	eigMaxR = maximum(eigvals(R))
	eigMinR = minimum(eigvals(R))
	L, Ld = 10, 10
	β = 1 #β1
	L1 = 10 # .9*(β*sqrt(eigMinM/eigMaxR))
	αa_UB = (8*eigMinR-4)/(eigMinR+2)
	αa = .25*αa_UB
	finalTime = zeros(0,)

	# x0 = [x1;0.0;0.0]
	# x0 = [x1; S.lastVelocity] # lastVelocity
	# xf = [x2; 0.0; 0.0]

    # p0 = (x0 - xf)'*M*(x0 - xf)

	Xhat = x0 - xf

    t_save = [0,]
	x_save = [[x0; Wc0; Wd0; 0],]
    u_save = [u0,] # u at every T
	uvec = [u0,] # all u
	d_save = [d0,] # d at every T
	dvec = [d0,] # all d

    u1Delay = interp2PWC(getindex.(u_save, 1), -1, 1) # return an interpolation function
	u2Delay = interp2PWC(getindex.(u_save, 2), -1, 1)
	x1Delay = interp2PWC(getindex.(x_save, 1), -1, 1)
    x2Delay = interp2PWC(getindex.(x_save, 2), -1, 1)
    x3Delay = interp2PWC(getindex.(x_save, 3), -1, 1)
	x4Delay = interp2PWC(getindex.(x_save, 4), -1, 1)
	dDelay = interp2PWC(d_save, -1, 1)
	pDelay = interp2PWC(getindex.(x_save, 37), -1, 1)

	# xdist = norm(x0[1:2] - xf[1:2])
	# error = 0.2
	# localKd = 0.0
	maxIter = 10000
	poseAndKd = Array{Tuple{Array{Float64,2},Float64}}(undef,0)

	trigCondvec = zeros(0)
	normESqvec = zeros(0)
	Wavec = zeros(8,)
	Wcvec = zeros(28,)
	Wdvec = zeros(4,)
	eventTimes = zeros(0,)

	numDataRecord = 0
	e_buff = Array{Float64, 1}(undef, 0)
	halfvec = [1; zeros(6); 1; zeros(5); 1; 0;0;0;0; 1; 0;0;0; .5; 0;0; .5; 0; -γ^2]
	Ω = Array{Float64,2}(undef, Int((n+m+q)*(n+m+q+1)/2), 0)

	# Basis parameters
	epsilon_basis = .1
	alpha_basis = 5
	time_basis = 20
	basis_on = 1 # 1-on, 2-off
	nu = zeros(Int((n+m+q)*(n+m+q+1)/2), Int((n+m+q)*(n+m+q+1)/2))
	nudelay = zeros(Int((n+m+q)*(n+m+q+1)/2), Int((n+m+q)*(n+m+q+1)/2))
	mu = zeros(n, n)
	xi = zeros(n, n)

    # solve ODEs
	for i = 1:maxIter
		t = ((i-1)*T, i*T)
		p = [xf, T, A, B, F, M, R, γ, Pt, αa, αc, αd, u1Delay, u2Delay, x1Delay, x2Delay, x3Delay, x4Delay, dDelay, pDelay]#,eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa10, Wa20, sizePlant, β, L, L1, Wavec, eventTimes]#, i,maxIter,start] percent, amplitude,
        sol = solve(ODEProblem(babyETC_robust_relaxedPE_journal20_A, x_save[end], t, p), DP5())#, reltol = 1e-4, abstol = 1e-4, dtmax = .05)
		t_save = [t_save; sol.t[2:end]] # vcat(t_save, sol.t) save time
		x_save = [x_save; sol.u[2:end]] # vcat(x_save, sol.u) save state
		u_save = [u_save; [uvec[end]]]
		d_save = [d_save; dvec[end]]
		# integral = x_save[end][end] - x_save[end-1][end]
		# vecInt = [vecInt; integral]
		Wcvec = [Wcvec sol.u[end][5:32]]
		Wdvec = [Wdvec sol.u[end][33:36]]
		# uvec1 = p[end-2]
		# uvec2 = p[end-1]
		# trigCondvec = p[32]
		# normESqvec = p[31]
		# println("$trigCondvec")
        u1Delay = interp2PWC(getindex.(u_save, 1), -1, i*T+.01) # interpolate control input
		u2Delay = interp2PWC(getindex.(u_save, 2), -1, i*T+.01)
		x1Delay = interp2PWC(getindex.(x_save, 1), -1, i*T+.01)
		x2Delay = interp2PWC(getindex.(x_save, 2), -1, i*T+.01)
		x3Delay = interp2PWC(getindex.(x_save, 3), -1, i*T+.01)
		x4Delay = interp2PWC(getindex.(x_save, 4), -1, i*T+.01)
		dDelay = interp2PWC(d_save, -1, i*T+.01)
		pDelay = interp2PWC(getindex.(x_save, 37), -1, i*T+.01)
		# trigCondDelay = interp2PWC(trigCondvec, -1, i*T+.01)

        # relaxed PE: record data
        if rank(Ω) < n+m+q
			∂nu∂t = zeros(Int((n+m+q)*(n+m+q+1)/2), Int((n+m+q)*(n+m+q+1)/2))
			if basis_on == 1
				r = time_basis - i*T
				for i = 1:Int((n+m+q)*(n+m+q+1)/2)
					for j = 1:Int((n+m+q)*(n+m+q+1)/2)
						if i+j < 23
							∂nu∂t[i,j] = alpha_basis*(i+j-1)*(1-tanh(epsilon_basis*r)^2)*(-epsilon_basis)*tanh(epsilon_basis*r)^(i+j-2) # alpha_basis*tanh(epsilon_basis*r)^(i+j-1); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-1); % cos(2*pi*(r)/10)^(i+j-1); %
						elseif i+j == 23
							∂nu∂t[i,j] = alpha_basis*(1-tanh(epsilon_basis*r)^2)*(-epsilon_basis) # alpha_basis*tanh(epsilon_basis*r); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2); % cos(2*pi*(r)/10); %
						else #if i+j > 23
							∂nu∂t[i,j] = alpha_basis*(i+j-(((n+m+q)*(n+m+q+1)/2)+1))*(1-tanh(epsilon_basis*r)^2)*(-epsilon_basis)*tanh(epsilon_basis*r)^(i+j-(((n+m+q)*(n+m+q+1)/2)+1)-1) # alpha_basis*tanh(epsilon_basis*r)^(i+j-22); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-22) ; % cos(2*pi*(r)/10)^(i+j-22); %
						end
					end
				end
			end
			x_dot = (x_save[end][1:4] - x_save[end-1][1:4])/T
			u_dot = (u_save[end] - u_save[end-1])/T
			d_dot = (d_save[end] - d_save[end-1])/T
			U = [x_save[end][1:4] - xf; u_save[end]; d_save[end]]
			ud = u_save[end]
			d = d_save[end]
			∂UkU∂t = [2*U[1]*x_dot[1]; x_dot[1]*U[2] + U[1]*x_dot[2]; x_dot[1]*U[3] + U[1]*x_dot[3]; x_dot[1]*U[4] + U[1]*x_dot[4]; x_dot[1]*ud[1] + U[1]*u_dot[1]; x_dot[1]*ud[2] + U[1]*u_dot[2]; x_dot[1]*d + U[1]*d_dot;
			 		  2*U[2]*x_dot[2]; x_dot[2]*U[3] + U[2]*x_dot[3]; x_dot[2]*U[4] + U[2]*x_dot[4]; x_dot[2]*ud[1] + U[2]*u_dot[1]; x_dot[2]*ud[2] + U[2]*u_dot[2]; x_dot[2]*d + U[2]*d_dot;
				      2*U[3]*x_dot[3]; x_dot[3]*U[4] + U[3]*x_dot[4]; x_dot[3]*ud[1] + U[3]*u_dot[1]; x_dot[3]*ud[2] + U[3]*u_dot[2]; x_dot[3]*d + U[3]*d_dot;
					  2*U[4]*x_dot[4]; x_dot[4]*ud[1] + U[4]*u_dot[1]; x_dot[4]*ud[2] + U[4]*u_dot[2]; x_dot[4]*d + U[4]*d_dot;
					  2*ud[1]*u_dot[1]; u_dot[1]*ud[2] + ud[1]*u_dot[2]; u_dot[1]*d + ud[1]*d_dot;
					  2*ud[2]*u_dot[2]; u_dot[2]*d + ud[2]*d_dot;
					  2*d*d_dot]
			ω = ∂nu∂t*UkU + nu*∂UkU∂t
			e_buff = [e_buff; (.5*halfvec'*UkU + Wc'*ω - .5*([Wa10 Wa20]'*((x_save[end][1:4] - xf) - Xhat))'*R*([Wa10 Wa20]'*((x_save[end][1:4] - xf) - Xhat)))];
			Ω = [Ω ω]
            numDataRecord += 1
        end

		# S.kino_dist is always the max kd till now
		# localKd = distPoint2Line(x_save[end][1:2], x1, x2)
		poseAndKd = [poseAndKd; ([x_save[end][1] x_save[end][2] x_save[end][3] x_save[end][4]], 0)]#copy(localKd))]
		# if norm(x_save[end][1:2] - xf[1:2]) < error*xdist break end
	end
	# S.lastVelocity = x_save[end][3:4]
	# lastVelocity = x_save[end][3:4]
	(poseAndKd, normESqvec, trigCondvec)#,t_save)#eventTimes)
	# println(length(normESqvec))
	# println(length(trigCondvec))
	# plot(1:length(normESqvec), normESqvec)
	# plot!(1:length(normESqvec), trigCondvec)
	# println(length(t_save))
	# plot(getindex.(getindex.(poseAndKd, 1), 1), getindex.(getindex.(poseAndKd, 1), 2))
	# maximum(getindex.(poseAndKd, 2))
	# println(distPoint2Line(poseAndKd[end][1][1:2], x1, x2))
	# plot(t_save, getindex.(x_save, 37)) # check cost
	# plot(T:T:T*length(vecInt), vecInt) # check integral
	# plot(eventTimes, Wavec[:,2:end]')
	# plot(1:size(Wcvec)[2], Wcvec',label=0)
	# plot(1:size(Wdvec)[2], Wdvec',label=0)
end

function babyETC_robust_relaxedPE_journal_A_trynormalization(dotx, x, p, t)

	global ud, eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa10, Wa20, sizePlant, β, L, L1, Ld, Wavec, Wcvec, Wdvec, numDataRecord, d_real, UkU, Wc, e_buff, Ω, mu, nu, nudelay, xi, basis_on, alpha_basis, epsilon_basis, time_basis, uvec, dvec, tvec

	# global xf, uf, Tf, T, A, B, M, R, Pt, per, amp, αa, αc, Quu, Qxu, Wcfinal, uDelay, x1Delay, x2Delay, x3Delay, x_save, t_save, uvec
	xf, T, A, B, F, M, R, γ, Pt, αa, αc, αd, u1Delay, u2Delay, x1Delay, x2Delay, x3Delay, x4Delay, dDelay, pDelay = p#,eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa10, Wa20, sizePlant, β, L, L1, Wavec, eventTimes=p#,i,maxIter,start = p # uu, vv useless here percent, amplitude,
    n, m, q = 4, 2, 1
	Wc = x[(n+1) : Int(n+(n+m+q)*(n+m+q+1)/2)] # 5-32
	Wd = x[Int(n+(n+m+q)*(n+m+q+1)/2+1) : Int(n+(n+m+q)*(n+m+q+1)/2+n)] # 33-36
	P = x[end] # 37
	d_real = 0#sin(t)#0#0.1

	# Basis
	if basis_on == 1
		r = time_basis - t
		rdelay = time_basis - (t - T)
		for i = 1:Int((n+m+q)*(n+m+q+1)/2)
			for j = 1:Int((n+m+q)*(n+m+q+1)/2)
				if i+j < 23
					nu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-1) # alpha_basis*tanh(epsilon_basis*r)^(i+j-1); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-1); % cos(2*pi*(r)/10)^(i+j-1); %
					nudelay[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-1)
				elseif i+j == 23
					nu[i,j] = alpha_basis*tanh(epsilon_basis*r) # alpha_basis*tanh(epsilon_basis*r); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2); % cos(2*pi*(r)/10); %
					nudelay[i,j] = alpha_basis*tanh(epsilon_basis*r)
				else #if i+j > 23
					nu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-(((n+m+q)*(n+m+q+1)/2)+1)) # alpha_basis*tanh(epsilon_basis*r)^(i+j-22); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-22) ; % cos(2*pi*(r)/10)^(i+j-22); %
					nudelay[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-(((n+m+q)*(n+m+q+1)/2)+1))
				end
			end
		end
		for i = 1:n
			for j = 1:n
				if i+j < 6
					mu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-1) # alpha_basis*tanh(epsilon_basis*r)^(i+j-1); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-1); %cos(2*pi*(r)/10)^(i+j-1); %
					xi[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-1)
				elseif i+j == 6
					mu[i,j] = alpha_basis*tanh(epsilon_basis*r) # alpha_basis*tanh(epsilon_basis*r); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2); % cos(2*pi*(r)/100); %
					xi[i,j] = alpha_basis*tanh(epsilon_basis*r)
				else #if i+j > 6
					mu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-(n+1)) # alpha_basis*tanh(epsilon_basis*r)^(i+j-5); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-5) ; % cos(2*pi*(r)/10)^(i+j-5); %
					xi[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-(n+1))
				end
			end
		end
	end

	# Update disturbance
	d = Wd'*xi*(x[1:4] - xf)

	# Delays
    uddelay = zeros(m,)
	uddelay[1] = u1Delay(t - T)
	uddelay[2] = u2Delay(t - T)
    xdelay = zeros(n,)
	xdelay[1] = x1Delay(t - T)
    xdelay[2] = x2Delay(t - T)
    xdelay[3] = x3Delay(t - T)
	xdelay[4] = x4Delay(t - T)
	ddelay = dDelay(t - T)
	pdelay = pDelay(t - T)

	# Kronecker products
	U = vcat(x[1:n] - xf, ud[1], ud[2], d) # augmented state
    UkU = vcat(U[1]^2, U[1]*U[2], U[1]*U[3], U[1]*U[4], U[1]*ud[1], U[1]*ud[2], U[1]*d,
               U[2]^2, U[2]*U[3], U[2]*U[4], U[2]*ud[1], U[2]*ud[2], U[2]*d,
               U[3]^2, U[3]*U[4], U[3]*ud[1], U[3]*ud[2], U[3]*d,
			   U[4]^2, U[4]*ud[1], U[4]*ud[2], U[4]*d,
               ud[1]^2, ud[1]*ud[2], ud[1]*d,
			   ud[2]^2, ud[2]*d,
			   d^2)
    UkUdelay = vcat(xdelay[1]^2, xdelay[1]*xdelay[2], xdelay[1]*xdelay[3], xdelay[1]*xdelay[4], xdelay[1]*uddelay[1], xdelay[1]*uddelay[2], xdelay[1]*ddelay,
                    xdelay[2]^2, xdelay[2]*xdelay[3], xdelay[2]*xdelay[4], xdelay[2]*uddelay[1], xdelay[2]*uddelay[2], xdelay[2]*ddelay,
                    xdelay[3]^2, xdelay[3]*xdelay[4], xdelay[3]*uddelay[1], xdelay[3]*uddelay[2], xdelay[3]*ddelay,
					xdelay[4]^2, xdelay[4]*uddelay[1], xdelay[4]*uddelay[2], xdelay[4]*ddelay,
                    uddelay[1]^2, uddelay[1]*uddelay[2], uddelay[1]*ddelay,
					uddelay[2]^2, uddelay[2]*ddelay,
					ddelay^2)
    Quu = [Wc[23] Wc[24]; Wc[24] Wc[26]] # m x m 24-27
    Quu_inv = inv(Quu)
    Qux = [Wc[5] Wc[11] Wc[16] Wc[20]; Wc[6] Wc[12] Wc[17] Wc[21]] #reshape(Wc[Int((n+m+q)*(n+m+q+1)/2-m^2-q^2-m*n+1):Int((n+m+q)*(n+m+q+1)/2-m^2-q^2)], (n, m)) # n x m 16-23
    Qxu = Qux' # m x n
    Qdx = [Wc[7] Wc[13] Wc[18] Wc[22]]#reshape(Wc[Int((n+m+q)*(n+m+q+1)/2-m^2-q^2-m*n-n*q+1):Int((n+m+q)*(n+m+q+1)/2-m^2-q^2-m*n)], (n, q)) # n x q 12-15
    Qxd = Qdx'
    # Qdd = Wc[end] # q x q
	Qdd_inv = 1.0/Wc[end]

	# Triggering condition
	normXSq = abs(U[1])^2 + abs(U[2])^2 + U[3]^2 + U[4]^2
	# trigCond = (((1 - β^2)*eigMinM - γ^2*Ld^2)/(2*L^2*eigMaxR))*normXSq + ((1-β^2)*eigMinR/(2*L^2*eigMaxR))*(abs(ud[1])^2+abs(ud[2])^2)
	trigCond = ((1 - β^2)*eigMinM*normXSq + (1-β^2)*eigMinR*(abs(ud[1])^2+abs(ud[2])^2) - γ^2*d^2)/(4*eigMaxR*(L^2+L1^2)) # d_real

	# Integral reinforcement dynamics
	# dP = 0.5*((x[1:n] - xf)'*M*(x[1:n] - xf) - (xdelay - xf)'*M*(xdelay - xf) + ud'*R*ud - uddelay'*R*uddelay - γ^2*d^2 + γ^2*ddelay^2)
	dP = 0.5*(U[1:n]'*M*U[1:n] + ud'*R*ud - γ^2*d_real^2) # d

	# approximation errors
	σ = nu*UkU - nudelay*UkUdelay 
	σ_f	= nu*UkU
    ec = P - pdelay + Wc'*σ #P + Wc'*σ
	ecfinal = 0.5*U[1:n]'*Pt*U[1:n] - Wc'*σ_f
	ed = Wd'*xi*U[1:n] + (Qdd_inv*Qdx*U[1:n])[1] # q x 1

	# gap
	e = Xhat - U[1:n]
	normESq = norm(e)^2

	# Critic dynamics
	relaxedPETerm = zeros(Int((n+m+q)*(n+m+q+1)/2),)
	for rrr in 1:numDataRecord
		relaxedPETerm = relaxedPETerm + e_buff[rrr]*Ω[:,rrr]./((1 + Ω[:,rrr]'*Ω[:,rrr]).^2)
	end
	dWc = -αc*((σ./(σ'*σ+1).^2)*ec + (σ_f./(σ_f'*σ_f+1).^2)*ecfinal + relaxedPETerm)

	# Triggering updates
	if normESq >= 0.95*trigCond # flows
	    Xhat = U[1:n]
		ea1 = Wa10'*mu*Xhat + [Quu_inv[1,:]'*Qux[:,1]; Quu_inv[1,:]'*Qux[:,2]; Quu_inv[1,:]'*Qux[:,3]; Quu_inv[1,:]'*Qux[:,4]]'*Xhat
		ea2 = Wa20'*mu*Xhat + [Quu_inv[2,:]'*Qux[:,1]; Quu_inv[2,:]'*Qux[:,2]; Quu_inv[2,:]'*Qux[:,3]; Quu_inv[2,:]'*Qux[:,4]]'*Xhat
		Wa10 = Wa10 - αa*mu*Xhat*ea1'/(1 + (mu*Xhat)'*(mu*Xhat))
		Wa20 = Wa20 - αa*mu*Xhat*ea2'/(1 + (mu*Xhat)'*(mu*Xhat))
		Wavec = [Wavec [Wa10; Wa20]]
		eventTimes = [eventTimes; t]
		# Update control
		ud = [Wa10'*mu*Xhat; Wa20'*mu*Xhat]
		# ud = zeros(m,)
		# ud[1] = Wa10'*mu*Xhat
		# ud[2] = Wa20'*mu*Xhat
	end

	# Disturbance dynamics
	dWd = -αd*xi*U[1:n]*ed'/(1 + (xi*U[1:n])'*(xi*U[1:n]))

	# System dynamics
	dx = A*U[1:n] + B*ud + F*d_real

	# save vecs
	# p[end] = unew
	# p[end-1:end] .= unew
	uvec = [uvec; [ud]]
	normESqvec = [normESqvec; normESq]
	# p[31] = [p[31]; normESq]
	trigCondvec = [trigCondvec; trigCond]
	tvec = [tvec; t]

	# p[32] = [p[32]; trigCond]
	# println("$trigCondvec")
	dvec = [dvec; d]
	# dd = d

	dotx .= vcat(dx, dWc, dWd, dP)
end

function sim_journal_A_trynormalization_ETC_robust_relaxedPE_Local(x1, x2)#, S) # lastVelocity # x1 -> x2, Journal 20

	global ud, eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa10, Wa20, sizePlant, β, L, L1, Ld, Wavec, Wcvec, Wdvec, numDataRecord, d_real, UkU, Wc, e_buff, Ω, mu, nu, nudelay, xi, basis_on, alpha_basis, epsilon_basis, time_basis, uvec, dvec, tvec

	# System Dynamics for feedback
    n, m, q = 4, 2, 1
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
    F = [0; 0; 1; 1] # n x q
    M = 5.0*Matrix(I, n, n) # n x n # 1, .1, 10
    R = 0.3*Matrix(I, m, m) # m x m # 1, 10, .1
	γ = 1
	Pt = 0.3*Matrix(I, n, n)

	# check controllability & observability
    Co = ctrb(A, B)
    unco = size(A, 1) - rank(Co)
    Ob = obsv(sqrt(M), A)
    unob = size(sqrt(M), 1) - rank(Ob)

	if unco + unob > 0	error("system uncontrollable and/or unobservable")	end

    # ODE parameters
	Tf, T, N = 10, 0.05, 300 # finite horizon
	αc, αd = 90, 3 # αa = 1.2
	
	# goal orientation
	# theta = atan((x2-x1)[2]/(x2-x1)[1])

	# Wc0 = 10*rand(28)[:]
	Wc0 = [5.0*ones(Int((n+m+q)*(n+m+q+1)/2 - m*(m+1)/2 - m*q - q^2), 1); .3; 0; 0; .3; 0; -γ^2][:] # (n+m+q)(n+m+q+1)/2 x 1; 28 x 1
	# Wc0 = 20*[0.328555; 0.756487; 0.273798; 0.309953; 0.438633; 0.307851; 0.335527; 0.0851001; 0.608546; 0.282202; 0.399089; 0.281325; 0.528826; 0.234641; 0.584131; 0.313394; 0.897723; 0.675292; 0.73439; 0.445563; 0.592254; 0.596938; 0.634189; 0.332015; 0.714059; 0.511094; 0.0298947; 0.0722483]
	# Wc0 = [10.0*ones(11); 6.0*ones(11); .525; 0; 0; .525; 0; -γ^2][:] # (n+m+q)(n+m+q+1)/2 x 1; 28 x 1
	# Wc0 = [10.0*ones(4); -5; 5; 10; 10.0*ones(3); -5; 5; 10; 10.0*ones(2); -5; 5; 10; 10; -5; 5; 10; .525; 0; 10; .525; 10; -1];
	# Wc0 = [10.0*ones(4); 5cos(theta); 5sin(theta); 0.5; 10.0*ones(3); 5cos(theta); 5sin(theta); 0.5; 10.0*ones(2); 5cos(theta); 5sin(theta); 0.5; 10.0; 5cos(theta); 5sin(theta); 0.5; .525; 0; 0; .525; 0; -γ^2];
	
	# Wa10 = 5*cos(theta)*ones(n,) # n x m; 4 x 2
	# Wa20 = 5*sin(theta)*ones(n,)
	Wa10 = 5*ones(n,) # n x m; 4 x 2 rand(n,)
	Wa20 = 5*ones(n,) #rand(n,)#.
	# Wa10 = -5*ones(n,)
	# Wa20 = 5*ones(n,)
	Wd0 = 0.5*ones(n,)
	# u0 = 0.5*ones(m,) # m x 1; 2 x 1
    d0 = 0
	# d_real = 1


	# Trigger parameters
	sizePlant, lixo = size(A)
	eigMinM = minimum(eigvals(M))
	eigMaxR = maximum(eigvals(R))
	eigMinR = minimum(eigvals(R))
	L, Ld = 10, 10
	β = .1 #β1
	L1 = 1 # .9*(β*sqrt(eigMinM/eigMaxR))
	αa_UB = (8*eigMinR-2)/(eigMinR+2)
	αa = .25*αa_UB
	finalTime = zeros(0,)

	x0 = [x1; -1.0; 1.0]
	# x0 = [x1; S.lastVelocity] # lastVelocity
	xf = [x2; 0.0; 0.0]

    # p0 = (x0 - xf)'*M*(x0 - xf)

	Xhat = x0 - xf

	# Basis parameters
	epsilon_basis = 1
	alpha_basis = 1
	time_basis = 20
	basis_on = 1 # 1-on, 2-off
	nu = zeros(Int((n+m+q)*(n+m+q+1)/2), Int((n+m+q)*(n+m+q+1)/2))
	nudelay = zeros(Int((n+m+q)*(n+m+q+1)/2), Int((n+m+q)*(n+m+q+1)/2))
	mu = zeros(n, n)
	xi = zeros(n, n)

	# if basis_on == 1
	r = time_basis - 0
	for i = 1:n
		for j = 1:n
			if i+j < 6
				mu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-1) # alpha_basis*tanh(epsilon_basis*r)^(i+j-1); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-1); %cos(2*pi*(r)/10)^(i+j-1); %
			elseif i+j == 6
				mu[i,j] = alpha_basis*tanh(epsilon_basis*r) # alpha_basis*tanh(epsilon_basis*r); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2); % cos(2*pi*(r)/100); %
			else #if i+j > 6
				mu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-(n+1)) # alpha_basis*tanh(epsilon_basis*r)^(i+j-5); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-5) ; % cos(2*pi*(r)/10)^(i+j-5); %
			end
		end
	end
	# end
	ud = [Wa10'*mu*Xhat; Wa20'*mu*Xhat]
	# println(Xhat)
	# println(mu)
	println(ud)

	t_save = [0,]
	tvec = [0,] # all t in baby, for plotting error gap
	x_save = [[x0; Wc0; Wd0; 0],]
    u_save = [ud,] # u at every T
	uvec = [ud,] # all u
	d_save = [d0,] # d at every T
	dvec = [d0,] # all d

    u1Delay = interp2PWC(getindex.(u_save, 1), -1, 1) # return an interpolation function
	u2Delay = interp2PWC(getindex.(u_save, 2), -1, 1)
	x1Delay = interp2PWC(getindex.(x_save, 1), -1, 1)
    x2Delay = interp2PWC(getindex.(x_save, 2), -1, 1)
    x3Delay = interp2PWC(getindex.(x_save, 3), -1, 1)
	x4Delay = interp2PWC(getindex.(x_save, 4), -1, 1)
	dDelay = interp2PWC(d_save, -1, 1)
	pDelay = interp2PWC(getindex.(x_save, 37), -1, 1)

	xdist = norm(x0[1:2] - xf[1:2])
	error = 0
	localKd = 0.0
	maxIter = 10000
	poseAndKd = Array{Tuple{Array{Float64,2},Float64}}(undef,0)

	trigCondvec = zeros(0)
	normESqvec = zeros(0)
	Wavec = zeros(8,)
	Wcvec = zeros(28,)
	Wdvec = zeros(4,)
	eventTimes = zeros(0,)

	numDataRecord = 0
	e_buff = Array{Float64, 1}(undef, 0)
	halfvec = [5; zeros(6); 5; zeros(5); 5; 0;0;0;0; 5; 0;0;0; .3; 0;0; .3; 0; -γ^2]
	Ω = Array{Float64,2}(undef, Int((n+m+q)*(n+m+q+1)/2), 0)

    # solve ODEs
	for i = 1:N#maxIter
		t = ((i-1)*T, i*T)
		p = [xf, T, A, B, F, M, R, γ, Pt, αa, αc, αd, u1Delay, u2Delay, x1Delay, x2Delay, x3Delay, x4Delay, dDelay, pDelay]#,eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa10, Wa20, sizePlant, β, L, L1, Wavec, eventTimes]#, i,maxIter,start] percent, amplitude,
        sol = solve(ODEProblem(babyETC_robust_relaxedPE_journal_A_trynormalization, x_save[end], t, p), DP5())#, reltol = 1e-4, abstol = 1e-4, dtmax = .05)
		t_save = [t_save; sol.t[2:end]] # vcat(t_save, sol.t) save time
		x_save = [x_save; sol.u[2:end]] # vcat(x_save, sol.u) save state
		u_save = [u_save; [uvec[end]]]
		d_save = [d_save; dvec[end]]

		Wcvec = [Wcvec sol.u[end][5:32]]
		Wdvec = [Wdvec sol.u[end][33:36]]
		# uvec1 = p[end-2]
		# uvec2 = p[end-1]
		# trigCondvec = p[32]
		# normESqvec = p[31]
		# println("$trigCondvec")
        u1Delay = interp2PWC(getindex.(u_save, 1), -1, i*T+.01) # interpolate control input
		u2Delay = interp2PWC(getindex.(u_save, 2), -1, i*T+.01)
		x1Delay = interp2PWC(getindex.(x_save, 1), -1, i*T+.01)
		x2Delay = interp2PWC(getindex.(x_save, 2), -1, i*T+.01)
		x3Delay = interp2PWC(getindex.(x_save, 3), -1, i*T+.01)
		x4Delay = interp2PWC(getindex.(x_save, 4), -1, i*T+.01)
		dDelay = interp2PWC(d_save, -1, i*T+.01)
		pDelay = interp2PWC(getindex.(x_save, 37), -1, i*T+.01)
		# trigCondDelay = interp2PWC(trigCondvec, -1, i*T+.01)

        # relaxed PE: record data
        if rank(Ω) < n+m+q
			∂nu∂t = zeros(Int((n+m+q)*(n+m+q+1)/2), Int((n+m+q)*(n+m+q+1)/2))
			if basis_on == 1
				r = time_basis - i*T
				for i = 1:Int((n+m+q)*(n+m+q+1)/2)
					for j = 1:Int((n+m+q)*(n+m+q+1)/2)
						if i+j < 23
							nu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-1)
							∂nu∂t[i,j] = alpha_basis*(i+j-1)*(1-tanh(epsilon_basis*r)^2)*(-epsilon_basis)*tanh(epsilon_basis*r)^(i+j-2) # alpha_basis*tanh(epsilon_basis*r)^(i+j-1); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-1); % cos(2*pi*(r)/10)^(i+j-1); %
						elseif i+j == 23
							nu[i,j] = alpha_basis*tanh(epsilon_basis*r)
							∂nu∂t[i,j] = alpha_basis*(1-tanh(epsilon_basis*r)^2)*(-epsilon_basis) # alpha_basis*tanh(epsilon_basis*r); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2); % cos(2*pi*(r)/10); %
						else #if i+j > 23
							nu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-(((n+m+q)*(n+m+q+1)/2)+1))
							∂nu∂t[i,j] = alpha_basis*(i+j-(((n+m+q)*(n+m+q+1)/2)+1))*(1-tanh(epsilon_basis*r)^2)*(-epsilon_basis)*tanh(epsilon_basis*r)^(i+j-(((n+m+q)*(n+m+q+1)/2)+1)-1) # alpha_basis*tanh(epsilon_basis*r)^(i+j-22); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-22) ; % cos(2*pi*(r)/10)^(i+j-22); %
						end
					end
				end
			end
			x_dot = (x_save[end][1:4] - x_save[end-1][1:4])/T
			u_dot = (u_save[end] - u_save[end-1])/T
			d_dot = (d_save[end] - d_save[end-1])/T
			U = [x_save[end][1:n] - xf; u_save[end]; d_save[end]]
			ud = u_save[end]
			d = d_save[end]
			∂UkU∂t = [2*U[1]*x_dot[1]; x_dot[1]*U[2] + U[1]*x_dot[2]; x_dot[1]*U[3] + U[1]*x_dot[3]; x_dot[1]*U[4] + U[1]*x_dot[4]; x_dot[1]*ud[1] + U[1]*u_dot[1]; x_dot[1]*ud[2] + U[1]*u_dot[2]; x_dot[1]*d + U[1]*d_dot;
			 		  2*U[2]*x_dot[2]; x_dot[2]*U[3] + U[2]*x_dot[3]; x_dot[2]*U[4] + U[2]*x_dot[4]; x_dot[2]*ud[1] + U[2]*u_dot[1]; x_dot[2]*ud[2] + U[2]*u_dot[2]; x_dot[2]*d + U[2]*d_dot;
				      2*U[3]*x_dot[3]; x_dot[3]*U[4] + U[3]*x_dot[4]; x_dot[3]*ud[1] + U[3]*u_dot[1]; x_dot[3]*ud[2] + U[3]*u_dot[2]; x_dot[3]*d + U[3]*d_dot;
					  2*U[4]*x_dot[4]; x_dot[4]*ud[1] + U[4]*u_dot[1]; x_dot[4]*ud[2] + U[4]*u_dot[2]; x_dot[4]*d + U[4]*d_dot;
					  2*ud[1]*u_dot[1]; u_dot[1]*ud[2] + ud[1]*u_dot[2]; u_dot[1]*d + ud[1]*d_dot;
					  2*ud[2]*u_dot[2]; u_dot[2]*d + ud[2]*d_dot;
					  2*d*d_dot]

			UkU = vcat(U[1]^2, U[1]*U[2], U[1]*U[3], U[1]*U[4], U[1]*ud[1], U[1]*ud[2], U[1]*d,
			U[2]^2, U[2]*U[3], U[2]*U[4], U[2]*ud[1], U[2]*ud[2], U[2]*d,
			U[3]^2, U[3]*U[4], U[3]*ud[1], U[3]*ud[2], U[3]*d,
			U[4]^2, U[4]*ud[1], U[4]*ud[2], U[4]*d,
			ud[1]^2, ud[1]*ud[2], ud[1]*d,
			ud[2]^2, ud[2]*d,
			d^2)
			∂UkU∂̄x	= hcat([2*U[1]; U[2:7]; zeros(21,)],  # 28 x 4
					[0; U[2]; zeros(5,); 2U[2]; U[3:7]; zeros(15,1)],
			        [0; 0; U[1]; zeros(5,); U[2]; zeros(4,); 2U[3]; U[4:7]; zeros(10,1)],
					[0; 0; 0; U[1]; zeros(5,); U[2]; zeros(4,); U[3]; zeros(3,1); 2U[4]; U[5:7]; zeros(6,)])					
			ω = ∂nu∂t*UkU + nu*∂UkU∂t + nu*∂UkU∂̄x*x_dot
			e_buff = [e_buff; (.5*halfvec'*UkU + Wc'*ω - .5*([Wa10 Wa20]'*mu*(U[1:n] - Xhat))'*R*([Wa10 Wa20]'*mu*(U[1:n] - Xhat)))];
			Ω = [Ω ω]
            numDataRecord += 1
        end

		# S.kino_dist is always the max kd till now
		localKd = distPoint2Line(x_save[end][1:2], x1, x2)
		poseAndKd = [poseAndKd; ([x_save[end][1] x_save[end][2] x_save[end][3] x_save[end][4]], copy(localKd))]
		if norm(x_save[end][1:2] - xf[1:2]) < error*xdist break end
	end
	# S.lastVelocity = x_save[end][3:4]
	# lastVelocity = x_save[end][3:4]
	(poseAndKd, normESqvec, trigCondvec)#,t_save)#eventTimes)
	# println(length(normESqvec))
	# println(length(trigCondvec))
	# println(length(poseAndKd))
	# plot(1:length(normESqvec), normESqvec)
	# plot!(1:length(normESqvec), trigCondvec)
	println(length(eventTimes))
	# println(eventTimes)
	println((eventTimes[2:end] - eventTimes[1:end-1]))

	println(getindex.(u_save,1))
	println(getindex.(u_save,2))

	# println(getindex.(x_save, 1))
	# println(getindex.(x_save, 2))
	# println(getindex.(x_save, 3))
	# println(getindex.(x_save, 4))
	println(x_save[end][1])
	println(x_save[end][2])
	println(x_save[end][3])
	println(x_save[end][4])

	# saveData(trigCondvec, "temp/Trig.txt")
	# saveData(normESqvec, "temp/Esq.txt")
	# saveData(getindex.(getindex.(poseAndKd, 1), 1), "temp/x1.txt")
	# saveData(getindex.(getindex.(poseAndKd, 1), 2), "temp/x2.txt")
	# saveData(getindex.(getindex.(poseAndKd, 1), 3), "temp/x3.txt")
	# saveData(getindex.(getindex.(poseAndKd, 1), 4), "temp/x4.txt")
	# saveData(getindex.(x_save, 1), "temp/x1.txt")
	# saveData(getindex.(x_save, 2), "temp/x2.txt")
	# saveData(getindex.(x_save, 3), "temp/x3.txt")
	# saveData(getindex.(x_save, 4), "temp/x4.txt")
	# saveData(tvec, "temp/tvec.txt")
	# saveData(t_save, "temp/tsave.txt")
	# println(t_save)
	# println(tvec)

	# println(trigCondvec)
	# println(normESqvec)

	# plot(tvec[2:end], normESqvec, lab = "squared norm error", lw = 2, grid = ":no")
	# plot!(tvec[2:end], trigCondvec, lab = "threshold", lw = 2)
	# plot(tvec[2:end], [normESqvec trigCondvec], lab = ["squared norm error" "threshold"], lc = [:pink :red], lw = 2, grid = :no, legend = :topright)
	# plot(range(0, stop = 20, length = length(normESqvec)), [normESqvec trigCondvec], lab = ["squared norm error" "threshold"], lc = [:pink :red], lw = 2, grid = :no, legend = :topright)
	plot(range(0, stop = 15, length = length(u_save)), [getindex.(u_save,1) getindex.(u_save,2)], lab = ["u1" "u2"], lc = [:pink :red], lw = 2, grid = :no, legend = :topright)

	# plot(t_save, [getindex.(x_save, 1) getindex.(x_save, 2) getindex.(x_save, 3) getindex.(x_save, 4)], lab = ["x1" "x2" "x3" "x4"], lw = 2, grid = :no, legend = :topright) # lc = [:pink :red],[getindex.(getindex.(poseAndKd, 1), 1) getindex.(getindex.(poseAndKd, 1), 2) getindex.(getindex.(poseAndKd, 1), 3) getindex.(getindex.(poseAndKd, 1), 4)]

	# plot(getindex.(getindex.(poseAndKd, 1), 1), getindex.(getindex.(poseAndKd, 1), 2))
	# maximum(getindex.(poseAndKd, 2))
	# println(distPoint2Line(poseAndKd[end][1][1:2], x1, x2))
	# plot(t_save, getindex.(x_save, 37)) # check cost
	# plot(T:T:T*length(vecInt), vecInt) # check integral
	# plot(eventTimes, Wavec[:,2:end]')
	# plot(1:size(Wcvec)[2], Wcvec',label=0)
	# plot(1:size(Wdvec)[2], Wdvec',label=0)
end


function babyETC_robust_relaxedPE_journal_A_nobasis(dotx, x, p, t)

	global ud, eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa10, Wa20, sizePlant, β, L, L1, Ld, Wavec, Wcvec, Wdvec, numDataRecord, d_real, UkU, Wc, e_buff, Ω, mu, nu, nudelay, xi, basis_on, alpha_basis, epsilon_basis, time_basis, uvec, dvec, tvec

	# global xf, uf, Tf, T, A, B, M, R, Pt, per, amp, αa, αc, Quu, Qxu, Wcfinal, uDelay, x1Delay, x2Delay, x3Delay, x_save, t_save, uvec
	xf, T, A, B, F, M, R, γ, Pt, αa, αc, αd, u1Delay, u2Delay, x1Delay, x2Delay, x3Delay, x4Delay, dDelay, pDelay = p#,eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa10, Wa20, sizePlant, β, L, L1, Wavec, eventTimes=p#,i,maxIter,start = p # uu, vv useless here percent, amplitude,
    n, m, q = 4, 2, 1
	Wc = x[(n+1) : Int(n+(n+m+q)*(n+m+q+1)/2)] # 5-32
	Wd = x[Int(n+(n+m+q)*(n+m+q+1)/2+1) : Int(n+(n+m+q)*(n+m+q+1)/2+n)] # 33-36
	P = x[end] # 37
	d_real = 0#sin(t)#0#0.1

	# Basis
	# if basis_on == 1
	# 	r = time_basis - t
	# 	rdelay = time_basis - (t - T)
	# 	for i = 1:Int((n+m+q)*(n+m+q+1)/2)
	# 		for j = 1:Int((n+m+q)*(n+m+q+1)/2)
	# 			if i+j < 23
	# 				nu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-1) # alpha_basis*tanh(epsilon_basis*r)^(i+j-1); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-1); % cos(2*pi*(r)/10)^(i+j-1); %
	# 				nudelay[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-1)
	# 			elseif i+j == 23
	# 				nu[i,j] = alpha_basis*tanh(epsilon_basis*r) # alpha_basis*tanh(epsilon_basis*r); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2); % cos(2*pi*(r)/10); %
	# 				nudelay[i,j] = alpha_basis*tanh(epsilon_basis*r)
	# 			else #if i+j > 23
	# 				nu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-(((n+m+q)*(n+m+q+1)/2)+1)) # alpha_basis*tanh(epsilon_basis*r)^(i+j-22); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-22) ; % cos(2*pi*(r)/10)^(i+j-22); %
	# 				nudelay[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-(((n+m+q)*(n+m+q+1)/2)+1))
	# 			end
	# 		end
	# 	end
	# 	for i = 1:n
	# 		for j = 1:n
	# 			if i+j < 6
	# 				mu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-1) # alpha_basis*tanh(epsilon_basis*r)^(i+j-1); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-1); %cos(2*pi*(r)/10)^(i+j-1); %
	# 				xi[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-1)
	# 			elseif i+j == 6
	# 				mu[i,j] = alpha_basis*tanh(epsilon_basis*r) # alpha_basis*tanh(epsilon_basis*r); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2); % cos(2*pi*(r)/100); %
	# 				xi[i,j] = alpha_basis*tanh(epsilon_basis*r)
	# 			else #if i+j > 6
	# 				mu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-(n+1)) # alpha_basis*tanh(epsilon_basis*r)^(i+j-5); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-5) ; % cos(2*pi*(r)/10)^(i+j-5); %
	# 				xi[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-(n+1))
	# 			end
	# 		end
	# 	end
	# end

	# Update disturbance
	d = Wd'*(x[1:4] - xf)

	# Delays
    uddelay = zeros(m,)
	uddelay[1] = u1Delay(t - T)
	uddelay[2] = u2Delay(t - T)
    xdelay = zeros(n,)
	xdelay[1] = x1Delay(t - T)
    xdelay[2] = x2Delay(t - T)
    xdelay[3] = x3Delay(t - T)
	xdelay[4] = x4Delay(t - T)
	ddelay = dDelay(t - T)
	pdelay = pDelay(t - T)

	# Kronecker products
	U = vcat(x[1:n] - xf, ud[1], ud[2], d) # augmented state
    UkU = vcat(U[1]^2, U[1]*U[2], U[1]*U[3], U[1]*U[4], U[1]*ud[1], U[1]*ud[2], U[1]*d,
               U[2]^2, U[2]*U[3], U[2]*U[4], U[2]*ud[1], U[2]*ud[2], U[2]*d,
               U[3]^2, U[3]*U[4], U[3]*ud[1], U[3]*ud[2], U[3]*d,
			   U[4]^2, U[4]*ud[1], U[4]*ud[2], U[4]*d,
               ud[1]^2, ud[1]*ud[2], ud[1]*d,
			   ud[2]^2, ud[2]*d,
			   d^2)
    UkUdelay = vcat(xdelay[1]^2, xdelay[1]*xdelay[2], xdelay[1]*xdelay[3], xdelay[1]*xdelay[4], xdelay[1]*uddelay[1], xdelay[1]*uddelay[2], xdelay[1]*ddelay,
                    xdelay[2]^2, xdelay[2]*xdelay[3], xdelay[2]*xdelay[4], xdelay[2]*uddelay[1], xdelay[2]*uddelay[2], xdelay[2]*ddelay,
                    xdelay[3]^2, xdelay[3]*xdelay[4], xdelay[3]*uddelay[1], xdelay[3]*uddelay[2], xdelay[3]*ddelay,
					xdelay[4]^2, xdelay[4]*uddelay[1], xdelay[4]*uddelay[2], xdelay[4]*ddelay,
                    uddelay[1]^2, uddelay[1]*uddelay[2], uddelay[1]*ddelay,
					uddelay[2]^2, uddelay[2]*ddelay,
					ddelay^2)
    Quu = [Wc[23] Wc[24]; Wc[24] Wc[26]] # m x m 24-27
    Quu_inv = inv(Quu)
    Qux = [Wc[5] Wc[11] Wc[16] Wc[20]; Wc[6] Wc[12] Wc[17] Wc[21]] #reshape(Wc[Int((n+m+q)*(n+m+q+1)/2-m^2-q^2-m*n+1):Int((n+m+q)*(n+m+q+1)/2-m^2-q^2)], (n, m)) # n x m 16-23
    Qxu = Qux' # m x n
    Qdx = [Wc[7] Wc[13] Wc[18] Wc[22]]#reshape(Wc[Int((n+m+q)*(n+m+q+1)/2-m^2-q^2-m*n-n*q+1):Int((n+m+q)*(n+m+q+1)/2-m^2-q^2-m*n)], (n, q)) # n x q 12-15
    Qxd = Qdx'
    # Qdd = Wc[end] # q x q
	Qdd_inv = 1.0/Wc[end]

	# Triggering condition
	normXSq = abs(U[1])^2 + abs(U[2])^2 + U[3]^2 + U[4]^2
	# trigCond = (((1 - β^2)*eigMinM - γ^2*Ld^2)/(2*L^2*eigMaxR))*normXSq + ((1-β^2)*eigMinR/(2*L^2*eigMaxR))*(abs(ud[1])^2+abs(ud[2])^2)
	trigCond = ((1 - β^2)*eigMinM*normXSq + (1-β^2)*eigMinR*(abs(ud[1])^2+abs(ud[2])^2) - γ^2*d_real^2)/(4*eigMaxR*(L^2+L1^2)) # when using d, threshold always < 0

	# Integral reinforcement dynamics
	# dP = 0.5*((x[1:n] - xf)'*M*(x[1:n] - xf) - (xdelay - xf)'*M*(xdelay - xf) + ud'*R*ud - uddelay'*R*uddelay - γ^2*d^2 + γ^2*ddelay^2)
	dP = 0.5*(U[1:n]'*M*U[1:n] + ud'*R*ud - γ^2*d_real^2) # d

	# approximation errors
	σ = UkU - UkUdelay 
	σ_f	= UkU
    ec = P - pdelay + Wc'*σ #P + Wc'*σ
	ecfinal = 0.5*U[1:n]'*Pt*U[1:n] - Wc'*σ_f
	ed = Wd'*U[1:n] + (Qdd_inv*Qdx*U[1:n])[1] # q x 1

	# gap
	e = Xhat - U[1:n]
	normESq = norm(e)^2

	# Critic dynamics
	relaxedPETerm = zeros(Int((n+m+q)*(n+m+q+1)/2),)
	for rrr in 1:numDataRecord
		relaxedPETerm = relaxedPETerm + e_buff[rrr]*Ω[:,rrr]./((1 + Ω[:,rrr]'*Ω[:,rrr]).^2)
	end
	dWc = -αc*((σ./(σ'*σ+1).^2)*ec + (σ_f./(σ_f'*σ_f+1).^2)*ecfinal + relaxedPETerm)

	# Triggering updates
	if normESq >= 0.95*trigCond # flows
	    Xhat = U[1:n]
		ea1 = Wa10'*Xhat + [Quu_inv[1,:]'*Qux[:,1]; Quu_inv[1,:]'*Qux[:,2]; Quu_inv[1,:]'*Qux[:,3]; Quu_inv[1,:]'*Qux[:,4]]'*Xhat
		ea2 = Wa20'*Xhat + [Quu_inv[2,:]'*Qux[:,1]; Quu_inv[2,:]'*Qux[:,2]; Quu_inv[2,:]'*Qux[:,3]; Quu_inv[2,:]'*Qux[:,4]]'*Xhat
		Wa10 = Wa10 - αa*Xhat*ea1'/(1 + (Xhat)'*(Xhat))
		Wa20 = Wa20 - αa*Xhat*ea2'/(1 + (Xhat)'*(Xhat))
		Wavec = [Wavec [Wa10; Wa20]]
		eventTimes = [eventTimes; t]
		# Update control
		ud = [Wa10'*Xhat; Wa20'*Xhat]
		# ud = zeros(m,)
		# ud[1] = Wa10'*mu*Xhat
		# ud[2] = Wa20'*mu*Xhat
	end

	# Disturbance dynamics
	dWd = -αd*U[1:n]*ed'/(1 + (U[1:n])'*(U[1:n]))

	# System dynamics
	dx = A*U[1:n] + B*ud + F*d_real

	# save vecs
	# p[end] = unew
	# p[end-1:end] .= unew
	uvec = [uvec; [ud]]
	normESqvec = [normESqvec; normESq]
	# p[31] = [p[31]; normESq]
	trigCondvec = [trigCondvec; trigCond]
	tvec = [tvec; t]

	# p[32] = [p[32]; trigCond]
	# println("$trigCondvec")
	dvec = [dvec; d]
	# dd = d

	dotx .= vcat(dx, dWc, dWd, dP)
end

function sim_journal_A_ETC_robust_relaxedPE_Local_nobasis(x1, x2)#, S) # lastVelocity # x1 -> x2, Journal 20

	global ud, eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa10, Wa20, sizePlant, β, L, L1, Ld, Wavec, Wcvec, Wdvec, numDataRecord, d_real, UkU, Wc, e_buff, Ω, mu, nu, nudelay, xi, basis_on, alpha_basis, epsilon_basis, time_basis, uvec, dvec, tvec

	# System Dynamics for feedback
    n, m, q = 4, 2, 1
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
    F = [0; 0; 1; 1] # n x q
    M = 5.0*Matrix(I, n, n) # n x n # 1, .1, 10
    R = 1.0*Matrix(I, m, m) # m x m # 1, 10, .1
	γ = .1
	Pt = 0.3*Matrix(I, n, n)

	# check controllability & observability
    Co = ctrb(A, B)
    unco = size(A, 1) - rank(Co)
    Ob = obsv(sqrt(M), A)
    unob = size(sqrt(M), 1) - rank(Ob)

	if unco + unob > 0	error("system uncontrollable and/or unobservable")	end

    # ODE parameters
	Tf, T, N = 15, 0.05, 300 # finite horizon
	αc, αd = 90, 3 # αa = 1.2
	
	# goal orientation
	# theta = atan((x2-x1)[2]/(x2-x1)[1])

	# Wc0 = 10*rand(28)[:]
	Wc0 = [5.0*ones(Int((n+m+q)*(n+m+q+1)/2 - m*(m+1)/2 - m*q - q^2), 1); 1; 0; 0; 1; 0; -γ^2][:] # (n+m+q)(n+m+q+1)/2 x 1; 28 x 1
	# Wc0 = 20*[0.328555; 0.756487; 0.273798; 0.309953; 0.438633; 0.307851; 0.335527; 0.0851001; 0.608546; 0.282202; 0.399089; 0.281325; 0.528826; 0.234641; 0.584131; 0.313394; 0.897723; 0.675292; 0.73439; 0.445563; 0.592254; 0.596938; 0.634189; 0.332015; 0.714059; 0.511094; 0.0298947; 0.0722483]
	# Wc0 = [10.0*ones(11); 6.0*ones(11); .525; 0; 0; .525; 0; -γ^2][:] # (n+m+q)(n+m+q+1)/2 x 1; 28 x 1
	# Wc0 = [10.0*ones(4); -5; 5; 10; 10.0*ones(3); -5; 5; 10; 10.0*ones(2); -5; 5; 10; 10; -5; 5; 10; .525; 0; 10; .525; 10; -1];
	# Wc0 = [10.0*ones(4); 5cos(theta); 5sin(theta); 0.5; 10.0*ones(3); 5cos(theta); 5sin(theta); 0.5; 10.0*ones(2); 5cos(theta); 5sin(theta); 0.5; 10.0; 5cos(theta); 5sin(theta); 0.5; .525; 0; 0; .525; 0; -γ^2];
	
	# Wa10 = 5*cos(theta)*ones(n,) # n x m; 4 x 2
	# Wa20 = 5*sin(theta)*ones(n,)
	Wa10 = -5*ones(n,) # n x m; 4 x 2 rand(n,)
	Wa20 = 5*ones(n,) #rand(n,)#.
	# Wa10 = -5*ones(n,)
	# Wa20 = 5*ones(n,)
	Wd0 = 0.5*ones(n,)
	# u0 = 0.5*ones(m,) # m x 1; 2 x 1
    d0 = 0
	# d_real = 1


	# Trigger parameters
	sizePlant, lixo = size(A)
	eigMinM = minimum(eigvals(M))
	eigMaxR = maximum(eigvals(R))
	eigMinR = minimum(eigvals(R))
	L, Ld = 10, 10
    β = .1 #β1
	L1 = 1 # .9*(β*sqrt(eigMinM/eigMaxR))
	αa_UB = (8*eigMinR-2)/(eigMinR+2)
	αa = .25*αa_UB
	finalTime = zeros(0,)

	x0 = [x1; 1.0; -1.0]
	# x0 = [x1; S.lastVelocity] # lastVelocity
	xf = [x2; 0.0; 0.0]

    # p0 = (x0 - xf)'*M*(x0 - xf)

	Xhat = x0 - xf

	# Basis parameters
	# epsilon_basis = 1
	# alpha_basis = 5
	# time_basis = 25
	# basis_on = 1 # 1-on, 2-off
	# nu = zeros(Int((n+m+q)*(n+m+q+1)/2), Int((n+m+q)*(n+m+q+1)/2))
	# nudelay = zeros(Int((n+m+q)*(n+m+q+1)/2), Int((n+m+q)*(n+m+q+1)/2))
	# mu = zeros(n, n)
	# xi = zeros(n, n)

	# if basis_on == 1
	# r = time_basis - 0
	# for i = 1:n
	# 	for j = 1:n
	# 		if i+j < 6
	# 			mu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-1) # alpha_basis*tanh(epsilon_basis*r)^(i+j-1); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-1); %cos(2*pi*(r)/10)^(i+j-1); %
	# 		elseif i+j == 6
	# 			mu[i,j] = alpha_basis*tanh(epsilon_basis*r) # alpha_basis*tanh(epsilon_basis*r); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2); % cos(2*pi*(r)/100); %
	# 		else #if i+j > 6
	# 			mu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-(n+1)) # alpha_basis*tanh(epsilon_basis*r)^(i+j-5); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-5) ; % cos(2*pi*(r)/10)^(i+j-5); %
	# 		end
	# 	end
	# end
	# end
	ud = [Wa10'*Xhat; Wa20'*Xhat]
	# println(Xhat)
	# println(mu)
	println(ud)

	t_save = [0,]
	tvec = [0,] # all t in baby, for plotting error gap
	x_save = [[x0; Wc0; Wd0; 0],]
    u_save = [ud,] # u at every T
	uvec = [ud,] # all u
	d_save = [d0,] # d at every T
	dvec = [d0,] # all d

    u1Delay = interp2PWC(getindex.(u_save, 1), -1, 1) # return an interpolation function
	u2Delay = interp2PWC(getindex.(u_save, 2), -1, 1)
	x1Delay = interp2PWC(getindex.(x_save, 1), -1, 1)
    x2Delay = interp2PWC(getindex.(x_save, 2), -1, 1)
    x3Delay = interp2PWC(getindex.(x_save, 3), -1, 1)
	x4Delay = interp2PWC(getindex.(x_save, 4), -1, 1)
	dDelay = interp2PWC(d_save, -1, 1)
	pDelay = interp2PWC(getindex.(x_save, 37), -1, 1)

	xdist = norm(x0[1:2] - xf[1:2])
	error = 0
	localKd = 0.0
	maxIter = 10000
	poseAndKd = Array{Tuple{Array{Float64,2},Float64}}(undef,0)

	trigCondvec = zeros(0)
	normESqvec = zeros(0)
	Wavec = zeros(8,)
	Wcvec = zeros(28,)
	Wdvec = zeros(4,)
	eventTimes = zeros(0,)

	numDataRecord = 0
	e_buff = Array{Float64, 1}(undef, 0)
	halfvec = [5; zeros(6); 5; zeros(5); 5; 0;0;0;0; 5; 0;0;0; 1; 0;0; 1; 0; -γ^2]
	Ω = Array{Float64,2}(undef, Int((n+m+q)*(n+m+q+1)/2), 0)

    # solve ODEs
	for i = 1:N#maxIter
		t = ((i-1)*T, i*T)
		p = [xf, T, A, B, F, M, R, γ, Pt, αa, αc, αd, u1Delay, u2Delay, x1Delay, x2Delay, x3Delay, x4Delay, dDelay, pDelay]#,eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa10, Wa20, sizePlant, β, L, L1, Wavec, eventTimes]#, i,maxIter,start] percent, amplitude,
        sol = solve(ODEProblem(babyETC_robust_relaxedPE_journal_A_nobasis, x_save[end], t, p), DP5())#, reltol = 1e-4, abstol = 1e-4, dtmax = .05)
		t_save = [t_save; sol.t[2:end]] # vcat(t_save, sol.t) save time
		x_save = [x_save; sol.u[2:end]] # vcat(x_save, sol.u) save state
		u_save = [u_save; [uvec[end]]]
		d_save = [d_save; dvec[end]]

		Wcvec = [Wcvec sol.u[end][5:32]]
		Wdvec = [Wdvec sol.u[end][33:36]]
		# uvec1 = p[end-2]
		# uvec2 = p[end-1]
		# trigCondvec = p[32]
		# normESqvec = p[31]
		# println("$trigCondvec")
        u1Delay = interp2PWC(getindex.(u_save, 1), -1, i*T+.01) # interpolate control input
		u2Delay = interp2PWC(getindex.(u_save, 2), -1, i*T+.01)
		x1Delay = interp2PWC(getindex.(x_save, 1), -1, i*T+.01)
		x2Delay = interp2PWC(getindex.(x_save, 2), -1, i*T+.01)
		x3Delay = interp2PWC(getindex.(x_save, 3), -1, i*T+.01)
		x4Delay = interp2PWC(getindex.(x_save, 4), -1, i*T+.01)
		dDelay = interp2PWC(d_save, -1, i*T+.01)
		pDelay = interp2PWC(getindex.(x_save, 37), -1, i*T+.01)
		# trigCondDelay = interp2PWC(trigCondvec, -1, i*T+.01)

        # relaxed PE: record data
        if rank(Ω) < n+m+q
			# ∂nu∂t = zeros(Int((n+m+q)*(n+m+q+1)/2), Int((n+m+q)*(n+m+q+1)/2))
			# if basis_on == 1
			# 	r = time_basis - i*T
			# 	for i = 1:Int((n+m+q)*(n+m+q+1)/2)
			# 		for j = 1:Int((n+m+q)*(n+m+q+1)/2)
			# 			if i+j < 23
			# 				nu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-1)
			# 				∂nu∂t[i,j] = alpha_basis*(i+j-1)*(1-tanh(epsilon_basis*r)^2)*(-epsilon_basis)*tanh(epsilon_basis*r)^(i+j-2) # alpha_basis*tanh(epsilon_basis*r)^(i+j-1); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-1); % cos(2*pi*(r)/10)^(i+j-1); %
			# 			elseif i+j == 23
			# 				nu[i,j] = alpha_basis*tanh(epsilon_basis*r)
			# 				∂nu∂t[i,j] = alpha_basis*(1-tanh(epsilon_basis*r)^2)*(-epsilon_basis) # alpha_basis*tanh(epsilon_basis*r); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2); % cos(2*pi*(r)/10); %
			# 			else #if i+j > 23
			# 				nu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-(((n+m+q)*(n+m+q+1)/2)+1))
			# 				∂nu∂t[i,j] = alpha_basis*(i+j-(((n+m+q)*(n+m+q+1)/2)+1))*(1-tanh(epsilon_basis*r)^2)*(-epsilon_basis)*tanh(epsilon_basis*r)^(i+j-(((n+m+q)*(n+m+q+1)/2)+1)-1) # alpha_basis*tanh(epsilon_basis*r)^(i+j-22); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-22) ; % cos(2*pi*(r)/10)^(i+j-22); %
			# 			end
			# 		end
			# 	end
			# end
			x_dot = (x_save[end][1:4] - x_save[end-1][1:4])/T
			u_dot = (u_save[end] - u_save[end-1])/T
			d_dot = (d_save[end] - d_save[end-1])/T
			U = [x_save[end][1:n] - xf; u_save[end]; d_save[end]]
			ud = u_save[end]
			d = d_save[end]
			∂UkU∂t = [2*U[1]*x_dot[1]; x_dot[1]*U[2] + U[1]*x_dot[2]; x_dot[1]*U[3] + U[1]*x_dot[3]; x_dot[1]*U[4] + U[1]*x_dot[4]; x_dot[1]*ud[1] + U[1]*u_dot[1]; x_dot[1]*ud[2] + U[1]*u_dot[2]; x_dot[1]*d + U[1]*d_dot;
			 		  2*U[2]*x_dot[2]; x_dot[2]*U[3] + U[2]*x_dot[3]; x_dot[2]*U[4] + U[2]*x_dot[4]; x_dot[2]*ud[1] + U[2]*u_dot[1]; x_dot[2]*ud[2] + U[2]*u_dot[2]; x_dot[2]*d + U[2]*d_dot;
				      2*U[3]*x_dot[3]; x_dot[3]*U[4] + U[3]*x_dot[4]; x_dot[3]*ud[1] + U[3]*u_dot[1]; x_dot[3]*ud[2] + U[3]*u_dot[2]; x_dot[3]*d + U[3]*d_dot;
					  2*U[4]*x_dot[4]; x_dot[4]*ud[1] + U[4]*u_dot[1]; x_dot[4]*ud[2] + U[4]*u_dot[2]; x_dot[4]*d + U[4]*d_dot;
					  2*ud[1]*u_dot[1]; u_dot[1]*ud[2] + ud[1]*u_dot[2]; u_dot[1]*d + ud[1]*d_dot;
					  2*ud[2]*u_dot[2]; u_dot[2]*d + ud[2]*d_dot;
					  2*d*d_dot]

			UkU = vcat(U[1]^2, U[1]*U[2], U[1]*U[3], U[1]*U[4], U[1]*ud[1], U[1]*ud[2], U[1]*d,
			U[2]^2, U[2]*U[3], U[2]*U[4], U[2]*ud[1], U[2]*ud[2], U[2]*d,
			U[3]^2, U[3]*U[4], U[3]*ud[1], U[3]*ud[2], U[3]*d,
			U[4]^2, U[4]*ud[1], U[4]*ud[2], U[4]*d,
			ud[1]^2, ud[1]*ud[2], ud[1]*d,
			ud[2]^2, ud[2]*d,
			d^2)
			∂UkU∂̄x	= hcat([2*U[1]; U[2:7]; zeros(21,)],  # 28 x 4
					[0; U[2]; zeros(5,); 2U[2]; U[3:7]; zeros(15,1)],
			        [0; 0; U[1]; zeros(5,); U[2]; zeros(4,); 2U[3]; U[4:7]; zeros(10,1)],
					[0; 0; 0; U[1]; zeros(5,); U[2]; zeros(4,); U[3]; zeros(3,1); 2U[4]; U[5:7]; zeros(6,)])					
			ω = ∂UkU∂t + ∂UkU∂̄x*x_dot
			e_buff = [e_buff; (.5*halfvec'*UkU + Wc'*ω - .5*([Wa10 Wa20]'*(U[1:n] - Xhat))'*R*([Wa10 Wa20]'*(U[1:n] - Xhat)))];
			Ω = [Ω ω]
            numDataRecord += 1
        end

		# S.kino_dist is always the max kd till now
		localKd = distPoint2Line(x_save[end][1:2], x1, x2)
		poseAndKd = [poseAndKd; ([x_save[end][1] x_save[end][2] x_save[end][3] x_save[end][4]], copy(localKd))]
		if norm(x_save[end][1:2] - xf[1:2]) < error*xdist break end
	end
	# S.lastVelocity = x_save[end][3:4]
	# lastVelocity = x_save[end][3:4]
	(poseAndKd, normESqvec, trigCondvec)#,t_save)#eventTimes)
	# println(length(normESqvec))
	# println(length(trigCondvec))
	# println(length(poseAndKd))
	# plot(1:length(normESqvec), normESqvec)
	# plot!(1:length(normESqvec), trigCondvec)
	println(length(eventTimes))
	# println(eventTimes)
	# println((eventTimes[2:end] - eventTimes[1:end-1]))
	scatter(eventTimes[2:end], (eventTimes[2:end] - eventTimes[1:end-1]), lab = false, m = :cross, lw = 2, grid = :no, legend = :topright)

	# println(getindex.(u_save,1))
	# println(getindex.(u_save,2))

	# println(getindex.(x_save, 1))
	# println(getindex.(x_save, 2))
	# println(getindex.(x_save, 3))
    # println(getindex.(x_save, 4))
    
	# println(x_save[end][1])
	# println(x_save[end][2])
	# println(x_save[end][3])
	# println(x_save[end][4])

	# println(t_save)
	# println(tvec)

	# println(trigCondvec)
	# println(normESqvec)

	# plot(tvec[2:end], normESqvec, lab = "squared norm error", lw = 2, grid = ":no")
	# plot!(tvec[2:end], trigCondvec, lab = "threshold", lw = 2)
	# plot(tvec[2:end], [normESqvec trigCondvec], lab = ["squared norm error" "threshold"], lc = [:pink :red], lw = 2, grid = :no, legend = :topright)
	# plot(range(0, stop = 15, length = length(normESqvec)), [normESqvec trigCondvec], lab = ["squared norm error" "threshold"], lc = [:pink :red], lw = 2, grid = :no, legend = :topright)
	# plot(range(0, stop = 15, length = length(u_save)), [getindex.(u_save,1) getindex.(u_save,2)], lab = ["u1" "u2"], lc = [:pink :red], lw = 2, grid = :no, legend = :topright)

	# plot(t_save, [getindex.(x_save, 1) getindex.(x_save, 2) getindex.(x_save, 3) getindex.(x_save, 4)], lab = ["x1" "x2" "x3" "x4"], lw = 2, grid = :no, legend = :topright) # lc = [:pink :red],[getindex.(getindex.(poseAndKd, 1), 1) getindex.(getindex.(poseAndKd, 1), 2) getindex.(getindex.(poseAndKd, 1), 3) getindex.(getindex.(poseAndKd, 1), 4)]

	# plot(getindex.(getindex.(poseAndKd, 1), 1), getindex.(getindex.(poseAndKd, 1), 2))
	# maximum(getindex.(poseAndKd, 2))
	# println(distPoint2Line(poseAndKd[end][1][1:2], x1, x2))
	# plot(t_save, getindex.(x_save, 37)) # check cost
	# plot(T:T:T*length(vecInt), vecInt) # check integral
	# plot(eventTimes, Wavec[:,2:end]')
	# plot(1:size(Wcvec)[2], Wcvec',label=0)
	# plot(1:size(Wdvec)[2], Wdvec',label=0)
end

# hard to tune basis
function babyETC_robust_relaxedPE_journal20_A_aircraft(dotx, x, p, t)

	global eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa0, sizePlant, β, L, L1, Ld, Wavec, Wcvec, Wdvec, numDataRecord, d_real, UkU, Wc, e_buff, Ω, mu, nu, nudelay, xi, basis_on, alpha_basis, epsilon_basis, time_basis, uvec, dvec

	# global xf, uf, Tf, T, A, B, M, R, Pt, per, amp, αa, αc, Quu, Qxu, Wcfinal, uDelay, x1Delay, x2Delay, x3Delay, x_save, t_save, uvec
	xf, T, A, B, F, M, R, γ, Pt, αa, αc, αd, uDelay, x1Delay, x2Delay, x3Delay, dDelay, pDelay = p#,eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa10, Wa20, sizePlant, β, L, L1, Wavec, eventTimes=p#,i,maxIter,start = p # uu, vv useless here percent, amplitude,
    n, m, q = 3, 1, 1
	Wc = x[(n+1) : Int(n+(n+m+q)*(n+m+q+1)/2)] # 4-18
	Wd = x[Int(n+(n+m+q)*(n+m+q+1)/2+1) : Int(n+(n+m+q)*(n+m+q+1)/2+n)] # 19-21
	P = x[end] # 22
	d_real = 0#.3*sin(t)

	# Basis
	if basis_on == 1
		r = time_basis - t
		rdelay = time_basis - (t - T)
		for i = 1:Int((n+m+q)*(n+m+q+1)/2)
			for j = 1:Int((n+m+q)*(n+m+q+1)/2)
				if i+j < 23
					nu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-1) # alpha_basis*tanh(epsilon_basis*r)^(i+j-1); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-1); % cos(2*pi*(r)/10)^(i+j-1); %
					nudelay[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-1)
				elseif i+j == 23
					nu[i,j] = alpha_basis*tanh(epsilon_basis*r) # alpha_basis*tanh(epsilon_basis*r); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2); % cos(2*pi*(r)/10); %
					nudelay[i,j] = alpha_basis*tanh(epsilon_basis*r)
				else #if i+j > 23
					nu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-(((n+m+q)*(n+m+q+1)/2)+1)) # alpha_basis*tanh(epsilon_basis*r)^(i+j-22); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-22) ; % cos(2*pi*(r)/10)^(i+j-22); %
					nudelay[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-(((n+m+q)*(n+m+q+1)/2)+1))
				end
			end
		end
		for i = 1:n
			for j = 1:n
				if i+j < 6
					mu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-1) # alpha_basis*tanh(epsilon_basis*r)^(i+j-1); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-1); %cos(2*pi*(r)/10)^(i+j-1); %
					xi[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-1)
				elseif i+j == 6
					mu[i,j] = alpha_basis*tanh(epsilon_basis*r) # alpha_basis*tanh(epsilon_basis*r); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2); % cos(2*pi*(r)/100); %
					xi[i,j] = alpha_basis*tanh(epsilon_basis*r)
				else #if i+j > 6
					mu[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-(n+1)) # alpha_basis*tanh(epsilon_basis*r)^(i+j-5); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-5) ; % cos(2*pi*(r)/10)^(i+j-5); %
					xi[i,j] = alpha_basis*tanh(epsilon_basis*r)^(i+j-(n+1))
				end
			end
		end
	end

	# Update control
    # ud = zeros(m,)
	# ud[1] = Wa10'*(x[1:n] - xf)
	# ud[2] = Wa20'*(x[1:n] - xf)
	ud = Wa0'*mu*Xhat

	# Update disturbance
	d = Wd'*xi*(x[1:n] - xf)

    # Delays
    # println(t-T)
	uddelay = uDelay(t - T)
    xdelay = zeros(n,)
	xdelay[1] = x1Delay(t - T)
    xdelay[2] = x2Delay(t - T)
    xdelay[3] = x3Delay(t - T)
	ddelay = dDelay(t - T)
	pdelay = pDelay(t - T)

	# Kronecker products
	U = vcat(x[1:n] - xf, ud, d) # augmented state
    UkU = vcat(U[1]^2, U[1]*U[2], U[1]*U[3], U[1]*ud, U[1]*d,
               U[2]^2, U[2]*U[3], U[2]*ud, U[2]*d,
               U[3]^2, U[3]*ud, U[3]*d,
               ud^2, ud*d,
			   d^2)
    UkUdelay = vcat(xdelay[1]^2, xdelay[1]*xdelay[2], xdelay[1]*xdelay[3], xdelay[1]*uddelay, xdelay[1]*ddelay,
                    xdelay[2]^2, xdelay[2]*xdelay[3], xdelay[2]*uddelay, xdelay[2]*ddelay,
                    xdelay[3]^2, xdelay[3]*uddelay, xdelay[3]*ddelay,
                    uddelay^2, uddelay*ddelay,
					ddelay^2)
    Quu = Wc[13] # m x m 
    Quu_inv = 1.0/Quu
    Qux = [Wc[4] Wc[8] Wc[11]] # m x n
    Qxu = Qux' # n x m
    Qdx = [Wc[5] Wc[9] Wc[12]] # q x n 
    Qxd = Qdx'
    # Qdd = Wc[end] # q x q
	Qdd_inv = 1.0/Wc[end]

	# Triggering condition
	normXSq = U[1]^2 + U[2]^2 + U[3]^2
	trigCond = ((1 - β^2)*eigMinM*normXSq + (1-β^2)*eigMinR*ud[1]^2 - γ^2*d[1]^2)/(4*eigMaxR*(L^2+L1^2)) # d_real

	# Integral reinforcement dynamics
	# dP = 0.5*((x[1:n] - xf)'*M*(x[1:n] - xf) - (xdelay - xf)'*M*(xdelay - xf) + ud'*R*ud - uddelay'*R*uddelay - γ^2*d^2 + γ^2*ddelay^2)
	dP = 0.5*(U[1:n]'*M*U[1:n] + ud'*R*ud - γ^2*d_real^2) # d_real

	# approximation errors
	σ = nu*UkU - nudelay*UkUdelay 
	σ_f	= nu*UkU
    ec = P - pdelay + Wc'*σ #P + Wc'*σ
	ecfinal = 0.5*U[1:n]'*Pt*U[1:n] - Wc'*σ_f
	ed = Wd'*xi*U[1:n] + (Qdd_inv*Qdx*U[1:n])[1] # q x 1

	# gap
	e = Xhat - U[1:n]
	normESq = norm(e)^2

	# Critic dynamics
	relaxedPETerm = zeros(Int((n+m+q)*(n+m+q+1)/2),)
	for rrr in 1:numDataRecord
		relaxedPETerm = relaxedPETerm + e_buff[rrr]*Ω[:,rrr]./((1 + Ω[:,rrr]'*Ω[:,rrr]).^2)
	end
	dWc = -αc*((σ./(σ'*σ+1).^2)*ec + (σ_f./(σ_f'*σ_f+1).^2)*ecfinal + relaxedPETerm)

	# Triggering updates
	if normESq >= 0.95*trigCond # flows
	    Xhat = U[1:n]
		ea = Wa0'*mu*Xhat + (Quu_inv*Qux*Xhat)[1]
		Wa0 = Wa0 - αa*mu*Xhat*ea'/(1 + (mu*Xhat)'*(mu*Xhat))
		# Wavec = [Wavec Wa0]
		eventTimes = [eventTimes; t]
	end

	# Disturbance dynamics
	dWd = -αd*xi*U[1:n]*ed'/(1 + (xi*U[1:n])'*(xi*U[1:n]))

	# System dynamics
	dx = A*U[1:n] + B*ud + F*d_real

	# save vecs
	# p[end] = unew
	# p[end-1:end] .= unew
	uvec = [uvec; ud]
	normESqvec = [normESqvec; normESq]
	# p[31] = [p[31]; normESq]
	trigCondvec = [trigCondvec; trigCond]
	# p[32] = [p[32]; trigCond]
	# println("$trigCondvec")
	dvec = [dvec; d]
    # dd = d
	# println(length(vcat(dx, dWc, dWd, dP)))
	# println(length(x))
    # dotx .= vcat(dx, dWc, dWd, dP)
    dotx[1:n] = dx
    dotx[n+1:18] = dWc
    dotx[19:21] = dWd
    dotx[end] = dP
end

function sim_journal20_A_ETC_robust_relaxedPE_aircraft(x0, xf)#, S) # lastVelocity # x1 -> x2, Journal 20

	global eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa0, sizePlant, β, L, L1, Ld, Wavec, Wcvec, Wdvec, numDataRecord, d_real, UkU, Wc, e_buff, Ω, mu, nu, nudelay, xi, basis_on, alpha_basis, epsilon_basis, time_basis, uvec, dvec

	# System Dynamics for feedback
    n, m, q = 3, 1, 1

    A = [-1.01887    0.90506   -0.00215; 0.8225   -1.07741   -0.17555; 0  0  -20.2] # n x n
    B = [0 0 20.2]'#[0.8527   -2.7619    2.5127]'#[0.0002978 1.012 -0.0448 -0.1]' # n x m
    F = [0 0 1]'#[-0.1024    -0.5441     1.4631]' #[-0.09663 0.06398 0.0006301 0]' # n x q
    M = 5.0*Matrix(I, n, n) # n x n # 1, .1, 10
    R = 0.5 # 0.5*Matrix(I, m, m) # m x m # 1, 10, .1
	γ = 1
	Pt = 0.5*Matrix(I, n, n)

	# check controllability & observability
    Co = ctrb(A, B)
    unco = size(A, 1) - rank(Co)
    Ob = obsv(sqrt(M), A)
    unob = size(sqrt(M), 1) - rank(Ob)

	if unco + unob > 0	error("system uncontrollable and/or unobservable")	end

    # ODE parameters
    Tf, T, N = 10, 0.05, 200 # finite horizon
	αc, αd = 90, 3 # αa = 1.2
	
	# goal orientation
	# theta = atan((x2-x1)[2]/(x2-x1)[1])

	# Wc0 = 10*rand(28)[:]
	Wc0 = [5.0*ones(Int((n+m+q)*(n+m+q+1)/2 - m*(m+1)/2 - m*q - q^2), 1); .5; 0; -γ^2][:] # (n+m+q)(n+m+q+1)/2 x 1; 28 x 1
	# Wc0 = 20*[0.328555; 0.756487; 0.273798; 0.309953; 0.438633; 0.307851; 0.335527; 0.0851001; 0.608546; 0.282202; 0.399089; 0.281325; 0.528826; 0.234641; 0.584131; 0.313394; 0.897723; 0.675292; 0.73439; 0.445563; 0.592254; 0.596938; 0.634189; 0.332015; 0.714059; 0.511094; 0.0298947; 0.0722483]
	# Wc0 = [10.0*ones(11); 6.0*ones(11); .525; 0; 0; .525; 0; -γ^2][:] # (n+m+q)(n+m+q+1)/2 x 1; 28 x 1
	# Wc0 = [10.0*ones(4); -5; 5; 10; 10.0*ones(3); -5; 5; 10; 10.0*ones(2); -5; 5; 10; 10; -5; 5; 10; .525; 0; 10; .525; 10; -1];
	# Wc0 = [10.0*ones(4); 5cos(theta); 5sin(theta); 0.5; 10.0*ones(3); 5cos(theta); 5sin(theta); 0.5; 10.0*ones(2); 5cos(theta); 5sin(theta); 0.5; 10.0; 5cos(theta); 5sin(theta); 0.5; .525; 0; 0; .525; 0; -γ^2];
	
	# Wa10 = 5*cos(theta)*ones(n,) # n x m; 4 x 2
	# Wa20 = 5*sin(theta)*ones(n,)
	Wa0 = .5*ones(n,) # n x m; 4 x 2
	# Wa20 = .5*ones(n,)
	Wd0 = 0.5*ones(n,)
    u0 = 0#zeros(m,) # m x 1; 2 x 1
    d0 = 0
    # d_real = 1

	# Trigger parameters
	sizePlant, lixo = size(A)
	eigMinM = minimum(eigvals(M))
	eigMaxR = maximum(eigvals(R))
	eigMinR = minimum(eigvals(R))
	L = 10
	β = 1 #β1
	L1 = 1 # .9*(β*sqrt(eigMinM/eigMaxR))
	αa_UB = (8*eigMinR-4)/(eigMinR+2)
	αa = .25*αa_UB
	finalTime = zeros(0,)

	# x0 = [x1;0.0;0.0]
	# x0 = [x1; S.lastVelocity] # lastVelocity
	# xf = [x2; 0.0; 0.0]

	Xhat = x0 - xf

    t_save = [0,]
    x_save = [[x0; Wc0; Wd0; 0],] # 22
    u_save = [u0,] # u at every T
	uvec = [u0,] # all u
	d_save = [d0,] # d at every T
	dvec = [d0,] # all d

    uDelay = interp2PWC(u_save[:], -1, 1) # return an interpolation function
	x1Delay = interp2PWC(getindex.(x_save, 1), -1, 1)
    x2Delay = interp2PWC(getindex.(x_save, 2), -1, 1)
    x3Delay = interp2PWC(getindex.(x_save, 3), -1, 1)
	dDelay = interp2PWC(d_save[:], -1, 1)
	pDelay = interp2PWC(getindex.(x_save, 22), -1, 1)

	# xdist = norm(x0[1:2] - xf[1:2])
	error = 0
	localKd = 0.0
	maxIter = 10000
	poseAndKd = Array{Tuple{Array{Float64,2},Float64}}(undef,0)

	trigCondvec = zeros(0)
	normESqvec = zeros(0)
	# Wavec = zeros(3,)
	# Wcvec = zeros(15,)
	# Wdvec = zeros(3,)
	eventTimes = zeros(0,)

	numDataRecord = 0
	e_buff = Array{Float64, 1}(undef, 0)
	halfvec = [10; zeros(4); 10; zeros(3); 5; 0;0; .5; 0; -γ^2]
	Ω = Array{Float64,2}(undef, Int((n+m+q)*(n+m+q+1)/2), 0)

	# Basis parameters
	epsilon_basis = 3
	alpha_basis = 1
	time_basis = 50
	basis_on = 1 # 1-on, 2-off
	nu = zeros(Int((n+m+q)*(n+m+q+1)/2), Int((n+m+q)*(n+m+q+1)/2))
	nudelay = zeros(Int((n+m+q)*(n+m+q+1)/2), Int((n+m+q)*(n+m+q+1)/2))
	mu = zeros(n, n)
	xi = zeros(n, n)

    # solve ODEs
	for i = 1:N
		t = ((i-1)*T, i*T)
		p = [xf, T, A, B, F, M, R, γ, Pt, αa, αc, αd, uDelay, x1Delay, x2Delay, x3Delay, dDelay, pDelay]#,eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa10, Wa20, sizePlant, β, L, L1, Wavec, eventTimes]#, i,maxIter,start] percent, amplitude,
        sol = solve(ODEProblem(babyETC_robust_relaxedPE_journal20_A_aircraft, x_save[end], t, p), DP5())#, reltol = 1e-4, abstol = 1e-4, dtmax = .05)
		t_save = [t_save; sol.t[2:end]] # vcat(t_save, sol.t) save time
		x_save = [x_save; sol.u[2:end]] # vcat(x_save, sol.u) save state
		u_save = [u_save; [uvec[end]]]
		d_save = [d_save; dvec[end]]

		# Wcvec = [Wcvec sol.u[end][4:18]]
		# Wdvec = [Wdvec sol.u[end][19:21]]
		# uvec1 = p[end-2]
		# uvec2 = p[end-1]
		# trigCondvec = p[32]
		# normESqvec = p[31]
		# println("$trigCondvec")
        uDelay = interp2PWC(u_save[:], -1, i*T+.01) # interpolate control input
		x1Delay = interp2PWC(getindex.(x_save, 1), -1, i*T+.01)
		x2Delay = interp2PWC(getindex.(x_save, 2), -1, i*T+.01)
		x3Delay = interp2PWC(getindex.(x_save, 3), -1, i*T+.01)
		dDelay = interp2PWC(d_save[:], -1, i*T+.01)
		pDelay = interp2PWC(getindex.(x_save, 22), -1, i*T+.01)
		# trigCondDelay = interp2PWC(trigCondvec, -1, i*T+.01)

        # relaxed PE: record data
        if rank(Ω) < n+m+q
			∂nu∂t = zeros(Int((n+m+q)*(n+m+q+1)/2), Int((n+m+q)*(n+m+q+1)/2))
			if basis_on == 1
				r = time_basis - i*T
				for i = 1:Int((n+m+q)*(n+m+q+1)/2)
					for j = 1:Int((n+m+q)*(n+m+q+1)/2)
						if i+j < 23
							∂nu∂t[i,j] = alpha_basis*(i+j-1)*(1-tanh(epsilon_basis*r)^2)*(-epsilon_basis)*tanh(epsilon_basis*r)^(i+j-2) # alpha_basis*tanh(epsilon_basis*r)^(i+j-1); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-1); % cos(2*pi*(r)/10)^(i+j-1); %
						elseif i+j == 23
							∂nu∂t[i,j] = alpha_basis*(1-tanh(epsilon_basis*r)^2)*(-epsilon_basis) # alpha_basis*tanh(epsilon_basis*r); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2); % cos(2*pi*(r)/10); %
						else #if i+j > 23
							∂nu∂t[i,j] = alpha_basis*(i+j-(((n+m+q)*(n+m+q+1)/2)+1))*(1-tanh(epsilon_basis*r)^2)*(-epsilon_basis)*tanh(epsilon_basis*r)^(i+j-(((n+m+q)*(n+m+q+1)/2)+1)-1) # alpha_basis*tanh(epsilon_basis*r)^(i+j-22); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-22) ; % cos(2*pi*(r)/10)^(i+j-22); %
						end
					end
				end
			end
			x_dot = (x_save[end][1:n] - x_save[end-1][1:n])/T
			u_dot = (u_save[end] - u_save[end-1])/T
			d_dot = (d_save[end] - d_save[end-1])/T
			U = [x_save[end][1:n] - xf; u_save[end]; d_save[end]]
			ud = u_save[end]
			d = d_save[end]
			∂UkU∂t = [2*U[1]*x_dot[1]; x_dot[1]*U[2] + U[1]*x_dot[2]; x_dot[1]*U[3] + U[1]*x_dot[3]; x_dot[1]*ud + U[1]*u_dot; x_dot[1]*d + U[1]*d_dot;
			 		2*U[2]*x_dot[2]; x_dot[2]*U[3] + U[2]*x_dot[3]; x_dot[2]*ud + U[2]*u_dot; x_dot[2]*d + U[2]*d_dot;
				    2*U[3]*x_dot[3]; x_dot[3]*ud + U[3]*u_dot; x_dot[3]*d + U[3]*d_dot;
					2*ud*u_dot; u_dot*d + ud*d_dot;
                    2*d*d_dot]
            UkU = vcat(U[1]^2, U[1]*U[2], U[1]*U[3], U[1]*ud[1], U[1]*d,
                    U[2]^2, U[2]*U[3], U[2]*ud[1], U[2]*d,
                    U[3]^2, U[3]*ud[1], U[3]*d,
                    ud[1]^2, ud[1]*d,
                    d^2)
            ∂UkU∂̄x	= hcat([2*U[1]; U[2:5]; zeros(10,)],  # 15 x 3
                    [0; U[1]; zeros(3,); 2U[2]; U[3:5]; zeros(6,1)],
                    [0; 0; U[1]; zeros(3,); U[2]; zeros(2,); 2U[3]; U[4:5]; zeros(3,1)])
            ω = ∂nu∂t*UkU + nu*∂UkU∂t + nu*∂UkU∂̄x*x_dot
            # println(Wa0'*((x_save[end][1:n] - xf) - Xhat))
            # println(typeof(.5*(Wa0'*((x_save[end][1:n] - xf) - Xhat))'*R*(Wa0'*((x_save[end][1:n] - xf) - Xhat))))
            # println(typeof(.5*halfvec'*UkU + Wc'*ω))
            # println(((Wa0'*mu*(U[1:n] - Xhat))'*R*(Wa0'*mu*(U[1:n] - Xhat)))[1][1])
            e_buff = [e_buff; .5*halfvec'*UkU + Wc'*ω - .5*(Wa0'*mu*(U[1:n] - Xhat))'*R*(Wa0'*mu*(U[1:n] - Xhat))];
			Ω = [Ω ω]
            numDataRecord += 1
        end

		# S.kino_dist is always the max kd till now
		# localKd = distPoint2Line(x_save[end][1:2], x1, x2)
		# poseAndKd = [poseAndKd; ([x_save[end][1] x_save[end][2] x_save[end][3]], 0)]#copy(localKd))]
		# if norm(x_save[end][1:2] - xf[1:2]) < error*xdist break end
	end
	# S.lastVelocity = x_save[end][3:4]
	# lastVelocity = x_save[end][3:4]
	# (poseAndKd, normESqvec, trigCondvec)#,t_save)#eventTimes)
	# println(length(normESqvec))
	# println(length(trigCondvec))
	# plot(1:length(normESqvec), normESqvec)
	# plot!(1:length(normESqvec), trigCondvec)
	# println(length(t_save))
	# plot(getindex.(getindex.(poseAndKd, 1), 1), getindex.(getindex.(poseAndKd, 1), 2))
	# maximum(getindex.(poseAndKd, 2))
	# println(distPoint2Line(poseAndKd[end][1][1:2], x1, x2))
	# plot(t_save, getindex.(x_save, 37)) # check cost
	# plot(T:T:T*length(vecInt), vecInt) # check integral
	# plot(eventTimes, Wavec[:,2:end]')
	# plot(1:size(Wcvec)[2], Wcvec',label=0)
    # plot(1:size(Wdvec)[2], Wdvec',label=0)
	plot(t_save, [getindex.(x_save, 1) getindex.(x_save, 2) getindex.(x_save, 3)], lab = ["x1" "x2" "x3"], lw = 2, grid = :no, legend = :topright) # lc = [:pink :red],[getindex.(getindex.(poseAndKd, 1), 1) getindex.(getindex.(poseAndKd, 1), 2) getindex.(getindex.(poseAndKd, 1), 3) getindex.(getindex.(poseAndKd, 1), 4)]

end


# actually we are using this nobasis version
function babyETC_robust_relaxedPE_journal20_A_nobasis_aircraft(dotx, x, p, t)

	global eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa0, sizePlant, β, L, L1, Ld, Wavec, Wcvec, Wdvec, numDataRecord, d_real, UkU, Wc, e_buff, Ω, mu, nu, nudelay, xi, basis_on, alpha_basis, epsilon_basis, time_basis, uvec, dvec

	# global xf, uf, Tf, T, A, B, M, R, Pt, per, amp, αa, αc, Quu, Qxu, Wcfinal, uDelay, x1Delay, x2Delay, x3Delay, x_save, t_save, uvec
	xf, T, A, B, F, M, R, γ, Pt, αa, αc, αd, uDelay, x1Delay, x2Delay, x3Delay, dDelay, pDelay = p#,eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa10, Wa20, sizePlant, β, L, L1, Wavec, eventTimes=p#,i,maxIter,start = p # uu, vv useless here percent, amplitude,
    n, m, q = 3, 1, 1
	Wc = x[(n+1) : Int(n+(n+m+q)*(n+m+q+1)/2)] # 4-18
	Wd = x[Int(n+(n+m+q)*(n+m+q+1)/2+1) : Int(n+(n+m+q)*(n+m+q+1)/2+n)] # 19-21
	P = x[end] # 22
	d_real = 0#.3*sin(t)

	# Update control
    # ud = zeros(m,)
	# ud[1] = Wa10'*(x[1:n] - xf)
	# ud[2] = Wa20'*(x[1:n] - xf)
	ud = Wa0'*Xhat

	# Update disturbance
	d = Wd'*(x[1:n] - xf)

    # Delays
    # println(t-T)
	uddelay = uDelay(t - T)
    xdelay = zeros(n,)
	xdelay[1] = x1Delay(t - T)
    xdelay[2] = x2Delay(t - T)
    xdelay[3] = x3Delay(t - T)
	ddelay = dDelay(t - T)
	pdelay = pDelay(t - T)

	# Kronecker products
	U = vcat(x[1:n] - xf, ud, d) # augmented state
    UkU = vcat(U[1]^2, U[1]*U[2], U[1]*U[3], U[1]*ud, U[1]*d,
               U[2]^2, U[2]*U[3], U[2]*ud, U[2]*d,
               U[3]^2, U[3]*ud, U[3]*d,
               ud^2, ud*d,
			   d^2)
    UkUdelay = vcat(xdelay[1]^2, xdelay[1]*xdelay[2], xdelay[1]*xdelay[3], xdelay[1]*uddelay, xdelay[1]*ddelay,
                    xdelay[2]^2, xdelay[2]*xdelay[3], xdelay[2]*uddelay, xdelay[2]*ddelay,
                    xdelay[3]^2, xdelay[3]*uddelay, xdelay[3]*ddelay,
                    uddelay^2, uddelay*ddelay,
					ddelay^2)
    Quu = Wc[13] # m x m 
    Quu_inv = 1.0/Quu
    Qux = [Wc[4] Wc[8] Wc[11]] # m x n
    Qxu = Qux' # n x m
    Qdx = [Wc[5] Wc[9] Wc[12]] # q x n 
    Qxd = Qdx'
    # Qdd = Wc[end] # q x q
	Qdd_inv = 1.0/Wc[end]

	# Triggering condition
	normXSq = U[1]^2 + U[2]^2 + U[3]^2
	trigCond = ((1 - β^2)*eigMinM*normXSq + (1-β^2)*eigMinR*ud[1]^2 - γ^2*d_real^2)/(4*eigMaxR*(L^2+L1^2)) # d_real

	# Integral reinforcement dynamics
	# dP = 0.5*((x[1:n] - xf)'*M*(x[1:n] - xf) - (xdelay - xf)'*M*(xdelay - xf) + ud'*R*ud - uddelay'*R*uddelay - γ^2*d^2 + γ^2*ddelay^2)
	dP = 0.5*(U[1:n]'*M*U[1:n] + ud'*R*ud - γ^2*d_real^2) # d_real

	# approximation errors
	σ = UkU - UkUdelay 
	σ_f	= UkU
    ec = P - pdelay + Wc'*σ #P + Wc'*σ
	ecfinal = 0.5*U[1:n]'*Pt*U[1:n] - Wc'*σ_f
	ed = Wd'*U[1:n] + (Qdd_inv*Qdx*U[1:n])[1] # q x 1

	# gap
	e = Xhat - U[1:n]
	normESq = norm(e)^2

	# Critic dynamics
	relaxedPETerm = zeros(Int((n+m+q)*(n+m+q+1)/2),)
	for rrr in 1:numDataRecord
		relaxedPETerm = relaxedPETerm + e_buff[rrr]*Ω[:,rrr]./((1 + Ω[:,rrr]'*Ω[:,rrr]).^2)
	end
	dWc = -αc*((σ./(σ'*σ+1).^2)*ec + (σ_f./(σ_f'*σ_f+1).^2)*ecfinal + relaxedPETerm)

	# Triggering updates
	if normESq >= trigCond # flows
	    Xhat = U[1:n]
		ea = Wa0'*Xhat + (Quu_inv*Qux*Xhat)[1]
		Wa0 = Wa0 - αa*Xhat*ea'/(1 + (Xhat)'*(Xhat))
		# Wavec = [Wavec Wa0]
		eventTimes = [eventTimes; t]
	end

	# Disturbance dynamics
	dWd = -αd*U[1:n]*ed'/(1 + (U[1:n])'*(U[1:n]))

	# System dynamics
	dx = A*U[1:n] + B*ud + F*d_real

	# save vecs
	# p[end] = unew
	# p[end-1:end] .= unew
	uvec = [uvec; ud]
	normESqvec = [normESqvec; normESq]
	# p[31] = [p[31]; normESq]
	trigCondvec = [trigCondvec; trigCond]
	# p[32] = [p[32]; trigCond]
	# println("$trigCondvec")
	dvec = [dvec; d]
    # dd = d
	# println(length(vcat(dx, dWc, dWd, dP)))
	# println(length(x))
	# dotx .= vcat(dx, dWc, dWd, dP)
    dotx[1:n] = dx
    dotx[n+1:18] = dWc
    dotx[19:21] = dWd
    dotx[end] = dP
end

function sim_journal20_A_ETC_robust_relaxedPE_nobasis_aircraft(x0, xf)#, S) # lastVelocity # x1 -> x2, Journal 20

	global eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa0, sizePlant, β, L, L1, Ld, Wavec, Wcvec, Wdvec, numDataRecord, d_real, UkU, Wc, e_buff, Ω, mu, nu, nudelay, xi, basis_on, alpha_basis, epsilon_basis, time_basis, uvec, dvec

	# System Dynamics for feedback
    n, m, q = 3, 1, 1

    A = [-1.01887    0.90506   -0.00215; 0.8225   -1.07741   -0.17555; 0  0  -20.2] # n x n
    B = [0 0 20.2]'#[0.8527   -2.7619    2.5127]'#[0.0002978 1.012 -0.0448 -0.1]' # n x m
    F = [0 0 1]'#[-0.1024    -0.5441     1.4631]' #[-0.09663 0.06398 0.0006301 0]' # n x q
    M = 5.0*Matrix(I, n, n) # n x n # 1, .1, 10
    R = 0.5 # 0.5*Matrix(I, m, m) # m x m # 1, 10, .1
	γ = 1
	Pt = 0.5*Matrix(I, n, n)

	# check controllability & observability
    Co = ctrb(A, B)
    unco = size(A, 1) - rank(Co)
    Ob = obsv(sqrt(M), A)
    unob = size(sqrt(M), 1) - rank(Ob)

	if unco + unob > 0	error("system uncontrollable and/or unobservable")	end

    # ODE parameters
    Tf, T, N = 15, 0.05, 300 # finite horizon
	αc, αd = 90, 3 # αa = 1.2
	
	# goal orientation
	# theta = atan((x2-x1)[2]/(x2-x1)[1])

	# Wc0 = 10*rand(28)[:]
	Wc0 = [5.0*ones(Int((n+m+q)*(n+m+q+1)/2 - m*(m+1)/2 - m*q - q^2), 1); .5; 0; -γ^2][:] # (n+m+q)(n+m+q+1)/2 x 1; 28 x 1
	# Wc0 = 20*[0.328555; 0.756487; 0.273798; 0.309953; 0.438633; 0.307851; 0.335527; 0.0851001; 0.608546; 0.282202; 0.399089; 0.281325; 0.528826; 0.234641; 0.584131; 0.313394; 0.897723; 0.675292; 0.73439; 0.445563; 0.592254; 0.596938; 0.634189; 0.332015; 0.714059; 0.511094; 0.0298947; 0.0722483]
	# Wc0 = [10.0*ones(11); 6.0*ones(11); .525; 0; 0; .525; 0; -γ^2][:] # (n+m+q)(n+m+q+1)/2 x 1; 28 x 1
	# Wc0 = [10.0*ones(4); -5; 5; 10; 10.0*ones(3); -5; 5; 10; 10.0*ones(2); -5; 5; 10; 10; -5; 5; 10; .525; 0; 10; .525; 10; -1];
	# Wc0 = [10.0*ones(4); 5cos(theta); 5sin(theta); 0.5; 10.0*ones(3); 5cos(theta); 5sin(theta); 0.5; 10.0*ones(2); 5cos(theta); 5sin(theta); 0.5; 10.0; 5cos(theta); 5sin(theta); 0.5; .525; 0; 0; .525; 0; -γ^2];
	
	# Wa10 = 5*cos(theta)*ones(n,) # n x m; 4 x 2
	# Wa20 = 5*sin(theta)*ones(n,)
	Wa0 = [0.18856, 0.438721, 0.121578]#0.5*rand(n,)#ones(n,) # n x m; 4 x 2
	# println("Wa0 = $Wa0")
	# Wa20 = .5*ones(n,)
	Wd0 = [0.0745647, 0.24397, 0.120925]#0.5*rand(n,)#ones(n,)
	# println("Wd0 = $Wd0")
    u0 = Wa0'*(x0 - xf)#zeros(m,) # m x 1; 2 x 1
    d0 = Wd0'*(x0 - xf)
    # d_real = 1

	# Trigger parameters
	sizePlant, lixo = size(A)
	eigMinM = minimum(eigvals(M))
	eigMaxR = maximum(eigvals(R))
	eigMinR = minimum(eigvals(R))
	L = 10
	β = .1 #β1
	L1 = 1 # .9*(β*sqrt(eigMinM/eigMaxR))
	αa_UB = (8*eigMinR-2)/(eigMinR+2)
	αa = .25*αa_UB
	finalTime = zeros(0,)

	# x0 = [x1;0.0;0.0]
	# x0 = [x1; S.lastVelocity] # lastVelocity
	# xf = [x2; 0.0; 0.0]

	Xhat = x0 - xf

    t_save = [0,]
    x_save = [[x0; Wc0; Wd0; 0],] # 22
    u_save = [u0,] # u at every T
	uvec = [u0,] # all u
	d_save = [d0,] # d at every T
	dvec = [d0,] # all d

    uDelay = interp2PWC(u_save[:], -1, 1) # return an interpolation function
	x1Delay = interp2PWC(getindex.(x_save, 1), -1, 1)
    x2Delay = interp2PWC(getindex.(x_save, 2), -1, 1)
    x3Delay = interp2PWC(getindex.(x_save, 3), -1, 1)
	dDelay = interp2PWC(d_save[:], -1, 1)
	pDelay = interp2PWC(getindex.(x_save, 22), -1, 1)

	# xdist = norm(x0[1:2] - xf[1:2])
	error = 0
	localKd = 0.0
	maxIter = 10000
	poseAndKd = Array{Tuple{Array{Float64,2},Float64}}(undef,0)

	trigCondvec = zeros(0)
	normESqvec = zeros(0)
	# Wavec = zeros(3,)
	# Wcvec = zeros(15,)
	# Wdvec = zeros(3,)
	eventTimes = zeros(0,)

	numDataRecord = 0
	e_buff = Array{Float64, 1}(undef, 0)
	halfvec = [10; zeros(4); 10; zeros(3); 10; 0;0; .5; 0; -γ^2]
	Ω = Array{Float64,2}(undef, Int((n+m+q)*(n+m+q+1)/2), 0)

	# Basis parameters
	# epsilon_basis = 3
	# alpha_basis = 1
	# time_basis = 50
	# basis_on = 1 # 1-on, 2-off
	# nu = zeros(Int((n+m+q)*(n+m+q+1)/2), Int((n+m+q)*(n+m+q+1)/2))
	# nudelay = zeros(Int((n+m+q)*(n+m+q+1)/2), Int((n+m+q)*(n+m+q+1)/2))
	# mu = zeros(n, n)
	# xi = zeros(n, n)

    # solve ODEs
	for i = 1:N
		t = ((i-1)*T, i*T)
		p = [xf, T, A, B, F, M, R, γ, Pt, αa, αc, αd, uDelay, x1Delay, x2Delay, x3Delay, dDelay, pDelay]#,eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa10, Wa20, sizePlant, β, L, L1, Wavec, eventTimes]#, i,maxIter,start] percent, amplitude,
        sol = solve(ODEProblem(babyETC_robust_relaxedPE_journal20_A_nobasis_aircraft, x_save[end], t, p), DP5())#, reltol = 1e-4, abstol = 1e-4, dtmax = .05)
		t_save = [t_save; sol.t[2:end]] # vcat(t_save, sol.t) save time
		x_save = [x_save; sol.u[2:end]] # vcat(x_save, sol.u) save state
		u_save = [u_save; uvec[end]]
		d_save = [d_save; dvec[end]]

		# Wcvec = [Wcvec sol.u[end][4:18]]
		# Wdvec = [Wdvec sol.u[end][19:21]]
		# uvec1 = p[end-2]
		# uvec2 = p[end-1]
		# trigCondvec = p[32]
		# normESqvec = p[31]
		# println("$trigCondvec")
        uDelay = interp2PWC(u_save[:], -1, i*T+.01) # interpolate control input
		x1Delay = interp2PWC(getindex.(x_save, 1), -1, i*T+.01)
		x2Delay = interp2PWC(getindex.(x_save, 2), -1, i*T+.01)
		x3Delay = interp2PWC(getindex.(x_save, 3), -1, i*T+.01)
		dDelay = interp2PWC(d_save[:], -1, i*T+.01)
		pDelay = interp2PWC(getindex.(x_save, 22), -1, i*T+.01)
		# trigCondDelay = interp2PWC(trigCondvec, -1, i*T+.01)

        # relaxed PE: record data
        if rank(Ω) < n+m+q
			# ∂nu∂t = zeros(Int((n+m+q)*(n+m+q+1)/2), Int((n+m+q)*(n+m+q+1)/2))
			# if basis_on == 1
			# 	r = time_basis - i*T
			# 	for i = 1:Int((n+m+q)*(n+m+q+1)/2)
			# 		for j = 1:Int((n+m+q)*(n+m+q+1)/2)
			# 			if i+j < 23
			# 				∂nu∂t[i,j] = alpha_basis*(i+j-1)*(1-tanh(epsilon_basis*r)^2)*(-epsilon_basis)*tanh(epsilon_basis*r)^(i+j-2) # alpha_basis*tanh(epsilon_basis*r)^(i+j-1); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-1); % cos(2*pi*(r)/10)^(i+j-1); %
			# 			elseif i+j == 23
			# 				∂nu∂t[i,j] = alpha_basis*(1-tanh(epsilon_basis*r)^2)*(-epsilon_basis) # alpha_basis*tanh(epsilon_basis*r); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2); % cos(2*pi*(r)/10); %
			# 			else #if i+j > 23
			# 				∂nu∂t[i,j] = alpha_basis*(i+j-(((n+m+q)*(n+m+q+1)/2)+1))*(1-tanh(epsilon_basis*r)^2)*(-epsilon_basis)*tanh(epsilon_basis*r)^(i+j-(((n+m+q)*(n+m+q+1)/2)+1)-1) # alpha_basis*tanh(epsilon_basis*r)^(i+j-22); % sigmaRBF*(exp(- ((r)/ sigma_sq))^2)^(i+j-22) ; % cos(2*pi*(r)/10)^(i+j-22); %
			# 			end
			# 		end
			# 	end
			# end
			x_dot = (x_save[end][1:n] - x_save[end-1][1:n])/T
			u_dot = (u_save[end] - u_save[end-1])/T
			d_dot = (d_save[end] - d_save[end-1])/T
			U = [x_save[end][1:n] - xf; u_save[end]; d_save[end]]
			ud = u_save[end]
			d = d_save[end]
			∂UkU∂t = [2*U[1]*x_dot[1]; x_dot[1]*U[2] + U[1]*x_dot[2]; x_dot[1]*U[3] + U[1]*x_dot[3]; x_dot[1]*ud + U[1]*u_dot; x_dot[1]*d + U[1]*d_dot;
			 		2*U[2]*x_dot[2]; x_dot[2]*U[3] + U[2]*x_dot[3]; x_dot[2]*ud + U[2]*u_dot; x_dot[2]*d + U[2]*d_dot;
				    2*U[3]*x_dot[3]; x_dot[3]*ud + U[3]*u_dot; x_dot[3]*d + U[3]*d_dot;
					2*ud*u_dot; u_dot*d + ud*d_dot;
                    2*d*d_dot]
            UkU = vcat(U[1]^2, U[1]*U[2], U[1]*U[3], U[1]*ud[1], U[1]*d,
                    U[2]^2, U[2]*U[3], U[2]*ud[1], U[2]*d,
                    U[3]^2, U[3]*ud[1], U[3]*d,
                    ud[1]^2, ud[1]*d,
                    d^2)
            ∂UkU∂̄x	= hcat([2*U[1]; U[2:5]; zeros(10,)],  # 15 x 3
                    [0; U[1]; zeros(3,); 2U[2]; U[3:5]; zeros(6,1)],
                    [0; 0; U[1]; zeros(3,); U[2]; zeros(2,); 2U[3]; U[4:5]; zeros(3,1)])
            ω = ∂UkU∂t + ∂UkU∂̄x*x_dot
            # println(Wa0'*((x_save[end][1:n] - xf) - Xhat))
            # println(typeof(.5*(Wa0'*((x_save[end][1:n] - xf) - Xhat))'*R*(Wa0'*((x_save[end][1:n] - xf) - Xhat))))
            # println(typeof(.5*halfvec'*UkU + Wc'*ω))
            # println(((Wa0'*mu*(U[1:n] - Xhat))'*R*(Wa0'*mu*(U[1:n] - Xhat)))[1][1])
            e_buff = [e_buff; .5*halfvec'*UkU + Wc'*ω - .5*(Wa0'*(U[1:n] - Xhat))'*R*(Wa0'*(U[1:n] - Xhat))];
			Ω = [Ω ω]
            numDataRecord += 1
        end

		# S.kino_dist is always the max kd till now
		# localKd = distPoint2Line(x_save[end][1:2], x1, x2)
		# poseAndKd = [poseAndKd; ([x_save[end][1] x_save[end][2] x_save[end][3]], 0)]#copy(localKd))]
		# if norm(x_save[end][1:2] - xf[1:2]) < error*xdist break end
	end
		# S.lastVelocity = x_save[end][3:4]
	# lastVelocity = x_save[end][3:4]
	# (poseAndKd, normESqvec, trigCondvec)#,t_save)#eventTimes)
	# println(length(normESqvec))
	# println(length(trigCondvec))
	# println(length(poseAndKd))
	# plot(1:length(normESqvec), normESqvec)
	# plot!(1:length(normESqvec), trigCondvec)
	println(length(eventTimes))
	# println(eventTimes)
	# println((eventTimes[2:end] - eventTimes[1:end-1]))
	scatter(eventTimes[2:end], (eventTimes[2:end] - eventTimes[1:end-1]), lab = false, m = :cross, lw = 2, grid = :no, legend = :topright)

	# println(getindex.(u_save,1))

	# println(getindex.(x_save, 1))
	# println(getindex.(x_save, 2))
	# println(getindex.(x_save, 3))
    
	# println(x_save[end][1])
	# println(x_save[end][2])
	# println(x_save[end][3])

	# println(t_save)

	# println(trigCondvec)
	# println(normESqvec)

	# plot(range(0, stop = 15, length = length(normESqvec)), [normESqvec trigCondvec], lab = ["squared norm error" "threshold"], lc = [:pink :red], lw = 2, grid = :no, legend = :topright)
	# plot(range(0, stop = 15, length = length(u_save)), u_save, lab = "u", lc = :red, lw = 2, grid = :no, legend = :topright)

	# plot(t_save, [getindex.(x_save, 1) getindex.(x_save, 2) getindex.(x_save, 3)], lab = ["x1" "x2" "x3"], lw = 2, grid = :no, legend = :topright) # lc = [:pink :red],[getindex.(getindex.(poseAndKd, 1), 1) getindex.(getindex.(poseAndKd, 1), 2) getindex.(getindex.(poseAndKd, 1), 3) getindex.(getindex.(poseAndKd, 1), 4)]

	# plot(getindex.(getindex.(poseAndKd, 1), 1), getindex.(getindex.(poseAndKd, 1), 2))
	# maximum(getindex.(poseAndKd, 2))
	# println(distPoint2Line(poseAndKd[end][1][1:2], x1, x2))
	# plot(t_save, getindex.(x_save, 37)) # check cost
	# plot(T:T:T*length(vecInt), vecInt) # check integral
	# plot(eventTimes, Wavec[:,2:end]')
	# plot(1:size(Wcvec)[2], Wcvec',label=0)
	# plot(1:size(Wdvec)[2], Wdvec',label=0)

end