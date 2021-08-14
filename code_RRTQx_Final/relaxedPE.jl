# ode fcn for simulation B in TNNLS with ETC, H-infinity, and relaxed PE
function babyETC_robust_relaxedPE_B(dotx, x, p, t)

	global eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa10, Wa20, sizePlant, β, L, L1, Ld, Wavec, Wcvec, Wdvec, numDataRecord, d_real, UkU, Wc, e_buff, Ω, mu, nu, nudelay, xi, basis_on, alpha_basis, epsilon_basis, time_basis, uvec, dvec

	# global xf, uf, Tf, T, A, B, M, R, Pt, per, amp, αa, αc, Quu, Qxu, Wcfinal, uDelay, x1Delay, x2Delay, x3Delay, x_save, t_save, uvec
	xf, T, A, B, F, M, R, γ, Pt, αa, αc, αd, u1Delay, u2Delay, x1Delay, x2Delay, x3Delay, x4Delay, dDelay, pDelay = p#,eigMinM, eigMaxR, eigMinR, Xhat, normESqvec, trigCondvec, eventTimes, Wa10, Wa20, sizePlant, β, L, L1, Wavec, eventTimes=p#,i,maxIter,start = p # uu, vv useless here percent, amplitude,
    n, m, q = 4, 2, 1
	Wc = x[(n+1) : Int(n+(n+m+q)*(n+m+q+1)/2)] # 5-32
	Wd = x[Int(n+(n+m+q)*(n+m+q+1)/2+1) : Int(n+(n+m+q)*(n+m+q+1)/2+n)] # 33-36
	P = x[end] # 37
	d_real = 0#3sin(t)

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
	Quu = [Wc[23] Wc[24]; Wc[24] Wc[26]] # m x m 
	Quu_inv = inv(Quu)
	Qux = [Wc[5] Wc[11] Wc[16] Wc[20]; Wc[6] Wc[12] Wc[17] Wc[21]] #reshape(Wc[Int((n+m+q)*(n+m+q+1)/2-m^2-q^2-m*n+1):Int((n+m+q)*(n+m+q+1)/2-m^2-q^2)], (n, m)) # n x m 16-23
	Qxu = Qux' # m x n
	Qdx = [Wc[7] Wc[13] Wc[18] Wc[22]]#reshape(Wc[Int((n+m+q)*(n+m+q+1)/2-m^2-q^2-m*n-n*q+1):Int((n+m+q)*(n+m+q+1)/2-m^2-q^2-m*n)], (n, q)) # n x q 12-15
	Qxd = Qdx'
    # Qdd = Wc[end] # q x q
	Qdd_inv = 1.0/Wc[end]

	# Triggering condition
	normXSq = abs(U[1])^2 + abs(U[2])^2 + U[3]^2 + U[4]^2
	trigCond = (((1 - β^2)*eigMinM - γ^2*Ld^2)/(4*(L^2+L1^2)*eigMaxR))*normXSq + (eigMinR/(4*(L^2+L1^2)*eigMaxR))*(abs(ud[1])^2+abs(ud[2])^2) - γ^2*d^2/(4*(L^2+L1^2)*eigMaxR)

	# Integral reinforcement dynamics
	# dP = 0.5*((x[1:n] - xf)'*M*(x[1:n] - xf) - (xdelay - xf)'*M*(xdelay - xf) + ud'*R*ud - uddelay'*R*uddelay - γ^2*d^2 + γ^2*ddelay^2)
	dP = 0.5*(U[1:n]'*M*U[1:n] + ud'*R*ud - γ^2*d^2)

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
		Wa10 = Wa10 - αa*Xhat*ea1'/(1 + Xhat'*Xhat)
		Wa20 = Wa20 - αa*Xhat*ea2'/(1 + Xhat'*Xhat)
		Wavec = [Wavec [Wa10; Wa20]]
		eventTimes = [eventTimes; t]
	end

	# Disturbance dynamics
	dWd = -αd*U[1:n]*ed'

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



function sim_TNNLS_B_ETC_robust_relaxedPE_Local(x1, x2)#, S) # lastVelocity # x1 -> x2, Journal 20

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

    A = [0 0 1 0; 0 0 0 1; -kx/m_total 0 -cx/m_total 0; 0 -ky/m_total 0 -cx/m_total] # n x n
    B = [0 0; 0 0; 1/m_total 0; 0 1/m_total] # n x m
    F = [0; 0; 1; 1] # n x q
    M = 10.0*Matrix(I, n, n) # n x n # 1, .1, 10
    R = 0.525*Matrix(I, m, m) # m x m # 1, 10, .1
	γ = 10
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

    Wc0 = [10.0*ones(Int((n+m+q)*(n+m+q+1)/2 - m^2 - q^2), 1); reshape(R, (m^2, 1)); -γ^2][:] # (n+m+q)(n+m+q+1)/2 x 1; 28 x 1
	Wa10 = 0.5*ones(n,) # n x m; 4 x 2
	Wa20 = 0.5*ones(n,)
	Wd0 = 0.5*ones(n,)
    u0 = zeros(m,) # m x 1; 2 x 1
    d0 = 0.5
    # d_real = 1

	# Trigger parameters
	sizePlant, lixo = size(A)
	eigMinM = minimum(eigvals(M))
	eigMaxR = maximum(eigvals(R))
	eigMinR = minimum(eigvals(R))
	L, β, Ld = 30, .6, 30
	L1 = 30 # .9*(β*sqrt(eigMinM/eigMaxR))
	αa_UB = (8*eigMinR-4)/(eigMinR+2)
	αa = .25*αa_UB
	finalTime = zeros(0,)

	x0 = [x1;0.0;0.0]
	# x0 = [x1; S.lastVelocity] # lastVelocity
	xf = [x2; 0.0; 0.0]

    # p0 = (x0 - xf)'*M*(x0 - xf)
	# integral = 0.0
	# vecInt = [0.0,]

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

	xdist = norm(x0[1:2] - xf[1:2])
	error = 0.1
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
	halfvec = [10; zeros(6); 10; zeros(5); 10; 0;0;0;0; 10; 0;0;0; .525; 0;0; .525; 0; -γ^2]
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
        sol = solve(ODEProblem(babyETC_robust_relaxedPE_B, x_save[end], t, p), DP5())#, reltol = 1e-4, abstol = 1e-4, dtmax = .05)
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
		localKd = distPoint2Line(x_save[end][1:2], x1, x2)
		poseAndKd = [poseAndKd; ([x_save[end][1] x_save[end][2]], copy(localKd))]
		if norm(x_save[end][1:2] - xf[1:2]) < error*xdist break end
	end
	# S.lastVelocity = x_save[end][3:4]
	(poseAndKd, normESqvec, trigCondvec)#,t_save)#eventTimes)
	plot(getindex.(getindex.(poseAndKd, 1), 1), getindex.(getindex.(poseAndKd, 1), 2))
	# plot(t_save, getindex.(x_save, 37)) # check cost
	# plot(T:T:T*length(vecInt), vecInt) # check integral
	# plot(eventTimes, Wavec[:,2:end]')
	# plot(1:size(Wcvec)[2], Wcvec',label=0)
	# plot(1:size(Wdvec)[2], Wdvec',label=0)
end
