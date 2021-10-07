"""
This module defines the deriitives and jacobian for the model
described in Tufto Am. Nats. 2001.
"""
module tufto_2001

include("utils.jl")
using LinearAlgebra
using NLsolve
using Plots
using DifferentialEquations


params = [1,1,1,1,3.5,0.0]

function derivitives(x,params)
    N = x[1]
    z = x[2]
    
    r = params[1]
    K = params[2]
    V = params[3]
    s = params[4]
    z1 = params[5]
    m = params[6]
    
    
    dN = N*(r*(1-N/K) - s/2 * (V+z^2)) + m
    dz = -V*s*z - m/N*(z-z1)
    
    return [dN,dz]
end 


function jacobian(x,params)
    N = x[1]
    z = x[2]
    
    r = params[1]
    K = params[2]
    V = params[3]
    s = params[4]
    z1 = params[5]
    m = params[6]
    
    NdN = r*(1-N/K) - r*N/K - s*(z^2+V)/2
    Ndz = -s*N*z
    zdz = -V*s-m/N
    zdN = m*(z-z1)/N^2
    return [NdN Ndz; zdN zdz]
end 




function get_zeros(params, n; tol = 0.00001)
    
    function f!(f,x)
        y = derivitives(x,params)
        f[1] = y[1]
        f[2] = y[2]
        #return tufto_2001.derivitives(x,par)
    end

    function J!(J,x)
        J1 =  jacobian(x,params)
        J[1,1] = J1[1,1]
        J[2,1] = J1[2,1]
        J[1,2] = J1[1,2]
        J[2,2] = J1[2,2]
    end 
    
    eqs = zeros(n,2)
    for i in 1:n
        init_x = (rand(2) .- [0.0,0.25]) .* [0.75,4.0]
        opt = NLsolve.nlsolve(f!,J!,init_x, method = :newton, iterations = 10000)
        y = opt.zero
        check = mapslices(x -> sum((x.-y).^2), eqs, dims = 2)
        if !any(check .< tol)
            eqs[i,:] = opt.zero
        end
    end 
    
    inds = mapslices(x -> !any(x.<=0), eqs, dims = 2)
    return eqs[reshape(inds,n),:]
end


function max_eigen(x,params)
    J = jacobian(x,params)
    v = eigen(J).values
    return v[argmax(v)]
end 

function get_stability(eqs, params)
    n = size(eqs)[1]
    vals = mapslices(x -> max_eigen(x,params) < 0, eqs, dims = 2)
    return reshape(vals,n)
end



function phase_portrait(params)
    ###### plot a phase portrait
    Ns = 0.01:0.1:1.01
    vs = -1:0.5:4

    df(x) = tufto_2001.derivitives(x,params)./20
    xxs = [N for N in Ns for v in vs]
    yys = [v for N in Ns for v in vs]
    uv = mapslices(df, hcat(xxs,yys), dims = 2)
    
    # plot points and arrows with 10% head sizes
    p1 = Plots.scatter(ylim = [-1.5,4.5],xlim = [-0.1,1.0])
    for (x,y,u,v) in zip(xxs,yys,uv[:,1],uv[:,2])
        utils.arrow0!( x, y, u, v; as=0.1, lc=:blue, la=1)
    end
    
    eqs = get_zeros(params, 100)

    
    inds = get_stability(eqs, params)
    
#     
    if any(inds .!= 1)
        Plots.scatter!(p1,eqs[inds.!= 1,1],eqs[inds .!= 1,2], color = "red")
        
        point = eqs[inds.!= 1,:]
        point = point[1,:]
        J = jacobian(point,params)
        v = eigen(J).values
        V = eigen(J).vectors
        
        v1 = eigen(J).vectors[:,argmax(v)]
        v2 = eigen(J).vectors[:,argmin(v)]
        
        u0 = point .+ 0.05.*v2
        f(u,p,t) = -1 .*tufto_2001.derivitives(u,params)
        tspan = (0.0,1.5)
        prob = ODEProblem(f,u0,tspan)
        sol = solve(prob)
        Plots.plot!(p1, sol,vars = (1,2),color = "red")
        
        u0 = point .- 0.05.*v2
        f(u,p,t) = -1 .*tufto_2001.derivitives(u,params)
        tspan = (0.0,1.25)
        prob = ODEProblem(f,u0,tspan)
        sol = solve(prob)
        Plots.plot!(p1, sol,vars = (1,2),color = "red")
        
    end
    
    Plots.scatter!(p1,eqs[inds,1],eqs[inds,2], color = "black",ylim = [-1.5,4.5],xlim = [-0.1,1.0])
    
    return p1
end








param_nd = [1.0,8.0, 0.1]

function derivitives_nd(x,params)
    nu = x[1]
    omega = x[2]
    
    lambda = params[1]
    sigma = params[2]
    mu = params[3]
    
    
    dnu = nu*((1-(1-lambda)*nu)-lambda*(1+omega^2))+mu
    domega = -2*lambda*omega - mu*nu^-1*(omega - sigma)
    
    return [dnu,domega]
end 


function jacobian_nd(x,params)
    nu = x[1]
    omega = x[2]
    
    lambda = params[1]
    sigma = params[2]
    mu = params[3]
    

    nudnu = 1 + -2*(1-lambda)*nu - lambda*(1+omega^2)
    nudomega = -2*nu*lambda*omega
    omegadomega = -2*lambda - mu/nu
    omegadnu = mu*(omega-sigma)/nu^2
    return [nudnu nudomega; omegadnu omegadomega]
end 







function get_zeros_nd(params, n; tol = 0.00001)
    
    function f!(f,x)
        y = derivitives_nd(x,params)
        f[1] = y[1]
        f[2] = y[2]
        #return tufto_2001.derivitives(x,par)
    end

    function J!(J,x)
        J1 =  jacobian_nd(x,params)
        J[1,1] = J1[1,1]
        J[2,1] = J1[2,1]
        J[1,2] = J1[1,2]
        J[2,2] = J1[2,2]
    end 
    
    eqs = zeros(n,2)
    for i in 1:n
        init_x = (rand(2) .- [0.0,0.25]) .* [2.0,4.0]
        opt = NLsolve.nlsolve(f!,J!,init_x, method = :newton, iterations = 10000,ftol = 10^-10)
        y = opt.zero
        check = mapslices(x -> sum((x.-y).^2), eqs, dims = 2)
        if !all(check .< tol) && converged(opt)
            eqs[i,:] = opt.zero
        end

    end 
    
    inds = mapslices(x -> !any(x.<=0), eqs, dims = 2)
    return eqs[reshape(inds,n),:]
end


function max_eigen_nd(x,params)
    J = jacobian_nd(x,params)
    v = eigen(J).values

    return v[argmax(v)]
end 

function get_stability_nd(eqs, params)
    n = size(eqs)[1]
    vals = mapslices(x -> max_eigen_nd(x,params) < 0, eqs, dims = 2)
    return reshape(vals,n)
end



function phase_portrait_nd(params)
    ###### plot a phase portrait
    sigma = params[2]
    Ns = 0.01:0.2:2.01
    vs = -(sigma/10):(sigma/10):(sigma + sigma/10)

    df(x) = tufto_2001.derivitives_nd(x,params)./20
    xxs = [N for N in Ns for v in vs]
    yys = [v for N in Ns for v in vs]
    uv = mapslices(df, hcat(xxs,yys), dims = 2)
    
    # plot points and arrows with 10% head sizes
    p1 = Plots.scatter(ylim = [-1,1.4],xlim = [-0.1,2.01])
    for (x,y,u,v) in zip(xxs,yys,uv[:,1],uv[:,2])
        utils.arrow0!( x, y, u, v; as=0.1, lc=:blue, la=1)
    end
    
    eqs = get_zeros_nd(params, 500)


    inds = get_stability_nd(eqs, params)
    
#     
    if any(inds .!= 1)
        Plots.scatter!(p1,eqs[inds.!= 1,1],eqs[inds .!= 1,2], color = "red")
        
        point = eqs[inds.!= 1,:]
        point = point[1,:]
        J = jacobian_nd(point,params)
        v = eigen(J).values
        V = eigen(J).vectors
        
        v1 = eigen(J).vectors[:,argmax(v)]
        v2 = eigen(J).vectors[:,argmin(v)]
        
        u0 = point .+ 0.05.*v2
        f(u,p,t) = -1 .*tufto_2001.derivitives_nd(u,params)
        tspan = (0.0,1.5)
        prob = ODEProblem(f,u0,tspan)
        sol = solve(prob)
        Plots.plot!(p1, sol,vars = (1,2),color = "red")
        
        u0 = point .- 0.05.*v2
        f(u,p,t) = -1 .*tufto_2001.derivitives_nd(u,params)
        tspan = (0.0,1.25)
        prob = ODEProblem(f,u0,tspan)
        sol = solve(prob)
        Plots.plot!(p1, sol,vars = (1,2),color = "red")
        
    end
    
    Plots.scatter!(p1,eqs[inds,1],eqs[inds,2], color = "black",ylim = [-(sigma/10),(sigma + sigma/10)],xlim = [-0.1,2.01])
    
    lambda = params[1]
    sigma = params[2]
    mu = params[3]
    
    omgea_0(x) = mu*sigma/(2*lambda+mu*x^-1)/x
    function nu_0(x) 
        v = 1/lambda - (1-lambda)*x/lambda - 1 + mu/(x*lambda)
        if v > 0
            return sqrt(v)
        end 
        return NaN
    end 
    
    x = 0.01:0.01:2.01
    y_omega = omgea_0.(x)
    y_nu = nu_0.(x)
    
    Plots.plot!(p1,x,y_omega, color = "black", label = "omega = 0")
    Plots.plot!(p1,x,y_nu, color = "grey", label = "nu = 0")
    return p1
end

## non autonomous


function derivitives_na(x,m,params,t)
    N = x[1]
    z = x[2]
    
    r = params[1]
    K = params[2]
    V = params[3]
    s = params[4]
    z1 = params[5]
    
    
    dN = N*(r*(1-N/K) - s/2 * (V+z^2)) + m(t)
    dz = -V*s*z - m(t)/N*(z-z1)
    
    return [dN,dz]
end 

end # module 