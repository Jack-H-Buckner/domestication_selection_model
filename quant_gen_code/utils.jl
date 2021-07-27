module utils

normal_density(mu, sigma, x) = 1/(sigma*sqrt(2*pi))*exp(-1/2*((x - mu)/sigma)^2)

squared_exponential(opt,scale,x) = exp(-1/2*((x - opt)/scale)^2)

dsigma_normal_density(mu,sigma,x) = - ((sigma^2 - x^2+ 2*mu*x-mu^2)*exp(-1/2*((x-mu)/sigma)^2))/(sigma^4*sqrt(2*pi))


function conv(step_function, differnce, product, n)
    vals = abs.(differnce.- step_function.nodes[n])
    tf = vals .< log(step_function.n)*step_function.interval
    return sum( tf .* product)/(2*log(step_function.n))
end

function conv_2(f,g,product,difference,n)
    @assert f.interval == g.interval
    return sum((abs.(difference .- f.nodes[n]) .< f.interval) .* product )/2
end 

end # module