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



using Plots; gr(legend=false)
print("!!")
# as: arrow head size 0-1 (fraction of arrow length)
# la: arrow alpha transparency 0-1
function arrow0!(x, y, u, v; as=0.07, lc=:black, la=1)
    nuv = sqrt(u^2 + v^2)
    v1, v2 = [u;v] / nuv,  [-v;u] / nuv
    v4 = (3*v1 + v2)/3.1623  # sqrt(10) to get unit vector
    v5 = v4 - 2*(v4'*v2)*v2
    v4, v5 = as*nuv*v4, as*nuv*v5
    plot!([x,x+u], [y,y+v], lc=lc,la=la)
    plot!([x+u,x+u-v5[1]], [y+v,y+v-v5[2]], lc=lc, la=la)
    plot!([x+u,x+u-v4[1]], [y+v,y+v-v4[2]], lc=lc, la=la)
end


end # module