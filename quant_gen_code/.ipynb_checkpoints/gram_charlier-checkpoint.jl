module gram_charlier


function factorial(i::Int)
    if i == 0
        return 1
    else
        return i*factorial(i-1)
    end    
end 

function phi(x::Float64)::Float64
    return exp(-0.5*x^2)/sqrt(2*pi)
end 

function H(x::Float64,i::Int)::Float64
    @assert i in [3,4,5]
    if i == 3
        return x^3 - 3*x
    elseif i == 4
        return x^4 - 6*x^2 + 3
    else
        return x^5 - 10*x^3 + 15*x
    end
end 

mutable struct aproximation
    degree::Int
    moments::AbstractVector{Float64}
end 

function evaluate(x::Float64, aproximation)::Float64
    n = aproximation.degree
    mean = aproximation.moments[1]
    var = aproximation.moments[2]
    x = (x - mean)/sqrt(var)
    d = phi(x)
    for i in 3:n
        d += d*H(x,i)*aproximation.moments[i]/factorial(i) 
    end 
    return d
end 

end # module 