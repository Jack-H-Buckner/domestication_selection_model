module utils

function collapse_v(v)
    return v[1:2:length(v)] .+ v[2:2:(length(v)+1)]
end 

function colapse_M(M)
    M = M[:,1:2:(size(M)[2])] .+ M[:,2:2:(size(M)[2]+1)]
    M = M[1:2:(size(M)[1]),:] .+ M[2:2:size(M)[1],:]
    return M
end 

end 