module Decompose

function decompose(n,nfactors)
    if nfactors == 1
        return [[n]]
    end

    result = Vector{Vector{Int32}}()
    for factor in 1:n
        if n%factor == 0
            subdecompositions = decompose(n/factor,nfactors-1)
            for subdecomposition in subdecompositions
                push!(result,[factor,subdecomposition...])
            end
        end
    end
    return result
end

end
