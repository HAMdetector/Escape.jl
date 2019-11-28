struct EpitopeMap 
    name::String
    epitope::Vector{String}
    start::Vector{Int}
    stop::Vector{Int}
    allele::Vector{HLAAllele}

    function EpitopeMap(name, epitope, start, stop)
        if !(length(epitope) == length(start) == length(stop))
            error("epitope sequence and start and stop positions" * 
                "must have the same length.")
        elseif !all(start .<= stop)
            error("epitope start must be less than epitope stop.")
        elseif !(all((start .> 0) .& (stop .> 0)))
            error("start and stop positions must be positive.")
        elseif !all(length.(epitope) .== (stop .- start .+ 1))
            error("epitope length must match the given start and stop positions.")
        end

        new(name , epitope, start, stop)
    end
end

