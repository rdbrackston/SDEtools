module gMAM

#  - Evaluate gMAM, the geometric, time-independent action function
#  - Include non-gradient optimisation methods, e.g. GA, simulated annealing


"""
Evaluate the time-independent geometric action for a discretised path φ.
"""
function geometricaction(φ::AbstractArray,
           f::Function, g::Function,
           X₀::AbstractArray, Xₑ::AbstractArray,
           N::Int64, n::Int64)

    # Loop over the elements in φ


end


end
