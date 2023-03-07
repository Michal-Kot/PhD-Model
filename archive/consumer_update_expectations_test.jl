unit_possessed = [false, false]
iter = 1
unit_first_possessed_time = 1
λ_ind = 0.5
λ_wom = 0.25
quality_expectation = [0.66, 0.53]
quality_of_unit_possessed = [0.55, 0]
buying_neighbours = [0,0]
average_quality_neighbours = [0.59, 0]

if any(unit_possessed)

    if unit_first_possessed_time == iter

        println(quality_expectation .+ 
            λ_ind .* unit_possessed .* (quality_of_unit_possessed .- quality_expectation) .+ 
            λ_wom .* sign.(buying_neighbours) .* (average_quality_neighbours .- quality_expectation))
        
    else

        println(quality_expectation .+ 
            λ_wom .* sign.(buying_neighbours) .* (average_quality_neighbours .- quality_expectation))
        
    end

else

    println(quality_expectation .+ 
        λ_wom .* sign.(buying_neighbours) .* (average_quality_neighbours .- quality_expectation))

end