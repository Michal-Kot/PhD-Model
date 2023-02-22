function calculate_state_profit(K::Float64, eK_dist::Vector{Float64}, D::Float64, eD_dist::Vector{Float64}, M::Float64, Q::Float64, o_K::Float64, o_D::Float64, o_P::Float64, N::Int64, μ_c::Float64, cc::Float64, ρ_dist::Vector{Float64}, β_dist::Vector{Float64}, return_type::String, product_life::Int64)
    """

    Funkcja licząca oczekiwany zysk z danego stanu. Wykorzystywana przez firmę do szacowania efektu zmiany stanu.

    """

    s = β_dist
    eK = eK_dist
    eD = eD_dist
    eρ = ρ_dist

    @assert all(0 .<= s .<= 1)
    @assert all(0 .<= eK .<= 1)
    @assert all(0 .<= eD .<= 1)
    @assert all(0 .<= eρ .<= 1)

    o_U = s .* sum_of_geom_series_finite(o_K, eρ * o_D; t = product_life) .- o_P # użyteczność dobra konkurencji, jeśli liczba konkurentów > 1, to o_k, o_D i o_P są średnimi

    U = s .* sum_of_geom_series_finite.(eK, eρ .* eD; t = product_life)  .- cost_coefficient(K, D, cc) .* sum_of_geom_series_infinite(K, D) .* M # użyteczność mojego dobra przy parametrach K, D, M

    demand = sum((U .> 0) .& (U .> o_U) .& (rand(N) .< 1/product_life)) # szacowany popyt. warunek 1: moja użyteczność > 0, warunek 2: moja użyteczność wyższa niż użyteczność dobra konkurencyjnego, warunek 3: oczekiwana liczba klientów poszukujących dobra - skalowanie dla dóbr trwałych > 1 okres

    @assert demand >= 0

    price = cost_coefficient(K, D, cc) * sum_of_geom_series_finite(K, D; t = product_life) * M # marża na 1 sprzedanym produkcie

    @assert price >= 0

    profit = min(demand,Q) .* price .+ min(0, Q - demand) .* (1 - μ_c) .* cost_coefficient(K, D, cc) .* sum_of_geom_series_finite(K, D; t = product_life) - Q * cost_coefficient(K, D, cc) .* sum_of_geom_series_finite(K, D; t = product_life)  # oczekiwany zysk firmy

    @assert min(demand, Q) >= 0
    @assert 1 - μ_c >= 0
    @assert cost_coefficient(K, D, cc) >= 0
    @assert sum_of_geom_series_finite(K, D; t = product_life) >= 0
    
    if return_type == "profit"
        return profit
    elseif return_type == "demand"
        return demand
    end

end