#### AUX functions

function sum_of_geom_series_infinite(a0,q)
    return a0 ./ (1 .- q.^2)
end

function sum_of_geom_series_finite(a0,q)
    t = 1 ./ (1 .- q)
    return a0 .* (1 .- q.^t) ./ (1 .- q)
end

function in_boundaries(x::Float64,lb::Float64,ub::Float64)::Float64
    return max(min(x,ub),lb)
end

function in_boundaries(x::Float64,lb::Int64,ub::Int64)::Float64
    return max(min(x,ub),lb)
end

function probability_still_valid(d,t)
    return (1-d)^(t-1)*d
end

function cost_coefficient(k,d,cc)
    return cc
end

function softmax(x::Vector, scale::Bool = true)

    y = copy(x)

    if scale

        y = y / std(y)

    end



end

function create_weights(x::Vector, method::String, scale::Bool = true)

    y = copy(x)

    @assert (method == "softmax") | (method == "relu") | (method == "argmax")

    if method == "softmax"

        if scale

            y = y / std(y)

        end

        return exp.(y) ./ sum(exp.(y)) 

    elseif method == "relu"

        return max.(0,y)

    elseif method == "argmax"

        w = fill(0, maximum(axes(y,1)))
        w[argmax(y)] = 1

        return w

    end

end

function calculate_state_profit(K::Float64, eK_dist::Vector{Float64}, D::Float64, eD_dist::Vector{Float64}, M::Float64, Q::Float64, o_K::Float64, o_D::Float64, o_P::Float64, N::Int64, μ_c::Float64, cc::Float64, ρ_dist::Vector{Float64}, return_type::String)
    """

    Funkcja licząca oczekiwany zysk z danego stanu. Wykorzystywana przez firmę do szacowania efektu zmiany stanu.

    """

    s = rand(Uniform(0, 1), N) # standard reservation price, założone dla N klientów

    eK = sample(eK_dist, N)
    eD = sample(eD_dist, N)
    eρ = sample(ρ_dist, N)

    o_U = s .* sum_of_geom_series_infinite(o_K, eρ * o_D) .- o_P # użyteczność dobra konkurencji, jeśli liczba konkurentów > 1, to o_k, o_D i o_P są średnimi

    Us = s .* sum_of_geom_series_infinite.(eK, eρ .* eD)  .- cost_coefficient(K, D, cc) .* sum_of_geom_series_infinite(K, D) .* M # użyteczność mojego dobra przy parametrach K, D, M

    Ul = s .* sum_of_geom_series_infinite.(eK, eρ .* eD) .- cost_coefficient(K, D, cc) .* sum_of_geom_series_infinite(K, D) .* M .* (1.1) .* (1-D) .* (1 .- eρ .^ (1 / (1-D))) ./ (1 .- eρ)
    
    U = max(Us, Ul)

    demand = sum((U .> 0) .& (U .> o_U) .& (rand(N) .< (1-D))) # szacowany popyt. warunek 1: moja użyteczność > 0, warunek 2: moja użyteczność wyższa niż użyteczność dobra konkurencyjnego, warunek 3: oczekiwana liczba klientów poszukujących dobra - skalowanie dla dóbr trwałych > 1 okres

    price = cost_coefficient(K, D, cc) * sum_of_geom_series_infinite(K, D) * M # marża na 1 sprzedanym produkcie

    profit = min(demand,Q) .* price .+ max.(0, Q .- demand) .* (1 - μ_c) .* cost_coefficient(K, D, cc) .* sum_of_geom_series_infinite(K, D) - Q * cost_coefficient(K, D, cc) .* sum_of_geom_series_infinite(K, D)  # oczekiwany zysk firmy

    if return_type == "profit"
        return profit
    elseif return_type == "demand"
        return demand
    end

end

function calculate_cost(_seller::seller)::Float64
    cost = sum_of_geom_series_infinite(_seller.quality, _seller.durability) * cost_coefficient(_seller.quality, _seller.durability, _seller.cost_coefficient) #_seller.cost_coefficient
    return cost
end


function calculate_cost_history(_seller::seller)::Vector{Float64}
    cost = sum_of_geom_series_infinite.(_seller.quality_history, _seller.durability_history) .* cost_coefficient.(_seller.quality_history, _seller.durability_history, _seller.cost_coefficient) # _seller.cost_coefficient
    return cost
end

function calculate_price(_seller::seller)::Float64
    price = calculate_cost(_seller) * _seller.margin
    return price
end

function calculate_price_history(_seller::seller)::Vector{Float64}
    price = calculate_cost_history(_seller) .* _seller.margin_history
    return price
end

function calculate_lease_single(_seller::seller, interest_rate::Float64)::Float64
    lease = calculate_cost(_seller) * _seller.margin * interest_rate * (1 - _seller.durability)
    return lease
end

function calculate_lease_total(_seller::seller, interest_rate::Float64)::Float64
    lease = calculate_cost(_seller) * _seller.margin * interest_rate
    return lease
end




function u2w(u::Vector{Float64}, p_min::Float64 = 0.1)::Vector{Float64}

    if p_min > 0

        u_min = count(u .== minimum(u))

        ut = u .- minimum(u)
        wgt = ut ./ sum(ut) .* (1 - u_min * p_min)
        wgt[argmin(u)] = p_min

    else

        wgt = u ./ sum(u)

    end

    return wgt
end

function calculate_profit_history(_seller::seller; trim=1)::Vector{Float64}
    profit_history = _seller.selling_income_history .+ _seller.leasing_income_history .- _seller.cost_of_production_history .+ _seller.utilization_cost_history
    return profit_history[trim:end]
end

function create_bool_purchase(n::Int64,k::Int64)::Vector{Bool}
    x = fill(false, n)
    x[k] = true
    return x
end

function calculate_surplus(sim_single, type::String, cumulated::Bool)

    if type == "consumer,pm"

        if cumulated

            return sum(sum(getfield.(sim_single.buyers, :realized_surplus_pm_history)) .- sum(getfield.(sim_single.buyers, :leasing_interest_history)))

        else

            return sum(getfield.(sim_single.buyers, :realized_surplus_pm_history) .- getfield.(sim_single.buyers, :leasing_interest_history))

        end

    elseif type == "consumer,sm,b"

        if cumulated

            return sum(sum(getfield.(sim_single.buyers, :realized_surplus_sm_b_history)))

        else

            return sum(getfield.(sim_single.buyers, :realized_surplus_sm_b_history))

        end


    elseif type == "consumer,sm,s"

        if cumulated

            return sum(sum(getfield.(sim_single.buyers, :realized_surplus_sm_s_history)))

        else

            return sum(getfield.(sim_single.buyers, :realized_surplus_sm_s_history))

        end

    elseif type == "consumer,total"

        if cumulated

            return sum(sum(getfield.(sim_single.buyers, :realized_surplus_pm_history))) + sum(sum(getfield.(sim_single.buyers, :realized_surplus_sm_b_history))) + sum(sum(getfield.(sim_single.buyers, :realized_surplus_sm_s_history))) - sum(sum(getfield.(sim_single.buyers, :leasing_interest_history))) + sum(sum(getfield.(sim_single.buyers, :realized_surplus_breakage_history)))

        else

            return sum(getfield.(sim_single.buyers, :realized_surplus_pm_history)) + sum(getfield.(sim_single.buyers, :realized_surplus_sm_b_history)) + sum(getfield.(sim_single.buyers, :realized_surplus_sm_s_history)) - sum(getfield.(sim_single.buyers, :leasing_interest_history)) + sum(getfield.(sim_single.buyers, :realized_surplus_breakage_history))

        end

    elseif type == "consumer,damage loss"

        if cumulated

            return sum(sum(getfield.(sim_single.buyers, :realized_surplus_breakage_history)))

        else

            return sum(getfield.(sim_single.buyers, :realized_surplus_breakage_history))


        end

    elseif type == "producer"

        if cumulated

            return sum(sum(calculate_profit_history.(sim_single.sellers)))

        else

            return sum(calculate_profit_history.(sim_single.sellers))

        end


    elseif type == "total"

        if cumulated

            return sum(sum(getfield.(sim_single.buyers, :realized_surplus_pm_history))) + sum(sum(getfield.(sim_single.buyers, :realized_surplus_sm_b_history))) + sum(sum(getfield.(sim_single.buyers, :realized_surplus_sm_s_history))) + sum(sum(calculate_profit_history.(sim_single.sellers))) - sum(sum(getfield.(sim_single.buyers, :leasing_interest_history))) + sum(sum(getfield.(sim_single.buyers, :realized_surplus_breakage_history)))

        else

            return sum(getfield.(sim_single.buyers, :realized_surplus_pm_history)) + sum(getfield.(sim_single.buyers, :realized_surplus_sm_b_history)) + sum(getfield.(sim_single.buyers, :realized_surplus_sm_s_history)) + sum(calculate_profit_history.(sim_single.sellers)) - sum(getfield.(sim_single.buyers, :leasing_interest_history)) + sum(getfield.(sim_single.buyers, :realized_surplus_breakage_history))

        end

    end

end

function cut_integer(x::Vector{Float64},k::Int64)

    bin_width = 1 / k * (maximum(x) - minimum(x)) 
    bin_up_bounds = minimum(x) .+ collect(1:k) .* bin_width
    
    x_categorical = fill(k, length(x))
    
    for i in reverse(1:(k-1))
        x_index = x .<= bin_up_bounds[i]
        x_categorical[x_index] .= i
    end

    return x_categorical, bin_up_bounds

end

function calculate_expectation(sim_res, metric, cumulated = false)

    if metric == "quality"
        if cumulated
            return mean.([getindex.(mean(getfield.(sim_res.buyers, :quality_expectation_history)),x) for x in 1:sim_res.function_args.num_sellers])
        else
            return [getindex.(mean(getfield.(sim_res.buyers, :quality_expectation_history)),x) for x in 1:sim_res.function_args.num_sellers]
        end
    elseif metric == "durability"
        if cumulated
            return mean.([getindex.(mean(getfield.(sim_res.buyers, :durability_expectation_history)),x) for x in 1:sim_res.function_args.num_sellers])
        else
            return [getindex.(mean(getfield.(sim_res.buyers, :durability_expectation_history)),x) for x in 1:sim_res.function_args.num_sellers]
        end
    end

end