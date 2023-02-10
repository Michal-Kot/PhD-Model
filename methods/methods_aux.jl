#### AUX functions

function sum_of_geom_series_infinite(a0,q)
    return a0 ./ (1 .- q.^2)
end

function sum_of_geom_series_finite(a0,q;t)
    return a0 .* (1 .- q.^t) ./ (1 .- q)
end

function trim_first(x; trimmed = 3)
    return x[trimmed:end]
end

function trim_extremes(x, t = 10)

    return x[(x .>= percentile(x, t)) .& (x .<= percentile(x, 100-t))]

end

function in_boundaries(x::Float64,lb::Float64,ub::Float64)::Float64
    return max(min(x,ub),lb)
end

function cost_coefficient(k,d,cc)
    return cc
end

average_nan(v,d) = d == 0 ? 0 : v/d

function average_signal(signal_value, signal_volume)
    signal_average = average_nan.(signal_value, signal_volume)
    return signal_average
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


function calculate_state_profit(K::Float64, eK_dist::Vector{Float64}, D::Float64, eD_dist::Vector{Float64}, M::Float64, Q::Float64, o_K::Float64, o_D::Float64, o_P::Float64, N::Int64, μ_c::Float64, cc::Float64, ρ_dist::Vector{Float64}, β_dist::Vector{Float64}, return_type::String, product_life::Int64)
    """

    Funkcja licząca oczekiwany zysk z danego stanu. Wykorzystywana przez firmę do szacowania efektu zmiany stanu.

    """

    s = β_dist
    eK = eK_dist
    eD = eD_dist
    eρ = ρ_dist

    @assert all(0 .<= s .<= 1)
    #@assert all(0 .<= eK .<= 1)
    #@assert all(0 .<= eD .<= 1)
    @assert all(0 .<= eρ .<= 1)

    o_U = s .* sum_of_geom_series_finite(o_K, eρ * o_D; t = product_life) .- o_P # użyteczność dobra konkurencji, jeśli liczba konkurentów > 1, to o_k, o_D i o_P są średnimi

    U = s .* sum_of_geom_series_finite.(eK, eρ .* eD; t = product_life)  .- cost_coefficient(K, D, cc) .* sum_of_geom_series_infinite(K, D) .* M # użyteczność mojego dobra przy parametrach K, D, M

    demand = Int(round(mean([sum((U .> 0) .& (U .> o_U) .& (rand(N) .< 1/product_life)) for i in 1:10]))) # szacowany popyt. warunek 1: moja użyteczność > 0, warunek 2: moja użyteczność wyższa niż użyteczność dobra konkurencyjnego, warunek 3: oczekiwana liczba klientów poszukujących dobra - skalowanie dla dóbr trwałych > 1 okres

    @assert demand >= 0

    price = cost_coefficient(K, D, cc) * sum_of_geom_series_finite(K, D; t = product_life) * M  # marża na 1 sprzedanym produkcie

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

function calculate_cost(_seller::seller; product_life::Int64)::Float64
    cost = sum_of_geom_series_finite(_seller.quality, _seller.durability; t=product_life) * cost_coefficient(_seller.quality, _seller.durability, _seller.cost_coefficient) #_seller.cost_coefficient
    return cost
end


function calculate_cost_history(_seller::seller; product_life::Int64)::Vector{Float64}
    cost = sum_of_geom_series_finite.(_seller.quality_history, _seller.durability_history; t = product_life) .* cost_coefficient.(_seller.quality_history, _seller.durability_history, _seller.cost_coefficient)
    return cost
end


function calculate_price(_seller::seller; product_life::Int64)::Float64
    price = calculate_cost(_seller; product_life =  product_life) * _seller.margin
    return price
end

function calculate_price_history(_seller::seller; product_life::Int64)::Vector{Float64}
    price = calculate_cost_history(_seller; product_life =  product_life) .* _seller.margin_history
    return price
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

function calculate_profit_history(_seller::seller)::Vector{Float64}
    profit_history = _seller.selling_income_history .- _seller.cost_of_production_history .+ _seller.utilization_cost_history
    return profit_history
end

function create_bool_purchase(n::Int64,k::Int64)::Vector{Bool}
    x = fill(false, n)
    x[k] = true
    return x
end

function calculate_surplus(sim_single, type::String, cumulated::Bool)

    if type == "consumer,pm"

        if cumulated

            return sum(sum(getfield.(sim_single.buyers, :realized_surplus_pm_history)))

        else

            return sum(getfield.(sim_single.buyers, :realized_surplus_pm_history))

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

            return sum(sum(getfield.(sim_single.buyers, :realized_surplus_pm_history))) + sum(sum(getfield.(sim_single.buyers, :realized_surplus_sm_b_history))) + sum(sum(getfield.(sim_single.buyers, :realized_surplus_sm_s_history)))

        else

            return sum(getfield.(sim_single.buyers, :realized_surplus_pm_history)) + sum(getfield.(sim_single.buyers, :realized_surplus_sm_b_history)) + sum(getfield.(sim_single.buyers, :realized_surplus_sm_s_history))

        end

    elseif type == "producer"

        if cumulated

            return sum(sum(calculate_profit_history.(sim_single.sellers)))

        else

            return sum(calculate_profit_history.(sim_single.sellers))

        end


    elseif type == "total"

        if cumulated

            return sum(sum(getfield.(sim_single.buyers, :realized_surplus_pm_history))) + sum(sum(getfield.(sim_single.buyers, :realized_surplus_sm_b_history))) + sum(sum(getfield.(sim_single.buyers, :realized_surplus_sm_s_history))) + sum(sum(calculate_profit_history.(sim_single.sellers)))

        else

            return sum(getfield.(sim_single.buyers, :realized_surplus_pm_history)) + sum(getfield.(sim_single.buyers, :realized_surplus_sm_b_history)) + sum(getfield.(sim_single.buyers, :realized_surplus_sm_s_history)) + sum(calculate_profit_history.(sim_single.sellers))

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

function test_if_buyer_has_no_product(_buyer::buyer)
    x1 = _buyer.expected_surplus == 0
    x2 = all(_buyer.unit_possessed .== false) 
    x3 = _buyer.unit_possessed_is_new == false
    x4 = isnothing(_buyer.unit_possessed_time)
    x5 = isnothing(_buyer.unit_first_possessed_time)
    x6 = isnothing(_buyer.quality_of_unit_possessed_history)
    x7 = isnothing(_buyer.durability_of_unit_possessed)
    x8 = isnothing(_buyer.current_quality_of_unit_possessed)

    is_cleared = x1 & x2 & x3 & x4 & x5 & x6 & x7 & x8

    return is_cleared

end

mean_nothing(x) = length(x) == 0 ? missing : mean(x)
maximum_nothing(x) = length(x) == 0 ? missing : maximum(x)
minimum_nothing(x) = length(x) == 0 ? missing : minimum(x)
std_nothing(x) = length(x) == 0 ? missing : std(x)

function savefigs(plt, target)
    Plots.savefig(plt, pwd() * target * ".svg") # for pptx
    Plots.savefig(plt, pwd() * target * ".pdf") # for thesis, latex
end