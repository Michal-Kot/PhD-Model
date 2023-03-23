#######################################################################################################

# Auxiliary functions used in the process of  modelling

#######################################################################################################

function sum_of_geom_series_infinite(a0,q)

    """
    Calculates a sum of infinite geometric series for initial value of a0 and difference q
    Inputs:
        a0 - initial value
        q - difference between subsequent elements
    Returns:
    sum of inifinite series
    """

    return a0 ./ (1 .- q)

end

function sum_of_geom_series_finite(a0,q;t)

    """
    Calculates a sum of t elements of geometric series for initial value of a0 and difference q
    Inputs:
        a0 - initial value
        q - difference between subsequent elements
        t - number of elements
    Returns:
    sum of finite series
    """

    return a0 .* (1 .- q.^t) ./ (1 .- q)
end

function trim_first(x; trimmed = 3)

    """
    Trims first elements of a vector
    Inputs:
        x - vector
        trimmed - number of elements to be trimmed - 1
    Returns:
    vector without first trimmed-1 elements
    """

    @assert trimmed <= lastindex(x)

    return x[trimmed:end]

end

function trim_extremes(x, t = 10)

    """
    Trims first and last elements of a vector, based on percentile t
    Inputs:
        x - vector
        t - percentile to be trimmed at each extreme, 100-2t middle elements are left
    Returns:
    vector without trimmed elements

    """

    @assert t < 50

    return x[(x .>= percentile(x, t)) .& (x .<= percentile(x, 100-t))]

end

function in_boundaries(x::Float64,lb::Float64,ub::Float64)::Float64

    """
    Ensure that all elements of x are between lb (lower bound) and ub (upper bound)
    Inputs:
        x - vector
        lb - lower bound
        ub - upper bound
    Returns:
    vector in which all elements are between lb and ub, if in initial vector any element was outside feasible region (lb, ub) it is changed to the boundary value
    """

    return max(min(x,ub),lb)
end

function cost_coefficient(k,d,cc)

    """
    Calculate cost coefficient
    Inputs:
        k - quality
        d - durability
        cc - cost coefficient of production
    Returns:
    cost coefficient, if k & d are not used then production function is constant
    """

    @assert 0 < cc < 1

    return cc
end

average_nan(v,d) = d == 0 ? 0 : v/d

function average_signal(signal_value, signal_volume)

    """
    Calculate average signal from neighbours
    Inputs:
        signal_value - average of signals
        signal_volume - number of signals
    Returns:
    vector of average signals, elements being means of 0, if signal_volume = 0
    """

    signal_average = average_nan.(signal_value, signal_volume)
    return signal_average
end

function create_weights(x::Vector, method::String, scale::Bool = true)

    """
    Calculate transformation of weights
    Inputs:
        x - weights to be transformed
        method - method to transform weights: softmax, relu or argmax
        scale - if weights shall be scaled prior to transformation
    Returns:
    vector of weights transformed in line with softmax
    """

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


function calculate_state_profit(K::Float64, eK_dist::Vector{Float64}, D::Float64, eD_dist::Vector{Float64}, M::Float64, Q::Float64, o_K::Float64, o_D::Float64, o_P::Float64, N::Int64, μ_c::Float64, cc::Float64, ρ_dist::Vector{Float64}, β_dist::Vector{Float64}, ps_dist::Vector{Float64}, return_type::String, product_life::Int64)
    """
    Calculate expected profit or demand from given state, used by firms to estimate effects of strategy change
    Inputs:
        K - average quality, known to producer
        eK_dist - expected quality, based on research or assumptions <- most important difference between firms of two types
        D - average durability, known to producer
        eD_dist - expected durablity, based on research or assumptions <- most important difference between firms of two types
        M - margin
        Q - quantity of production
        o_K - quality of competitors products, averaged
        o_D - durability of competitors products, averaged
        o_P - price of competitors products, averaged
        N - population size
        μ_c - discount in the secondary B2C market
        cc - cost coefficient
        ρ_dist - future discount of consumers based on research or assumptions <- most important difference between firms of two types
        β_dist - taste for quality of consumers based on research or assumptions <- most important difference between firms of two types
        ps_dist - price sensitivity of consumers based on research or assumptions <- most important difference between firms of two types
        return_type - type of return, demand or profit
        product_life - products' ability to serve, in steps
    Returns:
    estimated profit or demand
    """

    s = β_dist
    eK = eK_dist
    eD = eD_dist
    eρ = ρ_dist

    @assert all(0 .<= s .<= 1)
    @assert all(-0.1 .<= eK .<= 1.1)
    @assert all(-0.1 .<= eD .<= 1.1)
    @assert all(0 .<= eρ .<= 1)

    o_U = s .* sum_of_geom_series_finite(o_K, eρ * o_D; t = product_life) .- o_P # użyteczność dobra konkurencji, jeśli liczba konkurentów > 1, to o_k, o_D i o_P są średnimi

    price = cost_coefficient(K, D, cc) * sum_of_geom_series_finite(K, D; t = product_life) * M  # marża na 1 sprzedanym produkcie

    U = s .* sum_of_geom_series_finite.(eK, eρ .* eD; t = product_life)  .- price # użyteczność mojego dobra przy parametrach K, D, M

    upper_price_limit = β_dist ./ ps_dist

    demand = sum((U .> 0) .& (U .> o_U) .& (rand(N) .< 1/product_life) .& (price .<= upper_price_limit)) # szacowany popyt. warunek 1: moja użyteczność > 0, warunek 2: moja użyteczność wyższa niż użyteczność dobra konkurencyjnego, warunek 3: oczekiwana liczba klientów poszukujących dobra - skalowanie dla dóbr trwałych > 1 okres

    @assert demand >= 0
    #@assert price >= 0

    profit = min(demand,Q) .* price .+ max(0, Q - demand) .* (1 - μ_c) .* cost_coefficient(K, D, cc) .* sum_of_geom_series_finite(K, D; t = product_life) - Q * cost_coefficient(K, D, cc) .* sum_of_geom_series_finite(K, D; t = product_life)  # oczekiwany zysk firmy

    @assert max(demand, Q) >= 0
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

    """
    Calculate cost of production
    Inputs:
        _seller - producer object
        product_life - products' ability to serve, in steps
    Returns:
    cost of good production for given parameters
    """

    cost = sum_of_geom_series_finite(_seller.quality, _seller.durability; t=product_life) * cost_coefficient(_seller.quality, _seller.durability, _seller.cost_coefficient) #_seller.cost_coefficient
    return cost
end


function calculate_cost_history(_seller::seller; product_life::Int64)::Vector{Float64}

    """
    Calculate history of cost of production
    Inputs:
        _seller - producer object
        product_life - products' ability to serve, in steps
    Returns:
    historical cost of good production for given parameters
    """

    cost = sum_of_geom_series_finite.(_seller.quality_history, _seller.durability_history; t = product_life) .* cost_coefficient.(_seller.quality_history, _seller.durability_history, _seller.cost_coefficient)
    return cost
end


function calculate_price(_seller::seller; product_life::Int64)::Float64

    """
    Calculate price of product
    Inputs:
        _seller - producer object
        product_life - products' ability to serve, in steps
    Returns:
    price of product for given parameters
    """

    price = calculate_cost(_seller; product_life =  product_life) * _seller.margin
    return price
end

function calculate_price_history(_seller::seller; product_life::Int64)::Vector{Float64}

    """
    Calculate history of prices of products
    Inputs:
        _seller - producer object
        product_life - products' ability to serve, in steps
    Returns:
    historical prices of goods for given parameters
    """

    price = calculate_cost_history(_seller; product_life =  product_life) .* _seller.margin_history
    return price
end

function calculate_profit_history(_seller::seller)::Vector{Float64}

    """
    Calculate history of profit
    Inputs:
        _seller - producer object
    Returns:
    historical profit if producer
    """

    profit_history = _seller.selling_income_history .- _seller.cost_of_production_history .+ _seller.utilization_cost_history
    return profit_history
end

function create_bool_purchase(n::Int64,k::Int64)::Vector{Bool}

    """
    Present purchase decision as vector of bool instead of product's ID = k
    Inputs:
        n - number of available products
        k - chosen product
    Returns:
    vector of bool, true if product is chosen, false otherwise
    """

    x = fill(false, n)
    x[k] = true

    @assert sum(x) == 1

    return x
end

function calculate_surplus(sim_single, type::String, cumulated::Bool)

    """
    Calculate history of surpluses
    Inputs:
        sim_single - simulation results
        type - type of surplus to be calculated
        cumulated - if cumulated surplus be returned, or vector of surpluses in each step
    Returns:
    surplus of given type, cumulated or not
    """

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

    """
    Cut numeric vector into k groups
    Inputs:
        x - vector to be cut
        k - number of groups
    Returns:
    group assignment, boundaries between groups
    """

    bin_width = 1 / k * (maximum(x) - minimum(x)) 
    bin_up_bounds = minimum(x) .+ collect(1:k) .* bin_width
    
    x_categorical = fill(k, length(x))
    
    for i in reverse(1:(k-1))
        x_index = x .<= bin_up_bounds[i]
        x_categorical[x_index] .= i
    end

    return x_categorical, bin_up_bounds

end

function initial_m(β, k1, d1, k2, d2, H, c, m2)

    """
    Calculate initial m to ensure that for customer with β expected utility of both goods is equal
    Inputs:
        β - taste for quality
        k1, d1 - quality and durability of good #1
        k2, d2 - quality and durability of good #2
        H - products' ability to serve, in steps
        c - cost coefficient
        m2 - margin of the competitor's product
    Returns:
    margin of firm #1
    """
    
    m1= (β*sum_of_geom_series_finite(k1,d1;t=H)-sum_of_geom_series_finite(k2,d2;t=H) * (β - m2 * c)) / (c * sum_of_geom_series_finite(k1,d1;t=H))
    return m1
end

function transition_percentile(x)

    """
    calculate probabiity of loss
    Inputs:
       x - vector of profits
    Returns:
    percentile in which profit becomes positive
    """
    
    trans_perc = 0.0
    for i in 100.0:-0.01:0.01
        if (percentile(x, i) >= 0) & (percentile(x, i - 0.01) < 0)
            trans_perc = i
        end
    end
    return trans_perc
end


function find_nash_eq_pure(payoff_matrix::Matrix)

    """
    find pure nash equilibria in (2,2) Matrix
    Inputs:
        payoff_matrix - matrix of payoffs, each element being a tuple (payoff_player_1, payoff_player_2)
    Returns:
    bool matrix with true if given strategy is NE, false otherwise
    """

    # player 1
    is_nash_eq_1 = zeros(Int64, size(payoff_matrix))
    for col in 1:size(payoff_matrix, 2)
        max_val = maximum(getindex.(payoff_matrix[:,col], 1))
        for row in 1:size(payoff_matrix, 1)
            if payoff_matrix[row, col][1] == max_val
                is_nash_eq_1[row, col] = 1
            end
        end
    end
    # player 2
    is_nash_eq_2 = zeros(Int64, size(payoff_matrix))
    for row in 1:size(payoff_matrix, 1)
        max_val = maximum(getindex.(payoff_matrix[row,:], 2))
        for col in 1:size(payoff_matrix, 2)
            if payoff_matrix[row, col][2] == max_val
                is_nash_eq_2[row, col] = 1
            end
        end
    end
    is_nash_eq = is_nash_eq_1 .* is_nash_eq_2
    return is_nash_eq
end

function find_nash_eq_mixed(payoff_matrix::Matrix)

    """
    find mixed nash equilibria in (2,2) Matrix
    Inputs:
        payoff_matrix - matrix of payoffs, each element being a tuple (payoff_player_1, payoff_player_2)
    Returns:
    (p,q) - probabilities of playing strategy one by each of players
    """

    # player 1
    payoffs1 = getindex.(payoff_matrix, 1)
    q = (payoffs1[2,2] - payoffs1[1,2]) / (payoffs1[1,1] + payoffs1[2,2] - payoffs1[1,2] - payoffs1[2,1])
    payoffs2 = getindex.(payoff_matrix, 2)
    p = (payoffs2[2,2] - payoffs2[2,1]) / (payoffs2[1,1] + payoffs2[2,2] - payoffs2[1,2] - payoffs2[2,1])
    #if (p >= 0) & (p <= 1) & (q >= 0) & (q <= 1)
        return (p,q)
    #else
    #    return missing
    #end 
end

function simulate_ne_for_costs(payoff1, payoff2, costs, ret)

    """
    find all equilibria for given payoffs and cost of research
    Inputs:
        payoff1 - payoffs of player 1
        payoff2 - payoffs of player 2
        costs - cost of research
        ret - type of return
    Returns:
    pure / mixed equlibria or check if prisoners' dilemma
    """

    equilibriums = []

    for c in costs
        payoff1_cost = payoff1 .- [0, c, 0, c]
        payoff2_cost = payoff2 .- [0, 0, c, c]
        PM = construct_payoff_matrix(payoff1_cost, payoff2_cost, (2,2))
        isPD = check_if_prisoners_dilemma(PM)
        pNE = find_nash_eq_pure(PM)
        mNE = find_nash_eq_mixed(PM)
        push!(equilibriums, (PM, pNE, mNE, isPD))
    end

    nash_equilibriums_pure = getindex.(equilibriums, 2)
    positions = [1 3
    2 4]

    nash_equilibriums_pure = [ne .* positions for ne in nash_equilibriums_pure]

    if ret == "pNE"

        return nash_equilibriums_pure

    elseif ret == "mNE"

        return getindex.(equilibriums, 3)

    elseif ret == "PD"

        return getindex.(equilibriums, 4)

    end

end

function construct_payoff_matrix(payoff1, payoff2, shape = (2,2))

    """
    construct payoff matrix for given payoffs
    Inputs:
        payoff1 - payoffs of player 1
        payoff2 - payoffs of player 2
        shape - shape of payoff matrix
    Returns:
    payoff matrix
    """

    payoff_matrix = [(x,y) for (x,y) in zip(payoff1, payoff2)]
    payoff_matrix = reshape(payoff_matrix, shape)
    return payoff_matrix
end

function check_if_prisoners_dilemma(pm, detailed = false)

    """
    check if given payoffs implicate problem being a prisoner's dilemma
    Inputs:
        pm - payoff matrix
    Returns:
    if prisoner's dilemma
    """

# T>R>P>S: betray > cooperate > punishment > fool

    pm1 = getindex.(pm, 1)

    TR1 = pm1[1,2] > pm1[2,2]
    RP1 = pm1[2,2] > pm1[1,1]
    PS1 = pm1[1,1] > pm1[2,1]

    pd1 = TR1 & RP1 & PS1

    pm2 = getindex.(pm, 2)

    TR2 = pm2[2,1] > pm2[2,2]
    RP2 = pm2[2,2] > pm2[1,1]
    PS2 = pm2[1,1] > pm2[1,2]

    pd2 = TR2 & RP2 & PS2

    if detailed
        return (pd1&pd2, pd1, (TR1, RP1, PS1), pd2, (TR2, RP2, PS2))
    else
        return pd1&pd2
    end
end

function get_expectation_buyers(buyers, metric; s, T)

    """
    calculate expected features ofor buyers only
    Inputs:
        buyers - population of consumers
        metric - feature of product to calculate expectation upon
        s - number of sellers
        T - number of periods
    Returns:
    vector of expected features
    """

    buyers_expectations = [[getindex.([multi_with_missing.(x,y) for (x,y) in zip(getfield(b, :received_signal_history), getfield(b, metric))], x) for x in 1:s] for b in buyers]

    return [mean_na.([getindex.(getindex.(buyers_expectations, s), t) for t in 1:T]) for s in 1:2]

end

######### Additional methods ###################

true_or_missing(x; k) = x == 1 ? k : missing
fTC(Q, K, D, H, c) = Q .* c .* K .* (1 .- D .^ H) ./ (1 .- D)
fATC(K, D, H, c) = c .* K .* (1 .- D .^ H) ./ (1 .- D)
any_vec(x, s) = any.(getindex.(x,s) .> 0)
multi_with_missing(x,y) = x == 0 ? missing : x * y
mean_na(x) = all(ismissing.(x)) ? missing : mean(x[.!ismissing.(x)])
RMSE(x,y) = sqrt(mean((x .- y).^2))
xydiff(x,y) = mean(x .- y)
cv(x) = std(x) / mean(x)
multi(x,y) = x .* y
divide(x,y) = x ./ y
mean_nothing(x) = length(x) == 0 ? missing : mean(x)
maximum_nothing(x) = length(x) == 0 ? missing : maximum(x)
minimum_nothing(x) = length(x) == 0 ? missing : minimum(x)
std_nothing(x) = length(x) == 0 ? missing : std(x)
sum_na(x) = lastindex(x) == 0 ? 0 : sum(x)