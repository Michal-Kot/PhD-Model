using Distributions, StatsBase
using OrderedCollections
using Plots
using StatsPlots
using LinearAlgebra
using InvertedIndices
using Graphs
using GraphPlot, Compose, Cairo, Fontconfig
using HypothesisTests
using SmoothingSplines
using GLM
using DataFrames
using SmoothingSplines
using Random

### TO DO: add duration impact on buyers decisions


################################# STRUCTS ####################################################################

"""

Sellers form a supply side of the model

Sellers' main goal is to maximize profit π
π = pq - c(q)p
p - price of product, q - quantity of sales, c - function of cost, linear on q, a - advertising intensity

"""

mutable struct seller
    id::Int64
    quality::Float64
    durability::Float64
    average_cost::Float64
    margin::Float64
    quantity_history::Vector{Float64}
    margin_history::Vector{Float64}
    advertising::Float64
    advertising_history::Vector{Float64}
    risk_decisions::Float64
    persuasiveness::Float64
    price_change::Bool
end

function create_sellers(num_sellers::Int64,q::Vector{Float64},c::Vector{Float64},m::Vector{Float64},a::Vector{Float64},r::Vector{Float64},p::Vector{Float64},d::Vector{Int64})::Vector{seller}
    @assert length(q) == num_sellers
    @assert length(c) == num_sellers
    @assert length(m) == num_sellers
    @assert length(a) == num_sellers
    @assert length(r) == num_sellers
    @assert length(p) == num_sellers
    sellers_vector = []
    for s in 1:num_sellers
        new_seller = seller(s, q[s], d[s], c[s], m[s], [], [], a[s], [], r[s], p[s], true)
        push!(sellers_vector, new_seller)
    end
    return sellers_vector
end

"""

Buyers form a demand side of the model

Buyers main goal is to maximize utility of consumption υ(q) = s*eq - p, subject to budget constraint ρ*eq >= p
s - seeking for quality parameter, eq = expected quality of the product, p - price of product, ρ - willingness to pay for unit of quality

"""

mutable struct buyer
    id::Int64
    neighbours::Vector{Int64}
    std_reservation_price::Float64
    future_discount::Float64
    quality_seeking::Float64
    wealth::Float64
    quality_expectation::Vector{Float64}
    quality_expectation_history::Vector{Vector{Float64}}
    durability_expectation::Vector{Float64}
    durability_expectation_history::Vector{Vector{Float64}}
    unit_possessed::Vector{Bool}
    unit_possessed_history::Vector{Vector{Bool}}
    unit_possessed_time::Int64
    unit_possessed_durability::Int64
    quality_of_unit_possessed::Float64
    quality_of_unit_possessed_history::Vector{Vector{Float64}}
    ad_received::Vector{Bool}
    ad_received_history::Vector{Vector{Bool}}
    surplus_history::Vector{Float64}
end

function create_buyers(num_buyers::Int64, num_sellers::Int64, network::SimpleGraph{Int64}, initial_quality_expectation::Vector{Float64}=fill(1.0, num_sellers), initial_duration_expectation::Vector{Float64}=fill(5.0, num_sellers))::Vector{buyer}
    buyers_vector = []
    for b in 1:num_buyers
        my_neighbours = neighbors(network, b)
        qs_srp = rand(Uniform(0,1))
        new_buyer = buyer(b, #id
            my_neighbours, #neighbours
            qs_srp, #std_reservation_price
            rand(Uniform(0.5,0.8)), #future_discount
            qs_srp, #quality_seeking
            qs_srp, #wealth
            initial_quality_expectation, #quality_expectation
            [], #quality_expectation_history
            initial_duration_expectation, #duration_expectation
            [], #duration_expectation_history
            fill(false, num_sellers), #unit_possessed
            [], #unit_possessed_history,
            0, #unit_possessed_time
            0, #unit_possessed_durability
            0.0, #quality_of_unit_possessed
            [], #quality_of_unit_possessed_history
            [], #ad_received
            [], #ad_received_history
            []) #surplus_history
        push!(buyers_vector, new_buyer)
    end
    return buyers_vector
end

"""

Create network of a given type - vertices are buyers

if random type, then use Erdos Renyi algorithm, specify number of vertices and edges
if preferential attachment type , then use Albert Barabasi algorithm, specify number of vertices and number of links for a single vertex (constant)

"""

function create_network(type::String; num_buyers::Int64, num_links::Int64=0, pref_attachment_links::Int64=0)

    if type == "random"
        return Graphs.erdos_renyi(num_buyers, num_links)
    elseif type == "preferential_attachment"
        return Graphs.barabasi_albert(num_buyers, pref_attachment_links)
    end
end

#################################### AUX FUNCTIONS ##############################################################

function in_boundaries(x::Float64,lb::Float64,ub::Float64)::Float64
    return max(min(x,ub),lb)
end

function calculate_cost(_seller::seller)::Float64
    cost = _seller.durability * _seller.quality * _seller.average_cost
    return cost
end

function calculate_price(_seller::seller)::Float64
    price = calculate_cost(_seller) * (1 + _seller.margin)
    return price
end

function calculate_price_history(_seller::seller)::Vector{Float64}
    price = calculate_cost(_seller) .* (1 .+ _seller.margin_history)
    return price
end

function u2w(u::Vector{Float64}, p_min::Float64 = 0.1)::Vector{Float64}
    @assert p_min < 1.0
    ut = u .- minimum(u)
    wgt = ut ./ sum(ut) .* (1 - p_min)
    wgt[argmin(u)] = p_min
    return wgt
end

function calculate_profit_history(_seller::seller)::Vector{Float64}
    (calculate_price_history(_seller) .- calculate_cost(_seller)) .* _seller.quantity_history
end

function plot_phasediagram(_seller::seller, function_args::Vector)
    x = calculate_price_history(_seller)
    y = calculate_profit_history(_seller, function_args.γ, function_args.num_buyers)
    d_x = diff(x)
    push!(d_x, 0)
    d_y = diff(y)
    push!(d_y,0)
    labs = string.(1:length(x))
    points_labels = ([mod(x,100) == 0 ? labs[x] : "" for x in 1:length(labs)])
    points_labels[1] = "1"
    quiver(x,y,quiver=(d_x, d_y), color = :green)
    scatter!(x,y, xlabel = "Price", ylabel = "Profit", markershape = :none, markercolor = :white, markerstrokecolor = :white, series_annotations = points_labels, markersize = 0)
end

function add_smoothing_spline(x,y,clr,λ=0.05)
    spl = fit(SmoothingSpline, x, y, λ)
    y_hat = predict(spl)
    plot!(sort(x), y_hat[sortperm(x)], label = "", color = clr)
end

function calculate_average_elasticity(q,p;trim = 10)
    q = q[(trim+1):end]
    p = p[(trim+1):end]
    dq = diff(q) ./ mean(q)
    dp = diff(p) ./ mean(p)
    eqp = dq ./ dp
    eqp = eqp[eqp .!= Inf]
    eqp = eqp[eqp .!= -Inf]
    eqp = eqp[.!isnan.(eqp)]
    return mean(eqp)
end

function trim_outliers(x, trim = 5)
    lb = percentile(x, trim)
    ub = percentile(x, 100-trim)
    y = copy(x)
    y = y[(y .>= lb) .& (y .<= ub)]
    return y
end

function create_bool_purchase(n::Int64,k::Int64)::Vector{Bool}
    x = fill(false, n)
    x[k] = true
    return x
end

function calculate_total_surplus(sim_res, return_type::String = "total", cumulated::Bool = true)

    producer_surplus = sum(calculate_profit_history.(sim_res.sellers, sim_res.function_args.γ, sim_res.function_args.num_buyers))
    consumer_surplus = sum(getfield.(sim_res.buyers, :surplus_history))
    total_surplus = producer_surplus .+ consumer_surplus

    if return_type == "total"

        if cumulated 

            return sum(total_surplus)

        else

            return total_surplus

        end

    elseif return_type == "producer"

        if cumulated 

            return sum(producer_surplus)

        else

            return producer_surplus

        end

    elseif return_type == "consumer"

        if cumulated 

            return sum(consumer_surplus)

        else

            return consumer_surplus

        end

    end

end

function calculate_profit(my_m, my_dq, comp_p, comp_dq, my_c, δ_m)

    my_p_current = my_c * my_dq * (1 + my_m)
    my_p_priceup = my_c * my_dq * (1 + my_m + δ_m)
    my_p_pricedown = my_c * my_dq * (1 + my_m - δ_m)

    s = LinRange(0:0.01:1.0)

    my_w = s .* my_dq
    comp_w = [s .* cdq for cdq in comp_dq]

    my_spl_current = my_w .- my_p_current
    my_spl_priceup = my_w .- my_p_priceup
    my_spl_pricedown = my_w .- my_p_pricedown

    comp_spl = [cw .- cp for (cw,cp) in zip(comp_w, comp_p)]

    my_spl_diff_current = reduce(.*,[my_spl_current .> cspl for cspl in comp_spl])
    my_spl_diff_priceup = reduce(.*,[my_spl_priceup .> cspl for cspl in comp_spl])
    my_spl_diff_pricedown = reduce(.*,[my_spl_pricedown .> cspl for cspl in comp_spl])

    my_demand_current = count((my_spl_current .>= 0) .& (my_spl_diff_current))
    my_demand_priceup = count((my_spl_priceup .>= 0) .& (my_spl_diff_priceup))
    my_demand_pricedown = count((my_spl_pricedown .>= 0) .& (my_spl_diff_pricedown))

    my_margin_current = my_p_current - my_dq * my_c
    my_margin_priceup = my_p_priceup - my_dq * my_c
    my_margin_pricedown = my_p_pricedown - my_dq * my_c

    my_π_current = my_demand_current * my_margin_current
    my_π_priceup = my_demand_priceup * my_margin_priceup
    my_π_pricedown = my_demand_pricedown * my_margin_pricedown

    my_decision = argmax([my_π_pricedown, my_π_current, my_π_priceup])

    if my_decision == 1

        return my_m - δ_m

    elseif my_decision == 2

        return my_m

    elseif my_decision == 3

        return my_m + δ_m

    end

end

function seller_price_adjustment_observing_competitors(sellers, iter, allow_negative_margin)

    if iter == 1

        # no change to price in iter 1

        for _seller in sellers

            push!(_seller.margin_history, _seller.margin)
            #push!(_seller.advertising_history, _seller.advertising)

        end

    elseif iter >= 2

        # change based on P(price) and K(quality-duration)
            # if P/K > P̂/K̂ -> test lower price, to increase profit
            # if P/K = P̂/K̂ -> test new, random price, to increase profit
            # if P/K < P̂/K̂ -> test higher price, to increase profit
        # price changes by δ
        # advertising changes by α

        for _seller in sellers

            if _seller.price_change == true

                # seek for a new price

                my_margin = _seller.margin

                my_dq = _seller.quality * _seller.durability

                comp = copy(sellers)
                comp = comp[Not(_seller.id)]

                comp_p = calculate_price.(comp)

                comp_dq = getfield.(comp, :quality) .* getfield.(comp, :durability)

                my_c = _seller.average_cost

                my_r = _seller.risk_decisions

                new_margin = calculate_profit(my_margin, my_dq, comp_p, comp_dq, my_c, my_r)

                _seller.price_change = false

            else

                # verify if last move was profitable

                profit_history = calculate_profit_history(_seller)

                profit_change = profit_history[end] - profit_history[end-1]

                if profit_change <= 0

                    # if not profitable, then return to the previous price

                    new_margin = _seller.margin_history[end-1]

                else

                    new_margin = _seller.margin

                end

                _seller.price_change = true

            end

            if !allow_negative_margin
                new_margin = in_boundaries(new_margin, 0, 2)
            end

            _seller.margin = new_margin
            push!(_seller.margin_history, _seller.margin)

            push!(_seller.advertising_history, _seller.advertising)

        end

    end

    return sellers

end

function cut_integer(x::Vector{Float64},k::Int64)::Vector{Int64}

    bin_width = 1 / k * (maximum(x) - minimum(x)) 
    bin_up_bounds = collect(1:k) .* bin_width
    
    x_categorical = fill(k, length(x))
    
    for i in reverse(1:(k-1))
        x_index = x .<= bin_up_bounds[i]
        x_categorical[x_index] .= i
    end

    return x_categorical

end

function plot_quantity(sim_res)
    groupedbar(reduce(hcat,[getfield.(sim_res.sellers, :quantity_history)[x] for x in 1:length(sim_res.sellers)]),linecolor=nothing, bar_position = :stack, xlabel = "Time", ylabel = "Quantity of sales", label = reshape("Producer " .* string.(1:sim_res.function_args.num_sellers), 1, sim_res.function_args.num_sellers))
end

function plot_quality_expectation(sim_res)
    plot([getindex.(mean(getfield.(sim_res.buyers, :quality_expectation_history)),x) for x in 1:sim_res.function_args.num_sellers], label = reshape("Producer " .* string.(1:sim_res.function_args.num_sellers), 1, sim_res.function_args.num_sellers))
end

function plot_profit_history(sim_res)
    plot(calculate_profit_history.(sim_res.sellers))
end

function plot_advertising(sim_res)
    plot(getfield.(sim_res.sellers, :advertising_history), xlabel = "Time", ylabel = "Advertising intensity [share of population targeted]", label = reshape("Producer " .* string.(1:sim_res.function_args.num_sellers), 1, sim_res.function_args.num_sellers))
end

function plot_margin(sim_res)
    Plots.plot([getfield.(sim_res.sellers, :margin_history)[x] for x in 1:length(sim_res.sellers)], xlabel = "Time", ylabel = "Margin", label = reshape("Producer " .* string.(1:sim_res.function_args.num_sellers), 1, sim_res.function_args.num_sellers))
end

function plot_ecdf(metric, label, xlabel, ylabel, title, is_new)
    if is_new
        Plots.plot(sort(metric), (1:length(metric))./length(metric), xlabel = xlabel, ylabel = ylabel, title = title, label = label, legend=:bottomright)
    else
        Plots.plot!(sort(metric), (1:length(metric))./length(metric), xlabel = xlabel, ylabel = ylabel, title = title, label = label, legend=:bottomright)
    end
end

function expected_duration(d_max)
    pi = 0.5 .* collect(1:(2*d_max)) / d_max
    return sum([prod(1 .- (pi[1:(i-1)]))*pi[i]*i for i in 1:length(pi)])
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

function sum_of_geom_series(a1,q,n)
    a1 .* (1 .- q .^ n) ./ (1 .- q)
end

function consumer_choice(_buyer, prices, p, consumer_behaviour, iter)

    wtp_advertising = 1 .+ (_buyer.ad_received .* p)

    buyer_wtp = _buyer.wealth .* _buyer.quality_expectation .* sum_of_geom_series(1, _buyer.future_discount, _buyer.durability_expectation) .* wtp_advertising

    diff_p_w = prices .- buyer_wtp

    min_time_to_purchase = in_boundaries.(Float64.(floor.(diff_p_w ./ (_buyer.std_reservation_price .* _buyer.quality_expectation .* sum_of_geom_series(1, _buyer.future_discount, _buyer.durability_expectation)))), 0.0, 100000.0)

    expected_discounted_utility = _buyer.future_discount .^ min_time_to_purchase .* _buyer.std_reservation_price .* _buyer.quality_expectation .* sum_of_geom_series(1, _buyer.future_discount, _buyer.durability_expectation)

    if consumer_behaviour == "deterministic"

        # consumer choses product that maximizes utility

        chosen_product = argmax(expected_discounted_utility)

    elseif consumer_behaviour == "stochastic"

        # consumer randomly choses products, accordingly to utility (may choose not the best one)

        weight = u2w(expected_discounted_utility)

        chosen_product = sample(1:num_sellers, Weights(weight))

    end

    chosen_time_to_purchase = min_time_to_purchase[chosen_product]

    if chosen_time_to_purchase == 0

        # buy now

        return chosen_product

    else 

        # postpone purchase

        return nothing

    end

end

############################################ MAIN SIM FUNCTION #####################################################

function TO_GO(num_sellers::Int64, num_buyers::Int64, max_iter::Int64, λ_ind::Float64, λ_wom::Float64, consumer_behaviour::String; q::Vector{Float64} = fill(1.0, num_sellers), c::Vector{Float64} = fill(0.8, num_sellers), m::Vector{Float64} = fill(0.2, num_sellers), a::Vector{Float64} = fill(0.05, num_sellers), r::Vector{Float64} = fill(0.05, num_sellers), ϵ::Vector{Float64} = fill(1/3, num_sellers), q_init::Vector{Float64} = fill(1.0, num_sellers), p::Vector{Float64} = fill(0.05, num_sellers), d::Vector{Int64} = fill(5, num_sellers), num_links::Int64 = 200, δ::Float64 = 0.05, γ::Float64 = 0.050, α::Float64 = 0.005, variant_advertising::Bool = true, allow_negative_margin::Bool = true)


    function_args = (num_sellers = num_sellers, num_buyers = num_buyers, max_iter = max_iter, λ_ind = λ_ind, λ_wom = λ_wom, q = q, c = c, m = m, a = a, r = r, ϵ = ϵ, q_init = q_init, p=p, d=d, num_links = num_links, δ = δ, γ = γ, α=α)

    # create sellers and buyers

    sellers = create_sellers(num_sellers, q, c, m, a, r, p, d)

    buyers_network = create_network("random", num_buyers = num_buyers, num_links = num_links)
    buyers = create_buyers(num_buyers, num_sellers, buyers_network, q_init, Float64.(d))

    # utlity function definition

    u(k,d,p;β) = ((0.5 + β) .* k .* d) .- p

    # TEMP 🍎

    utility_history = []
    durability_history = []

    # sellers are choosing new price

    #margin_lb = -2
    #margin_ub = 2

    for iter in 1:max_iter

        sellers = seller_price_adjustment_observing_competitors(sellers, iter, allow_negative_margin)

        #### RYNEK ####

        demand = Int64[]

        # REKLAMY

        advertising_intensity = getfield.(sellers, :advertising)

        for _buyer in buyers
            # Advertising

            if rand() < sum(advertising_intensity)
                firm_ad_seen = sample(1:num_sellers, Weights(advertising_intensity))
                advertising_bool = fill(false, num_sellers)
                advertising_bool[firm_ad_seen] = true
                _buyer.ad_received = advertising_bool
                push!(_buyer.ad_received_history, _buyer.ad_received)
            else
                _buyer.ad_received = fill(false, num_sellers)
                push!(_buyer.ad_received_history, _buyer.ad_received)
            end
            
        end 

        for _buyer in buyers

            # Product usage -> if broken

            if iter >= 2

                if any(_buyer.unit_possessed)

                    product_amortization = iter - _buyer.unit_possessed_time
                    product_durability = _buyer.unit_possessed_durability

                    if product_amortization >= product_durability

                        #_buyer.broken_product = true
                        _buyer.unit_possessed = fill(false, num_sellers)

                        push!(durability_history, (argmax(_buyer.unit_possessed), product_amortization))

                        _buyer.durability_expectation = _buyer.durability_expectation .+ λ_ind .* _buyer.unit_possessed .* (product_amortization .- _buyer.durability_expectation)

                        for idn in _buyer.neighbours
                            buyers[idn].durability_expectation = buyers[idn].durability_expectation .+ λ_wom .* _buyer.unit_possessed .* (product_amortization .- buyers[idn].durability_expectation)
                        end
                    end

                end

            end

        end

        prices = calculate_price.(sellers)

        for _buyer in buyers

            # buying new product in iter == 1, or replacing product in iter >= 2

            if all(_buyer.unit_possessed .== false)

                chosen_product = consumer_choice(_buyer, prices, p, consumer_behaviour, iter)

                if !isnothing(chosen_product)

                    _buyer.unit_possessed = create_bool_purchase(num_sellers, chosen_product)
                    _buyer.unit_possessed_time = iter
                    _buyer.unit_possessed_durability = rand(Poisson(sellers[chosen_product].durability))
                    push!(_buyer.unit_possessed_history, _buyer.unit_possessed)
                    push!(demand, chosen_product)
                    _buyer.wealth = 0    

                else

                    _buyer.unit_possessed = fill(false, num_sellers)
                    _buyer.unit_possessed_durability = 0
                    push!(_buyer.unit_possessed_history, _buyer.unit_possessed)

                end

            else

                push!(_buyer.unit_possessed_history, _buyer.unit_possessed)

            end

            _buyer.wealth += _buyer.std_reservation_price

        end

        for _buyer in buyers

            # calculate quality driven by the product

            if any(_buyer.unit_possessed) # if product is possessed

                chosen_product = argmax(_buyer.unit_possessed)

                _buyer.quality_of_unit_possessed = rand(Uniform(sellers[chosen_product].quality - ϵ[chosen_product], sellers[chosen_product].quality + ϵ[chosen_product]))

                push!(_buyer.quality_of_unit_possessed_history, _buyer.unit_possessed .*  _buyer.quality_of_unit_possessed)

                _buyer.quality_expectation = _buyer.quality_expectation .+ λ_ind .* _buyer.unit_possessed .* (_buyer.quality_of_unit_possessed .- _buyer.quality_expectation)

                for idn in _buyer.neighbours
                    buyers[idn].quality_expectation = buyers[idn].quality_expectation .+ λ_wom .* _buyer.unit_possessed .* (_buyer.quality_of_unit_possessed .- buyers[idn].quality_expectation)
                end

                if iter == _buyer.unit_possessed_time

                    push!(_buyer.surplus_history, _buyer.std_reservation_price * _buyer.quality_of_unit_possessed - prices[chosen_product])

                else

                    push!(_buyer.surplus_history, _buyer.std_reservation_price * _buyer.quality_of_unit_possessed)

                end
            
            else

                _buyer.quality_of_unit_possessed = 0.0
                push!(_buyer.quality_of_unit_possessed_history, fill(0.0, num_sellers))
                push!(_buyer.surplus_history, 0.0)

            end

        end

        # kupujący weryfikują oczekiwaną jakość produktów

        for _buyer in buyers

            push!(_buyer.quality_expectation_history, _buyer.quality_expectation)
            push!(_buyer.durability_expectation_history, _buyer.durability_expectation)

        end

        # firma weryfikuje Q sprzedaży

        for _seller in sellers

            push!(_seller.quantity_history, count(demand .== _seller.id))
        end

    end

    return (sellers = sellers, buyers = buyers, function_args = function_args, ut_his = utility_history, dur_his = durability_history)

end

__precompile__()