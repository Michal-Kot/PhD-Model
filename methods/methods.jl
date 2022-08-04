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

function calculate_optimal_margin(my_m, my_dq, comp_p, comp_dq, my_c, δ_m)

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

function calculate_optimal_margin_and_advertising(my_m, my_dq, comp_p, comp_dq, my_c, δ_m, my_a, δ_a, my_i, comp_a, comp_i)

    my_p_current = my_c * my_dq * (1 + my_m)
    my_p_priceup = my_c * my_dq * (1 + my_m + δ_m)
    my_p_pricedown = my_c * my_dq * (1 + my_m - δ_m)

    my_a_current = my_a
    my_a_adsup = my_a + δ_a
    my_a_adsdown = my_a - δ_a
    my_a_none = 0

    s = LinRange(0:0.01:1.0)

    comp_w = [s .* cdq .* (1 + ca * ci) for (cdq, ca, ci) in zip(comp_dq, comp_a, comp_i)]

    m_a_profit = []

    for a in [my_a_none, my_a_adsdown, my_a_current, my_a_adsup]

        my_w = s .* my_dq .* (1 + a * my_i)

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

        my_π_current = my_demand_current * my_margin_current - 0.05 * 100 * a
        my_π_priceup = my_demand_priceup * my_margin_priceup - 0.05 * 100 * a
        my_π_pricedown = my_demand_pricedown * my_margin_pricedown - 0.05 * 100 * a

        push!(m_a_profit, [my_π_pricedown, my_π_current, my_π_priceup])

    end

    profit_matrix = hcat(m_a_profit...)

    my_decision = argmax(profit_matrix)

    m_decision = my_decision[1]
    a_decision = my_decision[2]

    if m_decision == 1

        new_m = 0
    
    elseif m_decision == 2

        new_m = my_m - δ_m

    elseif m_decision == 3

        new_m = my_m

    elseif m_decision == 4

        new_m = my_m + δ_m

    end

    if a_decision == 1

        new_a = my_a - δ_a

    elseif a_decision == 2

        new_a = my_a

    elseif a_decision == 3

        new_a = my_a + δ_a

    end

    return new_m, new_a

end

function seller_price_and_advertising_adjustment_observing_competitors(sellers, iter, allow_negative_margin)

    if iter == 1

        # no change to price in iter 1

        for _seller in sellers

            push!(_seller.margin_history, _seller.margin)
            push!(_seller.advertising_history, _seller.advertising)

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

                # seek for a new price and advertising

                my_margin = _seller.margin

                my_dq = _seller.quality * _seller.durability

                comp = copy(sellers)
                comp = comp[Not(_seller.id)]

                comp_p = calculate_price.(comp)

                comp_dq = getfield.(comp, :quality) .* getfield.(comp, :durability)

                my_c = _seller.average_cost

                my_a = _seller.advertising

                my_i = _seller.persuasiveness

                comp_a = getfield.(comp, :advertising)

                comp_i = getfield.(comp, :persuasiveness)

                new_margin, new_advertising = calculate_optimal_margin_and_advertising(my_margin, my_dq, comp_p, comp_dq, my_c, 0.01, my_a, 0.005, my_i, comp_a, comp_i)

                _seller.price_change = false

            else

                # verify if the last move was profitable

                profit_history = calculate_profit_history(_seller)

                profit_change = profit_history[end] - profit_history[end-1]

                if profit_change <= 0

                    # if not profitable, then return to the previous price and advertising

                    rr = rand()

                    if rr < 1/3

                        new_margin = _seller.margin_history[end-1]
                        new_advertising = _seller.advertising

                    elseif rr > 2/3

                        new_margin = _seller.margin
                        new_advertising = _seller.advertising_history[end-1]

                    else

                        new_margin = _seller.margin_history[end-1]
                        new_advertising = _seller.advertising_history[end-1]

                    end

                else

                    new_margin = _seller.margin
                    new_advertising = _seller.advertising

                end

                _seller.price_change = true

            end

            if !allow_negative_margin
                new_margin = in_boundaries(new_margin, 0, 2)
            end

            new_advertising = in_boundaries(new_advertising, 0, 1)

            _seller.margin = new_margin
            _seller.advertising = new_advertising

            push!(_seller.margin_history, _seller.margin)
            push!(_seller.advertising_history, _seller.advertising)

        end

    end

    return sellers

end

function consumer_choice(_buyer, prices, p, consumer_behaviour)

    num_sellers = length(_buyer.quality_expectation)

    wtp_advertising = 1 .+ (_buyer.ad_received .* p)

    buyer_wtp = _buyer.wealth .* _buyer.quality_expectation .* sum_of_geom_series(1, _buyer.future_discount, _buyer.durability_expectation)

    diff_p_w = prices .- buyer_wtp

    min_time_to_purchase = in_boundaries.(Float64.(floor.(diff_p_w ./ (_buyer.std_reservation_price .* _buyer.quality_expectation .* sum_of_geom_series(1, _buyer.future_discount, _buyer.durability_expectation)))), 0.0, 100000.0)

    expected_discounted_utility = (_buyer.future_discount .^ min_time_to_purchase) .* (_buyer.std_reservation_price .* _buyer.quality_expectation .* sum_of_geom_series(1, _buyer.future_discount, _buyer.durability_expectation) .* wtp_advertising .- prices)

    if consumer_behaviour == "deterministic"

        # consumer choses product that maximizes utility

        chosen_product = argmax(expected_discounted_utility)

    elseif consumer_behaviour == "stochastic"

        # consumer randomly choses products, accordingly to utility (may choose not the best one)

        weight = u2w(expected_discounted_utility)

        chosen_product = sample(1:num_sellers, Weights(weight))

    end

    chosen_time_to_purchase = min_time_to_purchase[chosen_product]

    #println((fd = _buyer.future_discount, mttp = min_time_to_purchase, fda = _buyer.future_discount .^ min_time_to_purchase,edu = expected_discounted_utility, cp = chosen_product, cttp = chosen_time_to_purchase))

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

        sellers = seller_price_and_advertising_adjustment_observing_competitors(sellers, iter, allow_negative_margin)

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

                chosen_product = consumer_choice(_buyer, prices, p, consumer_behaviour)

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

###### Load AUX

include(pwd() * "\\methods\\methods_visualization.jl")
include(pwd() * "\\methods\\methods_aux.jl")

__precompile__()