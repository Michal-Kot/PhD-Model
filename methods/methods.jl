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
using Serialization
using PlotlyJS
using CSV

################################# STRUCTS ####################################################################

"""

Sellers form a supply side of the model

Sellers' main goal is to maximize profit π
π = pq - c(q)p
p - price of product, q - quantity of sales, c - function of cost, linear on q

"""

mutable struct seller
    id::Int64
    quality::Float64
    durability::Float64

    cost_coefficient::Float64
    margin::Float64

    cost_of_production::Float64
    cost_of_production_history::Vector{Float64}
    
    quantity_produced::Float64
    quantity_sold::Float64

    quantity_produced_history::Vector{Float64}
    quantity_sold_history::Vector{Float64}

    quality_history::Vector{Float64}
    durability_history::Vector{Float64}
    margin_history::Vector{Float64}

    reselling_history::Vector{Int64}

    selling_income::Float64

    selling_income_history::Vector{Float64}

    expected_income::Float64
    expected_income_history::Vector{Float64}

    quality_range::Vector{Float64}
    durability_range::Vector{Float64}
    margin_range::Vector{Float64}

    utilization_cost_history::Vector{Float64}

    consumer_research::Bool
    sample_size::Float64

    consumers_asked::Vector{Vector{Int64}}

end

function create_sellers(num_sellers::Int64,c::Vector{Float64},m::Vector{Float64}, kr::Vector{Vector{Float64}}, dr::Vector{Vector{Float64}}, mr::Vector{Vector{Float64}}, cr::Vector{Bool}, ss::Vector{Float64})::Vector{seller}
    
    @assert length(c) == num_sellers
    @assert length(m) == num_sellers
    @assert length(kr) == num_sellers
    @assert length(dr) == num_sellers
    @assert length(mr) == num_sellers

    sellers_vector = []

    for s in 1:num_sellers

        k = rand(Uniform(kr[s][1], kr[s][2]))
        d = rand(Uniform(dr[s][1], dr[s][2]))

        new_seller = seller(s, #id
        k, #quality
        d, #durability
        c[s], #cost_coefficient
        m[s], #margin
        0.0, #cost_of_production
        [], #cost_of_production_history
        0.0, #quantity produced
        0.0, #quantity_sold
        [], #quantity_produced_history
        [], #quantity_sold_history
        [], #quality_history
        [], #durability_history
        [], # margin_history
        [], #reselling_history
        0.0,
        [],
        0.0,
        [],
        kr[s], #quality_range
        dr[s], #durability_range
        mr[s], #margin_range
        [],
        cr[s],
        ss[s],
        []) #utilization_cost_history
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
    quality_expectation::Vector{Float64}
    quality_expectation_history::Vector{Vector{Float64}}
    durability_expectation::Vector{Float64}
    durability_expectation_history::Vector{Vector{Float64}}
    expected_surplus::Float64
    unit_possessed::Vector{Bool}
    unit_possessed_is_new::Bool
    unit_buying_selling_history::Vector{}
    unit_possessed_time::Int64
    unit_first_possessed_time::Int64
    quality_of_unit_possessed::Float64
    quality_of_unit_possessed_history::Vector{Vector{Float64}}
    durability_of_unit_possessed::Float64
    durability_of_unit_possessed_history::Vector{Vector{Float64}}
    current_quality_of_unit_possessed::Float64
    unit_for_sale::Vector{Bool}
    quality_for_sale::Float64
    durability_for_sale::Float64
    price_for_sale::Float64
    age_for_sale::Int64
    realized_surplus_pm_history::Vector{Float64}
    realized_surplus_sm_b_history::Vector{Float64}
    realized_surplus_sm_s_history::Vector{Float64}
    choosing_surplus::Vector{Vector{Float64}}
    reselling_probability::Float64
    signal_volume::Int64
end

function create_buyers(num_buyers::Int64, sellers::Vector{seller}, num_sellers::Int64, network::SimpleGraph{Int64}, future_discount_range::Vector{Float64}, product_life::Int64)::Vector{buyer}
    buyers_vector = []
    for b in 1:num_buyers
        my_neighbours = neighbors(network, b)

        qs_srp = rand(Uniform(0,1))
        future_discount = rand(Uniform(future_discount_range[1],future_discount_range[2]))

        initial_quality_expectation = getfield.(sellers, :quality)
        initial_durability_expectation = getfield.(sellers, :durability)

        # initial random assignment of product #

        probability_to_posess_product = qs_srp .> (getfield.(sellers, :cost_coefficient) .* getfield.(sellers, :margin))

        if any(probability_to_posess_product .> 0)

            possessed_product = sample(getfield.(sellers, :id), Weights(probability_to_posess_product))

            unit_possessed = fill(false, num_sellers)
            unit_possessed[possessed_product] = true

            unit_possessed_time = -1 * sample(1:product_life)
            unit_first_possessed_time = unit_possessed_time
            quality_of_unit_possessed = initial_quality_expectation[possessed_product]
            durability_of_unit_possessed = initial_durability_expectation[possessed_product]
            current_quality_of_unit_possessed = quality_of_unit_possessed * durability_of_unit_possessed ^ (-1 * unit_possessed_time)

        else

            unit_possessed = fill(false, num_sellers)
            unit_possessed_time = 0
            unit_first_possessed_time = 0
            quality_of_unit_possessed = 0.0
            durability_of_unit_possessed = 0.0
            current_quality_of_unit_possessed = 0.0

        end

        # #

        new_buyer = buyer(b, #id
            my_neighbours, #neighbours
            qs_srp, #std_reservation_price
            future_discount, #future_discount
            initial_quality_expectation, #quality_expectation
            [], #quality_expectation_history
            initial_durability_expectation, #duration_expectation
            [], #duration_expectation_history
            0.0, #expected_surplus_buying
            unit_possessed, #unit_possessed
            false, #unit_possessed_is_new
            [], #unit_possessed_history,
            unit_possessed_time, #unit_possessed_time
            unit_first_possessed_time, # unit first possessed time
            quality_of_unit_possessed, #quality_of_unit_possessed
            [], #quality_of_unit_possessed_history
            durability_of_unit_possessed, #durability_of_unit_possessed
            [], #durability_of_unit_possessed_history
            current_quality_of_unit_possessed, #current_quality_of_unit_possessed
            fill(false, num_sellers), # unit for sale
            0.0, # quality_for_sale
            0.0, # durability_for_sale
            0.0, # price for sale
            0,
            [],
            [],
            [],
            [],
            rand(Uniform(0,1)),
            0) # age_for_sale
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

###### Load AUX

include(pwd() * "\\methods\\methods_visualization.jl")
include(pwd() * "\\methods\\methods_aux.jl")

__precompile__()

function sellers_choose_qp_k_d_m_states(buyers::Vector{buyer}, sellers::Vector{seller}, randPeriod::Int64, iter::Int64, num_buyers::Int64, profit_expected::Vector, μ_c::Float64, method_weight::String, product_life::Int64, samples_fixed::Bool)::Vector{seller}

    """

    Producers make decisions about changes in production process

    """

    for _seller in sellers

        δ_k = round(sample(0:0.01:0.05), digits = 2)
        δ_d = round(sample(0:0.01:0.05), digits = 2)
        δ_m = round(sample(0:0.01:0.05), digits = 2)
        δ_q = sample(1:5)

        if iter >= (randPeriod + 1)
            
                #change phase

                M_range = [max(_seller.margin_range[1], _seller.margin - δ_m), _seller.margin, min(_seller.margin_range[2], _seller.margin + δ_m)]

                Q_range = [max(0, _seller.quantity_produced - δ_q), _seller.quantity_produced, min(num_buyers, _seller.quantity_produced + δ_q)]

                K_range = [max(_seller.quality_range[1], _seller.quality - δ_k), _seller.quality, min(_seller.quality_range[2], _seller.quality + δ_k)]

                D_range = [max(_seller.durability_range[1], _seller.durability - δ_d), _seller.durability, min(_seller.durability_range[2], _seller.durability + δ_d)]

                o_p = mean(calculate_price_history.(sellers[Not(_seller.id)]; product_life =  product_life))[end]

                if _seller.consumer_research

                    sample_size = _seller.sample_size

                    e_k_all = getindex.(getfield.(buyers, :quality_expectation), _seller.id)
                    e_d_all = getindex.(getfield.(buyers, :durability_expectation), _seller.id)
                    ρ_all = getfield.(buyers, :future_discount)
                    β_all = getfield.(buyers, :std_reservation_price)
                    e_k_o_all = [mean(x[Not(_seller.id)]) for x in getfield.(buyers, :quality_expectation)]
                    e_d_o_all = [mean(x[Not(_seller.id)]) for x in getfield.(buyers, :durability_expectation)]

                    kdρ = [x for x in zip(e_k_all, e_d_all, ρ_all, β_all, e_k_o_all, e_d_o_all)]

                    if samples_fixed

                        if iter == (randPeriod + 1)

                            kdρ_idx = sample(1:length(kdρ), Int(round(sample_size * num_buyers)), replace = false)

                            push!(_seller.consumers_asked, kdρ_idx)

                        else

                            kdρ_idx = _seller.consumers_asked[end]
                            push!(_seller.consumers_asked, kdρ_idx)

                        end

                    else

                        kdρ_idx = sample(1:length(kdρ), Int(round(sample_size * num_buyers)), replace = false)
                        push!(_seller.consumers_asked, kdρ_idx)

                    end

                    observed_kdρ = kdρ[kdρ_idx]

                    observed_kdρ = sample(observed_kdρ, num_buyers, replace = true)

                    e_k = getindex.(observed_kdρ, 1)
                    e_d = getindex.(observed_kdρ, 2)
                    ρ_mean = getindex.(observed_kdρ, 3)
                    β_mean = getindex.(observed_kdρ, 4)

                    o_k = mean(getindex.(kdρ, 5))
                    o_d = mean(getindex.(kdρ, 6))

                else

                    e_k = fill(_seller.quality, num_buyers)
                    e_d = fill(_seller.durability, num_buyers)

                    o_k = mean(getfield.(sellers[Not(_seller.id)], :quality_history))[end]
                    o_d = mean(getfield.(sellers[Not(_seller.id)], :durability_history))[end]

                    ρ_mean = rand(Uniform(0.7, 1.0), num_buyers)
                    β_mean = rand(Uniform(0,1), num_buyers)

                end

                expected_profit_around = [calculate_state_profit(k, e_k, d, e_d, m, q, o_k, o_d, o_p, num_buyers, μ_c, _seller.cost_coefficient, ρ_mean, β_mean, "profit", product_life) for k in K_range, d in D_range, m in M_range, q in Q_range] # prior

                h_K = _seller.quality_history
                h_D = _seller.durability_history
                h_M = _seller.margin_history
                h_P = calculate_profit_history(_seller)
                h_Q = _seller.quantity_produced_history

                h_o_K = mean(getfield.(sellers[Not(_seller.id)], :quality_history))
                h_o_D = mean(getfield.(sellers[Not(_seller.id)], :durability_history))
                h_o_P = mean(calculate_price_history.(sellers[Not(_seller.id)]; product_life =  product_life))

                known_profits_around = [sum(h_P[(h_K .== k) .& (h_D .== d) .& (h_M .== m) .& (h_Q .== q) .& (h_o_K .== o_k).& (h_o_D .== o_d).& (h_o_P .== o_p)]) for k in K_range, d in D_range, m in M_range, q in Q_range]

                known_demands_around = [mean(h_Q[(h_K .== k) .& (h_D .== d) .& (h_M .== m) .& (h_Q .== q) .& (h_o_K .== o_k).& (h_o_D .== o_d).& (h_o_P .== o_p)]) for k in K_range, d in D_range, m in M_range, q in Q_range]

                known_items_around = [count((h_K .== k) .& (h_D .== d) .& (h_M .== m).& (h_Q .== q) .& (h_o_K .== o_k).& (h_o_D .== o_d).& (h_o_P .== o_p)) for k in K_range, d in D_range, m in M_range, q in Q_range]

                known_profits_around[isnan.(known_profits_around)] .= 0
                known_demands_around[isnan.(known_demands_around)] .= 0

                posterior_expected_profit_around = (expected_profit_around .+ known_profits_around) ./ (1 .+ known_items_around)

                weights = create_weights(vec(posterior_expected_profit_around), method_weight)

                optimal_profit_args = sample(vec(CartesianIndices(posterior_expected_profit_around)), Weights(weights))

                push!(profit_expected, (_seller.id, "p", posterior_expected_profit_around[optimal_profit_args]))

                new_quality = K_range[optimal_profit_args[1]]
                new_durability = D_range[optimal_profit_args[2]]
                new_margin = M_range[optimal_profit_args[3]]
                new_quantity = Q_range[optimal_profit_args[4]]

                expected_demand = calculate_state_profit(new_quality, e_k,new_durability, e_d, new_margin, new_quantity, o_k, o_d, o_p, num_buyers, μ_c, _seller.cost_coefficient, ρ_mean, β_mean, "demand", product_life)

                posterior_expected_demand = Int(ceil((expected_demand + known_demands_around[optimal_profit_args]) / (1 + known_items_around[optimal_profit_args])))

                if iter == (randPeriod + 1)

                    new_quantity = posterior_expected_demand

                end

                push!(profit_expected, (_seller.id, "qp", new_quantity))

                push!(profit_expected, (_seller.id, "ed", posterior_expected_demand))

                _seller.quality = new_quality
                _seller.durability = new_durability
                _seller.margin = new_margin
                _seller.quantity_produced = new_quantity

                @assert _seller.quality_range[1] <= _seller.quality <= _seller.quality_range[2]
                @assert _seller.durability_range[1] <= _seller.durability <= _seller.durability_range[2]
                @assert _seller.margin_range[1] <= _seller.margin <= _seller.margin_range[2]
                @assert _seller.quality_range[1] <= _seller.quality <= _seller.quality_range[2]

                @assert _seller.quantity_produced >= 0

        end

    end

    for _seller in sellers

        push!(_seller.quantity_produced_history, _seller.quantity_produced)
        push!(_seller.quality_history, _seller.quality)
        push!(_seller.durability_history, _seller.durability)
        push!(_seller.margin_history, _seller.margin)
        _seller.cost_of_production = calculate_cost(_seller; product_life = product_life) * _seller.quantity_produced
        push!(_seller.cost_of_production_history, _seller.cost_of_production)

    end

    return sellers

end


function consumers_compare_offers(buyers::Vector{buyer}, sellers::Vector{seller}, iter::Int64, product_life::Int64)::Vector{buyer}

    """

    Calculate further stream of quality, if unit is possessed
    For those who buy, expected surplus is stream of quality
    For others, expected_surplus is 0

    """

    for _buyer in buyers
        if any(_buyer.unit_possessed)
            _buyer.expected_surplus =  _buyer.std_reservation_price .* sum_of_geom_series_finite(_buyer.durability_of_unit_possessed^(iter - _buyer.unit_first_possessed_time) * _buyer.quality_of_unit_possessed, _buyer.future_discount * _buyer.durability_of_unit_possessed; t = _buyer.unit_first_possessed_time + product_life - iter)
        else
            _buyer.expected_surplus = 0.0
        end
    end

    return buyers

end

function consumers_make_decision(buyers::Vector{buyer}, sellers::Vector{seller}, num_sellers::Int64, iter::Int64, buyer_behaviour::String, secondary_market_exists::Bool, product_life::Int64)

    for _seller in sellers
        _seller.selling_income = 0.0
        _seller.expected_income = 0.0
    end

    shuffle!(buyers)

    for _buyer in buyers

        supply = getfield.(sellers, :quantity_produced) .- getfield.(sellers, :quantity_sold)

        supply_requirement = supply .> 0

        buy_prices = calculate_price.(sellers; product_life = product_life)

        utility = _buyer.std_reservation_price .* sum_of_geom_series_finite.(_buyer.quality_expectation, _buyer.future_discount * _buyer.durability_expectation; t = product_life) 
        reselling_price = 0.0

        if iter >= 2
            if any(_buyer.unit_possessed)
                reselling_price = secondary_market_exists .* _buyer.reselling_probability .* _buyer.expected_surplus
            end
        end

        buy_surplus = utility .+ reselling_price .- buy_prices

        buy_vs_keep = buy_surplus .- _buyer.expected_surplus

        push!(_buyer.choosing_surplus, buy_vs_keep)

        buy_requirement = buy_vs_keep .> 0

        requirements = vcat(supply_requirement .* buy_requirement .* buy_surplus)

        chosen_product = -1
        pm_surplus = 0.0

        if any(requirements .> 0)

            if buyer_behaviour == "deterministic"

                chosen_product = argmax(requirements)

            elseif buyer_behaviour == "stochastic"

                weight = requirements
                chosen_product = sample(1:length(requirements), Weights(weight))

            end

        end

        if chosen_product > 0

            if any(_buyer.unit_possessed)

                    # if previous product was bought -> buyer goes to secondary market

                _buyer.unit_for_sale = _buyer.unit_possessed
                _buyer.quality_for_sale = _buyer.current_quality_of_unit_possessed
                _buyer.durability_of_unit_possessed = _buyer.durability_of_unit_possessed
                _buyer.price_for_sale = _buyer.expected_surplus
                _buyer.age_for_sale = iter - _buyer.unit_first_possessed_time

            end

            push!(_buyer.unit_buying_selling_history, (i = _buyer.id, t=iter, d="buy, primary market", p=chosen_product, s=buy_surplus))

            pm_surplus = (utility .- buy_prices)[chosen_product]

            _buyer.unit_possessed_is_new = true
            _buyer.unit_possessed = create_bool_purchase(num_sellers, chosen_product)
            _buyer.unit_first_possessed_time = iter
            _buyer.unit_possessed_time = iter
            sellers[chosen_product].quantity_sold += 1
            push!(_buyer.realized_surplus_pm_history, pm_surplus)

            sellers[chosen_product].selling_income = sellers[chosen_product].selling_income + calculate_price(sellers[chosen_product]; product_life = product_life)

            sellers[chosen_product].expected_income = sellers[chosen_product].expected_income + calculate_price(sellers[chosen_product]; product_life = product_life)

        else

            push!(_buyer.realized_surplus_pm_history, 0.0)

        end

    end

    @assert all(getfield.(sellers, :quantity_produced) .>= getfield.(sellers, :quantity_sold))

    for _seller in sellers

        push!(_seller.quantity_sold_history, _seller.quantity_sold)
        push!(_seller.selling_income_history, _seller.selling_income)
        push!(_seller.expected_income_history, _seller.expected_income)

    end

    buyers = buyers[sortperm(getfield.(buyers, :id))]

    return buyers, sellers

end

function consumers_discover_q_d(buyers::Vector{buyer}, sellers::Vector{seller}, iter::Int64, ϵ_q::Float64=0.1, ϵ_d::Float64=0.05)::Vector{buyer}
    for _buyer in buyers
        if any(_buyer.unit_possessed)
            if _buyer.unit_first_possessed_time == iter
                chosen_product = argmax(_buyer.unit_possessed)
                q_lb = in_boundaries(sellers[chosen_product].quality - ϵ_q, 0.,1.)
                q_ub = in_boundaries(sellers[chosen_product].quality + ϵ_q, 0.,1.)
                if q_lb == q_ub
                    _buyer.quality_of_unit_possessed = q_lb
                    _buyer.current_quality_of_unit_possessed = q_lb
                else
                    _buyer.quality_of_unit_possessed = rand(Uniform(q_lb, q_ub))
                    _buyer.current_quality_of_unit_possessed = _buyer.quality_of_unit_possessed
                end

                d_lb = in_boundaries(sellers[chosen_product].durability - ϵ_d, 0., 1.)
                d_ub = in_boundaries(sellers[chosen_product].durability + ϵ_d, 0., 1.)

                if d_lb == d_ub
                    _buyer.durability_of_unit_possessed = d_lb
                else
                    _buyer.durability_of_unit_possessed = rand(Uniform(d_lb, d_ub))
                end
            end
        end
    end

    return buyers

end

function consumers_update_expectations(buyers, iter, λ_ind, λ_wom)

    for _buyer in buyers

        if any(_buyer.unit_possessed)

            if _buyer.unit_first_possessed_time == iter

                _buyer.quality_expectation = _buyer.quality_expectation .+ λ_ind .* _buyer.unit_possessed .* (_buyer.quality_of_unit_possessed .- _buyer.quality_expectation)

                _buyer.durability_expectation = _buyer.durability_expectation .+ λ_ind .* _buyer.unit_possessed .* (_buyer.durability_of_unit_possessed .- _buyer.durability_expectation)            

                """for idn in _buyer.neighbours                            

                    buyers[idn].quality_expectation = buyers[idn].quality_expectation .+ λ_wom .* _buyer.unit_possessed .* (_buyer.quality_of_unit_possessed .- buyers[idn].quality_expectation)

                    buyers[idn].durability_expectation = buyers[idn].durability_expectation .+ λ_wom .* _buyer.unit_possessed .* (_buyer.durability_of_unit_possessed .- buyers[idn].durability_expectation)

                end"""

                
            end

        end

        buying_neighbours = zeros(Int64, length(_buyer.unit_possessed))
        quality_neighbours = zeros(Float64, length(_buyer.unit_possessed))
        durability_neighbours = zeros(Float64, length(_buyer.unit_possessed))

        for idn in _buyer.neighbours

            if any(buyers[idn].unit_possessed)

                if buyers[idn].unit_first_possessed_time == iter

                    buying_neighbours = buying_neighbours .+ buyers[idn].unit_possessed
                    quality_neighbours[argmax(buyers[idn].unit_possessed)] += buyers[idn].quality_of_unit_possessed
                    durability_neighbours[argmax(buyers[idn].unit_possessed)] += buyers[idn].durability_of_unit_possessed

                end

            end

        end

        if any(buying_neighbours .> 0)

            _buyer.signal_volume += 1

            for i in eachindex(buying_neighbours)

                if buying_neighbours[i] > 0

                    average_quality_neighbours = quality_neighbours[i] / buying_neighbours[i]
                    average_durability_neighbours = durability_neighbours[i] / buying_neighbours[i]

                    _buyer.quality_expectation[i] = _buyer.quality_expectation[i] .+ λ_wom .* (average_quality_neighbours .- _buyer.quality_expectation[i])

                    _buyer.durability_expectation[i] = _buyer.durability_expectation[i] .+ λ_wom .* (average_durability_neighbours .- _buyer.durability_expectation[i]) 

                end
            
            end

        end

    end

    for _buyer in buyers

        push!(_buyer.quality_expectation_history, _buyer.quality_expectation)
        push!(_buyer.durability_expectation_history, _buyer.durability_expectation)

    end

    return buyers

end

function sellers_utilize_not_sold_products(sellers::Vector{seller}, μ_c::Float64, product_life::Int64)

    for _seller in sellers

        push!(_seller.utilization_cost_history, (1 - μ_c) * calculate_cost(_seller; product_life = product_life) * (_seller.quantity_produced - _seller.quantity_sold))
        _seller.quantity_sold = 0

    end

    return sellers

end

function buyers_products_age(buyers::Vector{buyer}, sellers::Vector{seller}, iter::Int64, product_life::Int64)

    for _buyer in buyers

        if any(_buyer.unit_possessed)

            if iter >= 2

                if (iter - _buyer.unit_first_possessed_time) > product_life

                    # breakage

                    # consumer loses further utility

                    push!(_buyer.unit_buying_selling_history, (i = _buyer.id, t=iter, d="destroyed", p=argmax(_buyer.unit_possessed)))

                    _buyer.unit_possessed = fill(false, length(_buyer.unit_possessed))
                    _buyer.quality_of_unit_possessed = 0.0
                    _buyer.current_quality_of_unit_possessed = 0.0
                    _buyer.durability_of_unit_possessed = 0.0
                    _buyer.unit_possessed_is_new = false
                    _buyer.unit_possessed_time = 0
                    _buyer.unit_first_possessed_time = 0
                    

                else

                    _buyer.current_quality_of_unit_possessed = _buyer.current_quality_of_unit_possessed * _buyer.durability_of_unit_possessed

                end

            end

        end

    end

    return buyers, sellers

end

function buyers_choose_secondary_market(buyers::Vector{buyer}, sellers::Vector{seller}, iter::Int64, buyer_behaviour::String, secondary_market_exists::Bool, method_weight::String, product_life::Int64)

    products_for_resale = any.(getfield.(buyers, :unit_for_sale))

    reselling_quantity = fill(0, length(sellers))
    
    if secondary_market_exists

        if any(products_for_resale)

            shuffle!(buyers)

            for _buyer in buyers

                if any(_buyer.unit_for_sale)

                    possessed_product = argmax(_buyer.unit_for_sale)

                    current_quality = getindex.(getfield.(buyers, :durability_expectation), possessed_product) .^ _buyer.age_for_sale .* getindex.(getfield.(buyers, :quality_expectation), possessed_product)

                    product_iters_of_life = product_life - _buyer.age_for_sale

                    q = getfield.(buyers, :future_discount) .* getindex.(getfield.(buyers, :durability_expectation), possessed_product)

                    offers = getfield.(buyers, :std_reservation_price) .* sum_of_geom_series_finite.(current_quality, q; t = product_iters_of_life) .* .!any.(getfield.(buyers, :unit_possessed))

                    if any(offers .> _buyer.price_for_sale)

                        if buyer_behaviour == "deterministic"

                            chosen_buyer = argmax(offers .- _buyer.price_for_sale)

                        elseif buyer_behaviour == "stochastic"

                            u = max.(offers .- _buyer.price_for_sale, 0)

                            weight = create_weights(u, method_weight)

                            chosen_buyer = sample(1:length(buyers), Weights(weight))

                        end

                        push!(_buyer.unit_buying_selling_history, (i = _buyer.id, t=iter, d="sell, secondary market", p=argmax(_buyer.unit_for_sale)))
                        push!(buyers[chosen_buyer].unit_buying_selling_history, (i = _buyer.id, t = iter, d="buy, secondary market", p = argmax(_buyer.unit_for_sale)))

                        reselling_quantity = reselling_quantity .+ _buyer.unit_for_sale

                        buyers[chosen_buyer].unit_possessed = _buyer.unit_for_sale
                        buyers[chosen_buyer].unit_possessed_time = iter 
                        buyers[chosen_buyer].unit_first_possessed_time = iter - _buyer.age_for_sale
                        buyers[chosen_buyer].quality_of_unit_possessed = _buyer.quality_for_sale
                        buyers[chosen_buyer].durability_of_unit_possessed = _buyer.durability_of_unit_possessed
                        buyers[chosen_buyer].current_quality_of_unit_possessed = _buyer.quality_for_sale
                        buyers[chosen_buyer].unit_possessed_is_new = false

                        if length(_buyer.realized_surplus_sm_s_history) == iter       

                            deleteat!(_buyer.realized_surplus_sm_b_history, iter)
                            deleteat!(_buyer.realized_surplus_sm_s_history, iter)
        
                        end

                        push!(_buyer.realized_surplus_sm_s_history, _buyer.price_for_sale)
                        push!(_buyer.realized_surplus_sm_b_history, 0.0)

                        if length(buyers[chosen_buyer].realized_surplus_sm_s_history) == iter       
                                        
                            deleteat!(buyers[chosen_buyer].realized_surplus_sm_b_history, iter)
                            deleteat!(buyers[chosen_buyer].realized_surplus_sm_s_history, iter)
        
                        end

                        push!(buyers[chosen_buyer].realized_surplus_sm_s_history, 0.0)
                        push!(buyers[chosen_buyer].realized_surplus_sm_b_history, offers[chosen_buyer] - _buyer.price_for_sale)

                        _buyer.reselling_probability = 1/2 * (_buyer.reselling_probability + 1)

                    else

                        _buyer.reselling_probability = 1/2 * (_buyer.reselling_probability + 0)

                    end

                    _buyer.unit_for_sale = fill(false, length(_buyer.unit_for_sale))
                    _buyer.price_for_sale = 0.0
                    _buyer.age_for_sale = 0

                else

                    if length(_buyer.realized_surplus_sm_b_history) != iter

                        push!(_buyer.realized_surplus_sm_s_history, 0.0)
                        push!(_buyer.realized_surplus_sm_b_history, 0.0)

                    end

                end

            end

            buyers = buyers[sortperm(getfield.(buyers, :id))]

        else

            for _buyer in buyers

                push!(_buyer.realized_surplus_sm_s_history, 0.0)
                push!(_buyer.realized_surplus_sm_b_history, 0.0)

            end

        end

    else

        for _buyer in buyers

            push!(_buyer.realized_surplus_sm_s_history, 0.0)
            push!(_buyer.realized_surplus_sm_b_history, 0.0)

        end

    end

    for _seller in sellers

        push!(_seller.reselling_history, reselling_quantity[_seller.id])

    end

    for _buyer in buyers
        if length(_buyer.realized_surplus_sm_b_history) < iter
            push!(_buyer.realized_surplus_sm_b_history, 0.0)
            push!(_buyer.realized_surplus_sm_s_history, 0.0)
        end
    end

    return buyers, sellers

end

function TO_GO(maxIter, num_sellers, num_buyers, num_links, c, m, network_type, λ_ind, λ_wom, buyer_behaviour, kr, dr, mr, μ_c, secondary_market_exists, rand_period, future_discount_range, method_weight, consumer_research, sample_size, product_life, samples_fixed)

    sellers = create_sellers(num_sellers, c, m, kr, dr, mr, consumer_research, sample_size)

    buyers_network = create_network(network_type, num_buyers = num_buyers, num_links = num_links)
    buyers = create_buyers(num_buyers, sellers, num_sellers, buyers_network,future_discount_range, product_life)

    profit_expected = []
 
    for iter in 0:maxIter

        if iter == 0

            sellers = sellers_choose_qp_k_d_m_states(buyers, sellers, rand_period, iter, num_buyers, profit_expected, μ_c, method_weight, product_life, samples_fixed)

        elseif iter == 1

            buyers = consumers_compare_offers(buyers, sellers, iter, product_life)

            buyers, sellers = consumers_make_decision(buyers, sellers, num_sellers, iter, buyer_behaviour, secondary_market_exists, product_life)

            buyers = consumers_discover_q_d(buyers, sellers, iter)

            buyers = consumers_update_expectations(buyers, iter, λ_ind, λ_wom)

            sellers = sellers_utilize_not_sold_products(sellers, μ_c, product_life)

            buyers, sellers = buyers_products_age(buyers, sellers, iter, product_life)

            sellers = sellers_choose_qp_k_d_m_states(buyers, sellers, rand_period, iter, num_buyers, profit_expected, μ_c, method_weight, product_life, samples_fixed)


        elseif iter >= 2

            buyers, sellers = buyers_products_age(buyers, sellers, iter, product_life)

            buyers = consumers_compare_offers(buyers, sellers, iter, product_life)

            buyers, sellers = consumers_make_decision(buyers, sellers, num_sellers, iter, buyer_behaviour, secondary_market_exists, product_life)

            buyers, sellers = buyers_choose_secondary_market(buyers, sellers, iter, buyer_behaviour, secondary_market_exists, method_weight, product_life)

            buyers = consumers_discover_q_d(buyers, sellers, iter)

            buyers = consumers_update_expectations(buyers, iter, λ_ind, λ_wom)

            sellers = sellers_utilize_not_sold_products(sellers, μ_c, product_life)

            if iter < maxIter

                sellers = sellers_choose_qp_k_d_m_states(buyers, sellers, rand_period, iter, num_buyers, profit_expected, μ_c, method_weight, product_life, samples_fixed)

            end

        end

    end

    return (buyers = buyers, sellers = sellers, profit_expected = profit_expected)

end
