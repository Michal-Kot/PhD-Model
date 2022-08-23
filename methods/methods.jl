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
p - price of product, q - quantity of sales, c - function of cost, linear on q

"""

mutable struct seller
    id::Int64
    quality::Float64
    durability::Float64
    cost_coefficient::Float64
    margin::Float64
    
    quantity_produced::Float64
    quantity_sold::Float64

    quantity_produced_history::Vector{Float64}
    quantity_sold_history::Vector{Float64}
    quality_history::Vector{Float64}
    durability_history::Vector{Float64}
    margin_history::Vector{Float64}

    reselling_history::Vector{Int64}

    quality_range::Vector{Float64}
    durability_range::Vector{Float64}
    margin_range::Vector{Float64}

    determinism::String

end

function create_sellers(num_sellers::Int64,k::Vector{Float64},c::Vector{Float64},m::Vector{Float64},d::Vector{Float64}, qr::Vector{Vector{Float64}}, dr::Vector{Vector{Float64}}, mr::Vector{Vector{Float64}}, ω::Vector{String})::Vector{seller}
    
    @assert length(k) == num_sellers
    @assert length(c) == num_sellers
    @assert length(m) == num_sellers
    @assert length(d) == num_sellers
    @assert length(qr) == num_sellers
    @assert length(dr) == num_sellers
    @assert length(mr) == num_sellers

    sellers_vector = []

    for s in 1:num_sellers
        new_seller = seller(s, #id
        k[s], #quality
        d[s], #durability
        c[s], #cost_coefficient
        m[s], #margin
        0.0, #quantity produced
        0.0, #quantity_sold
        [], #quantity_produced_history
        [], #quantity_sold_history
        [], #quality_history
        [], #durability_history
        [], # margin_history
        [],
        qr[s],
        dr[s],
        mr[s],
        ω[s]) #reselling_history
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
    std_reservation_price::Float64 #s
    quality_expectation::Vector{Float64}
    quality_expectation_history::Vector{Vector{Float64}}
    durability_expectation::Vector{Float64}
    durability_expectation_history::Vector{Vector{Float64}}
    surplus::Vector{Float64}
    unit_possessed::Vector{Bool}
    unit_possessed_is_new::Bool
    unit_possessed_history::Vector{Vector{Bool}}
    unit_possessed_time::Int64
    quality_of_unit_possessed::Float64
    quality_of_unit_possessed_history::Vector{Vector{Float64}}
    durability_of_unit_possessed::Float64
    durability_of_unit_possessed_history::Vector{Vector{Float64}}
    current_quality_of_unit_possessed::Float64
    surplus_history::Vector{Vector{Float64}}
    unit_for_sale::Vector{Bool}
    quality_for_sale::Float64
    durability_for_sale::Float64
    price_for_sale::Float64
    age_for_sale::Int64
    realized_surplus_pm_history::Vector{Float64}
    realized_surplus_sm_b_history::Vector{Float64}
    realized_surplus_sm_s_history::Vector{Float64}
end

function create_buyers(num_buyers::Int64, num_sellers::Int64, network::SimpleGraph{Int64}, initial_quality_expectation::Vector{Float64}=fill(1.0, num_sellers), initial_duration_expectation::Vector{Float64}=fill(0.5, num_sellers))::Vector{buyer}
    buyers_vector = []
    for b in 1:num_buyers
        my_neighbours = neighbors(network, b)
        qs_srp = rand(Uniform(0,1))
        new_buyer = buyer(b, #id
            my_neighbours, #neighbours
            qs_srp, #std_reservation_price
            initial_quality_expectation, #quality_expectation
            [], #quality_expectation_history
            initial_duration_expectation, #duration_expectation
            [], #duration_expectation_history
            [], #surplus
            fill(false, num_sellers), #unit_possessed
            false, #unit_possessed_is_new
            [], #unit_possessed_history,
            0, #unit_possessed_time
            0.0, #quality_of_unit_possessed
            [], #quality_of_unit_possessed_history
            0.0, #durability_of_unit_possessed
            [], #durability_of_unit_possessed_history
            0.0, #current_quality_of_unit_possessed
            [], #surplus_history
            fill(false, num_sellers), # unit for sale
            0.0, # quality_for_sale
            0.0, # durability_for_sale
            0.0, # price for sale
            0,
            [],
            [],
            []) # age_for_sale
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

function sellers_choose_qp_k_d_m_deterministic(sellers::Vector{seller}, iter::Int64, num_buyers::Int64)::Vector{seller}
    for _seller in sellers
        if iter == 0

            _seller.quantity_produced = sample(1:num_buyers)

        elseif iter >= 1

            # 🍎 tu skończone, dopisać 

            δ_k = 0.01
            δ_d = 0.01
            δ_m = 0.01

            my_k = _seller.quality
            my_k = my_k .+ [-δ_k, 0, δ_k]
            my_k = my_k[(my_k .>= _seller.quality_range[1]) .& (my_k .<= _seller.quality_range[2])]

            my_d = _seller.durability
            my_d = my_d .+ [-δ_d, 0, δ_d]
            my_d = my_d[(my_d .>= _seller.durability_range[1]) .& (my_d .<= _seller.durability_range[2])]

            my_m = _seller.margin
            my_m = my_m .+ [-δ_m, 0, δ_m]
            my_m = my_m[(my_m .>= _seller.margin_range[1]) .& (my_m .<= _seller.margin_range[2])]

            my_c = _seller.cost_coefficient

            combinations = vec(collect(Iterators.product(my_k, my_d, my_m)))
                
            comp = copy(sellers)
            comp = comp[Not(_seller.id)]

            comp_k = getfield.(comp, :quality)
            comp_d = getfield.(comp, :durability)
            comp_p = calculate_price.(comp)

            strategy_outcomes = []

            for comb in combinations

                my_k = comb[1]
                my_d = comb[2]
                my_m = comb[3]
                
                y = simulate_outcome_of_strategy(my_k, my_d, my_m, my_c, comp_k, comp_d, comp_p)
                
                push!(strategy_outcomes, (comb, y[1], y[2]))
                    
            end

            if _seller.determinism == "deterministic"

            best_strategy = getindex.(strategy_outcomes,1)[argmax(getindex.(strategy_outcomes,2))]    
                
            elseif _seller.determinism == "stochastic"

                weights = getindex.(strategy_outcomes,2)
                weights = max.(0, weights)
                best_strategy = sample(getindex.(strategy_outcomes,1), Weights(weights))

            end

            new_quality = best_strategy[1]
            new_durability = best_strategy[2]
            new_margin = best_strategy[3]

            # correction

            produced_q = _seller.quantity_produced_history[end]
            sold_q = _seller.quantity_sold_history[end]
            
            if produced_q > sold_q

                planned_q = ceil(produced_q - rand() * (produced_q - sold_q))

            elseif produced_q == sold_q

                planned_q = ceil(produced_q * (1 + rand()))

            end

            _seller.quality = new_quality
            _seller.durability = new_durability
            _seller.margin = new_margin
            _seller.quantity_produced = planned_q

        end

        push!(_seller.quantity_produced_history, _seller.quantity_produced)
        push!(_seller.quality_history, _seller.quality)
        push!(_seller.durability_history, _seller.durability)
        push!(_seller.margin_history, _seller.margin)

    end

    return sellers

end

function sellers_choose_qp_k_d_m_stochastic(sellers::Vector{seller}, iter::Int64, num_buyers::Int64)::Vector{seller}

    δ_k = 0.01
    δ_d = 0.01
    δ_m = 0.01

    for _seller in sellers

        if iter <= 50

            _seller.quality = _seller.quality + sample([-δ_k, 0, δ_k], Weights([1/3, 1/3, 1/3]))
            _seller.durability = _seller.durability + sample([-δ_d, 0, δ_d], Weights([1/3, 1/3, 1/3]))
            _seller.margin = _seller.margin + sample([-δ_m, 0, δ_m], Weights([1/3, 1/3, 1/3]))
            _seller.quantity_produced = sample(1:num_buyers)

        elseif iter >= 51

            δ_k_history = sign.(diff(_seller.quality_history))
            δ_d_history = sign.(diff(_seller.durability_history))
            δ_m_history = sign.(diff(_seller.margin_history))

            δ_π_history = diff(calculate_profit_history(_seller))

            possible_k_d_m = vec(collect(Iterators.product([-1, 0, 1], [-1, 0, 1], [-1, 0, 1])))

            expected_π = []

            for kdm in possible_k_d_m

                requirements = (δ_k_history .== kdm[1]) .& (δ_d_history .== kdm[2]) .& (δ_m_history .== kdm[3])

                if any(requirements)

                    π_e = mean(δ_π_history[requirements])
                    push!(expected_π, π_e)

                else

                    push!(expected_π, 0)

                end

            end

            if _seller.determinism == "deterministic"

                best_strategy = possible_k_d_m[argmax(expected_π)]    
                
            elseif _seller.determinism == "stochastic"

                weights = expected_π
                weights = max.(0, weights)

                best_strategy = sample(possible_k_d_m, Weights(weights))

            end

            new_quality = _seller.quality + best_strategy[1] * δ_k
            new_durability = _seller.durability + best_strategy[2] * δ_d
            new_margin = _seller.margin + best_strategy[3] * δ_m

            # correction

            produced_q = _seller.quantity_produced_history[end]
            sold_q = _seller.quantity_sold_history[end]
            
            if produced_q > sold_q

                planned_q = ceil(produced_q - rand() * (produced_q - sold_q))

            elseif produced_q == sold_q

                planned_q = ceil(produced_q * (1 + rand()))

            end

            _seller.quality = new_quality
            _seller.durability = new_durability
            _seller.margin = new_margin
            _seller.quantity_produced = planned_q

        end

        push!(_seller.quantity_produced_history, _seller.quantity_produced)
        push!(_seller.quality_history, _seller.quality)
        push!(_seller.durability_history, _seller.durability)
        push!(_seller.margin_history, _seller.margin)

    end

    return sellers

end

function sum_of_geom_series(a0,q)
    return a0 ./ (1 .- q)
end

function consumers_compare_offers(buyers::Vector{buyer}, sellers::Vector{seller})::Vector{buyer}

    prices = calculate_price.(sellers)
    for _buyer in buyers
        _buyer.surplus = _buyer.std_reservation_price .* sum_of_geom_series(_buyer.quality_expectation, _buyer.durability_expectation) .- prices
        push!(_buyer.surplus_history, _buyer.surplus)
    end

    return buyers

end

function consumers_make_decision(buyers::Vector{buyer}, sellers::Vector{seller}, num_sellers::Int64, iter::Int64, buyer_behaviour::String)

    shuffle!(buyers)

    for _buyer in buyers
        if iter == 1

            supply = getfield.(sellers, :quantity_produced) .- getfield.(sellers, :quantity_sold)
            surplus_requirement = _buyer.surplus .> 0
            supply_requirement = supply .> 0
            if any(surplus_requirement .& supply_requirement)
                requirements = surplus_requirement .* supply_requirement .* _buyer.surplus

                if buyer_behaviour == "deterministic"

                    chosen_product = argmax(requirements)

                elseif buyer_behaviour == "stochastic"

                    weight = u2w(requirements, 0.0)
                    chosen_product = sample(1:length(sellers), Weights(weight))

                end

                _buyer.unit_possessed_is_new = true
                _buyer.unit_possessed = create_bool_purchase(num_sellers, chosen_product)
                _buyer.unit_possessed_time = iter
                sellers[chosen_product].quantity_sold += 1
                push!(_buyer.realized_surplus_pm_history, requirements[chosen_product])

            else

                push!(_buyer.realized_surplus_pm_history, 0.0)
            
            end

        elseif iter >= 2

            supply = getfield.(sellers, :quantity_produced) .- getfield.(sellers, :quantity_sold)
            surplus_requirement = _buyer.surplus .> 0
            supply_requirement = supply .> 0
            if any(_buyer.unit_possessed)

                expected_surplus_from_possessed_product = _buyer.std_reservation_price * _buyer.current_quality_of_unit_possessed / (1 - _buyer.durability_of_unit_possessed)
                reselling_requirement = surplus_requirement .> expected_surplus_from_possessed_product

                if any(surplus_requirement .& supply_requirement) & any(reselling_requirement) & _buyer.unit_possessed_is_new

                    requirements = surplus_requirement .* supply_requirement .* _buyer.surplus

                    if buyer_behaviour == "deterministic"

                        chosen_product = argmax(requirements)
    
                    elseif buyer_behaviour == "stochastic"
    
                        weight = u2w(requirements, 0.0)
                        chosen_product = sample(1:length(sellers), Weights(weight))
    
                    end
    

                    _buyer.unit_for_sale = _buyer.unit_possessed
                    _buyer.quality_for_sale = _buyer.current_quality_of_unit_possessed
                    _buyer.durability_of_unit_possessed = _buyer.durability_of_unit_possessed
                    _buyer.price_for_sale = expected_surplus_from_possessed_product
                    _buyer.age_for_sale = iter - _buyer.unit_possessed_time

                    _buyer.unit_possessed = create_bool_purchase(num_sellers, chosen_product)
                    _buyer.unit_possessed_time = iter
                    _buyer.unit_possessed_is_new = true
                    sellers[chosen_product].quantity_sold += 1
                    push!(_buyer.realized_surplus_pm_history, requirements[chosen_product])

                else

                    push!(_buyer.realized_surplus_pm_history, 0.0)

                end

            else

                if any(surplus_requirement .& supply_requirement)

                    requirements = surplus_requirement .* supply_requirement .* _buyer.surplus

                    if buyer_behaviour == "deterministic"

                        chosen_product = argmax(requirements)
    
                    elseif buyer_behaviour == "stochastic"
    
                        weight = u2w(requirements, 0.0)
                        chosen_product = sample(1:length(sellers), Weights(weight))
    
                    end

                    _buyer.unit_possessed = create_bool_purchase(num_sellers, chosen_product)
                    _buyer.unit_possessed_time = iter
                    _buyer.unit_possessed_is_new = true
                    sellers[chosen_product].quantity_sold += 1
                    push!(_buyer.realized_surplus_pm_history, requirements[chosen_product])

                else

                    push!(_buyer.realized_surplus_pm_history, 0.0)

                end
            end
        end

        push!(_buyer.unit_possessed_history, _buyer.unit_possessed)

    end

    for _seller in sellers
        push!(_seller.quantity_sold_history, _seller.quantity_sold)
    end

    buyers = buyers[sortperm(getfield.(buyers, :id))]

    return buyers, sellers

end

function consumers_discover_q_d(buyers::Vector{buyer}, sellers::Vector{seller}, iter::Int64, ϵ_q::Float64=0.33, ϵ_d::Float64=0.10)::Vector{buyer}
    for _buyer in buyers
        if any(_buyer.unit_possessed)
            if _buyer.unit_possessed_time == iter
                chosen_product = argmax(_buyer.unit_possessed)
                q_lb = in_boundaries(sellers[chosen_product].quality - ϵ_q, 0., 2.)
                q_ub = in_boundaries(sellers[chosen_product].quality + ϵ_q, 0., 2.)
                if q_lb == q_ub
                    _buyer.quality_of_unit_possessed = q_lb
                    _buyer.current_quality_of_unit_possessed = q_lb
                else
                    _buyer.quality_of_unit_possessed = rand(Uniform(q_lb, q_ub))
                    _buyer.current_quality_of_unit_possessed = _buyer.quality_of_unit_possessed
                end

                d_lb = in_boundaries(sellers[chosen_product].durability - ϵ_d, 0., 1.)
                d_ub = in_boundaries(sellers[chosen_product].durability + ϵ_d, 0., 1.
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

function simulate_outcome_of_strategy(my_k, my_d, my_m, my_c, comp_k, comp_d, comp_p)

    s = LinRange(0:0.01:1.0)

    my_p = my_c * my_k / (1 - my_d) * my_m
    my_w = s .* my_k / (1 - my_d) .- my_p
    my_w_d = s .* my_k / (1 - my_d) * (1 - my_d^2) .- my_p

    comp_w = [s .* ck / (1 - cd) .- cp for (ck, cd, cp) in zip(comp_k, comp_d, comp_p)]

    is_my_w_positive = my_w .>= 0
    is_my_w_d_positive = my_w_d .>= 0
    is_my_w_best = reduce(.*, [my_w .> cw for cw in comp_w])

    my_demand = count(is_my_w_best .& is_my_w_positive .& is_my_w_d_positive)
    my_profit = my_demand * (my_p - my_c * my_k / (1 - my_d))

    return my_profit, my_demand

end


function consumers_update_expectations(buyers, iter, λ_ind, λ_wom)

    for _buyer in buyers
        if any(_buyer.unit_possessed)
            if _buyer.unit_possessed_time == iter

                _buyer.quality_expectation = _buyer.quality_expectation .+ λ_ind .* _buyer.unit_possessed .* (_buyer.quality_of_unit_possessed .- _buyer.quality_expectation)        
                _buyer.durability_expectation = _buyer.durability_expectation .+ λ_ind .* _buyer.unit_possessed .* (_buyer.durability_of_unit_possessed .- _buyer.durability_expectation)

                for idn in _buyer.neighbours
                    buyers[idn].quality_expectation = buyers[idn].quality_expectation .+ λ_wom .* _buyer.unit_possessed .* (_buyer.quality_of_unit_possessed .- buyers[idn].quality_expectation)
                    buyers[idn].durability_expectation = buyers[idn].durability_expectation .+ λ_wom .* _buyer.unit_possessed .* (_buyer.durability_of_unit_possessed .- buyers[idn].durability_expectation)
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

function sellers_utilize_not_sold_products(sellers::Vector{seller})

    for _seller in sellers

        _seller.quantity_sold = 0

    end

    return sellers

end

function buyers_products_age(buyers::Vector{buyer}, iter::Int64)

    for _buyer in buyers

        if any(_buyer.unit_possessed)

            if rand() > (_buyer.durability_of_unit_possessed ^ (iter - _buyer.unit_possessed_time))

                # breakage

                _buyer.unit_possessed = fill(false, length(_buyer.unit_possessed))
                _buyer.quality_of_unit_possessed = 0.0
                _buyer.current_quality_of_unit_possessed = 0.0
                _buyer.durability_of_unit_possessed = 0.0
                _buyer.unit_possessed_is_new = false
                _buyer.unit_possessed_time = 0

            else

                _buyer.current_quality_of_unit_possessed = _buyer.current_quality_of_unit_possessed * _buyer.durability_of_unit_possessed

            end

        end

    end

    return buyers

end

function buyers_choose_secondary_market(buyers::Vector{buyer}, sellers::Vector{seller}, iter::Int64, buyer_behaviour::String)

    products_for_resale = any.(getfield.(buyers, :unit_for_sale))

    reselling_quantity = fill(0, length(sellers))

    if any(products_for_resale)

        shuffle!(buyers)

        for _buyer in buyers

            if any(_buyer.unit_for_sale)

                possessed_product = argmax(_buyer.unit_for_sale)

                current_quality = getindex.(getfield.(buyers, :durability_expectation), possessed_product) .^ _buyer.age_for_sale .* getindex.(getfield.(buyers, :quality_expectation), possessed_product)

                offers = getfield.(buyers, :std_reservation_price) .* sum_of_geom_series.(current_quality, getindex.(getfield.(buyers, :durability_expectation), possessed_product)) .* .!any.(getfield.(buyers, :unit_possessed))

                if any(offers .> _buyer.price_for_sale)

                    if buyer_behaviour == "deterministic"

                        chosen_buyer = argmax(offers .- _buyer.price_for_sale)

                    elseif buyer_behaviour == "stochastic"

                        u = max.(offers .- _buyer.price_for_sale, 0)

                        weight = u2w(u, 0.0)

                        chosen_buyer = sample(1:length(buyers), Weights(weight))

                    end

                    reselling_quantity = reselling_quantity .+ _buyer.unit_for_sale

                    buyers[chosen_buyer].unit_possessed = _buyer.unit_for_sale
                    buyers[chosen_buyer].unit_possessed_time = iter - _buyer.age_for_sale
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

function TO_GO(maxIter, num_sellers, num_buyers, num_links, k, c, m, d, network_type, λ_ind, λ_wom, buyer_behaviour, determinism, qr, dr, mr)

    sellers = create_sellers(num_sellers, k, c, m, d, qr, dr, mr, determinism)

    buyers_network = create_network(network_type, num_buyers = num_buyers, num_links = num_links)
    buyers = create_buyers(num_buyers, num_sellers, buyers_network, k, d)

    for iter in 0:maxIter

        if iter == 0

            # iter = 0

            sellers = sellers_choose_qp_k_d_m_stochastic(sellers, iter, num_buyers)

        elseif iter == 1

            # iter = 1

            buyers = consumers_compare_offers(buyers, sellers)

            buyers, sellers = consumers_make_decision(buyers, sellers, num_sellers, iter, buyer_behaviour)

            buyers = consumers_discover_q_d(buyers, sellers, iter)

            buyers = consumers_update_expectations(buyers, iter, λ_ind, λ_wom)

            sellers = sellers_choose_qp_k_d_m_stochastic(sellers, iter, num_buyers)

            sellers = sellers_utilize_not_sold_products(sellers)

        elseif iter >= 2

            # iter = 2,3...,T

            buyers = buyers_products_age(buyers, iter)

            buyers = consumers_compare_offers(buyers, sellers)

            buyers, sellers = consumers_make_decision(buyers, sellers, num_sellers, iter, buyer_behaviour)

            #getfield.(buyers[any.(getfield.(buyers, :unit_for_sale))], :price_for_sale)

            buyers, sellers = buyers_choose_secondary_market(buyers, sellers, iter, buyer_behaviour)

            buyers = consumers_discover_q_d(buyers, sellers, iter)

            buyers = consumers_update_expectations(buyers, iter, λ_ind, λ_wom)

            sellers = sellers_utilize_not_sold_products(sellers)

            if iter < maxIter

                sellers = sellers_choose_qp_k_d_m_stochastic(sellers, iter, num_buyers)

            end

        end

    end

    return (buyers = buyers, sellers = sellers)

end
