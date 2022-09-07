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

    adaptation_rate::Float64
    increase_q::Float64

    utlization_cost_history::Vector{Float64}

end

function create_sellers(num_sellers::Int64,k::Vector{Float64},c::Vector{Float64},m::Vector{Float64},d::Vector{Float64}, qr::Vector{Vector{Float64}}, dr::Vector{Vector{Float64}}, mr::Vector{Vector{Float64}}, a::Vector{Float64}, i::Vector{Float64})::Vector{seller}
    
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
        a[s],
        i[s],
        []) #reselling_history
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
    reselling_probability::Float64
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
            [],
            rand(Uniform(0,1))) # age_for_sale
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

function sellers_choose_qp_k_d_m_stochastic(sellers::Vector{seller}, iter::Int64, num_buyers::Int64, randPeriod::Int64)::Vector{seller}

    δ_k = 0.01
    δ_d = 0.01
    δ_m = 0.00

    for _seller in sellers

        if iter <= randPeriod

            _seller.quality = in_boundaries(_seller.quality + sample([-δ_k, 0, δ_k], Weights([1/3, 1/3, 1/3])), _seller.quality_range[1], _seller.quality_range[2])
            _seller.durability = in_boundaries(_seller.durability + sample([-δ_d, 0, δ_d], Weights([1/3, 1/3, 1/3])), _seller.durability_range[1], _seller.durability_range[2])
            _seller.margin = in_boundaries(_seller.margin + sample([-δ_m, 0, δ_m], Weights([1/3, 1/3, 1/3])), _seller.margin_range[1], _seller.margin_range[2])
            _seller.quantity_produced = sample(1:floor(0.2*num_buyers))

        elseif iter >= (randPeriod + 1)

            δ_k_history = sign.(diff(_seller.quality_history))
            δ_d_history = sign.(diff(_seller.durability_history))
            δ_m_history = sign.(diff(_seller.margin_history))

            δ_π_history = diff(calculate_profit_history(_seller))

            steps_k = [-1, 0, 1]
            steps_d = [-1, 0, 1]
            steps_m = [-1, 0, 1]

            steps_k = steps_k[((_seller.quality .+ [-1,0,1] .* δ_k) .>= _seller.quality_range[1]) .& ((_seller.quality .+ [-1,0,1] .* δ_k) .<= _seller.quality_range[2])]

            steps_d = steps_d[((_seller.durability .+ [-1,0,1] .* δ_d) .>= _seller.durability_range[1]) .& ((_seller.durability .+ [-1,0,1] .* δ_d) .<= _seller.durability_range[2])]

            steps_m = steps_m[((_seller.margin .+ [-1,0,1] .* δ_m) .>= _seller.margin_range[1]) .& ((_seller.margin .+ [-1,0,1] .* δ_m) .<= _seller.margin_range[2])]

            possible_k_d_m = vec(collect(Iterators.product(steps_k, steps_d, steps_m)))

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

            #if _seller.determinism == "deterministic"

            #    best_strategy = possible_k_d_m[argmax(expected_π)]    
                
            #elseif _seller.determinism == "stochastic"

            #if _seller.id == 1

            #    println(round.(expected_π, digits=2))

            #end

            #possible_k_d_m = possible_k_d_m[expected_π .>= 0]
            #expected_π = expected_π[expected_π .>= 0]

            weights = expected_π
            #weights = max.(0, weights)
            weights = (weights .- minimum(weights)) ./ (maximum(weights) - minimum(weights))
            best_strategy = sample(possible_k_d_m, Weights(weights))

            #end

            new_quality = _seller.quality + best_strategy[1] * δ_k
            new_durability = _seller.durability + best_strategy[2] * δ_d
            new_margin = _seller.margin + best_strategy[3] * δ_m

            # correction

            produced_q = _seller.quantity_produced_history[end]
            sold_q = _seller.quantity_sold_history[end]
            
            if produced_q > sold_q

                planned_q = ceil(produced_q - _seller.adaptation_rate * (produced_q - sold_q))

            elseif produced_q == sold_q

                planned_q = ceil(produced_q * (1 + _seller.increase_q))

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

buyers = sim_single.buyers
secondary_market_exists = true

function consumers_make_decision2(buyers::Vector{buyer}, sellers::Vector{seller}, iter::Int64, secondary_market_exists::Bool)

    num_sellers = length(sellers)
    num_buyers = length(buyers)

    already_bougth = falses(num_buyers)

    if iter == 1

        ladder = collect(zip(vcat([fill(x, num_sellers) for x in 1:num_buyers]...), vcat(fill(collect(1:num_sellers), num_buyers)...), vcat(getfield.(buyers, :surplus)...)))

        ordered_ladder = ladder[sortperm(getindex.(ladder,3), rev = true)]

    elseif iter >= 2

        surplus = getfield.(buyers, :surplus)
        surplus_selling = [secondary_market_exists .* _buyer.reselling_probability .* _buyer.std_reservation_price * _buyer.unit_possessed_is_new * _buyer.current_quality_of_unit_possessed / (1 - _buyer.durability_of_unit_possessed) for _buyer in buyers]
        surplus_total = [sur .+ sur_sel for (sur, sur_sel) in zip(surplus, surplus_selling)]

        ladder = collect(zip(vcat([fill(x, num_sellers) for x in 1:num_buyers]...), vcat(fill(collect(1:num_sellers), num_buyers)...), vcat(surplus_total...)))

        ordered_ladder = ladder[sortperm(getindex.(ladder,3), rev = true)]

    end

    for ladder_item in ordered_ladder

        supply = getfield.(sellers, :quantity_produced) .- getfield.(sellers, :quantity_sold)

        if (ladder_item[3] > 0) & (supply[ladder_item[2]] > 0) & !(already_bougth[ladder_item[1]])

            chosen_product = ladder_item[2]
            choosing_buyer = ladder_item[1]

            if (iter >= 2) & any(buyers[choosing_buyer].unit_possessed)

                buyers[choosing_buyer].unit_for_sale = buyers[choosing_buyer].unit_possessed
                buyers[choosing_buyer].quality_for_sale = buyers[choosing_buyer].current_quality_of_unit_possessed
                buyers[choosing_buyer].durability_of_unit_possessed = buyers[choosing_buyer].durability_of_unit_possessed
                buyers[choosing_buyer].price_for_sale = buyers[choosing_buyer].std_reservation_price * buyers[choosing_buyer].current_quality_of_unit_possessed / (1 - buyers[choosing_buyer].durability_of_unit_possessed)
                buyers[choosing_buyer].age_for_sale = iter - buyers[choosing_buyer].unit_possessed_time

            end

            buyers[choosing_buyer].unit_possessed_is_new = true
            buyers[choosing_buyer].unit_possessed = create_bool_purchase(num_sellers, chosen_product)
            buyers[choosing_buyer].unit_possessed_time = iter
            sellers[chosen_product].quantity_sold += 1

            already_bougth[ladder_item[1]] = true

        end

    end

    for _seller in sellers
        push!(_seller.quantity_sold_history, _seller.quantity_sold)
    end

    return buyers, sellers

end





function consumers_make_decision(buyers::Vector{buyer}, sellers::Vector{seller}, num_sellers::Int64, iter::Int64, buyer_behaviour::String, secondary_market_exists::Bool)

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
            supply_requirement = supply .> 0

            surplus_requirement = _buyer.surplus .> 0

            if any(_buyer.unit_possessed)

                #if secondary_market_exists

                    expected_surplus_from_possessed_product = _buyer.std_reservation_price * _buyer.current_quality_of_unit_possessed / (1 - _buyer.durability_of_unit_possessed)

                    reselling_requirement = (_buyer.surplus .+ (secondary_market_exists .* _buyer.reselling_probability .* expected_surplus_from_possessed_product))  .> expected_surplus_from_possessed_product

                    if any(surplus_requirement .& supply_requirement) & any(reselling_requirement) & _buyer.unit_possessed_is_new

                        requirements = surplus_requirement .* supply_requirement .* (_buyer.surplus .+ secondary_market_exists * _buyer.reselling_probability .* expected_surplus_from_possessed_product)

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

                #else

                #    push!(_buyer.realized_surplus_pm_history, 0.0)                  

                #end

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

function sellers_utilize_not_sold_products(sellers::Vector{seller}, μ_c::Float64)

    for _seller in sellers

        push!(_seller.utlization_cost_history, μ_c * (_seller.quantity_produced - _seller.quantity_sold))
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

function buyers_choose_secondary_market(buyers::Vector{buyer}, sellers::Vector{seller}, iter::Int64, buyer_behaviour::String, secondary_market_exists::Bool)

    products_for_resale = any.(getfield.(buyers, :unit_for_sale))

    reselling_quantity = fill(0, length(sellers))
    
    if secondary_market_exists

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

function TO_GO(maxIter, num_sellers, num_buyers, num_links, k, c, m, d, network_type, λ_ind, λ_wom, buyer_behaviour, adaptation, increase_q, qr, dr, mr, μ_c, secondary_market_exists, rand_period = 20)

    sellers = create_sellers(num_sellers, k, c, m, d, qr, dr, mr, adaptation, increase_q)

    buyers_network = create_network(network_type, num_buyers = num_buyers, num_links = num_links)
    buyers = create_buyers(num_buyers, num_sellers, buyers_network, k, d)

    for iter in 0:maxIter

        if iter == 0

            # iter = 0

            sellers = sellers_choose_qp_k_d_m_stochastic(sellers, iter, num_buyers, rand_period)

        elseif iter == 1

            # iter = 1

            buyers = consumers_compare_offers(buyers, sellers)

            buyers, sellers = consumers_make_decision2(buyers, sellers, iter, secondary_market_exists)

            buyers = consumers_discover_q_d(buyers, sellers, iter)

            buyers = consumers_update_expectations(buyers, iter, λ_ind, λ_wom)

            sellers = sellers_choose_qp_k_d_m_stochastic(sellers, iter, num_buyers, rand_period)

            sellers = sellers_utilize_not_sold_products(sellers, μ_c)

        elseif iter >= 2

            # iter = 2,3...,T

            buyers = buyers_products_age(buyers, iter)

            buyers = consumers_compare_offers(buyers, sellers)

            buyers, sellers = consumers_make_decision2(buyers, sellers, iter, secondary_market_exists)

            #getfield.(buyers[any.(getfield.(buyers, :unit_for_sale))], :price_for_sale)

            buyers, sellers = buyers_choose_secondary_market(buyers, sellers, iter, buyer_behaviour, secondary_market_exists)

            buyers = consumers_discover_q_d(buyers, sellers, iter)

            buyers = consumers_update_expectations(buyers, iter, λ_ind, λ_wom)

            sellers = sellers_utilize_not_sold_products(sellers, μ_c)

            if iter < maxIter

                sellers = sellers_choose_qp_k_d_m_stochastic(sellers, iter, num_buyers, rand_period)

            end

        end

    end

    return (buyers = buyers, sellers = sellers)

end
