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
    interest_rate::Float64

    cost_of_production::Float64
    cost_of_production_history::Vector{Float64}
    
    quantity_produced::Float64
    quantity_sold::Float64
    quantity_leased::Float64

    quantity_produced_history::Vector{Float64}
    quantity_sold_history::Vector{Float64}
    quantity_leased_history::Vector{Float64}

    quality_history::Vector{Float64}
    durability_history::Vector{Float64}
    margin_history::Vector{Float64}

    reselling_history::Vector{Int64}

    selling_income::Float64
    leasing_income::Float64

    selling_income_history::Vector{Float64}
    leasing_income_history::Vector{Float64}

    expected_income::Float64
    expected_income_history::Vector{Float64}

    adaptation_rate::Float64

    quality_range::Vector{Float64}
    durability_range::Vector{Float64}
    margin_range::Vector{Float64}

    utilization_cost_history::Vector{Float64}

end

function create_sellers(num_sellers::Int64,c::Vector{Float64},m::Vector{Float64}, kr::Vector{Vector{Float64}}, dr::Vector{Vector{Float64}}, mr::Vector{Vector{Float64}}, a::Vector{Float64}, inrt::Vector{Float64})::Vector{seller}
    
    @assert length(c) == num_sellers
    @assert length(m) == num_sellers
    @assert length(kr) == num_sellers
    @assert length(dr) == num_sellers
    @assert length(mr) == num_sellers

    sellers_vector = []

    for s in 1:num_sellers
        k = mean(kr[s])
        d = mean(dr[s])
        new_seller = seller(s, #id
        k, #quality
        d, #durability
        c[s], #cost_coefficient
        m[s], #margin
        inrt[s],
        0.0, #cost_of_production
        [], #cost_of_production_history
        0.0, #quantity produced
        0.0, #quantity_sold
        0.0, #quantity_leased
        [], #quantity_produced_history
        [], #quantity_sold_history
        [], #quantity_leased_history
        [], #quality_history
        [], #durability_history
        [], # margin_history
        [], #reselling_history
        0.0,
        0.0,
        [],
        [],
        0.0,
        [],
        a[s], #adaptation_rate
        kr[s], #quality_range
        dr[s], #durability_range
        mr[s], #margin_range
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
    unit_possessed_is_leased::Bool
    unit_possessed_history::Vector{Vector{Bool}}
    unit_possessed_time::Int64
    quality_of_unit_possessed::Float64
    quality_of_unit_possessed_history::Vector{Vector{Float64}}
    durability_of_unit_possessed::Float64
    durability_of_unit_possessed_history::Vector{Vector{Float64}}
    current_quality_of_unit_possessed::Float64
    leasing_interest::Float64
    leasing_total::Float64
    leasing_interest_history::Vector{Float64}
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

[rand(Uniform(x[1], x[2])) for x in [[1,2], [4,5]]]

function create_buyers(num_buyers::Int64, num_sellers::Int64, network::SimpleGraph{Int64}, initial_quality_expectation_range::Vector{Vector{Float64}}, initial_durability_expectation_range::Vector{Vector{Float64}})::Vector{buyer}
    buyers_vector = []
    for b in 1:num_buyers
        my_neighbours = neighbors(network, b)
        qs_srp = rand(Uniform(0,1))
        future_discount = rand(Uniform(0,1))
        initial_quality_expectation = mean.(initial_quality_expectation_range)
        initial_durability_expectation = mean.(initial_durability_expectation_range)
        new_buyer = buyer(b, #id
            my_neighbours, #neighbours
            qs_srp, #std_reservation_price
            future_discount, #future_discount
            initial_quality_expectation, #quality_expectation
            [], #quality_expectation_history
            initial_durability_expectation, #duration_expectation
            [], #duration_expectation_history
            0.0, #expected_surplus_buying
            fill(false, num_sellers), #unit_possessed
            false, #unit_possessed_is_new
            false, #unit_possessed_is_leased
            [], #unit_possessed_history,
            0, #unit_possessed_time
            0.0, #quality_of_unit_possessed
            [], #quality_of_unit_possessed_history
            0.0, #durability_of_unit_possessed
            [], #durability_of_unit_possessed_history
            0.0, #current_quality_of_unit_possessed
            0.0, #leasing interest
            0.0, # leasing total
            [], # leasing interest history
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

    """

    Producers make decisions about changes in production process

    """

    δ_k = rand() * 0.05
    δ_d = rand() * 0.05
    δ_m = rand() * 0.05
    δ_q = sample(1:5)

    for _seller in sellers

        if iter <= randPeriod

            _seller.quality = in_boundaries(_seller.quality + sample([-δ_k, 0, δ_k], Weights([1/3, 1/3, 1/3])), _seller.quality_range[1], _seller.quality_range[2])
            _seller.durability = in_boundaries(_seller.durability + sample([-δ_d, 0, δ_d], Weights([1/3, 1/3, 1/3])), _seller.durability_range[1], _seller.durability_range[2])
            _seller.margin = in_boundaries(_seller.margin + sample([-δ_m, 0, δ_m], Weights([1/3, 1/3, 1/3])), _seller.margin_range[1], _seller.margin_range[2])
            _seller.quantity_produced = sample(1:floor(0.2*num_buyers))

        elseif iter >= (randPeriod + 1)

            if (rand() < 1/3) & (iter >= 5)
                #eval phase

                if diff(_seller.expected_income_history .- _seller.cost_of_production_history)[end] < 0

                    _seller.quality = _seller.quality_history[end-1]
                    _seller.durability = _seller.durability_history[end-1]
                    _seller.margin = _seller.margin_history[end-1]

                    produced_q = _seller.quantity_produced_history[end]
                    leased_q = _seller.quantity_leased_history[end]
                    sold_q = _seller.quantity_sold_history[end]

                    new_quantity = _seller.quantity_produced_history[end-1]
                
                    if produced_q > (sold_q + leased_q)
    
                        new_quantity = ceil(new_quantity - _seller.adaptation_rate * (produced_q - sold_q - leased_q))
    
                    end

                    _seller.quantity_produced = new_quantity

                end
                
            else
                    #change phase
                δ_k_history = sign.(diff(_seller.quality_history))
                δ_d_history = sign.(diff(_seller.durability_history))
                δ_m_history = sign.(diff(_seller.margin_history))
                δ_q_history = sign.(diff(_seller.quantity_produced_history))

                δ_π_history = diff(_seller.expected_income_history .- _seller.cost_of_production_history)

                steps_k = [-1, 0, 1]
                steps_d = [-1, 0, 1]
                steps_m = [-1, 0, 1]
                steps_q = [-1, 0, 1]

                steps_k = steps_k[((_seller.quality .+ [-1,0,1] .* δ_k) .>= _seller.quality_range[1]) .& ((_seller.quality .+ [-1,0,1] .* δ_k) .<= _seller.quality_range[2])]

                steps_d = steps_d[((_seller.durability .+ [-1,0,1] .* δ_d) .>= _seller.durability_range[1]) .& ((_seller.durability .+ [-1,0,1] .* δ_d) .<= _seller.durability_range[2])]

                steps_m = steps_m[((_seller.margin .+ [-1,0,1] .* δ_m) .>= _seller.margin_range[1]) .& ((_seller.margin .+ [-1,0,1] .* δ_m) .<= _seller.margin_range[2])]

                steps_q = steps_q[(_seller.quantity_produced .+ [-1,0,1] .* δ_q) .>= 0]

                possible_k_d_m_q = vec(collect(Iterators.product(steps_k, steps_d, steps_m, steps_q)))

                expected_π = []

                for kdmq in possible_k_d_m_q

                    requirements = (δ_k_history .== kdmq[1]) .& (δ_d_history .== kdmq[2]) .& (δ_m_history .== kdmq[3]) .& (δ_q_history .== kdmq[4])

                    if any(requirements)

                        π_e = mean(δ_π_history[requirements])

                        push!(expected_π, π_e)

                    else

                        push!(expected_π, 0)

                    end

                end

                weights = expected_π .- minimum(expected_π) ./ (maximum(expected_π) - minimum(expected_π))
                best_strategy = sample(possible_k_d_m_q, Weights(weights))

                new_quality = _seller.quality + best_strategy[1] * δ_k
                new_durability = _seller.durability + best_strategy[2] * δ_d
                new_margin = _seller.margin + best_strategy[3] * δ_m
                new_quantity = _seller.quantity_produced + best_strategy[4] * δ_q

                produced_q = _seller.quantity_produced_history[end]
                leased_q = _seller.quantity_leased_history[end]
                sold_q = _seller.quantity_sold_history[end]
                
                if produced_q > (sold_q + leased_q)

                    new_quantity = ceil(new_quantity - _seller.adaptation_rate * (produced_q - sold_q - leased_q))

                end

                _seller.quality = new_quality
                _seller.durability = new_durability
                _seller.margin = new_margin
                _seller.quantity_produced = new_quantity

            end

        end

        push!(_seller.quantity_produced_history, _seller.quantity_produced)
        push!(_seller.quality_history, _seller.quality)
        push!(_seller.durability_history, _seller.durability)
        push!(_seller.margin_history, _seller.margin)
        _seller.cost_of_production = calculate_cost(_seller) * _seller.quantity_produced
        push!(_seller.cost_of_production_history, _seller.cost_of_production)

    end

    return sellers

end

function consumers_compare_offers(buyers::Vector{buyer}, sellers::Vector{seller}, iter::Int64)::Vector{buyer}

    """

    Calculate further stream of quality, if unit is possessed
    For those who lease, expected surplus is stream of quality - price of leasing
    For those who buy, expected surplus is stream of quality
    For others, expected_surplus is 0

    """

    for _buyer in buyers
        if any(_buyer.unit_possessed)
            _buyer.expected_surplus =  _buyer.std_reservation_price .* sum_of_geom_series(_buyer.durability_of_unit_possessed^(iter - _buyer.unit_possessed_time) * _buyer.quality_of_unit_possessed, _buyer.durability_of_unit_possessed)
        else
            _buyer.expected_surplus = 0.0
        end
    end

    return buyers

end

function consumers_make_decision(buyers::Vector{buyer}, sellers::Vector{seller}, num_sellers::Int64, iter::Int64, buyer_behaviour::String, secondary_market_exists::Bool, lease_available::Bool)

    for _seller in sellers
        _seller.selling_income = 0.0
        _seller.leasing_income = 0.0
        _seller.expected_income = 0.0
    end

    shuffle!(buyers)

    for _buyer in buyers

        supply = getfield.(sellers, :quantity_produced) .- getfield.(sellers, :quantity_sold) .- getfield.(sellers, :quantity_leased)

        supply_requirement = supply .> 0

        if any(supply_requirement)

            buy_prices = calculate_price.(sellers)
            lease_prices = calculate_lease_single.(sellers) .* 1 ./ (1 .- _buyer.future_discount) #pv of leasing

            utility = _buyer.std_reservation_price .* sum_of_geom_series.(_buyer.quality_expectation, _buyer.durability_expectation) # buying/leasing expected utility

            reselling_price = 0.0

            if iter >= 2
                if any(_buyer.unit_possessed)
                    if !_buyer.unit_possessed_is_leased
                        expected_surplus_from_possessed_product = _buyer.std_reservation_price * _buyer.current_quality_of_unit_possessed / (1 - _buyer.durability_of_unit_possessed)
                        reselling_price = secondary_market_exists .* _buyer.reselling_probability .* expected_surplus_from_possessed_product
                    end
                end
            end

            buy_surplus = utility .+ reselling_price .- buy_prices
            lease_surplus = utility .+ reselling_price .- lease_prices

            buy_vs_keep = buy_surplus .> (_buyer.expected_surplus .+ _buyer.leasing_total)
            lease_vs_keep = lease_surplus .> (_buyer.expected_surplus .+ _buyer.leasing_total)

            buy_requirement = buy_vs_keep .> 0
            lease_requirement = lease_vs_keep .> 0

            if !lease_available
                lease_requirement = falses(length(lease_requirement))
            end

            requirements = vcat(supply_requirement .* buy_requirement .* buy_surplus, supply_requirement .* lease_requirement .* lease_surplus)

            chosen_product = -1
            decision = ""
            pm_surplus = 0.0

            if any(requirements .> 0)

                decisions = repeat(["buy", "lease"], inner = length(sellers))

                if buyer_behaviour == "deterministic"

                    chosen_product = argmax(requirements)
                    decision = decisions[chosen_product] 

                    if chosen_product > length(sellers)
                        chosen_product = chosen_product - length(sellers)
                    end

                elseif buyer_behaviour == "stochastic"

                    weight = u2w(requirements, 0.0)
                    chosen_product = sample(1:length(requirements), Weights(weight))
                    decision = decisions[chosen_product]

                    if chosen_product > length(sellers)
                        chosen_product = chosen_product - length(sellers)
                    end

                end

            end

            if (chosen_product > 0) & (decision == "buy")

                if any(_buyer.unit_possessed)

                    if _buyer.unit_possessed_is_leased

                        # if previous product was leased -> seller sells it immidiately

                        sellers[argmax(_buyer.unit_possessed)].leasing_income = sellers[argmax(_buyer.unit_possessed)].leasing_income + _buyer.leasing_total

                    else

                        # if previous product was bought -> buyer goes to secondary market

                        _buyer.unit_for_sale = _buyer.unit_possessed
                        _buyer.quality_for_sale = _buyer.current_quality_of_unit_possessed
                        _buyer.durability_of_unit_possessed = _buyer.durability_of_unit_possessed
                        _buyer.price_for_sale = expected_surplus_from_possessed_product
                        _buyer.age_for_sale = iter - _buyer.unit_possessed_time

                    end

                end

                pm_surplus = (utility .- buy_prices)[chosen_product]

                _buyer.leasing_interest = 0.0
                _buyer.leasing_total = 0.0
                _buyer.unit_possessed_is_new = true
                _buyer.unit_possessed_is_leased = false
                _buyer.unit_possessed = create_bool_purchase(num_sellers, chosen_product)
                _buyer.unit_possessed_time = iter
                sellers[chosen_product].quantity_sold += 1
                push!(_buyer.realized_surplus_pm_history, pm_surplus)

                sellers[chosen_product].selling_income = sellers[chosen_product].selling_income + calculate_price(sellers[chosen_product])
                sellers[chosen_product].expected_income = sellers[chosen_product].expected_income + calculate_price(sellers[chosen_product])

            elseif (chosen_product > 0) & (decision == "lease")

                #lease new

                if any(_buyer.unit_possessed)

                    if _buyer.unit_possessed_is_leased

                        # if previous product was leased -> seller sells it immidiately

                        sellers[argmax(_buyer.unit_possessed)].leasing_income = sellers[argmax(_buyer.unit_possessed)].leasing_income + + _buyer.leasing_total

                    else

                        # if previous product was bought -> buyer goes to secondary market

                        _buyer.unit_for_sale = _buyer.unit_possessed
                        _buyer.quality_for_sale = _buyer.current_quality_of_unit_possessed
                        _buyer.durability_of_unit_possessed = _buyer.durability_of_unit_possessed
                        _buyer.price_for_sale = expected_surplus_from_possessed_product
                        _buyer.age_for_sale = iter - _buyer.unit_possessed_time

                    end
    
                end

                pm_surplus = utility[chosen_product]

                _buyer.leasing_interest = calculate_lease_single(sellers[chosen_product])
                _buyer.leasing_total = calculate_lease_total(sellers[chosen_product])
                _buyer.unit_possessed_is_new = true
                _buyer.unit_possessed_is_leased = true
                _buyer.unit_possessed = create_bool_purchase(num_sellers, chosen_product)
                _buyer.unit_possessed_time = iter
                sellers[chosen_product].quantity_leased += 1
                push!(_buyer.realized_surplus_pm_history, pm_surplus)
                sellers[chosen_product].expected_income = sellers[chosen_product].expected_income + calculate_lease_total(sellers[chosen_product])

            else

                push!(_buyer.realized_surplus_pm_history, 0.0)

            end

        else

            push!(_buyer.realized_surplus_pm_history, 0.0)

        end

    end

    for _seller in sellers

        push!(_seller.quantity_sold_history, _seller.quantity_sold)
        push!(_seller.quantity_leased_history, _seller.quantity_leased)
        push!(_seller.selling_income_history, _seller.selling_income)
        push!(_seller.expected_income_history, _seller.expected_income)

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

        push!(_seller.utilization_cost_history, μ_c * calculate_cost(_seller) * (_seller.quantity_produced - _seller.quantity_sold - _seller.quantity_leased))
        _seller.quantity_sold = 0
        _seller.quantity_leased = 0

    end

    return sellers

end

function buyers_products_age_and_lease(buyers::Vector{buyer}, sellers::Vector{seller}, iter::Int64)

    for _buyer in buyers

        if any(_buyer.unit_possessed)

            if iter >= 2

                if rand() > _buyer.durability_of_unit_possessed

                    # breakage

                    if _buyer.unit_possessed_is_leased

                        # pays total left interests

                        sellers[argmax(_buyer.unit_possessed)].leasing_income = sellers[argmax(_buyer.unit_possessed)].leasing_income + _buyer.leasing_total
                        push!(_buyer.leasing_interest_history, _buyer.leasing_total)

                    else

                        #does nothing

                        push!(_buyer.leasing_interest_history, 0.0)

                    end

                    _buyer.unit_possessed = fill(false, length(_buyer.unit_possessed))
                    _buyer.quality_of_unit_possessed = 0.0
                    _buyer.current_quality_of_unit_possessed = 0.0
                    _buyer.durability_of_unit_possessed = 0.0
                    _buyer.unit_possessed_is_new = false
                    _buyer.unit_possessed_time = 0
                    _buyer.unit_possessed_is_leased = false
                    _buyer.leasing_interest = 0.0
                    _buyer.leasing_total = 0.0
                    

                else

                    _buyer.current_quality_of_unit_possessed = _buyer.current_quality_of_unit_possessed * _buyer.durability_of_unit_possessed

                    if _buyer.unit_possessed_is_leased

                        # pays next interest

                        interest_paid = max(min(_buyer.leasing_interest, _buyer.leasing_total), 0)

                        sellers[argmax(_buyer.unit_possessed)].leasing_income = sellers[argmax(_buyer.unit_possessed)].leasing_income + interest_paid

                        push!(_buyer.leasing_interest_history, interest_paid)

                        _buyer.leasing_total = _buyer.leasing_total - interest_paid

                        if _buyer.leasing_total <= 0

                            _buyer.leasing_interest = 0.0

                        end

                    else

                        # does nothing

                        push!(_buyer.leasing_interest_history, 0.0)

                    end       

                end

            elseif iter == 1

                if _buyer.unit_possessed_is_leased

                    # pays next interest

                    interest_paid = max(min(_buyer.leasing_interest, _buyer.leasing_total), 0)

                    sellers[argmax(_buyer.unit_possessed)].leasing_income = sellers[argmax(_buyer.unit_possessed)].leasing_income + interest_paid

                    push!(_buyer.leasing_interest_history, interest_paid)

                    _buyer.leasing_total = _buyer.leasing_total - interest_paid

                    if _buyer.leasing_total <= 0

                        _buyer.leasing_interest = 0.0

                    end

                else

                    # does nothing

                    push!(_buyer.leasing_interest_history, 0.0)

                end       

            end

        else

            push!(_buyer.leasing_interest_history, 0.0)

        end

    end

    for _seller in sellers

        push!(_seller.leasing_income_history, _seller.leasing_income)

    end

    return buyers, sellers

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
                        buyers[chosen_buyer].unit_possessed_is_leased = false

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

function TO_GO(maxIter, num_sellers, num_buyers, num_links, c, m, network_type, λ_ind, λ_wom, buyer_behaviour, adaptation, inrt, kr, dr, mr, μ_c, secondary_market_exists, lease_available, rand_period = 10)

    sellers = create_sellers(num_sellers, c, m, kr, dr, mr, adaptation, inrt)

    buyers_network = create_network(network_type, num_buyers = num_buyers, num_links = num_links)
    buyers = create_buyers(num_buyers, num_sellers, buyers_network, kr, dr)

    for iter in 0:maxIter

        if iter == 0

            # iter = 0

            sellers = sellers_choose_qp_k_d_m_stochastic(sellers, iter, num_buyers, rand_period)

        elseif iter == 1

            # iter = 1

            buyers = consumers_compare_offers(buyers, sellers, iter)

            buyers, sellers = consumers_make_decision(buyers, sellers, num_sellers, iter, buyer_behaviour, secondary_market_exists, lease_available)

            buyers = consumers_discover_q_d(buyers, sellers, iter)

            buyers = consumers_update_expectations(buyers, iter, λ_ind, λ_wom)

            sellers = sellers_choose_qp_k_d_m_stochastic(sellers, iter, num_buyers, rand_period)

            sellers = sellers_utilize_not_sold_products(sellers, μ_c)

            buyers, sellers = buyers_products_age_and_lease(buyers, sellers, iter) # only pays interest of leasing

        elseif iter >= 2

            # iter = 2,3...,T

            buyers, sellers = buyers_products_age_and_lease(buyers, sellers, iter)

            buyers = consumers_compare_offers(buyers, sellers, iter)

            buyers, sellers = consumers_make_decision(buyers, sellers, num_sellers, iter, buyer_behaviour, secondary_market_exists, lease_available)

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
