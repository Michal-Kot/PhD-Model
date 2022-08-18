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

end

function create_sellers(num_sellers::Int64,k::Vector{Float64},c::Vector{Float64},m::Vector{Float64},d::Vector{Float64})::Vector{seller}
    @assert length(k) == num_sellers
    @assert length(c) == num_sellers
    @assert length(m) == num_sellers
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

function sellers_choose_qp_k_d_m(sellers::Vector{seller}, iter::Int64, num_buyers::Int64)::Vector{seller}
    for _seller in sellers
        if iter == 0

            #_seller.quality = rand(Uniform(0,2))
            #_seller.durability = rand(Uniform(0.5,1))
            _seller.quantity_produced = sample(1:num_buyers)
            #_seller.margin = rand(Uniform(0.5,1.5))

        elseif iter >= 1

            # 🍎 tu skończone, dopisać 

            #if _seller.phase == "setting"

                δ_k = 0.01
                δ_d = 0.01
                δ_m = 0.01

                my_k = _seller.quality
                my_k = my_k .+ [-δ_k, 0, δ_k]
                my_d = _seller.durability
                my_d = my_d .+ [-δ_d, 0, δ_d]
                my_m = _seller.margin
                my_m = my_m .+ [-δ_m, 0, δ_m]
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

                best_strategy = getindex.(strategy_outcomes,1)[argmax(getindex.(strategy_outcomes,2))]
                expected_marketshare = getindex.(strategy_outcomes,3)[argmax(getindex.(strategy_outcomes,2))]       
                
                new_quality = in_boundaries(best_strategy[1], 0.50, 1.50)
                new_durability = in_boundaries(best_strategy[2], 0.20, 0.50)
                new_margin = in_boundaries(best_strategy[3], 0.8, 2.)
                
                _seller.quality = new_quality
                _seller.durability = new_durability
                _seller.margin = new_margin
                _seller.quantity_produced = ceil(expected_marketshare / 100 * num_buyers)

                #_seller.phase = "evaluating"

            #elseif (iter >= 2) & (_seller.phase == "evaluating")



            #end

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

function consumers_make_decision(buyers::Vector{buyer}, sellers::Vector{seller}, num_sellers::Int64, iter::Int64)

    shuffle!(buyers)

    for _buyer in buyers
        if iter == 1

            supply = getfield.(sellers, :quantity_produced) .- getfield.(sellers, :quantity_sold)
            surplus_requirement = _buyer.surplus .> 0
            supply_requirement = supply .> 0
            if any(surplus_requirement .& supply_requirement)
                requirements = surplus_requirement .* supply_requirement .* _buyer.surplus
                chosen_product = argmax(requirements)

                if _buyer.id == 1

                    println("Buyer " * string(_buyer.id) * " bought product " * string(chosen_product) * " in iteration " * string(iter))

                end

                _buyer.unit_possessed_is_new = true
                _buyer.unit_possessed = create_bool_purchase(num_sellers, chosen_product)
                _buyer.unit_possessed_time = iter
                sellers[chosen_product].quantity_sold += 1
            end

        elseif iter >= 2

            supply = getfield.(sellers, :quantity_produced) .- getfield.(sellers, :quantity_sold)
            surplus_requirement = _buyer.surplus .> 0
            supply_requirement = supply .> 0
            if any(_buyer.unit_possessed)

                expected_surplus_from_possessed_product = _buyer.std_reservation_price * _buyer.current_quality_of_unit_possessed / (1 - _buyer.durability_of_unit_possessed)
                reselling_requirement = surplus_requirement .> expected_surplus_from_possessed_product

                x = ""

                if any(surplus_requirement .& supply_requirement)
                    x = x * "A"
                end

                if any(reselling_requirement)
                    x = x * "B"
                end

                if _buyer.unit_possessed_is_new
                    x = x * "C"
                end

                println(x)


                if any(surplus_requirement .& supply_requirement) & any(reselling_requirement) & _buyer.unit_possessed_is_new

                    println(string(_buyer.id) * " wants to resell!")

                    requirements = surplus_requirement .* supply_requirement .* _buyer.surplus
                    chosen_product = argmax(requirements)

                    if _buyer.id == 1

                        println("Buyer " * string(_buyer.id) * " bought product " * string(chosen_product) * " in iteration " * string(iter) * " and reselling")

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

                end

            else

                if any(surplus_requirement .& supply_requirement)

                    requirements = surplus_requirement .* supply_requirement .* _buyer.surplus
                    chosen_product = argmax(requirements)

                    if _buyer.id == 1

                        println("Buyer " * string(_buyer.id) * " bought product " * string(chosen_product) * " in iteration " * string(iter))

                    end

                    _buyer.unit_possessed = create_bool_purchase(num_sellers, chosen_product)
                    _buyer.unit_possessed_time = iter
                    _buyer.unit_possessed_is_new = true
                    sellers[chosen_product].quantity_sold += 1

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
                _buyer.durability_of_unit_possessed = rand(Uniform(d_lb, d_ub))
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

function buyers_choose_secondary_market(buyers::Vector{buyer}, sellers::Vector{seller}, iter)

    products_for_resale = any.(getfield.(buyers, :unit_for_sale))

    reselling_quantity = fill(0, length(sellers))

    if any(products_for_resale)

        for _buyer in buyers

            if any(_buyer.unit_for_sale)

                possessed_product = argmax(_buyer.unit_for_sale)

                #println("Secondary market exchange!")

                current_quality = getindex.(getfield.(buyers, :durability_expectation), possessed_product) .^ _buyer.age_for_sale .* getindex.(getfield.(buyers, :quality_expectation), possessed_product)

                offers = getfield.(buyers, :std_reservation_price) .* sum_of_geom_series.(current_quality, getindex.(getfield.(buyers, :durability_expectation), possessed_product)) .* .!any.(getfield.(buyers, :unit_possessed))
                
                println((maximum(offers), _buyer.price_for_sale))

                if any(offers .> _buyer.price_for_sale)

                    chosen_buyer = argmax(offers .- _buyer.price_for_sale)

                    reselling_quantity = reselling_quantity .+ _buyer.unit_for_sale

                    buyers[chosen_buyer].unit_possessed = _buyer.unit_for_sale
                    buyers[chosen_buyer].unit_possessed_time = iter - _buyer.age_for_sale
                    buyers[chosen_buyer].quality_of_unit_possessed = _buyer.quality_for_sale
                    buyers[chosen_buyer].durability_of_unit_possessed = _buyer.durability_of_unit_possessed
                    buyers[chosen_buyer].current_quality_of_unit_possessed = _buyer.quality_for_sale
                    buyers[chosen_buyer].unit_possessed_is_new = false

                end

                _buyer.unit_for_sale = fill(false, length(_buyer.unit_for_sale))
                _buyer.price_for_sale = 0.0
                _buyer.age_for_sale = 0

            end

        end

    end

    for _seller in sellers

        push!(_seller.reselling_history, reselling_quantity[_seller.id])

    end

    return buyers, sellers

end


#################################### AUX FUNCTIONS ##############################################################

num_sellers = 2
k = [1.2, 0.8]
c = [0.4, 0.5]
m = [1.0, 1.0]
d = [0.3, 0.3]

num_buyers = 200
num_links = 150

sellers = create_sellers(num_sellers, k, c, m, d)

buyers_network = create_network("random", num_buyers = num_buyers, num_links = num_links)
buyers = create_buyers(num_buyers, num_sellers, buyers_network, k, d)

getfield.(sellers, :quantity_sold)

maxIter = 150

for iter in 0:maxIter

    #println(iter)

    if iter == 0

        # iter = 0

        sellers = sellers_choose_qp_k_d_m(sellers, iter, num_buyers)

    elseif iter == 1

        # iter = 1

        buyers = consumers_compare_offers(buyers, sellers)

        buyers, sellers = consumers_make_decision(buyers, sellers, num_sellers, iter)

        #countmap(map(x -> all(diff(x) .== 0) ? -1 : argmax(x), getfield.(buyers, :unit_possessed)))

        buyers = consumers_discover_q_d(buyers, sellers, iter)

        #getfield.(buyers[any.(getfield.(buyers, :unit_possessed))], :durability_of_unit_possessed)

        buyers = consumers_update_expectations(buyers, iter, 0.25, 0.25)

        sellers = sellers_choose_qp_k_d_m(sellers, iter, num_buyers)

        sellers = sellers_utilize_not_sold_products(sellers)

    elseif iter >= 2

        # iter = 2,3...,T

        buyers = buyers_products_age(buyers, iter)

        buyers = consumers_compare_offers(buyers, sellers)

        buyers, sellers = consumers_make_decision(buyers, sellers, num_sellers, iter)

        #getfield.(buyers[any.(getfield.(buyers, :unit_for_sale))], :price_for_sale)

        buyers, sellers = buyers_choose_secondary_market(buyers, sellers, iter)

        buyers = consumers_discover_q_d(buyers, sellers, iter)

        buyers = consumers_update_expectations(buyers, iter, 0.25, 0.25)

        sellers = sellers_utilize_not_sold_products(sellers)

        if iter < maxIter

            sellers = sellers_choose_qp_k_d_m(sellers, iter, num_buyers)

        end

    end

end


plot(getfield.(sellers, :quantity_produced_history))
plot!(getfield.(sellers, :quantity_sold_history))

plot(getfield.(sellers, :reselling_history))

sum.([any.(x) for x in getfield.(buyers, :unit_possessed_history)])

plot(getindex.(getindex(getfield.(buyers, :surplus_history), 1), 1))

plot(calculate_profit_history.(sellers))

plot(calculate_price_history.(sellers))
plot(getfield.(sellers, :durability_history))
plot(getfield.(sellers, :quality_history))
plot(getfield.(sellers, :margin_history))
