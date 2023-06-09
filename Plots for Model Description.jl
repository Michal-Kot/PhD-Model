############## użyteczność vs. k i d

using Plots
using LaTeXStrings
pyplot()

function utility_stream(k,d;h)
    return k * (1-d^h)/(1-d)
end

K = LinRange(0.010:0.010:0.990)
D = LinRange(0.010:0.010:0.990)

plotattr(:countour)

ac_p1 = Plots.contour(K,D, (K,D)->utility_stream(K,D;h=1), levels=5, xlabel = "Jakość dobra, " * L"K_{jt}", ylabel = "Trwałość dobra, " * L"D_{jt}", xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6, clabels=true, color=:turbo, colorbar_title = L"W_{jt}", title = "H = 1", titlefontsize = 10, colorbar_titlefontsize = 8, colorbar_tickfontsize = 6)
ac_p2 = Plots.contour(K,D, (K,D)->utility_stream(K,D;h=2), levels=5, xlabel = "Jakość dobra, " * L"K_{jt}", ylabel = "Trwałość dobra, " * L"D_{jt}", xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6, clabels=true, color=:turbo, colorbar_title = L"W_{jt}", title = "H = 2", titlefontsize = 10, colorbar_titlefontsize = 8, colorbar_tickfontsize = 6)
ac_p5 = Plots.contour(K,D, (K,D)->utility_stream(K,D;h=5), levels=5, xlabel = "Jakość dobra, " * L"K_{jt}", ylabel = "Trwałość dobra, " * L"D_{jt}", xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6, clabels=true, color=:turbo, colorbar_title = L"W_{jt}", title = "H = 5", titlefontsize = 10, colorbar_titlefontsize = 8, colorbar_tickfontsize = 6)
ac_p10 = Plots.contour(K,D, (K,D)->utility_stream(K,D;h=10), levels=5, xlabel = "Jakość dobra, " * L"K_{jt}", ylabel = "Trwałość dobra, " * L"D_{jt}", xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6, clabels=true, color=:turbo, colorbar_title = L"W_{jt}", title = "H = 10", titlefontsize = 10, colorbar_titlefontsize = 8, colorbar_tickfontsize = 6)

ac_p = Plots.plot(ac_p1, ac_p2, ac_p5, ac_p10, layout = (2,2))

Plots.savefig(ac_p, pwd() * "\\plots_to_export\\contour_utility_stream.pdf")

############## sieć powiązań

using Graphs
using GraphPlot
using Compose, Cairo, Fontconfig

function rectangle_graph(a::Int,b::Int,T=SimpleWeightedGraph; seed=0, rrand=0.9, scale=1.0)
    # utils
    mt = MersenneTwister(seed)
    add_edge2!(g,a,b,w) = add_edge!.(Ref(g),[a,b],[b,a],w) 
    r() = rand(mt)*rrand-rrand/2
  
    locs = DataFrame(vec([(;x=j+r(), y=i+r()) for i=1:a for j=1:b]))  # pozycje węzłów
    locs.x .*= scale
    locs.y .*= scale
    dist(n1,n2) = √((locs.x[n1]-locs.x[n2])^2+(locs.y[n1]-locs.y[n2])^2)  # funkcja licząca dystans

    g = T(nrow(locs))
    for i=1:a, j=1:b
        n = b*(i-1) + j
        j <= b-1 && add_edge2!(g,n,n+1,dist(n,n+1))
        i <= a-1 && add_edge2!(g,n,n+b,dist(n,n+b))
    end
    (;g,locs)
end

function neighbourhood(g, i, k)
    n_v = []
    for kk in 1:k
        n_v_l = []
        if kk == 1
            n_v_l = Graphs.neighbors(g,i)
        else
            for ii in n_v
                n_v_l = unique(vcat(n_v_l, Graphs.neighbors(g,ii)...))
            end
        end
        n_v = unique(vcat(n_v, n_v_l...))
    end
    n_v = n_v[n_v .!= i]
    return unique(n_v)
end

G = Graphs.erdos_renyi(30, 45)
typeof(G)

using Random
using DataFrames
using SimpleWeightedGraphs
G_l = rectangle_graph(5,5)

# n = 1

println(cgrad(:matter, 5, categorical = true)[1])

colors = fill(cgrad(:matter, 5, categorical = true)[5], Graphs.nv(G_l.g))
colors[neighbourhood(G_l.g, 1, 4)] .= cgrad(:matter, 5, categorical = true)[4]
colors[neighbourhood(G_l.g, 1, 3)] .= cgrad(:matter, 5, categorical = true)[3]
colors[neighbourhood(G_l.g, 1, 2)] .= cgrad(:matter, 5, categorical = true)[2]
colors[neighbourhood(G_l.g, 1, 1)] .= cgrad(:matter, 5, categorical = true)[1]
colors[1] = RGBA(1.,1.,1.)
gr1 = gplot(G_l.g, G_l.locs.x, G_l.locs.y, nodelabel = 1:Graphs.nv(G_l.g), nodefillc = colors, nodestrokec = "black", nodestrokelw  = 1)
draw(PDF(pwd() * "\\plots_to_export\\network_graph.pdf", 16cm, 16cm), gr1)

############## oczekiwana dalsza użyteczność

using Plots
using LaTeXStrings

H = 10
ts = 1
k = 0.5
d = 0.5
β = 0.75
ρ = 0.95

function utility_next(t, ts, H, β, k, d, ρ)
    @assert t >= ts
    @assert t <= ts + H
    β * k * d^(t-ts) * (1 - (ρ * d)^(H+ts-t)) / (1 - ρ * d)
end

[utility_next(t, ts, H, β, k, d, ρ) for t in 1:10]

eu_p = plot([utility_next(t, ts, H, β, k, 0.25, ρ) for t in 1:H], markershape = :circle, xlim = (0,H+1), xticks = collect(0:(H+1)), xlabel = "Czas, zakup w momencie t=1", ylabel = "Oczekiwana użyteczność z konsumpcji dobra", xlabelfontsize = 10, ylabelfontsize = 10, xtickfontsize = 6, ytickfontsize = 6, label = L"\hat D_{ijt} = 0.25")

plot!([utility_next(t, ts, H, β, k, 0.50, ρ) for t in 1:H], markershape = :circle, xlim = (0,H+1), xticks = collect(0:(H+1)), xlabel = "Czas, zakup w momencie t=1", ylabel = "Oczekiwana użyteczność z konsumpcji dobra", xlabelfontsize = 10, ylabelfontsize = 10, xtickfontsize = 6, ytickfontsize = 6, label = L"\hat D_{ijt} = 0.50")

plot!([utility_next(t, ts, H, β, k, 0.75, ρ) for t in 1:H], markershape = :circle, xlim = (0,H+1), xticks = collect(0:(H+1)), xlabel = "Czas, zakup w momencie t=1", ylabel = "Oczekiwana użyteczność z konsumpcji dobra", xlabelfontsize = 10, ylabelfontsize = 10, xtickfontsize = 6, ytickfontsize = 6, label = L"\hat D_{ijt} = 0.75")

savefig(eu_p, pwd() * "\\plots_to_export\\expected_further_utility.pdf")


############## cechy dobra a przekonania

using Plots
using LaTeXStrings
pyplot()

ρ = 0.95

function q(t, ts, H, k, d)
    @assert t >= ts
    @assert t <= ts + H
    k * d^(t-ts) * (1 - d^(H+ts-t)) / (1 - d)
end

function beta_req(ke,de;t,ts,H,k,d,c,m,ρ)
    return c * m * q(t,ts,H,k,d) / q(t,ts,H, ke*ρ, de*ρ)
end

c = 0.4
m = 1.1

q(1,1,10,0.5,0.5)
q(1,1,10,0.5,0.5)
beta_req(0.5,0.5;t=1,ts=1,H=10,k=0.5,d=0.5,c=c,m=m,ρ=ρ)

K = LinRange(0.30:0.010:0.70)
D = LinRange(0.30:0.010:0.70)

beta_min_p1 = contour(K,D, (K,D)->beta_req(K,D;t=1,ts=1,H=10,k=0.5,d=0.5,c=c,m=m,ρ=ρ), levels=20, xlabel = "Przekonanie o jakości dobra, " * L"K_{ijt}", ylabel = "Przekonanie o trwałości dobra, " * L"D_{ijt}", titlefontsize = 10, xlabelfontsize = 10, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6, clabels=true, color=:turbo, colorbar_title = L"\beta_i")
Plots.scatter!([0.5], [0.5], color = :black, legend = false)

savefig(beta_min_p1, pwd() * "\\plots_to_export\\beta_min_to_purchase.pdf")

#### probability of choosing

using Plots
using LaTeXStrings

β = 0.8
c = 0.5
m = 1.1

k1 = 0.85
k2 = 0.15

d1 = 0.85
d2 = 0.45

ρ = 0.85
H = 6

betas = LinRange(0.001:0.001:0.999)

surplus(β, k, d, H, c, m, ρ) = β * k * (1-(ρ*d)^H)/(1-ρ*d) - c*m * k * (1-(d)^H)/(1-d)

s1 = [surplus(β, k1, d1, H, c, m, ρ) for β in betas]
s2 = [surplus(β, k2, d2, H, c, m, ρ) for β in betas]

Plots.plot(betas, s1)
Plots.plot!(betas, s2)

s1 = max.(0, s1)
s2 = max.(0, s2)

calc_prob(x1,x2) = (x1 == 0) & (x2 == 0) ? 0 : x1 / (x1+x2)

beta_prob = Plots.plot(betas, calc_prob.(s1, s2), xlabel = L"\beta_i", ylabel = "Prawdopodobieństwo wyboru", legend = :topleft, label = "Firma 1, " * L"(K_{jt} = 0.85, D_{jt} = 0.85)", titlefontsize = 8, xlabelfontsize = 10, ylabelfontsize = 10, xtickfontsize = 6, ytickfontsize = 6, linewidth = 2, ylim = (-0.01, 1.01))
Plots.plot!(betas, calc_prob.(s2, s1), label = "Firma 2, " * L"(K_{jt} = 0.15, D_{jt} = 0.45)", linewidth = 2)

calc_prob.(s2, s1)[739]

savefig(beta_prob, pwd() * "\\plots_to_export\\beta_vs_probability.pdf")
savefig(beta_prob, pwd() * "\\thesis_plots\\beta_vs_probability.svg")
# kupowanie

using StatsBase

repeat(1:10, 5)

[sample(1:2, Weights([calc_prob.(s1, s2)[x], calc_prob.(s2, s1)[x]])) for x in repeat(500:999,5)]

#### lost opportunity

using Plots
using LaTeXStrings

β = 0.8
c = 0.4
m = 1.2

k1 = 0.85
k2 = 0.15

d1 = 0.85
d2 = 0.45

ρ = 0.95
H = 10

betas = LinRange(0.001:0.001:0.999)

surplus(β, k, d, H, c, m, ρ) = β * k * (1-(ρ*d)^H)/(1-ρ*d) - c*m * k * (1-(d)^H)/(1-d)

s1 = [surplus(β, k1, d1, H, c, m, ρ) for β in betas]
s2 = [surplus(β, k2, d2, H, c, m, ρ) for β in betas]

plot(s1)
plot!(s2)

s1 = max.(0, s1)
s2 = max.(0, s2)

calc_prob(x1,x2) = (x1 == 0) & (x2 == 0) ? 0 : x1 / (x1+x2)

beta_prob = Plots.plot(betas, calc_prob.(s1, s2), xlabel = "Cena rezerwacji, " * L"\beta_i", ylabel = "Prawdopodobieństwo wyboru", ylim = (0,1), legend = :topleft, label = "Firma 1, " * L"(K_{jt} = 0.85, D_{jt} = 0.85)", titlefontsize = 8, xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6, title = "Prawdopodobieństwo zakupu a cena rezerwacji")
Plots.plot!(betas, calc_prob.(s2, s1), label = "Firma 2, " * L"(K_{jt} = 0.15, D_{jt} = 0.45)")

savefig(beta_prob, pwd() * "\\Documents\\GitHub\\PhD-Model\\thesis_plots\\beta_vs_probability.pdf")


#### diff between producer resarch and lack of

using Plots
using LaTeXStrings

β = 0.8
c = 0.4
m = 1.2
ke = 0.5
de = 0.5
ρ = 0.95
H = 5

K = LinRange(0.45:0.001:0.55)
D = LinRange(0.45:0.001:0.55)

utility_diff(k, d; ke, de, H, ρ) = (β * ke * (1-(ρ*de)^H)/(1-ρ*de) - β * k * (1-(ρ*d)^H)/(1-ρ*d)) / β * ke * (1-(ρ*de)^H)/(1-ρ*de)

utility_diff(0.5, 0.5; ke = 0.5, de = 0.5, H = 5, ρ = 0.95)

kd_vs_kede = contour(K,D, (K,D)->utility_diff(K,D;ke = 0.5, de = 0.5, H=5 ,ρ=ρ), levels=20, xlabel = "Zakładana przez producenta jakość dobra, " * L"K_{jt}", ylabel = "Zakładana przez producenta trwałość dobra, " * L"D_{jt}", titlefontsize = 10, xlabelfontsize = 10, ylabelfontsize = 10, xtickfontsize = 6, ytickfontsize = 6, clabels=true, color=:turbo, colorbar_title = "Różnica [%] oczekiwanej użyteczności \n i użyteczności zakładanej przez producenta")
Plots.scatter!([0.5], [0.5], color = :black, legend = false)

savefig(kd_vs_kede, pwd() * "\\plots_to_export\\kd_vs_kede.pdf")

#### diff between producer resarch and lack of

H = 7
N = 1000
T = 100
Q = 50

(H-1)*(1-Q/N)^(H-1)*Q/N

US = zeros(Int64, N, T)
BN = zeros(Int64, N, T)

for i in 1:N
    for t in 1:T
        if t == 1
            if rand() < Q/N
                BN[i,t] = 1
                US[i,t] = 1
            end
        elseif 2 <= t <= H
            if US[i,t-1] == 1
                US[i,t] = 1
            else
                if rand() < Q/N
                    BN[i,t] = 1
                    US[i,t] = 1
                end
            end
        else
            if US[i,t-1] == 1
                if BN[i,t-H] == 1
                    if rand() < Q/N
                        BN[i,t] = 1
                        US[i,t] = 1
                    end
                else
                    US[i,t] = 1
                end
            else
                if rand() < Q/N
                    BN[i,t] = 1
                    US[i,t] = 1
                end
            end
        end            
    end
end



Plots.plot(vec(sum(US, dims=1))/N)

#####

using Plots
using LaTeXStrings
pyplot()


K1 = 0.50
K2 = 0.60

if_i = 1
if_w = 1



calculate_posterior(λ_i, λ_w; K1, K2) = K1*(1 - λ_i - λ_w) + K2 * (λ_i + λ_w)

li = LinRange(0.01:0.01:0.99)
lw = LinRange(0.01:0.01:0.99)

ac_p1 = contour(li, lw, (li,lw)->calculate_posterior(li,lw; K1=K1, K2=K2), levels=50, xlabel = "Szybkość uczenia się, własne doświadczenie, " * L"K_{jt}", ylabel = "Trwałość dobra, " * L"D_{jt}",title = "Całkowita użyteczność dla H = 1", titlefontsize = 8, xlabelfontsize = 8, ylabelfontsize = 8, xtickfontsize = 6, ytickfontsize = 6)

#### Clip trans

using Plots
using LaTeXStrings
pyplot()
binomial(5,2)
x = [-1,0,1,2]
clip(x) = min(max(0, x), 1)
y = clip.(x)

clip_p = Plots.plot(x,y, xlabel = "x", ylabel = "f(x)", legend = nothing, xlim = (-1,2))

savefig(clip_p, pwd() * "\\plots_to_export\\clip_transform.pdf")

#### Bates Distributions

using Plots
using LaTeXStrings

bates(x;n) = n/(2*factorial(n-1)) * sum((-1)^k*binomial(n,k)*(n*x-k)^(n-1)*sign(n*x-k) for k in 0:n)
moved_bates(x;n,a,b) = 1/(b-a) * bates((x-a)/(b-a);n)

a = 0.40
b = 0.60

bates_dist_p = Plots.plot(LinRange((a+0.001):0.001:(b-0.001)), [moved_bates(x;n=1, a=a, b=b) for x in LinRange((a+0.001):0.001:(b-0.001))], label = L"|v_k(i)| = 1", xlabel = "x", ylabel = "g(x; " * L"|v_k(i)|, K_{jt} - \epsilon^K, K_{jt} + \epsilon^K" *")")
Plots.plot!(LinRange((a+0.001):0.001:(b-0.001)), [moved_bates(x;n=2, a=a, b=b) for x in LinRange((a+0.001):0.001:(b-0.001))], label = L"|v_k(i)| = 2")
Plots.plot!(LinRange((a+0.001):0.001:(b-0.001)), [moved_bates(x;n=5, a=a, b=b) for x in LinRange((a+0.001):0.001:(b-0.001))], label = L"|v_k(i)| = 5")
Plots.plot!(LinRange((a+0.001):0.001:(b-0.001)), [moved_bates(x;n=10, a=a, b=b) for x in LinRange((a+0.001):0.001:(b-0.001))], label = L"|v_k(i)| = 10")

savefig(bates_dist_p, pwd() * "\\thesis_plots\\bates_dist.pdf")

#### Expectation walk


#### beta and ps

using Distributions
using Plots
using LaTeXStrings

beta = rand(TriangularDist(0,1,0.5), 100000)
ps = rand(TriangularDist(1,2,1.5), 100000)

p_beta = Plots.plot([0,1], [0,1], color = "black", legend = nothing, ylabel = L"\beta_i", xlabel = L"\frac{\beta_i}{\gamma_i}", margin=2.5Plots.mm)
Plots.plot!([0,0.5], [0,1.], color = "black")
Plots.annotate!(0.2, 0.55, L"\gamma_i = 2")
Plots.annotate!(0.45, 0.55, L"\gamma_i = 1")

Plots.savefig(p_beta, pwd() * "\\plots_to_export\\beta_ps.pdf")

####

include(pwd() * "\\methods\\methods.jl")

Random.seed!(12345)

@time sim_single = TO_GO(200, 2, 200, 400, [0.6, 0.6], [1.1, 1.1], "random", .8, .8, "stochastic", [[0.05, 0.95], [0.05, 0.95]], [[0.05, 0.95], [0.05, 0.95]], [[1., 2.], [1., 2.]], 0.10, true, 1, [0.85, 1.], "softmax", ["market research", "market research"], [0.1, 0.1], 3, "dist", TriangularDist(0,1,0.5))


quality_expectation_buyers = [[] for _ = 1:200] # to T
sim_single.buyers[1].unit_buying_selling_history
for b in sim_single.buyers
    ubsh = b.unit_buying_selling_history
    for item in ubsh
        if item.decision == "buy, primary market"
            qe = getindex.(b.quality_expectation_history, 1)[item.t]
            push!(quality_expectation_buyers[item.t], qe)
        end
    end
end
p = Plots.plot(xlabel = "T", ylabel = "Jakość / przekonanie o jakości", legend = :bottomleft, ylim = (0.2, 0.7))
for i in 1:200
    p = Plots.plot!(getindex.(sim_single.buyers[i].quality_expectation_history, 1), color = "grey", linealpha = 0.10, label = nothing)
end
p
Plots.plot!(sim_single.sellers[1].quality_history, label = "Średnia jakość, producent", linewidth = 2, color = "blue")
Plots.plot!(mean([getindex.(x,1) for x in getfield.(sim_single.buyers, :quality_expectation_history)]), color = "red", linewidth = 2, label = "Przekonanie o jakości, cała populacja")
Plots.plot!(mean_nothing.(quality_expectation_buyers), label = "Przekonanie o jakości, kupujący w t", linewidth = 2, color = "orange")

Plots.savefig(p, pwd() * "\\plots_to_export\\quality_expected_average.pdf")