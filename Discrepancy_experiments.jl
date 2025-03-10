import Pkg
#Pkg.add(url="https://github.com/Vishnupriya-Anupindi/ReducedDigitalNets.jl")
using ReducedDigitalNets, LinearAlgebra, CairoMakie, DataFrames, CSV

b = 2
s = 2
m = 2

vol(z) = prod(z)
vol_h(z, pts::Vector) = count( all(p .< z) for p in pts ) / length(pts)
vol_h_cl(z, pts::Vector) = count( all(p .<= z) for p in pts ) / length(pts)
δ(z,pts) = vol_h(z, pts) - vol(z)
δ_cl(z,pts) = vol_h_cl(z, pts) - vol(z)

function discr(pts::Vector)
    N = length(pts) #sort.by=()
    disc = Vector{Float64}()
    pts = sort(pts;by=x->x[1]) #sort!(pts)
    for i in eachindex(pts)
        sort_pts = sort(pts[1:i];by=x->x[2])
        for j in 1:i
            z = (i/N,sort_pts[j][2])
            push!(disc, abs(δ(z,pts)), abs(δ_cl(z,pts)))
        end
    end
    return disc
    #d = maximum(disc)
    #return d
end



#Identity matrix
C_1 = diagm(0 => fill(1,m))

################ Use Ctrl + / to uncomment for LP and Hammersley net ######
C_2_1 = fill(1,(m,m))
C_2a = triu(C_2_1)
C_2 = C_2a[:, end:-1:1]

C_1a = C_1[:, end:-1:1]
#Inverse of C_1a = C_1a

# begin
#     C_2_1 = fill(1,(m,m))
#     C_2_t = triu(C_2_1)
#     C_2 = zeros(Int64, m,m)
#     for j in 1:m
#         C_2[:,j] = C_2_t[:,m-j+1] 
#     end
#     C_2
# end

#LP net
C = [C_1,C_2]
P = DigitalNetGenerator(b,m,s,C) 
pts = genpoints(P)


#Hammersley
H = [C_1,C_1a]
P_h = DigitalNetGenerator(b,m,s,H) 
pts_h = genpoints(P_h)

# For LP-net and Hammersley

d = discr(pts)
maximum(d)
d_h = discr(pts_h)
maximum(d_h)

begin
    x = LinRange(0,1,100)
    y = LinRange(0,1,100)
    fig = Figure(resolution = (1200, 800))

    ax = Axis(fig[1,1], title = "volume_interval")
    B = [a*b for a in x, b in y]
    heatmap!(x,y, B,colormap = "thermometer" )

    ax = Axis(fig[1,2], title="Upper-1-mat_volume_h(z)")
    A = [vol_h((a,b),pts[1:end]) for a in x, b in y]
    heatmap!(x,y, A, colormap = "thermometer" )

    ax = Axis(fig[1,3], title = "Upper-1-mat_delta(z)")
    hm = heatmap!(x,y, A-B, colormap = "thermometer" )
    scatter!( Point2.(pts), color = :black, markersize = 5 )
    Colorbar(fig[1,4], hm)

    ax = Axis(fig[2,1], title = "volume_interval")
    B = [a*b for a in x, b in y]
    heatmap!(x,y, B,colormap = "thermometer" )

    ax = Axis(fig[2,2], title="Hammersley_volume_h(z)")
    A = [vol_h((a,b),pts_h[1:end]) for a in x, b in y]
    heatmap!(x,y, A, colormap = "thermometer" )

    ax = Axis(fig[2,3], title = "Hammersley_delta(z)")
    hm = heatmap!(x,y, A-B, colormap = "thermometer" )
    scatter!( Point2.(pts_h), color = :black, markersize = 5 )
    Colorbar(fig[2,4], hm)
    fig
end


# begin
#     fig = Figure(resolution = (800, 400))

#     ax = Axis(fig[1,1], title = "Diagonal Upper-1-mat", limits = (nothing,nothing, nothing, 1), xticks = [0.0, 0.25, 0.5, 0.75,1.0], yticks = [0.0, 0.25, 0.5, 0.75,1.0],xminorgridvisible = false, yminorticksvisible = false)
#     scatter!( Point2.(pts) )

#     ax = Axis(fig[1,2], title = "Hammersley", limits = (nothing,nothing, nothing, 1), xticks = [0.0, 0.25, 0.5, 0.75,1.0], yticks = [0.0, 0.25, 0.5, 0.75,1.0],xminorgridvisible = false, yminorticksvisible = false)
#     scatter!( Point2.(pts_h) )

#     fig
# end


# # Helper function for loading Sobol matrices
# function row_to_mat(row,b,m)
#     return stack(first(digits(x, base = b, pad = m),m) for x in first(row,m))
# end

# function load_seq_mat(filename,b,m,s)
#     df = CSV.read(filename, DataFrame, header = false, delim = ' ')
#     K = Matrix(df[1:s,1:m])
#     C = [row_to_mat(K[i,:],b,m) for i in 1:s]
#     return C
# end

# ################################# Pascal net ##################################
# CP = [C_1a, C_3]
# P_pa = DigitalNetGenerator(b,m,s,CP) 
# pts_pa = sort(genpoints(P_pa))

# #Reduced Pascal net
# w = (0,1)
# C3_red = zeros(Int64,m,m)
# for j in 1:m-w[2]
#     C3_red[:,j] = C_3[:,j]
# end
# C3_red
# CP_red = [C_1a,C3_red]
# P_pa_red = DigitalNetGenerator(b,m,s,CP_red)
# pts_pa_red = sort(genpoints(P_pa_red))


# C_a = [C_1a,C_2a]
# P_a = DigitalNetGenerator(b,m,s,C_a) 
# pts_a = genpoints(P_a)

#pts == pts_a

# C_sob = load_seq_mat("sobol_Cs.txt",b, m, s)

# C_sob[2]

# # Pascal matrix
# C_3 = zeros(Int64,m,m)
# for i in 1:m
#     for j in 1:m
#         if i <= j
#             C_3[i,j] = binomial(j-1, i-1) % b
#         end
#     end
# end
# return C_3


# C_3*C_1a

# M = zeros(Int64, m, m)
#     for i in 1:m
#         for j in 1:m
#             M[i, j] = binomial(m-j, i-1) % 2
#         end
#     end
# return M

# C_1a*C_3

# begin
#     fig = Figure(resolution = (400, 400))

#     ax = Axis(fig[1,1], title = "LP-net", limits = (nothing,nothing, nothing, 1), xticks = [0.0, 0.25, 0.5, 0.75,1.0], yticks = [0.0, 0.25, 0.5, 0.75,1.0],xminorgridvisible = false, yminorticksvisible = false)
#     scatter!( Point2.(pts) )

#     # ax = Axis(fig[1,2], title = "Reduced net, w = $w", limits = (nothing,nothing, nothing, 1), xticks = [0.0, 0.25, 0.5, 0.75,1.0], yticks = [0.0, 0.25, 0.5, 0.75,1.0],xminorgridvisible = false, yminorticksvisible = false)
#     # scatter!( Point2.(pts_pa_red) )

#     fig
# end