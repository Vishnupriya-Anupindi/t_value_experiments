import Pkg
#Pkg.add(url="https://github.com/Vishnupriya-Anupindi/ReducedDigitalNets.jl")
using ReducedDigitalNets, LinearAlgebra, CairoMakie


b = 2
s = 2
m = 8

#Identity matrix
C_1 = diagm(0 => fill(1,m))

C_2_1 = fill(1,(m,m))
C_2a = triu(C_2_1)
C_2 = C_2a[:, end:-1:1]

C_1a = C_1[:, end:-1:1]


# begin
#     C_2_1 = fill(1,(m,m))
#     C_2_t = triu(C_2_1)
#     C_2 = zeros(Int64, m,m)
#     for j in 1:m
#         C_2[:,j] = C_2_t[:,m-j+1] 
#     end
#     C_2
# end

C = [C_1,C_2]
P = DigitalNetGenerator(b,m,s,C) 
pts = genpoints(P)


#Hammersley
H = [C_1,C_1a]
P_h = DigitalNetGenerator(b,m,s,H) 
pts_h = genpoints(P_h)

# C_a = [C_1a,C_2a]
# P_a = DigitalNetGenerator(b,m,s,C_a) 
# pts_a = genpoints(P_a)

#pts == pts_a

begin
    fig = Figure(resolution = (800, 400))

    ax = Axis(fig[1,1], title = "Diagonal Upper-1-mat", limits = (nothing,nothing, nothing, 1), xticks = [0.0, 0.25, 0.5, 0.75,1.0], yticks = [0.0, 0.25, 0.5, 0.75,1.0],xminorgridvisible = false, yminorticksvisible = false)
    scatter!( Point2.(pts) )

    ax = Axis(fig[1,2], title = "Hammersley", limits = (nothing,nothing, nothing, 1), xticks = [0.0, 0.25, 0.5, 0.75,1.0], yticks = [0.0, 0.25, 0.5, 0.75,1.0],xminorgridvisible = false, yminorticksvisible = false)
    scatter!( Point2.(pts_h) )

    fig
end