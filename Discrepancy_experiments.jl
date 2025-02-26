import Pkg
#Pkg.add(url="https://github.com/Vishnupriya-Anupindi/ReducedDigitalNets.jl")
using ReducedDigitalNets, LinearAlgebra, CairoMakie


b = 2
s = 2
m = 4

#Identity matrix
C_1 = diagm(0 => fill(1,m))

################ Use Ctrl + / to uncomment for LP and Hammersley net ######
# C_2_1 = fill(1,(m,m))
# C_2a = triu(C_2_1)
# C_2 = C_2a[:, end:-1:1]

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

# #LP net
# C = [C_1,C_2]
# P = DigitalNetGenerator(b,m,s,C) 
# pts = genpoints(P)


# #Hammersley
# H = [C_1,C_1a]
# P_h = DigitalNetGenerator(b,m,s,H) 
# pts_h = genpoints(P_h)

# C_a = [C_1a,C_2a]
# P_a = DigitalNetGenerator(b,m,s,C_a) 
# pts_a = genpoints(P_a)

#pts == pts_a


# Pascal matrix
C_3 = zeros(Int64,m,m)
for i in 1:m
    for j in 1:m
        if i <= j
            C_3[i,j] = binomial(j-1, i-1) % b
        end
    end
end
return C_3


C_3*C_1a

M = zeros(Int64, m, m)
    for i in 1:m
        for j in 1:m
            M[i, j] = binomial(m-j, i-1) % 2
        end
    end
return M

C_1a*C_3

################################# Pascal net ##################################
CP = [C_1, C_3]
P_pa = DigitalNetGenerator(b,m,s,CP) 
pts_pa = sort(genpoints(P_pa))

#Reduced Pascal net
w = (0,1)
C3_red = zeros(Int64,m,m)
for j in 1:m-w[2]
    C3_red[:,j] = C_3[:,j]
end
C3_red
CP_red = [C_1,C3_red]
P_pa_red = DigitalNetGenerator(b,m,s,CP_red)
pts_pa_red = sort(genpoints(P_pa_red))

begin
    fig = Figure(resolution = (800, 400))

    ax = Axis(fig[1,1], title = "Some (0,$m,$s)-net", limits = (nothing,nothing, nothing, 1), xticks = [0.0, 0.25, 0.5, 0.75,1.0], yticks = [0.0, 0.25, 0.5, 0.75,1.0],xminorgridvisible = false, yminorticksvisible = false)
    scatter!( Point2.(pts_pa) )

    ax = Axis(fig[1,2], title = "Reduced net, w = $w", limits = (nothing,nothing, nothing, 1), xticks = [0.0, 0.25, 0.5, 0.75,1.0], yticks = [0.0, 0.25, 0.5, 0.75,1.0],xminorgridvisible = false, yminorticksvisible = false)
    scatter!( Point2.(pts_pa_red) )

    fig
end

# For LP-net and Hammersley

# begin
#     fig = Figure(resolution = (800, 400))

#     ax = Axis(fig[1,1], title = "Diagonal Upper-1-mat", limits = (nothing,nothing, nothing, 1), xticks = [0.0, 0.25, 0.5, 0.75,1.0], yticks = [0.0, 0.25, 0.5, 0.75,1.0],xminorgridvisible = false, yminorticksvisible = false)
#     scatter!( Point2.(pts) )

#     ax = Axis(fig[1,2], title = "Hammersley", limits = (nothing,nothing, nothing, 1), xticks = [0.0, 0.25, 0.5, 0.75,1.0], yticks = [0.0, 0.25, 0.5, 0.75,1.0],xminorgridvisible = false, yminorticksvisible = false)
#     scatter!( Point2.(pts_h) )

#     fig
# end