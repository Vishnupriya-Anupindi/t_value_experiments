using Test


get_badic(b, m) = collect.(Iterators.product(fill(0:b-1, m)...))[:]
get_matrix(i, b, m) = reshape(int_to_matrix(i, b, m), m, m)


vol(z) = prod(z)
vol_h(z, pts::Vector) = count( all(p .< z) for p in pts ) / length(pts)
vol_h_cl(z, pts::Vector) = count( all(p .<= z) for p in pts ) / length(pts)
δ(z,pts) = vol_h(z, pts) - vol(z)
δ_cl(z,pts) = vol_h_cl(z, pts) - vol(z)

# function discr(pts::Vector)
#     N = length(pts) #sort.by=()
#     disc = Vector{Float64}()
#     pts = sort(pts;by=x->x[1]) #sort!(pts)
#     for i in eachindex(pts)
#         sort_pts = sort(pts[1:i];by=x->x[2])
#         for j in 1:i
#             z = ((i-1)/N,sort_pts[j][2])
#             push!(disc, abs(δ(z,pts)), abs(δ_cl(z,pts)))
#             z_1 = ((i-1)/N, 1)
#             push!(disc, abs(δ(z_1,pts)), abs(δ_cl(z_1,pts)))
#         end
#     end
#     return maximum(disc)
# end

function discr(pts::Vector)
    N = length(pts) #sort.by=()
    disc = Vector{Float64}()
    pts = sort(pts;by=x->x[1]) #sort!(pts)
    for i in eachindex(pts)
        for j in 1:N
            z = ((i-1)/N,(j-1)/N)
            push!(disc, abs(δ(z,pts)), abs(δ_cl(z,pts)))
            z_1 = ((i-1)/N, 1)
            push!(disc, abs(δ(z_1,pts)), abs(δ_cl(z_1,pts)))
        end
    end
    return maximum(disc)
end

# begin
#     b = 2
#     s = 2
#     m = 3
#     N = b^m

#     i = 3
#     j = N
#     z = ((i-1)/N,(j-1)/N)

# end


# generates the ith matrix (returned as a vector) 
# in base b with m rows and m columns 
function int_to_matrix!(C, i, b, m) 
    C .= 0
    if i > 0
        digits!(C, i, base=b)
    end
    return C
end 

function int_to_matrix(i, b, m) 
    C = zeros(Int, m*m)
    return int_to_matrix!(C, i, b, m)
end 


function compare_rows_matrix(C1,C2)
    C = (C1,C2)
    if C[1][1,:] == C[2][1,:]
        return true
    end

    if C[1][1,:] == C[2][2,:]
        return true
    end

    if C[1][2,:] == C[2][1,:]
        return true
    end
    return false
end

function compare_rows_lin_ind(C1,C2)
    if m ==3
        M_1 = zeros(Int64, m,m)
        M_1[1,:] = C[1][1,:]
        M_1[2,:] = C[1][2,:]
        M_1[3,:] = C[2][1,:]
        M_1
        det(M_1) % b == 0
    
        M_2 = zeros(Int64, m,m)
        M_2[1,:] = C[2][1,:]
        M_2[2,:] = C[2][2,:]
        M_2[3,:] = C[1][1,:]
        M_2
        det(M_2) % b == 0
    
        if det(M_1) % b == 0
            return true
        end
    
        if det(M_2) % b == 0
            return true
        end
        return false
    end
end


# for partitions
# d = collect(partitions(m, s))
# d[1][1]
# d[1][2]