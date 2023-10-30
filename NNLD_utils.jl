using Test

vol(z) = prod(z)
vol_h(z, pts::Vector) = count( all(p .< z) for p in pts ) / length(pts)
vol_h(z, pts::Matrix) = count( all(pts[j,i] < z[j] for j in axes(pts,1)) for i in axes(pts,2) ) / size(pts,2)
δ(z,pts) = vol_h(z, pts) - vol(z)

@test vol_h( [0.5,0.5], [0.4 0.6; 0.4 0.4] ) == 0.5
@test vol_h( [0.5,0.5], [[0.4,0.4],[0.6,0.4]] ) == 0.5

function is_NNLD(c_z, s, pts)
    NNLD = true
    z = zeros(s)
    for i in 1:c_z 
        rand!(z)
        if δ(z,pts) < 0
            NNLD = false
            break
        end
    end
    return NNLD
end

function is_NNLD_d(pts::Vector, tol = eps(1.0))
    NNLD = true
    z = similar(pts[1])
    for i in eachindex(pts)
        @. z = pts[i] - tol 
        @. z = max(z, 0.0)

        if δ(z,pts) < 0
            NNLD = false
            break
        end
    end
    return NNLD
end

@inline function norm_coord(v, b, bf = float(b))
    v_1 = 0.0
    for i in eachindex(v)
        v_1 += v[i] * bf^(-i)
    end
    return v_1
end

@test norm_coord([1 1 0 1],2) == 13/16

# generates the ith matrix (returned as a vector) 
# in base b with m rows and m columns 
function int_2_matrix!(C, i, b, m) 
    C .= 0
    if i > 0
        digits!(C, i, base=b)
    end
    return C
end 

function int_2_matrix(i, b, m) 
    C = zeros(Int, m*m)
    return int_2_matrix!(C, i, b, m)
end 

@test int_2_matrix(0, 3, 2) == [0, 0, 0, 0]
@test int_2_matrix(5, 3, 2) == [2, 1, 0, 0]
@test int_2_matrix(3^(2*2)-1, 3, 2) == [2, 2, 2, 2]


# pts needs to be of size (s, N) where N is length(badic) = b^m
function get_points!(pts, C, badic, m, b, bf = float(b))  
    Cn = zeros(Int64, m) 
    @inbounds for j in eachindex(pts)  # 1:N
        n = badic[j]

        for i in eachindex(C)
            mul!(Cn, C[i], n)
            for k in eachindex(Cn)
                Cn[k] = Cn[k] % b
            end 
            pts[j][i] = norm_coord(Cn, b, bf)
        end
    end
end

function get_points(C, badic, m, b, bf = float(b))  
    pts = [zeros(s) for i in 1:length(badic)]
    get_points!(pts, C, badic, m, b, bf)
    return pts
end

get_badic(b, m) = collect.(Iterators.product(fill(0:b-1, m)...))[:]

@testset begin "get points"
    C = ( [1 0; 0 1], [0 1; 1 0] )
    b = 2
    m = 2
    s = 2
    badic = get_badic(b,m)
    bf = float(b)
    pts = [zeros(s) for i in 1:length(badic)]
    get_points!(pts, C, badic, m, b, bf)
    @test pts == [[0.0, 0.0], [0.5, 0.25],[0.25; 0.5], [0.75, 0.75]]

    @test is_NNLD(100, 2, pts) == true
    @test is_NNLD_d(pts) == true
    #TODO add example with no NNLD
end

@testset begin "get points"
    C = ( [2 1 1; 1 1 0; 1 0 0], [0 1 1; 1 1 0; 1 0 0] )
    b = 3
    m = 3
    s = 2
    badic = get_badic(b,m)
    bf = float(b)
    pts = [zeros(s) for i in 1:length(badic)]
    get_points!(pts, C, badic, m, b, bf)
    #@test pts == [0 0.25 0.5 0.75; 0 0.5 0.25 0.75]

    @test is_NNLD(50, 2, pts) == false
    @test is_NNLD_d(pts) == false
    #TODO add example with no NNLD
end


# # t_value = 0 
# # for t_test in 1:m
# #     test_matrix = zeros(t_test, m)
# #     for d_1 in 0:m 
# #         for d_2 in 0:t_test-d_1 
# #             # d_1 + d_2 = m - t_test

# #             test_matrix[1:d_1,:] .= C1[1:d_1, :]
# #             test_matrix[d_1:d_1+d_2,:] .= C1[1:d_2, :]

# #             if rank(test_matrix) != m - t_test
# #                 break
# #             end
# #         end
# #     end
# #     t_value += 1
# # end