using Test

vol(z) = prod(z)
vol_h(z, pts::Vector) = count( all(p .< z) for p in pts ) / length(pts)
#vol_h(z, pts::Matrix) = count( all(pts[j,i] < z[j] for j in axes(pts,1)) for i in axes(pts,2) ) / size(pts,2)
δ(z,pts) = vol_h(z, pts) - vol(z)

#@test vol_h( [0.5,0.5], [0.4 0.6; 0.4 0.4] ) == 0.5
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

@test int_to_matrix(0, 3, 2) == [0, 0, 0, 0]
@test int_to_matrix(5, 3, 2) == [2, 1, 0, 0]
@test int_to_matrix(3^(2*2)-1, 3, 2) == [2, 2, 2, 2]


# pts needs to be of size (s, N) where N is length(badic) = b^m
function get_points!(pts, C, badic, b, m,bf = float(b))  
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

function get_points(C, badic, b, m, s, bf = float(b))  
    pts = [zeros(s) for i in 1:length(badic)]
    get_points!(pts, C, badic, b, m, bf)
    return pts
end

get_badic(b, m) = collect.(Iterators.product(fill(0:b-1, m)...))[:]

@testset "get points" begin
    C = ( [1 0; 0 1], [0 1; 1 0] )
    b = 2
    m = 2
    s = 2
    badic = get_badic(b,m)
    bf = float(b)
    pts = [zeros(s) for i in 1:length(badic)]
    get_points!(pts, C, badic, b, m, bf)
    @test pts == [[0.0, 0.0], [0.5, 0.25],[0.25; 0.5], [0.75, 0.75]]

    @test is_NNLD(100, 2, pts) == true
    
    #TODO add example with no NNLD
end

@testset "get points" begin
    C = ([2 1 1; 1 1 0; 1 0 0], [0 1 1; 1 1 0; 1 0 0])
    b = 3
    m = 3
    s = 2
    badic = get_badic(b,m)
    bf = float(b)
    pts = [zeros(s) for i in 1:length(badic)]
    get_points!(pts, C, badic, b, m, bf)
    #@test pts == [0 0.25 0.5 0.75; 0 0.5 0.25 0.75]

    @test is_NNLD(100, 2, pts) == false
    
    #TODO add example with no NNLD
end

get_matrix(i, b, m) = reshape(int_to_matrix(i, b, m), m, m)
#get_matrices(i1,i2,b,m) = [ get_matrix(i1,b,m) get_matrix(i2,b,m) ]

function validate_NNLD(idxs,c_z,b,m,s)
    C = [get_matrix(i,b,m) for i in idxs]
    pts = get_points(C, get_badic(b,m), b, m,s)

    return is_NNLD(c_z,s,pts)
end

@testset "validate NNLD" begin
    @test validate_NNLD((3354,3062),10000,3,3,2) == true
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

function compare_rows_matrix(C1,C2,C3)
    C = (C1,C2,C3)
    for d1 in 1:3
        for d2 in 1:3
            if d1<d2 && C[d1][1,:] == C[d2][1,:]
                return true
            end

            #if d1!=d2 && C[d1][1,:] == C[d2][2,:]
            #    return true
            #end
        end
    end
    return false
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

function is_d_NNLD(pts::Vector,s)
    N = length(pts) #sort.by=()
    pts = sort(pts;by=x->x[1]) #sort!(pts)
    for i in eachindex(pts)
        sort_pts = sort(pts[1:i];by=x->x[2])
        for j in 1:i
            if s==2
                #@show (j-1) i*sort_pts[j][2] 
                if (j-1) < i*sort_pts[j][2]  # prod(sort_pts[j][k] for k in 2:s)
                    return false
                end

            elseif s==3
                for k in 1:j
                    #@show (k-1) i*prod(sort_pts[k][d] for d in 2:s)
                    #there are bugs
                    sort_pts_2 = sort(sort_pts[1:j];by=x->x[3])
                    if (k-1) < i*prod(sort_pts_2[k][d] for d in 2:s)
                        @show (k-1) i*prod(sort_pts_2[k][d] for d in 2:s)
                        return false
                    end 
                end
                
            else 
                @error "not implemented"
            end
        end
    end 
    return true

end

#@testset "det_NNLD" 
begin
    (b,m,s) = (2,3,2)
    idxs = (273,140)
    C = [get_matrix(i,b,m) for i in idxs]
    pts = get_points(C, get_badic(b,m), b, m,s)
    is_d_NNLD(pts,s),is_NNLD(10000,s,pts)
end

begin
    (b,m,s) = (2,3,3)
    idxs = (273,106,84)
    C = [get_matrix(i,b,m) for i in idxs]
    pts = get_points(C, get_badic(b,m), b, m,s)
    is_d_NNLD(pts,s),is_NNLD(10000,s,pts)
end
