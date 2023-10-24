using Random, Permutations, DataFrames, CSV, LinearAlgebra, ProgressMeter
include("NNLD_utils.jl")



nnld_counter = 0

#@profview let 
let 
    b = 3
    m = 3
    s = 2
    ρ = m 
    t = m - ρ 
    N = b^m

    C1 = zeros(Int64, m*m)
    C2 = zeros(Int64, m*m)

    bf = float(b)
    c_z = 60
    pts = zeros(s, N)
    badic = collect.(Iterators.product(fill(0:b-1, m)...))[:]


    df = DataFrame(i1 = Int64[], i2 = Int64[])
    CSV.write("output.csv", df)
    
    matrix_range = 0:1000 # 0:b^(m*m)-1

    @showprogress for i1 in matrix_range
        for i2 in 0:i1-1
            int_2_matrix!(C1, i1, b, m)
            int_2_matrix!(C2, i2, b, m)

            C = (reshape(C1, m, m), reshape(C2, m, m))
            if det(C[1]) % b == 0 || det(C[2]) % b == 0
                continue
            end


            get_points!(pts, C, badic, m, b, bf)  
            
            # L = 0 
            # for j in axes(pts, 2)
            #     if j > 1 && @views pts[:,j] == pts[:,1]
            #         break 
            #     else 
            #         L += 1
            #     end 
            # end 

            nnld = is_NNLD(c_z, s, pts)
            if nnld == true
                push!(df, (i1, i2))
            end
        end

        if i1 % 10 == 0 
            CSV.write("output.csv", df, append = true)
            empty!(df)
        end
    end
end

# df_result = CSV.read("output.csv", DataFrame)

# use this to get the corresponding matrices
get_matrix(i, b, m) = reshape(int_2_matrix(i, b, m), m, m)
get_matrices(i1,i2,b,m) = [ get_matrix(i1,b,m) get_matrix(i2,b,m) ]
