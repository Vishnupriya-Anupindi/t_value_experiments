using Random, DataFrames, CSV, LinearAlgebra, ProgressMeter
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
    C = (reshape(C1, m, m), reshape(C2, m, m))

    bf = float(b)
    c_z = 60
    pts = [zeros(s) for i in 1:N ]
    badic = get_badic(b,m)


    df = DataFrame(i1 = Int64[], i2 = Int64[])
    CSV.write("output.csv", df)
    
    matrix_range = 0:1000 # 0:b^(m*m)-1

    @showprogress for i1 in matrix_range

        # write in the beginning, to ensure that we don't accidently shit the writing part
        if i1 % 10 == 0 
            CSV.write("output.csv", df, append = true)
            empty!(df)
        end

        int_2_matrix!(C1, i1, b, m)
        C[1] .= reshape(C1, m, m)

        if det(C[1]) % b == 0
            continue
        end

        for i2 in 0:i1-1
            int_2_matrix!(C2, i2, b, m)
            C[2] .= reshape(C2, m, m)

            if det(C[2]) % b == 0
                continue
            end

            if C[1][1,:] == C[2][1,:]
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

            nnld = is_NNLD_d(pts)
            if nnld == true
                push!(df, (i1, i2))
            end
        end
    end
end



# df_result = CSV.read("output.csv", DataFrame)

# use this to get the corresponding matrices
get_matrix(i, b, m) = reshape(int_2_matrix(i, b, m), m, m)
get_matrices(i1,i2,b,m) = [ get_matrix(i1,b,m) get_matrix(i2,b,m) ]


C1 = get_matrix(828, 3, 3)
C2 = get_matrix(820, 3, 3)

pts = get_points((C1, C2), badic, m, b)

is_NNLD_d(pts)



include("NNLD_plots.jl")

plot_points(pts)