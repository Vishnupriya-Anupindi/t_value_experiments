import Pkg
#Pkg.add(url="https://github.com/Vishnupriya-Anupindi/ReducedDigitalNets.jl")
using ReducedDigitalNets, LinearAlgebra, CairoMakie, DataFrames, CSV, ProgressMeter
include("Discrepancy_utils.jl")
using Combinatorics

mkpath("data")
fn_postfix = "id"


begin
    b = 2
    s = 2
    m = 5
    N = b^m

    bf = float(b)
    pts = [zeros(s) for i in 1:N ]
    badic = get_badic(b,m)

    println("b = $b, m = $m, s = $s, # matrices = $(b^(m*m)), # zero row matrices = $(b^(m*(m-1)))" )
end

begin
    matrix_range = 0:b^(m*m)-1
    @assert length(matrix_range) <= b^(m*m) "matrix range large"
    @assert last(matrix_range) >= b^(m*(m-1)) "zero rows"

    if s==2
        C = [diagm(0 => fill(1,m)),zeros(BigInt,m,m)]
    
        df = DataFrame(i1 = Int64[], d_s = Float64[])
        CSV.write("data/disc_output_b$(b)_m$(m)_s$(s)$(fn_postfix).csv ", df)
    
        prog = Progress(Int(length(matrix_range)*(length(matrix_range)-1)/2))

        for i1 in matrix_range

            # write in the beginning, to ensure that we don't accidently shift the writing part
            if i1 % 100 == 0 
                CSV.write("data/disc_output_b$(b)_m$(m)_s$(s)$(fn_postfix).csv ", df, append = true)
                empty!(df)
            end

            
                next!(prog)
                
                C[2] .= reshape(int_to_matrix(i1, b, m), m, m)

                if det(C[2]) % b == 0
                    continue
                end

                if compare_rows_matrix(C...)
                    continue
                end

                if compare_rows_lin_ind(C...)
                    continue
                end
                    
                P = DigitalNetGenerator(b,m,s,C) 
                pts = genpoints(P) 
                
                d_s = discr(pts)
                push!(df, (i1, d_s))
        end
        CSV.write("data/disc_output_b$(b)_m$(m)_s$(s)$(fn_postfix).csv ", df, append = true)
    else 
        @warn "not implemented"
    end

end

df_result = CSV.read("data/disc_output_b$(b)_m$(m)_s$(s)$(fn_postfix).csv ", DataFrame)

ma = maximum(df_result[!,2])
mi = minimum(df_result[!, 2])

filter!(row -> row.d_s == mi, df_result)

CSV.write("data/disc_filter_b$(b)_m$(m)_s$(s)$(fn_postfix).csv ",df_result)

begin
    idxs = values(df_result[12,:])[1] # Use Vector() if number of row entries are large
    C_2 = [get_matrix(i,b,m) for i in idxs]
    J = hcat(C_2...)
    display(J)
    @show idxs

    C = [diagm(0 => fill(1,m)),zeros(BigInt,m,m)]
    C[2] .= reshape(int_to_matrix(idxs, b, m), m, m)
    P = DigitalNetGenerator(b,m,s,C) 
    pts = genpoints(P)
    #display(pts)

    include("NNLD_plots.jl")
    plot_points(pts)

end
det(C[2])

begin
    #guess_4 = 1 + 2^5 + 2^9 + 2^10 + 2^12 + 2^15
    #guess = 1 + 2^6 + 2^12 + 2^16 + 2^18 + 2^20 + 2^24
    guess = 1 + 2^6 + 2^12 + 2^18 + 2^20 + 2^24
    C = [diagm(0 => fill(1,m)),zeros(BigInt,m,m)]
    C[2] .= reshape(int_to_matrix(guess, b, m), m, m)
    display(C[2])
    P = DigitalNetGenerator(b,m,s,C) 
    pts = genpoints(P)
    d_s = discr(pts)
    display(d_s)
    plot_points(pts)
end
get_matrix(guess,b,m)


df_result_2 = CSV.read("data/disc_output_b$(b)_m$(m)_s$(s)$(fn_postfix).csv ", DataFrame)

# lp_disc = 0.25
lp_disc = ma

filter!(row -> row.d_s == lp_disc, df_result_2)

CSV.write("data/disc_filter_lp_b$(b)_m$(m)_s$(s)$(fn_postfix).csv ",df_result_2)

begin
    idxs = values(df_result_2[1,:])[1] # Use Vector() if number of row entries are large
    C_2 = [get_matrix(i,b,m) for i in idxs]
    J = hcat(C_2...)
    display(J)
    @show idxs

    C = [diagm(0 => fill(1,m)),zeros(BigInt,m,m)]
    C[2] .= reshape(int_to_matrix(idxs, b, m), m, m)
    P = DigitalNetGenerator(b,m,s,C) 
    pts = genpoints(P)
    #display(pts)

    include("NNLD_plots.jl")
    plot_points(pts)

end



idxs = values(df_result[1,:])[1] # Use Vector() if number of row entries are large
C_2 = get_matrix(idxs,b,m)
C_2 != [0 1 1; 1 1 1; 0 1 0]
i_fil = df[!,1]
i_fil[1]


for i in 1:length(i_fil)
    up_mat = Vector{Int64}()
    C_2 = get_matrix(i_fil[i],b,m)
    if C_2 == [1 1 1; 1 1 0; 1 0 0]
        push!(up_mat,i)
    end
    return up_mat
end

int_to_matrix(95, b, m) 
m