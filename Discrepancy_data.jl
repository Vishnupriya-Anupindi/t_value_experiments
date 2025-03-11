import Pkg
#Pkg.add(url="https://github.com/Vishnupriya-Anupindi/ReducedDigitalNets.jl")
using ReducedDigitalNets, LinearAlgebra, CairoMakie, DataFrames, CSV, ProgressMeter

mkpath("data")
fn_postfix = "id"

begin
    b = 2
    s = 2
    m = 3
    N = b^m

    bf = float(b)
    pts = [zeros(s) for i in 1:N ]
    badic = get_badic(b,m)

    println("b = $b, m = $m, s = $s, # matrices = $(b^(m*m)), # zero row matrices = $(b^(m*(m-1)))" )
end


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
    return maximum(disc)
end

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

get_badic(b, m) = collect.(Iterators.product(fill(0:b-1, m)...))[:]
get_matrix(i, b, m) = reshape(int_to_matrix(i, b, m), m, m)

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

begin
    matrix_range = 0:b^(m*m)-1
    @assert length(matrix_range) <= b^(m*m) "matrix range large"
    @assert last(matrix_range) >= b^(m*(m-1)) "zero rows"

    if s==2
        C = [diagm(0 => fill(1,m)),zeros(Int,m,m)]
    
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

mi = minimum(df_result[!, 2])

filter!(row -> row.d_s == mi, df_result)

CSV.write("data/disc_filter_b$(b)_m$(m)_s$(s)$(fn_postfix).csv ",df_result)

begin
    idxs = values(df_result[8,:])[1] # Use Vector() if number of row entries are large
    C_2 = [get_matrix(i,b,m) for i in idxs]
    J = hcat(C_2...)
    display(J)
    @show idxs

    C = [diagm(0 => fill(1,m)),zeros(Int,m,m)]
    C[2] .= reshape(int_to_matrix(idxs, b, m), m, m)
    P = DigitalNetGenerator(b,m,s,C) 
    pts = genpoints(P)
    display(pts)

    include("NNLD_plots.jl")
    plot_points(pts)

end

idxs = values(df_result[13,:])[1] # Use Vector() if number of row entries are large
C_2 = get_matrix(idxs,b,m)
C_2 != [0 1 1; 1 1 1; 0 1 0]
i_fil = df_result[!,1]
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