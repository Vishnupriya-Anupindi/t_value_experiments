using CairoMakie

function plot_points(pts)
    x = LinRange(0,1,100)
    y = LinRange(0,1,100)
    fig = Figure(resolution = (1200, 400))
    
    ax = Axis(fig[1,1], title = "volume_interval")
    B = [a*b for a in x, b in y]
    heatmap!(x,y, B, colormap = "thermometer" )

    ax = Axis(fig[1,2], title="volume_h(z)")
    A = [vol_h((a,b),pts[1:end]) for a in x, b in y]
    heatmap!(x,y, A, colormap = "thermometer" )
    
    ax = Axis(fig[1,3], title = "delta(z)")
    hm = heatmap!(x,y, A-B, colormap = "thermometer" )
    scatter!( Point2.(pts), color = :black, markersize = 5 )
    Colorbar(fig[1,4], hm)
    fig    
end