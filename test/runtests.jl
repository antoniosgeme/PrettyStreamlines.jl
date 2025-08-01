using Test
using PrettyStreamlines  
using Statistics

@testset "Uniform flow streamlines" begin
    x = y = -2:0.5:2
    U = fill(1.0, length(y), length(x))
    V = fill(0.0, length(y), length(x))

    
    streams = compute_streamlines(x, y, U, V)
    @test !isempty(streams) 

    br = vcat(0,findall(isnan.(streams[:,1])))

    streamdiff = zeros(length(br)-1)
    for i = firstindex(br):lastindex(br)-1
        ys = streams[br[i]+1:br[i+1]-1,2]
        streamdiff[i]  = isapprox(sum(diff(ys)), 0.0; atol=1e-8)
    end
    @test all(streamdiff .== 1)
end

@testset "Rotation field preserves radius" begin
    x = y = -1:0.1:1
    X = [j for i in y, j in x]
    Y = [i for i in y, j in x]
    U = -Y
    V = X

    streams = compute_streamlines(x, y, U, V)
    br = vcat(0,findall(isnan.(streams[:,1])))

    radii = zeros(length(br)-1)
    for i = firstindex(br):lastindex(br)-1
        xs = streams[br[i]+1:br[i+1]-1,1]
        ys = streams[br[i]+1:br[i+1]-1,2]
        rs = sqrt.(xs.^2 .+ ys.^2)
        radii[i] = isapprox(sum(diff(rs)), 0.0; atol=1e-3)
    end
    @test all(radii .== 1)
end


@testset "Seed injection" begin
    x = y = -1:0.2:1
    X = [j for i in y, j in x]
    Y = [i for i in y, j in x]
    U = -Y
    V = X

    rc = [0.5,0.25]
    seed = [(rc[1], 0.0),(0.0,rc[2])]
    streams = compute_streamlines(x, y, U, V; seeds=seed)
    @test !isempty(streams)
    br = vcat(0,findall(isnan.(streams[:,1])))

    for i = firstindex(br):lastindex(br)-1
        xs = streams[br[i]+1:br[i+1]-1,1]
        ys = streams[br[i]+1:br[i+1]-1,2]
        rs = sqrt.(xs.^2 .+ ys.^2)
        rtest = rc[argmin(abs.(mean(rs) .- rc))]
        @test all(isapprox.(rs.- rtest, 0.0; atol=1e-5) )
    end 
end