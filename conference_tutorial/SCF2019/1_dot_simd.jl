using BenchmarkTools

function dot1(x, y)
    s = 0.0
    for i in 1:length(x)
        @inbounds s += x[i]*y[i]
    end
    s
end

function dot2(x, y)
    s = 0.0
    @simd for i in 1:length(x)
        @inbounds s += x[i]*y[i]
    end
    s
end


x = rand(100)
y = rand(100)

@btime dot1($x, $y)
@btime dot2($x, $y)

