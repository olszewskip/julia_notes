# © Przemyslaw Szufel, 2018, run under Julia 1.0.0

using PyPlot


function plot_schelling(aggField::Matrix{Int8},stepNo::Int,img)
    if img == nothing
        cla()
        ion()
        img = imshow(aggField)
    else
        img[:set_data](aggField)
    end
    title("Step $stepNo")
    show()
    sleep(0.1)
    return img
end
