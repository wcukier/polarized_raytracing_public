

using Plots
using DelimitedFiles
using Statistics
using ImageFiltering
ENV["GKSwstype"]="nul"

dir = "sim"
pth = "/scratch/gpfs/wcukier/Polarized_Ray_Tracing/0_stored_data/$(dir)/image"
suffix = "_phase"
# suffix = ""

# n_phi = 1
# n_theta = 40

n_phi = 24 # Number of azimuthal bins
n_theta = 10 # Number of polar bins


PHI_AVG = false


function mask_data(aI, aQ, aU, j)
    nx = size(aI,2)

    if (~PHI_AVG)
        i = Int((j%n_phi)*10 + floor(j/n_phi))
        theta = (90 / n_theta) * (i % n_theta)
        phi = (360 / n_phi) * Int(floor(i / n_theta))
        mask = Int(1 + i*nx):Int((i+1)*nx)
        bI = aI[mask, :]
        bQ = aQ[mask, :]
        bU = aU[mask, :]
        bI = bI .+ 10^(-45)
    else # Phi average a 3d simulation
        bI = zeros(nx,nx) .+ 10^(-45)
        bQ = zeros(nx,nx)
        bU = zeros(nx,nx)

        for k in 0:n_phi-1
            i = Int(k*n_theta + floor(j/n_phi))
            mask = Int(1 + i*nx):Int((i+1)*nx)

            bI += aI[mask, :]
            bQ += aQ[mask, :]
            bU += aU[mask, :]

            phi = 0
            theta = (90 / n_theta) * (i % n_theta)
        end
    end
    
    return bI, bQ, bU, phi, theta, i
end

function plot_location!(theta, phi)
    θs = LinRange(0, π/2, n_theta+1)
    ϕs = LinRange(0,2π,n_phi+1)

    for θ in θs
        for ϕ in ϕs
            scatter!([7+2.5*sin(θ)* sin(ϕ)],[7+2.5*cos(θ)], color="blue", 
                        label=nothing)
        end
    end
    90 < phi < 270 ? c="magenta" : c="cyan"
    scatter!([7+2.5*sin((theta+8)/180*π)* sin(phi/180*π)],
                [7+2.5*cos((theta+8)/180*π)], color=c, label=nothing)



    annotate!(-9.5, 8, 
            text("$(Int(round(phi)))-$(Int(round(phi+(360 / n_phi))))", 
                :cyan, :left, 24))
    annotate!(-9.5, -8, 
                text("$(Int(round(theta)))-$(Int(round(theta+(90 / n_theta))))", 
                :red, :left, 24))
end

function median_filter(bI)
    patch_size = (3,3)
    bI = mapwindow(median, bI, patch_size)
    return bI
end


function stokes_I_plot(aI, aQ, aU, j, clims)
    nx = size(aI,2)
    n_polar = 20
    x_polar = LinRange(-10,10, Int(n_polar))
    ratio = Int(nx/n_polar)

    nx = size(aI,2)
    x = LinRange(-10,10, Int(nx))

    bI, bQ, bU, phi, theta, i = mask_data(aI, aQ, aU, j)
    P = .√(bQ.^2 + bU.^2)

    m = median_filter(P./bI)
    bI = median_filter(bI)


    # Plot the intensity
    heatmap(x, x, bI, size=(650,600), clims=clims, 
            background_color="transparent")

    plot_location!(theta, phi)

    Pₘ = min(maximum(P)*1/100, 10^(-5.5))
    plot_ticks!(bI, bQ, bU, m, Pₘ, x_polar, ratio, n_polar)
    return phi, theta, i
end

function mag_P_plot(aI, aQ, aU, j, clims)
    nx = size(aI,2)
    n_polar = 20
    x_polar = LinRange(-10,10, Int(n_polar))
    ratio = Int(nx/n_polar)

    nx = size(aI,2)
    x = LinRange(-10,10, Int(nx))

    bI, bQ, bU, phi, theta, i = mask_data(aI, aQ, aU, j)
    P = .√(bQ.^2 + bU.^2)

    m = median_filter(P./bI)
    bI = median_filter(bI)
    P = median_filter(P)


    # Plot the intensity
    heatmap(x, x, P, size=(650,600), clims=clims, 
            background_color="transparent")

    plot_location!(theta, phi)

    Pₘ = min(maximum(P)*1/100, 10^(-5.5))
    plot_ticks!(bI, bQ, bU, m, Pₘ, x_polar, ratio, n_polar)
    return phi, theta, i
end




function raytraced_image(aI, aQ, aU, nalpha, type, clims=(10^(-8), 10^(-5)))

    anim = @animate for j ∈ 0:(nalpha-1)
        phi, theta, i = stokes_I_plot(aI, aQ, aU, j, clims)
        saveplot("$(pth)/$(type)", "$(i)_theta-$(Int(round(theta)))_phi-$(Int(round(phi)))")
    end

    savegif(anim, pth, type)
end

function saveplot(path, name)
    try
        savefig("$(path)/$(name)")
    catch
        mkdir(path)
        savefig("$(path)/$(name)")
    end
end

function savegif(anim, path, name, fps=15)
    try
        gif(anim, "$(path)/gifs/$(name).gif", fps=fps)
    catch
        mkdir("$(path)/gifs")
        gif(anim, "$(path)/gifs/$(name).gif", fps=fps)
    end
end

function plot_ticks!(bI, bQ, bU, m, Pₘ, x_polar, ratio, n_polar)
    for j in (1):n_polar
        for k in 1:n_polar
            mask = [(k-1)*ratio+1:k*ratio,(j-1)*ratio+1:j*ratio]
            if (mean(bI[mask...])> Pₘ) && (sum(bI[mask...] .> 0)) > ratio
                s = max(mean(m[mask...]), 0.3)
            else
                 s = 0
            end

            χ= 1/2*atan(sum(bU[mask...]), sum(bQ[mask...]))
            plot!([x_polar[j] - sin(χ)/2*s, x_polar[j] + sin(χ)/2*s],
            [x_polar[k] + cos(χ)/2*s, x_polar[k] - cos(χ)/2*s],
            color="cyan",
            label=nothing,
            lw=3);
        end
    end
end


function Q_U_loops(aI, aQ, aU, nalpha, type)
    type *= "_QU"
    Qs = zeros(n_phi)
    Us = zeros(n_phi)
    phis = zeros(n_phi)
    thetas = zeros(n_phi)
    iϕ = 1
    for j ∈ 0:nalpha-1
        println(iϕ)
        I, bQ, bU, phi, theta, i = mask_data(aI, aQ, aU, j)
        Qs[iϕ] = sum(bQ) / sum(I)
        Us[iϕ] = sum(bU) / sum(I)
        phis[iϕ] = phi
        thetas[iϕ] = theta

        if iϕ==n_phi
            iϕ = 0
            print("phis: ")
            println(phis)

            print("thetas: ")
            println(thetas)
            
            plot(Qs, Us, label=nothing, lw=3, xlabel="Q/I", ylabel="U/I")
            println("saving $(pth)/$(type)/theta-$(Int(round(theta))).png")
            saveplot("$(pth)/$(type)", "theta-$(Int(round(theta))).png")
        end
        iϕ += 1

    end
end


function depolarization_hist(aI, aQ, aU, name)
    P = .√(aU.^2 + aQ.^2)
    m = (P[aI.!=0]./aI[aI.!=0])[:]
    if size(m)[1]>0
        histogram(m, yaxis = ((1, Inf)), label=nothing, bins=51)
        # histogram(m, yaxis = (:log10, (1, Inf)), label=nothing, bins=51)
        saveplot("$(pth)/hist", name)
    else
        println("$(pth)/hist not saved--no data")
    end
end


function EVPA_plot(aI, aQ, aU, j)
    nx = size(aI,2)
    x = LinRange(-10,10, Int(nx))
    n_polar = 20
    x_polar = LinRange(-10,10, Int(n_polar))
    ratio = Int(nx/n_polar)
    bI, bQ, bU, phi, theta, i = mask_data(aI, aQ, aU, j)
        
        
    m = median_filter((bU.^2 + bQ .^2)./bI)
    bI = median_filter(bI)
    bQ = median_filter(bQ)
    bU = median_filter(bU)

    P = bQ .+ im.*bU
    P[bI.< 1e-40] .= NaN
    χ = 1/2 .* atan.(bU , bQ)
    χ[bI .< 1e-40] .= NaN
    heatmap(x, x, χ, c=:phase, size=(650,600), background_color="transparent")
    plot_location!(theta, phi)
    Pₘ = 1e-20
    plot_ticks!(bI, bQ, bU, m, Pₘ, x_polar, ratio, 20)
    return phi, theta, i
end

function polarization_arg(aI, aQ, aU, nalpha, type)
    type = type * "_arg"




    anim = @animate for j ∈ 0:(nalpha-1)

        phi, theta, i = EVPA_plot(aI, aQ, aU, j)        
        saveplot("$(pth)/$(type)", "$(i)_theta-$(Int(round(theta)))_phi-$(Int(round(phi)))")

    end
    savegif(anim, pth, type)
end

function P_over_I_plot(aI, aQ, aU, j)
    nx = size(aI,2)
    x = LinRange(-10,10, Int(nx))
    bI, bQ, bU, phi, theta, i = mask_data(aI, aQ, aU, j)
    P = .√(bQ.^2 + bU.^2)

    m = median_filter(P./bI)
    bI = median_filter(bI)

    heatmap(x, x, m, size=(650,600), background_color="transparent")
    plot_location!(theta, phi)
    return phi, theta, i
end

function polar_degree_image(aI, aQ, aU, nalpha, type)
    type = type * "_degree"


    anim = @animate for j ∈ 0:(nalpha-1)

        phi, theta, i = P_over_I_plot(aI, aQ, aU, j)
        saveplot("$(pth)/$(type)", "$(i)_theta-$(Int(round(theta)))_phi-$(Int(round(phi)))")
    end

    savegif(anim, pth, type)
end




function plotimage(n::Int64, side::String, nalpha::Int64, f="1")

    a = readdlm("$(pth)/im_E_$(side)_$(n)$(suffix).dat")
    aQ = readdlm("$(pth)/im_Q_$(side)_$(n)$(suffix).dat")
    aU = readdlm("$(pth)/im_U_$(side)_$(n)$(suffix).dat")

    # Q_U_loops(a, f.*aQ, f.*aU, nalpha, "$(side)$(n)$(suffix)")
    polar_degree_image(a, f.*aQ, f.*aU, nalpha, "$(side)$(n)$(suffix)_01")
    raytraced_image(a, f.*aQ, f.*aU, nalpha,"$(side)$(n)$(suffix)_01")
    # depolarization_hist(a, aQ, aU, "010_$(n)")
    polarization_arg(a, f.*aQ, f.*aU, nalpha,"$(side)$(n)$(suffix)")
    print(n)
end

function plotall(nalpha)

    side = "cs"
    n = 1
    a = readdlm("$(pth)/im_E_$(side)_$(n)$(suffix).dat")
    aQ = readdlm("$(pth)/im_Q_$(side)_$(n)$(suffix).dat")
    aU = readdlm("$(pth)/im_U_$(side)_$(n)$(suffix).dat")

    for n ∈ 2:3
        a = a .+ readdlm("$(pth)/im_E_$(side)_$(n)$(suffix).dat")
        aQ = aQ .+ readdlm("$(pth)/im_Q_$(side)_$(n)$(suffix).dat")
        aU = aU .+ readdlm("$(pth)/im_U_$(side)_$(n)$(suffix).dat")
    end

    side = "pc"

    for n ∈ 1:3
        a = a .+ readdlm("$(pth)/im_E_$(side)_$(n)$(suffix).dat")
        aQ = aQ .+ readdlm("$(pth)/im_Q_$(side)_$(n)$(suffix).dat")
        aU = aU .+ readdlm("$(pth)/im_U_$(side)_$(n)$(suffix).dat")
    end


    raytraced_image(a, aQ, aU, nalpha, "all$(suffix)")
end

function plotturningpts(nalpha)
    a = readdlm("$(pth)/im_E_d$(suffix).dat")
    aQ = readdlm("$(pth)/im_Q_d$(suffix).dat")
    aU = readdlm("$(pth)/im_U_d$(suffix).dat")

    polar_degree_image(a, aQ, aU, nalpha, "direct$(suffix)")
    raytraced_image(a, aQ, aU, nalpha, "direct$(suffix)", (10^(-9), 10^(-6)))


    a = readdlm("$(pth)/im_E_id$(suffix).dat")
    aQ = readdlm("$(pth)/im_Q_id$(suffix).dat")
    aU = readdlm("$(pth)/im_U_id$(suffix).dat")

    polar_degree_image(a, aQ, aU, nalpha, "indirect$(suffix)")
    raytraced_image(a, aQ, aU, nalpha, "indirect$(suffix)", (10^(-9), 10^(-6)))
end

function plot_nturns(nalpha)
    a = readdlm("$(pth)/im_T$(suffix).dat")
    aQ = zeros(size(a))
    aU = zeros(size(a))

    polar_degree_image(a, aQ, aU, nalpha, "turns")
    raytraced_image(a, aQ, aU, nalpha, "turns", (0.0,4.0))

end


function julia_main()::Cint
    # plotall(n_theta*n_phi)
    # plotturningpts(n_theta*n_phi)
    # plot_nturns(n_theta*n_phi)
    plotimage(1, "cs", n_theta*n_phi, -1)
    plotimage(2, "cs", n_theta*n_phi, -1)
    plotimage(3, "cs", n_theta*n_phi, -1)
    # plotimage(1, "pc", n_theta*n_phi)
    # plotimage(2, "pc", n_theta*n_phi)
    # plotimage(3, "pc", n_theta*n_phi)
    return 0
end

# julia_main()





