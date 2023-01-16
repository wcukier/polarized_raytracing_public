using Plots
using DelimitedFiles
using Statistics
using ImageFiltering
using OffsetArrays
using LinearAlgebra
ENV["GKSwstype"]="nul"

# dir = "monopole_a0"
# pth = "/scratch/gpfs/wcukier/Polarized_Ray_Tracing/0_stored_data/$(dir)/image"
# suffix = "_phase"
# suffix = "id"

function azimuthalfourier(aI, aQ, aU, nx, nφ=100, mₘᵢₙ=-4, mₘₐₓ=4, ρₘᵢₙ=-10, ρₘₐₓ=100)
    x = LinRange(-10, 10, nx)
    yy = (x * ones(nx)')
    xx = (x * ones(nx)')'

    φφ = atan.(xx, yy) # φ is defined "East of North" (Palumbo 2020)
    ρρ = .√(xx.^2 + yy.^2)

    ρ_mask = ρₘᵢₙ .< .√(xx.^2 + yy.^2) .< ρₘₐₓ

    dφ = 2π/nφ
    Iₐₙₙ = 0
    
    β = offset_beta_array(mₘᵢₙ, mₘₐₓ)
    # for i in 1:nφ
    for i in 1:1
        # mask =  ((i-1)*dφ .< φφ .< i*dφ) .* ρ_mask
        mask = ρ_mask
        Iₐₙₙ += sum(aI[mask])
        for m ∈ mₘᵢₙ:mₘₐₓ
            # summation = (sum(aQ[mask]) + im*sum(aU[mask]))
            # beta = exp(im * 2 * (i-1/2)*dφ)*exp(-im * m * (i-1/2)*dφ)
            beta = sum((aQ[mask] .+ im.*aU[mask]) .* exp.(-im * m * φφ[mask]))
            # beta = sum(exp.(im .* 2 .* φφ[mask]) .* exp.(-im .* m .* φφ[mask])) 

            
            # beta = phase
            β[m] += beta
            # β[-m] += conj(beta/im)

        end
    end
    println(Iₐₙₙ)
    β ./= Iₐₙₙ

    β[norm.(β) .< 1e-7] .= 0 + im*0 
    return β
end

function offset_beta_array(mₘᵢₙ, mₘₐₓ)
    return OffsetArray(zeros(mₘₐₓ-mₘᵢₙ+1)+im*zeros(mₘₐₓ-mₘᵢₙ+1), mₘᵢₙ:mₘₐₓ)
end

function offset_beta_array(arr)
    mₘₐₓ = 4
    mₘᵢₙ = -4
    β = OffsetArray(zeros(mₘₐₓ-mₘᵢₙ+1)+im*zeros(mₘₐₓ-mₘᵢₙ+1), mₘᵢₙ:mₘₐₓ)
    for i in -4:4
        β[i] = arr[i+5]
    end
    # println(β)
    return β
end

function generate_fakedata_arr(arr)
    return generate_fakedata(offset_beta_array(arr), 1000)
end

function test_fakedata(arr)
    β = azimuthalfourier(generate_fakedata_arr(arr)..., 1000, 1000)
    plot_beta_ring(β, "test.png")
    return β
end


function generate_fakedata(β, nx=200, ρₘᵢₙ=5.5, ρₘₐₓ=5.75)
    mₘᵢₙ = β.offsets[1]+1
    mₘₐₓ = mₘᵢₙ + size(β)[1] - 1

    x = LinRange(-10, 10, nx)
    yy = (x * ones(nx)')
    xx = (x * ones(nx)')'

    φφ = atan.(xx, yy) # φ is defined "East of North" (Palumbo 2020)
    ρ_mask = ρₘᵢₙ .< .√(xx.^2 + yy.^2) .< ρₘₐₓ

    P = zeros((nx,nx)) + im*zeros((nx,nx))
    for m ∈ mₘᵢₙ:mₘₐₓ
        P[ρ_mask] += β[m] .* exp.(im * m * φφ[ρ_mask])
    end

    I = norm.(P)
    Q = real.(P)
    U = imag.(P)

    return I, Q, U
end


function bar_plot(β, file_path)
    mₘᵢₙ = β.offsets[1]+1
    mₘₐₓ = mₘᵢₙ + size(β)[1] - 1
    bar(mₘᵢₙ:mₘₐₓ, collect(abs.(β)))
    savefig("$(file_path)_abs.png")
    bar(mₘᵢₙ:mₘₐₓ, collect(angle.(β)))
    savefig("$(file_path)_arg.png")
end

function plot_beta_ring(β, file_path, n_ticks=17)
    mₘᵢₙ = β.offsets[1]+1
    mₘₐₓ = mₘᵢₙ + size(β)[1] - 1

    r = 3√3

    φs = LinRange(0, 2π, n_ticks)
    P = zeros(n_ticks) + im * zeros(n_ticks)
    # ϕs = atan.(r .* cos.(φs), r .* sin.(φs))
    # ϕs = (ϕs .+ 2π) .% 2π
    # println(ϕs)
    for m in mₘᵢₙ:mₘₐₓ
        P .+= β[m] .* exp.(im .* m .* (π/2 .- φs)) 
    end
    # P .+= β[2] .* exp.(im .* 2 .* φs) 

    plot(size=(600,600), xlims=(-1.3*r, 1.3*r), ylims=(-1.3*r, 1.3*r))


    s = 1
    for i ∈ 1:n_ticks
        φ = φs[i]
        # χ = χs[i]


        # Proper Calculation
        x₀ = r * sin(φ)
        y₀ = r * cos(φ)
        

        χ = 1/2 * atan(imag(P[i]), real(P[i]))

        # if (x₀ > 1)*(y₀>4) println(φs[i]/π, P[i], χ/π) end
        plot!([x₀ - sin(χ)/2*s, x₀ + sin(χ)/2*s],
        [y₀ + cos(χ)/2*s, y₀ - cos(χ)/2*s],
        color="black",
        label=nothing,
        lw=3);

  # plot!([x₀ + sin(χ)/2*s, x₀ - sin(χ)/2*s],
        # [y₀ - cos(χ)/2*s, y₀ + cos(χ)/2*s],
        # color="black",
        # label=nothing,
        # lw=3);

        # Kludged Calculation
        # x₀ = r * sin(φ)
        # y₀ = r * cos(φ)
        # χ = 1/2 * atan(imag(P[i]), real(P[i]))
        # plot!([x₀ - cos(χ)/2*s, x₀ + cos(χ)/2*s],
        # [y₀ - sin(χ)/2*s, y₀ + sin(χ)/2*s],
        # color="black",
        # label=nothing,
        # lw=3);
    end
    savefig(file_path)

end

function QU_loop(β, n, path)
    mₘᵢₙ = β.offsets[1]+1
    mₘₐₓ = mₘᵢₙ + size(β)[1] - 1

    ϕ = LinRange(0, 2π, n)
    P = zeros(n) .+ 0*im

    for m ∈ mₘᵢₙ:mₘₐₓ
        P .+= β[m] .* exp.(im * m * ϕ)
    end
    Q = real.(P)
    U = imag.(P)
    plot(Q, U, label=nothing, lw=3)


    savefig(path)
end


function main()

    dir = "monopole_a99"
    pth = "/scratch/gpfs/wcukier/Polarized_Ray_Tracing/0_stored_data/$(dir)/image"
    # suffix = "_phase"
    suffix = ""
    side="cs_"
    n="1"

    # β = offset_beta_array(-4, 4)
    # β[2] = 1
    # β[0] = 0
    # β[-1] = 0
    # n_phi = 1
    # n_theta = 1
    # nx = 1000
    # println(β)
    # aI, aQ, aU = generate_fakedata(β, nx)

    aI = readdlm("$(pth)/im_E_$(side)$(n)$(suffix).dat")
    aQ = readdlm("$(pth)/im_Q_$(side)$(n)$(suffix).dat")
    aU = readdlm("$(pth)/im_U_$(side)$(n)$(suffix).dat")

    for n ∈ 2:3
        aI = aI .+ readdlm("$(pth)/im_E_$(side)$(n)$(suffix).dat")
        aQ = aQ .+ readdlm("$(pth)/im_Q_$(side)$(n)$(suffix).dat")
        aU = aU .+ readdlm("$(pth)/im_U_$(side)$(n)$(suffix).dat")
    end

    nx = size(aI, 2)
    n_phi = 1
    n_theta = 40


    n_alpha =  n_phi * n_theta
    for path ∈ ["$(pth)/ring/", "$(pth)/QU_loop/", "$(pth)/bar/"]
        try
            mkdir(path)
        catch
            1==1
        end
    end

    for j ∈ 0:(n_theta-1)
        for k in 0:n_phi-1
            i = Int(k*n_theta + floor(j/n_phi))
            bI = aI[Int(1 + i*nx):Int((i+1)*nx), :]
            bQ = aQ[Int(1 + i*nx):Int((i+1)*nx), :]
            bU = aU[Int(1 + i*nx):Int((i+1)*nx), :]
            # β = (azimuthalfourier(bI, bQ, bU, nx, 1000))
            β = (azimuthalfourier(bI', bQ', bU', nx, 1000))

            # println(norm.(β))
            plot_beta_ring(β, "$(pth)/ring/$(suffix)_$(j)-$(k).png")
            QU_loop(β, 50, "$(pth)/QU_loop/$(suffix)_$(j)-$(k).png")
            bar_plot(β, "$(pth)/bar/$(suffix)_$(j)-$(k).png")

        end
    end
    
end


# main()