using Plots
include("Visualize_v2.jl")
include("azimuthal_fourier.jl")

# pth = "/scratch/gpfs/wcukier/Polarized_Ray_Tracing/0_stored_data/$(dir)/image"
# suffix = ""

dirs=["monopole_a0", "monopole_a99", "toroidal_a0", "toroidal_a99", "wald_a0", "sim", "wald_a99"]
# suffixes = ["cs_1", "cs_1", "cs_2", "cs_2", "cs_2", "cs_2", "cs_2"]
# suffixes = ["d", "d", "d", "d", "d", "d", "d"]
suffixes = ["id", "id", "id", "id", "id", "id", "id"]

# dir = "monopole_a0"
# side = "cs"
# n="1"

ndirs = size(dirs)[1]
angle_index = [8, 32]
angle_label = ["face", "edge"]

for i ∈ 1:ndirs
    dir = dirs[i]
    pth = "/scratch/gpfs/wcukier/Polarized_Ray_Tracing/0_stored_data/$(dir)/image"
    suffix = suffixes[i]

    aI = readdlm("$(pth)/im_E_$(suffix).dat")
    aQ = readdlm("$(pth)/im_Q_$(suffix).dat")
    aU = readdlm("$(pth)/im_U_$(suffix).dat")


    nx = size(aI, 2)
    n_phi = 1
    n_theta = 40

    for j in 1:2
        index = angle_index[j]
        label = angle_label[j]
        stokes_I_plot(aI, aQ, aU, index, :default)
        savefig("../../../../figures/intenisty_$(dir)_$(label).png")
        mag_P_plot(aI, aQ, aU, index, :default)
        savefig("../../../../figures/polarized_intensity_$(dir)_$(label).png")
        EVPA_plot(aI, aQ, aU, index)
        savefig("../../../../figures/EVPA_$(dir)_$(label).png")
        P_over_I_plot(aI, aQ, aU, index)
        savefig("../../../../figures/degree_$(dir)_$(label).png")


        bI, bQ, bU, phi, theta, i = mask_data(aI, aQ, aU, index)


        β = (azimuthalfourier(bI', bQ', bU', nx, 1000))

        QU_loop(β, 50, "../../../../figures/QUloop_$(dir)_$(label).png")
        bar_plot(β, "../../../../figures/betabar_$(dir)_$(label).png")
    
    end
end
