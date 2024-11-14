using Plots; theme(:ggplot2, linewidth=2)
using LaTeXStrings
using Revise
using XLSX

includet("LiftingLineBaselineEdited.jl")

# Define the data array to hold aerofoil data
aerofoil_data = Matrix{Any}(undef, 4, 4)
# Define the aerofoils and related data
aerofoils = ["SD7037", "AG24", "NACA 2410", "RG15"]
m_values = [6.143602236, 6.269267229, 6.061452736, 6.066685483]
Cd0_values = [0.0114917, 0.0057377, 0.0090146, 0.0061735]
alpha0_values = [-3.414194915, -2.609649123, -2.104885057, -2.61732852]
aerofoil_data[:, 1] = aerofoils
aerofoil_data[:, 2] = m_values
aerofoil_data[:, 3] = Cd0_values
aerofoil_data[:, 4] = alpha0_values
println(aerofoil_data)

# Define ranges for wing properties
span_values = 1.5:0.5:4.5
chord_root_values = [0.220872446, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]
alpha_root_values = 2:0.5:6
alpha_tip_values = 2:0.5:6

# Propulsion elements
engine_no = 1
mass_engine = 2
mass_aux=mass_engine*0.5
engine_power = 3.7e3
prop_efficiency = 0.8
SFC = 1e-5


# Define ranges for mass and fuel mass
#mass = 150
mass_payload = 5
mass_equipment = 1.8 * mass_payload
#mass_payload = (mass/5.312168)^(1/0.972)
mass = (5.312168*mass_equipment^0.972)
mass_empty =  (0.8105*mass^0.9425)
mass_propulsion = mass_engine*engine_no+mass_aux
mass_airframe = mass_empty - mass_propulsion

mass_fuel = mass-mass_empty
# mass_fuel = 30

mass_initial=mass
mass_final=mass-mass_fuel

# Define constants for calculations
g = 9.81
rho = 1.225

# Track the best configuration for endurance factor
global_best_endurance_factor = 0
global_best_endurance_configuration = Dict{String, Any}()

extracted_data = Dict(
    "chord_root" => Float64[], 
    "alpha_root" => Float64[], 
    "endurance" => Float64[],
    "AR" => Float64[],
)

# Nested loop over each combination of wing properties, mass, and fuel mass
for span in span_values
    for chord_root in chord_root_values
        for alpha_root in alpha_root_values
            for alpha_tip in alpha_tip_values
                        for i in 1:length(m_values)
                            m_const = m_values[i]
                            Cd0 = Cd0_values[i]
                            alpha0 = alpha0_values[i]

                            y = uniform_grid(span, 35)
                            ecc = eccentricity(span, chord_root)

                            θ(y) = theta(y, span)
                            c(y) = 2 * chordhalf.(y, span, ecc)
                            a(y) = linear_function(y, alpha_root, alpha_tip, span)
                            a0(y) = alpha0
                            m(y) = m_const

                            C, D = build_linear_system(y, span, θ, c, m, a, a0)
                            A = C \ D

                            G0 = gamma0(y, A, θ)
                            area = wing_area(span, chord_root)
                            AR = aspect_ratio(span, area)
                            δ = delta(A)
                            E = efficiency_factor(δ)
                            Cl = lift_coefficient(A, AR)
                            Cd = Cd0 + drag_coefficient(Cl, AR, E)
                            MAC = c((span / 2) * 0.4244)

                            if Cl >= 0
                                V = sqrt((2 * mass * g) / (rho * area * Cl))

                                if V < 40
                                    endurance = (prop_efficiency/SFC)*sqrt(2 * rho * area) * (Cl^(3/2) / Cd) * (((mass_final*9.81)^(-1/2))-((mass_initial*9.81)^(-1/2)))/3600
                                    range = (prop_efficiency/SFC)*(Cl/Cd)*log((mass_initial*9.81)/(mass_final*9.81))/1000
                                    #endurance = (Cl^(3/2) / Cd)
                                    D = 0.5 * rho * V^2 * area * Cd
                                    P = D * V
                                    T = engine_no * engine_power * prop_efficiency / V
                                    
                                    #println("Checking endurance factor: ", endurance)
                                    if T >= D
                                    if span == 4.5 && aerofoils[i] == "AG24"
                                        # Store extracted data
                                         push!(extracted_data["chord_root"], chord_root)
                                         push!(extracted_data["alpha_root"], alpha_root)
                                         push!(extracted_data["endurance"], endurance)
                                         push!(extracted_data["AR"], AR)
                                    end
                                        if endurance > global_best_endurance_factor
                                            #println("Updating best configuration with endurance factor: ", endurance)
                                            global global_best_endurance_factor = endurance
                                            global global_best_endurance_configuration = Dict(
                                                "Aerofoil" => aerofoils[i],
                                                "m_value" => m_const,
                                                "Cd0_value" => Cd0,
                                                "Cl" => Cl,
                                                "Cd" => Cd,
                                                "Cl/Cd" => Cl / Cd,
                                                "Endurance Factor" => endurance,
                                                "Freestream Velocity (V)" => V,
                                                "AR" => AR,
                                                "span" => span,
                                                "chord_root" => chord_root,
                                                "alpha_root" => alpha_root,
                                                "alpha_tip" => alpha_tip,
                                                "alpha0" => alpha0,
                                                "G0" => G0,
                                                "MAC" => MAC,
                                                "area" => area,
                                                "efficiency_factor" => E,
                                                "delta" => δ,
                                                "mass" => mass,
                                                "mass_fuel" => mass_fuel,
                                                "payload_mass" => mass_payload,
                                                "empty_weight" => mass_empty,
                                                "Lift" => mass * g,  # Lift equals weight in steady level flight
                                                "Drag" => D,
                                                "Power" => P,
                                                "y" => y,
                                                "θ" => θ,
                                                "c" => c,
                                                "a" => a,
                                                "a0" => a0,
                                                "range" => range
                                            )
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

# Output the overall best configuration for endurance factor
println("\nBest endurance factor configuration across all configurations:")
println("Aerofoil: $(global_best_endurance_configuration["Aerofoil"])")
println("With m_value: $(round(global_best_endurance_configuration["m_value"], sigdigits=3)), Cd0_value: $(round(global_best_endurance_configuration["Cd0_value"], sigdigits=3)), Alpha0: $(round(global_best_endurance_configuration["alpha0"], sigdigits=3))\n")
println("Cl: $(round(global_best_endurance_configuration["Cl"], sigdigits=3)), Cd: $(round(global_best_endurance_configuration["Cd"], sigdigits=5))")
println("Cl/Cd: $(round(global_best_endurance_configuration["Cl/Cd"], sigdigits=3))")
println("Endurance [h]: $(round(global_best_endurance_configuration["Endurance Factor"], sigdigits=3)), Range [km]: $(round(global_best_endurance_configuration["range"], sigdigits=3))")

println("Freestream Velocity (V): $(round(global_best_endurance_configuration["Freestream Velocity (V)"], sigdigits=3)) m/s\n")
println("Lift: $(round(global_best_endurance_configuration["Lift"], sigdigits=3)) N")
println("Drag: $(round(global_best_endurance_configuration["Drag"], sigdigits=3)) N")
println("Power: $(round(global_best_endurance_configuration["Power"], sigdigits=3)) W\n")
println("Span: $(round(global_best_endurance_configuration["span"], sigdigits=3)), Chord Root: $(round(global_best_endurance_configuration["chord_root"], sigdigits=3))")
println("Aspect Ratio (AR): $(round(global_best_endurance_configuration["AR"], sigdigits=3))")
println("Alpha Root: $(global_best_endurance_configuration["alpha_root"]), Alpha Tip: $(global_best_endurance_configuration["alpha_tip"])\n")
println("Mass: $(round(global_best_endurance_configuration["mass"], sigdigits=3)), Fuel Mass: $(round(global_best_endurance_configuration["mass_fuel"], sigdigits=3))")
println("Payload Mass: $(round(global_best_endurance_configuration["payload_mass"], sigdigits=3))")
println("Empty Mass: $(round(global_best_endurance_configuration["empty_weight"], sigdigits=3))")
println("Airframe Mass: $(round(mass_airframe, sigdigits=3))")
# Post-processing results for the best endurance factor configuration
y_best = global_best_endurance_configuration["y"]
G0_best = global_best_endurance_configuration["G0"]
span_best = global_best_endurance_configuration["span"]
chord_root_best = global_best_endurance_configuration["chord_root"]

# # Plot Circulation and Chord Length
# p1 = plot(
#     y_best, G0_best, label="Circulation", legend=false,
#     xlabel="Spanwise coordinate [m]", ylabel=L"\Gamma")

# p2 = plot(
#     y_best, global_best_endurance_configuration["c"](y_best), label="Chord Length",
#     aspect_ratio=:equal, xlabel="Spanwise coordinate [m]", ylabel="Chord Length [m]",
#     title="Chord Length vs Spanwise Coordinate in Upper-Right Quadrant"
# )
# plot!(xlims=(0, span_best / 2), ylims=(0, chord_root_best + 1))

# Plot wing platform as an ellipse
t = range(0, 2π, length=100)
x_ellipse = (span_best / 2) * cos.(t)
y_ellipse = (chord_root_best / 2) * sin.(t)

p3 = plot(
    x_ellipse, y_ellipse,
    label="Wing Planform",
    xlabel="Spanwise coordinate [m]", ylabel="Chord [m]",
    aspect_ratio=:equal,
    title="Wing Planform"
)
plot!(xlims=(-span_best / 2 -0.5, span_best / 2 + 0.5), ylims=(-chord_root_best / 2 -0.5, chord_root_best / 2 + 0.5))

# # Combine plots
# plot(p1, p2, p3, layout=(1, 3), size=(1500, 500))

# chord_roots = extracted_data["chord_root"]
# alpha_roots = extracted_data["alpha_root"]
# endurances = extracted_data["endurance"]
# ARs = extracted_data["AR"]

# surface(
#     ARs, alpha_roots, endurances, 
#     xlabel = "Aspect Ratio [-]", 
#     ylabel = "Alpha Root [deg]", 
#     zlabel = "Endurance [h]", 
#     title = "Endurance Variation for Optimal Configuration"
# )

# surface(
#     chord_roots, alpha_roots, endurances, 
#     xlabel = "Root Chord [m]", 
#     ylabel = "Alpha Root [deg]", 
#     zlabel = "Endurance [h]", 
#     title = "Endurance Variation for Optimal Configuration"
# )
