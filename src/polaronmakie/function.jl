# function.jl

"""
Uses the AbstractPlotting.jl and MakieLayout.jl plotting packages to produce 3D interactive graphs of complex functions. The Real and Imaginary parts of the function are displayed side-by-side in two different 3D plots and can be altered by any finite number of sliders to change the value of any arguments they may depend on. Axis sliders are also present to either widen or focus the area of the graphs under investigation.
"""

mutable struct FunctionView
    params::Any
    sliders::Any
    annotations::Any
    limits::Any
    ℜ_matrix::Any
    ℑ_matrix::Any
    scenes::Any
    layout::Any
end
FunctionView(xs::AbstractArray) = FunctionView(xs..., [], [])

function functionview(
    ℜ,
    ℑ,
    params::Tuple,
    sliders::Tuple,
    annotations::Dict{String,String};
    axis_limits = ((-20, 20), (-20, 20), (-20, 20)),
    len = 100,
)

    x_lim = lift(x -> range(axis_limits[1][1], stop = x, length = len), sliders[1].value)
    y_lim = lift(x -> range(axis_limits[2][1], stop = x, length = len), sliders[2].value)
    limits = (x_lim, y_lim)

    ℜ_matrix = lift(map(s -> s.value, sliders)...) do x_max, y_max, z_max, params...
        ℜf(x, y) = ℜ(x, y, params...)
        [
            ifelse(abs(ℜf(x, y)) < z_max, ℜf(x, y), NaN)
            for
            x in range(axis_limits[1][1], stop = x_max, length = len),
            y in range(axis_limits[2][1], stop = y_max, length = len)
        ]
    end

    ℑ_matrix = lift(map(s -> s.value, sliders)...) do x_max, y_max, z_max, params...
        ℑf(x, y) = ℑ(x, y, params...)
        [
            ifelse(abs(ℑf(x, y)) < z_max, ℑf(x, y), NaN)
            for
            x in range(axis_limits[1][1], stop = x_max, length = len),
            y in range(axis_limits[2][1], stop = y_max, length = len)
        ]
    end

    return FunctionView(map(
        x -> x,
        [params, sliders, annotations, limits, ℜ_matrix, ℑ_matrix],
    ))
end

function viewfunction(
    ℜ,
    ℑ,
    params::Tuple,
    param_initial_values::Tuple,
    param_limits,
    annotations::Dict{String,String};
    axis_limits = ((-20, 20), (-20, 20), (-20, 20)),
    len = 100,
    resolution = (1500, 800),
)
    # Set the scene
    outer_padding = 30
    scene, layout = layoutscene(
        outer_padding,
        resolution = resolution,
        backgroundcolor = RGBf0(0.99, 0.99, 0.99),
    )
    makescene() = LScene(
        scene,
        camera = cam3d!,
        raw = true,
        scenekw = (backgroundcolor = RGBf0(0.96, 0.96, 0.96), clear = false),
    )

    real_scene = layout[2, 2] = makescene()
    imag_scene = layout[2, 3] = makescene()

    colsize!(layout, 2, Relative(0.45))
    colsize!(layout, 3, Relative(0.45))

    real_scene_title =
        layout[1, 2] = LText(
            scene,
            annotations["ℜ_title"],
            textsize = 20,
            padding = (0, 0, 10, 0),
            font = "Noto Sans Bold",
            tellwidth = false,
        )
    imag_scene_title =
        layout[1, 3] = LText(
            scene,
            annotations["ℑ_title"],
            textsize = 20,
            padding = (0, 0, 10, 0),
            font = "Noto Sans Bold",
            tellwidth = false,
        )

    # Axis sliders

    sliders = []
    coordinates = ("x", "y", "z")

    for c = 1:length(coordinates)
        coord = coordinates[c]
        slider_label = Symbol(coord, '_', "slider", '_', "label")
        slider = Symbol(coord, '_', "slider")
        slider_value = Symbol(coord, '_', "slider", '_', "value")

        @eval $slider_label = $(layout[3+c, 1] = LText(scene, coord, textsize = 20))

        @eval $slider =
            $(
                layout[3+c, 2:3] = LSlider(
                    scene,
                    range = axis_limits[c][1]:0.01:axis_limits[c][2],
                    horizontal = true,
                    tellheight = true,
                    width = nothing,
                    height = Auto(),
                    startvalue = axis_limits[c][2],
                )
            )
        @eval value = $slider.value
        @eval $slider_value =
            $(layout[3+c, 4] = LText(scene, lift(x -> "$(x)", value), textsize = 20))

        push!(sliders, eval(slider))
    end

    # Parameter sliders

    para_num = length(params)

    for p = 1:para_num
        param = params[p]
        slider_label = Symbol(param, '_', "slider", '_', "label")
        slider = Symbol(param, '_', "slider")
        slider_value = Symbol(param, '_', "slider", '_', "value")

        @eval $slider_label = $(layout[1, 4+p] = LText(scene, "$(param)", textsize = 20))
        @eval $slider =
            $(
                layout[2, 4+p] = LSlider(
                    scene,
                    range = param_limits[p][1]:0.01:param_limits[p][2],
                    horizontal = false,
                    tellwidth = true,
                    height = nothing,
                    width = Auto(),
                    startvalue = param_initial_values[p],
                )
            )
        @eval value = $slider.value
        @eval $slider_value =
            $(layout[3, 4+p] = LText(scene, lift(x -> "$(x)", value), textsize = 20))

        push!(sliders, eval(slider))
        colsize!(layout, 4 + p, Relative(0.07 / para_num))
    end

    sliders = tuple(sliders...)

    # Title
    title =
        layout[0, :] = LText(
            scene,
            "PolaronMakie.jl",
            textsize = 30,
            padding = (0, 0, 50, 0),
            font = "Noto Sans Bold",
            color = (:black, 0.25),
        )

    # Function data
    fv = functionview(ℜ, ℑ, params, sliders, annotations, axis_limits = axis_limits)

    # Plotting
    sℜ = AbstractPlotting.surface!(real_scene, fv.limits[1], fv.limits[2], fv.ℜ_matrix)
    cℜ = AbstractPlotting.contour!(
        real_scene,
        fv.limits[1],
        fv.limits[2],
        fv.ℜ_matrix,
        levels = 30,
        linewidth = 1,
        transformation = (:xy, -5.0 + minimum(filter(!isnan, [to_value(fv.ℜ_matrix)...]))),
    )
    AbstractPlotting.wireframe!(
        real_scene,
        fv.limits[1],
        fv.limits[2],
        fv.ℜ_matrix,
        overdraw = true,
        transparency = true,
        color = (:black, 0.15),
    )

    sℑ = AbstractPlotting.surface!(
        imag_scene,
        fv.limits[1],
        fv.limits[2],
        fv.ℑ_matrix,
        scale_plot = false,
    )
    cℑ = AbstractPlotting.contour!(
        imag_scene,
        fv.limits[1],
        fv.limits[2],
        fv.ℑ_matrix,
        levels = 30,
        linewidth = 1,
        transformation = (:xy, -5.0 + minimum(filter(!isnan, [to_value(fv.ℑ_matrix)...]))),
    )

    AbstractPlotting.wireframe!(
        imag_scene,
        fv.limits[1],
        fv.limits[2],
        fv.ℑ_matrix,
        overdraw = true,
        transparency = true,
        color = (:black, 0.15),
    )

    lift(fv.ℜ_matrix, fv.ℑ_matrix, z_slider.value) do ℜz, ℑz, z
        AbstractPlotting.translate!(cℜ, 0, 0, -5.0 + minimum(filter(!isnan, [to_value(ℜz)...])))
        AbstractPlotting.translate!(cℑ, 0, 0, -5.0 + minimum(filter(!isnan, [to_value(ℑz)...])))
    end

    AbstractPlotting.display(scene)

    fv.scenes = [scene, real_scene, imag_scene]
    fv.layout = layout

    return fv
end
