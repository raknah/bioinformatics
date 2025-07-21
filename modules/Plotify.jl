function plotifylightax(fig, row, col, title, xlabel, ylabel, grid = true)
	
    ax = Axis(
		fig[row, col],

        # TITLE
		title = title,
		titlesize = 49,
		titlegap = 21,

        # X-LABEL
		xlabel = xlabel,
		xlabelsize = 30,
		xticklabelsize = 21,
		xlabelpadding = 21,

        # Y-LABEL
		ylabel = ylabel,
		ylabelsize = 30,
		yticklabelsize = 21,
		ylabelpadding = 21,

        # SPINES
		spinewidth = 1.5,
		leftspinevisible = true,
		bottomspinevisible = true,
		rightspinevisible = false,
		topspinevisible = false,
	)

    # GRIDS
	if grid
		ax.xgridvisible[] = true
		ax.ygridvisible[] = true
		ax.xminorgridvisible[] = true
		ax.yminorgridvisible[] = true
		ax.xgridstyle[] = :dash
		ax.ygridstyle[] = :dash
		ax.xgridwidth[] = 2.1
		ax.ygridwidth[] = 2.1
		ax.xgridcolor[] = RGBAf(0.7, 0.7, 0.7, 0.7)
		ax.ygridcolor[] = RGBAf(0.7, 0.7, 0.7, 0.7)
	else
		ax.xgridvisible[] = false
		ax.ygridvisible[] = false
		ax.xminorgridvisible[] = false
		ax.yminorgridvisible[] = false
	end

	return ax
end
