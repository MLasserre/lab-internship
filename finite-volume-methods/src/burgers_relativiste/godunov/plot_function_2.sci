function [] = plot_function_2(filename, graphic_name)

	// Recuperation of data in the file
    unit1 = file('open', filename, 'old')
    infos = read(unit1, 1, 7)

    n = infos(1)                   // nb of points for x
    m = infos(2)                   // nb of points for time
    x_left = infos(3)
    x_right = infos(4)
    t_f = infos(5)
    x_step = infos(6)
    t_step = infos(7)

    approx= read(unit1, -1, n)';
    x=linspace(x_left,x_right,n)'; // Domain for x
    m = size(approx, 2)            // In case it stops prematurely
    pech = m/4                     // Rate of sampling in time
    file('close', unit1)


    // Setting and plot of curves
    xset("window", 1)
    //xtitle(graph_title, "position (x)", "u(x,t)") // Titles
    xtitle('${\Huge\mathrm{Shock\ Wave}}$')
    plot2d(x, approx(:, 1:pech:m)) // Plot

	a=get("current_axes") // Get the handle of the newly created axes
	a.axes_visible="on"; // Makes the axes visible

	a.children.children(1).polyline_style=1;
	a.children.children(1).foreground=27;
	a.children.children(1).thickness=1;
	a.children.children(1).line_style=8;

	a.children.children(2).polyline_style=1;
	a.children.children(2).foreground=19;
	a.children.children(2).thickness=1;
	a.children.children(2).line_style=5;

	a.children.children(3).polyline_style=1;
	a.children.children(3).foreground=13;
	a.children.children(3).thickness=1;
	a.children.children(3).line_style=2;

	a.children.children(4).polyline_style=1;
	a.children.children(4).foreground=1;
	a.children.children(4).thickness=1;
	a.children.children(1).line_style=8;

	a.font_size=3; // Set the tics label font size
	a.x_location="bottom"; // Set the x axis position
	a.data_bounds=[x_left,-0.6;
                   x_right,0.6]; //Set the boundary values for x
	a.sub_tics=[0,0];
	a.labels_font_color=1;
	a.labels_font_size=3;
	//a.grid=[2,2];
	a.box="on";
	//a.thickness=1;
	//a.background=2;
	//a.foreground=2;

	// Graphic export
	eps_graphic_name = graphic_name + ".eps"
	//png_graphic_name = graphic_name + ".png"

	// EPS export
	xs2eps(gcf(), eps_graphic_name)
	// PNG export
	//xs2png(gcf(), png_graphic_name)

endfunction

