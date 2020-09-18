function [] = plot_function(filename, graph_title)
    unit1 = file('open', filename, 'old')
    infos = read(unit1, 1, 1)
    n = infos(1) // nb de points en x
    x_left = infos(3)
    x_right = infos(4)
    t_f = infos(5)
    x_step = infos(6)
    t_step = infos(7)
    approx= read(unit1, -1, n)';
    m = size(approx, 2) // en cas d'arret premature
    // taux d'echantillonnage en temps
    pech = m/10
    x=linspace(x_left,x_right,n)';
    file('close', unit1)
    xset("window", 1)
    //xbasc(1)
    xtitle(graph_title, "position (x)", "u(x,t)")
    plot2d(x, approx(:, 1))
endfunction

