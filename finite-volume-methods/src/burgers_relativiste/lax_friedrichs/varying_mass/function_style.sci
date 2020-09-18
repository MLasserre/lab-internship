function function_style()
	sda()      // remet a plat les proprietes des axes par defaut
	a = gda(); // a est le handle sur ces axes par defaut
	a.title.font_size = 1;
	a.title.font_style = %helvetica_bold;
	a.x_label.font_size = %14pts;
	a.x_label.font_style = %helvetica_italic;
	a.y_label.font_size = %14pts;
	a.y_label.font_style = %helvetica_italic;
	a.z_label.font_size = %14pts;
	a.z_label.font_style = %helvetica_italic;
	a.font_size = %14pts;
endfunction

