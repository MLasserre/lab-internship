unit1 = file('open', "lax_friedrichs.dat", 'old')
dims = read(unit1, 1, 2)
n = dims(1) // nb de points en x
m = dims(2) // nb de points calcules en temps
approx= read(unit1, -1, n)';
m = size(approx, 2) // en cas d'arret premature
// taux d'echantillonnage en temps
pech = m/10
x=linspace(0,n-1,n)';
file('close', unit1)
xset("window", 1)
//xbasc(1)
xtitle("Burgers : schema de Lax-Friedrichs", "position (x)", "u(x,t)")
plot2d(x, approx(:, 1:pech:m))

