mu = 0.823
sigma = 0.6455


f(x) = 1./(2*sigma*sqrt(pi)) * exp( -((x-mu)**2)/(2*sigma*sigma) ) + 1./(2*sigma*sqrt(pi)) * exp( -((x+mu)**2)/(2*sigma*sigma) )
fteo(x) = f(x)*f(x)

plot "Psi.out" w l
replot fteo(x) w l ls 2
