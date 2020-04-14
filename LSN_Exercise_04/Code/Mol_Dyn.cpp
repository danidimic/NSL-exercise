#include <vector>
#include <cmath>
using namespace std;

double position_Verlet(double r0, double rminus, double a, double dt){
	return 2*r0-rminus+pow(dt,2)*a;
}

double velocity_Verlet(double rminus, double rplus, double dt){
	return (rplus-rminus)/(2*dt);
}
