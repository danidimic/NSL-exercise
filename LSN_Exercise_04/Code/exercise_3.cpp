#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"
#include "statistics.h"
using namespace std;


int main (int argc, char *argv[]){

	Random rnd = RandomGenerator();		//imposto il generatore di numeri casuali

	cout<<"eserc 3"<<endl;

	rnd.SaveSeed();
	return 0;
}
