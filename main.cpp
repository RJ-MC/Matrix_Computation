#include <iostream>
using namespace std;

#include "includes/matrix.h"

int main() {
	
	cout<<(5*eye(5)/(8.0*ones(matrix<int>(2,1,{5,5}))));

	return 0;
}

