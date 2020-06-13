#include <iostream>
using namespace std;

#include "includes/matrix.h"

int main() {
	
	cout<<(matrix<int>(2,1,{2,2})^matrix<int>(2,1,{2,3}));

	return 0;
}

