#ifndef CONSTANTS_h
#define CONSTANTS_h

#define ROUND(x) (int(x<0?(x-0.5):(x+0.5)))
#define MAX(a,b) (a>b?a:b)
#define MIN(a,b) (a>b?b:a)
#define SIGN(x)  (x>0?1:-1)

constexpr auto EPS = 1e-6;
constexpr auto PI = 3.14159265358979323846;



#endif

