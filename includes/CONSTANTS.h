#ifndef CONSTANTS_h
#define CONSTANTS_h

#define ROUND(x) (int((x)<0?(x-0.5):(x+0.5)))
#define SIGN(x)  (int((x)<0?-1:((x)>0?1:0)))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))
#define IS_ZERO(a) (abs(a)<EPS)
#define IS_NON_ZERO(a) (abs(a)>EPS)

constexpr auto EPS 	= 1e-6;
constexpr auto PI 	= 3.14159265358979323846;


#endif

