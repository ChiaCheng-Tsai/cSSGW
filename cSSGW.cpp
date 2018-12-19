
#include"SSGW.m.hpp"
#include<armadillo>
#include<vector>

int main()
{
 cx_vec ws;
 cx_vec zs;
 vec PP;
 SSGW(datum::inf,0.3,zs,ws,PP);
// SSGW(2.,0.3,zs,ws,PP);

 cout<<std::setprecision(15);
 //(zs).raw_print(std::cout);
 cout<<endl;
 //(ws).raw_print(std::cout);
 cout<<endl;
 (PP).raw_print(std::cout);
 system("pause");
 return 0;
}
