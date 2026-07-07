#include <stdio.h>
#include <math.h>

main()
{
 int i,j,nn=100;
 double alpha=0.08333333,test1,test2,arg;

 for(i=1;i<=nn;i++)
 {
  for(j=i;j<=nn;j++)
  {
   test1=erfc(alpha*(j-i));
   arg=(j-i)*alpha;
   arg*=arg;
   test2=exp(-arg/2);
   printf("delta,test1,test2= %d %e %e\n",(j-i),test1,test2);
  }
 }

}
