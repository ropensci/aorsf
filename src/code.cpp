#include <Rcpp.h>
using namespace Rcpp;


bool is_constant(NumericVector x){
 for(int i=1; i<x.length(); i++){
  if(x[i] != x[0]) return(false);
 }
 return(true);
}

double lrtestC(NumericVector time,
               NumericVector status,
               NumericVector grp){

 int n = time.length();
 double Y=n;
 double Y1=sum(grp);

 int lwr=0, upr=0;

 for(int i=0; i<n; i++){
  if(status[i]==0){
   upr++;
  } else {
   break;
  }
 }

 IntegerVector indx = seq(lwr,upr);
 double d = sum(as<NumericVector>(status[indx]));
 double d1= sum(as<NumericVector>(status[indx]) * as<NumericVector>(grp[indx]));

 double e1=Y1*d/Y;
 double e0=(Y-Y1)*d/Y;
 double o1=d1;
 double o0=d-d1;

 double V = (Y-Y1)*Y1*d*(Y-d) / (pow(Y,2)*(Y-1));
 Y -= indx.length();
 Y1 -= sum(as<NumericVector>(grp[indx]));

 lwr=upr+1;
 int counter=1;

 for( ; ; ){

  upr=lwr;
  while((status[upr]==0) & (upr<n-1)){
   upr++;
  }

  if(upr==n-1){
   if(status[upr]==0){
    break;
   } else {

    indx=seq(lwr,upr);

    d=sum(as<NumericVector>(status[indx]));
    d1=sum(as<NumericVector>(status[indx]) * as<NumericVector>(grp[indx])  );

    e1+=(Y1*d/Y);
    e0+=((Y-Y1)*d/Y);
    o1+=d1;
    o0+=d-d1;

    V += (Y-Y1)*Y1*d*(Y-d) / (pow(Y,2)*(Y-1));
    Y -= indx.length();
    Y1 -= sum(as<NumericVector>(grp[indx]));

    break;
   }
  }

  indx=seq(lwr,upr);

  d=sum(as<NumericVector>(status[indx]));
  d1=sum(as<NumericVector>(status[indx]) * as<NumericVector>(grp[indx])  );

  e1+=(Y1*d/Y);
  e0+=((Y-Y1)*d/Y);
  o1+=d1;
  o0+=d-d1;

  V += (Y-Y1)*Y1*d*(Y-d) / (pow(Y,2)*(Y-1));
  Y -= indx.length();
  Y1 -= sum(as<NumericVector>(grp[indx]));
  counter++;
  lwr=upr+1;

  if(Y==1) break;

 }

 return pow(o1-e1,2) / V;

}
