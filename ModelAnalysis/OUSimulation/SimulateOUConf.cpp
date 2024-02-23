#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericMatrix SimulateOUConfidence(int n, NumericVector params, double delta=0.01,
                                   double maxT=9, bool stop_on_error=true) {
  NumericMatrix out(n, 3);
  
  double a   =    params[0];
  double zr  =    params[1];
  double szr =    params[2];
  double v   =    params[3];
  double sv  =    params[4];
  double k   =    params[5];
  double tau =    params[6];
  double t0  =    params[7];
  double st0 =    params[8];
  double s   =    params[9];
  double lambda = params[10];
  
  bool valid = true;
  
  if (a <= 0)                         { valid = false; Rcpp::Rcout << "error: invalid parameter a = " << a << std::endl;  }
  if (szr < 0 || szr > 1)             { valid = false; Rcpp::Rcout << "error: invalid parameter szr = " << szr << std::endl; }
  if (st0 < 0)                        { valid = false; Rcpp::Rcout << "error: invalid parameter st0 = " << st0 << std::endl; }
  if (sv < 0)                         { valid = false; Rcpp::Rcout << "error: invalid parameter sv = " << sv << std::endl; }
  if (t0 - 0.5*st0 < 0) { valid = false; Rcpp::Rcout << "error: invalid parameter combination t0 = " << t0 << ", st0 =" << st0 << std::endl; }
  if (zr - 0.5*szr <= 0)              { valid = false; Rcpp::Rcout << "error: invalid parameter combination zr = " << zr << ", szr = " << szr << std::endl;}
  if (zr + 0.5*szr >= 1)              { valid = false; Rcpp::Rcout << "error: invalid parameter combination zr = " << zr << ", szr = " << szr << std::endl;}
  if (tau < 0)                        { valid = false; Rcpp::Rcout << "error: invalid parameter tau = " << tau << std::endl;}
  if (k < 0 || k > 1)             { valid = false; Rcpp::Rcout << "error: invalid parameter k = " << szr << std::endl; }
  
  if (!valid) {
    if (stop_on_error) { Rcpp::stop("Error validating parameters.\n"); }
    else {return out;}
  }
  
  
  double mu, x0, t, conf;
  int resp;
  
  for (int i=0; i < n; i++) {
    mu = R::rnorm(v, sv);
    x0 = a* (R::runif(zr-szr/2, zr+szr/2)-0.5);
    t = 0;
    while ((x0 > -a/2) && (x0 < a/2) && (t < maxT)) {
      x0 = x0 - delta*k*x0 + R::rnorm(delta*mu, sqrt(delta)*s);
      t = t + delta;
    }
    if (x0 >= a/2) {
      resp = 1;
    } else {
      if (x0 <= -a/2) {
        resp = -1;
      } else {
        resp = 0;
      }
    }
    //if (w == 1) {
    if (tau > 0) {
      conf = resp*(exp(-k*tau)*x0 + mu/k*(1-exp(-k*tau))+ 
        R::rnorm(0, sqrt((1-exp(-2*k*tau))/(2*k))));
    } else {
      conf = resp*x0;
    }
    // save response, response time and state of evidence accumulator
    out( i , 0 ) = std::max(0.0, t )  + R::runif(t0-st0/2, t0+st0/2);;
    out( i , 1 ) = resp;
    out( i , 2 ) = conf/ pow(t+tau, lambda); // evidence term

    if (i % 200 ==0 ) Rcpp::checkUserInterrupt();
    
  }
  
  return out;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
// 
// /*** R
// timesTwo(42)
// */






// [[Rcpp::export]]
NumericMatrix SimulateMeanOUConfidence(int n, NumericVector params, double delta=0.01,
                                   double maxT=9, bool stop_on_error=true) {
  NumericMatrix out(2,3);
  
  double a   =    params[0];
  double zr  =    params[1];
  double szr =    params[2];
  double v   =    params[3];
  double sv  =    params[4];
  double k   =    params[5];
  double tau =    params[6];
  double t0  =    params[7];
  double st0 =    params[8];
  double s   =    params[9];
  double lambda = params[10];
  
  
  bool valid = true;
  
  if (a <= 0)                         { valid = false; Rcpp::Rcout << "error: invalid parameter a = " << a << std::endl;  }
  if (szr < 0 || szr > 1)             { valid = false; Rcpp::Rcout << "error: invalid parameter szr = " << szr << std::endl; }
  if (st0 < 0)                        { valid = false; Rcpp::Rcout << "error: invalid parameter st0 = " << st0 << std::endl; }
  if (sv < 0)                         { valid = false; Rcpp::Rcout << "error: invalid parameter sv = " << sv << std::endl; }
  if (t0 - 0.5*st0 < 0) { valid = false; Rcpp::Rcout << "error: invalid parameter combination t0 = " << t0 << ", st0 =" << st0 << std::endl; }
  if (zr - 0.5*szr <= 0)              { valid = false; Rcpp::Rcout << "error: invalid parameter combination zr = " << zr << ", szr = " << szr << std::endl;}
  if (zr + 0.5*szr >= 1)              { valid = false; Rcpp::Rcout << "error: invalid parameter combination zr = " << zr << ", szr = " << szr << std::endl;}
  if (tau < 0)                        { valid = false; Rcpp::Rcout << "error: invalid parameter tau = " << tau << std::endl;}
  if (k < 0 || k > 1)             { valid = false; Rcpp::Rcout << "error: invalid parameter k = " << szr << std::endl; }
  
  if (!valid) {
    if (stop_on_error) { Rcpp::stop("Error validating parameters.\n"); }
    else {return out;}
  }
  
  
  double mu, x0, t, conf;
  int resp;
  double conf1 = 0;
  double conf2 = 0;
  double count1 = 0;
  double count2 = 0;
  
  for (int i=0; i < n; i++) {
    mu = R::rnorm(v, sv);
    x0 = a* (R::runif(zr-szr/2, zr+szr/2)-0.5);
    t = 0;
    while ((x0 > -a/2) && (x0 < a/2) && (t < maxT)) {
      x0 = x0 - delta*k*x0 + R::rnorm(delta*mu, sqrt(delta)*s);
      t = t + delta;
    }
    if (x0 >= a/2) {
      resp = 1;
      count1 = count1+1;
      if (tau > 0) {
        conf = resp*(exp(-k*tau)*x0 + mu/k*(1-exp(-k*tau))+ 
          R::rnorm(0, sqrt((1-exp(-2*k*tau))/(2*k))));
      } else {
        conf = resp*x0;
      }
      conf1 = conf1 + conf / pow(t+tau, lambda);
    } else {
      if (x0 <= -a/2) {
        resp = -1;
        count2 = count2+1;
        if (tau > 0) {
          conf = resp*(exp(-k*tau)*x0 + mu/k*(1-exp(-k*tau))+ 
            R::rnorm(0, sqrt((1-exp(-2*k*tau))/(2*k))));
        } else {
          conf = resp*x0;
        }
        conf2 = conf2 + conf/ pow(t+tau, lambda);
      }
    }

    if (i % 200 ==0 ) Rcpp::checkUserInterrupt();
    
  }
  conf1 = conf1/count1;
  conf2 = conf2/count2;
  
  out(0,0) = 1;
  out(1,0) = -1;
  out(0,1) = conf1;
  out(1,1) = conf2;
  out(0,2) = count1;
  out(1,2) = count2;
  return out;
}
