#include <Rcpp.h>
using namespace Rcpp;

// source this function into an R session using the Rcpp::sourceCpp 

// [[Rcpp::export]]
NumericMatrix RNG_WEV_SubjTime (int n, NumericVector params, double delta=0.01,
                                double maxT=9, bool stop_on_error=true)
{
  NumericMatrix out(n, 6);
  
  double a   =    params[0];
  double v   =    params[1];
  double t0  =    params[2];
  double d   =    params[3];
  double szr =    params[4];
  double sv  =    params[5];
  double st0 =    params[6];
  double zr  =    params[7];
  double tau =    params[8];
  double lambda =  params[9];
  double w   =    params[10];
  double muvis =  params[11];
  double sigvis = params[12];
  double svis =   params[13];
  double Timerate = params[14];
  double Timediff = params[15];
  
  bool valid = true;
  
  if (a <= 0)                         { valid = false; Rcpp::Rcout << "error: invalid parameter a = " << a << std::endl;  }
  if (szr < 0 || szr > 1)             { valid = false; Rcpp::Rcout << "error: invalid parameter szr = " << szr << std::endl; }
  if (st0 < 0)                        { valid = false; Rcpp::Rcout << "error: invalid parameter st0 = " << st0 << std::endl; }
  if (sv < 0)                         { valid = false; Rcpp::Rcout << "error: invalid parameter sv = " << sv << std::endl; }
  if (t0 - fabs(0.5*d) - 0.5*st0 < 0) { valid = false; Rcpp::Rcout << "error: invalid parameter combination t0 = " << t0 << ", d = " << d << ", st0 =" << st0 << std::endl; }
  if (zr - 0.5*szr <= 0)              { valid = false; Rcpp::Rcout << "error: invalid parameter combination zr = " << zr << ", szr = " << szr << std::endl;}
  if (zr + 0.5*szr >= 1)              { valid = false; Rcpp::Rcout << "error: invalid parameter combination zr = " << zr << ", szr = " << szr << std::endl;}
  if (tau < 0)                        { valid = false; Rcpp::Rcout << "error: invalid parameter tau = " << tau << std::endl;}
  if (w<0 || w > 1)                   { valid = false; Rcpp::Rcout << "error: invalid parameter w = " << w << ", allowed: w in [0,1]" <<  std::endl; }
  if (sigvis < 0)                      { valid = false; Rcpp::Rcout << "error: invalid parameter sigvis = " << sigvis <<  std::endl; }
  if (svis <= 0)                        { valid = false; Rcpp::Rcout << "error: invalid parameter svis = " << svis <<  std::endl; }
  if (lambda < 0)                        { valid = false; Rcpp::Rcout << "error: invalid parameter lambda = " << lambda <<  std::endl; }
  if (Timediff <= 0)                        { valid = false; Rcpp::Rcout << "error: invalid parameter Timediff = " << Timediff <<  std::endl; }
  
  if (!valid) {
    if (stop_on_error) { Rcpp::stop("Error validating parameters.\n"); }
    else {return out;}
  }
  
  
  double mu, x0, t, conf, vis, subj_time;
  int resp;
  
  for (int i=0; i < n; i++) {
    mu = R::rnorm(v, sv);
    x0 = a* R::runif(zr-szr/2, zr+szr/2);
    t = 0;
    while ((x0 > 0) && (x0 < a) && (t < maxT)) {
      x0 = x0 + R::rnorm(delta*mu, sqrt(delta));
      t = t + delta;
    }
    if (x0 >= a) {
      resp = 1;
    } else {
      if (x0 <= 0) {
        resp = -1;
      } else {
        resp = 0;
      }
    }
    //if (w == 1) {
    if (tau > 0) {
      conf = resp*(x0 + R::rnorm(tau*mu, sqrt(tau)) - a*zr);
    } else {
      conf = resp*(x0 - a*zr);
    }
    // save response, response time and state of evidence accumulator
    out( i , 0 ) = std::max(0.0, t - resp*d/2)  + R::runif(t0-st0/2, t0+st0/2);;
    out( i , 1 ) = resp;
    out( i , 3 ) = conf; // evidence term
    
    vis = R::rnorm((tau+t)*muvis, sqrt(svis*svis*(tau+t)+(t+tau)*(t+tau)*sigvis*sigvis));
    
    if (lambda >0) {
      // subj_time = R::rnorm((tau+t)*Timerate, sqrt(Timediff*Timediff*(tau+t)));
      // if (subj_time <=0) {
      //   conf = 9999999999999;
      // } else {
      //   conf = (w*conf + (1-w)*vis)/pow(subj_time, lambda);
      // }
      subj_time = R::rgamma((tau+t)*Timerate*Timerate/(Timediff*Timediff), Timediff*Timediff/Timerate);
      conf = (w*conf + (1-w)*vis)/pow(subj_time, lambda);
    } else {
      conf = (w*conf + (1-w)*vis);
    }
    // save final confidence and value of visibility process
    out( i , 2 ) = conf;
    out( i , 4 ) = vis;
    out( i , 5 ) = mu;
    
    if (i % 200 ==0 ) Rcpp::checkUserInterrupt();
    
  }
  
  return out;
}


// [[Rcpp::export]]
NumericMatrix RNG_WEV_SAT_SubjTime (int n, double a, double v, 
                           NumericVector RT2s,
                           NumericVector params, 
                           double delta=0.01,double maxT=9, 
                           bool stop_on_error=true)
{
  
  
  NumericMatrix out(n, 4);
  
  //double a   =    params[0];
  //double v   =    params[0];
  double t0  =    params[0];
  //double d   =    params[1];
  double szr =    params[1];
  double sv  =    params[2];
  double st0 =    params[3];
  //double zr  =    params[4];
  double tau0 =    params[4];
  double lambda =  params[5];
  double w   =    params[6];
  //double muvis =  params[9];
  double sigvis = params[7];
  double svis =   params[8];
  double Timerate = params[9];
  double Timediff = params[10];
  
  bool valid = true;
  
  //if (a <= 0)                         { valid = false; Rcpp::Rcout << "error: invalid parameter a = " << a << std::endl;  }
  if (szr < 0 || szr > 1)             { valid = false; Rcpp::Rcout << "error: invalid parameter szr = " << szr << std::endl; }
  if (st0 < 0)                        { valid = false; Rcpp::Rcout << "error: invalid parameter st0 = " << st0 << std::endl; }
  if (sv < 0)                         { valid = false; Rcpp::Rcout << "error: invalid parameter sv = " << sv << std::endl; }
  //if (t0 - fabs(0.5*d) - 0.5*st0 < 0) { valid = false; Rcpp::Rcout << "error: invalid parameter combination t0 = " << t0 << ", d = " << d << ", st0 =" << st0 << std::endl; }
  //if (zr - 0.5*szr <= 0)              { valid = false; Rcpp::Rcout << "error: invalid parameter combination zr = " << zr << ", szr = " << szr << std::endl;}
  //if (zr + 0.5*szr >= 1)              { valid = false; Rcpp::Rcout << "error: invalid parameter combination zr = " << zr << ", szr = " << szr << std::endl;}
  if (tau0 < 0)                        { valid = false; Rcpp::Rcout << "error: invalid parameter tau0 = " << tau0 << std::endl;}
  if (w<0 || w > 1)                   { valid = false; Rcpp::Rcout << "error: invalid parameter w = " << w << ", allowed: w in [0,1]" <<  std::endl; }
  if (sigvis < 0)                      { valid = false; Rcpp::Rcout << "error: invalid parameter sigvis = " << sigvis <<  std::endl; }
  if (svis <= 0)                        { valid = false; Rcpp::Rcout << "error: invalid parameter svis = " << svis <<  std::endl; }
  if (lambda < 0)                        { valid = false; Rcpp::Rcout << "error: invalid parameter lambda = " << lambda <<  std::endl; }
  
  if (!valid) {
    if (stop_on_error) { Rcpp::stop("Error validating parameters.\n"); }
    else {return out;}
  }
  
  
  double mu, x0, t, conf, vis, tau, subj_time;
  int resp;
  if (RT2s.length() < n) {
    Rcpp::Rcout << "error: RT2s is not long enough; length of RT2s: " << RT2s.length() << "but required are n=" << n <<  std::endl; 
    Rcpp::stop("RT2s length error\n");
  }
  for (int i=0; i < n; i++) {
    tau = RT2s[i] - tau0 + t0;
    
    mu = R::rnorm(v, sv);
    x0 = a* (0.5+R::runif(-szr/2, szr/2));//R::runif(zr-szr/2, zr+szr/2);
    t = 0;
    while ((x0 > 0) && (x0 < a) && (t < maxT)) {
      x0 = x0 + R::rnorm(delta*mu, sqrt(delta));
      t = t + delta;
    }
    if (x0 >= a) {
      resp = 1;
    } else {
      if (x0 <= 0) {
        resp = -1;
      } else {
        resp = 0;
      }
    }
    //if (w == 1) {
    
    
    // save experimental condition, response, response time
    out( i , 0) = RT2s[i];
    
    out( i , 1) = t  + R::runif(t0-st0/2, t0+st0/2); //std::max(0.0, t - resp*d/2)+
    out( i , 2) = resp;
    
    // state of evidence accumulator
    if (tau <=0 ) tau = 0.00001;
    tau = tau;//+t0
    conf = resp*(x0 + R::rnorm(tau*mu, sqrt(tau))- a/2);//  //use zr here
    // state of visibility accumulator
    vis = R::rnorm((tau+t)*v, sqrt(svis*svis*(tau+t)+(t+tau)*(t+tau)*sigvis*sigvis)); // maybe use muvis here
    if (lambda >0) {
      subj_time = R::rgamma((tau+t)*Timerate*Timerate/(Timediff*Timediff), Timediff*Timediff/Timerate);
      conf = (w*conf + (1-w)*vis)/pow(subj_time, lambda);
    } else {
      conf = (w*conf + (1-w)*vis);
    }
    // save final confidence
    out( i , 3) = conf;
    
    if (i % 200 ==0 ) Rcpp::checkUserInterrupt();
    
  }
  return out;
}

NumericMatrix RNG_WEV_SAT_SubjTime (int n, NumericVector params, NumericVector tau,
                                double delta=0.01,
                                double maxT=9, bool stop_on_error=true)
{
  NumericMatrix out(n, 6);
  
  double a   =    params[0];
  double v   =    params[1];
  double t0  =    params[2];
  double d   =    params[3];
  double szr =    params[4];
  double sv  =    params[5];
  double st0 =    params[6];
  double zr  =    params[7];
  //double tau =    params[8];
  double lambda =  params[9];
  double w   =    params[10];
  double muvis =  params[11];
  double sigvis = params[12];
  double svis =   params[13];
  double Timerate = params[14];
  double Timediff = params[15];
  
  bool valid = true;
  
  if (a <= 0)                         { valid = false; Rcpp::Rcout << "error: invalid parameter a = " << a << std::endl;  }
  if (szr < 0 || szr > 1)             { valid = false; Rcpp::Rcout << "error: invalid parameter szr = " << szr << std::endl; }
  if (st0 < 0)                        { valid = false; Rcpp::Rcout << "error: invalid parameter st0 = " << st0 << std::endl; }
  if (sv < 0)                         { valid = false; Rcpp::Rcout << "error: invalid parameter sv = " << sv << std::endl; }
  if (t0 - fabs(0.5*d) - 0.5*st0 < 0) { valid = false; Rcpp::Rcout << "error: invalid parameter combination t0 = " << t0 << ", d = " << d << ", st0 =" << st0 << std::endl; }
  if (zr - 0.5*szr <= 0)              { valid = false; Rcpp::Rcout << "error: invalid parameter combination zr = " << zr << ", szr = " << szr << std::endl;}
  if (zr + 0.5*szr >= 1)              { valid = false; Rcpp::Rcout << "error: invalid parameter combination zr = " << zr << ", szr = " << szr << std::endl;}
  if (w<0 || w > 1)                   { valid = false; Rcpp::Rcout << "error: invalid parameter w = " << w << ", allowed: w in [0,1]" <<  std::endl; }
  if (sigvis < 0)                      { valid = false; Rcpp::Rcout << "error: invalid parameter sigvis = " << sigvis <<  std::endl; }
  if (svis <= 0)                        { valid = false; Rcpp::Rcout << "error: invalid parameter svis = " << svis <<  std::endl; }
  if (lambda < 0)                        { valid = false; Rcpp::Rcout << "error: invalid parameter lambda = " << lambda <<  std::endl; }
  if (Timediff <= 0)                        { valid = false; Rcpp::Rcout << "error: invalid parameter Timediff = " << Timediff <<  std::endl; }
  
  if (!valid) {
    if (stop_on_error) { Rcpp::stop("Error validating parameters.\n"); }
    else {return out;}
  }
  
  
  double mu, x0, t, conf, vis, subj_time;
  int resp;
  
  for (int i=0; i < n; i++) {
    mu = R::rnorm(v, sv);
    x0 = a* R::runif(zr-szr/2, zr+szr/2);
    t = 0;
    while ((x0 > 0) && (x0 < a) && (t < maxT)) {
      x0 = x0 + R::rnorm(delta*mu, sqrt(delta));
      t = t + delta;
    }
    if (x0 >= a) {
      resp = 1;
    } else {
      if (x0 <= 0) {
        resp = -1;
      } else {
        resp = 0;
      }
    }
    //if (w == 1) {
    if (tau[i] > 0) {
      conf = resp*(x0 + R::rnorm(tau[i]*mu, sqrt(tau[i])) - a*zr);
    } else {
      conf = resp*(x0 - a*zr);
    }
    // save response, response time and state of evidence accumulator
    out( i , 0 ) = std::max(0.0, t - resp*d/2)  + R::runif(t0-st0/2, t0+st0/2);;
    out( i , 1 ) = resp;
    out( i , 3 ) = conf; // evidence term
    
    vis = R::rnorm((tau[i]+t)*muvis, sqrt(svis*svis*(tau[i]+t)+(t+tau[i])*(t+tau[i])*sigvis*sigvis));
    
    if (lambda >0) {
      // subj_time = R::rnorm((tau+t)*Timerate, sqrt(Timediff*Timediff*(tau+t)));
      // if (subj_time <=0) {
      //   conf = 9999999999999;
      // } else {
      //   conf = (w*conf + (1-w)*vis)/pow(subj_time, lambda);
      // }
      subj_time = R::rgamma((tau[i]+t)*Timerate*Timerate/(Timediff*Timediff), Timediff*Timediff/Timerate);
      conf = (w*conf + (1-w)*vis)/pow(subj_time, lambda);
    } else {
      conf = (w*conf + (1-w)*vis);
    }
    // save final confidence and value of visibility process
    out( i , 2 ) = conf;
    out( i , 4 ) = vis;
    out( i , 5 ) = mu;
    
    if (i % 200 ==0 ) Rcpp::checkUserInterrupt();
    
  }
  
  return out;
}

