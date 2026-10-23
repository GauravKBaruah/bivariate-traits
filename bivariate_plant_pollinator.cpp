// Bivariate plant-pollinator ODE: C++ translation of the supplied R model.
// Compile from R with Rcpp::sourceCpp(). RcppArmadillo is NOT required.
// [[Rcpp::plugins(cpp11)]]
#include <vector>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#ifndef BIVARIATE_STANDALONE
#include <Rcpp.h>
#endif

// Four entries are stored explicitly so the 2 x 2 calculations are visible.
struct Matrix2 {
  double a, b, c, d;
};
Matrix2 add(Matrix2 x, Matrix2 y) {
  return {x.a+y.a, x.b+y.b, x.c+y.c, x.d+y.d};
}
double determinant(Matrix2 x) { return x.a*x.d-x.b*x.c; }
Matrix2 inverse(Matrix2 x) {
  double det=determinant(x);
  if (!(det>0) || !std::isfinite(det))
    throw std::runtime_error("Nonpositive covariance determinant");
  return {x.d/det,-x.b/det,-x.c/det,x.a/det};
}
Matrix2 multiply(Matrix2 x, Matrix2 y) {
  return {x.a*y.a+x.b*y.c, x.a*y.b+x.b*y.d,
          x.c*y.a+x.d*y.c, x.c*y.b+x.d*y.d};
}
double cutoff_cpp(double x) {
  if(x<=0) return 0;
  if(x>=1) return 1;
  return x*x*x*(10+x*(-15+6*x));
}

struct CompetitionPair {
  Matrix2 invD;
  double amplitude;
};

struct BivariateModel {
  int animals, plants, species, quadrature_size;
  double gi, ki, bw, a, H, mutualistic_strength, alpha_ii, trait_evolve;
  double inv_w1_squared, inv_w2_squared;
  std::vector<Matrix2> phenotype, genetic_response;
  std::vector<double> link, degree, offset1, offset2, weight;
  std::vector<CompetitionPair> competition;

  // State: Na,Np,muA1,muP1,muA2,muP2 (all species within each block).
  std::vector<double> rhs(const std::vector<double>& state, double temperature) const {
    const int S=species, K=quadrature_size;
    if(static_cast<int>(state.size())!=3*S)
      throw std::runtime_error("State must have length 3*(animals+plants)");
    std::vector<double> N(S), mu1(S), mu2(S), z1(S*K), z2(S*K);
    for(int i=0;i<S;++i) {
      N[i]=std::max(0.0,state[i]);
      mu1[i]=state[S+i]; mu2[i]=state[2*S+i];
      for(int q=0;q<K;++q) {
        z1[i*K+q]=mu1[i]+offset1[i*K+q];
        z2[i*K+q]=mu2[i]+offset2[i*K+q];
      }
    }
    std::vector<double> derivative(3*S,0.0);
    for(int i=0;i<S;++i) {
      // 1. Exact intrinsic-growth average and its centered trait moments.
      double thermal_width=bw-a*mu1[i];
      double total_width=thermal_width*thermal_width+phenotype[i].a;
      double lag=temperature-mu1[i];
      double birth=gi/std::sqrt(total_width)*std::exp(-lag*lag/(2*total_width));
      double mean_growth=birth-ki-alpha_ii*N[i];
      double moment1=phenotype[i].a*birth*lag/total_width;
      double moment2=phenotype[i].c*birth*lag/total_width;

      // 2. Exact interspecific competition averages and centered moments.
      int first=(i<animals) ? 0 : animals;
      int last=(i<animals) ? animals : S;
      for(int j=first;j<last;++j) {
        if(i==j) continue; // self-regulation was included above
        const CompetitionPair& pair=competition[i*S+j];
        double difference1=mu1[i]-mu1[j], difference2=mu2[i]-mu2[j];
        double u1=pair.invD.a*difference1+pair.invD.b*difference2;
        double u2=pair.invD.c*difference1+pair.invD.d*difference2;
        double alpha=pair.amplitude*std::exp(-0.5*(difference1*u1+difference2*u2));
        double beta1=-(phenotype[i].a*u1+phenotype[i].b*u2)*alpha;
        double beta2=-(phenotype[i].c*u1+phenotype[i].d*u2)*alpha;
        mean_growth-=alpha*N[j];
        moment1-=beta1*N[j]; moment2-=beta2*N[j];
      }

      // 3. Type-II mutualism: EXACTLY the supplied nested quadrature rule.
      // q = focal node, j = partner species, r = partner node.
      double mutualistic_growth=0, mutualistic_moment1=0, mutualistic_moment2=0;
      int partner_first=(i<animals) ? animals : 0;
      int partner_last=(i<animals) ? S : animals;
      for(int q=0;q<K;++q) {
        double input=0;
        for(int j=partner_first;j<partner_last;++j) {
          double coefficient=link[i*S+j]*N[j]*mutualistic_strength;
          if(coefficient==0) continue;
          double average_matching=0;
          for(int r=0;r<K;++r) {
            double difference1=z1[i*K+q]-z1[j*K+r];
            double difference2=z2[i*K+q]-z2[j*K+r];
            double matching=std::exp(-difference1*difference1*inv_w1_squared
                                    -difference2*difference2*inv_w2_squared);
            average_matching+=matching*weight[r];
          }
          input+=coefficient*average_matching;
        }
        // Preserve degree scaling in the numerator ONLY.
        double benefit=(input/degree[i])/(1+H*input);
        mutualistic_growth+=benefit*weight[q];
        mutualistic_moment1+=offset1[i*K+q]*benefit*weight[q];
        mutualistic_moment2+=offset2[i*K+q]*benefit*weight[q];
      }
      mean_growth+=mutualistic_growth;
      moment1+=mutualistic_moment1; moment2+=mutualistic_moment2;

      // 4. Abundance and inherited mean-trait responses: G * inverse(P) * moments.
      derivative[i]=N[i]*mean_growth*cutoff_cpp(N[i]/1e-6);
      if(N[i]>0) {
        Matrix2 response=genetic_response[i];
        derivative[S+i]=trait_evolve*(response.a*moment1+response.b*moment2);
        derivative[2*S+i]=trait_evolve*(response.c*moment1+response.d*moment2);
      }
    }
    return derivative;
  }
};

// Cache only quantities that are fixed during a simulation.
void prepare_fixed_calculations(BivariateModel& model,
                                const std::vector<Matrix2>& G,
                                const std::vector<Matrix2>& E,
                                const std::vector<double>& nodes,
                                const std::vector<double>& weights,
                                double wc1,double wc2,double alpha0) {
  int S=model.species, points=nodes.size();
  int K=points*points; model.quadrature_size=K;
  model.phenotype.resize(S); model.genetic_response.resize(S);
  model.offset1.resize(S*K); model.offset2.resize(S*K); model.weight.resize(K);
  for(int i=0;i<S;++i) {
    Matrix2 P=add(G[i],E[i]);
    model.phenotype[i]=P;
    model.genetic_response[i]=multiply(G[i],inverse(P));
    double sigma1=std::sqrt(P.a), sigma2=std::sqrt(P.d);
    double rho=P.b/(sigma1*sigma2);
    double remaining=std::sqrt(std::max(0.0,1-rho*rho));
    int q=0;
    for(int u=0;u<points;++u) {
      for(int v=0;v<points;++v,++q) {
        model.offset1[i*K+q]=sigma1*nodes[u];
        model.offset2[i*K+q]=sigma2*(rho*nodes[u]+remaining*nodes[v]);
        model.weight[q]=weights[u]*weights[v];
      }
    }
  }
  Matrix2 Wc={wc1*wc1/2,0,0,wc2*wc2/2};
  model.competition.resize(S*S);
  for(int i=0;i<S;++i) {
    for(int j=0;j<S;++j) {
      if((i<model.animals)!=(j<model.animals) || i==j) continue;
      Matrix2 D=add(Wc,add(model.phenotype[i],model.phenotype[j]));
      model.competition[i*S+j]={inverse(D),alpha0*std::sqrt(determinant(Wc)/determinant(D))};
    }
  }
}

#ifndef BIVARIATE_STANDALONE
// R stores covariance arrays as [species, row, column], column-major.
Matrix2 covariance_from_R(Rcpp::NumericVector x,int i,int S) {
  return {x[i],x[i+2*S],x[i+S],x[i+3*S]};
}

// [[Rcpp::export]]
SEXP prepare_bivariate_cpp(Rcpp::List pars,Rcpp::NumericVector nodes,
                           Rcpp::NumericVector weights) {
  Rcpp::XPtr<BivariateModel> model(new BivariateModel(),true);
  Rcpp::NumericMatrix net=pars["matrix"], preference=pars["pref_matrix"];
  model->animals=net.ncol(); model->plants=net.nrow();
  model->species=model->animals+model->plants;
  int A=model->animals, P=model->plants, S=model->species;
  model->gi=Rcpp::as<double>(pars["gi"]);
  model->ki=Rcpp::as<double>(pars["ki"]);
  model->bw=Rcpp::as<double>(pars["bw"]);
  model->a=Rcpp::as<double>(pars["a"]);
  model->H=Rcpp::as<double>(pars["H"]);
  model->mutualistic_strength=Rcpp::as<double>(pars["mut.strength"]);
  model->alpha_ii=Rcpp::as<double>(pars["alpha_ii"]);
  model->trait_evolve=Rcpp::as<double>(pars["trait_evolve"]);
  Rcpp::NumericVector w=pars["w"], wc=pars["w_c"];
  model->inv_w1_squared=1/(w[0]*w[0]);
  model->inv_w2_squared=1/(w[1]*w[1]);
  model->link.assign(S*S,0.0); model->degree.assign(S,1.0);
  for(int a=0;a<A;++a) for(int p=0;p<P;++p) {
    double value=net(p,a)*preference(p,a);
    model->link[a*S+(A+p)]=value;
    model->link[(A+p)*S+a]=value;
  }
  if(Rcpp::as<bool>(pars["degree_scale"])) {
    Rcpp::NumericVector da=pars["dganimals"],dp=pars["dgplants"];
    for(int a=0;a<A;++a) model->degree[a]=std::max(1.0,da[a]);
    for(int p=0;p<P;++p) model->degree[A+p]=std::max(1.0,dp[p]);
  }
  Rcpp::NumericVector ga=pars["G_A"],gp=pars["G_P"],ea=pars["E_A"],ep=pars["E_P"];
  std::vector<Matrix2> G(S),E(S);
  for(int a=0;a<A;++a) {G[a]=covariance_from_R(ga,a,A);E[a]=covariance_from_R(ea,a,A);}
  for(int p=0;p<P;++p) {G[A+p]=covariance_from_R(gp,p,P);E[A+p]=covariance_from_R(ep,p,P);}
  prepare_fixed_calculations(*model,G,E,Rcpp::as<std::vector<double> >(nodes),
                            Rcpp::as<std::vector<double> >(weights),wc[0],wc[1],
                            Rcpp::as<double>(pars["alpha0"]));
  return model;
}

// [[Rcpp::export]]
Rcpp::NumericVector bivariate_rhs_cpp(SEXP pointer,Rcpp::NumericVector state,
                                     double temperature) {
  Rcpp::XPtr<BivariateModel> model(pointer);
  if(!std::isfinite(temperature)) Rcpp::stop("Temperature must be finite");
  for(double value:state) if(!std::isfinite(value)) Rcpp::stop("Nonfinite ODE state");
  return Rcpp::wrap(model->rhs(Rcpp::as<std::vector<double> >(state),temperature));
}
#endif
