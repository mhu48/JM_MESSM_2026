# =============================================================================
# Joint Log-Likelihoods: Cubic-Cubic and Cubic-Local Models
# =============================================================================
# Rcpp/inline implementations of the full joint log-likelihood for the two
# alternative individual trajectory models used in the sensitivity analysis.
# Adapted from llh_survival_x_OU.R; only the state-space matrices differ.
#
# Defines two R functions:
#   llh_survival_x_cs_cs(theta_trans, data_split, z, D_in, X_cube, A_list, AB_list)
#     Joint likelihood for cubic spline population + cubic spline individual.
#     Individual state alpha_v is 2q-dimensional (position + slope per biomarker).
#     D_in = [I_q x (1,0), I_q x (1,0)] selecting position components from both processes.
#
#   llh_survival_x_cs_local(theta_trans, data_split, z, D_in, X_cube, A_list, AB_list)
#     Joint likelihood for cubic spline population + local level individual.
#     Individual state alpha_v is q-dimensional (position only).
#     D_in = [I_q x (1,0), I_q] same as OU model.
#
# Full parameter vector theta_trans: see llh_survival_x_OU.R header.
# Key difference from OU: no xi (mean-reversion) parameter; lambda replaces nu2.
#
# Inputs (R-level): same as llh_survival_x_OU.R
# =============================================================================

library(Rcpp)
library(inline)

## ============================================================
## CUBIC-CUBIC JOINT MODEL
## ============================================================

# First, source the longitudinal-only likelihood
source("llh_cs_cs.R")  # Creates llh_cs function

# Now create joint version by adapting llh_survival_x_OU
# Key changes: d_v = 2*q instead of q, T_v is cubic spline transition
rcpp_inc_cs_cs_joint = '
#include<iostream>
using namespace Rcpp;
using namespace arma;
using namespace std;

mat my_inv(mat& Q);
mat my_chol(mat& Q);
vec delta_est_cs_cs(vec& theta_trans, cube& data_split_cube, List& A_list);
double neg2llh_joint_cs_cs(vec& theta_trans, cube& data_split_cube, mat& z_mat, mat& D_mat, cube& X00_cube, List& A_list, List& AB_list);

// Helper functions (same as OU version)
mat my_inv(mat& Q){
  vec eigval, sgn;
  mat eigvec;
  Q = (Q + Q.t()) / 2.0;
  eig_sym(eigval, eigvec, Q);
  sgn = sign(eigval);
  for (int i = 0; i < eigval.n_elem; ++i){
    if (abs(eigval(i)) < 0.0000000000000001) eigval(i) = sgn(i) * 0.0000000001;
  }
  return eigvec * diagmat(1 / eigval) * eigvec.t();
}

mat my_chol(mat& Q){
  vec eigval;
  mat eigvec;
  Q = (Q + Q.t()) / 2.0;
  eig_sym(eigval, eigvec, Q);
  for (int i = 0; i < eigval.n_elem; ++i){
    if (eigval(i) < 0) eigval(i) = 0;
  }
  return eigvec * diagmat(sqrt(eigval));
}

// Delta estimation for cubic-cubic model
vec delta_est_cs_cs(vec& theta_trans, cube& data_split_cube, List& A_list){
  
  vec theta, zeta_vec, lambda_vec, d_vec, l_vec, sum_v_a, delta_term2, delta_hat, sum_V_A,
      a_u_a, a_u_a_tt;
  
  mat y, D, L, tmp, tmp_lower, Z00, T00, P0, P1, P2, P3, Z_u, Z_v, T_u, T_v,
      Sigma, Lambda, Sig_u, Sig_v,
      F1mF2, F1, F2, A, B, AmB, C0, C1, C2, M0, M1, M2, M3, 
      a_a_tt, P0_tt, P1_tt, P2_tt, P3_tt,
      Za_a, ZA_A, a_a, v_a, a_v_a, a_v_a_tt, A_u_A, A_u_A_tt,
      mat_1m, delta_term1, sum_V_d, sum_VFV, sum_VFv, V_temp;
  
  cube A_A, A_A_tt, V_A, A_v_A, A_v_A_tt;
  
  double dt, kappa;
  
  int q, npt, m, mt, d_u, d_v, from, to;
  
  uvec idx, A_vec;
  
  // Initialize
  q = data_split_cube.n_rows;
  npt = data_split_cube.n_cols;
  m = data_split_cube.n_slices;
  d_u = 2 * q;
  d_v = 2 * q;  // CUBIC SPLINE for individuals
  mat_1m = mat(1, m, fill::ones);
  
  // Transform parameters
  theta = theta_trans;
  from = 0;
  to = 3 * q - 1;
  theta.subvec(from, to) = exp(theta_trans.subvec(from, to));
  from = 0.5 * q * q + 2.5 * q;
  theta(from) = exp(theta_trans(from));
  
  // Extract parameters
  from = 0; to = q - 1;
  zeta_vec = theta.subvec(from, to);
  
  from = q; to = 2 * q - 1;
  lambda_vec = theta.subvec(from, to);  // smoothing for individuals
  
  from = 2 * q; to = 3 * q - 1;
  d_vec = theta.subvec(from, to);
  D = diagmat(d_vec);
  
  if (q > 1){
    from = 3 * q;
    to = 0.5 * q * q + 2.5 * q - 1;
    l_vec = theta.subvec(from, to);
    L = eye(q, q);
    tmp = ones(size(L));
    tmp_lower = trimatl(tmp);
    tmp_lower.diag().zeros();
    idx = find(tmp_lower > 0);
    L.elem(idx) = l_vec;
  } else L = eye(1, 1);
  
  from = 0.5 * q * q + 2.5 * q;
  kappa = theta(from);
  
  // System matrices
  Z00 = zeros(1, 2);
  Z00(0, 0) = 1;
  Z_u = kron(eye(q, q), Z00);
  Z_v = kron(eye(q, q), Z00);  // CUBIC SPLINE
  
  dt = 1.0;
  T00 = eye(2, 2);
  T00(0, 1) = dt;
  T_u = kron(eye(q, q), T00);
  T_v = kron(eye(q, q), T00);  // CUBIC SPLINE transition
  
  Lambda = zeros(2, 2);
  Lambda(0, 0) = dt * dt * dt / 3;
  Lambda(0, 1) = Lambda(1, 0) = dt * dt / 2;
  Lambda(1, 1) = dt;
  Sig_u = kron(diagmat(zeta_vec), Lambda);
  Sig_v = kron(diagmat(lambda_vec), Lambda);  // CUBIC SPLINE noise
  
  // Initialize
  P0 = zeros(d_u, d_u);
  P1 = zeros(d_u, d_v);
  P3 = zeros(d_v, d_v);
  P2 = kappa * eye(d_v, d_v);
  
  a_u_a = zeros(d_u, 1);
  a_v_a = zeros(d_v, m);
  
  A_u_A = eye(d_u, d_u);
  A_v_A = zeros(d_v, m, d_u);
  A_u_A_tt = zeros(d_u, d_u);
  A_v_A_tt = zeros(d_v, m, d_u);
  V_A = zeros(q, m, d_u);
  
  Sigma = L * D * L.t();
  
  delta_term1 = mat(d_u, d_u, fill::zeros);
  delta_term2 = vec(d_u, fill::zeros);
  
  // Filtering loop for delta estimation
  for (int j = 0; j < npt; ++j){
    A_vec = as<uvec>(A_list[j]);
    A_vec -= 1;
    mt = A_vec.n_elem;
    
    y = data_split_cube.subcube(0, j, 0, q - 1, j, m - 1);
    Za_a = kron(mat_1m, Z_u * a_u_a) + Z_v * a_v_a;
    v_a = y - Za_a;
    
    for (int k = 0; k < d_u; ++k){
      ZA_A = kron(mat_1m, Z_u * A_u_A.col(k)) + Z_v * A_v_A.slice(k);
      V_A.slice(k) = - ZA_A;
    }
    
    // F components
    A = Z_v * P2 * Z_v.t() + Sigma;
    B = Z_u * P0 * Z_u.t() + Z_u * P1 * Z_v.t() + Z_v * P1.t() * Z_u.t() + Z_v * P3 * Z_v.t();
    
    F1 = inv(A);
    AmB = A + mt * B;
    F2 = - F1 * B * inv(AmB);
    F1 = (F1 + F1.t()) / 2.0;
    F2 = (F2 + F2.t()) / 2.0;
    F1mF2 = F1 + mt * F2;
    
    sum_v_a = sum(v_a.cols(A_vec), 1);
    
    C0 = (P0 * Z_u.t() + P1 * Z_v.t()) * F1mF2;
    C1 = P2 * Z_v.t() * F1;
    C2 = (P1.t() * Z_u.t() + P3 * Z_v.t()) * F1mF2 + P2 * Z_v.t() * F2;
    
    a_u_a_tt = a_u_a + C0 * sum_v_a;
    a_v_a_tt = a_v_a + C1 * v_a + kron(mat_1m, C2 * sum_v_a);
    
    for (int k = 0; k < d_u; ++k){
      sum_V_A = sum(V_A.slice(k).cols(A_vec), 1);
      A_u_A_tt.col(k) = A_u_A.col(k) + C0 * sum_V_A;
      A_v_A_tt.slice(k) = A_v_A.slice(k) + C1 * V_A.slice(k) + kron(mat_1m, C2 * sum_V_A);
    }
    
    M0 = mt * (P0 * Z_u.t() + P1 * Z_v.t()) * F1mF2 * (Z_u * P0 + Z_v * P1.t());
    M1 = (P0 * Z_u.t() + P1 * Z_v.t()) * F1mF2 * (mt * Z_u * P1 + Z_v * (P2 + mt * P3));
    M2 = P2 * Z_v.t() * F1 * Z_v * P2;
    M3 = P2 * Z_v.t() * F2 * Z_v * P2 + 
         P2 * Z_v.t() * F1mF2 * (Z_u * P1 + Z_v * P3) + 
         (P1.t() * Z_u.t() + P3 * Z_v.t()) * F1mF2 * Z_v * P2 +
         mt * (P1.t() * Z_u.t() + P3 * Z_v.t()) * F1mF2 * (Z_u * P1 + Z_v * P3);
    
    P0_tt = P0 - M0;
    P1_tt = P1 - M1;
    P2_tt = P2 - M2;
    P3_tt = P3 - M3;
    
    // Delta terms
    sum_V_d = mat(q, d_u, fill::zeros);
    sum_VFV = mat(d_u, d_u, fill::zeros);
    sum_VFv = vec(d_u, fill::zeros);
    for (int i = 0; i < mt; ++i){
      V_temp = mat(V_A.subcube(0, A_vec[i], 0, q - 1, A_vec[i], d_u - 1));
      sum_V_d += V_temp;
      sum_VFV += V_temp.t() * F1 * V_temp;
      sum_VFv += V_temp.t() * F1 * v_a.col(A_vec[i]);
    }
    
    delta_term1 += sum_VFV + sum_V_d.t() * F2 * sum_V_d;
    delta_term2 += sum_VFv + sum_V_d.t() * F2 * sum_v_a;
    
    a_u_a = T_u * a_u_a_tt;
    a_v_a = T_v * a_v_a_tt;
    
    for (int k = 0; k < d_u; ++k){
      A_u_A.col(k) = T_u * A_u_A_tt.col(k);
      A_v_A.slice(k) = T_v * A_v_A_tt.slice(k);
    }
    
    P0 = T_u * P0_tt * T_u.t() + Sig_u;
    P1 = T_u * P1_tt * T_v.t();
    P2 = T_v * P2_tt * T_v.t() + Sig_v;
    P3 = T_v * P3_tt * T_v.t();
    
    P0 = (P0 + P0.t()) / 2.0;
    P2 = (P2 + P2.t()) / 2.0;
    P3 = (P3 + P3.t()) / 2.0;
  }
  
  delta_hat = - inv(delta_term1) * delta_term2;
  return delta_hat;
}

// Joint likelihood for cubic-cubic model
double neg2llh_joint_cs_cs(vec& theta_trans, cube& data_split_cube, mat& z_mat, mat& D_mat, cube& X00_cube, List& A_list, List& AB_list){
  
  vec theta, gamma_vec, gamma_x, delta, zeta_vec, lambda_vec, d_vec, l_vec, sum_v, Z_nonNA, p_z_sim, cond_mean, sum_diff,
      theta_t_plus_1, z_term1, z_term2, z_term, a_u, a_u_tt, a_u_tt_nonNA, alpha_u_gen;
  
  mat data_split_j, D, L, tmp, tmp_lower, Z00, T00, P0, P1, P2, P3, Z_u, Z_v, T_u, T_v,
    Sigma, Lambda, Sig_u, Sig_v, P2_tt_inv, temp, var_mat, a_v, a_v_tt, a_v_tt_nonNA, alpha_v_gen,
    F1mF2, F1, F2, A, B, AmB, C0, C1, C2, M0, M1, M2, M3, 
    P0_tt, P1_tt, P2_tt, P3_tt, cond_chol, mat_1mt_next, stacked_mat, stacked_mat2, Za, v_split, mat_1m, S_u, chol_terms_u;
  
  int q, p, npt, m, mt, mt_next, d_u, d_v, from, to, nsim = 100, dim_gamma;
  
  double dt, llh, llh_const, llh_term1, llh_term2, p_z, p_hat, kappa;
  
  uvec idx, A_vec, AB_vec_next;
  
  cube S_v, chol_terms_v;
  
  // Initialize
  q = data_split_cube.n_rows;
  npt = data_split_cube.n_cols;
  m = data_split_cube.n_slices;
  p = X00_cube.n_cols;
  d_u = 2 * q;
  d_v = 2 * q;  // KEY DIFFERENCE: cubic spline for individuals
  mat_1m = mat(1, m, fill::ones);
  dim_gamma = D_mat.n_rows;
  
  delta = delta_est_cs_cs(theta_trans, data_split_cube, A_list);
  
  // Transform parameters
  theta = theta_trans;
  from = 0;
  to = 3 * q - 1;
  theta.subvec(from, to) = exp(theta_trans.subvec(from, to));
  from = 0.5 * q * q + 2.5 * q;
  theta(from) = exp(theta_trans(from));
  
  // Extract parameters
  from = 0; to = q - 1;
  zeta_vec = theta.subvec(from, to);
  
  from = q; to = 2 * q - 1;
  lambda_vec = theta.subvec(from, to);  // smoothing params for individuals
  
  from = 2 * q; to = 3 * q - 1;
  d_vec = theta.subvec(from, to);
  D = diagmat(d_vec);
  
  if (q > 1){
    from = 3 * q;
    to = 0.5 * q * q + 2.5 * q - 1;
    l_vec = theta.subvec(from, to);
    L = eye(q, q);
    tmp = ones(size(L));
    tmp_lower = trimatl(tmp);
    tmp_lower.diag().zeros();
    idx = find(tmp_lower > 0);
    L.elem(idx) = l_vec;
  } else L = eye(1, 1);
  
  from = 0.5 * q * q + 2.5 * q;
  kappa = theta(from);
  
  // Gamma parameters
  from = 0.5 * q * q + 2.5 * q + 1;
  to = 0.5 * q * q + 2.5 * q + dim_gamma;
  gamma_vec = theta.subvec(from, to);
  
  from = 0.5 * q * q + 2.5 * q + dim_gamma + 1;
  to = 0.5 * q * q + 2.5 * q + dim_gamma + p + 1;
  gamma_x = theta.subvec(from, to);
  
  // System matrices
  Z00 = zeros(1, 2);
  Z00(0, 0) = 1;
  Z_u = Z_v = kron(eye(q, q), Z00);  // Both cubic
  
  dt = 1.0;
  T00 = eye(2, 2);
  T00(0, 1) = dt;
  T_u = T_v = kron(eye(q, q), T00);  // KEY DIFFERENCE: T_v is cubic
  
  Lambda = zeros(2, 2);
  Lambda(0, 0) = dt * dt * dt / 3;
  Lambda(0, 1) = Lambda(1, 0) = dt * dt / 2;
  Lambda(1, 1) = dt;
  Sig_u = kron(diagmat(zeta_vec), Lambda);
  Sig_v = kron(diagmat(lambda_vec), Lambda);  // KEY DIFFERENCE: cubic noise
  
  // Initialize
  P0 = zeros(d_u, d_u);
  P1 = zeros(d_u, d_v);
  P3 = zeros(d_v, d_v);
  P2 = kappa * eye(d_v, d_v);
  
  a_u = delta;
  a_v = zeros(d_v, m);
  
  Sigma = L * D * L.t();
  llh = 0.0;
  
  // FILTERING LOOP (same structure as OU)
  for (int j = 0; j < npt; ++j){
    A_vec = as<uvec>(A_list[j]);
    A_vec -= 1;
    mt = A_vec.n_elem;
    
    data_split_j = data_split_cube.subcube(0, j, 0, q - 1, j, m - 1);
    Za = kron(mat_1m, Z_u * a_u) + Z_v * a_v;
    v_split = data_split_j - Za;
    
    // F components
    A = Z_v * P2 * Z_v.t() + Sigma;
    B = Z_u * P0 * Z_u.t() + Z_u * P1 * Z_v.t() + Z_v * P1.t() * Z_u.t() + Z_v * P3 * Z_v.t();
    
    F1 = inv(A);
    AmB = A + mt * B;
    F2 = - F1 * B * inv(AmB);
    F1 = (F1 + F1.t()) / 2.0;
    F2 = (F2 + F2.t()) / 2.0;
    F1mF2 = F1 + mt * F2;
    
    sum_v = sum(v_split.cols(A_vec), 1);
    
    // Update
    C0 = (P0 * Z_u.t() + P1 * Z_v.t()) * F1mF2;
    C1 = P2 * Z_v.t() * F1;
    C2 = (P1.t() * Z_u.t() + P3 * Z_v.t()) * F1mF2 + P2 * Z_v.t() * F2;
    
    a_u_tt = a_u + C0 * sum_v;
    a_v_tt = a_v + C1 * v_split + kron(mat_1m, C2 * sum_v);
    
    // P update
    M0 = mt * (P0 * Z_u.t() + P1 * Z_v.t()) * F1mF2 * (Z_u * P0 + Z_v * P1.t());
    M1 = (P0 * Z_u.t() + P1 * Z_v.t()) * F1mF2 * (mt * Z_u * P1 + Z_v * (P2 + mt * P3));
    M2 = P2 * Z_v.t() * F1 * Z_v * P2;
    M3 = P2 * Z_v.t() * F2 * Z_v * P2 + 
         P2 * Z_v.t() * F1mF2 * (Z_u * P1 + Z_v * P3) + 
         (P1.t() * Z_u.t() + P3 * Z_v.t()) * F1mF2 * Z_v * P2 +
         mt * (P1.t() * Z_u.t() + P3 * Z_v.t()) * F1mF2 * (Z_u * P1 + Z_v * P3);
    
    P0_tt = P0 - M0;
    P1_tt = P1 - M1;
    P2_tt = P2 - M2;
    P3_tt = P3 - M3;
    
    // Longitudinal likelihood
    llh_term1 = 0.0;
    for (int k = 0; k < mt; ++k){
      llh_term1 += as_scalar(v_split.col(A_vec(k)).t() * F1 * v_split.col(A_vec(k)));
    }
    llh_term2 = as_scalar(sum_v.t() * F2 * sum_v);
    llh_const = - (q * mt) / 2.0 * log(2.0 * datum::pi);
    
    if (det(F1mF2) > 0 & det(F1) > 0)
      llh += llh_const + 0.5 * (log(det(F1mF2)) + (mt - 1.0) * log(det(F1))) - 0.5 * (llh_term1 + llh_term2);
    
    // SURVIVAL LIKELIHOOD (same as OU, just with d_v = 2*q)
    if ((j < npt - 1) & (j > 0)){
      AB_vec_next = as<uvec>(AB_list[j + 1]);
      AB_vec_next -= 1;
      mt_next = AB_vec_next.n_elem;
      if (mt_next <= 0) break;
      
      S_v = zeros(d_v, d_v, mt_next);
      chol_terms_v = zeros(d_v, d_v, mt_next);
      
      P2_tt_inv = inv(P2_tt);
      P2_tt_inv = (P2_tt_inv + P2_tt_inv.t()) / 2.0;
      
      temp = P2_tt + mt_next * P3_tt;
      S_u = P1_tt * (P2_tt_inv - mt_next * P2_tt_inv * P3_tt * inv(temp));
      var_mat = P0_tt - mt_next * S_u * P1_tt.t();
      var_mat = (var_mat + var_mat.t()) / 2.0;
      chol_terms_u = my_chol(var_mat);
      
      if (mt_next > 2){
        for (int k = 1; k <= mt_next - 2; ++k){
          temp = P2_tt + (mt_next - k) * P3_tt;
          S_v.slice(k - 1) = P3_tt * (P2_tt_inv - (mt_next - k) * P2_tt_inv * P3_tt * inv(temp));
          var_mat = (P2_tt + P3_tt) - (mt_next - k) * S_v.slice(k - 1) * P3_tt;
          var_mat = (var_mat + var_mat.t()) / 2.0;
          chol_terms_v.slice(k - 1) = my_chol(var_mat);
        }
      }
      
      if (mt_next > 1){
        temp = P2_tt + P3_tt;
        S_v.slice(mt_next - 2) = P3_tt * inv(temp);
        var_mat = (P2_tt + P3_tt) - S_v.slice(mt_next - 2) * P3_tt;
        var_mat = (var_mat + var_mat.t()) / 2.0;
        chol_terms_v.slice(mt_next - 2) = my_chol(var_mat);
      }
      
      var_mat = P2_tt + P3_tt;
      var_mat = (var_mat + var_mat.t()) / 2.0;
      chol_terms_v.slice(mt_next - 1) = my_chol(var_mat);
      
      a_u_tt_nonNA = a_u_tt;
      a_v_tt_nonNA = a_v_tt.cols(AB_vec_next);
      
      Z_nonNA = vec(mt_next);
      for (int i = 0; i < mt_next; ++i){
        Z_nonNA(i) = as_scalar(z_mat(AB_vec_next(i), j + 1));
      }
      
      alpha_v_gen = mat(d_v, mt_next);
      p_z_sim = vec(2 * nsim, fill::zeros);
      
      for (int k = 0; k < nsim; ++k){
        cond_mean = a_v_tt_nonNA.col(mt_next - 1);
        cond_chol = chol_terms_v.slice(mt_next - 1);
        alpha_v_gen.col(mt_next - 1) = cond_mean + cond_chol * randn(d_v);
        
        sum_diff = alpha_v_gen.col(mt_next - 1) - a_v_tt_nonNA.col(mt_next - 1);
        for (int i = mt_next - 1; i >= 1; --i){
          sum_diff += alpha_v_gen.col(i) - a_v_tt_nonNA.col(i);
          cond_mean = a_v_tt_nonNA.col(i - 1) + S_v.slice(i - 1) * sum_diff;
          cond_chol = chol_terms_v.slice(i - 1);
          alpha_v_gen.col(i - 1) = cond_mean + cond_chol * randn(d_v);
        }
        
        sum_diff += alpha_v_gen.col(0) - a_v_tt_nonNA.col(0);
        cond_mean = a_u_tt_nonNA + S_u * sum_diff;
        cond_chol = chol_terms_u;
        alpha_u_gen = cond_mean + cond_chol * randn(d_u);
        
        mat_1mt_next = mat(1, mt_next, fill::ones);
        stacked_mat = zeros(d_u + d_v, mt_next);
        stacked_mat.submat(0, 0, d_u - 1, mt_next - 1) = kron(mat_1mt_next, alpha_u_gen);
        stacked_mat.submat(d_u, 0, d_u + d_v - 1, mt_next - 1) = alpha_v_gen;
        
        from = 1; to = p;
        theta_t_plus_1 = stacked_mat.t() * D_mat.t() * gamma_vec + 
                         X00_cube.slice(j).rows(AB_vec_next) * gamma_x.subvec(from, to) + gamma_x(0);
        
        z_term1 = (Z_nonNA - 1) % theta_t_plus_1;
        z_term2 = log(1 + exp(- theta_t_plus_1));
        z_term = z_term1 - z_term2;
        p_z = exp(sum(z_term));
        p_z_sim(2 * k) = p_z;
        
        stacked_mat2 = zeros(d_u + d_v, mt_next);
        stacked_mat2.submat(0, 0, d_u - 1, mt_next - 1) = kron(mat_1mt_next, 2 * a_u_tt_nonNA - alpha_u_gen);
        stacked_mat2.submat(d_u, 0, d_u + d_v - 1, mt_next - 1) = 2 * a_v_tt_nonNA - alpha_v_gen;
        
        theta_t_plus_1 = stacked_mat2.t() * D_mat.t() * gamma_vec + 
                         X00_cube.slice(j).rows(AB_vec_next) * gamma_x.subvec(from, to) + gamma_x(0);
        
        z_term1 = (Z_nonNA - 1) % theta_t_plus_1;
        z_term2 = log(1 + exp(- theta_t_plus_1));
        z_term = z_term1 - z_term2;
        p_z = exp(sum(z_term));
        p_z_sim(2 * k + 1) = p_z;
      }
      
      p_hat = mean(p_z_sim);
      if (p_hat > 0) llh += log(p_hat);
    }
    
    // Predict next state
    a_u = T_u * a_u_tt;
    a_v = T_v * a_v_tt;
    
    P0 = T_u * P0_tt * T_u.t() + Sig_u;
    P1 = T_u * P1_tt * T_v.t();
    P2 = T_v * P2_tt * T_v.t() + Sig_v;
    P3 = T_v * P3_tt * T_v.t();
    
    P0 = (P0 + P0.t()) / 2.0;
    P2 = (P2 + P2.t()) / 2.0;
    P3 = (P3 + P3.t()) / 2.0;
  }
  
  return -2.0 * llh;
}
'

src_cs_cs = '
vec theta_trans;
mat z_mat, D_mat;
cube data_split_cube, X00_cube;
List A_list, AB_list;

theta_trans = as<vec>(theta_trans_in);
z_mat = as<mat>(z_in);
D_mat = as<mat>(D_in);
data_split_cube = as<cube>(data_split_in);
X00_cube = as<cube>(X00_cube_in);
A_list = as<List>(A_list_in);
AB_list = as<List>(AB_list_in);

return wrap(neg2llh_joint_cs_cs(theta_trans, data_split_cube, z_mat, D_mat, X00_cube, A_list, AB_list));
'

llh_survival_x_cs_cs = cxxfunction(
  signature(theta_trans_in = "numeric", data_split_in = "numeric", z_in = "numeric", 
            D_in = "numeric", X00_cube_in = "numeric", A_list_in = "numeric", AB_list_in = "numeric"),
  includes = rcpp_inc_cs_cs_joint, body = src_cs_cs, plugin = "RcppArmadillo"
)

## ============================================================
## CUBIC-LOCAL JOINT MODEL
## ============================================================

# Source the longitudinal-only likelihood
source("llh_cs_local.R")  # Creates llh_cs function

# Key changes: d_v = q, T_v = I_q (no dynamics), Sig_v = diag(lambda)
rcpp_inc_cs_local_joint = '
#include<iostream>
using namespace Rcpp;
using namespace arma;
using namespace std;

mat my_inv(mat& Q);
mat my_chol(mat& Q);
vec delta_est_cs_local(vec& theta_trans, cube& data_split_cube, List& A_list);
double neg2llh_cs_local(vec& theta_trans, cube& data_split_cube, mat& z_mat, mat& D_mat, cube& X00_cube, List& A_list, List& AB_list);

// implementation of my_inv() ---------------------------------------------------------
  mat my_inv(mat& Q){
    vec eigval, sgn;
    mat eigvec;
    
    Q = (Q + Q.t()) / 2.0;
    
    eig_sym(eigval, eigvec, Q);
    sgn = sign(eigval);
    for (int i = 0; i < eigval.n_elem; ++i){
      if (abs(eigval(i)) < 0.0000000000000001) eigval(i) = sgn(i) * 0.0000000001;
    }
    return eigvec * diagmat(1 / eigval) * eigvec.t();
  }

// implementation of my_chol() ----------------------------------------------------------------------------
  mat my_chol(mat& Q){
    vec eigval;
    mat eigvec;
    Q = (Q + Q.t()) / 2.0;
    eig_sym(eigval, eigvec, Q);
    for (int i = 0; i < eigval.n_elem; ++i){
      if (eigval(i) < 0) eigval(i) = 0;
    }
    return eigvec * diagmat(sqrt(eigval));
  }

// implementation of function delta_est_cs_local ---------------------------------------------------------------------------
  vec delta_est_cs_local(vec& theta_trans, cube& data_split_cube, List& A_list){
    
    vec theta, zeta_vec, lambda_vec, d_vec, l_vec, sum_v_a, delta_term2, delta_hat, sum_V_A,
    a_u_a, a_u_a_tt;
    
    mat y, D, L, tmp, tmp_lower, Z00, T00, P0, P1, P2, P3, Z_u, Z_v, T_u, T_v,
    Sigma, Lambda, Sig_u, Sig_v,
    F1mF2, F1, F2, A, B, AmB, C0, C1, C2, M0, M1, M2, M3, 
    a_a_tt, P0_tt, P1_tt, P2_tt, P3_tt,
    Za_a, ZA_A, a_a, v_a, a_v_a, a_v_a_tt, A_u_A, A_u_A_tt,
    mat_1m, delta_term1, sum_V_d, sum_VFV, sum_VFv, V_temp;
    
    cube A_A, A_A_tt, V_A, A_v_A, A_v_A_tt;
    
    double dt, kappa;
    
    int q, npt, m, mt, d_u, d_v, from, to;
    
    uvec idx, A_vec;
    
    // Initialize
    q = data_split_cube.n_rows;
    npt = data_split_cube.n_cols;
    m = data_split_cube.n_slices;
    d_u = 2 * q;
    d_v = q;  // LOCAL LEVEL: just q, not 2*q
    mat_1m = mat(1, m, fill::ones);
    
    // Transform parameters
    theta = theta_trans;
    from = 0;
    to = 3 * q - 1;
    theta.subvec(from, to) = exp(theta_trans.subvec(from, to));
    from = 0.5 * q * q + 2.5 * q;
    theta(from) = exp(theta_trans(from));
    
    // Extract parameters
    from = 0; to = q - 1;
    zeta_vec = theta.subvec(from, to);
    
    from = q; to = 2 * q - 1;
    lambda_vec = theta.subvec(from, to);  // variances for local level
    
    from = 2 * q; to = 3 * q - 1;
    d_vec = theta.subvec(from, to);
    D = diagmat(d_vec);
    
    if (q > 1){
      from = 3 * q;
      to = 0.5 * q * q + 2.5 * q - 1;
      l_vec = theta.subvec(from, to);
      L = eye(q, q);
      tmp = ones(size(L));
      tmp_lower = trimatl(tmp);
      tmp_lower.diag().zeros();
      idx = find(tmp_lower > 0);
      L.elem(idx) = l_vec;
    } else L = eye(1, 1);
    
    from = 0.5 * q * q + 2.5 * q;
    kappa = theta(from);
    
    // System matrices
    Z00 = zeros(1, 2);
    Z00(0, 0) = 1;
    Z_u = kron(eye(q, q), Z00);
    Z_v = eye(q, q);  // LOCAL LEVEL: just identity
    
    dt = 1.0;
    T00 = eye(2, 2);
    T00(0, 1) = dt;
    T_u = kron(eye(q, q), T00);
    T_v = eye(q, q);  // LOCAL LEVEL: random walk, no dynamics
    
    Lambda = zeros(2, 2);
    Lambda(0, 0) = dt * dt * dt / 3;
    Lambda(0, 1) = Lambda(1, 0) = dt * dt / 2;
    Lambda(1, 1) = dt;
    Sig_u = kron(diagmat(zeta_vec), Lambda);
    Sig_v = diagmat(lambda_vec);  // LOCAL LEVEL: simple variance
    
    // Initialize
    P0 = zeros(d_u, d_u);
    P1 = zeros(d_u, d_v);
    P3 = zeros(d_v, d_v);
    P2 = kappa * eye(d_v, d_v);
    
    a_u_a = zeros(d_u, 1);
    a_v_a = zeros(d_v, m);
    
    A_u_A = eye(d_u, d_u);
    A_v_A = zeros(d_v, m, d_u);
    A_u_A_tt = zeros(d_u, d_u);
    A_v_A_tt = zeros(d_v, m, d_u);
    V_A = zeros(q, m, d_u);
    
    Sigma = L * D * L.t();
    
    delta_term1 = mat(d_u, d_u, fill::zeros);
    delta_term2 = vec(d_u, fill::zeros);
    
    // Filtering loop for delta estimation
    for (int j = 0; j < npt; ++j){
      A_vec = as<uvec>(A_list[j]);
      A_vec -= 1;
      mt = A_vec.n_elem;
      
      y = data_split_cube.subcube(0, j, 0, q - 1, j, m - 1);
      Za_a = kron(mat_1m, Z_u * a_u_a) + Z_v * a_v_a;
      v_a = y - Za_a;
      
      for (int k = 0; k < d_u; ++k){
        ZA_A = kron(mat_1m, Z_u * A_u_A.col(k)) + Z_v * A_v_A.slice(k);
        V_A.slice(k) = - ZA_A;
      }
      
      // F components
      A = Z_v * P2 * Z_v.t() + Sigma;
      B = Z_u * P0 * Z_u.t() + Z_u * P1 * Z_v.t() + Z_v * P1.t() * Z_u.t() + Z_v * P3 * Z_v.t();
      
      F1 = inv(A);
      AmB = A + mt * B;
      F2 = - F1 * B * inv(AmB);
      F1 = (F1 + F1.t()) / 2.0;
      F2 = (F2 + F2.t()) / 2.0;
      F1mF2 = F1 + mt * F2;
      
      sum_v_a = sum(v_a.cols(A_vec), 1);
      
      C0 = (P0 * Z_u.t() + P1 * Z_v.t()) * F1mF2;
      C1 = P2 * Z_v.t() * F1;
      C2 = (P1.t() * Z_u.t() + P3 * Z_v.t()) * F1mF2 + P2 * Z_v.t() * F2;
      
      a_u_a_tt = a_u_a + C0 * sum_v_a;
      a_v_a_tt = a_v_a + C1 * v_a + kron(mat_1m, C2 * sum_v_a);
      
      for (int k = 0; k < d_u; ++k){
        sum_V_A = sum(V_A.slice(k).cols(A_vec), 1);
        A_u_A_tt.col(k) = A_u_A.col(k) + C0 * sum_V_A;
        A_v_A_tt.slice(k) = A_v_A.slice(k) + C1 * V_A.slice(k) + kron(mat_1m, C2 * sum_V_A);
      }
      
      M0 = mt * (P0 * Z_u.t() + P1 * Z_v.t()) * F1mF2 * (Z_u * P0 + Z_v * P1.t());
      M1 = (P0 * Z_u.t() + P1 * Z_v.t()) * F1mF2 * (mt * Z_u * P1 + Z_v * (P2 + mt * P3));
      M2 = P2 * Z_v.t() * F1 * Z_v * P2;
      M3 = P2 * Z_v.t() * F2 * Z_v * P2 + 
        P2 * Z_v.t() * F1mF2 * (Z_u * P1 + Z_v * P3) + 
        (P1.t() * Z_u.t() + P3 * Z_v.t()) * F1mF2 * Z_v * P2 +
        mt * (P1.t() * Z_u.t() + P3 * Z_v.t()) * F1mF2 * (Z_u * P1 + Z_v * P3);
      
      P0_tt = P0 - M0;
      P1_tt = P1 - M1;
      P2_tt = P2 - M2;
      P3_tt = P3 - M3;
      
      // Delta terms
      sum_V_d = mat(q, d_u, fill::zeros);
      sum_VFV = mat(d_u, d_u, fill::zeros);
      sum_VFv = vec(d_u, fill::zeros);
      for (int i = 0; i < mt; ++i){
        V_temp = mat(V_A.subcube(0, A_vec[i], 0, q - 1, A_vec[i], d_u - 1));
        sum_V_d += V_temp;
        sum_VFV += V_temp.t() * F1 * V_temp;
        sum_VFv += V_temp.t() * F1 * v_a.col(A_vec[i]);
      }
      
      delta_term1 += sum_VFV + sum_V_d.t() * F2 * sum_V_d;
      delta_term2 += sum_VFv + sum_V_d.t() * F2 * sum_v_a;
      
      a_u_a = T_u * a_u_a_tt;
      a_v_a = T_v * a_v_a_tt;  // T_v = I, so this is just a_v_a_tt
      
      for (int k = 0; k < d_u; ++k){
        A_u_A.col(k) = T_u * A_u_A_tt.col(k);
        A_v_A.slice(k) = T_v * A_v_A_tt.slice(k);  // T_v = I
      }
      
      P0 = T_u * P0_tt * T_u.t() + Sig_u;
      P1 = T_u * P1_tt * T_v.t();
      P2 = T_v * P2_tt * T_v.t() + Sig_v;  // T_v = I, so P2_tt + Sig_v
      P3 = T_v * P3_tt * T_v.t();  // T_v = I, so just P3_tt
      
      P0 = (P0 + P0.t()) / 2.0;
      P2 = (P2 + P2.t()) / 2.0;
      P3 = (P3 + P3.t()) / 2.0;
    }
    
    delta_hat = - inv(delta_term1) * delta_term2;
    return delta_hat;
  }

// Implementation of function neg2llh_cs_local ----------------------------------------------------------------------
  double neg2llh_cs_local(vec& theta_trans, cube& data_split_cube, mat& z_mat, mat& D_mat, cube& X00_cube, List& A_list, List& AB_list){
    
    vec theta, gamma_vec, gamma_x, delta, zeta_vec, lambda_vec, d_vec, l_vec, sum_v, Z_nonNA, p_z_sim, cond_mean, sum_diff,
    theta_t_plus_1, z_term1, z_term2, z_term, a_u, a_u_tt, a_u_tt_nonNA, alpha_u_gen;
    
    mat data_split_j, D, L, tmp, tmp_lower, Z00, T00, P0, P1, P2, P3, Z_u, Z_v, T_u, T_v,
    Sigma, Lambda, Sig_u, Sig_v, P2_tt_inv, temp, var_mat, a_v, a_v_tt, a_v_tt_nonNA, alpha_v_gen,
    F1mF2, F1, F2, C0, C1, C2, M0, M1, M2, M3, 
    P0_tt, P1_tt, P2_tt, P3_tt, cond_chol, mat_1mt_next, stacked_mat, stacked_mat2,
    Za, v_split, mat_1m, S_u, chol_terms_u, A, B, AmB;
    
    int q, p, npt, m, mt, mt_next, d_u, d_v, from, to, nsim = 100, dim_gamma;
    
    double dt, llh, llh_const, llh_term1, llh_term2, p_z, p_hat, kappa;
    
    uvec idx, A_vec, AB_vec_next;
    
    cube S_v, chol_terms_v;
    
    // Estimate delta
    delta = delta_est_cs_local(theta_trans, data_split_cube, A_list);
    
    // Initialize
    q = data_split_cube.n_rows;
    npt = data_split_cube.n_cols;
    m = data_split_cube.n_slices;
    p = X00_cube.n_cols;
    d_u = 2 * q;
    d_v = q;  // LOCAL LEVEL
    mat_1m = mat(1, m, fill::ones);
    dim_gamma = D_mat.n_rows;
    
    // Transform parameters
    theta = theta_trans;
    from = 0;
    to = 3 * q - 1;
    theta.subvec(from, to) = exp(theta_trans.subvec(from, to));
    from = 0.5 * q * q + 2.5 * q;
    theta(from) = exp(theta_trans(from));
    
    // Extract parameters
    from = 0; to = q - 1;
    zeta_vec = theta.subvec(from, to);
    
    from = q; to = 2 * q - 1;
    lambda_vec = theta.subvec(from, to);
    
    from = 2 * q; to = 3 * q - 1;
    d_vec = theta.subvec(from, to);
    D = diagmat(d_vec);
    
    if (q > 1){
      from = 3 * q;
      to = 0.5 * q * q + 2.5 * q - 1;
      l_vec = theta.subvec(from, to);
      L = eye(q, q);
      tmp = ones(size(L));
      tmp_lower = trimatl(tmp);
      tmp_lower.diag().zeros();
      idx = find(tmp_lower > 0);
      L.elem(idx) = l_vec;
    } else L = eye(1, 1);
    
    from = 0.5 * q * q + 2.5 * q;
    kappa = theta(from);
    
    from = 0.5 * q * q + 2.5 * q + 1;
    to = 0.5 * q * q + 2.5 * q + dim_gamma;
    gamma_vec = theta.subvec(from, to);
    
    from = 0.5 * q * q + 2.5 * q + dim_gamma + 1;
    to = 0.5 * q * q + 2.5 * q + dim_gamma + p + 1;
    gamma_x = theta.subvec(from, to);
    
    // System matrices
    Z00 = zeros(1, 2);
    Z00(0, 0) = 1;
    Z_u = kron(eye(q, q), Z00);
    Z_v = eye(q, q);  // LOCAL LEVEL
    
    dt = 1.0;
    T00 = eye(2, 2);
    T00(0, 1) = dt;
    T_u = kron(eye(q, q), T00);
    T_v = eye(q, q);  // LOCAL LEVEL: random walk
    
    Lambda = zeros(2, 2);
    Lambda(0, 0) = dt * dt * dt / 3;
    Lambda(0, 1) = Lambda(1, 0) = dt * dt / 2;
    Lambda(1, 1) = dt;
    Sig_u = kron(diagmat(zeta_vec), Lambda);
    Sig_v = diagmat(lambda_vec);  // LOCAL LEVEL
    
    // Initialize
    P0 = zeros(d_u, d_u);
    P1 = zeros(d_u, d_v);
    P3 = zeros(d_v, d_v);
    P2 = kappa * eye(d_v, d_v);
    
    a_u = delta;
    a_v = zeros(d_v, m);
    
    Sigma = L * D * L.t();
    llh = 0.0;
    
    // FILTERING LOOP
    for (int j = 0; j < npt; ++j){
      A_vec = as<uvec>(A_list[j]);
      A_vec -= 1;
      mt = A_vec.n_elem;
      
      data_split_j = data_split_cube.subcube(0, j, 0, q - 1, j, m - 1);
      Za = kron(mat_1m, Z_u * a_u) + Z_v * a_v;
      v_split = data_split_j - Za;
      
      // F components
      A = Z_v * P2 * Z_v.t() + Sigma;
      B = Z_u * P0 * Z_u.t() + Z_u * P1 * Z_v.t() + Z_v * P1.t() * Z_u.t() + Z_v * P3 * Z_v.t();
      
      F1 = inv(A);
      AmB = A + mt * B;
      F2 = - F1 * B * inv(AmB);
      F1 = (F1 + F1.t()) / 2.0;
      F2 = (F2 + F2.t()) / 2.0;
      F1mF2 = F1 + mt * F2;
      
      sum_v = sum(v_split.cols(A_vec), 1);
      
      // Update
      C0 = (P0 * Z_u.t() + P1 * Z_v.t()) * F1mF2;
      C1 = P2 * Z_v.t() * F1;
      C2 = (P1.t() * Z_u.t() + P3 * Z_v.t()) * F1mF2 + P2 * Z_v.t() * F2;
      
      a_u_tt = a_u + C0 * sum_v;
      a_v_tt = a_v + C1 * v_split + kron(mat_1m, C2 * sum_v);
      
      // P update
      M0 = mt * (P0 * Z_u.t() + P1 * Z_v.t()) * F1mF2 * (Z_u * P0 + Z_v * P1.t());
      M1 = (P0 * Z_u.t() + P1 * Z_v.t()) * F1mF2 * (mt * Z_u * P1 + Z_v * (P2 + mt * P3));
      M2 = P2 * Z_v.t() * F1 * Z_v * P2;
      M3 = P2 * Z_v.t() * F2 * Z_v * P2 + 
        P2 * Z_v.t() * F1mF2 * (Z_u * P1 + Z_v * P3) + 
        (P1.t() * Z_u.t() + P3 * Z_v.t()) * F1mF2 * Z_v * P2 +
        mt * (P1.t() * Z_u.t() + P3 * Z_v.t()) * F1mF2 * (Z_u * P1 + Z_v * P3);
      
      P0_tt = P0 - M0;
      P1_tt = P1 - M1;
      P2_tt = P2 - M2;
      P3_tt = P3 - M3;
      
      // Longitudinal likelihood
      llh_term1 = 0.0;
      for (int k = 0; k < mt; ++k){
        llh_term1 += as_scalar(v_split.col(A_vec(k)).t() * F1 * v_split.col(A_vec(k)));
      }
      
      llh_term2 = as_scalar(sum_v.t() * F2 * sum_v);
      llh_const = - (q * mt) / 2.0 * log(2.0 * datum::pi);
      
      if (det(F1mF2) > 0 & det(F1) > 0)
        llh += llh_const + 0.5 * (log(det(F1mF2)) + (mt - 1.0) * log(det(F1))) - 0.5 * (llh_term1 + llh_term2);
      
      // SURVIVAL LIKELIHOOD (same structure as OU)
      if ((j < npt - 1) & (j > 0)){
        AB_vec_next = as<uvec>(AB_list[j + 1]);
        AB_vec_next -= 1;
        mt_next = AB_vec_next.n_elem;
        if (mt_next <= 0) break;
        
        S_v = zeros(d_v, d_v, mt_next);
        chol_terms_v = zeros(d_v, d_v, mt_next);
        
        P2_tt_inv = inv(P2_tt);
        P2_tt_inv = (P2_tt_inv + P2_tt_inv.t()) / 2.0;
        
        temp = P2_tt + mt_next * P3_tt;
        S_u = P1_tt * (P2_tt_inv - mt_next * P2_tt_inv * P3_tt * inv(temp));
        var_mat = P0_tt - mt_next * S_u * P1_tt.t();
        var_mat = (var_mat + var_mat.t()) / 2.0;
        chol_terms_u = my_chol(var_mat);
        
        if (mt_next > 2){
          for (int k = 1; k <= mt_next - 2; ++k){
            temp = P2_tt + (mt_next - k) * P3_tt;
            S_v.slice(k - 1) = P3_tt * (P2_tt_inv - (mt_next - k) * P2_tt_inv * P3_tt * inv(temp));
            var_mat = (P2_tt + P3_tt) - (mt_next - k) * S_v.slice(k - 1) * P3_tt;
            var_mat = (var_mat + var_mat.t()) / 2.0;
            chol_terms_v.slice(k - 1) = my_chol(var_mat);
          }
        }
        
        if (mt_next > 1){
          temp = P2_tt + P3_tt;
          S_v.slice(mt_next - 2) = P3_tt * inv(temp);
          var_mat = (P2_tt + P3_tt) - S_v.slice(mt_next - 2) * P3_tt;
          var_mat = (var_mat + var_mat.t()) / 2.0;
          chol_terms_v.slice(mt_next - 2) = my_chol(var_mat);
        }
        
        var_mat = P2_tt + P3_tt;
        var_mat = (var_mat + var_mat.t()) / 2.0;
        chol_terms_v.slice(mt_next - 1) = my_chol(var_mat);
        
        a_u_tt_nonNA = a_u_tt;
        a_v_tt_nonNA = a_v_tt.cols(AB_vec_next);
        
        Z_nonNA = vec(mt_next);
        for (int i = 0; i < mt_next; ++i){
          Z_nonNA(i) = as_scalar(z_mat(AB_vec_next(i), j + 1));
        }
        
        alpha_v_gen = mat(d_v, mt_next);
        p_z_sim = vec(2 * nsim, fill::zeros);
        
        for (int k = 0; k < nsim; ++k){
          cond_mean = a_v_tt_nonNA.col(mt_next - 1);
          cond_chol = chol_terms_v.slice(mt_next - 1);
          alpha_v_gen.col(mt_next - 1) = cond_mean + cond_chol * randn(d_v);
          
          sum_diff = alpha_v_gen.col(mt_next - 1) - a_v_tt_nonNA.col(mt_next - 1);
          for (int i = mt_next - 1; i >= 1; --i){
            sum_diff += alpha_v_gen.col(i) - a_v_tt_nonNA.col(i);
            cond_mean = a_v_tt_nonNA.col(i - 1) + S_v.slice(i - 1) * sum_diff;
            cond_chol = chol_terms_v.slice(i - 1);
            alpha_v_gen.col(i - 1) = cond_mean + cond_chol * randn(d_v);
          }
          
          sum_diff += alpha_v_gen.col(0) - a_v_tt_nonNA.col(0);
          cond_mean = a_u_tt_nonNA + S_u * sum_diff;
          cond_chol = chol_terms_u;
          alpha_u_gen = cond_mean + cond_chol * randn(d_u);
          
          mat_1mt_next = mat(1, mt_next, fill::ones);
          stacked_mat = zeros(d_u + d_v, mt_next);
          stacked_mat.submat(0, 0, d_u - 1, mt_next - 1) = kron(mat_1mt_next, alpha_u_gen);
          stacked_mat.submat(d_u, 0, d_u + d_v - 1, mt_next - 1) = alpha_v_gen;
          
          from = 1; to = p;
          theta_t_plus_1 = stacked_mat.t() * D_mat.t() * gamma_vec + 
            X00_cube.slice(j).rows(AB_vec_next) * gamma_x.subvec(from, to) + gamma_x(0);
          
          z_term1 = (Z_nonNA - 1) % theta_t_plus_1;
          z_term2 = log(1 + exp(- theta_t_plus_1));
          z_term = z_term1 - z_term2;
          p_z = exp(sum(z_term));
          p_z_sim(2 * k) = p_z;
          
          stacked_mat2 = zeros(d_u + d_v, mt_next);
          stacked_mat2.submat(0, 0, d_u - 1, mt_next - 1) = kron(mat_1mt_next, 2 * a_u_tt_nonNA - alpha_u_gen);
          stacked_mat2.submat(d_u, 0, d_u + d_v - 1, mt_next - 1) = 2 * a_v_tt_nonNA - alpha_v_gen;
          
          theta_t_plus_1 = stacked_mat2.t() * D_mat.t() * gamma_vec + 
            X00_cube.slice(j).rows(AB_vec_next) * gamma_x.subvec(from, to) + gamma_x(0);
          
          z_term1 = (Z_nonNA - 1) % theta_t_plus_1;
          z_term2 = log(1 + exp(- theta_t_plus_1));
          z_term = z_term1 - z_term2;
          p_z = exp(sum(z_term));
          p_z_sim(2 * k + 1) = p_z;
        }
        
        p_hat = mean(p_z_sim);
        if (p_hat > 0) llh += log(p_hat);
      }
      
      // Predict next state
      a_u = T_u * a_u_tt;
      a_v = T_v * a_v_tt;  // T_v = I, so just a_v_tt
      
      P0 = T_u * P0_tt * T_u.t() + Sig_u;
      P1 = T_u * P1_tt * T_v.t();
      P2 = T_v * P2_tt * T_v.t() + Sig_v;  // T_v = I, so P2_tt + Sig_v
      P3 = T_v * P3_tt * T_v.t();  // T_v = I, so P3_tt
      
      P0 = (P0 + P0.t()) / 2.0;
      P2 = (P2 + P2.t()) / 2.0;
      P3 = (P3 + P3.t()) / 2.0;
    }
    
    return -2.0 * llh;
  }
'

src = '
vec theta_trans;
mat z_mat, D_mat;
cube data_split_cube, X00_cube;
List A_list, AB_list;

theta_trans = as<vec>(theta_trans_in);
z_mat = as<mat>(z_in);
D_mat = as<mat>(D_in);
data_split_cube = as<cube>(data_split_in);
X00_cube = as<cube>(X00_cube_in);
A_list = as<List>(A_list_in);
AB_list = as<List>(AB_list_in);

return wrap(neg2llh_cs_local(theta_trans, data_split_cube, z_mat, D_mat, X00_cube, A_list, AB_list));
'

llh_survival_x_cs_local = cxxfunction(
  signature(theta_trans_in = "numeric", data_split_in = "numeric", z_in = "numeric", 
            D_in = "numeric", X00_cube_in = "numeric", A_list_in = "numeric", AB_list_in = "numeric"),
  includes = rcpp_inc_cs_local_joint, body = src, plugin = "RcppArmadillo"
)

# Parameter structure for cubic-local:
# theta_trans = (log(zeta_vec[q]), log(lambda_vec[q]), log(d_vec[q]), l_vec[q*(q-1)/2], log(kappa), gamma_vec[q], gamma_x[p+1])
# Total length: q + q + q + q*(q-1)/2 + 1 + q + (p+1) = 0.5*q^2 + 3.5*q + p + 2
