# =============================================================================
# Joint Log-Likelihood: Cubic Spline Population + OU Individual + Survival
# =============================================================================
# Rcpp/inline implementation of the full joint log-likelihood (Section 3.4).
# Combines the longitudinal OU model with a discrete-time logistic survival model.
#
# Returns: -2 * joint log-likelihood (scalar, to be minimized by optim)
#
# Full parameter vector theta_trans (all log-transformed except l_vec, gamma):
#   [1:q]            log(zeta_vec)   population smoothing variances
#   [q+1:2q]         log(nu2_vec)    OU innovation variances
#   [2q+1:3q]        log(xi_vec)     OU mean-reversion rates
#   [3q+1:4q]        log(d_vec)      Cholesky diagonal of cross-biomarker Sigma
#   [4q+1:n_long-1]  l_vec           Cholesky lower triangle of Sigma
#   [n_long]         log(kappa)      survival baseline hazard scale
#   [n_long+1:       gamma_vec       q-vector linking individual OU state to survival log-odds
#    n_long+q]
#   [n_long+q+1:     gamma_x         (p+1)-vector: intercept + p survival covariate coefficients
#    n_long+q+p+1]
#   where n_long = 0.5*q^2 + 3.5*q + 1
#
# Inputs (R-level):
#   theta_trans_in:  numeric vector (length n_long + q + p + 1)
#   data_split_in:   q x npt x m array of standardized biomarker observations
#   z_in:            m x npt matrix; z_in[i,t] = 1 if subject i died at month t
#   D_in:            q x (2q+q) matrix mapping state vector to observations
#                    For OU: D_in = [I_q x (1,0), I_q] selecting position components
#   X00_cube_in:     m x p x npt covariate array
#   A_list_in:       list of length npt; alive-at-t subject indices
#   AB_list_in:      list of length npt; present-at-t subject indices (alive or dying)
# =============================================================================

library(Rcpp)
library(inline)
rcpp_inc = '
#include<iostream>
using namespace Rcpp;
using namespace arma;
using namespace std;

mat my_inv(mat& Q);
// inverse using eigendecomposition. If abs(eigenvalue) < 10^(-16), set eigenvalue = 10^(-8)

mat my_chol(mat& Q);
// eigen decomposition of Q, setting negative eigenvalues to 0

vec delta_est(vec& theta_trans, cube& data_split_cube, List& A_list);
// theta_trans = (log(zeta_vec, lambda_vec, d_vec), l_vec, log(kappa))
// returns a vector, delta_hat

double neg2llh(vec& theta_trans, cube& data_split_cube, mat& z_mat, mat& D_mat, cube& X00_cube, List& A_list, List& AB_list);
// returns -2 * llh

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


// implementation of function delta_est ---------------------------------------------------------------------------
vec delta_est(vec& theta_trans, cube& data_split_cube, List& A_list){
// -----------------------------------------
// Declare variables
// -----------------------------------------

vec theta, zeta_vec, nu2_vec, xi_vec, d_vec, l_vec, sum_v_a, delta_term2, delta_hat, sum_V_A,
a_u_a, a_u_a_tt;


mat y, D, L, L_inv, tmp, tmp_lower, Z00, T00, P0, P1, P2, P3, Z_u, Z_v, T_u, T_v, Q00,
Sigma, Lambda, Sig_u, Sig_v,
F1mF2, F1, F2, K0, K1, K2, 
A, B, AmB, C0, C1, C2, M0, M1, M2, M3, 
a_a_tt, P0_tt, P1_tt, P2_tt, P3_tt,
Za_a, ZA_A, a_a, v_a, a_v_a, a_v_a_tt, A_u_A, A_u_A_tt,
mat_1m,  
delta_term1, sum_V_d, sum_VFV, sum_VFv, V_temp;

cube A_A, A_A_tt, V_A, A_v_A, A_v_A_tt;  

double dt;

int q, npt, m, mt, d_u, d_v, from, to;

uvec idx, A_vec; // for creating L from l_vec

// -----------------------------------------
// Initializing necessary quantities
// -----------------------------------------
// q, npt, m, p, mat_1m
q = data_split_cube.n_rows;      // data_split_cube is q x npt x m
npt = data_split_cube.n_cols;
m = data_split_cube.n_slices;
d_u = 2 * q; 
d_v = q;
mat_1m = mat(1, m, fill::ones);


// transform parameters to original scale
theta = theta_trans;
from = 0;
to = 4 * q - 1;
theta.subvec(from, to) = exp(theta_trans.subvec(from, to));

// smoothing parameters, population level
from = 0;
to = q - 1;
zeta_vec = theta.subvec(from, to);

// nu2_vec
from = q;
to = 2 * q - 1;
nu2_vec = theta.subvec(from, to);

// xi_vec
from = 2 * q;
to = 3 * q - 1;
xi_vec = theta.subvec(from, to);

// diagonal elements of D, where Sigma = LDt(L)
from = 3 * q;
to = 4 * q - 1;
d_vec = theta.subvec(from, to);
D = diagmat(d_vec);

// L
if (q > 1){
from = 4 * q;
to = 0.5 * q * q + 3.5 * q - 1;
l_vec = theta.subvec(from, to);
L = eye(q, q);
// Find the indices for the lower triangular part of L
tmp = ones(size(L));
tmp_lower = trimatl(tmp);
tmp_lower.diag().zeros();
idx = find(tmp_lower > 0);
L.elem(idx) = l_vec;
} else L = eye(1, 1);


// -----------------------------------
// system matrices
// -----------------------------------
// components of Z
Z00 = zeros(1, 2);
Z00(0, 0) = 1;
Z_u = kron(eye(q, q), Z00);
Z_v = eye(q, q);

// components of T
dt = 1.0;
T00 = eye(2, 2);
T00(0, 1) = dt;
T_u = kron(eye(q, q), T00);
T_v = diagmat(exp(- xi_vec * dt));

// components of Q
Lambda = zeros(2, 2);
Lambda(0, 0) = dt * dt * dt / 3;
Lambda(0, 1) = Lambda(1, 0) = dt * dt / 2;
Lambda(1, 1) = dt;
Sig_u = kron(diagmat(zeta_vec), Lambda);
Sig_v = diagmat(0.5 * nu2_vec % (1 - exp(-2 * xi_vec * dt)) / xi_vec );
// the diagonal block of Q0
Q00 = diagmat(0.5 * nu2_vec / xi_vec);  

// Initialize P(t_1) and a(t_1)
P0 = zeros(d_u, d_u);
P1 = zeros(d_u, d_v);
P3 = zeros(d_v, d_v);
P2 = Q00;

// initialize a_{a, t}, A_{A, t}, A_{A, tt}, V_{A, t} 
// initial state: fixed effects and random effects components
a_u_a = zeros(d_u, 1);
a_v_a = zeros(d_v, m);
// augmented part
A_u_A = eye(d_u, d_u);
A_v_A = zeros(d_v, m, d_u);

A_u_A_tt = zeros(d_u, d_u);
A_v_A_tt = zeros(d_v, m, d_u);

V_A = zeros(q, m, d_u);

// Sigma
Sigma = L * D * L.t();

// initialize delta_term1 and delta_term2
delta_term1 = mat(d_u, d_u, fill::zeros);
delta_term2 = vec(d_u, fill::zeros);


// -------------------------------------------------------------------------
// Filtering
// -------------------------------------------------------------------------

// loop through all the time points
for (int j = 0; j < npt; ++j){
A_vec = as<uvec>(A_list[j]); // IDs start from 1
A_vec -= 1; // make IDs start from 0
mt = A_vec.n_elem;

// filtering
// data of all subjects at time point j
y = data_split_cube.subcube(0, j, 0, q - 1, j, m - 1); // q x m, originally qm x 1

// Z_t %*% a_{a, t}
Za_a = kron(mat_1m, Z_u * a_u_a) + Z_v * a_v_a;
// v_{a, t}
v_a = y - Za_a; // v_a = (v_{a, 1}, ..., v_{a, m}) at time point j, with dimension q x m
// augmented part
for (int k = 0; k < d_u; ++k){
// Z_t * A_{A, t, k}, where A_{A, t, k} is the kth column of A_{A, t}, split into a d_u x 1 vector (fixed effects) and a d_v x m matrix (random effects)
ZA_A = kron(mat_1m, Z_u * A_u_A.col(k)) + Z_v * A_v_A.slice(k);
V_A.slice(k) = - ZA_A;
}

// Covariance part of KF
// components of F
A = Z_v * P2 * Z_v.t() + Sigma;
B = Z_u * P0 * Z_u.t() + Z_u * P1 * Z_v.t() + Z_v * P1.t() * Z_u.t() + Z_v * P3 * Z_v.t();

// Components of F_inv
F1 = inv(A);
AmB = A + mt * B;
F2 = - F1 * B * inv(AmB);
// Symmetrize F1 and F2
F1 = (F1 + F1.t()) / 2.0;
F2 = (F2 + F2.t()) / 2.0;
F1mF2 = F1 + mt * F2;  

// v_1 + ... + v_m, sum of nonNA v_i
sum_v_a = sum(v_a.cols(A_vec), 1);


// a_{t|t}
C0 = (P0 * Z_u.t() + P1 * Z_v.t()) * F1mF2;
C1 = P2 * Z_v.t() * F1;
C2 = (P1.t() * Z_u.t() + P3 * Z_v.t()) * F1mF2 + P2 * Z_v.t() * F2;

// v_{a, 1} + ... + v_{a, m}, sum of nonNA v_i
sum_v_a = sum(v_a.cols(A_vec), 1);
// a_{a, tt}
a_u_a_tt = a_u_a;
a_v_a_tt = a_v_a;
a_u_a_tt += C0 * sum_v_a;
a_v_a_tt += C1 * v_a + kron(mat_1m, C2 * sum_v_a);



// augmented part
for (int k = 0; k < d_u; ++k){
sum_V_A = sum(V_A.slice(k).cols(A_vec), 1);
A_u_A_tt.col(k) = A_u_A.col(k);
A_v_A_tt.slice(k) = A_v_A.slice(k);
A_u_A_tt.col(k) += C0 * sum_V_A;
A_v_A_tt.slice(k) += C1 * V_A.slice(k) + kron(mat_1m, C2 * sum_V_A);
}




// P_{t|t}
M0 = mt * (P0 * Z_u.t() + P1 * Z_v.t()) * F1mF2 * (Z_u * P0 + Z_v * P1.t());
M1 = (P0 * Z_u.t() + P1 * Z_v.t()) * F1mF2 * (mt * Z_u * P1 + Z_v * (P2 + mt * P3));
M2 = P2 * Z_v.t() * F1 * Z_v * P2;
M3 = P2 * Z_v.t() * F2 * Z_v * P2 + 
P2 * Z_v.t() * F1mF2 * (Z_u * P1 + Z_v * P3) + (P1.t() * Z_u.t() + P3 * Z_v.t()) * F1mF2 * Z_v * P2 +
mt * (P1.t() * Z_u.t() + P3 * Z_v.t()) * F1mF2 * (Z_u * P1 + Z_v * P3);


P0_tt = P0 - M0;
P1_tt = P1 - M1;
P2_tt = P2 - M2;
P3_tt = P3 - M3;



// delta_terms ------------------------------------------------------------------------------
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


// a_{a, t + 1}
a_u_a = T_u * a_u_a_tt;
a_v_a = T_v * a_v_a_tt;
// augmented part, A_{A, t + 1}
for (int k = 0; k < d_u; ++k){
A_u_A.col(k) = T_u * A_u_A_tt.col(k);
A_v_A.slice(k) = T_v * A_v_A_tt.slice(k);
}


// P_{t + 1}
P0 = T_u * P0_tt * T_u.t() + Sig_u;
P1 = T_u * P1_tt * T_v.t();
P2 = T_v * P2_tt * T_v.t() + Sig_v;
P3 = T_v * P3_tt * T_v.t();

P0 = (P0 + P0.t()) / 2.0;
P2 = (P2 + P2.t()) / 2.0;
P3 = (P3 + P3.t()) / 2.0;
}

//cout << delta_term1 << endl;
delta_hat = - inv(delta_term1) * delta_term2;
return delta_hat;
}


// Implementation of function neg2llh ----------------------------------------------------------------------
double neg2llh(vec& theta_trans, cube& data_split_cube, mat& z_mat, mat& D_mat, cube& X00_cube, List& A_list, List& AB_list){

vec theta, gamma_vec, gamma_x, delta, zeta_vec, nu2_vec, xi_vec, d_vec, l_vec, sum_v, Z_nonNA, p_z_sim, cond_mean, sum_diff,
theta_t_plus_1, z_term1, z_term2, z_term, a_u, a_u_tt, a_u_tt_nonNA, alpha_u_gen;

mat data_split_j, D, L, L_inv, tmp, tmp_lower, Z00, T00, P0, P1, P2, P3, Z_u, Z_v, T_u, T_v,
Sigma, Lambda, Q00, Sig_u, Sig_v, P2_tt_inv, temp, var_mat, a_v, a_v_tt, a_v_tt_nonNA, alpha_v_gen,
F1mF2, F1, F2, K0, K1, K2, 
A, B, AmB, C0, C1, C2, M0, M1, M2, M3, 
P0_tt, P1_tt, P2_tt, P3_tt,
cond_chol, mat_1mt_next, stacked_mat, stacked_mat2,
Za, v_split, 
mat_1m, S_u, chol_terms_u;

int q, p, npt, m, mt, mt_next, d_u, d_v, from, to, nsim = 100, dim_gamma;

double dt, llh, llh_const, llh_term1, llh_term2, p_z, p_hat; 

uvec idx, A_vec, AB_vec_next; // for creating L from l_vec

cube S_v, chol_terms_v;


// -----------------------------------------
// Convert input quantities to C++ types
// -----------------------------------------

// estimate delta using delta_est function
delta = delta_est(theta_trans, data_split_cube, A_list);

// -----------------------------------------
// Initializing necessary quantities
// -----------------------------------------
// q, npt, m, d, mat_1m
q = data_split_cube.n_rows;        // number of longitudinal variables
npt = data_split_cube.n_cols;
m = data_split_cube.n_slices;
p = X00_cube.n_cols;               // X00_cube is m x p x npt
d_u = 2 * q; 
d_v = q;
mat_1m = mat(1, m, fill::ones);  
dim_gamma = D_mat.n_rows;

// transform parameters to original scale
theta = theta_trans;
from = 0;
to = 4 * q - 1;
theta.subvec(from, to) = exp(theta_trans.subvec(from, to));

// smoothing parameters, population level
from = 0;
to = q - 1;
zeta_vec = theta.subvec(from, to);

// nu2_vec
from = q;
to = 2 * q - 1;
nu2_vec = theta.subvec(from, to);

// xi_vec
from = 2 * q;
to = 3 * q - 1;
xi_vec = theta.subvec(from, to);

// diagonal elements of D, where Sigma = LDt(L)
from = 3 * q;
to = 4 * q - 1;
d_vec = theta.subvec(from, to);
D = diagmat(d_vec);

// L
if (q > 1){
from = 4 * q;
to = 0.5 * q * q + 3.5 * q - 1;
l_vec = theta.subvec(from, to);
L = eye(q, q);
// Find the indices for the lower triangular part of L
tmp = ones(size(L));
tmp_lower = trimatl(tmp);
tmp_lower.diag().zeros();
idx = find(tmp_lower > 0);
L.elem(idx) = l_vec;
} else L = eye(1, 1);

// gamma_vec
from = 0.5 * q * q + 3.5 * q;
to = 0.5 * q * q + 3.5 * q + dim_gamma - 1;
gamma_vec = theta.subvec(from, to);


// gamma_x = (intercept, for x), length p + 1
from = 0.5 * q * q + 3.5 * q + dim_gamma;
to = 0.5 * q * q + 3.5 * q + dim_gamma + p;
gamma_x = theta.subvec(from, to);


// -----------------------------------
// system matrices
// -----------------------------------
// components of Z
Z00 = zeros(1, 2);
Z00(0, 0) = 1;
Z_u = kron(eye(q, q), Z00);
Z_v = eye(q, q);

// components of T
dt = 1.0;
T00 = eye(2, 2);
T00(0, 1) = dt;
T_u = kron(eye(q, q), T00);
T_v = diagmat(exp(- xi_vec * dt));

// components of Q
Lambda = zeros(2, 2);
Lambda(0, 0) = dt * dt * dt / 3;
Lambda(0, 1) = Lambda(1, 0) = dt * dt / 2;
Lambda(1, 1) = dt;
Sig_u = kron(diagmat(zeta_vec), Lambda);
Sig_v = diagmat(0.5 * nu2_vec % (1 - exp(-2 * xi_vec * dt)) / xi_vec );
// the diagonal block of Q0
Q00 = diagmat(0.5 * nu2_vec / xi_vec);  

// Initialize P(t_1) and a(t_1)
P0 = zeros(d_u, d_u);
P1 = zeros(d_u, d_v);
P3 = zeros(d_v, d_v);
P2 = Q00;

// initialize a1
a_u = delta;
a_v = zeros(d_v, m);

// Sigma
Sigma = L * D * L.t();


// -------------------------------------------------------------------------
// Filtering. And storing matrices for use in the smoothing step
// -------------------------------------------------------------------------
// initialize log-likelihood
llh = 0.0;

// loop through all the time points
for (int j = 0; j < npt; ++j){
A_vec = as<uvec>(A_list[j]); // IDs start from 1
A_vec -= 1; // make IDs start from 0
mt = A_vec.n_elem;

// filtering
// data of all subjects at time point j
data_split_j = data_split_cube.subcube(0, j, 0, q - 1, j, m - 1);

// Z_t %*% a_t
Za = kron(mat_1m, Z_u * a_u) + Z_v * a_v;
v_split = data_split_j - Za; // v_split = (v_1, ..., v_m) at time point j, with dimension q x m


// components of F
A = Z_v * P2 * Z_v.t() + Sigma;
B = Z_u * P0 * Z_u.t() + Z_u * P1 * Z_v.t() + Z_v * P1.t() * Z_u.t() + Z_v * P3 * Z_v.t();

// Components of F_inv

F1 = inv(A);

AmB = A + mt * B;
F2 = - F1 * B * inv(AmB);
// Symmetrize F1 and F2
F1 = (F1 + F1.t()) / 2.0;
F2 = (F2 + F2.t()) / 2.0;
F1mF2 = F1 + mt * F2; 

// v_1 + ... + v_m, sum of nonNA v_i
sum_v = sum(v_split.cols(A_vec), 1);


// a_{t|t}
C0 = (P0 * Z_u.t() + P1 * Z_v.t()) * F1mF2;
C1 = P2 * Z_v.t() * F1;
C2 = (P1.t() * Z_u.t() + P3 * Z_v.t()) * F1mF2 + P2 * Z_v.t() * F2;
a_u_tt = a_u;
a_v_tt = a_v;
a_u_tt += C0 * sum_v;
a_v_tt += C1 * v_split + kron(mat_1m, C2 * sum_v);


// P_{t|t}
M0 = mt * (P0 * Z_u.t() + P1 * Z_v.t()) * F1mF2 * (Z_u * P0 + Z_v * P1.t());
M1 = (P0 * Z_u.t() + P1 * Z_v.t()) * F1mF2 * (mt * Z_u * P1 + Z_v * (P2 + mt * P3));
M2 = P2 * Z_v.t() * F1 * Z_v * P2;
M3 = P2 * Z_v.t() * F2 * Z_v * P2 + 
P2 * Z_v.t() * F1mF2 * (Z_u * P1 + Z_v * P3) + (P1.t() * Z_u.t() + P3 * Z_v.t()) * F1mF2 * Z_v * P2 +
mt * (P1.t() * Z_u.t() + P3 * Z_v.t()) * F1mF2 * (Z_u * P1 + Z_v * P3);


P0_tt = P0 - M0;
P1_tt = P1 - M1;
P2_tt = P2 - M2;
P3_tt = P3 - M3;


// log-likelihood ------------------------------------------------------------------------------
// llh_term1
llh_term1 = 0.0;
for (int k = 0; k < mt; ++k){
llh_term1 += as_scalar(v_split.col(A_vec(k)).t() * F1 * v_split.col(A_vec(k)));
}


// llh_term2
llh_term2 = as_scalar(sum_v.t() * F2 * sum_v);
// constant part of llh
llh_const = - (q * mt) / 2.0 * log(2.0 * datum::pi);
if (det(F1mF2) > 0 & det(F1) > 0)
llh += llh_const + 0.5 * (log(det(F1mF2)) + (mt - 1.0) * log(det(F1))) - 0.5 * (llh_term1 + llh_term2);



// survival llh ---------------------------------------
if ((j < npt - 1) & (j > 0)){
AB_vec_next = as<uvec>(AB_list[j + 1]); 
// ids of subjects whose z_i(t_{j + 1}) is 0 or 1 (at the next time point). 
// That is, subjects who are alive at the current time point.
AB_vec_next -= 1;
mt_next = AB_vec_next.n_elem;
if (mt_next <= 0) break;

S_v = zeros(d_v, d_v, mt_next);
chol_terms_v = zeros(d_v, d_v, mt_next);

P2_tt_inv = inv(P2_tt);

P2_tt_inv = (P2_tt_inv + P2_tt_inv.t()) / 2.0;

// S_0, V_0
temp = P2_tt + mt_next * P3_tt;
S_u = P1_tt * (P2_tt_inv - mt_next * P2_tt_inv * P3_tt * inv(temp));
var_mat = P0_tt - mt_next * S_u * P1_tt.t();
var_mat = (var_mat + var_mat.t()) / 2.0;
chol_terms_u = my_chol(var_mat);

// S_k, V_k for k = 1, ..., mt_next - 2
if (mt_next > 2){
for (int k = 1; k <= mt_next - 2; ++k){
temp = P2_tt + (mt_next - k) * P3_tt;
S_v.slice(k - 1) = P3_tt * (P2_tt_inv - (mt_next - k) * P2_tt_inv * P3_tt * inv(temp));
var_mat = (P2_tt + P3_tt) - (mt_next - k) * S_v.slice(k - 1) * P3_tt;
var_mat = (var_mat + var_mat.t()) / 2.0;
chol_terms_v.slice(k - 1) = my_chol(var_mat);
}
}


// S_{mt_next - 1}, V_{mt_next - 1}
if (mt_next > 1){
temp = P2_tt + P3_tt;
S_v.slice(mt_next - 2) = P3_tt * inv(temp);
var_mat = (P2_tt + P3_tt) - S_v.slice(mt_next - 2) * P3_tt;
var_mat = (var_mat + var_mat.t()) / 2.0;
chol_terms_v.slice(mt_next - 2) = my_chol(var_mat);
}


// S_{mt_next} = 0, V_{mt_next}
var_mat = P2_tt + P3_tt;
var_mat = (var_mat + var_mat.t()) / 2.0;
chol_terms_v.slice(mt_next - 1) = my_chol(var_mat);


// Generating \alpha_t
// non-missing a_{t|t}, only consider those whose z_i missing at t + 1, z_i(t + 1) can be 0 or 1

a_u_tt_nonNA = a_u_tt;
a_v_tt_nonNA = a_v_tt.cols(AB_vec_next);

Z_nonNA = vec(mt_next);
for (int i = 0; i < mt_next; ++i){
Z_nonNA(i) = as_scalar(z_mat(AB_vec_next(i), j + 1));
}
alpha_v_gen = mat(d_v, mt_next);
p_z_sim = vec(2 * nsim, fill::zeros);

// simulate alpha_t
for (int k = 0; k < nsim; ++k){

// term mt_next
cond_mean = a_v_tt_nonNA.col(mt_next - 1);
cond_chol = chol_terms_v.slice(mt_next - 1);
alpha_v_gen.col(mt_next - 1) = cond_mean + cond_chol * randn(d_v);

// terms mt_next - 1, ..., 1 --------------------------------------------

// sum_diff = sum_{k = i + 1}^{mt_next}, where i = mt - 1, ..., 1 is the subject index, starting from 1
// initialize sum_diff
sum_diff = alpha_v_gen.col(mt_next - 1) - a_v_tt_nonNA.col(mt_next - 1);
for (int i = mt_next - 1; i >= 1; --i){
sum_diff += alpha_v_gen.col(i) - a_v_tt_nonNA.col(i);
//sum_diff = sum(alpha_v_gen.cols(i, mt_next - 1) - a_v_tt_nonNA.cols(i, mt_next - 1), 1);
cond_mean = a_v_tt_nonNA.col(i - 1) + S_v.slice(i - 1) * sum_diff;
cond_chol = chol_terms_v.slice(i - 1);
alpha_v_gen.col(i - 1) = cond_mean + cond_chol * randn(d_v);
}

sum_diff += alpha_v_gen.col(0) - a_v_tt_nonNA.col(0);
//sum_diff = sum(alpha_v_gen - a_v_tt_nonNA, 1);
cond_mean = a_u_tt_nonNA + S_u * sum_diff;
cond_chol = chol_terms_u;
alpha_u_gen = cond_mean + cond_chol * randn(d_u);

mat_1mt_next = mat(1, mt_next, fill::ones);

// stacked_mat = ((u^T, v_1^T)^T, ..., (u^T, v_{mt_next}^T)^T)
stacked_mat = zeros(d_u + d_v, mt_next);
stacked_mat.submat(0, 0, d_u - 1, mt_next - 1) = kron(mat_1mt_next, alpha_u_gen);
stacked_mat.submat(d_u, 0, d_u + d_v - 1, mt_next - 1) = alpha_v_gen;
from = 1;
to = p;
theta_t_plus_1 = stacked_mat.t() * D_mat.t() * gamma_vec + X00_cube.slice(j).rows(AB_vec_next) * gamma_x.subvec(from, to) + gamma_x(0);

// p(z_{t + 1} | alpha_t)
z_term1 = (Z_nonNA - 1) % theta_t_plus_1;
z_term2 = log(1 + exp(- theta_t_plus_1));
z_term = z_term1 - z_term2;
p_z = exp(sum(z_term));
from = 2 * k;
p_z_sim(from) = p_z;

// antithetic variable: 2a_tt - \alpha
stacked_mat2 = zeros(d_u + d_v, mt_next);
stacked_mat2.submat(0, 0, d_u - 1, mt_next - 1) = kron(mat_1mt_next, 2 * a_u_tt_nonNA - alpha_u_gen);
stacked_mat2.submat(d_u, 0, d_u + d_v - 1, mt_next - 1) = 2 * a_v_tt_nonNA - alpha_v_gen;
from = 1;
to = p;
theta_t_plus_1 = stacked_mat2.t() * D_mat.t() * gamma_vec + X00_cube.slice(j).rows(AB_vec_next) * gamma_x.subvec(from, to) + gamma_x(0);



// p(z_{t + 1} | alpha_t)
z_term1 = (Z_nonNA - 1) % theta_t_plus_1;
z_term2 = log(1 + exp(- theta_t_plus_1));
z_term = z_term1 - z_term2;
p_z = exp(sum(z_term));
from = 2 * k + 1;
p_z_sim(from) = p_z;

} 

p_hat = mean(p_z_sim);
//cout << p_hat << endl;
if (p_hat > 0) 
llh += log(p_hat); 

}




// ------------------------------------------------------------------------------------------------------


// a_{t + 1}
a_u = T_u * a_u_tt;
a_v = T_v * a_v_tt;

// P_{t + 1}
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

src = '
vec theta_trans;
mat z_mat, D_mat;
cube data_split_cube, X00_cube;

List A_list, AB_list;

theta_trans = as<vec>(theta_trans_in);
z_mat = as<mat>(z_in);
D_mat = as<mat>(D_in);
data_split_cube = as<cube>(data_split_in); // q x npt x m. data_split_in is a 3D array, data_split_in[,, i] is q x npt, the data of subject i at all time points.
X00_cube = as<cube>(X00_cube_in);
A_list = as<List>(A_list_in); // length npt. Storing ids in set A_j at each time point, in ascending order.
AB_list = as<List>(AB_list_in); // length npt - 1, storing ids in set A_j U B_j at time j, j = 2, ..., npt.


return wrap(neg2llh(theta_trans, data_split_cube, z_mat, D_mat, X00_cube, A_list, AB_list));
'



llh_survival_x_OU = cxxfunction(signature(theta_trans_in = "numeric", data_split_in = "numeric", z_in = "numeric", D_in = "numeric", X00_cube_in = "numeric", A_list_in = "numeric", AB_list_in = "numeric"), includes = rcpp_inc, body = src, plugin = "RcppArmadillo")
# parameters:
# theta_trans_in: (log(zeta_vec), log(nu2_vec), log(xi_vec), log(d_vec), l_vec, gamma_vec, gamma_x), gamma_vec: dim(D_in)[1] x 1, gamma_x: (p + 1) x 1
# data_split_in: q x npt x m
# z_in: m x npt: survival status, contains NA. Note: Although z_in[, 1] is input, it is not used by the function.
# X00_cube_in: m x p x npt: covariates of the m subjects at each time point
# A_list_in: a list of A_j, j = 1, ..., npt. Indices subjects whose z_i(t_j) = 0.
# AB_list_in: a list of (A_j U B_j). j = 2, ..., npt. Indices of subjects whose z_i(t_j) = 0 or 1.




















