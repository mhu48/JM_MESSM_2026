# local level model, longitudinal part only, \Sigma unstructured
# estimate initial delta using Rosenberg approach

###################################################
# mistakes I made:
# (1) integer division: 1 / m should be 1.0 / m
###################################################

library(Rcpp)
library(inline)
rcpp_inc = '
#include<iostream>
using namespace Rcpp;
using namespace arma;
using namespace std;


vec delta_est(vec& theta_trans, cube& data_split_cube, List& nonNAid_list);
// theta_trans = (log(zeta_vec, nu2_vec, xi_vec, d_vec), l_vec)
// returns a vector, delta_hat

List filter_est(vec& theta_trans, cube& data_split_cube, List& nonNAid_list);
// returns -2 * llh

// implementation of funciton delta_est ---------------------------------------------------------------------------
vec delta_est(vec& theta_trans, cube& data_split_cube, List& nonNAid_list){
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

uvec idx, nonNAid_vec; // for creating L from l_vec

// -----------------------------------------
// Initializing necessary quantities
// -----------------------------------------
// q, npt, m, p, mat_1m
q = data_split_cube.n_rows;   // data_split_cube is q x npt x m
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
nonNAid_vec = as<uvec>(nonNAid_list[j]); // IDs start from 1
nonNAid_vec -= 1; // make IDs start from 0
mt = nonNAid_vec.n_elem;

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
sum_v_a = sum(v_a.cols(nonNAid_vec), 1);


// a_{t|t}
C0 = (P0 * Z_u.t() + P1 * Z_v.t()) * F1mF2;
C1 = P2 * Z_v.t() * F1;
C2 = (P1.t() * Z_u.t() + P3 * Z_v.t()) * F1mF2 + P2 * Z_v.t() * F2;

// v_{a, 1} + ... + v_{a, m}, sum of nonNA v_i
sum_v_a = sum(v_a.cols(nonNAid_vec), 1);
// a_{a, tt}
a_u_a_tt = a_u_a;
a_v_a_tt = a_v_a;
a_u_a_tt += C0 * sum_v_a;
a_v_a_tt += C1 * v_a + kron(mat_1m, C2 * sum_v_a);



// augmented part
for (int k = 0; k < d_u; ++k){
sum_V_A = sum(V_A.slice(k).cols(nonNAid_vec), 1);
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
V_temp = mat(V_A.subcube(0, nonNAid_vec[i], 0, q - 1, nonNAid_vec[i], d_u - 1));
sum_V_d += V_temp;
sum_VFV += V_temp.t() * F1 * V_temp;
sum_VFv += V_temp.t() * F1 * v_a.col(nonNAid_vec[i]);
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


delta_hat = - inv(delta_term1) * delta_term2;
return delta_hat;
}

// Implementation of funciton neg2llh ----------------------------------------------------------------------
List filter_est(vec& theta_trans, cube& data_split_cube, List& nonNAid_list){

vec theta, delta, zeta_vec, nu2_vec, xi_vec, d_vec, l_vec, sum_v, sum_r, sum_u, a_u, a_u_tt;

mat data_split_j, D, L, L_inv, tmp, tmp_lower, Z00, T00, P0, P1, P2, P3, Z_u, Z_v, T_u, T_v, Q00,
Sigma, Lambda, Sig_u, Sig_v, a_u_mat, f_pred_mat, a_v_i, a_v_tt_i, a_u_tt_mat, f_filter_mat,
F1mF2, F1, F2, K0, K1, K2, 
A, B, AmB, C0, C1, C2, M0, M1, M2, M3, 
a_tt_split, P0_tt, P1_tt, P2_tt, P3_tt,
Za, a_v, a_v_tt, v_split, Ta_split, 
mat_1m, X_t;

cube a_v_cube, b_pred_cube, a_v_tt_cube, b_filter_cube;

int q, npt, m, mt, d_u, d_v, from, to;

double dt, llh, llh_const, llh_term1, llh_term2; // val, sign, log_det_F1mF2, log_det_F1;

uvec idx, nonNAid_vec; // for creating L from l_vec




// -----------------------------------------
// Convert input quantities to C++ types
// -----------------------------------------

// estimate delta using delta_est function
delta = delta_est(theta_trans, data_split_cube, nonNAid_list);

// -----------------------------------------
// Initializing necessary quantities
// -----------------------------------------
// q, npt, m, mat_1m
q = data_split_cube.n_rows;        // number of longitudinal variables
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

// initialize a1
a_u = delta;
a_v = zeros(d_v, m);

// Sigma
Sigma = L * D * L.t();


// -------------------------------------------------------------------------
// Filtering.
// -------------------------------------------------------------------------
// initialize log-likelihood
llh = 0.0;
a_u_mat = zeros(d_u, npt);
a_v_cube = zeros(d_v, m, npt);
a_u_tt_mat = zeros(d_u, npt);
a_v_tt_cube = zeros(d_v, m, npt);

// loop through all the time points
for (int j = 0; j < npt; ++j){

a_u_mat.col(j) = a_u;
a_v_cube.slice(j) = a_v;

nonNAid_vec = as<uvec>(nonNAid_list[j]); // IDs start from 1
nonNAid_vec -= 1; // make IDs start from 0
mt = nonNAid_vec.n_elem;

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
sum_v = sum(v_split.cols(nonNAid_vec), 1);


// a_{t|t}
C0 = (P0 * Z_u.t() + P1 * Z_v.t()) * F1mF2;
C1 = P2 * Z_v.t() * F1;
C2 = (P1.t() * Z_u.t() + P3 * Z_v.t()) * F1mF2 + P2 * Z_v.t() * F2;
a_u_tt = a_u;
a_v_tt = a_v;
a_u_tt += C0 * sum_v;
a_v_tt += C1 * v_split + kron(mat_1m, C2 * sum_v);

// store the filtering estimates
a_u_tt_mat.col(j) = a_u_tt;
a_v_tt_cube.slice(j) = a_v_tt;


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
llh_term1 += as_scalar(v_split.col(nonNAid_vec(k)).t() * F1 * v_split.col(nonNAid_vec(k)));
}


// llh_term2
llh_term2 = as_scalar(sum_v.t() * F2 * sum_v);
// constant part of llh
llh_const = - (q * mt) / 2.0 * log(2.0 * datum::pi);
llh += llh_const + 0.5 * (log(det(F1mF2)) + (mt - 1.0) * log(det(F1))) - 0.5 * (llh_term1 + llh_term2);

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

// one step ahead prediction for fixed and random effects
f_pred_mat = Z_u * a_u_mat;
b_pred_cube = zeros(q, m, npt);

// filtering estimate for fixed and random effects
f_filter_mat = Z_u * a_u_tt_mat;
b_filter_cube = zeros(q, m, npt);


for (int i = 0; i < m; ++i){
  a_v_i = a_v_cube.subcube(0, i, 0, d_v - 1, i, npt - 1);
  b_pred_cube.subcube(0, i, 0, q - 1, i, npt - 1) = Z_v * a_v_i;
  a_v_tt_i = a_v_tt_cube.subcube(0, i, 0, d_v - 1, i, npt - 1);
  b_filter_cube.subcube(0, i, 0, q - 1, i, npt - 1) = Z_v * a_v_tt_i;
}
return List::create(
Named("a_u_mat") = a_u_mat,
Named("a_v_cube") = a_v_cube,
Named("a_u_tt_mat") = a_u_tt_mat,
Named("a_v_tt_cube") = a_v_tt_cube,
Named("f_pred_mat") = f_pred_mat,
Named("b_pred_cube") = b_pred_cube,
Named("f_filter_mat") = f_filter_mat,
Named("b_filter_cube") = b_filter_cube);
}

'

src = '
vec theta_trans;
mat Z_mat;
cube data_split_cube;

List nonNAid_list;

theta_trans = as<vec>(theta_trans_in);
data_split_cube = as<cube>(data_split_in); // q x npt x m. data_split_in is a 3D array, data_split_in[,, i] is q x npt, the data of subject i at all time points.
nonNAid_list = as<List>(nonNAid_list_in); // length npt. Storing nonNA ids at each time point, in ascending order.

return wrap(filter_est(theta_trans, data_split_cube, nonNAid_list));
'



filter_estimate_OU = cxxfunction(signature(theta_trans_in = "numeric", data_split_in = "numeric", nonNAid_list_in = "numeric"), includes = rcpp_inc, body = src, plugin = "RcppArmadillo")

