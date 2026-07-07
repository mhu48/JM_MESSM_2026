# local level model, longitudinal part only, \Sigma unstructured
# setting delta to 0, carry out one round of filtering

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



double neg2llh(vec& theta_trans, cube& data_split_cube, List& nonNAid_list);
// returns -2 * llh


// Implementation of funciton neg2llh ----------------------------------------------------------------------
double neg2llh(vec& theta_trans, cube& data_split_cube, List& nonNAid_list){

vec theta, delta, zeta_vec, nu2_vec, xi_vec, d_vec, l_vec, sum_v, sum_r, sum_u, a_u, a_u_tt;

mat data_split_j, D, L, L_inv, tmp, tmp_lower, Z00, T00, P0, P1, P2, P3, Z_u, Z_v, T_u, T_v, Q00,
Sigma, Lambda, Sig_u, Sig_v,
F1mF2, F1, F2, K0, K1, K2, 
A, B, AmB, C0, C1, C2, M0, M1, M2, M3, 
a_tt_split, P0_tt, P1_tt, P2_tt, P3_tt,
eps_hat_mat, y_hat, f_hat, b_hat,
Za, a_v, a_v_tt, v_split, Ta_split, Kv_split, Talpha, alpha_hat_split, r_split, eta_hat_split, 
u, u_term1, u_term2, r_term1, r_term2, Pr, eps_hat_split, alpha0_hat, b_hat_split, y_hat_split,
mat_1m, X_t;


int q, npt, m, mt, d_u, d_v, from, to;

double dt, llh, llh_const, llh_term1, llh_term2; // val, sign, log_det_F1mF2, log_det_F1;

uvec idx, nonNAid_vec; // for creating L from l_vec




// -----------------------------------------
// Convert input quantities to C++ types
// -----------------------------------------


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
delta = zeros(d_u, 1);

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
// Filtering. And storing matrices for use in the smoothing step
// -------------------------------------------------------------------------
// initialize log-likelihood
llh = 0.0;

// loop through all the time points
for (int j = 0; j < npt; ++j){
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


return -2.0 * llh;
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

return wrap(neg2llh(theta_trans, data_split_cube, nonNAid_list));
'



llh_OU = cxxfunction(signature(theta_trans_in = "numeric", data_split_in = "numeric", nonNAid_list_in = "numeric"), includes = rcpp_inc, body = src, plugin = "RcppArmadillo")

