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



double neg2llh(vec& theta_trans, mat& data, int& q);
// returns -2 * llh


// Implementation of funciton neg2llh ----------------------------------------------------------------------
double neg2llh(vec& theta_trans, mat& data, int& q){

vec theta, zeta_vec, lambda_vec, d_vec, l_vec, a_tilde, delta, nu2_vec, xi_vec;

mat data_trans, L, tmp_lower, L_inv, A, A_lower, N, N_to_k, C_inv, Z0, Z_f, Z_b, one_m, Z,
T0, T_u, T_v, T, Z00, T00, Q00,
Lambda, Sig_u, Sig_v, Q, P_tilde, Sig_tilde, H_tilde,
Z_tilde_here, K_tilde, tmp;

int d_u, d_v, npt, dim_alpha, from, to, row_from, row_to, col_from, col_to;

double dt, llh, v_tilde, F_tilde, m;

uvec idx, idx1, idx2;


// some quantities
npt = data.n_cols; // data is qm x n
d_u = 2 * q;
d_v = q;
dt = 1.0;
m = data.n_rows / q * 1.0;
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

// OU parameters
from = q;
to = 2 * q - 1;
nu2_vec = theta.subvec(from, to);

// OU parameters
from = 2 * q;
to = 3 * q - 1;
xi_vec = theta.subvec(from, to);

// diagonal elements of D, where Sigma = LDt(L)
from = 3 * q;
to = 4 * q - 1;
d_vec = theta.subvec(from, to);
//D = diagmat(d_vec);

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



// compute C_inv
L_inv = inv(L);
C_inv = kron(eye(m, m), L_inv);

// transform the original data to univariate model (with missing values) ------------------------------------
data_trans = C_inv * data;

// system matrices-------------------------------------------------------------------
dim_alpha = d_u + m * d_v;

// observation matrix Z_tilde
Z00 = zeros(1, 2);
Z00(0, 0) = 1;
Z_f = kron(eye(q, q), Z00);
Z_b = eye(q, q);

// m Z_f matrices (population level) in the first 2q columns
one_m = ones(m, 1);
Z = join_rows(kron(one_m, Z_f), kron(eye(m, m), Z_b));

// transformed Z (for univariate model)
mat Z_tilde;
Z_tilde = C_inv * Z;

// state transition matrix T ----------------------------------------------
dt = 1.0;
T00 = eye(2, 2);
T00(0, 1) = dt;
T_u = kron(eye(q, q), T00);
T_v = diagmat(exp(- xi_vec * dt));

row_from = col_from = 0;
row_to = col_to = d_u - 1;
T = zeros(dim_alpha, dim_alpha);
T.submat(row_from, col_from, row_to, col_to) = T_u;
row_from = col_from = d_u;
row_to = col_to = dim_alpha - 1;
T.submat(row_from, col_from, row_to, col_to) = kron(eye(m, m), T_v);

// innovation covariance matrix Q --------------------------------------------
Lambda = zeros(2, 2);
Lambda(0, 0) = dt * dt * dt / 3;
Lambda(0, 1) = Lambda(1, 0) = dt * dt / 2;
Lambda(1, 1) = dt;
Sig_u = kron(diagmat(zeta_vec), Lambda);
Sig_v = diagmat(0.5 * nu2_vec % (1 - exp(-2 * xi_vec * dt)) / xi_vec );

// Q
Q = zeros(dim_alpha, dim_alpha);
row_from = col_from = 0;
row_to = col_to = d_u - 1;
Q.submat(row_from, col_from, row_to, col_to) = Sig_u;
row_from = col_from = d_u;
row_to = col_to = dim_alpha - 1;
Q.submat(row_from, col_from, row_to, col_to) = kron(eye(m, m), Sig_v);


// initial variance matrix P^*(t_1) (denoted as P.tilde) of state vector \alpha^*(t_1) for univariate model -----------------------
P_tilde = mat(dim_alpha, dim_alpha, fill::zeros);
row_from = col_from = d_u;
row_to = col_to = dim_alpha - 1;
// the diagonal block of Q0
Q00 = diagmat(0.5 * nu2_vec / xi_vec);  
P_tilde.submat(row_from, col_from, row_to, col_to) = kron(eye(m, m), Q00);


// initial mean a^*(t_1) (denoted as a.tilde) of state vector \alpha^*(t_1) for univariate model------------------------------
a_tilde = zeros<vec>(dim_alpha);
from = 0;
to = d_u - 1;
a_tilde.subvec(from, to) = delta;

// observational error variance matrix H.tilde (transformed for univariate model)
Sig_tilde = diagmat(d_vec);
H_tilde = kron(eye(m, m), Sig_tilde);


//--------------------------------------------
// Filtering algorithm and log-likelihood
//--------------------------------------------

llh = 0.0; // initialize log-likelihood

// loop through all the timepoints and data.trans points

for (int j = 0; j < npt; j++){
for (int i = 0; i < q * m; i++){
// if missing value, go to next observation
// if (isnan(data_trans(i, j))) continue;
if (data_trans(i, j) != data_trans(i, j)) continue;
// else filtering
Z_tilde_here  = Z_tilde.row(i);

// v_{t, i}
v_tilde = data_trans(i, j) - as_scalar(Z_tilde_here * a_tilde);
F_tilde = as_scalar(Z_tilde_here * P_tilde * Z_tilde_here.t()) + H_tilde(i, i);

llh = llh - 0.5 * log(2 * datum::pi) - 0.5 * log(F_tilde) - 0.5 * v_tilde * v_tilde / F_tilde;

K_tilde = P_tilde * Z_tilde_here.t() / F_tilde;


a_tilde += K_tilde * v_tilde;

P_tilde -= K_tilde * F_tilde * K_tilde.t();
P_tilde = 0.5 * (P_tilde + P_tilde.t());


}

// transition from {t, qm + 1} to {t + 1, 1}
a_tilde = T * a_tilde;
P_tilde = T * P_tilde * T.t() + Q;
}

return -2.0 * llh;
}

'

src = '
vec theta_trans;
mat data;
int q;

// convert input quantities to C++ types
theta_trans = as<vec>(theta_trans_in);
data = as<mat>(data_in);
q = as<int>(q_in);

return wrap(neg2llh(theta_trans, data, q));
'



llh_OU_unitrt = cxxfunction(signature(theta_trans_in = "numeric", data_in = "numeric", q_in = "numeric"), includes = rcpp_inc, body = src, plugin = "RcppArmadillo")


