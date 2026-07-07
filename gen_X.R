# Calculate -2 * log-likelihood
library(Rcpp)
library(inline)
rcpp_inc = '
using namespace Rcpp;
using namespace arma;

List generate_X0(int& m, int& p, int& q);

List generate_X0(int& m, int& p, int& q){

mat X00, X0;
int row_from, row_to, col_from, col_to;

X00 = randn(m, p);
X0 = zeros(q * m, p * q);
for (int i = 0; i < m; ++i){
row_from = i * q;
col_from = 0;
row_to = i * q + q - 1;
col_to = p * q - 1;
X0.submat(row_from, col_from, row_to, col_to) = kron(eye(q, q), X00.row(i));
}

return List::create(
Named("X00") = X00, // m x p
Named("X0") = X0); // qm x pq
}

'



src = '
int m, p, q;

m = as<int>(m_in);
p = as<int>(p_in);
q = as<int>(q_in);

return wrap(generate_X0(m, p, q));
'

gen_X = cxxfunction(signature(m_in = "numeric", p_in = "numeric", q_in = "numeric"), includes = rcpp_inc, body = src, plugin = "RcppArmadillo")
