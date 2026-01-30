//Changed the objective function to step by step estimation of sigmaB and sigmaE

#include <RcppArmadillo.h>
#include <iostream>
#include <armadillo>
#include <math.h>


// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace std;
using namespace arma;

// Function to compute Kronecker product between identity matrix and vector
// [[Rcpp::export]]
arma::sp_mat kron_iv(int n, const arma::mat& u) {
    arma::vec v = vectorise(u);
    int m = v.n_elem;  // Length of vector v

    // Initialize the result matrix with size (n * m) by n
    arma::mat result(n * m, n, fill::zeros);

    // Fill the result matrix: Each block row corresponds to the identity matrix row multiplied by v
    for (int i = 0; i < n; ++i) {
        result(arma::span(i * m, (i + 1) * m - 1),i) = v;
    }

    return arma::sp_mat(result);
}

// Utility function to compute (i, j)-th index of y_ij
// [[Rcpp::export]]
int y_ij(int i, int j, const arma::ivec& nis) {
    return (i == 1) ? j : accu(nis.subvec(0, i - 2)) + j;
}

// Utility function to compute (i, j)-th indices for Z
// [[Rcpp::export]]
List z_ij(int i, int j, const arma::ivec& nis, int q) {
    int result_i = (i == 1) ? j : accu(nis.subvec(0, i - 2)) + j;
    int start_j = q * (i - 1) + 1;
    int end_j = q * i;
    return List::create(Named("result_i") = result_i,
                        Named("result_j_start") = start_j,
                        Named("result_j_end") = end_j);
}

// Compute diagonal blocks of W_ik kronecker W_ij
// [[Rcpp::export]]
arma::mat diag_ww(const arma::mat& Z, int i, int j, int k, const arma::ivec& nis, int q){
    List z_indices_j = z_ij(i, j, nis, q);
    int result_i_j = as<int>(z_indices_j["result_i"]);
    int result_j_start = as<int>(z_indices_j["result_j_start"]) - 1;
    int result_j_end = as<int>(z_indices_j["result_j_end"]) - 1;

    List z_indices_k = z_ij(i, k, nis, q);
    int result_i_k = as<int>(z_indices_k["result_i"]);
    int result_k_start = as<int>(z_indices_k["result_j_start"]) - 1;
    int result_k_end = as<int>(z_indices_k["result_j_end"]) - 1;

    arma::mat zik = Z(result_i_k - 1, arma::span(result_k_start,result_k_end));
    arma::mat zij = Z(result_i_j - 1, arma::span(result_j_start,result_j_end));

    return kron(zik.t() * zik, zij.t() * zij);
}

// Compute diagonal blocks of Hessian (triple sum of WW)
// [[Rcpp::export]]
arma::mat triple_sum_diag_ww(const arma::mat& Z, const arma::ivec& nis, int m, int q){
    arma::mat result(q * q, q * q);
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= nis[i - 1]; j++) {
            for (int k = 1; k <= nis[i - 1]; k++) {
                if (j != k){
                    result += diag_ww(Z,i,j,k,nis,q);
                }
            }
        }
    }
    return result;
}

// make a dense diagonal sparse
// [[Rcpp::export]]
arma::sp_mat triple_sum_transform(const arma::mat& Z, const arma::ivec& nis, int m, int d, int q){
    arma::mat block = triple_sum_diag_ww(Z, nis, m, q);
    arma::sp_mat I = arma::speye(d, d);
    arma::sp_mat big_matrix(d*q*q, d*q*q);
    for(int i = 1; i <= q; i++){
        for(int j = 1; j <= q; j++){
            arma::mat submatrix = block(arma::span((i - 1) * q, i * q - 1),arma::span((j - 1) * q, j * q - 1));
            arma::sp_mat new_block = arma::kron(I,arma::sp_mat(submatrix));
            big_matrix(arma::span((i - 1) * d * q, i * d * q - 1),
                       arma::span((j - 1) * d * q, j * d * q - 1)) = new_block;
        }
    }
    return big_matrix;
}

// Compute triple.sum.ww
// [[Rcpp::export]]
arma::sp_mat triple_sum_ww(const arma::mat& Z, const arma::ivec& nis, int m, int d, int q) {
    arma::sp_mat I = arma::speye(d, d);
    arma::sp_mat block = triple_sum_transform(Z, nis, m, d, q);

    return arma::kron(I, block);
}


// Compute double.sum.y
// [[Rcpp::export]]
arma::mat double_sum_y(const arma::mat& Y, const arma::mat& X, const arma::mat& beta, const arma::ivec& nis, int m, int d) {
    arma::mat result = arma::zeros(d, d);
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= nis[i - 1]; j++) {
            arma::rowvec yij = Y.row(y_ij(i, j, nis) - 1);
            arma::rowvec xij = X.row(y_ij(i, j, nis) - 1);
            arma::colvec residual = yij.t() - beta.t() * xij.t();
            result += residual * residual.t();
        }
    }
    return result;
}

// Compute double.sum.z.b
// [[Rcpp::export]]
arma::sp_mat double_sum_z_b(const arma::mat& Z, const arma::mat& sigmaB, const arma::ivec& nis, int m, int d, int q) {
    arma::sp_mat result(d, d);
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= nis[i - 1]; j++) {
            List z_indices = z_ij(i, j, nis, q);
            int result_i = as<int>(z_indices["result_i"]);
            int result_j_start = as<int>(z_indices["result_j_start"]) - 1;
            int result_j_end = as<int>(z_indices["result_j_end"]) - 1;

            arma::mat submatrix = Z(result_i - 1, arma::span(result_j_start,result_j_end));
            arma::sp_mat kron_sub = kron_iv(d, submatrix);
            result += kron_sub.t() * sigmaB * kron_sub;
        }
    }
    return result;
}

// Compute triple.sum.y
// [[Rcpp::export]]
arma::sp_mat triple_sum_y(const arma::mat& Y, const arma::mat& Z, const arma::mat& X, const arma::mat& beta, const arma::ivec& nis, int m, int d, int q) {
    arma::sp_mat result(d * q, d * q);
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= nis[i - 1]; j++) {
            for (int k = 1; k <= nis[i - 1]; k++) {
                if (j != k){
                    List z_indices_j = z_ij(i, j, nis, q);
                    List z_indices_k = z_ij(i, k, nis, q);

                    int result_i_j = as<int>(z_indices_j["result_i"]) - 1;
                    int result_i_k = as<int>(z_indices_k["result_i"]) - 1;

                    int result_j_start = as<int>(z_indices_j["result_j_start"]) - 1;
                    int result_j_end = as<int>(z_indices_j["result_j_end"]) - 1;

                    int result_k_start = as<int>(z_indices_k["result_j_start"]) - 1;
                    int result_k_end = as<int>(z_indices_k["result_j_end"]) - 1;

                    arma::mat submatrix_j = Z(result_i_j, arma::span(result_j_start,result_j_end));
                    arma::mat submatrix_k = Z(result_i_k, arma::span(result_k_start,result_k_end));

                    arma::sp_mat kron_j = kron_iv(d, submatrix_j);
                    arma::sp_mat kron_k = kron_iv(d, submatrix_k);

                    arma::rowvec yij = Y.row(result_i_j);
                    arma::rowvec xij = X.row(result_i_j);
                    arma::colvec residual_ij = yij.t() - beta.t() * xij.t();

                    arma::rowvec yik = Y.row(result_i_k);
                    arma::rowvec xik = X.row(result_i_k);
                    arma::colvec residual_ik = yik.t() - beta.t() * xik.t();

                    result += kron_j * (-residual_ij * residual_ik.t()) * kron_k.t();
                }
            }
        }
    }
    return result;
}

// Compute the covariance matrix \Sigma for vec(Y_i^T)
// [[Rcpp::export]]
arma::mat cov_yi(const arma::mat& Z, const arma::mat& sigmaB, const arma::mat& sigmaE,
                 const arma::ivec& nis, int i, int d, int q){
    int ni = nis(i-1);
    arma::mat result(ni * d, ni * d);
    for (int j = 1; j <= ni; j++){
        for (int k = 1; k <= ni; k++){
            List z_indices_j = z_ij(i, j, nis, q);
            int result_i_j = as<int>(z_indices_j["result_i"]);
            int result_j_start = as<int>(z_indices_j["result_j_start"]) - 1;
            int result_j_end = as<int>(z_indices_j["result_j_end"]) - 1;

            List z_indices_k = z_ij(i, k, nis, q);
            int result_i_k = as<int>(z_indices_k["result_i"]);
            int result_k_start = as<int>(z_indices_k["result_j_start"]) - 1;
            int result_k_end = as<int>(z_indices_k["result_j_end"]) - 1;

            arma::mat zik = Z(result_i_k - 1, arma::span(result_k_start,result_k_end));
            arma::mat zij = Z(result_i_j - 1, arma::span(result_j_start,result_j_end));

            arma::mat block = kron_iv(d,zij).t() * sigmaB * kron_iv(d,zik);
            result(arma::span((j - 1) * d, j * d - 1),
                   arma::span((k - 1) * d, k * d - 1)) = block;
        }
    }
    arma::sp_mat I = arma::speye(ni, ni);
    result = result + arma::kron(I, arma::sp_mat(sigmaE));
    return result;
}

// [[Rcpp::export]]
List triple_sum_for_estimate(const arma::mat& Y, const arma::mat& X, const arma::mat& Z, const arma::mat& beta,
                           const arma::ivec& nis, int m, int d, int q,
                           const List &active_set_mat, int n){
    arma::sp_mat z4(n * n, n * n);
    arma::sp_mat zzy(n, n);
    for (int i = 1; i <= m; i++) {
        for (int j = 1; j <= nis[i - 1]; j++) {
            for (int k = 1; k <= nis[i - 1]; k++) {
                if (j != k){
                    List z_indices_j = z_ij(i, j, nis, q);
                    List z_indices_k = z_ij(i, k, nis, q);

                    int result_i_j = as<int>(z_indices_j["result_i"]) - 1;
                    int result_i_k = as<int>(z_indices_k["result_i"]) - 1;

                    int result_j_start = as<int>(z_indices_j["result_j_start"]) - 1;
                    int result_j_end = as<int>(z_indices_j["result_j_end"]) - 1;

                    int result_k_start = as<int>(z_indices_k["result_j_start"]) - 1;
                    int result_k_end = as<int>(z_indices_k["result_j_end"]) - 1;

                    arma::mat Zij = Z(result_i_j, arma::span(result_j_start,result_j_end));
                    arma::mat Zik = Z(result_i_k, arma::span(result_k_start,result_k_end));

                    arma::sp_mat Zij_tilde(n, d);
                    arma::sp_mat Zik_tilde(n, d);

                    int A = 0;
                    for (int l = 0; l < d; ++l) {
                        arma::uvec indices_l = active_set_mat[l];
                        int A_l = indices_l.n_elem;
                        if (A_l == 0){
                            continue;
                        }

                        Zij_tilde(arma::span(A, A + A_l - 1), l) = Zij.elem(indices_l);
                        Zik_tilde(arma::span(A, A + A_l - 1), l) = Zik.elem(indices_l);

                        A = A + A_l;
                    }

                    //For triple.sum.z4
                    arma::sp_mat ZZk = Zik_tilde * Zik_tilde.t();
                    arma::sp_mat ZZj = Zij_tilde * Zij_tilde.t();

                    z4 += arma::kron(ZZk, ZZj);

                    //For triple.sum.zzy
                    arma::rowvec yij = Y.row(result_i_j);
                    arma::rowvec xij = X.row(result_i_j);
                    arma::colvec residual_ij = yij.t() - beta.t() * xij.t();

                    arma::rowvec yik = Y.row(result_i_k);
                    arma::rowvec xik = X.row(result_i_k);
                    arma::colvec residual_ik = yik.t() - beta.t() * xik.t();

                    zzy += Zij_tilde * (-residual_ij * residual_ik.t()) * Zik_tilde.t();
                }
            }
        }
    }

    List output = List::create(
            Named("z4") = z4,
            Named("zzy") = zzy
            );

    return output;
}
