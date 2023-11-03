/*
 * This file uses the Catch unit testing library, alongside
 * testthat's simple bindings, to test a C++ function.
 *
 * For your own packages, ensure that your test files are
 * placed within the `src/` folder, and that you include
 * `LinkingTo: testthat` within your DESCRIPTION file.
 */

// All test files should include the <testthat.h>
// header file.
#include <testthat.h>
#include "RcppArmadillo.h"
#include "multinethist.h"

using namespace Rcpp;

context("GroupSize Check") {
  GroupSize G(200,13), F(200, 10);
  test_that("GroupSize constructor?") {
    expect_true(G.number_group == 15);
    expect_true(F.number_group == 20);
    expect_true(G.group_number.n_elem == 15);
    expect_true(F.group_number.n_elem == 20);
    expect_true(G.group_number.at(0) == 13);
    expect_true(F.group_number.at(0) == 10);
    expect_true(G.group_number.at(14) == 18);
    expect_true(F.group_number.at(19) == 10);
  }
  test_that("length()") {
    expect_true(G.length() == 15);
    expect_true(F.length() == 20);
  }
}

context("Assignment check") {
  arma::icube A(10,10,1);
  GroupSize G(10,5), F(10, 3);
  arma::uvec label_F = {0,0,0,1,1,1,2,2,2,2};
  arma::uvec label_G = {0,0,0,0,0,1,1,1,1,1};
  const arma::mat cell_size_F = {{3.0, 9.0, 12.0}, {9.0, 3.0, 12.0}, {12.0, 12.0, 6.0}};
  const arma::mat cell_size_G = {{10.0, 25.0}, {25.0, 10.0}};
  
  A.at(0,1,0) = 1;
  A.at(1,0,0) = 1;
  A.at(0,3,0) = 1;
  A.at(3,0,0) = 1;
  
  Assignment assignF(A, label_F, F);
  Assignment assignG(A, label_G, G);
  Assignment assignF2 = assignF;
  AssignCommonF assignCF_F(A, label_F, F);
  AssignCommonF assignCF_G(A, label_G, G);
  AssignCommonF assignCF_F2 = assignCF_F;
  AssignSingleLayer assignSL_F(A.slice(0), label_F, F);
  AssignSingleLayer assignSL_G(A.slice(0), label_G, G);
  AssignSingleLayer assignSL_F2 = assignSL_F;
  
  test_that("Assignment Constructor") {
    expect_true(assignF.group_size.length() == F.length());
    expect_true(assignG.group_size.length() == G.length());
    expect_true(all(vectorise(assignF.bin_cell_counts == cell_size_F)));
    expect_true(all(vectorise(assignG.bin_cell_counts == cell_size_G)));
  }
  test_that("Assignment comparison") {
    // assignF.likelihood == -5.049031
    // assignG.likelihood == -5.004024
    // mulitnethist
    expect_false(assignF > assignG);
    expect_true(assignF < assignG);
    expect_false(assignF == assignG);
    // homogeneous multinethist
    expect_false(assignCF_F > assignCF_G);
    expect_true(assignCF_F < assignCF_G);
    expect_false(assignCF_F == assignCF_G);
    // single layer nethist
    expect_false(assignSL_F > assignSL_G);
    expect_true(assignSL_F < assignSL_G);
    expect_false(assignSL_F == assignSL_G);
  }
  test_that("Assignment method") {
    arma::uvec expected_label = {0,1,0,0,1,1,2,2,2,2};
    arma::mat expected_theta_hat = {{1.0/3.0, 1.0/9.0, 0},{1.0/9.0,0,0},{0,0,0}};

    //swap_nodes
    assignF2.swap_nodes(assignF, {1, 3});
    expect_true(all(assignF2.node_labels == expected_label));
    assignCF_F2.swap_nodes(assignCF_F, {1, 3});
    expect_true(all(assignCF_F2.node_labels == expected_label));
    assignSL_F2.swap_nodes(assignSL_F, {1, 3});
    expect_true(all(assignSL_F2.node_labels == expected_label));
    //update_thetahat
    assignF2.update_thetahat({1,3}, A);
    expect_true(approx_equal(assignF2.estimated_theta.slice(0), expected_theta_hat, "absdiff", 2*eps));
    assignCF_F2.update_thetahat({1,3}, A);
    expect_true(approx_equal(assignCF_F2.estimated_theta.slice(0), expected_theta_hat, "absdiff", 2*eps));
    assignSL_F2.update_thetahat({1,3}, A);
    expect_true(approx_equal(assignSL_F2.estimated_theta, expected_theta_hat, "absdiff", 2*eps));
    //updateLL
    assignF2.updateLL();
    expect_true(abs(assignF2.likelihood - (-5.049031)) < 1e-6);
    assignCF_F2.updateLL();
    expect_true(abs(assignCF_F2.likelihood - (-5.049031)) < 1e-6);
    assignSL_F2.updateLL();
    expect_true(abs(assignSL_F2.likelihood - (-5.049031)) < 1e-6);
    //copy_labels_theta_LL
    assignF2.copy_labels_theta_LL(assignF);
    expect_true(assignF2==assignF);
    assignCF_F2.copy_labels_theta_LL(assignCF_F);
    expect_true(assignCF_F2==assignCF_F);
    assignSL_F2.copy_labels_theta_LL(assignSL_F);
    expect_true(assignSL_F2==assignSL_F);
  }
}


