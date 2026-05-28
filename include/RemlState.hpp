/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * RemlState — shared post-REML descriptor used by both MLMA_stream and
 * MLMA_loco.  Holds either a dense V^{-1} (GOBY path) or the Woodbury
 * low-rank factors (TUNA path), plus the fixed-effect vector b.
 *
 * This header intentionally has no heavy includes so it can be included
 * by both legacy (main/) and v2 (src/) translation units.
 */
#pragma once

#include "mlma_woodbury.hpp"
#include <Eigen/Dense>
#include <cstdint>

struct RemlState {
    int32_t n           = 0;
    int32_t x_c         = 0;
    bool    is_woodbury = false;
    // LLT path: is_llt=true means Vi_L_f holds lower Cholesky of V (not V^{-1}).
    // Apply V^{-1} via two triangular solves: L^{-T}(L^{-1}·v), avoiding dpotri.
    bool    is_llt      = false;
    Eigen::MatrixXf Vi;      // GOBY: full n×n V^{-1}  (is_woodbury=false, is_llt=false)
    Eigen::MatrixXf Vi_L_f;  // LLT:  lower Cholesky of V (is_llt=true)
    WoodburyMLMACache wb;    // TUNA: Woodbury low-rank factors (is_woodbury=true)
    Eigen::VectorXf b;       // fixed-effect coefficients (x_c elements)
};
