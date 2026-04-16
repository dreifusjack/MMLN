# Shell configuration for git bash compatibility
ifeq ($(OS),Windows_NT)
    SHELL := C:/Program Files/Git/bin/bash.exe
    ifeq ($(wildcard $(SHELL)),)
        SHELL := C:/Program Files/Git/usr/bin/bash.exe
    endif
    ifeq ($(wildcard $(SHELL)),)
        SHELL := C:/Program Files (x86)/Git/bin/bash.exe
    endif
    ifeq ($(wildcard $(SHELL)),)
        SHELL := C:/Program Files (x86)/Git/usr/bin/bash.exe
    endif
    export PATH := C:/Program Files/Git/usr/bin:$(PATH)
else
    SHELL := /bin/bash
endif
.SHELLFLAGS := -euo pipefail -c

.PHONY: help generate-data clone-upstream reference test compare clean clean-all clean-upstream real-data-reference real-data-test

# Colors for output
GREEN  := \033[0;32m
RED    := \033[0;31m
YELLOW := \033[0;33m
BLUE   := \033[0;34m
CYAN   := \033[0;36m
NC     := \033[0m
OK     := [OK]
FAIL   := [FAIL]

# Defaults
FUNC         ?= all
SCENARIO     ?= all
REAL_N_ITER  ?= 2000
REAL_BURN_IN ?= 500

# Paths
SCRIPTS  := $(CURDIR)/scripts
DATA_DIR := $(CURDIR)/data
FIX_DIR  := $(CURDIR)/fixtures
RES_DIR  := $(CURDIR)/results
REP_DIR  := $(CURDIR)/report
REPO     := $(CURDIR)/..
UPSTREAM := $(CURDIR)/upstream

help: ## Show available commands
	@printf "$(CYAN)====================================================\n$(NC)"
	@printf "$(CYAN)  MMLN Optimization Test Suite\n$(NC)"
	@printf "$(CYAN)====================================================\n\n$(NC)"
	@printf "  $(GREEN)%-25s$(NC) %s\n" "help"      "Show available commands"
	@printf "  $(GREEN)%-25s$(NC) %s\n" "check-deps" "Check required R packages are installed"
	@printf "  $(GREEN)%-25s$(NC) %s\n" "test"      "Run fork tests and compare (use FUNC= and SCENARIO= to filter)"
	@printf "  $(GREEN)%-25s$(NC) %s\n" "compare"   "Re-run comparison only (no model runs)"
	@printf "  $(GREEN)%-25s$(NC) %s\n" "clean"     "Remove fork results and reports (keep fixtures and data)"
	@printf "  $(GREEN)%-25s$(NC) %s\n" "clean-upstream" "Remove the local upstream/ clone"
	@printf "\n$(YELLOW)  Variables:$(NC)\n"
	@printf "    FUNC=fmln|mmln|mdres|all    Function to test (default: all)\n"
	@printf "    SCENARIO=name|\"n1 n2\"|all   Scenario(s) to test (default: all)\n"
	@printf "\n$(YELLOW)  Setup (generally one-time commands)$(NC)\n\n"
	@printf "  $(GREEN)%-25s$(NC) %s\n" "clone-upstream" "Clone upstream MMLN into upstream/ for local editing"
	@printf "  $(GREEN)%-25s$(NC) %s\n" "generate-data" "Generate synthetic test datasets (no package install needed)"
	@printf "  $(GREEN)%-25s$(NC) %s\n" "reference"     "Install upstream from local copy and generate reference fixtures"
	@printf "  $(GREEN)%-25s$(NC) %s\n" "clean-all"     "Remove everything (including fixtures and data)"
	@printf "\n$(YELLOW)  Scenarios (all baseball-scale, N~1200-2000, n_iter=1000):$(NC)\n"
	@printf "    $(GREEN)baseline$(NC)        MLB-mirrored (N=2000, m=50, p=10, d=3)\n"
	@printf "    $(GREEN)zero_counts$(NC)     Many zero cells across 8 categories (N=1500, m=30, d=7, n_mean=10)\n"
	@printf "    $(GREEN)tiny_counts$(NC)     Extreme ALR values, wide W variance (N=1200, m=40, d=3, n_mean=5)\n"
	@printf "    $(GREEN)large_counts$(NC)    Precision in exp/log1p accumulates (N=1200, m=30, d=3, n_mean=5000)\n"
	@printf "    $(GREEN)high_dim$(NC)        9 categories, 8x8 covariance (N=1200, m=30, p=8, d=8)\n"
	@printf "    $(GREEN)near_singular$(NC)   Sigma/Phi with rho=0.92 (N=1400, m=40, d=4)\n"
	@printf "    $(GREEN)large_coefs$(NC)     beta in [-4,4], exp(W) overflow risk (N=1200, m=30, p=8, d=3)\n"
	@printf "    $(GREEN)unbalanced$(NC)      Group sizes 1-80, single-obs groups (N~2000, m=50, d=3)\n"
	@printf "    $(GREEN)many_groups$(NC)     500 random intercepts like MLB players (N=2000, m=500, d=2)\n"
	@printf "    $(GREEN)dominant_cat$(NC)    One category ~90%% prob, ALR near -inf (N=1200, m=30, d=2)\n"
	@printf "\n$(YELLOW)  Examples:$(NC)\n"
	@printf "    make test FUNC=fmln SCENARIO=baseline\n"
	@printf "    make test FUNC=mdres\n"
	@printf "    make test FUNC=mmln SCENARIO=\"baseline zero_counts\"\n"
	@printf "\n$(YELLOW)  Real-data benchmarks (upstream vs fork on actual pollen + MLB data):$(NC)\n\n"
	@printf "  $(GREEN)%-25s$(NC) %s\n" "real-data-reference" "Install upstream; fit FMLN (pollen) + MMLN (MLB), save to fixtures/"
	@printf "  $(GREEN)%-25s$(NC) %s\n" "real-data-test"      "Install fork;     fit FMLN (pollen) + MMLN (MLB), compare timing"
	@printf "\n$(YELLOW)  Real-data variables:$(NC)\n"
	@printf "    REAL_N_ITER=<n>   MCMC iterations for real-data runs (default: 2000)\n"
	@printf "    REAL_BURN_IN=<n>  Burn-in for real-data runs           (default: 500)\n"
	@printf "\n$(YELLOW)  Real-data examples:$(NC)\n"
	@printf "    make real-data-reference\n"
	@printf "    make real-data-test\n"
	@printf "    make real-data-test REAL_N_ITER=3000\n"
	@printf "\n"

check-deps: ## Check required dependencies (system + R packages)
	@printf "$(BLUE)Checking dependencies...$(NC)\n\n"
	@printf "$(YELLOW)  System tools:$(NC)\n"
	@command -v Rscript >/dev/null 2>&1 && \
		printf "  $(GREEN)$(OK)$(NC) Rscript:     $$(Rscript --version 2>&1 | head -1)\n" || \
		printf "  $(RED)$(FAIL)$(NC) Rscript:     Not found (sudo apt-get install r-base)\n"
	@command -v g++ >/dev/null 2>&1 && \
		printf "  $(GREEN)$(OK)$(NC) g++:         $$(g++ --version 2>&1 | head -1)\n" || \
		printf "  $(RED)$(FAIL)$(NC) g++:         Not found (sudo apt-get install build-essential)\n"
	@command -v make >/dev/null 2>&1 && \
		printf "  $(GREEN)$(OK)$(NC) make:        $$(make --version 2>&1 | head -1)\n" || \
		printf "  $(RED)$(FAIL)$(NC) make:        Not found (sudo apt-get install build-essential)\n"
	@command -v cmake >/dev/null 2>&1 && \
		printf "  $(GREEN)$(OK)$(NC) cmake:       $$(cmake --version 2>&1 | head -1)\n" || \
		printf "  $(RED)$(FAIL)$(NC) cmake:       Not found (sudo apt-get install cmake)\n"
	@printf "\n$(YELLOW)  R packages:$(NC)\n"
	@Rscript -e 'if (requireNamespace("devtools", quietly=TRUE)) cat("  ok\n") else cat("  missing\n")' 2>/dev/null | \
		grep -q "ok" && \
		printf "  $(GREEN)$(OK)$(NC) devtools\n" || \
		printf "  $(RED)$(FAIL)$(NC) devtools     (install.packages('devtools'))\n"
	@Rscript -e 'if (requireNamespace("mvnfast", quietly=TRUE)) cat("  ok\n") else cat("  missing\n")' 2>/dev/null | \
		grep -q "ok" && \
		printf "  $(GREEN)$(OK)$(NC) mvnfast\n" || \
		printf "  $(RED)$(FAIL)$(NC) mvnfast      (install.packages('mvnfast'))\n"
	@Rscript -e 'if (requireNamespace("Rcpp", quietly=TRUE)) cat("  ok\n") else cat("  missing\n")' 2>/dev/null | \
		grep -q "ok" && \
		printf "  $(GREEN)$(OK)$(NC) Rcpp\n" || \
		printf "  $(RED)$(FAIL)$(NC) Rcpp         (install.packages('Rcpp'))\n"
	@Rscript -e 'if (requireNamespace("RcppArmadillo", quietly=TRUE)) cat("  ok\n") else cat("  missing\n")' 2>/dev/null | \
		grep -q "ok" && \
		printf "  $(GREEN)$(OK)$(NC) RcppArmadillo\n" || \
		printf "  $(RED)$(FAIL)$(NC) RcppArmadillo (install.packages('RcppArmadillo'))\n"
	@Rscript -e 'if (requireNamespace("RcppEigen", quietly=TRUE)) cat("  ok\n") else cat("  missing\n")' 2>/dev/null | \
		grep -q "ok" && \
		printf "  $(GREEN)$(OK)$(NC) RcppEigen\n" || \
		printf "  $(RED)$(FAIL)$(NC) RcppEigen    (install.packages('RcppEigen'))\n"
	@printf "\n$(YELLOW)  R packages (real-data benchmarks):$(NC)\n"
	@Rscript -e 'if (requireNamespace("MM", quietly=TRUE)) cat("  ok\n") else cat("  missing\n")' 2>/dev/null | \
		grep -q "ok" && \
		printf "  $(GREEN)$(OK)$(NC) MM           (pollen data)\n" || \
		printf "  $(RED)$(FAIL)$(NC) MM           (needs libgmp-dev, then install.packages('MM'))\n"
	@Rscript -e 'if (requireNamespace("Lahman", quietly=TRUE)) cat("  ok\n") else cat("  missing\n")' 2>/dev/null | \
		grep -q "ok" && \
		printf "  $(GREEN)$(OK)$(NC) Lahman       (MLB data)\n" || \
		printf "  $(RED)$(FAIL)$(NC) Lahman       (install.packages('Lahman'))\n"
	@Rscript -e 'if (requireNamespace("dplyr", quietly=TRUE)) cat("  ok\n") else cat("  missing\n")' 2>/dev/null | \
		grep -q "ok" && \
		printf "  $(GREEN)$(OK)$(NC) dplyr        (MLB data cleaning)\n" || \
		printf "  $(RED)$(FAIL)$(NC) dplyr        (install.packages('dplyr'))\n"
	@printf "\n"

generate-data: ## Generate synthetic test datasets (no package install needed)
	@printf "$(BLUE)Generating synthetic test datasets...$(NC)\n"
	@Rscript "$(SCRIPTS)/generate_data.R" "$(REPO)"
	@printf "  $(GREEN)$(OK)$(NC) Datasets saved to data/\n\n"

reference: ## Install upstream from local copy and generate reference fixtures (use FUNC= and SCENARIO= to filter)
	@printf "$(BLUE)Installing upstream MMLN from upstream/...$(NC)\n"
	@Rscript -e 'Rcpp::compileAttributes("$(UPSTREAM)")'
	@Rscript -e 'devtools::install("$(UPSTREAM)", force = TRUE, quiet = TRUE)'
	@printf "  $(GREEN)$(OK)$(NC) Upstream installed\n\n"
	@printf "$(BLUE)Generating reference fixtures (FUNC=$(FUNC), SCENARIO=$(SCENARIO))...$(NC)\n"
	@Rscript "$(SCRIPTS)/generate_reference.R" "$(REPO)" "$(FUNC)" "$(SCENARIO)"
	@printf "\n  $(GREEN)$(OK)$(NC) Reference fixtures saved to fixtures/\n\n"

test: ## Run fork tests and compare (use FUNC= and SCENARIO= to filter)
	@printf "$(BLUE)Installing fork MMLN (local)...$(NC)\n"
	@Rscript -e 'Rcpp::compileAttributes("$(REPO)")'
	@Rscript -e 'devtools::install("$(REPO)", force = TRUE, quiet = TRUE)'
	@printf "  $(GREEN)$(OK)$(NC) Fork installed\n\n"
	@printf "$(BLUE)Running fork tests (FUNC=$(FUNC), SCENARIO=$(SCENARIO))...$(NC)\n"
	@Rscript "$(SCRIPTS)/run_test.R" "$(REPO)" "$(FUNC)" "$(SCENARIO)"
	@printf "\n$(BLUE)Comparing results...$(NC)\n"
	@Rscript "$(SCRIPTS)/compare.R" "$(REPO)" "$(FUNC)" "$(SCENARIO)"

compare: ## Re-run comparison only (no model runs)
	@printf "$(BLUE)Comparing results (FUNC=$(FUNC), SCENARIO=$(SCENARIO))...$(NC)\n"
	@Rscript "$(SCRIPTS)/compare.R" "$(REPO)" "$(FUNC)" "$(SCENARIO)"

clean: ## Remove fork results and reports (keep fixtures and data)
	@printf "$(YELLOW)Cleaning fork results and reports...$(NC)\n"
	@rm -rf "$(RES_DIR)"/* "$(REP_DIR)"/*
	@printf "  $(GREEN)$(OK)$(NC) Cleaned results/ and report/\n\n"

clean-all: clean ## Remove everything (including fixtures and data)
	@printf "$(YELLOW)Cleaning fixtures and data...$(NC)\n"
	@rm -rf "$(FIX_DIR)"/* "$(DATA_DIR)"/*
	@printf "  $(GREEN)$(OK)$(NC) Cleaned fixtures/ and data/\n\n"

clean-upstream: ## Remove the local upstream/ clone
	@printf "$(YELLOW)Removing upstream/ clone...$(NC)\n"
	@rm -rf "$(UPSTREAM)"
	@printf "  $(GREEN)$(OK)$(NC) Removed upstream/\n\n"

real-data-reference: ## Install upstream and run FMLN (pollen) + MMLN (MLB), save timing to fixtures/
	@printf "$(BLUE)Installing upstream MMLN from upstream/...$(NC)\n"
	@Rscript -e 'Rcpp::compileAttributes("$(UPSTREAM)")'
	@Rscript -e 'devtools::install("$(UPSTREAM)", force = TRUE, quiet = TRUE)'
	@printf "  $(GREEN)$(OK)$(NC) Upstream installed\n\n"
	@printf "$(BLUE)Running real-data reference (REAL_N_ITER=$(REAL_N_ITER), REAL_BURN_IN=$(REAL_BURN_IN))...$(NC)\n"
	@Rscript "$(SCRIPTS)/run_real_data_reference.R" "$(REPO)" "$(REAL_N_ITER)" "$(REAL_BURN_IN)"
	@printf "\n  $(GREEN)$(OK)$(NC) Reference fixtures saved to fixtures/\n\n"

real-data-test: ## Install fork and run FMLN (pollen) + MMLN (MLB), print timing vs reference
	@printf "$(BLUE)Installing fork MMLN (local)...$(NC)\n"
	@Rscript -e 'Rcpp::compileAttributes("$(REPO)")'
	@Rscript -e 'devtools::install("$(REPO)", force = TRUE, quiet = TRUE)'
	@printf "  $(GREEN)$(OK)$(NC) Fork installed\n\n"
	@printf "$(BLUE)Running real-data fork test (REAL_N_ITER=$(REAL_N_ITER), REAL_BURN_IN=$(REAL_BURN_IN))...$(NC)\n"
	@Rscript "$(SCRIPTS)/run_real_data_test.R" "$(REPO)" "$(REAL_N_ITER)" "$(REAL_BURN_IN)"
