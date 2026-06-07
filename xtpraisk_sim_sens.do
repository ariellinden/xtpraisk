*==============================================================================*
*  xtpraisk simulation study -- Do-file 4: Sensitivity (cross-panel correlation)
*
*  Tests xtpraisk and xtscc when panels are contemporaneously correlated.
*  Cross-panel correlation via common factor (Hoechle 2007):
*    eps_it = lambda * f_t + eta_it
*    f_t   ~ iid N(0,1)  same draw for all N panels at time t
*    eta_it ~ iid N(0,1) idiosyncratic
*
*  Parsimonious design -- AR(2) DGP only (middle case); two scenarios
*  (mild positive and high persistent); reduced T grid; b1 in {0, 1.0} only.
*  Varies N in {10,15,20} and lambda in {0.5,1.0} to examine whether
*  larger N improves PCSE accuracy and how lambda strength affects coverage.
*
*  Design: 2 scenarios x 3 N x 5 T x 2 lambda x 2 b1 = 120 conditions
*          x 2 methods x 2,000 reps = 480,000 datasets
*
*  Rho structures follow Linden (2026b):
*    S1: rho = (0.4,  0.2)  mild positive
*    S3: rho = (0.7,  0.2)  high persistent
*  (oscillatory omitted; cross-panel correlation story does not change
*   qualitatively across the oscillatory scenario)
*
*  Output column "reject" = proportion rejecting H0: b1=0
*    => type I error when b1_true==0, power when b1_true!=0
*
*  Output: xtpraisk_sim_sensitivity.dta
*==============================================================================*

clear all
set more off
set seed 20260204

local b0    = 10
local sigma = 1
local reps  = 2000
local alpha = 0.05

*--- AR(2) scenarios used: S1 (mild positive) and S3 (high persistent)
*    Linden (2026b)
local rho1_s1 =  0.4
        local rho2_s1 =  0.2   // mild positive
local rho1_s3 =  0.7
        local rho2_s3 =  0.2   // high persistent

local N_list   "10 15 20"
local T_list   "10 20 30 50 100"
local lam_list "0.5 1.0"
local b1_list  "0 0.20"

postfile sim_res          ///
    int    scenario       ///  1=mild positive (S1)  3=high persistent (S3)
    int    N_panels       ///
    int    T_periods      ///
    double lambda         ///  common factor loading
    double b1_true        ///
    int    method         ///  1=xtpraisk  2=xtscc
    double reject         ///  rejection rate (type I if b1=0; power if b1!=0)
    double coverage       ///  95% CI coverage
    double pct_bias       ///  percentage bias
    double rmse           ///  root mean squared error
    double se_ratio       ///  mean(SE_hat) / SD(b1hat)
    using "C:\Users\Ariel\Desktop\ITSA_stuff\xtprais\xtpraisk_sim_sensitivity.dta", replace

foreach scen in 1 3 {

    local rho1 = `rho1_s`scen''
    local rho2 = `rho2_s`scen''

    * Galbraith & Galbraith (1974) stationary variance matrix elements
    local g0  = (`sigma'^2)*(1-`rho2') / ///
                ((1+`rho2')*((1-`rho2')^2 - `rho1'^2))
    local g1  = `rho1'*`g0'/(1-`rho2')
    local L11 = sqrt(`g0')
    local L21 = `g1'/`L11'
    local L22 = sqrt(max(`g0' - (`g1'^2)/`g0', 1e-12))

    foreach N of numlist `N_list' {
    foreach T of numlist `T_list' {
    foreach lam of numlist `lam_list' {
    foreach b1 of numlist `b1_list' {

        di as txt _newline ///
            "Sens  Scen=`scen' rho=(`rho1',`rho2')  N=`N'  T=`T'" ///
            "  lambda=`lam'  b1=`b1'"

        local sum_b1_pw=0
            local sum_b1_dk=0
        local sum_sq_pw=0
        local sum_sq_dk=0
        local sum_se_pw=0
        local sum_se_dk=0
        local n_rej_pw=0
        local n_rej_dk=0
        local n_cov_pw=0
        local n_cov_dk=0
        local n_ok=0

        _dots 0, title(Replications) reps(`reps')
        forvalues rep = 1/`reps' {
            quietly {
                *--- generate common factor: one draw per time period
                clear
                set obs `T'
                gen int t = _n
                gen double f_t = rnormal(0,1)
                tempfile factor
                save `factor'

                *--- generate full panel dataset
                clear
                set obs `= `N'*`T''
                gen int panid = ceil(_n/`T')
                gen int t     = mod(_n-1,`T')+1
                merge m:1 t using `factor', nogenerate
                sort panid t

                *--- covariate (independent of factor: strict exogeneity holds)
                gen double x = rnormal(0,1)

                *--- correlated innovations
                gen double eps = `lam'*f_t + rnormal(0,1)
                drop f_t

                *--- AR(2) errors: Galbraith initialization for first 2 obs
                *    then AR(2) filter using correlated innovations
                gen double u = .
                by panid: replace u = `L11'*rnormal() if _n==1
                by panid: replace u = `L21'*(u[_n-1]/`L11') + `L22'*rnormal() ///
                    if _n==2
                by panid: replace u = `rho1'*u[_n-1] + `rho2'*u[_n-2] + eps ///
                    if _n>=3

                xtset panid t
                gen double y = `b0' + `b1'*x + u

                *--- xtpraisk: PCSE sandwich accounts for cross-panel correlation
                            }

            capture quietly xtpraisk y x, lag(2) nolog
            if _rc != 0 {
                _dots `rep' 1
                continue
            }
            local bpw = _b[x]
            local spw = _se[x]
            if missing(`bpw') | missing(`spw') {
                _dots `rep' 1
                continue
            }

            capture quietly xtscc y x, lag(2)
            if _rc != 0 {
                _dots `rep' 1
                continue
            }
            local bdk = _b[x]
            local sdk = _se[x]
            if missing(`bdk') | missing(`sdk') {
                _dots `rep' 1
                continue
            }

            local ppw  = 2*(1-normal(abs(`bpw'/`spw')))
            local lopw = `bpw' - invnormal(0.975)*`spw'
            local hipw = `bpw' + invnormal(0.975)*`spw'
            local df   = `T' - 1
            local pdk  = 2*ttail(`df',abs(`bdk'/`sdk'))
            local lodk = `bdk' - invttail(`df',0.025)*`sdk'
            local hidk = `bdk' + invttail(`df',0.025)*`sdk'

            local sum_b1_pw = `sum_b1_pw' + `bpw'
            local sum_b1_dk = `sum_b1_dk' + `bdk'
            local sum_sq_pw = `sum_sq_pw' + (`bpw'-`b1')^2
            local sum_sq_dk = `sum_sq_dk' + (`bdk'-`b1')^2
            local sum_se_pw = `sum_se_pw' + `spw'
            local sum_se_dk = `sum_se_dk' + `sdk'
            if `ppw'<`alpha'                  local n_rej_pw = `n_rej_pw'+1
            if `lopw'<=`b1'&`b1'<=`hipw'     local n_cov_pw = `n_cov_pw'+1
            if `pdk'<`alpha'                  local n_rej_dk = `n_rej_dk'+1
            if `lodk'<=`b1'&`b1'<=`hidk'     local n_cov_dk = `n_cov_dk'+1
            local n_ok = `n_ok' + 1
            _dots `rep' 0
        }

        *--- performance measures
        local rmse_pw = sqrt(`sum_sq_pw'/`reps')
        local rmse_dk = sqrt(`sum_sq_dk'/`reps')
        local mse_pw  = `sum_se_pw'/`reps'
        local mse_dk  = `sum_se_dk'/`reps'
        local mb1_pw  = `sum_b1_pw'/`reps'
        local mb1_dk  = `sum_b1_dk'/`reps'

        if `b1'!=0 {
            local bias_pw = 100*(`mb1_pw'-`b1')/`b1'
            local bias_dk = 100*(`mb1_dk'-`b1')/`b1'
        }
        else {
            local bias_pw = 100*`mb1_pw'
            local bias_dk = 100*`mb1_dk'
        }

        post sim_res (`scen') (`N') (`T') (`lam') (`b1') (1) ///
            (`n_rej_pw'/`n_ok') (`n_cov_pw'/`n_ok') ///
            (`bias_pw') (`rmse_pw') (`mse_pw'/`rmse_pw')
        post sim_res (`scen') (`N') (`T') (`lam') (`b1') (2) ///
            (`n_rej_dk'/`n_ok') (`n_cov_dk'/`n_ok') ///
            (`bias_dk') (`rmse_dk') (`mse_dk'/`rmse_dk')

    }   // b1
    }   // lam
    }   // T
    }   // N
}   // scenario

postclose sim_res

*--- label output
use "C:\Users\Ariel\Desktop\ITSA_stuff\xtprais\xtpraisk_sim_sensitivity.dta", clear

label define scen_lbl 1 "Mild positive (0.4, 0.2)"   ///
                      3 "High persistent (0.7, 0.2)"
label define meth_lbl 1 "xtpraisk" 2 "xtscc"
label values scenario scen_lbl
label values method   meth_lbl

label var scenario  "Autocorrelation scenario (AR(2) DGP)"
label var N_panels  "Number of panels (N)"
label var T_periods "Number of time periods (T)"
label var lambda    "Common factor loading (cross-panel correlation strength)"
label var b1_true   "True value of b1 (0=null)"
label var method    "Estimator"
label var reject    "Rejection rate (type I error if b1=0; power if b1!=0)"
label var coverage  "95% CI coverage rate"
label var pct_bias  "Percentage bias in b1"
label var rmse      "Root mean squared error of b1"
label var se_ratio  "SE ratio: mean(SE_hat)/SD(b1hat)"

save "C:\Users\Ariel\Desktop\ITSA_stuff\xtprais\xtpraisk_sim_sensitivity.dta", replace

di as txt _newline "Do-file 4 (sensitivity: cross-panel correlation) complete."
