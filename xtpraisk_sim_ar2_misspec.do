*==============================================================================*
*  xtpraisk simulation study -- Do-file 2b: AR(2) misspecification
*
*  MISSPECIFICATION design follows Linden (2026b):
*    DGP=AR(2); both estimators fitted with lag(1) [underspecified]
*    xtpraisk lag(1)  vs  xtscc lag(1)
*
*  xtpraisk and xtscc always use the same lag -- complete alignment.
*  3 scenarios x 3 N x 6 T x 2 b1 = 108 conditions
*  2,000 replications each = 432,000 datasets
*
*  DGP (Beck & Katz 1995; Hoechle 2007):
*    y_it = b0 + b1*x_it + u_it
*    u_it = rho1*u_i,t-1 + rho2*u_i,t-2 + eps_it,  eps_it ~ iid N(0,1)
*    x_it ~ iid N(0,1) redrawn each replication
*    First 2 obs per panel initialized from stationary AR(2) distribution
*    using Galbraith & Galbraith (1974) closed-form V_2 matrix
*
*  Rho structures follow Linden (2026b):
*    S1: rho = (0.4,  0.2)  mild positive
*    S2: rho = (0.5, -0.4)  oscillatory
*    S3: rho = (0.7,  0.2)  high persistent
*
*  Output column "reject" = proportion rejecting H0: b1=0
*    => type I error when b1_true==0, power when b1_true!=0
*
*  Output: xtpraisk_sim_ar2_misspec.dta
*==============================================================================*

clear all
set more off
set seed 20260205

local b0    = 10
local sigma = 1
local reps  = 2000
local alpha = 0.05
local path  "C:\Users\Ariel\Desktop\ITSA_stuff\xtprais"

*--- AR(2) scenarios: Linden (2026b)
local rho1_s1 =  0.4
local rho2_s1 =  0.2
local rho1_s2 =  0.5
local rho2_s2 = -0.4
local rho1_s3 =  0.7
local rho2_s3 =  0.2

local N_list  "10 15 20"
local T_list  "10 20 30 50 75 100"

postfile mis_res          ///
    int    scenario       ///
    int    N_panels       ///
    int    T_periods      ///
    double b1_true        ///
    int    method         ///  1=xtpraisk lag(1)  2=xtscc lag(1)
    double reject         ///
    double coverage       ///
    double pct_bias       ///
    double rmse           ///
    double se_ratio       ///
    using "`path'\xtpraisk_sim_ar2_misspec.dta", replace

foreach scen of numlist 1/3 {

    local rho1 = `rho1_s`scen''
    local rho2 = `rho2_s`scen''

    local g0  = (`sigma'^2)*(1-`rho2') / ///
                ((1+`rho2')*((1-`rho2')^2 - `rho1'^2))
    local g1  = `rho1'*`g0'/(1-`rho2')
    local L11 = sqrt(`g0')
    local L21 = `g1'/`L11'
    local L22 = sqrt(max(`g0' - (`g1'^2)/`g0', 1e-12))

    foreach N of numlist `N_list' {
    foreach T of numlist `T_list' {
    foreach b1 of numlist 0 0.20 {

        di as txt _newline ///
            "AR(2) Misspec  Scen=`scen' rho=(`rho1',`rho2')  N=`N'  T=`T'  b1=`b1'"

        local sum_b1_pw = 0
        local sum_b1_dk = 0
        local sum_sq_pw = 0
        local sum_sq_dk = 0
        local sum_se_pw = 0
        local sum_se_dk = 0
        local n_rej_pw  = 0
        local n_rej_dk  = 0
        local n_cov_pw  = 0
        local n_cov_dk  = 0
        local n_ok      = 0   // count successful reps

        _dots 0, title(Replications) reps(`reps')
        forvalues rep = 1/`reps' {

            local ok = 0

            quietly {
                clear
                set obs `= `N'*`T''
                gen int panid = ceil(_n/`T')
                gen int t     = mod(_n-1,`T')+1
                xtset panid t
                gen double x = rnormal(0,1)
                gen double u = .
                by panid: replace u = `L11'*rnormal() if _n==1
                by panid: replace u = `L21'*(u[_n-1]/`L11') + ///
                    `L22'*rnormal() if _n==2
                by panid: replace u = `rho1'*u[_n-1] + `rho2'*u[_n-2] + ///
                    rnormal(0,`sigma') if _n>=3
                gen double y = `b0' + `b1'*x + u
            }

            * capture estimation errors (e.g. missing values in small samples)
            capture quietly xtpraisk y x, lag(1) nolog
            if _rc != 0 {
                _dots `rep' 1
                continue
            }
            local bpw  = _b[x]
            local spw  = _se[x]
            if missing(`bpw') | missing(`spw') {
                _dots `rep' 1
                continue
            }

            capture quietly xtscc y x, lag(1)
            if _rc != 0 {
                _dots `rep' 1
                continue
            }
            local bdk  = _b[x]
            local sdk  = _se[x]
            if missing(`bdk') | missing(`sdk') {
                _dots `rep' 1
                continue
            }

            local ppw  = 2*(1-normal(abs(`bpw'/`spw')))
            local lopw = `bpw' - invnormal(0.975)*`spw'
            local hipw = `bpw' + invnormal(0.975)*`spw'
            local df   = `T' - 1
            local pdk  = 2*ttail(`df', abs(`bdk'/`sdk'))
            local lodk = `bdk' - invttail(`df', 0.025)*`sdk'
            local hidk = `bdk' + invttail(`df', 0.025)*`sdk'

            local sum_b1_pw = `sum_b1_pw' + `bpw'
            local sum_b1_dk = `sum_b1_dk' + `bdk'
            local sum_sq_pw = `sum_sq_pw' + (`bpw'-`b1')^2
            local sum_sq_dk = `sum_sq_dk' + (`bdk'-`b1')^2
            local sum_se_pw = `sum_se_pw' + `spw'
            local sum_se_dk = `sum_se_dk' + `sdk'
            if `ppw' < `alpha'                  local n_rej_pw = `n_rej_pw' + 1
            if `lopw'<=`b1' & `b1'<=`hipw'     local n_cov_pw = `n_cov_pw' + 1
            if `pdk' < `alpha'                  local n_rej_dk = `n_rej_dk' + 1
            if `lodk'<=`b1' & `b1'<=`hidk'     local n_cov_dk = `n_cov_dk' + 1
            local n_ok = `n_ok' + 1
            _dots `rep' 0
        }

        * use n_ok as denominator to handle any skipped reps
        if `n_ok' > 0 {
            local rmse_pw = sqrt(`sum_sq_pw'/`n_ok')
            local rmse_dk = sqrt(`sum_sq_dk'/`n_ok')
            local mse_pw  = `sum_se_pw'/`n_ok'
            local mse_dk  = `sum_se_dk'/`n_ok'
            local mb1_pw  = `sum_b1_pw'/`n_ok'
            local mb1_dk  = `sum_b1_dk'/`n_ok'
            if `b1' != 0 {
                local bias_pw = 100*(`mb1_pw'-`b1')/`b1'
                local bias_dk = 100*(`mb1_dk'-`b1')/`b1'
            }
            else {
                local bias_pw = 100*`mb1_pw'
                local bias_dk = 100*`mb1_dk'
            }
            post mis_res (`scen') (`N') (`T') (`b1') (1) ///
                (`n_rej_pw'/`n_ok') (`n_cov_pw'/`n_ok') ///
                (`bias_pw') (`rmse_pw') (`mse_pw'/`rmse_pw')
            post mis_res (`scen') (`N') (`T') (`b1') (2) ///
                (`n_rej_dk'/`n_ok') (`n_cov_dk'/`n_ok') ///
                (`bias_dk') (`rmse_dk') (`mse_dk'/`rmse_dk')
        }

    }
    }
    }
}

postclose mis_res

*--- label output
use "`path'\xtpraisk_sim_ar2_misspec.dta", clear
label define scen_lbl 1 "Mild positive (0.4, 0.2)"  ///
                      2 "Oscillatory (0.5, -0.4)"   ///
                      3 "High persistent (0.7, 0.2)"
label define meth_lbl 1 "xtpraisk lag(1) underspecified" ///
                      2 "xtscc lag(1) underspecified"
label values scenario scen_lbl
label values method   meth_lbl
label var scenario  "Autocorrelation scenario"
label var N_panels  "Number of panels (N)"
label var T_periods "Number of time periods (T)"
label var b1_true   "True value of b1 (0=null)"
label var method    "Estimator (both underspecified: DGP=AR(2), fit=AR(1))"
label var reject    "Rejection rate (type I error if b1=0; power if b1!=0)"
label var coverage  "95% CI coverage rate"
label var pct_bias  "Percentage bias in b1"
label var rmse      "Root mean squared error of b1"
label var se_ratio  "SE ratio: mean(SE_hat)/SD(b1hat)"
save "`path'\xtpraisk_sim_ar2_misspec.dta", replace

di as txt _newline "Do-file 2b (AR(2) misspecification) complete."
