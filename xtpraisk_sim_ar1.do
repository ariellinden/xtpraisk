*==============================================================================*
*  xtpraisk simulation study -- Do-file 1: AR(1) primary
*
*  Estimators : xtpraisk lag(1)  vs  xtscc lag(1)
*  Design     : 3 scenarios x 3 N x 6 T x 4 b1 = 216 conditions
*               2,000 replications each = 432,000 datasets
*
*  DGP (Beck & Katz 1995; Hoechle 2007):
*    y_it = b0 + b1*x_it + u_it
*    u_it = rho1*u_i,t-1 + eps_it,  eps_it ~ iid N(0,1)
*    x_it ~ iid N(0,1) redrawn each replication
*    u_i1 ~ N(0, sigma^2/(1-rho1^2))  [stationary initialization]
*
*  Rho structures follow Linden (2026b):
*    S1: rho = 0.4  (mild positive)
*    S2: rho = -0.4 (oscillatory)
*    S3: rho = 0.7  (high persistent)
*
*  Output columns:
*    scenario  N_panels  T_periods  b1_true  method
*    reject    coverage  pct_bias   rmse     se_ratio
*  where reject = proportion of replications rejecting H0: b1=0
*    => type I error when b1_true==0, power when b1_true!=0
*
*  Output file: xtpraisk_sim_ar1.dta
*==============================================================================*

clear all
set more off
set seed 20260201

local b0    = 10
local sigma = 1
local reps  = 2000
local alpha = 0.05

*--- AR(1) scenarios: Linden (2026a)
local rho_s1 =  0.4
local rho_s2 = -0.4
local rho_s3 =  0.7

local N_list  "10 15 20"
local T_list  "10 20 30 50 75 100"
local b1_list "0 0.05 0.10 0.20 0.30"

postfile sim_res          ///
    int    scenario       ///
    int    N_panels       ///
    int    T_periods      ///
    double b1_true        ///
    int    method         ///  1=xtpraisk  2=xtscc
    double reject         ///  rejection rate (type I when b1=0, power when b1!=0)
    double coverage       ///  95% CI coverage
    double pct_bias       ///  percentage bias
    double rmse           ///  root mean squared error
    double se_ratio       ///  mean(estimated SE) / SD(b1hat)
    using "C:\Users\Ariel\Desktop\ITSA_stuff\xtprais\xtpraisk_sim_ar1.dta", replace

foreach scen of numlist 1/3 {

    local rho1    = `rho_s`scen''
    local var_u1  = (`sigma'^2) / (1 - `rho1'^2)

    foreach N of numlist `N_list' {
    foreach T of numlist `T_list' {
    foreach b1 of numlist `b1_list' {

        di as txt _newline "AR(1) Scen=`scen' rho=`rho1'  N=`N'  T=`T'  b1=`b1'"

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

        _dots 0, title(Replications) reps(`reps')
        forvalues rep = 1/`reps' {
            quietly {
                clear
                set obs `= `N'*`T''
                gen int panid = ceil(_n/`T')
                gen int t     = mod(_n-1,`T')+1
                xtset panid t

                gen double x = rnormal(0,1)

                gen double u = .
                by panid: replace u = rnormal(0,sqrt(`var_u1')) if _n==1
                by panid: replace u = `rho1'*u[_n-1] + rnormal(0,`sigma') if _n>1

                gen double y = `b0' + `b1'*x + u

                * --- xtpraisk
                xtpraisk y x, lag(1) nolog
                local bpw = _b[x]
                local spw = _se[x]
                local ppw = 2*(1-normal(abs(`bpw'/`spw')))
                local lopw = `bpw' - invnormal(0.975)*`spw'
                local hipw = `bpw' + invnormal(0.975)*`spw'

                * --- xtscc (lag matches xtpraisk specification)
                xtscc y x, lag(1)
                local bdk = _b[x]
                local sdk = _se[x]
                local df  = `T'-1
                local pdk = 2*ttail(`df',abs(`bdk'/`sdk'))
                local lodk = `bdk' - invttail(`df',0.025)*`sdk'
                local hidk = `bdk' + invttail(`df',0.025)*`sdk'
            }

            local sum_b1_pw = `sum_b1_pw' + `bpw'
            local sum_b1_dk = `sum_b1_dk' + `bdk'
            local sum_sq_pw = `sum_sq_pw' + (`bpw'-`b1')^2
            local sum_sq_dk = `sum_sq_dk' + (`bdk'-`b1')^2
            local sum_se_pw = `sum_se_pw' + `spw'
            local sum_se_dk = `sum_se_dk' + `sdk'
            if `ppw'<`alpha'                          local n_rej_pw=`n_rej_pw'+1
            if `lopw'<=`b1' & `b1'<=`hipw'           local n_cov_pw=`n_cov_pw'+1
            if `pdk'<`alpha'                          local n_rej_dk=`n_rej_dk'+1
            if `lodk'<=`b1' & `b1'<=`hidk'           local n_cov_dk=`n_cov_dk'+1
            _dots `rep' 0
        }

        * --- performance measures
        local mb1_pw  = `sum_b1_pw'/`reps'
        local mb1_dk  = `sum_b1_dk'/`reps'
        local rmse_pw = sqrt(`sum_sq_pw'/`reps')
        local rmse_dk = sqrt(`sum_sq_dk'/`reps')
        local mse_pw  = `sum_se_pw'/`reps'
        local mse_dk  = `sum_se_dk'/`reps'

        if `b1'!=0 {
            local bias_pw = 100*(`mb1_pw'-`b1')/`b1'
            local bias_dk = 100*(`mb1_dk'-`b1')/`b1'
        }
        else {
            local bias_pw = 100*`mb1_pw'
            local bias_dk = 100*`mb1_dk'
        }

        * se_ratio = mean(SE_hat) / empirical SD of b1hat
        local ser_pw = `mse_pw'/`rmse_pw'
        local ser_dk = `mse_dk'/`rmse_dk'

        post sim_res (`scen') (`N') (`T') (`b1') (1) ///
            (`n_rej_pw'/`reps') (`n_cov_pw'/`reps') ///
            (`bias_pw') (`rmse_pw') (`ser_pw')
        post sim_res (`scen') (`N') (`T') (`b1') (2) ///
            (`n_rej_dk'/`reps') (`n_cov_dk'/`reps') ///
            (`bias_dk') (`rmse_dk') (`ser_dk')

    }  // b1
    }  // T
    }  // N
}  // scenario

postclose sim_res

*--- label dataset
use "C:\Users\Ariel\Desktop\ITSA_stuff\xtprais\xtpraisk_sim_ar1.dta", clear
label define scen_lbl 1 "Mild positive (rho=0.4)"   ///
                      2 "Oscillatory (rho=-0.4)"     ///
                      3 "High persistent (rho=0.7)"
label define meth_lbl 1 "xtpraisk" 2 "xtscc"
label values scenario scen_lbl
label values method   meth_lbl
label var scenario  "Autocorrelation scenario"
label var N_panels  "Number of panels (N)"
label var T_periods "Number of time periods (T)"
label var b1_true   "True value of b1 (0=null)"
label var method    "Estimator"
label var reject    "Rejection rate (type I error if b1=0; power if b1!=0)"
label var coverage  "95% CI coverage rate"
label var pct_bias  "Percentage bias in b1"
label var rmse      "Root mean squared error of b1"
label var se_ratio  "SE ratio: mean(SE_hat)/SD(b1hat)"
save "C:\Users\Ariel\Desktop\ITSA_stuff\xtprais\xtpraisk_sim_ar1.dta", replace

di as txt _newline "Do-file 1 (AR(1) primary) complete."
