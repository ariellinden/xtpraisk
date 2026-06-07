*==============================================================================
* xtpraisk Applied Example
* Prediabetes disease management: panel of 10 regional health systems
* Monthly average fasting blood glucose (mg/dL), T=30 months, N=10 units
* Intervention at month 16 (common to all units)
*
* DGP: y_it = b0 + b1*x_it + u_it
*      b0=108, b1=-2.5, SD=3, tpost=16
*      AR(1): rho=0.7
*      AR(2): rho=(0.7, 0.2)
*==============================================================================

clear all
set more off

local N     = 10
local T     = 30
local b0    = 108
local b1    = -2.5
local sigma = 3
local tpost = 16
set seed 29

* Initialise results matrix: 4 rows x 5 cols (2 AR orders x 2 methods)
matrix results_mat = J(4, 5, .)
matrix colnames results_mat = Coefficient SE p-value CI_lower CI_upper

*==============================================================================
* AR(1): rho = 0.7
*==============================================================================
clear
quietly set obs `=`N'*`T''
quietly gen id   = ceil(_n/`T')
quietly gen time = mod(_n-1,`T')+1
quietly xtset id time
quietly gen x = (time >= `tpost')
quietly gen u = 0
forvalues i = 1/`N' {
    quietly replace u = rnormal(0,`sigma') if id==`i'
    forvalues t = 2/`T' {
        quietly replace u = 0.7*u[_n-1] + rnormal(0,`sigma') ///
            if id==`i' & time==`t'
    }
}
quietly gen y = `b0' + `b1'*x + u

* xtpraisk
quietly xtpraisk y x, lag(1)
quietly lincom x
matrix results_mat[1,1] = r(estimate)
matrix results_mat[1,2] = r(se)
matrix results_mat[1,3] = r(p)
matrix results_mat[1,4] = r(lb)
matrix results_mat[1,5] = r(ub)

* xtscc
quietly xtscc y x, lag(1)
local df1 = e(df_r)
quietly lincom x, df(`df1')
matrix results_mat[2,1] = r(estimate)
matrix results_mat[2,2] = r(se)
matrix results_mat[2,3] = r(p)
matrix results_mat[2,4] = r(lb)
matrix results_mat[2,5] = r(ub)

*==============================================================================
* AR(2): rho = (0.7, 0.2)
*==============================================================================
clear
quietly set obs `=`N'*`T''
quietly gen id   = ceil(_n/`T')
quietly gen time = mod(_n-1,`T')+1
quietly xtset id time
quietly gen x = (time >= `tpost')
quietly gen u = 0
forvalues i = 1/`N' {
    quietly replace u = rnormal(0,`sigma') if id==`i'
    quietly replace u = 0.7*u[_n-1] + rnormal(0,`sigma') ///
        if id==`i' & time==2
    forvalues t = 3/`T' {
        quietly replace u = 0.7*u[_n-1] + 0.2*u[_n-2] + rnormal(0,`sigma') ///
            if id==`i' & time==`t'
    }
}
quietly gen y = `b0' + `b1'*x + u

* xtpraisk
quietly xtpraisk y x, lag(2)
quietly lincom x
matrix results_mat[3,1] = r(estimate)
matrix results_mat[3,2] = r(se)
matrix results_mat[3,3] = r(p)
matrix results_mat[3,4] = r(lb)
matrix results_mat[3,5] = r(ub)

* xtscc
quietly xtscc y x, lag(2)
local df2 = e(df_r)
quietly lincom x, df(`df2')
matrix results_mat[4,1] = r(estimate)
matrix results_mat[4,2] = r(se)
matrix results_mat[4,3] = r(p)
matrix results_mat[4,4] = r(lb)
matrix results_mat[4,5] = r(ub)

*==============================================================================
* Display formatted table
*==============================================================================
di _newline(2)
di as text "Table 3. Applied example: intervention effect (mg/dL) by AR order and method."
di as text "{hline 80}"
di as text %12s "AR Order" %14s "Method" %12s "Coeff" %10s "SE" ///
           %12s "p-value" %22s "95% CI"
di as text "{hline 80}"

forvalues ar = 1/2 {
    if `ar' == 1 local arstr "AR(1)"
    if `ar' == 2 local arstr "AR(2)"
    forvalues m = 1/2 {
        if `m' == 1 local mstr "xtpraisk"
        if `m' == 2 local mstr "xtscc"
        local row = (`ar'-1)*2 + `m'
        local coef = results_mat[`row',1]
        local se   = results_mat[`row',2]
        local pval = results_mat[`row',3]
        local cilo = results_mat[`row',4]
        local cihi = results_mat[`row',5]
        if `pval' < 0.001 local pstr "<0.001"
        else               local pstr = string(`pval',"%6.4f")
        if `m' == 1 {
            di as text %12s "`arstr'" %14s "`mstr'" ///
               as result %12.4f `coef' %10.4f `se' ///
               %12s "`pstr'" "   (" %7.4f `cilo' ", " %7.4f `cihi' ")"
        }
        else {
            di as text %12s ""        %14s "`mstr'" ///
               as result %12.4f `coef' %10.4f `se' ///
               %12s "`pstr'" "   (" %7.4f `cilo' ", " %7.4f `cihi' ")"
        }
    }
    if `ar' < 2 di as text "{hline 80}"
}

di as text "{hline 80}"
di as text "Note: Both methods applied to the same dataset at each AR order."
di as text "      xtpraisk: z-based inference; xtscc: t-based inference (df=e(df_r))."

di _newline as text "Raw matrix:"
matrix list results_mat, format(%9.4f)

*==============================================================================
* End of file
*==============================================================================
