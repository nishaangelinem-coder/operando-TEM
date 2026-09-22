"""Verify that every number quoted in the manuscript matches its source table.

A manuscript and its tables drift apart silently: a study is re-run, a table
updates, and the prose keeps the old value.  This script pins the
load-bearing numbers in manuscript/manuscript.md to the CSV files in
results/tables and fails if any has moved.

Run from the repository root:  python tests/check_manuscript_numbers.py
"""
import csv
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
def rows(name):
    with open(f'results/tables/table_{name}.csv') as fh:
        lines=[l for l in fh if not l.startswith('#')]
    return list(csv.DictReader(lines))
def index(name, key):
    return {str(float(x[key])) if _num(x[key]) else x[key]: x for x in rows(name)}
def _num(v):
    try: float(v); return True
    except Exception: return False

checks=[]
def chk(label, tv, qv, tol=0.02):
    try:
        t=float(tv); q=float(qv)
        ok = abs(t-q) <= tol*max(abs(t),1e-12)+1e-12
    except Exception:
        ok=False; t=tv; q=qv
    checks.append((ok,label,tv,qv))

r=index('S1_transport_calibration','flow_sccm')
chk('S1 tau_fit @2sccm', r['2.0']['tau_fit_mean'],'3.203')
chk('S1 sigma_fit @2sccm', r['2.0']['sigma_fit_mean'],'0.8386')
r={x['descriptor']:x for x in rows('S2a_descriptor_accuracy')}
chk('S2a N1 rel_bias', r['N1']['rel_bias'],'-0.412')
chk('S2a tau_dimer rel_bias', r['tau_dimer']['rel_bias'],'-0.1225',0.03)
r={x['descriptor_set']:x for x in rows('S2b_variance_explained')}
chk('S2b census_only r2', r['census_only']['r2_mean'],'0.9860')
chk('S2b census+dyn r2', r['census_plus_dynamic']['r2_mean'],'0.9866')
r=index('S3a_event_density_sweep','n_pt_atoms')
chk('S3a lamW @8', r['8.0']['event_density_lambda_W'],'0.33')
chk('S3a dec/true @8', r['8.0']['ratio_dec_over_true'],'0.76')
chk('S3a dec/true @6', r['6.0']['ratio_dec_over_true'],'0.85')
r={x['correction']:x for x in rows('S4a_transport_correction')}
chk('S4a deconv ratio', r['deconvolve']['recovery_ratio'],'0.49')
chk('S4a shift ratio', r['shift']['recovery_ratio'],'0.29')
r={x['correction']:x for x in rows('S4b_lead_lag')}
chk('S4b raw peak lag', r['raw']['peak_lag_s'],'3.05')
r={x['instrument']:x for x in rows('S4c_observational_vs_interventional')}
chk('S4c multi ratio weak', r['weak']['multivariable_over_target'],'1.069')
chk('S4c uni ratio weak', r['weak']['univariate_over_target'],'2.13')
chk('S4c first_stage strong', r['strong']['first_stage_r2'],'0.051',0.05)
r=index('S5_projection_bias','corrug_amp_nm')
chk('S5 bias @0.70', r['0.7']['projection_bias_frac'],'-0.155')
r={x['model'][:6]:x for x in rows('S6b_dose_models')}
chk('S6b linear zero-dose', r['linear']['zero_dose_value'],'1.713')
chk('S6b site zero-dose', r['site-a']['zero_dose_value'],'1.169')
r={x['quantity']:x['value'] for x in rows('S7a_representativeness_budget')}
chk('S7a f_rep', r['f_rep'],'4.571e-12')
chk('S7a shortfall', r['field_detector_shortfall'],'6.602e11')
r=index('S7d_patch_vs_chip','dose_rate')
chk('S7d rel diff @1e4', r['10000.0']['rel_difference'],'-0.7413')
r=rows('S7f_position_heterogeneity')[0]
chk('S7f anova p', r['anova_p'],'0.0002123',0.05)
r=rows('S8b_control_paired_test')[0]
chk('S8b held-out ratio', r['held_out_yield_ratio_vs_fixed'],'0.9807')
chk('S8b p vs fixed', r['paired_p_vs_fixed'],'0.6722')
chk('S8b mde relative', r['mde_80pct_relative'],'0.1248')

for ok,l,tv,qv in checks:
    print(('PASS ' if ok else 'FAIL '), f'{l:28s} table={tv} quoted={qv}')
bad=sum(1 for c in checks if not c[0])
print(f'\n{len(checks)-bad}/{len(checks)} numeric cross-checks passed')
raise SystemExit(1 if bad else 0)
