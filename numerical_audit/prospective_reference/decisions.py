"""Predeclared descriptive contrasts and practical correction decision rule."""
import statistics

FLOOR=1e-12
ARMS=('current','pair_reference','optimized_JD','optimized_J0','pauli_JD',
      'pauli_JD_large','jw_signed','jw_magnitude','symmetric_JD')


def contrast(pairs):
    """Pairs are (numerator error, denominator error), with absolute ties."""
    wins=sum(a-b < -FLOOR for a,b in pairs)
    losses=sum(a-b > FLOOR for a,b in pairs)
    ratios=[a/b for a,b in pairs if b>FLOOR]
    return dict(n=len(pairs),wins=wins,ties=len(pairs)-wins-losses,losses=losses,
                ratios_retained=len(ratios),floor_exclusions=len(pairs)-len(ratios),
                median_ratio=statistics.median(ratios) if ratios else None,
                max_ratio=max(ratios) if ratios else None)


def compute(rows,contexts):
    """contexts maps (case,state) to the frozen per-method construction costs."""
    lookup={}
    for r in rows:
        key=(r['case'],r['state'],r['arm'],r['budget_units'])
        if key in lookup:raise ValueError('Duplicate decision input row')
        lookup[key]=r['one_minus_overlap']
    budgets=sorted({r['budget_units'] for r in rows})
    required={(case,state,arm,b) for case,state in contexts for arm in ARMS for b in budgets}
    if set(lookup)!=required or 10 not in budgets:raise ValueError('Incomplete declared comparisons')
    contrasts={}
    for budget in budgets:
        block={}
        for numerator in ('current','pair_reference'):
            for denominator in ARMS:
                if numerator==denominator:continue
                block[numerator+'_over_'+denominator]=contrast([
                    (lookup[case,state,numerator,budget],lookup[case,state,denominator,budget])
                    for case,state in contexts])
        contrasts[str(budget)]=block
    costs=[]
    for (case,state),methods in sorted(contexts.items()):
        totals={arm:methods[arm]['setup_cpu_seconds']+methods[arm]['search_cpu_seconds']
                for arm in ('current','pair_reference')}
        if totals['current']<=0 or totals['pair_reference']<0:raise ValueError('Invalid construction cost')
        costs.append(dict(case=case,state=state,current_cpu=totals['current'],
                          corrected_cpu=totals['pair_reference'],ratio=totals['pair_reference']/totals['current']))
    loss_ratios=[]
    for case,state in sorted(contexts):
        corrected=lookup[case,state,'pair_reference',10];current=lookup[case,state,'current',10]
        if corrected-current>FLOOR:
            loss_ratios.append(dict(case=case,state=state,ratio=corrected/max(current,FLOOR)))
    primary=contrasts['10']['pair_reference_over_current']
    median_cpu=statistics.median(x['ratio'] for x in costs)
    maximum_loss=max((x['ratio'] for x in loss_ratios),default=1.)
    conditions=dict(at_least_four_wins=primary['wins']>=4,
                    no_loss_exceeds_twofold=maximum_loss<=2,
                    median_construction_CPU_at_most_twofold=median_cpu<=2)
    return dict(contrasts=contrasts,costs=costs,primary=primary,strict_loss_ratios=loss_ratios,
                median_construction_CPU_ratio=median_cpu,maximum_strict_loss_ratio=maximum_loss,
                conditions=conditions,correction_passes=all(conditions.values()),tie_and_denominator_floor=FLOOR,
                loss_guard='strict losses use corrected/max(current,1e-12); descriptive ratios exclude denominator<=1e-12')
