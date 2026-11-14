from app.study_schema import rows_by_id, probability


def validate(case):
    specimens = rows_by_id(case['specimens'])
    pools = rows_by_id(case['pools'])
    for specimen in specimens.values():
        if probability(specimen['inclusion_probability']) == 0:
            raise ValueError('入样概率必须大于零')
        if specimen['lot'] not in case['lot_population']:
            raise ValueError('未知批次')
    for pool in pools.values():
        if not pool['members'] or len(set(pool['members'])) != len(pool['members']):
            raise ValueError('混样成员必须非空且不重复')
        if not set(pool['members']) <= specimens.keys():
            raise ValueError('混样引用未知果实')
        if not isinstance(pool['positive'], bool):
            raise ValueError('混样结果必须为布尔值')
        for parent in pool['selected_after']:
            if pools[parent]['stage'] >= pool['stage']:
                raise ValueError('后续混样不得依赖未来结果')
    for calibration in case['calibration']:
        lower, upper = calibration['positive_probability']
        if probability(lower) > probability(upper):
            raise ValueError('校准区间逆序')
    return specimens, pools


def preview(case):
    specimens, pools = validate(case)
    appearances = {identifier: 0 for identifier in specimens}
    missing = []
    for pool in pools.values():
        for member in pool['members']:
            appearances[member] += 1
        present = any(row['assay'] == pool['assay'] and row['pool_size'] == len(pool['members'])
                      for row in case['calibration'])
        if not present:
            missing.append(pool['id'])
    return {'snapshot': case['_snapshot'], 'specimens': len(specimens),
            'tests': len(pools), 'positive_tests': sum(pool['positive'] for pool in pools.values()),
            'membership_counts': appearances, 'uncalibrated_pools': missing,
            'descriptive_only': True}
