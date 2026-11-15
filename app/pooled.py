"""混样筛查联合推断：共享果实、抽样机会与稀释/校准误差。

推断语义（对应 docs/contracts.md）：
- 同一果实的缺陷状态同时进入其所有混管的似然，管结果不按独立伯努利处理；
  重复入管不增加样本量，每颗果实只计一次。
- 校准区间是系统性不确定（同一未知阳性概率的上下界，不是逐行独立随机源）：
  推断在三个角点（全部下界 / 中点 / 全部上界）上分别精确计算，报告三点包络，
  并在 calibration_sensitivity 中给出各角点明细，说明校准缺口对采购依据的影响。
- 缺少某 (assay, pool_size) 校准的混管不借用其他大小的校准，其读数从似然中剔除
  并列入 diagnostics。
- 追加拆样（stage/selected_after）与最初抽样不是同一选择机制：它作为同一批果实上
  的联合似然因子进入，不视为新的随机样本，也不再次使用入样概率。
- 批次层面采用有限总体模型：批次缺陷总数 D 取 0..N 的均匀先验；抽样与缺陷状态
  独立时入样概率在似然中约去（诊断中给出设计核算与 HT 交叉核对）。
"""

import math
from itertools import product

from app.study import validate
from app.study_schema import probability

CREDIBLE_LEVEL = 0.9       # 批次比例可信区间的后验概率水平
FEASIBLE_CUTOFF = 0.01     # 样本缺陷数“仍相容”的后验概率阈值
ASSESSABLE_WIDTH = 0.5     # 批次比例区间宽度超过该值视为不足以支撑采购判断
GRID_TOLERANCE = 1e-9      # 后验可忽略的网格点裁剪阈值（相对最大质量）
MAX_FRONTIER = 16          # 精确推断允许的最大前沿果实数
MAX_STEP_WORK = 2000000    # 单管状态数 × 赋值数的上限
CORNERS = ('low', 'mid', 'high')  # 校准盒角点：全部下界 / 中点 / 全部上界
ROUND = 6


# ---------------------------------------------------------------- 数据准备

def _calibration_table(case):
    table = {}
    for row in case['calibration']:
        key = (row['assay'], row['pool_size'], row['bad_count'])
        lower, upper = row['positive_probability']
        table[key] = (probability(lower), probability(upper))
    return table


def _pool_factors(pool, table):
    """一管在三个校准角点下的似然因子表（按实际坏果数索引）；缺校准返回 None。"""
    size = len(pool['members'])
    rows = []
    for bad in range(size + 1):
        entry = table.get((pool['assay'], size, bad))
        if entry is None:
            return None
        rows.append(entry)
    factors = {corner: [] for corner in CORNERS}
    for bad, (qlo, qhi) in enumerate(rows):
        for corner, q in (('low', qlo), ('mid', (qlo + qhi) / 2.0), ('high', qhi)):
            factors[corner].append(q if pool['positive'] else 1.0 - q)
    return factors


def _build_steps(member_lists):
    """按给定管序构造变量消除计划：每管的新果实、保留集合、状态域与赋值模板。"""
    count = len(member_lists)
    future = [set() for _ in range(count + 1)]
    for j in range(count - 1, -1, -1):
        future[j] = future[j + 1] | set(member_lists[j])
    seen = set()
    steps = []
    for j, members in enumerate(member_lists):
        new = [fruit for fruit in members if fruit not in seen]
        seen.update(members)
        keep = future[j + 1]
        active_before = sorted((seen - set(new)) & future[j])
        if len(active_before) > MAX_FRONTIER:
            raise ValueError('混样依赖图过宽，超出精确推断容量')
        domain = []
        for mask in range(1 << len(active_before)):
            domain.append(frozenset(active_before[k] for k in range(len(active_before)) if mask >> k & 1))
        assignments = []
        for mask in range(1 << len(new)):
            defective = frozenset(new[k] for k in range(len(new)) if mask >> k & 1)
            bits = tuple((fruit, fruit in defective) for fruit in new)
            assignments.append((defective, len(defective), bits))
        if max(1, len(domain)) * max(1, len(assignments)) > MAX_STEP_WORK:
            raise ValueError('单管状态空间过大，超出精确推断容量')
        steps.append({
            'member_set': frozenset(members),
            'new': new,
            'keep': keep,
            'assignments': assignments,
            'domain': domain,
            'frontier': len(active_before),
        })
    return steps


# ---------------------------------------------------------------- 精确推断核

def _forward(steps, factors_list, weight0, weight1):
    alpha = [{frozenset(): 1.0}]
    for step, factors in zip(steps, factors_list):
        current = alpha[-1]
        nxt = {}
        mset, keep = step['member_set'], step['keep']
        for state, weight in current.items():
            shared = len(state & mset)
            for defective, added, bits in step['assignments']:
                factor = factors[shared + added]
                if factor == 0.0:
                    continue
                w = weight * factor
                for fruit, bit in bits:
                    w *= weight1[fruit] if bit else weight0[fruit]
                if w == 0.0:
                    continue
                out = (state | defective) & keep
                nxt[out] = nxt.get(out, 0.0) + w
        alpha.append(nxt)
    return alpha


def _backward(steps, factors_list, weight0, weight1):
    count = len(steps)
    beta = [None] * (count + 1)
    beta[count] = {frozenset(): 1.0}
    for j in range(count - 1, -1, -1):
        step, factors = steps[j], factors_list[j]
        nxt = beta[j + 1]
        mset, keep = step['member_set'], step['keep']
        current = {}
        for state in step['domain']:
            shared = len(state & mset)
            total = 0.0
            for defective, added, bits in step['assignments']:
                factor = factors[shared + added]
                if factor == 0.0:
                    continue
                w = factor
                for fruit, bit in bits:
                    w *= weight1[fruit] if bit else weight0[fruit]
                if w == 0.0:
                    continue
                total += w * nxt.get((state | defective) & keep, 0.0)
            if total != 0.0:
                current[state] = total
        beta[j] = current
    return beta


def _marginals(steps, factors_list, weight0, weight1, alpha, beta):
    """每颗果实在其首次入管步的质量分解：fruit -> [x=0 质量, x=1 质量]。"""
    mass = {}
    for j, (step, factors) in enumerate(zip(steps, factors_list)):
        before = alpha[j]
        after = beta[j + 1]
        mset, keep = step['member_set'], step['keep']
        for state, weight in before.items():
            shared = len(state & mset)
            for defective, added, bits in step['assignments']:
                factor = factors[shared + added]
                if factor == 0.0:
                    continue
                w = weight * factor
                for fruit, bit in bits:
                    w *= weight1[fruit] if bit else weight0[fruit]
                if w == 0.0:
                    continue
                tail = after.get((state | defective) & keep, 0.0)
                if tail == 0.0:
                    continue
                contribution = w * tail
                for fruit, bit in bits:
                    entry = mass.setdefault(fruit, [0.0, 0.0])
                    entry[1 if bit else 0] += contribution
    return mass


def _count_likelihood(steps, factors_list, fruit_lot, lot_count):
    """按批次缺陷计数聚合的似然：{计数向量: 质量}（仅含已校准管中的果实）。"""
    zero = (0,) * lot_count
    alpha = {(frozenset(), zero): 1.0}
    for step, factors in zip(steps, factors_list):
        nxt = {}
        mset, keep = step['member_set'], step['keep']
        for (state, counts), weight in alpha.items():
            shared = len(state & mset)
            for defective, added, bits in step['assignments']:
                factor = factors[shared + added]
                if factor == 0.0:
                    continue
                updated = list(counts)
                for fruit, bit in bits:
                    if bit:
                        updated[fruit_lot[fruit]] += 1
                key = ((state | defective) & keep, tuple(updated))
                nxt[key] = nxt.get(key, 0.0) + weight * factor
        if len(nxt) > MAX_STEP_WORK * 4:
            raise ValueError('计数状态空间过大，超出精确推断容量')
        alpha = nxt
    totals = {}
    for (state, counts), weight in alpha.items():
        totals[counts] = totals.get(counts, 0.0) + weight
    return totals


# ---------------------------------------------------------------- 后验工具

def _credible_interval(values, mass, level):
    """最小宽度网格区间，其后验概率不低于 level。values 为升序列表。"""
    count = len(values)
    prefix = [0.0] * (count + 1)
    for i, value in enumerate(values):
        prefix[i + 1] = prefix[i] + mass.get(value, 0.0)
    total = prefix[count]
    best = None
    for i in range(count):
        for j in range(i, count):
            inside = prefix[j + 1] - prefix[i]
            if total > 0.0 and inside / total >= level:
                width = values[j] - values[i]
                if best is None or width < best[2] - 1e-12:
                    best = (values[i], values[j], width)
                break  # 同一 i 继续扩大 j 只会更宽
    if best is None:
        return [values[0], values[-1]]
    return [best[0], best[1]]


def _moments(mass_by_value):
    total = sum(mass_by_value.values())
    if total <= 0.0:
        return 0.0, 0.0
    mean = sum(value * mass for value, mass in mass_by_value.items()) / total
    second = sum(value * value * mass for value, mass in mass_by_value.items()) / total
    return mean, max(0.0, second - mean * mean)


# ---------------------------------------------------------------- 主入口

def analyze(case, audit_budget=None):
    specimens, pools = validate(case)
    for lot, size in case['lot_population'].items():
        if isinstance(size, bool) or not isinstance(size, int) or size <= 0:
            raise ValueError('批次总体必须为正整数')
    if audit_budget is None:
        audit_budget = case.get('audit_budget', 0)
    if isinstance(audit_budget, bool) or not isinstance(audit_budget, int) or audit_budget < 0:
        raise ValueError('复检额度必须是非负整数')

    lots = sorted(case['lot_population'])
    lot_index = {lot: i for i, lot in enumerate(lots)}
    population = {lot: case['lot_population'][lot] for lot in lots}
    sampled = {lot: [s['id'] for s in case['specimens'] if s['lot'] == lot] for lot in lots}
    for lot in lots:
        if len(sampled[lot]) > population[lot]:
            raise ValueError('抽样数超过批次总体：' + lot)
    inclusion = {lot: sorted({s['inclusion_probability'] for s in case['specimens'] if s['lot'] == lot})
                 for lot in lots}

    table = _calibration_table(case)
    order = list(pools)  # 保持资料给定顺序（追加管排在其触发管之后）
    calibrated, uncalibrated = [], []
    factors = {corner: [] for corner in CORNERS}
    for pid in order:
        pool = pools[pid]
        pool_factors = _pool_factors(pool, table)
        if pool_factors is None:
            uncalibrated.append(pool)
        else:
            calibrated.append(pool)
            for corner in CORNERS:
                factors[corner].append(pool_factors[corner])

    tested_ids = sorted({m for pool in calibrated for m in pool['members']})
    tested_set = set(tested_ids)
    untested = {lot: [s for s in sampled[lot] if s not in tested_set] for lot in lots}
    fruits = tested_ids
    fruit_index = {f: i for i, f in enumerate(fruits)}
    fruit_lot = [lot_index[specimens[f]['lot']] for f in fruits]
    member_lists = [[fruit_index[m] for m in pool['members']] for pool in calibrated]
    steps = _build_steps(member_lists)
    max_frontier = max((step['frontier'] for step in steps), default=0)

    # ---- 批次比例网格（批次缺陷总数均匀先验）----
    grids = {lot: [d / population[lot] for d in range(population[lot] + 1)] for lot in lots}
    grid_points = list(product(*(grids[lot] for lot in lots)))
    prior = 1.0 / len(grid_points)

    def weights(theta):
        w1 = [theta[fruit_lot[i]] for i in range(len(fruits))]
        w0 = [1.0 - t for t in w1]
        return w0, w1

    # ---- 计数所需的先验预测权重 w(k) = E_theta[theta^k (1-theta)^(n-k)] ----
    prior_predictive = {}
    for lot in lots:
        n = len(sampled[lot])
        size = population[lot]
        weights_k = []
        for k in range(n + 1):
            acc = 0.0
            for d in range(size + 1):
                theta = d / size
                acc += theta ** k * (1.0 - theta) ** (n - k)
            weights_k.append(acc / (size + 1))
        prior_predictive[lot] = weights_k

    # 未检果实按组合数并入计数（其似然为 1，只增加缺陷计数的多重性）
    def fold_untested(counts):
        totals = {}
        ranges = [range(len(untested[lot]) + 1) for lot in lots]
        for base, weight in counts.items():
            for extra in product(*ranges):
                key = tuple(base[i] + extra[i] for i in range(len(lots)))
                mult = 1
                for i, lot in enumerate(lots):
                    mult *= math.comb(len(untested[lot]), extra[i])
                totals[key] = totals.get(key, 0.0) + weight * mult
        return totals

    # ---- 逐角点精确计算 ----
    corner_mass = {}       # corner -> {theta: 后验未归一化质量}
    corner_marginals = {}  # corner -> {theta: {fruit: [m0, m1]}}
    corner_counts = {}     # corner -> {计数向量: 后验未归一化质量}
    corner_retained = {}
    for corner in CORNERS:
        corner_factors = factors[corner]
        if steps:
            mass = {}
            for theta in grid_points:
                w0, w1 = weights(theta)
                mass[theta] = sum(_forward(steps, corner_factors, w0, w1)[-1].values())
        else:
            mass = {theta: 1.0 for theta in grid_points}
        peak = max(mass.values()) if mass else 0.0
        retained = ([theta for theta in grid_points if mass[theta] >= GRID_TOLERANCE * peak]
                    if peak > 0.0 else list(grid_points))
        corner_retained[corner] = retained
        marginals = {}
        for theta in retained:
            if not steps:
                marginals[theta] = {}
                continue
            w0, w1 = weights(theta)
            alpha = _forward(steps, corner_factors, w0, w1)
            beta = _backward(steps, corner_factors, w0, w1)
            marginals[theta] = _marginals(steps, corner_factors, w0, w1, alpha, beta)
        if steps:
            counts = _count_likelihood(steps, corner_factors, fruit_lot, len(lots))
        else:
            counts = {tuple(0 for _ in lots): 1.0}
        counts = fold_untested(counts)
        count_mass = {}
        for count_vec, lik in counts.items():
            w = 1.0
            for i, lot in enumerate(lots):
                w *= prior_predictive[lot][count_vec[i]]
            if w > 0.0:
                count_mass[count_vec] = lik * w
        corner_mass[corner] = {theta: m * prior for theta, m in mass.items()}
        corner_marginals[corner] = marginals
        corner_counts[corner] = count_mass

    def lot_value_masses(corner, lot):
        idx = lot_index[lot]
        per_value = {}
        for theta, m in corner_mass[corner].items():
            per_value[theta[idx]] = per_value.get(theta[idx], 0.0) + m
        return per_value

    # ---- 单果后验（逐角点，再取包络）----
    fruit_prob_corners = {corner: {} for corner in CORNERS}
    for corner in CORNERS:
        mass = corner_mass[corner]
        marginals = corner_marginals[corner]
        retained = corner_retained[corner]
        total = sum(mass.values())
        for fruit in fruits:
            i = fruit_index[fruit]
            a = sum(marginals[theta].get(i, [0.0, 0.0])[1] for theta in retained)
            b = sum(marginals[theta].get(i, [0.0, 0.0])[0] for theta in retained)
            fruit_prob_corners[corner][fruit] = a / (a + b) if (a + b) > 0.0 else 0.0
        for lot in lots:
            idx = lot_index[lot]
            mean = sum(theta[idx] * m for theta, m in mass.items()) / total if total > 0 else 0.0
            for fruit in untested[lot]:
                fruit_prob_corners[corner][fruit] = mean
    fruit_prob = {}
    for lot in lots:
        for fruit in sampled[lot]:
            values = [fruit_prob_corners[corner][fruit] for corner in CORNERS]
            fruit_prob[fruit] = (min(values), max(values))

    # ---- 批次层面输出（逐角点，再取包络）----
    lot_rows = []
    sensitivity = {}
    design = {}
    unassessable = []
    prevalence_corners = {}
    for lot in lots:
        idx = lot_index[lot]
        values = grids[lot]
        n_sampled = len(sampled[lot])
        n_tested = sum(1 for f in sampled[lot] if f in tested_set)

        # 样本缺陷数可行范围：各角点下后验 >= FEASIBLE_CUTOFF 的并集
        feasible_union = set()
        for corner in CORNERS:
            count_mass = corner_counts[corner]
            total = sum(count_mass.values())
            if total <= 0.0:
                feasible_union.update(range(n_sampled + 1))
                continue
            for k in range(n_sampled + 1):
                inside = sum(m for counts, m in count_mass.items() if counts[idx] == k)
                if inside / total >= FEASIBLE_CUTOFF:
                    feasible_union.add(k)
        if n_sampled == 0:
            sample_defects = None
        elif feasible_union:
            ordered = sorted(feasible_union)
            sample_defects = [ordered[0], ordered[-1]]
        else:
            sample_defects = [0, n_sampled]

        # 批次比例可信区间：各角点最小宽度区间的端点包络
        per_corner = {}
        if n_tested == 0:
            prevalence = None
        else:
            for corner in CORNERS:
                per_corner[corner] = _credible_interval(
                    values, lot_value_masses(corner, lot), CREDIBLE_LEVEL)
            prevalence = [min(per_corner[c][0] for c in CORNERS),
                          max(per_corner[c][1] for c in CORNERS)]
            prevalence_corners[lot] = per_corner

        if sample_defects is None:
            pop_feasible = None
        else:
            pop_feasible = [sample_defects[0], sample_defects[1] + (population[lot] - n_sampled)]

        width = None if prevalence is None else prevalence[1] - prevalence[0]
        assessable = prevalence is not None and width <= ASSESSABLE_WIDTH
        if not assessable:
            unassessable.append(lot)

        if prevalence is not None:
            sensitivity[lot] = {
                'corner_intervals': {corner: per_corner[corner] for corner in CORNERS},
                'envelope': prevalence,
                'width_spread_from_calibration': max(
                    0.0, max(per_corner[c][1] for c in CORNERS)
                    - min(per_corner[c][1] for c in CORNERS)),
            }
        else:
            sensitivity[lot] = {
                'corner_intervals': None,
                'envelope': None,
                'width_spread_from_calibration': None,
            }

        # HT 交叉核对：D_hat = sum_i x_i / pi_i（以单果后验包络代替未知 x_i）
        ht_lo = sum(fruit_prob[f][0] / specimens[f]['inclusion_probability'] for f in sampled[lot])
        ht_hi = sum(fruit_prob[f][1] / specimens[f]['inclusion_probability'] for f in sampled[lot])
        design[lot] = {
            'unique_fruits': n_sampled,
            'tested_fruits': n_tested,
            'pool_memberships': sum(1 for pool in pools.values() for m in pool['members']
                                    if specimens[m]['lot'] == lot),
            'inclusion_probabilities': inclusion[lot],
            'ht_defects_crosscheck': [min(ht_lo, float(population[lot])),
                                      min(ht_hi, float(population[lot]))],
        }

        lot_rows.append({
            'lot': lot,
            'population': population[lot],
            'sampled': n_sampled,
            'tested': n_tested,
            'sample_defects': sample_defects,
            'sample_defects_basis': '可行：任一校准角点下后验概率不低于 %.2f 的样本缺陷数范围（角点并集）' % FEASIBLE_CUTOFF,
            'prevalence': prevalence,
            'prevalence_basis': '可信：各校准角点 %.1f 最小宽度后验区间端点的包络（校准区间按系统性误差同时取界）' % CREDIBLE_LEVEL,
            'population_defects_feasible': pop_feasible,
            'population_defects_basis': '可行：样本缺陷可行范围加上未抽果实全部缺陷的极端情形',
            'assessable': assessable,
        })

    fruit_rows = []
    for lot in lots:
        for fruit in sampled[lot]:
            lower, upper = fruit_prob[fruit]
            if fruit not in fruit_index:
                status = 'untested'
            elif upper <= 0.05:
                status = 'cleared'
            elif lower >= 0.5:
                status = 'suspect'
            else:
                status = 'undetermined'
            fruit_rows.append({
                'specimen': fruit,
                'lot': lot,
                'tested': fruit in fruit_index,
                'defective_probability': [lower, upper],
                'status': status,
            })

    # ---- 仍相容的计数配置（任一角点后验不低于阈值）----
    compatible = []
    all_counts = set()
    for corner in CORNERS:
        all_counts.update(corner_counts[corner])
    for counts in sorted(all_counts):
        per_corner = []
        for corner in CORNERS:
            count_mass = corner_counts[corner]
            total = sum(count_mass.values())
            if total > 0.0:
                per_corner.append(count_mass.get(counts, 0.0) / total)
        if per_corner and max(per_corner) >= FEASIBLE_CUTOFF:
            compatible.append({
                'counts': {lot: counts[lot_index[lot]] for lot in lots},
                'posterior': [min(per_corner), max(per_corner)],
            })
    compatible.sort(key=lambda row: (-row['posterior'][1], tuple(row['counts'][lot] for lot in lots)))
    truncated = len(compatible) > 50
    compatible = compatible[:50]

    # ---- 复检规划：期望后验方差下降（短视、按悲观角点排序）----
    def lot_variance(corner, lot):
        idx = lot_index[lot]
        per = {}
        for theta, m in corner_mass[corner].items():
            per[theta[idx]] = per.get(theta[idx], 0.0) + m
        return _moments(per)[1]

    variance = {(corner, lot): lot_variance(corner, lot) for corner in CORNERS for lot in lots}

    def score_fruit(fruit):
        """返回 (悲观评分, 乐观评分, 各批次下降范围)：评分=期望后验方差下降。"""
        reductions = {lot: [] for lot in lots}
        for corner in CORNERS:
            mass = corner_mass[corner]
            marginals = corner_marginals[corner]
            retained = corner_retained[corner]
            total = sum(mass.values())
            if total <= 0.0:
                continue
            if fruit in fruit_index:
                i = fruit_index[fruit]
                given = {v: {theta: marginals[theta].get(i, [0.0, 0.0])[v] * prior
                             for theta in retained}
                         for v in (0, 1)}
            else:
                idx = lot_index[specimens[fruit]['lot']]
                given = {v: {theta: m * (theta[idx] if v == 1 else 1.0 - theta[idx])
                             for theta, m in mass.items()}
                         for v in (0, 1)}
            for lot in lots:
                idx = lot_index[lot]
                expected = 0.0
                for v in (0, 1):
                    per = {}
                    for theta, m in given[v].items():
                        per[theta[idx]] = per.get(theta[idx], 0.0) + m
                    expected += (sum(given[v].values()) / total) * _moments(per)[1]
                reductions[lot].append(max(0.0, variance[(corner, lot)] - expected))
        per_lot = {}
        for lot in lots:
            values = reductions[lot] or [0.0]
            per_lot[lot] = [min(values), max(values)]
        lower = sum(per_lot[lot][0] for lot in lots)
        upper = sum(per_lot[lot][1] for lot in lots)
        return lower, upper, per_lot

    scores = {}
    for lot in lots:
        for fruit in sampled[lot]:
            scores[fruit] = score_fruit(fruit)
    ranking = sorted(scores, key=lambda f: (-scores[f][0], -scores[f][1], f))
    selected = [f for f in ranking if scores[f][1] > 0.0][:audit_budget]
    candidates = []
    for fruit in selected:
        lower, upper, per_lot = scores[fruit]
        candidates.append({
            'specimen': fruit,
            'lot': specimens[fruit]['lot'],
            'tested': fruit in fruit_index,
            'defective_probability': list(fruit_prob[fruit]),
            'expected_variance_reduction': [lower, upper],
            'reduces': {lot: per_lot[lot] for lot in lots if per_lot[lot][1] > 0.0},
        })

    # ---- 诊断 ----
    membership = {}
    for pid in order:
        for m in pools[pid]['members']:
            membership.setdefault(m, []).append(pid)
    shared = [{'specimen': fruit, 'pools': pids}
              for fruit, pids in sorted(membership.items()) if len(pids) > 1]

    assumptions = [
        '同一果实的缺陷状态联合进入其所有混管的似然，管结果不按独立伯努利处理',
        '重复入管不增加样本量：每批次样本量按唯一果实计，不按入管次数计',
        '追加拆样（stage>=2）是最初样本上的联合似然因子，不是新的随机抽样，也不再次使用入样概率',
        '校准区间按系统性误差处理：所有校准行同时取下界/中点/上界三个角点，报告三点包络，'
        '不把每行区间当作独立随机来源',
        '缺少校准的混管读数已从似然剔除，未借用其他混样大小的校准',
        '批次缺陷总数取 0..N 的均匀先验；抽样与缺陷状态独立时入样概率在似然中约去',
        '复检规划假设单果复检可确定该果状态，评分为短视（非自适应）期望方差下降，按悲观角点排序',
    ]

    result = {
        'snapshot': case.get('_snapshot'),
        'lots': lot_rows,
        'identification': {
            'compatible_count_configurations': compatible,
            'compatible_configurations_truncated': truncated,
            'fruit_posteriors': fruit_rows,
            'notes': [
                '共享果实使各管读数相互依赖，不能按管独立解释；仍相容的批次缺陷计数配置见 compatible_count_configurations',
                'fruit_posteriors 给出每果缺陷概率的校准角点包络；untested 表示该果不在任何已校准管中',
                '追加管与初抽管在联合似然中一并处理，追加管不构成独立样本',
            ],
        },
        'followup': {
            'budget': audit_budget,
            'selected_specimens': selected,
            'candidates': candidates,
            'method': '短视期望后验方差下降（单果复检视为确定该果状态），按校准角点悲观值排序，取前 budget 个不同果实',
            'note': '同一混管内的果实相关，联合复检的下降不超过各果评分之和；评分仅用于排序，不改变已报告的后验区间',
        },
        'diagnostics': {
            'uncalibrated_pools': [
                {'pool': p['id'], 'assay': p['assay'], 'pool_size': len(p['members']),
                 'members': p['members']} for p in uncalibrated],
            'excluded_pool_readings': len(uncalibrated),
            'unassessable_lots': unassessable,
            'untested_specimens': {lot: untested[lot] for lot in lots if untested[lot]},
            'shared_specimens': shared,
            'design': design,
            'calibration_sensitivity': sensitivity,
            'assumptions': assumptions,
            'parameters': {
                'credible_level': CREDIBLE_LEVEL,
                'feasible_cutoff': FEASIBLE_CUTOFF,
                'assessable_width': ASSESSABLE_WIDTH,
                'grid_tolerance': GRID_TOLERANCE,
                'calibration_corners': list(CORNERS),
                'prior': '批次缺陷总数 0..N 均匀',
            },
            'compute': {
                'grid_points': len(grid_points),
                'grid_retained': {corner: len(corner_retained[corner]) for corner in CORNERS},
                'max_frontier': max_frontier,
                'calibrated_pools': len(calibrated),
            },
        },
    }
    return _clean(result)


def _clean(obj):
    if isinstance(obj, bool) or obj is None or isinstance(obj, (int, str)):
        return obj
    if isinstance(obj, float):
        if not math.isfinite(obj):
            raise ValueError('输出含非有限数值')
        return round(obj, ROUND)
    if isinstance(obj, dict):
        return {key: _clean(value) for key, value in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(value) for value in obj]
    raise ValueError('输出含不支持类型：' + type(obj).__name__)
