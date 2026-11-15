import json
import threading
import unittest

from app.client import analysis
from app.pooled import analyze
from app.server import make_server
from app.study_store import load_case


def tiny_case():
    """两果一正一拆：确定性仪器下结论可手算。"""
    return {
        'kind': 'pools',
        'specimens': [
            {'id': 'a', 'lot': 'L', 'inclusion_probability': 1.0},
            {'id': 'b', 'lot': 'L', 'inclusion_probability': 1.0},
        ],
        'pools': [
            {'id': 'p1', 'members': ['a', 'b'], 'positive': True,
             'assay': 'det', 'stage': 1, 'selected_after': []},
            {'id': 'p2', 'members': ['a'], 'positive': False,
             'assay': 'det', 'stage': 2, 'selected_after': ['p1']},
        ],
        'calibration': [
            {'assay': 'det', 'pool_size': 1, 'bad_count': 0, 'positive_probability': [0.0, 0.0]},
            {'assay': 'det', 'pool_size': 1, 'bad_count': 1, 'positive_probability': [1.0, 1.0]},
            {'assay': 'det', 'pool_size': 2, 'bad_count': 0, 'positive_probability': [0.0, 0.0]},
            {'assay': 'det', 'pool_size': 2, 'bad_count': 1, 'positive_probability': [1.0, 1.0]},
            {'assay': 'det', 'pool_size': 2, 'bad_count': 2, 'positive_probability': [1.0, 1.0]},
        ],
        'lot_population': {'L': 2},
        'audit_budget': 1,
        '_snapshot': 'tiny',
    }


def exact_case():
    """两果一管阳性，点校准（区间退化为点），后验可手算。"""
    return {
        'kind': 'pools',
        'specimens': [
            {'id': 'a', 'lot': 'L', 'inclusion_probability': 1.0},
            {'id': 'b', 'lot': 'L', 'inclusion_probability': 1.0},
        ],
        'pools': [
            {'id': 'p1', 'members': ['a', 'b'], 'positive': True,
             'assay': 'ex', 'stage': 1, 'selected_after': []},
        ],
        'calibration': [
            {'assay': 'ex', 'pool_size': 2, 'bad_count': 0, 'positive_probability': [0.02, 0.02]},
            {'assay': 'ex', 'pool_size': 2, 'bad_count': 1, 'positive_probability': [0.9, 0.9]},
            {'assay': 'ex', 'pool_size': 2, 'bad_count': 2, 'positive_probability': [0.9, 0.9]},
        ],
        'lot_population': {'L': 4},
        'audit_budget': 0,
        '_snapshot': 'exact',
    }


class TinyCaseTests(unittest.TestCase):
    def test_shared_fruit_determined_by_joint_likelihood(self):
        result = analyze(tiny_case())
        probs = {row['specimen']: row['defective_probability']
                 for row in result['identification']['fruit_posteriors']}
        # p1 阳性且 p2（拆样 a）阴性 ⇒ a 清白、b 携带；追加管不视为新样本
        self.assertEqual(probs['a'], [0.0, 0.0])
        self.assertEqual(probs['b'], [1.0, 1.0])
        lot = result['lots'][0]
        self.assertEqual(lot['sample_defects'], [1, 1])
        self.assertEqual(lot['prevalence'], [0.5, 0.5])

    def test_nothing_left_to_learn_selects_nobody(self):
        result = analyze(tiny_case())
        self.assertEqual(result['followup']['budget'], 1)
        self.assertEqual(result['followup']['selected_specimens'], [])


class ExactCaseTests(unittest.TestCase):
    def test_hand_computed_posterior(self):
        result = analyze(exact_case())
        probs = {row['specimen']: row['defective_probability']
                 for row in result['identification']['fruit_posteriors']}
        # 手算：w(0)=0.375, w(1)=0.125, w(2)=0.375；P(x=1|pos)=0.45/0.57
        for fruit in ('a', 'b'):
            self.assertAlmostEqual(probs[fruit][0], 0.45 / 0.57, places=5)
            self.assertAlmostEqual(probs[fruit][1], 0.45 / 0.57, places=5)
        lot = result['lots'][0]
        # 后验：theta=0:0.007, 0.25:0.142, 0.5:0.239, 0.75:0.296, 1:0.316
        self.assertEqual(lot['prevalence'], [0.25, 1.0])
        # k=0 后验 0.0132 >= 0.01 仍相容
        self.assertEqual(lot['sample_defects'], [0, 2])

    def test_compatible_configurations_match_hand_count(self):
        result = analyze(exact_case())
        configs = result['identification']['compatible_count_configurations']
        by_count = {row['counts']['L']: row['posterior'] for row in configs}
        self.assertAlmostEqual(by_count[0][1], 0.0075 / 0.57, places=5)
        self.assertAlmostEqual(by_count[1][1], 0.225 / 0.57, places=5)
        self.assertAlmostEqual(by_count[2][1], 0.3375 / 0.57, places=5)


class DesignAccountingTests(unittest.TestCase):
    def test_repeated_membership_is_not_new_sample(self):
        result = analyze(load_case('standard'))
        design = result['diagnostics']['design']
        self.assertEqual(design['north']['unique_fruits'], 12)
        self.assertEqual(design['north']['pool_memberships'], 21)
        self.assertEqual(design['south']['unique_fruits'], 12)
        self.assertEqual(design['south']['pool_memberships'], 16)
        shared = {row['specimen'] for row in result['diagnostics']['shared_specimens']}
        self.assertIn('fruit02', shared)  # fruit02 同时进入 pool0 与 pool1

    def test_stage2_pool_is_not_a_fresh_sample(self):
        result = analyze(load_case('standard'))
        north = result['diagnostics']['design']['north']
        # 追加管 followup 的两颗果实已计入唯一果实数，不额外增加样本量
        self.assertEqual(north['unique_fruits'], 12)
        text = ' '.join(result['diagnostics']['assumptions'])
        self.assertIn('不是新的随机抽样', text)


class CalibrationGapTests(unittest.TestCase):
    def test_missing_calibration_excludes_pool_without_borrowing(self):
        base = exact_case()
        base['calibration'] = [
            {'assay': 'ex', 'pool_size': 3, 'bad_count': 0, 'positive_probability': [0.02, 0.02]},
            {'assay': 'ex', 'pool_size': 3, 'bad_count': 1, 'positive_probability': [0.9, 0.9]},
            {'assay': 'ex', 'pool_size': 3, 'bad_count': 2, 'positive_probability': [0.9, 0.9]},
            {'assay': 'ex', 'pool_size': 3, 'bad_count': 3, 'positive_probability': [0.9, 0.9]},
        ]
        with_pool = analyze(base)
        without_pool = analyze({**base, 'pools': []})
        self.assertEqual(with_pool['lots'], without_pool['lots'])
        self.assertEqual(with_pool['identification']['fruit_posteriors'],
                         without_pool['identification']['fruit_posteriors'])
        self.assertEqual(len(with_pool['diagnostics']['uncalibrated_pools']), 1)

    def test_edge_case_surfaces_information_limits(self):
        result = analyze(load_case('edge'))
        diagnostics = result['diagnostics']
        self.assertTrue(diagnostics)
        self.assertEqual(len(diagnostics['uncalibrated_pools']), 5)
        self.assertIn('south', diagnostics['unassessable_lots'])
        self.assertTrue(diagnostics['untested_specimens'])
        south = next(row for row in result['lots'] if row['lot'] == 'south')
        self.assertGreater(south['prevalence'][1] - south['prevalence'][0], 0.5)

    def test_midpoint_interval_inside_envelope(self):
        result = analyze(load_case('standard'))
        for lot, row in result['diagnostics']['calibration_sensitivity'].items():
            envelope = row['envelope']
            mid = row['corner_intervals']['mid']
            self.assertLessEqual(envelope[0], mid[0])
            self.assertLessEqual(mid[1], envelope[1])
            self.assertGreaterEqual(row['width_spread_from_calibration'], 0.0)


class FollowupContractTests(unittest.TestCase):
    def test_budget_zero_selects_nobody(self):
        result = analyze(load_case('standard'), audit_budget=0)
        self.assertEqual(result['followup']['selected_specimens'], [])

    def test_selection_is_unique_within_budget_and_valid(self):
        case = load_case('standard')
        result = analyze(case, audit_budget=3)
        selected = result['followup']['selected_specimens']
        self.assertEqual(len(selected), len(set(selected)))
        self.assertLessEqual(len(selected), 3)
        valid = {row['id'] for row in case['specimens']}
        self.assertTrue(set(selected) <= valid)
        for candidate in result['followup']['candidates']:
            lo, hi = candidate['expected_variance_reduction']
            self.assertLessEqual(lo, hi)
            self.assertGreater(hi, 0.0)

    def test_default_budget_comes_from_case(self):
        result = analyze(load_case('standard'))
        self.assertEqual(result['followup']['budget'], 5)

    def test_invalid_budget_rejected(self):
        for bad in (-1, 1.5, True, '2'):
            with self.assertRaises(ValueError):
                analyze(load_case('standard'), audit_budget=bad)


class OutputContractTests(unittest.TestCase):
    def test_snapshot_and_finite_json(self):
        for name in ('standard', 'edge'):
            case = load_case(name)
            result = analyze(case)
            self.assertEqual(result['snapshot'], case['_snapshot'])
            json.dumps(result, allow_nan=False)

    def test_lots_ranges_ordered_and_bounded(self):
        case = load_case('standard')
        result = analyze(case)
        counts = {}
        for row in case['specimens']:
            counts[row['lot']] = counts.get(row['lot'], 0) + 1
        for row in result['lots']:
            if row['sample_defects'] is not None:
                self.assertLessEqual(0, row['sample_defects'][0])
                self.assertLessEqual(row['sample_defects'][0], row['sample_defects'][1])
                self.assertLessEqual(row['sample_defects'][1], counts[row['lot']])
            if row['prevalence'] is not None:
                self.assertLessEqual(0, row['prevalence'][0])
                self.assertLessEqual(row['prevalence'][0], row['prevalence'][1])
                self.assertLessEqual(row['prevalence'][1], 1)
            self.assertIn('可行', row['sample_defects_basis'])
            self.assertIn('可信', row['prevalence_basis'])

    def test_analysis_over_real_http(self):
        server = make_server()
        worker = threading.Thread(target=server.serve_forever, daemon=True)
        worker.start()
        try:
            base = 'http://127.0.0.1:' + str(server.server_port)
            result = analysis(base)
            self.assertEqual(result['snapshot'], load_case('standard')['_snapshot'])
            for field in ('lots', 'identification', 'followup', 'diagnostics'):
                self.assertIn(field, result)
            result = analysis(base, 'edge', audit_budget=2)
            self.assertLessEqual(len(result['followup']['selected_specimens']), 2)
        finally:
            server.shutdown()
            worker.join(timeout=3)
            server.server_close()


if __name__ == '__main__':
    unittest.main()
