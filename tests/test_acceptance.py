import json
import os
import threading
import unittest
from app.client import analysis
from app.server import make_server
from app.study_store import load_case


@unittest.skipUnless(os.environ.get('RUN_ACCEPTANCE') == '1', '分析消费接口尚未接入')
class AnalysisContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.server = make_server()
        cls.worker = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.worker.start()
        cls.base = 'http://127.0.0.1:' + str(cls.server.server_port)

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.worker.join(timeout=3)
        cls.server.server_close()

    def test_consumer_contract_and_provenance(self):
        result = analysis(self.base)
        self.assertEqual(result['snapshot'], load_case('standard')['_snapshot'])
        for field in ['lots', 'identification', 'followup', 'diagnostics']:
            self.assertIn(field, result)
        json.dumps(result, allow_nan=False)

    def test_information_limits_are_visible(self):
        result = analysis(self.base, 'edge')
        self.assertEqual(result['snapshot'], load_case('edge')['_snapshot'])
        self.assertTrue(result['diagnostics'])

    def test_sample_and_population_bounds_use_their_own_units(self):
        case = load_case('standard')
        result = analysis(self.base)
        for row in result['lots']:
            count = sum(fruit['lot'] == row['lot'] for fruit in case['specimens'])
            if row['sample_defects'] is not None:
                self.assertLessEqual(0, row['sample_defects'][0])
                self.assertLessEqual(row['sample_defects'][0], row['sample_defects'][1])
                self.assertLessEqual(row['sample_defects'][1], count)
            if row['prevalence'] is not None:
                self.assertLessEqual(0, row['prevalence'][0])
                self.assertLessEqual(row['prevalence'][0], row['prevalence'][1])
                self.assertLessEqual(row['prevalence'][1], 1)

    def test_followup_does_not_reuse_the_same_fruit_or_exceed_budget(self):
        result = analysis(self.base, audit_budget=2)
        selected = result['followup']['selected_specimens']
        self.assertEqual(len(selected), len(set(selected)))
        self.assertLessEqual(len(selected), 2)
        self.assertTrue(set(selected) <= {fruit['id'] for fruit in load_case('standard')['specimens']})
