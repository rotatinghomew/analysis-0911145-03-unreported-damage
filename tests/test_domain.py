import copy
import unittest
from app.claims import capture_patterns
from app.outreach import contact_summary
from app.store import load_case


class CoverageTests(unittest.TestCase):
    def test_duplicate_list_entry_does_not_duplicate_person(self):
        case = load_case('standard')
        previous = capture_patterns(case)
        case['confirmed_links'].append(dict(case['confirmed_links'][0]))
        self.assertEqual(capture_patterns(case), previous)

    def test_unanswered_callback_is_not_a_negative_label(self):
        summary = contact_summary(load_case('edge')['outreach'])
        self.assertEqual(summary['connected'], 0)
        self.assertGreater(summary['sampled'], 0)

    def test_conflicting_identity_is_visible(self):
        case = load_case('standard')
        case['confirmed_links'].append({**case['confirmed_links'][0], 'order': 'p239'})
        with self.assertRaises(ValueError):
            capture_patterns(case)
