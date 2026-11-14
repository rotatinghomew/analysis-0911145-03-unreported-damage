from collections import Counter, defaultdict
from app.outreach import contact_summary
from app.store import unique


def capture_patterns(case):
    population = unique(case['population'], 'id')
    records = unique(case['lists'], 'record')
    linked = {}
    by_order = defaultdict(set)
    for link in case['confirmed_links']:
        if link['record'] not in records or link['order'] not in population:
            raise ValueError('link outside observed records or population')
        previous = linked.get(link['record'])
        if previous is not None and previous != link['order']:
            raise ValueError('conflicting confirmed identity')
        linked[link['record']] = link['order']
        by_order[link['order']].add(records[link['record']]['source'])
    patterns = Counter()
    for order, sources in by_order.items():
        patterns[(population[order]['stratum'], tuple(sorted(sources)))] += 1
    return [{'stratum': stratum, 'sources': list(sources), 'count': count}
            for (stratum, sources), count in sorted(patterns.items())]


def preview(case):
    linked_records = {link['record'] for link in case['confirmed_links']}
    patterns = capture_patterns(case)
    return {'snapshot': case['_snapshot'], 'population': len(case['population']),
            'confirmed_patterns': patterns,
            'unlinked_records': sum(row['record'] not in linked_records for row in case['lists']),
            'outreach': contact_summary(case['outreach']),
            'observed_only': True}
