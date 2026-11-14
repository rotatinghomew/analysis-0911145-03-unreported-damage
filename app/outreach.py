def contact_summary(rows):
    sampled = 0
    connected = 0
    reported = 0
    for row in rows:
        probability = row['inclusion_probability']
        if not 0 < probability <= 1:
            raise ValueError('invalid inclusion probability')
        sampled += int(row['sampled'])
        if row['connected']:
            connected += 1
            reported += int(row['damage_reported'])
        elif row['damage_reported'] is not None:
            raise ValueError('unanswered contact cannot have a label')
    return {'sampled': sampled, 'connected': connected, 'reported': reported}
