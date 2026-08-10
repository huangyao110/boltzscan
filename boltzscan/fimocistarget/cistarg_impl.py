"""FIMO hit de-duplication helpers used by the promoter/model pipeline."""
import numpy as np


def filter_fimo_by_seq_and_overlap(df, overlap_threshold=0.0, extend_range=5):
    """Deduplicate FIMO hits while preserving distinct motif core sites.

    MEME Suite FIMO reports 1-based, inclusive coordinates. ``real_start``
    and ``real_stop`` are the corresponding 0-based, half-open DNA-window
    coordinates used for Python slicing.
    """
    if extend_range < 0:
        raise ValueError("extend_range must be non-negative")
    if not 0.0 <= overlap_threshold <= 1.0:
        raise ValueError("overlap_threshold must be between 0 and 1")

    print("正在初始化并重置索引...")
    data = df.copy().reset_index(drop=True)

    # FIMO coordinates are 1-based and inclusive, unlike Python slices.
    core_starts = data[['start', 'stop']].min(axis=1).astype(int).values - 1
    core_stops = data[['start', 'stop']].max(axis=1).astype(int).values
    starts = np.maximum(core_starts - extend_range, 0)
    stops = np.minimum(core_stops + extend_range, data['promoter_seq'].str.len().values)
    full_seqs = data['promoter_seq'].values

    data['extracted_motif_seq'] = [
        sequence[start:stop]
        for sequence, start, stop in zip(full_seqs, starts, stops)
    ]
    data['real_start'] = starts
    data['real_stop'] = stops
    data['length'] = stops - starts
    data['_core_start'] = core_starts
    data['_core_stop'] = core_stops

    # p-values, unlike raw PWM scores, are comparable across motifs. The
    # remaining keys make equal-score choices reproducible.
    data = data.sort_values(
        by=['pvalue', 'score', 'sequence_name', 'start', 'stop', 'strand', 'motif_id'],
        ascending=[True, False, True, True, True, True, True],
        kind='mergesort',
    )
    print(f"原始数据行数: {len(df)}")

    # Same promoter, TF, and DNA bait is one structural/Y1H input.
    data['tf_name'] = data['tf_name'].fillna('Unknown_TF')
    before_seq_filter = len(data)
    data = data.drop_duplicates(
        subset=['sequence_name', 'tf_name', 'extracted_motif_seq'],
        keep='first',
    )
    print(f"剔除相同序列内容后: {len(data)} (减少 {before_seq_filter - len(data)})")

    # Expanded flanks define the modeled DNA window, not whether two FIMO
    # cores are separate binding sites.
    # ponytail: O(n²) per group — fine for typical FIMO sets; for tens of millions of
    # rows, upgrade to interval tree / sweep-line.
    print("正在进行物理重叠去重 (此步骤处理千万级数据较慢，请稍候)...")
    keep_indices = []
    grouped = data.groupby(['sequence_name', 'tf_name'], sort=False)

    for count, (_, group) in enumerate(grouped, start=1):
        g_data = group[['_core_start', '_core_stop']].values
        g_indices = group.index.values
        accepted_intervals = []

        for i, (curr_start, curr_stop) in enumerate(g_data):
            curr_len = curr_stop - curr_start
            is_redundant = False
            for acc_start, acc_stop in accepted_intervals:
                overlap_len = min(curr_stop, acc_stop) - max(curr_start, acc_start)
                if overlap_len > 0 and overlap_len / curr_len > overlap_threshold:
                    is_redundant = True
                    break
            if not is_redundant:
                accepted_intervals.append((curr_start, curr_stop))
                keep_indices.append(g_indices[i])

        if count % 10000 == 0:
            print(f"已处理 {count}/{len(grouped)} 个 Group...", end='\r')

    print("\n物理去重完成。")
    result = data.loc[keep_indices].sort_values(
        by=['sequence_name', 'start', 'stop', 'motif_id'], kind='mergesort'
    ).drop(columns=['_core_start', '_core_stop'])
    print(f"最终结果行数: {len(result)}")
    return result
