"""Finding 2 capture analysis; writes only to an explicit analysis directory."""
import argparse
import csv
import hashlib
import importlib.util
import json
import re
import statistics
from collections import defaultdict
from pathlib import Path

PERIOD = 62500
SHAPE = re.compile(r'^(Spawning|Built) Shape: (\S+) \(V=(\d+), E=(\d+), F=(\d+), I=(\d+)\)')
ISR = re.compile(r'^(isr_\w+)\s+n=(\d+) cyc min/avg/max=(\d+)/(\d+)/(\d+) ns min/avg/max=(\d+)/(\d+)/(\d+) total=(\d+) us cpu=([\d.]+)%')


def stats(rows):
    result = {'frames': len(rows)}
    for field in ('render_us', 'wall_us'):
        values = sorted(r[field] for r in rows)
        result[field] = dict(min=min(values), mean=statistics.mean(values),
                             median=statistics.median(values), max=max(values),
                             p95=values[min(len(values)-1, int(.95*len(values)))],
                             p99=values[min(len(values)-1, int(.99*len(values)))])
    result['peak_frames'] = [r['frame'] for r in rows if r['render_us'] == result['render_us']['max']]
    result['spilled'] = sum(r['render_us'] > PERIOD for r in rows)
    result['spill_pct'] = result['spilled'] / len(rows) * 100
    return result


def analyze(config, parser, expected):
    path = Path(config['log'])
    raw = path.read_text(encoding='utf-8-sig', errors='replace')
    windows, effect, _ = parser.parse_capture(path)
    errors = []
    if effect != 'IslamicStars':
        errors.append(f'Wrong effect: {effect}')
    rows, raw_windows, spawns = {}, {}, []
    shape, geometry, current, last_frame = None, None, None, 0
    for line in raw.splitlines():
        m = SHAPE.match(line)
        if m:
            kind, name, *sizes = m.groups()
            if kind == 'Spawning':
                shape, geometry = name, None
                spawns.append({'name': name, 'next_frame': last_frame + 1})
            else:
                if shape != name:
                    errors.append(f'Built shape mismatch {name}/{shape}')
                geometry = dict(zip(('V', 'E', 'F', 'I'), map(int, sizes)))
        m = parser.FRAME_RE.match(line)
        if m:
            n, wall, render = map(int, m.groups())
            if n in rows or n <= last_frame:
                errors.append(f'Repeated/nonmonotonic frame {n}')
            rows[n] = dict(frame=n, wall_us=wall, render_us=render, shape=shape, geometry=geometry)
            last_frame = n
        m = parser.HEADER_RE.match(line)
        if m:
            current = {'counters': [], 'isr': {}}
            raw_windows[int(m.group(4))] = current
        m = parser.COUNTER_RE.match(line)
        if m and current is not None:
            indent, label, us, pct, calls, cyc, tags = m.groups()
            current['counters'].append(dict(label=label, depth=len(indent)//2, us=int(us),
                                           parent_pct=int(pct), calls=int(calls), cyc=int(cyc), tags=tags.strip()))
        m = ISR.match(line)
        if m and current is not None:
            label, *vals = m.groups()
            current['isr'][label] = dict(zip(('n', 'cyc_min', 'cyc_avg', 'cyc_max', 'ns_min', 'ns_avg', 'ns_max', 'total_us'), map(int, vals[:-1])))
    names = [s['name'] for s in spawns]
    if len(set(names)) != expected:
        errors.append(f'Expected {expected} unique shapes, found {len(set(names))}')
    if len(names) <= expected or names[expected] != names[0] or len(set(names[:expected])) != expected:
        errors.append('No complete shape cycle followed by wrap to first shape')
    complete, enriched, holds = [], [], defaultdict(list)
    previous = 0
    for w in windows:
        if w.f_start != previous + 1:
            errors.append(f'Window gap before {w.f_start}')
        previous = w.f_end
        if (w.w, w.h) != (288, 144) or w.effect != effect:
            errors.append(f'Window dimensions/effect mismatch {w.f_start}')
        ns = [r[0] for r in w.frame_rows]
        if ns != list(range(w.f_start, w.f_end + 1)):
            errors.append(f'Incomplete/noncontiguous closed window {w.f_start}')
            continue
        wr = [rows[n] for n in ns]
        complete += wr
        detail = raw_windows[w.f_start]
        cs = detail['counters']
        def count(label):
            return sum(c['calls'] for c in cs if c['label'] == label)
        def cycles(label):
            return sum(c['cyc'] for c in cs if c['label'] == label)
        if count('frame') != w.frames:
            errors.append(f'Root call count mismatch {w.f_start}')
        if not any(c['label'].endswith('buffer_wait') for c in cs):
            errors.append(f'Missing sync wait scope {w.f_start}')
        root_error = abs(cycles('frame') / 600 - w.wall[3]) / w.wall[3] * 1e6
        detail.update(start=w.f_start, end=w.f_end, frames=w.frames, window_us=w.window_us,
                      metrics=stats(wr), wall_header=w.wall, root_wall_error_ppm=root_error,
                      mesh_scan_ms=cycles('is_mesh_scan') / 600000 / w.frames,
                      filter_calls=count('filter_blend'), filter_cycles=cycles('filter_blend'))
        own = {r['shape'] for r in wr}
        geo = wr[0]['geometry']
        clean = (w.f_start > 1 and len(own) == 1 and geo is not None
                 and all(r['geometry'] == geo for r in wr)
                 and count('is_build_draw') == 0 and count('is_draw_shape') == w.frames
                 and count('scan_mesh_raster') == w.frames * geo['F'])
        detail['clean_hold'] = clean
        detail['build'] = count('is_build_draw') > 0
        if clean:
            detail['shape'], detail['geometry'] = wr[0]['shape'], geo
            holds[wr[0]['shape']].append(detail)
        enriched.append(detail)
    richest = max(enriched, key=lambda w: w['mesh_scan_ms'])
    if richest['root_wall_error_ppm'] > 5:
        errors.append(f"Richest window root/wall discrepancy {richest['root_wall_error_ppm']:.2f} ppm")
    runtime = [r for r in complete if r['frame'] != 1]
    if not runtime:
        raise ValueError(f'No runtime data: {path}')
    by_shape = {name: stats([r for r in runtime if r['shape'] == name]) for name in dict.fromkeys(names)}
    for data in by_shape.values():
        data['bucket'] = 'green' if data['spilled'] == 0 else ('yellow' if data['spill_pct'] < 25 else 'red')
    after_setup = [w for w in enriched if w['start'] > 1]
    isr = {}
    elapsed = sum(w['window_us'] for w in after_setup)
    for label in ('isr_wake', 'isr_pack', 'isr_dma_submit'):
        items = [w['isr'][label] for w in after_setup if label in w['isr']]
        if len(items) != len(after_setup):
            errors.append(f'Missing {label} rows')
        if items:
            n = sum(i['n'] for i in items)
            total = sum(i['total_us'] for i in items)
            isr[label] = dict(n=n, total_us=total, elapsed_us=elapsed, cpu_pct=100*total/elapsed,
                              rate_hz=n*1e6/elapsed, mean_ns=total*1000/n,
                              min_ns=min(i['ns_min'] for i in items), max_ns=max(i['ns_max'] for i in items))
    provenance = {}
    sidecar = path.with_name(path.name.removesuffix('.log.txt') + '.provenance.txt') if path.name.endswith('.log.txt') else path.with_suffix('.provenance')
    provenance_text = raw + ('\n' + sidecar.read_text(errors='replace') if sidecar.exists() else '')
    for line in provenance_text.splitlines():
        line = line.removeprefix('profile provenance: ')
        if re.match(r'^(source_sha|compiler|.*_sha256|artifact_\w+)=', line):
            key, val = line.split('=', 1)
            provenance[key] = val
    if provenance.get('source_sha') and not provenance['source_sha'].startswith(config['source_sha']):
        errors.append('Source provenance does not match expected source SHA')
    artifacts = {}
    for kind in ('profile', 'phantasm'):
        key = 'artifact_' + kind + '_elf'
        if key in provenance:
            artifact = Path(config['tree']) / provenance[key]
            actual = hashlib.sha256(artifact.read_bytes()).hexdigest() if artifact.exists() else None
            artifacts[kind] = dict(path=str(artifact), sha256=actual,
                                   verified=actual == provenance.get(kind+'_elf_sha256'))
            if not artifacts[kind]['verified']:
                errors.append(f'{kind} ELF hash verification failed')
    total_filter = sum(w['filter_calls'] for w in after_setup)
    total_filter_cyc = sum(w['filter_cycles'] for w in after_setup)
    wrap_frame = spawns[expected]['next_frame'] if len(spawns) > expected else None
    result = dict(config=config, log_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  provenance=provenance, artifacts=artifacts, validation_errors=errors, valid=not errors,
                  setup_frame=rows.get(1), runtime=stats(runtime), shape_order=names,
                  first_cycle=stats([r for r in runtime if r['frame'] < wrap_frame]) if wrap_frame else None,
                  first_wrap_frame=wrap_frame,
                  spawns=spawns, by_shape=by_shape, rows=runtime, windows=enriched,
                  clean_holds=dict(holds), isr=isr, trailing_unclosed_frames=len(rows)-len(complete),
                  perpixel=dict(filter_calls=total_filter, filter_cycles=total_filter_cyc,
                                cycles_per_blend=total_filter_cyc/total_filter if total_filter else None))
    return result


def compare(before, after):
    a, b = ({r['frame']: r for r in run['rows']} for run in (before, after))
    common = sorted(a.keys() & b.keys())
    mismatches = [n for n in common if (a[n]['shape'], a[n]['geometry']) != (b[n]['shape'], b[n]['geometry'])]
    paired = [dict(frame=n, shape=a[n]['shape'], before_render_us=a[n]['render_us'],
                   after_render_us=b[n]['render_us'], delta_render_us=b[n]['render_us']-a[n]['render_us'],
                   before_wall_us=a[n]['wall_us'], after_wall_us=b[n]['wall_us']) for n in common]
    deltas = [r['delta_render_us'] for r in paired]
    per_shape = {}
    for shape in dict.fromkeys(r['shape'] for r in paired):
        ns = [r['frame'] for r in paired if r['shape'] == shape]
        am, bm = stats([a[n] for n in ns]), stats([b[n] for n in ns])
        per_shape[shape] = dict(before=am, after=bm,
                               mean_delta_us=bm['render_us']['mean']-am['render_us']['mean'])
    cycle_end = min(before['first_wrap_frame'], after['first_wrap_frame'])
    cycle_ns = [n for n in common if n < cycle_end]
    return dict(matched_frames=len(common), schedule_mismatch_frames=mismatches, by_shape=per_shape,
                first_cycle_before=stats([a[n] for n in cycle_ns]),
                first_cycle_after=stats([b[n] for n in cycle_ns]),
                before=stats([a[n] for n in common]), after=stats([b[n] for n in common]),
                mean_delta_us=statistics.mean(deltas), median_delta_us=statistics.median(deltas),
                min_delta_us=min(deltas), max_delta_us=max(deltas),
                positive_delta_frames=sum(d > 0 for d in deltas), rows=paired)


def report(run, title):
    rt = run['runtime']
    lines = [f'# {title}', '', 'Analysis draft; candidate is unlanded. Runtime excludes setup frame 1.', '',
             '## Setup', '', f"Source: `{run['config']['source_sha']}`; log: `{run['config']['log']}`.",
             f"Capture SHA256: `{run['log_sha256']}`.", '',
             'Teensy 4.0, 600 MHz, segmented 288×144, four segments, 480 RPM, TS=4.',
             'TS=4 changes build/ripple sampling as well as hold duration; this is a matched TS=4 workload.', '',
             '## Frame cadence', '',
             f"{rt['frames']} runtime frames; mean render {rt['render_us']['mean']/1000:.3f} ms; "
             f"peak {rt['render_us']['max']/1000:.3f} ms; spilled {rt['spilled']}/{rt['frames']} ({rt['spill_pct']:.3f}%).",
             f"Wall min/mean/max: {rt['wall_us']['min']/1000:.3f}/{rt['wall_us']['mean']/1000:.3f}/{rt['wall_us']['max']/1000:.3f} ms.", '',
             '## Phase readouts', '']
    windows = [w for w in run['windows'] if w['start'] > 1]
    candidates = [('Build', [w for w in windows if w['build']]),
                  ('Finished/ripple', [w for w in windows if w['clean_hold']]),
                  ('Peak runtime window', [w for w in windows if any(w['start'] <= n <= w['end'] for n in rt['peak_frames'])])]
    seen = set()
    for name, group in candidates:
        if not group:
            continue
        w = max(group, key=lambda x: x['metrics']['render_us']['max'])
        if w['start'] in seen:
            lines += [f"The {name.lower()} is the already shown frames {w['start']}–{w['end']}; its tree is not repeated.", '']
            continue
        seen.add(w['start'])
        root = next(c['cyc'] for c in w['counters'] if c['label'] == 'frame')
        lines += [f"### {name}: frames {w['start']}–{w['end']}", '', '```text']
        for i, c in enumerate(w['counters']):
            label = '  '*c['depth'] + c['label']
            us = c['cyc']/600/w['frames']
            cyc = c['cyc']/w['frames']
            time = f'{us/1000:.2f}ms' if us >= 1000 else f'{us:.2f}us'
            cycles = f'{cyc/1e6:.2f}M' if cyc >= 1e6 else (f'{cyc/1000:.2f}k' if cyc >= 1000 else f'{cyc:.0f}c')
            leaf = i+1 == len(w['counters']) or w['counters'][i+1]['depth'] <= c['depth']
            suffix = f" x{c['calls']/w['frames']:.1f} {c['cyc']/600/c['calls']:.2f}us" if leaf and c['calls'] else ''
            lines.append(f"{label:<28} {time:>8} {cycles:>7} {100*c['cyc']/root:5.1f}%{suffix}")
        lines += ['```', '']
        tagged = sorted({c['label'] + ' (' + c['tags'] + ')' for c in w['counters'] if c['tags']})
        if tagged:
            lines += ['Shared/nonexclusive counters: ' + ', '.join(tagged) + '.', '']
        m = w['metrics']['wall_us']
        lines += [f"Wall min/mean/max: {m['min']/1000:.3f}/{m['mean']/1000:.3f}/{m['max']/1000:.3f} ms. Percentages use the root frame counter; tagged shared or mixed-parent nodes are inclusive costs.", '']
    lines += ['### Per-preset table', '', 'Ranked by worst clean-hold is_mesh_scan window. Window counts use the strict finished-geometry predicate below. Source indices follow the first 23 spawn markers; the cycle wraps back to its first entry. Window fps reflects measured wall time and boundary jitter; zero-spill shapes retain the observed 16 fps cadence.', '', '| # | Shape | V/E/F/I | Windows | Scan ms/frame | Blends/frame | Render ms/frame | Window fps |', '|---:|---|---|---:|---:|---:|---:|---:|']
    for name, holds in sorted(run['clean_holds'].items(), key=lambda kv: max(w['mesh_scan_ms'] for w in kv[1]), reverse=True):
        w = max(holds, key=lambda x: x['mesh_scan_ms'])
        geo = '/'.join(str(w['geometry'][k]) for k in ('V','E','F','I'))
        lines.append(f"| {run['shape_order'].index(name)} | {name} | {geo} | {len(holds)} | {w['mesh_scan_ms']:.3f} | {w['filter_calls']/w['frames']:.1f} | {w['metrics']['render_us']['mean']/1000:.3f} | {1e6/w['metrics']['wall_us']['mean']:.2f} |")
    missing = set(run['by_shape']) - set(run['clean_holds'])
    lines += ['', f"No strict clean-hold window: {', '.join(sorted(missing)) or 'none'}.", '', '### Per-pixel figures', '',
              f"Filter blend: {run['perpixel']['cycles_per_blend']} cycles/call, including instrumentation.", '',
              '## Column ISR', '', '| Scope | Calls | Rate Hz | Min ns | Mean ns | Max ns | CPU % |', '|---|---:|---:|---:|---:|---:|---:|']
    for name, v in run['isr'].items():
        lines.append(f"| {name} | {v['n']} | {v['rate_hz']:.2f} | {v['min_ns']} | {v['mean_ns']:.1f} | {v['max_ns']} | {v['cpu_pct']:.3f} |")
    lines += ['', '## Caveats', '', 'Scope times include interrupts. Wake includes pack and DMA submit; these are not additive. Counter trees preserve duplicate labels and mixed parents; shared scopes do not support exclusive phase attribution. ISR means use rounded window totals and header elapsed time. Asynchronous 600-byte SPI wire transfer is 200 µs at 24 MHz.', '',
              'Clean holds require one finished geometry, no build draw, one shape draw per frame and exactly F raster calls per frame. Geometry comes from Built Shape, not seed geometry from Spawning Shape.', '',
              f"Trailing unclosed frames excluded: {run['trailing_unclosed_frames']}. Validation errors: {run['validation_errors']}."]
    return '\n'.join(lines)+'\n'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--manifest', default=str(Path(__file__).with_name('capture_manifest.json')))
    ap.add_argument('--out', required=True)
    ap.add_argument('--runs', nargs='*')
    ap.add_argument('--drafts', action='store_true')
    args = ap.parse_args()
    manifest = json.loads(Path(args.manifest).read_text(encoding='utf-8-sig'))
    spec = importlib.util.spec_from_file_location('profile_parser', manifest['parser'])
    parser = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(parser)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    runs = {key: analyze(cfg, parser, manifest['expected_shapes']) for key, cfg in manifest['runs'].items() if not args.runs or key in args.runs}
    comparisons = {}
    for config in ('ship', 'o3'):
        if 'before_'+config in runs and 'after_'+config in runs:
            comparisons[config] = compare(runs['before_'+config], runs['after_'+config])
            with (out/f'matched-{config}.csv').open('w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=list(comparisons[config]['rows'][0]))
                writer.writeheader()
                writer.writerows(comparisons[config]['rows'])
    (out/'analysis.json').write_text(json.dumps(dict(manifest=manifest, runs=runs, comparisons=comparisons), indent=2))
    if args.drafts:
        if len(runs) != 4 or any(not r['valid'] for r in runs.values()) or any(c['schedule_mismatch_frames'] for c in comparisons.values()):
            raise SystemExit('Drafts require all four valid matched captures')
        for key, run in runs.items():
            (out/f'{key}-draft.md').write_text(report(run, key), encoding='utf-8')
    print(json.dumps({key: dict(valid=r['valid'], errors=r['validation_errors'], runtime=r['runtime'], holds=len(r['clean_holds'])) for key, r in runs.items()}, indent=2))


if __name__ == '__main__':
    main()
