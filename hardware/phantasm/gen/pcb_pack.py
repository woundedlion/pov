"""Pure placement packing and board reservation geometry."""
import math
from pcb_model import FAR_CONNS
from pcb_model import HUB_CONNS
from pcb_model import MOUNTING_HOLE_INSET
from pcb_model import MOUNTING_KEEPOUT_HALF_SIDE
from pcb_model import MOUNTING_RESERVE_HALF_SIDE
from pcb_model import PCB_W


def mounting_holes(L):
    """R-MECH-2: four M2.5 clearance holes, centres 3.5 mm in from each corner."""
    d = MOUNTING_HOLE_INSET
    return {"H1": (d, d), "H2": (d, PCB_W - d),
            "H3": (L - d, d), "H4": (L - d, PCB_W - d)}


def keepout_rects(L):
    """All-copper keepout square (screw head) around each mounting hole, as
    ref -> (min x, min y, max x, max y)."""
    half = MOUNTING_KEEPOUT_HALF_SIDE
    return {ref: (x - half, y - half, x + half, y + half)
            for ref, (x, y) in mounting_holes(L).items()}


def mounting_reserve_rects(L):
    """Square each mounting hole holds other footprints out of, as
    ref -> (min x, min y, max x, max y)."""
    half = MOUNTING_RESERVE_HALF_SIDE
    return {ref: (x - half, y - half, x + half, y + half)
            for ref, (x, y) in mounting_holes(L).items()}


def _boxes_overlap(a, b, tol=1e-6):
    return (a[0] < b[2] - tol and a[2] > b[0] + tol
            and a[1] < b[3] - tol and a[3] > b[1] + tol)


def keepout_clashes(place, crt_bxs, L):
    """Refs reaching into a mounting-hole reservation, as "REF/HOLE".

    `crt_bxs` are courtyard boxes, which overrun the pads: a part whose copper
    clears the screw head can still touch the hole footprint's courtyard, and
    KiCad reads that as an error-severity overlap.
    """
    reserved = mounting_reserve_rects(L)
    out = []
    for ref, (x, y, rot) in sorted(place.items()):
        mnx, mny, mxx, mxy = _rot_bb(crt_bxs[ref], rot)
        box = (x + mnx, y + mny, x + mxx, y + mxy)
        out += [f"{ref}/{hole}" for hole in sorted(reserved)
                if _boxes_overlap(box, reserved[hole])]
    return out


def outline_overflows(place, pad_bxs, L, refs):
    """Refs of `refs` whose copper reaches past the L x PCB_W outline."""
    out = []
    for ref in sorted(refs):
        x, y, rot = place[ref]
        mnx, mny, mxx, mxy = _rot_bb(pad_bxs[ref], rot)
        if x + mnx < 0 or y + mny < 0 or x + mxx > L or y + mxy > PCB_W:
            out.append(f"{ref} [{x+mnx:.2f},{y+mny:.2f}]-"
                       f"[{x+mxx:.2f},{y+mxy:.2f}]")
    return out


def _rot_bb(bb, rot):
    """Bounding box of `bb` rotated by a right angle.

    Any other angle is an error, not the unrotated box.
    """
    mnx, mny, mxx, mxy = bb
    pts = [(mnx, mny), (mxx, mny), (mxx, mxy), (mnx, mxy)]
    if rot == 90:
        pts = [(y, -x) for x, y in pts]
    elif rot == 180:
        pts = [(-x, -y) for x, y in pts]
    elif rot == 270:
        pts = [(-y, x) for x, y in pts]
    elif rot != 0:
        raise ValueError(f"placement rotation must be 0/90/180/270: {rot!r}")
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    return (min(xs), min(ys), max(xs), max(ys))


def _rotatable(ref):
    # Placement aligns the rotated bbox corner to the cell, so any origin packs
    # correctly. Rotate chip passives (R*/C*). Keep solder jumpers (JP*, custom
    # pads + clearance-outline trip DRC at 90), C_IN (a through-hole radial
    # electrolytic), diodes, fuse, bead and ICs at rot 0.
    if ref.startswith("JP") or ref == "C_IN":
        return False
    return ref[0] in "RC"


def _column_height(refs, bxs, gap):
    """Y extent of a rot-0 column of `refs` separated by `gap`."""
    heights = [_rot_bb(bxs[r], 0)[3] - _rot_bb(bxs[r], 0)[1] for r in refs]
    return sum(heights) + gap * max(0, len(heights) - 1)


def _stack(refs, bxs, x0, edge, gap, place, y0=0.0):
    """Place `refs` as a column stacked along Y from y0 at left edge x0 (rot 0).
    Returns (column right-edge x, total stacked height)."""
    y = y0
    w = 0.0
    for ref in refs:
        rb = _rot_bb(bxs[ref], 0)
        place[ref] = (round(edge + x0 - rb[0], 3), round(edge + y - rb[1], 3), 0)
        y += (rb[3] - rb[1]) + gap
        w = max(w, rb[2] - rb[0])
    return x0 + w, y - gap - y0


def _skyline_free_x(sky, yb, h):            # leftmost x a [yb,yb+h] slot can start at
    x = 0.0
    for s0, s1, sx in sky:
        if s1 > yb + 1e-9 and s0 < yb + h - 1e-9:
            x = max(x, sx)
    return x


def _skyline_best_position(sky, usable, h):             # (x, yb) bottom-left slot for a band of height h
    cands = {0.0}
    for s0, s1, _ in sky:
        cands.add(s0); cands.add(s1 - h)
    best = None
    for yb in sorted(cands):
        if yb < -1e-6 or yb + h > usable + 1e-6:
            continue
        r = (_skyline_free_x(sky, yb, h), yb)
        if best is None or r < best:
            best = r
    return best


def _skyline_reserve(sky, yb, w, h, x):    # raise the skyline over [yb,yb+h] to x+w
    out = []
    for s0, s1, sx in sky:
        if s1 <= yb or s0 >= yb + h:
            out.append((s0, s1, sx)); continue
        if s0 < yb:
            out.append((s0, yb, sx))
        out.append((max(s0, yb), min(s1, yb + h), x + w))
        if s1 > yb + h:
            out.append((yb + h, s1, sx))
    out.sort()
    merged = []
    for seg in out:
        if merged and abs(merged[-1][2] - seg[2]) < 1e-6 \
                and abs(merged[-1][1] - seg[0]) < 1e-6:
            merged[-1] = (merged[-1][0], seg[1], seg[2])
        else:
            merged.append(seg)
    return merged


def _pack_interior(bxs, pinned, edge, gap, usable, sky, place):
    interior = sorted((r for r in bxs if r not in pinned),
                      key=lambda r: -max(bxs[r][2] - bxs[r][0], bxs[r][3] - bxs[r][1]))
    for ref in interior:
        bb = bxs[ref]
        rots = (0, 90) if _rotatable(ref) else (0,)
        choice = None
        for rot in rots:
            rb = _rot_bb(bb, rot)
            wg, hg = rb[2] - rb[0] + gap, rb[3] - rb[1] + gap
            if hg > usable + 1e-6:
                continue
            bp = _skyline_best_position(sky, usable, hg)
            if bp is None:
                continue
            x, yb = bp
            if choice is None or x < choice[0] - 1e-9:
                choice = (x, yb, rot, rb, wg, hg)
        if choice is None:        # too tall even rotated — force rot 0 at the end
            rb = _rot_bb(bb, 0)
            wg, hg = rb[2] - rb[0] + gap, rb[3] - rb[1] + gap
            x = max(s[2] for s in sky); yb = 0.0
            choice = (x, yb, 0, rb, wg, hg)
        x, yb, rot, rb, wg, hg = choice
        place[ref] = (round(edge + x - rb[0], 3), round(edge + yb - rb[1], 3), rot)
        sky = _skyline_reserve(sky, yb, wg, hg, x)
    return sky


def pack(bxs, width, edge=1.0, gap=1.2):
    """Connectors pinned to the ends + skyline (bottom-left-fill) interior pack.
    Board WIDTH is fixed; minimise LENGTH. Returns ({ref:(x,y,rot)}, length mm).

    The hub-end mounting-hole keepouts sit at a fixed x, so they seed the skyline;
    the far-end pair follows the length, so it is reserved as tail length past the
    rightmost part instead. Either way no part lands under a screw head."""
    usable = width - 2 * edge
    hub = [r for r in HUB_CONNS if r in bxs]
    far = [r for r in FAR_CONNS if r in bxs]
    pinned = set(hub) | set(far)
    place = {}
    sky = [(0.0, usable, 0.0)]    # skyline over the width: (y0, y1, x_right)

    # hub-end screw heads: block their squares before anything is placed
    half = MOUNTING_RESERVE_HALF_SIDE
    for cy in (MOUNTING_HOLE_INSET, width - MOUNTING_HOLE_INSET):
        yb = max(0.0, cy - half - edge)
        yt = min(usable, cy + half - edge)
        if yt > yb:
            sky = _skyline_reserve(sky, yb, 2 * half, yt - yb, MOUNTING_HOLE_INSET - half - edge)

    # hub connectors: column at the left edge; reserve it so the interior packs right
    if hub:
        hubh = _column_height(hub, bxs, gap) + gap
        x0, y0 = _skyline_best_position(sky, usable, hubh) or (_skyline_free_x(sky, 0.0, hubh), 0.0)
        hubw, _ = _stack(hub, bxs, x0, edge, gap, place, y0)
        sky = _skyline_reserve(sky, y0, hubw + gap, hubh, 0.0)

    sky = _pack_interior(bxs, pinned, edge, gap, usable, sky, place)

    # Strip and sync connectors occupy the far edge; power occupies the hub.
    right = max(s[2] for s in sky)
    if far:
        sumh = _column_height(far, bxs, 0.0)
        fgap = gap
        if len(far) > 1:
            fgap = max(0.6, min(gap, (usable - sumh) / (len(far) - 1)))
        right, _ = _stack(far, bxs, right + gap, edge, fgap, place)

    # length: enough tail for the far-end screw heads to clear every part
    extent = max((x + _rot_bb(bxs[ref], rot)[2] for ref, (x, _, rot) in place.items()),
                 default=edge + right)
    tail = max(edge, MOUNTING_HOLE_INSET + MOUNTING_RESERVE_HALF_SIDE)
    return place, math.ceil((extent + tail) * 100 - 1e-9) / 100


def unplaced_layout(bxs, L, width, margin=2.0, gap=2.0):
    """Stage every footprint in a grid BELOW the board outline (rot 0), so an
    autoplacer (Quilter) sees them as unplaced. Returns {ref:(x,y,rot)}."""
    place = {}
    x = margin
    y = width + 6.0           # start clear of the outline (which is 0..width in Y)
    rowh = 0.0
    wrap = max(L, 40.0)
    for ref in sorted(bxs):
        mnx, mny, mxx, mxy = bxs[ref]
        w, h = mxx - mnx, mxy - mny
        if x + w > wrap - margin:
            x = margin; y += rowh + gap; rowh = 0.0
        place[ref] = (round(x - mnx, 3), round(y - mny, 3), 0)
        x += w + gap; rowh = max(rowh, h)
    return place
