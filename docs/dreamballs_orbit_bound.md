# DreamBalls orbit bounds

DreamBalls evaluates a tangent-plane orbit with Bhaskara sine and cosine,
then normalizes it onto the unit sphere. Its motion planner uses a bound on
the normalized path's second derivative to choose endpoint chords or
three-capture parabolas. The bounds include the approximation's variation in
radius and angular rate.

## Normalized path

Let `b`, `u`, and `v` be an orthonormal frame, `r` the tangent radius in
`[0,1]`, and `w` the phase change over a frame. In coordinates along that
frame, write `q = (r c, r s, 1)` and `p = q / N`, where
`N² = 1 + r²(c²+s²)`. The woven head uses a fixed parallel transport, which
preserves the tangent frame's inner products and therefore the same bounds.

On a quarter turn, set `z = 2 theta / pi` in `[0,1]` and
`f(z) = 4 z (2-z) / (z²-2z+5)`. Then
`g = (c,s) = (f(1-z), f(z))`. The remaining quarters follow by signed
coordinate swaps. Derivatives below are with respect to `theta`.

The exact rational certificate below proves these conservative inequalities:

| Quantity | Bound |
|---|---:|
| `g·g` | at least `0.996` |
| norm of `g'` | at most `1.02` |
| norm of `g''` | at most `1.026` |
| absolute `h = g·g'` | at most `0.005` |
| absolute `k = g'·g' + g·g''` | at most `0.025` |

Since `N² >= 0.996(1+r²)`, put `rho = r / sqrt(1+r²)`. Projection onto the
plane perpendicular to `p` gives `|p'| <= r |g'| / N`, so
`|p'| <= 1.02 rho / sqrt(0.996) < 1.025 rho`.

Differentiating the normalization twice gives

```text
p'' = q''/N - (2 r² h q' + r² k q)/N³ + 3 r⁴ h² q/N⁵.
```

Consequently,

```text
|p''| <= r |g''|/N + 2 r³ |h| |g'|/N³
         + r² |k|/N² + 3 r⁴ h²/N⁴.
```

Set `a = r/N` and `A = a/rho`. For `r <= 1`, the preceding norm bound
implies `A <= 1.002`, `a² <= 0.502`, and `a <= 0.708`. Dividing the four
terms by `rho` bounds their sum by

```text
A (|g''| + 2 a² |h| |g'| + a |k| + 3 a³ h²)
<= 1.002 (1.026 + 2*0.502*0.005*1.02 + 0.708*0.025
          + 3*0.708*0.502*0.005²)
< 1.051 < 1.06.
```

Frame derivatives multiply these phase bounds by `|w|` and `w²`.
`DreamBalls::orbit_bounds` therefore supplies motion `1.025 rho |w|` and
second derivative `1.06 rho w²`. `Spacetime::PathBound::chord_error` divides
the second derivative by eight. The chord bound is about 7.55 to 10.67 times
tighter than `r w²` over the supported nonzero radii.

Bhaskara sine is continuously differentiable across its half-period folds,
but its second derivative jumps there. The first derivative remains
Lipschitz, so the chord bound holds across folds. A moving path has no finite
global third-derivative bound; its third bound is infinity. Stationary paths
have zero bounds. The mesh planner uses the second bound to choose a chord
or capture a midpoint; it does not consume the third bound.

The certificate treats the real rational kernel. Native C++ tests separately
check the float evaluator, normalized derivative formulas, stationary paths,
fold crossings, and midpoint work. Their capture comparisons allow ordinary
float rounding; the analytic derivative checks use double precision without
an additive error allowance.

## Reproducible exact certificate

Run this standalone Python program. It clears positive rational denominators
and proves the resulting polynomials nonnegative on `[0,1]`: all Bernstein
coefficients are nonnegative on each dyadic subinterval. Every operation uses
`Fraction`. The rational `200000/314159` exceeds `2/pi`, so upper derivative
bounds also cover the exact scale. The five checks need at most four
subintervals each.

```python
from fractions import Fraction as F
from math import comb

def add(a,b):
 r=[F(0)]*max(len(a),len(b))
 for i,v in enumerate(a):r[i]+=v
 for i,v in enumerate(b):r[i]+=v
 return r

def mul(a,b):
 r=[F(0)]*(len(a)+len(b)-1)
 for i,v in enumerate(a):
  for j,w in enumerate(b):r[i+j]+=v*w
 return r

def scale(a,s):return [v*s for v in a]
def power(a,n):
 r=[F(1)]
 for _ in range(n):r=mul(r,a)
 return r

def difference(a,b):return add(a,scale(b,-1))

def bernstein(a,lo,hi):
 n=len(a)-1
 shifted=[F(0)]*(n+1)
 for k,v in enumerate(a):
  for j in range(k+1):shifted[j]+=v*comb(k,j)*lo**(k-j)*(hi-lo)**j
 return [sum((shifted[j]*F(comb(i,j),comb(n,j)) for j in range(i+1)),F(0)) for i in range(n+1)]

def certify(a,lo=F(0),hi=F(1),depth=0):
 coeff=bernstein(a,lo,hi)
 if min(coeff)>=0:return (1,depth)
 if depth>=14:raise ValueError((float(lo),float(hi),float(min(coeff))))
 mid=(lo+hi)/2
 left=certify(a,lo,mid,depth+1);right=certify(a,mid,hi,depth+1)
 return (left[0]+right[0],max(left[1],right[1]))

s=[F(0),F(8),F(-4)];c=[F(4),F(0),F(-4)]
q=[F(5),F(-2),F(1)];p=[F(4),F(0),F(1)]
ds=[F(40),F(-40)];dc=[F(0),F(-40)]
dds=[F(-40),F(-240),F(120)];ddc=[F(-160),F(0),F(120)]
K=F(200000,314159)
# K is a rational upper bound on 2/pi.
norm_num=add(mul(power(s,2),power(p,2)),mul(power(c,2),power(q,2)))
norm_den=mul(power(q,2),power(p,2))
speed_num=add(mul(power(ds,2),power(p,4)),mul(power(dc,2),power(q,4)))
speed_den=mul(power(q,4),power(p,4))
acc_num=add(mul(power(dds,2),power(p,6)),mul(power(ddc,2),power(q,6)))
acc_den=mul(power(q,6),power(p,6))
h_num=add(mul(mul(s,ds),power(p,3)),mul(mul(c,dc),power(q,3)))
h_den=mul(power(q,3),power(p,3))
# g_u'^2 + g*g_u'' common denominator q^4*p^4.
k_num=add(speed_num,add(mul(mul(s,dds),power(p,4)),mul(mul(c,ddc),power(q,4))))
checks={
 'norm2 >= .996':difference(norm_num,scale(norm_den,F(249,250))),
 'speed2 <= 1.02^2':difference(scale(speed_den,F(51,50)**2),scale(speed_num,K*K)),
 'acc2 <= 1.026^2':difference(scale(acc_den,F(513,500)**2),scale(acc_num,K**4)),
 'h^2 <= .005^2':difference(scale(power(h_den,2),F(1,200)**2),scale(power(h_num,2),K*K)),
 'k^2 <= .025^2':difference(scale(power(speed_den,2),F(1,40)**2),scale(power(k_num,2),K**4)),
}
for name,poly in checks.items():
 print(name, 'degree',len(poly)-1,'certificate(leaves,max_depth)',certify(poly),flush=True)

assert F(500, 499) < F(1002, 1000)**2
assert F(250, 499) < F(502, 1000)
assert F(250, 499) < F(708, 1000)**2
A, a, a2 = F(1002, 1000), F(708, 1000), F(502, 1000)
normalized_second = A * (F(1026, 1000) + 2*a2*F(5, 1000)*F(102, 100)
                         + a*F(25, 1000) + 3*a*a2*F(5, 1000)**2)
assert normalized_second < F(106, 100)
assert A * F(102, 100) < F(1025, 1000)
print('normalized second < 1.06; normalized speed < 1.025')
```
