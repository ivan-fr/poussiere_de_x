import Foundation
import Metal

struct ComplexD {
    var re: Double
    var im: Double
}

struct Complex2 {
    var re: Float
    var im: Float
}

func +(lhs: ComplexD, rhs: ComplexD) -> ComplexD {
    ComplexD(re: lhs.re + rhs.re, im: lhs.im + rhs.im)
}

func -(lhs: ComplexD, rhs: ComplexD) -> ComplexD {
    ComplexD(re: lhs.re - rhs.re, im: lhs.im - rhs.im)
}

prefix func -(value: ComplexD) -> ComplexD {
    ComplexD(re: -value.re, im: -value.im)
}

func *(lhs: ComplexD, rhs: ComplexD) -> ComplexD {
    ComplexD(re: lhs.re * rhs.re - lhs.im * rhs.im, im: lhs.re * rhs.im + lhs.im * rhs.re)
}

func *(lhs: Double, rhs: ComplexD) -> ComplexD {
    ComplexD(re: lhs * rhs.re, im: lhs * rhs.im)
}

func *(lhs: ComplexD, rhs: Double) -> ComplexD {
    ComplexD(re: lhs.re * rhs, im: lhs.im * rhs)
}

func /(lhs: ComplexD, rhs: ComplexD) -> ComplexD {
    let den = max(rhs.re * rhs.re + rhs.im * rhs.im, 1.0e-300)
    return ComplexD(re: (lhs.re * rhs.re + lhs.im * rhs.im) / den, im: (lhs.im * rhs.re - lhs.re * rhs.im) / den)
}

func cabs2(_ z: ComplexD) -> Double {
    z.re * z.re + z.im * z.im
}

func isFinite(_ z: ComplexD) -> Bool {
    z.re.isFinite && z.im.isFinite
}

func phase(_ theta: Double) -> ComplexD {
    ComplexD(re: cos(theta), im: sin(theta))
}

func vectorNorm(_ y: [ComplexD]) -> Double {
    sqrt(y.reduce(0.0) { $0 + cabs2($1) })
}

func vectorDistance(_ a: [ComplexD], _ b: [ComplexD]) -> Double {
    var acc = 0.0
    for i in 0..<min(a.count, b.count) {
        acc += cabs2(a[i] - b[i])
    }
    return sqrt(acc)
}

func scaleVector(_ y: [ComplexD], _ scale: Double) -> [ComplexD] {
    y.map { $0 * scale }
}

func addVectors(_ a: [ComplexD], _ b: [ComplexD]) -> [ComplexD] {
    var out: [ComplexD] = []
    out.reserveCapacity(min(a.count, b.count))
    for i in 0..<min(a.count, b.count) {
        out.append(a[i] + b[i])
    }
    return out
}

func splitmix64(_ input: UInt64) -> UInt64 {
    var x = input &+ 0x9E3779B97F4A7C15
    x = ((x ^ (x >> 30)) &* 0xBF58476D1CE4E5B9)
    x = ((x ^ (x >> 27)) &* 0x94D049BB133111EB)
    return x ^ (x >> 31)
}

func u01Hash(_ input: UInt64) -> Double {
    let y = (splitmix64(input) >> 11) & ((UInt64(1) << 53) - 1)
    return Double(y) / Double(UInt64(1) << 53)
}

struct SplitMixRNG {
    var state: UInt64
    var hasSpare = false
    var spare = 0.0

    mutating func nextUInt64() -> UInt64 {
        state = state &+ 0x9E3779B97F4A7C15
        var z = state
        z = ((z ^ (z >> 30)) &* 0xBF58476D1CE4E5B9)
        z = ((z ^ (z >> 27)) &* 0x94D049BB133111EB)
        return z ^ (z >> 31)
    }

    mutating func uniform() -> Double {
        let y = (nextUInt64() >> 11) & ((UInt64(1) << 53) - 1)
        return Double(y) / Double(UInt64(1) << 53)
    }

    mutating func gaussian() -> Double {
        if hasSpare {
            hasSpare = false
            return spare
        }
        let u1 = max(1.0e-300, uniform())
        let u2 = uniform()
        let r = sqrt(-2.0 * log(u1))
        let theta = 2.0 * Double.pi * u2
        spare = r * sin(theta)
        hasSpare = true
        return r * cos(theta)
    }
}

func stableSeed(n: Int, d: Int, seedIndex: Int = 0, salt: Int = 0) -> UInt64 {
    let base = UInt64(0x50414E44524F5349)
        &+ UInt64(1_000_003 * n)
        &+ UInt64(9_176 * d)
        &+ UInt64(97 * seedIndex)
        &+ UInt64(salt)
    return splitmix64(base) & UInt64(0x7fffffff)
}

func compositionsLeq(_ d: Int, _ n: Int) -> [Int32] {
    var out: [Int32] = []
    var cur = Array(repeating: 0, count: n)
    func rec(_ pos: Int, _ remaining: Int) {
        if pos == n - 1 {
            for k in 0...remaining {
                cur[pos] = k
                for v in cur {
                    out.append(Int32(v))
                }
            }
            return
        }
        for k in 0...remaining {
            cur[pos] = k
            rec(pos + 1, remaining - k)
        }
    }
    rec(0, d)
    return out
}

func logFactorials(_ d: Int) -> [Double] {
    var out = Array(repeating: 0.0, count: d + 1)
    if d >= 1 {
        for k in 1...d {
            out[k] = out[k - 1] + log(Double(k))
        }
    }
    return out
}

func kostlanWeights(exps: [Int32], n: Int, d: Int) -> [Double] {
    let terms = exps.count / n
    let logfac = logFactorials(d)
    var weights: [Double] = []
    weights.reserveCapacity(terms)
    for t in 0..<terms {
        var total = 0
        for j in 0..<n {
            total += Int(exps[t * n + j])
        }
        var logs = logfac[d] - logfac[d - total]
        for j in 0..<n {
            logs -= logfac[Int(exps[t * n + j])]
        }
        weights.append(exp(0.5 * logs))
    }
    return weights
}

struct DenseKostlanSystem {
    let n: Int
    let d: Int
    let seed: UInt64
    let exps: [Int32]
    let coeff: [ComplexD]
    let weights: [Double]
    let generationSeconds: Double

    var terms: Int { exps.count / n }
    var totalTerms: Int { n * terms }
    var bezout: Int {
        var out = 1
        for _ in 0..<n {
            out *= d
        }
        return out
    }

    static func make(n: Int, d: Int, seedIndex: Int, equationNormalize: Bool) -> DenseKostlanSystem {
        let t0 = CFAbsoluteTimeGetCurrent()
        let exps = compositionsLeq(d, n)
        let weights = kostlanWeights(exps: exps, n: n, d: d)
        let seed = stableSeed(n: n, d: d, seedIndex: seedIndex)
        var rng = SplitMixRNG(state: seed)
        let terms = weights.count
        var coeff = Array(repeating: ComplexD(re: 0.0, im: 0.0), count: n * terms)
        let invSqrt2 = 1.0 / sqrt(2.0)
        for eq in 0..<n {
            for t in 0..<terms {
                let z = ComplexD(re: rng.gaussian() * invSqrt2 * weights[t], im: rng.gaussian() * invSqrt2 * weights[t])
                coeff[eq * terms + t] = z
            }
        }
        if equationNormalize {
            for eq in 0..<n {
                var norm2 = 0.0
                for t in 0..<terms {
                    norm2 += cabs2(coeff[eq * terms + t])
                }
                let norm = sqrt(norm2)
                if norm > 0.0 {
                    for t in 0..<terms {
                        coeff[eq * terms + t] = coeff[eq * terms + t] * (1.0 / norm)
                    }
                }
            }
        }
        return DenseKostlanSystem(
            n: n,
            d: d,
            seed: seed,
            exps: exps,
            coeff: coeff,
            weights: weights,
            generationSeconds: CFAbsoluteTimeGetCurrent() - t0
        )
    }
}

final class SolveContext {
    let system: DenseKostlanSystem
    var evalCount = 0
    var slopeCount = 0
    var secondsEval = 0.0
    var secondsSlope = 0.0

    init(system: DenseKostlanSystem) {
        self.system = system
    }

    func powers(_ z: ComplexD) -> [ComplexD] {
        var out = Array(repeating: ComplexD(re: 1.0, im: 0.0), count: max(1, system.d + 1))
        if system.d >= 1 {
            out[1] = z
            if system.d >= 2 {
                for k in 2...system.d {
                    out[k] = out[k - 1] * z
                }
            }
        }
        return out
    }

    func eval(_ y: [ComplexD]) -> [ComplexD] {
        let t0 = CFAbsoluteTimeGetCurrent()
        let pows = y.map { powers($0) }
        var out = Array(repeating: ComplexD(re: 0.0, im: 0.0), count: system.n)
        for eq in 0..<system.n {
            var sum = ComplexD(re: 0.0, im: 0.0)
            for t in 0..<system.terms {
                var mon = ComplexD(re: 1.0, im: 0.0)
                for j in 0..<system.n {
                    mon = mon * pows[j][Int(system.exps[t * system.n + j])]
                }
                sum = sum + system.coeff[eq * system.terms + t] * mon
            }
            out[eq] = sum
        }
        evalCount += 1
        secondsEval += CFAbsoluteTimeGetCurrent() - t0
        return out
    }

    func residual(_ y: [ComplexD]) -> Double {
        let f = eval(y)
        let r = sqrt(f.reduce(0.0) { $0 + cabs2($1) })
        return r.isFinite ? r : Double.infinity
    }

    func slopeTable(a: ComplexD, b: ComplexD) -> [ComplexD] {
        var out = Array(repeating: ComplexD(re: 0.0, im: 0.0), count: max(1, system.d + 1))
        if system.d >= 1 {
            var acc = ComplexD(re: 0.0, im: 0.0)
            let powsB = powers(b)
            for m in 1...system.d {
                acc = powsB[m - 1] + a * acc
                out[m] = acc
            }
        }
        return out
    }

    func slopeMatrix(_ a: [ComplexD], _ b: [ComplexD]) -> [ComplexD] {
        let t0 = CFAbsoluteTimeGetCurrent()
        let powsA = a.map { powers($0) }
        let powsB = b.map { powers($0) }
        let slopes = (0..<system.n).map { slopeTable(a: a[$0], b: b[$0]) }
        var q = Array(repeating: ComplexD(re: 0.0, im: 0.0), count: system.n * system.n)
        for eq in 0..<system.n {
            for col in 0..<system.n {
                var sum = ComplexD(re: 0.0, im: 0.0)
                for t in 0..<system.terms {
                    var term = slopes[col][Int(system.exps[t * system.n + col])]
                    if term.re == 0.0 && term.im == 0.0 {
                        continue
                    }
                    if col > 0 {
                        for k in 0..<col {
                            term = term * powsB[k][Int(system.exps[t * system.n + k])]
                        }
                    }
                    if col + 1 < system.n {
                        for k in (col + 1)..<system.n {
                            term = term * powsA[k][Int(system.exps[t * system.n + k])]
                        }
                    }
                    sum = sum + system.coeff[eq * system.terms + t] * term
                }
                q[eq * system.n + col] = sum
            }
        }
        slopeCount += 1
        secondsSlope += CFAbsoluteTimeGetCurrent() - t0
        return q
    }

    func stats() -> [String: Any] {
        [
            "eval_count": evalCount,
            "slope_count": slopeCount,
            "seconds_eval": secondsEval,
            "seconds_slope": secondsSlope,
            "terms_per_poly": system.terms,
            "total_terms": system.totalTerms
        ]
    }
}

func solveLinear(_ matrix: [ComplexD], _ rhs: [ComplexD], n: Int) -> [ComplexD]? {
    var a = matrix
    var b = rhs
    for k in 0..<n {
        var pivot = k
        var pivotNorm = cabs2(a[k * n + k])
        if k + 1 < n {
            for r in (k + 1)..<n {
                let norm = cabs2(a[r * n + k])
                if norm > pivotNorm {
                    pivot = r
                    pivotNorm = norm
                }
            }
        }
        if pivotNorm < 1.0e-300 || !pivotNorm.isFinite {
            return nil
        }
        if pivot != k {
            for c in k..<n {
                a.swapAt(k * n + c, pivot * n + c)
            }
            b.swapAt(k, pivot)
        }
        let diag = a[k * n + k]
        if k + 1 < n {
            for r in (k + 1)..<n {
                let factor = a[r * n + k] / diag
                a[r * n + k] = ComplexD(re: 0.0, im: 0.0)
                for c in (k + 1)..<n {
                    a[r * n + c] = a[r * n + c] - factor * a[k * n + c]
                }
                b[r] = b[r] - factor * b[k]
            }
        }
    }
    var x = Array(repeating: ComplexD(re: 0.0, im: 0.0), count: n)
    for ii in 0..<n {
        let i = n - 1 - ii
        var sum = b[i]
        if i + 1 < n {
            for c in (i + 1)..<n {
                sum = sum - a[i * n + c] * x[c]
            }
        }
        x[i] = sum / a[i * n + i]
        if !isFinite(x[i]) {
            return nil
        }
    }
    return x
}

func finiteResidual(_ ctx: SolveContext, _ y: [ComplexD]) -> Double {
    let r = ctx.residual(y)
    return r.isFinite ? r : Double.infinity
}

func rawDirection(n: Int, trial: UInt64, seed: UInt64, normalize: Bool = true) -> [ComplexD] {
    var vals: [ComplexD] = []
    vals.reserveCapacity(n)
    var norm2 = 0.0
    for j in 0..<n {
        let jj = UInt64(j + 1)
        let h1 = seed &+ 0xD1A5E &+ 0x1000003 &* trial &+ 0x9E37 &* jj
        let h2 = seed &+ 0xBADC0DE &+ 0x1000033 &* trial &+ 0xC2B2 &* jj
        let ang = 2.0 * Double.pi * u01Hash(h1)
        let amp = exp(0.45 * (2.0 * u01Hash(h2) - 1.0))
        let v = amp * phase(ang)
        vals.append(v)
        norm2 += cabs2(v)
    }
    if normalize && norm2 > 0.0 {
        let scale = sqrt(Double(max(1, n))) / sqrt(norm2)
        vals = vals.map { $0 * scale }
    }
    return vals
}

func defaultPowers() -> [Double] {
    var out: [Double] = []
    for k in -20...24 { out.append(pow(2.0, Double(k))) }
    for k in -16...21 { out.append(3.0 * pow(2.0, Double(k))) }
    for k in -14...19 { out.append(5.0 * pow(2.0, Double(k))) }
    for k in -6...6 { out.append(pow(10.0, Double(k))) }
    return out
}

let defaultAnglesDeg: [Double] = [0, 6, 12, 18, 24, 32, 40, 48, 56, 64, 72, 80, 86, 89, 90, 91, 94, 100, 108, 116, 128, 140, 152, 164, 172]
let defaultRadii: [Double] = [0.025, 0.04, 0.06, 0.08, 0.12, 0.18, 0.27, 0.40, 0.60, 0.85, 1.15, 1.55, 2.05, 2.75, 3.60, 4.80, 6.40]
let defaultGains: [Double] = [0.035, 0.055, 0.085, 0.12, 0.18, 0.27, 0.40, 0.58, 0.78, 1.0, 1.28, 1.65, 2.2, 3.0, 4.2, 6.0, 8.5, 12.0]

func parseFloatList(_ raw: String?, default defaults: [Double], positive: Bool = false) -> [Double] {
    guard let raw, !raw.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty else {
        return defaults
    }
    var out: [Double] = []
    for part in raw.replacingOccurrences(of: ";", with: ",").split(separator: ",") {
        if let x = Double(part.trimmingCharacters(in: .whitespacesAndNewlines)), x.isFinite, (!positive || x > 0.0) {
            out.append(x)
        }
    }
    return out.isEmpty ? defaults : out
}

func mobiusHomothetyStart(
    n: Int,
    trial: Int,
    seed: UInt64,
    powers: [Double],
    angles: [Double],
    radii: [Double],
    cap: Double,
    rootsFound: Int = 0,
    duplicates: Int = 0,
    failures: Int = 0,
    targetCount: Int = 1
) -> [ComplexD] {
    let powers2 = Array(Set(powers.filter { $0 > 0.0 }.map { min(max($0, 1.0e-300), cap) })).sorted()
    let pows = powers2.isEmpty ? [1.0] : powers2
    let lp = pows.count
    let la = max(1, angles.count)
    let lr = max(1, radii.count)
    let phi = 0.6180339887498948482
    let q = (Double(trial) * phi + 0.071 * Double(rootsFound) + 0.013 * Double(duplicates)).truncatingRemainder(dividingBy: 1.0)
    let powerIndex = (Int(q * Double(lp)) + 37 * trial + 11 * (trial / 7) + 5 * rootsFound) % lp
    let basePower = pows[powerIndex]
    let dupPressure = (Double(duplicates) + 1.0) / (Double(rootsFound) + 1.0)
    let failPressure = (Double(failures) + 1.0) / (Double(trial) + 1.0)
    let progress = min(1.0, max(0.0, Double(rootsFound) / max(1.0, Double(targetCount))))
    let thrustLadder: [Double] = [1.0, 1.6, 2.5, 4.0, 6.5, 10.0, 16.0, 25.0, 40.0, 64.0, 100.0, 160.0, 256.0]
    let thrust = thrustLadder[(trial * 17 + rootsFound * 3 + duplicates) % thrustLadder.count]
    let amp = pow(thrust, 0.18 + 0.82 * progress) * pow(1.0 + dupPressure, 0.42) / pow(1.0 + 0.25 * failPressure, 0.15)
    let pwr = min(cap, max(1.0e-300, basePower * amp))
    let theta0 = angles[(trial * 19 + rootsFound * 7 + duplicates * 3) % la]
    let thetaJitter = (4.0 * Double.pi / 180.0) * sin(1.324717957244746 * Double(trial + 1) + 0.31 * Double(rootsFound))
    let radius0 = radii[(trial * 13 + failures * 5 + rootsFound * 2) % lr]
    let radius = max(1.0e-300, radius0 * exp(0.22 * sin(0.754877666 * Double(trial + 1) + 0.17 * Double(duplicates))))
    let dir = rawDirection(n: n, trial: UInt64(trial), seed: seed, normalize: true)
    var out: [ComplexD] = []
    out.reserveCapacity(n)
    for j in 0..<n {
        let u = radius * dir[j]
        let hj = seed &+ 0xA11CE &+ UInt64(982451653 * trial) &+ UInt64(1009 * (j + 1))
        let pole = phase(2.0 * Double.pi * u01Hash(hj))
        let tj = theta0 + thetaJitter * cos(0.5 + Double(j)) + (2.0 * Double.pi / 180.0) * sin(Double((j + 1) * (trial + 1)) * 0.38196601125)
        let c = cos(tj)
        let s = sin(tj)
        let w = u / pole
        var denom = ComplexD(re: -s * w.re + c, im: -s * w.im)
        if cabs2(denom) < 1.0e-24 {
            denom = denom + 1.0e-12 * phase(0.37 + Double(j))
        }
        let num = ComplexD(re: c * w.re + s, im: c * w.im)
        out.append(pwr * (pole * (num / denom)))
    }
    return out
}

struct CandidateRow {
    let index: Int
    let trial: Int
    let residual: Double
    let y: [ComplexD]
}

struct MetalParams {
    var n: UInt32
    var equations: UInt32
    var terms: UInt32
    var candidates: UInt32
    var targetCount: UInt32
    var powersCount: UInt32
    var anglesCount: UInt32
    var radiiCount: UInt32
    var seedLo: UInt32
    var seedHi: UInt32
    var powerCap: Float
}

let metalSource = """
#include <metal_stdlib>
using namespace metal;

struct Params {
    uint n;
    uint equations;
    uint terms;
    uint candidates;
    uint targetCount;
    uint powersCount;
    uint anglesCount;
    uint radiiCount;
    uint seedLo;
    uint seedHi;
    float powerCap;
};

static inline float2 cmul(float2 a, float2 b) {
    return float2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}

static inline float2 cdiv(float2 a, float2 b) {
    float den = max(b.x * b.x + b.y * b.y, 1.0e-30f);
    return float2((a.x * b.x + a.y * b.y) / den, (a.y * b.x - a.x * b.y) / den);
}

static inline float2 cphase(float theta) {
    return float2(cos(theta), sin(theta));
}

static inline float2 cpow_int(float2 z, int e) {
    float2 out = float2(1.0f, 0.0f);
    for (int k = 0; k < e; ++k) {
        out = cmul(out, z);
    }
    return out;
}

static inline ulong splitmix64(ulong x) {
    x = (x + 0x9E3779B97F4A7C15UL);
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9UL);
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EBUL);
    return (x ^ (x >> 31));
}

static inline float u01_hash(ulong x) {
    ulong y = splitmix64(x);
    return float((y >> 40) & 0xFFFFFFUL) / 16777216.0f;
}

static inline ulong make_seed(constant Params &p) {
    return (ulong(p.seedHi) << 32) | ulong(p.seedLo);
}

static inline void raw_direction(float2 outv[8], uint n, uint trial, ulong seed, bool normalize) {
    float norm2 = 0.0f;
    for (uint j = 0; j < n; ++j) {
        ulong h1 = splitmix64(seed + 0xD1A5EUL + 0x1000003UL * ulong(trial) + 0x9E37UL * ulong(j + 1));
        ulong h2 = splitmix64(seed + 0xBADC0DEUL + 0x1000033UL * ulong(trial) + 0xC2B2UL * ulong(j + 1));
        float ang = 6.283185307179586f * u01_hash(h1);
        float amp = exp(0.45f * (2.0f * u01_hash(h2) - 1.0f));
        outv[j] = amp * cphase(ang);
        norm2 += dot(outv[j], outv[j]);
    }
    if (normalize && norm2 > 0.0f) {
        float scale = sqrt(float(n)) / sqrt(norm2);
        for (uint j = 0; j < n; ++j) {
            outv[j] *= scale;
        }
    }
}

static inline void mobius_start(float2 y[8], uint trial, constant Params &p, device const float *powers, device const float *angles, device const float *radii) {
    uint pressureSpan = max(1u, p.candidates / max(1u, p.targetCount));
    uint rootsFound = min(max(0u, p.targetCount - 1u), trial / pressureSpan);
    uint duplicates = (trial / 17u + 3u * rootsFound) % max(2u, 2u * p.targetCount + 3u);
    uint failures = (trial / 29u) % max(2u, p.targetCount + 5u);
    uint lp = max(1u, p.powersCount);
    uint la = max(1u, p.anglesCount);
    uint lr = max(1u, p.radiiCount);
    float q = fract(float(trial) * 0.6180339887498948f + 0.071f * float(rootsFound) + 0.013f * float(duplicates));
    uint powerIndex = (uint(q * float(lp)) + 37u * trial + 11u * (trial / 7u) + 5u * rootsFound) % lp;
    float basePower = powers[powerIndex];
    float dupPressure = (float(duplicates) + 1.0f) / (float(rootsFound) + 1.0f);
    float failPressure = (float(failures) + 1.0f) / (float(trial) + 1.0f);
    float progress = clamp(float(rootsFound) / max(1.0f, float(p.targetCount)), 0.0f, 1.0f);
    constexpr float thrustLadder[13] = {1.0f, 1.6f, 2.5f, 4.0f, 6.5f, 10.0f, 16.0f, 25.0f, 40.0f, 64.0f, 100.0f, 160.0f, 256.0f};
    float thrust = thrustLadder[(trial * 17u + rootsFound * 3u + duplicates) % 13u];
    float amp = pow(thrust, 0.18f + 0.82f * progress) * pow(1.0f + dupPressure, 0.42f) / pow(1.0f + 0.25f * failPressure, 0.15f);
    float pwr = min(p.powerCap, max(1.0e-30f, basePower * amp));
    float theta0 = angles[(trial * 19u + rootsFound * 7u + duplicates * 3u) % la];
    float thetaJitter = 0.06981317007977318f * sin(1.324717957244746f * float(trial + 1u) + 0.31f * float(rootsFound));
    float radius0 = radii[(trial * 13u + failures * 5u + rootsFound * 2u) % lr];
    float radius = max(1.0e-30f, radius0 * exp(0.22f * sin(0.754877666f * float(trial + 1u) + 0.17f * float(duplicates))));
    float2 dir[8];
    ulong seed = make_seed(p);
    raw_direction(dir, p.n, trial, seed, true);
    for (uint j = 0; j < p.n; ++j) {
        float2 u = radius * dir[j];
        ulong hj = splitmix64(seed + 0xA11CEUL + 982451653UL * ulong(trial) + 1009UL * ulong(j + 1));
        float2 pole = cphase(6.283185307179586f * u01_hash(hj));
        float tj = theta0 + thetaJitter * cos(0.5f + float(j)) + 0.03490658503988659f * sin(float(j + 1u) * float(trial + 1u) * 0.38196601125f);
        float c = cos(tj);
        float s = sin(tj);
        float2 w = cdiv(u, pole);
        float2 denom = float2(-s * w.x + c, -s * w.y);
        if (dot(denom, denom) < 1.0e-24f) {
            denom += 1.0e-12f * cphase(0.37f + float(j));
        }
        float2 num = float2(c * w.x + s, c * w.y);
        y[j] = pwr * cmul(pole, cdiv(num, denom));
    }
}

kernel void start_eval_dense_ks(
    device const int *exps [[buffer(0)]],
    device const float2 *coeff [[buffer(1)]],
    device const float *powers [[buffer(2)]],
    device const float *angles [[buffer(3)]],
    device const float *radii [[buffer(4)]],
    device float2 *out [[buffer(5)]],
    device float2 *pointsOut [[buffer(6)]],
    constant Params &p [[buffer(7)]],
    uint gid [[thread_position_in_grid]]
) {
    uint total = p.candidates * p.equations;
    if (gid >= total || p.n > 8u) {
        return;
    }
    uint pointIndex = gid / p.equations;
    uint equation = gid - pointIndex * p.equations;
    float2 y[8];
    mobius_start(y, pointIndex, p, powers, angles, radii);
    if (equation == 0u) {
        for (uint j = 0; j < p.n; ++j) {
            pointsOut[pointIndex * p.n + j] = y[j];
        }
    }
    float2 sum = float2(0.0f, 0.0f);
    for (uint t = 0; t < p.terms; ++t) {
        float2 mon = float2(1.0f, 0.0f);
        for (uint j = 0; j < p.n; ++j) {
            int e = exps[t * p.n + j];
            mon = cmul(mon, cpow_int(y[j], e));
        }
        float2 c = coeff[equation * p.terms + t];
        sum += cmul(c, mon);
    }
    out[gid] = sum;
}
"""

final class MetalStartSelector {
    let device: MTLDevice
    let queue: MTLCommandQueue
    let pipeline: MTLComputePipelineState

    init?() {
        guard let device = MTLCreateSystemDefaultDevice(),
              let queue = device.makeCommandQueue() else {
            return nil
        }
        do {
            let library = try device.makeLibrary(source: metalSource, options: nil)
            guard let fn = library.makeFunction(name: "start_eval_dense_ks") else {
                return nil
            }
            self.pipeline = try device.makeComputePipelineState(function: fn)
            self.device = device
            self.queue = queue
        } catch {
            return nil
        }
    }

    func select(system: DenseKostlanSystem, candidates: Int, targetCount: Int, topK: Int, seed: UInt64, powerCap: Double, powers: [Double], angles: [Double], radii: [Double]) throws -> ([CandidateRow], [String: Any]) {
        if system.n > 8 {
            throw NSError(domain: "MetalStartSelector", code: 1, userInfo: [NSLocalizedDescriptionKey: "Metal start selector supports n <= 8"])
        }
        let coeff2 = system.coeff.map { Complex2(re: Float($0.re), im: Float($0.im)) }
        let powersF = powers.map { Float($0) }
        let anglesF = angles.map { Float($0) }
        let radiiF = radii.map { Float($0) }
        let outCount = candidates * system.n
        let pointsCount = candidates * system.n
        guard let expsBuffer = device.makeBuffer(bytes: system.exps, length: system.exps.count * MemoryLayout<Int32>.stride, options: .storageModeShared),
              let coeffBuffer = device.makeBuffer(bytes: coeff2, length: coeff2.count * MemoryLayout<Complex2>.stride, options: .storageModeShared),
              let powersBuffer = device.makeBuffer(bytes: powersF, length: powersF.count * MemoryLayout<Float>.stride, options: .storageModeShared),
              let anglesBuffer = device.makeBuffer(bytes: anglesF, length: anglesF.count * MemoryLayout<Float>.stride, options: .storageModeShared),
              let radiiBuffer = device.makeBuffer(bytes: radiiF, length: radiiF.count * MemoryLayout<Float>.stride, options: .storageModeShared),
              let outBuffer = device.makeBuffer(length: outCount * MemoryLayout<Complex2>.stride, options: .storageModeShared),
              let pointsBuffer = device.makeBuffer(length: pointsCount * MemoryLayout<Complex2>.stride, options: .storageModeShared) else {
            throw NSError(domain: "MetalStartSelector", code: 2, userInfo: [NSLocalizedDescriptionKey: "cannot create Metal buffers"])
        }
        var params = MetalParams(
            n: UInt32(system.n),
            equations: UInt32(system.n),
            terms: UInt32(system.terms),
            candidates: UInt32(candidates),
            targetCount: UInt32(targetCount),
            powersCount: UInt32(powersF.count),
            anglesCount: UInt32(anglesF.count),
            radiiCount: UInt32(radiiF.count),
            seedLo: UInt32(seed & 0xffffffff),
            seedHi: UInt32((seed >> 32) & 0xffffffff),
            powerCap: Float(powerCap)
        )
        guard let paramsBuffer = device.makeBuffer(bytes: &params, length: MemoryLayout<MetalParams>.stride, options: .storageModeShared),
              let commandBuffer = queue.makeCommandBuffer(),
              let encoder = commandBuffer.makeComputeCommandEncoder() else {
            throw NSError(domain: "MetalStartSelector", code: 3, userInfo: [NSLocalizedDescriptionKey: "cannot create command buffer"])
        }
        let t0 = CFAbsoluteTimeGetCurrent()
        encoder.setComputePipelineState(pipeline)
        encoder.setBuffer(expsBuffer, offset: 0, index: 0)
        encoder.setBuffer(coeffBuffer, offset: 0, index: 1)
        encoder.setBuffer(powersBuffer, offset: 0, index: 2)
        encoder.setBuffer(anglesBuffer, offset: 0, index: 3)
        encoder.setBuffer(radiiBuffer, offset: 0, index: 4)
        encoder.setBuffer(outBuffer, offset: 0, index: 5)
        encoder.setBuffer(pointsBuffer, offset: 0, index: 6)
        encoder.setBuffer(paramsBuffer, offset: 0, index: 7)
        let threadsPerGroup = min(pipeline.maxTotalThreadsPerThreadgroup, 256)
        encoder.dispatchThreads(MTLSize(width: outCount, height: 1, depth: 1), threadsPerThreadgroup: MTLSize(width: threadsPerGroup, height: 1, depth: 1))
        encoder.endEncoding()
        commandBuffer.commit()
        commandBuffer.waitUntilCompleted()
        let kernelSeconds = CFAbsoluteTimeGetCurrent() - t0
        let output = outBuffer.contents().bindMemory(to: Complex2.self, capacity: outCount)
        let points = pointsBuffer.contents().bindMemory(to: Complex2.self, capacity: pointsCount)
        var rows: [CandidateRow] = []
        rows.reserveCapacity(candidates)
        for point in 0..<candidates {
            var norm2 = 0.0
            for eq in 0..<system.n {
                let v = output[point * system.n + eq]
                norm2 += Double(v.re) * Double(v.re) + Double(v.im) * Double(v.im)
            }
            let residual = sqrt(norm2)
            if residual.isFinite {
                var y: [ComplexD] = []
                y.reserveCapacity(system.n)
                for j in 0..<system.n {
                    let p = points[point * system.n + j]
                    y.append(ComplexD(re: Double(p.re), im: Double(p.im)))
                }
                rows.append(CandidateRow(index: point, trial: point, residual: residual, y: y))
            }
        }
        rows.sort {
            if $0.residual == $1.residual {
                return $0.index < $1.index
            }
            return $0.residual < $1.residual
        }
        return (
            Array(rows.prefix(min(max(1, topK), rows.count))),
            [
                "device": device.name,
                "kernel_seconds": kernelSeconds,
                "candidates": candidates,
                "top_k": topK,
                "mode": "metal-start-selector"
            ]
        )
    }
}

func cpuSelectStarts(ctx: SolveContext, args: Args, powers: [Double], angles: [Double], radii: [Double]) -> ([CandidateRow], [String: Any]) {
    let t0 = CFAbsoluteTimeGetCurrent()
    var rows: [CandidateRow] = []
    rows.reserveCapacity(args.metalCandidates)
    let pressureSpan = max(1, args.metalCandidates / max(1, args.count))
    for trial in 0..<args.metalCandidates {
        let rootsFound = min(max(0, args.count - 1), trial / pressureSpan)
        let duplicates = (trial / 17 + 3 * rootsFound) % max(2, 2 * args.count + 3)
        let failures = (trial / 29) % max(2, args.count + 5)
        let y = mobiusHomothetyStart(n: ctx.system.n, trial: trial, seed: ctx.system.seed + 0x113000, powers: powers, angles: angles, radii: radii, cap: args.powerCap, rootsFound: rootsFound, duplicates: duplicates, failures: failures, targetCount: args.count)
        let r = finiteResidual(ctx, y)
        if r.isFinite {
            rows.append(CandidateRow(index: trial, trial: trial, residual: r, y: y))
        }
    }
    rows.sort {
        if $0.residual == $1.residual {
            return $0.index < $1.index
        }
        return $0.residual < $1.residual
    }
    let topK = args.selectorTop > 0 ? args.selectorTop : args.refineTop * 8
    return (
        Array(rows.prefix(min(topK, rows.count))),
        [
            "kernel_seconds": 0.0,
            "total_seconds": CFAbsoluteTimeGetCurrent() - t0,
            "candidates": args.metalCandidates,
            "top_k": topK,
            "mode": "cpu-start-selector"
        ]
    )
}

func diversify(_ selected: [CandidateRow], limit: Int, sep: Double) -> [CandidateRow] {
    let wanted = max(1, limit)
    if sep <= 0.0 {
        return Array(selected.prefix(wanted))
    }
    var chosen: [CandidateRow] = []
    var chosenPoints: [[ComplexD]] = []
    var deferred: [CandidateRow] = []
    for cand in selected {
        let zn = max(1.0, vectorNorm(cand.y))
        var near = false
        for prev in chosenPoints {
            let dist = vectorDistance(cand.y, prev) / max(zn, vectorNorm(prev), 1.0)
            if dist < sep {
                near = true
                break
            }
        }
        if near {
            deferred.append(cand)
        } else {
            chosen.append(cand)
            chosenPoints.append(cand.y)
            if chosen.count >= wanted {
                return chosen
            }
        }
    }
    var seen = Set(chosen.map { $0.index })
    for cand in deferred {
        if seen.contains(cand.index) {
            continue
        }
        chosen.append(cand)
        seen.insert(cand.index)
        if chosen.count >= wanted {
            break
        }
    }
    return chosen
}

func probeEndpoint(ctx: SolveContext, y: [ComplexD], residual: Double, prevDelta: [ComplexD]?, epoch: Int, directionSeed: UInt64, probeScale: Double, probeCandidates: Int, probeRadii: [Double], includeSelf: Bool) throws -> ([ComplexD], [String: Any]) {
    let n = y.count
    let ynorm = max(1.0, vectorNorm(y))
    let usableRadii = probeRadii.filter { $0 >= 0.0 }.isEmpty ? [1.0] : probeRadii.filter { $0 >= 0.0 }
    var candidates: [(String, [ComplexD])] = []
    if includeSelf {
        candidates.append(("self", y))
    }
    if let prevDelta {
        let pdn = max(1.0e-300, vectorNorm(prevDelta))
        let base = scaleVector(prevDelta, min(max(pdn, probeScale * ynorm), 2.5 * ynorm) / pdn)
        candidates.append(("inertial", addVectors(y, base)))
    }
    let budget = max(1, probeCandidates)
    var k = 0
    while candidates.count < budget {
        let rad = probeScale * ynorm * usableRadii[k % usableRadii.count]
        var qdir = rawDirection(n: n, trial: UInt64(epoch + 1) &* 104729 &+ UInt64(k + 1) &* 7919 &+ directionSeed, seed: directionSeed ^ (0x116116 &+ UInt64(17 * k)), normalize: true)
        let qnorm = max(1.0e-300, vectorNorm(qdir))
        qdir = scaleVector(qdir, sqrt(Double(max(1, n))) / qnorm)
        let ph = phase(0.6180339887498948 * Double(epoch + 1) + 2.399963229728653 * Double(k + 1))
        var step: [ComplexD] = []
        step.reserveCapacity(n)
        for j in 0..<n {
            var s = rad * (ph * qdir[j])
            if ynorm > 0.0 {
                s = s + (0.12 * rad / ynorm) * (y[j] * phase(0.38196601125 * Double(k + 1)))
            }
            let tiny = 1.0e-12 * ynorm
            if sqrt(cabs2(s)) < tiny {
                s = s + tiny * phase(0.17 + Double(j) + Double(epoch) + Double(k))
            }
            step.append(s)
        }
        candidates.append(("geom-\(k)", addVectors(y, step)))
        k += 1
    }
    var bestName = ""
    var bestB: [ComplexD]? = nil
    var bestRes = Double.infinity
    var bestDistance = 0.0
    var evals = 0
    for (name, b) in candidates.prefix(budget) {
        let rb = finiteResidual(ctx, b)
        evals += 1
        let dist = vectorDistance(b, y)
        let score = rb.isFinite ? log1p(max(0.0, rb)) + 1.0e-14 * log1p(dist) : Double.infinity
        let old = bestRes.isFinite ? log1p(max(0.0, bestRes)) + 1.0e-14 * log1p(bestDistance) : Double.infinity
        if score.isFinite && score < old {
            bestName = name
            bestB = b
            bestRes = rb
            bestDistance = dist
        }
    }
    guard let selected = bestB else {
        throw NSError(domain: "Pandrosion", code: 1, userInfo: [NSLocalizedDescriptionKey: "no finite probe"])
    }
    return (
        selected,
        [
            "probe_name": bestName,
            "probe_candidates": min(budget, candidates.count),
            "probe_evals": evals,
            "probe_residual": bestRes,
            "probe_distance": bestDistance
        ]
    )
}

func pandrosionCorrector(ctx: SolveContext, y0: [ComplexD], args: Args, directionSeed: UInt64, maxEpochs: Int? = nil, acceptOverride: Double? = nil, lineSearchOverride: Int? = nil) -> [String: Any] {
    let t0 = CFAbsoluteTimeGetCurrent()
    let epochsLimit = max(1, maxEpochs ?? args.epochs)
    let accept = acceptOverride ?? args.accept
    let lineSearch = lineSearchOverride ?? args.lineSearch
    var y = y0
    var bestY = y
    var bestR = finiteResidual(ctx, y)
    var ok = false
    var status = "started"
    var epochs = 0
    var prevDelta: [ComplexD]? = nil
    var totalProbeEvals = 0
    for ep in 0..<epochsLimit {
        let f = ctx.eval(y)
        let r = sqrt(f.reduce(0.0) { $0 + cabs2($1) })
        if r.isFinite && r < bestR {
            bestR = r
            bestY = y
        }
        if r <= max(args.tol, accept) && (accept <= 0.0 || r < accept) {
            ok = true
            status = "converged"
            break
        }
        let b: [ComplexD]
        do {
            let probe = try probeEndpoint(ctx: ctx, y: y, residual: r, prevDelta: prevDelta, epoch: ep, directionSeed: directionSeed, probeScale: args.probeScale, probeCandidates: args.probeCandidates, probeRadii: args.probeRadiiList, includeSelf: args.probeSelf)
            b = probe.0
            if let pe = probe.1["probe_evals"] as? Int {
                totalProbeEvals += pe
            }
        } catch {
            status = "probe-error"
            break
        }
        let q = ctx.slopeMatrix(y, b)
        guard var delta = solveLinear(q, f.map { -$0 }, n: y.count) else {
            status = "slope-solve-error"
            break
        }
        let ynorm = max(1.0, vectorNorm(y))
        let dnorm = vectorNorm(delta)
        if dnorm > 18.0 * ynorm {
            delta = scaleVector(delta, (18.0 * ynorm) / max(dnorm, 1.0e-300))
        }
        var accepted = false
        let baseR = r
        for k in 0..<max(1, lineSearch) {
            let lambda = 1.0 / pow(2.0, Double(k))
            let yy = addVectors(y, scaleVector(delta, lambda))
            let rr = finiteResidual(ctx, yy)
            if rr.isFinite && (rr < baseR || rr < bestR) {
                prevDelta = scaleVector(delta, lambda)
                y = yy
                if rr < bestR {
                    bestY = yy
                    bestR = rr
                }
                accepted = true
                break
            }
        }
        epochs = ep + 1
        if !accepted {
            status = "no-decrease"
            break
        }
    }
    let finalR = finiteResidual(ctx, bestY)
    if finalR <= max(args.tol, accept) && (accept <= 0.0 || finalR < accept) {
        ok = true
        status = "converged"
    }
    return [
        "accepted": accept <= 0.0 ? ok : (finalR.isFinite && finalR < accept),
        "ok": ok,
        "status": status,
        "epochs": epochs,
        "residual": finalR,
        "y": bestY,
        "seconds": CFAbsoluteTimeGetCurrent() - t0,
        "probe_total_evals": totalProbeEvals
    ]
}

func startopt(ctx: SolveContext, y0: [ComplexD], trial: Int, seed: UInt64, args: Args) -> ([ComplexD], [String: Any]) {
    var best = y0
    var bestR = finiteResidual(ctx, best)
    let initial = bestR
    var evals = 1
    var microTotal = 0
    var chosenGain = 1.0
    let gains = args.startoptGainsList.isEmpty ? [1.0] : args.startoptGainsList
    if args.startoptSteps > 0 {
        for step in 0..<args.startoptSteps {
            let base = best
            var pool: [(Double, [ComplexD])] = [(1.0, base)]
            if args.startoptCandidates > 1 {
                for c in 0..<(args.startoptCandidates - 1) {
                    let gain = gains[(trial + 3 * step + c) % gains.count]
                    var cand = scaleVector(base, gain)
                    if c % 3 == 1 {
                        var pert: [ComplexD] = []
                        for j in 0..<cand.count {
                            let h1 = seed &+ 0x5157A47 &+ UInt64(65537 * trial) &+ UInt64(4099 * c) &+ UInt64(193 * (j + 1))
                            let h2 = seed &+ 0x7150F7 &+ UInt64(104729 * trial) &+ UInt64(8191 * c) &+ UInt64(389 * (j + 1))
                            let ph = 0.23 * (2.0 * u01Hash(h1) - 1.0)
                            let amp = exp(0.28 * (2.0 * u01Hash(h2) - 1.0))
                            pert.append(cand[j] * amp * phase(ph))
                        }
                        cand = pert
                    } else if c % 3 == 2 {
                        let fresh = rawDirection(n: base.count, trial: UInt64(trial + 31 * (step + 1) + c), seed: seed, normalize: true)
                        let nm = max(1.0e-300, vectorNorm(base))
                        let freshNorm = max(1.0e-300, vectorNorm(fresh))
                        cand = addVectors(scaleVector(cand, 0.70), scaleVector(fresh, 0.30 * gain * nm / freshNorm))
                    }
                    pool.append((gain, cand))
                }
            }
            for (gain, cand0) in pool {
                evals += 1
                var yscore = cand0
                var r = finiteResidual(ctx, yscore)
                if args.startoptMicroEpochs > 0 {
                    let loc = pandrosionCorrector(ctx: ctx, y0: yscore, args: args, directionSeed: seed, maxEpochs: args.startoptMicroEpochs, acceptOverride: 0.0, lineSearchOverride: 6)
                    microTotal += loc["epochs"] as? Int ?? 0
                    if let yy = loc["y"] as? [ComplexD], let rr = loc["residual"] as? Double, rr < r {
                        yscore = yy
                        r = rr
                    }
                }
                let score = r.isFinite ? log1p(max(0.0, r)) + 1.0e-15 * log1p(vectorNorm(yscore)) : Double.infinity
                let old = bestR.isFinite ? log1p(max(0.0, bestR)) + 1.0e-15 * log1p(vectorNorm(best)) : Double.infinity
                if score < old {
                    best = yscore
                    bestR = r
                    chosenGain = gain
                }
            }
        }
    }
    return (
        best,
        [
            "startopt_r0": initial,
            "startopt_r1": bestR,
            "startopt_ratio": (initial.isFinite && bestR.isFinite && initial > 0.0) ? bestR / initial : NSNull(),
            "startopt_evals": evals,
            "startopt_micro_epochs": microTotal,
            "startopt_gain": chosenGain
        ]
    )
}

func clusterIndex(_ roots: [[String: Any]], _ z: [ComplexD], sep: Double) -> Int? {
    let zn = max(1.0, vectorNorm(z))
    for (i, root) in roots.enumerated() {
        guard let rzJSON = root["z"] as? [[Double]] else {
            continue
        }
        let rz = rzJSON.map { ComplexD(re: $0[0], im: $0[1]) }
        let dist = vectorDistance(z, rz) / max(zn, vectorNorm(rz), 1.0)
        if dist <= sep {
            return i
        }
    }
    return nil
}

func rootJSON(_ y: [ComplexD]) -> [[Double]] {
    y.map { [$0.re, $0.im] }
}

func parseCase(_ raw: String) -> (Int, Int) {
    let clean = raw.replacingOccurrences(of: "x", with: ",").replacingOccurrences(of: ":", with: ",")
    let parts = clean.split(separator: ",").map { String($0).trimmingCharacters(in: .whitespacesAndNewlines) }
    if parts.count != 2 {
        fatalError("case must be n,d, got \(raw)")
    }
    guard let n = Int(parts[0]), let d = Int(parts[1]), n > 0, d > 0 else {
        fatalError("n,d must be positive")
    }
    return (n, d)
}

struct Args {
    var cases = "2,34"
    var seedIndex = 0
    var equationNormalize = false
    var count = 8
    var pool = 512
    var epochs = 24
    var tol = 1.0e-12
    var accept = 1.0e-8
    var clusterSep = 1.0e-8
    var lineSearch = 12
    var probeScale = 0.035
    var probeCandidates = 8
    var probeRadii = "0,0.35,0.7,1,1.6,2.6,4.2"
    var probeSelf = true
    var powers: String? = nil
    var powerCap = 1_048_576.0
    var angles: String? = nil
    var rays: String? = nil
    var startoptSteps = 1
    var startoptCandidates = 12
    var startoptGains: String? = nil
    var startoptMicroEpochs = 0
    var metalCandidates = 2048
    var selectorTop = 0
    var refineTop = 96
    var diversitySep = 0.45
    var out: String? = nil
    var outdir = "verification/129_swift_metal_full"
    var verboseTrials = false
    var keepTrials = 160

    var powersList: [Double] { parseFloatList(powers, default: defaultPowers(), positive: true).map { min(max($0, 1.0e-300), powerCap) }.sorted() }
    var anglesList: [Double] { parseFloatList(angles, default: defaultAnglesDeg).map { $0 * Double.pi / 180.0 } }
    var radiiList: [Double] { parseFloatList(rays, default: defaultRadii, positive: true) }
    var probeRadiiList: [Double] { parseFloatList(probeRadii, default: [0.0, 0.35, 0.7, 1.0, 1.6, 2.6, 4.2]).filter { $0 >= 0.0 } }
    var startoptGainsList: [Double] { parseFloatList(startoptGains, default: defaultGains, positive: true) }
}

func parseArgs() -> Args {
    var args = Args()
    let argv = CommandLine.arguments
    var i = 1
    func value(_ name: String) -> String {
        if i + 1 >= argv.count {
            fatalError("missing value for \(name)")
        }
        i += 1
        return argv[i]
    }
    while i < argv.count {
        let key = argv[i]
        switch key {
        case "--cases": args.cases = value(key)
        case "--seed-index": args.seedIndex = Int(value(key)) ?? args.seedIndex
        case "--equation-normalize": args.equationNormalize = true
        case "--no-equation-normalize": args.equationNormalize = false
        case "--count", "--thales-count", "--useful-count": args.count = Int(value(key)) ?? args.count
        case "--pool", "--thales-pool", "--useful-pool": args.pool = Int(value(key)) ?? args.pool; args.metalCandidates = args.pool
        case "--epochs", "--thales-epochs", "--useful-epochs": args.epochs = Int(value(key)) ?? args.epochs
        case "--tol": args.tol = Double(value(key)) ?? args.tol
        case "--accept", "--residual-accept": args.accept = Double(value(key)) ?? args.accept
        case "--cluster-sep": args.clusterSep = Double(value(key)) ?? args.clusterSep
        case "--line-search": args.lineSearch = Int(value(key)) ?? args.lineSearch
        case "--probe-scale": args.probeScale = Double(value(key)) ?? args.probeScale
        case "--probe-candidates": args.probeCandidates = Int(value(key)) ?? args.probeCandidates
        case "--probe-radii": args.probeRadii = value(key)
        case "--probe-self": args.probeSelf = true
        case "--no-probe-self": args.probeSelf = false
        case "--powers", "--thales-powers": args.powers = value(key)
        case "--power-cap", "--thales-power-cap": args.powerCap = Double(value(key)) ?? args.powerCap
        case "--angles", "--thales-angles": args.angles = value(key)
        case "--rays", "--thales-rays": args.rays = value(key)
        case "--startopt-steps": args.startoptSteps = Int(value(key)) ?? args.startoptSteps
        case "--startopt-candidates": args.startoptCandidates = Int(value(key)) ?? args.startoptCandidates
        case "--startopt-gains": args.startoptGains = value(key)
        case "--startopt-micro-epochs": args.startoptMicroEpochs = Int(value(key)) ?? args.startoptMicroEpochs
        case "--metal-candidates": args.metalCandidates = Int(value(key)) ?? args.metalCandidates
        case "--selector-top": args.selectorTop = Int(value(key)) ?? args.selectorTop
        case "--refine-top": args.refineTop = Int(value(key)) ?? args.refineTop
        case "--diversity-sep": args.diversitySep = Double(value(key)) ?? args.diversitySep
        case "--out", "--thales-out", "--useful-out": args.out = value(key)
        case "--outdir": args.outdir = value(key)
        case "--verbose-trials": args.verboseTrials = true
        case "--keep-trials": args.keepTrials = Int(value(key)) ?? args.keepTrials
        default:
            fatalError("unknown arg \(key)")
        }
        i += 1
    }
    return args
}

func runCase(_ args: Args, _ caseRaw: String, selector: MetalStartSelector?) -> [String: Any] {
    let tCase = CFAbsoluteTimeGetCurrent()
    let (n, d) = parseCase(caseRaw)
    let system = DenseKostlanSystem.make(n: n, d: d, seedIndex: args.seedIndex, equationNormalize: args.equationNormalize)
    let ctx = SolveContext(system: system)
    let powers = args.powersList
    let angles = args.anglesList
    let radii = args.radiiList
    let selectorTop = args.selectorTop > 0 ? args.selectorTop : args.refineTop * 8
    let selectedAll: [CandidateRow]
    let selectorStats: [String: Any]
    if let selector, n <= 8 {
        do {
            let result = try selector.select(system: system, candidates: args.metalCandidates, targetCount: args.count, topK: selectorTop, seed: system.seed + 0x113000, powerCap: args.powerCap, powers: powers, angles: angles, radii: radii)
            selectedAll = result.0
            selectorStats = result.1
        } catch {
            let result = cpuSelectStarts(ctx: ctx, args: args, powers: powers, angles: angles, radii: radii)
            selectedAll = result.0
            selectorStats = result.1
        }
    } else {
        let result = cpuSelectStarts(ctx: ctx, args: args, powers: powers, angles: angles, radii: radii)
        selectedAll = result.0
        selectorStats = result.1
    }
    let selected = diversify(selectedAll, limit: args.refineTop, sep: args.diversitySep)
    var roots: [[String: Any]] = []
    var trials: [[String: Any]] = []
    var duplicates = 0
    var failures = 0
    let tExtract = CFAbsoluteTimeGetCurrent()
    for cand in selected {
        if roots.count >= args.count {
            break
        }
        let start = startopt(ctx: ctx, y0: cand.y, trial: cand.trial, seed: system.seed + 0x112555, args: args)
        let loc = pandrosionCorrector(ctx: ctx, y0: start.0, args: args, directionSeed: system.seed + UInt64(7919 * cand.trial))
        let y = loc["y"] as? [ComplexD] ?? start.0
        let r = finiteResidual(ctx, y)
        let accepted = r.isFinite && r < args.accept
        var rec: [String: Any] = [
            "trial": cand.trial,
            "selector_rank": cand.index,
            "selector_residual": cand.residual,
            "accepted": accepted,
            "status": loc["status"] as? String ?? "",
            "r1": r,
            "epochs": loc["epochs"] as? Int ?? 0,
            "seconds": loc["seconds"] as? Double ?? 0.0
        ]
        for (k, v) in start.1 {
            rec[k] = v
        }
        if !accepted {
            failures += 1
            trials.append(rec)
            continue
        }
        if let dup = clusterIndex(roots, y, sep: args.clusterSep) {
            duplicates += 1
            rec["status"] = "duplicate"
            rec["cluster"] = dup
            trials.append(rec)
            continue
        }
        let rid = roots.count
        let root: [String: Any] = [
            "id": rid,
            "source": "129-swift-metal-full",
            "trial": cand.trial,
            "selector_rank": cand.index,
            "selector_residual": cand.residual,
            "residual": r,
            "epochs": loc["epochs"] as? Int ?? 0,
            "seconds": loc["seconds"] as? Double ?? 0.0,
            "z": rootJSON(y),
            "y": rootJSON(y)
        ]
        roots.append(root)
        rec["status"] = "new-root"
        rec["root_id"] = rid
        trials.append(rec)
    }
    let extractSeconds = CFAbsoluteTimeGetCurrent() - tExtract
    return [
        "script": "129_pandrosion_swift_metal_full.swift",
        "autonomous": true,
        "mode": "swift-metal-full-general-kostlan/metal-start-selector/swift-double-pandrosion",
        "case": "\(n),\(d)",
        "family": "ks",
        "seed_index": args.seedIndex,
        "seed": Int(system.seed),
        "n": n,
        "degree": d,
        "terms_per_poly": system.terms,
        "terms": system.totalTerms,
        "bezout": system.bezout,
        "roots": roots,
        "trials": args.verboseTrials ? trials : Array(trials.prefix(min(trials.count, args.keepTrials))),
        "selector": selectorStats,
        "summary": [
            "requested_roots": args.count,
            "unique_roots": roots.count,
            "success": roots.count >= args.count,
            "trials_used": trials.count,
            "duplicates": duplicates,
            "failures": failures,
            "generation_seconds": system.generationSeconds,
            "extract_seconds": extractSeconds,
            "total_seconds": CFAbsoluteTimeGetCurrent() - tCase,
            "eval_stats": ctx.stats()
        ]
    ]
}

let args = parseArgs()
let selector = MetalStartSelector()
let cases = args.cases.replacingOccurrences(of: "|", with: ";").split(separator: ";").map { String($0).trimmingCharacters(in: .whitespacesAndNewlines) }.filter { !$0.isEmpty }
let outputs = cases.map { runCase(args, $0, selector: selector) }
let final: [String: Any] = outputs.count == 1 ? outputs[0] : ["script": "129_pandrosion_swift_metal_full.swift", "autonomous": true, "cases": outputs]
let outPath: String
if let out = args.out {
    outPath = out
} else {
    let first = (cases.first ?? "case").replacingOccurrences(of: ",", with: "x")
    outPath = "\(args.outdir)/129_swift_metal_full_\(first).json"
}
let outURL = URL(fileURLWithPath: outPath)
try FileManager.default.createDirectory(at: outURL.deletingLastPathComponent(), withIntermediateDirectories: true)
let jsonData = try JSONSerialization.data(withJSONObject: final, options: [.prettyPrinted, .sortedKeys])
try jsonData.write(to: outURL)

print(String(repeating: "=", count: 120))
print("129 autonomous Swift + Metal full general Pandrosion")
print("Kostlan generation in Swift; Metal batch start selector; Swift Double generic finite-slope corrector")
print(String(repeating: "=", count: 120))
for result in outputs {
    let summary = result["summary"] as? [String: Any] ?? [:]
    let n = result["n"] as? Int ?? 0
    let d = result["degree"] as? Int ?? 0
    let seed = result["seed"] as? Int ?? 0
    let terms = result["terms"] as? Int ?? 0
    print("case=ks(\(n),\(d)), seed=\(seed), terms=\(terms), roots=\(summary["unique_roots"] ?? 0)/\(summary["requested_roots"] ?? 0), success=\(summary["success"] ?? false)")
    let gen = summary["generation_seconds"] as? Double ?? 0.0
    let ext = summary["extract_seconds"] as? Double ?? 0.0
    let total = summary["total_seconds"] as? Double ?? 0.0
    print(String(format: "seconds: generation=%.4f, extract=%.4f, total=%.4f", gen, ext, total))
}
print("out=\(outPath)")
