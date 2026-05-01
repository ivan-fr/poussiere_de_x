import Foundation
import Metal

struct StartSelectJob: Decodable {
    let n: Int
    let equations: Int
    let terms: Int
    let candidates: Int
    let targetCount: Int
    let topK: Int
    let seed: UInt64
    let powerCap: Float
    let exps: [Int32]
    let coeffRe: [Float]
    let coeffIm: [Float]
    let powers: [Float]
    let angles: [Float]
    let radii: [Float]
}

struct LoadSystemJob: Decodable {
    let op: String
    let n: Int
    let equations: Int
    let terms: Int
    let exps: [Int32]
    let coeffRe: [Float]
    let coeffIm: [Float]
    let coeffRe64: [Double]?
    let coeffIm64: [Double]?
}

struct StoredStartSelectJob: Decodable {
    let op: String
    let candidates: Int
    let targetCount: Int
    let topK: Int
    let seed: UInt64
    let powerCap: Float
    let powers: [Float]
    let angles: [Float]
    let radii: [Float]
}

struct EvalPointsJob: Decodable {
    let op: String
    let points: Int
    let pointsRe: [Float]
    let pointsIm: [Float]
}

struct Polish2Job: Decodable {
    let op: String
    let points: Int
    let epochs: Int
    let probeCandidates: Int
    let lineSearch: Int
    let seed: UInt64
    let probeScale: Float
    let pointsRe: [Float]
    let pointsIm: [Float]
}

struct Solve2Job: Decodable {
    let op: String
    let candidates: Int
    let targetCount: Int
    let topK: Int
    let refineTop: Int
    let count: Int
    let seed: UInt64
    let powerCap: Float
    let powers: [Float]
    let angles: [Float]
    let radii: [Float]
    let diversitySep: Double
    let epochs: Int
    let tol: Double
    let accept: Double
    let clusterSep: Double
    let lineSearch: Int
    let probeScale: Double
    let probeCandidates: Int
    let probeRadii: [Double]
    let probeSelf: Bool
    let startoptSteps: Int
    let startoptCandidates: Int
    let startoptGains: [Double]
    let startoptMicroEpochs: Int
    let metalPolish2: Bool
    let metalPolishEpochs: Int
    let metalPolishProbes: Int
    let metalPolishLineSearch: Int
    let metalPolishProbeScale: Float
}

struct Params {
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

struct PolishParams {
    var n: UInt32
    var equations: UInt32
    var terms: UInt32
    var points: UInt32
    var epochs: UInt32
    var probeCandidates: UInt32
    var lineSearch: UInt32
    var seedLo: UInt32
    var seedHi: UInt32
    var probeScale: Float
}

struct Complex2 {
    var re: Float
    var im: Float
}

struct ComplexD {
    var re: Double
    var im: Double
}

struct Candidate {
    let index: Int
    let trial: Int
    let residual: Double
}

final class LoadedSystem {
    let n: Int
    let equations: Int
    let terms: Int
    let degree: Int
    let exps: [Int32]
    let coeffD: [ComplexD]
    let expsBuffer: MTLBuffer
    let coeffBuffer: MTLBuffer

    init(n: Int, equations: Int, terms: Int, degree: Int, exps: [Int32], coeffD: [ComplexD], expsBuffer: MTLBuffer, coeffBuffer: MTLBuffer) {
        self.n = n
        self.equations = equations
        self.terms = terms
        self.degree = degree
        self.exps = exps
        self.coeffD = coeffD
        self.expsBuffer = expsBuffer
        self.coeffBuffer = coeffBuffer
    }
}

func fail(_ message: String) -> Never {
    FileHandle.standardError.write((message + "\n").data(using: .utf8)!)
    exit(1)
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

struct PolishParams {
    uint n;
    uint equations;
    uint terms;
    uint points;
    uint epochs;
    uint probeCandidates;
    uint lineSearch;
    uint seedLo;
    uint seedHi;
    float probeScale;
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

static inline ulong make_seed(constant Params &p) {
    return (ulong(p.seedHi) << 32) | ulong(p.seedLo);
}

static inline ulong splitmix64(ulong x) {
    x = (x + 0x9E3779B97F4A7C15UL);
    x = ((x ^ (x >> 30)) * 0xBF58476D1CE4E5B9UL);
    x = ((x ^ (x >> 27)) * 0x94D049BB133111EBUL);
    return (x ^ (x >> 31));
}

static inline ulong make_seed_polish(constant PolishParams &p) {
    return (ulong(p.seedHi) << 32) | ulong(p.seedLo);
}

static inline float u01_hash(ulong x) {
    ulong y = splitmix64(x);
    return float((y >> 40) & 0xFFFFFFUL) / 16777216.0f;
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

static inline void mobius_start(
    float2 y[8],
    uint trial,
    constant Params &p,
    device const float *powers,
    device const float *angles,
    device const float *radii
) {
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

kernel void eval_points_residuals(
    device const int *exps [[buffer(0)]],
    device const float2 *coeff [[buffer(1)]],
    device const float2 *points [[buffer(2)]],
    device float *residuals [[buffer(3)]],
    constant Params &p [[buffer(4)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid >= p.candidates || p.n > 8u) {
        return;
    }

    float norm2 = 0.0f;
    for (uint equation = 0; equation < p.equations; ++equation) {
        float2 sum = float2(0.0f, 0.0f);
        for (uint t = 0; t < p.terms; ++t) {
            float2 mon = float2(1.0f, 0.0f);
            for (uint j = 0; j < p.n; ++j) {
                int e = exps[t * p.n + j];
                float2 z = points[gid * p.n + j];
                mon = cmul(mon, cpow_int(z, e));
            }
            float2 c = coeff[equation * p.terms + t];
            sum += cmul(c, mon);
        }
        norm2 += dot(sum, sum);
    }
    residuals[gid] = sqrt(norm2);
}

static inline float eval_system_residual2(
    float2 y0,
    float2 y1,
    device const int *exps,
    device const float2 *coeff,
    constant PolishParams &p,
    thread float2 f[2]
) {
    f[0] = float2(0.0f, 0.0f);
    f[1] = float2(0.0f, 0.0f);
    for (uint equation = 0; equation < 2u; ++equation) {
        float2 sum = float2(0.0f, 0.0f);
        for (uint t = 0; t < p.terms; ++t) {
            int e0 = exps[t * 2u];
            int e1 = exps[t * 2u + 1u];
            float2 mon = cmul(cpow_int(y0, e0), cpow_int(y1, e1));
            float2 c = coeff[equation * p.terms + t];
            sum += cmul(c, mon);
        }
        f[equation] = sum;
    }
    return sqrt(dot(f[0], f[0]) + dot(f[1], f[1]));
}

static inline void slope_sum_pow(float2 a, float2 b, int e, thread float2 &out) {
    float2 acc = float2(0.0f, 0.0f);
    for (int r = 0; r < e; ++r) {
        float2 left = cpow_int(b, e - 1 - r);
        float2 right = cpow_int(a, r);
        acc += cmul(left, right);
    }
    out = acc;
}

static inline void slope_matrix2(
    float2 a0,
    float2 a1,
    float2 b0,
    float2 b1,
    device const int *exps,
    device const float2 *coeff,
    constant PolishParams &p,
    thread float2 q[4]
) {
    q[0] = float2(0.0f, 0.0f);
    q[1] = float2(0.0f, 0.0f);
    q[2] = float2(0.0f, 0.0f);
    q[3] = float2(0.0f, 0.0f);
    for (uint equation = 0; equation < 2u; ++equation) {
        float2 col0 = float2(0.0f, 0.0f);
        float2 col1 = float2(0.0f, 0.0f);
        for (uint t = 0; t < p.terms; ++t) {
            int e0 = exps[t * 2u];
            int e1 = exps[t * 2u + 1u];
            float2 s0;
            float2 s1;
            slope_sum_pow(a0, b0, e0, s0);
            slope_sum_pow(a1, b1, e1, s1);
            float2 term0 = cmul(s0, cpow_int(a1, e1));
            float2 term1 = cmul(cpow_int(b0, e0), s1);
            float2 c = coeff[equation * p.terms + t];
            col0 += cmul(c, term0);
            col1 += cmul(c, term1);
        }
        q[equation * 2u] = col0;
        q[equation * 2u + 1u] = col1;
    }
}

static inline bool solve2(thread float2 q[4], float2 f0, float2 f1, thread float2 &d0, thread float2 &d1) {
    float2 a = q[0];
    float2 b = q[1];
    float2 c = q[2];
    float2 d = q[3];
    float2 det = cmul(a, d) - cmul(b, c);
    float den = dot(det, det);
    if (!isfinite(den) || den < 1.0e-30f) {
        return false;
    }
    float2 nf0 = -f0;
    float2 nf1 = -f1;
    d0 = cdiv(cmul(d, nf0) - cmul(b, nf1), det);
    d1 = cdiv(cmul(a, nf1) - cmul(c, nf0), det);
    return isfinite(d0.x) && isfinite(d0.y) && isfinite(d1.x) && isfinite(d1.y);
}

static inline float probe_radius(uint k) {
    switch (k % 7u) {
        case 0u: return 0.0f;
        case 1u: return 0.35f;
        case 2u: return 0.7f;
        case 3u: return 1.0f;
        case 4u: return 1.6f;
        case 5u: return 2.6f;
        default: return 4.2f;
    }
}

kernel void polish2_points(
    device const int *exps [[buffer(0)]],
    device const float2 *coeff [[buffer(1)]],
    device const float2 *pointsIn [[buffer(2)]],
    device float2 *pointsOut [[buffer(3)]],
    device float *residualsOut [[buffer(4)]],
    constant PolishParams &p [[buffer(5)]],
    uint gid [[thread_position_in_grid]]
) {
    if (gid >= p.points || p.n != 2u || p.equations != 2u) {
        return;
    }
    float2 y0 = pointsIn[gid * 2u];
    float2 y1 = pointsIn[gid * 2u + 1u];
    float2 best0 = y0;
    float2 best1 = y1;
    float2 f[2];
    float bestR = eval_system_residual2(y0, y1, exps, coeff, p, f);
    ulong seed = make_seed_polish(p) + 7919UL * ulong(gid);

    for (uint ep = 0; ep < p.epochs; ++ep) {
        float r = eval_system_residual2(y0, y1, exps, coeff, p, f);
        if (isfinite(r) && r < bestR) {
            bestR = r;
            best0 = y0;
            best1 = y1;
        }

        float ynorm = max(1.0f, sqrt(dot(y0, y0) + dot(y1, y1)));
        float2 bestB0 = y0;
        float2 bestB1 = y1;
        float bestProbeR = r;
        uint budget = max(1u, p.probeCandidates);
        for (uint k = 1; k < budget; ++k) {
            float rad = p.probeScale * ynorm * probe_radius(k - 1u);
            float2 dir[8];
            raw_direction(dir, 2u, uint(seed + 104729UL * ulong(ep + 1u) + 7919UL * ulong(k + 1u)), seed ^ (0x116116UL + 17UL * ulong(k)), true);
            float ph = 0.6180339887498948f * float(ep + 1u) + 2.399963229728653f * float(k + 1u);
            float2 phasev = cphase(ph);
            float2 step0 = rad * cmul(phasev, dir[0]);
            float2 step1 = rad * cmul(phasev, dir[1]);
            if (ynorm > 0.0f) {
                float2 radialPhase = cphase(0.38196601125f * float(k + 1u));
                step0 += (0.12f * rad / ynorm) * cmul(y0, radialPhase);
                step1 += (0.12f * rad / ynorm) * cmul(y1, radialPhase);
            }
            float2 b0 = y0 + step0;
            float2 b1 = y1 + step1;
            float2 pf[2];
            float pr = eval_system_residual2(b0, b1, exps, coeff, p, pf);
            if (isfinite(pr) && pr < bestProbeR) {
                bestProbeR = pr;
                bestB0 = b0;
                bestB1 = b1;
            }
        }

        float2 q[4];
        slope_matrix2(y0, y1, bestB0, bestB1, exps, coeff, p, q);
        float2 d0;
        float2 d1;
        if (!solve2(q, f[0], f[1], d0, d1)) {
            break;
        }
        float dnorm = sqrt(dot(d0, d0) + dot(d1, d1));
        if (dnorm > 18.0f * ynorm) {
            float scale = (18.0f * ynorm) / max(dnorm, 1.0e-30f);
            d0 *= scale;
            d1 *= scale;
        }

        bool accepted = false;
        uint lineSteps = max(1u, p.lineSearch);
        for (uint ls = 0; ls < lineSteps; ++ls) {
            float lam = exp2(-float(ls));
            float2 yy0 = y0 + lam * d0;
            float2 yy1 = y1 + lam * d1;
            float2 lf[2];
            float rr = eval_system_residual2(yy0, yy1, exps, coeff, p, lf);
            if (isfinite(rr) && (rr < r || rr < bestR)) {
                y0 = yy0;
                y1 = yy1;
                if (rr < bestR) {
                    bestR = rr;
                    best0 = yy0;
                    best1 = yy1;
                }
                accepted = true;
                break;
            }
        }
        if (!accepted) {
            break;
        }
    }
    pointsOut[gid * 2u] = best0;
    pointsOut[gid * 2u + 1u] = best1;
    residualsOut[gid] = bestR;
}
"""

guard let device = MTLCreateSystemDefaultDevice() else {
    fail("no Metal device")
}
guard let queue = device.makeCommandQueue() else {
    fail("cannot create Metal command queue")
}
let library = try device.makeLibrary(source: metalSource, options: nil)
guard let function = library.makeFunction(name: "start_eval_dense_ks") else {
    fail("missing start_eval_dense_ks function")
}
let pipeline = try device.makeComputePipelineState(function: function)
guard let evalFunction = library.makeFunction(name: "eval_points_residuals") else {
    fail("missing eval_points_residuals function")
}
let evalPipeline = try device.makeComputePipelineState(function: evalFunction)
guard let polishFunction = library.makeFunction(name: "polish2_points") else {
    fail("missing polish2_points function")
}
let polishPipeline = try device.makeComputePipelineState(function: polishFunction)

func runJob(_ job: StartSelectJob) throws -> [String: Any] {
    if job.n > 8 {
        throw NSError(domain: "StartSelect", code: 2, userInfo: [NSLocalizedDescriptionKey: "n > 8 is not supported"])
    }
    if job.exps.count != job.terms * job.n {
        throw NSError(domain: "StartSelect", code: 3, userInfo: [NSLocalizedDescriptionKey: "bad exps length"])
    }
    if job.coeffRe.count != job.equations * job.terms || job.coeffIm.count != job.equations * job.terms {
        throw NSError(domain: "StartSelect", code: 4, userInfo: [NSLocalizedDescriptionKey: "bad coeff length"])
    }
    if job.powers.isEmpty || job.angles.isEmpty || job.radii.isEmpty {
        throw NSError(domain: "StartSelect", code: 5, userInfo: [NSLocalizedDescriptionKey: "empty powers/angles/radii"])
    }

    let totalStart = CFAbsoluteTimeGetCurrent()
    let coeff = zip(job.coeffRe, job.coeffIm).map { Complex2(re: $0.0, im: $0.1) }
    let outCount = job.candidates * job.equations
    let pointsCount = job.candidates * job.n

    guard let expsBuffer = device.makeBuffer(bytes: job.exps, length: job.exps.count * MemoryLayout<Int32>.stride, options: .storageModeShared),
          let coeffBuffer = device.makeBuffer(bytes: coeff, length: coeff.count * MemoryLayout<Complex2>.stride, options: .storageModeShared),
          let powersBuffer = device.makeBuffer(bytes: job.powers, length: job.powers.count * MemoryLayout<Float>.stride, options: .storageModeShared),
          let anglesBuffer = device.makeBuffer(bytes: job.angles, length: job.angles.count * MemoryLayout<Float>.stride, options: .storageModeShared),
          let radiiBuffer = device.makeBuffer(bytes: job.radii, length: job.radii.count * MemoryLayout<Float>.stride, options: .storageModeShared),
          let outBuffer = device.makeBuffer(length: outCount * MemoryLayout<Complex2>.stride, options: .storageModeShared),
          let pointsBuffer = device.makeBuffer(length: pointsCount * MemoryLayout<Complex2>.stride, options: .storageModeShared) else {
        throw NSError(domain: "StartSelect", code: 6, userInfo: [NSLocalizedDescriptionKey: "cannot create Metal buffers"])
    }

    var params = Params(
        n: UInt32(job.n),
        equations: UInt32(job.equations),
        terms: UInt32(job.terms),
        candidates: UInt32(job.candidates),
        targetCount: UInt32(job.targetCount),
        powersCount: UInt32(job.powers.count),
        anglesCount: UInt32(job.angles.count),
        radiiCount: UInt32(job.radii.count),
        seedLo: UInt32(job.seed & 0xFFFFFFFF),
        seedHi: UInt32((job.seed >> 32) & 0xFFFFFFFF),
        powerCap: job.powerCap
    )
    guard let paramsBuffer = device.makeBuffer(bytes: &params, length: MemoryLayout<Params>.stride, options: .storageModeShared) else {
        throw NSError(domain: "StartSelect", code: 7, userInfo: [NSLocalizedDescriptionKey: "cannot create params buffer"])
    }

    let kernelStart = CFAbsoluteTimeGetCurrent()
    guard let commandBuffer = queue.makeCommandBuffer(),
          let encoder = commandBuffer.makeComputeCommandEncoder() else {
        throw NSError(domain: "StartSelect", code: 8, userInfo: [NSLocalizedDescriptionKey: "cannot create command buffer"])
    }
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
    let grid = MTLSize(width: outCount, height: 1, depth: 1)
    let group = MTLSize(width: threadsPerGroup, height: 1, depth: 1)
    encoder.dispatchThreads(grid, threadsPerThreadgroup: group)
    encoder.endEncoding()
    commandBuffer.commit()
    commandBuffer.waitUntilCompleted()
    let kernelSeconds = CFAbsoluteTimeGetCurrent() - kernelStart

    let selectStart = CFAbsoluteTimeGetCurrent()
    let output = outBuffer.contents().bindMemory(to: Complex2.self, capacity: outCount)
    let points = pointsBuffer.contents().bindMemory(to: Complex2.self, capacity: pointsCount)
    var candidates: [Candidate] = []
    candidates.reserveCapacity(job.candidates)
    for point in 0..<job.candidates {
        var norm = 0.0
        for eq in 0..<job.equations {
            let v = output[point * job.equations + eq]
            let re = Double(v.re)
            let im = Double(v.im)
            norm += re * re + im * im
        }
        let residual = sqrt(norm)
        if residual.isFinite {
            candidates.append(Candidate(index: point, trial: point, residual: residual))
        }
    }
    candidates.sort {
        if $0.residual == $1.residual {
            return $0.index < $1.index
        }
        return $0.residual < $1.residual
    }
    let selected = Array(candidates.prefix(min(max(1, job.topK), candidates.count)))
    let selectSeconds = CFAbsoluteTimeGetCurrent() - selectStart

    var rows: [[String: Any]] = []
    rows.reserveCapacity(selected.count)
    for cand in selected {
        var zr: [Double] = []
        var zi: [Double] = []
        zr.reserveCapacity(job.n)
        zi.reserveCapacity(job.n)
        for j in 0..<job.n {
            let p = points[cand.index * job.n + j]
            zr.append(Double(p.re))
            zi.append(Double(p.im))
        }
        rows.append([
            "rank": rows.count,
            "index": cand.index,
            "trial": cand.trial,
            "residual": cand.residual,
            "re": zr,
            "im": zi
        ])
    }

    return [
        "device": device.name,
        "candidates": job.candidates,
        "terms": job.terms,
        "equations": job.equations,
        "top_k": job.topK,
        "selected": rows,
        "kernel_seconds": kernelSeconds,
        "select_seconds": selectSeconds,
        "total_seconds": CFAbsoluteTimeGetCurrent() - totalStart
    ]
}

func loadSystem(_ job: LoadSystemJob) throws -> LoadedSystem {
    if job.n > 8 {
        throw NSError(domain: "StartSelect", code: 20, userInfo: [NSLocalizedDescriptionKey: "n > 8 is not supported"])
    }
    if job.exps.count != job.terms * job.n {
        throw NSError(domain: "StartSelect", code: 21, userInfo: [NSLocalizedDescriptionKey: "bad exps length"])
    }
    if job.coeffRe.count != job.equations * job.terms || job.coeffIm.count != job.equations * job.terms {
        throw NSError(domain: "StartSelect", code: 22, userInfo: [NSLocalizedDescriptionKey: "bad coeff length"])
    }
    if let coeffRe64 = job.coeffRe64, coeffRe64.count != job.equations * job.terms {
        throw NSError(domain: "StartSelect", code: 24, userInfo: [NSLocalizedDescriptionKey: "bad coeffRe64 length"])
    }
    if let coeffIm64 = job.coeffIm64, coeffIm64.count != job.equations * job.terms {
        throw NSError(domain: "StartSelect", code: 25, userInfo: [NSLocalizedDescriptionKey: "bad coeffIm64 length"])
    }
    let coeff = zip(job.coeffRe, job.coeffIm).map { Complex2(re: $0.0, im: $0.1) }
    var coeffD: [ComplexD] = []
    coeffD.reserveCapacity(job.equations * job.terms)
    if let coeffRe64 = job.coeffRe64, let coeffIm64 = job.coeffIm64 {
        for i in 0..<(job.equations * job.terms) {
            coeffD.append(ComplexD(re: coeffRe64[i], im: coeffIm64[i]))
        }
    } else {
        for i in 0..<(job.equations * job.terms) {
            coeffD.append(ComplexD(re: Double(job.coeffRe[i]), im: Double(job.coeffIm[i])))
        }
    }
    let degree = Int(job.exps.max() ?? 0)
    guard let expsBuffer = device.makeBuffer(bytes: job.exps, length: job.exps.count * MemoryLayout<Int32>.stride, options: .storageModeShared),
          let coeffBuffer = device.makeBuffer(bytes: coeff, length: coeff.count * MemoryLayout<Complex2>.stride, options: .storageModeShared) else {
        throw NSError(domain: "StartSelect", code: 23, userInfo: [NSLocalizedDescriptionKey: "cannot create loaded buffers"])
    }
    return LoadedSystem(n: job.n, equations: job.equations, terms: job.terms, degree: degree, exps: job.exps, coeffD: coeffD, expsBuffer: expsBuffer, coeffBuffer: coeffBuffer)
}

func runStoredJob(_ job: StoredStartSelectJob, system: LoadedSystem) throws -> [String: Any] {
    if job.powers.isEmpty || job.angles.isEmpty || job.radii.isEmpty {
        throw NSError(domain: "StartSelect", code: 30, userInfo: [NSLocalizedDescriptionKey: "empty powers/angles/radii"])
    }

    let totalStart = CFAbsoluteTimeGetCurrent()
    let outCount = job.candidates * system.equations
    let pointsCount = job.candidates * system.n

    guard let powersBuffer = device.makeBuffer(bytes: job.powers, length: job.powers.count * MemoryLayout<Float>.stride, options: .storageModeShared),
          let anglesBuffer = device.makeBuffer(bytes: job.angles, length: job.angles.count * MemoryLayout<Float>.stride, options: .storageModeShared),
          let radiiBuffer = device.makeBuffer(bytes: job.radii, length: job.radii.count * MemoryLayout<Float>.stride, options: .storageModeShared),
          let outBuffer = device.makeBuffer(length: outCount * MemoryLayout<Complex2>.stride, options: .storageModeShared),
          let pointsBuffer = device.makeBuffer(length: pointsCount * MemoryLayout<Complex2>.stride, options: .storageModeShared) else {
        throw NSError(domain: "StartSelect", code: 31, userInfo: [NSLocalizedDescriptionKey: "cannot create Metal buffers"])
    }

    var params = Params(
        n: UInt32(system.n),
        equations: UInt32(system.equations),
        terms: UInt32(system.terms),
        candidates: UInt32(job.candidates),
        targetCount: UInt32(job.targetCount),
        powersCount: UInt32(job.powers.count),
        anglesCount: UInt32(job.angles.count),
        radiiCount: UInt32(job.radii.count),
        seedLo: UInt32(job.seed & 0xFFFFFFFF),
        seedHi: UInt32((job.seed >> 32) & 0xFFFFFFFF),
        powerCap: job.powerCap
    )
    guard let paramsBuffer = device.makeBuffer(bytes: &params, length: MemoryLayout<Params>.stride, options: .storageModeShared) else {
        throw NSError(domain: "StartSelect", code: 32, userInfo: [NSLocalizedDescriptionKey: "cannot create params buffer"])
    }

    let kernelStart = CFAbsoluteTimeGetCurrent()
    guard let commandBuffer = queue.makeCommandBuffer(),
          let encoder = commandBuffer.makeComputeCommandEncoder() else {
        throw NSError(domain: "StartSelect", code: 33, userInfo: [NSLocalizedDescriptionKey: "cannot create command buffer"])
    }
    encoder.setComputePipelineState(pipeline)
    encoder.setBuffer(system.expsBuffer, offset: 0, index: 0)
    encoder.setBuffer(system.coeffBuffer, offset: 0, index: 1)
    encoder.setBuffer(powersBuffer, offset: 0, index: 2)
    encoder.setBuffer(anglesBuffer, offset: 0, index: 3)
    encoder.setBuffer(radiiBuffer, offset: 0, index: 4)
    encoder.setBuffer(outBuffer, offset: 0, index: 5)
    encoder.setBuffer(pointsBuffer, offset: 0, index: 6)
    encoder.setBuffer(paramsBuffer, offset: 0, index: 7)

    let threadsPerGroup = min(pipeline.maxTotalThreadsPerThreadgroup, 256)
    let grid = MTLSize(width: outCount, height: 1, depth: 1)
    let group = MTLSize(width: threadsPerGroup, height: 1, depth: 1)
    encoder.dispatchThreads(grid, threadsPerThreadgroup: group)
    encoder.endEncoding()
    commandBuffer.commit()
    commandBuffer.waitUntilCompleted()
    let kernelSeconds = CFAbsoluteTimeGetCurrent() - kernelStart

    let selectStart = CFAbsoluteTimeGetCurrent()
    let output = outBuffer.contents().bindMemory(to: Complex2.self, capacity: outCount)
    let points = pointsBuffer.contents().bindMemory(to: Complex2.self, capacity: pointsCount)
    var candidates: [Candidate] = []
    candidates.reserveCapacity(job.candidates)
    for point in 0..<job.candidates {
        var norm = 0.0
        for eq in 0..<system.equations {
            let v = output[point * system.equations + eq]
            let re = Double(v.re)
            let im = Double(v.im)
            norm += re * re + im * im
        }
        let residual = sqrt(norm)
        if residual.isFinite {
            candidates.append(Candidate(index: point, trial: point, residual: residual))
        }
    }
    candidates.sort {
        if $0.residual == $1.residual {
            return $0.index < $1.index
        }
        return $0.residual < $1.residual
    }
    let selected = Array(candidates.prefix(min(max(1, job.topK), candidates.count)))
    let selectSeconds = CFAbsoluteTimeGetCurrent() - selectStart

    var rows: [[String: Any]] = []
    rows.reserveCapacity(selected.count)
    for cand in selected {
        var zr: [Double] = []
        var zi: [Double] = []
        zr.reserveCapacity(system.n)
        zi.reserveCapacity(system.n)
        for j in 0..<system.n {
            let p = points[cand.index * system.n + j]
            zr.append(Double(p.re))
            zi.append(Double(p.im))
        }
        rows.append([
            "rank": rows.count,
            "index": cand.index,
            "trial": cand.trial,
            "residual": cand.residual,
            "re": zr,
            "im": zi
        ])
    }

    return [
        "device": device.name,
        "candidates": job.candidates,
        "terms": system.terms,
        "equations": system.equations,
        "top_k": job.topK,
        "selected": rows,
        "kernel_seconds": kernelSeconds,
        "select_seconds": selectSeconds,
        "total_seconds": CFAbsoluteTimeGetCurrent() - totalStart,
        "stateful": true
    ]
}

func runEvalPoints(_ job: EvalPointsJob, system: LoadedSystem) throws -> [String: Any] {
    if job.pointsRe.count != job.points * system.n || job.pointsIm.count != job.points * system.n {
        throw NSError(domain: "StartSelect", code: 50, userInfo: [NSLocalizedDescriptionKey: "bad points length"])
    }
    let totalStart = CFAbsoluteTimeGetCurrent()
    let points = zip(job.pointsRe, job.pointsIm).map { Complex2(re: $0.0, im: $0.1) }
    guard let pointsBuffer = device.makeBuffer(bytes: points, length: points.count * MemoryLayout<Complex2>.stride, options: .storageModeShared),
          let residualsBuffer = device.makeBuffer(length: job.points * MemoryLayout<Float>.stride, options: .storageModeShared) else {
        throw NSError(domain: "StartSelect", code: 51, userInfo: [NSLocalizedDescriptionKey: "cannot create eval buffers"])
    }
    var params = Params(
        n: UInt32(system.n),
        equations: UInt32(system.equations),
        terms: UInt32(system.terms),
        candidates: UInt32(job.points),
        targetCount: 1,
        powersCount: 0,
        anglesCount: 0,
        radiiCount: 0,
        seedLo: 0,
        seedHi: 0,
        powerCap: 0.0
    )
    guard let paramsBuffer = device.makeBuffer(bytes: &params, length: MemoryLayout<Params>.stride, options: .storageModeShared) else {
        throw NSError(domain: "StartSelect", code: 52, userInfo: [NSLocalizedDescriptionKey: "cannot create eval params"])
    }
    let kernelStart = CFAbsoluteTimeGetCurrent()
    guard let commandBuffer = queue.makeCommandBuffer(),
          let encoder = commandBuffer.makeComputeCommandEncoder() else {
        throw NSError(domain: "StartSelect", code: 53, userInfo: [NSLocalizedDescriptionKey: "cannot create eval command buffer"])
    }
    encoder.setComputePipelineState(evalPipeline)
    encoder.setBuffer(system.expsBuffer, offset: 0, index: 0)
    encoder.setBuffer(system.coeffBuffer, offset: 0, index: 1)
    encoder.setBuffer(pointsBuffer, offset: 0, index: 2)
    encoder.setBuffer(residualsBuffer, offset: 0, index: 3)
    encoder.setBuffer(paramsBuffer, offset: 0, index: 4)
    let threadsPerGroup = min(evalPipeline.maxTotalThreadsPerThreadgroup, 256)
    let grid = MTLSize(width: job.points, height: 1, depth: 1)
    let group = MTLSize(width: threadsPerGroup, height: 1, depth: 1)
    encoder.dispatchThreads(grid, threadsPerThreadgroup: group)
    encoder.endEncoding()
    commandBuffer.commit()
    commandBuffer.waitUntilCompleted()
    let kernelSeconds = CFAbsoluteTimeGetCurrent() - kernelStart

    let residualPtr = residualsBuffer.contents().bindMemory(to: Float.self, capacity: job.points)
    var residuals: [Double] = []
    residuals.reserveCapacity(job.points)
    for i in 0..<job.points {
        residuals.append(Double(residualPtr[i]))
    }
    return [
        "op": "eval_points",
        "points": job.points,
        "residuals": residuals,
        "kernel_seconds": kernelSeconds,
        "total_seconds": CFAbsoluteTimeGetCurrent() - totalStart,
        "stateful": true
    ]
}

func runPolish2(_ job: Polish2Job, system: LoadedSystem) throws -> [String: Any] {
    if system.n != 2 || system.equations != 2 {
        throw NSError(domain: "StartSelect", code: 60, userInfo: [NSLocalizedDescriptionKey: "polish2 requires n=2 and equations=2"])
    }
    if job.pointsRe.count != job.points * system.n || job.pointsIm.count != job.points * system.n {
        throw NSError(domain: "StartSelect", code: 61, userInfo: [NSLocalizedDescriptionKey: "bad polish points length"])
    }

    let totalStart = CFAbsoluteTimeGetCurrent()
    let points = zip(job.pointsRe, job.pointsIm).map { Complex2(re: $0.0, im: $0.1) }
    guard let pointsInBuffer = device.makeBuffer(bytes: points, length: points.count * MemoryLayout<Complex2>.stride, options: .storageModeShared),
          let pointsOutBuffer = device.makeBuffer(length: points.count * MemoryLayout<Complex2>.stride, options: .storageModeShared),
          let residualsBuffer = device.makeBuffer(length: job.points * MemoryLayout<Float>.stride, options: .storageModeShared) else {
        throw NSError(domain: "StartSelect", code: 62, userInfo: [NSLocalizedDescriptionKey: "cannot create polish buffers"])
    }
    var params = PolishParams(
        n: UInt32(system.n),
        equations: UInt32(system.equations),
        terms: UInt32(system.terms),
        points: UInt32(job.points),
        epochs: UInt32(max(0, job.epochs)),
        probeCandidates: UInt32(max(1, job.probeCandidates)),
        lineSearch: UInt32(max(1, job.lineSearch)),
        seedLo: UInt32(job.seed & 0xFFFFFFFF),
        seedHi: UInt32((job.seed >> 32) & 0xFFFFFFFF),
        probeScale: job.probeScale
    )
    guard let paramsBuffer = device.makeBuffer(bytes: &params, length: MemoryLayout<PolishParams>.stride, options: .storageModeShared) else {
        throw NSError(domain: "StartSelect", code: 63, userInfo: [NSLocalizedDescriptionKey: "cannot create polish params"])
    }

    let kernelStart = CFAbsoluteTimeGetCurrent()
    guard let commandBuffer = queue.makeCommandBuffer(),
          let encoder = commandBuffer.makeComputeCommandEncoder() else {
        throw NSError(domain: "StartSelect", code: 64, userInfo: [NSLocalizedDescriptionKey: "cannot create polish command buffer"])
    }
    encoder.setComputePipelineState(polishPipeline)
    encoder.setBuffer(system.expsBuffer, offset: 0, index: 0)
    encoder.setBuffer(system.coeffBuffer, offset: 0, index: 1)
    encoder.setBuffer(pointsInBuffer, offset: 0, index: 2)
    encoder.setBuffer(pointsOutBuffer, offset: 0, index: 3)
    encoder.setBuffer(residualsBuffer, offset: 0, index: 4)
    encoder.setBuffer(paramsBuffer, offset: 0, index: 5)
    let threadsPerGroup = min(polishPipeline.maxTotalThreadsPerThreadgroup, 256)
    let grid = MTLSize(width: job.points, height: 1, depth: 1)
    let group = MTLSize(width: threadsPerGroup, height: 1, depth: 1)
    encoder.dispatchThreads(grid, threadsPerThreadgroup: group)
    encoder.endEncoding()
    commandBuffer.commit()
    commandBuffer.waitUntilCompleted()
    let kernelSeconds = CFAbsoluteTimeGetCurrent() - kernelStart

    let pointsOut = pointsOutBuffer.contents().bindMemory(to: Complex2.self, capacity: points.count)
    let residualsOut = residualsBuffer.contents().bindMemory(to: Float.self, capacity: job.points)
    var candidates: [Candidate] = []
    candidates.reserveCapacity(job.points)
    for point in 0..<job.points {
        let residual = Double(residualsOut[point])
        if residual.isFinite {
            candidates.append(Candidate(index: point, trial: point, residual: residual))
        }
    }
    candidates.sort {
        if $0.residual == $1.residual {
            return $0.index < $1.index
        }
        return $0.residual < $1.residual
    }

    var rows: [[String: Any]] = []
    rows.reserveCapacity(candidates.count)
    for cand in candidates {
        var zr: [Double] = []
        var zi: [Double] = []
        zr.reserveCapacity(system.n)
        zi.reserveCapacity(system.n)
        for j in 0..<system.n {
            let p = pointsOut[cand.index * system.n + j]
            zr.append(Double(p.re))
            zi.append(Double(p.im))
        }
        rows.append([
            "rank": rows.count,
            "index": cand.index,
            "trial": cand.trial,
            "residual": cand.residual,
            "re": zr,
            "im": zi
        ])
    }

    return [
        "op": "polish2",
        "device": device.name,
        "points": job.points,
        "epochs": job.epochs,
        "probe_candidates": job.probeCandidates,
        "line_search": job.lineSearch,
        "selected": rows,
        "kernel_seconds": kernelSeconds,
        "total_seconds": CFAbsoluteTimeGetCurrent() - totalStart,
        "stateful": true
    ]
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

func phaseD(_ theta: Double) -> ComplexD {
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

func subVectors(_ a: [ComplexD], _ b: [ComplexD]) -> [ComplexD] {
    var out: [ComplexD] = []
    out.reserveCapacity(min(a.count, b.count))
    for i in 0..<min(a.count, b.count) {
        out.append(a[i] - b[i])
    }
    return out
}

func splitmix64D(_ input: UInt64) -> UInt64 {
    var x = input &+ 0x9E3779B97F4A7C15
    x = ((x ^ (x >> 30)) &* 0xBF58476D1CE4E5B9)
    x = ((x ^ (x >> 27)) &* 0x94D049BB133111EB)
    return x ^ (x >> 31)
}

func u01D(_ input: UInt64) -> Double {
    let y = (splitmix64D(input) >> 11) & ((UInt64(1) << 53) - 1)
    return Double(y) / Double(UInt64(1) << 53)
}

func rawDirectionD(n: Int, trial: UInt64, seed: UInt64, normalize: Bool = true) -> [ComplexD] {
    var vals: [ComplexD] = []
    vals.reserveCapacity(n)
    var norm2 = 0.0
    for j in 0..<n {
        let jj = UInt64(j + 1)
        let h1 = seed &+ 0xD1A5E &+ 0x1000003 &* trial &+ 0x9E37 &* jj
        let h2 = seed &+ 0xBADC0DE &+ 0x1000033 &* trial &+ 0xC2B2 &* jj
        let angle = 2.0 * Double.pi * u01D(h1)
        let amp = exp(0.45 * (2.0 * u01D(h2) - 1.0))
        let v = amp * phaseD(angle)
        vals.append(v)
        norm2 += cabs2(v)
    }
    if normalize && norm2 > 0.0 {
        let scale = sqrt(Double(max(1, n))) / sqrt(norm2)
        vals = vals.map { $0 * scale }
    }
    return vals
}

final class SolveContext {
    let system: LoadedSystem
    var evalCount = 0
    var slopeCount = 0
    var secondsEval = 0.0
    var secondsSlope = 0.0

    init(system: LoadedSystem) {
        self.system = system
    }

    func powers(_ z: ComplexD) -> [ComplexD] {
        var out = Array(repeating: ComplexD(re: 1.0, im: 0.0), count: max(1, system.degree + 1))
        if system.degree >= 1 {
            out[1] = z
            if system.degree >= 2 {
                for k in 2...system.degree {
                    out[k] = out[k - 1] * z
                }
            }
        }
        return out
    }

    func eval(_ y: [ComplexD]) -> [ComplexD] {
        let t0 = CFAbsoluteTimeGetCurrent()
        let p0 = powers(y[0])
        let p1 = powers(y[1])
        var out = Array(repeating: ComplexD(re: 0.0, im: 0.0), count: system.equations)
        for eq in 0..<system.equations {
            var sum = ComplexD(re: 0.0, im: 0.0)
            for t in 0..<system.terms {
                let e0 = Int(system.exps[t * system.n])
                let e1 = Int(system.exps[t * system.n + 1])
                let mon = p0[e0] * p1[e1]
                sum = sum + system.coeffD[eq * system.terms + t] * mon
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
        var out = Array(repeating: ComplexD(re: 0.0, im: 0.0), count: max(1, system.degree + 1))
        if system.degree >= 1 {
            var acc = ComplexD(re: 0.0, im: 0.0)
            let powsB = powers(b)
            for m in 1...system.degree {
                acc = powsB[m - 1] + a * acc
                out[m] = acc
            }
        }
        return out
    }

    func slopeMatrix(_ a: [ComplexD], _ b: [ComplexD]) -> [ComplexD] {
        let t0 = CFAbsoluteTimeGetCurrent()
        let powsA1 = powers(a[1])
        let powsB0 = powers(b[0])
        let s0 = slopeTable(a: a[0], b: b[0])
        let s1 = slopeTable(a: a[1], b: b[1])
        var q = Array(repeating: ComplexD(re: 0.0, im: 0.0), count: 4)
        for eq in 0..<system.equations {
            var col0 = ComplexD(re: 0.0, im: 0.0)
            var col1 = ComplexD(re: 0.0, im: 0.0)
            for t in 0..<system.terms {
                let e0 = Int(system.exps[t * system.n])
                let e1 = Int(system.exps[t * system.n + 1])
                let c = system.coeffD[eq * system.terms + t]
                col0 = col0 + c * (s0[e0] * powsA1[e1])
                col1 = col1 + c * (powsB0[e0] * s1[e1])
            }
            q[eq * 2] = col0
            q[eq * 2 + 1] = col1
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
            "total_terms": system.terms * system.equations
        ]
    }
}

func solve2x2(_ q: [ComplexD], _ f: [ComplexD]) -> [ComplexD]? {
    let det = q[0] * q[3] - q[1] * q[2]
    if !det.re.isFinite || !det.im.isFinite || cabs2(det) < 1.0e-300 {
        return nil
    }
    let nf0 = -f[0]
    let nf1 = -f[1]
    let d0 = (q[3] * nf0 - q[1] * nf1) / det
    let d1 = (q[0] * nf1 - q[2] * nf0) / det
    if !isFinite(d0) || !isFinite(d1) {
        return nil
    }
    return [d0, d1]
}

func rowPoint(_ row: [String: Any]) -> [ComplexD] {
    let re = row["re"] as? [Double] ?? []
    let im = row["im"] as? [Double] ?? []
    var out: [ComplexD] = []
    out.reserveCapacity(min(re.count, im.count))
    for i in 0..<min(re.count, im.count) {
        out.append(ComplexD(re: re[i], im: im[i]))
    }
    return out
}

func rootJSON(_ y: [ComplexD]) -> [[Double]] {
    y.map { [$0.re, $0.im] }
}

func doubleValue(_ row: [String: Any], _ key: String, _ fallback: Double = 0.0) -> Double {
    if let x = row[key] as? Double {
        return x
    }
    if let x = row[key] as? Int {
        return Double(x)
    }
    return fallback
}

func intValue(_ row: [String: Any], _ key: String, _ fallback: Int = 0) -> Int {
    if let x = row[key] as? Int {
        return x
    }
    if let x = row[key] as? Double {
        return Int(x)
    }
    return fallback
}

func diversifyRows(_ selected: [[String: Any]], limit: Int, sep: Double) -> [[String: Any]] {
    let wanted = max(1, limit)
    if sep <= 0 {
        return Array(selected.prefix(wanted))
    }
    var chosen: [[String: Any]] = []
    var chosenPoints: [[ComplexD]] = []
    var deferred: [[String: Any]] = []
    for row in selected {
        let z = rowPoint(row)
        let zn = max(1.0, vectorNorm(z))
        var near = false
        for prev in chosenPoints {
            let dist = vectorDistance(z, prev) / max(zn, vectorNorm(prev), 1.0)
            if dist < sep {
                near = true
                break
            }
        }
        if near {
            deferred.append(row)
            continue
        }
        chosen.append(row)
        chosenPoints.append(z)
        if chosen.count >= wanted {
            return chosen
        }
    }
    var seen = Set<Int>(chosen.map { intValue($0, "index", -1) })
    for row in deferred {
        let idx = intValue(row, "index", -1)
        if seen.contains(idx) {
            continue
        }
        chosen.append(row)
        seen.insert(idx)
        if chosen.count >= wanted {
            break
        }
    }
    return chosen
}

func finiteResidual(_ ctx: SolveContext, _ y: [ComplexD]) -> Double {
    let r = ctx.residual(y)
    return r.isFinite ? r : Double.infinity
}

func probeEndpoint(
    ctx: SolveContext,
    y: [ComplexD],
    residual: Double,
    prevDelta: [ComplexD]?,
    epoch: Int,
    directionSeed: UInt64,
    probeScale: Double,
    probeCandidates: Int,
    probeRadii: [Double],
    includeSelf: Bool
) throws -> ([ComplexD], [String: Any]) {
    let n = y.count
    let ynorm = max(1.0, vectorNorm(y))
    let radii = probeRadii.filter { $0 >= 0.0 }
    let usableRadii = radii.isEmpty ? [1.0] : radii
    var candidates: [(String, [ComplexD])] = []
    if includeSelf {
        candidates.append(("self", y))
    }
    if let prev = prevDelta {
        let pdn = max(1.0e-300, vectorNorm(prev))
        let base = scaleVector(prev, min(max(pdn, probeScale * ynorm), 2.5 * ynorm) / pdn)
        candidates.append(("inertial", addVectors(y, base)))
    }
    let budget = max(1, probeCandidates)
    var k = 0
    while candidates.count < budget {
        let rad = probeScale * ynorm * usableRadii[k % usableRadii.count]
        var qdir = rawDirectionD(n: n, trial: UInt64(epoch + 1) &* 104729 &+ UInt64(k + 1) &* 7919 &+ directionSeed, seed: directionSeed ^ (0x116116 &+ UInt64(17 * k)), normalize: true)
        let qnorm = max(1.0e-300, vectorNorm(qdir))
        qdir = scaleVector(qdir, sqrt(Double(max(1, n))) / qnorm)
        let ph = phaseD(0.6180339887498948 * Double(epoch + 1) + 2.399963229728653 * Double(k + 1))
        var step: [ComplexD] = []
        step.reserveCapacity(n)
        for j in 0..<n {
            var s = rad * (ph * qdir[j])
            if ynorm > 0.0 {
                s = s + (0.12 * rad / ynorm) * (y[j] * phaseD(0.38196601125 * Double(k + 1)))
            }
            let tiny = 1.0e-12 * ynorm
            if sqrt(cabs2(s)) < tiny {
                s = s + tiny * phaseD(0.17 + Double(j) + Double(epoch) + Double(k))
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
        throw NSError(domain: "StartSelect", code: 70, userInfo: [NSLocalizedDescriptionKey: "no finite probe"])
    }
    return (
        selected,
        [
            "probe_mode": "swift-double-theorem-guided-residual-min",
            "probe_name": bestName,
            "probe_candidates": min(budget, candidates.count),
            "probe_evals": evals,
            "probe_residual": bestRes,
            "probe_distance": bestDistance,
            "probe_improvement_proxy": (bestRes.isFinite && bestRes > 0.0 && residual.isFinite) ? residual / bestRes : NSNull(),
            "probe_self_enabled": includeSelf
        ]
    )
}

func pandrosionCorrectorSwift(
    ctx: SolveContext,
    y0: [ComplexD],
    maxEpochs: Int,
    tol: Double,
    accept: Double,
    lineSearch: Int,
    probeScale: Double,
    directionSeed: UInt64,
    probeCandidates: Int,
    probeRadii: [Double],
    includeSelf: Bool
) -> [String: Any] {
    let t0 = CFAbsoluteTimeGetCurrent()
    var y = y0
    var bestY = y
    var bestR = finiteResidual(ctx, y)
    var ok = false
    var status = "started"
    var epochs = 0
    var prevDelta: [ComplexD]? = nil
    var lastProbeMeta: [String: Any] = [:]
    var totalProbeEvals = 0

    for ep in 0..<max(1, maxEpochs) {
        let f = ctx.eval(y)
        let r = sqrt(f.reduce(0.0) { $0 + cabs2($1) })
        if r.isFinite && r < bestR {
            bestR = r
            bestY = y
        }
        if r <= max(tol, accept) && (accept <= 0.0 || r < accept) {
            ok = true
            status = "converged"
            break
        }
        let b: [ComplexD]
        do {
            let probe = try probeEndpoint(
                ctx: ctx,
                y: y,
                residual: r,
                prevDelta: prevDelta,
                epoch: ep,
                directionSeed: directionSeed,
                probeScale: probeScale,
                probeCandidates: probeCandidates,
                probeRadii: probeRadii,
                includeSelf: includeSelf
            )
            b = probe.0
            lastProbeMeta = probe.1
            totalProbeEvals += intValue(probe.1, "probe_evals", 0)
        } catch {
            status = "probe-error:\(type(of: error))"
            break
        }
        let q = ctx.slopeMatrix(y, b)
        guard var delta = solve2x2(q, f) else {
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
    if finalR <= max(tol, accept) && (accept <= 0.0 || finalR < accept) {
        ok = true
        status = "converged"
    }
    var out: [String: Any] = [
        "accepted": accept <= 0.0 ? ok : (finalR.isFinite && finalR < accept),
        "ok": ok,
        "status": status,
        "epochs": epochs,
        "residual": finalR,
        "y": bestY,
        "seconds": CFAbsoluteTimeGetCurrent() - t0,
        "corrector": "swift-double-probe-aware-pure-pandrosion-exact-telescopic-slope",
        "probe_total_evals": totalProbeEvals
    ]
    for (k, v) in lastProbeMeta {
        out[k] = v
    }
    return out
}

func startoptSwift(
    ctx: SolveContext,
    y0: [ComplexD],
    trial: Int,
    seed: UInt64,
    steps: Int,
    candidates: Int,
    gains: [Double],
    microEpochs: Int,
    job: Solve2Job
) -> ([ComplexD], [String: Any]) {
    var best = y0
    var bestR = finiteResidual(ctx, best)
    let initial = bestR
    var evals = 1
    var microTotal = 0
    var chosenGain = 1.0
    let usableGains = gains.isEmpty ? [1.0] : gains

    if steps > 0 {
        for step in 0..<steps {
            let base = best
            var pool: [(Double, [ComplexD])] = [(1.0, base)]
            if candidates > 1 {
                for c in 0..<(candidates - 1) {
                    let gain = usableGains[(trial + 3 * step + c) % usableGains.count]
                    var cand = scaleVector(base, gain)
                    if c % 3 == 1 {
                        var pert: [ComplexD] = []
                        pert.reserveCapacity(cand.count)
                        for j in 0..<cand.count {
                            let h1 = seed &+ 0x5157A47 &+ UInt64(65537 * trial) &+ UInt64(4099 * c) &+ UInt64(193 * (j + 1))
                            let h2 = seed &+ 0x7150F7 &+ UInt64(104729 * trial) &+ UInt64(8191 * c) &+ UInt64(389 * (j + 1))
                            let ph = 0.23 * (2.0 * u01D(h1) - 1.0)
                            let amp = exp(0.28 * (2.0 * u01D(h2) - 1.0))
                            pert.append(cand[j] * amp * phaseD(ph))
                        }
                        cand = pert
                    } else if c % 3 == 2 {
                        let fresh = rawDirectionD(n: base.count, trial: UInt64(trial + 31 * (step + 1) + c), seed: seed, normalize: true)
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
                if microEpochs > 0 {
                    let loc = pandrosionCorrectorSwift(
                        ctx: ctx,
                        y0: yscore,
                        maxEpochs: microEpochs,
                        tol: 1.0e-12,
                        accept: 0.0,
                        lineSearch: 6,
                        probeScale: job.probeScale,
                        directionSeed: seed,
                        probeCandidates: job.probeCandidates,
                        probeRadii: job.probeRadii,
                        includeSelf: job.probeSelf
                    )
                    microTotal += intValue(loc, "epochs", 0)
                    if let yy = loc["y"] as? [ComplexD] {
                        let rr = doubleValue(loc, "residual", Double.infinity)
                        if rr < r {
                            yscore = yy
                            r = rr
                        }
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
            "startopt_enabled": steps > 0,
            "startopt_r0": initial,
            "startopt_r1": bestR,
            "startopt_ratio": (initial.isFinite && bestR.isFinite && initial > 0.0) ? bestR / initial : NSNull(),
            "startopt_steps": max(0, steps),
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

func polishedRows(_ result: [String: Any], originals: [[String: Any]]) -> [[String: Any]] {
    guard let selected = result["selected"] as? [[String: Any]] else {
        return originals
    }
    var rows: [[String: Any]] = []
    rows.reserveCapacity(selected.count)
    for row in selected {
        let localIndex = intValue(row, "index", -1)
        if localIndex < 0 || localIndex >= originals.count {
            continue
        }
        var out = originals[localIndex]
        out["pre_polish_rank"] = intValue(out, "rank", -1)
        out["pre_polish_residual"] = doubleValue(out, "residual", Double.infinity)
        out["polish2_rank"] = intValue(row, "rank", -1)
        out["polish2_input_index"] = localIndex
        out["polish2_residual"] = doubleValue(row, "residual", Double.infinity)
        out["rank"] = intValue(row, "rank", rows.count)
        out["index"] = intValue(out, "index", localIndex)
        out["residual"] = doubleValue(row, "residual", Double.infinity)
        out["re"] = row["re"] as? [Double] ?? []
        out["im"] = row["im"] as? [Double] ?? []
        out["metal_polish2"] = true
        rows.append(out)
    }
    return rows
}

func refineRowsSwift(ctx: SolveContext, selected: [[String: Any]], job: Solve2Job) -> [String: Any] {
    var roots: [[String: Any]] = []
    var trials: [[String: Any]] = []
    var failures = 0
    var duplicates = 0
    let t0 = CFAbsoluteTimeGetCurrent()

    for cand in selected {
        if roots.count >= job.count {
            break
        }
        let trial = intValue(cand, "trial", intValue(cand, "index", 0))
        let yRaw = rowPoint(cand)
        if yRaw.count != 2 {
            failures += 1
            continue
        }
        let start = startoptSwift(
            ctx: ctx,
            y0: yRaw,
            trial: trial,
            seed: job.seed &+ 0x112555,
            steps: job.startoptSteps,
            candidates: job.startoptCandidates,
            gains: job.startoptGains,
            microEpochs: job.startoptMicroEpochs,
            job: job
        )
        let y0 = start.0
        let loc = pandrosionCorrectorSwift(
            ctx: ctx,
            y0: y0,
            maxEpochs: job.epochs,
            tol: job.tol,
            accept: job.accept,
            lineSearch: job.lineSearch,
            probeScale: job.probeScale,
            directionSeed: job.seed &+ UInt64(7919 * trial),
            probeCandidates: job.probeCandidates,
            probeRadii: job.probeRadii,
            includeSelf: job.probeSelf
        )
        let y = loc["y"] as? [ComplexD] ?? y0
        let rOrig = finiteResidual(ctx, y)
        let accepted = rOrig.isFinite && rOrig < job.accept
        var rec: [String: Any] = [
            "trial": trial,
            "selector_rank": intValue(cand, "rank", -1),
            "selector_residual": doubleValue(cand, "residual", Double.infinity),
            "accepted": accepted,
            "status": loc["status"] as? String ?? "",
            "r1": rOrig,
            "epochs": intValue(loc, "epochs", 0),
            "seconds": doubleValue(loc, "seconds", 0.0),
            "probe_total_evals": intValue(loc, "probe_total_evals", 0)
        ]
        for (k, v) in start.1 {
            rec[k] = v
        }
        if !accepted {
            failures += 1
            trials.append(rec)
            continue
        }
        if let dup = clusterIndex(roots, y, sep: job.clusterSep) {
            duplicates += 1
            rec["status"] = "duplicate"
            rec["cluster"] = dup
            trials.append(rec)
            continue
        }
        let rid = roots.count
        roots.append([
            "id": rid,
            "source": "128-swift-metal-selector-polish2/swift-double-refine",
            "trial": trial,
            "selector_rank": intValue(cand, "rank", -1),
            "selector_residual": doubleValue(cand, "residual", Double.infinity),
            "residual": rOrig,
            "epochs": intValue(loc, "epochs", 0),
            "seconds": doubleValue(loc, "seconds", 0.0),
            "z": rootJSON(y),
            "y": rootJSON(y)
        ])
        rec["status"] = "new-root"
        rec["root_id"] = rid
        trials.append(rec)
    }

    return [
        "roots": roots,
        "trials": trials,
        "summary": [
            "requested_roots": job.count,
            "unique_roots": roots.count,
            "success": roots.count >= job.count,
            "trials_used": trials.count,
            "duplicates": duplicates,
            "failures": failures,
            "refine_seconds": CFAbsoluteTimeGetCurrent() - t0,
            "eval_stats": ctx.stats()
        ]
    ]
}

func runSolve2(_ job: Solve2Job, system: LoadedSystem) throws -> [String: Any] {
    if system.n != 2 || system.equations != 2 {
        throw NSError(domain: "StartSelect", code: 80, userInfo: [NSLocalizedDescriptionKey: "solve2 requires n=2 and equations=2"])
    }
    let totalStart = CFAbsoluteTimeGetCurrent()
    let startJob = StoredStartSelectJob(
        op: "start_select",
        candidates: job.candidates,
        targetCount: job.targetCount,
        topK: job.topK,
        seed: job.seed,
        powerCap: job.powerCap,
        powers: job.powers,
        angles: job.angles,
        radii: job.radii
    )
    let selector = try runStoredJob(startJob, system: system)
    let selectorRows = selector["selected"] as? [[String: Any]] ?? []
    var selected = diversifyRows(selectorRows, limit: job.refineTop, sep: job.diversitySep)
    let selectedBeforePolish = selected
    var polishResult: [String: Any] = [
        "enabled": false,
        "points": 0,
        "kernel_seconds": 0.0,
        "total_seconds": 0.0
    ]
    if job.metalPolish2 {
        var pointsRe: [Float] = []
        var pointsIm: [Float] = []
        for row in selected {
            for z in rowPoint(row) {
                pointsRe.append(Float(z.re))
                pointsIm.append(Float(z.im))
            }
        }
        let polishJob = Polish2Job(
            op: "polish2",
            points: selected.count,
            epochs: job.metalPolishEpochs,
            probeCandidates: job.metalPolishProbes,
            lineSearch: job.metalPolishLineSearch,
            seed: job.seed &+ 0x125000,
            probeScale: job.metalPolishProbeScale,
            pointsRe: pointsRe,
            pointsIm: pointsIm
        )
        polishResult = try runPolish2(polishJob, system: system)
        polishResult["enabled"] = true
        selected = diversifyRows(polishedRows(polishResult, originals: selected), limit: job.refineTop, sep: job.diversitySep)
    }

    let ctx = SolveContext(system: system)
    let refined = refineRowsSwift(ctx: ctx, selected: selected, job: job)
    let refinedSummary = refined["summary"] as? [String: Any] ?? [:]
    var summary = refinedSummary
    summary["generation_seconds"] = 0.0
    summary["metal_select_kernel_seconds"] = selector["kernel_seconds"] as? Double ?? 0.0
    summary["metal_select_process_seconds"] = selector["total_seconds"] as? Double ?? 0.0
    summary["metal_polish2_seconds"] = polishResult["total_seconds"] as? Double ?? 0.0
    summary["metal_polish2_kernel_seconds"] = polishResult["kernel_seconds"] as? Double ?? 0.0
    summary["total_seconds"] = CFAbsoluteTimeGetCurrent() - totalStart

    return [
        "script": "125_swift_metal_start_select_server.swift",
        "mode": "128-swift-metal-selector-polish2/swift-double-refine",
        "op": "solve2",
        "n": system.n,
        "equations": system.equations,
        "terms_per_poly": system.terms,
        "terms": system.terms * system.equations,
        "device": device.name,
        "selector": selector,
        "selected_before_polish_count": selectedBeforePolish.count,
        "selected_for_refine_count": selected.count,
        "polish2": polishResult,
        "roots": refined["roots"] as? [[String: Any]] ?? [],
        "trials": refined["trials"] as? [[String: Any]] ?? [],
        "summary": summary
    ]
}

let decoder = JSONDecoder()
var loadedSystem: LoadedSystem? = nil
while let line = readLine() {
    if line == "__quit__" {
        break
    }
    do {
        let data = Data(line.utf8)
        let raw = try JSONSerialization.jsonObject(with: data, options: []) as? [String: Any]
        let op = raw?["op"] as? String
        let payload: [String: Any]
        if op == "load" {
            let job = try decoder.decode(LoadSystemJob.self, from: data)
            loadedSystem = try loadSystem(job)
            payload = [
                "ok": true,
                "op": "load",
                "device": device.name,
                "n": job.n,
                "equations": job.equations,
                "terms": job.terms
            ]
        } else if op == "start_select" {
            guard let system = loadedSystem else {
                throw NSError(domain: "StartSelect", code: 40, userInfo: [NSLocalizedDescriptionKey: "no loaded system"])
            }
            let job = try decoder.decode(StoredStartSelectJob.self, from: data)
            payload = try runStoredJob(job, system: system)
        } else if op == "eval_points" {
            guard let system = loadedSystem else {
                throw NSError(domain: "StartSelect", code: 41, userInfo: [NSLocalizedDescriptionKey: "no loaded system"])
            }
            let job = try decoder.decode(EvalPointsJob.self, from: data)
            payload = try runEvalPoints(job, system: system)
        } else if op == "polish2" {
            guard let system = loadedSystem else {
                throw NSError(domain: "StartSelect", code: 42, userInfo: [NSLocalizedDescriptionKey: "no loaded system"])
            }
            let job = try decoder.decode(Polish2Job.self, from: data)
            payload = try runPolish2(job, system: system)
        } else if op == "solve2" {
            guard let system = loadedSystem else {
                throw NSError(domain: "StartSelect", code: 43, userInfo: [NSLocalizedDescriptionKey: "no loaded system"])
            }
            let job = try decoder.decode(Solve2Job.self, from: data)
            payload = try runSolve2(job, system: system)
        } else {
            let job = try decoder.decode(StartSelectJob.self, from: data)
            payload = try runJob(job)
        }
        let json = try JSONSerialization.data(withJSONObject: payload, options: [.sortedKeys])
        FileHandle.standardOutput.write(json)
        FileHandle.standardOutput.write("\n".data(using: .utf8)!)
        fflush(stdout)
    } catch {
        let payload: [String: Any] = ["error": "\(error)"]
        let json = try JSONSerialization.data(withJSONObject: payload, options: [.sortedKeys])
        FileHandle.standardOutput.write(json)
        FileHandle.standardOutput.write("\n".data(using: .utf8)!)
        fflush(stdout)
    }
}
