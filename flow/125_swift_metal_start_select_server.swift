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

struct Candidate {
    let index: Int
    let trial: Int
    let residual: Double
}

final class LoadedSystem {
    let n: Int
    let equations: Int
    let terms: Int
    let expsBuffer: MTLBuffer
    let coeffBuffer: MTLBuffer

    init(n: Int, equations: Int, terms: Int, expsBuffer: MTLBuffer, coeffBuffer: MTLBuffer) {
        self.n = n
        self.equations = equations
        self.terms = terms
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
    let coeff = zip(job.coeffRe, job.coeffIm).map { Complex2(re: $0.0, im: $0.1) }
    guard let expsBuffer = device.makeBuffer(bytes: job.exps, length: job.exps.count * MemoryLayout<Int32>.stride, options: .storageModeShared),
          let coeffBuffer = device.makeBuffer(bytes: coeff, length: coeff.count * MemoryLayout<Complex2>.stride, options: .storageModeShared) else {
        throw NSError(domain: "StartSelect", code: 23, userInfo: [NSLocalizedDescriptionKey: "cannot create loaded buffers"])
    }
    return LoadedSystem(n: job.n, equations: job.equations, terms: job.terms, expsBuffer: expsBuffer, coeffBuffer: coeffBuffer)
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
