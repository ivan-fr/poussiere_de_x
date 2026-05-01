import Foundation
import Metal

struct StartSelectInput: Decodable {
    let n: Int
    let equations: Int
    let terms: Int
    let candidates: Int
    let targetCount: Int
    let seed: UInt64
    let powerCap: Float
    let exps: [Int32]
    let coeffRe: [Float]
    let coeffIm: [Float]
    let powers: [Float]
    let angles: [Float]
    let radii: [Float]
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

struct Complex2 {
    var re: Float
    var im: Float
}

struct Candidate {
    let index: Int
    let trial: Int
    let residual: Double
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
"""

if CommandLine.arguments.count < 4 {
    fail("usage: 124_swift_metal_start_select <input.json> <top_k> <output.json>")
}

let inputPath = CommandLine.arguments[1]
let topK = max(1, Int(CommandLine.arguments[2]) ?? 1)
let outputPath = CommandLine.arguments[3]
let totalStart = CFAbsoluteTimeGetCurrent()

let inputData = try Data(contentsOf: URL(fileURLWithPath: inputPath))
let input = try JSONDecoder().decode(StartSelectInput.self, from: inputData)

if input.n > 8 {
    fail("n > 8 is not supported by this Metal prototype")
}
if input.exps.count != input.terms * input.n {
    fail("bad exps length")
}
if input.coeffRe.count != input.equations * input.terms || input.coeffIm.count != input.equations * input.terms {
    fail("bad coeff length")
}
if input.powers.isEmpty || input.angles.isEmpty || input.radii.isEmpty {
    fail("powers/angles/radii must be nonempty")
}

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

let coeff = zip(input.coeffRe, input.coeffIm).map { Complex2(re: $0.0, im: $0.1) }
let outCount = input.candidates * input.equations
let pointsCount = input.candidates * input.n

guard let expsBuffer = device.makeBuffer(bytes: input.exps, length: input.exps.count * MemoryLayout<Int32>.stride, options: .storageModeShared),
      let coeffBuffer = device.makeBuffer(bytes: coeff, length: coeff.count * MemoryLayout<Complex2>.stride, options: .storageModeShared),
      let powersBuffer = device.makeBuffer(bytes: input.powers, length: input.powers.count * MemoryLayout<Float>.stride, options: .storageModeShared),
      let anglesBuffer = device.makeBuffer(bytes: input.angles, length: input.angles.count * MemoryLayout<Float>.stride, options: .storageModeShared),
      let radiiBuffer = device.makeBuffer(bytes: input.radii, length: input.radii.count * MemoryLayout<Float>.stride, options: .storageModeShared),
      let outBuffer = device.makeBuffer(length: outCount * MemoryLayout<Complex2>.stride, options: .storageModeShared),
      let pointsBuffer = device.makeBuffer(length: pointsCount * MemoryLayout<Complex2>.stride, options: .storageModeShared) else {
    fail("cannot create Metal buffers")
}

var params = Params(
    n: UInt32(input.n),
    equations: UInt32(input.equations),
    terms: UInt32(input.terms),
    candidates: UInt32(input.candidates),
    targetCount: UInt32(input.targetCount),
    powersCount: UInt32(input.powers.count),
    anglesCount: UInt32(input.angles.count),
    radiiCount: UInt32(input.radii.count),
    seedLo: UInt32(input.seed & 0xFFFFFFFF),
    seedHi: UInt32((input.seed >> 32) & 0xFFFFFFFF),
    powerCap: input.powerCap
)
guard let paramsBuffer = device.makeBuffer(bytes: &params, length: MemoryLayout<Params>.stride, options: .storageModeShared) else {
    fail("cannot create params buffer")
}

let kernelStart = CFAbsoluteTimeGetCurrent()
guard let commandBuffer = queue.makeCommandBuffer(),
      let encoder = commandBuffer.makeComputeCommandEncoder() else {
    fail("cannot create Metal command buffer")
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
candidates.reserveCapacity(input.candidates)
for point in 0..<input.candidates {
    var norm = 0.0
    for eq in 0..<input.equations {
        let v = output[point * input.equations + eq]
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
let selected = Array(candidates.prefix(min(topK, candidates.count)))
let selectSeconds = CFAbsoluteTimeGetCurrent() - selectStart

var rows: [[String: Any]] = []
rows.reserveCapacity(selected.count)
for cand in selected {
    var zr: [Double] = []
    var zi: [Double] = []
    zr.reserveCapacity(input.n)
    zi.reserveCapacity(input.n)
    for j in 0..<input.n {
        let p = points[cand.index * input.n + j]
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

let totalSeconds = CFAbsoluteTimeGetCurrent() - totalStart
let payload: [String: Any] = [
    "device": device.name,
    "candidates": input.candidates,
    "terms": input.terms,
    "equations": input.equations,
    "top_k": topK,
    "selected": rows,
    "kernel_seconds": kernelSeconds,
    "select_seconds": selectSeconds,
    "total_seconds": totalSeconds
]
let json = try JSONSerialization.data(withJSONObject: payload, options: [.sortedKeys])
try json.write(to: URL(fileURLWithPath: outputPath))
FileHandle.standardOutput.write(json)
FileHandle.standardOutput.write("\n".data(using: .utf8)!)
