import Foundation
import Metal

struct SelectInput: Decodable {
    let n: Int
    let equations: Int
    let terms: Int
    let points: Int
    let exps: [Int32]
    let coeffRe: [Float]
    let coeffIm: [Float]
    let pointsRe: [Float]
    let pointsIm: [Float]
    let trialIds: [Int32]
}

struct Params {
    var n: UInt32
    var equations: UInt32
    var terms: UInt32
    var points: UInt32
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
    uint points;
};

static inline float2 cmul(float2 a, float2 b) {
    return float2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}

static inline float2 cpow_int(float2 z, int e) {
    float2 out = float2(1.0f, 0.0f);
    for (int k = 0; k < e; ++k) {
        out = cmul(out, z);
    }
    return out;
}

kernel void eval_dense_ks(
    device const int *exps [[buffer(0)]],
    device const float2 *coeff [[buffer(1)]],
    device const float2 *points [[buffer(2)]],
    device float2 *out [[buffer(3)]],
    constant Params &p [[buffer(4)]],
    uint gid [[thread_position_in_grid]]
) {
    uint total = p.points * p.equations;
    if (gid >= total) {
        return;
    }

    uint pointIndex = gid / p.equations;
    uint equation = gid - pointIndex * p.equations;
    float2 sum = float2(0.0f, 0.0f);

    for (uint t = 0; t < p.terms; ++t) {
        float2 mon = float2(1.0f, 0.0f);
        for (uint j = 0; j < p.n; ++j) {
            int e = exps[t * p.n + j];
            float2 z = points[pointIndex * p.n + j];
            mon = cmul(mon, cpow_int(z, e));
        }
        float2 c = coeff[equation * p.terms + t];
        sum += cmul(c, mon);
    }

    out[gid] = sum;
}
"""

if CommandLine.arguments.count < 4 {
    fail("usage: 123_swift_metal_select <input.json> <top_k> <output.json>")
}

let inputPath = CommandLine.arguments[1]
let topK = max(1, Int(CommandLine.arguments[2]) ?? 1)
let outputPath = CommandLine.arguments[3]
let totalStart = CFAbsoluteTimeGetCurrent()

let inputData = try Data(contentsOf: URL(fileURLWithPath: inputPath))
let input = try JSONDecoder().decode(SelectInput.self, from: inputData)

if input.exps.count != input.terms * input.n {
    fail("bad exps length")
}
if input.coeffRe.count != input.equations * input.terms || input.coeffIm.count != input.equations * input.terms {
    fail("bad coeff length")
}
if input.pointsRe.count != input.points * input.n || input.pointsIm.count != input.points * input.n {
    fail("bad points length")
}
if input.trialIds.count != input.points {
    fail("bad trialIds length")
}

guard let device = MTLCreateSystemDefaultDevice() else {
    fail("no Metal device")
}
guard let queue = device.makeCommandQueue() else {
    fail("cannot create Metal command queue")
}

let library = try device.makeLibrary(source: metalSource, options: nil)
guard let function = library.makeFunction(name: "eval_dense_ks") else {
    fail("missing eval_dense_ks function")
}
let pipeline = try device.makeComputePipelineState(function: function)

let coeff = zip(input.coeffRe, input.coeffIm).map { Complex2(re: $0.0, im: $0.1) }
let points = zip(input.pointsRe, input.pointsIm).map { Complex2(re: $0.0, im: $0.1) }
let outCount = input.points * input.equations

guard let expsBuffer = device.makeBuffer(
    bytes: input.exps,
    length: input.exps.count * MemoryLayout<Int32>.stride,
    options: .storageModeShared
) else {
    fail("cannot create exps buffer")
}
guard let coeffBuffer = device.makeBuffer(
    bytes: coeff,
    length: coeff.count * MemoryLayout<Complex2>.stride,
    options: .storageModeShared
) else {
    fail("cannot create coeff buffer")
}
guard let pointsBuffer = device.makeBuffer(
    bytes: points,
    length: points.count * MemoryLayout<Complex2>.stride,
    options: .storageModeShared
) else {
    fail("cannot create points buffer")
}
guard let outBuffer = device.makeBuffer(
    length: outCount * MemoryLayout<Complex2>.stride,
    options: .storageModeShared
) else {
    fail("cannot create output buffer")
}

var params = Params(
    n: UInt32(input.n),
    equations: UInt32(input.equations),
    terms: UInt32(input.terms),
    points: UInt32(input.points)
)
guard let paramsBuffer = device.makeBuffer(
    bytes: &params,
    length: MemoryLayout<Params>.stride,
    options: .storageModeShared
) else {
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
encoder.setBuffer(pointsBuffer, offset: 0, index: 2)
encoder.setBuffer(outBuffer, offset: 0, index: 3)
encoder.setBuffer(paramsBuffer, offset: 0, index: 4)

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
var candidates: [Candidate] = []
candidates.reserveCapacity(input.points)
for point in 0..<input.points {
    var norm = 0.0
    for eq in 0..<input.equations {
        let v = output[point * input.equations + eq]
        let re = Double(v.re)
        let im = Double(v.im)
        norm += re * re + im * im
    }
    let residual = sqrt(norm)
    if residual.isFinite {
        candidates.append(Candidate(index: point, trial: Int(input.trialIds[point]), residual: residual))
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
    "points": input.points,
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
