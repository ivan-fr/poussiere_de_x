import Foundation
import Metal

struct BenchInput: Decodable {
    let n: Int
    let equations: Int
    let terms: Int
    let points: Int
    let exps: [Int32]
    let coeffRe: [Float]
    let coeffIm: [Float]
    let pointsRe: [Float]
    let pointsIm: [Float]
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

if CommandLine.arguments.count < 3 {
    fail("usage: 122_swift_metal_eval <input.json> <loops>")
}

let inputPath = CommandLine.arguments[1]
let loops = max(1, Int(CommandLine.arguments[2]) ?? 1)
let totalStart = CFAbsoluteTimeGetCurrent()

let inputData = try Data(contentsOf: URL(fileURLWithPath: inputPath))
let input = try JSONDecoder().decode(BenchInput.self, from: inputData)

if input.exps.count != input.terms * input.n {
    fail("bad exps length")
}
if input.coeffRe.count != input.equations * input.terms || input.coeffIm.count != input.equations * input.terms {
    fail("bad coeff length")
}
if input.pointsRe.count != input.points * input.n || input.pointsIm.count != input.points * input.n {
    fail("bad points length")
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

func dispatchOnce() {
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
}

dispatchOnce()
let kernelStart = CFAbsoluteTimeGetCurrent()
for _ in 0..<loops {
    dispatchOnce()
}
let kernelSeconds = CFAbsoluteTimeGetCurrent() - kernelStart

let output = outBuffer.contents().bindMemory(to: Complex2.self, capacity: outCount)
var checksumRe = 0.0
var checksumIm = 0.0
var checksumNorm = 0.0
for i in 0..<min(outCount, 4096) {
    let re = Double(output[i].re)
    let im = Double(output[i].im)
    checksumRe += re
    checksumIm += im
    checksumNorm += re * re + im * im
}

let totalSeconds = CFAbsoluteTimeGetCurrent() - totalStart
let payload: [String: Any] = [
    "device": device.name,
    "points": input.points,
    "equations": input.equations,
    "terms": input.terms,
    "loops": loops,
    "kernel_seconds": kernelSeconds,
    "total_seconds": totalSeconds,
    "point_evals": input.points * input.equations * loops,
    "point_evals_per_second": Double(input.points * input.equations * loops) / max(kernelSeconds, 1.0e-300),
    "checksum_re": checksumRe,
    "checksum_im": checksumIm,
    "checksum_norm": checksumNorm
]
let json = try JSONSerialization.data(withJSONObject: payload, options: [.sortedKeys])
FileHandle.standardOutput.write(json)
FileHandle.standardOutput.write("\n".data(using: .utf8)!)
