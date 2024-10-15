// Define the uniform for window size
struct WindowSize {
    windowSize: vec2<f32>,
}

// Bind the uniform to the shader
@group(0) @binding(0) var<uniform> window_size: WindowSize;

struct Vertex {
    @location(0) position: vec3<f32>,
    @location(1) colour: vec3<f32>,
}

struct VertexPayload {
    @builtin(position) position: vec4<f32>,
    @location(0) colour: vec3<f32>,
}

@vertex
fn vs_main(vertex: Vertex) -> VertexPayload {
    var out: VertexPayload;

    // Transform the vertex position from normalized space to pixel space
    let ndc_position = vertex.position.xy;

    // Convert NDC (-1 to 1) to screen space (0 to windowSize)
    let pixel_position = ndc_position / window_size.windowSize;

    // // Map pixel position back to NDC space (-1 to 1)
    // let adjusted_position = (pixel_position / window_size.windowSize) * 2.0 - vec2<f32>(1.0);

    out.position = vec4<f32>(pixel_position, vertex.position.z, 1.0);

    out.colour = vertex.colour;

    return out;
}

@fragment
fn fs_main(input: VertexPayload) -> @location(0) vec4<f32> {
    return vec4<f32>(input.colour, 1.0);
}
