use super::util;
// use rand::seq::SliceRandom;
// use rand::Rng;
use wgpu::util::DeviceExt;

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Vertex {
    pub position: [f32; 2],
    pub colour: [f32; 3],
}

impl Vertex {
    pub fn get_layout() -> wgpu::VertexBufferLayout<'static> {
        // No idea why this has to be a const
        const ATTRIBUTES: [wgpu::VertexAttribute; 2] =
            wgpu::vertex_attr_array![0 => Float32x2, 1 => Float32x3];

        wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<Vertex>() as u64,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &ATTRIBUTES,
        }
    }
}

#[rustfmt::skip] // For nicer formatting of the vertices array
pub fn make_triangle(device: &wgpu::Device) -> wgpu::Buffer {
    let size = 1200.0; // Size
    let vertices: [Vertex; 3] = [
        Vertex {position: [1.0*size,  -1.0*size], colour: [1.0, 0.0, 0.0]},
        Vertex {position: [0.0*size,   1.0*size], colour: [0.0, 1.0, 0.0]},
        Vertex {position: [-1.0*size, -1.0*size], colour: [0.0, 0.0, 1.0]},
    ];

    let buffer_descriptor = wgpu::util::BufferInitDescriptor {
        label: Some("Vertex buffer"),
        contents: unsafe { util::any_as_u8_slice(&vertices) },
        usage: wgpu::BufferUsages::VERTEX,
    };

    device.create_buffer_init(&buffer_descriptor)
}

pub fn from_vertices(
    device: &wgpu::Device,
    vertices: &[Vertex],
) -> wgpu::Buffer {
    let vertex_data_size = std::mem::size_of::<Vertex>() * vertices.len();
    let vertex_data_ptr = vertices.as_ptr() as *const u8;
    let vertex_data_slice = unsafe {
        std::slice::from_raw_parts(vertex_data_ptr, vertex_data_size)
    };

    let buffer_descriptor = wgpu::util::BufferInitDescriptor {
        label: Some("Vertex buffer"),
        contents: vertex_data_slice,
        usage: wgpu::BufferUsages::VERTEX,
    };

    device.create_buffer_init(&buffer_descriptor)
}
