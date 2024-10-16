use super::util;
use rand::seq::SliceRandom;
use rand::Rng;
use wgpu::util::DeviceExt;

#[repr(C)]
#[derive(Debug, Clone)]
pub struct Vertex {
    position: [f32; 2],
    colour: [f32; 3],
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

//  Algorithm to generate random convex polygons is converted to Rust from my
//  Lua implementation, which was converted from a Java implementation.
//
//  The original Java implementation is here:
//  https://cglab.ca/~sander/misc/ConvexGeneration/ValtrAlgorithm.java
//  The algorithm is explained in a blog post here:
//  https://cglab.ca/~sander/misc/ConvexGeneration/convex.html
//
//  And lists the following steps:
//  1. Generate two lists of random X and Y coordinates
//  2. Sort them
//  3. Isolate the extreme points
//  4. Randomly divide the interior points into two chains
//  5. Extract the vector components
//  6. Randomly pair up the X- and Y-components
//  7. Combine the paired up components into vectors
//  8. Sort the vectors by angle
//  9. Lay them end-to-end to form a polygon
//  10. Move the polygon to the original min and max coordinates

pub fn generate_random_convex_polygon(
    device: &wgpu::Device,
    n: usize,
    colour: [f32; 3],
) -> wgpu::Buffer {
    let mut rng = rand::thread_rng();

    // Step 1: Generate two lists of random X and Y coordinates
    let mut x_pool: Vec<f32> = (0..n).map(|_| rng.gen::<f32>()).collect();
    let mut y_pool: Vec<f32> = (0..n).map(|_| rng.gen::<f32>()).collect();

    // Step 2: Sort them
    x_pool.sort_by(|a, b| a.partial_cmp(b).unwrap());
    y_pool.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Step 3: Isolate the extreme points
    let min_x = x_pool[0];
    let max_x = x_pool[n - 1];
    let min_y = y_pool[0];
    let max_y = y_pool[n - 1];

    // Step 4: Randomly divide the interior points into two chains
    let mut x_vec = Vec::with_capacity(n);
    let mut last_top = min_x;
    let mut last_bot = min_x;

    for i in 1..n - 1 {
        let x = x_pool[i];
        if rng.gen_bool(0.5) {
            x_vec.push(x - last_top);
            last_top = x;
        } else {
            x_vec.push(last_bot - x);
            last_bot = x;
        }
    }

    // Add the extreme points to the chains
    x_vec.push(max_x - last_top);
    x_vec.push(last_bot - max_x);

    let mut y_vec = Vec::with_capacity(n);
    let mut last_left = min_y;
    let mut last_right = min_y;

    for i in 1..n - 1 {
        let y = y_pool[i];
        if rng.gen_bool(0.5) {
            y_vec.push(y - last_left);
            last_left = y;
        } else {
            y_vec.push(last_right - y);
            last_right = y;
        }
    }

    // Add the extreme points to the chains
    y_vec.push(max_y - last_left);
    y_vec.push(last_right - max_y);

    // Step 6: Shuffle the Y-components
    y_vec.shuffle(&mut rng);

    // Step 7: Combine the paired up components into vectors
    let mut vec = Vec::with_capacity(n);
    for i in 0..n {
        vec.push((x_vec[i], y_vec[i]));
    }

    // Step 8: Sort the vectors by angle
    vec.sort_by(|a, b| {
        let angle_a = (a.1).atan2(a.0);
        let angle_b = (b.1).atan2(b.0);
        angle_a.partial_cmp(&angle_b).unwrap()
    });

    // Step 9: Lay them end-to-end to form a polygon
    let mut vertices = Vec::with_capacity(n);
    let mut x = 0.0;
    let mut y = 0.0;
    let mut min_polygon_x: f32 = 0.0;
    let mut min_polygon_y: f32 = 0.0;

    for point in &vec {
        vertices.push(Vertex {
            position: [x as f32, y as f32],
            colour,
        });

        x += point.0;
        y += point.1;

        min_polygon_x = min_polygon_x.min(x);
        min_polygon_y = min_polygon_y.min(y);
    }

    // Step 10: Move the polygon
    let center_vertex = calculate_center(&vertices);
    let x_shift = -center_vertex.position[0];
    let y_shift = -center_vertex.position[1];

    let size = 600.0;

    for vertex in &mut vertices {
        vertex.position[0] += x_shift as f32;
        vertex.position[1] += y_shift as f32;

        vertex.position[0] *= size;
        vertex.position[1] *= size;
    }

    let triangulated_mesh = triangulate_polygon(&vertices);
    // println!("{:#?}", &triangulated_mesh);
    // println!(
    //     "Number of vertices after triangulation: {}",
    //     triangulated_mesh.len()
    // );

    let vertex_data_size =
        std::mem::size_of::<Vertex>() * triangulated_mesh.len();
    let vertex_data_ptr = triangulated_mesh.as_ptr() as *const u8;
    let vertex_data_slice = unsafe {
        std::slice::from_raw_parts(vertex_data_ptr, vertex_data_size)
    };

    let buffer_descriptor = wgpu::util::BufferInitDescriptor {
        label: Some("Polygon"),
        contents: vertex_data_slice,
        usage: wgpu::BufferUsages::VERTEX,
    };

    let buffer = device.create_buffer_init(&buffer_descriptor);

    println!("Created buffer with size: {} bytes", buffer.size());
    println!("Expected size: {} bytes", vertex_data_size);

    // println!("Created buffer with size: {} bytes", buffer.size());

    buffer
}

fn triangulate_polygon(outline: &[Vertex]) -> Vec<Vertex> {
    if outline.len() < 3 {
        return Vec::new(); // Not enough vertices to form a polygon
    }

    let mut triangulated = Vec::with_capacity(3 * (outline.len() - 2));

    // Calculate the center point
    let center = calculate_center(outline);

    // Triangulate by creating triangles from the center to each edge
    for i in 1..outline.len() - 1 {
        triangulated.push(center.clone());
        triangulated.push(outline[i].clone());
        triangulated.push(outline[i + 1].clone());
    }

    // Add the last triangle to close the fan
    triangulated.push(center.clone());
    triangulated.push(outline[outline.len() - 1].clone());
    triangulated.push(outline[1].clone());

    triangulated
}

fn calculate_center(vertices: &[Vertex]) -> Vertex {
    let count = vertices.len() as f32;
    let sum = vertices.iter().fold([0.0f32; 2], |acc, v| {
        [acc[0] + v.position[0], acc[1] + v.position[1]]
    });

    Vertex {
        position: [sum[0] / count, sum[1] / count],
        colour: [1.0, 1.0, 1.0], // I might want to adjust this
    }
}
