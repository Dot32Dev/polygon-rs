use super::renderer_backend::mesh::Vertex;
use rand::seq::SliceRandom;
use rand::Rng;

// 704 bytes for a single 10 sided polygon :sob:
pub struct Polygon {
    pub vertices: Vec<Vertex>,
    pub vertex_count: usize,
    position: [f32; 2],
    linear_velocity: [f32; 2],
    angle: f32,
    angular_velocity: f32,
    area: f32,
    mass: f32,
    center: [f32; 2],
    density: f32,
    inertia: f32,
}

impl Polygon {
    pub fn generate_with_position(
        n: usize,
        colour: [f32; 3],
        position: [f32; 2],
    ) -> Self {
        let mut polygon = Self::generate(n, colour);
        polygon.position = position;
        polygon
    }

    // No real idea how this works. The original Java implementation is here:
    // https://cglab.ca/~sander/misc/ConvexGeneration/ValtrAlgorithm.java
    pub fn generate(n: usize, colour: [f32; 3]) -> Self {
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
        // For the X:
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

        // For the Y:
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

        // Step 5 is a lie.

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

        let vertex_count = vertices.len();
        // Define it
        let mut polygon = Self {
            vertices,
            vertex_count,
            position: [0.0; 2],
            linear_velocity: [0.0; 2],
            angle: 0.0,
            angular_velocity: 0.0,
            area: 0.0,
            mass: 0.0,
            center: [0.0; 2],
            density: 1.0,
            inertia: 0.0,
        };
        // Calculate its properties
        polygon.calculate_properties_from_edge_vertices();
        polygon.triangulate();

        let x_shift = -polygon.center[0];
        let y_shift = -polygon.center[1];

        let size = 600.0;

        // Shift it
        for vertex in &mut polygon.vertices {
            vertex.position[0] += x_shift as f32;
            vertex.position[1] += y_shift as f32;

            vertex.position[0] *= size;
            vertex.position[1] *= size;
        }

        // Now that we shifted all the vertices, this is the new center.
        polygon.center = [0.0, 0.0];

        polygon
    }

    // Should be ran before the mesh is triangulated.
    // Don't ask me how it works.
    fn calculate_properties_from_edge_vertices(&mut self) {
        self.area = 0.0;
        self.center = [0.0, 0.0];
        self.inertia = 0.0;

        for (vertex, next_vertex) in self
            .vertices
            .iter()
            .zip(self.vertices.iter().cycle().skip(1))
        {
            let cross = cross_2d(&vertex.position, &next_vertex.position);
            self.area += cross;

            // Center calculation
            self.center[0] +=
                (vertex.position[0] + next_vertex.position[0]) * cross;
            self.center[1] +=
                (vertex.position[1] + next_vertex.position[1]) * cross;

            // Inertia calculation
            // let factor = vertex.position[0] * vertex.position[0]
            //     + vertex.position[0] * next_vertex.position[0]
            //     + next_vertex.position[0] * next_vertex.position[0]
            //     + vertex.position[1] * vertex.position[1]
            //     + vertex.position[1] * next_vertex.position[1]
            //     + next_vertex.position[1] * next_vertex.position[1];
            let factor = dot_product(&vertex.position, &vertex.position)
                + dot_product(&vertex.position, &next_vertex.position)
                + dot_product(&next_vertex.position, &next_vertex.position);

            self.inertia += cross * factor;
        }

        // Finalize calculations
        self.area = (self.area * 0.5).abs();
        self.mass = self.area * self.density;
        self.inertia = (self.density / 12.0) * self.inertia.abs();

        // Correct the center calculation
        self.center[0] /= 6.0 * self.area;
        self.center[1] /= 6.0 * self.area;
    }

    // Returns the vertices transformed by the position and angle
    pub fn transformed_vertices(&self) -> Vec<Vertex> {
        let cos_theta = self.angle.cos();
        let sin_theta = self.angle.sin();

        let matrix = [
            [cos_theta, -sin_theta, self.position[0]],
            [sin_theta, cos_theta, self.position[1]],
            [0.0, 0.0, 1.0],
        ];

        let mut transform_vertices = self.vertices.clone();

        for vertex in &mut transform_vertices {
            vertex.position[0] = vertex.position[0] * matrix[0][0]
                + vertex.position[1] * matrix[0][1]
                + matrix[0][2];

            vertex.position[1] = vertex.position[0] * matrix[1][0]
                + vertex.position[1] * matrix[1][1]
                + matrix[1][2];
        }

        transform_vertices
    }

    fn triangulate(&mut self) {
        if self.vertices.len() < 3 {
            return; // Not enough vertices to form a polygon
        }

        let mut triangulated =
            Vec::with_capacity(3 * (self.vertices.len() - 2));

        // Calculate the center point
        let center = Vertex {
            position: [self.center[0], self.center[1]],
            colour: self.vertices[0].colour,
            // colour: [1.0, 1.0, 1.0],
        };

        // Triangulate by creating triangles from the center to each edge
        for i in 1..self.vertices.len() - 1 {
            triangulated.push(center.clone());
            triangulated.push(self.vertices[i].clone());
            triangulated.push(self.vertices[i + 1].clone());
        }

        // Add the last triangle to close the fan
        triangulated.push(center.clone());
        triangulated.push(self.vertices[self.vertices.len() - 1].clone());
        triangulated.push(self.vertices[1].clone());

        self.vertices = triangulated;
    }
}

// 2D cross product? idk man I dont make the rules.
fn cross_2d(a: &[f32; 2], b: &[f32; 2]) -> f32 {
    a[0] * b[1] - a[1] * b[0]
}

fn dot_product(a: &[f32; 2], b: &[f32; 2]) -> f32 {
    a[0] * b[0] + a[1] * b[1]
}
