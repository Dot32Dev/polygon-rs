use std::env::current_dir;
use std::fs;

pub struct Pipleline<'a> {
    shader_filename: String,
    vertex_entry: String,
    fragment_entry: String,
    pixel_format: wgpu::TextureFormat,
    vertex_buffer_layouts: Vec<wgpu::VertexBufferLayout<'a>>,
    bind_group_layouts: Vec<&'a wgpu::BindGroupLayout>,
}

impl<'a> Pipleline<'a> {
    pub fn new(
        shader_filename: &str,
        vertex_entry: &str,
        fragment_entry: &str,
        pixel_format: wgpu::TextureFormat,
    ) -> Self {
        Pipleline {
            shader_filename: shader_filename.to_owned(),
            vertex_entry: vertex_entry.to_owned(),
            fragment_entry: fragment_entry.to_owned(),
            pixel_format,
            vertex_buffer_layouts: Vec::new(),
            bind_group_layouts: Vec::new(),
        }
    }

    fn reset(&mut self) {
        self.vertex_buffer_layouts.clear();
        self.bind_group_layouts.clear();
    }

    pub fn add_vertex_buffer_layout(
        &mut self,
        layout: wgpu::VertexBufferLayout<'static>,
    ) {
        self.vertex_buffer_layouts.push(layout);
    }

    pub fn add_bind_group_layout(&mut self, layout: &'a wgpu::BindGroupLayout) {
        self.bind_group_layouts.push(layout);
    }

    pub fn build(&mut self, device: &wgpu::Device) -> wgpu::RenderPipeline {
        let mut filepath = current_dir().unwrap();
        filepath.push("src/shaders/");
        filepath.push(self.shader_filename.as_str());
        let filepath = filepath.into_os_string().into_string().unwrap();
        println!("{}", filepath);
        let source_code =
            fs::read_to_string(filepath).expect("Cant read shader file");

        let shader_module =
            device.create_shader_module(wgpu::ShaderModuleDescriptor {
                label: Some("Shader module"),
                source: wgpu::ShaderSource::Wgsl(source_code.into()),
            });

        let pipleline_layout =
            device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
                label: Some("Render Pipleline layout"),
                bind_group_layouts: &self.bind_group_layouts,
                push_constant_ranges: &[],
            });

        let render_targets = [Some(wgpu::ColorTargetState {
            format: self.pixel_format,
            blend: Some(wgpu::BlendState::REPLACE),
            write_mask: wgpu::ColorWrites::ALL,
        })];

        let render_pipleline_descriptor = wgpu::RenderPipelineDescriptor {
            label: Some("Render pipeline"),
            layout: Some(&pipleline_layout),
            vertex: wgpu::VertexState {
                module: &shader_module,
                entry_point: &self.vertex_entry,
                buffers: &self.vertex_buffer_layouts,
            },
            primitive: wgpu::PrimitiveState {
                topology: wgpu::PrimitiveTopology::TriangleList,
                strip_index_format: None,
                front_face: wgpu::FrontFace::Ccw,
                cull_mode: Some(wgpu::Face::Back),
                polygon_mode: wgpu::PolygonMode::Fill,
                unclipped_depth: false,
                conservative: false,
            },
            fragment: Some(wgpu::FragmentState {
                module: &shader_module,
                entry_point: &self.fragment_entry,
                targets: &render_targets,
            }),
            depth_stencil: None,
            multisample: wgpu::MultisampleState {
                count: 1,
                mask: !0,
                alpha_to_coverage_enabled: false,
            },
            multiview: None,
        };

        device.create_render_pipeline(&render_pipleline_descriptor)
    }
}
