use winit::{
    event::{Event, WindowEvent},
    event_loop::{ControlFlow, EventLoop},
    window::{Window, WindowBuilder},
};

mod renderer_backend;
use renderer_backend::mesh;
use renderer_backend::pipeline::Pipleline;
use renderer_backend::uniform;
use renderer_backend::util;

fn main() {
    // WGPU logs via this crate. We must call init to enable it
    env_logger::init();

    // Requires eventloop::EventLoop
    let event_loop = EventLoop::new().unwrap();
    // POLL Runs event loop continuously, even when no events were registered.
    // Good for games, bad for UI
    event_loop.set_control_flow(ControlFlow::Poll);

    // Input the event loop and this is where the window will send its events
    let window = WindowBuilder::new().build(&event_loop).unwrap();

    let mut state = pollster::block_on(State::new(&window));

    event_loop
        // Run supply's the event and event_loop_control_target inputs
        .run(move |event, event_loop_window_target| {
            match event {
                Event::WindowEvent {
                    event: WindowEvent::CloseRequested,
                    // We don't care what the window ID is, so we dont put it in
                    // the pattern
                    ..
                } => {
                    // The close button was pressed, close the window
                    event_loop_window_target.exit();
                }
                Event::WindowEvent {
                    event: WindowEvent::Resized(new_size),
                    ..
                } => {
                    state.resize(new_size);
                }
                Event::AboutToWait => {
                    // Update
                    state.update();
                }
                Event::WindowEvent {
                    event: WindowEvent::RedrawRequested,
                    ..
                } => {
                    // Draw
                    state.render().unwrap();
                }
                _ => (),
            }
        })
        .unwrap();
}

struct State<'a> {
    // The surface is the texture we draw to
    surface: wgpu::Surface<'a>,
    // The device is how we access the GPU
    device: wgpu::Device,
    queue: wgpu::Queue,
    config: wgpu::SurfaceConfiguration,
    size: winit::dpi::PhysicalSize<u32>,
    window: &'a Window,
    render_pipeline: wgpu::RenderPipeline,
    triangle_mesh: wgpu::Buffer,
    uniform_buffer: wgpu::Buffer,
    uniform_bind_group: wgpu::BindGroup,
}

impl<'a> State<'a> {
    // Creating some of the wgpu types requires async code
    async fn new(window: &'a Window) -> Self {
        let size = window.inner_size();

        // Handle to our GPU
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends: wgpu::Backends::all(),
            ..Default::default()
        });

        // The surface needs to live as long as the window that created it.
        // State owns the window, so this should be safe.
        // let surface =
        //     unsafe { instance.create_surface(window.clone()) }.unwrap();
        let surface = instance.create_surface(window).unwrap();

        // Create the adaptor!
        // We also use this to query the capabilities of the system.
        // We use this to create the device and queue
        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions {
                power_preference: wgpu::PowerPreference::default(),
                compatible_surface: Some(&surface),
                force_fallback_adapter: false,
            })
            .await
            .unwrap();

        let (device, queue) = adapter
            .request_device(
                &wgpu::DeviceDescriptor {
                    required_features: wgpu::Features::empty(),
                    // If we want to build to WebGL, this must be more
                    // restrictive
                    required_limits: wgpu::Limits::default(),
                    label: None,
                },
                None, // Trace parth ?
            )
            .await
            .unwrap();

        // Loop through all the possible surface formats and find one that
        // supports sRGB
        let surface_caps = surface.get_capabilities(&adapter);
        let surface_format = surface_caps
            .formats
            .iter()
            .copied()
            .filter(|f| f.is_srgb())
            .next()
            .unwrap_or(surface_caps.formats[0]);

        let config = wgpu::SurfaceConfiguration {
            // The only supported usage is RENDER_ATTACHMENT lmao
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            // The format we found that supports sRGB
            format: surface_format,
            width: size.width,
            height: size.height,
            present_mode: wgpu::PresentMode::Fifo,
            desired_maximum_frame_latency: 2,
            alpha_mode: surface_caps.alpha_modes[0],
            view_formats: vec![],
        };
        surface.configure(&device, &config);

        let triangle_mesh = mesh::make_triangle(&device);

        let mut render_pipeline_builder =
            Pipleline::new("shader.wgsl", "vs_main", "fs_main", config.format);

        render_pipeline_builder
            .add_vertex_buffer_layout(mesh::Vertex::get_layout());

        // let aspect_ratio = get_aspect_ratio(size.width, size.height);

        let uniforms =
            uniform::Uniforms::new(size.width as f32, size.height as f32);

        let uniform_buffer = uniforms.to_buffer(&device);

        let (uniform_bind_group, uniform_layout) =
            uniform::bind_group_and_layout(&device, &uniform_buffer);

        render_pipeline_builder.add_bind_group_layout(&uniform_layout);

        let render_pipeline = render_pipeline_builder.build(&device);

        Self {
            window,
            surface,
            device,
            queue,
            config,
            size,
            render_pipeline,
            triangle_mesh,
            uniform_buffer,
            uniform_bind_group,
        }
    }

    fn resize(&mut self, new_size: winit::dpi::PhysicalSize<u32>) {
        if new_size.width > 0 && new_size.height > 0 {
            self.size = new_size;
            self.config.width = new_size.width;
            self.config.height = new_size.height;
            self.surface.configure(&self.device, &self.config);

            let uniforms = uniform::Uniforms::new(
                new_size.width as f32,
                new_size.height as f32,
            );
            self.queue.write_buffer(&self.uniform_buffer, 0, unsafe {
                util::any_as_u8_slice(&uniforms)
            });
        }
    }

    fn update(&mut self) {
        self.window.request_redraw();
    }

    fn render(&mut self) -> Result<(), wgpu::SurfaceError> {
        // Grabbing an image for us to draw to.
        let drawable = self.surface.get_current_texture()?;
        let image_view_descriptor = wgpu::TextureViewDescriptor::default();
        let image_view = drawable.texture.create_view(&image_view_descriptor);

        let command_encoder_descriptor = wgpu::CommandEncoderDescriptor {
            label: Some("Render encoder"),
        };
        let mut command_encoder = self
            .device
            .create_command_encoder(&command_encoder_descriptor);

        let colour_attachment = wgpu::RenderPassColorAttachment {
            view: &image_view,
            resolve_target: None,
            ops: wgpu::Operations {
                load: wgpu::LoadOp::Clear(wgpu::Color {
                    r: 0.05,
                    g: 0.05,
                    b: 0.05,
                    a: 0.0,
                }),
                store: wgpu::StoreOp::Store,
            },
        };

        // Which resources will be modified by rendering
        let render_pass_descriptor = wgpu::RenderPassDescriptor {
            label: Some("Renderpass"),
            color_attachments: &[Some(colour_attachment)],
            depth_stencil_attachment: None,
            occlusion_query_set: None,
            timestamp_writes: None,
        };

        {
            let mut renderpass =
                command_encoder.begin_render_pass(&render_pass_descriptor);
            renderpass.set_pipeline(&self.render_pipeline);
            renderpass.set_bind_group(0, &self.uniform_bind_group, &[]);
            renderpass.set_vertex_buffer(0, self.triangle_mesh.slice(..));
            renderpass.draw(0..3, 0..1);
        }

        self.queue.submit(std::iter::once(command_encoder.finish()));

        drawable.present();

        Ok(())
    }
}

// fn get_aspect_ratio(width: u32, height: u32) -> f32 {
//     return width as f32 / height as f32;
// }
